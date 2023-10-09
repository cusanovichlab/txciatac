#!/usr/bin/env python
import argparse
import subprocess
import sys
import os
import gzip
import io
import itertools
import time
import json
import collections
import pickle
from Bio.SeqIO.QualityIO import FastqGeneralIterator
#import barcode_to_well
#import barcode_constants as bc


READ_BUFFER_SIZE = 1000000

def correct_barcode(barcode, mismatch_map):
    """
    Correct an observed raw barcode to one of a list of whitelists of mismatches.
    Args:
            barcode (string): barcode sequence to be corrected
            mismatch_map (list of dict dict): list of dict of mismatched sequences to real sequences
    Returns:
            string: corrected barcodes or None if barcode not correctable.
    """
    for mismatch_whitelist in mismatch_map:
        corrected = mismatch_whitelist.get(barcode, None)

        if corrected:
            return corrected

    return None


def generate_mismatches(sequence, num_mismatches, allow_n=True):
    """
    Generate a list of mimatched sequences to a given sequence. Must only contain ATGC.
    This is heavily based on a biostars answer.
    Args:
        sequence (str): The sequence must contain only A, T, G, and C
        num_mismatches (int): number of mismatches to generate sequences for
        allow_n (bool): True to allow N bases and False if not
    Yield:
    """
    letters = 'ACGT'

    if allow_n:
        letters += 'N'

    sequence = sequence.upper()
    mismatches = []

    for locs in itertools.combinations(range(len(sequence)), num_mismatches):
        sequence_list = [[char] for char in sequence]
        for loc in locs:
            orig_char = sequence[loc]
            sequence_list[loc] = [l for l in letters if l != orig_char]

        for poss in itertools.product(*sequence_list):
            mismatches.append(''.join(poss))

    return mismatches


def construct_mismatch_to_whitelist_map(whitelist, edit_distance, allow_n=True):
    """
    Constructs a precomputed set of all mimatches within a specified edit distance and the barcode whitelist.
    Args:
        whitelist (set of str): set of whitelist sequences
        edit_distance (int): max edit distance to consider
        allow_n (bool): True to allow N bases and False if not
    Returns:
        dict: mapping of mismatched sequences to their whitelist sequences
    """

    mismatch_to_whitelist_map = [None] * (edit_distance + 1)

    mismatch_to_whitelist_map[0] = {k: k for k in whitelist}

    conflicting_mismatches = []  # tracks conflicts where mismatches map to different sequences

    # Doesn't really matter as correction function will never see it,
    # but exclude any perfect matches to actual seqs by mismatches
    conflicting_mismatches.extend(list(whitelist))

    for mismatch_count in range(1, edit_distance + 1):
        mismatch_to_whitelist_map[mismatch_count] = {}

        for sequence in whitelist:
            sequence = sequence.upper()

            # Generate all possible mismatches in range
            mismatches = generate_mismatches(sequence, num_mismatches=mismatch_count, allow_n=allow_n)

            # Construct a mapping to the intended sequences
            for mismatch in mismatches:
                # Check for conflict with existing sequence and track if so
                if mismatch in mismatch_to_whitelist_map[mismatch_count]:
                    conflicting_mismatches.append(mismatch)
                mismatch_to_whitelist_map[mismatch_count][mismatch] = sequence

        # Go back and remove any conflicting mismatches
        for mismatch in set(conflicting_mismatches):
            if mismatch in mismatch_to_whitelist_map[mismatch_count]:
                del mismatch_to_whitelist_map[mismatch_count][mismatch]

    return mismatch_to_whitelist_map


def reverse_complement(x):
    complements = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    xrev = x[::-1]
    xrevcomp = ''.join([complements[z] for z in xrev])
    return xrevcomp

def get_barcode_seqs(r1_name):
    """
    Extract the correct sequences from the R1 name.
    """
    # The barcodes include 8 bp sample index, 16 bp 10x index and 8 bp tn5 index with a "+" in between
    barcodes = r1_name[-34:]
    sample_idx_seq = barcodes[0:8]
    gem_idx_seq = barcodes[9:25]
    tn5_idx_seq = barcodes[26:34]

    return sample_idx_seq, gem_idx_seq, tn5_idx_seq

def indexsplitter(indexrange):
	if len(indexrange) < 3:
		indexout = [int(indexrange)-1]
	elif "-" in indexrange or ',' in indexrange:
		range_list = [x for x in indexrange.split(",")]
		indexout = []
		for myrange in range_list:
			index_range = myrange.split('-')
			
			if len(index_range) == 1:
				start = int(index_range[0]) - 1
				end = start + 1
			elif len(index_range) == 2:
				start = int(index_range[0]) - 1
				end = int(index_range[1])
			else:
				raise ValueError('Invalid index range %s' % myrange)

			indexout.extend(range(start, end))
	else:
		raise ValueError('Invalid format for index range: %s' % indexrange)
	return indexout


def get_first_round_sample_lookup(samplesheet, sample_idx_whitelist, tn5_idx_whitelist):
    first_round_sample_lookup = {}
    for line in samplesheet:
        if line.startswith('sample_id\tranges'):
            continue

        entries = line.strip().split('\t')
        sample, indices = entries
        indexsplit = indices.split(':')

        sample_indices = indexsplitter(indexsplit[0])
        tn5_indices = indexsplitter(indexsplit[1])

        for i in sample_indices:
			for j in tn5_indices:
				first_round_sample_lookup[sample_idx_whitelist[i] + tn5_idx_whitelist[j]] = sample

    return first_round_sample_lookup

cell_bc = []
with open("/whitelist/737K-cratac-v1.txt", 'r') as bc_10x:
	for ele in bc_10x:
		cell_bc.append(ele.strip())

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='A program to fix erroneous barcodes in scATAC data.')
    parser.add_argument('-1', '--input1', nargs='?', type=argparse.FileType('r'), default=sys.stdin, required=True, help='Text piped in from stdin for R1.')
    parser.add_argument('-2', '--input2', nargs='?', type=argparse.FileType('r'), default=sys.stdin, required=True, help='Text piped in from stdin for R2.')
    parser.add_argument('--samplesheet', required=True, help='Samplesheet describing the layout of the samples.')
    parser.add_argument('-o', '--output_prefix', help='Prefix to add to output files.')
    parser.add_argument('--stats_out', required=True, help='JSON file with output stats about processed reads and correction rates.')
    parser.add_argument('-X', '--nextseq', help='NextSeq run indicator', dest='nextseq', action="store_true")
    args = parser.parse_args()

    # Set up the right index set depending on the indices
    sample_idx_whitelist = set(['TAAGGCGA', 'CGTACTAG', 'AGGCAGAA', 'TCCTGAGC', 'GGACTCCT', 'TAGGCATG', 'CTCTCTAC', 'CAGAGAGG', 'GCTACGCT', 'CGAGGCTG', 'AAGAGGCA', 'GTAGAGGA'])
    gem_idx_whitelist = set(cell_bc)
    tn5_idx_whitelist = set(['GAACCGCG', 'AGGTTATA', 'TCATCCTT', 'CTGCTTCC', 'GGTCACGA', 'AACTGTAG', 'GTGAATAT', 'ACAGGCGC', 'CATAGAGT', 'TGCGAGAC', 'GACGTCTT', 'AGTACTCC', 'TGGCCGGT', 'CAATTAAC', 'ATAATGTG', 'GCGGCACA', 'CTAGCGCT', 'TCGATATC', 'CGTCTGCG', 'TACTCATA', 'ACGCACCT', 'GTATGTTC', 'CGCTATGT', 'TATCGCAC', 'TCTGTTGG', 'CTCACCAA', 'TATTAGCT', 'CGCCGATC', 'TCTCTACT', 'CTCTCGTC', 'CCAAGTCT', 'TTGGACTC', 'GGCTTAAG', 'AATCCGGA', 'TAATACAG', 'CGGCGTGA', 'ATGTAAGT', 'GCACGGAC', 'GGTACCTT', 'AACGTTCC', 'GCAGAATT', 'ATGAGGCC', 'ACTAAGAT', 'GTCGGAGC', 'CCGCGGTT', 'TTATAACC', 'GGACTTGG', 'AAGTCCAA', 'ATCCACTG', 'GCTTGTCA', 'CAAGCTAG', 'TGGATCGA', 'AGTTCAGG', 'GACCTGAA', 'TGACGAAT', 'CAGTAGGC', 'AGCCTCAT', 'GATTCTGC', 'TCGTAGTG', 'CTACGACA', 'TAAGTGGT', 'CGGACAAC', 'ATATGGAT', 'GCGCAAGC', 'AAGATACT', 'GGAGCGTC', 'ATGGCATG', 'GCAATGCA', 'GTTCCAAT', 'ACCTTGGC', 'CTTATCGG', 'TCCGCTAA', 'GCTCATTG', 'ATCTGCCA', 'CTTGGTAT', 'TCCAACGC', 'CCGTGAAG', 'TTACAGGA', 'GGCATTCT', 'AATGCCTC', 'TACCGAGG', 'CGTTAGAA', 'CACGAGCG', 'TGTAGATA', 'GATCTATC', 'AGCTCGCT', 'CGGAACTG', 'TAAGGTCA', 'TTGCCTAG', 'CCATTCGA', 'ACACTAAG', 'GTGTCGGA', 'TTCCTGTT', 'CCTTCACC', 'GCCACAGG', 'ATTGTGAA'])
        
    # Build up sample mapping from first round index to sample
    first_round_sample_lookup = get_first_round_sample_lookup(open(args.samplesheet), sample_idx_whitelist, tn5_idx_whitelist)

    if args.nextseq:
        gem_whitelist_rc = {reverse_complement(k):k for k in gem_idx_whitelist}
        gem_idx_whitelist = set([reverse_complement(x) for x in gem_idx_whitelist])

    sample_correction_map = construct_mismatch_to_whitelist_map(sample_idx_whitelist, 2)
    gem_correction_map = construct_mismatch_to_whitelist_map(gem_idx_whitelist, 2)
    tn5_correction_map = construct_mismatch_to_whitelist_map(tn5_idx_whitelist, 2)

    # Set up all input/output files
    output_files = {}
    for sample in list(set(first_round_sample_lookup.values())):
        output_file_1 = '%s.%s_R1.fastq' % (args.output_prefix, sample)
        output_file_2 = '%s.%s_R2.fastq' % (args.output_prefix, sample)
        output_files[sample] = {}
        output_files[sample]['r1'] = open(output_file_1, 'w')
        output_files[sample]['r1_name'] = output_file_1
        output_files[sample]['r2'] = open(output_file_2, 'w')
        output_files[sample]['r2_name'] = output_file_2

    output_files['Unknown'] = {}
    output_files['Unknown']['r1'] = open(args.output_prefix + ".Unknown_R1.fastq", 'w')
    output_files['Unknown']['r1_name'] = args.output_prefix + ".Unknown_R1.fastq"
    output_files['Unknown']['r2'] = open(args.output_prefix + ".Unknown_R2.fastq", 'w')
    output_files['Unknown']['r2_name'] = args.output_prefix + ".Unknown_R2.fastq"

    if1 = FastqGeneralIterator(args.input1)
    if2 = FastqGeneralIterator(args.input2)

    totreads = 0
    validreads = {}
    validreads['gem_idx'] = 0
    validreads['tn5_idx'] = 0
    validreads['sample_idx'] = 0
    validreads['all_barcodes'] = 0

    start = time.time()

    for (r1_name, r1_seq, r1_qual),(r2_name, r2_seq, r2_qual) in zip(if1, if2):

        totreads += 1

        r1_rename = r1_name.split()[0]
        r2_rename = r2_name.split()[0]
	
        # Get barcodes and correct
        sample_idx_seq, gem_idx_seq, tn5_idx_seq = get_barcode_seqs(r1_name)
        sample_idx_seq = correct_barcode(sample_idx_seq, sample_correction_map)
        gem_idx_seq = correct_barcode(gem_idx_seq, gem_correction_map)
        #pcr_i5_seq = correct_barcode(pcr_i5_seq, pcr_i5_correction_map)
        tn5_idx_seq = correct_barcode(tn5_idx_seq, tn5_correction_map)

        # Skip invalid reads and track valid read count for error checking
        if sample_idx_seq is not None:
            validreads['sample_idx'] += 1
        if tn5_idx_seq is not None:
            validreads['tn5_idx'] += 1
        if gem_idx_seq is not None:
            validreads['gem_idx'] += 1
        
        if sample_idx_seq is None or gem_idx_seq is None or tn5_idx_seq is None:
            continue

        validreads['all_barcodes'] += 1

        # Map back to original whitelist if on nextseq so barcodes are always same on every sequencer
        if args.nextseq:
            gem_idx_seq = gem_whitelist_rc[gem_idx_seq]

        barcodes_string = sample_idx_seq + gem_idx_seq + tn5_idx_seq
        first_round_index = '%s%s' % (sample_idx_seq, tn5_idx_seq)
		
        try:
			sample = first_round_sample_lookup[first_round_index]
        except KeyError:
			sample = 'Unknown'
        output_files[sample]['r1'].write(''.join(['@', barcodes_string, ':', r1_rename, ' 1', '\n', r1_seq, '\n+\n', r1_qual, '\n']))
        output_files[sample]['r2'].write(''.join(['@', barcodes_string, ':', r2_rename, ' 2', '\n', r2_seq, '\n+\n', r2_qual, '\n']))

    if totreads == 0:
        raise ValueError('No reads found in fastq input.')

    # Output basic stats
    for stat in validreads:
        validreads[stat] = validreads[stat] / float(totreads)
    validreads['total_input_reads'] = totreads

    with open(args.stats_out, 'w') as f:
        f.write(json.dumps(validreads, f, indent=4))

    # Error checking and compress output
    if validreads['all_barcodes'] < 0.05:
        raise ValueError('Warning, you had less than 5 percent of all reads pass index correction. Something may have gone wrong here w.r.t. index sets or the expected library configuration not matching the data...')

    print('Done correcting barcodes in %s minutes. Starting compression...' % ((time.time() - start) / 60.0))
    start = time.time()
    for sample in output_files:
        output_files[sample]['r1'].close()
        output_files[sample]['r2'].close()
        subprocess.check_call('gzip %s' % output_files[sample]['r1_name'], shell=True)
        subprocess.check_call('gzip %s' % output_files[sample]['r2_name'], shell=True)
    print('Done compressing with gzip in %s minutes.' % ((time.time() - start) / 60.0))
