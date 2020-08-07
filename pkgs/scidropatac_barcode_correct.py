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

def get_barcode_seqs(r1_name, nextseq, two_level_indexed_tn5):
    """
    Extract the correct sequences from the R1 name.
    """
    # In 3LV runs, this is the the P5 + P7 index seq with a + in between
    # which is 20 + 1 + 20 (so have to skip "+" position)
    # Similar for two-level indexed Tn5, but Tn5 barcodes are 8bp
    if not two_level_indexed_tn5:
		barcodes = r1_name[-34:]
		#tagmentation_i7_seq = sample index
		tagmentation_i7_seq = barcodes[0:8]
		#pcr_i7_seq = 10x cell barcode
		pcr_i7_seq = barcodes[9:25]
		#pcr_i5_seq = barcodes[21:31]
		#tagmentation_i5_seq = tn5 barcodes
		tagmentation_i5_seq = barcodes[26:34]
    else:
        barcodes = r1_name[-37:]
        tagmentation_i7_seq = barcodes[0:8]
        pcr_i7_seq = barcodes[8:18]

        if nextseq:
            pcr_i5_seq = barcodes[27:37]
            tagmentation_i5_seq = barcodes[19:27]
        else:
            pcr_i5_seq = barcodes[19:29]
            tagmentation_i5_seq = barcodes[29:37]

    return tagmentation_i7_seq, pcr_i7_seq, tagmentation_i5_seq

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


def get_first_round_sample_lookup(samplesheet, tagi5):
    first_round_sample_lookup = {}
    for line in samplesheet:
        if line.startswith('sample_id\tranges'):
            continue

        entries = line.strip().split('\t')
        sample, indices = entries
        indexsplit = indices.split(':')

        # tagi7_indices = indexsplitter(indexsplit[0])
        tagi5_indices = indexsplitter(indexsplit[1])

        for i in tagi5_indices:
            first_round_sample_lookup[tagi5[i]] = sample

    return first_round_sample_lookup

cell_bc = []
with open("/xdisk/darrenc/mig2020/rsgrps/cusanovichlab/sbin/sciATAC/737K-cratac-v1.txt", 'r') as bc_10x:
	for ele in bc_10x:
		cell_bc.append(ele.strip())

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='A program to fix erroneous barcodes in scATAC data.')
    parser.add_argument('-1', '--input1', nargs='?', type=argparse.FileType('r'), default=sys.stdin, required=True, help='Text piped in from stdin for R1.')
    parser.add_argument('-2', '--input2', nargs='?', type=argparse.FileType('r'), default=sys.stdin, required=True, help='Text piped in from stdin for R2.')
    parser.add_argument('--samplesheet', required=True, help='Samplesheet describing the layout of the samples.')
    parser.add_argument('-o', '--output_prefix', help='Prefix to add to output files.')
    parser.add_argument('--stats_out', required=True, help='JSON file with output stats about processed reads and correction rates.')
    parser.add_argument('--two_level_indexed_tn5', action='store_true', help='Flag to run assuming that the library is a two-level indexed-TN5 sample.')
    parser.add_argument('--wells_384', action='store_true', help='Flag to run assuming That the known barcode set is the 384 well set.')
    parser.add_argument('--well_ids', action='store_true', help='Flag to output cell IDs that are composed of well IDs rather than the actual sequences.')
    parser.add_argument('-X', '--nextseq', help='NextSeq run indicator', dest='nextseq', action="store_true")
    args = parser.parse_args()

    if args.two_level_indexed_tn5 and args.wells_384:
        raise ValueError('There is no 384 well barcode set for indexed Tn5, may not specify both --two_level_indexed_tn5 and --wells_384.')

	# Set up the right index set depending on the indices
    if args.two_level_indexed_tn5:
        tagi7 = bc.nex_i7_two_level_indexed_tn5_list
        pcri7 = bc.pcr_i7_two_level_indexed_tn5_list
        pcri5 = bc.pcr_i5_two_level_indexed_tn5_list
        tagi5 = bc.nex_i5_two_level_indexed_tn5_list
    else:
        if args.wells_384:
            tagi7 = bc.lig_i7_list_384
            pcri7 = bc.pcr_i7_list_384
            pcri5 = bc.pcr_i5_list_384
            tagi5 = bc.lig_i5_list_384
            lig_i7_to_well = bc.lig_i7_to_well_384
            lig_i5_to_well = bc.lig_i5_to_well_384
            pcr_to_well = bc.pcr_to_well_384
        else:
            #tagi7 = bc.lig_i7_list
			#tagi7 = ['TAAGGCGA', 'CGTACTAG', 'AGGCAGAA', 'TCCTGAGC', 'GGACTCCT', 'TAGGCATG', 'CTCTCTAC', 'CAGAGAGG', 'GCTACGCT', 'CGAGGCTG', 'AAGAGGCA', 'GTAGAGGA']
            #pcri7 = bc.pcr_i7_list
			pcri7 = cell_bc
            #pcri5 = bc.pcr_i5_list
            #tagi5 = bc.lig_i5_list
			tagi5 = ['GAACCGCG', 'AGGTTATA', 'TCATCCTT', 'CTGCTTCC', 'GGTCACGA', 'AACTGTAG', 'GTGAATAT', 'ACAGGCGC', 'CATAGAGT', 'TGCGAGAC', 'GACGTCTT', 'AGTACTCC', 'TGGCCGGT', 'CAATTAAC', 'ATAATGTG', 'GCGGCACA', 'CTAGCGCT', 'TCGATATC', 'CGTCTGCG', 'TACTCATA', 'ACGCACCT', 'GTATGTTC', 'CGCTATGT', 'TATCGCAC', 'TCTGTTGG', 'CTCACCAA', 'TATTAGCT', 'CGCCGATC', 'TCTCTACT', 'CTCTCGTC', 'CCAAGTCT', 'TTGGACTC', 'GGCTTAAG', 'AATCCGGA', 'TAATACAG', 'CGGCGTGA', 'ATGTAAGT', 'GCACGGAC', 'GGTACCTT', 'AACGTTCC', 'GCAGAATT', 'ATGAGGCC', 'ACTAAGAT', 'GTCGGAGC', 'CCGCGGTT', 'TTATAACC', 'GGACTTGG', 'AAGTCCAA', 'ATCCACTG', 'GCTTGTCA', 'CAAGCTAG', 'TGGATCGA', 'AGTTCAGG', 'GACCTGAA', 'TGACGAAT', 'CAGTAGGC', 'AGCCTCAT', 'GATTCTGC', 'TCGTAGTG', 'CTACGACA', 'TAAGTGGT', 'CGGACAAC', 'ATATGGAT', 'GCGCAAGC', 'AAGATACT', 'GGAGCGTC', 'ATGGCATG', 'GCAATGCA', 'GTTCCAAT', 'ACCTTGGC', 'CTTATCGG', 'TCCGCTAA', 'GCTCATTG', 'ATCTGCCA', 'CTTGGTAT', 'TCCAACGC', 'CCGTGAAG', 'TTACAGGA', 'GGCATTCT', 'AATGCCTC', 'TACCGAGG', 'CGTTAGAA', 'CACGAGCG', 'TGTAGATA', 'GATCTATC', 'AGCTCGCT', 'CGGAACTG', 'TAAGGTCA', 'TTGCCTAG', 'CCATTCGA', 'ACACTAAG', 'GTGTCGGA', 'TTCCTGTT', 'CCTTCACC', 'GCCACAGG', 'ATTGTGAA']
            #lig_i7_to_well = bc.lig_i7_to_well
            #lig_i5_to_well = bc.lig_i5_to_well
            #pcr_to_well = bc.pcr_to_well

    # Build up sample mapping from first round index to sample
    first_round_sample_lookup = get_first_round_sample_lookup(open(args.samplesheet), tagi5)

    # TODO would eventually like to restrict these to the set that we use but no big deal for now
    # tagmentation_i7_whitelist = tagi7
    tagmentation_i5_whitelist = tagi5

    if not args.two_level_indexed_tn5:

        if args.nextseq:
            # p5_pcr_rc_map = {reverse_complement(k):k for k in pcri5}
            # p5_tagmentation_rc_map = {reverse_complement(k):k for k in tagi5}
			p7_pcr_rc_map = {reverse_complement(k):k for k in pcri7}
            # pcr_i5_whitelist = set([reverse_complement(x) for x in pcri5])
            # tagmentation_i5_whitelist = set([reverse_complement(x) for x in tagi5])
			pcr_i7_whitelist = set([reverse_complement(x) for x in pcri7])
        else:
            #pcr_i5_whitelist = pcri5
            pcr_i7_whitelist = pcri7

    else:
        tagmentation_i7_whitelist = bc.nex_i7_two_level_indexed_tn5
        pcr_i7_whitelist = bc.pcr_i7_two_level_indexed_tn5

        if args.nextseq:
            p5_pcr_rc_map = {reverse_complement(k):k for k in bc.pcr_i5_two_level_indexed_tn5}
            p5_tagmentation_rc_map = {reverse_complement(k):k for k in bc.nex_i5_two_level_indexed_tn5}

            pcr_i5_whitelist = set([reverse_complement(x) for x in bc.pcr_i5_two_level_indexed_tn5])
            tagmentation_i5_whitelist = set([reverse_complement(x) for x in bc.nex_i5_two_level_indexed_tn5])
        else:
            pcr_i5_whitelist = bc.pcr_i5_two_level_indexed_tn5
            tagmentation_i5_whitelist = bc.nex_i5_two_level_indexed_tn5

    # tagmentation_i7_correction_map = construct_mismatch_to_whitelist_map(tagmentation_i7_whitelist, 2)
    pcr_i7_correction_map = construct_mismatch_to_whitelist_map(pcr_i7_whitelist, 2)
    #pcr_i5_correction_map = construct_mismatch_to_whitelist_map(pcr_i5_whitelist, 2)
    tagmentation_i5_correction_map = construct_mismatch_to_whitelist_map(tagmentation_i5_whitelist, 2)

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
    #validreads['pcr_i5'] = 0
    validreads['pcr_i7'] = 0
    validreads['tagmentation_i5'] = 0
    #validreads['tagmentation_i7'] = 0
    validreads['all_barcodes'] = 0

    start = time.time()

    for (r1_name, r1_seq, r1_qual),(r2_name, r2_seq, r2_qual) in zip(if1, if2):

        totreads += 1

        # Get barcodes and correct
        tagmentation_i7_seq, pcr_i7_seq, tagmentation_i5_seq = get_barcode_seqs(r1_name, args.nextseq, args.two_level_indexed_tn5)
        #tagmentation_i7_seq = correct_barcode(tagmentation_i7_seq, tagmentation_i7_correction_map)
        pcr_i7_seq = correct_barcode(pcr_i7_seq, pcr_i7_correction_map)
        #pcr_i5_seq = correct_barcode(pcr_i5_seq, pcr_i5_correction_map)
        tagmentation_i5_seq = correct_barcode(tagmentation_i5_seq, tagmentation_i5_correction_map)

        # Skip invalid reads and track valid read count for error checking
        #if tagmentation_i7_seq is not None:
        #    validreads['tagmentation_i7'] += 1
        if tagmentation_i5_seq is not None:
            validreads['tagmentation_i5'] += 1
        if pcr_i7_seq is not None:
            validreads['pcr_i7'] += 1
        #if pcr_i5_seq is not None:
            #validreads['pcr_i5'] += 1
        
        if pcr_i7_seq is None or tagmentation_i5_seq is None:
            continue

        validreads['all_barcodes'] += 1

        # Map back to original whitelist if on nextseq so barcodes are always same on every sequencer
        if args.nextseq:
			pcr_i7_seq = p7_pcr_rc_map[pcr_i7_seq]

        # Convert to well IDs if requested
        # Note that for two level Tn5 barcode well comes first then PCR,
        # for three-level will be Tn5 N7, Tn5 N5, PCR WELL ID
        if args.two_level_indexed_tn5:
            barcodes_string = barcode_to_well.get_two_level_barcode_string(tagmentation_i7_seq, pcr_i7_seq, pcr_i5_seq, tagmentation_i5_seq, bc.nex_two_level_indexed_tn5_to_well, bc.pcr_two_level_indexed_tn5_to_well, args.well_ids)
            first_round_index = '%s%s' % (tagmentation_i5_seq, tagmentation_i7_seq)
        else:
			barcodes_string = tagmentation_i7_seq + pcr_i7_seq + tagmentation_i5_seq
			first_round_index = tagmentation_i5_seq
        try:
			sample = first_round_sample_lookup[first_round_index]
        except KeyError:
			sample = 'Unknown'
        output_files[sample]['r1'].write(''.join(['@', barcodes_string, ':', str(totreads), '#0000/1', '\n', r1_seq, '\n+\n', r1_qual, '\n']))
        output_files[sample]['r2'].write(''.join(['@', barcodes_string, ':', str(totreads), '#0000/2', '\n', r2_seq, '\n+\n', r2_qual, '\n']))

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
        subprocess.check_call('pigz %s' % output_files[sample]['r1_name'], shell=True)
        subprocess.check_call('pigz %s' % output_files[sample]['r2_name'], shell=True)
    print('Done compressing with pigz in %s minutes.' % ((time.time() - start) / 60.0))
