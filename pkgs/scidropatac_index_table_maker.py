from __future__ import print_function

#Script to generate a table of indexes and sample names
#For each sample name, input the nextera i7 index numbers,
#pcr i7 index numbers, pcr i5 index numbers and nextera i5
#index numbers in the format lig_i7:pcr_i7:pcr_i5:lig_i5.
#To assign multiple indexes, separate numbers between colons
#with commas (e.g. '1,3,5').  To input a range, separate
#start and end with '-' (e.g. '1-6').
import sys
import argparse
#import barcode_to_well
#import barcode_constants as bc

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

cell_bc = []
with open("/groups/darrenc/sbin/sciATAC/737K-cratac-v1.txt", 'r') as bc_10x:
	for ele in bc_10x:
		cell_bc.append(ele.strip())

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Script to generate an index table for sci-ATAC runs.')
	parser.add_argument('indices', help='Index string. See overall pipeline docs for how these work.')
	parser.add_argument('name', help='A name for the sample who fit this index range.')
	parser.add_argument('--output', required=True, help='full path to output file')
	args = parser.parse_args()

	indices = args.indices
	name = args.name
	indexsplit = indices.strip().split(':')

	tagi7_indices = indexsplitter(indexsplit[0])
	pcri7_indices = indexsplitter(indexsplit[1])
	tagi5_indices = indexsplitter(indexsplit[2])

	tagi7 = ['TAAGGCGA', 'CGTACTAG', 'AGGCAGAA', 'TCCTGAGC', 'GGACTCCT', 'TAGGCATG', 'CTCTCTAC', 'CAGAGAGG', 'GCTACGCT', 'CGAGGCTG', 'AAGAGGCA', 'GTAGAGGA']
	pcri7 = cell_bc
	tagi5 = ['GAACCGCG', 'AGGTTATA', 'TCATCCTT', 'CTGCTTCC', 'GGTCACGA', 'AACTGTAG', 'GTGAATAT', 'ACAGGCGC', 'CATAGAGT', 'TGCGAGAC', 'GACGTCTT', 'AGTACTCC', 'TGGCCGGT', 'CAATTAAC', 'ATAATGTG', 'GCGGCACA', 'CTAGCGCT', 'TCGATATC', 'CGTCTGCG', 'TACTCATA', 'ACGCACCT', 'GTATGTTC', 'CGCTATGT', 'TATCGCAC', 'TCTGTTGG', 'CTCACCAA', 'TATTAGCT', 'CGCCGATC', 'TCTCTACT', 'CTCTCGTC', 'CCAAGTCT', 'TTGGACTC', 'GGCTTAAG', 'AATCCGGA', 'TAATACAG', 'CGGCGTGA', 'ATGTAAGT', 'GCACGGAC', 'GGTACCTT', 'AACGTTCC', 'GCAGAATT', 'ATGAGGCC', 'ACTAAGAT', 'GTCGGAGC', 'CCGCGGTT', 'TTATAACC', 'GGACTTGG', 'AAGTCCAA', 'ATCCACTG', 'GCTTGTCA', 'CAAGCTAG', 'TGGATCGA', 'AGTTCAGG', 'GACCTGAA', 'TGACGAAT', 'CAGTAGGC', 'AGCCTCAT', 'GATTCTGC', 'TCGTAGTG', 'CTACGACA', 'TAAGTGGT', 'CGGACAAC', 'ATATGGAT', 'GCGCAAGC', 'AAGATACT', 'GGAGCGTC', 'ATGGCATG', 'GCAATGCA', 'GTTCCAAT', 'ACCTTGGC', 'CTTATCGG', 'TCCGCTAA', 'GCTCATTG', 'ATCTGCCA', 'CTTGGTAT', 'TCCAACGC', 'CCGTGAAG', 'TTACAGGA', 'GGCATTCT', 'AATGCCTC', 'TACCGAGG', 'CGTTAGAA', 'CACGAGCG', 'TGTAGATA', 'GATCTATC', 'AGCTCGCT', 'CGGAACTG', 'TAAGGTCA', 'TTGCCTAG', 'CCATTCGA', 'ACACTAAG', 'GTGTCGGA', 'TTCCTGTT', 'CCTTCACC', 'GCCACAGG', 'ATTGTGAA']
	idx_tbl = open(args.output, 'w')
	for tagi7_id in tagi7_indices:
		for pcri7_id in pcri7_indices:
			for tagi5_id in tagi5_indices:
				try:
					tagmentation_i7_seq = tagi7[tagi7_id]
					pcr_i7_seq = pcri7[pcri7_id]
					tagmentation_i5_seq = tagi5[tagi5_id]
					barcodes_string = tagmentation_i7_seq + pcr_i7_seq + tagmentation_i5_seq
					print(barcodes_string + '\t' + name, file = idx_tbl)
				except IndexError:
					raise ValueError('One of your specified indices did not fall within the bounds of expected barcodes.')
	idx_tbl.close()
