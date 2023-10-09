#!/usr/bin/env python
import argparse
import subprocess
import sys
import gzip
import time
from Bio.SeqIO.QualityIO import FastqGeneralIterator

parser = argparse.ArgumentParser(description='A program to append barcode to fastq for scidrop data.')
parser.add_argument('-1', '--input1', nargs='?', type=argparse.FileType('r'), default=sys.stdin, required=True, help='Text piped in from stdin for R1.')
parser.add_argument('-2', '--input2', nargs='?', type=argparse.FileType('r'), default=sys.stdin, required=True, help='Text piped in from stdin for R2.')
parser.add_argument('-I1', '--idx1', nargs='?', type=argparse.FileType('r'), default=sys.stdin, required=True, help='Text piped in from stdin for R1 index.')
parser.add_argument('-I2', '--idx2', nargs='?', type=argparse.FileType('r'), default=sys.stdin, required=True, help='Text piped in from stdin for R2 index.')
parser.add_argument('-o', '--output_prefix', help='Prefix to add to output files.')
args = parser.parse_args()

# set up input and output files
output_files = {}
output_file_1 = '%s_R1.fastq' % (args.output_prefix)
output_file_2 = '%s_R2.fastq' % (args.output_prefix)
		
output_files['r1'] = open(output_file_1, 'w')
output_files['r2'] = open(output_file_2, 'w')

if1 = FastqGeneralIterator(args.input1)
if2 = FastqGeneralIterator(args.input2)
index1 = FastqGeneralIterator(args.idx1)
index2 = FastqGeneralIterator(args.idx2)

start = time.time()

for (r1_name, r1_seq, r1_qual),(r2_name, r2_seq, r2_qual),(i1_name, i1_seq, i1_qual),(i2_name, i2_seq, i2_qual) in zip(if1, if2, index1, index2):
	p7_bc = i1_seq
	cell_bc = i2_seq
	tn5_bc = r2_seq[0:8]
	bc = ''.join([p7_bc + "+" + cell_bc + "+" + tn5_bc])
	r2_seq_short = r2_seq[27:]
	r2_qual_short = r2_qual[27:]
	output_files['r1'].write(''.join(['@', r1_name, ':', bc, '\n', r1_seq, '\n+\n', r1_qual, '\n']))
	output_files['r2'].write(''.join(['@', r2_name, ':', bc, '\n', r2_seq_short, '\n+\n', r2_qual_short, '\n']))

print('Done adding barcodes in %s minutes. Starting compression...' % ((time.time() - start) / 60.0))
start = time.time()

output_files['r1'].close()
output_files['r2'].close()
subprocess.check_call('gzip %s' % output_file_1, shell=True)
subprocess.check_call('gzip %s' % output_file_2, shell=True)
print('Done compressing with gzip in %s minutes.' % ((time.time() - start) / 60.0))
