from __future__ import print_function
import argparse
import os
import subprocess
import pysam

def submitter(commander):
        """Submits commands directly to the command line and waits for the process to finish."""
        submiting = subprocess.Popen(commander,shell=True)
        submiting.wait()

parser = argparse.ArgumentParser(description='A program to count duplicate reads for each cell for scATAC-seq analysis.')
parser.add_argument('-B','--inbam', help='Input bam file (not deduplicated)',dest='inbam',required=True)
parser.add_argument('-D','--indedup', help='Input deduplicated bam file',dest='indedup',required=True)
parser.add_argument('-O','--output',help='Output file name',required=True)
args = parser.parse_args()

print("Counting mapped reads...")
indexdic = {}
bamfile = pysam.Samfile(args.inbam,'rb')
for read in bamfile:
	cellname = read.qname.split(':')[0]
	if cellname not in indexdic:
		indexdic[cellname] = [0, 0]
	indexdic[cellname][0] += 1

print("Counting deduplicated reads...")

dedupfile = pysam.Samfile(args.indedup,'rb')
for read in dedupfile:
	cellname = read.qname.split(':')[0]
	if cellname not in indexdic:
		indexdic[cellname] = [0, 0]
	indexdic[cellname][1] += 1

print("Writing results for %s cells to %s..." % (len(indexdic), args.output))

outter = open(args.output, 'w')
outter.write("Barcode\tMapped\tDeduplicated\n")
for barcode in indexdic:
	if indexdic[barcode][0] > 0:
		outter.write(barcode + "\t" + str(indexdic[barcode][0]) + "\t" + str(indexdic[barcode][1]) + '\t' + '\n')

outter.close()
