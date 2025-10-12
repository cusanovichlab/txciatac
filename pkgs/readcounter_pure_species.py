# py3
import sys
import os
import pysam
import pybedtools
from pybedtools import BedTool

#print len(sys.argv)
if len(sys.argv) != 5:
	sys.exit('Usage: python sc_atac_readcounter.py [Input Bam file] [Input Index table (NA if no table)] [peak BED] [Output file]')

inputbam = sys.argv[1]
indextable = sys.argv[2]
type1bed = sys.argv[3]
#type2bed = sys.argv[4]
outfile = sys.argv[4]


type1desc = os.path.basename(type1bed)
#type2desc = os.path.basename(type2bed)

totalct = {}
#species1ct = {}
#species2ct = {}
species1dhsct = {}
#species2dhsct = {}
descriptor = {}
descdic = {}
if indextable != 'NA':
	with open(indextable, 'r') as descer:
		for line in descer:
			liner = line.strip().split()
			if len(liner) >= 2:
				descdic[liner[0]] = liner[1]

print("Counting total reads...")
bamfile = pysam.AlignmentFile(inputbam, 'rb')

# for read in  bamfile.fetch():
	# tagger = read.query_name.split(':')[0]
	# try:
		# totalct[tagger] += 1
	# except KeyError:
		# totalct[tagger] = 1
		# species1dhsct[tagger] = 0
		# try:
			# descriptor[tagger] = descdic[tagger]
		# except KeyError:
			# descriptor[tagger] = 'bkgd'

for read in bamfile.fetch():
    tagger = read.query_name.split(':')[0]
    totalct[tagger] = totalct.get(tagger, 0) + 1
    species1dhsct.setdefault(tagger, 0);
    if tagger not in descriptor:
        descriptor[tagger] = descdic.get(tagger, 'bkgd')

bamfile.close()

# def lister(bedfile):
#         currfile = open(bedfile,'r')
#         currrecout = [line.strip().split()[0:3] for line in currfile]
#         currfile.close()
#         return currrecout

# print "Counting celltype-specific reads..."
# print "Building " + type1desc + " map..."
# rec1list = lister(type1bed)
# ##print "Building " + type2desc + " map..."
# ##rec2list = lister(type2bed)
# bamfile = pysam.Samfile(inputbam,'rb')

bamfile = BedTool(inputbam)
print("Counting celltype-specific reads...")
print("Building " + type1desc + " map...")
rec1list = BedTool(type1bed)
inter1list = bamfile.intersect(rec1list, u=True, bed=True)

# def counter(bedtuple,countdic,labeler):
# 	for rec in bedtuple:
# 		reads = bamfile.fetch(rec[0], int(rec[1]), int(rec[2]))
# 		for read in reads:
# 			readname = read.qname.split(':')[0]
# 			if 'CTF' in readname or 'AMBIG' in readname:
# 				continue
# 			countdic[readname] += 1

def counter(bedinter,countdic,labeler):
	for read in bedinter:
		readname = read[3].split(':')[0]
		if 'CTF' in readname or 'AMBIG' in readname:
			continue
		countdic[readname] += 1

print("Counting " + type1desc + " reads...")
counter(inter1list,species1dhsct,type1desc)
#print "Counting " + type2desc + " reads..."
#counter(rec2list,species2dhsct,type2desc)

with open(outfile, 'w') as outter:
	outter.write('Tag\tTotal\t' + type1desc + '\n')
	for tag in sorted(totalct.keys()):
		outter.write(
			f"{tag}\t{descriptor[tag]}\t{totalct[tag]}\t{species1dhsct[tag]}\n"
		)
# outter = open(outfile,'w')
# print >> outter, 'Tag\tTotal\t' + type1desc
# for tag in sorted(totalct.keys()):
	# print >> outter, tag + '\t' + descriptor[tag] + '\t' + str(totalct[tag]) + '\t' + str(species1dhsct[tag])

# outter.close()
