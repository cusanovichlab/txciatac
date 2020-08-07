import sys
import os
import pysam

#print len(sys.argv)
if len(sys.argv) != 6:
	sys.exit('Usage: python sc_atac_readcounter.py [Input Bam file] [Input Index table (NA if no table)] [Cell Type 1 BED] [Cell Type 2 BED] [Output file]')

inputbam = sys.argv[1]
indextable = sys.argv[2]
type1bed = sys.argv[3]
type2bed = sys.argv[4]
outfile = sys.argv[5]


type1desc = os.path.basename(type1bed)
type2desc = os.path.basename(type2bed)

totalct = {}
species1ct = {}
species2ct = {}
species1dhsct = {}
species2dhsct = {}
descriptor = {}
descdic = {}
if indextable != 'NA':
	descer = open(indextable,'r')
	for line in descer:
		liner = line.strip().split()
		descdic[liner[0]] = liner[1]

print "Counting total reads..."
bamfile = pysam.Samfile(inputbam,'rb')
for read in  bamfile.fetch():
	tagger = read.qname.split(':')[0]
	try:
		totalct[tagger] += 1
	except KeyError:
		totalct[tagger] = 1
		species1ct[tagger] = 0
		species2ct[tagger] = 0
		species1dhsct[tagger] = 0
		species2dhsct[tagger] = 0
		try:
			descriptor[tagger] = descdic[tagger]
		except KeyError:
			descriptor[tagger] = 'bkgd'
	if 'hg19' in bamfile.getrname(read.reference_id):
		species1ct[tagger] += 1
	if 'mm9' in bamfile.getrname(read.reference_id):
		species2ct[tagger] += 1

bamfile.close()

def lister(bedfile):
        currfile = open(bedfile,'r')
        currrecout = [line.strip().split()[0:3] for line in currfile]
        currfile.close()
        return currrecout

print "Counting celltype-specific reads..."
print "Building " + type1desc + " map..."
rec1list = lister(type1bed)
print "Building " + type2desc + " map..."
rec2list = lister(type2bed)
bamfile = pysam.Samfile(inputbam,'rb')

def counter(bedtuple,countdic,labeler):
	for rec in bedtuple:
		reads = bamfile.fetch(rec[0], int(rec[1]), int(rec[2]))
		for read in reads:
			readname = read.qname.split(':')[0]
			if 'CTF' in readname or 'AMBIG' in readname:
				continue
			countdic[readname] += 1

print "Counting " + type1desc + " reads..."
counter(rec1list,species1dhsct,type1desc)
print "Counting " + type2desc + " reads..."
counter(rec2list,species2dhsct,type2desc)

outter = open(outfile,'w')
print >> outter, 'Tag\tTotal\thg19\tmm9\t' + type1desc + '_DHS\t' + type2desc + '_DHS'
for tag in sorted(totalct.keys()):
	print >> outter, tag + '\t' + descriptor[tag] + '\t' + str(totalct[tag]) + "\t" + str(species1ct[tag]) + '\t' + str(species2ct[tag]) + '\t' + str(species1dhsct[tag]) + '\t' + str(species2dhsct[tag])

outter.close()
