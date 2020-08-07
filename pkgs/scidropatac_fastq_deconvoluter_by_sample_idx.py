# python 3
import sys

infastq1 = sys.argv[1]
infastq2 = sys.argv[2]
inbarcodes = sys.argv[3]
outprefix = sys.argv[4]

inindex = open(inbarcodes,'r')
#inindex[0] is barcode sequence, inindex[1] is sample name

libdic = {}
outdic1 = {}
outdic2 = {}
for line in inindex:
	liner = line.strip().split()
	libdic[liner[0]] = liner[1]
	try:
		outdic1[liner[1]]
	except KeyError:
		outdic1[liner[1]] = open(outprefix + "." + liner[1] + ".R1.fastq",'w')
		outdic2[liner[1]] = open(outprefix + "." + liner[1] + ".R2.fastq",'w')

outdic1['Unknown'] = open(outprefix + ".Unknown.R1.fastq",'w')
outdic2['Unknown'] = open(outprefix + ".Unknown.R2.fastq",'w')
inindex.close()

# with open(infastq1, 'r') as R1:
    # line1 = R1.readline().strip().split(':')
    # sequencer_id = line1[0]


readsin1 = open(infastq1,"r")
readsin2 = open(infastq2,"r")
barcodesout = open(outprefix + '.barcode_report.txt','w')
print("Library\tsci_index", file = barcodesout)

for line in readsin1:
	liner = line.strip().split()[1].split(':')[-1]
	matcher = liner[0:8]
	try:
		currlib = libdic[matcher]
	except KeyError:
		currlib = 'Unknown'
	print(currlib + '\t' + matcher, file = barcodesout)
	print(line.strip(), file = outdic1[currlib])
	print(readsin1.readline().strip(), file = outdic1[currlib])
	print(readsin1.readline().strip(), file = outdic1[currlib])
	print(readsin1.readline().strip(), file = outdic1[currlib])
	print(readsin2.readline().strip(), file = outdic2[currlib])
	print(readsin2.readline().strip(), file = outdic2[currlib])
	print(readsin2.readline().strip(), file = outdic2[currlib])
	print(readsin2.readline().strip(), file = outdic2[currlib])

barcodesout.close()

for fileh in outdic1.values():
	fileh.close()

for fileh in outdic2.values():
	fileh.close()

readsin1.close()
readsin2.close()
