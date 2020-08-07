# python 3
import sys

infastq1 = sys.argv[1]
infastq2 = sys.argv[2]
outprefix = sys.argv[3]


outdic1 = open(outprefix + ".R1.fastq",'w')
outdic2 = open(outprefix + ".R2.fastq",'w')

readsin1 = open(infastq1,"r")
readsin2 = open(infastq2,"r")

for line in readsin2:
	second_line = readsin2.readline().strip()
	tn5_bc = second_line[0:8]
	first_line = line.strip() + "+" + tn5_bc
	print(first_line, file = outdic2)
	print(second_line[27:], file = outdic2)
	print(readsin2.readline().strip(), file = outdic2)
	fourth_line = readsin2.readline().strip()
	print(fourth_line[27:], file = outdic2)
	print(readsin1.readline().strip() + "+" + tn5_bc, file = outdic1)
	print(readsin1.readline().strip(), file = outdic1)
	print(readsin1.readline().strip(), file = outdic1)
	print(readsin1.readline().strip(), file = outdic1)


outdic1.close()
outdic2.close()
readsin1.close()
readsin2.close()
