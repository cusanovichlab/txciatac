setwd("C:\\Users\\hzhan\\Documents\\Dr.Cusanovich\\data\\sci_atac_drop\\test_blocking_conditions\\cell_bc_frequency")
library(tidyverse)
library(pheatmap)
library(plyr)

# check number of cells per well
filename = "Blocking_oligo.readcounts.readdepth.cells.indextable.txt"
base = "blocking_oligo\\Blocking_oligo_cells.per.well.pdf"

cell_seq = 10000
cell_well = 20000
well_num = 64
est_cell_num = cell_seq/(cell_well*well_num)*cell_well
break_max = est_cell_num

cellbc = read_tsv(filename, col_names = F)

table_transform = transform(cellbc, sample_bc = substr(X1, 1, 8), cell_bc = substr(X1, 9, 24), tn5_bc = substr(X1, 25, 32))
counts_tn5 = count(table_transform, "tn5_bc")

tagi5 = c('GAACCGCG', 'AGGTTATA', 'TCATCCTT', 'CTGCTTCC', 'GGTCACGA', 'AACTGTAG', 'GTGAATAT', 'ACAGGCGC', 'CATAGAGT', 'TGCGAGAC', 'GACGTCTT', 'AGTACTCC', 'TGGCCGGT', 'CAATTAAC', 'ATAATGTG', 'GCGGCACA', 'CTAGCGCT', 'TCGATATC', 'CGTCTGCG', 'TACTCATA', 'ACGCACCT', 'GTATGTTC', 'CGCTATGT', 'TATCGCAC', 'TCTGTTGG', 'CTCACCAA', 'TATTAGCT', 'CGCCGATC', 'TCTCTACT', 'CTCTCGTC', 'CCAAGTCT', 'TTGGACTC', 'GGCTTAAG', 'AATCCGGA', 'TAATACAG', 'CGGCGTGA', 'ATGTAAGT', 'GCACGGAC', 'GGTACCTT', 'AACGTTCC', 'GCAGAATT', 'ATGAGGCC', 'ACTAAGAT', 'GTCGGAGC', 'CCGCGGTT', 'TTATAACC', 'GGACTTGG', 'AAGTCCAA', 'ATCCACTG', 'GCTTGTCA', 'CAAGCTAG', 'TGGATCGA', 'AGTTCAGG', 'GACCTGAA', 'TGACGAAT', 'CAGTAGGC', 'AGCCTCAT', 'GATTCTGC', 'TCGTAGTG', 'CTACGACA', 'TAAGTGGT', 'CGGACAAC', 'ATATGGAT', 'GCGCAAGC', 'AAGATACT', 'GGAGCGTC', 'ATGGCATG', 'GCAATGCA', 'GTTCCAAT', 'ACCTTGGC', 'CTTATCGG', 'TCCGCTAA', 'GCTCATTG', 'ATCTGCCA', 'CTTGGTAT', 'TCCAACGC', 'CCGTGAAG', 'TTACAGGA', 'GGCATTCT', 'AATGCCTC', 'TACCGAGG', 'CGTTAGAA', 'CACGAGCG', 'TGTAGATA', 'GATCTATC', 'AGCTCGCT', 'CGGAACTG', 'TAAGGTCA', 'TTGCCTAG', 'CCATTCGA', 'ACACTAAG', 'GTGTCGGA', 'TTCCTGTT', 'CCTTCACC', 'GCCACAGG', 'ATTGTGAA')

count_tagi5 = rep(NA, 96)
names(count_tagi5) = tagi5

# counts_tn5_sum$tn5_bc[i] is factor, need to convert to character
for (i in 1:dim(counts_tn5)[1]) {
  tb5_bc = as.character(counts_tn5$tn5_bc[i])
  count_tagi5[tb5_bc] = counts_tn5$freq[i]
}

rowid = c("A","B","C","D","E","F","G","H")
colnid = c(1:12)
plate = matrix(count_tagi5, nrow = 8, ncol = 12, 
               byrow = T, dimnames = list(rowid, colnid))

maxnum = max(plate[!is.na(plate)])
mimnum = min(plate[!is.na(plate)])
mediannum = median(plate[!is.na(plate)])

pdf(base)
pheatmap(plate, na_col = "black", scale = "none",
         breaks = seq(0, break_max, length = 101),
         cluster_rows = F, cluster_cols = F, angle_col = 45,
         # border_color = "grey", 
         cellwidth = 30, cellheight = 30,
         show_rownames = T, show_colnames = T,
         gaps_row = c(1:7), gaps_col = c(1:11))
dev.off()

###################################################
###################################################
#
# check number of reads per well
#
###################################################
###################################################
filename = "DecoyDNA.readcounts.report.txt"
base = "DecoyDNA_reads.per.well.pdf"
complexity= rep(c(30220, 30717, 49420, 1, 25795, 30220), time = c(1, 1, 2, 4, 2, 2))

counts = read_tsv(filename, skip = 1, col_names = F)

table_transform = transform(counts, sample_bc = substr(X1, 1, 8), cell_bc = substr(X1, 9, 24), tn5_bc = substr(X1, 25, 32))
counts_tn5 = table_transform[,c(10,3)]

counts_tn5_sum = counts_tn5 %>%
  group_by(tn5_bc) %>%
  summarise(sum = sum(X3))

tagi5 = c('GAACCGCG', 'AGGTTATA', 'TCATCCTT', 'CTGCTTCC', 'GGTCACGA', 'AACTGTAG', 'GTGAATAT', 'ACAGGCGC', 'CATAGAGT', 'TGCGAGAC', 'GACGTCTT', 'AGTACTCC', 'TGGCCGGT', 'CAATTAAC', 'ATAATGTG', 'GCGGCACA', 'CTAGCGCT', 'TCGATATC', 'CGTCTGCG', 'TACTCATA', 'ACGCACCT', 'GTATGTTC', 'CGCTATGT', 'TATCGCAC', 'TCTGTTGG', 'CTCACCAA', 'TATTAGCT', 'CGCCGATC', 'TCTCTACT', 'CTCTCGTC', 'CCAAGTCT', 'TTGGACTC', 'GGCTTAAG', 'AATCCGGA', 'TAATACAG', 'CGGCGTGA', 'ATGTAAGT', 'GCACGGAC', 'GGTACCTT', 'AACGTTCC', 'GCAGAATT', 'ATGAGGCC', 'ACTAAGAT', 'GTCGGAGC', 'CCGCGGTT', 'TTATAACC', 'GGACTTGG', 'AAGTCCAA', 'ATCCACTG', 'GCTTGTCA', 'CAAGCTAG', 'TGGATCGA', 'AGTTCAGG', 'GACCTGAA', 'TGACGAAT', 'CAGTAGGC', 'AGCCTCAT', 'GATTCTGC', 'TCGTAGTG', 'CTACGACA', 'TAAGTGGT', 'CGGACAAC', 'ATATGGAT', 'GCGCAAGC', 'AAGATACT', 'GGAGCGTC', 'ATGGCATG', 'GCAATGCA', 'GTTCCAAT', 'ACCTTGGC', 'CTTATCGG', 'TCCGCTAA', 'GCTCATTG', 'ATCTGCCA', 'CTTGGTAT', 'TCCAACGC', 'CCGTGAAG', 'TTACAGGA', 'GGCATTCT', 'AATGCCTC', 'TACCGAGG', 'CGTTAGAA', 'CACGAGCG', 'TGTAGATA', 'GATCTATC', 'AGCTCGCT', 'CGGAACTG', 'TAAGGTCA', 'TTGCCTAG', 'CCATTCGA', 'ACACTAAG', 'GTGTCGGA', 'TTCCTGTT', 'CCTTCACC', 'GCCACAGG', 'ATTGTGAA')

count_tagi5 = rep(NA, 96)
names(count_tagi5) = tagi5

# counts_tn5_sum$tn5_bc[i] is factor, need to convert to character
for (i in 1:dim(counts_tn5_sum)[1]) {
  tb5_bc = as.character(counts_tn5_sum$tn5_bc[i])
  count_tagi5[tb5_bc] = counts_tn5_sum$sum[i]
}

rowid = c("A","B","C","D","E","F","G","H")
colnid = c(1:12)
plate = matrix(count_tagi5, nrow = 8, ncol = 12, 
                         byrow = T, dimnames = list(rowid, colnid))
plate_norm = t(t(plate)/complexity)

maxnum = max(plate_norm[!is.na(plate_norm)])
mimnum = min(plate_norm[!is.na(plate_norm)])
mediannum = median(plate_norm[!is.na(plate_norm)])

pdf(base)
pheatmap(plate_norm, na_col = "black", scale = "none",
         breaks = seq(0, 1.2*maxnum, length = 101),
         cluster_rows = F, cluster_cols = F, angle_col = 45,
         # border_color = "grey", 
         cellwidth = 30, cellheight = 30,
         show_rownames = T, show_colnames = T,
         gaps_row = c(1:7), gaps_col = c(1:11))
dev.off()