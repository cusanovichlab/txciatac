setwd("C:\\Users\\hzhan\\Documents\\Dr.Cusanovich\\script\\scATAC\\20200716_sci-drop-atac-100k\\samplesheets")
library(tidyverse)

indextb = read_tsv("make_indextable_samplesheet.txt", col_names = F)
# samplesheet = indextb[1:4, ]
samplesheet = separate(indextb, col = "X3", into = c("p7", "beads", "tn5"), sep = ":")
samplesheet$tn5 = gsub("-", ":", samplesheet$tn5, fixed = T)

# create list
# sp_list = vector(mode = "list", length = length(unique(samplesheet$X1)))
sp_vec = rep(NA, 96)
sp_list = rep(list(sp_vec), length(unique(samplesheet$X1)))
names(sp_list) = unique(samplesheet$X1)

# sp_list$P7_linear = samplesheet$X2
tn5bc = list()
# split by "," , strsplit output is a list
# when index the slot in list, it have to double quotes
for (i in 1:dim(samplesheet)[1]) {
  tn5bc[[i]] = unlist(strsplit(samplesheet$tn5[i], ","))
  for (j in tn5bc[[i]]) {
    if (grepl(":",  j, fixed = TRUE)) {
      num = unlist(strsplit( j, ":"))
      idx = c(as.integer(num[1]):as.integer(num[2]))
      samplename = samplesheet$X1[i]
      sp_list[[samplename]][idx] = samplesheet$X2[i]
    } else {
      idx = c(as.integer(j))
      sp_list[[samplename]][idx] = samplesheet$X2[i]
    }
  }
}

# make 96 well plate
rowid = c("A","B","C","D","E","F","G","H")
colnid = c(1:12)
P7_linear_plate = matrix(sp_list$P7_linear, nrow = 8, ncol = 12, 
                         byrow = T, dimnames = list(rowid, colnid))
P7_expo_plate = matrix(sp_list$P7_expo, nrow = 8, ncol = 12, 
                       byrow = T, dimnames = list(rowid, colnid))
DecoyDNA_plate = matrix(sp_list$DecoyDNA, nrow = 8, ncol = 12, 
                        byrow = T, dimnames = list(rowid, colnid))
Blocking_olig_plate = matrix(sp_list$Blocking_olig, nrow = 8, ncol = 12, 
                             byrow = T, dimnames = list(rowid, colnid))

