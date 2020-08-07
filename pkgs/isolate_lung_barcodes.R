# use plyr pkg to count charactor frequence in a column
# isolate lung barcode
library(tidyverse)
library(plyr)

setwd("C:\\Users\\hzhan\\Documents\\Dr.Cusanovich\\data\\sci_atac_drop\\test_blocking_conditions\\cell_bc_frequency")
# table = read_tsv("all_cells.indextable.txt", col_names = FALSE)
index_block = "Blocking_oligo.readcounts.readdepth.cells.indextable.txt"
index_decoy = "DecoyDNA.readcounts.readdepth.cells.indextable.txt"
index_p7linear = "P7_linear.readcounts.readdepth.cells.indextable.txt"
index_p7expo = "P7_expo.readcounts.readdepth.cells.indextable.txt"

# load cell index table
tb_block = read_tsv(index_block, col_names = FALSE)
tb_decoy = read_tsv(index_decoy, col_names = FALSE)
tb_p7linear = read_tsv(index_p7linear, col_names = FALSE)
tb_p7expo = read_tsv(index_p7expo, col_names = FALSE)

# merge p7 linear and p7 expo
tb_p7 = rbind(tb_p7linear, tb_p7expo)

# create list not used
# tb_list = list(block = tb_block, decoy = tb_decoy, 
#                p7linear = tb_p7linear, p7expo = tb_p7expo)

# isolate mm lung barcodes
block_lung = tb_block %>%
  filter(X2 %in% c("lung_mix", "mm_lung") & X3 == "Mouse")

decoy_lung = tb_decoy %>%
  filter(X2 %in% c("lung_mix", "mm_lung") & X3 == "Mouse")

p7_lung = tb_p7 %>%
  filter(X2 %in% c("lung_mix", "mm_lung") & X3 == "Mouse")

# replace sample name by species
block_lung = block_lung[, c(1, 3, 2)]
decoy_lung = decoy_lung[, c(1, 3, 2)]
p7_lung = p7_lung[, c(1, 3, 2)]

write_tsv(block_lung, "Blocking_oligo_mm_lung.indextable.txt", col_names = FALSE)
write_tsv(decoy_lung, "DecoyDNA_mm_lung.indextable.txt", col_names = FALSE)
write_tsv(p7_lung, "P7_mm_lung.indextable.txt", col_names = FALSE)

