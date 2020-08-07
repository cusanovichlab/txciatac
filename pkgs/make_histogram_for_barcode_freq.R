# use plyr pkg to count charactor frequence in a column
library(tidyverse)
library(plyr)

setwd("C:\\Users\\hzhan\\Documents\\Dr.Cusanovich\\data\\sci_atac_drop\\test_blocking_conditions\\cell_bc_frequency")
# table = read_tsv("all_cells.indextable.txt", col_names = FALSE)
table = read_tsv("all.cells_scidrop.Blocking_oligo.readcounts.readdepth.cells.indextable.txt", col_names = FALSE)

# saparate cell barcodes [9:24] from tn5 barcodes [25,32], define columne name use "cell_bc" and "tn5_bc"
table_transform = transform(table, sample_bc = substr(X1, 1, 8), cell_bc = substr(X1, 9, 24), tn5_bc = substr(X1, 25, 32))

# count character frequence in column cell_bc
table_count = count(table_transform, "cell_bc")

pdf("cell_bc_freq.pdf", 5, 5)
ggplot(table_count, aes(x = freq)) +
  geom_histogram(binwidth = 1, color = "black", fill="dark grey") +
  theme_classic()+
  scale_x_continuous(breaks = seq(0, 39, by=5)) +
  scale_y_continuous(breaks = seq(0, 1700, by = 250)) +
  geom_hline(yintercept = 28, size = 1, linetype="dashed")
dev.off()

tn5_count = count(table_transform, "tn5_bc")
hist(tn5_count$freq, col="dark green")

collision_tb = table_transform %>%
	filter(X3 == "Collision")

# filter collision tn5 bc
tn5_collision = tn5_count %>%
	filter(tn5_bc %in% collision_tb$tn5_bc)
hist(tn5_collision$freq, col="dark blue")

# filter collision cell bc
cell_collision = table_count %>%
	filter(cell_bc %in% collision_tb$cell_bc)

# hist(cell_collision$freq, col="red", breaks = seq(1, 5, by = 1))
ggplot(cell_collision, aes(x = freq)) +
  geom_histogram(binwidth = 1, color = "black", fill="#FF6666") +
  theme_classic()+
  scale_x_continuous(breaks = seq(0, 39, by=5)) +
  scale_y_continuous(breaks = seq(0, 30, by = 5))

# frequency of collision in each number of cell barcode
single_cell.bc = table_count %>%
  filter(freq == 1)
singlet_tb = table_transform %>%
  filter(cell_bc %in% single_cell.bc$cell_bc) %>%
  mutate(Number_tn5 = "singlet")

duo_cell.bc = table_count %>%
  filter(freq == 2)
doublet_tb = table_transform %>%
  filter(cell_bc %in% duo_cell.bc$cell_bc) %>%
  mutate(Number_tn5 = "doublet")

tri_cell.bc = table_count %>%
  filter(freq == 3)
triplet_tb = table_transform %>%
  filter(cell_bc %in% tri_cell.bc$cell_bc) %>%
  mutate(Number_tn5 = "triplet")

quad_cell.bc = table_count %>%
  filter(freq == 4)
quadruplet_tb = table_transform %>%
  filter(cell_bc %in% quad_cell.bc$cell_bc) %>%
  mutate(Number_tn5 = "quadruplet")

quin_cell.bc = table_count %>%
  filter(freq == 5)
quintuple_tb = table_transform %>%
  filter(cell_bc %in% quin_cell.bc$cell_bc) %>%
  mutate(Number_tn5 = "quintuple")

collision = rbind(singlet_tb, doublet_tb, triplet_tb, quadruplet_tb, quintuple_tb)

collision_matrix = table(collision$Number_tn5, collision$X3)
collision_df = as.data.frame.matrix(t(collision_matrix))

collision_df = collision_df %>%
  rownames_to_column()

collision_df = collision_df %>%
  gather(key = "collision_type", value = "Number_of_cells", -rowname)

# play with collision dataframe to add percentage to plot
# Calculate the percentages
collision_df = ddply(collision_df, .(collision_type), transform, percent = Number_of_cells/sum(Number_of_cells) * 100)

# Format the labels and calculate their positions
# collision_df = ddply(collision_df, .(collision_type), transform, pos = (cumsum(Number_of_cells) - 0.5 * Number_of_cells))
collision_df$label = paste0(sprintf("%.2f", collision_df$percent), " %")

collision_df$collision_type = factor(collision_df$collision_type, 
                                     levels = c("singlet", "doublet", "triplet", "quadruplet", "quintuple"))

pdf("Collision_percet_in_GEM.pdf", 6, 7)
ggplot(collision_df, aes(collision_type, percent, fill = rowname)) +
  geom_bar(position= position_stack(), stat="identity", width = 0.6) +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 5) +
  theme(panel.background = element_blank(),
          axis.line = element_line(size = 1.5, colour = "black"),
        axis.ticks.length.y = unit(.15, "cm"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(face="bold", size=12),
        axis.text.y = element_text(face="bold", size=12), 
        axis.ticks = element_line(colour = "black", size = 1.5),
        axis.title.x = element_blank(),
        legend.key = element_rect(fill = "white"),
        legend.key.height = unit(1, 'cm'),
        legend.text = element_text(size = 12, face = "bold"))
dev.off()

# extract droplet contain no more than 3 droplet

tri_all_cell.bc = table_count %>%
  filter(freq <= 3)
tri_all_tb = table_transform %>%
  filter(cell_bc %in% tri_all_cell.bc$cell_bc) %>%
  select(X1, X2, X3)

write_tsv(tri_all_tb, "triplet.indextable.txt", col_names = FALSE)

######################################################
######################################################
# count the species number that can be separate by doublelet GEM
double_cell.bc = table_count %>%
	filter(freq == 2)
doublelet_sub = table_transform %>%
	filter(cell_bc %in% double_cell.bc$cell_bc)
write.csv(doublelet_sub, file = "cell.bc_with_multi_tn5.csv")

# cell_bc was saved as integer, so need to convert to string
doublelet_df = data.frame("type" = doublelet_sub$X3, cell_bc = as.character(doublelet_sub$cell_bc))
doublelet_matrix = table(doublelet_df$cell_bc, doublelet_df$type)

write.csv(doublelet_matrix, file = "doublelet_species_separation.csv")

doublelet_table = as.data.frame.matrix(doublelet_matrix)

doublelet_species_count = doublelet_table %>%
	mutate(Species = ifelse(Human == 2, "Human", ifelse(Mouse==2, "Mouse", ifelse(Collision==2, "Collision", ifelse(Collision==1, "ambig", "both")))))
write.csv(doublelet_species_count, file = "doublelet_species_count.csv")

species_count = count(doublelet_species_count, "Species")

barplot(species_count$freq, names.arg = c("Ambiguity", "Both", "Collision", "Human", "Mouse"), main = "Efficiency of Species separation in doublelet", ylab = "Cell barcode amount")

#########################################################
# pull out the cell barcode that only associated with one tn5 barcode
unique.bc = table_count %>%
	filter(freq == 1)

unique_table = table_transform %>%
	filter(cell_bc %in% unique.bc$cell_bc)

unique_output = unique_table %>%
	select(X1, X2, X3)

write_tsv(unique_output, "unique_tn5_beads.indexable.txt", col_names = FALSE)

################################################################
# plot cell recovery rate
s10k = c(10000, 7236, 5539)
s75k = c(75000, 53518, 19416)
name = c("input", "recovered cells", "unique tn5 cells")
barplot(c(s10k, s75k), names.arg = c(name, name), col = rep(c("lavender", "coral"), each=3), ylim = c(0,80000))

#################################################################
# isolate tn5 unique barcode and merge 10k and 75k
all_index = read_tsv("unique_tn5_beads.indexable.txt", col_names = FALSE)
lung_mix = all_index %>%
	filter(X2 == "Lung_mix")
lung_75k = lung_mix %>%
	filter(X3 != "Collision")
write_tsv(lung_mix_clear, "lung_mix_75k_unique_tn5.indextable.txt", col_names = FALSE)
unique_tn5_10k = read_tsv("C:\\Users\\hzhan\\Documents\\Dr.Cusanovich\\data\\sci_atac_drop_10k\\unique_tn5_beads.indexable.txt", col_names = FALSE)
lung_10k = unique_tn5_10k %>%
	filter(X2 == "lung_mix" & X3 != "Collision")

#replace "Lung_mix" with "lung_mix"
lung_75k$X2 = recode(lung_75k$X2, Lung_mix = "lung_mix")

# combine two dataframe vertically using rbind
lung_merged = rbind(lung_10k, lung_75k)

# combine two column into one
lung_merged = unite(lung_merged, sample.name, c(X2, X3), remove=FALSE)

write_tsv(lung_merged, "lung_merged_unique_tn5.indextable.txt", col_names = FALSE)

#############################################################
# pull out the cell barcode that associated with 2 tn5 barcode
double.tn5 = table_count %>%
  filter(freq == 2)

double_tn5tb = table_transform %>%
  filter(cell_bc %in% double.tn5$cell_bc)

# cell_bc was saved as integer, so need to convert to string
double_df = data.frame("type" = double_tn5tb$X3, cell_bc = as.character(double_tn5tb$cell_bc))
double_matrix = table(double_df$cell_bc, double_df$type)

double_table = as.data.frame.matrix(double_matrix)

double_species_count = double_table %>%
  mutate(Species = ifelse(Human == 2, "Human", ifelse(Mouse==2, "Mouse", ifelse(Collision==2, "Collision", ifelse(Collision==1, "ambig", "both")))))

double_species_count = double_species_count %>% add_rownames()

double_species_count = double_species_count %>%
  filter(Species != "Human" & Species != "Mouse")

double_species_tb = table_transform %>%
  filter(cell_bc %in% double_species_count$rowname)

double_species_tb = double_species_tb[,1:3]

write_tsv(double_species_tb, "double_tn5.indexable.txt", col_names = FALSE)
