# txci-ATAC-seq
This pipeline is used for processing txci-ATAC-seq data, including removing Tn5 barcode (8bp) from Read 2 and appending it to the header, correcting barcodes, removing sequence adapters, and aligning, filtering, and deduplicating reads, and generating species-specific read count report for mixed species assay.  

The scripts used in this repo are saved in the 'pkgs' folder, and example samplesheets are provided in the 'samplesheet_example' folder

## Packages required to install
### Linux packages
bcl2fastq, trimmomatic, bowtie2, samtools, pypy, macs2 and bedtools
### Python packages
argparse, pysam, pybedtools, Bio.SeqIO.QualityIO
### R packages
mclust tidyverse Matrix
## Define the output directory and the path for this repo and Rscript
```
OUTDIR=/path/to/outdir
path=/path/to/this/repo
Rscript=/path/to/Rscript
```
## Generate following folders
```
mkdir $OUTDIR/fastqs
mkdir $OUTDIR/fastqs/add_bc
mkdir $OUTDIR/fastqs/bc_correct
mkdir $OUTDIR/bams
mkdir $OUTDIR/samplesheet
mkdir $OUTDIR/reports
mkdir $OUTDIR/scripts
mkdir $OUTDIR/scripts/mk_idx_tb
```
## Generate FASTQ files
```
bcl_dir=/PATH/To/BCL/Directory
####################################

bcl2fastq --runfolder-dir $bcl_dir \
-o $OUTDIR/fastqs \
--sample-sheet $path/samplesheet_example/SampleSheet.csv \
--ignore-missing-bcls \
--create-fastq-for-index-reads \
--no-lane-splitting
```
## Format fastq files
Remove the first 27 bp (8 bp Tn5 barcode + 19 bp Tn5 mosaic-end adapter) from Read 2 and append the Tn5 barcode to fastq header
```
fastq1=$OUTDIR/fastqs/Undetermined_S0_R1_001.fastq.gz
fastq2=$OUTDIR/fastqs/Undetermined_S0_R2_001.fastq.gz
idx1=$OUTDIR/fastqs/Undetermined_S0_I1_001.fastq.gz
idx2=$OUTDIR/fastqs/Undetermined_S0_I2_001.fastq.gz
base=test
##########################

echo 'formating fastqs...'
python $path/pkgs/scidropatac_add_bc.py \
-1 <(zcat $fastq1) \
-2 <(zcat $fastq2) \
-I1 <(zcat $idx1) \
-I2 <(zcat $idx2) \
-o $OUTDIR/fastqs/add_bc/$base
```
## Correct barcode
In the sample sheet, the value/range before ":" is the index position (One-based indexing) of the P7 sample barcode, and the value/range after ":" is the index position (One-based indexing) of the Tn5 barcode (row-wisely ordered on a 96-well plate).

During barcode correction, samples can be demultiplexed using the combination of the P7 sample barcode and Tn5 barcode if the barcode indices of each sample are provided on separate lines within the sample sheet. However, as a standard practice, We do not perform sample demultiplexing at this step. Therefore, we specify all possible barcode indices for all samples on a single line.  

Use the "-X" option for data sequenced on the NextSeq platform.  

The output file name is in a format of "${base}.${sample id in samplesheet}_R1.fastq.gz".
```
samplesheet=$path/samplesheet_example/fixbc_samplesheet.txt
############################

echo 'correcting barcode...'
pypy $path/pkgs/scidropatac_barcode_correct.py \
--samplesheet $samplesheet \
-1 <(zcat $OUTDIR/fastqs/add_bc/${base}_R1.fastq.gz) -2 <(zcat $OUTDIR/fastqs/add_bc/${base}_R2.fastq.gz) \
-o $OUTDIR/fastqs/bc_correct/$base \
--stats_out $OUTDIR/reports/${base}_bc.fixed.stats -X
```
## Read alignment, filtering, and deduplication
```
base=test.bc_correct
fastq1=$OUTDIR/fastqs/bc_correct/${base}_R1.fastq.gz
fastq2=$OUTDIR/fastqs/bc_correct/${base}_R2.fastq.gz
genome_dir=/path/to/bowtie/hg38_mm10/hg38_mm10
adapter_path=$path/Trimmomatic-0.36/adapters
#############################

echo 'trimming adaptors...'
java -Xmx1G -jar $path/Trimmomatic-0.36/trimmomatic-0.36.jar \
PE \
-threads 24 \
$fastq1 $fastq2 \
$OUTDIR/fastqs/${base}_trimmed_R1_paired.fq.gz $OUTDIR/fastqs/${base}_trimmed_R1_unpaired.fq.gz \
$OUTDIR/fastqs/${base}_trimmed_R2_paired.fq.gz $OUTDIR/fastqs/${base}_trimmed_R2_unpaired.fq.gz \
ILLUMINACLIP:$adapter_path/iTSM-PE.fa:2:30:10:1:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:10 MINLEN:20 \
2> $OUTDIR/reports/${base}_trimmomatic.log;

echo 'mapping reads...'
bowtie2 -p 8 -X 2000 -3 1 -x $genome_dir \
-1 $OUTDIR/fastqs/${base}_trimmed_R1_paired.fq.gz \
-2 $OUTDIR/fastqs/${base}_trimmed_R2_paired.fq.gz \
2> $OUTDIR/reports/${base}.bowtie2.log | samtools view -bS - > $OUTDIR/bams/${base}.bam

echo 'filtering reads...'
samtools view -h -f3 -F12 -q10 $OUTDIR/bams/${base}.bam \
| grep -v '[0-9]'$'\t'hg38chrM | grep -v '[0-9]'$'\t'mm10chrM \
| grep -v '[0-9]'$'\t'GL | grep -v '[0-9]'$'\t'KI | grep -v '[0-9]'$'\t'JH \
| samtools view -Su - | samtools sort -@ 16 -T $OUTDIR/bams/${base}.sorttemp \
-o $OUTDIR/bams/${base}.q10.sort.bam

echo 'deduplicating...'
python $path/pkgs/sc_atac_true_dedup.py $OUTDIR/bams/${base}.q10.sort.bam $OUTDIR/bams/${base}.q10.sort.dedup.bam

echo 'indexing dedup.bam...'
samtools index $OUTDIR/bams/${base}.q10.sort.dedup.bam

echo 'generating deduplicate counter report...'
python $path/pkgs/duplicate_counter.py \
-B $OUTDIR/bams/${base}.q10.sort.bam \
-D $OUTDIR/bams/${base}.q10.sort.dedup.bam \
-O $OUTDIR/reports/${base}.deduplicate_counter_results.txt
```
## Create index table 
Generate all possible barcode combinations for each sample that were pooled together.  

The resulting index table will be used to generate the read count report. This allows us to calculate the collision rate for each sample separately.  

In the sample sheet, the first column is the ${base} variable used in the above section, the second column is the sample name, and the third column is the index position for each barcode.  

For barcode indices, the value/range before the first ":" is the index position () of the P7 sample barcode, the value/range between the first ":" and second ":" is the index position () of the 10x GEM barcode (this is usually a fixed range from 1-737280), and the last value/range is the index position of the Tn5 barcode (row-wisely ordered on a 96-well plate).

We submit an array job using SLURM to make an index table for each sample in parallel, and then merge all index tables.
```
# Variable N is the number of samples that will be generated after demultiplexing (equal to the number of lines in the demultiplexing samplesheet)
N=4
samplesheet=$path/samplesheet_example/make_indextable_samplesheet.txt
#####################################################################

echo 'Making index table...'
sed "2 s/make_index_table/${base}_idxtb/" $path/pkgs/scidropatac_make_index.tables.sbatch > $OUTDIR/scripts/mk_idx_tb/scidropatac_make_index_${base}.sbatch
sed -i "10 s/N/$N/" $OUTDIR/scripts/mk_idx_tb/scidropatac_make_index_${base}.sbatch
sed -i "25 s~OUTDIR~$OUTDIR~2" $OUTDIR/scripts/mk_idx_tb/scidropatac_make_index_${base}.sbatch
sed -i "26 s~path~$path~2" $OUTDIR/scripts/mk_idx_tb/scidropatac_make_index_${base}.sbatch
sed -i "29 s~samplesheet~$samplesheet~" $OUTDIR/scripts/mk_idx_tb/scidropatac_make_index_${base}.sbatch
sed -i "30 s~samplesheet~$samplesheet~" $OUTDIR/scripts/mk_idx_tb/scidropatac_make_index_${base}.sbatch
sed -i "31 s~samplesheet~$samplesheet~" $OUTDIR/scripts/mk_idx_tb/scidropatac_make_index_${base}.sbatch

# submit jobs
cd $OUTDIR/scripts/mk_idx_tb/
sbatch scidropatac_make_index_${base}.sbatch

date
echo "waiting for making indextable......"
sleep 5m
echo "counting index tables generated......"
date

x=$(find $OUTDIR/reports/indices_table/*_mk_indextable_log.txt |wc -l)
while [ $x -lt $N ]
do
echo "generated indextable for $x samples"
sleep 15s
x=$(find $OUTDIR/reports/indices_table/*_mk_indextable_log.txt |wc -l)
done

while [ $x == $N ]
do
echo "generated indextable for all samples"
x=$(( $x + 1 ))
done

echo 'combining index table...'
cat $OUTDIR/reports/indices_table/*_indices_table.txt \
>$OUTDIR/reports/${base}.indextable.txt
```
## Generate read count report
Quantify the species-specific reads for each cell barcode, and calculate the Fraction of Reads in Peaks (FRiP). 
```
index_table=$OUTDIR/reports/${base}.indextable.txt
hs_peaks=/path/to/human/peak.bed
mm_peaks=/path/to/mouse/peak.bed
######################################

echo 'generating readcounts report...'
python $path/pkgs/readcounter.py $OUTDIR/bams/${base}.q10.sort.dedup.bam \
$index_table \
$hs_peaks \
$mm_peaks \
$OUTDIR/reports/${base}.readcounts.report.txt

echo 'calling cells passing a given threshold...'
# In the second argument, either specify a number or "mclust" to calculate the cell threshold automatically
$Rscript $path/pkgs/sc_atac_twospecies_cellcall_plots.R $OUTDIR/reports/${base}.readcounts 1000

echo 'plotting deduplication report...'
$Rscript $path/pkgs/sc_atac_dupreport_plotter.R \
$OUTDIR/reports/${base}.deduplicate_counter_results.txt \
$OUTDIR/reports/${base}.readcounts.readdepth.cells.indextable.txt \
$OUTDIR/reports/${base}
```
