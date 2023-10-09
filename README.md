# txci-ATAC-seq
This pipeline is used for processing txci-ATAC-seq data, including removing Tn5 barcode (8bp) from Read 2 and appending it to the header, correcting barcodes, removing sequence adapters, and aligning, filtering, and deduplicating reads, and generating species-specific read count report for mixed species assay.  
The scripts used in this repo are saved in the 'pkgs' folder, and example samplesheets are provided in the 'samplesheet_example' folder

## Packages required to install
### Linux packages
bcl2fastq, trimmomatic, bowtie2, samtools
### Python packages
argparse, pysam, pybedtools, Bio.SeqIO.QualityIO
### R packages
mclust
## Define the Output directory and the path to the pkgs folder and Rscript
```
OUTDIR=/path/to/outdir
path=/path/to/pkgs
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
OUTDIR=/PATH/To/Output/Directory
bcl_dir=/PATH/To/BCL/Directory

mkdir $OUTDIR/fastqs
mkdir $OUTDIR/fastqc

bcl2fastq --runfolder-dir $bcl_dir \
-o $OUTDIR/fastqs \
--sample-sheet $OUTDIR/bcl2fastq_SampleSheet.csv \
--ignore-missing-bcls \
--no-lane-splitting
```
## Demultiplexing FASTQ using P7 index
### pbs file: scidropatac_fastq_deconvoluter_by_sample_idx.py; use Python 3.6
```
mkdir $OUTDIR/fastqs/temp

index=$OUTDIR/samplesheet/p7_8bp.index_demultiplex_table.txt
fastq1=$OUTDIR/fastqs/Undetermined_S0_R1_001.fastq.gz
fastq2=$OUTDIR/fastqs/Undetermined_S0_R2_001.fastq.gz
output_prefix=$OUTDIR/fastqs/temp/scidrop

python3.6 scidropatac_fastq_deconvoluter_by_sample_idx.py <(zcat $fastq1) <(zcat $fastq2) $index $output_prefix
```
## Append tn5 barcode
### Tn5 barcodes are attached to Read2. This code will move the Tn5 bc from Read2 to the index line in FASTQ file.
### pbs file: scidropatac_add_tn5.bc.py; use Python 3.6
```
# mv Unknown fastqs to unknown folder
mkdir $OUTDIR/fastqs/temp/unknown
for file in $OUTDIR/fastqs/temp/*.Unknown.*.fastq
do
mv ${file} $OUTDIR/fastqs/temp/unknown
done

# append tn5 barcode
for file in $OUTDIR/fastqs/temp/*.R1.fastq
do
base=$(basename ${file} .R1.fastq)
fastq1=${base}.R1.fastq
fastq2=${base}.R2.fastq
outprefix=$OUTDIR/fastqs/temp/${base}_tn5.bc
python3.6 scidropatac_add_tn5.bc.py $OUTDIR/fastqs/temp/$fastq1 $OUTDIR/fastqs/temp/$fastq2 $outprefix
done
```
## Correct barcodes
### DO NOT demultiplex samples labeled by Tn5 barcodes in this step
### Need to check the samplesheet to make sure there is no quotation mark in it.
### Samplesheet requires beads barcode and tn5 barcode (row-wise), eg.
```
sample_id	ranges
bcfixed	1-737280:1-96
```
### Use scidropatac_barcode_correct.pbs as a template to generate the pbs files for each sample and run them parallelly
### Use pyhton 2.
```
mkdir $OUTDIR/scripts

for file in $OUTDIR/fastqs/temp/*_tn5.bc.R1.fastq
do
base=$(basename ${file} _tn5.bc.R1.fastq)
fastq1=${base}_tn5.bc.R1.fastq
fastq2=${base}_tn5.bc.R2.fastq
samplesheet=fixbc_samplesheet.txt
stats=${base}_bc.fixed.stats
sed "4 s/scidrop/${base}_fixbc/" scidropatac_barcode_correct.pbs > $OUTDIR/scripts/${base}_bc_correct.sh
sed -i "44 s~OUTDIR~$OUTDIR~2" $OUTDIR/scripts/${base}_bc_correct.sh
sed -i "45 s/samplesheet/$samplesheet/2" $OUTDIR/scripts/${base}_bc_correct.sh
sed -i "46 s/fastq1/$fastq1/2" $OUTDIR/scripts/${base}_bc_correct.sh
sed -i "47 s/fastq2/$fastq2/2" $OUTDIR/scripts/${base}_bc_correct.sh
sed -i "48 s/output/$base/2" $OUTDIR/scripts/${base}_bc_correct.sh
sed -i "50 s/stats/$stats/2" $OUTDIR/scripts/${base}_bc_correct.sh
qsub $OUTDIR/scripts/bc_correct_${base}.sh
done
```
### mv unknown fastq to unknown folder
```
mkdir $OUTDIR/fastqs/unknown
mv $OUTDIR/fastqs/*.Unknown_R*.fastq $OUTDIR/fastqs/unknown/
```
## Deduplicate
### pbs file: scidropatac_deduplicate.pbs; use Python 2.7.16 
```
for file in $OUTDIR/fastqs/*.bcfixed_R1.fastq
do
base=$(basename ${file} .bcfixed_R1.fastq)
fastq1=${base}.bcfixed_R1.fastq
fastq2=${base}.bcfixed_R2.fastq
sed "4 s/scidrop/${base}_duplicate/" scidropatac_deduplicate.pbs > $OUTDIR/scripts/${base}_deduplicate.sh
sed -i "42 s~OUTDIR~$OUTDIR~2" $OUTDIR/scripts/${base}_deduplicate.sh
sed -i "43 s/base/$base/2" $OUTDIR/scripts/${base}_deduplicate.sh
sed -i "44 s/fastq1/$fastq1/2" $OUTDIR/scripts/${base}_deduplicate.sh
sed -i "45 s/fastq2/$fastq2/2" $OUTDIR/scripts/${base}_deduplicate.sh
qsub $OUTDIR/scripts/${base}_deduplicate.sh
done
```
## Demultiplex samples labeled by Tn5 barcodes
### Create a samplesheet to assign the indices for each sample
### an example samplesheet is saved in 'samplesheet' named as make_indextable_samplesheet.txt
### column1: sample name labeled by P7 index; column2: sample name labeled by Tn5 barcodes; column3: indices
### the order of indices is: P7 barcode:10x beads barcodes:Tn5 barcodes(row-wise)
```
scidrop.SBS100	decoy_true	1:1-737280:1-8,13-20,25-32,37-44
scidrop.SBS100	decoy_pseudo	1:1-737280:9-12,21-24,33-36,45-48
scidrop.SBS100	ice_true	1:1-737280:49-56,61-68,73-80,85-92
scidrop.SBS100	ice_pseudo	1:1-737280:57-60,69-72,81-84,93-96
scidrop.SBS800	decoy_true	2:1-737280:1-8,13-20,25-32,37-44
scidrop.SBS800	decoy_pseudo	2:1-737280:9-12,21-24,33-36,45-48
scidrop.SBS800	ice_true	2:1-737280:49-56,61-68,73-80,85-92
scidrop.SBS800	ice_pseudo	2:1-737280:57-60,69-72,81-84,93-96
```
### Need to verify the Tn5 barcodes in the samplesheet using check_sample_well_id.R in 'pkg' folder
#### This R script will create a 8x12 matrix to resemble a 96-well plate. The sample name in each well should be consistent with the experimental setting.
### Make index table using scidropatac_make_index.table.pbs
#### Need to change the number of jobs in line 13, which is equal to the number of lines in make_indextable_samplesheet.txt
#### Need to change the output dir in line 26
#### After changing the above lines, submit the pbs file using qsub
```
qsub $OUTDIR/scripts/scidropatac_make_index.tables.pbs
```
### combine all index tables
```
find $OUTDIR/reports/indices_table/*.txt  -printf "%f\n" |\
sed 's/_/ /g' | awk '{print $1}'| sort | uniq >$OUTDIR/reports/indices_table/idx_base

for base in `cat $OUTDIR/reports/indices_table/idx_base`
do
cat $OUTDIR/reports/indices_table/${base}*indices_table.txt >$OUTDIR/reports/indices_table/${base}_allsample.indextable.txt
done

# remove individual index table
rm $OUTDIR/reports/indices_table/*indices_table.txt
```
## Create barnyard plots
### pbs file: scidropatac_readcounter.pbs; use Python 2
```
hs_peaks=$hs_lung_peaks
mm_peaks=$mm_lung_peaks

for file in $OUTDIR/reports/indices_table/*_allsample.indextable.txt
do
base=$(basename ${file} _allsample.indextable.txt)
sed "4 s/scidrop/${base}_readscount/" scidropatac_readcounter.pbs > $OUTDIR/scripts/${base}_readcounter.sh
sed -i "41 s~OUTDIR~$OUTDIR~2" $OUTDIR/scripts/${base}_readcounter.sh
sed -i "42 s/base/$base/2" $OUTDIR/scripts/${base}_readcounter.sh
sed -i "43 s~hs_peaks~$hs_peaks~2" $OUTDIR/scripts/${base}_readcounter.sh
sed -i "44 s~mm_peaks~$mm_peaks~2" $OUTDIR/scripts/${base}_readcounter.sh
sed -i "45 s/index_table/${base}_allsample.indextable.txt/2" $OUTDIR/scripts/${base}_readcounter.sh
qsub $OUTDIR/scripts/${base}_readcounter.sh
done
```
## For clustering, please see R codes in harmony_clustering_anno_lung.R in 'pkg' folder
