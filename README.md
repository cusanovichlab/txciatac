# scidropatac
## Generate FASTQ files
```
OUTDIR=/PATH/To/OUTDIR
bcl_dir=/PATH/To/BCL DIR

mkdir $OUTDIR/fastqs
mkdir $OUTDIR/fastqc

bcl2fastq --runfolder-dir $bcl_dir \
-o $OUTDIR/fastqs \
--sample-sheet $OUTDIR/SampleSheet.csv \
--ignore-missing-bcls \
--no-lane-splitting
```
## Demultiplexing FASTQ using P7 sample index
```
mkdir $OUTDIR/fastqs/temp

index=$OUTDIR/samplesheet/p7_8bp.index_demultiplex_table.txt
fastq1=$OUTDIR/fastqs/Undetermined_S0_R1_001.fastq.gz
fastq2=$OUTDIR/fastqs/Undetermined_S0_R2_001.fastq.gz
output_prefix=$OUTDIR/fastqs/temp/scidrop

python3.6 scidropatac_fastq_deconvoluter_by_sample_idx.py <(zcat $fastq1) <(zcat $fastq2) $index $output_prefix
```





