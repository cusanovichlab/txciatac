# scidropatac
### The scripts used in this repo are saved in 'pkgs', and samplesheets are saved in 'samplesheets'
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
## Append tn5 barcode to the index line in FASTQ and remove the first 27 bases which are ME (19bp) and tn5 bc (8bp) from Read2
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



