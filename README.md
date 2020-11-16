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
