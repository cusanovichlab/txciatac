#!/bin/bash
#SBATCH --job-name=scidrop
#SBATCH --ntasks=28
#SBATCH --ntasks-per-node=28
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=5gb
#SBATCH --time=5:00:00
#SBATCH --partition=standard
#SBATCH --account=cusanovichlab
#SBATCH -o %x.out

#############################
##EDIT THESE AS NEEDED:
#############################

date
#module load R/3.6.1

#############################
##EDIT THESE AS NEEDED:
#############################
OUTDIR=OUTDIR
base=base
hs_peaks=hs_peaks
mm_peaks=mm_peaks
index_table=index_table
#bam_dir=bam_dir
short=short
#############################
#mkdir $OUTDIR/bams
#mkdir $OUTDIR/reports
#mkdir $OUTDIR/reports/indices_table

date
# generate reads count report, if not install pysam, use "conda config --add channels defaults"; "conda config --add channels conda-forge"; "conda config --add channels bioconda"; "conda install pysam"
# pysam requires "libffi.so.6" so add this line to .bashrc file: export LD_LIBRARY_PATH="/home/u23/haozhang1/miniconda2/envs/atac/lib"
echo 'generating readcounts report...'
python /groups/darrenc/sbin/scidropatac/readcounter.py $OUTDIR/bams/${base}.q10.sort.dedup.bam \
$OUTDIR/reports/indices_table/${index_table} \
$hs_peaks \
$mm_peaks \
$OUTDIR/reports/${short}.readcounts.report.txt

# call cell plots using R, need to install mclust
echo 'calling cells...'
/home/u23/haozhang1/miniconda2/envs/r_env/bin/Rscript /groups/darrenc/sbin/scidropatac/sc_atac_twospecies_cellcall_plots.R $OUTDIR/reports/${short}.readcounts 1000

echo 'plotting dedupreport...'
/home/u23/haozhang1/miniconda2/envs/r_env/bin/Rscript /groups/darrenc/sbin/scidropatac/sc_atac_dupreport_plotter.R \
$OUTDIR/reports/dedup/${base}.deduplicate_counter_results.txt \
$OUTDIR/reports/${short}.readcounts.readdepth.cells.indextable.txt \
$OUTDIR/reports/${short}

echo 'done'
date
