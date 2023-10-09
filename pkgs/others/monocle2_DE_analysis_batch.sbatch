#!/bin/bash
#SBATCH --job-name=monocle
#SBATCH --ntasks=72
#SBATCH --ntasks-per-node=72
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=5gb
#SBATCH --time=24:00:00
#SBATCH --partition=standard
#SBATCH --account=cusanovichlab
#SBATCH -o %x.out

#############################
##EDIT THESE AS NEEDED:
#############################
date
# module load R/3.6.1
# module load samtools/1.10

#############################
##EDIT THESE AS NEEDED:
#############################
OUTPUT=OUTPUT
cluster=cluster
seurat_obj=seurat_obj
group=group
ctl=ctl
treat=treat
min_cell=min_cell
#############################
echo 'Call Monocle...'
/home/u23/haozhang1/miniconda2/envs/seurat/bin/Rscript /groups/darrenc/sbin/scidropatac/monocle2_by_treatment_batch.R \
--monocle_patch /groups/darrenc/sbin/scidropatac/monocle_patch.R \
--seurat_obj $seurat_obj \
--cluster $cluster \
--group $group \
--ctl $ctl \
--treat $treat \
--min_cell $min_cell \
--cores 1 \
--output_directory $OUTPUT

echo 'done'
date
