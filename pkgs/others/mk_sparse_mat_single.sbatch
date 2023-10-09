#!/bin/bash
#SBATCH --job-name=sparse_mat
#SBATCH --ntasks=28
#SBATCH --ntasks-per-node=28
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=5gb
#SBATCH --time=24:00:00
#SBATCH --partition=standard
#SBATCH --account=cusanovichlab
#SBATCH -o %x.out

### run your excutable program with begin and end date and time output
date
#############################
##EDIT THESE AS NEEDED:
#############################
OUTDIR=OUTDIR
binary_matrix=binary_matrix
out_prefix=out_prefix
#############################
mkdir $OUTDIR/reports
mkdir $OUTDIR/reports/rdata
echo 'making sparse matrix ${out_prefix}...'
# arg1 = path to binary matrix
# arg2 = output prefix
/home/u23/haozhang1/miniconda2/envs/r_env/bin/Rscript /groups/darrenc/sbin/scidropatac/mk_sparse_mat.R \
$binary_matrix \
$OUTDIR/reports/rdata/$out_prefix

echo 'done'
date
