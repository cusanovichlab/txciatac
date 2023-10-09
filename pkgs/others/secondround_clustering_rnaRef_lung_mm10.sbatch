#!/bin/bash
#SBATCH --job-name=second_clustering_lung_mm10
#SBATCH --ntasks=72
#SBATCH --ntasks-per-node=72
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
OUTDIR=/xdisk/darrenc/haozhang1/20221117_scidrop_cc16ko_combined
base=combined_cc16_mm10
fragpath=$OUTDIR/fragments/${base}_fragments.tsv.gz
outfile=$OUTDIR/reports/rdata/$base/${base}_cluster_specific
dmat_path=$OUTDIR/reports/rdata/${base}_cluster.specific_c
meta_path=$OUTDIR/reports/rdata/$base/${base}_firstround.seurat.obj.filtered_singlet_harmony.rds
peak_path=$OUTDIR/reports/readcount_matrix/${base}_cluster.specific_readcounter.matrix.binary.rownames.txt
ref_path=/xdisk/darrenc/haozhang1/mouse_lung_snRNAseq/lung_snrna_seurat.obj.integrated.rds
hvar=Sample
group_var=Genotype
second_var=Input
pc=30
res=0.8
n=$(find $OUTDIR/reports/rdata/${base}_cluster.specific_c* |wc -l)

#############################
echo 'second round clustering...'
/home/u23/haozhang1/miniconda2/envs/seurat/bin/Rscript /groups/darrenc/sbin/scidropatac/secondround_clustering_rnaRef_mm.R \
--fragpath=$fragpath \
--outfile=$outfile \
--dmat_path=$dmat_path \
--meta_path=$meta_path \
--peak_path=$peak_path \
--ref_path=$ref_path \
--hvar=$hvar \
--group_var=$group_var \
--second_var=$second_var \
--pc=$pc \
--res=$res \
-n=$n

echo 'done'
date
