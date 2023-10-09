#!/bin/bash
#SBATCH --job-name=first_clustering_lung_mm10
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
tss=/groups/darrenc/references/annotations/mm10/mm10.v23.tss.bed.gz
indextable=$OUTDIR/reports/rdata/$base/${base}_firstround.cluster_specific_cells.indextable.txt
blacklist=/groups/darrenc/Hao/reference/blacklist/mm10-blacklist.v2.sorted.bed
#hg38_blacklist=$OUTDIR/barnyard_hg38-blacklist.v3.bed
dedupbam=$OUTDIR/bams/${base}.q10.sort.dedup.bam
genome=mm
hvar=Sample
vln_group=Sample
cluster_var=Genotype
lane_var=Input
pc=30
res=0.8
estimated_doublet_rate=0.05

mkdir $OUTDIR/reports/rdata
mkdir $OUTDIR/reports/rdata/$base
mkdir $OUTDIR/scripts/cluster_peaks

#############################

echo 'First round clustering...'
# arg1 = path to binary matrix
# arg2 = output prefix
/home/u23/haozhang1/miniconda2/envs/seurat/bin/Rscript /groups/darrenc/sbin/scidropatac/firstround_clustering_mm_wt.R \
--fragpath=$OUTDIR/fragments/${base}_fragments.tsv.gz \
--outfile=$OUTDIR/reports/rdata/$base/${base}_firstround \
--dmat_path=$OUTDIR/reports/rdata/${base}.dgcmat.rds \
--meta_path=$OUTDIR/reports/${base}.qc_metadata.txt \
--tss_path=$tss \
--hvar=$hvar \
--vln_group=$vln_group \
--cluster_var=$cluster_var \
--lane_var=$lane_var \
--pc=$pc \
--res=$res \
--estimated_doublet_rate=$estimated_doublet_rate

# # use call_cluster.peaks.sbatch to call peaks on each cluster

# echo 'Call cluster specific peaks...'

# sed "2 s/cluster_peaks/${base}_cluster_peaks/" /groups/darrenc/sbin/scidropatac/call_cluster.peaks.sbatch > $OUTDIR/scripts/cluster_peaks/${base}_cluster_peaks.sbatch
# sed -i "23 s~OUTDIR~$OUTDIR~2" $OUTDIR/scripts/cluster_peaks/${base}_cluster_peaks.sbatch
# sed -i "24 s/base/$base/2" $OUTDIR/scripts/cluster_peaks/${base}_cluster_peaks.sbatch
# sed -i "25 s~indextable~$indextable~2" $OUTDIR/scripts/cluster_peaks/${base}_cluster_peaks.sbatch
# sed -i "26 s~blacklist~$blacklist~2" $OUTDIR/scripts/cluster_peaks/${base}_cluster_peaks.sbatch
# sed -i "27 s/genome/$genome/2" $OUTDIR/scripts/cluster_peaks/${base}_cluster_peaks.sbatch
# sed -i "28 s~dedupbam~$dedupbam~2" $OUTDIR/scripts/cluster_peaks/${base}_cluster_peaks.sbatch
# # # determine how many matrix needed to split based on total cell numbers
# # # and delete the lines exceeding total cell number and delete the last number, see below e.g.
# sed -i "87 s/50000//" $OUTDIR/scripts/cluster_peaks/${base}_cluster_peaks.sbatch
# sed -i '88,89d' $OUTDIR/scripts/cluster_peaks/${base}_cluster_peaks.sbatch

# cd $OUTDIR/scripts/cluster_peaks
# sbatch ${base}_cluster_peaks.sbatch

echo 'done'
date
