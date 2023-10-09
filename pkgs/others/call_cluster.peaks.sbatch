#!/bin/bash
#SBATCH --job-name=cluster_peaks
#SBATCH --ntasks=72
#SBATCH --ntasks-per-node=72
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=5gb
#SBATCH --time=36:00:00
#SBATCH --partition=standard
#SBATCH --account=cusanovichlab
#SBATCH -o %x.out

#############################
##EDIT THESE AS NEEDED:
#############################

date
# module load R/3.6.1
module load samtools/1.10
#############################
##EDIT THESE AS NEEDED:
#############################
#
OUTDIR=OUTDIR
base=base
indextable=indextable
blacklist=blacklist
genome=genome
dedupbam=dedupbam
#############################
# isolate cells from each clusters

mkdir $OUTDIR/bams/clusters
mkdir $OUTDIR/macs/clusters

echo "Isolating cluster-specific cells..."
python /home/u23/haozhang1/scatac/process/sc_atac_library_deconvoluter.py $dedupbam $indextable $OUTDIR/bams/clusters/${base} .q10.sort.dedup.bam

# make a name list for cluster
cut -f2 ${indextable} | sort | uniq >$OUTDIR/${base}.cluster.name.txt

# call peaks on cluster specific beds
echo "Calling cluster-specific peaks using $genome genome..."

for cluster in `cat $OUTDIR/${base}.cluster.name.txt`
do
bedtools bamtobed -i $OUTDIR/bams/clusters/${base}.${cluster}.q10.sort.dedup.bam > $OUTDIR/beds/${base}.${cluster}.q10.sort.dedup.bed;
macs2 callpeak -t $OUTDIR/beds/${base}.${cluster}.q10.sort.dedup.bed -f BED -g $genome --nomodel --keep-dup all --extsize 200 --shift -100 --call-summits -n $OUTDIR/macs/clusters/${base}.${cluster}.nodups;
sort -k 8gr,8gr $OUTDIR/macs/clusters/${base}.${cluster}.nodups_peaks.narrowPeak | awk 'BEGIN{OFS="\t"}{$4="Peak_"NR ; print $0}' |\
sort -k1,1 -k2,2n -k3,3n | bedtools merge -i stdin > $OUTDIR/macs/clusters/${base}.${cluster}.nodups.narrowPeak_merged.bed
done

cat $OUTDIR/macs/clusters/${base}.*.nodups.narrowPeak_merged.bed \
| sort -k1,1 -k2,2n -k3,3n | bedtools merge -i stdin > $OUTDIR/macs/${base}_cluster.merged.peaks.bed

echo "removing ${blacklist} overlapping peaks..."
bedtools intersect -v -a $OUTDIR/macs/${base}_cluster.merged.peaks.bed -b $blacklist > $OUTDIR/peak_ref/${base}_cluster.whitelist.bed

echo "All peaks identified..."
cat $OUTDIR/macs/${base}_cluster.merged.peaks.bed |wc -l
echo "Peaks after removing blacklist..."
cat $OUTDIR/peak_ref/${base}_cluster.whitelist.bed |wc -l

cut -f1 $OUTDIR/peak_ref/${base}_cluster.whitelist.bed | uniq

# making peak count matrix using cluster specific peaks...
echo "generating readcount matrix..."

python /home/u23/haozhang1/scatac/preprocess/sc_atac_individual_dhs_readcounter.py \
$dedupbam \
$indextable \
$OUTDIR/peak_ref/${base}_cluster.whitelist.bed \
$OUTDIR/reports/readcount_matrix/${base}_cluster.specific_readcounter.matrix.txt \
True

echo "making sparse matrix for ${base}..."
# # arg1 = path to binary matrix
# # arg2 = output prefix 
# /home/u23/haozhang1/miniconda2/envs/r_env/bin/Rscript /groups/darrenc/sbin/scidropatac/mk_sparse_mat.nosplit.R \
# $OUTDIR/reports/readcount_matrix/${base}_cluster.specific_readcounter.matrix.binary.txt \
# $OUTDIR/reports/rdata/${base}_cluster
#############################
# split binary matrix
cut -f5-10000 $OUTDIR/reports/readcount_matrix/${base}_cluster.specific_readcounter.matrix.binary.txt >$OUTDIR/reports/readcount_matrix/${base}_cluster.specific_readcounter.matrix.binary_c1.txt
cut -f10001-20000 $OUTDIR/reports/readcount_matrix/${base}_cluster.specific_readcounter.matrix.binary.txt >$OUTDIR/reports/readcount_matrix/${base}_cluster.specific_readcounter.matrix.binary_c2.txt
cut -f20001-30000 $OUTDIR/reports/readcount_matrix/${base}_cluster.specific_readcounter.matrix.binary.txt >$OUTDIR/reports/readcount_matrix/${base}_cluster.specific_readcounter.matrix.binary_c3.txt
cut -f30001-40000 $OUTDIR/reports/readcount_matrix/${base}_cluster.specific_readcounter.matrix.binary.txt >$OUTDIR/reports/readcount_matrix/${base}_cluster.specific_readcounter.matrix.binary_c4.txt
cut -f40001-50000 $OUTDIR/reports/readcount_matrix/${base}_cluster.specific_readcounter.matrix.binary.txt >$OUTDIR/reports/readcount_matrix/${base}_cluster.specific_readcounter.matrix.binary_c5.txt
cut -f50001-60000 $OUTDIR/reports/readcount_matrix/${base}_cluster.specific_readcounter.matrix.binary.txt >$OUTDIR/reports/readcount_matrix/${base}_cluster.specific_readcounter.matrix.binary_c6.txt
cut -f60001- $OUTDIR/reports/readcount_matrix/${base}_cluster.specific_readcounter.matrix.binary.txt >$OUTDIR/reports/readcount_matrix/${base}_cluster.specific_readcounter.matrix.binary_c7.txt

# extract row and column from sparse matrix
cut -f4 $OUTDIR/reports/readcount_matrix/${base}_cluster.specific_readcounter.matrix.binary.txt | tail -n +2 >$OUTDIR/reports/readcount_matrix/${base}_cluster.specific_readcounter.matrix.binary.rownames.txt
#head -1 $OUTDIR/reports/readcount_matrix/lung_hg38_readcounter.matrix.binary.txt | cut -f5- >$OUTDIR/reports/readcount_matrix/lung_hg38_readcounter.matrix.binary.colnames.txt

mkdir $OUTDIR/scripts/mk_sparse_mat
cd $OUTDIR/scripts/mk_sparse_mat

echo "making sparse matrix scripts for ${base}..."
for file in $OUTDIR/reports/readcount_matrix/${base}_cluster.specific_readcounter.matrix.binary_c*.txt
do
sample=$(basename ${file} .txt)
out_prefix=$(echo $sample | sed 's/_readcounter.matrix.binary//')
sed "2 s/sparse_mat/${out_prefix}_sparse_mat/" /groups/darrenc/sbin/scidropatac/mk_sparse_mat_single.sbatch >$OUTDIR/scripts/mk_sparse_mat/${sample}_mk_dgcmat.sbatch
sed -i "17 s~OUTDIR~$OUTDIR~2" $OUTDIR/scripts/mk_sparse_mat/${sample}_mk_dgcmat.sbatch
sed -i "18 s~binary_matrix~$file~2" $OUTDIR/scripts/mk_sparse_mat/${sample}_mk_dgcmat.sbatch
sed -i "19 s/out_prefix/$out_prefix/2" $OUTDIR/scripts/mk_sparse_mat/${sample}_mk_dgcmat.sbatch
sbatch $OUTDIR/scripts/mk_sparse_mat/${sample}_mk_dgcmat.sbatch
done

echo 'done'
date
