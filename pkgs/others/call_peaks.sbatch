#!/bin/bash
#SBATCH --job-name=call_peaks
#SBATCH --ntasks=72
#SBATCH --ntasks-per-node=72
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=5gb
#SBATCH --time=48:00:00
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
# sample1 and 2 match the sample names in indextable.txt
OUTDIR=/xdisk/darrenc/haozhang1/20221117_scidrop_cc16ko_combined
base=combined_cc16
base1=lung100k
base2=lung200k
sample1=wt
sample2=ko
idx1=/xdisk/darrenc/haozhang1/20220909_scidrop_cc16ko/reports/${base1}_lung_dhs_peak.readcounts.readdepth.cells.indextable.txt
idx2=/xdisk/darrenc/haozhang1/20221114_scidrop_cc16ko_200k/reports/${base2}_lung_dhs_peak.readcounts.readdepth.cells.indextable.txt
dedup_bam1=/xdisk/darrenc/haozhang1/20220909_scidrop_cc16ko/bams/${base1}.q10.sort.dedup.bam
dedup_bam2=/xdisk/darrenc/haozhang1/20221114_scidrop_cc16ko_200k/bams/${base2}.q10.sort.dedup.bam
# barnyard_dedup_bam=$OUTDIR/bams/barnyard.q10.sort.dedup.bam
q10_bam1=/xdisk/darrenc/haozhang1/20220909_scidrop_cc16ko/bams/${base1}_merge.q10.sort.bam
q10_bam2=/xdisk/darrenc/haozhang1/20221114_scidrop_cc16ko_200k/bams/${base2}_merge.q10.sort.bam
# barnyard_q10_bam=$OUTDIR/bams/barnyard_merge.q10.sort.bam
#genome=hs
genome=mm
mm10_blacklist=/groups/darrenc/Hao/reference/blacklist/mm10-blacklist.v2.sorted.bed
# hg38_blacklist=/groups/darrenc/Hao/reference/blacklist/GRCh38_unified_blacklist_ENCFF356LFX.bed
#############################
#

mkdir $OUTDIR/bams
mkdir $OUTDIR/peak_ref
mkdir $OUTDIR/reports
mkdir $OUTDIR/beds
mkdir $OUTDIR/macs
mkdir $OUTDIR/scripts
mkdir $OUTDIR/bams/samples
mkdir $OUTDIR/reports/readcount_matrix
mkdir $OUTDIR/reports/rdata
mkdir $OUTDIR/fragments
mkdir $OUTDIR/scripts/mk_sparse_mat
# sparse matrix rds file will be saved in $OUTDIR/reports/rdata
mkdir $OUTDIR/reports/rdata

## call peaks
# sort -k1,1 -k2,2n -k3,3n /groups/darrenc/Hao/reference/blacklist/mm10-blacklist.v2.bed >/groups/darrenc/Hao/reference/blacklist/mm10-blacklist.v2.sorted.bed

echo "Isolating species specific barcodes from ${base1}..."
/home/u23/haozhang1/miniconda2/envs/r_env/bin/Rscript /groups/darrenc/sbin/scidropatac/isolate_sample_barcodes.R \
--output $OUTDIR/reports/${base1} \
--idx $idx1

echo "Isolating species specific barcodes from ${base2}..."
/home/u23/haozhang1/miniconda2/envs/r_env/bin/Rscript /groups/darrenc/sbin/scidropatac/isolate_sample_barcodes.R \
--output $OUTDIR/reports/${base2} \
--idx $idx2

echo 'merge index table...'
cat $OUTDIR/reports/${base1}_mm10.indextable.txt $OUTDIR/reports/${base2}_mm10.indextable.txt >$OUTDIR/reports/${base}_mm10.indextable.txt

echo 'Merging dedup bams...'
samtools merge $OUTDIR/bams/${base}.q10.sort.dedup.bam $dedup_bam1 $dedup_bam2

echo 'Isolate mouse reads from dedup bams...'
samtools view -h $OUTDIR/bams/${base}.q10.sort.dedup.bam | grep mm10 | sed "s/mm10//" | samtools view -Su - | samtools sort -@ 8 - -o $OUTDIR/bams/${base}_mm10.q10.sort.dedup.bam
samtools index $OUTDIR/bams/${base}_mm10.q10.sort.dedup.bam

# echo 'Isolate human reads from dedup bams...'
# samtools view -h $OUTDIR/bams/${base}.q10.sort.dedup.bam | grep hg38 | sed "s/hg38//" | samtools view -Su - | samtools sort -@ 8 - -o $OUTDIR/bams/${base}_hg38.q10.sort.dedup.bam
# samtools index $OUTDIR/bams/${base}_hg38.q10.sort.dedup.bam

echo 'demultiplex mouse dedup bams...'
python /home/u23/haozhang1/scatac/process/sc_atac_library_deconvoluter.py $OUTDIR/bams/${base}_mm10.q10.sort.dedup.bam $OUTDIR/reports/${base}_mm10.indextable.txt $OUTDIR/bams/samples/${base} _mm10.q10.sort.dedup.bam
samtools index $OUTDIR/bams/samples/${base}.${sample1}_mm10.q10.sort.dedup.bam
samtools index $OUTDIR/bams/samples/${base}.${sample2}_mm10.q10.sort.dedup.bam

# echo 'demultiplex human dedup bams...'
# python /home/u23/haozhang1/scatac/process/sc_atac_library_deconvoluter.py $OUTDIR/bams/${base}_hg38.q10.sort.dedup.bam $OUTDIR/reports/${base}.lung_hg38.indextable.txt $OUTDIR/bams/samples/${base} _hg38.q10.sort.dedup.bam
# samtools index $OUTDIR/bams/samples/${base}.${sample1}_hg38.q10.sort.dedup.bam

echo 'Call peaks for mouse samples...'
for file in $OUTDIR/bams/samples/${base}.*_mm10.q10.sort.dedup.bam
do
short=$(basename ${file} _mm10.q10.sort.dedup.bam)
bedtools bamtobed -i $OUTDIR/bams/samples/${short}_mm10.q10.sort.dedup.bam > $OUTDIR/beds/${short}_mm10.q10.sort.dedup.bed;
echo "Calling peaks on mouse ${short}......"
macs2 callpeak -t $OUTDIR/beds/${short}_mm10.q10.sort.dedup.bed -f BED -g $genome --nomodel --keep-dup all --extsize 200 --shift -100 --call-summits -n $OUTDIR/macs/${short}_mm10.nodups;
sort -k 8gr,8gr $OUTDIR/macs/${short}_mm10.nodups_peaks.narrowPeak | awk 'BEGIN{OFS="\t"}{$4="Peak_"NR ; print $0}' |\
sort -k1,1 -k2,2n -k3,3n | gzip -c > $OUTDIR/macs/${short}_mm10.nodups.narrowPeak.gz;
zcat $OUTDIR/macs/${short}_mm10.nodups.narrowPeak.gz | bedtools merge -i stdin > $OUTDIR/macs/${short}_mm10.nodups.narrowPeak_merged.bed
bedtools intersect -v -a $OUTDIR/macs/${short}_mm10.nodups.narrowPeak_merged.bed -b $mm10_blacklist > $OUTDIR/peak_ref/${short}_mm10_whitelist.bed
cut -f1 $OUTDIR/peak_ref/${short}_mm10_whitelist.bed | uniq
before=$(cat $OUTDIR/macs/${short}_mm10.nodups.narrowPeak_merged.bed |wc -l)
after=$(cat $OUTDIR/peak_ref/${short}_mm10_whitelist.bed |wc -l)
echo "For ${short}, identify ${before} peak; ${after} peaks in whitelist..."
done

# echo 'Call peaks for human samples...'
# for file in $OUTDIR/bams/samples/${base}.*_hg38.q10.sort.dedup.bam
# do
# short=$(basename ${file} _hg38.q10.sort.dedup.bam)
# bedtools bamtobed -i $OUTDIR/bams/samples/${short}_hg38.q10.sort.dedup.bam > $OUTDIR/beds/${short}_hg38.q10.sort.dedup.bed;
# echo "Calling peaks on human ${short}......"
# macs2 callpeak -t $OUTDIR/beds/${short}_hg38.q10.sort.dedup.bed -f BED -g hs --nomodel --keep-dup all --extsize 200 --shift -100 --call-summits -n $OUTDIR/macs/${short}_hg38.nodups;
# sort -k 8gr,8gr $OUTDIR/macs/${short}_hg38.nodups_peaks.narrowPeak | awk 'BEGIN{OFS="\t"}{$4="Peak_"NR ; print $0}' |\
# sort -k1,1 -k2,2n -k3,3n | gzip -c > $OUTDIR/macs/${short}_hg38.nodups.narrowPeak.gz;
# zcat $OUTDIR/macs/${short}_hg38.nodups.narrowPeak.gz | bedtools merge -i stdin > $OUTDIR/macs/${short}_hg38.nodups.narrowPeak_merged.bed
# bedtools intersect -v -a $OUTDIR/macs/${short}_hg38.nodups.narrowPeak_merged.bed -b $hg38_blacklist > $OUTDIR/peak_ref/${short}_hg38_whitelist.bed
# cut -f1 $OUTDIR/peak_ref/${short}_hg38_whitelist.bed | uniq
# before=$(cat $OUTDIR/macs/${short}_hg38.nodups.narrowPeak_merged.bed |wc -l)
# after=$(cat $OUTDIR/peak_ref/${short}_hg38_whitelist.bed |wc -l)
# echo "For ${short}, identify ${before} peak; ${after} peaks in whitelist..."
# done

# merge lung peaks called from different conditions
echo 'Merging peak whitelists...'
cat $OUTDIR/peak_ref/${base}.${sample1}_mm10_whitelist.bed $OUTDIR/peak_ref/${base}.${sample2}_mm10_whitelist.bed |\
sort -k1,1 -k2,2n -k3,3n | bedtools merge -i stdin >$OUTDIR/peak_ref/${base}_mm10.merged_whitelist.bed

allpeaks=$(cat $OUTDIR/peak_ref/${base}_mm10.merged_whitelist.bed |wc -l)
echo "Resulting in ${allpeaks} peaks after merging..."

# # merge cell-only bam files
# echo "merging cell-only bams..."
# samtools merge $OUTDIR/bams/samples/${base}.${merge_bam}.q10.sort.dedup.bam $OUTDIR/bams/samples/${base}.${sample1}.q10.sort.dedup.bam $OUTDIR/bams/samples/${base}.${sample2}.q10.sort.dedup.bam
# samtools index $OUTDIR/bams/samples/${base}.${merge_bam}.q10.sort.dedup.bam

# ################################
echo 'Merging q10.sort.bam...'
samtools merge $OUTDIR/bams/${base}.q10.sort.bam $q10_bam1 $q10_bam2

echo 'Isolate mouse reads from q10 bams...'
samtools view -h $OUTDIR/bams/${base}.q10.sort.bam | grep mm10 | sed "s/mm10//" | samtools view -Su - | samtools sort -@ 8 - -o $OUTDIR/bams/${base}_mm10.q10.sort.bam
samtools index $OUTDIR/bams/${base}_mm10.q10.sort.bam

# echo 'Isolate human reads from q10 bams...'
# samtools view -h $OUTDIR/bams/${base}.q10.sort.bam | grep hg38 | sed "s/hg38//" | samtools view -Su - | samtools sort -@ 8 - -o $OUTDIR/bams/${base}_hg38.q10.sort.bam
# samtools index $OUTDIR/bams/${base}_hg38.q10.sort.bam

# echo 'demultiplex mouse q10 bams...'
# python /home/u23/haozhang1/scatac/process/sc_atac_library_deconvoluter.py $OUTDIR/bams/${base}_mm10.q10.sort.bam $OUTDIR/reports/${base}_mm10.indextable.txt $OUTDIR/bams/samples/${base} _mm10.q10.sort.bam
# samtools index $OUTDIR/bams/samples/${base}.${sample1}_mm10.q10.sort.bam
# samtools index $OUTDIR/bams/samples/${base}.${sample2}_mm10.q10.sort.bam

# echo 'demultiplex human q10 bams...'
# python /home/u23/haozhang1/scatac/process/sc_atac_library_deconvoluter.py $OUTDIR/bams/${base}_hg38.q10.sort.bam $OUTDIR/reports/${base}.lung_hg38.indextable.txt $OUTDIR/bams/samples/${base} _hg38.q10.sort.bam
# samtools index $OUTDIR/bams/samples/${base}.${sample1}_hg38.q10.sort.bam

## Make read count matrix...
## args for mk_sparse_mat
# arg1 = path to binary matrix
# arg2 = output prefix

for file in $OUTDIR/bams/*_mm10.q10.sort.dedup.bam
do
short=$(basename ${file} .q10.sort.dedup.bam)
peak_whitelist=$OUTDIR/peak_ref/${short}.merged_whitelist.bed
echo "generating readcount matrix for ${short} using ${peak_whitelist}..."
python /home/u23/haozhang1/scatac/preprocess/sc_atac_individual_dhs_readcounter.py \
$OUTDIR/bams/${short}.q10.sort.dedup.bam \
$OUTDIR/reports/${short}.indextable.txt \
$peak_whitelist \
$OUTDIR/reports/readcount_matrix/${short}_readcounter.matrix.txt \
True;
echo "Make fragment files for ${short}..."
samtools view -h $OUTDIR/bams/${short}.q10.sort.bam \
| samtools sort -n -T $OUTDIR/bams/${short}.q10.sorted -@ 8 -\
| bedtools bamtobed -bedpe \
| awk '{print $1,$2-4,$6+5,$7}' | cut -f1 -d":" | sort -k1,1 -k2,2n -k3,3n -k4,4 | uniq -c \
| sed -e 's/ *//' -e 's/ /\t/' | awk -v OFS='\t' '{print $2,$3,$4,$5,$1}' \
| bgzip > $OUTDIR/fragments/${short}_fragments.tsv.gz;
tabix --preset=bed $OUTDIR/fragments/${short}_fragments.tsv.gz;
echo "Split readcount matrix..."
cut -f4 $OUTDIR/reports/readcount_matrix/${short}_readcounter.matrix.binary.txt | tail -n +2 >$OUTDIR/reports/readcount_matrix/${short}_readcounter.matrix.binary.rownames.txt
cut -f5-20000 $OUTDIR/reports/readcount_matrix/${short}_readcounter.matrix.binary.txt >$OUTDIR/reports/readcount_matrix/${short}_readcounter.matrix.binary_c1.txt
cut -f20001-40000 $OUTDIR/reports/readcount_matrix/${short}_readcounter.matrix.binary.txt >$OUTDIR/reports/readcount_matrix/${short}_readcounter.matrix.binary_c2.txt
cut -f40001-60000 $OUTDIR/reports/readcount_matrix/${short}_readcounter.matrix.binary.txt >$OUTDIR/reports/readcount_matrix/${short}_readcounter.matrix.binary_c3.txt
cut -f60001-80000 $OUTDIR/reports/readcount_matrix/${short}_readcounter.matrix.binary.txt >$OUTDIR/reports/readcount_matrix/${short}_readcounter.matrix.binary_c4.txt
cut -f80001-100000 $OUTDIR/reports/readcount_matrix/${short}_readcounter.matrix.binary.txt >$OUTDIR/reports/readcount_matrix/${short}_readcounter.matrix.binary_c5.txt
cut -f100001-120000 $OUTDIR/reports/readcount_matrix/${short}_readcounter.matrix.binary.txt >$OUTDIR/reports/readcount_matrix/${short}_readcounter.matrix.binary_c6.txt
cut -f120001-140000 $OUTDIR/reports/readcount_matrix/${short}_readcounter.matrix.binary.txt >$OUTDIR/reports/readcount_matrix/${short}_readcounter.matrix.binary_c7.txt
cut -f140001-160000 $OUTDIR/reports/readcount_matrix/${short}_readcounter.matrix.binary.txt >$OUTDIR/reports/readcount_matrix/${short}_readcounter.matrix.binary_c8.txt
cut -f160001-180000 $OUTDIR/reports/readcount_matrix/${short}_readcounter.matrix.binary.txt >$OUTDIR/reports/readcount_matrix/${short}_readcounter.matrix.binary_c9.txt
cut -f180001-200000 $OUTDIR/reports/readcount_matrix/${short}_readcounter.matrix.binary.txt >$OUTDIR/reports/readcount_matrix/${short}_readcounter.matrix.binary_c10.txt
cut -f200001- $OUTDIR/reports/readcount_matrix/${short}_readcounter.matrix.binary.txt >$OUTDIR/reports/readcount_matrix/${short}_readcounter.matrix.binary_c11.txt
done

cd $OUTDIR/scripts/mk_sparse_mat

for file in $OUTDIR/reports/readcount_matrix/*_readcounter.matrix.binary_c*.txt
do
sample=$(basename ${file} .txt)
out_prefix=$(echo $sample | sed 's/_readcounter.matrix.binary//')
echo "making sparse matrix scripts for ${sample}..."
sed "2 s/sparse_mat/${out_prefix}_sparse_mat/" /groups/darrenc/sbin/scidropatac/mk_sparse_mat_single.sbatch >$OUTDIR/scripts/mk_sparse_mat/${sample}_mk_dgcmat.sbatch
sed -i "17 s~OUTDIR~$OUTDIR~2" $OUTDIR/scripts/mk_sparse_mat/${sample}_mk_dgcmat.sbatch
sed -i "18 s~binary_matrix~$file~2" $OUTDIR/scripts/mk_sparse_mat/${sample}_mk_dgcmat.sbatch
sed -i "19 s/out_prefix/$out_prefix/2" $OUTDIR/scripts/mk_sparse_mat/${sample}_mk_dgcmat.sbatch
sbatch $OUTDIR/scripts/mk_sparse_mat/${sample}_mk_dgcmat.sbatch
done

echo 'done'
date
