#!/bin/bash
#SBATCH --job-name=make_index_table
#SBATCH --ntasks=7
#SBATCH --ntasks-per-node=7
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=5gb
#SBATCH --time=5:00:00
#SBATCH --partition=standard
#SBATCH --account=cusanovichlab
#SBATCH --array=1-N
#SBATCH -o %x_%a.out

echo 'running'

#############################
##EDIT THESE AS NEEDED:
#############################

date
# module load R/3.6.1

#############################
##EDIT THESE AS NEEDED:
#############################
OUTDIR=OUTDIR
path=path
declare -i x='SLURM_ARRAY_TASK_ID'

base=$(head -$x samplesheet |tail -1 | cut -f 1)
SAMPLE=$(head -$x samplesheet |tail -1 | cut -f 2)
INDICES=$(head -$x samplesheet |tail -1 | cut -f 3)

echo "base = $base, SAMPLE = $SAMPLE, INDICES = $INDICES"

#############################
mkdir $OUTDIR/reports
mkdir $OUTDIR/reports/indices_table

date
echo 'generating indices table...'
# make indice table-generate a table with all the combination of index used in your experiment
python $path/scidropatac_index_table_maker.py $INDICES $SAMPLE --output $OUTDIR/reports/indices_table/${base}_${SAMPLE}_indices_table.txt

echo "done indextable for $SAMPLE" >$OUTDIR/reports/indices_table/${base}_${SAMPLE}_mk_indextable_log.txt

echo 'done'
date
