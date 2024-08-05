#!/bin/bash

#SBATCH --array=1-2%2 --cpus-per-task=20 --mem=100G --time=10-00:00:00

sampleArray=("BM_rep2" "SP_rep2")
sample=${sampleArray[${SLURM_ARRAY_TASK_ID}-1]}
fastqDir=/scratch/twlab/fanc/spbm2/fastq/

cellranger-arc count --id $sample --reference /scratch/twlab/fanc/software/cellrangerArc_2.0/refdata-cellranger-arc-mm10-2020-A-2.0.0 \
--localcores 18 --localmem 90 \
--gex-exclude-introns \
--libraries ${fastqDir}/${sample}/${sample}.csv 1>count.log.${sample} 2>&1