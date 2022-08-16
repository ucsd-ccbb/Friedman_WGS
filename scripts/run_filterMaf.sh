#!/bin/bash
#$ -cwd
#$ -pe smp 4

export PATH=/shared/workspace/software/bcbio/anaconda/bin:$PATH

workspace=/scratch/maf/$sample
mkdir -p $workspace

aws s3 cp s3://menieres-analysis-results/joint_genotype/$sample/vep/"$sample".vep.maf $workspace/

Rscript /shared/workspace/projects/friedman/scripts/filterMaf.R $workspace/"$sample".vep.maf

aws s3 cp $workspace/"$sample".filt.maf s3://menieres-analysis-results/joint_genotype/$sample/vep/
