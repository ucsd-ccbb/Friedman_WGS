#!/bin/bash
#$ -pe smp 32

workspace=/scratch

aws s3 cp s3://menieres-analysis-results/joint_genotype/2022-04-30_friedman_joint_metadata_04182022/ $workspace/ --recursive --exclude "*" --include "*norm.vep.maf"

for maf in $(ls $workspace/*/*/*maf)
    do out=$(echo $maf | sed 's/vep/vep.coding/'); awk -F '\t' '{ if($1 != "Unknown") { print }}' $maf > $out
done

RScript filterMaf.R
