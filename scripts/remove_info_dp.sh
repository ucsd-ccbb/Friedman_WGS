#!/bin/bash
#$ -cwd
#$ -pe smp 4

export PATH=/shared/workspace/software/bcbio/anaconda/bin:$PATH

workspace=/scratch/$sample
mkdir -p $workspace

echo "downloading "$sample"-joint-gatk-haplotype-annotated.vcf.gz"
aws s3 cp s3://menieres-analysis-results/$sample/ $workspace --recursive --exclude "*" --include "*gatk-haplotype-annotated.vcf.gz*"
mv $workspace/*/*gatk-haplotype-annotated.vcf.gz* $workspace

echo "removing INFO/DP"
bcftools annotate -x ^INFO/DP $workspace/"$sample"-joint-gatk-haplotype-annotated.vcf.gz -o $workspace/tmp.vcf.gz -Oz --threads 4
mv $workspace/tmp.vcf.gz $workspace/"$sample"-joint-gatk-haplotype-annotated-fixed.vcf.gz
echo "indexing $workspace/"$sample"-joint-gatk-haplotype-annotated-fixed.vcf.gz"
tabix -p vcf $workspace/"$sample"-joint-gatk-haplotype-annotated-fixed.vcf.gz

aws s3 cp $workspace/"$sample"-joint-gatk-haplotype-annotated-fixed.vcf.gz s3://menieres-analysis-results/$sample/$sample/
aws s3 cp $workspace/"$sample"-joint-gatk-haplotype-annotated-fixed.vcf.gz.tbi s3://menieres-analysis-results/$sample/$sample/
