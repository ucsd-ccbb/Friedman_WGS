#!/bin/bash
#$ -pe smp 32

export PATH=/shared/workspace/software/bcbio/anaconda/envs/bcftools1.15/bin:/shared/workspace/software/bcbio/anaconda/bin:$PATH
bcftools=/shared/workspace/software/bcbio/anaconda/envs/bcftools1.15/bin/bcftools
s3=s3://menieres-analysis-results/joint_genotype/2022-04-30_friedman_joint_metadata_04182022
genome_fasta=/shared/workspace/software/bcbio/genomes/Hsapiens/hg38/seq/hg38.fa
workspace=/scratch
vcf=menieres.norm.vep.vcf.gz
mkdir -p $workspace/$sample/maf
mkdir -p $workspace/$sample/vcf
echo $sample

if [ ! $(find $workspace/$vcf -type f -size +0c 2>/dev/null) ]; then
    aws s3 cp $s3/$vcf $workspace/
fi

$bcftools view -f.,PASS -c1 -s $sample -o $workspace/$sample/vcf/$sample.norm.vep.vcf $workspace/$vcf

#if [ ! $(find $workspace/$sample/maf/$sample.norm.vep.maf -type f -size +0c 2>/dev/null) ]; then
    perl `which vcf2maf.pl` \
        --tumor-id $sample \
        --normal-id $sample \
        --ref-fasta $genome_fasta \
        --ncbi-build GRCh38 \
        --inhibit-vep \
        --filter-vcf 0 \
	--retain-fmt GT,DP,GQ \
        --input-vcf $workspace/$sample/vcf/$sample.norm.vep.vcf \
        --output-maf $workspace/$sample/maf/$sample.norm.vep.maf
#fi

aws s3 cp $workspace/$sample/ $s3/$sample/ --recursive

