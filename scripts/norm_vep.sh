#!/bin/bash
#$ -pe smp 32

export PATH=/shared/workspace/software/bcbio/anaconda/bin:/shared/workspace/software/bcbio/tools/bin:$PATH

s3=s3://menieres-analysis-results/joint_genotype/2022-04-30_friedman_joint_metadata_04182022
vcf=menieres-joint-gatk-haplotype-joint-annotated-vqsr.vcf.gz
fixed_vcf=$(echo $vcf | sed 's/vcf.gz/fixed.vcf.gz/')
genome_fasta=/shared/workspace/software/bcbio/genomes/Hsapiens/hg38/seq/hg38.fa
vep_cache=/shared/workspace/software/bcbio/genomes/Hsapiens/hg38/vep
workspace=/scratch/$chr
mkdir -p $workspace/

aws s3 cp $s3/$vcf $workspace/
aws s3 cp $s3/$vcf.tbi $workspace/

bcftools view -h $workspace/$vcf > $workspace/hdr.txt
sed -i -E -e 's/Number=1|Number=R|Number=A/Number=./g' $workspace/hdr.txt
bcftools reheader -h $workspace/hdr.txt -o $workspace/$fixed_vcf $workspace/$vcf
tabix -p vcf $workspace/$fixed_vcf

bcftools norm -m - -c s --fasta-ref $genome_fasta -Oz \
-o $workspace/menieres-joint-gatk-haplotype-joint-annotated-vqsr-normed.vcf.gz \
$workspace/$fixed_vcf
tabix -p vcf $workspace/menieres-joint-gatk-haplotype-joint-annotated-vqsr-normed.vcf.gz

`which vep` \
	-i $workspace/menieres-joint-gatk-haplotype-joint-annotated-vqsr-normed.vcf.gz \
	-o $workspace/menieres-joint-gatk-haplotype-joint-annotated-vqsr-normed-vep.vcf \
	--vcf \
	--verbose \
	--cache \
	--cache_version 100 \
	--dir_cache $vep_cache \
	--symbol \
	--everything \
	--fork 4 \
	--no_intergenic \
	--coding_only \
	--pick \
	--no_stats \
	--fasta $genome_fasta \
	--sift b --polyphen b --ccds \
	--symbol --numbers --domains \
	--regulatory --canonical --protein \
	--biotype --uniprot --tsl \
	--appris --gene_phenotype --af \
	--af_1kg --af_esp --af_gnomad \
	--max_af --pubmed --variant_class --format vcf
bgzip $workspace/$workspace/menieres-joint-gatk-haplotype-joint-annotated-vqsr-normed-vep.vcf
tabix -p vcf $workspace/$workspace/menieres-joint-gatk-haplotype-joint-annotated-vqsr-normed-vep.vcf.gz
aws s3 cp $workspace/$workspace/menieres-joint-gatk-haplotype-joint-annotated-vqsr-normed-vep.vcf.gz $s3/
aws s3 cp $workspace/$workspace/menieres-joint-gatk-haplotype-joint-annotated-vqsr-normed-vep.vcf.gz.tbi $s3/
