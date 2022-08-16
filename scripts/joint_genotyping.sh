#!/bin/bash
#$ -N joint_wgs
#$ -cwd
#$ -pe smp 24

export PATH=/shared/workspace/software/bcbio/anaconda/bin:/shared/workspace/software/bcbio/tools/bin:$PATH

template_yaml=/shared/workspace/projects/friedman/config/joint_genotyping_from_gvcf_template.yaml
s3Upload=s3://menieres-analysis-results/joint_genotype
mkdir -p $workspace

metadata=$metadata
base_meta=$(basename $metadata)

#for sample in $(awk -F ',' '{print $2}' $metadata | grep -v description); do
#	if [ ! -f $workspace/"$sample"-joint-gatk-haplotype-annotated-fixed.vcf.gz ]; then
#                aws s3 cp s3://menieres-analysis-results/$sample/$sample/"$sample"-joint-gatk-haplotype-annotated-fixed.vcf.gz $workspace/
#                aws s3 cp s3://menieres-analysis-results/$sample/$sample/"$sample"-joint-gatk-haplotype-annotated-fixed.vcf.gz.tbi $workspace/
#	fi
#done

cd $workspace &&  bcbio_nextgen.py \
	-w template $template_yaml \
	$metadata \
	$workspace/*vcf.gz \
	--workdir $workspace

bcbio_nextgen.py \
	$workspace/"${base_meta%.*}"/config/"${base_meta%.*}".yaml \
	-n 24 \
	--workdir $workspace

aws s3 cp $workspace/final/ $s3Upload/ --recursive
