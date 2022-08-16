#!/bin/bash
#$ -cwd
#$ -pe smp 16

export PATH=/shared/workspace/software/bcbio/anaconda/bin:/shared/workspace/software/bcbio/tools/bin:$PATH

sample=$sample
workspace=/scratch/germline/$sample # Edit this
template_yaml=/shared/workspace/projects/friedman/config/germline_template.yaml
s3Download=$s3Download # Edit this
s3Upload=$s3Upload # Edit this
mkdir -p $workspace/fastq

echo $metadata
metadata=$metadata
base_meta=$(basename $metadata)

if [ ! -f $workspace/"$sample"_R1.fastq.gz ]; then
	aws s3 cp $s3Download/ $workspace/fastq/ --recursive --exclude "*" --include "$sample*gz"
	cat $workspace/fastq/"$sample"*R1_001.fastq*gz > $workspace/"$sample"_R1.fastq.gz
	cat $workspace/fastq/"$sample"*R2_001.fastq*gz > $workspace/"$sample"_R2.fastq.gz
	rm -r $workspace/fastq/
fi

# Create wgs config file from template yaml
# If your metadata file is named wgs_metadata.csv,
# then a directory in your workspace will be created called wgs_metadata
# with a subdirectory called config, where wgs_metadata.yaml will be made.
cd $workspace && bcbio_nextgen.py \
	-w template $template_yaml \
	$metadata \
	$workspace \
	--workdir $workspace

bcbio_nextgen.py \
	$workspace/"${base_meta%.*}"/config/"${base_meta%.*}".yaml \
	-n 16 \
	--workdir $workspace

aws s3 cp $workspace/final/ $s3Upload/ --recursive
#rm -r $workspace
