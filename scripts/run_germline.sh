#!/bin/bash
batch=$1
s3Download=s3://menieres-raw-data/batch_$batch
s3Upload=s3://menieres-analysis-results
complete=$(aws s3 ls $s3Upload | grep PRE | awk '{print $NF}')
inprogress=$(qstat | grep ubuntu | awk '{print $3}')

for sample in $(cat /shared/workspace/projects/friedman/metadata/friedman_wgs_samples_batch$batch.txt); do
    if ! echo $complete $inprogress | grep -w -q $sample; then
        metadata=/shared/workspace/projects/friedman/metadata/20220125_friedman_wgs_"$sample".csv
        echo -e "samplename,description" > $metadata
        echo "/scratch/germline/$sample/"$sample"_R1.fastq.gz,$sample" >> $metadata
        echo "/scratch/germline/$sample/"$sample"_R2.fastq.gz,$sample" >> $metadata
        qsub \
        -N "$sample" \
        -v sample=$sample \
        -v s3Download=$s3Download \
        -v s3Upload=$s3Upload/$sample \
        -v metadata=$metadata \
        /shared/workspace/projects/friedman/scripts/germline.sh
    else
        echo "Sample $sample is complete. Skipping."
    fi
done
