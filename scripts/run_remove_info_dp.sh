#!/bin/bash
s3Upload=s3://menieres-analysis-results
complete=$(aws s3 ls $s3Upload --recursive | grep fixed.vcf.gz$ | awk -F '/' '{print $(NF-1)}')

for sample in $(aws s3 ls s3://menieres-analysis-results/ | grep RF | awk '{print $NF}' | sed 's/\///g'); do
    if ! echo $complete | grep -w -q $sample; then
        qsub \
	-v sample=$sample \
        /shared/workspace/projects/friedman/scripts/remove_info_dp.sh
    else
        echo "Sample $sample is complete. Skipping."
    fi
done
