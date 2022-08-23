#sample=$1
for sample in $(cat /shared/workspace/projects/friedman/metadata/531_samples.txt); do
	qsub -v sample=$sample \
	-N $sample \
	/shared/workspace/projects/friedman/scripts/vcf2maf.sh
done
