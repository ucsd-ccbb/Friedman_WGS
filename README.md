# Friedman_WGS

This repository should demonstrate how to transform raw WGS fastq files into actionable, rare, deleterious variants in the context of Meniere's disease.

##### Steps

1. Call variants on each sample using bcbio, implementing GATK's best practices.
2. Fix header of VCF files for joint genotyping (possibly a bcbio fault).
3. Joint genotype all samples.
4. Variant Quality Score Recalibration (VQSR).
5. Normalize variants - store variants at one record per line in VCF.
6. Annotate with VEP.
7. Convert VEP VCF to individual sample MAF.
8. Filter variants. `for maf in $(ls */*/*maf); do out=$(echo $maf | sed 's/vep/vep.coding/'); awk -F '\t' '{ if($1 != "Unknown") { print }}' $maf > $out; done`
9. Annotate with gnomAD GENOME.
10. Further filter for rare and deleterious variants.

### 1. Variant Calling

You will first have to create a list of samples to run, and then execute `run_germline.sh`. This script submits `germline.sh` for each sample.

```
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
```
### 2. Fix VCF header

For some reason, bcbio creates both INFO/DP and FORMAT/DP in the VCF header. The joint genotyping code doesn't like this, so we have code to fix it. Execute `run_remove_info_dp.sh`, which runs `remove_info_dp.sh` for each sample and looks like this:
```
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
```
All it does is remove the INFO/DP part of the VCF header.

### 3. Joint Genotyping
Make a metadata file named `friedman_joint_metadata_04182022.csv` containing the path of the fixed VCF, sample name, and batch name ("menieres-joint"). There should be one row per sample, formatted as follows:
```
samplename,description,batch
/scratch/joint/RF_0001-joint-gatk-haplotype-annotated-fixed.vcf.gz,RF_0001,menieres-joint
/scratch/joint/RF_0002-joint-gatk-haplotype-annotated-fixed.vcf.gz,RF_0002,menieres-joint
```
Add the path of the metadata file to `run_joint.sh` and execute. It will run bcbio's joint genotyping pipeline.
```
#!/bin/bash

metadata=/shared/workspace/projects/friedman/metadata/friedman_joint_metadata_04182022.csv
workspace=/scratch/joint

qsub \
-v workspace=$workspace \
-v metadata=$metadata \
/shared/workspace/projects/friedman/scripts/joint_genotyping.sh
```
### 4. VQSR

We have to run Variant Quality Score Recalibration after joint genotyping. It may be a parameter you can add in bcbio for joint genotyping, but that was realized after the fact, so we are running it manually.

### 5. Normalize Variants and 6. Annotate with VEP

The purpose of this is to normalize the VCF, so only one variant allele is represented per line so that when the variants are annotated, each gets their own true annotation because VEP only assigns one annotation per line.  Execute `norm_vep.sh` to normalize and annotate with VEP.
```
chr1    15484   .       G       A,T
```
becomes
```
chr1    15484   .       G       A
chr1    15484   .       G       T
```
### 7. Convert VEP VCF to individual sample MAF.
Run `vcf2maf.sh` for each sample.

### 8. Filter Variants
Copy MAF from S3 to local directory and run:
```
for maf in $(ls */*/*maf)
    do out=$(echo $maf | sed 's/vep/vep.coding/'); awk -F '\t' '{ if($1 != "Unknown") { print }}' $maf > $out
done
```

### 9. Annotate with gnomAD GENOME
We have to annotate the variants with gnomAD Genome because VEP only annotates with gnomAD Exome. We need the WGS allele frequencies because there are huge discrepencies. To do this we have to:

    1. Combine all sample MAF files into one cohort MAF. Can do it in R with do.call(rbind, lapply(maf_files, fread))
   
    2. Convert MAF back to VCF so we can annotate.

    3. Annotate using bcftools to query gnomAD's remote DB. This takes a while so thats why we are doing it on filtered variants.
    
    4. Concatenate and sort variants back together.
    
    5. Reannotate with VEP.
    
Steps 2-5 are done in `maf2vcf_gnomadv3.sh`

### 10. Further filter for rare and deleterious variants.
Filter for final set and visualize. This is done in RStudio with `maftools` in `menierees_variant_summary.R`

NOTE: You can see this isn't fully automated and there are some manual linker steps required.
