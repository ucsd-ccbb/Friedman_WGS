#! /bin/bash

genome=/shared/software/genomes/hg38/seq/hg38.fa
variation_dir=/shared/software/genomes/hg38/variation
workspace=/shared/projects/friedman/20210120_friedman_harris_Meniueres_WGS/workspace

gatk --java-options "-Djava.io.tmpdir=/tmp -Xms15G -Xmx30G -XX:ParallelGCThreads=2" VariantRecalibrator \
  -tranche 100.0 -tranche 99.95 -tranche 99.9 \
  -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 \
  -tranche 95.0 -tranche 94.0 \
  -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
  -R $genome \
  -V $workspace/menieres-joint-gatk-haplotype-joint-annotated.vcf.gz \
  --resource:hapmap,known=false,training=true,truth=true,prior=15.0 \
  $variation_dir/hapmap_3.3.vcf.gz  \
  --resource:omni,known=false,training=true,truth=false,prior=12.0 \
  $variation_dir/1000G_omni2.5.vcf.gz \
  --resource:1000G,known=false,training=true,truth=false,prior=10.0 \
  $variation_dir/1000G_phase1.snps.high_confidence.vcf.gz \
  -an QD -an FS -an SOR -an DP  \
  -mode SNP -O $workspace/merged_SNP1.recal --tranches-file $workspace/output_SNP1.tranches \
  --rscript-file $workspace/output_SNP1.plots.R

gatk --java-options "-Djava.io.tmpdir=/tmp -Xms15G -Xmx30G -XX:ParallelGCThreads=2" VariantRecalibrator \
  -tranche 100.0 -tranche 99.95 -tranche 99.9 \
  -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 \
  -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 \
  -tranche 92.0 -tranche 91.0 -tranche 90.0 \
  -R $genome \
  -V $workspace/menieres-joint-gatk-haplotype-joint-annotated.vcf.gz \
  --resource:mills,known=false,training=true,truth=true,prior=12.0 \
  $variation_dir/Mills_and_1000G_gold_standard.indels.vcf.gz \
  --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 \
  $variation_dir/dbsnp-153.vcf.gz \
  -an QD -an FS -an SOR -an DP \
  -mode INDEL -O $workspace/merged_indel1.recal --tranches-file $workspace/output_indel1.tranches \
  --rscript-file $workspace/output_indel1.plots.R

# run ApplyVQSR on SNP then INDEL

gatk --java-options "-Djava.io.tmpdir=/tmp -Xms15G -Xmx30G -XX:ParallelGCThreads=2" ApplyVQSR \
  -V $workspace/menieres-joint-gatk-haplotype-joint-annotated.vcf.gz \
  --recal-file $workspace/merged_SNP1.recal \
  -mode SNP \
  --tranches-file $workspace/output_SNP1.tranches \
  --truth-sensitivity-filter-level 99.9 \
  --create-output-variant-index true \
  -O $workspace/SNP.recalibrated_99.9.vcf.gz

gatk --java-options "-Djava.io.tmpdir=/tmp -Xms30G -Xmx60G -XX:ParallelGCThreads=2" ApplyVQSR \
  -V $workspace/SNP.recalibrated_99.9.vcf.gz \
  -mode INDEL \
  --recal-file $workspace/merged_indel1.recal \
  --tranches-file $workspace/output_indel1.tranches \
  --truth-sensitivity-filter-level 99.9 \
  --create-output-variant-index true \
  -O $workspace/indel.SNP.recalibrated_99.9.vcf.gz

