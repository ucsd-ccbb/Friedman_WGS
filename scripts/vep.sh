#!/bin/bash

export PATH=/shared/workspace/software/bcbio/anaconda/bin:/shared/workspace/software/bcbio/tools/bin:$PATH
vcf=/scratch/joint/menieres-joint-gatk-haplotype-joint-annotated-vqsr-normed.vcf.gz.vcf.gz
genome_fasta=/shared/workspace/software/bcbio/genomes/Hsapiens/hg38/seq/hg38.fa
vep_cache=/shared/workspace/software/bcbio/genomes/Hsapiens/hg38/vep

export PATH=/shared/workspace/software/bcbio/anaconda/bin:/shared/workspace/software/bcbio/tools/bin:$PATH

`which vep` \
--vcf --verbose --cache \
--cache_version 100 \
--dir_cache $vep_cache \
--symbol --everything --fork 4 \
--no_intergenic --coding_only \
--pick --filter_common --no_stats \
--fasta $genome_fasta \
--sift b --polyphen b \
--ccds --symbol --numbers \
--domains --regulatory --canonical \
--protein --biotype --uniprot \
--tsl --appris --gene_phenotype \
--af --af_1kg --af_esp --af_gnomad \
--max_af --pubmed --variant_class \
--format vcf --compress_output bgzip \
-o menieres-joint-gatk-haplotype-joint-annotated-vqsr-normed-vep.vcf.gz -i $vcf
