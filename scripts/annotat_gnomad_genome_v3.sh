# Annotate with gnomad genome v3
export PATH=/shared/workspace/software/bcbio/anaconda/envs/bcftools1.15/bin:/shared/workspace/software/bcbio/anaconda/bin:$PATH

genome_fasta=/shared/workspace/software/bcbio/genomes/Hsapiens/hg38/seq/hg38.fa
vep_cache=/shared/workspace/software/bcbio/genomes/Hsapiens/hg38/vep
workspace=/scratch

# Convert filtered variants back to VCF
maf2vcf.pl --input-maf $workspace/menieres.531.filtered.tsv --output-dir $workspace --output-vcf $workspace/menieres.531.filtered.vcf --ref-fasta $genome_fasta
bgzip $workspace/menieres.531.filtered.vcf
bcftools sort -Oz -o $workspace/menieres.531.filtered.sorted.vcf.gz -T $workspace $workspace/menieres.531.filtered.vcf.gz
tabix -p vcf $workspace/menieres.531.filtered.sorted.vcf.gz

# Annotate VCF with gnomad genome v3.1.2 https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr"$chr".vcf.bgz
parallel "
	cd $workspace && wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr{}.vcf.bgz.tbi
	bcftools view --regions chr{} $workspace/menieres.531.filtered.sorted.vcf.gz \
	-Oz -o $workspace/menieres.531.filtered.sorted.chr{}.vcf.gz
	tabix -p vcf $workspace/menieres.531.filtered.sorted.chr{}.vcf.gz
	bcftools annotate \
	-a https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr{}.vcf.bgz \
	-c gnomad_v3_af_popmax:=AF_popmax \
	-Oz -o $workspace/menieres.531.gnomadv3.chr{}.vcf.gz \
	$workspace/menieres.531.filtered.sorted.chr{}.vcf.gz" ::: {1..22} X Y


# Concatenate the chr back to one VCF
bcftools concat \
$workspace/menieres.531.gnomadv3.chr*.vcf.gz | \
bcftools sort \
-Oz -o $workspace/menieres.531.gnomadv3.vcf.gz
tabix -p vcf $workspace/menieres.531.gnomadv3.vcf.gz

# Reannotate with VEP
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
-o $workspace/menieres.531.gnomadv3.vep.vcf.gz \
-i $workspace/menieres.531.gnomadv3.vcf.gz 
tabix -p vcf $workspace/menieres.531.gnomadv3.vep.vcf.gz

# Separate each sample out and convert to MAF again
for sample in $(bcftools query -l $workspace/menieres.531.gnomadv3.vep.vcf.gz); do
	bcftools view $workspace/menieres.531.gnomadv3.vep.vcf.gz -c1 -f.,PASS -s $sample \
	-o $workspace/$sample.gnomadv3.vep.vcf
    perl `which vcf2maf.pl` \
        --tumor-id $sample \
        --normal-id $sample \
        --ref-fasta $genome_fasta \
        --ncbi-build GRCh38 \
        --inhibit-vep \
        --filter-vcf 0 \
        --retain-info gnomad_v3_af_popmax \
        --retain-fmt GT,DP \
        --input-vcf $workspace/$sample.gnomadv3.vep.vcf \
        --output-maf $workspace/$sample.gnomadv3.vep.maf
done



