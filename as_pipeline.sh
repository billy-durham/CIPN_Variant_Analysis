#!/bin/bash

# Usage: ./ampseq_pipeline.sh SAMPLE_ID
sampleID=$1

# Paths
SAMPLE_DIR="/pvol/sample/${sampleID}"
REF="/pvol/HG/hg38.fa"
DBSNP="/pvol/HG/00-Allref.vcf.gz"
GATK="/pvol/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar"
LOG="${SAMPLE_DIR}/${sampleID}_pipeline.log"

# Sort hg19 BAM by read name
samtools sort -n -@ 8 -o ${SAMPLE_DIR}/${sampleID}_hg19_sorted.bam ${SAMPLE_DIR}/${sampleID}_hg19.bam
samtools index ${SAMPLE_DIR}/${sampleID}_hg19_sorted.bam

# Convert to FASTQ
samtools fastq -@ 8 ${SAMPLE_DIR}/${sampleID}_hg19_sorted.bam | gzip > ${SAMPLE_DIR}/${sampleID}.fastq.gz
fastqc -o ${SAMPLE_DIR} -t 4 ${SAMPLE_DIR}/${sampleID}.fastq.gz

# Align to hg38 / Sort
bwa mem -t 8 -R "@RG\tID:${sampleID}\tSM:${sampleID}\tPL:IONTORRENT" /pvol/HG/hg38.fa ${SAMPLE_DIR}/${sampleID}.fastq.gz | \
samtools sort -@ 8 -o ${SAMPLE_DIR}/${sampleID}_hg38_sorted.bam

samtools index ${SAMPLE_DIR}/${sampleID}_hg38_sorted.bam

# HaploTypeCaller
java -jar $GATK HaplotypeCaller \
  -R $REF \
  -I ${SAMPLE_DIR}/${sampleID}_hg38_sorted.bam \
  -O ${SAMPLE_DIR}/${sampleID}.g.vcf.gz \
  -ERC GVCF \
  --native-pair-hmm-threads 8


# Create Sample Map
java -jar /pvol/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar GenomicsDBImport \
  --genomicsdb-workspace-path /pvol/sample/all_gvcf/genomicsdb_workspace \
  --batch-size 50 \
  --sample-name-map /pvol/sample/all_gvcf/sample_map.txt \
  --reader-threads 8 \
  --intervals /pvol/sample/all_gvcf/hg38_exome_ccds.bed
  --merge-input-intervals


## Genotype GVCFS Together
java -jar /pvol/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar GenotypeGVCFs \
  -R /pvol/HG/hg38.fa \
  -V gendb:///pvol/sample/all_gvcf/genomicsdb_workspace \
  -O /pvol/sample/all_gvcf/ampliseq_jointcalled.vcf.gz

# Check File Integrity
bcftools view -h ampliseq_jointcalled.vcf.gz | head
bcftools index ampliseq_jointcalled.vcf.gz

### View Stats
bcftools stats ampliseq_jointcalled.vcf.gz > ampliseq_jointcalled.stats.txt
grep -E "^SN" ampliseq_jointcalled.stats.txt | column -t
grep TSTV ampliseq_jointcalled.stats.txt 

# Hard Filtering
gatk VariantFiltration   
  -R /pvol/HG/hg38.fa    
  -V ampliseq_jointcalled.vcf.gz    
  -O ampliseq_jointcalled_filtered.vcf.gz    
  --filter-name "QD_filter" --filter-expression "QD < 2.0"    
  --filter-name "MQ_filter" --filter-expression "MQ < 40.0"    
  --filter-name "DP_filter" --filter-expression "DP < 3"    
  --filter-name "ReadPos_filter" --filter-expression "ReadPosRankSum < -8.0"

### View Stats
bcftools view -f 'PASS,.' ampliseq_jointcalled_filtered.vcf.gz | grep -v "^#" | wc -l

### View Failed
bcftools query -f '%FILTER\n' ampliseq_jointcalled_filtered.vcf.gz | sort | uniq -c

# Sort and Index
bcftools sort ampliseq_jointcalled_filtered.vcf.gz -Oz -o ampliseq_jointcalled_sorted.vcf.gz
bcftools index ampliseq_jointcalled_sorted.vcf.gz 

# VEP Annotation
perl /pvol/ensembl-vep/vep \
  -i /pvol/sample/all_gvcf/ampliseq_jointcalled.sorted.vcf.gz \
  --plugin dbNSFP,/pvol/vep_cache/dbNSFP/dbNSFP4.8a_grch38.gz,ALL \
  --cache \
  --dir_cache /pvol/vep_cache/ \
  --dir_plugins ~/.vep/Plugins \
  --sift s --polyphen s --hgvs --symbol --canonical \
  --af_gnomade --af_gnomadg --pick \
  -o /pvol/sample/all_gvcf/ampliseq_jointcalled_vep.vcf \
  --force_overwrite --fork 8