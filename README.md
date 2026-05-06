# CIPN Variant Analysis Pipeline

## Overview
This repository contains scripts used for rare variant identification and pathway analysis in chemotherapy-induced peripheral neuropathy (CIPN).

## Workflow Summary
1. Alignment (BWA-MEM)
2. BAM processing (samtools, GATK)
3. Variant calling (HaplotypeCaller)
4. Filtering (GATK VQSR, HardFiltering)
5. Annotation (VEP)

## Requirements
- BWA
- Samtools
- GATK
- VEP
- bcftools

## Notes
- No raw sequencing data is included due to privacy restrictions.
- Two separate pipelines were designed to accommodate for BGI and IonTorrent AmpliSeq Data and are labeled accordingly.
