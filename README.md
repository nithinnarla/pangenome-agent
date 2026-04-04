# Pangenome Agent 🌱

## Overview
A computational pipeline for soybean pangenome variant analysis targeting SCN (Soybean Cyst Nematode) resistance detection across 9 soybean genomes. The pipeline integrates bioinformatics tools with AI-assisted variant interpretation to accelerate genomic discovery in agricultural research.

## Research Question
Can multi-agent computational frameworks improve the efficiency and interpretability of pangenome variant analysis for SCN resistance detection in soybean breeding programs?

## Background
Soybean Cyst Nematode is the most economically damaging soybean pathogen in the United States. Identifying resistance-associated variants across multiple soybean genomes is a computationally intensive task requiring integration of variant calling, gene annotation, and
functional interpretation. This pipeline addresses that challenge by combining established bioinformatics workflows with AI-assisted analysis.

## Input Files

| Format | Extension | Purpose |
|--------|-----------|---------|
| Variant Call Format | `.vcf` | Genomic variants |
| Gene Feature Format | `.gff` | Gene annotations |
| Graph Format | `.gfa` | Pangenome graph structure |
| FASTA sequence | `.fasta` | DNA sequences |

## Pipeline

1. Parse and integrate VCF, GFF, GFA and FASTA files
2. Filter variants by chromosomal position and quality
3. Cross-reference variants against GWAS databases
4. Classify variant significance using ACMG criteria
5. Map variants to gene features using GFF annotations
6. Generate structured variant report with SCN
   resistance findings

## Output

output/
├── report.md
├── result.json
├── variants.csv
└── reproducibility/
└── commands.sh

## Tech Stack
Python, biopython, pandas, ClawBio, bioinformatics tools via Galaxy bridge, numpy, bcftools, VCF tools, GWAS databases, Galaxy bioinformatics suite

## Research Timeline
- March 2026: Project initiated
- April 2026: Pipeline development and testing
- May-December 2026: Optimization and validation
- Spring 2027: Target submission

## Status
🔬 Active development
Target venue: Bioinformatics / ISMB / PLOS Computational Biology 2027

## Authors
Nithin Narla, Lucas Borges dos Santos
University of Illinois Urbana-Champaign

## License
MIT
