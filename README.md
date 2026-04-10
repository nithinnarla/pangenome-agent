# The Pangenome Agent

## Overview
An automated pipeline for pangenome graph analysis to facilitate the extraction of structural variants of regions of interest from graph structures. 

This tool levarages AI agents as bridges between the user and their data, where minimal coding is needed to extract desired information from pangenome graphs. 

Simply ask the agent "Which variants are present between region X and Z?" and the AI will orchestrate the rest. Mostly importantly, the data never leaves your computer. 

All the bioinformatic pipelines included as skills have been created by humans (Not AI) and are fully auditable. All the running logs and raw outputs are available to the user, as well as summarized and user-friendly versions for a quicker first look. 

## Research Question
Can agentic computational frameworks improve the efficiency and interpretability of pangenome variant analysis?

## Background
For decades, the "linear reference" paradigm dominated genomics. However, a single reference genome often fails to capture the entire genetic diversity present in a population. Pangenomes address this by representing genetic diversity as a variation graphs.

Variation graphs encode genomic information into a format called GFA (Graphical Fragment Assembly), which is less user-friendly than other genomic data formats (GFF,FASTA,VCF) because it's not human-readable. 

Therefore, parsing information from GFAs is not a simple task. It involves a series of steps from converting linear genome coordinates to graph positions, sorting and parsing files, to filtering relevant information to the user.

Agentic genomics has recently emerged as a promising field in bioinformatics, filling the knowledge gap of a fast-evolving computational landscape where technologies advance faster than scientist can learn them. 

In this framework, AI Agents act as "orchestrators". Unlike traditional "black-box" AI that might attempt to hallucinate sequence data, these agents utilize Large Language Models (LLMs) to reason from the user needs to triggering trusted bioinformatic tools that will generate the data.

The development of the Pangenome Agent is built upon the architectural foundations of ClawBio, the first bioinformatics-native AI agent skill library. While ClawBio provides a library of general pharmacological skills, our agent introduces specialized skills focused on pangenome analysis.

## Input Files

| Format | Extension | Purpose |
|--------|-----------|---------|
| Variant Call Format | `.vcf` | Genomic variants |
| Gene Feature Format | `.gff` | Gene annotations |
| Graph Format | `.gfa` | Pangenome graph structure |
| FASTA sequence | `.fasta` | DNA sequences |

## Pipeline

1. Parse and integrate GFA, VCF, GFF and FASTA files
2. Extract variants in the graph from linear chromosomal position
3. Cross-reference variants against GWAS databases (FUTURE)
4. Classify variant significance using ACMG criteria (FUTURE)
5. Map variants to gene features using GFF annotations (FUTURE)

## Output

output/
├── report.md
├── result.json
├── variants.csv
└── reproducibility/
└── commands.sh

## Tech Stack
ClawBio, Python, biopython, pandas, bioinformatics tools via Galaxy bridge, numpy, bcftools, VCF tools, GWAS databases, Galaxy bioinformatics suite

## Research Timeline
- Late March 2026: Project initiated
- April 2026: Prototype development and testing
- May-December 2026: Optimization and validation
- Spring 2027: Target submission

## Status
🔬 Active development
Target venue: Bioinformatics / ISMB / PLOS Computational Biology 2027

## Authors
Lucas Borges dos Santos - Conceptualization, Pipeline development, Methodology design
Nithin Narla - Agent implementation, Code optimization
University of Illinois Urbana-Champaign

## License
MIT