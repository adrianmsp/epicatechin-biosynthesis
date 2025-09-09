# Evolutionary Insigths Into (-)-Epicatechin Biosynthesis
This repository contains Python and R scripts for the evolutionary and functional analysis of (-)-epicatechin biosynthesis, including orthogroup extraction from OrthoFinder, protein sequence counting, post-BLAST hit extraction, species presence-absence matrix per orthogroup, and gene expression analysis in Arabidopsis thaliana. 
## Features

### Python Scripts
1. **Orthogroup Sequence Extraction** (´orthogroup_extraction.py´)
- Identifies orthogroups containing target protein IDs from *Arabidopsis thaliana*.
- Saves one FASTA file per orthogroup with standarised headers.
2. **Counter** (´orthogroup_counter.py´)
- Process output Fasta files from orthogroup_extraction.py.
- Counts protein sequence per orthogroup.
- Renames orthogroups with enzyme annotations.
- Creates a bar plot highlighting the enzymes of the phenylpropanoid pathway and flavonoid biosynthesis.
3. **Extract BLAST** (´extract_blast.py´)
- Process BLAST.tsv results
- Extracts the best-scoring hit per species.
4. **Species Matrix** (´species_matrix.py´)
- Builds a species presence-absence matrix from the orthogroups FASTA files.
- Visualises orthogroup distribution as a heatmap/matrix.
- Creates Clade Matrix and saves it as an excel file.

### R Scripts
1. **Gene Expression Analysis** (´gene_expression.R´)
- Analyses gene expression patterns of (-)-epicatechin biosynthesis genes in *Arabidopsis thaliana*.
- Filters and extracts baseline expression data for key tissues: Leaf, Flower, Fruit, and Root.
- Normalises expression values with log2 transformation.
- Generates a heatmap visualisation.
  
