# Usage for alien_hunter_to_proteins.py
This script is intended to match the Alien_Hunter HGT regions to an annotated genome.  

## Setup & Requirements 
Run the python script in the directory containing your Alien_Hunter output files and the required files mentioned below.  
  
**Requirements:**
1. A .gbff file for your assembly (e.g., downloaded from NCBI) **named "{prefix}_genomic.gbff"**
2. A .faa file for your assembly **named "{prefix}_proteins.faa"**

## Output
1. *`{prefix}_HGT_gene_hits.tsv`*: A file containing information for the genes found within putative HGT regions
2. *`{prefix}_HGT_proteins.faa`*: A file containing protein sequences for genes found within putative HGT regions
