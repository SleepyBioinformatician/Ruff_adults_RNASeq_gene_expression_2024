[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11032423.svg)](https://doi.org/10.5281/zenodo.11032423)

#  A single gene orchestrates androgen variation underlying male mating morphs
![ruffs](https://github.com/azemella/Ruff_adults_RNASeq_gene_expression_2024/assets/160619704/58acbe54-f214-4a29-aae5-bc000ec03dbd)

This repository contains data, metadata, data processing scripts, data analysis scripts and plotting scripts <br /> used in 
J. L. Loveland, A. Zemella *et al*., A single gene orchestrates androgen regulation underlying alternative <br /> mating morphs (2024)

## Before you start
Due to size limitations of Git repositories, some large files exceeding 500 MB are not included here. To ensure the complete reproducibility of all data analyses performed in this work, 
download the following external files and place them in the correct folders as indicated below:

### GitHub repository download:
Clone this GitHub repository by typing in the command line: <br /> 
- git clone https:<n/>//github.com/azemella/Ruff_adults_RNASeq_gene_expression_2024
- Rename the folder as 'ruff_adults_gene_expression'. This will be the working directory.

### *Calidris pugnax* (ruff) reference genome assembly [GCF_001431845.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001431845.1/) and related files:
- Genome assembly fasta file. Download the file in ./snakemake/genome/ncbi_genome/ and rename it as 'GCF_001431845.1_ASM143184v1_genomic.fna'
- Genome assembly fasta index file. Download the file in ./snakemake/genome/ncbi_genome/ and rename it as 'GCF_001431845.1_ASM143184v1_genomic.fna.fai'
- Genome annotation GTF file. Download the file in ./snakemake/genome/ncbi_genome/ and rename it as 'GCF_001431845.1_ASM143184v1_genomic.gtf'
- Save another copy of the same GTF file in ./metadata/ and rename it as 'NCBI_RefSeq_annotation.gtf'
- Genome annotation GFF3 file. Download the file in ./snakemake/genome/ncbi_genome/ and rename it as 'GCF_001431845.1_ASM143184v1_genomic.gff'
- Proteome fasta file. Download the file in ./eggNOG-mapper/ and rename it as 'GCF_001431845.1_ASM143184v1_proteome.faa'

### RNA-Sequencing libraries from [NCBI BioProject 1099138](https://www.ncbi.nlm.nih.gov/bioproject/1099138):
- Download all the 238 RNA-Seq libraries used in this work in the folder ./snakemake/data/
- Make sure the file names match those in the first column of the data frame ./metadata/metadata.csv

### Ruff BSgenome to be imported in karyoploteR:
There is no available *Calidris pugnax* full genome representation on [BSGenome](https://kasperdanielhansen.github.io/genbioconductor/html/BSgenome.html), a R package used to generate <br /> 
Biostrings-based genome data packages. This is required to replicate the plots displayed in Figure 1C-E made <br /> 
using the R package karyoploteR. To overcome this issue, we made a custom ruff BSgenome using the function <br /> 
forgeBSgenomeDataPkg() following the indications [in this tutorial](https://www.bioconductor.org/packages/release/bioc/vignettes/BSgenome/inst/doc/BSgenomeForge.pdf). Briefly, two files are needed to be created <br /> 
and stored in the same folder:
- A 2Bit ruff genome sequence data file named 'CalidrisPugnax.dna.2bit'
- A BSgenome data package seed file named 'BSgenome.Cpugnax.NCBI.GCF_001431845.1-seed' <br />
You can find this file in ./others/
- Follow the indications on page 9 of the tutorial to generate the custom BSgenome package and <br />
install it in R as a package

### Ruff Gene Ontology (GO) annotations for functional enrichment analysis:
Currently, there are no functional annotations available for *Calidris pugnax* genes in public databases. However, functional gene annotations can be inferred from orthologs using the online 
web tool [eggNOG-mapper](http://eggnog-mapper.embl.de/):
- Upload the file 'GCF_001431845.1_ASM143184v1_proteome.faa' that you can find in ./eggNOG-mapper/ to <br />
eggNOG-mapper (Proteins) and write your email address. Change the taxonomic scope to 'Vertebrata - 7742' <br />
in the annotations options
- Move all eggNOG-mapper output files to the folder eggNOG-mapper/eggNOG_functional_annotations/
- The main eggNOG-mapper output file used downstream is named 'out.emapper.annotations.csv'
- Download from [NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets/gene/GCF_001431845.1/) the gene annotations in ./metadata/ and rename the file as 'ncbi-gene-dataset.csv'
- Download from [NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets/gene/GCF_001431845.1/) the transcript annotations in ./metadata/ and rename the file as 'ncbi-transcript-dataset.csv'
  
The following steps are utilized to create a custom ruff Org.Db object and store all the previously downloaded information within it. To build the object, use the makeOrgPackage() function of the 
R package AnnotationForge as indicated in [this tutorial](https://bioconductor.org/packages/release/bioc/vignettes/AnnotationForge/inst/doc/MakingNewOrganismPackages.html):
- Run the R script 'create_files_for_makeorgdb.R' that you can find in ./makeorgdb/scripts/ to create all the necessary input files
- Run the R script 'makeorgdb.R' that you can find in ./makeorgdb/scripts/ to generate the *Calidris pugnax* Org.Db object necessary to perform GO functional enrichment analysis with clusterProfiler

### Note:
To comply with the journal's guidelines, in the final figures we edited some aspects (e.g., font size for axis titles or legend labels) for a few plots 
after generating them in R. These edits were conducted in Inkscape. We provide SVG files for both the original and edited versions. 
