# snvworkflow

[![DOI](https://zenodo.org/badge/151792411.svg)](https://zenodo.org/badge/latestdoi/151792411)


A pipeline to find allelic ratios in regions of clone-specific LOH in 10X scRNA-seq data. This uses a modified version of the [scAlleleCount](https://github.com/barkasn/scAlleleCount) script.

# Usage

To run `snvworkflow`, edit config.yaml with the required files. Briefly, you'll need:

- A BAM file from CellRanger of the scRNAseq
- A tsv file with sites of germline heterozygous SNPs
- A csv file detailing clone specific copy number regions

The `Dockerfile` is a template for creating a Docker container with all required software.