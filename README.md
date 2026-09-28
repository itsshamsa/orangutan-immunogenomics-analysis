# Orangutan Immunogenomics Analysis

Python scripts developed for the analysis and visualization of immunoglobulin genomic variation in orangutans.

This repository contains research scripts used for sequence-based analyses, switch-region characterization, haplotype comparison, gene-density analysis, and small sequence-processing utilities.

## Repository structure

### `switch_regions/`
Scripts for analysis of immunoglobulin switch regions, including:
- identification of candidate switch-like regions from FASTA sequences
- counting degenerate sequence motifs
- comparison of motif densities between species

### `haplotype_analysis/`
Scripts for comparing immunoglobulin loci across haplotypes and visualizing gene-level haplotype differences.

### `locus_analysis/`
Scripts for quantitative comparison and visualization of immunoglobulin gene densities across loci, genomes, and species.

### `utilities/`
Small sequence-analysis utilities, including extraction and generation of k-mers around specified nucleotide positions.

## Technologies

- Python
- Pandas
- NumPy
- Matplotlib
- Biological sequence analysis

## Research context

These scripts were developed as part of comparative immunogenomics analyses of *Pongo abelii* and *Pongo pygmaeus*, with a focus on immunoglobulin loci, genomic variation, and haplotype-resolved genome assemblies.

## Status

This repository is being organized and documented for reproducible research use.