# Homology Service - Command Line Tools

The homology service uses a set of underling service-side CLI tools to implement its operations. This document describes the tools and the data sets that back them up.

We support the following searches. The current PATRIC service offers them via BLAST and names them as BLAST services, but here we provide potentially multiple implementations. In particular, we use diamond for any searches against large data sets as it can be several orders of magnitude faster than BLAST.

## Searches Types Supported

The service supports all combinations of the standard homology searches.

For protein input:

- Search protein database
- Search translated nucleotide database

For nucleotide input:

- Search nucleotide database
- Search protein database using translated query
- Search translated nucleotide database using translated query

## Search Databases Supported

It also supports the following predefined databases. Genomic database are provided in protein feature, protein nucleotide, and contigs nucleotide formats. 

- NCBI reference genomes
- NCBI representative genomes
- NCBI representative and reference genomes
- Proteins with transcriptomic data (amino acid and nucleotide database)
- Phage genomes
- Plasmid genomes
- Proteins with specialty gene reference
- PATRIC 16s RNA genes (nucleotide database)

Users may also choose to create custom databases on the fly:

- Search a selected genome
- Search within a selected genome group
- Search within a taxon

