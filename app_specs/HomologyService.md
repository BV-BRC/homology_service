
# Application specification: Homology

This is the application specification for service with identifier Homology.

The backend script implementing the application is [App-Homology.pl](../service-scripts/App-Homology.pl).

The raw JSON file for this specification is [Homology.json](Homology.json).

This service performs the following task:   Perform homology searches on sequence data

It takes the following parameters:

| id | label | type | required | default value |
| -- | ----- | ---- | :------: | ------------ |
| input_type | Type of input (dna or aa) | enum  | :heavy_check_mark: |  |
| input_source | Source of input (id_list, fasta_data, fasta_file, feature_group, genome_group_features, genome_group_sequences) | enum  | :heavy_check_mark: |  |
| input_fasta_data | Input sequence in fasta formats | string  |  |  |
| input_id_list | Input sequence as a list of sequence identifiers | array  |  |  |
| input_fasta_file | Input sequence as a workspace file of fasta data | wsid  |  |  |
| input_feature_group | Input sequence as a workspace feature group | wsid  |  |  |
| input_genome_group | Input sequence as a workspace genome group | wsid  |  |  |
| db_type | Database type to search (protein / DNA / RNA / contigs) | enum  | :heavy_check_mark: |  |
| db_source | Source of database (id_list, fasta_data, fasta_file, genome_list, taxon_list, feature_group, genome_group, precomputed_database) | enum  | :heavy_check_mark: |  |
| db_fasta_data | Database sequences as fasta | string  |  |  |
| db_fasta_file | Database fasta file | wsid  |  |  |
| db_id_list | Database as a list of sequence identifiers | array  |  |  |
| db_feature_group | Database feature group | wsid  |  |  |
| db_genome_group | Database genome group | wsid  |  |  |
| db_genome_list | Database genome list | array  |  |  |
| db_taxon_list | Database taxon list | array  |  |  |
| db_precomputed_database | Precomputed database name | string  |  |  |
| blast_program | BLAST program to use | enum  |  |  |
| blast_evalue_cutoff | BLAST  E-value cutoff | float  |  | 1e-05 |
| blast_max_hits | BLAST max-hits | int  |  | 300 |
| blast_min_coverage | BLAST mininum coverage | int  |  |  |
| output_path | Output Folder | folder  | :heavy_check_mark: |  |
| output_file | File Basename | wsid  | :heavy_check_mark: |  |

