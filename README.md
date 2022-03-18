# homology_service
Callable service that provides access to a homolog database that covers prokaryotes, eukaryotes and metagenomes.  It can be used for navigation between sequences as well as for computation on protein families defined on the homologs.

Notes on installation:

This service requires per-genome and NR blast databases to be installed. See deploy.cfg for
the settings used to define their locations.

The default deployment will download the specific version of the BLAST+  command line tools 
required (this code uses the JSON output format which was not intially correct in BLAST).

====


Regarding the creation of non-redundant feature databases.

Especially in the larger genera, we can see a large performance boost if we
only save unique proteins. This has been done in the SEED for many years and
we are quite familiar with the technology to make it efficient.

The twist here is that we wish to use these databases with the support BLAST
has for requesting subset searches based on a set of taxonomy ids. When we create
the nonredundant databases, we erase the taxonomic identity of the proteins.

Here, we bring this back by computing a pseudo taxon ID that is stored in the
BLAST database.

It is constructed by computing for each unique protein the set of taxon IDs
that this protein occurs in. Call them T1, T2, ..., Tn.

We sort these IDs numerically, and concatenate them to create an aggregate
identifier string of the form "T1,T2,T3,...,Tn".

Each member of the set of aggregate identifiers is then assigned a unique
identifier S. The proteins in the BLAST database are each assigned the
appropriate S value in the field that BLAST uses to store taxon IDs (BLAST
itself does not interpret these values so we can assign meaning to them at will).

In order to construct the appropriate set of taxids to use to perform a subset
search with BLAST, we create mapping tables when we construct the BLAST databases.

For each taxon ID that appears in the database, we maintain a mapping table that
translates the taxon ID to the set of pseudo IDs that represent aggregate
IDs in which the taxon ID appears.

We also maintain a mapping table from the sequence ID ("lcl|md5-string") to the
set of protein feature IDs that have that sequence. This table is used to
expand the BLAST result.

We store these lookup tables in Berkeley DB files

    tax_to_pseudo  stores the taxon id to pseudo-id set (stored as records with multiple duplicate keys)

    md5_to_feature stores the mapping from md5 ID to feature ID set (stored as records with multiple duplicate keys)




===

create viral reference

p3x-create-blast-db --no-quality-check --title Viral\ Refs --parallel 2 aa features vout/tmpref --no-check-files --batch-size 500  --reference --representative  --taxon 10239 -f --viral

===

Loading sqlite database from the individual genome files

Set up db file

rm db.sqlite3
sqlite3 db.sqlite
> .read db-schema.sql

Create nodes.tsv in the NCBI taxdump directory (need notes on that)
Load taxonomy data

sqlite3 db.sqlite
sqlite> .mode csv
sqlite> .import nodes.tsv TaxNode


Load blast data
p3x-create-databases-lookup --curated-directory ref --sqlite db.sqlite /vol/blastdb/bvbrc-service blast.db

Create lineage lookup
perl ~/P3/dev-slurm/dev_container/modules/homology_service/mk-lineage /disks/tmp/blast.sqlite 
