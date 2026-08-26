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

Loading the sqlite catalog database (db.sqlite) from the individual genome files.

Everything referenced here lives in this module's scripts/ directory. The order
below is load-bearing: mk-lineage reads TaxonInDatabase, which does not exist
until p3x-create-databases-lookup has walked the database tree.

1. Download the NCBI taxdump and build the node table.

   p3x-create-taxonomy-nodes

   Creates taxdump-<YYYY-MM-DD>/ under the current directory and writes
   nodes.tsv into it, with columns tax_id,parent,rank. The file is
   comma-separated despite the name, which is why the import below needs
   .mode csv. The script refuses to run if the dated directory exists.

2. Create the database and load the taxonomy.

   sqlite3 db.sqlite
   sqlite> .read scripts/db-schema.sql
   sqlite> .mode csv
   sqlite> .import taxdump-<YYYY-MM-DD>/nodes.tsv TaxNode

3. Load the blast database catalog.

   p3x-create-databases-lookup --curated-directory ref --sqlite db.sqlite \
       /vol/blastdb/bvbrc-service blast.db

4. Build the lineage lookup.

   perl scripts/mk-lineage db.sqlite

Check afterwards that TaxNode has as many rows as nodes.tsv has lines. Note
that sqlite3 is not necessarily on PATH; perl -MDBI works for inspection.

Two things to be aware of when reading the result:

  - --curated-directory matches dirname() of the path relative to the database
    directory, so the argument has to be exactly that relative name. In
    data.2022-0916 it did not match and all 4058 GenomeGroup rows came out
    curated=0, which makes the "NOT g.curated" filter in
    BlastDatabasesSQL::search_taxa inert.

  - mk-lineage inserts a (taxon_id, NULL) row for any taxon in
    TaxonInDatabase that is absent from TaxNode, and those taxa then match
    nothing. data.2022-0916 has 63 of them. Do step 1 at catalog-build time rather
    than reusing a taxdump fetched at the start of the rebuild, so the
    taxonomy is never older than the genome set.
