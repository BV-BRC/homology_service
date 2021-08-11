# homology_service
Callable service that provides access to a homolog database that covers prokaryotes, eukaryotes and metagenomes.  It can be used for navigation between sequences as well as for computation on protein families defined on the homologs.

Notes on installation:

This service requires per-genome and NR blast databases to be installed. See deploy.cfg for
the settings used to define their locations.

The default deployment will download the specific version of the BLAST+  command line tools 
required (this code uses the JSON output format which was not intially correct in BLAST).

