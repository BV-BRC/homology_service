module HomologyService
{
    typedef structure {
	string id;
	string accession;
	string title;
	int taxid;
	string sciname;
    } HitDescr;

    typedef structure {
	int num;
	float bit_score;
	float score;
	float evalue;
	int identity;
	int positive;
	int density;
	int pattern_from;
	int pattern_to;
	int query_from;
	int query_to;
	string query_strand;
	int query_frame;
	int hit_from;
	int hit_to;
	string hit_strand;
	int hit_frame;
	int align_len;
	int gaps;
	string qseq;
	string hseq;
	string midline;
    } Hsp;

    typedef structure {
	int num;
	list<HitDescr> description;
	int len;
	list<Hsp> hsps;
    } Hit;

    typedef structure {
	int db_num;
	int db_len;
	int hsp_len;
	int eff_space;
	float kappa;
	float lambda;
	float entropy;
    } Statistics;

    typedef structure {
	string query_id;
	string query_title;
	int query_len;
	/* need: query-masking */
	list<Hit> hits;
	Statistics stat;
    } Search;

    typedef structure {
	string program;
	string version;
	string reference;
	structure {
	    string db;
	    string subjects;
	} search_target;
	structure {
	    string matrix;
	    float expect;
	    float include;
	    int sc_match;
	    int sc_mismatch;
	    int gap_open;
	    int gap_extend;
	    string filter;
	    string pattern;
	    string entrez_query;
	    int cbs;
	    int query_gencode;
	    int db_gencode;
	    string bl2seq_mode;
	} params;
	structure {
/*	    structure {
	    } iterations;
*/
	    Search search;
/*	    structure {
	    } bl2seq;
*/	    
	} result;
    } Report;

    typedef string genome_id;

    typedef structure
    {
	string function;
	string genome_name;
	string genome_id;
	string md5;
	string locus_tag;
	string alt_locus_tag;
	int match_count;
    } FeatureMetadata;

    typedef structure
    {
	float evalue_cutoff;
	int max_hits;
	float min_coverage;
    } BlastParameters;

    funcdef blast_fasta_to_genomes(string fasta_data,
				   string program,
				   list<genome_id> genomes,
				   /* subject_type is "contigs" or "features" */
				   string subject_type,
				   /* Post demo we will slot this in here.
				      BlastParameters blast_parameters */
				   float evalue_cutoff,
				   int max_hits,
				   float min_coverage)
	returns(list<Report> reports, mapping<string, FeatureMetadata> metadata) authentication optional;

    funcdef blast_fasta_to_taxon(string fasta_data,
				   string program,
				   string taxon_id,
				   /* subject_type is "contigs" or "features" */
				   string subject_type,
				   /* Post demo we will slot this in here.
				      BlastParameters blast_parameters */
				   float evalue_cutoff,
				   int max_hits,
				   float min_coverage)
	returns(list<Report> reports, mapping<string, FeatureMetadata> metadata);

    typedef structure
    {
	string name;
	string key;
	/* db_type is either "dna" or "protein" */
	string db_type;
	int seq_count;
    } DatabaseDescription;

    funcdef enumerate_databases() returns (list<DatabaseDescription>);

    funcdef blast_fasta_to_database(string fasta_data, string program, string database_key,
				   /* Post demo we will slot this in here.
				      BlastParameters blast_parameters */
				    float evalue_cutoff,
				    int max_hits,
				    float min_coverage)
	returns (list<Report> reports,
		 mapping <string, FeatureMetadata> metadata,
		 mapping <string, list<tuple<string, FeatureMetadata>>> identical_proteins);
				    
};
