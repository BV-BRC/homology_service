package Bio::KBase::HomologyService::HomologyServiceImpl;
use strict;
# Use Semantic Versioning (2.0.0-rc.1)
# http://semver.org 
our $VERSION = "0.1.0";

=head1 NAME

HomologyService

=head1 DESCRIPTION



=cut

#BEGIN_HEADER

use Data::Dumper;
use Bio::KBase::HomologyService::Util;
use Bio::P3::DeploymentConfig;
use P3DataAPI;

#END_HEADER

sub new
{
    my($class, @args) = @_;
    my $self = {
    };
    bless $self, $class;
    #BEGIN_CONSTRUCTOR

    my $cfg = Bio::P3::DeploymentConfig->new($ENV{KB_SERVICE_NAME} || "HomologyService");

    $self->{_blast_db_genomes} = $cfg->setting('blast-db-genomes');
    $self->{_blast_db_private_genomes} = $cfg->setting('blast-db-private-genomes');
    $self->{_blast_db_databases} = $cfg->setting('blast-db-databases');

    $self->{_data_api} = P3DataAPI->new();
    delete $self->{data_api}->{token};

    my $thr = $cfg->setting('blast-threads');
    if ($thr =~ /^\s*(\d+)\s*$/)
    {
	$self->{_blast_threads} = $1;
    }
    elsif ($thr)
    {
	warn "Invalid format for blast-threads: '$thr'. Ignoring this value.\n";
    }

    #
    # Determine if there is a custom blast installed in the service
    # directory. If so default the prefix to that.
    #
    my $svcdir = $ENV{KB_SERVICE_DIR};
    if (-x "$svcdir/bin/blastp")
    {
	$self->{_blast_program_prefix} = "$svcdir/bin/";
    }
    
    $self->{_blast_program_suffix} = $cfg->setting('blast-program-suffix') if $cfg->setting('blast-program-suffix');
    $self->{_blast_program_prefix} = $cfg->setting('blast-program-prefix') if $cfg->setting('blast-program-prefix');

    my $util = Bio::KBase::HomologyService::Util->new($self);
    $self->{_util} = $util;

    #END_CONSTRUCTOR

    if ($self->can('_init_instance'))
    {
	$self->_init_instance();
    }
    return $self;
}
=head1 METHODS
=head2 blast_fasta_to_genomes

  $reports, $metadata = $obj->blast_fasta_to_genomes($fasta_data, $program, $genomes, $subject_type, $evalue_cutoff, $max_hits, $min_coverage)

=over 4


=item Parameter and return types

=begin html

<pre>
$fasta_data is a string
$program is a string
$genomes is a reference to a list where each element is a genome_id
$subject_type is a string
$evalue_cutoff is a float
$max_hits is an int
$min_coverage is a float
$reports is a reference to a list where each element is a Report
$metadata is a reference to a hash where the key is a string and the value is a FeatureMetadata
genome_id is a string
Report is a reference to a hash where the following keys are defined:
	program has a value which is a string
	version has a value which is a string
	reference has a value which is a string
	search_target has a value which is a reference to a hash where the following keys are defined:
	db has a value which is a string
	subjects has a value which is a string

	params has a value which is a reference to a hash where the following keys are defined:
	matrix has a value which is a string
	expect has a value which is a float
	include has a value which is a float
	sc_match has a value which is an int
	sc_mismatch has a value which is an int
	gap_open has a value which is an int
	gap_extend has a value which is an int
	filter has a value which is a string
	pattern has a value which is a string
	entrez_query has a value which is a string
	cbs has a value which is an int
	query_gencode has a value which is an int
	db_gencode has a value which is an int
	bl2seq_mode has a value which is a string

	result has a value which is a reference to a hash where the following keys are defined:
	search has a value which is a Search

Search is a reference to a hash where the following keys are defined:
	query_id has a value which is a string
	query_title has a value which is a string
	query_len has a value which is an int
	hits has a value which is a reference to a list where each element is a Hit
	stat has a value which is a Statistics
Hit is a reference to a hash where the following keys are defined:
	num has a value which is an int
	description has a value which is a reference to a list where each element is a HitDescr
	len has a value which is an int
	hsps has a value which is a reference to a list where each element is a Hsp
HitDescr is a reference to a hash where the following keys are defined:
	id has a value which is a string
	accession has a value which is a string
	title has a value which is a string
	taxid has a value which is an int
	sciname has a value which is a string
Hsp is a reference to a hash where the following keys are defined:
	num has a value which is an int
	bit_score has a value which is a float
	score has a value which is a float
	evalue has a value which is a float
	identity has a value which is an int
	positive has a value which is an int
	density has a value which is an int
	pattern_from has a value which is an int
	pattern_to has a value which is an int
	query_from has a value which is an int
	query_to has a value which is an int
	query_strand has a value which is a string
	query_frame has a value which is an int
	hit_from has a value which is an int
	hit_to has a value which is an int
	hit_strand has a value which is a string
	hit_frame has a value which is an int
	align_len has a value which is an int
	gaps has a value which is an int
	qseq has a value which is a string
	hseq has a value which is a string
	midline has a value which is a string
Statistics is a reference to a hash where the following keys are defined:
	db_num has a value which is an int
	db_len has a value which is an int
	hsp_len has a value which is an int
	eff_space has a value which is an int
	kappa has a value which is a float
	lambda has a value which is a float
	entropy has a value which is a float
FeatureMetadata is a reference to a hash where the following keys are defined:
	function has a value which is a string
	genome_name has a value which is a string
	genome_id has a value which is a string
	md5 has a value which is a string
	locus_tag has a value which is a string
	alt_locus_tag has a value which is a string
	match_count has a value which is an int
</pre>

=end html

=begin text

$fasta_data is a string
$program is a string
$genomes is a reference to a list where each element is a genome_id
$subject_type is a string
$evalue_cutoff is a float
$max_hits is an int
$min_coverage is a float
$reports is a reference to a list where each element is a Report
$metadata is a reference to a hash where the key is a string and the value is a FeatureMetadata
genome_id is a string
Report is a reference to a hash where the following keys are defined:
	program has a value which is a string
	version has a value which is a string
	reference has a value which is a string
	search_target has a value which is a reference to a hash where the following keys are defined:
	db has a value which is a string
	subjects has a value which is a string

	params has a value which is a reference to a hash where the following keys are defined:
	matrix has a value which is a string
	expect has a value which is a float
	include has a value which is a float
	sc_match has a value which is an int
	sc_mismatch has a value which is an int
	gap_open has a value which is an int
	gap_extend has a value which is an int
	filter has a value which is a string
	pattern has a value which is a string
	entrez_query has a value which is a string
	cbs has a value which is an int
	query_gencode has a value which is an int
	db_gencode has a value which is an int
	bl2seq_mode has a value which is a string

	result has a value which is a reference to a hash where the following keys are defined:
	search has a value which is a Search

Search is a reference to a hash where the following keys are defined:
	query_id has a value which is a string
	query_title has a value which is a string
	query_len has a value which is an int
	hits has a value which is a reference to a list where each element is a Hit
	stat has a value which is a Statistics
Hit is a reference to a hash where the following keys are defined:
	num has a value which is an int
	description has a value which is a reference to a list where each element is a HitDescr
	len has a value which is an int
	hsps has a value which is a reference to a list where each element is a Hsp
HitDescr is a reference to a hash where the following keys are defined:
	id has a value which is a string
	accession has a value which is a string
	title has a value which is a string
	taxid has a value which is an int
	sciname has a value which is a string
Hsp is a reference to a hash where the following keys are defined:
	num has a value which is an int
	bit_score has a value which is a float
	score has a value which is a float
	evalue has a value which is a float
	identity has a value which is an int
	positive has a value which is an int
	density has a value which is an int
	pattern_from has a value which is an int
	pattern_to has a value which is an int
	query_from has a value which is an int
	query_to has a value which is an int
	query_strand has a value which is a string
	query_frame has a value which is an int
	hit_from has a value which is an int
	hit_to has a value which is an int
	hit_strand has a value which is a string
	hit_frame has a value which is an int
	align_len has a value which is an int
	gaps has a value which is an int
	qseq has a value which is a string
	hseq has a value which is a string
	midline has a value which is a string
Statistics is a reference to a hash where the following keys are defined:
	db_num has a value which is an int
	db_len has a value which is an int
	hsp_len has a value which is an int
	eff_space has a value which is an int
	kappa has a value which is a float
	lambda has a value which is a float
	entropy has a value which is a float
FeatureMetadata is a reference to a hash where the following keys are defined:
	function has a value which is a string
	genome_name has a value which is a string
	genome_id has a value which is a string
	md5 has a value which is a string
	locus_tag has a value which is a string
	alt_locus_tag has a value which is a string
	match_count has a value which is an int

=end text



=item Description


=back

=cut

sub blast_fasta_to_genomes
{
    my $self = shift;
    my($fasta_data, $program, $genomes, $subject_type, $evalue_cutoff, $max_hits, $min_coverage) = @_;

    my @_bad_arguments;
    (!ref($fasta_data)) or push(@_bad_arguments, "Invalid type for argument \"fasta_data\" (value was \"$fasta_data\")");
    (!ref($program)) or push(@_bad_arguments, "Invalid type for argument \"program\" (value was \"$program\")");
    (ref($genomes) eq 'ARRAY') or push(@_bad_arguments, "Invalid type for argument \"genomes\" (value was \"$genomes\")");
    (!ref($subject_type)) or push(@_bad_arguments, "Invalid type for argument \"subject_type\" (value was \"$subject_type\")");
    (!ref($evalue_cutoff)) or push(@_bad_arguments, "Invalid type for argument \"evalue_cutoff\" (value was \"$evalue_cutoff\")");
    (!ref($max_hits)) or push(@_bad_arguments, "Invalid type for argument \"max_hits\" (value was \"$max_hits\")");
    (!ref($min_coverage)) or push(@_bad_arguments, "Invalid type for argument \"min_coverage\" (value was \"$min_coverage\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to blast_fasta_to_genomes:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::HomologyService::Service::CallContext;
    my($reports, $metadata);
    #BEGIN blast_fasta_to_genomes

    #
    # If we have a token, temporarily override the token in our P3DataAPI instance.
    #
    
    local $self->{_data_api}->{token} = $ctx->token;
    print STDERR "Data api: " . Dumper($self->{_data_api});
    
    ($reports, $metadata) = $self->{_util}->blast_fasta_to_genomes($fasta_data, $program, $genomes, $subject_type, $evalue_cutoff, $max_hits, $min_coverage);
    
    #END blast_fasta_to_genomes
    my @_bad_returns;
    (ref($reports) eq 'ARRAY') or push(@_bad_returns, "Invalid type for return variable \"reports\" (value was \"$reports\")");
    (ref($metadata) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"metadata\" (value was \"$metadata\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to blast_fasta_to_genomes:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($reports, $metadata);
}


=head2 blast_fasta_to_taxon

  $reports, $metadata = $obj->blast_fasta_to_taxon($fasta_data, $program, $taxon_id, $subject_type, $evalue_cutoff, $max_hits, $min_coverage)

=over 4


=item Parameter and return types

=begin html

<pre>
$fasta_data is a string
$program is a string
$taxon_id is a string
$subject_type is a string
$evalue_cutoff is a float
$max_hits is an int
$min_coverage is a float
$reports is a reference to a list where each element is a Report
$metadata is a reference to a hash where the key is a string and the value is a FeatureMetadata
Report is a reference to a hash where the following keys are defined:
	program has a value which is a string
	version has a value which is a string
	reference has a value which is a string
	search_target has a value which is a reference to a hash where the following keys are defined:
	db has a value which is a string
	subjects has a value which is a string

	params has a value which is a reference to a hash where the following keys are defined:
	matrix has a value which is a string
	expect has a value which is a float
	include has a value which is a float
	sc_match has a value which is an int
	sc_mismatch has a value which is an int
	gap_open has a value which is an int
	gap_extend has a value which is an int
	filter has a value which is a string
	pattern has a value which is a string
	entrez_query has a value which is a string
	cbs has a value which is an int
	query_gencode has a value which is an int
	db_gencode has a value which is an int
	bl2seq_mode has a value which is a string

	result has a value which is a reference to a hash where the following keys are defined:
	search has a value which is a Search

Search is a reference to a hash where the following keys are defined:
	query_id has a value which is a string
	query_title has a value which is a string
	query_len has a value which is an int
	hits has a value which is a reference to a list where each element is a Hit
	stat has a value which is a Statistics
Hit is a reference to a hash where the following keys are defined:
	num has a value which is an int
	description has a value which is a reference to a list where each element is a HitDescr
	len has a value which is an int
	hsps has a value which is a reference to a list where each element is a Hsp
HitDescr is a reference to a hash where the following keys are defined:
	id has a value which is a string
	accession has a value which is a string
	title has a value which is a string
	taxid has a value which is an int
	sciname has a value which is a string
Hsp is a reference to a hash where the following keys are defined:
	num has a value which is an int
	bit_score has a value which is a float
	score has a value which is a float
	evalue has a value which is a float
	identity has a value which is an int
	positive has a value which is an int
	density has a value which is an int
	pattern_from has a value which is an int
	pattern_to has a value which is an int
	query_from has a value which is an int
	query_to has a value which is an int
	query_strand has a value which is a string
	query_frame has a value which is an int
	hit_from has a value which is an int
	hit_to has a value which is an int
	hit_strand has a value which is a string
	hit_frame has a value which is an int
	align_len has a value which is an int
	gaps has a value which is an int
	qseq has a value which is a string
	hseq has a value which is a string
	midline has a value which is a string
Statistics is a reference to a hash where the following keys are defined:
	db_num has a value which is an int
	db_len has a value which is an int
	hsp_len has a value which is an int
	eff_space has a value which is an int
	kappa has a value which is a float
	lambda has a value which is a float
	entropy has a value which is a float
FeatureMetadata is a reference to a hash where the following keys are defined:
	function has a value which is a string
	genome_name has a value which is a string
	genome_id has a value which is a string
	md5 has a value which is a string
	locus_tag has a value which is a string
	alt_locus_tag has a value which is a string
	match_count has a value which is an int
</pre>

=end html

=begin text

$fasta_data is a string
$program is a string
$taxon_id is a string
$subject_type is a string
$evalue_cutoff is a float
$max_hits is an int
$min_coverage is a float
$reports is a reference to a list where each element is a Report
$metadata is a reference to a hash where the key is a string and the value is a FeatureMetadata
Report is a reference to a hash where the following keys are defined:
	program has a value which is a string
	version has a value which is a string
	reference has a value which is a string
	search_target has a value which is a reference to a hash where the following keys are defined:
	db has a value which is a string
	subjects has a value which is a string

	params has a value which is a reference to a hash where the following keys are defined:
	matrix has a value which is a string
	expect has a value which is a float
	include has a value which is a float
	sc_match has a value which is an int
	sc_mismatch has a value which is an int
	gap_open has a value which is an int
	gap_extend has a value which is an int
	filter has a value which is a string
	pattern has a value which is a string
	entrez_query has a value which is a string
	cbs has a value which is an int
	query_gencode has a value which is an int
	db_gencode has a value which is an int
	bl2seq_mode has a value which is a string

	result has a value which is a reference to a hash where the following keys are defined:
	search has a value which is a Search

Search is a reference to a hash where the following keys are defined:
	query_id has a value which is a string
	query_title has a value which is a string
	query_len has a value which is an int
	hits has a value which is a reference to a list where each element is a Hit
	stat has a value which is a Statistics
Hit is a reference to a hash where the following keys are defined:
	num has a value which is an int
	description has a value which is a reference to a list where each element is a HitDescr
	len has a value which is an int
	hsps has a value which is a reference to a list where each element is a Hsp
HitDescr is a reference to a hash where the following keys are defined:
	id has a value which is a string
	accession has a value which is a string
	title has a value which is a string
	taxid has a value which is an int
	sciname has a value which is a string
Hsp is a reference to a hash where the following keys are defined:
	num has a value which is an int
	bit_score has a value which is a float
	score has a value which is a float
	evalue has a value which is a float
	identity has a value which is an int
	positive has a value which is an int
	density has a value which is an int
	pattern_from has a value which is an int
	pattern_to has a value which is an int
	query_from has a value which is an int
	query_to has a value which is an int
	query_strand has a value which is a string
	query_frame has a value which is an int
	hit_from has a value which is an int
	hit_to has a value which is an int
	hit_strand has a value which is a string
	hit_frame has a value which is an int
	align_len has a value which is an int
	gaps has a value which is an int
	qseq has a value which is a string
	hseq has a value which is a string
	midline has a value which is a string
Statistics is a reference to a hash where the following keys are defined:
	db_num has a value which is an int
	db_len has a value which is an int
	hsp_len has a value which is an int
	eff_space has a value which is an int
	kappa has a value which is a float
	lambda has a value which is a float
	entropy has a value which is a float
FeatureMetadata is a reference to a hash where the following keys are defined:
	function has a value which is a string
	genome_name has a value which is a string
	genome_id has a value which is a string
	md5 has a value which is a string
	locus_tag has a value which is a string
	alt_locus_tag has a value which is a string
	match_count has a value which is an int

=end text



=item Description

Post demo we will slot this in here.
                                      BlastParameters blast_parameters
=back

=cut

sub blast_fasta_to_taxon
{
    my $self = shift;
    my($fasta_data, $program, $taxon_id, $subject_type, $evalue_cutoff, $max_hits, $min_coverage) = @_;

    my @_bad_arguments;
    (!ref($fasta_data)) or push(@_bad_arguments, "Invalid type for argument \"fasta_data\" (value was \"$fasta_data\")");
    (!ref($program)) or push(@_bad_arguments, "Invalid type for argument \"program\" (value was \"$program\")");
    (!ref($taxon_id)) or push(@_bad_arguments, "Invalid type for argument \"taxon_id\" (value was \"$taxon_id\")");
    (!ref($subject_type)) or push(@_bad_arguments, "Invalid type for argument \"subject_type\" (value was \"$subject_type\")");
    (!ref($evalue_cutoff)) or push(@_bad_arguments, "Invalid type for argument \"evalue_cutoff\" (value was \"$evalue_cutoff\")");
    (!ref($max_hits)) or push(@_bad_arguments, "Invalid type for argument \"max_hits\" (value was \"$max_hits\")");
    (!ref($min_coverage)) or push(@_bad_arguments, "Invalid type for argument \"min_coverage\" (value was \"$min_coverage\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to blast_fasta_to_taxon:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::HomologyService::Service::CallContext;
    my($reports, $metadata);
    #BEGIN blast_fasta_to_taxon

    ($reports, $metadata) = $self->{_util}->blast_fasta_to_taxon($fasta_data, $program, $taxon_id, $subject_type, $evalue_cutoff, $max_hits, $min_coverage);

    #END blast_fasta_to_taxon
    my @_bad_returns;
    (ref($reports) eq 'ARRAY') or push(@_bad_returns, "Invalid type for return variable \"reports\" (value was \"$reports\")");
    (ref($metadata) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"metadata\" (value was \"$metadata\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to blast_fasta_to_taxon:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($reports, $metadata);
}


=head2 enumerate_databases

  $return = $obj->enumerate_databases()

=over 4


=item Parameter and return types

=begin html

<pre>
$return is a reference to a list where each element is a DatabaseDescription
DatabaseDescription is a reference to a hash where the following keys are defined:
	name has a value which is a string
	key has a value which is a string
	db_type has a value which is a string
	seq_count has a value which is an int
</pre>

=end html

=begin text

$return is a reference to a list where each element is a DatabaseDescription
DatabaseDescription is a reference to a hash where the following keys are defined:
	name has a value which is a string
	key has a value which is a string
	db_type has a value which is a string
	seq_count has a value which is an int

=end text



=item Description


=back

=cut

sub enumerate_databases
{
    my $self = shift;

    my $ctx = $Bio::KBase::HomologyService::Service::CallContext;
    my($return);
    #BEGIN enumerate_databases

    $return = $self->{_util}->enumerate_databases;
    
    #END enumerate_databases
    my @_bad_returns;
    (ref($return) eq 'ARRAY') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to enumerate_databases:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($return);
}


=head2 blast_fasta_to_database

  $reports, $metadata, $identical_proteins = $obj->blast_fasta_to_database($fasta_data, $program, $database_key, $evalue_cutoff, $max_hits, $min_coverage)

=over 4


=item Parameter and return types

=begin html

<pre>
$fasta_data is a string
$program is a string
$database_key is a string
$evalue_cutoff is a float
$max_hits is an int
$min_coverage is a float
$reports is a reference to a list where each element is a Report
$metadata is a reference to a hash where the key is a string and the value is a FeatureMetadata
$identical_proteins is a reference to a hash where the key is a string and the value is a reference to a list where each element is a reference to a list containing 2 items:
	0: a string
	1: a FeatureMetadata
Report is a reference to a hash where the following keys are defined:
	program has a value which is a string
	version has a value which is a string
	reference has a value which is a string
	search_target has a value which is a reference to a hash where the following keys are defined:
	db has a value which is a string
	subjects has a value which is a string

	params has a value which is a reference to a hash where the following keys are defined:
	matrix has a value which is a string
	expect has a value which is a float
	include has a value which is a float
	sc_match has a value which is an int
	sc_mismatch has a value which is an int
	gap_open has a value which is an int
	gap_extend has a value which is an int
	filter has a value which is a string
	pattern has a value which is a string
	entrez_query has a value which is a string
	cbs has a value which is an int
	query_gencode has a value which is an int
	db_gencode has a value which is an int
	bl2seq_mode has a value which is a string

	result has a value which is a reference to a hash where the following keys are defined:
	search has a value which is a Search

Search is a reference to a hash where the following keys are defined:
	query_id has a value which is a string
	query_title has a value which is a string
	query_len has a value which is an int
	hits has a value which is a reference to a list where each element is a Hit
	stat has a value which is a Statistics
Hit is a reference to a hash where the following keys are defined:
	num has a value which is an int
	description has a value which is a reference to a list where each element is a HitDescr
	len has a value which is an int
	hsps has a value which is a reference to a list where each element is a Hsp
HitDescr is a reference to a hash where the following keys are defined:
	id has a value which is a string
	accession has a value which is a string
	title has a value which is a string
	taxid has a value which is an int
	sciname has a value which is a string
Hsp is a reference to a hash where the following keys are defined:
	num has a value which is an int
	bit_score has a value which is a float
	score has a value which is a float
	evalue has a value which is a float
	identity has a value which is an int
	positive has a value which is an int
	density has a value which is an int
	pattern_from has a value which is an int
	pattern_to has a value which is an int
	query_from has a value which is an int
	query_to has a value which is an int
	query_strand has a value which is a string
	query_frame has a value which is an int
	hit_from has a value which is an int
	hit_to has a value which is an int
	hit_strand has a value which is a string
	hit_frame has a value which is an int
	align_len has a value which is an int
	gaps has a value which is an int
	qseq has a value which is a string
	hseq has a value which is a string
	midline has a value which is a string
Statistics is a reference to a hash where the following keys are defined:
	db_num has a value which is an int
	db_len has a value which is an int
	hsp_len has a value which is an int
	eff_space has a value which is an int
	kappa has a value which is a float
	lambda has a value which is a float
	entropy has a value which is a float
FeatureMetadata is a reference to a hash where the following keys are defined:
	function has a value which is a string
	genome_name has a value which is a string
	genome_id has a value which is a string
	md5 has a value which is a string
	locus_tag has a value which is a string
	alt_locus_tag has a value which is a string
	match_count has a value which is an int
</pre>

=end html

=begin text

$fasta_data is a string
$program is a string
$database_key is a string
$evalue_cutoff is a float
$max_hits is an int
$min_coverage is a float
$reports is a reference to a list where each element is a Report
$metadata is a reference to a hash where the key is a string and the value is a FeatureMetadata
$identical_proteins is a reference to a hash where the key is a string and the value is a reference to a list where each element is a reference to a list containing 2 items:
	0: a string
	1: a FeatureMetadata
Report is a reference to a hash where the following keys are defined:
	program has a value which is a string
	version has a value which is a string
	reference has a value which is a string
	search_target has a value which is a reference to a hash where the following keys are defined:
	db has a value which is a string
	subjects has a value which is a string

	params has a value which is a reference to a hash where the following keys are defined:
	matrix has a value which is a string
	expect has a value which is a float
	include has a value which is a float
	sc_match has a value which is an int
	sc_mismatch has a value which is an int
	gap_open has a value which is an int
	gap_extend has a value which is an int
	filter has a value which is a string
	pattern has a value which is a string
	entrez_query has a value which is a string
	cbs has a value which is an int
	query_gencode has a value which is an int
	db_gencode has a value which is an int
	bl2seq_mode has a value which is a string

	result has a value which is a reference to a hash where the following keys are defined:
	search has a value which is a Search

Search is a reference to a hash where the following keys are defined:
	query_id has a value which is a string
	query_title has a value which is a string
	query_len has a value which is an int
	hits has a value which is a reference to a list where each element is a Hit
	stat has a value which is a Statistics
Hit is a reference to a hash where the following keys are defined:
	num has a value which is an int
	description has a value which is a reference to a list where each element is a HitDescr
	len has a value which is an int
	hsps has a value which is a reference to a list where each element is a Hsp
HitDescr is a reference to a hash where the following keys are defined:
	id has a value which is a string
	accession has a value which is a string
	title has a value which is a string
	taxid has a value which is an int
	sciname has a value which is a string
Hsp is a reference to a hash where the following keys are defined:
	num has a value which is an int
	bit_score has a value which is a float
	score has a value which is a float
	evalue has a value which is a float
	identity has a value which is an int
	positive has a value which is an int
	density has a value which is an int
	pattern_from has a value which is an int
	pattern_to has a value which is an int
	query_from has a value which is an int
	query_to has a value which is an int
	query_strand has a value which is a string
	query_frame has a value which is an int
	hit_from has a value which is an int
	hit_to has a value which is an int
	hit_strand has a value which is a string
	hit_frame has a value which is an int
	align_len has a value which is an int
	gaps has a value which is an int
	qseq has a value which is a string
	hseq has a value which is a string
	midline has a value which is a string
Statistics is a reference to a hash where the following keys are defined:
	db_num has a value which is an int
	db_len has a value which is an int
	hsp_len has a value which is an int
	eff_space has a value which is an int
	kappa has a value which is a float
	lambda has a value which is a float
	entropy has a value which is a float
FeatureMetadata is a reference to a hash where the following keys are defined:
	function has a value which is a string
	genome_name has a value which is a string
	genome_id has a value which is a string
	md5 has a value which is a string
	locus_tag has a value which is a string
	alt_locus_tag has a value which is a string
	match_count has a value which is an int

=end text



=item Description


=back

=cut

sub blast_fasta_to_database
{
    my $self = shift;
    my($fasta_data, $program, $database_key, $evalue_cutoff, $max_hits, $min_coverage) = @_;

    my @_bad_arguments;
    (!ref($fasta_data)) or push(@_bad_arguments, "Invalid type for argument \"fasta_data\" (value was \"$fasta_data\")");
    (!ref($program)) or push(@_bad_arguments, "Invalid type for argument \"program\" (value was \"$program\")");
    (!ref($database_key)) or push(@_bad_arguments, "Invalid type for argument \"database_key\" (value was \"$database_key\")");
    (!ref($evalue_cutoff)) or push(@_bad_arguments, "Invalid type for argument \"evalue_cutoff\" (value was \"$evalue_cutoff\")");
    (!ref($max_hits)) or push(@_bad_arguments, "Invalid type for argument \"max_hits\" (value was \"$max_hits\")");
    (!ref($min_coverage)) or push(@_bad_arguments, "Invalid type for argument \"min_coverage\" (value was \"$min_coverage\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to blast_fasta_to_database:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	die $msg;
    }

    my $ctx = $Bio::KBase::HomologyService::Service::CallContext;
    my($reports, $metadata, $identical_proteins);
    #BEGIN blast_fasta_to_database

    ($reports, $metadata, $identical_proteins) = $self->{_util}->blast_fasta_to_database($fasta_data, $program, $database_key, $evalue_cutoff,
											 $max_hits, $min_coverage);
    #END blast_fasta_to_database
    my @_bad_returns;
    (ref($reports) eq 'ARRAY') or push(@_bad_returns, "Invalid type for return variable \"reports\" (value was \"$reports\")");
    (ref($metadata) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"metadata\" (value was \"$metadata\")");
    (ref($identical_proteins) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"identical_proteins\" (value was \"$identical_proteins\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to blast_fasta_to_database:\n" . join("", map { "\t$_\n" } @_bad_returns);
	die $msg;
    }
    return($reports, $metadata, $identical_proteins);
}





=head2 version 

  $return = $obj->version()

=over 4

=item Parameter and return types

=begin html

<pre>
$return is a string
</pre>

=end html

=begin text

$return is a string

=end text

=item Description

Return the module version. This is a Semantic Versioning number.

=back

=cut

sub version {
    return $VERSION;
}



=head1 TYPES



=head2 HitDescr

=over 4


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
id has a value which is a string
accession has a value which is a string
title has a value which is a string
taxid has a value which is an int
sciname has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
id has a value which is a string
accession has a value which is a string
title has a value which is a string
taxid has a value which is an int
sciname has a value which is a string


=end text

=back



=head2 Hsp

=over 4


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
num has a value which is an int
bit_score has a value which is a float
score has a value which is a float
evalue has a value which is a float
identity has a value which is an int
positive has a value which is an int
density has a value which is an int
pattern_from has a value which is an int
pattern_to has a value which is an int
query_from has a value which is an int
query_to has a value which is an int
query_strand has a value which is a string
query_frame has a value which is an int
hit_from has a value which is an int
hit_to has a value which is an int
hit_strand has a value which is a string
hit_frame has a value which is an int
align_len has a value which is an int
gaps has a value which is an int
qseq has a value which is a string
hseq has a value which is a string
midline has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
num has a value which is an int
bit_score has a value which is a float
score has a value which is a float
evalue has a value which is a float
identity has a value which is an int
positive has a value which is an int
density has a value which is an int
pattern_from has a value which is an int
pattern_to has a value which is an int
query_from has a value which is an int
query_to has a value which is an int
query_strand has a value which is a string
query_frame has a value which is an int
hit_from has a value which is an int
hit_to has a value which is an int
hit_strand has a value which is a string
hit_frame has a value which is an int
align_len has a value which is an int
gaps has a value which is an int
qseq has a value which is a string
hseq has a value which is a string
midline has a value which is a string


=end text

=back



=head2 Hit

=over 4


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
num has a value which is an int
description has a value which is a reference to a list where each element is a HitDescr
len has a value which is an int
hsps has a value which is a reference to a list where each element is a Hsp

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
num has a value which is an int
description has a value which is a reference to a list where each element is a HitDescr
len has a value which is an int
hsps has a value which is a reference to a list where each element is a Hsp


=end text

=back



=head2 Statistics

=over 4


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
db_num has a value which is an int
db_len has a value which is an int
hsp_len has a value which is an int
eff_space has a value which is an int
kappa has a value which is a float
lambda has a value which is a float
entropy has a value which is a float

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
db_num has a value which is an int
db_len has a value which is an int
hsp_len has a value which is an int
eff_space has a value which is an int
kappa has a value which is a float
lambda has a value which is a float
entropy has a value which is a float


=end text

=back



=head2 Search

=over 4


=item Description

need: query-masking

=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
query_id has a value which is a string
query_title has a value which is a string
query_len has a value which is an int
hits has a value which is a reference to a list where each element is a Hit
stat has a value which is a Statistics

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
query_id has a value which is a string
query_title has a value which is a string
query_len has a value which is an int
hits has a value which is a reference to a list where each element is a Hit
stat has a value which is a Statistics


=end text

=back



=head2 Report

=over 4


=item Description

structure {
            } bl2seq;

=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
program has a value which is a string
version has a value which is a string
reference has a value which is a string
search_target has a value which is a reference to a hash where the following keys are defined:
db has a value which is a string
subjects has a value which is a string

params has a value which is a reference to a hash where the following keys are defined:
matrix has a value which is a string
expect has a value which is a float
include has a value which is a float
sc_match has a value which is an int
sc_mismatch has a value which is an int
gap_open has a value which is an int
gap_extend has a value which is an int
filter has a value which is a string
pattern has a value which is a string
entrez_query has a value which is a string
cbs has a value which is an int
query_gencode has a value which is an int
db_gencode has a value which is an int
bl2seq_mode has a value which is a string

result has a value which is a reference to a hash where the following keys are defined:
search has a value which is a Search


</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
program has a value which is a string
version has a value which is a string
reference has a value which is a string
search_target has a value which is a reference to a hash where the following keys are defined:
db has a value which is a string
subjects has a value which is a string

params has a value which is a reference to a hash where the following keys are defined:
matrix has a value which is a string
expect has a value which is a float
include has a value which is a float
sc_match has a value which is an int
sc_mismatch has a value which is an int
gap_open has a value which is an int
gap_extend has a value which is an int
filter has a value which is a string
pattern has a value which is a string
entrez_query has a value which is a string
cbs has a value which is an int
query_gencode has a value which is an int
db_gencode has a value which is an int
bl2seq_mode has a value which is a string

result has a value which is a reference to a hash where the following keys are defined:
search has a value which is a Search



=end text

=back



=head2 genome_id

=over 4


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 FeatureMetadata

=over 4


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
function has a value which is a string
genome_name has a value which is a string
genome_id has a value which is a string
md5 has a value which is a string
locus_tag has a value which is a string
alt_locus_tag has a value which is a string
match_count has a value which is an int

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
function has a value which is a string
genome_name has a value which is a string
genome_id has a value which is a string
md5 has a value which is a string
locus_tag has a value which is a string
alt_locus_tag has a value which is a string
match_count has a value which is an int


=end text

=back



=head2 BlastParameters

=over 4


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
evalue_cutoff has a value which is a float
max_hits has a value which is an int
min_coverage has a value which is a float

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
evalue_cutoff has a value which is a float
max_hits has a value which is an int
min_coverage has a value which is a float


=end text

=back



=head2 DatabaseDescription

=over 4


=item Description

db_type is either "dna" or "protein"

=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
name has a value which is a string
key has a value which is a string
db_type has a value which is a string
seq_count has a value which is an int

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
name has a value which is a string
key has a value which is a string
db_type has a value which is a string
seq_count has a value which is an int


=end text

=back


=cut

1;
