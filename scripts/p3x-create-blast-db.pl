
=head1 NAME

p3x-create-blast-db

=head1 SYNOPSIS

p3x-create-blast-db [params] blast-db-file

=head1 DESCRIPTION

Create a blast database. 

If --reference and/or --representative is selected, limit to only those genomes.

If --taxon is selected, limit to the genomes at or below that level of the taxonomy.

=cut

use strict;
use gjoseqlib;
use Getopt::Long::Descriptive;
use P3DataAPI;
use Data::Dumper;
use LPTScheduler;
use JSON::XS;
use File::Temp;
use File::Path qw(make_path remove_tree);
use IPC::Run qw(run);
use PerlIO::via::Blockwise;
use BerkeleyDB;

#
# Import this to properly intialize PATH for blast binaries.
#
use Bio::P3::HomologySearch::HomologySearch;

my($opt, $usage) = describe_options("%c %o dbtype ftype blast-db-file",
				    ['exclude-taxon=i@', "Do not include genomes in this taxon", { default => [] }],
				    ["reference", "Include reference genomes"],
				    ["representative", "Include representative genomes"],
				    ["create-nr", "Create nonredundant database"],
				    ["quality-check!", "Require quality of Good (invert with --no-quality-check)", { default => 1 }],
				    ["check-files!", "Check for download files (invert with --no-check-files)", { default => 1 }],
				    ["viral", "Create viral database. Causes CDS & mat_peptide features to be included"],
				    ["complete", "Include only genome_status=Complete genomes"],
				    ["blastdb-version=i", "Use this blastdb version", { default => 5 }],
				    ['taxon=i@', "Limit to this taxon", { default => [] }],
				    ["title=s", "Database title", { default => "blast db"}],
				    ["parallel|j=i", "Use this many threads", { default => 1 }],
				    ["batch-size=i", "Batch size for genome data lookups for data api", { default => 10 } ],
				    ["overwrite|f", "Overwrite existing database"],
				    ["dump", "Don't build, just dump the list of genomes"],
				    ["help|h", "Show this help message"]);

print($usage->text), exit 0 if $opt->help;
die($usage->text)  unless @ARGV == 3 || ($opt->dump && @ARGV == 0);

my $genome_base = "/vol/patric3/downloads/genomes";

my %type_suffix = (aa => "faa",
		   dna => "fna");
my %blast_type = (aa => "prot",
		  dna => "nucl");
my %blastdb_suffixes = (aa => ['pal', 'psq'],
			dna => ['nal', 'nsq']);

my $dbtype = shift;
my $ftype = shift;
my $db = shift;

if (!$opt->dump)
{
    if ($ftype ne 'features' && $ftype ne 'contigs')
    {
	die "Invalid ftype $ftype: must be features or contigs\n";
    }
    
    if ($dbtype ne 'aa' && $dbtype ne 'dna')
    {
	die "Invalid dbtype $dbtype: must be aa or dna\n";
    }
}

my $dbfile = "$db.$ftype.$type_suffix{$dbtype}";
my $taxids = "$dbfile.taxids";
my $taxlist = "$dbfile.taxlist";
my $glist = "$dbfile.glist";
my $nr_dir;
$nr_dir = "$db.$ftype.$type_suffix{$dbtype}.nrdb" if $opt->create_nr;

#
# Check for existing data.
#
my @blast_files = grep { -s $_ } map { "$dbfile.$_" } @{$blastdb_suffixes{$dbtype}};

if (-s $taxids &&
    @blast_files &&
    !$opt->overwrite)
{
    print STDERR "Blast data in $dbfile already exists, skipping build\n";
    exit 0;
}
print STDERR "Will build $dbfile\n";

my $api = P3DataAPI->new();

my @params;

my @rr;
push(@rr, "reference_genome:Reference") if $opt->reference;
push(@rr, "reference_genome:Representative") if $opt->representative;

#
# Require lineage.
#
push(@params, fq => 'taxon_lineage_ids:*');

for my $not (@{$opt->exclude_taxon})
{
    push(@params, fq => "-taxon_lineage_ids:$not");
}

push @params, fq => (join(" OR ", @rr)) if @rr;

if (@{$opt->taxon})
{
    push(@params, fq => "taxon_lineage_ids:" . join(" OR ",  @{$opt->taxon}));
}


push(@params, fq => "genome_status:Complete") if $opt->complete;

push(@params, fl => "genome_id,genome_name,taxon_id,genome_length,genus");

#
# Require good genomes if desired.
#
if ($opt->quality_check)
{
    push(@params, fq => 'genome_quality:Good');
}

push(@params, fq => 'genome_length:*');

push(@params, fq => "public:true");

push(@params, q => "*:*");

my $genomes = $api->solr_query_list("genome", \@params);

if ($opt->dump)
{
    print JSON::XS->new->pretty->canonical->encode($genomes);
    exit 0;
}

print STDERR "Building database with " . scalar(@$genomes) . " genomes\n";

if (@$genomes == 0)
{
    print STDERR "Skipping build of empty list\n";
    exit 0;
}

#open(DB, ">", $dbfile) or die "Cannot write $dbfile: $!";
#open(TI, ">", $taxids) or die "Cannot write $taxids: $!";

my $db_files;
my $cleanup = sub {};

if ($opt->check_files)
{
    ($db_files, $cleanup) = process_from_download_files();
}
else
{
    $db_files = process_from_data_api();
}

#
# Check for empty data.
#
if (! -s $taxids)
{
    print STDERR "Skipping creation of empty database for $dbfile\n";
    &$cleanup();
    exit 0;
}

my $ok = run(["cat", @$db_files],
	     "|",
	     ["makeblastdb",
	      "-in", "-",
	      "-out", $dbfile,
	      "-parse_seqids",
	      "-taxid_map", $taxids,
	      "-title", $opt->title,
	      "-blastdb_version", $opt->blastdb_version,
	      "-dbtype", $blast_type{$dbtype}]);

&$cleanup();

$ok or die "Error  creating blastdb: $?";
	   
#
# We write a single data file with results from the batched data api calls.
#

sub process_from_data_api
{
    my $dbfile = File::Temp->new(UNLINK => 0);
    print "create $dbfile\n";
    my $api = P3DataAPI->new();
    my $tax;
    open($tax, ">", $taxids) or die "cannot write $taxids: $!";
    open(GL, ">", $glist) or die "Cannot write $glist: $!";
    print GL "$_->{genome_id}\n" foreach sort {
	my($a1, $a2) = split(/\./, $a->{genome_id});
	my($b1, $b2) = split(/\./, $b->{genome_id});
	$a1 <=> $b1 or $a2 <=> $b2;
    } @$genomes;
    close(GL);

    my @todo = @$genomes;

    my $id_map = {};
    my $data_cb;
    my $db_env;
    my $tbl_tax_to_pseudo;
    my $tbl_pseudo_to_tax;
    my $tbl_md5_to_feature;
    if ($opt->create_nr)
    {
	$data_cb = \&data_callback_nr;
	if (-d $nr_dir)
	{
	    if ($opt->overwrite)
	    {
		remove_tree($nr_dir);
	    }
	    else
	    {
		die "NR dir $nr_dir already exists\n";
	    }
	}
	make_path($nr_dir);
	
	$db_env =  BerkeleyDB::Env->new(-Home => $nr_dir,
					-ErrFile => *STDERR,
					-Flags => DB_CREATE | DB_INIT_CDB | DB_INIT_MPOOL);
	$db_env or die "error creating NR db: $!";
	$tbl_tax_to_pseudo = BerkeleyDB::Btree->new(-Filename, "tax_to_pseudo",
						    -Env => $db_env,
						    -Flags => DB_CREATE,
						    -Property => DB_DUP);
	$tbl_pseudo_to_tax = BerkeleyDB::Btree->new(-Filename, "pseudo_to_tax",
						    -Env => $db_env,
						    -Flags => DB_CREATE,
						    -Property => DB_DUP);
	$tbl_md5_to_feature = BerkeleyDB::Btree->new(-Filename, "md5_to_feature",
						     -Env => $db_env,
						     -Flags => DB_CREATE,
						     -Property => DB_DUP);
    }
    else
    {
	$data_cb = \&data_callback;
    }

    my $ftype_cond;
    if ($opt->viral)
    {
	$ftype_cond = [ "in",  "feature_type", "(mat_peptide,CDS)" ];
    }
    else
    {
	$ftype_cond = [ "eq",  "feature_type", "CDS" ];
    }

    while (@todo)
    {
	my @batch = splice(@todo, 0, $opt->batch_size);

	my @genomes = map { $_->{genome_id} } @batch;

	# print "Batch @genomes\n";

	my $genome_cond = [ "in","genome_id", "(" . join(",", @genomes) . ")"];
	
	if ($ftype eq 'features')
	{
	    my $key;
	    if ($dbtype eq 'aa')
	    {
		$key = "aa_sequence_md5";
	    }		    
	    else
	    {
		$key = "na_sequence_md5";
	    }
	    $api->query_cb("genome_feature",
			   sub { $data_cb->($api, $dbfile, $tax, $key, $id_map, @_) },
			   
			   $ftype_cond,
			   [ "eq", "annotation", "PATRIC"],
			   $genome_cond,
			   [ "select", "patric_id,product,genome_id,taxon_id,$key" ]);

	}
	else
	{
	    $api->query_cb("genome_sequence",
			    sub {
				my ($data) = @_;
				for my $ent (@$data) {
				    my $id = "lcl|$ent->{sequence_id}";
				    print_alignment_as_fasta($dbfile,
							     [$id,
							      "$ent->{description} [ $ent->{genome_name} | $ent->{genome_id} ]",
							      $ent->{sequence}]);
				    print $tax "$id\t$ent->{taxon_id}\n";
				}
				return 1;
			    },
			   $genome_cond,
			   ["select", "accession,genome_id,taxon_id,description,genome_name,sequence,sequence_id"]);

	}
    }
    # print Dumper($id_map);
    open(TL, ">", $taxlist) or die "Cannot write $taxlist: $!";

    my %taxa_seen;
    if ($opt->create_nr && $ftype eq 'features')
    {
	#
	# create faux taxids for each unique set of taxids that appear in the id map.
	#
	my %mapid;
	my $nextid = 1;
	for my $md5 (sort keys %$id_map)
	{
	    my $ids = $id_map->{$md5};
	    my %tset;
	    
	    for my $ent (@$ids)
	    {
		my($fid, $tax) = @$ent;
		$tbl_md5_to_feature->db_put($md5, $fid);
		$taxa_seen{$tax} = $tset{$tax} = 1;
	    }

	    my @tlist = sort { $a <=> $b } keys %tset;
	    my $tid = join(",", @tlist);
	    my $mid = $mapid{$tid};
	    if (!$mid)
	    {
		$mid = $nextid++;
		$mapid{$tid} = $mid;

		$tbl_tax_to_pseudo->db_put($_, $mid) foreach @tlist;
		$tbl_pseudo_to_tax->db_put($mid, $_) foreach @tlist;
	    }
	    print $tax "lcl|$md5\t$mid\n";
	}
	# die Dumper(\%mapid,$taxlist);

	#
	# populate md5 table from id map
	#
    }
    else
    {
	my %mapid;
	for my $ids (values %$id_map)
	{
	    $taxa_seen{$_->[1]} = 1 foreach @$ids;
	}
    }
    print TL "$_\n" foreach sort { $a <=> $b } keys %taxa_seen;
    close(TL);
	 
    close($dbfile);
    close($tax);
    return [$dbfile];
}

sub data_callback
{
    my ($api, $fh, $tax, $key, $id_map, $data) = @_;
    my %by_md5;

    $api->lookup_sequence_data([map { $_->{$key} } @$data ], sub {
	my $ent = shift;
	$by_md5{$ent->{md5}} = $ent->{sequence};
    });
    
    for my $ent (@$data) {
	my $id = "gnl|$ent->{patric_id}";
	print_alignment_as_fasta($fh,
				 [
				  $id, $ent->{product},
				  $by_md5{$ent->{$key}}
				  ]
				);
	print $tax "$id\t$ent->{taxon_id}\n";
	push(@{$id_map->{$id}}, [$id, $ent->{taxon_id}]);
    }
    return 1;
}

sub data_callback_nr
{
    my ($api, $fh, $tax, $key, $id_map, $data) = @_;

    #
    # Check the idmap for any MD5s we haven't looked up sequence for yet. For those,
    # we will do a sequence lookup and write the data, keyed with the md5, to
    # our output file.
    #

    my %to_find = map { $_->{$key} => 1 } grep { !exists $id_map->{$_->{$key}} } @$data;

    my $nd = @$data;
    my $tf = keys %to_find;
    # print "Data size $nd looking up $tf\n";
   
    $api->lookup_sequence_data([keys %to_find ], sub {
	my $ent = shift;
	print_alignment_as_fasta($fh, ["lcl|$ent->{md5}", '', $ent->{sequence}]);
    });

    #
    # Update the idmap for all of the sequences.
    #
    
    for my $ent (@$data) {
	my $md5 = $ent->{$key};
	
	push(@{$id_map->{$md5}}, [$ent->{patric_id}, $ent->{taxon_id}]);
    }
    return 1;
}

sub process_from_download_files
{
    my $sched = LPTScheduler->new($opt->parallel);
    
    open(GL, ">", $glist) or die "Cannot write $glist: $!";
    my %tax;
    
    for my $g (@$genomes)
    {
	my($genome_id, $genome_name, $taxon_id, $genome_length) = @$g{qw(genome_id genome_name taxon_id genome_length)};
	my $path;
	if ($opt->check_files)
	{
	    $path = compute_path($genome_id, $ftype, $dbtype);
	    # next unless $path;
	}
	
	$g->{path} = $path;
	$sched->add_work($g, $genome_length);
	print GL "$genome_id\n";
	$tax{$taxon_id} = 1;
    }
    close(GL);
    open(TL, ">", $taxlist) or die "Cannot write $taxlist: $!";
    print TL "$_\n" foreach sort { $a <=> $b  } keys %tax;
    close(TL);

    #
    # Don't clean up. The forked processes will inherit the ref and
    # undef it at completion, triggering cleanup.
    #
    my $tmpdir = File::Temp->newdir(CLEANUP => 0);
    my $dbdir = "$tmpdir/db";
    my $taxdir = "$tmpdir/tax";
    
    make_path($dbdir, $taxdir);

    my $cleanup = sub {
	my $err;
	remove_tree("$tmpdir", { error => \$err });
	if ($err && @$err)
	{
	    for my $diag (@$err) {
		my ($file, $message) = %$diag;
		if ($file eq '') {
		    warn "general error: $message\n";
		}
		else {
		    warn "problem unlinking $file: $message\n";
		}
	    }
	}
    };

    my $boot = sub {
	my $dbfile = "$dbdir/$$";
	my $taxfile = "$taxdir/$$";
	my $db = IO::File->new($dbfile, "w");
	$db or die "Cannot write $dbfile: $!";
	my $tax = IO::File->new($taxfile, "w");
	$tax or die "Cannot write $taxfile: $!";
	return [$db, $tax, $api];
    };
    
    $sched->run($boot, sub {
	my($global, $work) = @_;
	my($db, $tax, $api) = @$global;
	
	my($path, $genome_id, $genome_name, $taxon_id, $genome_length) = @$work{qw(path genome_id genome_name taxon_id genome_length)};
	
	if (!$path || !open(P, "<:via(Blockwise)", $path))
	{
	    print STDERR "$path open failed, using data api\n";
	    #
	    # Load from data api.
	    #
	    
	    if ($ftype eq 'features')
	    {
		if ($dbtype eq 'aa')
		{
		    print STDERR "Loading $genome_id from api\n";
		    $path = $api->retrieve_protein_features_in_genomes_to_temp([$genome_id]);
		}
		else
		{
		    print STDERR "Loading $genome_id from api\n";
		    $path = $api->retrieve_dna_features_in_genomes_to_temp([$genome_id]);
		}
	    }
	    else		# contigs
	    {
		print STDERR "Loading $genome_id from api\n";
		$path = $api->retrieve_contigs_in_genomes_to_temp([$genome_id]);
	    }
	    
	    if (!open(P, "<:via(Blockwise)", $path))
	    {
		warn "Cannot open $path from data api build from $genome_id: $!\n";
		return;
	    }
	}
	my $skip;
	my %seen;
	if ($ftype eq 'features')
	{
	    while (my($rawid, $def, $seq) = read_next_fasta(\*P))
	    {
		if ($rawid =~ /^(fig\|\d+\.\d+\.[^.]+\.\d+)/)
		{
		    my $id = "gnl|$1";
		    
		    if ($seen{$id}++)
		    {
			warn "Skipping duplicate id $id\n";
		    }
		    else
		    {
			if ($seq !~ /^[a-z*]+$/i)
			{
			    warn "Bad sequence $id\n";
			}
			else
			{
			    print_alignment_as_fasta($db, [$id, undef, $seq]);
			    print $tax "$id\t$taxon_id\n";
			}
		    }
		}
	    }
	}
	else
	{
	    while (my($rawid, $def, $seq) = read_next_fasta(\*P))
	    {
		if ($rawid =~ /^(accn\|)?(\S+)/)
		{
		    my $id = "lcl|${genome_id}.$2";
		    
		    if ($seen{$id}++)
		    {
			warn "Skipping duplicate id $id\n";
		    }
		    else
		    {
			if ($seq !~ /^[a-z*]+$/i)
			{
			    warn "Bad sequence $id\n";
			}
			else
			{
			    print_alignment_as_fasta($db, [$id, undef, $seq]);
			    print $tax "$id\t$taxon_id\n";
			}
		    }
		}
	    }
	}
    });

    my @taxf = glob("$taxdir/*");
    my @dbf = glob("$dbdir/*");

    if (@taxf == 0)
    {
	warn "No taxids found for $dbfile\n";
	&$cleanup();
	exit 0;
    }

    my $ok = run(["cat", @taxf], ">", $taxids);

    if (!$ok)
    {
	&$cleanup();
	die "Failure creating $taxids from @taxf\n";
    }

    return(\@dbf, $cleanup);
}

sub compute_path
{
    my($genome, $ftype, $dbtype) = @_;
    
    my $dir = "$genome_base/$genome";

    # if (! -d $dir)
    # {
    # 	warn "No directory $dir\n";
    # 	return undef;
    # }

    if ($ftype eq 'contigs')
    {
	if ($dbtype eq 'aa')
	{
	    die "Invalid request: $ftype and $dbtype\n";
	}

	my $path = "$dir/$genome.fna";
	if (! -f $path)
	{
	    warn "Empty data for $path\n";
	    return undef;
	}
	# print STDERR "Contig path for $genome is $path\n";
	return $path;
    }
    else
    {
	my $path;
	if ($dbtype eq 'aa')
	{	
	    $path = "$dir/$genome.PATRIC.faa";
	}
	else
	{
	    $path = "$dir/$genome.PATRIC.ffn";
	}

	# if (! -f $path)
	# {
	#     warn "Empty data for $path\n";
	#     return undef;
	# }
	return $path;
    }
      
}

