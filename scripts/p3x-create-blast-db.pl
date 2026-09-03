
=head1 NAME

p3x-create-blast-db

=head1 SYNOPSIS

p3x-create-blast-db [params] blast-db-file

=head1 DESCRIPTION

Create a blast database. 

If --reference and/or --representative is selected, limit to only those genomes.

If --taxon is selected, limit to the genomes at or below that level of the taxonomy.

Sequence data comes from the per-genome files under the downloads tree when
they are present, and from the data API otherwise. Bacterial builds are mostly
cache hits (roughly 80%, though coverage is uneven across the genome id range);
the misses are fetched in batches, not one genome at a time.

Viral builds must pass --no-check-files. No viral genome is written to the
download site, so every cache probe would miss, and the probe itself is a stat
per genome against a busy shared filesystem.

=cut

use strict;
use gjoseqlib;
use Getopt::Long::Descriptive;
use P3DataAPI;
use Data::Dumper;
use LPTScheduler;
use JSON::XS;
use File::Temp;
use IO::File;
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
				    ['refrep-only-taxon=i@', "Within this taxon include only reference and representative genomes; genomes outside it are unaffected. Repeatable. Use for taxa whose clinical-surveillance depositions would otherwise swamp the database (SARS-CoV-2 in Coronaviridae)", { default => [] }],
				    ["reference", "Include reference genomes"],
				    ["representative", "Include representative genomes"],
				    ["create-nr", "Create nonredundant database"],
				    ["quality-check!", "Require quality of Good (invert with --no-quality-check)", { default => 1 }],
				    ["check-files!", "Build from the per-genome files under the downloads tree where they exist, fetching only the rest from the data API. Invert with --no-check-files to read everything from the data API -- which is what viral builds must do, as viral genomes are not written to the download site", { default => 1 }],
				    ["viral", "Create viral database. Causes CDS & mat_peptide features to be included"],
				    ["complete", "Include only genome_status=Complete genomes"],
				    ['taxon-filter=s%' => "Define a taxon filter file", { default => {} }],
				    ["blastdb-version=i", "Use this blastdb version", { default => 5 }],
				    ['taxon=i@', "Limit to this taxon", { default => [] }],
				    ["title=s", "Database title", { default => "blast db"}],
				    ["parallel|j=i", "Use this many threads", { default => 1 }],
				    ["batch-size=i", "Batch size for genome data lookups for data api", { default => 10 } ],
				    ["overwrite|f", "Overwrite existing database"],
				    ["max-missing-fraction=f", "Tolerated mismatch between database records and .taxids lines. Exact match is the norm for bacterial builds; viral sets need a few percent", { default => 0 }],
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

my %taxon_filter;
while (my($t, $f) = each(%{$opt->taxon_filter}))
{
    print "Loading taxon filter $f for $t\n";
    open(F, "<", $f) or die "cannot open taxon filter file $f for taxon $t: $!";
    while (<F>)
    {
	if (/(\d+\.\d+)/)
	{
	    $taxon_filter{$t}->{$1} = 1;
	}
    }
    close(F);
}


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

#
# Per-taxon reference/representative carve-out.
#
# --reference / --representative above are global: they constrain every genome
# in the query. This is the narrow form -- drop the genomes that are inside the
# named taxon AND are not reference or representative, leaving everything
# outside the taxon untouched. Written as a single negated clause rather than
# "-taxon OR refrep" because a bare negative as the left arm of an OR does not
# reliably match in Solr.
#
# The motivating case is Coronaviridae: 9,413,387 of its 9,487,138 genomes
# (99.2%) are SARS-CoV-2 clinical samples, of which exactly 22 are reference or
# representative. Without this the family enumerates ~9.5M genomes -- which the
# data API does not survive -- and would build an estimated 150 GB.
#
for my $rr_taxon (@{$opt->refrep_only_taxon})
{
    push(@params, fq => "-(taxon_lineage_ids:$rr_taxon -reference_genome:(Reference OR Representative))");
}

push @params, fq => (join(" OR ", @rr)) if @rr;

if (@{$opt->taxon})
{
    #
    # The parentheses are load bearing. "taxon_lineage_ids:A OR B" qualifies
    # only A; the bare B is resolved against Solr's default field, so every
    # --taxon after the first selected a different set than the one asked for
    # -- not a subset, a wrong set. Measured against the ref/ invocation
    # "--taxon 2 --taxon 2157": the correct union is 1,428,217 genomes
    # (1,390,544 Bacteria + 37,673 Archaea) and the unparenthesized form
    # returns 1,429,101, i.e. 884 genomes that match the string 2157 somewhere
    # rather than carrying it in their lineage. It went unnoticed because the
    # driver only ever passed one taxon -- until it had to pass two for the
    # genus names that denote more than one taxon.
    #
    push(@params, fq => "taxon_lineage_ids:(" . join(" OR ",  @{$opt->taxon}) . ")");
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


my @new;
for my $gid (@$genomes)
{
    if (my $tf = $taxon_filter{$gid->{taxon_id}})
    {
	next unless $tf->{$gid->{genome_id}};
    }
    push(@new, $gid);
}
@$genomes = @new;

#
# Guard against the same genome appearing twice.
#
# This was the shape of the 2026-08 corruption: deep paging in the data API
# returned some genomes twice and silently dropped an equal number, and a
# genome fetched twice produces duplicate FASTA seqids, which makeblastdb
# -parse_seqids rejects outright -- so the whole database fails to build, three
# stages downstream of the actual fault. P3DataAPI now uses cursor paging and
# should never hand us a duplicate, but this is the one point where the entire
# list is in hand, and the check is far cheaper than the failure it prevents.
#
{
    my %seen;
    my @uniq = grep { !$seen{$_->{genome_id}}++ } @$genomes;
    if (@uniq != @$genomes)
    {
	my @dups = sort grep { $seen{$_} > 1 } keys %seen;
	warn "Dropping " . (scalar(@$genomes) - scalar(@uniq)) .
	     " duplicate genome(s) from the query result: @dups\n";
	@$genomes = @uniq;
    }
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

verify_record_count($dbfile, $taxids);

#
# makeblastdb accepts a truncated input stream without complaint: if one of the
# child fasta files was cut short (out of space, I/O error, a child that died)
# we still get a database, just a smaller one. Nothing downstream notices --
# the catalog loader trusts the .taxids sidecar and never reads the database.
#
# The .pjs/.njs metadata records how many sequences actually landed, and
# .taxids has one line per sequence we intended to write. Here is the only
# point where both numbers are known, so compare them here.
#
# A small shortfall is normal: makeblastdb deterministically rejects a few
# records in some viral sets. A large one means truncation, so fail and remove
# .taxids -- that both stops the skip guard at the top of this script from
# accepting the database on a rerun and keeps p3x-create-databases-lookup from
# cataloging it.
#
sub verify_record_count
{
    my($dbfile, $taxids) = @_;

    my($meta) = grep { -f $_ } map { "$dbfile.$_" } qw(pjs njs);
    if (!$meta)
    {
	warn "Cannot verify $dbfile: no .pjs or .njs metadata file\n";
	return;
    }

    my $built;
    if (open(my $fh, "<", $meta))
    {
	local $/;
	my $txt = <$fh>;
	close($fh);
	$built = eval { decode_json($txt)->{"number-of-sequences"} };
    }
    if (!defined $built)
    {
	warn "Cannot verify $dbfile: no number-of-sequences in $meta\n";
	return;
    }

    my $expected = 0;
    if (open(my $fh, "<", $taxids))
    {
	$expected++ while <$fh>;
	close($fh);
    }
    else
    {
	warn "Cannot verify $dbfile: $taxids: $!\n";
	return;
    }

    return if $built == $expected;

    #
    # Either direction is a defect. A shortfall is truncated input; a surplus
    # means the taxid map does not describe what is in the database, which
    # breaks -taxidlist searches just as badly. Measured on data.2022-0916:
    # every clean bacterial database matches exactly, while the damaged ones
    # sit at 1.6% and 5.0% short and one runs 1270 records long -- so the
    # useful default is exact equality, not a percentage band.
    #
    my $delta = $expected - $built;
    my $frac = $expected ? abs($delta) / $expected : 1;

    if ($frac > $opt->max_missing_fraction)
    {
	unlink($taxids);
	die sprintf("Database %s does not match its taxid map: %d records built, " .
		    "%d expected (%d %s, %.3f%%, tolerance %.3f%%). Removed %s so " .
		    "this database is not reused or cataloged.\n",
		    $dbfile, $built, $expected, abs($delta),
		    $delta > 0 ? "missing" : "extra",
		    100 * $frac, 100 * $opt->max_missing_fraction, $taxids);
    }

    warn sprintf("Note: %s holds %d of %d records (%d %s, within tolerance)\n",
		 $dbfile, $built, $expected, abs($delta),
		 $delta > 0 ? "missing" : "extra");
}
	   
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
	print $tax "$id\t$ent->{taxon_id}\n" or die "Cannot write to tax file: $!";
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

    #
    # Don't clean up. The forked processes will inherit the ref and
    # undef it at completion, triggering cleanup.
    #
    my $tmpdir = File::Temp->newdir(CLEANUP => 0);
    my $dbdir = "$tmpdir/db";
    my $taxdir = "$tmpdir/tax";
    my $fetchdir = "$tmpdir/fetch";

    make_path($dbdir, $taxdir, $fetchdir);

    open(GL, ">", $glist) or die "Cannot write $glist: $!";
    my %tax;
    my @missing;

    for my $g (@$genomes)
    {
	my($genome_id, $taxon_id) = @$g{qw(genome_id taxon_id)};

	#
	# compute_path only checks for the file in the contigs case (the
	# feature branch has its -f test commented out), so test here rather
	# than relying on the worker's open to fail. A zero-length file counts
	# as a miss: it would contribute nothing and hide the genome.
	#
	$g->{path} = compute_path($genome_id, $ftype, $dbtype);
	push(@missing, $g) unless ($g->{path} && -s $g->{path});

	print GL "$genome_id\n";
	$tax{$taxon_id} = 1;
    }
    close(GL);

    my $fetched = prefetch_missing_genomes($api, \@missing, $fetchdir);
    $_->{path} = $fetched->{$_->{genome_id}} for @missing;

    $sched->add_work($_, $_->{genome_length}) for @$genomes;

    open(TL, ">", $taxlist) or die "Cannot write $taxlist: $!";
    print TL "$_\n" foreach sort { $a <=> $b  } keys %tax;
    close(TL);

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
			    $db->error
				and die "Write failed on $dbdir/$$ after $id: $!\n";
			    print $tax "$id\t$taxon_id\n"
				or die "Write failed on $taxdir/$$ after $id: $!\n";
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
			    $db->error
				and die "Write failed on $dbdir/$$ after $id: $!\n";
			    print $tax "$id\t$taxon_id\n"
				or die "Write failed on $taxdir/$$ after $id: $!\n";
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

#
# Fetch sequence data for the genomes that have no usable download file, and
# write it out as per-genome files shaped exactly like the download files, so
# the worker below cannot tell the difference.
#
# What this replaces: the worker used to notice the failed open and call
# retrieve_{protein,dna}_features_in_genomes_to_temp([$genome_id]) /
# retrieve_contigs_in_genomes_to_temp([$genome_id]) for that one genome. Those
# helpers loop over their genome-id list and issue a query per genome, so
# passing a list would not have batched anything -- it was one round trip per
# cache miss no matter how it was called. Viral builds run with
# --no-check-files and so never reach this path at all; bacterial builds miss
# on roughly a fifth of their genomes, and the misses clump by genome id
# range, so a single genus can be almost entirely uncached.
#
# Instead do what process_from_data_api already does: one in(genome_id,(...))
# query per --batch-size genomes, demultiplexed into per-genome files by the
# genome_id on each record.
#
sub prefetch_missing_genomes
{
    my($api, $missing, $fetchdir) = @_;

    my %path;
    return \%path unless @$missing;

    printf STDERR "Prefetching %d genome(s) with no download file, in batches of %d\n",
	scalar(@$missing), $opt->batch_size;

    my $ftype_cond = $opt->viral
	? [ "in", "feature_type", "(mat_peptide,CDS)" ]
	: [ "eq", "feature_type", "CDS" ];

    my @todo = @$missing;
    while (@todo)
    {
	my @batch = splice(@todo, 0, $opt->batch_size);
	my %fh;

	#
	# Open every file in the batch up front. A genome that turns out to
	# have no records then still gets an empty file, which is the honest
	# answer; leaving it absent would send the worker down its own
	# per-genome fallback and undo the batching.
	#
	for my $g (@batch)
	{
	    my $p = "$fetchdir/$g->{genome_id}";
	    my $f = IO::File->new($p, "w") or die "Cannot write $p: $!";
	    $fh{$g->{genome_id}} = $f;
	    $path{$g->{genome_id}} = $p;
	}

	my $genome_cond = [ "in", "genome_id",
			    "(" . join(",", map { $_->{genome_id} } @batch) . ")" ];

	if ($ftype eq 'features')
	{
	    my $key = ($dbtype eq 'aa') ? "aa_sequence_md5" : "na_sequence_md5";

	    $api->query_cb("genome_feature",
			   sub {
			       my($data) = @_;
			       my %by_md5;
			       $api->lookup_sequence_data([map { $_->{$key} } @$data], sub {
				   my $ent = shift;
				   $by_md5{$ent->{md5}} = $ent->{sequence};
			       });
			       for my $ent (@$data)
			       {
				   my $f = $fh{$ent->{genome_id}} or next;
				   my $seq = $by_md5{$ent->{$key}};
				   if (!defined($seq) || $seq eq '')
				   {
				       warn "No $key sequence for $ent->{patric_id}, skipping\n";
				       next;
				   }
				   print_alignment_as_fasta($f,
					[$ent->{patric_id}, $ent->{product}, $seq]);
				   $f->error
				       and die "Write failed on $path{$ent->{genome_id}}: $!\n";
			       }
			       return 1;
			   },
			   $ftype_cond,
			   [ "eq", "annotation", "PATRIC" ],
			   $genome_cond,
			   [ "select", "patric_id,product,genome_id,$key" ]);
	}
	else
	{
	    $api->query_cb("genome_sequence",
			   sub {
			       my($data) = @_;
			       for my $ent (@$data)
			       {
				   my $f = $fh{$ent->{genome_id}} or next;
				   #
				   # Match the download .fna header, whose id is
				   # the bare accession: the worker parses
				   # /^(accn\|)?(\S+)/ and builds
				   # lcl|<genome_id>.<accession> from it. A
				   # record with no accession would otherwise
				   # write a header the regex cannot key on.
				   #
				   my $id = $ent->{accession} || $ent->{sequence_id};
				   print_alignment_as_fasta($f,
					[$id,
					 "$ent->{description} [ $ent->{genome_name} | $ent->{genome_id} ]",
					 $ent->{sequence}]);
				   $f->error
				       and die "Write failed on $path{$ent->{genome_id}}: $!\n";
			       }
			       return 1;
			   },
			   $genome_cond,
			   [ "select", "accession,sequence_id,description,genome_name,genome_id,sequence" ]);
	}

	for my $gid (keys %fh)
	{
	    $fh{$gid}->close or die "Close failed on $path{$gid}: $!\n";
	}
    }

    return \%path;
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

