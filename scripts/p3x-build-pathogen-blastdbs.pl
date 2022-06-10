
use IPC::Run qw(run);
use JSON::XS;
use Data::Dumper;
use strict;
use File::Path qw(make_path);
use P3DataAPI;
use Getopt::Long::Descriptive;
use LPTScheduler;

#
# Import this to properly intialize PATH for blast binaries.
#
use Bio::P3::HomologySearch::HomologySearch;

my($opt, $usage) = describe_options("%c %o dbtype ftype output-dir",
				    ['dbtype is either "aa" or "dna"'],
				    ['ftype is either "features" or "contigs"'],
				    [],
				    ["reference", "Include reference genomes"],
				    ["representative", "Include representative genomes"],
				    ["complete", "Include only genome_status=Complete genomes"],
				    ["quality-check!", "Require quality of Good (invert with --no-quality-check)", { default => 1 }],
				    ["n-build-threads=i", "Number of parallel builds to run", { default => 6 }],
				    ["n-db-threads=i", "Number of database lookup threads to run inside each build", { default => 2}],
				    ["viral", "Build for viral families"],
				    ["build-tempdir=s", "Tmpdir for build tool", { default => "/dev/shm" }],
				    ["all-genera", "Build for all genera instead of selected pathogens"],
				    ["check-size!", "Check the genus size to schedule work", { default => 1 }],
				    ["log-output!", "Log output of individual blast runs to files", { default => 1 }],
				    ["log-dir=s", "Directory to write log data", { default => "." }],
				    ["help|h", "Show this help message"]);

print($usage->text), exit 0 if $opt->help;
die($usage->text) unless @ARGV == 3;

my $dbtype = shift;
my $ftype = shift;
my $out_dir = shift;

-d $out_dir or die "Output directory $out_dir does not exist\n";

$dbtype eq 'dna' or $dbtype eq 'aa' or die "Invalid dbtype $dbtype (valid values are 'aa' and 'dna')\n";
$ftype eq 'features' or $ftype eq 'contigs' or die "Invalid ftype $ftype (valid values are 'features' and 'contigs')\n";

make_path($opt->log_dir);

$ENV{TMPDIR} = "/disks/tmp";

my $build_other = 0;

my $n_procs = $opt->{n_build_threads};
my $n_build_procs = $opt->{n_db_threads};

my $sched = LPTScheduler->new($n_procs);

my $json = JSON::XS->new->pretty->canonical;

my $api = P3DataAPI->new();

my %genera = (
	      'Acinetobacter' => 469,
	      'Francisella' => 262,
	      'Bacillus' => 1386,
	      'Helicobacter' => 209,
	      'Bartonella' => 773,
	      'Listeria' => 1637,
	      'Borreliella' => 64895,
	      'Mycobacterium' => 1763,
	      'Brucella' => 234,
	      'Pseudomonas' => 286,
	      'Burkholderia' => 32008,
	      'Rickettsia' => 780,
	      'Campylobacter' => 194,
	      'Salmonella' => 590,
	      'Chlamydia' => 810,
	      'Shigella' => 620,
	      'Clostridium' => 1485,
	      'Staphylococcus' => 1279,
	      'Coxiella' => 776,
	      'Streptococcus' => 1301,
	      'Ehrlichia' => 943,
	      'Vibrio' => 662,
      'Escherichia' => 561,
	      'Yersinia' => 629,
	      #
	      # Plus these larger genera
	      #
	      'Neisseria' => 482,
	      'Klebsiella' => 570,
	      'Enterococcus' => 1350,
	      );


#
# We don't have quality data on viruses.
#
$opt->{quality_check} = 0 if $opt->viral;

#
# If viral, choose all genera.
#
$opt->{all_genera} = 1 if $opt->viral;


my @target_ids;
my %taxon_name;

if ($opt->all_genera)
{
    #
    # Clear the hardcoded list of genera and build for all.
    #
    # Note that if --viral was provided we build for viral families
    # (in which case %genera isn't the correct name but we are leaving it for now).
    #

    %genera = ();

    my @genera;
    if ($opt->viral)
    {
	@genera = $api->query('taxonomy',
			      ['eq', 'lineage_names', 'Viruses'],
			      ['eq', 'taxon_rank', 'family'],
			      ['select', 'taxon_id,taxon_name,lineage_names']);

	#
	# Make sure we didn't pick up any stray taxa with the nonspecific match for Viruses above
	#
	@genera = grep { $_->{lineage_names}->[0] eq 'Viruses'} @genera;
	
	#
	# Testing subset.
	#
	if (0)
	{
	    my @fams = qw(Flaviviridae Herpesviridae Drexlerviridae);
	    my %fams = map { $_ => 1 } @fams;
	    @genera = grep { $fams{$_->{taxon_name} } } @genera;
	}
    }
    else
    {
	@genera = $api->query('taxonomy',
			      ['eq', 'division', 'Bacteria'],
			      ['eq', 'taxon_rank', 'genus'],
			      ['select', 'taxon_id,taxon_name']);
    }
    for my $ent (@genera)
    {
	my($taxon_id, $taxon_name) = @$ent{qw(taxon_id taxon_name)};
	push(@target_ids, $taxon_id);
	$genera{$taxon_name} = $taxon_id;
	$taxon_name{$taxon_id} = $taxon_name;
    }
}
else
{
    #
    # Construct %taxon_name and @target_ids from the %genera hash.
    #
    while (my($name, $id)  = each(%genera))
    {
	$taxon_name{$id} = $name;
	push(@target_ids, $id);
    }
}

#
# Build a database for each taxon; also collect the list of taxa to build
# an exclude list for the "Other" set.
#

my @all;
my @exclude;

#
# Disable quality check for viral genomes.
# Also determine the catchall taxon id.
#
my @create_options;
my $catchall_taxon;
if ($opt->viral)
{
    push(@create_options,
	 "--no-quality-check",
	 "--no-check-files",
	 "--batch-size", 500,
	);
    $catchall_taxon = 10239;
}
else
{
    # push(@create_options, "--no-check-files");
    $catchall_taxon = 2;
}
    
@target_ids = sort { $a <=> $b } @target_ids;

my $genus_size = compute_genome_lists(\@target_ids);
for my $g (sort { $genus_size->{$a} <=> $genus_size->{$b} } keys %$genus_size)
{
    print "$g\t" . int($genus_size->{$g} / 1e6) . "\n";
}
#exit;
	

for my $tax (@target_ids)
{
    my $genus = $taxon_name{$tax};
    my $dat;
    my $create_params = ["--tax", $tax, @create_options];

    my $sz = 1;
    my $glist = [];

    my $sz = $genus_size->{$genus} / 1e6;

    $genus =~ s/\s/_/g;
    my $item = { genus => $genus, tax => $tax, size => $sz, list => $glist, params => $create_params };
    push(@all, $item);
    my $n = @$glist;
    print STDERR "$genus\t$tax\t$n\t$sz\n";
    push(@exclude, $tax);

    $sched->add_work($item, $sz);
}

my $dat;

if ($build_other)
{
    my $create_params = [(map { ("--exclude", $_) } @exclude),
			 @create_options,
			 "--tax", $catchall_taxon];
    
    my @cmd = ("p3x-create-blast-db",
	       "--dump",
	       @$create_params,
	      );
    # print STDERR "@cmd\n";
    run(\@cmd, ">", \$dat);
    my $glist = $json->decode($dat);
    my $sz = 0;
    $sz += $_->{genome_length} foreach @$glist;

    if ($build_other)
    {
	my $item = { genus => 'OTHER', tax => undef, size => $sz, list => $glist, params => $create_params };
	push(@all, $item);
	$sched->add_work($item, $sz);
    }
    my $n = @$glist;
    
    print STDERR "OTHER\t\t$n\t$sz\n";
}

# print $json->encode(\@all);

$sched->run(sub {}, sub {
    my($global, $item) = @_;

    my($genus, $tax, $size, $list, $params) = @$item{qw(genus tax size list params)};

    print STDERR "Run $$ $genus $tax @$params\n";

    my @nr;
    if ($opt->viral)
    {
	# @nr = ("--create-nr");
    }
    my $db = "$out_dir/$genus";
    my @cmd = ("p3x-create-blast-db",
	       @$params,
	       "--title", $genus,
	       "--parallel", $n_build_procs,
	       @nr,
	       $dbtype, $ftype, $db);
    print "@cmd\n";
    my $ok = run(\@cmd,
		 init => sub {
		     if ($size < 50_000)
		     {
			 $ENV{TMPDIR} = $opt->build_tempdir;
		     }
		 },
		 ($opt->log_output ? (">", $opt->log_dir . "/log.$genus.out",
				      "2>", $opt->log_dir . "/log.$genus.err",)
		  : ()),
		);
    $ok or die "Failure running @cmd: $! $?\n";
});

sub compute_genome_lists
{
    my($genera) = @_;
    my @params;
    
    my @rr;
    push(@rr, "reference_genome:Reference") if $opt->reference;
    push(@rr, "reference_genome:Representative") if $opt->representative;
    
    #
    # Require lineage.
    #
    push(@params, fq => 'taxon_lineage_ids:*');

    if (ref($genera) && @$genera)
    {
	push(@params, fq => "taxon_lineage_ids:(" . join(" OR ",  @$genera) . ")");
    }

    my $facet;
    if ($opt->viral)
    {
	push(@params, fq => 'superkingdom:Viruses');
	$facet = 'family';
	push(@params, fq => "-taxon_id:2697049");
    }
    else
    {
	push(@params, fq => '-superkingdom:Viruses');
	$facet = 'genus';
    }
    
    push(@params, fq => "(" . (join(" OR ", @rr)) . ")") if @rr;
    
    push(@params, fq => "genome_status:Complete") if $opt->complete;
    
    push(@params, fl => "genome_id,genome_name,taxon_id,genome_length,genus,family");

    
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

    # unshift(@params, facet => 'true', 'facet.field' => $facet);
    push(@params, "json.facet" => $json->encode({genera =>  { type => 'terms', field => $facet, limit => -1,
								  facet => { dna_size => 'sum(genome_length)' }
							      } }));

    my $buckets;
    my $genomes = $api->solr_query_list("genome", \@params, 1, sub {
	my ($resp) = @_;
	my $facets = $resp->{facets};
	$buckets = $facets->{genera}->{buckets};
	return ();
       });
    my %buckets = map { $_->{val} => $_->{dna_size} } @$buckets;
    return \%buckets;
}
    
