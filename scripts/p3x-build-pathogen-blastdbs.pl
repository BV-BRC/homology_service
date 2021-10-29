
use IPC::Run qw(run);
use JSON::XS;
use Data::Dumper;
use strict;
use File::Path qw(make_path);
use P3DataAPI;
use Getopt::Long::Descriptive;
use LPTScheduler;

my($opt, $usage) = describe_options("%c %o dbtype ftype output-dir",
				    ['dbtype is either "aa" or "dna"'],
				    ['ftype is either "features" or "contigs"'],
				    [],
				    ["n-build-threads=i", "Number of parallel builds to run", { default => 6 }],
				    ["n-db-threads=i", "Number of database lookup threads to run inside each build", { default => 2}],
				    ["viral", "Build for viral families"],
				    ["all-genera", "Build for all genera instead of selected pathogens"],
				    ["check-size!", "Check the genus size to schedule work", { default => 1 }], 
				    ["help|h", "Show this help message"]);

print($usage->text), exit 0 if $opt->help;
die($usage->text) unless @ARGV == 3;

my $dbtype = shift;
my $ftype = shift;
my $out_dir = shift;

-d $out_dir or die "Output directory $out_dir does not exist\n";

$dbtype eq 'dna' or $dbtype eq 'aa' or die "Invalid dbtype $dbtype (valid values are 'aa' and 'dna')\n";
$ftype eq 'features' or $ftype eq 'contigs' or die "Invalid ftype $ftype (valid values are 'features' and 'contigs')\n";

$ENV{TMPDIR} = "/disks/tmp";

my $build_other = 0;

my $n_procs = $opt->{n_build_threads};
my $n_build_procs = $opt->{n_db_threads};

my $sched = LPTScheduler->new($n_procs);

my $json = JSON::XS->new->pretty->canonical;

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


if ($opt->all_genera)
{
    #
    # Clear the hardcoded list of genera and build for all.
    #
    # Note that if --viral was provided we build for viral families
    # (in which case %genera isn't the correct name but we are leaving it for now).
    #

    %genera = ();

    my $api = P3DataAPI->new();

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
	$genera{$taxon_name} = $taxon_id;
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
    $catchall_taxon = 2;
}
    

for my $genus (sort keys %genera)
{
    my $tax = $genera{$genus};
    my $dat;
    my $create_params = ["--tax", $tax, @create_options];


    my $sz = 1;
    my $glist = [];

    if ($opt->check_size)
    {
	
	my @cmd = ("p3x-create-blast-db",
		   "--dump",
		   @$create_params);
	print STDERR "@cmd\n";
	run(\@cmd, ">", \$dat);
	eval {
	    $glist = $json->decode($dat);
	};
	if ($@)
	{
	    die "Failed to parse for $tax: $@\n$dat\n";
	}
	$sz = 0;
	$sz += $_->{genome_length} foreach @$glist;
	$sz /= 1e6;
    }

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
		">", "log.$genus.out",
		"2>", "log.$genus.err",
		);
    $ok or die "Failure running @cmd: $! $?\n";
});
