
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
				    ["exclude-deprecated!", "Exclude genome_status=Deprecated genomes (invert with --no-exclude-deprecated)", { default => 1 }],
				    ["n-build-threads=i", "Number of parallel builds to run", { default => 6 }],
				    ["n-db-threads=i", "Number of database lookup threads to run inside each build", { default => 2}],
				    ["viral", "Build for viral families"],
				    ['refrep-only-taxon=i@', "Within this taxon include only reference and representative genomes (passed through to p3x-create-blast-db). Repeatable", { default => [] }],
				    ["sars2-carveout!", "For viral builds, restrict SARS-CoV-2 (taxon 2697049) to its reference and representative genomes. On by default; --no-sars2-carveout includes all 9.4M clinical depositions", { default => 1 }],
				    ["build-tempdir=s", "Tmpdir for build tool", { default => "/dev/shm" }],
				    ["all-genera", "Build for all genera instead of selected pathogens"],
				    ["clade-ladder", "Build the genomes that the genus (or, with --viral, family) partition leaves out: one database per deepest taxonomic rank each genome actually carries. Replaces the genus enumeration entirely; use a separate output directory"],
				    ['taxon-filter=s%' => "Define a taxon filter file", { default => {} }],
				    ["check-size!", "Check the genus size to schedule work", { default => 1 }],
				    ["log-output!", "Log output of individual blast runs to files", { default => 1 }],
				    ["log-dir=s", "Directory to write log data", { default => "." }],
				    ["dump-sizes", "Run the size estimation only and dump data"],
				    ["overwrite|f", "Rebuild databases that already exist instead of skipping them"],
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

#
# Clear the previous run's failure markers, or a clean rerun still reports the
# failures it just repaired.
#
unlink($_) for glob($opt->log_dir . "/FAILED.*");

$ENV{TMPDIR} = "/disks/tmp";

#
# Superseded by --clade-ladder; left at 0, and do not turn it on.
#
# The "other" database is meant to catch the genomes the genus partition leaves
# out, but as written it cannot: it passes one --exclude per built genus, which
# is ~15,000 "-taxon_lineage_ids:<n>" clauses in a single Solr query (the API
# tier wedged on a 165 KB request body on 2026-09-03), and it is a *single* LPT
# work item, so 700 Gbp would build single-threaded while 11 of 12 workers idle.
#
# --clade-ladder covers the same genomes as one database per deepest rank each
# genome actually carries -- parallel, and reachable by real lineage.
#
my $build_other = 0;

my $n_procs = $opt->{n_build_threads};
my $n_build_procs = $opt->{n_db_threads};

#
# Create options to propagate taxon filters
#
my @taxon_opts;
while (my($k,$v) = each %{$opt->taxon_filter})
{
    push(@taxon_opts, "--taxon-filter", "$k=$v");
}

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

if ($opt->clade_ladder)
{
    #
    # The ladder does not enumerate taxa at all -- it selects on the genome
    # record's own rank fields, so there is no name-to-taxon-id map to build.
    #
    %genera = ();
}
elsif ($opt->all_genera)
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
# SARS-CoV-2 is 99.2% of Coronaviridae -- 9,413,387 of 9,487,138 genomes, all
# but 22 of them clinical surveillance depositions with no reference or
# representative flag. Enumerating that family is roughly 380 cursor pages and
# does not survive the data API (it 502'd the 2026-08 rebuild), and the database
# it would build is an estimated 150 GB, 30x what the rest of the viral set
# takes together. Restricting the one taxon to reference and representative
# genomes leaves the other 73,751 Coronaviridae untouched.
#
use constant SARS2_TAXON => 2697049;

my @refrep_only = @{$opt->refrep_only_taxon};
push(@refrep_only, SARS2_TAXON) if $opt->viral && $opt->sars2_carveout;

#
# Disable quality check for viral genomes.
# Also determine the catchall taxon id.
#
my @create_options;
my $catchall_taxon;
if ($opt->viral)
{
    #
    # makeblastdb deterministically rejects a small number of records in some
    # viral sets, so the record-count check in p3x-create-blast-db needs a
    # tolerance here. Bacterial builds match exactly and keep the default of 0.
    # Measured on data.2022-0916: 45 viral databases were short, all by 2.4% or
    # less and most by under 0.2%.
    #
    push(@create_options,
	 "--no-quality-check",
	 "--no-check-files",
	 "--batch-size", 500,
	 "--max-missing-fraction", 0.03,
	);
    $catchall_taxon = 10239;
}
else
{
    # push(@create_options, "--no-check-files");
    $catchall_taxon = 2;
}

#
# p3x-create-blast-db skips any database that already has a non-empty .taxids
# plus its blast files, so without this a rerun into a populated directory
# cannot rebuild anything -- which is how data.2022-0916 shipped damaged.
# Pass the flag through so a targeted repair is possible.
#
push(@create_options, "--overwrite") if $opt->overwrite;

#
# Both sides default to excluding Deprecated, so the common case needs no
# passthrough. Disabling it on the driver alone would leave base_query_params
# sizing a smaller set than the child builds -- propagate it explicitly.
#
push(@create_options, "--no-exclude-deprecated") unless $opt->exclude_deprecated;

push(@create_options, map { ("--refrep-only-taxon", $_) } @refrep_only);

push(@create_options, @taxon_opts);
    
@target_ids = sort { $a <=> $b } @target_ids;

my $genus_size;
my $ladder;

if ($opt->clade_ladder)
{
    $ladder = compute_clade_lists();
    for my $c (@$ladder)
    {
	print join("\t", $c->{tag}, ($c->{rank} // 'residue'),
		   int($c->{bp} / 1e6)) . "\n";
    }
    printf STDERR "clade ladder: %d databases, %.1f Gbp\n",
	scalar(@$ladder), (eval { my $t = 0; $t += $_->{bp} for @$ladder; $t } / 1e9);
}
else
{
    $genus_size = compute_genome_lists(\@target_ids);
    for my $g (sort { $genus_size->{$a} <=> $genus_size->{$b} } keys %$genus_size)
    {
	print "$g\t" . int($genus_size->{$g} / 1e6) . "\n";
    }
}
exit if $opt->dump_sizes;

#
# Ladder items. One database per clade, keyed on the clade *name* rather than
# on a taxon id: --field/--unset-field select on the genome record's own rank
# fields, so there is no taxon to pass, and the names are unique within the
# ladder and collision-free against every tag already built.
#
# The catalog still resolves these by taxon, because the .taxlist sidecar is
# derived from the genomes that land in the database, not from the selection.
#
if ($opt->clade_ladder)
{
    my $spec = clade_ladder_spec();
    my @roots = map { ("--taxon", $_) } @{$spec->{roots}};

    for my $c (@$ladder)
    {
	my $create_params = [ @roots,
			      ($c->{rank} ? ("--field", "$c->{rank}=$c->{name}") : ()),
			      (map { ("--unset-field", $_) } @{$c->{unset}}),
			      @create_options ];

	my $sz = $c->{bp} / 1e6;
	my $item = { genus => $c->{tag}, tax => ($c->{rank} // 'residue'),
		     taxa => [], size => $sz, list => [], params => $create_params };
	push(@all, $item);
	print STDERR "$c->{tag}\t$item->{tax}\t0\t$sz\n";

	$sched->add_work($item, $sz);
    }
}
	

#
# A genus name does not identify a genus. Taxonomy carries two distinct genera
# named Ponticoccus (519422, Actinobacteria; 983507, Alphaproteobacteria), and
# three further names collide the same way: "Candidatus Dwaynesavagella",
# Pseudoruminococcus, Venteria. Every output path here is keyed on the name, so
# iterating taxon ids built one database twice and the second build silently
# replaced the first -- every genome of the losing genus was dropped, and the
# only trace was a "already exists, skipping build" line in a child .err log.
#
# Build one database per *name*, passing every taxon that bears it.
# p3x-create-blast-db ORs its --taxon values into a single lineage filter, and
# the catalog maps taxa onto databases through the .taxlist sidecar rather than
# through the filename, so both genera remain individually addressable in
# db.sqlite. This also matches the size we schedule with: compute_genome_lists
# facets on the Solr `genus` field, which is the name, so its count already
# covers both genera and was never per-taxon to begin with.
#
my %ids_for_name;
push(@{$ids_for_name{$taxon_name{$_}}}, $_) for @target_ids;

for my $genus_name (sort keys %ids_for_name)
{
    my $taxa = [ sort { $a <=> $b } @{$ids_for_name{$genus_name}} ];
    my $genus = $genus_name;

    if (@$taxa > 1)
    {
	print STDERR "Note: '$genus_name' names " . scalar(@$taxa) .
	    " distinct taxa (" . join(", ", @$taxa) . "); " .
	    "building them into one database.\n";
    }

    my $create_params = [(map { ("--taxon", $_) } @$taxa), @create_options];

    my $glist = [];
    my $sz = $genus_size->{$genus} / 1e6;

    $genus =~ s/\s/_/g;
    my $tax_label = join(",", @$taxa);
    my $item = { genus => $genus, tax => $tax_label, taxa => $taxa,
		 size => $sz, list => $glist, params => $create_params };
    push(@all, $item);
    my $n = @$glist;
    print STDERR "$genus\t$tax_label\t$n\t$sz\n";
    push(@exclude, @$taxa);

    $sched->add_work($item, $sz);
}

my $dat;

if ($build_other)
{
    my $create_params = [(map { ("--exclude", $_) } @exclude),
			 @create_options,
			 "--taxon", $catchall_taxon];
    
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
print Dumper($sched);

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
    #
    # Do not die here. LPTScheduler::run forks one child per bin and each child
    # walks its bin sequentially, so a die abandons every item the bin has not
    # reached yet -- and the parent never learns, because the scheduler
    # discards the child exit status. The first viral clade run lost 182 of 901
    # databases to two failures that had nothing to do with them.
    #
    # Record the failure as a marker file instead and keep going; the parent
    # counts the markers after the run and exits nonzero, so a rerun (the skip
    # guard makes it cheap) repairs exactly what failed.
    #
    if (!$ok)
    {
	my $err = "Failure running @cmd: $! $?";
	warn "$err\n";
	if (open(my $fh, ">", $opt->log_dir . "/FAILED.$genus"))
	{
	    print $fh "$err\n";
	    close($fh);
	}
    }
});

my @failed = sort map { s/.*\/FAILED\.//r } glob($opt->log_dir . "/FAILED.*");
if (@failed)
{
    print STDERR "\n" . scalar(@failed) . " database(s) failed to build:\n";
    print STDERR "    $_\n" for @failed;
    print STDERR "Rerun the same command to repair them; databases that " .
	"succeeded are skipped.\n";
    exit 1;
}

#
# The filters every sizing query shares with the builds it is sizing. Kept in
# one place so compute_genome_lists and compute_clade_lists cannot drift: a
# sizing query that does not match the build's selection schedules work that
# does not exist (see the Coronaviridae note below).
#
sub base_query_params
{
    my @params;

    my @rr;
    push(@rr, "reference_genome:Reference") if $opt->reference;
    push(@rr, "reference_genome:Representative") if $opt->representative;

    #
    # Require lineage.
    #
    push(@params, fq => 'taxon_lineage_ids:*');

    push(@params, fq => ($opt->viral ? 'superkingdom:Viruses' : '-superkingdom:Viruses'));

    push(@params, fq => "(" . (join(" OR ", @rr)) . ")") if @rr;

    #
    # Mirror the per-taxon carve-out the builds will use, so the size this
    # returns is the size of the database that actually gets built. Without it
    # Coronaviridae is scheduled as ~9.5M genomes of sequence it will not
    # contain, and LPT puts it first and alone.
    #
    for my $rr_taxon (@refrep_only)
    {
	push(@params, fq => "-(taxon_lineage_ids:$rr_taxon -reference_genome:(Reference OR Representative))");
    }
    
    push(@params, fq => "genome_status:Complete") if $opt->complete;

    #
    # Must match p3x-create-blast-db's own selection, or the sizing schedules
    # work that will not be built. See the passthrough below.
    #
    push(@params, fq => '-genome_status:Deprecated') if $opt->exclude_deprecated;

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

    return @params;
}

#
# The deepest-available-rank ladder.
#
# The genus (bacterial) and family (viral) partitions leave out every genome
# that carries no value at that rank -- 270,225 bacterial/archaeal genomes
# (701.6 Gbp, almost all MAGs and uncultured environmental bins) and 191,648
# viral ones as of 2026-09. They are in no database at all, so no taxon_list
# search can reach them.
#
# The ladder covers them by walking outward from the missing rank and building
# one database per clade at the first rank each genome actually has. Rung i
# selects "has <rung i>, and has none of the ranks below it", which makes the
# rungs disjoint by construction; the residue database takes what is left.
#
# Viral rung 1 is *genus*, which is deeper than the family the viral partition
# is keyed on: 13,187 viral genomes carry a genus but no family, a taxonomy
# hole rather than an unclassified genome.
#
# species is deliberately not a rung. 270,214 of the 270,225 bacterial genomes
# have one, but they are placeholders -- "uncultured Lachnospiraceae bacterium"
# is 16,825 genomes -- so it would produce 9,078 databases whose names restate
# the clade above them.
#
sub clade_ladder_spec
{
    return $opt->viral
	? { missing  => 'family',
	    rungs    => [qw(genus order class phylum kingdom)],
	    roots    => [10239],
	    residue  => 'Viruses-unclassified' }
	: { missing  => 'genus',
	    rungs    => [qw(family order class phylum)],
	    roots    => [2, 2157],
	    residue  => 'Bacteria-Archaea-unclassified' };
}

#
# One entry per database to build, in ladder order:
#   { rank, name, tag, bp, unset => [ranks that must be absent] }
# The residue entry has rank and name undef.
#
sub compute_clade_lists
{
    my $spec = clade_ladder_spec();
    my @unset = ($spec->{missing});
    my @out;

    for my $rank (@{$spec->{rungs}})
    {
	my @params = base_query_params();
	push(@params, fq => "taxon_lineage_ids:(" . join(" OR ", @{$spec->{roots}}) . ")");
	push(@params, fq => "-$_:*") for @unset;
	push(@params, fq => "$rank:*");
	push(@params, "json.facet" =>
	     $json->encode({ b => { type => 'terms', field => $rank, limit => -1,
				    facet => { dna_size => 'sum(genome_length)' } } }));

	my $buckets;
	$api->solr_query_list("genome", \@params, 1, sub {
	    $buckets = $_[0]->{facets}->{b}->{buckets};
	    return ();
	});

	for my $b (sort { $a->{val} cmp $b->{val} } @$buckets)
	{
	    my $tag = $b->{val};
	    $tag =~ s/\s/_/g;
	    push(@out, { rank => $rank, name => $b->{val}, tag => $tag,
			 bp => $b->{dna_size}, unset => [@unset] });
	}

	push(@unset, $rank);
    }

    #
    # The residue: no value at any rung.
    #
    my @params = base_query_params();
    push(@params, fq => "taxon_lineage_ids:(" . join(" OR ", @{$spec->{roots}}) . ")");
    push(@params, fq => "-$_:*") for @unset;
    push(@params, "json.facet" => $json->encode({ dna_size => 'sum(genome_length)' }));

    my $bp;
    $api->solr_query_list("genome", \@params, 1, sub {
	$bp = $_[0]->{facets}->{dna_size};
	return ();
    });

    push(@out, { rank => undef, name => undef, tag => $spec->{residue},
		 bp => ($bp || 0), unset => [@unset] });

    return \@out;
}

sub compute_genome_lists
{
    my($genera) = @_;

    my @params = base_query_params();

    if (ref($genera) && @$genera)
    {
	push(@params, fq => "taxon_lineage_ids:(" . join(" OR ",  @$genera) . ")");
    }

    my $facet = $opt->viral ? 'family' : 'genus';

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
    
