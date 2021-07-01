
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
use File::Path qw(make_path);
use IPC::Run qw(run);
use PerlIO::via::Blockwise;

my($opt, $usage) = describe_options("%c %o dbtype ftype blast-db-file",
				    ['exclude-taxon=i@', "Do not include genomes in this taxon", { default => [] }],
				    ["reference", "Include reference genomes"],
				    ["representative", "Include representative genomes"],
				    ["complete", "Include only genome_status=Complete genomes"],
				    ['taxon=i@', "Limit to this taxon", { default => [] }],
				    ["title=s", "Database title", { default => "blast db"}],
				    ["parallel|j=i", "Use this many threads", { default => 1 }],
				    ["dump", "Don't build, just dump the list of genomes"],
				    ["help|h", "Show this help message"]);

print($usage->text), exit 0 if $opt->help;
die($usage->text)  unless @ARGV == 3 || ($opt->dump && @ARGV == 0);

my $genome_base = "/vol/patric3/downloads/genomes";

my %type_suffix = (aa => "faa",
		   dna => "fna");
my %blast_type = (aa => "prot",
		  dna => "nucl");

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
# Require good genomes.
#
push(@params, fq => 'genome_quality:Good');

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

my $dbfile = "$db.$ftype.$type_suffix{$dbtype}";
my $taxids = "$db.$ftype.$type_suffix{$dbtype}.taxids";

#open(DB, ">", $dbfile) or die "Cannot write $dbfile: $!";
#open(TI, ">", $taxids) or die "Cannot write $taxids: $!";

my $sched = LPTScheduler->new($opt->parallel);

for my $g (@$genomes)
{
    my($genome_id, $genome_name, $taxon_id, $genome_length) = @$g{qw(genome_id genome_name taxon_id genome_length)};
    my $path = compute_path($genome_id, $ftype, $dbtype);
    next unless $path;

    $g->{path} = $path;
    $sched->add_work($g, $genome_length);
}

my $tmpdir = File::Temp->newdir(CLEANUP => 0);;
my $dbdir = "$tmpdir/db";
my $taxdir = "$tmpdir/tax";

make_path($dbdir, $taxdir);

my $boot = sub {
    my $dbfile = "$dbdir/$$";
    my $taxfile = "$taxdir/$$";
    my $db = IO::File->new($dbfile, "w");
    $db or die "Cannot write $dbfile: $!";
    my $tax = IO::File->new($taxfile, "w");
    $tax or die "Cannot write $taxfile: $!";
    return [$db, $tax];
};
    
$sched->run($boot, sub {
    my($global, $work) = @_;
    my($db, $tax) = @$global;

    my($path, $genome_id, $genome_name, $taxon_id, $genome_length) = @$work{qw(path genome_id genome_name taxon_id genome_length)};

    if (!open(P, "<:via(Blockwise)", $path))
    {
	warn "Cannot open $path: $!\n";
	return;
    }
    my $skip;
    my %seen;
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
});

my @taxf = glob("$taxdir/*");
my @dbf = glob("$dbdir/*");
my $ok = run(["cat", @taxf], ">", $taxids);
$ok or die "Failure creating $taxids from @taxf\n";

$ok = run(["cat", @dbf],
	  "|",
	  ["makeblastdb",
	   "-in", "-",
	   "-out", $dbfile,
	   "-parse_seqids",
	   "-taxid_map", $taxids,
	   "-title", $opt->title,
	   "-dbtype", $blast_type{$dbtype}]);

$ok or die "Error  creating blastdb: $!";
	   

sub compute_path
{
    my($genome, $ftype, $dbtype) = @_;
    
    my $dir = "$genome_base/$genome";

    if (! -d $dir)
    {
	warn "No directory $dir\n";
	return undef;
    }

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

