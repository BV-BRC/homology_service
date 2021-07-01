
=head1 NAME

p3x-genomes-under-tax

=head1 SYNOPSIS

p3x-genomes-under-tax taxid

=head1 DESCRIPTION

List all the genome IDs present under the given tax id.

If --reference and/or --representative is selected, limit to only those genomes.

=cut

use strict;
use Getopt::Long::Descriptive;
use P3DataAPI;
use Data::Dumper;

my($opt, $usage) = describe_options("%c %o taxon",
				    ["commas", "Use comma as separator instead of newline"],
				    ["taxa", "Generate list of taxonomy ids instead of genomes"],
				    ["reference", "Include reference genomes"],
				    ["representative", "Include representative genomes"],
				    ["complete", "Include only genome_status=Complete genomes"],
				    ["alpha", "Use the alpha.bv-brc.org data api"],
				    ["help|h", "Show this help message"]);

print($usage->text), exit 0 if $opt->help;
die($usage->text)  unless @ARGV == 1;

my $taxon = shift;

my $url;
$url = "https://alpha.bv-brc.org/api" if $opt->alpha;

my $api = P3DataAPI->new($url);

my @params;

my @rr;
push(@rr, "reference_genome:Reference") if $opt->reference;
push(@rr, "reference_genome:Representative") if $opt->representative;

push @params, fq => (join(" OR ", @rr)) if @rr;

push(@params, fq => "genome_status:Complete") if $opt->complete;

push(@params, fq => "taxon_lineage_ids:$taxon");

push(@params, fl => "genome_id,taxon_id");

push(@params, fq => "public:true");

push(@params, q => "*:*");

my $genomes = $api->solr_query_list("genome", \@params);

my $key = $opt->taxa ? "taxon_id" : "genome_id";

my %out;
$out{$_->{$key}} = 1 foreach @$genomes;

my $sep = $opt->commas ? "," : "\n";

print join($sep, sort keys %out), "\n";
