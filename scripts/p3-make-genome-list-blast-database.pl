#
# Given a list of genome IDs and a database type, make a blast db.
#

use Getopt::Long::Descriptive;
use strict;
use Data::Dumper;
use P3DataAPI;

my($opt, $usage) = describe_options("%c %o db-file [genome-id genome-id...]",
				    ["ids-from=s", "Use the given file for genome IDs"],
				    ["db-type=s", "Create database of this type (faa, ffn, rn, fna)", {
					default => 'faa' }],
				    ["dbapi-url=s", "Data API URL"],
				    ["help|h", "Show this help message"]);

print($usage->text), exit 0 if $opt->help;
die($usage->text) if @ARGV < 1;

my $api = P3DataAPI->new($opt->dbapi_url);

my $db_file = shift;

my @genome_ids;
for my $g (@ARGV)
{
    if ($g =~ /^\d+\.\d+$/)
    {
	push(@genome_ids, $g);
    }
    else
    {
	die "Invalid genome ID on command line: $g\n";
    }
}

if ($opt->ids_from)
{
    open(F, "<", $opt->ids_from) or die "Cannot open ids file " . $opt->ids_from . ": $!\n";
    while(<F>)
    {
	if (/^(\d+\.\d+)$/)
	{
	    push(@genome_ids, $1);
	}
	else
	{
	    die "Invalid genome ID at line $. of " . $opt->ids_from . "\n";
	}
    }
}

if ($opt->db_type eq 'faa')
{
    build_prot_db(\@genome_ids, $db_file);
}
elsif ($opt->db_type eq 'ffn')
{
    build_dna_db(\@genome_ids, $db_file);
}
elsif ($opt->db_type eq 'frn')
{
    build_rna_db(\@genome_ids, $db_file);
}
elsif ($opt->db_type eq 'fna')
{
    build_contig_db(\@genome_ids, $db_file);
}
else
{
    die "Invalid database type. Valid values are faa, ffn, frn, fna.\n";
}

sub build_prot_db
{
    my($genome_ids, $db_file) = @_;

    my $tmp = $api->retrieve_protein_features_in_genomes_to_temp($genome_ids);
    my @params = ("-dbtype" => "prot",
		  "-in", "$tmp",
		  "-out", $db_file);

    my $rc = system("makeblastdb", @params);
    if ($rc != 0)
    {
	die "makeblastdb @params failed with rc=$rc\n";
    }
}

sub build_dna_db
{
    my($genome_ids, $db_file) = @_;

    my $tmp = $api->retrieve_dna_features_in_genomes_to_temp($genome_ids);
    my @params = ("-dbtype" => "nucl",
		  "-in", "$tmp",
		  "-out", $db_file);

    my $rc = system("makeblastdb", @params);
    if ($rc != 0)
    {
	die "makeblastdb @params failed with rc=$rc\n";
    }
}

sub build_rna_db
{
    my($genome_ids, $db_file) = @_;

    my $tmp = $api->retrieve_rna_features_in_genomes_to_temp($genome_ids);
    system("cat $tmp");
    my @params = ("-dbtype" => "nucl",
		  "-in", "$tmp",
		  "-out", $db_file);

    my $rc = system("makeblastdb", @params);
    if ($rc != 0)
    {
	die "makeblastdb @params failed with rc=$rc\n";
    }
}

sub build_contig_db
{
    my($genome_ids, $db_file) = @_;

    my $tmp = $api->retrieve_contigs_in_genomes_to_temp($genome_ids);
    my @params = ("-dbtype" => "nucl",
		  "-in", "$tmp",
		  "-out", $db_file);

    my $rc = system("makeblastdb", @params);
    if ($rc != 0)
    {
	die "makeblastdb @params failed with rc=$rc\n";
    }
}

