
=head1 NAME

p3x-create-databases-lookup

=head1 SYNOPSIS

p3x-create-databases-lookup database-directory databases-file

=head1 DESCRITPION

Given a hierarchy containing BLAST databases, create the databases.json file that describes them.

We can find the databases by inspecting the .taxids files as created by L<p3x-create-blast-db>.

=cut

use strict;
use Getopt::Long::Descriptive;
use Data::Dumper;
use File::Find;
use File::Slurp;
use File::Basename;
use Cwd qw(getcwd abs_path);
use JSON::XS;

my $json = JSON::XS->new->pretty->canonical;

my($opt, $usage) = describe_options("%c %o db-dir db-file",
				    ["help|h", "Show this help message"],
				    );
print($usage->text), exit 0 if $opt->help;
die($usage->text) if @ARGV != 2;

my $db_dir = shift;
my $db_file = shift;

my %by_tag;

find(\&process_path, $db_dir);

my @out;
for my $tag (sort keys %by_tag)
{
    push(@out, { name => $tag, db_list => $by_tag{$tag} });
}

write_file($db_file, $json->encode(\@out));


sub process_path
{
    return unless /\.taxids/;
    my $taxids_file = $_;
    my $dir = $File::Find::name;
    $dir =~ s,^$db_dir(/?),,;

    my $base = basename($_, '.taxids');

    my($tag, $ftype, $dbtype) = $base =~ /^(.*)\.([^.]+)\.([^.]+)$/;

    print "$dir: tag=$tag ftype=$ftype dbtype=$dbtype\n";
    if ($ftype eq 'features' && $dbtype eq 'fna')
    {
	$dbtype = 'ffn';
    }

    my $itempath = $dir;
    $itempath =~ s/\.taxids$//;

    my $item = {
	path => $itempath,
	type => $dbtype,
	ftype => $ftype,
    };

    my %tax;
    my %genome;
    #
    # Look for the shorter files created by later blast builder.
    #
    if (open(F, "<", "$base.glist"))
    {
	while (<F>)
	{
	    if (/^(\d+\.\d+)/)
	    {
		$genome{$1}++;
	    }
	}
	close(F);
    }
    if (open(F, "<", "$base.taxlist"))
    {
	while (<F>)
	{
	    if (/^(\d+)/)
	    {
		$tax{$1}++;
	    }
	}
	close(F);
    }

    #
    # If we didn't find one or the other, we need to scan taxids.
    # We will write the other lists as well when we are done.
    #
    if (%tax == 0 || %genome == 0)
    {
	%tax = ();
	%genome = ();

	open(F, "<", $taxids_file) or die "cannot open $taxids_file: $!";
	print "Scanning taxids $taxids_file\n";
	while (<F>)
	{
	    if (/fig\|(\d+\.\d+)[^\t]+\t(\d+)/)
	    {
		$genome{$1}++;
		$tax{$2}++;
	    }
	    elsif (/lcl\|(\d+\.\d+)[^\t]+\t(\d+)/)
	    {
		$genome{$1}++;
		$tax{$2}++;
	    }
	}
	close(F);
	if (open(F, ">", "$base.taxlist"))
	{
	    print F "$_\n" foreach sort { $a <=> $b } keys %tax;
	    close(F);
	}
	if (open(F, ">", "$base.glist"))
	{
	    print F "$_\n" foreach sort { $a <=> $b } keys %genome;
	    close(F);
	}
    }
    $item->{genome_counts} = \%genome;
    $item->{tax_counts} = \%tax;
    push(@{$by_tag{$tag}}, $item);
}

