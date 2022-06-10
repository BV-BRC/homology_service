
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
use DBI;

my $json = JSON::XS->new->pretty->canonical;

my($opt, $usage) = describe_options("%c %o db-dir db-file",
				    ["help|h", "Show this help message"],
				    ['curated-directory=s@', "Databases in this directory to be marked as curated. May be repeated", { default => []}],
				    ["process-prefix=s", "Only process files with this prefix"],
				    ["sqlite=s", "Use the given file as a SQLite database"],
				    );
print($usage->text), exit 0 if $opt->help;
die($usage->text) if @ARGV != 2;

my $db_dir = shift;
my $db_file = shift;

my $dbh;
if ($opt->sqlite)
{
    $dbh = DBI->connect('dbi:SQLite:' . $opt->sqlite, undef, undef, {
	AutoCommit => 1,
	RaiseError => 1,
    });
    $dbh->do("PRAGMA foreign_keys = ON");
}

my $process_re;
if ($opt->process_prefix)
{
    my $x = $opt->process_prefix;
    $process_re = qr/^$x/;
}

my %by_tag;

my %curated = map { $_ => 1 } @{$opt->curated_directory};

find(\&process_path, $db_dir);

my @out;
for my $tag (sort keys %by_tag)
{
    push(@out, { name => $tag, db_list => $by_tag{$tag} });
}

# write_file($db_file, $json->encode(\@out));


sub process_path
{
    return unless /\.taxids/;
    return if $process_re && ! /$process_re/;
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

    #
    # Look up database key for this genome if we are building a database.
    #
    my $sql_genome_group;
    if ($dbh)
    {
	my $curated = $curated{dirname($dir)} ? 1 : 0;
	print "$dir " . dirname($dir) . " curated=$curated\n";
	# return;
	$dbh->begin_work();
	$dbh->do(qq(INSERT OR IGNORE INTO GenomeGroup (name) VALUES (?)), undef, $tag);
	my $res = $dbh->selectcol_arrayref(qq(SELECT id FROM GenomeGroup WHERE name = ?), undef, $tag);
	$sql_genome_group = $res->[0];
	# $dbh->commit();
	print "group=$sql_genome_group\n";
    }

    my $item = {
	path => $itempath,
	type => $dbtype,
	ftype => $ftype,
    };

    #
    # Create database record for this database.
    #
    
    my($sql_dbtype, $sql_ftype, $sql_blastdb);
    if ($dbh)
    {
	my($t) = $dbh->selectcol_arrayref(qq(SELECT id FROM DbType WHERE suffix = ?), undef, $dbtype);
	$sql_dbtype = $t->[0];
	($t) = $dbh->selectcol_arrayref(qq(SELECT id FROM FeatureType WHERE name = ?), undef, $ftype);
	# warn Dumper($ftype, $t);
	$sql_ftype = $t->[0];

#	$dbh->begin_work();
	#
	# If we already have one, bail.
	#
	my $res = $dbh->selectcol_arrayref(qq(SELECT id
					      FROM BLastDatabase
					      WHERE
					      genome_group = ? AND
					      dbtype = ? AND
					      ftype = ?), undef, $sql_genome_group, $sql_dbtype, $sql_ftype);
	if (@$res > 0)
	{
	    die "Already have database loaded for $tag $dbtype $ftype\n";
	}
	$res = $dbh->do(qq(INSERT INTO BlastDatabase (genome_group, path, dbtype, ftype)
			   VALUES (?, ?, ?, ?)), undef,
			$sql_genome_group, $itempath, $sql_dbtype, $sql_ftype);
	$sql_blastdb = $dbh->sqlite_last_insert_rowid();
	# print Dumper($res, $sql_blastdb);
#	$dbh->commit();
    }



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
    if ($dbh)
    {
#	$dbh->begin_work();
	my $sth = $dbh->prepare(qq(INSERT INTO GenomeInDatabase (genome_id, blast_database) VALUES (?, $sql_blastdb)));
	for my $g (sort keys %genome)
	{
	    $sth->execute($g);
	}
	undef $sth;
	
	my $sth = $dbh->prepare(qq(INSERT INTO TaxonInDatabase (taxon_id, blast_database) VALUES (?, $sql_blastdb)));
	for my $t (sort keys %tax)
	{
	    $sth->execute($t);
	}
	undef $sth;
	$dbh->commit();
    }

    # $item->{genome_counts} = \%genome;
    # $item->{tax_counts} = \%tax;
    # push(@{$by_tag{$tag}}, $item);
}

