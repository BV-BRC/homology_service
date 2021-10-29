
=head1 NAME

BlastDatabasesSQL - wrapper for precomputed BLAST databases

=head1 DESCRIPTION

This class is used to wrap the SQLite3 database of the precomputed
BLAST databases as generated with L<p3x-create-blast-db> and summarized
in a SQLite3 database by L<p3x-create-databases-lookup>.

It provides methods to expand databases from name and type to
the database path, and to search the set of precomputed databases
for matches based on taxon id or genome id.

=head2 METHODS

=over 4

=cut


package Bio::P3::HomologySearch::BlastDatabasesSQL;

use strict;
use Data::Dumper;
use File::Basename;
use Module::Metadata;
use File::Slurp;
use JSON::XS;
use P3DataAPI;

use base 'Class::Accessor';
__PACKAGE__->mk_accessors(qw(api search_path db_file dbh));

=item B<Constructor>

    $db = Bio::P3::HomologySearch::BlastDatabasesSQL->new($search_path)

Create a new BlastDatabases object. C<$search_path> is a list of
base directories for the BLAST database files.

=cut

sub new
{
    my($class, $search_path, $db_file, $data_api) = @_;

    $data_api //= P3DataAPI->new();

    my $self = {
	search_path => (ref($search_path) ? $search_path : [$search_path]),,
	api => $data_api,
	db_file => $db_file,
    };

    my $dbh = DBI->connect("dbi:SQLite:$db_file", undef, undef, {
	AutoCommit => 1,	     
	RaiseError => 1,	     
    });	

    $dbh->do("PRAGMA foreign_keys = ON");

    $self->{dbh} = $dbh;
    
    bless $self, $class;
    return $self;
}

sub databases
{
    my($self) = @_;

    my $res = $self->dbh->selectcol_arrayref(qq(SELECT name FROM GenomeGroup WHERE curated));
    
    return wantarray ? @$res : $res;
}

=item B<search_taxa>

    ($matches, $missing) = $db->search_taxa($dbtype, $taxa)

Returns in C<$matches> a list of matching databases in the form of C<($db_desc, $taxon_list)>
and a list of taxa not found in C<$missing>.

=cut


sub search_taxa
{
    my($self, $dbtype, $taxa) = @_;

    my %desired_taxa = map { $_ => 1 } @$taxa;

    my @qry;
    my @params;

    push(@qry, "l.lineage_id IN (" . join(", ", map { "?" } @$taxa) . ")");
    push(@params, @$taxa);

    push(@qry, "t.suffix = ?");
    push(@params, $dbtype);

    #
    # Don't pull in curated databases here, since they will replicate
    # taxa (and the non-curated ones are more complete)
    #
    push(@qry, 'NOT g.curated');

    my $where = join(" AND ", @qry);

    my $qry = qq(SELECT DISTINCT b.path, d.taxon_id
		 FROM TaxonLineage l JOIN TaxonInDatabase d ON d.taxon_id = l.taxon_id
		 JOIN BlastDatabase b ON b.id = d.blast_database
		 JOIN DbType t ON t.id = b.dbtype
		 JOIN GenomeGroup g on g.id = b.genome_group
		 WHERE $where);

    print "$qry\n";
    my $res = $self->dbh->selectall_arrayref($qry, undef, @params);

    my %dbs;

    for my $ent (@$res)
    {
	my($path, $tax) = @$ent;
	push @{$dbs{$path}}, $tax;
    }

    my @to_search;
    while (my($path, $list) = each %dbs)
    {
	delete $desired_taxa{$_} foreach @$list;
	push(@to_search, [$path, $list]);
    }
    
    return (\@to_search, [keys %desired_taxa]);
}

=item B<find_databases_for_taxa>

    ($taxa, $files) = $db->find_databases_for_taxa($dbtype, $taxa)

Determine the set of precomputed databases that cover those taxa, and the ids within that
are to be specified in the blast invocation.

=cut
    
sub find_databases_for_taxa
{
    my($self, $dbtype, $taxa) = @_;

    my($to_search, $missing) = $self->search_taxa($dbtype, $taxa);

    if (@$to_search == 0)
    {
	return undef;
    }
    elsif (@$missing > 0)
    {
	warn "Could not find databases for taxa @$missing\n";
    }

    my %taxa;
    my @files;
    for my $to (@$to_search)
    {
	my($path, $taxa) = @$to;
	my $p = $self->find_db_in_path($path);
	if ($p)
	{
	    $taxa{$_} = 1 foreach @$taxa;
	    push(@files, $p);
	}
	else
	{
	    warn "Could not find database $path\n";
	}
    }

    return([keys %taxa], \@files);
}

=item B<expand_taxa>

    $expanded = $db->expand_taxa($taxa)

Given a list of taxa C<$taxa>, expand to include the taxa below each in the
taxonomy tree.

=cut

sub expand_taxa
{
    my($self, $taxa, @qry) = @_;
    my %desired_taxa;

    my $qry = '(' . join(",", @$taxa) . ')';
    
    my @res = $self->api->query('genome', ['in', 'taxon_lineage_ids', $qry], ['select', 'taxon_id'], @qry);
    # print Dumper(\@res);

    if (@res)
    {
	$desired_taxa{$_->{taxon_id}}++  foreach @res;
	my @desired_taxa = keys %desired_taxa;
	return \@desired_taxa;
    }
    else
    {
	return undef;
    }
}

=item B<find_precomputed_database>

    $path = $db->find_precomputed_database($name, $type)

=cut

sub find_precomputed_database
{
    my($self, $name, $type) = @_;

    my $res = $self->dbh->selectcol_arrayref(qq(SELECT b.path
						FROM GenomeGroup g JOIN BlastDatabase b ON g.id = b.genome_group
						JOIN DbType d ON d.id = b.dbtype
						WHERE g.name = ? AND d.suffix = ?), undef,
					     $name, $type);

    for my $path (@$res)
    {
	my $db_path = $self->find_db_in_path($path);
	return $db_path if $db_path;
    }

    return undef;
}
   
=item B<find_db_in_path>

    $file = $db->find_db_in_path($item_path)

Given an C<$item_path> from a database entry, return the full pathname to that item or
undef if it is missing.

=cut

sub find_db_in_path
{
    my($self, $item_path) = @_;
    
    for my $p (@{$self->{search_path}})
    {
	my $file = "$p/$item_path";
	for my $suf (qw(pal nal psq nsq))
	{
	    my $chk = "$file.$suf";
	    if (-f $chk)
	    {
		return ($file);
	    }
	}
    }
    return undef;
}



=pod	

=back

=cut

1;
