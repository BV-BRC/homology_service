
=head1 NAME

BlastDatabases - wrapper for precomputed BLAST databases

=head1 DESCRIPTION

This class is used to wrap the JSON description of the precomputed
BLAST databases as generated with L<p3x-create-blast-db> and summarized
in the databases.json file by L<p3x-create-databases-lookup>.

It provides methods to expand databases from name and type to
the database path, and to search the set of precomputed databases
for matches based on taxon id or genome id.

=head2 METHODS

=over 4

=cut


package Bio::P3::HomologySearch::BlastDatabases;

use strict;
use File::Basename;
use Module::Metadata;
use File::Slurp;
use JSON::XS;
use P3DataAPI;

use base 'Class::Accessor';
__PACKAGE__->mk_accessors(qw(api search_path));

=item B<Constructor>

    $db = Bio::P3::HomologySearch::BlastDatabases->new($search_path)

Create a new BlastDatabases object. C<$search_path> is a list of
base directories for the BLAST database files.

=cut

sub new
{
    my($class, $search_path, $data_api) = @_;

    $data_api //= P3DataAPI->new();

    my $self = {
	search_path => $search_path,
	api => $data_api
    };
    bless $self, $class;
    $self->load_databases();
    return $self;
}

sub load_databases
{
    my($self) = @_;
    #
    # Open databases file and load.
    #

    my $mpath = Module::Metadata->find_module_by_name(__PACKAGE__);
    $mpath = dirname($mpath);
    my $db_path = "$mpath/databases.json";

    eval {
	my $db_txt = read_file($db_path);
	$db_txt or die "Cannot load databases from $db_path: $!";
	$self->databases(decode_json($db_txt));
    };
    if ($@)
    {
	die "Failure loading database file $db_path: $@";
    }
}

sub databases
{
    my($self, $val) = @_;
    if (defined($val))
    {
	$self->{databases} = $val;
    }
    return wantarray ? @{$self->{databases}} : $self->{databases};
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

    my @to_search;

    for my $db ($self->databases)
    {
	for my $bdb (@{$db->{db_list}})
	{
	    if ($bdb->{type} eq $dbtype)
	    {
		if (my @c = grep { $bdb->{tax_counts}->{$_}} keys %desired_taxa)
		{
		    push(@to_search, [$bdb, \@c]);
		    delete $desired_taxa{$_} foreach @c;
		}
	    }
	}
    }
    return (\@to_search, [keys %desired_taxa]);
}

=item B<find_databases_for_taxa>

    ($taxa, $files) = $db->find_databases_for_taxa($dbtype, $taxa)

Given a list of taxa, search the Solr database for all taxa that include the given
taxa in their lineage (effectively enumerating the taxonomy tree beneath the requested
taxa).

Determine the set of precomputed databases that cover those taxa, and the ids within that
are to be specified in the blast invocation.

=cut
    
sub find_databases_for_taxa
{
    my($self, $dbtype, $taxa) = @_;

    my $expanded_taxa = $self->expand_taxa($taxa);

    my($to_search, $missing) = $self->search_taxa($dbtype, $expanded_taxa);

    if (@$to_search == 0)
    {
	return undef;
    }
    elsif (@$missing > 0)
    {
	warn "Could not find databases for taxa @$missing\n";
    }

    my @taxa;
    my @files;
    for my $to (@$to_search)
    {
	my($bdb, $taxa) = @$to;
	my $p = $self->find_db_in_path($bdb->{path});
	if ($p)
	{
	    push(@taxa, @$taxa);
	    push(@files, $p);
	}
	else
	{
	    warn "Could not find database $bdb->{path}\n";
	}
    }

    return(\@taxa, \@files);
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

    my @cands = grep { $_->{name} eq $name } $self->databases;
    for my $cand (@cands)
    {
	for my $db (grep { $_->{type} eq $type } @{$cand->{db_list}})
	{
	    for my $p (@{$self->{search_path}})
	    {
		my $file = "$p/$db->{path}";
		# print STDERR "CHECK $file\n";
		for my $suf (qw(pal nal psq nsq))
		{
		    my $chk = "$file.$suf";
		    if (-f $chk)
		    {
			return ($file);
		    }
		}
	    }
	}
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
