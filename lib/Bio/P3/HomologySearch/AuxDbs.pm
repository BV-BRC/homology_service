=head1 NAME

    AuxDbs - wrapper around the berekely db storing tax & md5 lookup info

=cut

package Bio::P3::HomologySearch::AuxDbs;

use strict;
use BerkeleyDB;

sub new
{
    my($class, $path) = @_;

    my $env = BerkeleyDB::Env->new(-Home => $path);

    my $self = {
	env => $env,
	tbl_tax_to_pseudo => BerkeleyDB::Btree->new(-Filename, "tax_to_pseudo",
						    -Env => $env,
						    -Property => DB_DUP),
	tbl_pseudo_to_tax => BerkeleyDB::Btree->new(-Filename, "pseudo_to_tax",
						    -Env => $env,
						    -Property => DB_DUP),
	tbl_md5_to_feature => BerkeleyDB::Btree->new(-Filename, "md5_to_feature",
						     -Env => $env,
						     -Property => DB_DUP),
    };

    return bless $self, $class;
}

sub get_values_for_all
{
    my($self, $table, $list) = @_;

    $list = [$list] unless ref($list);

    my %ret;
    for my $key (@$list)
    {
	my @vals = $self->get_values($table, $key);
	$ret{$_} = 1 foreach @vals;
    }
    return keys %ret;
}


sub get_values
{
    my($self, $table, $key) = @_;

    my($k, $v) = ($key, "");
    my $c = $table->db_cursor();

    my @ret;
    my $ok = $c->c_get($k, $v, DB_SET);
    while ($ok != DB_NOTFOUND)
    {
	push(@ret, $v);
	$ok = $c->c_get($k, $v, DB_NEXT_DUP);
    }
    return @ret;
}

sub taxon_id_to_pseudo
{
    my($self, $tax) = @_;

    return $self->get_values_for_all($self->{tbl_tax_to_pseudo}, $tax);
}

sub pseudo_to_taxon_id
{
    my($self, $pseudo) = @_;

    return $self->get_values_for_all($self->{tbl_pseudo_to_tax}, $pseudo);
}


sub md5_to_features
{
    my($self, $md5) = @_;

    return $self->get_values($self->{tbl_md5_to_feature}, $md5);
}

1;
