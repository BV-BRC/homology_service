use strict;
use Data::Dumper;
use DBI;

@ARGV == 1 or die "$0 dbfile\n";

my $db = shift;

my $dbh = DBI->connect("dbi:SQLite:$db", undef, undef, {
    AutoCommit => 1,	     
    RaiseError => 1,	     
});	

$dbh->do("PRAGMA foreign_keys = ON");

$dbh->do("DELETE FROM TaxonLineage");

my $sth = $dbh->prepare(qq(SELECT parent FROM TaxNode WHERE tax_id = ?));
my $isth = $dbh->prepare(qq(INSERT INTO TaxonLineage(taxon_id, lineage_id) VALUES (?, ?)));
my $res = $dbh->selectcol_arrayref(qq(SELECT DISTINCT taxon_id FROM TaxonInDatabase));
    $dbh->begin_work();
while (@$res)
{
    print scalar(@$res), "\n";
    my @batch = splice(@$res, 0, 10000);


    for my $ent (@batch)
    {
	my $cur = $ent;
	#
	# Our own taxon id is part of the lineage.
	#
	$isth->execute($cur, $cur);
	while (1)
	{
	    $sth->execute($cur);
	    my($p) = $sth->fetchrow_array();
	    last if $p eq $cur;
	    $isth->execute($ent, $p);
	    
	    $cur = $p;
	}
    }
    
}
    $dbh->commit();

