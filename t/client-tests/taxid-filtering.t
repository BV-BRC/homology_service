#
# Taxonomic filtering must actually filter.
#
# BLAST needs taxonomy4blast.sqlite3 to honour -taxids / -taxidlist, and it
# looks for that file only along BLASTDB -- not beside the database it is
# searching. When the file is missing BLAST prints a warning on stderr and
# then returns UNFILTERED results with exit status 0, so a broken deployment
# is indistinguishable from a working one unless something checks the hits.
#
# That is exactly why this test asserts on the staxids returned rather than on
# the exit status: a test that only checks rc=0 passes in the broken state,
# which is how this survived in production undetected.
#
# Point HOMOLOGY_TEST_GENERATION at a generation directory to test one that is
# not the configured search path.
#

use strict;
use Test::More;
use File::Temp 'tempdir';

# Loading this sets both PATH (blast binaries) and BLASTDB (taxonomy file).
use_ok('Bio::P3::HomologySearch::HomologySearch') or BAIL_OUT("cannot load module");
Bio::P3::HomologySearch::Config->import(qw(blast_db_search_path));

my @cfg = @{ Bio::P3::HomologySearch::Config::blast_db_search_path() };

ok(defined($ENV{BLASTDB}), "BLASTDB is set at module load")
    or diag("without it every -taxids is silently ignored");
like($ENV{BLASTDB} // '', qr/\Q$cfg[0]\E/, "BLASTDB includes the configured search path");

#
# The override has to be pushed into BLASTDB by hand: the module set it at load
# time from the configured path, and nothing re-reads it afterwards.
#
my @paths = @cfg;
if (my $g = $ENV{HOMOLOGY_TEST_GENERATION})
{
    @paths = ($g);
    $ENV{BLASTDB} = "$g:$ENV{BLASTDB}";
}

#
# Everything above is hermetic and pins the code change itself. Everything below
# needs a real generation and a real BLAST, so it is skipped rather than failed
# when either is absent -- a dev tree whose generated Config.pm still points at
# the Makefile default is indistinguishable here from a broken deployment, and
# the deployment case is covered by the warning the module emits at load.
#
my($taxfile) = grep { -f "$_/taxonomy4blast.sqlite3" } @paths;
my $have_blast = `which blastp 2>/dev/null` && `which blastdbcmd 2>/dev/null`;

if (!$taxfile || !$have_blast)
{
    diag(!$taxfile
	 ? "no taxonomy4blast.sqlite3 under @paths -- cannot exercise filtering; "
	   . "point HOMOLOGY_TEST_GENERATION at a published generation to run the rest"
	 : "blastp/blastdbcmd not on PATH -- cannot exercise filtering; "
	   . "source the BLAST runtime first");
    done_testing();
    exit 0;
}

#
# Find a small protein database carrying several taxa: with only one taxon a
# filter cannot be observed to do anything.
#
my($db, @taxa);
PATH: for my $p (@paths)
{
    for my $tl (sort glob("$p/*/*.features.faa.taxlist"))
    {
	my @t = do { open(my $fh, "<", $tl) or next; map { chomp; $_ } <$fh> };
	next if @t < 3;
	(my $base = $tl) =~ s/\.taxlist$//;
	my $sz = -s "$base.psq";
	next if !$sz || $sz > 200_000_000;
	($db, @taxa) = ($base, @t);
	last PATH;
    }
}

if (!$db)
{
    done_testing();
    exit 0;
}
diag("using $db (" . scalar(@taxa) . " taxa)");

my $dir = tempdir(CLEANUP => 1);
system("blastdbcmd -db $db -entry all -outfmt '%a' 2>/dev/null | sed -n 50p > $dir/acc");
chomp(my $acc = do { open(my $fh, "<", "$dir/acc"); scalar <$fh> } // '');
$acc or do { done_testing(); exit 0 };
system("blastdbcmd -db $db -entry '$acc' > $dir/q.faa 2>/dev/null");

sub staxids
{
    my(@extra) = @_;
    my @out = `blastp -db $db -query $dir/q.faa -outfmt '6 staxid' -max_target_seqs 200 -num_threads 4 @extra 2>/dev/null`;
    chomp @out;
    my %u; $u{$_} = 1 for grep { length } @out;
    return sort keys %u;
}

my @all = staxids();
cmp_ok(scalar(@all), '>', 1, "unfiltered search spans more than one taxon")
    or diag("database is not heterogeneous enough to detect a broken filter");

my($one) = grep { my $t = $_; grep { $_ eq $t } @all } @taxa;
SKIP: {
    skip "no taxon with hits", 2 unless $one;

    my @got = staxids("-taxids", $one);
    is_deeply(\@got, [$one], "-taxids returns only the requested taxon")
	or diag("got: @got -- filter is being ignored");

    my @three = grep { my $t = $_; grep { $_ eq $t } @all } @taxa[0..($#taxa < 5 ? $#taxa : 5)];
    skip "not enough taxa with hits", 1 if @three < 2;
    open(my $fh, ">", "$dir/tl") or die; print $fh "$_\n" for @three; close($fh);
    my @gotl = staxids("-taxidlist", "$dir/tl");
    my %want = map { $_ => 1 } @three;
    my @leak = grep { !$want{$_} } @gotl;
    is_deeply(\@leak, [], "-taxidlist returns nothing outside the requested taxa")
	or diag("leaked: @leak");
}

done_testing();
