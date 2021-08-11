#
# Using the contents of p3x-check-blast-dbs.pl to-run file,
# start running blast makedbs in parallel.
#

use strict;
use LPTScheduler;
use Data::Dumper;
use File::Basename;

my $parallel = 24;

my $sched = LPTScheduler->new($parallel);

my %dbtype = (faa => 'prot',
	      ffn => 'nucl',
	      frn => 'nucl',
	      fna => 'nucl');

while (<>)
{
    chomp;
    my($genome, $type, $size, $path) = split(/\t/);
    if ($size > 0)
    {
	$sched->add_work([$genome, $type, $path], $size);
    }
}

$sched->run(sub {}, sub {
    my($global, $work) = @_;
    my($genome, $type, $path) = @$work;

    my $base = basename($path);
    my $out = "/disks/tmp/blastdb/$type/$base";
    
    my @blast = ("makeblastdb",
		 "-blastdb_version", 4,
		 "-dbtype", $dbtype{$type},
		 "-in", $path,
		 "-out", $out);
    my $rc = system(@blast);
    if ($rc != 0)
    {
	print STDERR "Failed @$work: $rc\n";
    }
});
