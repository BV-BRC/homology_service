use strict;
use Getopt::Long::Descriptive;
use gjoseqlib;
use LPTScheduler;
use Data::Dumper;
use Time::HiRes qw(gettimeofday);
use IPC::Run qw(run);

#
# read fasta file
# sort by length
# run probes every N aa
#

@ARGV == 2 or die "usage: $0 fasta-file db-file\n";
my $fa_file = shift;
my $db_file = shift;

my $threads_per_blast = 4;
my $cpus = 20;
my $threads = int($cpus / $threads_per_blast);

my $sched = LPTScheduler->new($threads);

my @seqs;

open(FA, "<", $fa_file) or die "Cannot open $fa_file: $!";

while (my($id, $def, $seq) = read_next_fasta_seq(\*FA))
{
    push(@seqs, [$id, $seq]);
}

@seqs = sort { length($a->[1]) <=> length($b->[1]) } @seqs;

my $min_int = 10;

my $last_len = -($min_int + 1);

for my $ent (@seqs)
{
    my($id, $seq) = @$ent;
    my $len = length($seq);
    if ($len > $last_len + $min_int)
    {
	# print "Use $id $len\n";
	$last_len = $len;
	my @rpt;
	while ($seq =~ /(([A-Z])\2\2+)/g)
	{
	    push(@rpt, $1);
	}
	# print "$id @rpt\n";
	$sched->add_work([$id, $seq], $len);
	$sched->add_work([$id, $seq], $len);
    }
}

$sched->run(sub {}, sub {
    my($glob, $item) = @_;
    my($id, $seq) = @$item;
    print STDERR "run $id " . length($seq) . "\n";

    my $inp = ">$id\n$seq\n";
    my $start = gettimeofday;
    my @cmd = ("blastp",
	       "-query", "-",
	       "-db", $db_file,
	       "-outfmt", 6,
	       "-num_threads", $threads_per_blast,
	       "-out", "/dev/null",
	       );

    my $ok = run(\@cmd, "<", \$inp);
    my $end = gettimeofday;
    print join("\t", $id, length($seq), $end - $start), "\n";
    
    return 1;
});

