
use IPC::Run qw(run);
use JSON::XS;
use strict;
use File::Path qw(make_path);

use LPTScheduler;

my $out_dir = "/disks/tmp/all.dbs.2";
make_path($out_dir);

my $n_procs = 12;
my $n_build_procs = 2;

my $dbtype = 'aa';
my $ftype = 'features';

my $sched = LPTScheduler->new($n_procs);

my $json = JSON::XS->new->pretty->canonical;

my %genera = (
	      'Acinetobacter' => 469,
	      'Francisella' => 262,
	      'Bacillus' => 1386,
	      'Helicobacter' => 209,
	      'Bartonella' => 773,
	      'Listeria' => 1637,
	      'Borreliella' => 64895,
	      'Mycobacterium' => 1763,
	      'Brucella' => 234,
	      'Pseudomonas' => 286,
	      'Burkholderia' => 32008,
	      'Rickettsia' => 780,
	      'Campylobacter' => 194,
	      'Salmonella' => 590,
	      'Chlamydia' => 810,
	      'Shigella' => 620,
	      'Clostridium' => 1485,
	      'Staphylococcus' => 1279,
	      'Coxiella' => 776,
	      'Streptococcus' => 1301,
	      'Ehrlichia' => 943,
	      'Vibrio' => 662,
	      'Escherichia' => 561,
	      'Yersinia' => 629,
	      #
	      # Plus these larger genera
	      #
	      'Neisseria' => 482,
	      'Klebsiella' => 570,
	      'Enterococcus' => 1350,
	      );


#
# Build a database for each taxon; also collect the list of taxa to build
# an exclude list for the "Other" set.
#

my @all;
my @exclude;

for my $genus (sort keys %genera)
{
    my $tax = $genera{$genus};
    my $dat;
    my $create_params = ["--tax", $tax];
    my @cmd = ("p3x-create-blast-db",
	       "--dump",
	       @$create_params);
#    print STDERR "@cmd\n";
    run(\@cmd, ">", \$dat);
    my $glist;
    eval {
	$glist = $json->decode($dat);
    };
    if ($@)
    {
	die "Failed to parse for $tax: $@\n$dat\n";
    }
    my $sz = 0;
    $sz += $_->{genome_length} foreach @$glist;
    $sz /= 1e6;
    my $item = { genus => $genus, tax => $tax, size => $sz, list => $glist, params => $create_params };
    push(@all, $item);
    my $n = @$glist;
    print STDERR "$genus\t$tax\t$n\t$sz\n";
    push(@exclude, $tax);

    $sched->add_work($item, $sz);
}

my $dat;

my $create_params = [(map { ("--exclude", $_) } @exclude),
		     "--tax", 2];

my @cmd = ("p3x-create-blast-db",
	   "--dump",
	   @$create_params,
	   );
# print STDERR "@cmd\n";
run(\@cmd, ">", \$dat);
my $glist = $json->decode($dat);
my $sz = 0;
$sz += $_->{genome_length} foreach @$glist;
$sz /= 1e6;
my $item = { genus => 'OTHER', tax => undef, size => $sz, list => $glist, params => $create_params };
push(@all, $item);
$sched->add_work($item, $sz);
my $n = @$glist;

print STDERR "OTHER\t\t$n\t$sz\n";
# print $json->encode(\@all);

$sched->run(sub {}, sub {
    my($global, $item) = @_;

    my($genus, $tax, $size, $list, $params) = @$item{qw(genus tax size list params)};

    print STDERR "Run $$ $genus $tax @$params\n";

    my $db = "$out_dir/$genus";
    my @cmd = ("p3x-create-blast-db",
	       @$params,
	       "--title", $genus,
	       "--parallel", $n_build_procs,
	       $dbtype, $ftype, $db);
    print "@cmd\n";
    my $ok = run(\@cmd,
		">", "log.$genus.out",
		"2>", "log.$genus.err",
		);
    $ok or die "Failure running @cmd: $! $?\n";
});
