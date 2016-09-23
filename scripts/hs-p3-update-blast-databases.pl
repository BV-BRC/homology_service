use strict;
use Proc::ParallelLoop;

my $base = "/vol/patric3/downloads";
my $genomes = "$base/genomes";
my $blastdb = "$base/blastdb";

my $par = 12;

opendir(D, $genomes) or die "cannot opendir $genomes: $!";

my @work;
while (my $g = readdir(D))
{
    next unless $g =~ /^\d+\.\d+$/;

    my $g_faa = "$genomes/$g/$g.PATRIC.faa";
    my $b_faa = "$blastdb/$g.PATRIC.faa";

    my $g_ffn = "$genomes/$g/$g.PATRIC.ffn";
    my $b_ffn = "$blastdb/$g.PATRIC.ffn";

    my $g_fna = "$genomes/$g/$g.fna";
    my $b_fna = "$blastdb/$g.fna";

    if (-f $g_faa && ! -f "$b_faa.pin")
    {
	push(@work, ['prot', $g_faa, $b_faa]);
    }
    if (-f $g_ffn && ! -f "$b_ffn.nin")
    {
	push(@work, ['nucl', $g_ffn, $b_ffn]);
    }

    if (-f $g_fna && ! -f "$b_fna.nin")
    {
	push(@work, ['nucl', $g_fna, $b_fna]);
    }
}
my $n = @work;
print "$n items\n";

pareach \@work, sub {
    my($item) = @_;
    my($type, $in, $out) = @$item;
    my @cmd = ("/homes/olson/P3/blast/deployment/services/homology_service/bin/makeblastdb",
	       "-in", $in,
	       "-out", $out,
	       "-dbtype", $type);
    print "@cmd\n";
    my $rc = system(@cmd);
    if ($rc != 0)
    {
	die "Error $rc processing $type $in $out: @cmd\n";
    }
}, { Max_Workers => $par };
