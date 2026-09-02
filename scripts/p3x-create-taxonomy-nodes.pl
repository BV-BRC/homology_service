
use strict;
use POSIX;

my $now = strftime("%Y-%m-%d", localtime);
print "$now\n";
my $dir = "taxdump-$now";

-d $dir and die "$dir already exists\n";
mkdir($dir) or die "cannot mkdir $dir: $!";

chdir($dir);

my $tar = "taxdump-$now.tar.gz";
my $url = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz";

my $rc = system("curl", "-L", "-o", $tar, $url);
$rc == 0 or die "error downloading $url to $tar\n";

$rc = system("tar", "vxzfp", $tar, "nodes.dmp");
$rc == 0 or die "Error unpacking $tar\n";

#
# Create the CSV of node data.
#

open(N, "<", "nodes.dmp") or die "Cannot open nodes.dmp: $!";
open(T, ">", "nodes.tsv") or die "Cannot write nodes.tsv: $!";
while (<N>)
{
    chomp;
    my @f = split(/\t\|\t/);
    print T join(",", @f[0,1,2]), "\n";
}
close(N);
close(T);


    
