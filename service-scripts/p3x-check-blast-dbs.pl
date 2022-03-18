
# Check for existence of formatted blast dbs in
# /vol/patric3/downloads/blastdb
#
# Inputs are four files with lists of existing data
# eg /tmp/blast.faa  /tmp/blast.ffn  /tmp/blast.fna  /tmp/blast.frn
#
# These came from something like
# cd /vol/patric3/downloads/blastdb
# for i in faa ffn fna frn; do echo $i; /bin/ls -f $i > /tmp/blast.$i; done
# /bin/ls -f /vol/patric3/downloads/genomes/ > /tmp/blast.genomes
#
# Scan those files to mark which formats we have blast databases for.
#
# Then look up in data API for all public genomes. Check if we have
# blast db.
#
# If not, also check for download files. Make lists of genomes
# that we need download for and that we need BLAST dbs created for.
#

use strict;
use Data::Dumper;
use P3DataAPI;

@ARGV == 5 or die "Usage: $0 list.faa list.ffn list.fna list.frn list.genomes\n";

my @types = qw(faa ffn fna frn dl);

my @base_paths = qw(/disks/tmp/downloads /vol/patric3/downloads);

my $to_run = "blast.to.run";
print STDERR "Writing to run to $to_run\n";
open(TO_RUN, ">", $to_run) or die "Cannot write $to_run: $!";

#
# $have_blast{$type}->{$genome} = 1 if present
#
my %have_blast;
$have_blast{$_} = {} foreach @types;

#
# $missing->{$genome_id}->{$type} = 1 if missing
#
my %missing;

for my $tid (0..$#types)
{
    my $type = $types[$tid];
    my $file = $ARGV[$tid];

    print STDERR "$type $file\n";

    my $bhash = $have_blast{$type};

    open(F, "<", $file) or die "Cannot open $file: $!";
    while (<F>)
    {
	if (/^(\d+\.\d+)/)
	{
	    $bhash->{$1} = 1;
	}
    }
}

#
# Now we can process the genomes.
#

my $api = P3DataAPI->new;
my $genomes = $api->query_cb('genome', \&check, ['eq', 'public', 'true'], [select => 'genome_id,date_inserted']);

sub check
{
    my ($data) = @_;
    print STDERR "data\n";
    for my $ent (@$data)
    {
	my $id= $ent->{genome_id};
	my $date = $ent->{date_inserted};
	for my $type (@types)
	{
	    if (!$have_blast{$type}->{$id})
	    {
		$missing{$id}->{$type} = 1;
	    }
	}
	my @m = keys %{$missing{$id}};
	if (@m)
	{
	    print "Missing\t$id\t$date\t@m\n";

	    if (!$missing{$id}->{dl})
	    {
		for my $m (@m)
		{
		    my $fn;
		    if ($m eq 'fna')
		    {
			$fn = "$id.fna";
		    }
		    else
		    {
			$fn = "$id.PATRIC.$m";
		    }
		    for my $base (@base_paths)
		    {
			my $path = "$base/genomes/$id/$fn";
			if (-f $path)
			{
			    my $sz = -s _;
			    print TO_RUN "$id\t$m\t$sz\t$path\n";
			    last;
			}
		    }
			
		}
	    }
	}
    }
    return 1;
}
