#
# Remove genome_counts from databases.json and save without extra space.
# acts as a filter.
#

use strict;
use JSON::XS;
use File::Slurp;

my $in = read_file(\*STDIN);
my $dat = decode_json($in);

for my $e (@$dat)
{
    for my $x (@{$e->{db_list}})
    {
	delete $x->{genome_counts};
    }
}

print encode_json($dat);
