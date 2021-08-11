
#
# Homology search application
#
#
#

use Bio::KBase::AppService::AppScript;
use Bio::P3::HomologySearch::HomologySearch;
use strict;
use Data::Dumper;

my $search = Bio::P3::HomologySearch::HomologySearch->new();

my $script = Bio::KBase::AppService::AppScript->new(sub { $search->process(@_); }, sub { $search->preflight(@_); });

my $rc = $script->run(\@ARGV);

exit $rc;
