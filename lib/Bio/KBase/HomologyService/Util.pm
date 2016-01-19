package Bio::KBase::HomologyService::Util;

use strict;
use File::Temp;
use File::Basename;
use Data::Dumper;
use base 'Class::Accessor';
use IPC::Run 'run';
use JSON::XS;
use DB_File;

__PACKAGE__->mk_accessors(qw(blast_db_genomes impl json));

my %blast_command_subject_type = (blastp => 'a',
				  blastn => 'd',
				  blastx => 'a',
				  tblastn => 'd',
				  tblastx => 'd');
my %blast_command_exe_name = (blastp => 'blastp',
			      blastn => 'blastn',
			      blastx => 'blastx',
			      tblastn => 'tblastn',
			      tblastx => 'tblastx');

sub new
{
    my($class, $impl) = @_;

    my $self = {
	impl => $impl,
	blast_db_genomes => $impl->{_blast_db_genomes},
	json => JSON::XS->new->pretty(1),
    };
    return bless $self, $class;
}

sub find_genome_db
{
    my($self, $genome, $db_type, $type) = @_;

    my $base = $self->blast_db_genomes . "/$genome";
    my $file;
    my $check;
    if ($db_type =~ /^a/i)
    {
	$file = "$base.faa";
	$check = "$file.pin";
    }
    elsif ($db_type =~ /^d/i && $type =~ /^c/i)
    {
	$file = "$base.fna";
	$check = "$file.nin";
    }
    elsif ($db_type =~ /^d/i && $type =~ /^f/i)
    {
	$file = "$base.ffn";
	$check = "$file.nin";
    }
    else
    {
	die "find_genome_db: Invalid combination of db_type=$db_type and type=$type";
    }


    if (! -f $check)
    {
	die "find_genome_db: Could not find index file $check";
    }

    return $file;
}

#
# Build a blast alias database and return a File::Temp referring to it.
# If we just have a single genome we'll return the right db file. Don't
# delete the output from here (File::Temp will handle the deletion).
# 
sub build_alias_database
{
    my($self, $subj_genomes, $subj_db_type, $subj_type) = @_;

    if (@$subj_genomes == 0)
    {
	return;
    }
    elsif (@$subj_genomes == 1)
    {
	return $self->find_genome_db($subj_genomes->[0], $subj_db_type, $subj_type);
    }

    my @db_files;
    for my $g (@$subj_genomes)
    {
	my $f = $self->find_genome_db($g, $subj_db_type, $subj_type);
	push(@db_files, $f);
    }
    my $build_db;
    print STDERR Dumper(\@db_files);
    my $db_file = File::Temp->new(UNLINK => 1);
    close($db_file);
    $build_db = ["blastdb_aliastool",
		 "-dblist", join(" ", @db_files),
		 "-title", join(" ", @$subj_genomes),
		 "-dbtype", (($subj_db_type =~ /^a/i) ? 'prot' : 'nucl'),
			 "-out", $db_file];
    my $ok = run($build_db);
    $ok or die "Error running database build @$build_db\n";
    print STDERR "Built db $db_file\n";
    return $db_file;
}

sub construct_blast_command
{
    my($self, $program, $evalue_cutoff, $max_hits, $min_coverage) = @_;

    my $exe = $blast_command_exe_name{$program};

    if (!$exe)
    {
	warn "No blast executable found for $program\n";
	return undef;
    }

    my $suffix = $self->impl->{_blast_program_suffix};
    $exe .= $suffix if $suffix;
    my $prefix = $self->impl->{_blast_program_prefix};
    $exe = $prefix . $exe if $prefix;

    my @cmd = ($exe);
    
    if ($evalue_cutoff)
    {
	push(@cmd, "-evalue", $evalue_cutoff);
    }
    if ($max_hits)
    {
	push(@cmd, "-max_target_seqs", $max_hits);
#	push(@cmd, "-max_hsps", $max_hits);
    }
    if ($min_coverage)
    {
	push(@cmd, "-qcov_hsp_perc", $min_coverage);
    }

    if ($self->impl->{_blast_threads})
    {
	push(@cmd, "-num_threads", $self->impl->{_blast_threads});
    }

    return @cmd;
}

sub blast_fasta_to_genomes
{
    my ($self, $fasta_data, $program, $genomes, $subj_type, $evalue_cutoff, $max_hits, $min_coverage) = @_;

    my @cmd = $self->construct_blast_command($program, $evalue_cutoff, $max_hits, $min_coverage);

    my $subj_db_type = $blast_command_subject_type{$program};

    if (!$subj_db_type || !@cmd)
    {
	die "blast_fasta_to_genomes: Couldn't find blast program $program";
    }
    
    my $db_file = $self->build_alias_database($genomes, $subj_db_type, $subj_type);
    $db_file or die "Couldn't find db file for @$genomes with subj_db_type=$subj_db_type and subj_type=$subj_type";

    my $fmt = 15;		# JSON single file

    push(@cmd, "-db", "$db_file");
    push(@cmd, "-outfmt", $fmt);
    
    my $json;
    my $err;
    my $ok = run(\@cmd, "<", \$fasta_data, ">", \$json, "2>", \$err);

    if (!$ok)
    {
	die "Blast failed @cmd: $err";
    }

    my $doc = eval { $self->json->decode($json) };
    if ($@)
    {
	die "json parse failed: $@\nfor cmd @cmd\n$json\n";
    }
    $doc = $doc->{BlastOutput2};
    $doc or die "JSON output didn't have expected key BlastOutput2";

    my $metadata = {};
    for my $report (@$doc)
    {
	my $search = $report->{report}->{results}->{search};
	$search->{query_id} =~ s/^gnl\|//;
	if ($search->{query_id} =~ /^Query_\d+/)
	{
	    my($xid) = $search->{query_title} =~ /^(\S+)/;
	    $search->{query_id} = $xid if $xid;
	}
	for my $res (@{$search->{hits}})
	{
	    for my $desc (@{$res->{description}})
	    {
		my $md;
		if ($desc->{id} =~ /^gnl\|BL_ORD/)
		{
		    if ($desc->{title} =~ /^(\S+)\s+(.*)\s{3}\[(.*?)(\s*\|\s*(\S+))?\]\s*$/)
		    {
			$desc->{id} = $1;
			$md->{function} = $2;
			$md->{genome_name} = $3;
			$md->{genome_id} = $5 if $5;
		    }
		    elsif ($desc->{title} =~ /^(\S+)\s+\[(.*?)(\s*\|\s*(\S+))?\]\s*$/)
		    {
			$desc->{id} = $1;
			$md->{genome_name} = $2;
			$md->{genome_id} = $4 if $4;
		    }
		}
		else
		{
		    $desc->{id} =~ s/^gnl\|//;
		    if ($desc->{title} =~ /^\s*(.*)\s{3}\[(.*?)(\s*\|\s*(\S+))?\]\s*$/)
		    {
			$md->{function} = $1;
			$md->{genome_name} = $2;
			$md->{genome_id} = $4 if $4;
		    }
		    elsif ($desc->{title} =~ /^\s*\[(.*?)(\s*\|\s*(\S+))?\]\s*$/)
		    {
			$md->{genome_name} = $1;
			$md->{genome_id} = $3 if $3;
		    }
		}
		if ($desc->{id} =~ /^(kb\|g\.\d+)/)
		{
		    $md->{genome_id} = $1;
		}
		$metadata->{$desc->{id}} = $md if $md;
	    }
	}
    }
    return($doc, $metadata);
}

sub enumerate_databases
{
    my($self) = @_;

    my $dir = $self->impl->{_blast_db_databases};

    my %typemap = ('.faa' => 'protein', '.ffn' => 'dna');

    my $res = [];
    
    for my $db (<$dir/*.{faa,ffn}>)
    {
	my $key = basename($db);
	my($name, $path, $suffix) = fileparse($db, '.faa', '.ffn');

	my $descr = {
	    name => $name,
	    key => $key,
	    db_type => $typemap{$suffix},
	    seq_count => 0,
	};
	push(@$res, $descr);
    }
    return $res;
}

sub blast_fasta_to_database
{
    my($self, $fasta_data, $program, $database_key, $evalue_cutoff, $max_hits, $min_coverage) = @_;

    my @cmd = $self->construct_blast_command($program, $evalue_cutoff, $max_hits, $min_coverage);
    if (!@cmd)
    {
	die "blast_fasta_to_genomes: Couldn't find blast program $program";
    }

    my $db_file = $self->impl->{_blast_db_databases} . "/" . $database_key;
    -f $db_file or die "Couldn't find db file $db_file\n";

    my $map_file = $self->impl->{_blast_db_databases} . "/" . $database_key . ".map.btree";
    my %map;

    if (!tie %map, 'DB_File', $map_file, O_RDONLY, 0, $DB_BTREE)
    {
	warn "Could not map $map_file: $!";
    }

    my $fmt = 15;		# JSON single file

    push(@cmd, "-db", "$db_file");
    push(@cmd, "-outfmt", $fmt);
    
    my $json;
    my $err;
    my $ok = run(\@cmd, "<", \$fasta_data, ">", \$json, "2>", \$err);

#     my $ok = run(["cat", "$ENV{HOME}/nr.out"], ">", \$json);

    if (!$ok)
    {
	die "Blast failed @cmd: $err";
    }

    my $doc = eval { $self->json->decode($json) };
    if ($@)
    {
	die "json parse failed: $@\nfor cmd @cmd\n$json\n";
    }
    $doc = $doc->{BlastOutput2};
    $doc or die "JSON output didn't have expected key BlastOutput2";

    my $metadata = {};
    my $identical_proteins = {};
    for my $report (@$doc)
    {
	my $search = $report->{report}->{results}->{search};
	$search->{query_id} =~ s/^gnl\|//;
	if ($search->{query_id} =~ /^Query_\d+/)
	{
	    my($xid) = $search->{query_title} =~ /^(\S+)/;
	    $search->{query_id} = $xid if $xid;
	}
	for my $res (@{$search->{hits}})
	{
	    for my $desc (@{$res->{description}})
	    {
		my $md;

		#
		# We expect to get our form of the title since we're processing one of our MD5 NR databases.
		#

		# "title": "md5|b4dd51958f3c7a7a21e0f222dcbbd764|kb|g.1053.peg.3023   Threonine synthase (EC 4.2.3.1)   [Escherichia coli 101-1]   (1067 matches)"

		if ($desc->{title} =~ /^md5\|([a-z0-9]{32})\|((kb\|g\.\d+)\S+)\s{3}(.+)\s{3}\[(.*?)\]\s{3}\((\d+)\s+matches/)
		{
		    my $md5 = $1;
		    my $rep = $2;
		    my $rep_genome_id = $3;
		    my $rep_fn = $4;
		    my $rep_genome = $5;
		    my $matches = $6;

		    $desc->{id} = $rep;

		    $md->{function} = $rep_fn;
		    $md->{genome_name} = $rep_genome;
		    $md->{genome_id} = $rep_genome_id;
		    $md->{md5} = $md5;
		    $md->{match_count} = 0 + $matches;

		    $metadata->{$desc->{id}} = $md if $md;

		    my $identical = $map{$md5};
		    my @iden = split(/$;/, $identical);
		    for my $one (@iden)
		    {
			my($xid, $xfunc, $xgenome) = split(/\t/, $one);
			next if $xid eq $desc->{id};
			my($xgn) = $xid =~ /^(kb\|g\.\d+)/;
			push(@{$identical_proteins->{$desc->{id}}}, [$xid, {
			    function => $xfunc,
			    genome_name => $xgenome,
			    genome_id => $xgn,
			}]);
		    }
		}
	    }
	}
    }
    return($doc, $metadata, $identical_proteins);
}


