package Bio::KBase::HomologyService::Util;

use Bio::KBase::HomologyService::Bench;

use strict;
use File::Temp;
use File::Basename;
use Data::Dumper;
use base 'Class::Accessor';
use IPC::Run 'run';
use JSON::XS;
use DB_File;
use LWP::UserAgent;

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
	$file = "$base.PATRIC.faa";
	$check = "$file.pin";
    }
    elsif ($db_type =~ /^d/i && $type =~ /^c/i)
    {
	$file = "$base.fna";
	$check = "$file.nin";
    }
    elsif ($db_type =~ /^d/i && $type =~ /^f/i)
    {
	$file = "$base.PATRIC.ffn";
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
    my $db_file = File::Temp->new(UNLINK => 0);
    close($db_file);
    $build_db = ["blastdb_aliastool",
		 "-dblist", join(" ", @db_files),
		 "-title", join(" ", @$subj_genomes),
		 "-dbtype", (($subj_db_type =~ /^a/i) ? 'prot' : 'nucl'),
			 "-out", "$db_file"];
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

sub blast_fasta_to_taxon
{
    my ($self, $fasta_data, $program, $taxon_id, $subj_type, $evalue_cutoff, $max_hits, $min_coverage) = @_;

    #
    # Query the data API to determine the set of genomes that derive from taxon_id.
    #

    $taxon_id =~ /(\d+)/ or die "Invalid taxon ID";
    $taxon_id = $1;
    my $ua = LWP::UserAgent->new;
    my $res = $ua->post("https://www.patricbrc.org/api/genome",
		        {
			    fl => 'genome_id,genome_name',
			    q => "taxon_lineage_ids:$taxon_id",
			    start => 0,
			    rows => 1000
			},
			"Content-type" => "application/solrquery+x-www-form-urlencoded",
			"Accept", "application/solr+json");
    die "failure retrieving taxon data: " . $res->status_line . "\n" . $res->content if !$res->is_success;
    my $doc;
    eval { $doc = decode_json($res->content); };
    $doc or die "could not decode taxon data: $@";

    my $subj_db_type = $blast_command_subject_type{$program};
    
    my @genomes;
    for my $g (@{$doc->{response}->{docs}})
    {
	my $gid = $g->{genome_id};
	my $gname = $g->{genome_name};

	#
	# Try to find the db file and emit a warning if it is missing.
	#
	eval {
	    my $d = $self->find_genome_db($gid, $subj_db_type, $subj_type);
	    print "Found $d for $gid $gname\n";
	    push(@genomes, $gid);
	};
	if ($@)
	{
	    print STDERR "No blast database found for $gid $gname\n";
	    print STDERR $@;
	}
    }

    die "No blast databases found for $taxon_id" if @genomes == 0;

    return $self->blast_fasta_to_genomes($fasta_data, $program, \@genomes, $subj_type, $evalue_cutoff, $max_hits, $min_coverage);
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

    my $bench = Bio::KBase::HomologyService::Bench->new();
    print STDERR "cmd=@cmd\n";
    my $ok = run(\@cmd, "<", \$fasta_data, ">", \$json, "2>", \$err);
    $bench->finish();

    {
	my $ctx = $Bio::KBase::HomologyService::Service::CallContext;
    
	my $len = length($fasta_data);
	my $stat = $bench->stats_text;
	my $g = join(",", @$genomes);
	my $rsize = length($json);
	my $logstr = "blast_fasta_to_genomes len=$len $stat program=$program genomes=$g ok=$ok" .
		       ($ok ? " result_size=$rsize" : " err=$err");
	if ($ctx)
	{
	    $ctx->log_info($logstr);
	}
	else
	{
	    print STDERR $logstr, "\n";
	}
	    
    }

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
		    if ($desc->{title} =~ /^(\S+)\s{3}(.*)\s{3}(.*)$/)
		    {
			my $id = $1;
			my $fun = $2;
			my $ginfo = $3;
			$id =~ s/^((fig\|\d+\.\d+.[^.]+\.\d+)|(accn\|[^|]+))//;
			my $fid = $1;
			print "id=$id\n";
			$id =~ s/^\|//;
			$id =~ s/\|$//;
			my @rest = split(/\|/, $id);
			my $locus;
			my $alt;
			if (@rest == 2)
			{
			    ($locus, $alt) = @rest;
			    $md->{locus_tag} = $locus;
			    $md->{alt_locus_tag} = $alt;
			}
			elsif (@rest == 1)
			{
			    $alt = $rest[0];
			    $md->{alt_locus_tag} = $alt;
			}
			if ($ginfo =~ /^\[(.*) \| (\d+\.\d+)/)
			{
			    $md->{genome_name} = $1;
			    $md->{genome_id} = $2;
			}
			$desc->{id} = $fid;
			$md->{function} = $fun;
		    }
		    elsif ($desc->{title} =~ /^(\S+)\s+(.*)\s{3}\[(.*?)(\s*\|\s*(\S+))?\]\s*$/)
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
	die "blast_fasta_to_database: Couldn't find blast program $program";
    }

    my $db_file = $self->impl->{_blast_db_databases} . "/" . $database_key;
    -f $db_file or die "Couldn't find db file $db_file\n";

    my $map_file = $self->impl->{_blast_db_databases} . "/" . $database_key . ".map.btree";
    my %map;

    if (!tie %map, 'DB_File', $map_file, O_RDONLY, 0, $DB_BTREE)
    {
	# warn "Could not map $map_file: $!";
    }

    my $fmt = 15;		# JSON single file

    push(@cmd, "-db", "$db_file");
    push(@cmd, "-outfmt", $fmt);
    
    my $json;
    my $err;
    my $bench = Bio::KBase::HomologyService::Bench->new();
    my $ok = run(\@cmd, "<", \$fasta_data, ">", \$json, "2>", \$err);

    $bench->finish();

    {
	my $ctx = $Bio::KBase::HomologyService::Service::CallContext;
    
	my $len = length($fasta_data);
	my $stat = $bench->stats_text;
	my $rsize = length($json);
	my $logstr = "blast_fasta_to_database len=$len db_file=$db_file $stat program=$program ok=$ok" .
		       ($ok ? " result_size=$rsize" : " err=$err");
	if ($ctx)
	{
	    $ctx->log_info($logstr);
	}
	else
	{
	    print STDERR $logstr, "\n";
	}
	    
    }


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

		if ($desc->{title} =~ /^md5\|([a-z0-9]{32})\|((kb\|g\.\d+)\S+)\s{3}(.*)\s{3}\[(.*?)\]\s{3}\((\d+)\s+matches/)
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
		#
		# Contigs database
		#
		# >kb|g.0.c.1   [Escherichia coli K12]
		#
		elsif ($desc->{title} =~ /^((kb\|g\.\d+)\S+)\s+\[(.*)\]/)
		{
		    my $contig = $1;
		    my $gid = $2;
		    my $gname = $3;

		    $desc->{id} = $contig;

		    $md->{genome_name} = $gname;
		    $md->{genome_id} = $gid;

		    $metadata->{$desc->{id}} = $md if $md;
		}
	    }
	}
    }
    return($doc, $metadata, $identical_proteins);
}


