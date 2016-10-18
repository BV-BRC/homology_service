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

sub data_api
{
    my($self) = @_;
    return $self->{impl}->{_data_api};
}

sub compute_db_filename
{
    my($self, $genome, $db_type, $type) = @_;

    my @candidates;
    if ($db_type =~ /^a/i)
    {
	push(@candidates,
	     ["PATRIC.faa", "PATRIC.faa.pin"],
	     ["RefSeq.faa", "RefSeq.faa.pin"]);
    }
    elsif ($db_type =~ /^d/i && $type =~ /^c/i)
    {
	push(@candidates,
	     ["fna", "fna.nin"]);
    }
    elsif ($db_type =~ /^d/i && $type =~ /^f/i)
    {
	push(@candidates,
	     ["PATRIC.ffn", "PATRIC.ffn.nin"],
	     ["RefSeq.ffn", "RefSeq.ffn.nin"],
	     ["PATRIC.frn", "PATRIC.frn.nin"],
	     ["RefSeq.frn", "RefSeq.frn.nin"]);
    }
    else
    {
	die "find_genome_db: Invalid combination of db_type=$db_type and type=$type";
    }

    return @candidates;
}

sub find_genome_db
{
    my($self, $genome, $db_type, $type) = @_;

    my $base = $self->blast_db_genomes . "/$genome";

    my @candidates = $self->compute_db_filename($genome, $db_type, $type);

    my @paths;
    for my $cand (@candidates)
    {
	my($file, $check) = @$cand;
	my $path = "$base.$file";
	my $check_path = "$base.$check";
	if (-f $check_path)
	{
	    push(@paths, $path);
	}
    }
    return @paths;
}

sub find_private_genome_db
{
    my($self, $genome, $owner, $db_type, $type) = @_;

    my $base = $self->blast_db_private_genomes . "$owner/$genome";

    return ();
    
    my($file, $check) = $self->compute_db_filename($genome, $db_type, $type);

    if (! -f $check)
    {
	$self->create_private_genome_db($genome, $owner, $db_type, $type, $base, $file);
    }

    return $file;
}


sub create_private_genome_db
{
    my($self, $genome, $owner, $db_type, $type, $base, $file) = @_;

    
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

    #
    # We partition the genome list into private and public sets.
    #


    my @todo = @$subj_genomes;
    my $public = [];
    my $private = [];
    while (@todo)
    {
	my @block = splice(@todo, 0, 500);
	my $glist = join(",", @block);

	my @res = $self->data_api->query('genome',
					 ['in', 'genome_id', "($glist)"],
					 ['select', 'genome_id', 'public', 'owner']);

	for my $ent (@res)
	{
	    push(@{$ent->{public} ? $public : $private}, [$ent->{genome_id}, $ent->{owner}]);
	}
    }
    # print Dumper($public, $private);

    my @public_files = map { $self->find_genome_db($_->[0], $subj_db_type, $subj_type) } @$public;
    my @private_files = map { $self->find_private_genome_db($_->[0], $_->[1], $subj_db_type, $subj_type) } @$private;

    my @db_files = (@public_files, @private_files);
    
    if (@db_files == 1)
    {
	return $db_files[0];
    }

    my $build_db;
    # print STDERR Dumper(\@db_files);
    my $db_file = File::Temp->new(UNLINK => 0);
    close($db_file);
    my $dblist_tmp;
    my @dblist_arg;
    if (@db_files > 10)
    {
	$dblist_tmp = File::Temp->new();
	print $dblist_tmp "$_\n" foreach @db_files;
	close($dblist_tmp);
	@dblist_arg = ("-dblist_file", "$dblist_tmp");
    }
    else
    {
	@dblist_arg = ("-dblist", join(" ", @db_files));
    }

    my $title;
    if (@$subj_genomes > 10)
    {
	my $nmore = @$subj_genomes - 10;
	$title = "DB of " . join(",", @$subj_genomes[0..9]) . " and $nmore more";
    }
    else
    {
	$title = "DB of " . join(",", @$subj_genomes);
    }

    $build_db = ["blastdb_aliastool",
		 @dblist_arg,
		 "-title", $title,
		 "-dbtype", (($subj_db_type =~ /^a/i) ? 'prot' : 'nucl'),
		 "-out", "$db_file"];
    # print STDERR "@$build_db\n";
    my $ok = run($build_db);
    $ok or die "Error running database build @$build_db\n";
    print STDERR "Built db $db_file\n";
    return $db_file;
}

sub construct_blast_command
{
    my($self, $program, $evalue_cutoff, $max_hits, $min_coverage, $is_short) = @_;

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

    #
    # if we have hsort sequences, use the tasks optimized for that.
    #
    if ($program eq 'blastn' && $is_short)
    {
	push(@cmd, "-task", "blastn-short");
    }
    elsif ($program eq 'blastp' && $is_short)
    {
	push(@cmd, "-task", "blastp-short");
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

    my $is_short = $self->test_for_short_query($fasta_data);

    print STDERR "query data len=" . length($fasta_data) . "is_short=$is_short\n";
    print STDERR "!$fasta_data!\n";
    
    my @cmd = $self->construct_blast_command($program, $evalue_cutoff, $max_hits, $min_coverage, $is_short);

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
		my $md = $self->decode_title($desc);
		
		$metadata->{$desc->{id}} = $md if $md;
	    }
	}
    }
    return($doc, $metadata);
}

sub decode_title
{
    my($self, $desc) = @_;
	
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
    return $md;
}

sub enumerate_databases
{
    my($self) = @_;

    my $dir = $self->impl->{_blast_db_databases};

    my %typemap = ('.faa' => 'protein',
		   '.ffn' => 'dna',
		   '.faa.pal' => 'protein',
		   '.ffn.nal' => 'dna',
		   '.fna.nal' => 'dna',
		  );

    my $res = [];
    
    for my $db (<$dir/*.{faa,ffn,faa.pal,ffn.nal,fna.nal}>)
    {
	my $key = basename($db);
	my($name, $path, $suffix) = fileparse($db, qw(.faa .ffn .faa.pal  .ffn.nal .fna.nal));

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

    my $is_short = $self->test_for_short_query($fasta_data);

    my @cmd = $self->construct_blast_command($program, $evalue_cutoff, $max_hits, $min_coverage, $is_short);
    if (!@cmd)
    {
	die "blast_fasta_to_database: Couldn't find blast program $program";
    }

    #
    # Sanity check.
    #
    $database_key = basename($database_key);

    my $db_file = $self->impl->{_blast_db_databases} . "/" . $database_key;
    -f $db_file or warn "Couldn't find db file $db_file\n";

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
		my $md = $self->decode_title($desc);
		$metadata->{$desc->{id}} = $md if $md;
	    }
	}
    }
    return($doc, $metadata, $identical_proteins);
}

sub test_for_short_query
{
    my($self, $fasta_data) = @_;
    my $thresh = 30;

    my $cur;

    while ($fasta_data =~ /^(.*)/mg)
    {
        my $l = $1;
	    
        if ($l =~ /^>/)
        {
	    return 0 if $cur > $thresh;
            $cur = 0;
        }
        else
        {
            $cur += ($l =~ tr/[a-z][A-Z]//);
        }
    }
    return $cur > $thresh ? 0 : 1;
    
}
