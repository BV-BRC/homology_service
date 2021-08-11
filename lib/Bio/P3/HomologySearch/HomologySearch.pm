
=head1 NAME

HomologySearch - Homology search application.

=head1 DESCRIPTION

The homology search application searches one or more databases for
occurrences of sequences in its input set. The application parameters
contain all the values we describe below.

=head2 Database types

The search databases contain sequences of four different types. Each of these
types is given a unique suffix in the name for the databases containing
data of the given type. This suffix is also used as the key to select a 
database type in the C<db_type> parameter.

=over 4

=item * 

Protein features. These are the amino acid sequences of the protein encoding genes in 
the genome. The suffix is I<faa>.

=item *

DNA features. These are the DNA-based sequences of any genes as called in the genome. 
The suffix is I<ffn>.

=item * 

RNA features. The suffix is I<frn>.

=item *

Genome sequences (contigs). The suffix is I<fna>.

=back

=head2 Specifying input sequences

The type of input (DNA or amino acid) is specified by the C<input_type> parameter. 
It has one of the following values:

=over 4

=item dna

Input is DNA sequences.

=item aa

Input is amino acid sequences.

=back

The source of input is specified by the C<input_source> parameter. It has one of the
following values:

=over 4

=item id_list

Use a set of sequences from the PATRIC database named by the parameter C<input_id_list>. 

=item fasta_data

Use the sequence data from the parameter C<input_fasta_data>. 

=item fasta_file

Use the sequence data from the workspace file at path C<input_fasta_file>.

=back

=head2 Specifying search database

The source of the database to search is specified by the C<db_source> parameter. It
has one of the following values:

=over 4

=item fasta_data

Use the sequence data from the parameter C<db_fasta_data>.

=item fasta_file

Use the sequence data from the workspace file at path C<db_fasta_file>.

=item genome_list

Use the data from the set of PATRIC genome ids specified in C<db_genome_list>.

=item taxon_list

Use the data from genomes in PATRIC at or below the given taxon identifiers. given in C<db_taxon_list>.

=item precomputed_database

Use the precomputed database specified by C<db_precomputed_database>. The list of available databases
may be queried using the HomologySearch service.

=back

=head2 Specifying behavior

The BLAST program to use is defined by the C<blast_program> parameter. Valid values
are C<blastp>, C<blastn>, C<blastx>, C<tblastn>, C<tblastx>. A value is not required; if not
specified, the appropriate program is chosen for the given input and database types.
For DNA to DNA searches, C<blastn> is the default.

=head1 PRECOMPUTED DATABASES

To accelerate common searches, we maintain a number of precomputed BLAST databases:

=over 4

=item * 

Each of the 26 large and/or pathogen specific bacterial genera has a BLAST database containing all genomes in the genus.

=item * 

Each of the 21 viral families has a BLAST database containing all genomes in the family.

=item * 

A database of all the reference and representative genomes in the bacterial and archaeal taxonomies (~6700 genomes).

=back




=head1 METHODS

=cut

package Bio::P3::HomologySearch::HomologySearch;

use P3DataAPI;
use gjoseqlib;
use strict;
use Data::Dumper;
use POSIX;
use Cwd;
use base 'Class::Accessor';
use gjoseqlib;
use JSON::XS;
use Module::Metadata;
use List::Util qw(any none);
use IPC::Run qw(run);
use File::Basename;
use File::Copy qw(copy);
use File::Path qw(make_path remove_tree);
use Template;
use Text::CSV_XS qw(csv);
use Bio::KBase::AppService::ClientExt;
use Bio::KBase::AppService::AppConfig qw(data_api_url);
use Bio::KBase::AppService::FastaParser 'parse_fasta';
use Bio::P3::Workspace::WorkspaceClientExt;

use File::Slurp;

__PACKAGE__->mk_accessors(qw(app app_def params token task_id
			     blast_program blast_params
			     work_dir stage_dir
			     output_base output_folder 
			     contigs app_params  api
			     database_path
			     json
			    ));


#
# $search_type_map{input-type}->{db-type} = valid-blast-programs
#
our %search_type_map =
    (aa => {
	faa => ['blastp'],
	ffn => ['tblastn'],
	frn => ['tblastn'],
	fna => ['tblastn'],
    },
     dna => {
	 faa => ['blastx'],
	 ffn => ['blastn', 'tblastx'],
	 frn => ['blastn', 'tblastx'],
	 fna => ['blastn', 'tblastx'],
     });
				    
our %makeblastdb_dbtype = (faa => "prot",
			   ffn => "nucl",
			   frn => "nucl",
			   fna => "nucl");

sub new
{
    my($class) = @_;

    my $self = {
	api => P3DataAPI->new(data_api_url),
	database_path => ["/vol/bvbrc/blast"],
	json => JSON::XS->new->pretty,
    };

    bless $self, $class;

    $self->load_databases();

    return $self;
}

sub load_databases
{
    my($self) = @_;
    #
    # Open databases file and load.
    #

    my $mpath = Module::Metadata->find_module_by_name(__PACKAGE__);
    $mpath = dirname($mpath);
    my $db_path = "$mpath/databases.json";
    eval {
	my $db_txt = read_file($db_path);
	$db_txt or die "Cannot load databases from $db_path: $!";
	$self->databases(decode_json($db_txt));
    };
    if ($@)
    {
	die "Failure loading database file $db_path: $@";
    }
}

sub databases
{
    my($self, $val) = @_;
    if (defined($val))
    {
	$self->{databases} = $val;
    }
    return wantarray ? @{$self->{databases}} : $self->{databases};
}

sub preflight
{
    my($app, $app_def, $raw_params, $params) = @_;

    my $cpu = 8;

    my $pf = {
	cpu => $cpu,
	memory => "16G",
	runtime => 0,
	storage => 0,
	is_control_task => 0,
    };
    return $pf;
}

sub process
{
    my($self, $app, $app_def, $raw_params, $params) = @_;

    $self->app($app);
    $self->app_def($app_def);
    $self->params($params);
    $self->token($app->token);
    $self->task_id($app->task_id);

    print "Process homology search ", Dumper($app_def, $raw_params, $params);

    my $cwd = getcwd();
    my $work_dir = "$cwd/work";
    my $stage_dir = "$cwd/stage";
    
    -d $work_dir or mkdir $work_dir or die "Cannot mkdir $work_dir: $!";
    -d $stage_dir or mkdir $stage_dir or die "Cannot mkdir $stage_dir: $!";

    $self->work_dir($work_dir);
    $self->stage_dir($stage_dir);

    my $output_base = $self->params->{output_file};
    my $output_folder = $self->app->result_folder();

    $self->output_base($output_base);
    $self->output_folder($output_folder);

    my $input_file = $self->stage_input($params, $stage_dir);
    my($db_file, $blast_params) = $self->stage_database($params, $stage_dir);

    my $blast_program = $self->determine_blast_program($params);
    $self->blast_program($blast_program);

    $self->run_blast($params, $input_file, $db_file, $blast_program, $blast_params, $work_dir);

    # Write output

    my %typemap = (json => 'json',
		   archive => 'unspecified',
		   tsv => 'tsv');

    for my $out (glob("$work_dir/*"))
    {
	my $file = basename($out);
	my($suffix) = $out =~ /\.([^.]+)$/;
	my $type = $typemap{$suffix} // 'unspecified';
	$self->app->workspace->save_file_to_file($out, {}, "$output_folder/$file", $type, 1, 1);
    }
	
}

=head2 stage_input

    my $input_file = $self->stage_input($params, $stage_dir)

Stage the input as defined by C<$params> into C<$stage_dir>.

=cut

sub stage_input
{
    my($self, $params, $stage_dir) = @_;

    my $src = $params->{input_source};
    my $method = "stage_input_$src";
    my $file = $self->$method($params, $stage_dir);
    return $file;
}

sub stage_input_id_list
{
    my($self, $params, $stage_dir) = @_;

    my $ids = $self->app->validate_param_array('input_id_list');
    if (!defined($ids))
    {
	die "Input source defined as id_list, but no input_id_list parameter specified";
    }

    my $seqs;
    if ($params->{input_type} eq 'aa')
    {
	$seqs = $self->api->retrieve_protein_feature_sequence($ids);
    }
    else
    {
	$seqs = $self->api->retrieve_nucleotide_feature_sequence($ids);
    }

    my $file = "$stage_dir/input_$params->{input_type}.fa";
    if (open(my $fh, ">", $file))
    {
	for my $id (@$ids)
	{
	    print_alignment_as_fasta($fh, [$id, undef, $seqs->{$id}]);
	}
	close($fh);
    }
    else
    {
	die "Cannot open $file for writing: $!";
    }
    return $file;
}

sub stage_input_fasta_data
{
    my($self, $params, $stage_dir) = @_;

    my $file = "$stage_dir/input_$params->{input_type}.fa";
    if (open(my $fh, ">", $file))
    {
	print $fh $params->{input_fasta_data};
	close($fh);
    }
    else
    {
	die "Cannot open $file for writing: $!";
    }
    return $file;
}

#
# Stage input from a workspace file
sub stage_input_fasta_file
{
    my($self, $params, $stage_dir) = @_;

    my $file = "$stage_dir/input_$params->{input_type}.fa";
    my $ws = Bio::P3::Workspace::WorkspaceClientExt->new;

    if (open(my $fh, ">", $file))
    {

	$ws->copy_files_to_handles(1, undef, [[$params->{input_fasta_file}, $fh]]);

	close($fh);
    }
    else
    {
	die "Cannot open $file for writing: $!";
    }
    return $file;
}

=head2 stage_database

    ($db_file, $blast_params) $self->stage_database($params, $stage_dir)

Stage the database as defined by C<$params> into C<$stage_dir>.

We return the database file and an initial set of BLAST parameters. The BLAST parameters
may be required if the database specification in the input parameters includes
a subset of taxonomy IDs to be searched.

=cut

sub stage_database
{
    my($self, $params, $stage_dir) = @_;

    my $src = $params->{db_source};
    my $method = "stage_database_$src";

    my $blast_params = [];

    my $db_file = $self->$method($params, $stage_dir, $blast_params);
    return ($db_file, $blast_params);
}

=head2 stage_database_fasta_data

Stage the fasta data from C<$params->{db_fasta_data}> to a local file.

=cut

sub stage_database_fasta_data
{
    my($self, $params, $stage_dir, $blast_params) = @_;

    my $file = "$stage_dir/db_$params->{db_type}.fa";

    if (!exists $params->{db_fasta_data})
    {
	die "Database of type fasta_data requested but parameter db_fasta_data is missing\n";
    }

    open(my $fh, "<", \ $params->{db_fasta_data}) or die "Cannot open filehandle on data: $!";
    open(my $out, ">", $file) or die "Cannot open output file $file: $!";

    $self->read_and_validate_fasta($fh, $params->{db_type}, $out);

    close($fh);
    close($out);

    my $ok = run(["makeblastdb",
		  "-dbtype", $makeblastdb_dbtype{$params->{db_type}},
		  "-in", $file]);
    $ok or die "makeblastdb failed: $!";

    return $file;
}

sub stage_database_fasta_file
{
    my($self, $params, $stage_dir, $blast_params) = @_;
}

sub stage_database_genome_list
{
    my($self, $params, $stage_dir, $blast_params) = @_;
}

sub stage_database_taxon_list
{
    my($self, $params, $stage_dir, $blast_params) = @_;
}

sub stage_database_precomputed_database
{
    my($self, $params, $stage_dir, $blast_params) = @_;

    my $db = $params->{db_precomputed_database};
    defined($db) or die "Precomputed database was select but db_precomputed_database key was missing from parameters";
    my @cands = grep { $_->{name} eq $db } $self->databases;
    if (@cands == 0)
    {
	die "Cannot find precomputed database $db\n";
    }

    for my $cand (@cands)
    {
	for my $db (grep { $_->{type} eq $params->{db_type} } @{$cand->{db_list}})
	{
	    for my $p (@{$self->{database_path}})
	    {
		my $file = "$p/$db->{path}";
		for my $suf (qw(pal nal))
		{
		    my $chk = "$file.$suf";
		    if (-f $chk)
		    {
			return ($file);
		    }
		}
	    }
	}
    }
    return undef;
}

sub read_and_validate_fasta
{
    my($self, $in_fh, $type, $out_fh) = @_;
    
    my %val_re = (faa => qr/^[a-z]+$/i,
		  ffn => qr/^[actg]+$/i,
		  frn => qr/^[actg]+$/i,
		  fna => qr/^[actg]+$/i);

    my $is_prot = $type eq 'faa';

    my $re = $val_re{$type};
    $re or die "read_and_validate_fasta: invalid type $type\n";

    parse_fasta($in_fh, $out_fh, sub {
	my($id, $seq) = @_;
	if ($seq !~ $re)
	{
	    die "Invalid sequence $id in database file\n" . substr($seq, 0, 50), "\n";
	}
	return 1;
    }, $is_prot);
}

=head2 determine_blast_program
    
Determine the appropriate BLAST program based on the input and db type and
a setting for C<blast_program> from the parameters if present.
    
=cut

sub determine_blast_program
{
    my($self, $params) = @_;
    
    my $valid_programs = $search_type_map{$params->{input_type}}->{$params->{db_type}};
    if (!$valid_programs)
    {
	die "No valid program found for input type '$params->{input_type}' and db type '$params->{db_type}'\n";
    }
    my $prog = $params->{blast_program};
    if ($prog)
    {
	if (none { $_ eq $prog } @$valid_programs)
	{
	    die "$prog is not a valid blast program for input type '$params->{input_type}' and db type '$params->{db_type}'\n";
	}
    }
    else
    {
	$prog = $valid_programs->[0];
    }
    return $prog;
}

sub run_blast
{
    my($self, $params, $input_file, $db_file, $blast_program, $blast_params, $work_dir) = @_;

    #
    # Set output format to blast archive, then translate to desired formats.
    #

    my @opts = @$blast_params;

    my $cpus = $ENV{P3_ALLOCATED_CPU} // 1;

    push(@opts, "-num_threads", $cpus);

    my $out_file = "$work_dir/blast_out.archive";
    my $out_json_file = "$work_dir/blast_out.raw.json";
    my $out_proc_json_file = "$work_dir/blast_out.json";
    my $out_md_file = "$work_dir/blast_out.metadata.json";

    push(@opts,
	 "-query", $input_file,
	 "-db", $db_file,
	 "-outfmt", 11,
	 "-out", $out_file,
	);

    my $ok = run([$blast_program, @opts]);
    $ok or die "Error running $blast_program @opts: $!";

    #
    # We wrote blast archive output. Use blast_formatter to
    # create json output.
    #

    my $ok = run(["blast_formatter",
	       "-archive", $out_file,
	       "-outfmt", 15,
	       "-out", $out_json_file]);
    $ok or die "Error running blast formatter: $!";
    my $txt = read_file($out_json_file);
    my($doc, $metadata) = $self->massage_blast_json($txt);
    write_file($out_proc_json_file, $self->json->encode($doc));
    write_file($out_md_file, $self->json->encode($metadata));
}

sub massage_blast_json
{
    my($self, $json) = @_;
    my $doc = eval { $self->json->decode($json) };
    if ($@)
    {
	die "json parse failed: $@\n$json\n";
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


1;
