package PerlIO::via::Blockwise;

require 5.008;
use strict;
use warnings;
use IO::AIO;

our $VERSION = 0.01;

sub PUSHED {
    my ($class, $mode) = @_;
    return -1 if $mode ne 'r';
    # The following variables are updated / accessed via the closures below
    my $buffer = '';		# internal buffer for this layer
    my %inside = ();
    my $last_block = -1;
    bless {
	buffer => sub : lvalue { $buffer },
	last_block => sub : lvalue { $last_block },
    }, $class;
}

sub FILL {
    my ($self, $fh) = @_;

    my $blksz = 1048576;
    my $lb = $self->{last_block}->();
    if ($lb >= 0)
    {
#	print STDERR "Free blk $lb\n";
	my $r = IO::AIO::fadvise($fh, $lb * $blksz, $blksz, IO::AIO::FADV_DONTNEED);
	if ($r != 0)
	{
	    print STDERR "fadv err r='$r' '$!'\n";
	}
    }

    $self->{last_block}->()++;

    my($block, $len);
    $len = sysread($fh, $block, $blksz);
#    print STDERR "read $len\n";
    return undef if $len == 0;
#    $self->{buffer}->() .= $block;
#    return $self->{buffer}->();
    return $block;
}

1;
