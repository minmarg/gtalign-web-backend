package readopt;

## (C) 2009-2020 Mindaugas Margelevicius, Institute of Biotechnology, Vilnius University
## Read parameters along with their values

use strict;
use FindBin;
use lib "$FindBin::Bin";

## -------------------------------------------------------------------

BEGIN {
}

## ===================================================================

sub new {
    my $that = shift;
    my $parfile = shift;
    my $class = ref( $that ) || $that;
    my $self;

    $self->{DUMMY} = '';
    while( scalar( @_ )) {
        $self->{uc( $_[0] )} = $_[1];
        shift, shift;
    }

    bless( $self, $class );
    $self->Read( $parfile );
    return $self;
}

## -------------------------------------------------------------------
## read/write member methods
##

sub Exists     { my $self = shift; if( @_ ) { return exists $self->{$_[0]} }return 0; }
sub GetValue   { my $self = shift; if( @_ && $self->Exists( $_[0] )) { return $self->{$_[0]} }return 0; }
sub SetValue   { my $self = shift; if( @_ && defined( $_[0] )&& defined( $_[1] )) { $self->{$_[0]} = $_[1]; } }

## -------------------------------------------------------------------
## read parameters along with their values
##

sub Read
{
    my  $self = shift;
    my  $parfile = shift;
    my  $class = ref( $self ) ||
            die( "ERROR: ",__PACKAGE__,": 'Read' should be called by an object." );
    my ($key, $value );

    unless( open( F, $parfile )) {
        printf( STDERR "ERROR: ",__PACKAGE__,": Failed to open %s\n", $parfile );
        return( 0 );
    }
    while( <F> ) {
        chomp;
        next if /^$/ || /^\s*#/;
        s/\s*(.*)$/$1/;
        next unless /^\s*(\w+)\s*=\s*(.+)$/;
        $key = $1;
        $value = $2;

        if( $value =~ /^\s*['"]([^'"]*)['"]/ ) { $value = $1; }
        elsif( $value =~ /^\s*['"]?([^\s#]*)/ ) { $value = $1; }
        else {
            print(STDERR "WARNING: ",__PACKAGE__,": No value obtained for key $key\n");
            next;
        }
        $self->{$key} = $value;
    }
    close( F );
    return 1;
}

## -------------------------------------------------------------------

1;

