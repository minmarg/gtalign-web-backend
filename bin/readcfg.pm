package readcfg;

## (C) 2020 Mindaugas Margelevicius, Institute of Biotechnology, Vilnius University
## Perl package for reading backend configuration

use strict;
use FindBin;
use lib "$FindBin::Bin";

## -------------------------------------------------------------------

BEGIN {
    my  $MAILADDRESSEE = '';
    my  $MAILSENDER = '';
    my  $MAILSERVER = '';

    sub GetMailAddressee { return $MAILADDRESSEE; }
    sub GetMailSender  { return $MAILSENDER; }
    sub GetMailServer  { return $MAILSERVER; }
}

## ===================================================================

sub new {
    my $that = shift;
    my $confile = shift;
    my $class = ref( $that ) || $that;
    my $self;

    $self->{LOCPDBDIR} = '';
    $self->{LOCSCOPeDIR} = '';
    $self->{LOCECODDIR} = '';
    $self->{WEBSCOPePDB} = '';
    $self->{WEBECODPDB} = '';
    $self->{WEBAFPDB} = '';

    $self->{JOB_NUM_CPUS} = 1;
    $self->{JOB_MAX_NQUERIES} = 1;
    $self->{JOB_MAX_SZINPUT} = 1048576; ##1MB

    $self->{PATHSTRDB_PDB} = '';
    $self->{PATHSTRDB_SCOP} = '';
    $self->{PATHSTRDB_ECOD} = '';
    $self->{PATHSTRDB_SwissProt} = '';
    $self->{PATHSTRDB_Proteomes} = '';
    $self->{PATHSTRDB_UniRef30} = '';

    $self->{INSTALL_GTALIGN} = '';

    $self->{MAILADDRESSEE} = GetMailAddressee();
    $self->{MAILSENDER}  = GetMailSender();
    $self->{MAILSERVER}  = GetMailServer();

    while( scalar( @_ )) {
        $self->{uc( $_[0] )} = $_[1];
        shift, shift;
    }

    bless( $self, $class );
    $self->Read( $confile );
    return $self;
}

## -------------------------------------------------------------------
## read/write member methods
##

sub Path3dDb_PDB { my $self = shift; if (@_) { $self->{LOCPDBDIR} = shift } return $self->{LOCPDBDIR}; }
sub Path3dDb_SCOPe { my $self = shift; if (@_) { $self->{LOCSCOPeDIR} = shift } return $self->{LOCSCOPeDIR}; }
sub Path3dDb_ECOD { my $self = shift; if (@_) { $self->{LOCECODDIR} = shift } return $self->{LOCECODDIR}; }
sub WebAddr3dDb_SCOPe { my $self = shift; if (@_) { $self->{WEBSCOPePDB} = shift } return $self->{WEBSCOPePDB}; }
sub WebAddr3dDb_ECOD { my $self = shift; if (@_) { $self->{WEBECODPDB} = shift } return $self->{WEBECODPDB}; }
sub WebAddr3dDb_AF { my $self = shift; if (@_) { $self->{WEBAFPDB} = shift } return $self->{WEBAFPDB}; }

sub JobMaxNoQueries { my $self = shift; if (@_) { $self->{JOB_MAX_NQUERIES} = shift } return $self->{JOB_MAX_NQUERIES}; }
sub JobMaxSizeInput { my $self = shift; if (@_) { $self->{JOB_MAX_SZINPUT} = shift } return $self->{JOB_MAX_SZINPUT}; }

sub PathStrDb_PDB { my $self = shift; if (@_) { $self->{PATHSTRDB_PDB} = shift } return $self->{PATHSTRDB_PDB}; }
sub PathStrDb_SCOP { my $self = shift; if (@_) { $self->{PATHSTRDB_SCOP} = shift } return $self->{PATHSTRDB_SCOP}; }
sub PathStrDb_ECOD { my $self = shift; if (@_) { $self->{PATHSTRDB_ECOD} = shift } return $self->{PATHSTRDB_ECOD}; }
sub PathStrDb_SwissProt { my $self = shift; if (@_) { $self->{PATHSTRDB_SwissProt} = shift } return $self->{PATHSTRDB_SwissProt}; }
sub PathStrDb_Proteomes { my $self = shift; if (@_) { $self->{PATHSTRDB_Proteomes} = shift } return $self->{PATHSTRDB_Proteomes}; }
sub PathStrDb_UniRef30 { my $self = shift; if (@_) { $self->{PATHSTRDB_UniRef30} = shift } return $self->{PATHSTRDB_UniRef30}; }

sub InstallDir_GTalign { my $self = shift; if (@_) { $self->{INSTALL_GTALIGN} = shift } return $self->{INSTALL_GTALIGN}; }

sub MailAddressee { my $self = shift; if (@_) { $self->{MAILADDRESSEE} = shift } return $self->{MAILADDRESSEE}; }
sub MailSender  { my $self = shift; if (@_) { $self->{MAILSENDER} = shift }  return $self->{MAILSENDER}; }
sub MailServer  { my $self = shift; if (@_) { $self->{MAILSERVER} = shift }  return $self->{MAILSERVER}; }


sub Exists     { my $self = shift; if( @_ ) { return exists $self->{$_[0]} }return 0; }
sub GetValue   { my $self = shift; if( @_ && $self->Exists( $_[0] )) { return $self->{$_[0]} }return 0; }
sub SetValue   { my $self = shift; if( @_ && defined( $_[0] )&& defined( $_[1] )) { $self->{$_[0]} = $_[1]; } }

sub ExistsKeyByValue  { my $self = shift; if( @_ ) { my $value = shift; foreach( keys %$self ){ return 1 if $self->{$_} eq $value; } }return 0; }
sub GetKeyByValue     { my $self = shift; if( @_ ) { my $value = shift; foreach( keys %$self ){ return $_ if $self->{$_} eq $value; } }return 0; }

## -------------------------------------------------------------------
## read configuration variables
##

sub Read
{
    my  $self = shift;
    my  $confile = shift;
    my  $class = ref( $self ) || die( "ERROR: Read: Should be called by object." );
    my ($key, $value );

    unless( open( F, $confile )) {
        printf( STDERR "ERROR: Read: Failed to open %s\n", $confile );
        return( 0 );
    }

    while( <F> ) {
        chomp;
        next if /^$/ || /^\s*#/;
        s/\s*(.*)$/$1/;
        next unless /^\s*([\w\.]+)\s*=\s*(.+)$/;
        $key = $1;
        $value = $2;

        ##unless( grep {/^$key$/} keys %$self ) {
        ##    printf( STDERR "WARNING: Read: Unrecognized variable %s; skipped.\n", $key );
        ##    next;
        ##}

        if( $value =~ /^\s*['"]([^'"]*)['"]/ ) { $value = $1; }
        elsif( $value =~ /^\s*['"]?([^\s#]*)/ ) { $value = $1; }
        else {
            printf( STDERR "WARNING: Read: Unable to read value of variable.\n" );
            next;
        }
        $self->{$key} = $value;
    }

    close( F );

    return 1;
}

## -------------------------------------------------------------------

1;

