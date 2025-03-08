package tsmtp;

## A tiny smtp client package
##
## 2008 (C) Mindaugas Margelevicius
## Institute of Biotechnology
## Vilnius, Lithuania
##

use strict;
use Net::SMTP;

## -------------------------------------------------------------------

BEGIN {

    my  $MAILSERVER = 'mail.srv.lt';
    my  $SENDERADDR = 'moco@srv.lt';

    my  $hostname;

    chomp($hostname = `hostname`);

    my  $HOST = $hostname;
    $HOST.='.none.lt' unless $HOST =~ /\./;

    sub GetServer     { return $MAILSERVER; }
    sub GetSender     { return $SENDERADDR; }
    sub GetHostAddr   { return \$HOST; }
}

## ===================================================================

sub new {
    my $that  = shift;
    my $class = ref( $that ) || $that;
    my $self;

    $self->{SERVER}  = GetServer();   ## server address
    $self->{SENDER}  = GetSender();   ## sender address
    $self->{"_HOST"} = GetHostAddr(); ## sender host

    while( scalar( @_ )) {
        $self->{uc( $_[0] )} = $_[1];
        shift, shift;
    }

    return bless  $self, $class;
}

## -------------------------------------------------------------------
## read/write member methods
##

sub Server   { my $self = shift; if (@_) { $self->{SERVER} = shift } return $self->{SERVER}; }
sub Sender   { my $self = shift; if (@_) { $self->{SENDER} = shift } return $self->{SENDER}; }

## -------------------------------------------------------------------

sub get_day {
    my $num = shift;
    if ($num == 0)  { return 'Sun'; }
    elsif ($num == 1)  { return 'Mon'; }
    elsif ($num == 2)  { return 'Tue'; }
    elsif ($num == 3)  { return 'Wed'; }
    elsif ($num == 4)  { return 'Thu'; }
    elsif ($num == 5)  { return 'Fri'; }
    elsif ($num == 6)  { return 'Sat'; }
    return('');
}

sub get_month {
    my $num = shift;
    if ($num == 0)  { return 'Jan'; }
    elsif ($num == 1)  { return 'Feb'; }
    elsif ($num == 2)  { return 'Mar'; }
    elsif ($num == 3)  { return 'Apr'; }
    elsif ($num == 4)  { return 'May'; }
    elsif ($num == 5)  { return 'Jun'; }
    elsif ($num == 6)  { return 'Jul'; }
    elsif ($num == 7)  { return 'Aug'; }
    elsif ($num == 8)  { return 'Sep'; }
    elsif ($num == 9)  { return 'Oct'; }
    elsif ($num == 10) { return 'Nov'; }
    elsif ($num == 11) { return 'Dec'; }
    return('');
}

## -------------------------------------------------------------------
## determine full host name
## (moved to BEGIN)

#sub DetermineHost
#{
#    my  $refhost = shift;
#    my  $defaultname = shift;
#
#    my  $hostname;
#    my  $domainname;
#
#    $$refhost = $defaultname;
#
#    chomp( $hostname = `hostname` );
#    chomp( $domainname = `domainname` );
#
#    return if $hostname =~ /^\s*$/;
#
#    if( $domainname =~ /^\s*$/ ) {
#        $domainname = $defaultname; 
#        $domainname =~ s/^[^\.]+\.(.+)$/$1/;
#    }
#
#    $$refhost = "$hostname.$domainname";
#    $$refhost = 'localhost' unless $$refhost;
#}

## -------------------------------------------------------------------
## send message via SMTP
##

sub Send
{
    my  $self = shift;
    my  $class = ref( $self ) || die( "ERROR: Send: Should be called by object." );

    my  $sendto = shift;
    my  $subject = shift;
    my  $message = shift;

    my  $mailserver = $self->{SERVER};
    my  $senderaddr = $self->{SENDER};
    my  $host = ${$self->{"_HOST"}};

    my  $smtp;
    my  $recipient;
    my  @goodrecips;

    ##{{ RFC compliant time stamp
    my ($sec,$min,$hour,$mday,$mon,$year,$day) = gmtime();
    my  $timezone = '+0000';
    $year += 1900; $mon = get_month($mon); $day = get_day($day);
    my $date = sprintf("%s, %s %s %d %.2d:%.2d:%.2d %s",$day, $mday, $mon, $year, $hour, $min, $sec, $timezone);
    ##}}

    unless( $sendto ) {
        printf( STDERR "ERROR: Send: Address to send message to is empty.\n" );
        return 0;
    }

    if( ref( $message ) && ref( $message ) ne 'SCALAR' ) {
        printf( STDERR "ERROR: Send: Message refers to a non scalar object.\n" );
        return 0;
    }

    $smtp = Net::SMTP->new( $mailserver, Hello => $host );

    if( !ref( $smtp ) || $@ ) {
        printf( STDERR "ERROR: Send: Failed to connect to $mailserver: $@\n" );
        return 0;
    }


    unless( $smtp->mail( "<$senderaddr>" )) {
        printf( STDERR "ERROR: Send: Setting mail( $senderaddr ) failed.\n" );
        $smtp->quit();
        return 0;
    }

    $recipient = $sendto;
    $recipient =~ s/^[^<]*<([^>]+)>.*$/$1/;
    @goodrecips = $smtp->recipient( "<$recipient>", { SkipBad => 1 });

    if( $#goodrecips < 0 || $goodrecips[0] ne "<$recipient>" ) {
        printf( STDERR "ERROR: Send: Setting recipient $recipient failed.\n" );
        $smtp->quit();
        return 0;
    }

    unless( $smtp->data()) {
        printf( STDERR "ERROR: Send: Initilization of send to $sendto failed.\n" );
        $smtp->quit();
        return 0;
    }
    unless( $smtp->datasend( "From: $senderaddr\n" )) {
        printf( STDERR "ERROR: Send: Send of 'From: $senderaddr' failed.\n" );
        $smtp->quit();
        return 0;
    }
    unless( $smtp->datasend( "To: $sendto\n" )) {
        printf( STDERR "ERROR: Send: Send of 'To: $sendto' failed.\n" );
        $smtp->quit();
        return 0;
    }
    unless( $smtp->datasend( "Subject: $subject\n" )) {
        printf( STDERR "ERROR: Send: Send of 'Subject:' failed.\n" );
        $smtp->quit();
        return 0;
    }
    unless( $smtp->datasend( "Date: $date\n" )) {
        printf( STDERR "ERROR: Send: Send of 'Date:' failed.\n" );
        $smtp->quit();
        return 0;
    }
    unless( $smtp->datasend( "MIME-Version: 1.0\n" )) {
        printf( STDERR "ERROR: Send: Send of 'MIME-Version:' failed.\n" );
        $smtp->quit();
        return 0;
    }
    #unless( $smtp->datasend( "Content-Type: text/plain;\n  charset=\"us-ascii\"\n" )) {
    unless( $smtp->datasend( "Content-Type: text/plain; charset=utf-8\n" )) {
        printf( STDERR "ERROR: Send: Send of 'Content-Type:' failed.\n" );
        $smtp->quit();
        return 0;
    }
    unless( $smtp->datasend( "\n" )) {
        printf( STDERR "ERROR: Send: Send of newline failed.\n" );
        $smtp->quit();
        return 0;
    }
    unless( $smtp->datasend( ref($message)? ($$message? "$$message\n":""): ($message? "$message\n":"") )) {
        printf( STDERR "ERROR: Send: Send of message failed.\n" );
        $smtp->quit();
        return 0;
    }
    unless( $smtp->dataend()) {
        printf( STDERR "ERROR: Send: Finalization of send failed.\n" );
        $smtp->quit();
        return 0;
    }

    unless( $smtp->quit()) {
        printf( STDERR "ERROR: Send: quit failed: $@\n" );
        return 0;
    }

    return 1;
}

## -------------------------------------------------------------------

1;

