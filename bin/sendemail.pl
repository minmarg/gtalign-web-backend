#!/usr/bin/env perl
BEGIN {$^W=1}

## (C) 2020 Mindaugas Margelevicius, Institute of Biotechnology, Vilnius University
## Simple sending of a message by e-mail

use strict;
use FindBin;
use lib "$FindBin::Bin";
use readcfg;
use tsmtp;
use File::Basename;
use Getopt::Long;
use POSIX qw( strftime );

my  $MYPROGNAME = basename($0,('.pl'));
my  $CFGFILE = "$FindBin::Bin/../var/gtalign-ws-backend.conf";

my  $usage = <<EOIN;

Send an e-mail.
(C) 2020 Mindaugas Margelevicius, Institute of Biotechnology, Vilnius University

Usage:
$0 [OPTIONS]

Options:

--body <text_in_quotes> Message body.

--file <filename>       File to send as a message body 
                        if --body is not provided.

--sub <text>            Subject of the message.

--nochange              Do not change every \\n escaped in the 
                        message body with a newline.

--help                  Output this help text.

EOIN

#--to <email_address>    Addressee (mandatory).
#
#--from <email_address>  Sender (mandatory).
#
#--srv <server_name>     SMTP server name. If not provided,
#                        the message will not be sent.
#
#--help                  Output this help text.
#
#EOIN

my  ($Body,$File);
my  $Subject = '';
my  $Sender;
my  $nochange = 0;

my  $result = GetOptions(
               'body=s'   => \$Body,
               'file=s'   => \$File,
               'sub=s'    => \$Subject,
#               'to=s'     => \$Addressee,
#               'from=s'   => \$Sender,
#               'srv=s'    => \$Mailserver,
               'nochange' => sub {$nochange=1;},
               'help|h'   => sub {print $usage; exit(0);}
);

unless($result) {
    print STDERR "ERROR: $MYPROGNAME: Error in command-line arguments.\n";
    exit(1);
}
if(!$Body && $File && !-f $File) {
    print STDERR "ERROR: $MYPROGNAME: File not found: $File\n";
}
unless(-f $CFGFILE) {
    print STDERR "ERROR: $MYPROGNAME: Config file not found: $CFGFILE\n";
    exit(1);
}

## ===================================================================

my  $tsmtp = tsmtp->new();
my  $cfgvar = readcfg->new($CFGFILE);

unless($cfgvar->MailServer() && $cfgvar->MailAddressee()) {
    printf( STDERR "$MYPROGNAME: Mail unsent: Mail server and/or addressee is not provided.\n" );
    exit(0);
}

$Sender = $cfgvar->MailSender();
chomp($Sender = 'no-reply@'.`hostname`) unless($Sender);

$tsmtp->Server($cfgvar->MailServer());
$tsmtp->Sender($Sender);

printf(STDERR "$MYPROGNAME: Sending email to %s (from %s)...\n",
  $cfgvar->MailAddressee(),$tsmtp->Sender());

exit(1) unless SendMsg($tsmtp, $cfgvar->MailAddressee(), $Subject, \$Body, $File, $nochange);

printf(STDERR "$MYPROGNAME: ok.\n");

exit(0);

## ===================================================================

sub SendMsg
{
    my  $smtp = shift;       ## smtp object
    my  $sendto = shift;     ## recipient e-mail address
    my  $subject = shift;    ## subject of e-mail
    my  $rbody = shift;      ## reference to body text
    my  $file = shift;       ## name of file to send
    my  $nochange = shift;   ## do not change the body
    my  $nl = "\n";

    $$rbody =~ s/\\n/$nl/g if(!$nochange && $$rbody);

    if(!$$rbody && $file) {
        unless(open(FL, $file)) {
            printf(STDERR "ERROR: $MYPROGNAME: SendMsg: File open failed for $file: $!\n");
            return 0;
        }
        $$rbody .= $_ while(<FL>);
        close(FL);
    }

    return 0 unless $smtp->Send($sendto, $subject, $rbody);
    return 1;
}

## <<>>

