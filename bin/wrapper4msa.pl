#!/usr/bin/env perl
BEGIN {$^W=1}

## (C) 2024 Mindaugas Margelevicius, Institute of Biotechnology, Vilnius University
## Wrapper for constructing the multiple sequence alignment from a given set of
## pairwise alignments

use strict;
use FindBin;
use lib "$FindBin::Bin";
use readcfg;
use readopt;
use File::Spec;
use File::Basename;
use Getopt::Long;
use POSIX qw(strftime);

##constants:
my  $MYPROGNAME = basename($0);
my  $PWFA2MSAprog = File::Spec->catfile($FindBin::Bin,"pwfa2msa.pl");
my  $SENDMAILprog = File::Spec->catfile($FindBin::Bin,"sendemail.pl");
my  $CFGFILE = File::Spec->catfile($FindBin::Bin,File::Spec->updir(),"var","gtalign-ws-backend.conf");
my  $devnull = File::Spec->devnull();
my  $TARPROG = `which tar 2>$devnull`; chomp($TARPROG);
my  $GZPROG = `which gzip 2>$devnull`; chomp($GZPROG);
my  $MAXNCPUs = 6;##maximum number of CPU cores assigned for a job
##maximum number of sequences per search program permitted to be included in the results
my  $MAXNSEQS_perprog = 20000;

my ($AFAEXT) = ('afa');

my  $usage = <<EOIN;

Wrapper for constructing the multiple sequence alignment from a given set of
pairwise alignments.
(C) 2024 Mindaugas Margelevicius, Institute of Biotechnology, Vilnius University

Usage:
$0 <Options>

Options:

--in <input>     Input file of pairwise alignments. Format is provided below.

--opt <options>  Filename of options for this job.

--status <file>  Name of the output status file of computation progress
                 messages.

--log <logfile>  Name of the output log file generated during the execution of
                 this script.

--results <file> Name of the output compressed file containing the results
                 file(s). (Extension .gz will be added.)

--err <errfile>  Name of the file of high-level error messages to write to on
                 error.

--help           This help text.


Format of the input file of pairwise alignments:

><query_description> (ALN:<query_start>-<query_end>)
<aligned_query_sequence>
><db_sequence_description> (ALN:<dbseq_start>-<dbseq_end>)
<aligned_database_sequence>
//
...

EOIN


my  $INPFILENAME = '';
my  $OPTFILENAME = '';
my  $STAFILENAME = '';
my  $LOGFILENAME = '';
my  $RESULTSFILENAME = '';
my  $ERRFILENAME = '';
my  $Fail = 0;

my  $result = GetOptions(
               'in=s'       => \$INPFILENAME,
               'opt=s'      => \$OPTFILENAME,
               'status=s'   => \$STAFILENAME,
               'log=s'      => \$LOGFILENAME,
               'results=s'  => \$RESULTSFILENAME,
               'err=s'      => \$ERRFILENAME,
               'help|h'     => sub {print $usage; exit(0);}
);

print(STDERR "\n\n".GetDatetime()."\n\n");

my ($inpbasename,$inpdirname,$suffix) = fileparse($INPFILENAME, qr/\.[^.]*/);
$inpdirname = File::Spec->rel2abs($inpdirname);

unless($STAFILENAME) {
    $STAFILENAME = File::Spec->catfile(${inpdirname},"${inpbasename}.status");
    print(STDERR "WARNING: $MYPROGNAME: Filename for computation progress mesages not given; ".
          "it has been set to: '$STAFILENAME'\n");
}
unless($LOGFILENAME) {
    $LOGFILENAME = File::Spec->catfile(${inpdirname},"${inpbasename}__msa_out.log");
    print(STDERR "WARNING: $MYPROGNAME: Filename for log mesages not given; ".
          "it has been set to: '$LOGFILENAME'\n");
}
unless($RESULTSFILENAME) {
    $RESULTSFILENAME = File::Spec->catfile(${inpdirname},"${inpbasename}__msa_out.tar");
    print(STDERR "WARNING: $MYPROGNAME: Filename for results not given; ".
          "it has been set to: '$RESULTSFILENAME'\n");
}

my ($resbasename,$resdirname,$ressuffix) = fileparse($RESULTSFILENAME, qr/\.[^.]*/);
my $afafilename = File::Spec->catfile(${resdirname},"${resbasename}.${AFAEXT}");

unless($ERRFILENAME) {
    $ERRFILENAME = File::Spec->catfile(${inpdirname},"${inpbasename}.err");
    print(STDERR "WARNING: $MYPROGNAME: Filename for high-level error mesages not given; ".
          "it has been set to: '$ERRFILENAME'\n");
}
## truncate error file before using Error!
if(open(F,'>',$ERRFILENAME)){close(F);}
unless($result) {
    Error("ERROR: $MYPROGNAME: Error in command-line arguments.\n",
            "Command-line arguments error.\n");## h-l error message
    MyExit(1);
}
unless($INPFILENAME && -f $INPFILENAME) {
    Error("ERROR: $MYPROGNAME: Input file not found: '$INPFILENAME'\n",
            "Input file not found.\n");## h-l error message
    MyExit(1);
}
##unless($OPTFILENAME && -f $OPTFILENAME) {
##    Error("ERROR: $MYPROGNAME: Input job options file not found: '$OPTFILENAME'\n",
##            "Job options file not found.\n");## h-l error message
##    MyExit(1);
##}
unless(-f $CFGFILE) {
    Error("ERROR: $MYPROGNAME: Config file not found: '$CFGFILE'\n",
            "Configuration file not found.\n");## h-l error message
    MyExit(1);
}
unless(-f $PWFA2MSAprog) {
    Error("ERROR: $MYPROGNAME: Program file not found: '$PWFA2MSAprog'\n",
            "Some of the required program files not found.\n");## h-l error message
    MyExit(1);
}
unless($TARPROG && $GZPROG) {
    Error("ERROR: $MYPROGNAME: System programs 'tar' and/or 'gzip' not found.\n",
            "Some of the system programs not found.\n");## h-l error message
    MyExit(1);
}


## =============================================================================

my  $cfgvar = readcfg->new($CFGFILE);
my  %optionvalues;


##check particular programs
$optionvalues{prog_comer_comer2msa} = '';

$optionvalues{prog_comer_comer2msa} = File::Spec->catfile($FindBin::Bin,'comer2msa.pl');
unless(-f $optionvalues{prog_comer_comer2msa}) {
    Error("ERROR: $MYPROGNAME: comer2msa program not found: '".$optionvalues{prog_comer_comer2msa}."'\n",
        "Incomplete software installed on the system.\n");## h-l error message
    MyExit(1);
}

## =============================================================================
## MAIN
##

##my  $options = readopt->new($OPTFILENAME);
my  %optionkeys = (
    job_num_cpus => 'JOB_NUM_CPUS'
);
##my  $nsyscpus = GetNCPUs();
##my  $ncpus = $MAXNCPUs;

my  $errmsg = "ERROR: $MYPROGNAME: Invalid job options.\n";
my  $hlerrmsg = "Invalid job options.\n";

$errmsg = '';
$hlerrmsg = '';

ProgressMsg("Buiding MSA from pairwise alignments...\n");

my  $msasize = 0;##size of the results file
my  $nmsalines = 0;##number of lines in the results file

unless(RunBuilder(
    $INPFILENAME, $inpdirname, $inpbasename, $afafilename, $LOGFILENAME,
    \%optionkeys, \%optionvalues, \$msasize, \$nmsalines, \$errmsg, \$hlerrmsg)) { 
    Error($errmsg, $hlerrmsg);
    MyExit(1);
}

ProgressMsg(sprintf("  (MSA size: %.3fMB  #lines: %.3fK)\n",$msasize/1048576,$nmsalines/1000))
    ;##if($msasize && $nmsalines);

ProgressMsg("Finalizing results...\n");

my  @resfilelist;

unless(MakeResultsList(
    \@resfilelist, $inpdirname, $afafilename, $LOGFILENAME,
    \$errmsg, \$hlerrmsg)) {
    Error($errmsg, $hlerrmsg);
    MyExit(1);
}

unless(CompressResults(
    $TARPROG, $GZPROG, $inpdirname, $RESULTSFILENAME,
    \@resfilelist,
    $LOGFILENAME, $STAFILENAME, $ERRFILENAME,
    \$errmsg, \$hlerrmsg)) {
    Error($errmsg, $hlerrmsg);
    MyExit(1);
}

ProgressMsg("Finished.\n");

MyExit(0);

## =============================================================================
## -----------------------------------------------------------------------------
## GetTime: get time string
##
sub GetTime
{
    return strftime("%H:%M:%S ",localtime());
}

sub GetDatetime
{
    return strftime("%a %b %d %H:%M:%S %Z %Y",localtime());
}

## -----------------------------------------------------------------------------
## GetNCPUs: get the number of CPU cores in the system; return <0 on error;
##
sub GetNCPUs
{
    unless(open(H, "/proc/cpuinfo")) {
        print(STDERR "WARNING: Unable to determine #cpu cores.\n");
        return -1;##failed to open cpuinfo: $!
    }
    my $ncpus = scalar(grep(/^processor/,<H>)); 
    close(H);
    return $ncpus;
}

## -----------------------------------------------------------------------------
## ProgressMsg: Print progress message to file
##
sub ProgressMsg
{
    my  $msg = shift;
    return 0 unless(open(F,'>>',$STAFILENAME));
    print(F $msg);
    close(F);
    return 1;
}

## -----------------------------------------------------------------------------
## GetWarnings: extract warnings (if any) from file
##
sub GetWarnings
{
    my  $logfile = shift;##builder log file
    my  $rlist = shift;##ref to the output list of warnings
    unless(open(F, $logfile)) {
        print(STDERR "WARNING: Unable to open builder's log file: '$logfile'.\n");
        return 1;
    }
    @$rlist = grep(/WARNING/i,<F>);
    close(F);
    return 0;
}

## -----------------------------------------------------------------------------
## GetFileStats: get some stats for the given file
##
sub GetFileStats
{
    my  $resfile = shift;##file to get stats about
    my  $rsize = shift;##ref to the size
    my  $rnlines = shift;##ref to the number of lines
    $$rsize = $$rnlines = 0;
    unless(open(F, $resfile)) {
        print(STDERR "WARNING: Unable to open results file: '$resfile'.\n");
        return 1;
    }
    $$rsize = -s F;
    $$rnlines += tr/\n/\n/ while sysread(F, $_, 65536);
    close(F);
    return 0;
}

## -----------------------------------------------------------------------------
## AddFileToArchive: add a file to a given archive
##
sub AddFileToArchive
{
    my  $tarprog = shift;##full pathname to the tar program
    my  $archive = shift;##fullname of the resulting archive file
    my  $filepathname = shift;##full pathname to the file
    my  $dirname = shift;##name of directory where $filepathname is; may be empty
    my  $rerrmsg = shift;##ref to the error message string to be put in logs
    my  $rhlerrmsg = shift;##ref to the h-l error message string
    my  $create = shift;##flag of whether the archive is to be created

    my  $mysubname = (caller(0))[3];
    my  $preamb = "[ ${mysubname} ] ";
    my  $command = '';
    my  $opt = $create?'c':'r';
    my  $ret = 1;

    if($dirname) {
        $command = "${tarprog} -${opt}f \"${archive}\" -C \"$dirname\" \"${filepathname}\"";
    }
    else {
        my $filedirname = dirname($filepathname);
        my $filename = basename($filepathname);
        $command = "${tarprog} -${opt}f \"${archive}\" -C \"${filedirname}\" \"${filename}\"";
    }

    print(STDERR GetTime()."${preamb} ${command}\n");

    unless(ExecCommand($command)) {
        $$rerrmsg = "ERROR: $MYPROGNAME: $mysubname: Failed to add file to archive: '${filepathname}'\n";
        $$rhlerrmsg = "Failed to archive some of the results files.\n";
        return 0;
    }

    return $ret;
}

## -----------------------------------------------------------------------------
## RunBuilder: build the MSA for pairwise alignments
##
sub RunBuilder
{
    my  $inpfilename = shift;##full name of input
    my  $dirname = shift;##name of the job directory
    my  $basename = shift;##basename of the input transfered to the backend
    my  $ouputFILENAME = shift;##name of output file
    my  $logfilename = shift;##name of log file
    my  $roptionkeys = shift;##ref to the keys of options
    my  $roptionvalues = shift;## ref to the values of options
    my  $rsize = shift;##ref to the size of the results file
    my  $rnlines = shift;##ref to the number of lines in the results file
    my  $rerrmsg = shift;##ref to the error message string to be put in logs
    my  $rhlerrmsg = shift;##ref to the h-l error message string

    my  $mysubname = (caller(0))[3];
    my  $preamb = "[ ${mysubname} ] ";
    my  $command = '';
    my  @warnings;
    my  $ret = 1;

    ##run builder
    $command = "${PWFA2MSAprog} ".
            "-i \"${inpfilename}\" -o \"${ouputFILENAME}\" -a -f 0 -e 1.e19 ".
            "-N ${MAXNSEQS_perprog}  >\"${logfilename}\" 2>&1";

    print(STDERR GetTime()."${preamb} ${command}\n\n");

    unless(ExecCommand($command)) {
        $$rerrmsg = "ERROR: $MYPROGNAME: $mysubname: Unable to build the MSA.\n";
        $$rhlerrmsg = "Unable to build the MSA.\n";
        return 0;
    }

    ##extract warnings if any; no check of return code
    GetWarnings($logfilename, \@warnings);

    ##get some stats on the results file; no check of return code
    GetFileStats($ouputFILENAME, $rsize, $rnlines);

    if(0 <= $#warnings) {
        my $lst = join("", @warnings);
        my $text = "\nWarnings from building the MSA:\n\n$lst\n\n";
        Warning($text, $text, 1);##no e-mail
    }

    return $ret;
}

## -----------------------------------------------------------------------------
## MakeResultsList: make and write to a file the list of the results files
##
sub MakeResultsList
{
    my  $rfilelist = shift;##ref to the output list of files
    my  $dirname = shift;##name of the job directory
    my  $ouputFILENAME = shift;##name of the results file
    my  $logfilename = shift;##name of log file
    my  $rerrmsg = shift;##ref to the error message string to be put in logs
    my  $rhlerrmsg = shift;##ref to the h-l error message string

    my  $mysubname = (caller(0))[3];
    my  $preamb = "[ ${mysubname} ] ";
    my  $ret = 1;

    $ouputFILENAME =~ s/^${dirname}\/*(.+)$/$1/;
    push @$rfilelist, $ouputFILENAME;

    return $ret;
}

## -----------------------------------------------------------------------------
## CompressResults: compress all required results files
##
sub CompressResults
{
    my  $tarprog = shift;##full pathname to the tar program
    my  $gzprog = shift;##full pathname to the gzip program
    my  $dirname = shift;##name of the job directory
    my  $outarchivefile = shift;##name of the output archive file to be created
    my  $rfilelist = shift;##ref to the list of ($dirname-prefix-removed) files to include in the archive
    my  $logfilename = shift;##fullname of the log file
    my  $statusfile = shift;##fullname of the status messages file
    my  $errorfile = shift;##fullname of the high-level error messages file
    my  $rerrmsg = shift;##ref to the error message string to be put in logs
    my  $rhlerrmsg = shift;##ref to the h-l error message string

    my  $mysubname = (caller(0))[3];
    my  $preamb = "[ ${mysubname} ] ";
    my  $command = '';
    my  $ret = 1;

    for(my $i = 0; $i <= $#{$rfilelist}; $i++) {
        return 0 unless(AddFileToArchive($tarprog, 
            $outarchivefile, $$rfilelist[$i], $dirname, $rerrmsg, $rhlerrmsg, $i==0)); #dirname given; $i==0:create
    }

    unless(AddFileToArchive($tarprog, 
        $outarchivefile, $statusfile, '', $rerrmsg, $rhlerrmsg)){ #no dirname
        ##record an error, don't stop and  don't send an e-mail
        Error($$rerrmsg, $$rhlerrmsg, 1);
    }

    unless(AddFileToArchive($tarprog, 
        $outarchivefile, $errorfile, '', $rerrmsg, $rhlerrmsg)){ #no dirname
        ##record an error, don't stop and  don't send an e-mail
        Error($$rerrmsg, $$rhlerrmsg, 1);
    }

    $command = "${gzprog} -f \"${outarchivefile}\"";

    print(STDERR GetTime()."${preamb} ${command}\n\n");

    unless(ExecCommand($command)) {
        $$rerrmsg = "ERROR: $MYPROGNAME: $mysubname: Failed to gzip the archive: '${outarchivefile}'\n";
        $$rhlerrmsg = "Failed to compress the archive of results files.\n";
        return 0;
    }

    return $ret;
}

## General =====================================================================
##
##







## =============================================================================
## MyExit: print a code to file and exit
##
sub MyExit
{
    my $ecode = shift;##exit code
    if(open(F,">>",$ERRFILENAME)){print(F "\n${ecode}\n");close(F);}
    exit($ecode);
}

## -------------------------------------------------------------------
## Warning: output a warning message and optionally, send an email
##
sub Warning
{
    my ($msg, $hlerrmsg, $nosend) = @_;
    return Issue($msg, $hlerrmsg, $nosend,
        '',
        "Warning from gtalign-ws ($MYPROGNAME)");
}

## Error: output an error message and optionally, send an email
##
sub Error
{
    my ($msg, $hlerrmsg, $nosend) = @_;
    return Issue($msg, $hlerrmsg, $nosend,
        "ERROR: Issue on the server's backend side: ",
        "Error message from gtalign-ws ($MYPROGNAME)");
}

## Issue: output an issue message and send an email
##
sub Issue
{
    my ($msg, $hlerrmsg, $nosend,  $leadtxt1, $sbjct) = @_;
    print(STDERR GetTime().$msg);
    if($ERRFILENAME && $hlerrmsg) {
        if(open(EF,">>",$ERRFILENAME)) {
            print(EF $leadtxt1.$hlerrmsg);
            close(EF);
        } else {
            print(STDERR "ERROR: $MYPROGNAME: Write to error file skipped ".
                "due to the fail to open the file: '$ERRFILENAME'\n");
        }
    }
    return 1 if $nosend;
    unless(-f $SENDMAILprog) {
        print(STDERR "ERROR: $MYPROGNAME: Send-mail program not found: $SENDMAILprog\n");
        return 0;
    }
    my  $command = "$SENDMAILprog --sub \"$sbjct\" --body \"$msg\" 2>&1";
    print(STDERR "$command\n");
    return ExecCommand($command);
}

## -------------------------------------------------------------------
## try read lock files in a directory and return 1 if locks exist,
## 0 otherwise, -1 on error
##

sub ReadLocks
{
    return ReadFiles( shift, shift );
}

sub ReadFiles
{
    my  $dirname = shift;
    my  $pattern = shift;  ## pattern of files to look for as locks
    my  $reffiles = shift; ## reference to vector of files
    my  @files;
    my  $locref = defined( $reffiles )? $reffiles: \@files;

    unless( opendir( DIR, $dirname )) {
        printf( STDERR "ERROR: Cannot open directory $dirname.\n" );
        return -1;
    }
    @{$locref} = grep { /$pattern$/ && -f File::Spec->catfile($dirname,$_) } readdir( DIR );
    closedir( DIR );
    return 1 if 0 <= $#{$locref};
    return 0;
}

## -------------------------------------------------------------------
## ExecCommand: execute system command
##
sub CheckStatus
{
    ExecCommand();
}

sub ExecCommand
{
    my  $cmdline = shift;
    my  $returnstatus = shift;##if defined, return exit status of the command

    system($cmdline) if $cmdline;

    if($? == -1) {
        print(STDERR GetTime()."ERROR: Failed to execute command: $!\n");
        return($returnstatus? 1: 0);
    }
    if($? & 127) {
        printf(STDERR GetTime()."ERROR: Command terminated with signal %d, %s coredump.\n",
            ($? & 127), ($? & 128)? 'with': 'without');
        return($returnstatus? ($? & 127): 0);
    }
    else {
        if(($? >> 8) != 0) {
            printf( STDERR GetTime()."ERROR: Command failed and exited with status %d\n",
                    $? >> 8);
            return($returnstatus? ($? >> 8): 0);
        }
    }
    return($returnstatus? 0: 1);
}

## <<>>

