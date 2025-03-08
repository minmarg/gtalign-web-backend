#!/usr/bin/env perl
BEGIN {$^W=1}

## (C) 2024 Mindaugas Margelevicius, Institute of Biotechnology, Vilnius University
## Helper script for superimposing a pair of 3D structures aligned by GTalign

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
my  $SUPERPOSEprog = File::Spec->catfile($FindBin::Bin,"superpose1.py");
my  $SENDMAILprog = File::Spec->catfile($FindBin::Bin,"sendemail.pl");
my  $CFGFILE = File::Spec->catfile($FindBin::Bin,File::Spec->updir(),"var","gtalign-ws-backend.conf");
my  $devnull = File::Spec->devnull();
my  $TARPROG = `which tar 2>$devnull`; chomp($TARPROG);
my  $GZPROG = `which gzip 2>$devnull`; chomp($GZPROG);
my  $MAXNCPUs = 2;##maximum number of CPU cores assigned to a job

my ($ENTEXT, $PDBEXT) = ('ent','pdb');

my  $usage = <<EOIN;

Helper script for superimposing a pair of 3D structures aligned by GTalign.
(C) 2024 Mindaugas Margelevicius, Institute of Biotechnology, Vilnius University

Usage:
$0 <Options>

Options:

--in <input>     Input file of a GTalign alignment in JSON format. This file
                 should contain the "query" and a "hit_record" sections from
                 a GTalign .json output file.

--status <file>  Name of the output status file of computation progress
                 messages.

--log <logfile>  Name of the output log file generated during the process.

--reslst <file>  Name of the output file listing results files.

--results <file> Name of the compressed output file containing resulting
                 superimposition(s). (Extension .gz will be added.)

--err <errfile>  Name of the file of high-level error messages to write to on
                 error.

--help           This help text.

EOIN


my  $INPFILENAME = '';
my  $STAFILENAME = '';
my  $LOGFILENAME = '';
my  $RESLSTFILENAME = '';
my  $RESULTSFILENAME = '';
my  $ERRFILENAME = '';
my  $Fail = 0;

my  $result = GetOptions(
               'in=s'       => \$INPFILENAME,
               'status=s'   => \$STAFILENAME,
               'log=s'      => \$LOGFILENAME,
               'reslst=s'   => \$RESLSTFILENAME,
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
    $LOGFILENAME = File::Spec->catfile(${inpdirname},"${inpbasename}__3d_out.log");
    print(STDERR "WARNING: $MYPROGNAME: Filename for log mesages not given; ".
          "it has been set to: '$LOGFILENAME'\n");
}
unless($RESLSTFILENAME) {
    $RESLSTFILENAME = File::Spec->catfile(${inpdirname},"${inpbasename}__3d_out.lst");
    print(STDERR "WARNING: $MYPROGNAME: Filename for a results list not given; ".
          "it has been set to: '$RESLSTFILENAME'\n");
}
unless($RESULTSFILENAME) {
    $RESULTSFILENAME = File::Spec->catfile(${inpdirname},"${inpbasename}__3d_out.tar");
    print(STDERR "WARNING: $MYPROGNAME: Filename for results not given; ".
          "it has been set to: '$RESULTSFILENAME'\n");
}
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
unless(-f $CFGFILE) {
    Error("ERROR: $MYPROGNAME: Config file not found: '$CFGFILE'\n",
            "Configuration file not found.\n");## h-l error message
    MyExit(1);
}
unless($TARPROG && $GZPROG) {
    Error("ERROR: $MYPROGNAME: System programs 'tar' and/or 'gzip' not found.\n",
            "Some of the system programs not found.\n");## h-l error message
    MyExit(1);
}
unless(-f $SUPERPOSEprog) {
    Error("ERROR: $MYPROGNAME: Package program not found: '$SUPERPOSEprog'\n",
            "Some of the package programs not found.\n");## h-l error message
    MyExit(1);
}


## =============================================================================

my  $cfgvar = readcfg->new($CFGFILE);
my  %optionvalues;
my ($resbasename,$resdirname,$ressuffix) = fileparse($RESULTSFILENAME, qr/\.[^.]*/);
    
## =============================================================================
## MAIN
##

my  %optionkeys = (
    job_num_cpus => 'JOB_NUM_CPUS'
);
my  $nsyscpus = GetNCPUs();
my  $ncpus = $MAXNCPUs;

my  $errmsg = "ERROR: $MYPROGNAME: Invalid job options.\n";
my  $hlerrmsg = "Invalid job options.\n";

##if($Config{useithreads}) {
    if($cfgvar->Exists($optionkeys{job_num_cpus})) {
        $ncpus = $cfgvar->GetValue($optionkeys{job_num_cpus});
        $ncpus = $nsyscpus if($nsyscpus > 0 && $nsyscpus < $ncpus);
    }
    else {
        print(STDERR "WARNING: Job option $optionkeys{job_num_cpus} not specified: ".
            "Using default #cpus=${ncpus}\n");
    }
##}
##else {
##    $ncpus = 1;
##    print(STDERR "WARNING: Perl compiled WITHOUT thread support: ".
##        "Using only one cpu: #cpus=${ncpus}\n");
##}

my  $nqueries = 0;## :shared = 0;##number of individual queries/inputs
my  %inputs;## :shared;##all inputs divided

$errmsg = "ERROR: $MYPROGNAME: Preprocessing input failed.\n";
$hlerrmsg = "Preprocessing input failed.\n";

ProgressMsg("Preparing for superposition...\n");

unless(PreprocessInputFile(
        $INPFILENAME, $resbasename, $resdirname, 
        \$nqueries, \%inputs, \$errmsg, \$hlerrmsg)) {
    Error($errmsg, $hlerrmsg);## err msg and h-l error message
    MyExit(1);
}

ProgressMsg("Running superposition(s)...\n");

unless(RunSuperposition(
        $LOGFILENAME,
        \$nqueries, \%inputs, \$errmsg, \$hlerrmsg)) {
    Error($errmsg, $hlerrmsg);## err msg and h-l error message
    MyExit(1);
}

$errmsg = '';
$hlerrmsg = '';

ProgressMsg("Finalizing results...\n");

my  @resfilelist;

unless(MakeResultsList(
        \@resfilelist, $inpdirname, $RESLSTFILENAME,
        $nqueries, \%inputs, \$errmsg, \$hlerrmsg)) {
    Error($errmsg, $hlerrmsg);
    MyExit(1);
}

unless(CompressResults(
        $TARPROG, $GZPROG, $inpdirname, $RESULTSFILENAME, 
        \@resfilelist, 
        $LOGFILENAME, $RESLSTFILENAME, $STAFILENAME, $ERRFILENAME,
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
    my  $logfile = shift;##log file
    my  $rlist = shift;##ref to the output list of warnings
    unless(open(F, $logfile)) {
        print(STDERR "WARNING: Unable to open log file: '$logfile'.\n");
        return 0;
    }
    @$rlist = grep(/WARNING/i,<F>);
    close(F);
    return 1;
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
## PreprocessInputFile: parse the input file and change the name of each query 
## inline so that duplicate names are avoided (Modeller will crash if the 
## sequence name and a template name are the same)
##
sub PreprocessInputFile
{
    my  $inpfilename = shift;##the input to be parsed for multiple individual inputs
    my  $resbasename = shift;##basename of the results file
    my  $resdirname = shift;##dirname of the results file
    my  $rnqueries = shift;##ref to the number of queries
    my  $rinputs = shift;##ref to the hash of individual inputn
    my  $rerrmsg = shift;##ref to the error message string to be put in logs
    my  $rhlerrmsg = shift;##ref to the h-l error message string

    my  $mysubname = (caller(0))[3];
    my ($qf, $hf, $hsn, $rfnd, $rot, $tran) = (0,0,-1,0,'','');
    my ($qname, $qchn, $qmod) = ('','','');
    my ($rname, $rchn, $rmod) = ('','','');
    my  $outfile = '';
    my  $ret = 1;

    unless(open(F, $inpfilename)) {
        $$rerrmsg = "ERROR: $MYPROGNAME: $mysubname: Failed to open input file: '$inpfilename'\n";
        $$rhlerrmsg = "Input file not found.\n";
        return 0;
    }
    while(<F>) {
        chomp;

        do {$qf = 1; next} if(/\s*"query":\s*{/);
        do {$qf = 0; next} if($qf && /\s*},?$/);
        if($qf && /\s*"description":\s*"(.+)"$/) {
            $qname = $1;
            $qchn = $1 if($qname =~ /Chn:(\S+)/);
            $qmod = $1 if($qname =~ /\(M:(\S+)\)/);
            $qname =~ s/^(.*)\s+Chn:.*$/$1/;
            $qname =~ s/^(.*)\s+\(M:\S+\).*$/$1/;
            $qname =~ s/\\//g;
            next
        }

        do {$hf = 1; next} if(/\s*{"hit_record":\s*{/);
        do {$rfnd = $1; next} if($hf && /\s*"tfm_referenced":\s*(\d),?$/);
        do {$rot = $1; $rot =~ s/\s//g; next} if($hf && /\s*"rotation_matrix_rowmajor":\s*\[(.+)\],?$/);
        do {$tran = $1; $tran =~ s/\s//g; next} if($hf && /\s*"translation_vector":\s*\[(.+)\],?$/);
        if($hf && /\s*"reference_description":\s*"(.+)",?$/) {
            $rname = $1;
            $rchn = $1 if($rname =~ /Chn:(\S+)/);
            $rmod = $1 if($rname =~ /\(M:(\S+)\)/);
            $rname =~ s/^(.*)\s+Chn:.*$/$1/;
            $rname =~ s/^(.*)\s+\(M:\S+\).*$/$1/;
            $rname =~ s/\\//g;
            next
        }

        if($hf && /\s*}},?$/) {
            $hf = 0;
            $hsn++;

            $$rinputs{"${hsn}_rcode"} = 1;

            unless($qname) {
                $$rinputs{"${hsn}_error"} = "ERROR: $MYPROGNAME: $mysubname: Query undefined in Record No. $hsn. Skipped.\n";
                $$rinputs{"${hsn}_errhl"} = "Query undefined in Record No. $hsn. Skipped.\n";
            }
            elsif(!$rname) {
                $$rinputs{"${hsn}_error"} = "ERROR: $MYPROGNAME: $mysubname: Reference undefined in Record No. $hsn. Skipped.\n";
                $$rinputs{"${hsn}_errhl"} = "Reference undefined in Record No. $hsn. Skipped.\n";
            }
            elsif(!($rot && $tran)) {
                $$rinputs{"${hsn}_error"} = "ERROR: $MYPROGNAME: $mysubname: Rotation and/or translation undefined in Record No. $hsn. Skipped.\n";
                $$rinputs{"${hsn}_errhl"} = "Rotation and/or translation undefined in Record No. $hsn. Skipped.\n";
            }
            else {
                $outfile = File::Spec->catfile($resdirname,$resbasename) . "__${hsn}.$PDBEXT";

                $$rinputs{"${hsn}_cmdline"} = "$SUPERPOSEprog --i1=\"$qname\" --i2=\"$rname\"";
                $$rinputs{"${hsn}_cmdline"} .= " -o \"$outfile\"";
                $$rinputs{"${hsn}_cmdline"} .= " --c1=$qchn" if($qchn);
                $$rinputs{"${hsn}_cmdline"} .= " --m1=$qmod" if($qmod);
                $$rinputs{"${hsn}_cmdline"} .= " --c2=$rchn" if($rchn);
                $$rinputs{"${hsn}_cmdline"} .= " --m2=$rmod" if($rmod);
                $$rinputs{"${hsn}_cmdline"} .= " --rot=$rot --tran=$tran";
                $$rinputs{"${hsn}_cmdline"} .= " -2" if($rfnd);

                $$rinputs{"${hsn}_rcode"} = 0;
                $$rinputs{"${hsn}_error"} = '';##error message
                $$rinputs{"${hsn}_errhl"} = '';##high-level error message
                $$rinputs{"${hsn}_outfl"} = $outfile;
            }
            $rfnd = 0; $rot = $tran = '';
            $rname = $rchn = $rmod = '';
            next
        }
    }
    close(F);

    $$rnqueries = ($hsn + 1);

    return $ret;
}

## -----------------------------------------------------------------------------
## RunSuperposition: Run superpositions
##
sub RunSuperposition
{
    my  $logfilename = shift;##filename of log messages
    my  $rnqueries = shift;##ref to the number of queries
    my  $rinputs = shift;##ref to the hash of individual inputs
    my  $rerrmsg = shift;##ref to the error message string to be put in logs
    my  $rhlerrmsg = shift;##ref to the h-l error message string

    my  $mysubname = (caller(0))[3];
    my  $preamb = "[ ${mysubname} ] ";
    my  $command = '';
    my ($cnt, $ret) = (0, 1);

    for(my $q = 0; $q < $$rnqueries; $q++) {
        my $qnum = $q;
        if($$rinputs{"${qnum}_rcode"}) {
            if($$rinputs{"${qnum}_error"} || $$rinputs{"${qnum}_errhl"}) {
                Error($$rinputs{"${qnum}_error"}, $$rinputs{"${qnum}_errhl"}, 1);##no e-mail
            }
            next
        }
        unless($$rinputs{"${qnum}_cmdline"}) {
            $$rinputs{"${qnum}_rcode"} = 1;
            $$rinputs{"${qnum}_error"} = "ERROR: $MYPROGNAME: $mysubname: Command line not formed for Record No. ${qnum}\n";
            $$rinputs{"${qnum}_errhl"} = "Invalid command line for Record No. ${qnum}.\n";
            Error($$rinputs{"${qnum}_error"}, $$rinputs{"${qnum}_errhl"}, 1);##no e-mail
            next
        }

        $command = $$rinputs{"${qnum}_cmdline"};

        if(open(F,'>>',$logfilename)){print(F GetTime()."${preamb} ${command}\n\n");close(F);}

        unless(ExecCommand($command)) {
            $$rinputs{"${qnum}_rcode"} = 1;
            $$rinputs{"${qnum}_error"} = "ERROR: $MYPROGNAME: $mysubname: Failed to superimpose structures from Record No. ${qnum}\n";
            $$rinputs{"${qnum}_errhl"} = "Failed to superimpose structures from Record No. ${qnum}.\n";
            Error($$rinputs{"${qnum}_error"}, $$rinputs{"${qnum}_errhl"}, 1);##no e-mail
            next
        }

        $cnt++;
    }

    unless($cnt) {
        $$rerrmsg = "ERROR: $MYPROGNAME: $mysubname: Superposition for all pairs failed.\n";
        $$rhlerrmsg = "Superposition for all pairs failed.\n";
        return 0;
    }

    return $ret;
}

## -----------------------------------------------------------------------------
## MakeResultsList: make and write to a file the list of the results files over 
## all queries
##
sub MakeResultsList
{
    my  $rfilelist = shift;##ref to the list of files listed in the results list file
    my  $dirname = shift;##name of the job directory
    my  $reslstfile = shift;##name of the output results listing file to be generated
    my  $nqueries = shift;##number of queries
    my  $rinputs = shift;##ref to the hash of individual inputs
    my  $rerrmsg = shift;##ref to the error message string to be put in logs
    my  $rhlerrmsg = shift;##ref to the h-l error message string

    my  $mysubname = (caller(0))[3];
    my  $preamb = "[ ${mysubname} ] ";
    my  @qrynums = 0..$nqueries-1;
    my  $ret = 1;

    unless(open(F,'>',$reslstfile)){
        $$rerrmsg = "ERROR: $MYPROGNAME: $mysubname: Failed to open file for writing: '${reslstfile}'.\n";
        $$rhlerrmsg = "File open for writing failed.\n";
        return 0;
    }

    for(my $q = 0; $q <= $#qrynums; $q++) {
        my $qnum = $qrynums[$q];
        ##skip queries with pending errors
        next if($$rinputs{"${qnum}_rcode"});

        my @files = ( $$rinputs{"${qnum}_outfl"} );

        $_ =~ s/^${dirname}\/*(.+)$/$1/ foreach(@files);

        print(F "\"$files[0]\"\n");

        push @$rfilelist, @files;
    }

    close(F);
    return $ret;
}

## -----------------------------------------------------------------------------
## CompressResults: compress all required results files to a single archive
##
sub CompressResults
{
    my  $tarprog = shift;##full pathname to the tar program
    my  $gzprog = shift;##full pathname to the gzip program
    my  $dirname = shift;##name of the job directory
    my  $resultsfile = shift;##name of the output results archive file to be created
    my  $rfilelist = shift;##ref to the list of ($dirname-prefix-removed) files to include in the archive
    my  $logfile = shift;##fullname of the log file to include in the archive
    my  $reslstfile = shift;##fullname of the results listing file to include in the archive
    my  $statusfile = shift;##fullname of the status messages file to include in the archive
    my  $errorfile = shift;##fullname of the high-level error messages file to include in the archive
    my  $rerrmsg = shift;##ref to the error message string to be put in logs
    my  $rhlerrmsg = shift;##ref to the h-l error message string

    my  $mysubname = (caller(0))[3];
    my  $preamb = "[ ${mysubname} ] ";
    my  $command = '';
    my  $ret = 1;

    return 0 unless(AddFileToArchive($tarprog, 
        $resultsfile, $reslstfile, '', $rerrmsg, $rhlerrmsg, 1)); #no dirname; 1==create

    for(my $i = 0; $i <= $#{$rfilelist}; $i++) {
        return 0 unless(AddFileToArchive($tarprog, 
            $resultsfile, $$rfilelist[$i], $dirname, $rerrmsg, $rhlerrmsg)); #dirname given
    }

    unless(AddFileToArchive($tarprog, 
        $resultsfile, $statusfile, '', $rerrmsg, $rhlerrmsg)){ #no dirname
        ##record an error, don't stop and  don't send an e-mail
        Error($$rerrmsg, $$rhlerrmsg, 1);
    }

    unless(AddFileToArchive($tarprog, 
        $resultsfile, $errorfile, '', $rerrmsg, $rhlerrmsg)){ #no dirname
        ##record an error, don't stop and  don't send an e-mail
        Error($$rerrmsg, $$rhlerrmsg, 1);
    }

    $command = "${gzprog} -f \"${resultsfile}\"";

    print(STDERR GetTime()."${preamb} ${command}\n\n");

    unless(ExecCommand($command)) {
        $$rerrmsg = "ERROR: $MYPROGNAME: $mysubname: Failed to gzip the archive: '${resultsfile}'\n";
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

