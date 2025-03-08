package usearch;

## (C) 2024 Mindaugas Margelevicius, Institute of Biotechnology, Vilnius University
## Package for conducting GTalign engine-powered search given

use strict;
use Config;
use threads;
use threads::shared;
use FindBin;
use lib "$FindBin::Bin";
use readcfg;
use readopt;
use File::Spec;
use File::Copy;
use File::Basename;
use Scalar::Util qw(looks_like_number);
use POSIX qw(strftime);

## =============================================================================

sub new {
    my $that = shift;
    my $class = ref($that) || $that;
    my $self;

    ##constants:
    $self->{MYPROGNAME} = __PACKAGE__;
    $self->{GETCHAINprog} = File::Spec->catfile($FindBin::Bin,"getchainlengths.py");
    $self->{PWFA2MSAprog} = File::Spec->catfile($FindBin::Bin,"pwfa2msa.pl");
    $self->{SENDMAILprog} = File::Spec->catfile($FindBin::Bin,"sendemail.pl");
    $self->{CFGFILE} = File::Spec->catfile($FindBin::Bin,File::Spec->updir(),"var","gtalign-ws-backend.conf");
    $self->{devnull} = File::Spec->devnull();
    $self->{TARPROG} = `which tar 2>$self->{devnull}`; chomp($self->{TARPROG});
    $self->{GZPROG} = `which gzip 2>$self->{devnull}`; chomp($self->{GZPROG});
    $self->{MAXNCPUs} = 6;##maximum number of CPU cores assigned for a job
    ##threshold for the number of queries which, if exceeded, implies creating threads, each assigned 1 cpu:
    ##(use a large number for serialization with multiple cores)
    $self->{MAXNQUERIES_multicore} = 6;
    $self->{MAXNQUERIES} = 100;##maximum number of queries allowed in the input

    $self->{MAXSTRLENGTALIGN} = 9999;##maximum structure length for GTalign

    ($self->{QRYEXT}, $self->{FASEXT}, $self->{AFAEXT}, $self->{PWFEXT},
     $self->{PDBEXT}, $self->{CIFEXT}, $self->{OUTEXT}) = 
       ('qry','fa','afa', 'pwfa','pdb','cif','json');

    $self->{cmdline_arguments} = '';
    $self->{gtalign_cacheopt} = '';
    $self->{gtalign_def_options}->{dev_queries_per_chunk} = 2;
    $self->{gtalign_def_options}->{dev_queries_total_length_per_chunk} = 4000;
    $self->{gtalign_def_options}->{dev_max_length} = $self->{gtalign_set_options}->{dev_max_length} = 4000;
    $self->{gtalign_def_options}->{dev_min_length} = $self->{gtalign_set_options}->{dev_min_length} = 20;

    $self->{INPFILENAME} = '';
    $self->{OPTFILENAME} = '';
    $self->{STAFILENAME} = '';
    $self->{COMLOGFILENAME} = '';
    $self->{RESLSTFILENAME} = '';
    $self->{RESULTSFILENAME} = '';
    $self->{ERRFILENAME} = '';
    $self->{NORUN} = 0;

    while( scalar(@_)) {
        $self->{uc( $_[0] )} = $_[1];
        shift, shift;
    }

    bless( $self, $class );
    return $self;
}

## =============================================================================

sub SetOptions
{
    my  $self = shift;
    my  $mysubname = (caller(0))[3];
    my  $class = ref($self) || die("ERROR: $mysubname: Should be called by object.");
    $self->{INPFILENAME} = shift;
    $self->{OPTFILENAME} = shift;
    $self->{STAFILENAME} = shift;
    $self->{COMLOGFILENAME} = shift;
    $self->{RESLSTFILENAME} = shift;
    $self->{RESULTSFILENAME} = shift;
    $self->{ERRFILENAME} = shift;
    $self->{NORUN} = shift;
    $self->{METHOD} = shift;
}

## =============================================================================
##
sub Initialize
{
    my  $self = shift;
    my  $mysubname = (caller(0))[3];
    my  $class = ref($self) || die("ERROR: $mysubname: Should be called by object.");
    $self->Preinitialize();
    $self->ValidateOptionsFile();
    $self->CheckConfig();
    $self->InitializeOptions();
}

## =============================================================================
## preinitilalize and check paths and files
##
sub Preinitialize
{
    my  $self = shift;
    my  $mysubname = (caller(0))[3];
    my  $class = ref($self) || die("ERROR: $mysubname: Should be called by object.");

    print(STDERR "\n\n".GetDatetime()."\n\n");

    my $suffix;
    ($self->{inpbasename},$self->{inpdirname},$suffix) = fileparse($self->{INPFILENAME}, qr/\.[^.]*/);
    $self->{inpdirname} = File::Spec->rel2abs($self->{inpdirname});

    unless($self->{STAFILENAME}) {
        $self->{STAFILENAME} = File::Spec->catfile($self->{inpdirname},"$self->{inpbasename}.status");
        print(STDERR "WARNING: $self->{MYPROGNAME}: Filename for computation progress mesages not given; ".
              "it has been set to: '$self->{STAFILENAME}'\n");
    }

    unless($self->{NORUN}) {
        unless($self->{COMLOGFILENAME}) {
            $self->{COMLOGFILENAME} = File::Spec->catfile($self->{inpdirname},"$self->{inpbasename}__gtalign_out.log");
            print(STDERR "WARNING: $self->{MYPROGNAME}: Filename for GTalign log mesages not given; ".
                  "it has been set to: '$self->{COMLOGFILENAME}'\n");
        }
        unless($self->{RESLSTFILENAME}) {
            $self->{RESLSTFILENAME} = File::Spec->catfile($self->{inpdirname},"$self->{inpbasename}__gtalign_out.lst");
            print(STDERR "WARNING: $self->{MYPROGNAME}: Filename for a results list not given; ".
                  "it has been set to: '$self->{RESLSTFILENAME}'\n");
        }
        unless($self->{RESULTSFILENAME}) {
            $self->{RESULTSFILENAME} = File::Spec->catfile($self->{inpdirname},"$self->{inpbasename}__gtalign_out.tar");
            print(STDERR "WARNING: $self->{MYPROGNAME}: Filename for compressed results not given; ".
                  "it has been set to: '$self->{RESULTSFILENAME}'\n");
        }
    }
    unless($self->{ERRFILENAME}) {
        $self->{ERRFILENAME} = File::Spec->catfile($self->{inpdirname},"$self->{inpbasename}.err");
        print(STDERR "WARNING: $self->{MYPROGNAME}: Filename for high-level error mesages not given; ".
              "it has been set to: '$self->{ERRFILENAME}'\n");
    }
    ## truncate error file before using Error!
    if(open(F,'>',$self->{ERRFILENAME})){close(F);}
    ##unless($result) {
    ##    $self->Error("ERROR: $self->{MYPROGNAME}: Error in command-line arguments.\n",
    ##            "Command-line arguments error.\n");## h-l error message
    ##    $self->MyExit(1);
    ##}
    unless($self->{INPFILENAME} && -f $self->{INPFILENAME}) {
        $self->Error("ERROR: $self->{MYPROGNAME}: Input file not found: '$self->{INPFILENAME}'\n",
                "Input file not found.\n");## h-l error message
        $self->MyExit(1);
    }
    unless($self->{OPTFILENAME} && -f $self->{OPTFILENAME}) {
        $self->Error("ERROR: $self->{MYPROGNAME}: Input job options file not found: '$self->{OPTFILENAME}'\n",
                "Job options file not found.\n");## h-l error message
        $self->MyExit(1);
    }
    unless(-f $self->{CFGFILE}) {
        $self->Error("ERROR: $self->{MYPROGNAME}: Config file not found: '$self->{CFGFILE}'\n",
                "Configuration file not found.\n");## h-l error message
        $self->MyExit(1);
    }
    unless(-f $self->{GETCHAINprog}) {
        $self->Error("ERROR: $self->{GETCHAINprog}: Program file not found: '$self->{GETCHAINprog}'\n",
                "Some of the required program files not found.\n");## h-l error message
        $self->MyExit(1);
    }
    unless(-f $self->{PWFA2MSAprog}) {
        $self->Error("ERROR: $self->{MYPROGNAME}: Program file not found: '$self->{PWFA2MSAprog}'\n",
                "Some of the required program files not found.\n");## h-l error message
        $self->MyExit(1);
    }
    unless($self->{TARPROG} && $self->{GZPROG}) {
        $self->Error("ERROR: $self->{MYPROGNAME}: System programs 'tar' and/or 'gzip' not found.\n",
                "Some of the system programs not found.\n");## h-l error message
        $self->MyExit(1);
    }
}

## =============================================================================
## check configuration variables and initilalize corresponding paths
##
sub CheckConfig
{
    my  $self = shift;
    my  $mysubname = (caller(0))[3];
    my  $class = ref($self) || die("ERROR: $mysubname: Should be called by object.");

    $self->{cfgvar} = readcfg->new($self->{CFGFILE});##backend configuration object

    { ##if( $self->{METHOD} eq 'gtalign') {
        if($self->{cfgvar}->JobMaxNoQueries() < 1) {
            print(STDERR "WARNING: $self->{MYPROGNAME}: Maximum number of queries not given; ".
                  "it has been set to: '$self->{MAXNQUERIES}'\n");
        } else {
            $self->{MAXNQUERIES} = $self->{cfgvar}->JobMaxNoQueries();
        }
    }


    unless(-d $self->{cfgvar}->PathStrDb_PDB()) {
        $self->Error("ERROR: $self->{MYPROGNAME}: Database directory not found: '".$self->{cfgvar}->PathStrDb_PDB()."'\n",
                "Some of the database directories not found.\n");## h-l error message
        $self->MyExit(1);
    }
    unless(-d $self->{cfgvar}->PathStrDb_SCOP()) {
        $self->Error("ERROR: $self->{MYPROGNAME}: Database directory not found: '".$self->{cfgvar}->PathStrDb_SCOP()."'\n",
                "Some of the database directories not found.\n");## h-l error message
        $self->MyExit(1);
    }
    unless(-d $self->{cfgvar}->PathStrDb_ECOD()) {
        $self->Error("ERROR: $self->{MYPROGNAME}: Database directory not found: '".$self->{cfgvar}->PathStrDb_ECOD()."'\n",
                "Some of the database directories not found.\n");## h-l error message
        $self->MyExit(1);
    }
    unless(-d $self->{cfgvar}->PathStrDb_SwissProt()) {
        $self->Error("ERROR: $self->{MYPROGNAME}: Database directory not found: '".$self->{cfgvar}->PathStrDb_SwissProt()."'\n",
                "Some of the database directories not found.\n");## h-l error message
        $self->MyExit(1);
    }
    unless(-d $self->{cfgvar}->PathStrDb_Proteomes()) {
        $self->Error("ERROR: $self->{MYPROGNAME}: Database directory not found: '".$self->{cfgvar}->PathStrDb_Proteomes()."'\n",
                "Some of the database directories not found.\n");## h-l error message
        $self->MyExit(1);
    }
    unless(-d $self->{cfgvar}->PathStrDb_UniRef30()) {
        $self->Error("ERROR: $self->{MYPROGNAME}: Database directory not found: '".$self->{cfgvar}->PathStrDb_UniRef30()."'\n",
                "Some of the database directories not found.\n");## h-l error message
        $self->MyExit(1);
    }


    unless($self->{NORUN}) {
        unless(-d $self->{cfgvar}->InstallDir_GTalign()) {
            $self->Error("ERROR: $self->{MYPROGNAME}: GTalign installation directory not found: '".$self->{cfgvar}->InstallDir_GTalign()."'\n",
                "Some of the installation directories not found.\n");## h-l error message
            $self->MyExit(1);
        }
    }


    ##check particular programs
    $self->{optionvalues}->{prog_gtalign_gtalign} = '';


    unless($self->{NORUN}) {
        $self->{optionvalues}->{prog_gtalign_gtalign} = File::Spec->catfile($self->{cfgvar}->InstallDir_GTalign(),'bin','gtalign');
        unless(-f $self->{optionvalues}->{prog_gtalign_gtalign}) {
            $self->Error("ERROR: $self->{MYPROGNAME}: GTalign executable not found: '".$self->{optionvalues}->{prog_gtalign_gtalign}."'\n",
                "Incomplete software installed on the system.\n");## h-l error message
            $self->MyExit(1);
        }
    }
}

## =============================================================================
#$ helper function for validating values in the options file
sub ValidateHelper
{
    my  $self = shift;
    my  $optname = shift;
    my  $optcurvalue = shift;
    my  $optminvalue = shift;
    my  $optmaxvalue = shift;
    my  $optdefvalue = shift;
    my  $mysubname = (caller(0))[3];
    my  $class = ref($self) || die("ERROR: $mysubname: Should be called by object.");
    my  $warn = '';
    if(!looks_like_number($optcurvalue) || $optcurvalue < $optminvalue || $optcurvalue > $optmaxvalue) {
        $warn = "\nWARNING: Disallowed values: Option $optname changed: $optcurvalue -> $optdefvalue\n";
        $self->Warning($warn, $warn, 1);##no e-mail
        $_ = "$optname=$optdefvalue\n";
    }
}
## validate options file and modify illegal values
##
sub ValidateOptionsFile
{
    my  $self = shift;
    my  $mysubname = (caller(0))[3];
    my  $class = ref($self) || die("ERROR: $mysubname: Should be called by object.");
    my ($contents, $warn, $value) = ('','',0);

    return unless(open(F, $self->{OPTFILENAME}));
    while(<F>) {
        chomp;
        next if /^\s*#/;
        if(/^\s*(-s)\s*=\s*(\S+)/) { $self->ValidateHelper($1, $2, 0.0, 0.999, 0.5); }
        elsif(/^\s*(--sort)\s*=\s*(\S+)/) { $self->ValidateHelper($1, $2, 0, 8, 2); }
        elsif(/^\s*(--nhits)\s*=\s*(\S+)/) { $self->ValidateHelper($1, $2, 1, 10000, 2000); }
        elsif(/^\s*(--nalns)\s*=\s*(\S+)/) { $self->ValidateHelper($1, $2, 1, 10000, 2000); }
        elsif(/^\s*(--pre\-similarity)\s*=\s*(\S+)/) { $self->ValidateHelper($1, $2, 0.0, 100.0, 0.0); }
        elsif(/^\s*(--pre\-score)\s*=\s*(\S+)/) { $self->ValidateHelper($1, $2, 0.0, 0.999, 0.4); }
        elsif(/^\s*(--speed)\s*=\s*(\S+)/) { $self->ValidateHelper($1, $2, 0, 13, 13); }
        ##NOTE: Do not permit changing length limits, as this invalidates caching:
        ##elsif(/^\s*(--dev\-max\-length)\s*=\s*(\S+)/) { $self->ValidateHelper($1, $2, 100, 65535, 4000); $self->{gtalign_set_options}->{dev_max_length} = $2; }
        ##elsif(/^\s*(--dev\-min\-length)\s*=\s*(\S+)/) { $self->ValidateHelper($1, $2, 3, 32767, 20); $self->{gtalign_set_options}->{dev_min_length} = $2; }
        elsif(/^\s*(--dev\-max\-length)\s*=\s*(\S+)/) { $self->ValidateHelper($1, $2, 100, 65535, 4000); next; }
        elsif(/^\s*(--dev\-min\-length)\s*=\s*(\S+)/) { $self->ValidateHelper($1, $2, 3, 32767, 20); next; }
        s/^\s*(\S+)\s*=\s*(\S+).*$/ $1=$2/;
        s/^\s*(\S+).*$/ $1/;
        $contents .= $_ if /^\s*\-/;
    }
    close(F);

    $self->{cmdline_arguments} = $contents;

    ## no error check
    #move($self->{OPTFILENAME},"$self->{OPTFILENAME}.org");
    #unless(open(F, ">", $self->{OPTFILENAME})) {
    #    $self->Error("ERROR: $self->{MYPROGNAME}: $mysubname: ".
    #            "Failed to open job options file: '$self->{OPTFILENAME}'\n",
    #            "Faled to open job options file.\n");## h-l error message
    #    $self->MyExit(1);
    #}
    #print(F $contents);
    #close(F);
}

## =============================================================================
## initialize options for the job
##
sub InitializeOptions
{
    my  $self = shift;
    my  $mysubname = (caller(0))[3];
    my  $class = ref($self) || die("ERROR: $mysubname: Should be called by object.");

    $self->{options} = readopt->new($self->{OPTFILENAME});##options object
    $self->{optionkeys} = {
        job_num_cpus => 'JOB_NUM_CPUS',
        gtalign_db => 'gtalign_db'
    };
    my  $nsyscpus = $self->GetNCPUs();
    $self->{ncpus} = $self->{MAXNCPUs};

##    if($Config{useithreads}) {
        if($self->{cfgvar}->Exists($self->{optionkeys}->{job_num_cpus})) {
            $self->{ncpus} = $self->{cfgvar}->GetValue($self->{optionkeys}->{job_num_cpus});
            $self->{ncpus} = $nsyscpus if($nsyscpus > 0 && $nsyscpus < $self->{ncpus});
        }
        else {
            print(STDERR "WARNING: $self->{MYPROGNAME}: Job option $self->{optionkeys}->{job_num_cpus} not specified: ".
                "Using default #cpus=$self->{ncpus}\n");
        }
##    }
##    else {
##        $self->{ncpus} = 1;
##        print(STDERR "WARNING: $self->{MYPROGNAME}: Perl compiled WITHOUT thread support: ".
##            "Using only one cpu: #cpus=$self->{ncpus}\n");
##    }
}

## =============================================================================
##







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
    my  $self = shift;
    my  $mysubname = (caller(0))[3];
    my  $class = ref($self) || die("ERROR: $mysubname: Should be called by object.");

    unless(open(H, "/proc/cpuinfo")) {
        print(STDERR "WARNING: $self->{MYPROGNAME}: Unable to determine #cpu cores.\n");
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
    my  $self = shift;
    my  $msg = shift;
    my  $mysubname = (caller(0))[3];
    my  $class = ref($self) || die("ERROR: $mysubname: Should be called by object.");
    return 0 unless(open(F,'>>',$self->{STAFILENAME}));
    print(F $msg);
    close(F);
    return 1;
}

## -----------------------------------------------------------------------------
## GetWarnings: extract warnings from a log file
##
sub GetWarnings
{
    my  $self = shift;
    my  $logfile = shift;##search log file
    my  $rlist = shift;##ref to the output list of warnings
    my  $mysubname = (caller(0))[3];
    my  $class = ref($self) || die("ERROR: $mysubname: Should be called by object.");
    unless(open(F, $logfile)) {
        print(STDERR "WARNING: $self->{MYPROGNAME}: Unable to open a search tool's log file: '$logfile'.\n");
        return 0;
    }
    @$rlist = grep(/\sWARNING/,<F>);
    close(F);
    return 1;
}

## -----------------------------------------------------------------------------
## AddFileToArchive: add a file to a given archive
##
sub AddFileToArchive
{
    my  $self = shift;
    my  $filepathname = shift;##full pathname to the file
    my  $dirname = shift;##name of directory where $filepathname is; may be empty
    my  $rerrmsg = shift;##ref to the error message string to be put in logs
    my  $rhlerrmsg = shift;##ref to the h-l error message string
    my  $create = shift;##flag of whether the archive is to be created

    my  $mysubname = (caller(0))[3];
    my  $class = ref($self) || die("ERROR: $mysubname: Should be called by object.");
    my  $preamb = "[ ${mysubname} ] ";
    my  $archive = $self->{RESULTSFILENAME};
    my  $command = '';
    my  $opt = $create?'c':'r';
    my  $ret = 1;

    if($dirname) {
        $command = "$self->{TARPROG} -${opt}f \"${archive}\" -C \"${dirname}\" \"${filepathname}\"";
    }
    else {
        my $filedirname = dirname($filepathname);
        my $filename = basename($filepathname);
        $command = "$self->{TARPROG} -${opt}f \"${archive}\" -C \"${filedirname}\" \"${filename}\"";
    }

    print(STDERR GetTime()."${preamb} ${command}\n");

    unless($self->ExecCommand($command)) {
        $$rerrmsg = "ERROR: $self->{MYPROGNAME}: $mysubname: Failed to add file to archive: '${filepathname}'\n";
        $$rhlerrmsg = "Failed to archive some of the results files.\n";
        return 0;
    }

    return $ret;
}

## -----------------------------------------------------------------------------
## VerifyOptionValues: verify directory and filename information given in the 
## job options file
##
sub VerifyOptionValues 
{
    my  $self = shift;
    my  $rerrmsg = shift;##ref to the error message string
    my  $rhlerrmsg = shift;##ref to the h-l error message string

    my  $mysubname = (caller(0))[3];
    my  $class = ref($self) || die("ERROR: $mysubname: Should be called by object.");
    my ($filename, $fullname);
    my  $ret = 1;


    unless($self->{NORUN}) {
        ##if( $self->{METHOD} eq 'gtalign')
        {
            unless($self->{options}->Exists($self->{optionkeys}->{gtalign_db})) {
                $$rerrmsg = "ERROR: $self->{MYPROGNAME}: $mysubname: ".
                    "Option $self->{optionkeys}->{gtalign_db} not specified in the job options file.\n";
                $$rhlerrmsg = "Structure database not specified.\n";
                return 0;
            }

            $self->{optionvalues}->{$self->{optionkeys}->{gtalign_db}} = '';
            $filename = $self->{options}->GetValue($self->{optionkeys}->{gtalign_db});
            my @fnames = split(',', $filename);
            my(@fullnamelist, @cachenamelist);
            my @predetdbs = (
                $self->{cfgvar}->PathStrDb_PDB(), $self->{cfgvar}->PathStrDb_SCOP(),
                $self->{cfgvar}->PathStrDb_ECOD(), $self->{cfgvar}->PathStrDb_SwissProt(),
                $self->{cfgvar}->PathStrDb_Proteomes(), $self->{cfgvar}->PathStrDb_UniRef30(),
                $self->{cfgvar}->Exists('PATHSTRDB_PDB_SCOP_ECOD')? $self->{cfgvar}->GetValue('PATHSTRDB_PDB_SCOP_ECOD'): '',
                $self->{cfgvar}->Exists('PATHSTRDB_PDB_SCOP_ECOD_sw_prot')? $self->{cfgvar}->GetValue('PATHSTRDB_PDB_SCOP_ECOD_sw_prot'): ''
            );
            my @predetdbnames = (
                $self->{cfgvar}->Exists('STRDB_PDB_NAME')? $self->{cfgvar}->GetValue('STRDB_PDB_NAME'): '',
                $self->{cfgvar}->Exists('STRDB_SCOP_NAME')? $self->{cfgvar}->GetValue('STRDB_SCOP_NAME'): '',
                $self->{cfgvar}->Exists('STRDB_ECOD_NAME')? $self->{cfgvar}->GetValue('STRDB_ECOD_NAME'): '',
                $self->{cfgvar}->Exists('STRDB_SwissProt_NAME')? $self->{cfgvar}->GetValue('STRDB_SwissProt_NAME'): '',
                $self->{cfgvar}->Exists('STRDB_Proteomes_NAME')? $self->{cfgvar}->GetValue('STRDB_Proteomes_NAME'): '',
                $self->{cfgvar}->Exists('STRDB_UniRef30_NAME')? $self->{cfgvar}->GetValue('STRDB_UniRef30_NAME'): '',
                $self->{cfgvar}->Exists('STRDB_PDB_SCOP_ECOD_NAME')? $self->{cfgvar}->GetValue('STRDB_PDB_SCOP_ECOD_NAME'): '',
                $self->{cfgvar}->Exists('STRDB_PDB_SCOP_ECOD_sw_prot_NAME')? $self->{cfgvar}->GetValue('STRDB_PDB_SCOP_ECOD_sw_prot_NAME'): ''
            );
            my @predetdbcachedirs = (
                $self->{cfgvar}->Exists('PATHSTRDBcache_PDB')? $self->{cfgvar}->GetValue('PATHSTRDBcache_PDB'): '',
                $self->{cfgvar}->Exists('PATHSTRDBcache_SCOP')? $self->{cfgvar}->GetValue('PATHSTRDBcache_SCOP'): '',
                $self->{cfgvar}->Exists('PATHSTRDBcache_ECOD')? $self->{cfgvar}->GetValue('PATHSTRDBcache_ECOD'): '',
                $self->{cfgvar}->Exists('PATHSTRDBcache_SwissProt')? $self->{cfgvar}->GetValue('PATHSTRDBcache_SwissProt'): '',
                $self->{cfgvar}->Exists('PATHSTRDBcache_Proteomes')? $self->{cfgvar}->GetValue('PATHSTRDBcache_Proteomes'): '',
                $self->{cfgvar}->Exists('PATHSTRDBcache_UniRef30')? $self->{cfgvar}->GetValue('PATHSTRDBcache_UniRef30'): '',
                $self->{cfgvar}->Exists('PATHSTRDBcache_PDB_SCOP_ECOD')? $self->{cfgvar}->GetValue('PATHSTRDBcache_PDB_SCOP_ECOD'): '',
                $self->{cfgvar}->Exists('PATHSTRDBcache_PDB_SCOP_ECOD_sw_prot')? $self->{cfgvar}->GetValue('PATHSTRDBcache_PDB_SCOP_ECOD_sw_prot'): ''
            );
            foreach $filename(@fnames) {
                $fullname = '';

                ##foreach my $path(@predetdbs) {
                for(my $pi = 0; $pi <= $#predetdbs; $pi++) {
                    my $path = $predetdbs[$pi];
                    my $name = ($pi <= $#predetdbnames)? $predetdbnames[$pi]: '';
                    my $cachepath = ($pi <= $#predetdbcachedirs)? $predetdbcachedirs[$pi]: '';
                    next unless $path;
                    next unless $name eq $filename;
                    my @paths = split(',', $path);
                    my @filenamesplit = split(/\|/, $filename);
                    next unless $#paths == $#filenamesplit;
                    for(my $fi = 0; $fi <= $#paths; $fi++) {
                        my $tmpfname = File::Spec->catfile($paths[$fi], $filenamesplit[$fi]);
                        next unless -d $tmpfname;
                        $fullname .= ',' if $fullname;
                        $fullname .= $tmpfname;
                    }
                    next unless $fullname;
                    push @fullnamelist, $fullname;
                    push @cachenamelist, $cachepath if $cachepath;
                    last;
                }

                unless($fullname) {
                    my $text = "WARNING: Structure database not found: $filename\n";
                    $self->Warning($text, $text, 1);##no e-mail
                }
            }

            if($#fullnamelist < 0) {
                $$rerrmsg = "ERROR: $self->{MYPROGNAME}: $mysubname: ".
                  "All structure db(s) '${filename}' not found in local directories.\n";
                $$rhlerrmsg = "Structure database(s) not found.\n";
                return 0;
            }

            $self->{optionvalues}->{$self->{optionkeys}->{gtalign_db}} = join(',',@fullnamelist);

            if($#fullnamelist == 0 && $#cachenamelist == 0 && length($cachenamelist[0]) > 0) {
                $self->{gtalign_cacheopt} = "-c \"$cachenamelist[0]\"";
            }
        }
    }


    return $ret;
}

## -----------------------------------------------------------------------------
## GenGTalignHPCqueryOptions: generate query-specific part of the HPC options 
## for GTalign
##
sub GenGTalignHPCqueryOptions
{
    my  $self = shift;
    my  $rnqueries = shift;##ref to the number of queries
    my  $rinfmsg = shift;##ref to an information message
    my  $rerrmsg = shift;##ref to the error message string to be put in logs
    my  $rhlerrmsg = shift;##ref to the h-l error message string

    my  $mysubname = (caller(0))[3];
    my  $class = ref($self) || die("ERROR: $mysubname: Should be called by object.");
    my  $ret = 1;

    unless(-f $self->{INPFILENAME}) {
        $$rerrmsg = "ERROR: $self->{MYPROGNAME}: $mysubname: Input file not found: '$self->{INPFILENAME}'\n";
        $$rhlerrmsg = "Input file not found.\n";
        return 0;
    }

    my  $cress = `python3 $self->{GETCHAINprog} -i $self->{INPFILENAME}`;
    my ($maxlen0, $maxlen1, $maxlens, @lens) = (0,0,0);
    my ($ncif, $npdb) = (0,0);
    ##no error check
    foreach(split(/\n/, $cress)) {
        if(/^D:\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)/) {
            my ($lcfile, $lcname, $lcchain, $lclength, $lcmodel) = ($1,$2,$3,$4,$5);
            push(@lens, $lclength);
            if($lcname =~ /^None$/i) {$ncif++} else {$npdb++}
            if($lclength < 3) {
                my $warntext = "WARNING: Input file $lcfile: Chain $lcchain ";
                $warntext .= "(Model $lcmodel) " if 0 < $lcmodel;
                if($lclength < 1) {$warntext .= "contains no amino acids.\n";}
                else {$warntext .= "with <3 amino acids ignored.\n";}
                $self->Warning($warntext, $warntext, 1);##no e-mail
            }
        }
    }
    @lens = sort {$b <=> $a} @lens;
    $maxlen0 = $lens[0] if 0 <= $#lens;
    $maxlen1 = $lens[1] if 1 <= $#lens;
    if($maxlen0 + 100 < 50000) { $maxlen0 += 100; }##some error in interpreting structure files
    $maxlens = $maxlen0 + $maxlen1;

    my ($totlen, $nqrsperchunk) = (
        $self->{gtalign_def_options}->{dev_queries_total_length_per_chunk},
        $self->{gtalign_def_options}->{dev_queries_per_chunk}
    );

    my  $defdbsizeopts =
        ($self->{gtalign_set_options}->{dev_max_length} == $self->{gtalign_def_options}->{dev_max_length}) &&
        ($self->{gtalign_set_options}->{dev_min_length} == $self->{gtalign_def_options}->{dev_min_length});

    if(4000 < $maxlen0) {
        $totlen = $maxlen0;
        $nqrsperchunk = 1;
    } else {
        $nqrsperchunk = ($#lens < 1)? 1: 2;
        if(0 < $maxlens && $maxlens <= 4000) { $totlen = $maxlens; }
        else { $totlen = 4000; }
    }

    ## if cache can be used and db sizes are default values, use cache even if total length is less than default:
    ##if((length($self->{gtalign_cacheopt}) > 0) && $defdbsizeopts &&
    ##   $totlen <= $self->{gtalign_def_options}->{dev_queries_total_length_per_chunk})
    ##NOTE: always use cache if given
    if((length($self->{gtalign_cacheopt}) > 0) && $defdbsizeopts)
    {
        $totlen = $self->{gtalign_def_options}->{dev_queries_total_length_per_chunk};
        $nqrsperchunk = $self->{gtalign_def_options}->{dev_queries_per_chunk};
    }

    $self->{cmdline_arguments} .= " --dev-queries-total-length-per-chunk=${totlen} --dev-queries-per-chunk=${nqrsperchunk}";

    unless($totlen == $self->{gtalign_def_options}->{dev_queries_total_length_per_chunk} &&
           $nqrsperchunk == $self->{gtalign_def_options}->{dev_queries_per_chunk} && $defdbsizeopts)
    {
        $self->{gtalign_cacheopt} = '';
    }

    $$rnqueries = $#lens + 1;

    if($#lens < 0) {
        $$rerrmsg = $$rhlerrmsg = '' if $ret;
        $$rerrmsg .= "ERROR: $self->{MYPROGNAME}: $mysubname: No queries found in the input.\n";
        $$rhlerrmsg .= "Invalid input: No structure queries.\n";
        return 0;
    }

    $$rinfmsg = "$npdb structure(s) in PDB format" if 0 < $npdb;
    $$rinfmsg .= '; ' if(0 < $npdb && 0 < $ncif);
    $$rinfmsg = "$ncif structure(s) in mmCIF format" if 0 < $ncif;
    $$rinfmsg .= '.' if(0 < $npdb || 0 < $ncif);

    return $ret;
}

## -----------------------------------------------------------------------------
## PreprocessInput: preprocess input for specific information;
##
sub PreprocessInput
{
    my  $self = shift;
    my  $rnqueries = shift;##ref to the number of queries
    my  $rinfmsg = shift;##ref to an information message
    my  $rerrmsg = shift;##ref to the error message string to be put in logs
    my  $rhlerrmsg = shift;##ref to the h-l error message string

    my  $mysubname = (caller(0))[3];
    my  $class = ref($self) || die("ERROR: $mysubname: Should be called by object.");
    my  $ret = 1;

    unless($self->GenGTalignHPCqueryOptions($rnqueries, $rinfmsg, $rerrmsg, $rhlerrmsg)) {
        return 0;
    }

    return $ret;
}

## -----------------------------------------------------------------------------
## CheckForErrorsandWarnings: check for error codes and any messages that could 
## have been returned from threads
##
sub CheckForErrorsandWarnings
{
    my  $self = shift;
    my  $nqueries = shift;##number of queries
    my  $rinputs = shift;##ref to the hash of individual inputs

    my  @qrynums = 0..$nqueries-1;

    for(my $q = 0; $q <= $#qrynums; $q++) {
        my $qnum = $qrynums[$q];
        ##print errors issued by threads if any
        if($$rinputs{"${qnum}_rcode"}) { ##don't send an e-mail
            $self->Error($$rinputs{"${qnum}_error"}, $$rinputs{"${qnum}_errhl"}, 1);
            next;
        }
        elsif($$rinputs{"${qnum}_errhl"}) {
            ##no error, but messages can be present; print them
            $self->Warning($$rinputs{"${qnum}_error"}, $$rinputs{"${qnum}_errhl"}, 1);##no e-mail
        }
    }

}

## -----------------------------------------------------------------------------
## RunGTalign: Run GTalign search
##
sub RunGTalign
{
    my  $self = shift;
    my  $rnqueries = shift;##ref to the number of queries
    my  $rinputs = shift;##ref to the hash of individual inputs
    my  $rerrmsg = shift;##ref to the error message string to be put in logs
    my  $rhlerrmsg = shift;##ref to the h-l error message string

    my  $mysubname = (caller(0))[3];
    my  $class = ref($self) || die("ERROR: $mysubname: Should be called by object.");
    my  $preamb = "[ ${mysubname} ] ";
    my  $outputpfx = File::Spec->catfile($self->{inpdirname},"$self->{inpbasename}__gtalign_out");
    my  $outputsubdir = $outputpfx;
    ##no. trials to launch gtalign; NOTE: may fail when two or more processes request resources at the same time:
    my  $ntrials = 3;
    my  $command = '';
    my (@warnings, @resfiles);
    my  $ret = 1;

    my $configopts = '';
    if($self->{cfgvar}->Exists('JOB_GPU_MEM')) {
        my $valuemb = $self->{cfgvar}->GetValue('JOB_GPU_MEM') * 1024;
        $configopts = " --dev-mem=${valuemb}" if $valuemb > 256;
    }

    ##if($self->{cfgvar}->Exists('JOB_NUM_CPUS')) {
    ##    my $value = $self->{cfgvar}->GetValue('JOB_NUM_CPUS');
    ##    $configopts .= " --cpu-threads-reading=${value}";
    ##}
    ## NOTE: $self->{ncpus} has already been validated!
    $configopts .= " --cpu-threads-reading=$self->{ncpus}";

    $command = "$self->{optionvalues}->{prog_gtalign_gtalign} -v ".
            "--qrs=\"$self->{INPFILENAME}\" --rfs=\"".$self->{optionvalues}->{$self->{optionkeys}->{gtalign_db}}.
            "\" $self->{gtalign_cacheopt} -o \"${outputsubdir}\" --outfmt=1 ${configopts} --ter=0 --split=2 ".
            "$self->{cmdline_arguments} >\"$self->{COMLOGFILENAME}\" 2>&1";

    for(my $i=1;; $i++) {
        print(STDERR GetTime()."${preamb} ${command}\n\n");
        unless($self->ExecCommand($command)) {
            do{ sleep(2); next} if $i < $ntrials;
            $$rerrmsg = "ERROR: $self->{MYPROGNAME}: $mysubname: Unable to conduct a GTalign search.\n";
            $$rhlerrmsg = "Unable to conduct a GTalign search. Please make sure your input has correct format.\n";
            return 0;
        }
        last;
    }

    ##extract warnings if any; no check of return code
    $self->GetWarnings($self->{COMLOGFILENAME}, \@warnings);

    if(0 <= $#warnings) {
        my $lst = join("", @warnings);
        my $text = "\nWarnings from the GTalign search:\n\n$lst\n\n";
        $self->Warning($text, $text, 1);##no e-mail
    }

    ##associate gtalign results with particular queries
    $command = '';

    if($self->ReadFiles($outputsubdir, $self->{OUTEXT}, \@resfiles) <= 0) {
        ;##no files found in directory
    }

    $$rnqueries = ($#resfiles + 1);

    for(my $q = 0; $q <= $#resfiles; $q++) {
        my $qnum = $q;
        my $resfilename = $resfiles[$q];
        ##skip queries with pending errors and which have not been searched for
        next if(exists $$rinputs{"${qnum}_rcode"} && $$rinputs{"${qnum}_rcode"});

        my $outfile = File::Spec->catfile($outputsubdir, $resfilename);

        unless($resfilename && -f $outfile) {
            $$rinputs{"${qnum}_rcode"} = 1;
            $$rinputs{"${qnum}_error"} = "ERROR: $self->{MYPROGNAME}: $mysubname: GTalign results for query No.${qnum} not found: $outfile\n";
            $$rinputs{"${qnum}_errhl"} = "GTalign results for query No.${qnum} not found.\n";
            $self->Error($$rinputs{"${qnum}_error"}, $$rinputs{"${qnum}_errhl"}, 1);##no e-mail
            next
        }

        my ($qryfulldesc, $qryfullname, $qrybasename, $qrylength, $qsec) = ('','','',0,0);

        unless(open(F, $outfile)) {
            $$rinputs{"${qnum}_rcode"} = 1;
            $$rinputs{"${qnum}_error"} = "ERROR: $self->{MYPROGNAME}: $mysubname: Failed to open GTalign results for query No.${qnum}: $outfile\n";
            $$rinputs{"${qnum}_errhl"} = "Failed to open GTalign results for query No.${qnum}.\n";
            $self->Error($$rinputs{"${qnum}_error"}, $$rinputs{"${qnum}_errhl"}, 1);##no e-mail
            next
        }

        while(<F>) {
            next if /^\s*$/;
            $qsec = 1 if /\s+"query":\s*{/;
            next unless $qsec;
            $qrylength = $1 if /^\s+"length":\s*(\d+),/;
            $qryfulldesc = $1 if /^\s+"description":\s*"(.+)"$/;
            last if /^\s+"database":\s*{/;
        }

        close(F);

        ##$qryfulldesc =~ s/\n+/ /g;
        $qryfulldesc =~ s/\\//g;
        $qryfullname = $qrybasename = $qryfulldesc;
        my $chpos = index($qryfulldesc,' Chn:');
        my $mpos = index($qryfulldesc,' (M:');
        my $pos = (0 < $chpos && 0 < $mpos)? (($chpos < $mpos)? $chpos: $mpos): ((0 < $chpos)? $chpos: $mpos);
        $qryfullname = substr($qryfulldesc, 0, $pos) if 0 < $pos;
        $qrybasename = basename($qryfullname);

        $$rinputs{"${qnum}_rcode"} = 0;
        $$rinputs{"${qnum}_error"} = '';##error message
        $$rinputs{"${qnum}_errhl"} = '';##high-level error message
        $$rinputs{"${qnum}_descr"} = $qryfulldesc;
        $$rinputs{"${qnum}_fname"} = $qryfullname;
        $$rinputs{"${qnum}_bname"} = $qrybasename;
        $$rinputs{"${qnum}_lengt"} = $qrylength;
        $$rinputs{"${qnum}_outfl"} = $command = $outfile;
    }

    unless($command) {
        $$rerrmsg = "ERROR: $self->{MYPROGNAME}: $mysubname: No GTalign results produced for all queries.\n";
        $$rhlerrmsg = "No GTalign results produced for any of the queries.\n";
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
    my  $self = shift;
    my  $rfilelist = shift;##ref to the list of files listed in the results list file
    my  $nqueries = shift;##number of queries
    my  $rinputs = shift;##ref to the hash of individual inputs
    my  $rerrmsg = shift;##ref to the error message string to be put in logs
    my  $rhlerrmsg = shift;##ref to the h-l error message string

    my  $mysubname = (caller(0))[3];
    my  $class = ref($self) || die("ERROR: $mysubname: Should be called by object.");
    my  $preamb = "[ ${mysubname} ] ";
    my  @qrynums = 0..$nqueries-1;
    my  $ret = 1;

    unless(open(F,'>',$self->{RESLSTFILENAME})){
        $$rerrmsg = "ERROR: $self->{MYPROGNAME}: $mysubname: Failed to open file for writing: '$self->{RESLSTFILENAME}'.\n";
        $$rhlerrmsg = "Opening file for writing failed.\n";
        return 0;
    }

    print(F "# Search_results Query_full_description Query_length\n");

    for(my $q = 0; $q <= $#qrynums; $q++) {
        my $qnum = $qrynums[$q];
        ##skip queries with pending errors
        next if($$rinputs{"${qnum}_rcode"});

        my @files = ( $$rinputs{"${qnum}_outfl"},
                      $$rinputs{"${qnum}_descr"},
                      $$rinputs{"${qnum}_lengt"}
        );

        ##do{$_ = '' unless -f $_} foreach @files;
        $files[0] = '' unless(-f $files[0]);

        ##$_ =~ s/^$self->{inpdirname}\/*(.+)$/$1/ foreach @files;
        $files[0] =~ s/^$self->{inpdirname}\/*(.+)$/$1/;

        print(F "\"$files[0]\"\t\"$files[1]\"\t\"$files[2]\"\n");

        ##push @$rfilelist, @files;
        push @$rfilelist, $files[0];
    }

    close(F);
    return $ret;
}

## -----------------------------------------------------------------------------
## CompressResults: compress all required results files to a single archive
##
sub CompressResults
{
    my  $self = shift;
    my  $rfilelist = shift;##ref to the list of ($dirname-prefix-removed) files to include in the archive
    my  $rerrmsg = shift;##ref to the error message string to be put in logs
    my  $rhlerrmsg = shift;##ref to the h-l error message string

    my  $mysubname = (caller(0))[3];
    my  $class = ref($self) || die("ERROR: $mysubname: Should be called by object.");
    my  $preamb = "[ ${mysubname} ] ";
    my  $command = '';
    my  $ret = 1;

    return 0 unless($self->AddFileToArchive($self->{COMLOGFILENAME}, '', $rerrmsg, $rhlerrmsg, 1)); #no dirname; 1==create

    return 0 unless($self->AddFileToArchive($self->{RESLSTFILENAME}, '', $rerrmsg, $rhlerrmsg)); #no dirname

    for(my $i = 0; $i <= $#{$rfilelist}; $i++) {
        next unless $$rfilelist[$i];
        return 0 unless($self->AddFileToArchive($$rfilelist[$i], $self->{inpdirname}, $rerrmsg, $rhlerrmsg)); #dirname given
    }

    unless($self->AddFileToArchive($self->{STAFILENAME}, '', $rerrmsg, $rhlerrmsg)){ #no dirname
        ##record an error, don't stop and  don't send an e-mail
        $self->Error($$rerrmsg, $$rhlerrmsg, 1);
    }

    unless($self->AddFileToArchive($self->{ERRFILENAME}, '', $rerrmsg, $rhlerrmsg)){ #no dirname
        ##record an error, don't stop and  don't send an e-mail
        $self->Error($$rerrmsg, $$rhlerrmsg, 1);
    }

    $command = "$self->{GZPROG} -f \"$self->{RESULTSFILENAME}\"";

    print(STDERR GetTime()."${preamb} ${command}\n\n");

    unless($self->ExecCommand($command)) {
        $$rerrmsg = "ERROR: $self->{MYPROGNAME}: $mysubname: Failed to gzip the archive: '$self->{RESULTSFILENAME}'\n";
        $$rhlerrmsg = "Failed to compress the archive of results files.\n";
        return 0;
    }

    return $ret;
}



## -----------------------------------------------------------------------------
## CreateLink: create a symbolic link to file
##
sub CreateLink
{
    my  $self = shift;
    my  $oldfile = shift;##source file
    my  $newfile = shift;##destination/new file (link)
    my  $rerrmsg = shift;##ref to the error message string to be put in logs
    my  $rhlerrmsg = shift;##ref to the h-l error message string

    my  $mysubname = (caller(0))[3];
    my  $class = ref($self) || die("ERROR: $mysubname: Should be called by object.");
    my  $preamb = "[ ${mysubname} ] ";
    my  $command = '';

    my $symlink_exists = eval {symlink("",""); 1};

    ##create a link to the msa file
    unless(-f $newfile) {
        ##if($symlink_exists) {
            ##unless(symlink($oldfile, $newfile)) {
            $command = "ln -s $oldfile $newfile";
            print(STDERR GetTime()."${preamb} ${command}\n\n");
            unless($self->ExecCommand($command)) {
                $$rerrmsg = "ERROR: $self->{MYPROGNAME}: $mysubname: Failed to create a link to a file: '$oldfile' -> '$newfile'\n";
                $$rhlerrmsg = "Creating a link failed.\n";
                return 0;
            }
        ##}
    }
    return 1;
}







## THREADS =====================================================================
##







## General =====================================================================
##
##







## =============================================================================
## MyExit: print a code to file and exit
##
sub MyExit
{
    my $self = shift;
    my $ecode = shift;##exit code
    my $mysubname = (caller(0))[3];
    my $class = ref($self) || die("ERROR: $mysubname: Should be called by object.");
    if(open(F,">>",$self->{ERRFILENAME})){print(F "\n${ecode}\n");close(F);}
    exit($ecode);
}

## -------------------------------------------------------------------
## Warning: output a warning message and optionally, send an email
##
sub Warning
{
    my ($self, $msg, $hlerrmsg, $nosend) = @_;
    my $mysubname = (caller(0))[3];
    my $class = ref($self) || die("ERROR: $mysubname: Should be called by object.");
    return $self->Issue($msg, $hlerrmsg, $nosend,
        '',
        "Warning from gtalign-ws ($self->{MYPROGNAME})");
}

## Error: output an error message and optionally, send an email
##
sub Error
{
    my ($self, $msg, $hlerrmsg, $nosend) = @_;
    my $mysubname = (caller(0))[3];
    my $class = ref($self) || die("ERROR: $mysubname: Should be called by object.");
    return $self->Issue($msg, $hlerrmsg, $nosend,
        "ERROR: Issue on the server's backend side: ",
        "Error message from gtalign-ws ($self->{MYPROGNAME})");
}

## Issue: output an issue message and send an email
##
sub Issue
{
    my ($self, $msg, $hlerrmsg, $nosend,  $leadtxt1, $sbjct) = @_;
    my $mysubname = (caller(0))[3];
    my $class = ref($self) || die("ERROR: $mysubname: Should be called by object.");
    print(STDERR GetTime().$msg);
    if($self->{ERRFILENAME} && $hlerrmsg) {
        if(open(EF,">>",$self->{ERRFILENAME})) {
            print(EF $leadtxt1.$hlerrmsg);
            close(EF);
        } else {
            print(STDERR "ERROR: $self->{MYPROGNAME}: Write to error file skipped ".
                "due to the fail to open the file: '$self->{ERRFILENAME}'\n");
        }
    }
    return 1 if $nosend;
    unless(-f $self->{SENDMAILprog}) {
        print(STDERR "ERROR: $self->{MYPROGNAME}: Send-mail program not found: $self->{SENDMAILprog}\n");
        return 0;
    }
    my  $command = "$self->{SENDMAILprog} --sub \"$sbjct\" --body \"$msg\" 2>&1";
    print(STDERR "$command\n");
    return $self->ExecCommand($command);
}

## -------------------------------------------------------------------
## try read lock files in a directory and return 1 if locks exist,
## 0 otherwise, -1 on error
##

sub ReadLocks
{
    my  $self = shift;
    my  $mysubname = (caller(0))[3];
    my  $class = ref($self) || die("ERROR: $mysubname: Should be called by object.");
    return $self->ReadFiles(shift, shift);
}

sub ReadFiles
{
    my  $self = shift;
    my  $dirname = shift;
    my  $pattern = shift;  ## pattern of files to look for as locks
    my  $reffiles = shift; ## reference to vector of files
    my  $mysubname = (caller(0))[3];
    my  $class = ref($self) || die("ERROR: $mysubname: Should be called by object.");
    my  @files;
    my  $locref = defined($reffiles)? $reffiles: \@files;

    unless( opendir(DIR, $dirname)) {
        printf( STDERR "ERROR: Cannot open directory $dirname.\n" );
        return -1;
    }
    @{$locref} = grep { /$pattern$/ && -f File::Spec->catfile($dirname,$_) } readdir(DIR);
    closedir(DIR);
    return 1 if 0 <= $#{$locref};
    return 0;
}

## -------------------------------------------------------------------
## ExecCommand: execute system command
##
sub CheckStatus
{
    my  $self = shift;
    my  $mysubname = (caller(0))[3];
    my  $class = ref($self) || die("ERROR: $mysubname: Should be called by object.");
    $self->ExecCommand();
}

sub ExecCommand
{
    my  $self = shift;
    my  $cmdline = shift;
    my  $returnstatus = shift;##if defined, return exit status of the command
    my  $mysubname = (caller(0))[3];
    my  $class = ref($self) || die("ERROR: $mysubname: Should be called by object.");

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
## -----------------------------------------------------------------------------

1;

