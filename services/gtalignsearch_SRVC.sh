#!/bin/bash

## (C) 2024 Mindaugas Margelevicius, Institute of Biotechnology, Vilnius University
## Submit a gtalign search job to a workload manager's queue

dirname="$(dirname $0)"
[[ "${dirname:0:1}" != "/" ]] && dirname="$(pwd)/$dirname"
basename="$(basename $0)"
username="$(whoami)"
SUBPROG="$(which sbatch 2>/dev/null)"
SINFOPROG="$(which sinfo 2>/dev/null)"
SQUEUEPROG="$(which squeue 2>/dev/null)"

usage="
Submit a gtalign search job to a workload manager's queue.
(C)2024 Mindaugas Margelevicius, Institute of Biotechnology, Vilnius University

Usage:
$basename <Options>

Options:

-i <id>        Web service job identifier which also corresponds to the 
               name pattern of input file <id>.tar and 
               options file <id>.options
-p <directory> Directory where input files for this job can be found.
-b <program>   Pathname to a workload manager's job submitter.
       default=$SUBPROG
-h             short description.
"

## -------------------------------------------------------------------
## Functions
## =========
read -r -d '' MYFUNCS <<'EOF'
## SendEmail: send an email
## avoid single quotes for being able to use the same code in parent and child shells
##
function SendEmail() {
    sendmailprog="$1"
    errfl="$2"
    msg="$3"
    if [[ -f "${sendmailprog}" ]]; then
        cmd="${sendmailprog} --sub \"Error message from gtalign-ws\" --body \"${msg}\" >>${errfl} 2>&1"
        echo "${cmd}" >>${errfl}; eval "${cmd}"
    fi
}
EOF
eval "${MYFUNCS}"
## -------------------------------------------------------------------


while getopts "i:p:b:h" Option
do
    case $Option in
        i ) ID=${OPTARG} ;;
        p ) UPATH=${OPTARG} ;;
        b ) SUBPROG=${OPTARG} ;;
        h ) echo "$usage"; exit 0 ;;
        * ) echo ERROR: Unrecognized argument. >&2; exit 1 ;;
    esac
done
shift $(( $OPTIND - 1 ))

fail=""

if [[ -z "$fail" && -z "$ID" ]]; then fail+="ERROR: No ID given.\n"; fi
if [[ -z "$fail" && -z "$UPATH" ]]; then fail+="ERROR: $ID: No directory given.\n"; fi
if [[ -z "$fail" && ! -d "$UPATH" ]]; then fail+="ERROR: $ID: Directory does not exist: \"$UPATH\".\n"; fi
if [[ -z "$fail" && ( -z "$SUBPROG" || ! -f "$SUBPROG" ) ]]; then 
  fail+="ERROR: $ID: Job submitter program not found: \"$SUBPROG\".\n"
fi

[[ "${UPATH:0:1}" != "/" ]] && UPATH="$(pwd)/$UPATH"

filename="$UPATH/$ID"
tempfile="$UPATH/.$(date +%m%d%Y%H%M%S_$$)"

## bash stdout/stderr file of the job submitted to a job manager's queue:
submitoutfile="${filename}.submit"
shoutfile="${filename}.pbstdout"
errfile="${filename}.pbserr"
inpfile="${filename}.tar"
orgfile="${filename}.tar.org"
optfile="${filename}.options"
statusfile="${filename}.status"
gtalignlogfile="${filename}__gtalign_out.log"
reslstfile="${filename}__gtalign_out.lst"
resultsfile="${filename}__gtalign_out.tar"
errorfile="${filename}.err"
#outfile=${filename}.out

sendemail="${dirname}/../bin/sendemail.pl"
gtalignpl="${dirname}/../bin/gtalign_wrapper.pl"
#getids="${dirname}/../bin/getids.pl"
cfgfile="${dirname}/../var/gtalign-ws-backend.conf"

ncpus=1
ncpusyst=$(nproc)
lspci="/usr/sbin/lspci"
if [[ ! -f "${lspci}" && -n "$(which lspci)" ]]; then
  lspci="$(which lspci)"
fi
if [ -f "${lspci}" ]; then
  ngpusyst=$($lspci|grep -i 'nvidia'|grep -i 'vga\|3d\|2d'|wc -l)
fi

SLURM_OPTS='' ##additional SLURM options

if [[ ! -f "${gtalignpl}" ]]; then fail+="ERROR: $ID: Program not found: ${gtalignpl}\n"; fi
#if [[ ! -f "${getids}" ]]; then fail+="ERROR: $ID: Program not found: ${getids}\n"; fi
if [[ ! -f "${cfgfile}" ]]; then fail+="ERROR: $ID: Config file not found: ${cfgfile}\n"; fi

## determine the allowed number of cpus per job
if [ -z "$fail" ]; then
  ncpusperjob=$(perl -e 'while(<>){if(/^\s*JOB_NUM_CPUS\s*=\s*(\d+)/){print "$1";last}}' "${cfgfile}")
  ngpusystspec=$(perl -e 'while(<>){if(/^\s*SYST_NUM_GPUS\s*=\s*(\d+)/){print "$1";last}}' "${cfgfile}")
  partition=$(perl -e 'while(<>){if(/^\s*JOB_SLURM_PARTITION\s*=\s*(\w+)/){print "$1";last}}' "${cfgfile}")
  njobspergpu=$(perl -e 'while(<>){if(/^\s*JOB_NUM_PER_GPU\s*=\s*(\d+)/){print "$1";last}}' "${cfgfile}")
  #
  if [ -n "${ngpusyst}" ]; then
    ngpusyst_f=${ngpusyst}
  fi
  if [[ -n "${ngpusystspec}" && ${ngpusystspec} -gt 0 ]]; then
    if [[ -n "${ngpusyst_f}" && ${ngpusyst_f} -gt 0 ]]; then
      if [ ${ngpusystspec} -lt ${ngpusyst_f} ]; then ngpusyst_f=${ngpusystspec}; fi
    else
      ngpusyst_f=${ngpusystspec}
    fi
  fi
  #
  ncpusyst_f=1
  if [ -n "${ncpusyst}" ]; then
    ncpusyst_f=${ncpusyst} #$((ncpusyst/2 + 1))
  fi
  if [ -n "${partition}" ]; then
    SLURM_OPTS+=" -p ${partition}"
    if [ -f "${SINFOPROG}" ]; then
      ## get #cpus allocated for this particular SLURM partition
      ncpusyst_f=$(${SINFOPROG} -h -p "${partition}" -o%c)
    fi
  fi
  if [[ -n "${ngpusyst_f}" && ${ngpusyst_f} -gt 0 ]]; then
    ncpusyst_f_org=${ncpusyst_f}
    ## fraction of cpus is divided among GPUs
    if [[ -n "${njobspergpu}" && ${njobspergpu} -gt 0 && ${njobspergpu} -lt 10 ]]; then
      ngpunits=$((ngpusyst_f*njobspergpu))
    else
      ngpunits=${ngpusyst_f}
    fi
    ncpusyst_f=$(( ncpusyst_f/ngpunits))
    ## make sure that assigned #cpus will not exhaust system resources (avoid too small number)
    if [ $((ncpusyst_f_org/ncpusyst_f)) -gt ${ngpunits} ]; then ncpusyst_f=$((ncpusyst_f+1)); fi
  fi
  if [ ${ncpusyst_f} -lt 1 ]; then ncpusyst_f=1; fi
  #
  if [ -n "${ncpusperjob}" ]; then
    if [ ${ncpusperjob} -lt ${ncpusyst_f} ]; then
      ## avoid assigning too small #cpus for a job, as numerous jobs can 
      ## halt the system by requiring too much resources
      ncpus=${ncpusyst_f}
    else
      ncpus=${ncpusperjob}
    fi
  fi
fi

if [[ -n "${username}" && -n "${SQUEUEPROG}" ]]; then
  ##NOTE: this is a temporary solution when other users have in the queue jobs
  ##NOTE: running on one cpu (#total cpus is not divisble by #gpus)
  ntotjobs=$(${SQUEUEPROG} -h -O username|uniq|grep -vE "${username}"|wc -l)
  if [ ${ntotjobs} -gt 0 ]; then ((ncpus++)); fi
fi

if [ -n "$fail" ]; then
  echo -e "$fail" >$errfile
  echo -e "ERROR: The server's backend issue: Invalid data.\n\n1" >>"${errorfile}"
  SendEmail "${sendemail}" "${errfile}" "${UPATH} ${ID}\n\n${fail}"
  exit 1
fi


cmd="echo Submitting... >>\"${statusfile}\""
eval "${cmd}"

mem=10240 ##10GB
##ncpus=20 ##!!##

## set the maximum execution time to 1 day
##
$SUBPROG ${SLURM_OPTS} -n ${ncpus} -t 1440 -J gta-ws-search -o "${shoutfile}" \
  --mem-per-cpu=${mem} >"${submitoutfile}" 2>&1 \
             <<EOF
#!/bin/bash
echo "Job started: \$(date); Running on: \$(hostname)"
echo "JOB_ID: \${SLURM_JOB_ID}"
if [ -n "${partition}" ]; then echo "PARTITION: ${partition}"; fi
if [ -n "${ncpusperjob}" ]; then echo "JOB_NUM_CPUS (from config.): ${ncpusperjob}"; fi
echo "CPUs: ${ncpus} assigned (can be used as a limiting factor)"
if [ -n "${ncpusyst}" ]; then echo "  ${ncpusyst} CPU(s) on system"; fi
if [[ -n "${ngpusyst}" && ${ngpusyst} -gt 0 ]]; then echo "  ${ngpusyst} GPU(s) on system"; fi
if [[ -n "${ngpusyst_f}" && ${ngpusyst_f} -gt 0 ]]; then echo "  ${ngpusyst_f} GPU(s) dedicated"; fi
if [[ -n "${ngpunits}" && ${ngpunits} -gt 0 ]]; then echo "  ${ngpunits} GPU logical unit(s)"; fi
echo "## check the presence of the job in the SLURM queue by typing:"
echo "##   squeue -h -j <JOB_ID> -o%A"
echo

if [[ -n "${ngpunits}" && ${ngpunits} -gt 0 ]]; then
  res=\$((SLURM_JOB_ID - SLURM_JOB_ID/${ngpunits}*${ngpunits}))
  ## sleep before running each job so that the gpus are not requested at the same time
  ## (time+reserve required to allocate GPU memory)
  ## NOTE: this is not a general solution since jobs may terminate, making a 
  ## possibility for jobs to start simultaniously
  sleep \$((res*6))
fi

## same functional interface:
eval '${MYFUNCS}'

bfail=""

cmd="echo Running... >>\"${statusfile}\""
eval "\${cmd}"

if [ -z "\${bfail}" ]; then
  cmd="cp \"${inpfile}\" \"${orgfile}\""
  echo -e "\${cmd}" >>${errfile}; eval "\${cmd}"

  [ \$? -eq 0 ] || bfail+="ERROR: ${basename}: System error: Copy failed."
fi

#if [ -z "\${bfail}" ]; then
#  cmd="perl -e 'while(<>){s/\\\\x0d/\\\\n/g; print}' ${orgfile} >${inpfile}"
#  echo -e "\${cmd}\n" >>${errfile}; eval "\${cmd}"
#
#  [ \$? -eq 0 ] || bfail+="ERROR: ${basename}: System error: Removing CRs failed."
#fi

if [ -z "\${bfail}" ]; then
  cmd="$gtalignpl --in \"${inpfile}\" --opt \"${optfile}\" --status \"${statusfile}\""
  cmd+=" --logfile \"${gtalignlogfile}\" --reslst \"${reslstfile}\""
  cmd+=" --results \"${resultsfile}\" --err \"${errorfile}\" 2>>${errfile}"
  echo -e "\${cmd}\n" >>${errfile}; eval "\${cmd}"

  [ \$? -eq 0 ] || bfail+="ERROR: ${basename}: Command failed."
fi

#if [ -z "\${bfail}" ]; then
#  if [[ -f "${outfile}" ]]; then
#    cmd="$getids --in ${outfile} --opt ${optfile} 2>>${errfile}"
#    echo -e "\${cmd}\n" >>${errfile}; eval "\${cmd}"
#
#    [ \$? -eq 0 ] || bfail+="ERROR: ${basename}: Command failed."
#  fi
#fi


if [ -n "\${bfail}" ]; then
  echo -e "\${bfail}\n" >>${errfile}
  if [[ ! -f "${errorfile}" || -z "\$(tail -n3 '${errorfile}'|grep -e '^1$')" ]]; then
      echo -e "ERROR: The server's backend issue: System error.\n\n1" >>${errorfile}
  fi
  SendEmail "${sendemail}" "${errfile}" "${UPATH} ${ID}\n\n\${bfail}"
  exit 1
fi


echo
echo "Finished: \$(date)"
exit 0
EOF

if [ $? -ne 0 ]; then
  fail+="ERROR: $ID: System error: Failed to submit a job to a queue.\n"
  echo -e "$fail" >>${errfile}
  echo -e "ERROR: The server's backend issue: Job submission failed.\n\n1" >>${errorfile}
  SendEmail "${sendemail}" "${errfile}" "${UPATH} ${ID}\n\n${fail}"
  exit 1
fi

cmd="echo Queued. >>\"${statusfile}\""
eval "${cmd}"


## print job id
if [ -f "${submitoutfile}" ]; then
  perl -e 'print "$1\n" if <>=~/\s*(\d+)\s*/' "${submitoutfile}"
fi


exit 0

