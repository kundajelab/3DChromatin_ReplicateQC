#!/bin/bash

usage(){
cat <<EOF
usage: `basename $0` options
Installs HiCRep, HiC-Spector, QuASAR-Rep and QuASAR-QC.
OPTIONS
   -h               Show this message and exit
   --pathtopython   Path to python. DEFAULT: python
   --pathtor        Path to R. DEFAULT: R
   --rlib           Path to R libraries. DEFAULT=''
   --pathtobedtools Path to bedtools. DEFAULT=bedtools
   --modules        Names of modules to be loaded. Comma-delimited. This can be used on computing clusters with shared installations, and will be loaded as 'module load modulename'. DEFAULT=''

EOF
}

ARGS=`getopt -o "h" -l "pathtopython:,pathtor:,rlib:,pathtobedtools:,modules:" -- "$@"`
#eval set -- "$ARGS"

#DEFAULTS
PATHTOPYTHON="python"
PATHTOR="R"
RLIB=""
PATHTOBEDTOOLS=""
MODULES=""

while [ $# -gt 0 ]; do
    case $1 in
    -h) usage; exit 1;;
    --pathtopython) PATHTOPYTHON=$2; shift 2;;
    --pathtor) PATHTOR=$2; shift 2;;
    --rlib) RLIB=$2; shift 2;;
    --pathtobedtools) PATHTOBEDTOOLS=$2; shift 2;;
    --modules) MODULES=$2; shift 2;;
    *) usage; exit 1;;
    esac          
done

#============================
# install location
#============================
dir_of_script="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
repo_dir=$(dirname ${dir_of_script})

#============================
# install different software
#============================
mkdir -p ${repo_dir}/software

#if any modules should be loaded, load them now
for modulename in $(echo ${MODULES} | sed 's/,/ /g');
do
    module load ${modulename}
done

pythondir=$(dirname ${PATHTOPYTHON})

git clone http://github.com/kundajelab/genomedisco ${repo_dir}/software/genomedisco
rlibtext=""
if [[ ${RLIB} != "" ]];
then
    rlibtext="--rlib ${RLIB}"
fi

bedtoolstext=""
if [[ ${PATHTOBEDTOOLS} != "" ]];
then
    bedtoolstext="--pathtobedtools ${PATHTOBEDTOOLS}"
fi

modulestext=""
if [[ ${MODULES} != "" ]];
then
    modulestext="--modules ${MODULES}"
fi

${repo_dir}/software/genomedisco/install_scripts/install_genomedisco.sh --pathtopython ${PATHTOPYTHON} --pathtor ${PATHTOR} ${rlibtext} ${bedtoolstext} ${modulestext}
${repo_dir}/software/genomedisco/install_scripts/install_others.sh --pathtopython ${PATHTOPYTHON} --pathtor ${PATHTOR} ${rlibtext} ${bedtoolstext} ${modulestext}

#finally make a softlink for the code
ln -s ${repo_dir}/software/genomedisco/reproducibility_analysis/3DChromatin_ReplicateQC.py ${repo_dir}/3DChromatin_ReplicateQC.py
