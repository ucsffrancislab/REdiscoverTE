#!/usr/bin/env bash
#SBATCH --export=NONE   # required when using 'module'

hostname
echo "Slurm job id:${SLURM_JOBID}:"
date

#if [ $# -gt 0 ] ; then
#if [ $( basename ${0} ) == "slurm_script" ] ; then

	set -e	#	exit if any command fails
	set -u	#	Error on usage of unset variables
	set -o pipefail
	if [ -n "$( declare -F module )" ] ; then
		echo "Loading required modules"
		module load CBI r
	fi
	set -x	#	print expanded command before executing it

	#ARGS=$*
	ARGS=""
	#	given that these params are directed at an R script, they are separated by a =
	#echo "Parsing $* before calling rollup.R"
	while [ $# -gt 0 ] ; do
		case $1 in
#			--datadir)
#				shift; datadir=$1; shift;;
#			-i|--in|--indir)
#				shift; INDIR=$1; shift;;
			-o|--out|--outdir|--outbase)
				shift; OUTBASE=$1; shift;;
#			--outbase=*)
#				OUTBASE=$( echo $1 | cut -d= -f2 ); shift;;
			*)
				ARGS="${ARGS} $1"; shift;;
		esac
	done

#	#~/github/ucsffrancislab/REdiscoverTE/rollup.R $*
#	~/github/ucsffrancislab/REdiscoverTE/rollup.R $ARGS



#	ONLY MERGE ONCE, NOT AFTER EVERY SUBSET FINISHES


	if [ $( ls -1d ${OUTBASE}/rollup.* | wc -l ) -gt 1 ] ; then
		REdiscoverTE_rollup_merge.Rscript --outdir=${OUTBASE}
	else
		ln -s $( basename $( ls -1d ${OUTBASE}/rollup.* ) ) ${OUTBASE}/rollup.merged
	fi

	for f in ${OUTBASE}/rollup.merged/* ; do
		ln -s rollup.merged/$( basename ${f} ) ${OUTBASE}/$( basename ${f} )
	done

#else
#
#	INDIR="${PWD}/out"
#	k=15
#
#	date=$( date "+%Y%m%d%H%M%S" )
#
#	OUTBASE="${PWD}/rollup"
#
#	datadir=/francislab/data1/refs/REdiscoverTE/rollup_annotation/
#
#	while [ $# -gt 0 ] ; do
#		case $1 in
#			--datadir)
#				shift; datadir=$1; shift;;
#			-i|--in|--indir)
#				shift; INDIR=$1; shift;;
#			-o|--out|--outdir|--outbase)
#				shift; OUTBASE=$1; shift;;
##			--array)
##				shift; array=$1; shift;;
##			-t|--threads)
##				shift; threads=$1; shift;;
##			-o|--out|--outdir|-e|--extension)
##				array_options="${array_options} $1 $2"; shift; shift;;
##			--paired)
##				array_options="${array_options} $1"; shift;;
##			-h|--help)
##				usage;;
##			-*)
##				array_options="${array_options} $1"; shift;;
##			*)
##				echo "Unknown param :${1}: Assuming file"; 
##				realpath --no-symlinks $1 >> ${array_file}; shift;;
#		esac
#	done
#
#	mkdir -p ${OUTBASE}
#
#
#	ls -1 ${INDIR}/*.salmon.REdiscoverTE.k${k}/quant.sf.gz \
#		| awk -F/ '{split($8,a,".");print a[1]"\t"$0}' > ${OUTBASE}/REdiscoverTE.tsv
#
#	split --lines=1000 --numeric-suffixes --additional-suffix=.tsv ${OUTBASE}/REdiscoverTE.tsv ${OUTBASE}/REdiscoverTE.
#
#	#for f in ${OUTBASE}/REdiscoverTE.??.tsv ; do	#	always 2 digits?
#	for f in ${OUTBASE}/REdiscoverTE.*.tsv ; do
#		echo $f
#		d=$( basename $f .tsv )
#		d=${d#*.}
#
#		echo $d
#
#		OUTDIR="${OUTBASE}/rollup.${d}"
#
#		mkdir -p ${OUTDIR}
#
#		#echo -e "sample\tquant_sf_path" > ${OUTDIR}/REdiscoverTE.tsv
#
#		sed -i '1i sample\tquant_sf_path' $f
#
#		#	A time limit of zero requests that no time limit be imposed.  Acceptable time formats include "minutes", "minutes:seconds", 
#		#	"hours:minutes:sec‐onds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds".
#
#		sbatch --mail-user=$(tail -1 ~/.forward)  --mail-type=FAIL \
#			--job-name="Rollup${d}" \
#			--output="${PWD}/logs/REdiscoverTE.${date}.rollup.${d}.out" \
#			--time=14-0 --nodes=1 --ntasks=64 --mem=495G \
#			${0} \
#			--metadata=${OUTBASE}/REdiscoverTE.${d}.tsv \
#			--datadir=${datadir} \
#			--nozero --threads=64 --assembly=hg38 --outdir=${OUTDIR}/ --outbase=${OUTBASE}
#
#			#${PWD}/REdiscoverTE_rollup.bash \
##	Question or no question? That is the question.
##			--datadir=/francislab/data1/refs/REdiscoverTE/rollup_annotation.noquestion/ \
#
#	done
#
#fi

echo "Runtime : $((SECONDS/3600)) hrs $((SECONDS%3600/60)) mins $((SECONDS%3600%60)) secs"

