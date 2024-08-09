#!/usr/bin/env bash
#SBATCH --export=NONE   # required when using 'module'

hostname
echo "Slurm job id:${SLURM_JOBID}:"
date

set -e  #       exit if any command fails
set -u  #       Error on usage of unset variables
set -o pipefail
#set -x  #       print expanded command before executing it


function usage(){
	set +x
	echo
	echo "Usage:"
	echo
	echo $0 --threads 4
	echo --out /francislab/data1/working/20200609_costello_RNAseq_spatial/20230421-TEProF2/out
	echo --extension .STAR.hg38.Aligned.out.bam
	echo /francislab/data1/working/20200609_costello_RNAseq_spatial/20200615-STAR_hg38/out/*.STAR.hg38.Aligned.out.bam
	echo
	exit
}


if [ $( basename ${0} ) == "slurm_script" ] ; then

	echo "Running an individual array job"

	SALMON="/francislab/data1/refs/salmon"

	threads=${SLURM_NTASKS:-4}
	echo "threads :${threads}:"
	mem=${SLURM_MEM_PER_NODE:-30000M}
	echo "mem :${mem}:"

	extension="_R1.fastq.gz"
	paired=false
	k=15

	while [ $# -gt 0 ] ; do
		case $1 in
			--array_file)
				shift; array_file=$1; shift;;
			-k)
				shift; k=$1; shift;;
			--paired)
				paired=true; shift;;
			-o|--out)
				shift; OUT=$1; shift;;
			-e|--extension)
				shift; extension=$1; shift;;
			*)
				echo "Unknown param :${1}:"; usage ;;
		esac
	done
	
	if [ -n "$( declare -F module )" ] ; then
		echo "Loading required modules"
		#module load CBI samtools bwa bedtools2 cufflinks star/2.7.7a
	fi
	
	date
	
	mkdir -p ${OUT}

	line_number=${SLURM_ARRAY_TASK_ID:-1}
	echo "Running line_number :${line_number}:"

	#	Use a 1 based index since there is no line 0.

	echo "Using array_file :${array_file}:"

	line=$( sed -n "$line_number"p ${array_file} )
	echo $line

	if [ -z "${line}" ] ; then
		echo "No line at :${line_number}:"
		exit
	fi

	base=$( basename $line ${extension} )
	#bam=${line}
	R1=${line}

	if ${paired} ; then
		#R2=${line/_R1./_R2.}
		R2=${line/1.fast/2.fast}
	fi

	echo
	echo "base : ${base}"
	echo "ext : ${extension}"
	#echo "r1 : $R1"
	#echo "r2 : $R2"
	#echo "bam : $bam"
	echo 

	date=$( date "+%Y%m%d%H%M%S%N" )

	echo "Running"


	#f=${OUT}/${base}.salmon.REdiscoverTE.k15
	#if [ -d $f ] ; then
	#	echo "Directory $f exists. Skipping."

	f=${OUT}/${base}.salmon.REdiscoverTE.k${k}/quant.sf.gz
	if [ -f ${f} ] && [ ! -w ${f} ] ; then
		echo "Write-protected $f exists. Skipping."
	else

		d=$( dirname ${f} )

		cp ${R1} ${TMPDIR}/
		scratch_R1=${TMPDIR}/$( basename ${R1} )

		if ${paired} ; then
			cp ${R2} ${TMPDIR}/
			scratch_R2=${TMPDIR}/$( basename ${R2} )
			fastas="-1 ${scratch_R1} -2 ${scratch_R2}"
		else
			fastas="--unmatedReads ${scratch_R1}"
		fi

		#scratch_out=${TMPDIR}/$( basename ${f} )
		scratch_out=${TMPDIR}/$( basename ${d} )

		index=${SALMON}/REdiscoverTE.k${k}
		cp -r ${index} ${TMPDIR}/
		scratch_index=${TMPDIR}/$( basename ${index} )

		date

		~/.local/bin/salmon.bash quant --seqBias --gcBias --index ${scratch_index} \
			--no-version-check --libType A --validateMappings ${fastas} \
			-o ${scratch_out} --threads ${SLURM_NTASKS}

		date

		#	GOTTA move an existing dir or we'll move this INTO it.
		if [ -d ${d} ] ; then
			date=$( date "+%Y%m%d%H%M%S" --date="$( stat --printf '%z' ${d} )" )
			mv ${d} ${d}.${date}
		fi

		chmod -R +w ${scratch_out}	#	so script can move and delete the contents (not crucial but stops error messages)
		gzip ${scratch_out}/quant.sf

		mv ${scratch_out} ${d}
		chmod a-w ${f}

		/bin/rm -rf ${scratch_R1}
		if ${paired} ; then
			/bin/rm -rf ${scratch_R2}
		fi

	fi


#	+ salmonstatus=0
#	+ chmod a-w /scratch/gwendt/1475973/p323SF10750-v1_S9
#	+ '[' 0 -ne 0 ']'
#	Wed Jul  5 18:31:01 PDT 2023
#	chmod: cannot access ‘/francislab/data1/working/20230628-Costello/20230706-REdiscoverTE/out/p323SF10750-v1_S9.salmon.REdiscoverTE.k15’: No such file or directory




	date

	echo "Runtime : $((SECONDS/3600)) hrs $((SECONDS%3600/60)) mins $((SECONDS%3600%60)) secs"

else

	date=$( date "+%Y%m%d%H%M%S%N" )
	echo "Preparing array job :${date}:"
	array_file=${PWD}/$( basename $0 ).${date}
	array_options="--array_file ${array_file} "
	
	threads=8
	array=""
	time="1-0"

	while [ $# -gt 0 ] ; do
		case $1 in
			--array)
				shift; array=$1; shift;;
			--time)
				shift; time=$1; shift;;
			-t|--threads)
				shift; threads=$1; shift;;
			-o|--out|--outdir|-e|--extension|-k)
				array_options="${array_options} $1 $2"; shift; shift;;
			--paired)
				array_options="${array_options} $1"; shift;;
			-h|--help)
				usage;;
#			-*)
#				array_options="${array_options} $1"; shift;;
			*)
				echo "Unknown param :${1}: Assuming file"; 
				realpath --no-symlinks $1 >> ${array_file}; shift;;
		esac
	done

	#	True if file exists and has a size greater than zero.
	if [ -s ${array_file} ] ; then

		# using M so can be more precise-ish
		mem=$[threads*7500]M
		scratch_size=$[threads*28]G	#	not always necessary

		max=$( cat ${array_file} | wc -l )

		mkdir -p ${PWD}/logs

		set -x  #       print expanded command before executing it

		[ -z "${array}" ] && array="1-${max}"

		#	A time limit of zero requests that no time limit be imposed.  Acceptable time formats include "minutes", "minutes:seconds", 
		#	"hours:minutes:sec‐onds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds".

		array_id=$( sbatch --mail-user=$(tail -1 ~/.forward)  --mail-type=FAIL --array=${array}%1 \
			--parsable --job-name="$(basename $0)" \
			--time=${time} --nodes=1 --ntasks=${threads} --mem=${mem} --gres=scratch:${scratch_size} \
			--output=${PWD}/logs/$(basename $0).$( date "+%Y%m%d%H%M%S%N" )-%A_%a.out.log \
				$( realpath ${0} ) ${array_options} )

		echo "Throttle with ..."
		echo "scontrol update JobId=${array_id} ArrayTaskThrottle=8"

	else

		echo "No files given"
		usage

	fi

fi

