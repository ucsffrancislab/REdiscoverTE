#!/usr/bin/env bash

#BOX_BASE="https://dav.box.com/dav/Francis _Lab_Share"
BOX_BASE="ftps://ftp.box.com/Francis _Lab_Share"

PROJECT=$( basename ${PWD} )
DATA=$( basename $( dirname ${PWD} ) )

#BOX="${BOX_BASE}/${DATA}"
#curl -netrc -X MKCOL "${BOX}/"
#BOX="${BOX_BASE}/${DATA}/${PROJECT}"
#curl -netrc -X MKCOL "${BOX}/"

#curl -netrc -T metadata.csv "${BOX}/"
#BOX="${BOX_BASE}/${DATA}/${PROJECT}/rollup"
#curl -netrc -X MKCOL "${BOX}/"

BOX="${BOX_BASE}/${DATA}/${PROJECT}/rollup"
for f in rollup/REdiscoverTE.tsv rollup/rollup.merged/* ; do
	echo $f
	curl --silent --ftp-create-dirs -netrc -T ${f} "${BOX}/"
done

for d in rmarkdown_results* ; do
	echo $d
	BOX="${BOX_BASE}/${DATA}/${PROJECT}/${d}"
	#curl -netrc -X MKCOL "${BOX}/"
	for f in ${d}/* ; do
		echo $f
		curl  --silent --ftp-create-dirs -netrc -T ${f} "${BOX}/"
	done
done

echo "Runtime : $((SECONDS/3600)) hrs $((SECONDS%3600/60)) mins $((SECONDS%3600%60)) secs"

