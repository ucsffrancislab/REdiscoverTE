#	REdiscoverTE

IMPORTANT. THIS IS NOT MY SOFTWARE.

IT IS IMPORTED (and modified) FROM NON-GITHUB SOURCES.

https://www.nature.com/articles/s41467-019-13035-2

http://research-pub.gene.com/REdiscoverTEpaper/


```
module load r
R
```

```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.15")

BiocManager::install(c('tibble','readr','dplyr','Biobase','edgeR','parallel','EDASeq','ggplot2','RColorBrewer','pheatmap','gridExtra','grid','gtable','RColorBrewer','biomaRt'))

```


#	Create Salmon Index


```BASH
zcat original/REdiscoverTE/rollup_annotation/REdiscoverTE_whole_transcriptome_hg38-20/*.fa.gz > genome.fasta

salmon index \
	-t genome.fasta \
	--threads 64 \
	-i REdiscoverTE
```

#	Align Samples to Salmon Index

```BASH
for f in ${SAMPLE_DIR}/???.fastq.gz ; do

	echo $f

	base=${f%.fastq.gz}
	echo $base

	salmon quant --seqBias --gcBias \
		--index REdiscoverTE \
		--libType A --unmatedReads ${f} \
		--validateMappings \
		-o ${base}.salmon.REdiscoverTE \
		--threads 8

done
```

#	Rollup / Aggregate Alignments to RE repName

```BASH
echo -e "sample\tquant_sf_path" > ${SAMPLE_DIR}/REdiscoverTE.tsv
ls -1 ${SAMPLE_DIR}/*.sample.REdiscoverTE/quant.sf \
		| awk -F/ '{split($8,a,".");print a[1]"\t"$0}' \
		>> ${SAMPLE_DIR}/REdiscoverTE.tsv

REdiscoverTE/rollup.R \
		--metadata=${SAMPLE_DIR}/REdiscoverTE.tsv \
		--datadir=REdiscoverTE/rollup_annotation/ \
		--nozero --threads=64 --assembly=hg38 \
		--outdir=${SAMPLE_DIR}/REdiscoverTE_rollup/
```

#	Analyze Results with EdgeR













#	References

This software was initially from https://www.nature.com/articles/s41467-019-13035-2

http://research-pub.gene.com/REdiscoverTEpaper/

http://research-pub.gene.com/REdiscoverTEpaper/data/

http://research-pub.gene.com/REdiscoverTEpaper/data/REdiscoverTEdata_README.html

http://research-pub.gene.com/REdiscoverTEpaper/data/REdiscoverTEdata_1.0.1.tar.gz

http://research-pub.gene.com/REdiscoverTEpaper/software/

http://research-pub.gene.com/REdiscoverTEpaper/software/REdiscoverTE_README.html

http://research-pub.gene.com/REdiscoverTEpaper/software/REdiscoverTE_1.0.1.tar.gz

These files were downloaded and retained in the original/ directory.
They were untarred to minimize file size.


```BASH
gzip original/REdiscoverTE/EXPECTED_OUTPUT_FILES/Step_2_salmon_counts/quant.sf

mkdir -p original/REdiscoverTE/rollup_annotation/REdiscoverTE_whole_transcriptome_hg38-20
faSplit sequence original/REdiscoverTE/rollup_annotation/REdiscoverTE_whole_transcriptome_hg38.fa 20 original/REdiscoverTE/rollup_annotation/REdiscoverTE_whole_transcriptome_hg38-20/
gzip original/REdiscoverTE/rollup_annotation/REdiscoverTE_whole_transcriptome_hg38-20/*.fa
```


