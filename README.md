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























#	View TCGA matrix


```
dat <- readRDS("original/REdiscoverTEdata/inst/Fig4_data/eset_TCGA_TE_intergenic_logCPM.RDS")
head(dat)

ExpressionSet (storageMode: lockedEnvironment)
assayData: 6 features, 7353 samples 
  element names: exprs 
protocolData: none
phenoData
  sampleNames: TCGA-OR-A5KX-01A-11R-A29S-07
    TCGA-EJ-7784-01A-11R-2118-07 ... TCGA-CV-7416-01A-11R-2081-07 (7353
    total)
  varLabels: indication patient_ID ... paired (5 total)
  varMetadata: labelDescription
featureData
  featureNames: 7SK 7SLRNA ... AluJb (6 total)
  fvarLabels: repName repClass repFam
  fvarMetadata: labelDescription
experimentData: use 'experimentData(object)'
Annotation:  

> as.data.frame(dat)[1:3,1:3]
                                  X7SK    X7SLRNA      ACRO1
TCGA-OR-A5KX-01A-11R-A29S-07 -1.386203 -0.7561825  2.2834249
TCGA-EJ-7784-01A-11R-2118-07 -2.439085 -0.2542860 -0.2175945
TCGA-BW-A5NQ-01A-11R-A27V-07 -3.079350 -1.3472844 -1.5813776

#	I don't think that R data frame column names can begin with a number.


> dim(as.data.frame(dat))
[1] 7353 1209

#	7353 subjects
#	~1209 RE/TE counts

> fData(dat)[1:5,]
            repName  repClass repFam
7SK             7SK       RNA    RNA
7SLRNA       7SLRNA    srpRNA srpRNA
ACRO1         ACRO1 Satellite   acro
ALR/Alpha ALR/Alpha Satellite  centr
Alu             Alu      SINE    Alu

#	Perhaps the rows that begin with numbers should have an X added as well?

```



##	20240516

Getting rid of LFS (large files)


Not sure if this file is ever used, but if it is just concat the pieces.
```
cd original/REdiscoverTE/EXPECTED_OUTPUT_FILES/Step_2_salmon_counts/
split -d -b 20MB quant.sf quant.sf.
git add quant.sf.??
```

Also

```
cd original/REdiscoverTE/rollup_annotation
split -d -b 20MB rmsk_annotation.RDS rmsk_annotation.RDS.
git add rmsk_annotation.RDS.??
```







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


