#	REdiscoverTE



#	Create Salmon Index



#	Align Samples to Salmon Index




#	Rollup / Aggregate Alignments to RE repName





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
gzip original/REdiscoverTE/rollup_annotation/REdiscoverTE_whole_transcriptome_hg38-20/?.fa
```


