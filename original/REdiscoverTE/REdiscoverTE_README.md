# REdiscoverTE
Author: *Haiyin Chen*

Date: *2019-01*

# Overview
This software enables whole-transcriptome RNA-seq simultaneous quantification of repetitive elements and human transcripts. Full description can be found in the publication `Kong Y, ..., Chen-Harris H (2019)`. It is written in R.

*REdiscoverTE* takes RNA-seq `FASTQ` files as input, and produces both raw and normalized counts for the repetitive elements (REs) in subfamilies defined by [RepeatMasker](http://www.repeatmasker.org/).


# Pipeline overview
   * *REdiscoverTE* transcriptome (FASTA) sequences are provided for `hg38` (human). This file is used by the open source software `salmon` to generate reference index
   * Salmon takes input RNAseq `FASTQ` file -> `quant.sf`
   * The R script `rollup.R` combines the counts at members of a RE into a single value at the subfamily level and produces `.RDS` files that can be loaded into R and exported if desired.

# Installing & Running REdiscoverTE

Follow the steps below to run the *REdiscoverTE* pipeline on the included simulated RNAseq data. 

## Prerequisites
### 1. Required Software
  * Two programs must already be installed: `R` and `salmon`.
  * `R` (https://www.r-project.org/). Tested version is `3.4.3`.
  * `salmon`. Tested version is `0.8.2`, available at the developers' github page: https://github.com/COMBINE-lab/salmon/releases/tag/v0.8.2 .

### 2. Required RAM and disk space
  * `89 gigabytes` of disk space is required for the Salmon index.
  * `30 GB of RAM` is recommended for **Salmon**. Salmon requires a substantial amount of RAM with the included 5-million-entry FASTA reference.

### 3. Time for each step
   1. Index generation: 90 minutes on a 2017 16-core 2.6 GHz Xeon processor.
   2. Salmon alignment: ~30-90 minutes. ~30 minutes for a 64 million sequence 50bp SE input FASTQ file on the same Xeon CPU as above.
   3. Rollup: ~5 minutes per `quant.sf` file. (There will be one `quant.sf` file for each input sample: for single-end reads, there will be a 1:1 correspondence between `.fastq` and `quant.sf` files).

### 4. Running the *REdiscoverTE* pipeline on the included sample data
  * You should be able enter the project directory (`cd /path/to/your/rollup_dir/`) and then run `make all` to run the commands in the included Makefile.
  * This will run the two primary steps of this process: the `salmon` quantification step and the `Rollup` step.
  * You can include multiple FASTQ files at the same time: each sample will generate a single output column in the final results.
  * Both single-end and paired-end FASTQ files are accepted.

### 5. Interpreting the result files
  * The result files (referred to as `rollup` output files) R data files in `.RDS` format.
  * You can load and view (or re-export) them as follows:
    * `x <- readRDS("filename.RDS")                # load the data file`
    * The files are various data types; for example, `DGELists` (from the [EdgeR](https://www.rdocumentation.org/packages/edgeR/versions/3.14.0/topics/DGEList) package) and `ExpressionSets` (from the [BioBase](https://www.rdocumentation.org/packages/Biobase/versions/2.32.0/topics/ExpressionSet) package).

## Included scripts
  1. `rollup.R`
     * described above
  2. `Makefile`
     * a list of example commands for running a *REdiscoverTE* pipeline. Running `make all` will run the entire pipeline on the included sample data, if the prerequisite software is present.

## Included annotation files used by *Rollup*

These four required files are all for `hg38` (human). They are included with this distribution of *REdiscoverTE*.

  1. `rollup_annotation/GENCODE_v26_basic_transcript_rmsk_stdchr_hg38_introns_ov_REs.fa.xz` (human genome FASTA file, compressed with `xz`)
  2. `rollup_annotation/rmsk_annotation_Salmon.RDS` (R data object containing information about repetitive elements, loadable from within R with `readRDS`)
  3. `rollup_annotation/GENCODE.V26.Basic_Gene_Annotation_md5.RDS` (R data object containing information about *GENCODE v. 26 Basic* human transcripts, loadable from within R with `readRDS`)
  4. `rollup_annotation/repName_repFamily_repClass_map.tsv` (tab-delimited text file)

### File #1 of 4: `hg38` genome and associated elements
   * Path: `rollup_annotation/GENCODE_v26_basic_transcript_rmsk_stdchr_hg38_introns_ov_REs.fa.xz`
   * This is a standard FASTA file. The record names (e.g. `9bd267abb08d45904c1741db501b5bdd`) are MD5 checksums of each sequence.
   * This is the file used to build the 89 gigabyte Salmon index.
   * This FASTA file consists of three different types of sequences:
      1. `hg38` transcripts from GENCODE v. 26 (basic).
      2. Any introns from GENCODE v. 26 (basic) that *also* contained one or more repetitive elements.
      3. All repetitive elements from RepeatMasker.
   * This file is compressed with the `xz` utility, which is similar to `gzip` or `zip`. If desired, it can be uncompressed with the standard `unxz` utility. It can be viewed on the command line with the following command: `xzcat ./rollup_annotation/GENCODE_v26_basic_transcript_rmsk_stdchr_hg38_introns_ov_REs.fa.xz | less -S`

### Viewing `.RDS` files
The `.RDS` files can be viewed as follows:
  1. Begin a new `R` session.
  2. Within R: `x <- readRDS("PATH/TO/FILENAME")`
  3. Within R: `head(x) # view the top 6 rows`

### File #2 of 4: RepeatMasker annotation (hg38)
   * Path: `rollup_annotation/rmsk_annotation_Salmon.RDS`
   * This file contains the correspondences between the 32-character MD5s (see above) to each RMSK (RepeatMasker) sequence. It also contains RMSK annotations and *genomic context* (*intergenic*, *intron*, or *exon*) for each MD5.
   * Annotation is from `GENCODE v26 *Basic*`. snRNA sequences and all duplicate sequences were removed.
   * Number of unique MD5s: `5099056`
   * Number of unique RepeatMakser *repNames* (subfamilies): `17468`
   * Example of what this repetitive element annotation file looks like after being loaded into `R`:
```
md5                              n.loci   repName      repClass     repFamily selected_feature            idx
9bd267abb08d45904c1741db501b5bdd      1 (TAACCC)n Simple_repeat Simple_repeat       intergenic 1__10001_10468
e98ef316f41469ed5d18179d6f4c246f      1      TAR1     Satellite          telo       intergenic 1__10469_11447
72834a9eaede5c9cad25de2588a65486      1      MIR3          SINE           MIR           intron 1__15265_15355
accd2fc2b8c6fb63f05975ebeaf89439      1        L3          LINE           CR1           intron 1__19972_20405
...
```

### File #3 of 4: Human transcript annotation
   * Path: `rollup_annotation/GENCODE.V26.Basic_Gene_Annotation_snRNA_removed_md5.RDS`
   * Transcript (and gene) annotation for the hg38 human genome assembly
   * Mapping between MD5 of transcripts from `GENCODE v.26 BASIC`.
   * `98,029` *unique* ENSEMBL transcripts IDs (`ENST` IDs) in total
   * `57,374` *unique* ENSEMBL gene IDs (`ENSG` IDs).
   * `55,842` *unique* 'symbol' (gene symbols)
   * Example of what this transcript annotation file looks like:
```
md5                        n.genes transcript_ID    ensembl_ID      symbol     geneID
0001d490f52254a5ada6d26...       1 ENST000002630... ENSG00000083... ZNF264     9422
0002168cd59867bac9fbab4...       1 ENST000004428... ENSG00000109... CLNK       116449
0002d0449e073ab67c8a800...       1 ENST000003175... ENSG00000163... ZMYM6      9204
000400e08761c20fdf45ece...       1 ENST000006355... ENSG00000236... LINC011... NA
...                            ... ...              ...             ...        ...
```

### File #4 of 4: RE category (taxonomic hierarchy) mapping
   * Path: `rollup_annotation/repName_repFamily_repClass_map.tsv`
   * Three columns, from most-specific to least-specific annotation category: `repName`, `repFamily`, and `repClass`. The term `repName` is sometimes also referred to as `subfamily`.
   * Example of what this file looks like:
```
 repName    repFamily     repClass
 5S         rRNA          rRNA
 7SK        RNA           RNA
 7SLRNA     srpRNA        srpRNA
 (AAAAAAC)n Simple_repeat Simple_repeat
 ...        ...           ...
```
   * Note that this file differs slightly from `rmsk_annotation_Salmon.RDS`: the RDS file contains multiple categories joined with commas (e.g. `(AAAAACA)n,(CAAAAA)n`), which are not included in this mapping file.

## Output files in `.RDS` R data format.

To read these files into R, from within an R session:
   * `library(edgeR);`
   * `x <- readRDS("the_filename.RDS")`
   * `x` will be a `DGEList` data type, which is part of the `edgeR` package. See: https://www.rdocumentation.org/packages/edgeR/versions/3.14.0/topics/DGEList

   * `GENE`: Transcript-level expression values:
      * `GENE_1_raw_counts.RDS`: raw counts DGEList data object.
      * `GENE_2_counts_normalized.RDS`: normalized counts (by total gene counts) DGEList data object.
      * `GENE_3_TPM.RDS`: TPM (transcripts per million) for each gene/transcript. DGEList data object.

   * `RE all`: Repetitive elements found anywhere in the genome. Granularity for all RE reporting is at `repName` level (where the three levels we use are `name`, `family`, and `class`, in increasing order of generality).
      * `RE_all_1_raw_counts.RDS`: raw counts DGEList data object.
      * `RE_all_2_counts_normalized.RDS`: note that counts are normalized by total gene counts
      * `RE_all_3_TPM.RDS`: TPM (transcripts per million) for each repetitive element 'repName'. DGEList data object.

   * `RE exon`: Subset of `RE all`: only repetitive elements found at least partially within an annotated exon.
      * `RE_exon_1_raw_counts.RDS`: raw counts DGEList data object.
      * `RE_exon_2_counts_normalized.RDS`: counts are normalized by total gene counts
      * `RE_exon_3_TPM.RDS`: TPM (transcripts per million) for each repetitive element 'repName'. DGEList data object.

   * `RE intron`: Subset of `RE all`: only repetitive elements that do not have any overlap with an exon, and have some overlap with an intron.
      * `RE_intron_1_raw_counts.RDS`: raw counts DGEList data object.
      * `RE_intron_2_counts_normalized.RDS`:counts are normalized by total gene counts
      * `RE_intron_3_TPM.RDS`: TPM (transcripts per million) for each repetitive element 'repName'. DGEList data object.
      
   * `RE intergenic`: Subset of `RE all`: only repetitive elements that have no overlap with annotated introns or exons.
      * `RE_intergenic_1_raw_counts.RDS`: raw counts DGEList data object.
      * `RE_intergenic_2_counts_normalized.RDS`: counts are normalized by total gene counts
      * `RE_intergenic_3_TPM.RDS`: TPM (transcripts per million) for each repetitive element 'repName'. DGEList data object.



