#!/usr/bin/env Rscript
# This program requires 'R' to run. See https://www.r-project.org/ for installation details.
# For help or examples, run the following --> Rscript rollup.R --help

# The following five R libraries must also be installed:
#   * 1. tibble
#   * 2. readr    <-- These libraries are all available from
#   * 3. dplyr        <Bioconductor>: they are NOT on <CRAN>.
#   * 4. edgeR        See: https://www.bioconductor.org/
#   * 5. Biobase

options(stringsAsFactors=F, menu.graphics=F)
library(optparse)

option_list <- list(
	optparse::make_option(c("-m", "--metadata"), dest="metadata_tsv", type="character",
		help="Input metadata .TSV file with two columns, 'quant_sf_path' and 'sample'. See '--help' for details."),
	optparse::make_option(c("-t", "--threads")  , dest="num_threads", type="integer", default=2,
		help="Num processor threads to use. Default is 2."),
	optparse::make_option(c("-a", "--assembly") , dest="assembly", type="character", default='hg38',
		help="Genome assembly version. Must be 'hg38' for now."),
	optparse::make_option(c("-d", "--datadir")  , dest="datadir", type="character", default='rollup_annotation',
		help="Directory with reference transcriptome annotation data files. This data should have been included with your 'rollup' download."),
	optparse::make_option(c("-o", "--outdir"), dest="outdir", type="character", default='rollup_output_dir',
		help="The output directory for Rollup results, which will be created by this program."),
	optparse::make_option(c("-p", "--priorcount"),dest="priorcount", type="integer", default=5,
		help="Prior count used in edgeR::cpm(...) calculation. Default is 5."),
	optparse::make_option(c("-v", "--verbose")  , dest="verbose", action="store_true", default=FALSE,
		help="Prints verbose messages. Useful for debugging."),
	optparse::make_option(c("-z", "--nozero"), dest="remove_zero", action="store_true", default=FALSE,
		help="Removes any elements (rows) with zero counts. Note: this causes .RDS output files to have differing numbers of rows."))

opt<- optparse::parse_args(optparse::OptionParser(
	option_list=option_list,
	usage=paste("usage: Rscript rollup.R [options]",
		"Rollup takes 1-5 minutes per quant.sf file on a 2018 laptop for a quant.sf with 5+ million lines.",
		"USAGE EXAMPLE:",
		"     rollup.R --metadata=./METADATA.tsv    --datadir=./path/to/annotation_dir/  --assembly=hg38 --outdir=Rolled_up_out_dir",
		" or:",
		"     rollup.R --metadata=./METADATA.tsv    --datadir=./path/to/annotation_dir/  --priorcount=1 --nozero --threads=8 --assembly=hg38 --outdir=Rolled_up_out_dir",
		"",
		" The METADATA.tsv (tab-delimited) file needs two columns: 'sample' and 'quant_sf_path', as in this example:",
		"       sample         quant_sf_path                  <-- These column names are required. Capitalization is irrelevant.",
		"       SOME_ID_1      /full/path/to/ID1/quant.sf.gz",
		"       SOME_ID_2      /full/path/to/ID2/quant.sf.gz  <-- Rollup can handle gzipped files as well.",
		"       SOME_ID_3      /full/path/to/ID3/quant.sf         Gzipping your quant.sf files is recommended.",
		" ",
		" ABOUT THE OUTPUT FILES:",
		"     * 'genenorm' files: CPM was calculated using the gene-level counts (lib.sizes).",
		"         i.e., samples are normalized by gene-level counts, not the repetitive-element-level (RE-level) counts.",
		"         We make the underlying assumption that gene-level counts are more stable than RE-level counts.",
		"     * 'rle_normed_eset' have the same data as above, but are presented as ExpressionSet (R) objects.",
		" ", sep="\n"),
	description="RollUp takes one or more 'quant.sf' files from Salmon (each of which contains transcript-level counts for a single sample) and aggregates the results by repetitive element category ('repName').") 
)


library(tibble);  library(readr); library(dplyr) # <-- install via Bionconductor (part of the "Tidyverse")
library(Biobase); library(edgeR); # <-- install via Bioconductor
library(parallel); # <-- built-in to R

RAW_COUNTS_RDS_SUFFIX  = "_1_raw_counts.RDS"
NORM_COUNTS_RDS_SUFFIX = "_2_counts_normalized.RDS" # normalized by gene
TPM_RDS_SUFFIX         = "_3_TPM.RDS"

verboseMessage <- function(...) { if (opt$verbose) { message(...) } }

# similar to stopifnot
mandate <- function(condition, msg) {
	if (!is.null(condition) && !condition) { warning(msg); stop(msg); } 
}

#' read_rds_or_table
#'
#' Reads a file that can be in either .RDS (R data) format or tab-delimited format (.tsv).
#' @param path: File path to read from.
#' @param required_cols: Optional. For error checking. Columns to *require* in the input file. Exits with an error if these columns are not present.
#' @return Returns the data that from the RDS (if an RDS is read), or a tibble (if a TSV is read).
read_rds_or_table <- function(path, required_cols=c(), ...) {
	# Filename ends in RDS, so try to 'readRDS'
	if (grepl("[.]rds$", path, ignore.case=T, perl=T)) {
		dat = readRDS(path, ...)
	} else {
		dat = readr::read_tsv(path, ...) # note: we expect a TAB-DELIMITED file
	}
	stopifnot(all(required_cols %in% colnames(dat)))
	return(dat)
}

#' load_salmon_quant_sf
#'
#' Generates a three-column tibble from a 'quant.sf' file (which is the output of the aligner 'Salmon'). Note: this function is 'memoized'--it will only read the same filename from disk a single time. It then saves the results into GLOBAL_MEMOIZED_QUANT_SF_DATA.
#' @param filename: The salmon quant.sf file. Also can be quant.sf.gz (or any other compressed format that 'readr' can handle).
#' @return A tibble with the 'Name', 'TPM', and 'EffectiveLength' columns from the specified quant.sf. The "Length" and "EffectiveLength" columns are omitted.
load_salmon_quant_sf <- function(filename) { # Can read a quant.sf file or a quant.sf.gz file.
	mandate(is.character(filename) && grepl("[.]sf", filename, perl=T, ignore.case=T),
		paste0(
			"Warning: 'quant.sf' file did not have '.sf' anywhere in the name. The filename in question was --> ",
			filename)
	)
	if (filename %in% names(GLOBAL_MEMOIZED_QUANT_SF_DATA)) {
		qsf.tib <- GLOBAL_MEMOIZED_QUANT_SF_DATA[[filename]]
	} else {
		qsf.tib <- readr::read_tsv(filename, col_types="cdddd")
		# col_types="cdddd" means: "character", "double"x4 (Name, Length, Effectivelength, TPM, NumReads)
		mandate(ncol(qsf.tib)==5, 
			paste0("Error: input quant.sf file did not have the expected 5 columns. Instead it had the following column names: ",
			paste0(colnames(qsf.tib), collapse=", ")))
		mandate(all(colnames(qsf.tib)[4:5] == c("TPM", "NumReads")),
			paste0("Input quant.sf file lacked a 'TPM' column 4 and 'NumReads' column 5. Those exact column headers are required for a valid salmon 'quant.sf' file. Instead, we found these columns: ",
			paste0(colnames(qsf.tib), collapse=", ")))

		# get rid of Length and EffectiveLength
		qsf.tib <- qsf.tib %>% select(-matches("^(Length|EffectiveLength)$"))

		# "memoize" (save this quant.sf so we don't have to read it from disk again)
		GLOBAL_MEMOIZED_QUANT_SF_DATA[[filename]] <<- qsf.tib
	}
	stopifnot(colnames(qsf.tib) == c("Name","TPM","NumReads"))
	return(qsf.tib) # <-- 3-column tibble (Name, TPM, NumReads). We discard 'Length' and 'EffectiveLength'.
}

#' rollup_and_return
#'
#' This is the main component of RollUp: it runs (in parallel, using 'mclapply') Rollup on data loaded from various quant.sf files, and saves the output to an RDS.
#' @param filepaths.vec: Input character vector of quant.sf file paths. Each element should be a valid file path.
#' @param outdir: Location to save the output RDS file.
#' @param sample_ids.vec: Human-readable names for each quant.sf file. A character vector. Should be the same length as 'filepaths.vec'.
#' @param columns: A vector of characters. Each column name that will be used to roll up the individual elements. If this was set to (for example) repName, then the results would be aggregated at the 'repName' level. Example: if the input data were data.table(STATE=c('Arizona', 'Arizona', 'Colorado'), COUNTS=c(10, 10, 30)), then rolling up with columns='STATE' would result in: Arizona COUNTS=20 and Colorado COUNTS=30.
#' @param outname: Used as a prefix when generating an output filename.
#' @param dat.annot: An annotation object, usually a tibble or data frame, with information such as the repetitive element category for for each ID in each quant.sf file.
#' @return (none, but writes an RDS as a side effect)
rollup_and_return <- function(filepaths.vec, outdir, sample_ids.vec, columns, outname, dat.annot) {
	stopifnot(is.character(filepaths.vec))
	stopifnot(is.character(sample_ids.vec))
	stopifnot(length(sample_ids.vec) == length(filepaths.vec))
	summarize_single_quant_sf <- function(quant_sf_filepath, dat.annot, columns) {
		quant       <- load_salmon_quant_sf(quant_sf_filepath) #%>% rename("md5"="Name") # "Name" column --> "md5"
		merged.top  <- dat.annot %>% left_join(quant, by=c("md5"="Name")) # the md5 value is in the "Name" column in the quant.sf files
		dat.summary <- list()
		for (single_column_name in columns) {
			stopifnot(single_column_name %in% colnames(dat.annot))
			m.tib = merged.top %>% group_by_(.dots=single_column_name) %>% summarize(NumReads=sum(NumReads, na.rm=T), TPM=sum(TPM, na.rm=T))
			m.mat = m.tib %>% select(-1) %>% base::as.matrix() %>% `rownames<-`( m.tib %>% pull(1) )
			dat.summary[[single_column_name]] <- m.mat
		}
		return(dat.summary) # data for exactly one single quant.sf file
	}

	qsf_list = parallel::mclapply(X=filepaths.vec,
		FUN=function(fname) { summarize_single_quant_sf(fname, dat.annot, columns=columns) },
		mc.cores=opt$num_threads) %>% `names<-`( sample_ids.vec )

	return(qsf_list) # Combined version of the single quant.sf files
}

#' make_gene_DGEList
#'
#' Saves a DGEList to disk as a .RDS file, when given a path to a rolled-up RDS file. This function is for genes. Repetitive elements use a very similar function (below) named 'make_repetitive_element_DGEList'.
#' @param salmon_rollup_rds_path: The INPUT path to an RDS file to read.
#' @param metadata_table: Metadata associated with each sample, in data.table or tibble format. Not a filename.
#' @param fileprefix: Prefix for the output RDS filename.
#' @param outdir: The directory where we will write output files.
#' @return (The 'raw counts' and 'TPM' DGELists, in a list, and writes two RDS files as a side effect)
make_gene_DGEList <- function(salmon_data, metadata_table, outdir) {
	# DGEList of read counts mapped to genes (exonic and intronic). Not noramlized.
	# Part 1: obtain 'NumReads' column from the Salmon data.
	res                  = lapply(salmon_data, function(quant) { return(quant[[1]][, "NumReads", drop=F]) }) # [!]
	all_counts           = do.call(cbind, res)
	colnames(all_counts) = metadata_table[[METADATA_SAMPLE_COLNAME]]
	dgelist_counts       = edgeR::DGEList(counts=all_counts, samples=metadata_table, remove.zeros=opt$remove_zero)
	saveRDS(dgelist_counts, file=file.path(outdir, paste0("GENE", RAW_COUNTS_RDS_SUFFIX)))

	# Part 2: obtain 'TPM' column from the Salmon data.
	res_tpm           = lapply(salmon_data, function(quant) { return(quant[[1]][, "TPM", drop=F]) }) # [!]
	all_tpm           = do.call(cbind, res_tpm)
	colnames(all_tpm) = metadata_table[[METADATA_SAMPLE_COLNAME]]
	dgelist_tpm       = edgeR::DGEList(counts=all_tpm, samples=metadata_table, remove.zeros=opt$remove_zero)
	saveRDS(dgelist_tpm, file=file.path(outdir, paste0("GENE", TPM_RDS_SUFFIX)))

	return(list("DGE_RAW"=dgelist_counts, "DGE_TPM"=dgelist_tpm))
}

#' make_repetitive_element_DGEList
#'
#' Saves a DGEList to disk as a .RDS file, when given a path to a rolled-up RDS file. Used for REPETITIVE ELEMENT data only, not genes or introns.
#' @param salmon_rollup_rds_path: The INPUT path to an RDS file to read.
#' @param metadata_table: Metadata associated with each sample, in data.table or tibble format. Not a filename.
#' @param repeat_names_map: A data frame or tibble with three columns: repName, repFamily, and repClass. For human (hg38), we expect 15432 repNames in this file, plus one header line.
#' @param dge_fileprefix: Prefix for the output RDS filename.
#' @param outdir: The directory where we will write output files
#' @return (None, but writes two RDS files to disk.)
make_repetitive_element_DGEList <- function(salmon_data, metadata_table, repeat_names_map, dge_fileprefix, outdir) {
	# DGEList of read counts mapped to repeats, no normalization
	# Read in salmon quant.sf RE list. Length = number of samples, then each line is a list of length 3 (repName, repFamily, repClass)
	res.list <- lapply(salmon_data, function(quant) { return(list("repNameNumReads"=quant[["repName"]][, "NumReads", drop=F], "repNameTPM"     =quant[["repName"]][,      "TPM", drop=F] ))  })
	extract_data <- function(column_name) {
		x <- do.call(cbind, lapply(res.list, '[[', column_name)) %>% `colnames<-`( metadata_table[[METADATA_SAMPLE_COLNAME]] )
		return(x)
	}
	repName_counts   = extract_data("repNameNumReads")
	repName_TPM      = extract_data("repNameTPM")
	stopifnot(identical(rownames(repName_counts), rownames(repName_TPM)))
	repNames_to_keep.bool.vec <- (rownames(repName_counts) %in% repeat_names_map$repName)
	repName_counts_keep       <- repName_counts[repNames_to_keep.bool.vec, , drop=F]
	repName_TPM_keep          <-    repName_TPM[repNames_to_keep.bool.vec, , drop=F]
	# repNames.fdat appears to be repName, repFamily, repClass for each entry
	repNames.fdat             <- data.frame(
		repeat_names_map[match(rownames(repName_counts_keep), repeat_names_map$repName), ])
	stopifnot(colnames(repNames.fdat) == c("repName", "repFamily", "repClass"))
	stopifnot(all.equal(repNames.fdat$repName, rownames(repName_counts_keep)))
	stopifnot(all.equal(repNames.fdat$repName, rownames(repName_TPM_keep)))
	dgelist_repName_counts = edgeR::DGEList(
		counts=repName_counts_keep, samples=metadata_table, genes=repNames.fdat, remove.zeros=opt$remove_zero)
	saveRDS(dgelist_repName_counts, file=file.path(outdir, paste0(dge_fileprefix, RAW_COUNTS_RDS_SUFFIX)))
	dgelist_repName_TPM = DGEList(
		counts=repName_TPM_keep, samples=metadata_table, genes=repNames.fdat, remove.zeros=opt$remove_zero)
	saveRDS(dgelist_repName_TPM, file=file.path(outdir, paste0(dge_fileprefix, TPM_RDS_SUFFIX)))
}

#' make_repName_NormCountDGE_CPM
#'
#' @param srcdir: Directory to read RDS file data from.
#' @param destdir: Directory to save a processed ExpressionSet and DGEList file to, in RDS (R) format.
#' @param input_sample_annotated_data_frame: DGEList 'sample' data from an AnnotatedDataFrame.
#' @param rds_fileprefix: Prefix for the output RDS files.
#' @return (None, but saves two RDS files to disk.)
make_repName_NormCountDGE_CPM <- function(srcdir, destdir, input_sample_annotated_data_frame, rds_fileprefix) {
	# Calculates CPM for repetitive elements. Note: the denominator is NOT THE TOTAL of *RE* counts---it's the total of GENE counts instead.
	# This is because we replace the dge$samples for the REs with the data from the genes.
	# Here, we scan for specific for hard-coded filepaths with a certain format, e.g. 'RE_intergenic_repName_logCPM_rle_normed_eset.RDS'
	dge <- readRDS(file.path(srcdir, paste0(rds_fileprefix, RAW_COUNTS_RDS_SUFFIX)))
	# replace with the norm.factor and lib.size from gene DGE
	# <-- replaces the RE's lib.sizes with the GENE LEVEL counts.
	#		This is how we normalize the REs to gene-level totals instead of RE-level totals.
	dge$samples       = as(input_sample_annotated_data_frame, 'data.frame')
	cpm_matrix        = edgeR::cpm(dge, log=T, prior.count=opt$priorcount)
	cpm_eset          = ExpressionSet(
		assayData=cpm_matrix, phenoData=input_sample_annotated_data_frame, featureData=AnnotatedDataFrame(dge$genes))
	# <-- note the text 'genenorm' here, indicating that these REs are normalized by *gene* counts
	saveRDS(dge, file=file.path(destdir, paste0(rds_fileprefix, NORM_COUNTS_RDS_SUFFIX)))
}

#' assert_metadata_tibble_is_ok
#'
#' Checking for a few common ways that a Metadata file can be malformed. Exits if the metadata that was loaded in does not appear to be OK.
#' @param meta_tib: a tibble that was read in from the metadata file that specifies all the sample names (usually column 1) and paths to their corresponding quant.sf files (usually column 2).
#' @return (None, but terminates the program if the metadata file does not look OK
assert_metadata_tibble_is_ok <- function(meta_tib) {
	mandate(METADATA_QUANTSF_COL %in% colnames(meta_tib),
		paste0("ERROR: Input metadata file (", opt$metadata_tsv, ") was malformed; it needs a path in a column with this name: ", METADATA_QUANTSF_COL, ". But the only column names we found were: ",
		paste0(colnames(meta_tib), sep=", ")))
	mandate(METADATA_SAMPLE_COLNAME %in% colnames(meta_tib),
		paste0("ERROR: Input metadata file (", opt$metadata_tsv, ") was malformed; it needs an ID in a column with this name: ", METADATA_SAMPLE_COLNAME, ". But the only column names we found were: ",
		paste0(colnames(meta_tib), sep=", ")))
	mandate(!any(is.na(meta_tib[[METADATA_QUANTSF_COL]])),
		paste0("The quant_sf_path cannot be 'NA', but one of yours is. These must all be specified and be valid files."))
	mandate(all(file.exists(meta_tib[[METADATA_QUANTSF_COL]])),
		paste0("Every single file that is specified in the quant_sf_path must exist, but the following appear to be missing: ",
		paste0(meta_tib[[METADATA_QUANTSF_COL]][!file.exists(meta_tib[[METADATA_QUANTSF_COL]])], collapse=", ")))
}

mandate(!is.null(opt$metadata_tsv) && is.character(opt$metadata_tsv),
	"You must specify a metadata .TSV file with '--metadata=FILENAME'")
mandate(file.exists(opt$metadata_tsv),
	paste0("Could not read a metadata file at the specified path: ", opt$metadata_tsv))
mandate(!is.null(opt$assembly) && (opt$assembly %in% c("hg38")),
	paste0("Invalid assembly. Currently, we only accept 'hg38'. Your specified assembly was: ",
		as.character(opt$assembly)))
if (!dir.exists(opt$outdir)) { dir.create(opt$outdir) }

GLOBAL_MEMOIZED_QUANT_SF_DATA <- list() # Data structure used to 'memoize' (save for later) 'quant.sf' data to avoid multiple reads.
ALL_RE_LEVELS           <- c("repName", "repFamily", "repClass") # Hierarchical categories of repetitive element classifiction, from most specific (repName) to most general (repClass). Must be in this order.
METADATA_SAMPLE_COLNAME <- "sample"        # <-- usually the 1st column in the metadata.tsv. A human-readable sample identifier (e.g. 'TreatmentA').
METADATA_QUANTSF_COL    <- "quant_sf_path" # <-- usually the 2nd column in the metadata.tsv. Should be full paths to each quant.sf file (e.g. /file/path/to/TreatmentA/quant.sf.gz).

# Below: the main code is here, not encapsulated in a function.
metadata <- readr::read_tsv(opt$metadata_tsv, col_names=TRUE, quote="", comment="") %>% rename_all(tolower) # (lower-case column names)
assert_metadata_tibble_is_ok(metadata)
mandate(dir.exists(opt$datadir),
	paste0("Unable to find specified '--datadir' annotation directory: ",
		opt$datadir, ". You should specify a valid --datadir=/DIRECTORY/PATH/ "))

# Hard-coded file paths here
GENE_ANNOT_TABLE             <- read_rds_or_table(file.path(opt$datadir, 'GENCODE.V26.Basic_Gene_Annotation_md5.RDS')  , required_cols=c("md5", "ensembl_ID", "symbol", "geneID")) %>% select(-matches("^n[.]genes$"))
RE_ANNOT_TABLE               <- read_rds_or_table(file.path(opt$datadir, 'rmsk_annotation.RDS')                        , required_cols=c(ALL_RE_LEVELS, "selected_feature")) %>% select(-matches("^idx$"), -matches("^n[.]loci$"))
RE_INTERGENIC_ANNOT          <- RE_ANNOT_TABLE %>% filter(selected_feature == 'intergenic')
RE_INTRON_ANNOT              <- RE_ANNOT_TABLE %>% filter(selected_feature == 'intron')
RE_EXON_ANNOT                <- RE_ANNOT_TABLE %>% filter(selected_feature == 'exon')
GLOBAL_NAME_FAMILY_CLASS_MAP <- readr::read_tsv(file.path(opt$datadir, 'repName_repFamily_repClass_map.tsv'))
stopifnot(all(colnames(GLOBAL_NAME_FAMILY_CLASS_MAP) == ALL_RE_LEVELS)) # this order is also required

message("Rollup will generate output data in the following output directory: ", opt$outdir)

quant_raw_list = list()
quant_raw_list[["RE_all"]]        = rollup_and_return(metadata[[METADATA_QUANTSF_COL]], outdir=opt$outdir, sample_ids.vec=metadata[[METADATA_SAMPLE_COLNAME]], columns=ALL_RE_LEVELS, outname='RE_all'       , dat.annot=RE_ANNOT_TABLE)
quant_raw_list[["RE_exon"]]       = rollup_and_return(metadata[[METADATA_QUANTSF_COL]], outdir=opt$outdir, sample_ids.vec=metadata[[METADATA_SAMPLE_COLNAME]], columns=ALL_RE_LEVELS, outname='RE_exon'      , dat.annot=RE_EXON_ANNOT)
quant_raw_list[["RE_intron"]]     = rollup_and_return(metadata[[METADATA_QUANTSF_COL]], outdir=opt$outdir, sample_ids.vec=metadata[[METADATA_SAMPLE_COLNAME]], columns=ALL_RE_LEVELS, outname='RE_intron'    , dat.annot=RE_INTRON_ANNOT)
quant_raw_list[["RE_intergenic"]] = rollup_and_return(metadata[[METADATA_QUANTSF_COL]], outdir=opt$outdir, sample_ids.vec=metadata[[METADATA_SAMPLE_COLNAME]], columns=ALL_RE_LEVELS, outname='RE_intergenic', dat.annot=RE_INTERGENIC_ANNOT)
quant_raw_list[["GENE"]]          = rollup_and_return(metadata[[METADATA_QUANTSF_COL]], outdir=opt$outdir, sample_ids.vec=metadata[[METADATA_SAMPLE_COLNAME]], columns="ensembl_ID" , outname='GENE'         , dat.annot=GENE_ANNOT_TABLE)
for (re_category in c("RE_all", "RE_exon", "RE_intergenic", "RE_intron")) {
	# Saves two RDS files to disk for each category.
	make_repetitive_element_DGEList(
		salmon_data=quant_raw_list[[re_category]],
		metadata_table=metadata,
		repeat_names_map=GLOBAL_NAME_FAMILY_CLASS_MAP,
		dge_fileprefix=re_category,
		outdir=opt$outdir)
}
gene_dgelists.list = make_gene_DGEList(
	salmon_data=quant_raw_list[["GENE"]],
	metadata_table=metadata,
	outdir=opt$outdir)
# Note: we calculate the normalization factors based on the GENES, not the REs.
dgelist_gene_count_rle = edgeR::calcNormFactors(
	gene_dgelists.list[["DGE_RAW"]], method="RLE")
# Write the output normalized gene-level data
saveRDS(dgelist_gene_count_rle,
	file=file.path(opt$outdir, paste0("GENE", NORM_COUNTS_RDS_SUFFIX)))
# calculate logCPM
dgelist_gene_cpm       = edgeR::cpm(dgelist_gene_count_rle, log=T, prior.count=opt$priorcount)
# Used for normalization
geneDGEsamples         = Biobase::AnnotatedDataFrame(dgelist_gene_count_rle$samples)

for (fname in c("RE_all", "RE_intergenic", "RE_intron", "RE_exon")) {
	# geneDDEsamples is gene-level count data that has been RLE-normalized
	make_repName_NormCountDGE_CPM(
		srcdir=opt$outdir,
		destdir=opt$outdir,
		input_sample_annotated_data_frame=geneDGEsamples,
		rds_fileprefix=fname)
}
message("[DONE]")
