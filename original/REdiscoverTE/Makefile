SHELL     := /bin/bash

# These two executables: Rscript (R) and Salmon must already be installed on your system
RSCRIPT_EXE=Rscript
SALMON_0.8.2_EXE=salmon
ROLLUP_EXE=./rollup.R
ROLLUP_ANNOTATION_DATA_DIR=rollup_annotation
XZ_INPUT_GENOME=$(ROLLUP_ANNOTATION_DATA_DIR)/REdiscoverTE_whole_transcriptome_hg38.fa.xz

#FASTQ_PAIRED_TYPE=SINGLE
FASTQ_PAIRED_TYPE=PAIRED
FASTQ_READS_1=SIMULATED_FASTQS/input_R1.fq.gz
FASTQ_READS_2=SIMULATED_FASTQS/input_R2.fq.gz

SALMON_INDEX_DIR=Step_1_salmon_idx
SALMON_INDEX_STDERR_WARNINGS=Step_1b_salmon_index_construction_warnings.txt
SALMON_INDEX_HASH_BIN=$(SALMON_INDEX_DIR)/hash.bin

SALMON_COUNTS_DIR=Step_2_salmon_counts
SALMON_QUANT_SF=$(SALMON_COUNTS_DIR)/quant.sf

METADATA_FOR_ROLLUP=Step_3_metadata_for_rollup.tsv.txt
ROLLUP_RESULTS_DIR=Step_4_rollup
ROLLUP_RESULTS_TIME_FILE=$(ROLLUP_RESULTS_DIR)/done.time.txt

generated=$(SALMON_COUNTS_DIR) $(SALMON_QUANT_SF) $(ROLLUP_RESULTS_DIR) $(ROLLUP_RESULTS_TIME_FILE) $(SALMON_INDEX_DIR) $(SALMON_INDEX_HASH_BIN) genome.fasta $(METADATA_FOR_ROLLUP)
extras=$(SALMON_INDEX_STDERR_WARNINGS)

.PHONY: clean

all: $(generated)

clean:
	@echo "Deleting the Makefile-generated files, including the very-slow-to-regenerate Salmon index."
	${RM} $(generated)
	${RM} $(extras)

$(XZ_INPUT_GENOME):
	@echo "This file is included as with the base 'REdiscoverTE' download, and should already exist."

genome.fasta: $(XZ_INPUT_GENOME)
	@echo "Extracting the 'xz'-compressed input genome..."
	xz --keep --decompress --stdout --threads=0 $(XZ_INPUT_GENOME) > genome.fasta

$(SALMON_INDEX_HASH_BIN) $(SALMON_INDEX_DIR): genome.fasta
	@echo "Builds a Salmon index. Tested with Salmon version 0.8.2."
	@echo "Warning: Please note that you will need 90 gigabytes of free space for the generated index."
	@echo "This specific index is for ~90k human hg38 transcripts from:"
	@echo "                     1. Gencode BASIC (human transcripts)"
	@echo "                     2. ~5 million repetitive elements (individual loci)"
	@echo "                     3. Any introns containing 1+ REs (a large fraction of total introns)."
	@echo "Index building takes several hours to run and generates a 90 gigabyte output file."
	@echo "If you do not have at least 100+ GB of free disk space on the machine that is"
	@echo "generating this index, you should cancel this operation and free up more disk space."
	$(SALMON_0.8.2_EXE) index -t genome.fasta -i $(SALMON_INDEX_DIR) --type 'quasi' -s 1    2> $(SALMON_INDEX_STDERR_WARNINGS)

FASTQ_INPUT_ARGS="(Not set yet: this will be set after checking the 'FASTQ_PAIRED_TYPE' variable)"
ifeq ($(FASTQ_PAIRED_TYPE),PAIRED)
	FASTQ_INPUT_ARGS=-1 $(FASTQ_READS_1) -2 $(FASTQ_READS_2)
else ifeq ($(FASTQ_PAIRED_TYPE),SINGLE)
	FASTQ_INPUT_ARGS=-r $(FASTQ_READS_1)
else
	$(error The 'FASTQ_INPUT_ARGS' argument must be either PAIRED or SINGLE, capitalized.)
endif
$(SALMON_QUANT_SF) $(SALMON_COUNTS_DIR): $(SALMON_INDEX_DIR) $(FASTQ_READS_1)
	@echo "Running Salmon on input FASTQ data." # module load apps/salmon/0.8.2
	@echo "This code has been teste specifically with Salmon 0.8.2."
	@echo "Later versions of Salmon may also work, but have not been specifically verified."
	@echo "Additional Salmon documentation: https://salmon.readthedocs.io/en/latest/salmon.html "
	$(SALMON_0.8.2_EXE)  quant  --seqBias  --gcBias  -i $(SALMON_INDEX_DIR)  -l 'A'  $(FASTQ_INPUT_ARGS)  -o $(SALMON_COUNTS_DIR)

$(METADATA_FOR_ROLLUP): $(SALMON_QUANT_SF)
	echo -e "sample\tquant_sf_path" > $(METADATA_FOR_ROLLUP)
	echo -e "SalmonFile1\t$(SALMON_QUANT_SF)" >> $(METADATA_FOR_ROLLUP)

$(ROLLUP_RESULTS_TIME_FILE) $(ROLLUP_RESULTS_DIR): $(SALMON_QUANT_SF) $(METADATA_FOR_ROLLUP)
	$(RSCRIPT_EXE) $(ROLLUP_EXE) -m $(METADATA_FOR_ROLLUP) --datadir=$(ROLLUP_ANNOTATION_DATA_DIR) --threads=1 --outdir=$(ROLLUP_RESULTS_DIR) --assembly="hg38"
	date > $(ROLLUP_RESULTS_TIME_FILE)

