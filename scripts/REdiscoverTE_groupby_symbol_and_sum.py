#!/usr/bin/env python3
import sys
import pandas as pd
df = pd.read_csv(sys.argv[1],sep='\t')
df.groupby(['symbol']).sum().to_csv(sys.stdout,sep='\t')


#   ./groupby_symbol_and_sum.py REdiscoverTE_rollup_noquestion/GENE_1_raw_counts.symbols.tsv > REdiscoverTE_rollup_noquestion/GENE_1_raw_counts.symbols.summed.tsv

