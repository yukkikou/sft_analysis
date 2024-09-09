#!/usr/bin/env python

import pandas as pd
import anndata
import sys
import os
import argparse
import logging
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

# Configure the logging
logging.basicConfig(
    format='%(asctime)s - %(levelname)s - %(message)s',
    level=logging.INFO
)

logger = logging.getLogger()
logger.info("*** Start Analysis ***")

# environment
workdir = '/media/data5/hxy/SFT/RNAseq/4_Spatial/CCC/cellphoneDB/'
os.chdir(workdir)
print("Now working directory:", os.getcwd())

# user files:"A1-ZX-2","A-2401158-04","A-LSC-SM-REC","D1-GHL-3"
parser = argparse.ArgumentParser(description='Process sample prefix for CellPhoneDB analysis.')
parser.add_argument('--sample_prefix', type=str, required=True, help='Prefix for the sample to be analyzed')
args = parser.parse_args()
sample_prefix = args.sample_prefix

# data
cpdb_file_path = 'db/v5/cellphonedb.zip'
counts_file_path = os.path.join(sample_prefix, 'data_plus', 'normalised_log_counts.csv')
counts_matrix_path = os.path.join(sample_prefix, 'data_plus', 'normalised_log_matrix.csv')
meta_file_path = os.path.join(sample_prefix, 'data_plus', 'metadata.tsv')
matching_meta_file_path = os.path.join(sample_prefix, 'data_plus', 'matching_metadata.tsv')
microenvs_file_path = os.path.join(sample_prefix, 'data_plus', 'microenvironment.tsv')
active_tf_path = os.path.join(sample_prefix, 'data', 'active_TFs.tsv')
out_path = os.path.join(sample_prefix, 'results', 'method2_plus')

microenv = pd.read_csv(microenvs_file_path, sep = '\t')
counts_df = pd.read_csv(counts_file_path)
counts_matrix = counts_df.set_index('Gene')
counts_matrix.to_csv(counts_matrix_path)

metadata = pd.read_csv(meta_file_path, sep = '\t')
matching_metadata = metadata.set_index('barcode_sample').loc[counts_matrix.columns].reset_index()
matching_metadata.to_csv(matching_meta_file_path, sep='\t', index=False)

# main
logger.info("*** Starting CellphoneDB ***")

# Run statistical analysis
cpdb_results = cpdb_statistical_analysis_method.call(
    cpdb_file_path = cpdb_file_path,
    meta_file_path = matching_meta_file_path,
    counts_file_path = counts_matrix_path,
    counts_data = 'hgnc_symbol',
    microenvs_file_path = microenvs_file_path,
    score_interactions = True,
    iterations = 1000, # denotes the number of shufflings performed in the analysis.
    threshold = 0.1, # defines the min % of cells expressing a gene for this to be employed in the analysis.
    threads = 25, # number of threads to use in the analysis.
    debug_seed = 42, # debug randome seed. To disable >=0.
    result_precision = 3, # Sets the rounding for the mean values in significan_means.
    pvalue = 0.05, # P-value threshold to employ for significance.
    subsampling = False, # To enable subsampling the data (geometri sketching).
    subsampling_log = False, # (mandatory) enable subsampling log1p for non log-transformed data inputs.
    subsampling_num_pc = 100, # Number of componets to subsample via geometric skectching (dafault: 100).
    subsampling_num_cells = 1000, # Number of cells to subsample (integer) (default: 1/3 of the dataset).
    separator = '|', # "cellA|CellB".
    debug = False, # Saves all intermediate tables employed during the analysis in pkl format.
    output_path = out_path,
    output_suffix = sample_prefix
)

logger.info("*** Finished ***")
