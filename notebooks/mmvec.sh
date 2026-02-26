#!/bin/bash -l

#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 24:00:00
#SBATCH --mem=64G


conda activate qiime2-2020.6



# MMVEC -- Shiseido Project -- all microbial features and all metabolite features using all samples 05/23/2025

qiime mmvec paired-omics \
    --i-microbes ../data/mmvec/217045_rarefied_table_RefHit_FULL_mmvec_relative_abundance_subset_samples_2020.6.qza \
    --i-metabolites ../data/mmvec/Metabolomics_data_sample_05022025_FULL_subset_sample_2020.6.qza \
    --p-summary-interval 1 \
    --p-learning-rate 0.001 \
    --p-epochs 300  \
    --output-dir ../output/mmvec/null_summary_inv_scan_FULL_05232025 \
    --verbose

