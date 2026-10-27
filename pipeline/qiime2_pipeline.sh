#!/bin/bash
## q2pipe.sh -- Combined Script for Rapamycin Data using Qiime2 + DADA2
## Date: [10-12-2024]

# Set script to exit immediately on error and print commands as they execute
set -e  # Exit immediately on error
set -x  # Print commands as they execute

# Check if a configuration file is provided
if [ $# -ne 1 ]; then
    echo "USAGE: $0 configFile.txt"
    exit 1
fi

# Read configuration file
source $1  # The configuration file should define variables like SAMPLEMETADATA, MANIFEST, CLASSIFIERSEQS, THREADS, TRIM_LEFT_F, TRIM_LEFT_R, etc.

# Define a unique output directory for this run
RUN_TIMESTAMP=$(date +"%Y-%m-%d_%H-%M-%S")
OUTPUT_DIR="output_run_${RUN_TIMESTAMP}"
mkdir -p $OUTPUT_DIR/qiime_artifacts
mkdir -p $OUTPUT_DIR/exported_data
mkdir -p $OUTPUT_DIR/visualizations

# Step 1: Import data into QIIME2
echo "Importing data into QIIME2..."
singularity exec $QIIME2_SIF qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path $MANIFEST \
    --output-path $OUTPUT_DIR/qiime_artifacts/paired-end-demux.qza \
    --input-format PairedEndFastqManifestPhred33V2

# Step 2: Trim primers using cutadapt
echo "Trimming primers..."
singularity exec $QIIME2_SIF qiime cutadapt trim-paired \
    --i-demultiplexed-sequences $OUTPUT_DIR/qiime_artifacts/paired-end-demux.qza \
    --p-cores $THREADS \
    --p-front-f CCTACGGGNGGCWGCAG \
    --p-front-r GACTACHVGGGTATCTAATCC \
    --o-trimmed-sequences $OUTPUT_DIR/qiime_artifacts/trimmed-demux.qza \
    --verbose

# Step 3: Summarize the trimmed sequences
echo "Summarizing trimmed sequences..."
singularity exec $QIIME2_SIF qiime demux summarize \
    --i-data $OUTPUT_DIR/qiime_artifacts/trimmed-demux.qza \
    --o-visualization $OUTPUT_DIR/visualizations/trimmed-demux-summary.qzv

# Step 4: Denoise using DADA2
echo "Denoising sequences with DADA2..."
singularity exec $QIIME2_SIF qiime dada2 denoise-paired \
    --i-demultiplexed-seqs $OUTPUT_DIR/qiime_artifacts/trimmed-demux.qza \
    --p-trim-left-f $TRIM_LEFT_F \
    --p-trim-left-r $TRIM_LEFT_R \
    --p-trunc-len-f $TRUNC_LEN_F \
    --p-trunc-len-r $TRUNC_LEN_R \
    --p-n-threads $THREADS \
    --o-table $OUTPUT_DIR/qiime_artifacts/table.qza \
    --o-representative-sequences $OUTPUT_DIR/qiime_artifacts/rep-seqs.qza \
    --o-denoising-stats $OUTPUT_DIR/qiime_artifacts/denoising-stats.qza

# Step 5: Summarize feature table and representative sequences
echo "Summarizing feature table and representative sequences..."
singularity exec $QIIME2_SIF qiime feature-table summarize \
    --i-table $OUTPUT_DIR/qiime_artifacts/table.qza \
    --o-visualization $OUTPUT_DIR/visualizations/table.qzv \
    --m-sample-metadata-file $SAMPLEMETADATA

singularity exec $QIIME2_SIF qiime feature-table tabulate-seqs \
    --i-data $OUTPUT_DIR/qiime_artifacts/rep-seqs.qza \
    --o-visualization $OUTPUT_DIR/visualizations/rep-seqs.qzv

singularity exec $QIIME2_SIF qiime metadata tabulate \
    --m-input-file $OUTPUT_DIR/qiime_artifacts/denoising-stats.qza \
    --o-visualization $OUTPUT_DIR/visualizations/denoising-stats.qzv

# Step 6: Classify sequences
echo "Classifying sequences..."
singularity exec $QIIME2_SIF qiime feature-classifier classify-sklearn \
    --i-classifier $CLASSIFIERSEQS \
    --p-n-jobs $THREADS \
    --i-reads $OUTPUT_DIR/qiime_artifacts/rep-seqs.qza \
    --o-classification $OUTPUT_DIR/qiime_artifacts/taxonomy.qza

singularity exec $QIIME2_SIF qiime metadata tabulate \
    --m-input-file $OUTPUT_DIR/qiime_artifacts/taxonomy.qza \
    --o-visualization $OUTPUT_DIR/visualizations/taxonomy.qzv

# Step 7: Export results for downstream analysis
echo "Exporting results for downstream analysis..."
singularity exec $QIIME2_SIF qiime tools export \
    --input-path $OUTPUT_DIR/qiime_artifacts/table.qza \
    --output-path $OUTPUT_DIR/exported_data/

singularity exec $QIIME2_SIF qiime tools export \
    --input-path $OUTPUT_DIR/qiime_artifacts/taxonomy.qza \
    --output-path $OUTPUT_DIR/exported_data/

singularity exec $QIIME2_SIF qiime tools export \
    --input-path $OUTPUT_DIR/qiime_artifacts/rep-seqs.qza \
    --output-path $OUTPUT_DIR/exported_data/

# Convert BIOM table to TSV format
echo "Converting BIOM table to TSV format..."
biom convert \
    --to-tsv \
    --input-fp $OUTPUT_DIR/exported_data/feature-table.biom \
    --output-fp $OUTPUT_DIR/exported_data/feature-table.tsv

# Step 8: Generate Phylogenetic Tree
echo "Generating phylogenetic tree..."
singularity exec $QIIME2_SIF qiime phylogeny align-to-tree-mafft-fasttree \
    --i-sequences $OUTPUT_DIR/qiime_artifacts/rep-seqs.qza \
    --o-alignment $OUTPUT_DIR/qiime_artifacts/aligned-rep-seqs.qza \
    --o-masked-alignment $OUTPUT_DIR/qiime_artifacts/masked-aligned-rep-seqs.qza \
    --o-tree $OUTPUT_DIR/qiime_artifacts/unrooted-tree.qza \
    --o-rooted-tree $OUTPUT_DIR/qiime_artifacts/rooted-tree.qza

# Export the rooted tree for external use
singularity exec $QIIME2_SIF qiime tools export \
    --input-path $OUTPUT_DIR/qiime_artifacts/rooted-tree.qza \
    --output-path $OUTPUT_DIR/exported_data/

echo "Pipeline completed successfully!"
