#!/bin/bash
###############################################################################
# Script: Salmon.sh
# Description:
#   1. Run FastQC on all FASTQ files in the input directory.
#   2. Align RNA-seq reads to a genome or transcriptome using Salmon.
#   3. Generate count and TPM tables for downstream analysis.
#   4. Run MultiQC to summarize FastQC and Salmon outputs.
#
#   To use:
#     - Edit the variables below to point to your FASTQ files, reference
#       transcriptome, and output directory.
###############################################################################

# -------------------------
# === USER VARIABLES ===
# -------------------------
FASTQ_DIR="/path/to/your/fastq_files/*"
TRANSCRIPTOME="/path/to/your/reference_transcriptome.fa"
OUTPUT_DIR="/path/to/your/output_directory"
THREADS=10

# Create all necessary subdirectories inside OUTPUT_DIR
mkdir -p "$OUTPUT_DIR/fastqc"
mkdir -p "$OUTPUT_DIR/quant"

# -------------------------
# === RUN FASTQC ===
# -------------------------
echo "########## Running FastQC on all FASTQ files ##########"
for fq in $FASTQ_DIR/*_1.fq.gz $FASTQ_DIR/*_2.fq.gz; do
    fastqc -o "$OUTPUT_DIR/fastqc" -t "$THREADS" "$fq"
done
echo "FastQC completed. Results in $OUTPUT_DIR/fastqc"

# -------------------------
# === BUILD SALMON INDEX ===
# -------------------------
echo "########## Running Salmon index ##########"
salmon index -t "$TRANSCRIPTOME" -i "$OUTPUT_DIR/Reference_index"

# -------------------------
# === RUN SALMON QUANT ===
# -------------------------
echo "########## Starting Salmon quantification ##########"
for R1 in $FASTQ_DIR/*_1.fq.gz; do
    PREFIX=$(basename "$R1" _1.fq.gz)
    salmon quant \
        -i "$OUTPUT_DIR/Reference_index" \
        --libType A \
        -1 "$FASTQ_DIR/${PREFIX}_1.fq.gz" \
        -2 "$FASTQ_DIR/${PREFIX}_2.fq.gz" \
        -p "$THREADS" \
        --seqBias \
        --useVBOpt \
        --validateMappings \
        -o "$OUTPUT_DIR/quant/${PREFIX}"
done
echo "End of Salmon quantification"

# -------------------------
# === POST-PROCESSING QUANT FILES ===
# -------------------------
cd "$OUTPUT_DIR/quant"
for i in *; do
    cp "$i/quant.sf" "$i.sf"
done

FIRST_FILE=$(ls -1 *.sf | head -1)
awk '{print $1}' "$FIRST_FILE" > geneID

# Create TPM table
for i in *.sf; do
    COLNAME=$(basename "$i" .sf)
    awk '{print $4}' "$i" > "$i.temp"
    sed -i "s/TPM/$COLNAME/g" "$i.temp"
done
paste geneID *.temp > ../TPM.txt

# Create count table
for i in *.sf; do
    COLNAME=$(basename "$i" .sf)
    awk '{print $5}' "$i" > "$i.temp"
    sed -i "s/NumReads/$COLNAME/g" "$i.temp"
done
paste geneID *.temp > ../counts.txt

# Cleanup
rm *.temp geneID

# -------------------------
# === RUN MULTIQC ===
# -------------------------
echo "Launching MultiQC"
cd "$OUTPUT_DIR"
multiqc .

echo "########## Salmon pipeline completed ##########"
