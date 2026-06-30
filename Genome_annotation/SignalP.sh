#!/bin/bash
###############################################################################
# Script: SignalP.sh
# Description:
#   Run SignalP to predict secreted proteins from a proteome.
#   Generates a list of secreted protein IDs.
###############################################################################

# -------------------------
# === USER VARIABLES ===
# -------------------------
# Sample name
SAMPLE="sample_name"

# Path to your proteome (FASTA format)
PROTEOME="/path/to/${SAMPLE}_proteins.fa"

# Output directory for results
OUTPUT_DIR="./SignalP_results"

mkdir -p "$OUTPUT_DIR"

# -------------------------
# Run SignalP to detect secreted proteins
# -------------------------
echo "=== Running SignalP on $SAMPLE ==="

signalp -fasta "$PROTEOME" \
        -org euk \
        -format short \
        -prefix "$OUTPUT_DIR/Fungi_"

SUMMARY_FILE="$OUTPUT_DIR/Fungi__summary.signalp5"

# Extract secreted protein IDs (Signal peptide prediction)
awk '$2=="SP(Sec/SPI)" {print $1}' "$SUMMARY_FILE" > "$OUTPUT_DIR/secreted_ids.txt"

echo "SignalP prediction completed. Results are in $OUTPUT_DIR"
