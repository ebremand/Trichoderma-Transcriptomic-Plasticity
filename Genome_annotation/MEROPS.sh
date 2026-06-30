#!/bin/bash
###############################################################################
# Script: Merops.sh
# Description:
#   Identify peptidases in a proteome using MEROPS v12.4 via BLASTp.
#   Produces a final table with Query | PeptidaseFamily | CatalyticType.
#
#   To use:
#     - Edit the variables below to point to your proteome, output folder, and MEROPS database.
###############################################################################

# -------------------------
# === USER VARIABLES ===
# -------------------------
# Sample name
SAMPLE="sample_name"

# Path to your proteome (FASTA format)
PROTEOME="/path/to/${SAMPLE}_proteins.fa"

# Output directory for results
OUTPUT_DIR="./Merops_results"

# Path to MEROPS database (v12.4)
MEROPS_DB="/path/to/Merops_db/meropsscan.lib"

# Number of threads
THREADS=8

mkdir -p "$OUTPUT_DIR"

# -------------------------
# Run BLASTp against MEROPS
# -------------------------
echo "=== Running BLASTp for $SAMPLE against MEROPS v12.4 ==="
BLAST_OUT="$OUTPUT_DIR/${SAMPLE}_merops.blast.out"

blastp -query "$PROTEOME" \
       -db "$MEROPS_DB" \
       -out "$BLAST_OUT" \
       -evalue 1e-5 \
       -num_threads "$THREADS" \
       -outfmt 6

echo "BLASTp completed."

# -------------------------
# Annotate peptidase family and catalytic type, remove duplicates
# -------------------------
echo "=== Annotating peptidase families and catalytic types ==="
FINAL_OUT="$OUTPUT_DIR/${SAMPLE}_merops_final.tsv"

# Create header
echo -e "Query\tPeptidaseFamily\tCatalyticType" > "$FINAL_OUT"

# Process BLAST output
awk -v merops="$MEROPS_DB" '
{
    merid=$2
    query=$1
    # Search corresponding line in MEROPS database
    cmd="grep -m1 " merid " " merops
    cmd | getline header
    close(cmd)
    # Extract peptidase family
    split(header,a,"\\[|\\]"); family=a[2]
    # Extract catalytic type
    split(header,b,"#"); cat_type=b[2]
    # Avoid duplicates: one line per Query + Family
    key=query"\t"family
    if(!(key in seen)){
        print query"\t"family"\t"cat_type
        seen[key]=1
    }
}' "$BLAST_OUT" >> "$FINAL_OUT"

echo "Annotation completed. Final table is in $FINAL_OUT"
