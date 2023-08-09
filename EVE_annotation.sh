#!/bin/bash

# Input files
db1=$1 # Viral proteins from nr without retrovirus sequences <fasta/prot>
db2=$2 # NR-NCBI proteins <fasta/prot>
hostFasta=$3 # host genome <fasta/nucl>

# Filenames for intermediate files
host_dmd1=dmd1_$(basename "$hostFasta" | cut -d'.' -f1).dmd
host_bed1=dmd1_$(basename "$hostFasta" | cut -d'.' -f1).bed
host_dmdFasta1dmd1_=$(basename "$hostFasta" | cut -d'.' -f1).fasta

host_dmd2=dmd2_$(basename "$hostFasta" | cut -d'.' -f1).dmd
host_bed2=dmd2_$(basename "$hostFasta" | cut -d'.' -f1).bed

# Run Diamond for the first step
diamond blastx --more-sensitive --max-target-seqs 1000 --range-culling --min-score 40 --outfmt 6 -b 30 -c 1 -d "$db1" -q "$hostFasta" -o "$host_dmd1" 

# Run custom Python script
python3 bestHitsToFasta.py "$host_dmd1" "$host_bed1" 

# Extract genomic regions with viral best hits
bedtools getfasta -fi "$hostFasta" -bed "$host_bed1" -fo "$host_dmdFasta1"

# Run Diamond for the second step
diamond blastx --more-sensitive --max-target-seqs 1000  --range-culling  --min-score 40 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle -b 30 -c 1 -d "$db2" -q "$host_dmdFasta1" -o "$host_dmd2"


# Run taxonomizR
Rscript --vanilla absolutCoord.R "$host_dmd2" "$host_bed2"
Rscript --vanilla Get_EVE_annoation_summary.r # Adapt path in Rscript file
