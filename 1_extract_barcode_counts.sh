#!/bin/bash

#build index
subread-buildindex -o esx_rd2 barcodes.fa  #build index. barcodes.fa is fasta file with all strains and respective barcodes

# Loop through all fastq files in the directory
for fastq_file in *.fastq.gz; do
    # Extract the part of the name until the third '_'
    output_prefix=$(echo "$fastq_file" | cut -d '_' -f 1-3)
    
    # Output names
    aligned_bam="${output_prefix}_aligned.bam"
    sorted_bam="${output_prefix}_sorted.bam"
    idxstats_output="${output_prefix}_counts.txt"
    
    # Subread-align
    subread-align -i ../Broad_hypos -r "$fastq_file" -t 1 -o "$aligned_bam"
    
    # BAM file processing
    samtools sort -o "$sorted_bam" "$aligned_bam"
    samtools index "$sorted_bam"
    
    # Samtools counts
    samtools idxstats "$sorted_bam" > "$idxstats_output"
    
    echo "Aligned $fastq_file to $sorted_bam, indexed, and generated idxstats."
done



