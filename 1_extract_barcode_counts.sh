#!/bin/bash
#BSUB -n 5
#BSUB -R span[hosts=1]
#BSUB -R rusage[mem=10000]
#BSUB -J WGS
#BSUB -q long
#BSUB -W 96:00

#BSUB -o "/home/peter.oluoch-umw/Peter/out/"
#BSUB -e "/home/peter.oluoch-umw/Peter/err/"

module load subread/1.6.2
module load samtools/1.3
#build index
cd /home/peter.oluoch-umw/all_data/Peter_Hypomorphs/RIF_HMP23
#subread-buildindex -o esx_rd2 barcodes.fa  #build index. Include a line on how the barcode should look like

# Loop through all fastq files in the directory
for fastq_file in *.fastq.gz; do
    # Extract the part of the name until the third '_'
    output_prefix=$(echo "$fastq_file" | cut -d '_' -f 1-3)
    
    # Define the output file names
    aligned_bam="${output_prefix}_aligned.bam"
    sorted_bam="${output_prefix}_sorted.bam"
    idxstats_output="${output_prefix}_counts.txt"
    
    # Perform the subread-align command
    subread-align -i ../Broad_hypos -r "$fastq_file" -t 1 -o "$aligned_bam"
    
    # Sort the aligned BAM file
    samtools sort -o "$sorted_bam" "$aligned_bam"
    
    # Index the sorted BAM file
    samtools index "$sorted_bam"
    
    # Generate idxstats output
    samtools idxstats "$sorted_bam" > "$idxstats_output"
    
    echo "Aligned $fastq_file to $sorted_bam, indexed, and generated idxstats."
done



