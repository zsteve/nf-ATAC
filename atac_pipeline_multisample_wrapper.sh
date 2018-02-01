#!/bin/bash

# *.sh [sample_directory] [num_cpus] [ref_genome_name] [ref_genome_index] [ref_genome_fasta]

sample_dir=$1
num_cpus=$2
ref_genome_name=$3
ref_genome_index=$4
ref_genome_fasta=$5

ml load macs2
ml load fastqc
ml load cutadapt
ml load bowtie2
ml load picard/2.8.2
ml load samtools

nf_script_path="./atac_pipeline.nf"

for i in $sample_dir/*/; do
  echo $i
  nextflow $nf_script_path --input-dir $i\
      --num-cpus $num_cpus\
      --ref-genome-index $ref_genome_index\
      --ref-genome-fasta $ref_genome_fasta\
      --ref-genome-name $ref_genome_name &
done
