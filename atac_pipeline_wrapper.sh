# ATAC-seq pipeline
# Wrapper script
# Author: Stephen Zhang
# Date: 22 Jan 2018
#!/bin/bash

# Usage:
# atac_pipeline_wrapper $CPUS $output_dir $input1 $input2

EXPECTED_ARGS=4
ARGC=$#

ml load macs2
ml load fastqc
ml load cutadapt
ml load bowtie2
ml load picard/2.8.2
ml load samtools

if [ $ARGC -ne $EXPECTED_ARGS ]; then
  echo "atac_pipeline_wrapper"
  echo ""
  echo "Usage: atac_pipeline_wrapper \$CPUS \$OUTPUT_DIR"\
        "\$INPUT_R1.fastq.gz \$INPUT_R2.fastq.gz"
else
  CPUS=$1
  OUTPUT_DIR=$2
  INPUT1=$3
  INPUT2=$4
  bpipe run --dir $OUTPUT_DIR\
    -p num_cpus=$CPUS\
    -p output_dir=$OUTPUT_DIR\
    atac_pipeline.groovy $INPUT1 $INPUT2
fi
