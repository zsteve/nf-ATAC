**ATAC-seq pipeline**
Author:		Stephen Zhang (stephen.zhang@monash.edu)
Date:		5 Feb 2018


____

*Introduction*
A pipeline for processing ATAC-seq data written in Nextflow script (https://www.nextflow.io/).
Currently in early stages, this `README` will definitely be updated regularly (check often!)

Have an problem? Please log an issue on GitHub (https://github.com/zsteve/atac-seq-pipeline)

*Dependencies*
Please make sure these tools are installed before running the pipeline:

* `FastQC`
* `cutadapt`
* `bowtie2`
* `picard/2.8.2`
* `samtools`

*Installing Nextflow*

Nextflow can be downloaded by using the following command:

`curl -s https://get.nextflow.io | bash`

This will create a binary `nextflow` in the working directory. You can add this binary to your `PATH` for ease of use:

`export PATH=$PATH:[your path here]

The pipeline can be executed by running `nextflow`, specifying the script and relevant commandline arguments.

`nextflow <script>.nf <command line arguments>`

*Running the pipeline - single sample*

_Data preparation_

Paired-end read sample data in `.fastq.gz` format should be located in a directory with the desired sample name. Read pairs should be distinguishable in the format `*_R{1,2}*.fastq.gz`.

_Command_

Nextflow will create a `work` directory (containing pipeline data) in its working directory (i.e. `.`). Final pipeline output files will be output to a desired directory, however these will generally be _symlinks_ to the actual copy of the file within `work/**/your_file_here`. It is *very* important that `work` does *not* get deleted - otherwise your symlinks will mean nothing!

```

nextflow atac_pipeline.nf --num-cpus $NUM_CPUS 
			  --input-dir $INPUT_DIR
			  --output-dir $OUTPUT_DIR
			  --config-file $CONFIG_FILE
			  --ref-genome-name $GENOME_NAME
			  --ref-genome-index $GENOME_INDEX
			  --ref-genome-fasta $GENOME_FASTA
```

* `NUM_CPUS` - maximum number of CPUs to use for the _entire_ pipeline
* `INPUT_DIR` - path of the directory containing R1,R2 data
* `OUTPUT_DIR` - path of the directory to write outputs to (will be created if it doesn't already exist). This can be the same as INPUT_DIR.
* `CONFIG_FILE` (OPTIONAL) - path to `config.yaml` (in case one wants custom parameters for pipeline components). 
* `GENOME_NAME` - name of the reference genome (e.g. `danRer10`, `hg18`)
* `GENOME_INDEX` - path to `bowtie2` indexes for reference genome
* `GENOME_FASTA` - path to `FASTA` sequence of reference genome

Nextflow will output its data to your directory of choice.

*Running the pipeline - multiple samples*

_Data preparation_

For each sample, create a folder `SAMPLE_ID/` containing the paired-end read data in `fastq.gz` format. Create a *sample table* as a text file:

* Each line corresponds to *one* sample. Fields are as follows:

```
[Sample_ID] [path to sample input directory] [path to sample output directory]
```

_Command_

Pipeline will read in samples from the sample table `.txt` file and attempt to process those samples in _parallel_. 

```
nextflow atac_pipeline.nf --num-cpus $NUM_CPUS
			  --config-file $CONFIG_FILE
			  --multi-sample
			  --sample-table $SAMPLE_TABLE
			  --ref-genome-name $GENOME_NAME
			  --ref-genome-index $GENOME_INDEX
			  --ref-genome-fasta $GENOME_FASTA

```


