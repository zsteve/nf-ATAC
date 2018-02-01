**ATAC-seq pipeline**
Author:		Stephen Zhang
Date:		1 Feb 2018

____

*Introduction*
This is a pipeline for processing of ATAC-seq data written in Nextflow script. Currently in early stages, this `README` will definitely be updated regularly!

*Scripts*

Nextflow script can be used directly:
```
nextflow atac_pipeline.nf [--help] --input-dir <...>
			--num_cpus $NUMCPUS
            --ref-genome-index $REFGENOME_INDEX
            --ref-genome-fasta $REFGENOMEFA
            --ref-genome-name $REFGENOMENAME

```

* `--input-dir` specifies a directory containing paired-end `.fastq.gz` files matching the glob `*_R{1,2}*.fastq.gz`.
* `atac_pipeline.nf` will produce a set of folders containing outputs inside the input directory.

Alternatively, if you have a set of _N_ samples following the directory hierarchy, use `atac_pipeline_multisample_wrapper.sh`
```
sample_dir/
	sample1/
    	*_R{1,2}*.fastq.gz
    ...
    sampleN/
    	*R{1,2}*.fastq.gz
```

`atac_pipeline_multisample_wrapper.sh` takes arguments specified as follows:

```
atac_pipeline_multisample_wrapper.sh [sample_dir]
			[num_cpus]
            [ref_genome_name]
            [ref_genome_index]
            [ref_genome_fasta]
```

* In this instance, `[sample_dir]` points the _top level_ directory. Subdirectories represent samples, which in turn contain `.fastq.gz` read data.
