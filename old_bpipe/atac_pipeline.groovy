/* ATAC pipeline
  Uses bpipe (http://docs.bpipe.org/)
  Author: Stephen Zhang

  Date: 22 Jan 2018

  Usage (paired reads):
  bpipe run -p num_cpus=$CPU -p output_dir=$OUTPUT_DIR <input_R1> <input_R2>

 */

jvarkit_path="~/tools/jvarkit/dist"
homer_path="/home/szha0069/tools/homer/bin"

/* reference file paths */
ref_genome_index="/home/szha0069/reference_genomes/Danio_rerio/UCSC/danRer10/Sequence/Bowtie2Index/genome"
ref_genome_fa = "/home/szha0069/reference_genomes/Danio_rerio/UCSC/danRer10/Sequence/WholeGenomeFasta/genome.fa"

run_id = "default"

trim_reads = {
  /* trim reads for both R1 and R2 */
  output.dir = output_dir + '/' + 'trim_reads'
  transform('*.fastq.gz') to ('_trimmed.fastq.gz'){
    exec """
        cutadapt \$(cat params/cutadapt.params)
        -o $output1
        -p $output2
        $input1
        $input2
        | tee $output.dir/_trim_reads_output
        """
  }
}

fastqc = {
  /* run fastqc first on input files */
  output.dir = output_dir + '/fastqc'
  exec """
    fastqc \$(cat params/fastqc.params)
    -o $output.dir $input1 $input2
  """
}

map_to_reference = {
  /* map reads for both R1 and R2 */
  output.dir = output_dir +'/map_to_reference'
  produce(run_id + '.sam'){
    exec """
      bowtie2 \$(cat params/bowtie2.params)
      -p $num_cpus
      -x $ref_genome_index
      -1 $input1
      -2 $input2
      -S $output
      3>&1 1>&2 2>&3 | tee $output.dir/$run_id'_map_to_reference_output'
    """
  }
}

get_map_stats = {
  /* use samtools flagstat to retrieve mapping stats */
  output.dir = output_dir + '/map_to_reference'
  produce(run_id + '_samtools_flagstat.txt'){
      exec """
        samtools flagstat \$(cat params/samtools_flagstat.params) $input > $output
      """
  }
}

sam_to_bam = {
  /* convert sam (reads_mapped.sam) to bam */
  output.dir = output_dir + '/' + 'sam_to_bam'
  transform('*.sam') to ('.bam'){
    exec """ picard SortSam
    INPUT=$input
    OUTPUT=$output
    \$(cat params/picard_sortsam.params)
    """
  }
}


/* filter_reads_1 and filter_reads_2 provide read filtering */

filter_reads_1 = {
  /* discard mitoc reads and multimappers */
  output.dir = output_dir + '/' + 'filter_reads'
  transform('*.bam') to ('_filtered.bam'){
  exec """
    java -jar $jvarkit_path/samjs.jar
    --samoutputformat bam
    -o $output.bam
    -e '!(record.getReadUnmappedFlag() || record.getReferenceName().equals("chrM") || record.getMappingQuality() < 10)'
    $input.bam
    | tee $output.dir/_filter_reads_1_output
    """
  }
}

filter_reads_2 = {
  /* now discard duplicates */
  output.dir = output_dir + '/' + 'filter_reads'
  transform('*.bam') to ('_dedup.bam', '_dedup.metrics'){
    exec """
      picard MarkDuplicates
      INPUT=$input.bam
      OUTPUT=$output.bam
      METRICS_FILE=$output.metrics
      \$(cat params/picard_markduplicates.params)
      | tee $output.dir/_filter_reads_2_output
      """
  }
}

/* after filter_reads_1 and filter_reads_2 */

make_tags = {
  output.dir = output_dir + '/' + 'make_tags'
  produce(run_id+'_tags'){
    exec """
      $homer_path/makeTagDirectory
      $output
      -genome $ref_genome_fasta
      -checkGC $input.bam
      \$(params/homer_maketagdirectory.params)
      | tee $output.dir/_make_tags_output
  	"""
  }
}

make_ucsc_gb_track = {
  output.dir = output_dir + '/' + 'ucsc_track'
  out_track_name = output.dir+'/'+run_id+'_track'
	exec """
		$homer_path/makeUCSCfile
    $input
    -o $out_track_name
    -name $run_id
    \$(params/homer_makeucscfile.params)
      | tee $output.dir/_make_ucsc_gb_track_output
	"""
  out_track_name = out_track_name+'.gz'
  forward out_track_name
}

// call with output *.bam from make_tags
macs2_call_peaks = {
  output.dir = output_dir + '/' + 'macs2_call_peaks'
  exec """
    macs2 callpeak
    -t $input
    -n $run_id
    --outdir $output.dir
    \$(cat params/macs2_callpeaks.params)
  """
  forward(output.dir + '/'+run_id+'_summits.bed')
}

// call with output from make_tags
homer_call_peaks = {
  produce(run_id+'_homer_peaks.txt'){
    output.dir = output_dir + '/homer_call_peaks'
    exec """
      $homer_path/findPeaks
      $input.dir
      -o $output
      \$(cat params/homer_findpeaks.params)
    """
  }
}

// annotate peaks using homer annotatePeaks.pl
homer_annotate_peaks = {
  output.dir = output_dir + '/homer_annotate_peaks'
  exec """
    $homer_path/annotatePeaks.pl
    $input
    \$(cat params/homer_annotatepeaks.params) > $output
  """
}

get_id = {
  /* obtains the filename before _R1 or _R2 */
  run_id = "basename $input".execute().text
  run_id = run_id[0..run_id.indexOf('_R')-1]
}

multiqc = {
  output.dir = output_dir
  exec """
    multiqc --outdir $output.dir $output_dir
  """
}

/****
Pipeline schematic
  - trim reads using cutadapt
  - map to reference genome using bowtie2
    - get mapping stats
  - convert .sam to .bam file
  - filter reads (filter_reads_1/2)
    - make tags -> UCSC visualisation track
    - call peaks with macs2 -> annotate peaks
    - call peaks with homer -> annotate peaks
****/

run {
  /* get run ID */
  get_id + [trim_reads + map_to_reference + [get_map_stats, sam_to_bam] + filter_reads_1 + \
  filter_reads_2 + [make_tags + make_ucsc_gb_track, macs2_call_peaks + homer_annotate_peaks,\
  make_tags + homer_call_peaks + homer_annotate_peaks], fastqc] + multiqc
}
