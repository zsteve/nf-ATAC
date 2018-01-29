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
ref_genome_index="/home/szha0069/Danio_rerio/UCSC/danRer10/Sequence/Bowtie2Index/genome"
ref_genome_fasta="~/Danio_rerio/UCSC/danRer10/Sequence/WholeGenomeFasta/genome.fa"

run_id = "default"

trim_reads = {
  /* trim reads for both R1 and R2 */
  output.dir = output_dir + '/' + 'trim_reads'
  transform('*.fastq.gz') to ('_trimmed.fastq.gz'){
    exec """
        cutadapt -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT
        -q 20
        --minimum-length 36
        -o $output1
        -p $output2
        $input1
        $input2
        | tee $output.dir/_trim_reads_output
        """
  }
}

map_to_reference = {
  /* map reads for both R1 and R2 */
  output.dir = output_dir +'/map_to_reference'
  produce(run_id + '.sam'){
    exec """
      bowtie2 -X2000 --no-mixed --no-discordant -p $num_cpus
      -x $ref_genome_index
      -1 $input1
      -2 $input2
      -S $output
      3>&1 1>&2 2>&3 | tee $output.dir/_map_to_reference_output
    """
  }
}

get_map_stats = {
  /* use samtools flagstat to retrieve mapping stats */
  output.dir = output_dir + '/map_to_reference'
  produce(run_id + '_samtools_flagstat.txt'){
      exec """
        samtools flagstat $input > $output
      """
  }
}

sam_to_bam = {
  /* convert sam (reads_mapped.sam) to bam */
  output.dir = output_dir + '/' + 'sam_to_bam'
  transform('*.sam') to ('.bam'){
    exec """ picard SortSam SO=coordinate
		INPUT=$input
    OUTPUT=$output
    CREATE_INDEX=true
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
      REMOVE_DUPLICATES=true
      ASSUME_SORTED=true
      VALIDATION_STRINGENCY=LENIENT
      CREATE_INDEX=true
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
			-fsize 1e10
			-name $run_id
      | tee $output.dir/_make_ucsc_gb_track_output
	"""
  out_track_name = out_track_name+'.gz'
  forward out_track_name
}

// call with output *.bam from make_tags
macs2_call_peaks = {
  output.dir = output_dir + '/' + 'macs2_call_peaks'
  exec """
    macs2 callpeak --nomodel
    -t $input
    -n $run_id
    --keep-dup all
    --gsize 11250000000
    --shift -37
    --extsize 73
    --outdir $output.dir
    -B
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
      -style factor
      -o $output
    """
  }
}

// annotate peaks using homer annotatePeaks.pl
homer_annotate_peaks = {
  output.dir = output_dir + '/homer_annotate_peaks'
  exec """
    $homer_path/annotatePeaks.pl $input danRer10 > $output
  """
}

get_id = {
  /* obtains the filename before _R1 or _R2 */
  run_id = "basename $input".execute().text
  run_id = run_id[0..run_id.indexOf('_R')-1]
}

run {
  /* get run ID */
  get_id + trim_reads + map_to_reference + [get_map_stats, sam_to_bam] + filter_reads_1 + \
  filter_reads_2 + [make_tags + make_ucsc_gb_track, macs2_call_peaks + homer_annotate_peaks,\
  make_tags + homer_call_peaks + homer_annotate_peaks]
}
