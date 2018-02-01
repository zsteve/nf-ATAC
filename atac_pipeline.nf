params.help = false
params.inputDir = ''
params.numCpus = 8
params.refGenomeIndex = ''
params.refGenomeFasta = ''
params.refGenomeName=''

/*
nextflow atac_pipeline.nf --inpu
t-dir test_samples/EDM_TAGGCA_L001/ --num-cpus 8
--ref-genome-index /home/szha0069/reference_genomes/Danio_rerio/UCSC/danRer10/Sequence/Bowtie2Index/genome
--ref-genome-fasta /home/szha0069/reference_genomes/Danio_rerio/UCSC/danRer10/Sequence/WholeGenomeFasta/genome.fa
 --ref-genome-name danRer10 -resume
*/

params.jvarkitPath = "~/tools/jvarkit/dist"
params.homerPath = "/home/szha0069/tools/homer/bin"

inputDir = params.inputDir

if(params.help){
  // help here
}

if(inputDir.endsWith('/')){
  inputDir = inputDir.substring(0, inputDir.size()-1)
  println inputDir
}

/* call on directory containing read pairs */
Channel.fromFilePairs("$inputDir/*_R{1,2}*.fastq.gz")
      .into{readPairOut; readPairOut_FQC}

process sampleFastQC {
  //echo true
  publishDir inputDir
  input:
    set val(name), file(reads) from readPairOut_FQC
  output:
    file 'fastqc_output' into fastQCReport
  script:
    """mkdir -p fastqc_output;
    fastqc -t $params.numCpus -o fastqc_output $reads \
      > fastqc.stdout 2> fastqc.stderr
    """
}

process sampleCutadapt {
  //echo true
  input:
    set val(name), file(reads) from readPairOut
  output:
    set val(name), file('cutadapt_output/*_trimmed*') into trimmedPairOut
  // obtain read1 and read2 filenames
  script:
    read1_trimmed_name = file(reads[0]).getName().toString().replaceAll(/_R1/, '_R1_trimmed')
    read2_trimmed_name = file(reads[1]).getName().toString().replaceAll(/_R2/, '_R2_trimmed')
    """
    mkdir -p cutadapt_output; \
    cutadapt \
    -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT\
    --info-file=cutadapt_output/cutadapt_info_file\
    -q 20\
    --minimum-length 36\
    -o cutadapt_output/${read1_trimmed_name}\
    -p cutadapt_output/${read2_trimmed_name}\
    $reads > cutadapt_output/cutadapt.stdout 2> cutadapt_output/cutadapt.stderr
    """
}

process sampleMapToReference {
  /* use bowtie2 for mapping */
  //echo true
  publishDir inputDir
  input:
    set val(name), file(reads) from trimmedPairOut
  output:
    set val(name), file('bowtie2_output/*.sam') into mappedSamOut
  script:
    output_sam_name = name+'.sam'
    """
    mkdir -p bowtie2_output;
    bowtie2 -X2000 --no-mixed --no-discordant\
    -p $params.numCpus\
    -x $params.refGenomeIndex\
    -1 ${reads[0]}\
    -2 ${reads[1]}\
    -S bowtie2_output/$output_sam_name\
    > bowtie2_output/bowtie2.stdout 2> bowtie2_output/bowtie2.stderr
    """
}

process sampleSamToBam {
  //echo true
  input:
    set val(name), file(samFile) from mappedSamOut
  output:
    set val(name), file('*.bam') into mappedBamOut_getStat, mappedBamOut
  script:
    output_bam_name = name + '.bam'
    """
      picard SortSam SO=coordinate INPUT=$samFile OUTPUT=$output_bam_name\
      CREATE_INDEX=true\
      > picard_sortsam.stdout 2> picard_sortsam.stderr
    """
}

process sampleGetMapStats {
  //echo true
  publishDir inputDir
  input:
    set val(name), file(bamFile) from mappedBamOut_getStat
  output:
    file 'flagstat_output/*_flagstat.txt' into flagStatOut
  script:
    output_flagstat_name = name + '_flagstat.txt'
    """
      mkdir -p flagstat_output;
      samtools flagstat $bamFile > flagstat_output/$output_flagstat_name\
      2> samtools_flagstat.stderr
    """
}

process sampleFilterMMMR {
  // filter multimappers and mitochondrial reads
  //echo true
  input:
    set val(name), file(bamFile) from mappedBamOut
  output:
    set val(name), file('*_filtered_MMMR.bam') into filteredMMMROut
  script:
    output_bam_name = name + '_filtered_MMMR.bam'
    """
      java -jar $params.jvarkitPath/samjs.jar\
      --samoutputformat bam\
      -o $output_bam_name\
      -e '!(record.getReadUnmappedFlag() || record.getReferenceName().equals("chrM") || record.getMappingQuality() < 10)'\
      $bamFile > jvarkit_samjs.stdout 2> jvarkit_samjs.stderr
    """
}

process sampleFilterDedup {
  // discard duplicate reads
  input:
    set val(name), file(bamFile) from filteredMMMROut
  output:
    set val(name), file('*_dedup.bam') into filteredDedupOut_makeTags, \
                          filteredDedupOut_callPeaks_MACS
  script:
    output_bam_name = name + '_dedup.bam'
    output_metrics_name = name + '_dedup_metrics'
    """
      picard MarkDuplicates\
      INPUT=$bamFile\
      OUTPUT=$output_bam_name\
      METRICS_FILE=$output_metrics_name\
      REMOVE_DUPLICATES=true\
      ASSUME_SORTED=true\
      VALIDATION_STRINGENCY=LENIENT\
      CREATE_INDEX=true\
      > picard_markduplicates.stdout 2> picard_markduplicates.stderr
    """
}


process sampleMakeTags {
  publishDir inputDir
  input:
    set val(name), file(bamFile) from filteredDedupOut_makeTags
  output:
    set val(name), file('maketags_output') into makeTagsOut, makeTagsOut_callPeaks_homer
  script:
    """
    $params.homerPath/makeTagDirectory\
    maketags_output\
    -genome $params.refGenomeFasta\
    -checkGC $bamFile\
    > maketagdirectory.stdout 2> maketagdirectory.stderr
    """
}

process sampleMakeUCSCTrack {
  //echo true
  publishDir inputDir
  input:
    set val(name), file(tagDir) from makeTagsOut
  output:
    file "makeUCSCtrack_output/*.bedgraph.gz" into makeUCSCTrackOut
  script:
    output_track_name = name + '_UCSC.bedgraph'
    """
    mkdir -p makeUCSCtrack_output;
    $params.homerPath/makeUCSCfile $tagDir -o makeUCSCtrack_output/$output_track_name\
      -fsize 1e10\
      -name $output_track_name\
      > makeucsctrack.stdout 2> makeucsctrack.stderr
    """
}

process sampleMACS2CallPeaks {
  publishDir inputDir
  input:
    set val(name), file(bamFile) from filteredDedupOut_callPeaks_MACS
  output:
    set val(name), file("macs2_output/*.xls") into outputMACS2_xlsPeaks
  script:
    """
    mkdir -p macs2_output;
    macs2 callpeak --nomodel\
    -t $bamFile\
    -n $name\
    --keep-dup all\
    --gsize 11250000000\
    --shift -37\
    --extsize 73\
    --outdir macs2_output\
    -B > macs2_output/macs2.stdout 2> macs2_output/macs2.stderr
    """
}

process sampleHomerCallPeaks {
  publishDir inputDir
  errorStrategy 'ignore'
  input:
    set val(name), file(tagDir) from makeTagsOut_callPeaks_homer
  output:
    file 'homer_findpeaks_output/*.txt' into outputHomerFindpeaks
  script:
    output_name = name + '_homer_findpeaks.txt'
  """
    $params.homerPath/findPeaks\
    $tagDir
    -style factor\
    -o homer_findpeaks_output/$output_name
    > homer_findpeaks_output/homer_findpeaks.stdout 2> homer_findpeaks_output/homer_findpeaks.stderr
  """
}

process sampleAnnotatePeaks {
  publishDir inputDir
  input:
    set val(name), file(xlsFile) from outputMACS2_xlsPeaks
  output:
    set val(name), file("homer_annotatepeaks_output/$annotationFile") into outputAnnotatePeaks
  script:
    annotationFile = name + '_annotation.txt'
    """
    mkdir -p homer_annotatepeaks_output;
    $params.homerPath/annotatePeaks.pl\
      $xlsFile\
      $params.refGenomeName\
      > homer_annotatepeaks_output/$annotationFile 2> homer_annotatepeaks_output/homer_annotatepeaks_output.stderr
    """
}
