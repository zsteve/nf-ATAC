/*
 * nf-ATAC - ATAC-seq pipeline
 * Written with Nextflow (https://www.nextflow.io)
 * Author: Stephen Zhang (stephen.zhang@monash.edu)
 * Date: 5 Feb 2018
 *    
 */

import org.yaml.snakeyaml.Yaml

params.help = false
/* multisample flag - if set, then we will handle sample table */
params.multiSample = false
params.sampleTable = '' // must pass if we want to process multiple samples
params.inputDir = ''
params.outputDir = ''
params.numCpus = 8
params.refGenomeIndex = ''
params.refGenomeFasta = ''
params.refGenomeName=''
params.configFile = ''
params.jvarkitPath = ''

// default command-line parameters for tools we want to use in pipeline

def def_cmd_params = [:]
def_cmd_params["fastqc"] = [:]
def_cmd_params["cutadapt"] = ["-a":"CTGTCTCTTATACACATCT",
                              "-A":"CTGTCTCTTATACACATCT",
                              "-q":20,
                              "--minimum-length":36]
def_cmd_params["bowtie2"] = ["-p":params.numCpus,
                             "-x":params.refGenomeIndex,
                             "-X2000":true,
                             "--no-mixed":true,
                             "--no-discordant":true]

def_cmd_params["flagstat"] = [:]

def_cmd_params["samjs"] = ["--samoutputformat":"bam",
                           "-e": '!(record.getReadUnmappedFlag() || record.getReferenceName().equals("chrM") || record.getMappingQuality() < 10)']

def_cmd_params["markduplicates"] = ["REMOVE_DUPLICATES":"true",
                                    "ASSUME_SORTED":"true",
                                    "VALIDATION_STRINGENCY":"LENIENT",
                                    "CREATE_INDEX":"true"]

def_cmd_params["maketagdirectory"] = [:]

def_cmd_params["makeucscfile"] = ["-fsize":"1e10"]

def_cmd_params["macs2"] = ["--nomodel":true,
                           "--keep-dup":"all",
                           "--gsize":"11250000000",
                           "--shift":"-37",
                           "--extsize":"73",
                           "-B":true]

def_cmd_params["findpeaks"] = ["-style":"factor"]

def_cmd_params["annotatepeaks"] = [:]

if(params.configFile != ''){
    Yaml yaml = new Yaml()
    String configFile = new File(params.configFile).text
    Map configMap = (Map)yaml.load(configFile)

    /*
      Any parameters that are set iin config.yaml will overwrite def_cmd_params
      These will in turn be overwritten if one specifies custom sample-specific params
    */
    for(tool in def_cmd_params.keySet()){
        if(configMap.containsKey(tool)){
            for(option in def_cmd_params[tool].keySet()){
                if(configMap[tool].containsKey(option)){
                    def_cmd_params[tool][option] = configMap[tool][option]
                }
            }
        }else{
            // do nothing
        }
    }
}

/*
nextflow atac_pipeline.nf --input-dir test_samples/EDM_TAGGCA_L001/ --num-cpus 8 --ref-genome-index /home/szha0069/reference_genomes/Danio_rerio/UCSC/danRer10/Sequence/Bowtie2Index/genome --ref-genome-fasta /home/szha0069/reference_genomes/Danio_rerio/UCSC/danRer10/Sequence/WholeGenomeFasta/genome.fa --ref-genome-name danRer10 -resume
*/


def getFlagString(Map param_map, String tool, String option){
  return param_map[tool][option] ? option : ""
}

def getReadPairFromSampleDir(sampleDir){
  /* does what it says!
    sampleDir is a string specifying the path of the base directory
    Returns [read1_path, read2_path]
    */
  if(sampleDir.endsWith('/')){
    sampleDir = sampleDir.substring(0, sampleDir.length()-1)
  }

  def files = file(sampleDir).list()
  def read1_path = ''
  def read2_path = ''
  for( i in files ){
    if(i.contains('_R1')){
      read1_path += i
      break
    }
  }
  for( i in files ){
    if(i.contains('_R2')){
      read2_path+=i
      break
    }
  }
  // sanity check to see if R1 and R2 were both found!
  if(read1_path == '' || read2_path == ''){
    println "Error - didn't find R1/R2 pair in directory: $sampleDir"
  }
  read1_path = sampleDir + '/' + read1_path
  read2_path = sampleDir + '/' + read2_path
  return ["R1":read1_path, "R2":read2_path]
}

def regularizeDirPath(String path){
    /* regularise directory path - removes trailing '/' if there is one. */
    if(path.endsWith('/')){
        path = path.substring(0, path.size()-1)
    }
    return path
}

if(params.multiSample){
  /*
    Load sample table
    Format:
    [Sample_ID] [Sample_input_directory] [Sample_output_directory]
  */
  sampleTable = file(params.sampleTable)
  Channel.from(sampleTable.readLines())
        .map{
          line ->
          def entries = line.tokenize()
          read_paths = getReadPairFromSampleDir(entries[1])
      // create the output directory if it doesn't already exist
      file(entries[2]).mkdir()
          // format:
          // [[ID=sampleID, baseDir=sample_folder], [read1, read2]]
          [["ID":entries[0], "baseDir":entries[1], "baseDirOut":entries[2]], [file(read_paths["R1"]), file(read_paths["R2"])]]
        }
        .into{readPairOut; readPairOut_FQC}
}else{
  inputDir = regularizeDirPath(params.inputDir)
  outputDir = regularizeDirPath(params.outputDir)
  /* call on directory containing read pairs */
  // Channel.fromFilePairs("$inputDir/*_R{1,2}*.fastq.gz")
  //       .into{readPairOut; readPairOut_FQC}
  println([inputDir, outputDir])
  Channel.value([inputDir, outputDir])
    .map{dir ->
        def inputDir = dir[0]
    def outputDir = dir[1]
    def sampleID = file(inputDir).getBaseName().toString()
        file(outputDir).mkdir()
    read_paths = getReadPairFromSampleDir(inputDir)
        [["ID":sampleID, "baseDir":inputDir, "baseDirOut":outputDir], [file(read_paths["R1"]), file(read_paths["R2"])]]}
    .into{readPairOut; readPairOut_FQC}
}

if(params.help){
  log.info 'nf-ATAC - An integrated pipeline for ATAC-seq data' 
  log.info '(c) 2018 Stephen Zhang (stephen.zhang@monash.edu)'
  log.info 'Usage (single sample):'
  log.info 'nextflow nf_atac.nf'
  log.info '[--help]'
  log.info '--num-cpus \$NUM_CPUS'
  log.info '--input-dir \$INPUT_DIR'
  log.info '--output-dir \$OUTPUT_DIR'
  log.info '[--config-file \$CONFIG_FILE]'
  log.info '--ref-genome-name \$GENOME_NAME'
  log.info '--ref-genome-index \$GENOME_INDEX'
  log.info '--ref-genome-fasta \$GENOME_FASTA'
  log.info 'Usage (multiple samples):'
  log.info 'nextflow nf_atac.nf'
  tog.info '--num-cpus \$NUM_CPUS'
  log.info '--config-file \$CONFIG_FILE'
  log.info '--multi-sample'
  log.info '--sample-table \$SAMPLE_TABLE'
  log.info '--ref-genome-name \$GENOME_NAME'
  log.info '--ref-genome-index \$GENOME_INDEX'
  log.info '--ref-genome-fasta \$GENOME_FASTA'
  exit 0
}

process sampleFastQC {
  //echo true
  publishDir {sampleInfo["baseDirOut"]}
  input:
    set val(sampleInfo), file(reads) from readPairOut_FQC
  output:
    file '*'
    file 'fastqc_output' into fastQCReport
  script:
    """mkdir -p fastqc_output;
    fastqc -t $params.numCpus -o fastqc_output $reads \
      > fastqc.stdout 2> fastqc.stderr
    """
}

process sampleCutadapt {
  // echo true
  publishDir {sampleInfo["baseDirOut"]}
  input:
    set val(sampleInfo), file(reads) from readPairOut
  output:
    file 'cutadapt_output/*'
    set val(sampleInfo), file('cutadapt_output/*_trimmed*') into trimmedPairOut
  // obtain read1 and read2 filenames
  script:
    read1_trimmed_name = file(reads[0]).getName().toString().replaceAll(/_R1/, '_R1_trimmed')
    read2_trimmed_name = file(reads[1]).getName().toString().replaceAll(/_R2/, '_R2_trimmed')
    """
    mkdir -p cutadapt_output; \
    cutadapt \
    -a ${def_cmd_params["cutadapt"]["-a"]} -A ${def_cmd_params["cutadapt"]["-A"]}\
    --info-file=cutadapt_output/cutadapt_info_file\
    -q ${def_cmd_params["cutadapt"]["-q"]}\
    --minimum-length ${def_cmd_params["cutadapt"]["--minimum-length"]}\
    -o cutadapt_output/${read1_trimmed_name}\
    -p cutadapt_output/${read2_trimmed_name}\
    $reads > cutadapt_output/cutadapt.stdout 2> cutadapt_output/cutadapt.stderr
    """
}

process sampleMapToReference {
  /* use bowtie2 for mapping */
  //echo true
  publishDir {sampleInfo["baseDirOut"]}
  input:
    set val(sampleInfo), file(reads) from trimmedPairOut
  output:
    file 'bowtie2_output/*'
    set val(sampleInfo), file('bowtie2_output/*.sam') into mappedSamOut
  script:
    output_sam_name = sampleInfo["ID"]+'.sam'
    """
    mkdir -p bowtie2_output;
    bowtie2 ${getFlagString(def_cmd_params, "bowtie2", "-X2000")}\
    ${getFlagString(def_cmd_params, "bowtie2", "--no-mixed")}\
    ${getFlagString(def_cmd_params, "bowtie2", "--no-discordant")}\
    -p ${def_cmd_params["bowtie2"]["-p"]}\
    -x ${def_cmd_params["bowtie2"]["-x"]}\
    -1 ${reads[0]}\
    -2 ${reads[1]}\
    -S bowtie2_output/$output_sam_name\
    > bowtie2_output/bowtie2.stdout 2> bowtie2_output/${sampleInfo["ID"]}_bowtie2.out
    """
    // here we had to move bowtie2.stderr (where bowtie2 output stats are fed) to a named
    // output file because we want multiQC to pick up on the sample name
}

process sampleSamToBam {
  //echo true
  publishDir {sampleInfo["baseDirOut"]}
  echo true
  input:
    set val(sampleInfo), file(samFile) from mappedSamOut
  output:
    file 'samtobam_output/*'
		file 'samtobam_output/*.bam.bai'
    set val(sampleInfo), file('samtobam_output/*.bam') into mappedBamOut_getStat,mappedBamOut
    set val(sampleInfo), file('samtobam_output/*.bam'), file('samtobam_output/*.bai') into mappedBamOut_QC
  script:
    output_name = sampleInfo["ID"]
		/* we make a symlink fropm *.bai to *.bam.bai because ATACseqQC makes a few assumptions about our naming convention */
    """
      mkdir -p samtobam_output;
      picard SortSam SO=coordinate INPUT=$samFile OUTPUT=samtobam_output/${output_name}.bam CREATE_INDEX=true\
      > samtobam_output/picard_sortsam.stdout 2> samtobam_output/picard_sortsam.stderr
    	mv samtobam_output/*.bai samtobam_output/${output_name}.bam.bai
		"""
}

process sampleGetMapStats {
  //echo true
  /* our usage of flagstat may soon be superseded by qualimap */
  publishDir {sampleInfo["baseDirOut"]}
  input:
    set val(sampleInfo), file(bamFile) from mappedBamOut_getStat
  output:
    file 'flagstat_output/*'
		file 'qualimap_output*/*'
    file 'flagstat_output/*_flagstat.txt' into flagStatOut
  script:
    output_flagstat_name = sampleInfo["ID"] + '_flagstat.txt'
    """
      mkdir -p flagstat_output;
      samtools flagstat $bamFile > flagstat_output/$output_flagstat_name\
      2> samtools_flagstat.stderr
			
			mkdir -p qualimap_output_${sampleInfo["ID"]};
			qualimap bamqc -bam $bamFile -outdir qualimap_output_${sampleInfo["ID"]}
		"""	
}

process sampleFilterMMMR {
  // filter multimappers and mitochondrial reads
  //echo true
  publishDir {sampleInfo["baseDirOut"]}
  input:
    set val(sampleInfo), file(bamFile) from mappedBamOut
  output:
    file "filtering_output/*"
    set val(sampleInfo), file('filtering_output/*_filtered_MMMR.bam') into filteredMMMROut
  script:
    output_bam_name = sampleInfo["ID"] + '_filtered_MMMR.bam'
    """
      mkdir -p filtering_output;
      java -jar $params.jvarkitPath/samjs.jar\
      --samoutputformat ${def_cmd_params["samjs"]["--samoutputformat"]}\
      -o filtering_output/$output_bam_name\
      -e \'${def_cmd_params["samjs"]["-e"]}\'\
      $bamFile > filtering_output/jvarkit_samjs.stdout 2> filtering_output/jvarkit_samjs.stderr
    """
}

process sampleFilterDedup {
  // discard duplicate reads
  publishDir {sampleInfo["baseDirOut"]}
  input:
    set val(sampleInfo), file(bamFile) from filteredMMMROut
  output:
    file 'filtering_output/*'
    set val(sampleInfo), file('filtering_output/*_dedup.bam') into filteredDedupOut_makeTags, \
                          filteredDedupOut_callPeaks_MACS
  script:
    output_bam_name = sampleInfo["ID"] + '_dedup.bam'
    output_metrics_name = sampleInfo["ID"] + '_dedup_metrics'
    """
      mkdir -p filtering_output;
      picard MarkDuplicates\
      INPUT=$bamFile\
      OUTPUT=filtering_output/$output_bam_name\
      METRICS_FILE=filtering_output/$output_metrics_name\
      REMOVE_DUPLICATES=${def_cmd_params["markduplicates"]["REMOVE_DUPLICATES"]}\
      ASSUME_SORTED=${def_cmd_params["markduplicates"]["ASSUME_SORTED"]}\
      VALIDATION_STRINGENCY=${def_cmd_params["markduplicates"]["VALIDATION_STRINGENCY"]}\
      CREATE_INDEX=${def_cmd_params["markduplicates"]["CREATE_INDEX"]}\
      > filtering_output/picard_markduplicates.stdout 2> filtering_output/picard_markduplicates.stderr
    """
}


process sampleMakeTags {
  publishDir {sampleInfo["baseDirOut"]}
  input:
    set val(sampleInfo), file(bamFile) from filteredDedupOut_makeTags
  output:
    set val(sampleInfo), file('maketags_output') into makeTagsOut, makeTagsOut_callPeaks_homer
  script:
    """
    makeTagDirectory\
    maketags_output\
    -genome $params.refGenomeFasta\
    -checkGC $bamFile\
    > maketagdirectory.stdout 2> maketagdirectory.stderr
    """
}

process sampleMakeUCSCTrack {
  //echo true
  publishDir {sampleInfo["baseDirOut"]}
  input:
    set val(sampleInfo), file(tagDir) from makeTagsOut
  output:
    file "makeUCSCtrack_output/*"
    file "makeUCSCtrack_output/*.bedgraph.gz" into makeUCSCTrackOut 
  script:
    output_track_name = sampleInfo["ID"] + '_UCSC.bedgraph'
    """
    mkdir -p makeUCSCtrack_output;
    makeUCSCfile $tagDir -o makeUCSCtrack_output/$output_track_name\
      -fsize ${def_cmd_params["makeucscfile"]["-fsize"]}\
      -name $output_track_name\
      > makeucsctrack.stdout 2> makeucsctrack.stderr
    """
}

process sampleMACS2CallPeaks {
  publishDir {sampleInfo["baseDirOut"]}
  echo true
  input:
    set val(sampleInfo), file(bamFile) from filteredDedupOut_callPeaks_MACS
  output:
    file "macs2_output/*"
    set val(sampleInfo), file("macs2_output/*.xls") into outputMACS2_xlsPeaks
  script:
    """
    mkdir -p macs2_output;
    macs2 callpeak ${getFlagString(def_cmd_params, "macs2", "--nomodel")}\
    -t $bamFile\
    -n ${sampleInfo["ID"]}\
    --keep-dup ${def_cmd_params["macs2"]["--keep-dup"]}\
    --gsize ${def_cmd_params["macs2"]["--gsize"]}\
    --shift ${def_cmd_params["macs2"]["--shift"]}\
    --extsize ${def_cmd_params["macs2"]["--extsize"]}\
    --outdir macs2_output\
    -B ${getFlagString(def_cmd_params, "macs2", "-B")} > macs2_output/macs2.stdout 2> macs2_output/macs2.stderr
    """
}

process sampleHomerCallPeaks {
  publishDir {sampleInfo["baseDirOut"]}
  errorStrategy 'ignore'
  input:
    set val(sampleInfo), file(tagDir) from makeTagsOut_callPeaks_homer
  output:
    file 'homer_findpeaks_output/*'
    file 'homer_findpeaks_output/*.txt' into outputHomerFindpeaks
  script:
    output_name = sampleInfo["ID"] + '_homer_findpeaks.txt'
  """
    findPeaks\
    $tagDir\
    -style ${def_cmd_params["findpeaks"]["-style"]}\
    -o homer_findpeaks_output/$output_name\
    > homer_findpeaks_output/homer_findpeaks.stdout 2> homer_findpeaks_output/homer_findpeaks.stderr
  """
}

process sampleAnnotatePeaks {
  publishDir {sampleInfo["baseDirOut"]}
  input:
    set val(sampleInfo), file(xlsFile) from outputMACS2_xlsPeaks
  output:
    file 'homer_annotatepeaks_output/*'
    set val(sampleInfo), file("homer_annotatepeaks_output/$annotationFile") into outputAnnotatePeaks
  script:
    annotationFile = sampleInfo["ID"] + '_annotation.txt'
    """
    mkdir -p homer_annotatepeaks_output;
    annotatePeaks.pl\
      $xlsFile\
      $params.refGenomeName\
      > homer_annotatepeaks_output/$annotationFile 2> homer_annotatepeaks_output/homer_annotatepeaks_output.stderr
    """
}

process createQCReport {
  publishDir {sampleInfo["baseDirOut"]}
  echo true
  errorStrategy 'ignore'
  input:
    set val(sampleInfo), file(bamFile), file(bamIndex) from mappedBamOut_QC 
  output:
    file 'qc_output/*'
  exec: 
    rmdScriptDir = workflow.scriptFile.getParent().toString()
  shell:
    """	
    #mv !{bamIndex} !{bamFile}.bai
    mkdir -p qc_output;
    !{rmdScriptDir}/qc/create_qc.sh !{bamFile} !{sampleInfo["ID"]} BSgenome.Drerio.UCSC.danRer10 TxDb.Drerio.UCSC.danRer10.refGene !{rmdScriptDir}/qc/qc.rmd qc_output\
    > qc_output/create_qc_report_output.stdout 2> qc_output/create_qc_report_output.stderr
    """
}
