#!/bin/bash
# check for presence of dependencies

UNMET_DEP=0

check_dep(){
	if [ -z $(command -v $1) ]; then
		echo "$1 not found"; return 1;
	else
		echo "$1 found"; return 0;
	fi
}

check_dep fastqc
check_dep cutadapt
check_dep bowtie2
check_dep picard
check_dep samtools
check_dep homer
