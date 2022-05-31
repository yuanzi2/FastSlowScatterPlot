#!/bin/sh

fileIndex=$1

thread=8
inputFilePath=CCLE/data/
trimFilePath=CCLE/trimData/
outputPath=CCLE/results/mapping/

logFile=`basename ${inputFile} .fastq.gz`
trim_galore --gzip --length 20 --path_to_cutadapt /opt/Anaconda3/bin/cutadapt --illumina --rrbs ${inputFilePath}/${fileIndex}.fastq.gz -o ${trimFilePath} 2> ${logPath}/${logFile}.trim_galore.log
bismark --parallel ${thread} --genome /data1/reference/hg38/indexed_bismark/ --path_to_bowtie2 /opt/Anaconda3/bin --gzip -o ${outputPath} ${trimFilePath}/${fileIndex}_trimmed.fq.gz
samtools sort -@ ${thread} -O bam -o ${outputPath}/${fileIndex}.sort.bam ${outputPath}/${fileIndex}_trimmed_bismark_bt2.bam
samtools index ${outputPath}/${fileIndex}.sort.bam
Rscript runBismarkAln.R ${fileIndex}
