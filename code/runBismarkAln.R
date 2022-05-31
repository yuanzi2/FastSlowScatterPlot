library("methylKit")
args=commandArgs(T)

inputIndex=args[1]
inputPath="CCLE/results/mapping/"
outputPath="CCLE/results/aln/"

my.methRaw=processBismarkAln(location = paste0(inputPath, inputIndex, ".sort.bam"),
                         sample.id=inputIndex, assembly="hg38",
                         read.context="CpG", mincov = 5, minqual = 20, save.folder=outputPath)
