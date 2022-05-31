####download CCLE RRBS raw data from NCBI###
###analysis the RRBS data#####
##run the code/BismarkHg38.sh for each files####
library(readr)
library(valr)
library(ggplot2)
library(ggpubr)
#setwd("/Users/yueyuanzheng/Documents/GitHub/MethylomeScatterPlot")

##prepare the data
###download the SRR files from PRJNA523380 and transform them into fastq.gz file.
###put all fastq.gz files into plotdata/CCLE/data folder.
###for each file, run the runBismarkHg38.sh script.
##get the CpG methylation for each file.

###obtain the fastCpG and slowCpG mean methylation for each sample
inputPath="CCLE/results/aln/"
outputPath="plotdata/"
fileList=dir(inputPath, pattern="_CpG.txt", full.name=T)
blackListFile="meta/hg38_blacklist.bed"
blackListData=read_bed(blackListFile)
fast_CpG=read_bed("meta/zhou_NN_scores_PMD_prone_10ptile.hg38liftover.bed")
slow_CpG=read_bed("meta/zhou_NN_scores_PMD_resistant_10ptile.hg38liftover.bed")
result=as.data.frame(matrix(numeric(0),ncol=3))
for(file in fileList){
  print(basename(file))
  tmpdata=read.table(file, stringsAsFactors=F, header=T, sep="\t")
  tmpdata$start=tmpdata$base
  tmpdata[tmpdata$strand%in%"R",]$start=tmpdata[tmpdata$strand%in%"R",]$base-1
  tmpdata$end=tmpdata$start+1
  tmpdata$start=tmpdata$start-1
  tmpdata$methy=tmpdata$freqC/100
  tmpdata=tmpdata[,c(2,8,9,10,5,4)]
  tmpdata[tmpdata$strand%in%"R",]$strand="-"
  tmpdata[tmpdata$strand%in%"F",]$strand="+"
  tmpdata=tibble(chrom=tmpdata$chr, start=tmpdata$start, end=tmpdata$end, 
                 methylation=tmpdata$methy, coverage=tmpdata$coverage, strand=tmpdata$strand)
  tmpdata=bed_intersect(tmpdata, blackListData,invert=T)
  tmpdata_fast=bed_intersect(tmpdata, fast_CpG)
  tmpdata_fast=tmpdata_fast[tmpdata_fast$.overlap==1, 1:6]
  tmpdata_slow=bed_intersect(tmpdata, slow_CpG)
  tmpdata_slow=tmpdata_slow[tmpdata_slow$.overlap==1, 1:6]
  fast=mean(tmpdata_fast$methylation.x)
  slow=mean(tmpdata_slow$methylation.x) 
  tmpResult=data.frame(sample=gsub("_CpG.txt", "", basename(file)), fastCpG=fast, slowCpG=slow)
  result=rbind(result, tmpResult)
}
write.table(result, "plotdata/CCLE_fast_slow.txt", row.names=F, col.names=F, quote=F, sep="\t")

###the scatter plot####
my_theme=theme_bw()+theme(plot.background=element_blank(), panel.background=element_blank(),panel.grid.minor=element_blank(),
                          panel.grid.major=element_blank(),axis.title=element_text(color="black",size=14),
                          axis.text=element_text(size=12), legend.title = element_text(hjust=0.5, face="bold"))

data=read.table("plotdata/CCLE_fast_slow.txt", header=T, sep="\t")
colnames(data)[1]="sample"
cancerTypeColor=read.table("meta/cancerTypeColor.txt", sep="\t", stringsAsFactors = F, header=T, comment.char = "@")
CCLEsampleInfo=read_tsv("meta/Cell_lines_annotations_change.txt")
CCLEsampleInfo=CCLEsampleInfo[,colnames(CCLEsampleInfo)%in%c("CCLE_ID","depMapID", "Site_Primary",	"type", "tcga_code",
                                                             "Gender", "lineage_molecular_subtype", "lineage_subtype", "lineage_sub_subtype")]
CCLEsampleInfo=as.data.frame(CCLEsampleInfo[,c(1:3,5)])
data=merge(data, CCLEsampleInfo, by.x="sample", by.y="CCLE_ID")
cancerTypeColor=cancerTypeColor[cancerTypeColor$type%in%unique(data$type),]
data$type=factor(data$type, levels=cancerTypeColor$type)
data=data[order(data$type),]
p1 <- ggplot(data, aes(slowCpG, fastCpG))+geom_point(aes(color=type),size=2, alpha=0.8)+coord_equal(ratio=1)+xlab("Slow-loss PMD methylation")+ylab("Fast-loss PMD methylation")
p1 <- p1+scale_color_manual(values = cancerTypeColor$color)+xlim(0, 1)+ylim(0, 1)+my_theme
p1=p1+geom_abline(intercept = 0, slope = 1)
pdf("result/CCLE_allcells.pdf", width=10, height = 6)
print(p1)
dev.off()

part_data=data[data$type%in%c("upper_aerodigestive","esophagus","stomach","colorectal", "pancreas", "bile_duct", "liver"),]
cancerTypeColor2=cancerTypeColor[cancerTypeColor$type%in%unique(part_data$type),]
part_data$type=factor(part_data$type, levels=cancerTypeColor2$type)
part_data=part_data[order(part_data$type),]
p1 <- ggplot(part_data, aes(slowCpG, fastCpG))+geom_point(aes(color=type),size=2, alpha=0.8)+coord_equal(ratio=1)+xlab("Slow-loss PMD methylation")+ylab("Fast-loss PMD methylation")
p1 <- p1+scale_color_manual(values = cancerTypeColor2$color)+xlim(0, 1)+ylim(0, 1)+my_theme
p1=p1+geom_abline(intercept = 0, slope = 1)
pdf("result/CCLE_GI_cells.pdf", width=8, height = 6)
print(p1)
dev.off()

part_data=data[data$type%in%c("B-cell_ALL","T-cell_ALL","AML", "CML", "leukemia_other", "multiple_myeloma", "lymphoma_DLBCL", 
                              "lymphoma_Burkitt", "lymphoma_Hodgkin", "lymphoma_other"),]
cancerTypeColor2=cancerTypeColor[cancerTypeColor$type%in%unique(part_data$type),]
part_data$type=factor(part_data$type, levels=cancerTypeColor2$type)
part_data=part_data[order(part_data$type),]
p1 <- ggplot(part_data, aes(slowCpG, fastCpG))+geom_point(aes(color=type),size=2, alpha=0.8)+coord_equal(ratio=1)+xlab("Slow-loss PMD methylation")+ylab("Fast-loss PMD methylation")
p1 <- p1+scale_color_manual(values = cancerTypeColor2$color)+xlim(0, 1)+ylim(0, 1)+my_theme
p1=p1+geom_abline(intercept = 0, slope = 1)
pdf("result/CCLE_blood_cells.pdf", width=8, height = 6)
print(p1)
dev.off()

