setwd("Calculate_prone_resistant_Plot/")
library("parallel")
library(valr)
library(readr)
library(ggplot2)

###
fileList=dir("WangdingWGBS/", pattern = ".bed", full.names = T)
fast_CpG=read_bed("meta/zhou_NN_scores_PMD_prone_10ptile.hg19.bedgraph")
slow_CpG=read_bed("meta/zhou_NN_scores_PMD_resistant_10ptile.hg19.bedgraph")
result=mclapply(fileList, function(file){
   print(gsub(".bed.gz", "", basename(file)))
   tmpdata=read_bed(file)
   tmpdata_fast=bed_intersect(fast_CpG, tmpdata)
   tmpdata_fast=unique(tmpdata_fast[,c(1:4,7)])
   colnames(tmpdata_fast)=c("chr", "start", "end", "NNscore", "methylation")
   tmpdata_slow=bed_intersect(slow_CpG, tmpdata)
   tmpdata_slow=unique(tmpdata_slow[,c(1:4,7)])
   colnames(tmpdata_slow)=c("chr", "start", "end", "NNscore", "methylation")
   tmpRsult=data.frame(fileName=gsub(".bed.gz", "", basename(file)), fastCpG=mean(tmpdata_fast$methylation), slowCpG=mean(tmpdata_slow$methylation))
   return(tmpRsult)
  }, mc.cores = 20)
result=do.call(rbind, result)
write.table(plotdata, file="result/WangDing_prone_resistant.txt", row.names=F, col.names=T, sep="\t", quote=F)


my_theme=theme_bw()+theme(plot.background=element_blank(), panel.background=element_blank(),panel.grid.minor=element_blank(),
                          panel.grid.major=element_blank(),axis.title=element_text(color="black",size=12),
                          axis.text=element_text(size=12), legend.title = element_text(hjust=0.5, face="bold"))
bluePrintSamples=read.table("meta/BLUEPRINT_blood_sampleInfo.txt", sep="\t", header=T, stringsAsFactors = F)
bluePrintSamples=bluePrintSamples[order(bluePrintSamples$Type),]
bluePrintSamples_normal=bluePrintSamples[bluePrintSamples$Type%in%c("Lymphoid", "Myeloid", "Others"),]
bluePrintSamples_normal$Type=factor(bluePrintSamples_normal$Type, levels=unique(bluePrintSamples_normal$Type))
bluePrintSamples_normal$Subtype=factor(bluePrintSamples_normal$Subtype,
                                       levels=c("BCell", "memoryBCell", "TCell", "memoryTCell", "NKCell", "ProgenitorCell", "Plasma",
                                                "DendriticCell","Eosinophil", "Erythroblast", "Neutrophil", "Monocyte", "Macrophage", "Megakaryocyte", "Osteoclast"))
bluePrintSamples_tumor=bluePrintSamples[!bluePrintSamples$Type%in%c("Myeloid","Lymphoid","Others"),]
bluePrintSamples_tumor=bluePrintSamples_tumor[order(bluePrintSamples_tumor$Type),]
bluePrintSamples_tumor$Type=factor(bluePrintSamples_tumor$Type, levels=unique(bluePrintSamples_tumor$Type))
bluePrintSamples_tumor$Subtype=factor(bluePrintSamples_tumor$Subtype, levels=c("AcuteLymphocyticLeukemia", "TcellProlymphocyticLeukemia",
                                                                               "ChronicLymphocyticLeukemia", "MantleCellLymphoma", "MultipleMyeloma", "AcuteMyeloidLeukemia"))
data=read.table("result/WangDing_prone_resistant.txt", header=T, sep="\t", stringsAsFactors = F)
bluePrintSamples_tumor=merge(data, bluePrintSamples_tumor, by.x="fileName", by.y="Sample")
bluePrintSamples_tumor=bluePrintSamples_tumor[order(bluePrintSamples_tumor$Subtype),]
p1 <- ggplot(bluePrintSamples_tumor, aes(resistant, prone, shape=Type))+geom_point(aes(color=Subtype),size=4, alpha=0.8)+coord_equal(ratio=1)
p1=p1+xlab("mean PMD methylation(resistant)")+ylab("mean PMD methylation(prone)")
p1 <- p1+scale_shape_manual(values=c(15, 17))+scale_color_manual(values = as.character(unique(bluePrintSamples_tumor$Color2)))+xlim(0, 1)+ylim(0, 1)+my_theme
p1=p1+geom_abline(intercept = 0, slope = 1)
pdf("result/BLUEPRINT_Tumor_resistant_prone.pdf", width=7, height=5)
print(p1)
dev.off()

bluePrintSamples_normal=merge(data, bluePrintSamples_normal, by.x="fileName", by.y="Sample")
bluePrintSamples_normal=bluePrintSamples_normal[order(bluePrintSamples_normal$Subtype),]
p1 <- ggplot(bluePrintSamples_normal, aes(resistant, prone, shape=Type))+geom_point(aes(color=Subtype),size=4, alpha=0.8)+coord_equal(ratio=1)
p1=p1+xlab("mean PMD methylation(resistant)")+ylab("mean PMD methylation(prone)")
p1 <- p1+scale_shape_manual(values=c(15, 17, 18))+scale_color_manual(values = as.character(unique(bluePrintSamples_normal$Color2)))+xlim(0, 1)+ylim(0, 1)+my_theme
p1=p1+geom_abline(intercept = 0, slope = 1)
pdf("result/BLUEPRINT_Normal_resistant_prone.pdf", width=7, height=5)
print(p1)
dev.off()

######
bluePrintCelllineSamples=read.table("meta/BLUEPRINT_cellline_sampleInfo.txt", sep="\t", header=T,stringsAsFactors = F)
bluePrintCelllineSamples=merge(data, bluePrintCelllineSamples, by.x="fileName", by.y="Sample")
bluePrintCelllineSamples$subtype=factor(bluePrintCelllineSamples$subtype, levels=c("4star", "H1", "ESC.derived.endoderm","ESC.derived.ectoderm","ESC.derived.mesendoderm","ESC.derived.mesoderm","ESC.derived.trophoblast", "H9", "iPSC","HUES64",
                                           "Neurosphere","Fibroblast","ICF","ADS.derived", "ADS", "Breast", "BCell", "LiverCancer", "BreastCancer"))
bluePrintCelllineSamples=bluePrintCelllineSamples[order(bluePrintCelllineSamples$subtype),]
p1 <- ggplot(bluePrintCelllineSamples, aes(resistant, prone))+geom_point(aes(color=subtype),size=4, alpha=0.8, shape=15)+coord_equal(ratio=1)+xlab("mean PMD methylation(resistant)")+ylab("mean PMD methylation(prone)")
p1 <- p1+scale_color_manual(values = unique(bluePrintCelllineSamples$color))+xlim(0, 1)+ylim(0, 1)+my_theme
p1=p1+geom_abline(intercept = 0, slope = 1)
pdf("result/BP_Cellline_resistant_prone.pdf", width=7, height=5)
print(p1)
dev.off()
