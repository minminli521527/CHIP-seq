#1.download packages
options("repos"=c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN"))
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ChIPseeker")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
BiocManager::install("ReactomePA")
BiocManager::install("DOSE")
BiocManager::install('ChIPpeakAnno')
BiocManager::install('TxDb.Mmusculus.UCSC.mm10.knownGene',ask=F,suppressUpdates=T)
#安装clusterProfiler包
install.packages("udunits2",configure.args = '--with-udunits2-include=/usr/include/udunits2')
BiocManager::install('clusterProfiler')
#或者
devtools::install_github("guangchuangyu/clusterProfiler")


#2.读入bed文件
rm(list = ls())
getwd()
setwd("C:/Users/99039/Desktop/bed")
bedPeaksFile = 'Cbx7_trimmed_summits.bed'; 
bedPeaksFile


#3. loading packages
require(ChIPseeker)
require(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
require(org.Mm.eg.db)
require(clusterProfiler)


#4.可视化
peak <- readPeakFile( bedPeaksFile ) 
# 查看 Peak 文件中染色体信息
seqlevels(peak)
# 过滤掉带有_的质量不好的染色体
keepChr= !grepl('_',seqlevels(peak))      
seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
# 查看有多少个peak
cat(paste0('there are ',length(peak),' peaks for this data' ))
# 使用 annotatPeak 进行注释
peakAnno <- annotatePeak(peak, tssRegion = c(-3000, 3000), 
                         TxDb = txdb, annoDb = "org.Mm.eg.db") 
# 转变成 data.frame 格式文件，方便查看与后续操作
peakAnno_df <- as.data.frame(peakAnno)
# 获取样本名
sampleName = basename(strsplit(bedPeaksFile,'\\.')[[1]][1])
print(sampleName)

#4.1 对peak进行可视化
pdf("plotAnnoPie.pdf")
plotAnnoPie(peakAnno)
dev.off()

pdf("plotAnnoBare.pdf")
plotAnnoBar(peakAnno)
dev.off()

#4.2 peak在染色体上的分布
pdf("peak-chr.pdf")
covplot(peak)
dev.off()

covplot(peak, chr = c("chr1", "chr2"))   #指定染色体

#4.3 peak关联基因注释-得到注释文件
peakAnno <- annotatePeak(peak, tssRegion = c(-3000, 3000), TxDb = txdb, annoDb = "org.Mm.eg.db")
write.table(as.data.frame(peakAnno), "peak.annotation.tsv", sep="\t", row.names = F, quote = F)



