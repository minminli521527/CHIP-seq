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
#��װclusterProfiler��
install.packages("udunits2",configure.args = '--with-udunits2-include=/usr/include/udunits2')
BiocManager::install('clusterProfiler')
#����
devtools::install_github("guangchuangyu/clusterProfiler")


#2.����bed�ļ�
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


#4.���ӻ�
peak <- readPeakFile( bedPeaksFile ) 
# �鿴 Peak �ļ���Ⱦɫ����Ϣ
seqlevels(peak)
# ���˵�����_���������õ�Ⱦɫ��
keepChr= !grepl('_',seqlevels(peak))      
seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
# �鿴�ж��ٸ�peak
cat(paste0('there are ',length(peak),' peaks for this data' ))
# ʹ�� annotatPeak ����ע��
peakAnno <- annotatePeak(peak, tssRegion = c(-3000, 3000), 
                         TxDb = txdb, annoDb = "org.Mm.eg.db") 
# ת��� data.frame ��ʽ�ļ�������鿴���������
peakAnno_df <- as.data.frame(peakAnno)
# ��ȡ������
sampleName = basename(strsplit(bedPeaksFile,'\\.')[[1]][1])
print(sampleName)

#4.1 ��peak���п��ӻ�
pdf("plotAnnoPie.pdf")
plotAnnoPie(peakAnno)
dev.off()

pdf("plotAnnoBare.pdf")
plotAnnoBar(peakAnno)
dev.off()

#4.2 peak��Ⱦɫ���ϵķֲ�
pdf("peak-chr.pdf")
covplot(peak)
dev.off()

covplot(peak, chr = c("chr1", "chr2"))   #ָ��Ⱦɫ��

#4.3 peak��������ע��-�õ�ע���ļ�
peakAnno <- annotatePeak(peak, tssRegion = c(-3000, 3000), TxDb = txdb, annoDb = "org.Mm.eg.db")
write.table(as.data.frame(peakAnno), "peak.annotation.tsv", sep="\t", row.names = F, quote = F)


