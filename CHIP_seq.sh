#1. Install the software
conda create -n chipseq
source activate chipseq
conda install trim-galore
conda install sra-tools
conda install samtools
conda install deeptools
conda install homer
conda install meme
conda install macs2
conda install bowtie
conda install bowtie2
conda install gatk
conda install picard
conda install igv



#2. Download data
#2.1 Reference article: RYBP and Cbx7 define specific biological functions of polycomb complexes in mouse embryonic stem cells
#https://www.ncbi.nlm.nih.gov/pubmed/23273917
#RYBP and Cbx7 are both components of Polycomb repressive complex 1 (PRC1):
#Data are in: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42466
#2.2 Download the index of mouse reference genome
#The commonly used mm10 is used here, and the index size is 3.2GB. It is not recommended to download the genome construction by yourself.
wget -4 -q ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/mm10.zip
#2.3 Comment file download
#(1) First enter the Table Browser of UCSC: https://genome.ucsc.edu/cgi-bin/hgTables. (2) Set the parameters, Tips: Change the output format to "BED-browser extensible data", the output file must be filled in, if it is left blank, even if the output format is selected as the BED format, it is the output webpage format; file type returned is generally selected gzip compressed, the download is relatively faster. (3) Click Get output to download.



#3. Data analysis

#3.1 get fastq.gz
#sar转fastq
ls *sra |while read id; do fastq-dump --split-3 $id;done
#gzipCompression fastq-fq.gz
time gzip -c test.fq> test.fq.gz
#Multithreaded pigz compression fastq-fq.gz
time pigz -k -p 8 test.fq
#Download the SRA data interface information SraRunTable.txt according to NCBI, and change the SRR number to the sample name
#Replace comma into tab key
sed's/\,/\x09/g' SraRunTable.txt> SraRunTable1.txt
#View the first line (header) of the file, convert it into a column, and add the line number.
head -1 SraRunTable1.txt | tr'\t''\n'|cat -n
#View all information of the table
less -S SraRunTable1.txt
#Intercept the first and second columns of SraRunTable1.txt, change the space of the file to _, delete ", add 1/2 to the same sample name in the second column. Save it as a config file.
cut -f 1,2 SraRunTable1.txt | tr '''_' | sed's/"//g'|perl -alne'{$h{$F[1]}++ if exists $h{$F [1]}; $h{$F[1]}=1 unless exists $h{$F[1]};print "$F[0]\t$F[1]$h{$F[1] }\t$F[2]"}'> config
#Change SRR number to sample name



#3.2 Use trim_galore software for quality control
#3.2.1 Perform fastqc
cd qc
ls ../raw/*gz | xargs fastqc -t 10 -o ./
#3.2.2.1 Use trim_galore to control batch quality of single-ended *gz files and output to the current folder
cd clean
ls ../raw/*gz | while read fq1;
do
nohup trim_galore -q 25 --phred33 --length 35 -e 0.1 --stringency 4 -o ./ $fq1 &
done
#3.2.2.2 Use trim_galore for batch quality control of double-ended *gz files, and output to the current folder
cd raw
ls ./*_1.fastq.gz> ./1
ls ./*_2.fastq.gz> ./2
paste ./1 ./2> ./config
#New loop script qc.sh
cat> ./qc.sh
#Write the loop in the qc.sh script:
source chipseq
dir='/home/minminli/chipseq/clean'
cat /home/minminli/chipseq/raw/config |while read id
do
           arr=($id)
           fq1=${arr[0]}
           fq2=${arr[1]}
nohup trim_galore -q 25 --phred33 --length 36 -e 0.1 --stringency 3 --paired -o $dir $fq1 $fq2 &
done
source deactivate chipseq
#Run script
bash /home/minminli/chipseq/raw/qc.sh config
#Run in the background, the results are saved in /home/minminli/chipseq/clean.
#View background process
top
#3.2.3 Perform fastqc
cd multiqc
ls ../clean/*gz | xargs fastqc -t 10 -o ./
#Use multiqc to view all the results of fastqc processing in the multiqc file
conda activate python354
multiqc ./ -o ./
#or
multiqc *fastqc.zip --pdf



#3.3 Align
#3.3.1 Compare and generate bam file
cd align
cat> align.sh
#Write the following code into the align.sh file/You can also paste the code directly to run (-m 10M can be omitted, because the memory is small and cannot be sorted, so -m 10M is added)
bowtie2='/home/minminli/software/bowtie2/bowtie2-2.3.5.1-linux-x86_64/bowtie2'
bowtie2_index="/home/minminli/chip-seq/raw/mm10/mm10"
ls ../clean/*gz |while read id;
do
file=$(basename $id)
sample=${file%%.*}
echo $file $sample
$bowtie2 -p 2 -x $bowtie2_index -U $id | samtools sort -m 10M -O bam -@ 5 -o-> ${sample}.bam
done
#Run align.sh file
nohup bash align.sh &

#3.3.2 Combine bam files (optional)
#Applicable situation: Because a sample is divided into multiple lanes for sequencing, it is necessary to merge bam when performing peaks calling. Merge control_1.bam and control_2.bam into control.merge.bam file
samtools merge control.merge.bam control_1.bam control_2.bam
#Or batch processing
#sed's/_[0-9]_trimmed.bam//g' is to remove _[0-9]_trimmed.bam in the original file name
ls *.bam|sed's/_[0-9]_trimmed.bam//g' |sort -u |while read id;do samtools merge ./$id.merge.bam $id*.bam ;done

#3.3.3 Indexing of bam files + QC
cd align
#Create the index of the bam file to get the .bai index file. If samtools reports an error, reinstall the lower version conda install -c bioconda samtools openssl=1.0.
ls *.bam |xargs -i samtools index {}
#After getting the index file, perform QC on the bam file to get the statistical result of .stat
nohup ls *.bam | while read id ;do (samtools flagstat $id> $(basename $id ".bam").stat);done &
#View comparison success rate
grep% *.stat

#3.3.4 Remove PCR duplication from bam file + build index + QC
(1) Method 1: Use samtools software
#Remove PCR duplication from bam file*.bam to get .rmdup.bam file
ls *.bam | while read id ;do (nohup samtools rmdup $id $(basename $id ".bam").rmdup.bam & );done
#Create the index of the .rmdup.bam file to get the .bai index file
ls *.rmdup.bam |xargs -i samtools index {}
#After getting the index file, perform QC on the rmdup.bam file and get the .rmdup.stat statistical result
nohup ls *.rmdup.bam | while read id ;do (samtools flagstat $id> $(basename $id ".rmdup.bam").stat);done &
#View comparison success rate
grep% *.rmdup.stat

(2) Method 2: Use picard software
#Remove PCR duplication, Xmx4g: Set the maximum available memory of JVM to 4G.
picard -Xmx4g MarkDuplicates VALIDATION_STRINGENCY=LENIENT I=Cbx7_trimmed.bam O=Cbx7_trimmed.markdup.bam M=Cbx7_trimmed.markdup_metrics.txt
#Use picard BuildBamIndex alone to build index
picard BuildBamIndex -Xmx4g I=Cbx7_trimmed.markdup.bam
#Or remove PCR duplication + build index
picard -Xmx4g MarkDuplicates VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true I=Cbx7_trimmed.bam O=Cbx7_trimmed.markdup.bam M=Cbx7_trimmed.markdup_metrics.txt
#After getting the index file, perform QC on the markdup.bam file to get the .markdup.stat statistical result
samtools flagstat Cbx7_trimmed.markdup.bam> Cbx7_trimmed.markdup.stat
#Or batch processing
nohup ls *.markdup.bam | while read id ;do (samtools flagstat $id> $(basename $id ".markdup.bam").stat);done &
#View comparison success rate
grep% *.markdup.stat



#3.4 Use macs2 to find peaks
#3.4.1 Run macs2
cd peaks
#Basic parameters: -t: the output result of the experimental group; -c: the output result of the control group; -f: -t and -c provide the file format, which can be "ELAND", "BED", "ELANDMULTI", " ELANDEXPORT", "ELANDMULTIPET" (for pair-end tags), "SAM", "BAM", "BOWTIE", "BAMPE" "BEDPE" any one. If this option is not provided, it is the automatic detection option. Paired-end sequencing uses BAMPE. For single-ended sequencing, no parameters are required. The default is auto recognition, but BAMPE cannot be recognized. -g: Genome size. The hs, mm, ce, dm options are provided by default. If they are not included, such as Arabidopsis, you need to provide them yourself. -n: the prefix name of the output file; -B: will save more information in the bedGraph file, such as fragment pileup, control lambda, -log10pvalue and -log10qvalue scores; -q: q value, which is the minimum PDR threshold, The default is 0.05. The q value is calculated using BH based on the p value, which is the result of multiple trials after correction. -p: This is the p value. After specifying the p value, MACS2 will not use the q value. -m: It is related to MFOLD, and MFOLD is related to the MACS pre-built model. The default is 5:50. MACS will first look for more than 100 peak areas to build the model. Generally, there is no need to change it.
#Generate _model.r, _peaks.narrowPeak, _peaks.xls, _summits.bed, _pileup.bdg, lambda.bdg files
macs2 callpeak -c IgG_old_trimmed.markdup.bam -t Suz12_trimmed.markdup.bam -B -f BAM -g mm -n suz12 -q 0.00001
#Or use the following script to batch process
cd align
ls *markdup.bam |cut -d"." -f 1 |while read id;
do
    if [! -s ${id}_summits.bed ];
    then
echo $id
nohup macs2 callpeak -c IgG_old_trimmed.markdup.bam -t $id.markdup.bam -f BAM -B -g mm -n $id --outdir ../peaks 2> $id.log &
    fi
done
#View background process
top



#3.5 Count the distribution of gene characteristics of peaks in the whole genome
#3.5.1 Roughly check the number of peaks
#See how many peaks are found for each sample (bed format file, which stores the location information of each peak)
wc -l *.bed

#3.5.2 Deeptool result visualization
#(1)BAM to BW file
# Cbx7_trimmed.markdup.bam is the BAM file obtained from the previous comparison and indexed
bamCoverage -e 170 -bs 10 -b Cbx7_trimmed.markdup.bam -o Cbx7_trimmed.markdup.bw
#Or use the following script to batch process
ls *markdup.bam |while read id;do
nohup bamCoverage --normalizeUsing CPM -b $id -o ${id%%.*}.bw &
done

#(2)igv visualization
#First loading the reference genome mm10, then loading the control and processed .bw and .bed files

#(3)Check the signal strength near TSS: reference-point mode
cd tss
###(3.1) Single sample processing: use the following script (-p is thread; tss is near 2.5k)
computeMatrix reference-point --referencePoint TSS -p 2 \
-b 2500 -a 2500 \
-R /home/minminli/chip-seq/raw/reference.bed \
-S /home/minminli/chip-seq/peaks/Cbx7_trimmed.bw \
--skipZeros -o Cbx7_trimmed_TSS.gz \
--outFileSortedRegions Cbx7_trimmed_genes.bed
##plotHeatmap and plotProfile will both use the output of computeMatrix
plotHeatmap -m Cbx7_trimmed_TSS.gz -out Cbx7_Heatmap.png
plotHeatmap -m Cbx7_trimmed_TSS.gz -out Cbx7_Heatmap.pdf --plotFileFormat pdf --dpi 720
plotProfile -m Cbx7_trimmed_TSS.gz -out Cbx7_Profile.png
plotProfile -m Cbx7_trimmed_TSS.gz -out Cbx7_Profile.pdf --plotFileFormat pdf --perGroup --dpi 720

###(3.2) or batch processing: use the following script (-p is thread)
cd tss
cat> tss_2.5.sh
#Write the following code into the tss_2.5.sh script
bed=/home/minminli/chip-seq/raw/reference.bed
for id in /home/minminli/chip-seq/peaks/*bw;
do
echo $id
file=$(basename $id)
sample=${file%%.*}
echo $sample

computeMatrix reference-point --referencePoint TSS -p 2 \
-b 2500 -a 2500 \
-R $bed \
-S $id \
--skipZeros -o ${sample}_TSS.gz \
--outFileSortedRegions ${sample}_genes.bed

Both #plotHeatmap and plotProfile will use the output of computeMatrix
plotHeatmap -m ${sample}_TSS.gz -out ${sample}_Heatmap_2.5k.png
plotHeatmap -m ${sample}_TSS.gz -out ${sample}_Heatmap_2.5k.pdf --plotFileFormat pdf --dpi 720
plotProfile -m ${sample}_TSS.gz -out ${sample}_Profile_2.5k.png
plotProfile -m ${sample}_TSS.gz -out ${sample}_Profile_2.5k.pdf --plotFileFormat pdf --perGroup --dpi 720
done
#Use command batch submission:
nohup bash tss_2.5.sh 1>2.5k.log &
#View process
cat 2.5k.log

###(3.3) Finally, integrate all the bam files of chipseq, draw the profile and heatmap integration map near the TSS of the gene
computeMatrix reference-point -p 2 --referencePoint TSS -b 2500 -a 2500 -S ../*bw -R /home/minminli/chip-seq/raw/reference.bed --skipZeros -o all.merge.gz
plotHeatmap -m all.merge.gz -out all.merge_Heatmap.png
plotHeatmap --dpi 720 -m all.merge.gz -out all.merge_Heatmap.pdf --plotFileFormat pdf
plotProfile -m all.merge.gz -out all.merge_profile.png
plotProfile --dpi 720 -m all.merge.gz -out all.merge_profile.pdf --plotFileFormat pdf --perGroup


#(4) View the signal strength near genebody: scale-regions mode
###(4.1) Single sample processing: use the following script
computeMatrix scale-regions \ # select mode
       -b 3000 -a 5000 \ # The area of ​​interest, -b upstream, -a downstream
       -R ~/reference/gtf/TAIR10/TAIR10_GFF3_genes.bed \
       -S 03-read-coverage/ap2_chip_rep1_1.bw \
       --skipZeros \
       --outFileNameMatrix 03-read-coverage/matrix1_ap2_chip_rep1_1_scaled.tab \ # output as a file for plotHeatmap, plotProfile
       --outFileSortedRegions 03-read-coverage/regions1_ap2_chip_re1_1_genes.bed

###(4.2) or batch processing: use the following script (-p is thread)
cd tss
cat> genebody.sh
#Write the following code into the genebody.sh script
bed=/public/annotation/CHIPseq/mm10/ucsc.refseq.bed
for id in /home/llwu/jmzeng/align/*.rmdup.bw;
do
echo $id
file=$(basename $id)
sample=${file%%.*}
echo $sample

computeMatrix scale-regions -p 5 \
-b 3000 -a 3000 \
-R $bed \
-S $id \
--regionBodyLength 15000 --skipZeros \
-o matrix1_${sample}_genebody.gz \
--outFileSortedRegions regions1_${sample}_genebody.bed

##both plotHeatmap and plotProfile will use the output from computeMatrix
plotHeatmap -m matrix1_${sample}_genebody.gz -out ${sample}_Heatmap_genebody.png
plotHeatmap -m matrix1_${sample}_genebody.gz -out ${sample}_Heatmap_genebody.pdf --plotFileFormat pdf --dpi 720
plotProfile -m matrix1_${sample}_genebody.gz -out ${sample}_Profile_genebody.png
plotProfile -m matrix1_${sample}_genebody.gz -out ${sample}_Profile_genebody.pdf --plotFileFormat pdf --perGroup --dpi 720
done
#Use command batch submission:
nohup bash genebody.sh 1>genebody.log &
#View process
cat genebody.log

#(5) Use the R package to annotate the peaks file found
#Install R package
BiocInstaller::biocLite(c('airway','DESeq2','edgeR','limma'))
BiocInstaller::biocLite(c('ChIPpeakAnno','ChIPseeker'))
BiocInstaller::biocLite('TxDb.Hsapiens.UCSC.hg19.knownGene',ask=F,suppressUpdates=T)
BiocInstaller::biocLite('TxDb.Hsapiens.UCSC.hg38.knownGene',ask=F,suppressUpdates=T)
BiocInstaller::biocLite('TxDb.Mmusculus.UCSC.mm10.knownGene',ask=F,suppressUpdates=T)

#RCode
###(5.1).download packages###
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
#Install clusterProfiler package
install.packages("udunits2",configure.args ='--with-udunits2-include=/usr/include/udunits2')
BiocManager::install('clusterProfiler')
#or
devtools::install_github("guangchuangyu/clusterProfiler")

###(5.2).Read in bed file###
rm(list = ls())
getwd()
setwd("C:/Users/99039/Desktop/bed")
bedPeaksFile = 'Cbx7_trimmed_summits.bed'; 
bedPeaksFile

###(5.3). loading packages###
require(ChIPseeker)
require(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
require(org.Mm.eg.db)
require(clusterProfiler)

###(5.4).Visualization###
peak <- readPeakFile( bedPeaksFile)
# View chromosome information in Peak file
seqlevels(peak)
# Filter out poor quality chromosomes with _
keepChr = !grepl('_',seqlevels(peak))
seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
# See how many peaks there are
cat(paste0('there are',length(peak),' peaks for this data' ))
# Use annotatPeak for annotation
peakAnno <- annotatePeak(peak, tssRegion = c(-3000, 3000),
                         TxDb = txdb, annoDb = "org.Mm.eg.db")
# Convert to data.frame format file for easy viewing and subsequent operations
peakAnno_df <- as.data.frame(peakAnno)
# Get sample name
sampleName = basename(strsplit(bedPeaksFile,'\\.')[[1]][1])
print(sampleName)

###(5.4.1). Visualize peak###
pdf("plotAnnoPie.pdf")
plotAnnoPie(peakAnno)
dev.off()

pdf("plotAnnoBare.pdf")
plotAnnoBar(peakAnno)
dev.off()

###(5.4.2). Distribution of peak on chromosome###
pdf("peak-chr.pdf")
covplot(peak)
dev.off()

covplot(peak, chr = c("chr1", "chr2")) #specify chromosome

###(5.4.3). Peak associated gene annotation-get annotation file###
peakAnno <- annotatePeak(peak, tssRegion = c(-3000, 3000), TxDb = txdb, annoDb = "org.Mm.eg.db")
write.table(as.data.frame(peakAnno), "peak.annotation.tsv", sep="\t", row.names = F, quote = F)

#(6)homer software to find motif
###(6.1) Then use the configuration script attached to homer to download the database###
perl ~/miniconda3/envs/chipseq/share/homer-4.9.1-5/configureHomer.pl -install mm10
ls -lh ~/miniconda3/envs/chipseq/share/homer-4.9.1-5/data/
##After the download is successful, there will be more ~/miniconda3/envs/chipseq/share/homer-4.9.1-5/data/genomes/mm9/ folder, 4.9G in total
##This folder depends on where you installed the homer software.

## Or use the following code to install:
cd ~/biosoft
mkdir homer && cd homer
wget http://homer.salk.edu/homer/configureHomer.pl
perl configureHomer.pl -install
perl configureHomer.pl -install mm10

###(6.2) The homer software finds motif and integrates two methods, including database-dependent query, and de novo inference, both of which read the peaks file in bed format obtained from upstream analysis of ChIP-seq data ###
cd motif
#(6.2.1) Find motif for a single sample
#If it is the peaks record file found by MACS, extract the corresponding column to HOMER as the input file:
awk'{print $4"\t"$1"\t"$2"\t"$3"\t+"}' Cbx7_trimmed_summits.bed >Cbx7_trimmed_summits_homer.bed
findMotifsGenome.pl ./Cbx7_trimmed_summits_homer.bed mm10 ./output -size 200 -len 8,10,12
#(6.2.2) Batch processing samples to find motif
#Run the following script (note: use the full path here)
for id in /home/minminli/chipseq/peaks/*.bed;
do
echo $id
file=$(basename $id)
sample=${file%%.*}
echo $sample
awk'{print $4"\t"$1"\t"$2"\t"$3"\t+"}' $id >homer_peaks.tmp
findMotifsGenome.pl homer_peaks.tmp mm10 ${sample}_motifDir -len 8,10,12
annotatePeaks.pl homer_peaks.tmp mm10 1>${sample}.peakAnn.xls 2>${sample}.annLog.txt
done
#(6.2.3) Find differences between two/multiple samples motif
awk'{print $4"\t"$1"\t"$2"\t"$3"\t+"}' Cbx7_trimmed_summits.bed >Cbx7_trimmed_summits_homer.bed
awk'{print $4"\t"$1"\t"$2"\t"$3"\t+"}' IgG_trimmed_summits.bed >IgG_trimmed_summits_homer.bed
awk'{print $4"\t"$1"\t"$2"\t"$3"\t+"}' Ring1B_trimmed_summits.bed >Ring1B_trimmed_summits_homer.bed
awk'{print $4"\t"$1"\t"$2"\t"$3"\t+"}' RYBP_trimmed_summits.bed >RYBP_trimmed_summits_homer.bed
awk'{print $4"\t"$1"\t"$2"\t"$3"\t+"}' Suz12_trimmed_summits.bed >Suz12_trimmed_summits_homer.bed
#Two
findMotifsGenome.pl ./Cbx7_trimmed_summits_homer.bed mm10 ./output -size 200 -len 8,10,12 -bg IgG_trimmed_summits_homer.bed
#Multiple
findMotifsGenome.pl ./Cbx7_trimmed_summits_homer.bed mm10 ./output -size 200 -len 8,10,12 -bg ./Ring1B_trimmed_summits_homer.bed -bg ./RYBP_trimmed_summits_homer.bed -bg ./Suzmit12_trimmed_sumbed


#(7)meme software to find motif
###(7.1) Obtain fasta sequence by the coordinates of peaks in bed format
BiocManager::install('BSgenome.Mmusculus.UCSC.mm10')
BiocManager::install('regioneR')
library(BSgenome.Mmusculus.UCSC.mm10)
library(regioneR)
setwd("/home/minminli/chip-seq/peaks/dup/setwd")
bedPeaksFile='Cbx7_trimmed_summits.bed';
sampleName=strsplit(bedPeaksFile,'\\.')[[1]][1]
peak <- toGRanges(bedPeaksFile, format="BED", header=FALSE)
keepChr = !grepl('_',seqlevels(peak))
seqlevels(peak, force=TRUE) <- seqlevels(peak)[keepChr]

seq <- getAllPeakSequence(peak, upstream=20, downstream=20, genome=Hsapiens)

write2FASTA(seq, paste0(sampleName,'.fa'))
###(7.2)MEME, link: http://meme-suite.org/, submit the .fa file and download the result.


#(8) Difference peak analysis


