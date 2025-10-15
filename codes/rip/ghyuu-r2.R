#write by Sun Haiayng at 2024.07.09
#1、搭环境
{
  rm(list=ls())
  for (i in c("ReactomePA","ggplot2","stringr","ChIPseeker","DESeq2","clusterProfiler","topGO","Rgraphviz","pathview","org.Mm.eg.db","TxDb.Mmusculus.UCSC.mm10.knownGene")) {library(i, character.only = T)};rm(i)
}
#2、给索引
{
  path = "/ChIP_seq_2/aaSHY/ghyuu/ripseq/20231025/star/pairEND/"
  inx1 = "/Reference/aaSHY/Ensembl99/GRCm38/INDEX/bowtie2/mm10"
  inx2 = "/Reference/aaSHY/Ensembl99/GRCm38/INDEX/STARv279a"
  inx3 = "/Reference/aaSHY/Ensembl99/GRCm38/fa/GRCm38.dna.primary_assembly.fa"
  inx4 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/Mus_musculus.GRCm38.99.gtf"
  inx5 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf"
  inx6 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf.locInd"
  inx7 = "/Reference/aaSHY/Ensembl99/GRCm38/BED/TEsubFamily/MM::ERVL::LTR.bed"
  inx8 = "/Reference/aaSHY/Ensembl99/GRCm38/BED/special/MM_zyw.bed"
  inx9 = "/Reference/aaSHY/Ensembl99/GRCm38/BED/TEsubFamily/MM22::ERVL::LTR.bed"
  inxA = "/Reference/aaSHY/Ensembl99/GRCm38/BED/special/allGene.bed"
  inxB = "/Reference/aaSHY/Ensembl99/GRCm38/BED/special/twoC.bed"
  inxC = "/Reference/aaSHY/zOther/Rn45s/INDEX/STARv279a"
  inxD = "/Reference/aaSHY/zOther/TEconsensus/INDEX/MuERVL/STARv279a"
}
################################################################################
#3、跑命令
{
  for (i in c("src","fq","Log","peak","MM","rRNA","dumpROOM",
              "fisher","deeptools","bam","trim","PLOT","bw/bcv","bw/bcp","count")) {
    if (dir.exists(paths = paste0(path, i))==F) {
      dir.create(path = paste0(path,i), recursive = T)
    }
  }
  shell <- paste0(path,"src/run_a.sh")
  cat("#!/bin/bash\n", file = shell)
  prex <- c("IP","Input")
  cmd_01 <- paste0("trim_galore --phred33 --fastqc --illumina ",
                   "--clip_R1 10 --clip_R2 10 --three_prime_clip_R1 2 --three_prime_clip_R2 2 ",
                   "-j 2 --paired ", path,"fq/",prex,"_1.fq.gz ",path, 
                   "fq/", prex, "_2.fq.gz -o ", path, "trim", "\n")
  lay1 <- paste0(path,"trim/",prex,"_1_val_1.fq.gz")
  lay2 <- paste0(path,"trim/",prex,"_2_val_2.fq.gz")
  cmd_02 <- paste0("STAR --runThreadN 10 --readFilesCommand zcat ",
                   "--outMultimapperOrder Random --outSAMtype BAM SortedByCoordinate ",
                   "--outFileNamePrefix ",path,"bam/",prex," ",
                   "--outReadsUnmapped Fastx ",path,"dumpROOM/",prex,
                   "-mate1.fq ",path,"dumpROOM/",prex,"-mate2.fq ",
                   "--outSAMprimaryFlag AllBestScore --twopassMode Basic ",
                   "--genomeDir ",inx2," --readFilesIn ",lay1," ",lay2,"\n")
  cmd_03 <- paste0("mv ",path,"bam/",prex,"Aligned.sortedByCoord.out.bam ",path,
                   "bam/",prex,".bam","\n")
  cmd_04 <- paste0("samtools index ",path,"bam/",prex,".bam","\n")
  cmd_05 <- paste0("bamCoverage --normalizeUsing CPM -p 10 --minMappingQuality 1 ",
                   "--samFlagExclude 256 -of bigwig -bs 20 -b ",
                   path,"bam/",prex,".bam -o ",
                   path,"bw/bcv/",prex,".bw","\n")
  cmd_06 <- paste0("bamCompare -p 10 --scaleFactorsMethod SES --sampleLength 100000 -b1 ",
                   path,"bam/",prex[1],".bam -b2 ",
                   path,"bam/",prex[2],".bam -o ",
                   path,"bw/bcp/IP_Input.bw","\n")
  cmd_07 <- paste0("computeMatrix scale-regions -p 10 -a 2000 -b 2000 -m 5000 ",
                   "--missingDataAsZero ", "-R ", c(inx7,inx8,inx9,inxA,inxB)," -S ", path,
                   "bw/bcp/IP_Input.bw -o ", path,
                   "deeptools/", basename(c(inx7,inx8,inx9,inxA,inxB)), ".mat.gz", "\n")
  cmd_08 <- paste0("plotHeatmap --heatmapHeight 12 --colorList white,red -m ",path,
                   "deeptools/",basename(c(inx7,inx8,inx9,inxA,inxB)),
                   ".mat.gz -out ",path,"PLOT/",
                   basename(c(inx7,inx8,inx9,inxA,inxB)),".pdf","\n")
  cmd_09 <- paste0("macs2 callpeak -f BAM -g mm -p 0.001 --keep-dup 1 -t ",
                   path,"bam/",prex[1],".bam -c ",
                   path,"bam/",prex[2],".bam -n ghyuuN --outdir ",
                   path,"peak/", "\n")
  cmd_10 <- paste0("macs2 callpeak -f BAM -g mm -p 0.001 --broad --broad-cutoff 0.001 ",
                   "--keep-dup 1 -t ",path,"bam/", prex[1],
                   ".bam -c ", path,"bam/", prex[2],
                   ".bam -n ghyuuB --outdir ", path,"peak/", "\n")
  cmd_11 <- paste0("TEtranscripts -t ",
                   path,"bam/",prex[1],".bam"," -c ",
                   path,"bam/",prex[2],".bam"," ",
                   "--GTF ",inx4," --TE ",inx5," --sortByPos --mode multi ",
                   "--project ghyuu --outdir ",path,"count/","\n")
  cmd_12 <- paste0("TElocal -b ",
                   path,"bam/",prex,".bam ",
                   "--GTF ",inx4," --TE ",inx6," --sortByPos --project ",
                   path, "count/", prex,"\n")
  cmd_13 <- paste0("STAR --runThreadN 10 --readFilesCommand zcat ",
                   "--outMultimapperOrder Random --outSAMtype BAM SortedByCoordinate ",
                   "--outFileNamePrefix ",path,"rRNA/",prex," ",
                   "--outSAMprimaryFlag AllBestScore --twopassMode Basic ",
                   "--genomeDir ",inxC," --readFilesIn ",lay1," ",lay2,"\n")
  cmd_14 <- paste0("mv ",path,"rRNA/",prex,"Aligned.sortedByCoord.out.bam ",path,
                   "rRNA/",prex,".bam","\n")
  cmd_15 <- paste0("samtools index ",path,"rRNA/",prex,".bam","\n")
  cmd_16 <- paste0("bamCoverage --normalizeUsing CPM -p 10 --minMappingQuality 1 ",
                   "--samFlagExclude 256 -of bigwig -bs 20 -b ",path,
                   "rRNA/",prex,".bam -o ",path,
                   "rRNA/",prex,".rRNA.bw\n")
  cmd_17 <- paste0("STAR --runThreadN 10 --readFilesCommand zcat ",
                   "--outMultimapperOrder Random --outSAMtype BAM SortedByCoordinate ",
                   "--outFileNamePrefix ",path,"MM/",prex," ",
                   "--outSAMprimaryFlag AllBestScore --twopassMode Basic ",
                   "--genomeDir ",inxD," --readFilesIn ",lay1," ",lay2,"\n")
  cmd_18 <- paste0("mv ",path,"MM/",prex,"Aligned.sortedByCoord.out.bam ",path,
                   "MM/",prex,".bam","\n")
  cmd_19 <- paste0("samtools index ",path,"MM/",prex,".bam","\n")
  cmd_20 <- paste0("bamCoverage --normalizeUsing CPM -p 10 --minMappingQuality 1 ",
                   "--samFlagExclude 256 -of bigwig -bs 20 -b ",path,
                   "MM/",prex,".bam -o ",path,
                   "MM/",prex,".MM.bw\n")
  for (i in cmd_01) {cat(i, append = T, file = shell)}
  for (i in cmd_02) {cat(i, append = T, file = shell)}
  for (i in cmd_03) {cat(i, append = T, file = shell)}
  for (i in cmd_04) {cat(i, append = T, file = shell)}
  for (i in cmd_05) {cat(i, append = T, file = shell)}
  for (i in cmd_06) {cat(i, append = T, file = shell)}
  for (i in cmd_07) {cat(i, append = T, file = shell)}
  for (i in cmd_08) {cat(i, append = T, file = shell)}
  for (i in cmd_09) {cat(i, append = T, file = shell)}
  for (i in cmd_10) {cat(i, append = T, file = shell)}
  for (i in cmd_11) {cat(i, append = T, file = shell)}
  for (i in cmd_12) {cat(i, append = T, file = shell)}
  for (i in cmd_13) {cat(i, append = T, file = shell)}
  for (i in cmd_14) {cat(i, append = T, file = shell)}
  for (i in cmd_15) {cat(i, append = T, file = shell)}
  for (i in cmd_16) {cat(i, append = T, file = shell)}
  for (i in cmd_17) {cat(i, append = T, file = shell)}
  for (i in cmd_18) {cat(i, append = T, file = shell)}
  for (i in cmd_19) {cat(i, append = T, file = shell)}
  for (i in cmd_20) {cat(i, append = T, file = shell)}
  print(paste0("nohup bash ",shell," >",path,"Log/run_a.log"," 2>&1 &"))
}
################################################################################
#4、计算CPM Excel
{
  cnt1 <- read.table(paste0(path, "count/ghyuu.cntTable"), header = T, sep = "\t", row.names = 1)
  colnames(cnt1) <- c("IP","Input")
  tmp1 <- cnt1[grep(x = rownames(cnt1), pattern=":", invert = F),]#TE
  tmp2 <- cnt1[grep(x = rownames(cnt1), pattern=":", invert = T),]#Gene
  tmp3 <- cnt1#All
  tmp1 <- as.data.frame(apply(tmp1, 2, function(x){x/sum(x)*1000000}))
  tmp2 <- as.data.frame(apply(tmp2, 2, function(x){x/sum(x)*1000000}))
  tmp3 <- as.data.frame(apply(tmp3, 2, function(x){x/sum(x)*1000000}))
  tmp1$theName <- rownames(tmp1)
  tmp2$theName <- rownames(tmp2)
  tmp3$theName <- rownames(tmp3)
  openxlsx::write.xlsx(x = tmp1, file = "/ChIP_seq_2/aaSHY/ghyuu/ripseq/20231025/star/pairEND/count/ghyuu.cpm.onlyTE.xlsx")
  openxlsx::write.xlsx(x = tmp2, file = "/ChIP_seq_2/aaSHY/ghyuu/ripseq/20231025/star/pairEND/count/ghyuu.cpm.onlyGene.xlsx")
  openxlsx::write.xlsx(x = tmp3, file = "/ChIP_seq_2/aaSHY/ghyuu/ripseq/20231025/star/pairEND/count/ghyuu.cpm.allGeneTE.xlsx")
}
################################################################################
#5、CPM点图
{
  cnt1 <- openxlsx::read.xlsx(xlsxFile = "/ChIP_seq_2/aaSHY/ghyuu/ripseq/20231025/star/pairEND/count/ghyuu.cpm.onlyTE.xlsx")
  cnt1$mark <- "other"
  cnt1$mark[grep(x = cnt1$theName, pattern = "MM22")] <- "MM22"
  cnt1$mark[grep(x = cnt1$theName, pattern = "MM")] <- "MM"
  labe <- cnt1[which(cnt1$mark != "other"),]
  labe <- labe[which(labe$IP>labe$Input),]
  p_01 <- ggplot() +
    geom_point(data = cnt1, aes(x = log2(Input +1), y = log2(IP +1), color = mark)) +
    scale_x_continuous(limits = c(0,21)) +
    scale_y_continuous(limits = c(0,21)) +
    scale_color_manual(values = c("blue","red","gray")) +
    geom_text_repel(data = labe, aes(x = log2(Input +1), y = log2(IP +1), label=theName)) +
    labs(x="Log2(CPM+1)(Input)",y="Log2(CPM+1)(IP)",title = "TE: IP vs Input") +
    theme_few() +
    theme(text = element_text(size = 12,family = "serif",color = "black"),
          plot.title = element_text(size = 16,family = "serif",color = "black"),
          legend.position = c(0.8,0.2),
          axis.title = element_text(size = 12,family = "serif",color = "black"),
          axis.text  = element_text(size = 12,family = "serif",color = "black")) +
    geom_abline(slope = 1,intercept = 0,lty="dashed")
  ggsave(plot = p_01, 
         units = "cm", width = 12, height = 12,
         filename = "/ChIP_seq_2/aaSHY/ghyuu/ripseq/20231025/star/pairEND/count/cpm.point.pdf")
}
################################################################################
#5、peak Analysis (pvalue 0.001)
{
  gtf1 <- rtracklayer::import(con = inx4, format = "gtf")
  gtf1 <- as.data.frame(gtf1)
  colnames(gtf1)[c(10,12)] <- c("gene_name","gene_id")
  gtf1 <- GenomicRanges::makeGRangesFromDataFrame(df = gtf1, keep.extra.columns = T)
  txb1 <- suppressWarnings(GenomicFeatures::makeTxDbFromGRanges(gr = gtf1, taxonomyId = 10090))
  twoC <- read.table(inxB)
  gtf2 <- gtf1[which(gtf1$gene_id %in% unique(twoC$V4))]
  txb2 <- suppressWarnings(GenomicFeatures::makeTxDbFromGRanges(gr = gtf2, taxonomyId = 10090))
  #
  pk01 <- suppressWarnings(ChIPseeker::annotatePeak(peak = paste0(path, "peak/ghyuuN_peaks.narrowPeak"), 
                                                    tssRegion = c(-2000, 2000), TxDb = txb1))
  pdf(file = paste0(path, "PLOT/ghyuu-N-IP_chipseekerAnno_all.pdf"), width = 8, height = 8)
  ChIPseeker::plotAnnoPie(pk01, legend.position = "rightside")
  dev.off()
  write.csv(as.data.frame(pk01), paste0(path, "peak/ghyuu-N-IP_chipseekerAnno_all.csv"), quote = F, row.names = F)
  pk02 <- suppressWarnings(ChIPseeker::annotatePeak(peak = paste0(path, "peak/ghyuuN_peaks.narrowPeak"), 
                                                    tssRegion = c(-2000, 2000), TxDb = txb2))
  pdf(file = paste0(path, "PLOT/ghyuu-N-IP_chipseekerAnno_twoC.pdf"), width = 8, height = 8)
  ChIPseeker::plotAnnoPie(pk02, legend.position = "rightside")
  dev.off()
  write.csv(as.data.frame(pk02), paste0(path, "peak/ghyuu-N-IP_chipseekerAnno_twoC.csv"), quote = F, row.names = F)
  #
  pk03 <- suppressWarnings(ChIPseeker::annotatePeak(peak = paste0(path, "peak/ghyuuB_peaks.broadPeak"),
                                                    tssRegion=c(-2000, 2000), TxDb = txb1))
  pdf(file = paste0(path, "PLOT/ghyuu-B-IP_chipseekerAnno_all.pdf"), width = 8, height = 8)
  ChIPseeker::plotAnnoPie(pk03, legend.position = "rightside")
  dev.off()
  write.csv(as.data.frame(pk03), paste0(path, "peak/ghyuu-B-IP_chipseekerAnno_all.csv"), quote = F, row.names = F)
  pk04 <- suppressWarnings(ChIPseeker::annotatePeak(peak = paste0(path, "peak/ghyuuB_peaks.broadPeak"),
                                                    tssRegion=c(-2000, 2000), TxDb = txb2))
  pdf(file = paste0(path, "PLOT/ghyuu-B-IP_chipseekerAnno_twoC.pdf"), width = 8, height = 8)
  ChIPseeker::plotAnnoPie(pk04, legend.position = "rightside")
  dev.off()
  write.csv(as.data.frame(pk04), paste0(path, "peak/ghyuu-B-IP_chipseekerAnno_twoC.csv"), quote = F, row.names = F)
}
################################################################################
#6、TE Binding situation [RIP-seq Specific]
{
  #Observed
  red1 <- read.table(file = paste0(path, "count/ghyuu.cntTable"), sep = "\t", stringsAsFactors = F, header = T, row.names = 1)
  colnames(red1) <- sapply(str_split(colnames(red1),"\\."),"[",9)
  red2 <- red1[grep(rownames(red1), pattern = ":", invert = F),]
  #Expected
  mmu1 <- read.table(file = "/ChIP_seq_2/aaSHY/ghyuu/rnaSeq/count/ccghyuu.cntTable", sep = "\t", header = T, stringsAsFactors = F, row.names = 1)
  mmu1 <- mmu1[grep(rownames(mmu1), pattern = ":", invert = F),]
  mmu2 <- as.data.frame(apply(mmu1, 2, function(x){x/sum(x)*1000000}))[,c(4,5,6)]
  mmu2$cpm <- rowMeans(mmu2)
  mmu2 <- mmu2[which(mmu2$cpm>=1),]
  red2 <- red2[which(rownames(red2) %in% rownames(mmu2)),]
  mmu3 <- mmu2[which(rownames(mmu2) %in% rownames(red2)),]
  sumReads <- sum(rowMeans(mmu1[which(rownames(mmu1) %in% rownames(mmu3)),c(4,5,6)]))
  for (i in rownames(red2)){
    print(i)
    m = i
    red2[i,3] <- round(as.numeric(rowMeans(mmu1[m,c(4,5,6)])*sum(red2$IP))/sumReads,4)
    red2[i,4] <- round(as.numeric(rowMeans(mmu1[m,c(4,5,6)])*sum(red2$Input))/sumReads,4)
    red2[i,5] <- suppressWarnings(fisher.test(matrix(as.numeric(red2[i,1:4]),ncol = 2,byrow = F),alternative = "two.sided")$p.value)
    red2[i,6] <- suppressWarnings(fisher.test(matrix(as.numeric(red2[i,1:4]),ncol = 2,byrow = F),alternative = "greater")$p.value)
    red2[i,7] <- suppressWarnings(fisher.test(matrix(as.numeric(red2[i,1:4]),ncol = 2,byrow = F),alternative = "less")$p.value)
  }
  colnames(red2) <- c("ghyuuIP","ghyuuInput","expcIP","expcInput","pvalue","greate","lesser")
  red2$pMark[red2$pvalue<0.05] <- "arresting"
  red2$gMark[red2$greate<0.05] <- "arresting"
  red2$lMark[red2$lesser<0.05] <- "arresting"
  red2$geneName <- rownames(red2)
  red2 <- red2[order(red2$greate, decreasing = F),]
  openxlsx::write.xlsx(x = red2, asTable = T, file = paste0(path,"fisher/ghyuuFisher_TE_v1.xlsx"))
}
################################################################################
#7、Gene Binding situation [RIP-seq Specific]
{
  #Observed
  red1 <- read.table(file = paste0(path, "count/ghyuu.cntTable"), sep = "\t", stringsAsFactors = F, header = T, row.names = 1)
  colnames(red1) <- sapply(str_split(colnames(red1),"\\."),"[",9)
  red2 <- red1[grep(rownames(red1), pattern = ":", invert = T),]
  #Expected
  mmu1 <- read.table(file = "/ChIP_seq_2/aaSHY/ghyuu/rnaSeq/count/ccghyuu.cntTable", sep = "\t", header = T, stringsAsFactors = F, row.names = 1)
  mmu1 <- mmu1[grep(rownames(mmu1), pattern = ":", invert = T),]
  mmu2 <- as.data.frame(apply(mmu1, 2, function(x){x/sum(x)*1000000}))[,c(4,5,6)]
  mmu2$cpm <- rowMeans(mmu2)
  mmu2 <- mmu2[which(mmu2$cpm>=1),]
  red2 <- red2[which(rownames(red2) %in% rownames(mmu2)),]
  mmu3 <- mmu2[which(rownames(mmu2) %in% rownames(red2)),]
  sumReads <- sum(rowMeans(mmu1[which(rownames(mmu1) %in% rownames(mmu3)),c(4,5,6)]))
  for (i in rownames(red2)){
    print(i)
    m = i
    red2[i,3] <- round(as.numeric(rowMeans(mmu1[m,c(4,5,6)])*sum(red2$IP))/sumReads,4)
    red2[i,4] <- round(as.numeric(rowMeans(mmu1[m,c(4,5,6)])*sum(red2$Input))/sumReads,4)
    red2[i,5] <- suppressWarnings(fisher.test(matrix(as.numeric(red2[i,1:4]),ncol = 2,byrow = F),alternative = "two.sided")$p.value)
    red2[i,6] <- suppressWarnings(fisher.test(matrix(as.numeric(red2[i,1:4]),ncol = 2,byrow = F),alternative = "greater")$p.value)
    red2[i,7] <- suppressWarnings(fisher.test(matrix(as.numeric(red2[i,1:4]),ncol = 2,byrow = F),alternative = "less")$p.value)
  }
  colnames(red2) <- c("ghyuuIP","ghyuuInput","expcIP","expcInput","pvalue","greate","lesser")
  red2$pMark[red2$pvalue<0.05] <- "arresting"
  red2$gMark[red2$greate<0.05] <- "arresting"
  red2$lMark[red2$lesser<0.05] <- "arresting"
  red2$geneName <- rownames(red2)
  red2 <- red2[order(red2$greate, decreasing = F),]
  openxlsx::write.xlsx(x = red2, asTable = T, file = paste0(path,"fisher/ghyuuFisher_Gene_v1.xlsx"))
}









