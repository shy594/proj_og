#write by Sun Haiayng at 2024.08.29 [单端模式]
#01、搭环境
{
  rm(list=ls())
  for (i in c("ReactomePA","ggplot2","stringr","ChIPseeker","DESeq2",
              "clusterProfiler","topGO","Rgraphviz","pathview","org.Mm.eg.db","TxDb.Mmusculus.UCSC.mm10.knownGene")) {library(i, character.only = T)};rm(i)
}
#02、给索引
{
  path = "/ChIP_seq_2/aaSHY/ghyuu/ripSEQ/Aef/single/"
  inx1 = "/Reference/aaSHY/Ensembl99/GRCm38/INDEX/bowtie2/mm10"
  inx2 = "/Reference/aaSHY/Ensembl99/GRCm38/INDEX/STARv279a"
  inx3 = "/Reference/aaSHY/Ensembl99/GRCm38/fa/GRCm38.dna.primary_assembly.fa"
  inx4 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/Mus_musculus.GRCm38.99.gtf"
  inx5 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf"
  inx6 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf.locInd"
  inx7 = "/Reference/aaSHY/Ensembl99/GRCm38/BED/TEsubFamily/MM-int::ERVL::LTR.bed"
  inx8 = "/Reference/aaSHY/Ensembl99/GRCm38/BED/special/MM_zyw.bed"
  inx9 = "/Reference/aaSHY/Ensembl99/GRCm38/BED/TEsubFamily/MM22::ERVL::LTR.bed"
  inxA = "/Reference/aaSHY/Ensembl99/GRCm38/BED/special/allGene.bed"
  inxB = "/Reference/aaSHY/Ensembl99/GRCm38/BED/special/twoC.bed"
  inxC = "/Reference/aaSHY/zOther/Rn45s/INDEX/STARv279a"
  inxD = "/Reference/aaSHY/zOther/TEconsensus/INDEX/MuERVL/STARv279a"
  inxE = "/Reference/aaSHY/zOther/TEconsensus/INDEX/theAlus/STARv279a"
}
################################################################################
#03、跑命令
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
  lays <- paste0(path, "fq/", prex, ".fq.gz")
  cmd_01 <- paste0("trim_galore --phred33 --fastqc --illumina ",
                   "--clip_R1 10 --three_prime_clip_R1 2 ",
                   "-j 2 -o ",path,"trim ",lays,"\n")
  lays <- paste0(path, "trim/", prex, "_trimmed.fq.gz")
  cmd_02 <- paste0("STAR --runThreadN 10 --readFilesCommand zcat ",
                   "--outMultimapperOrder Random --outSAMtype BAM SortedByCoordinate ",
                   "--outFileNamePrefix ",path,"bam/",prex," ",
                   "--outSAMprimaryFlag AllBestScore --twopassMode Basic ",
                   "--genomeDir ",inx2," --readFilesIn ",lays,"\n")
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
                   path,"bam/",prex[2],".bam -n AefN --outdir ",
                   path,"peak/", "\n")
  cmd_10 <- paste0("macs2 callpeak -f BAM -g mm -p 0.001 --broad --broad-cutoff 0.001 ",
                   "--keep-dup 1 -t ",path,"bam/", prex[1],
                   ".bam -c ", path,"bam/", prex[2],
                   ".bam -n AefB --outdir ", path,"peak/", "\n")
  cmd_11 <- paste0("TEtranscripts -t ",
                   path,"bam/",prex[1],".bam"," -c ",
                   path,"bam/",prex[2],".bam"," ",
                   "--GTF ",inx4," --TE ",inx5," --sortByPos --mode multi ",
                   "--project Aef --outdir ",path,"count/","\n")
  cmd_12 <- paste0("TElocal -b ",
                   path,"bam/",prex,".bam ",
                   "--GTF ",inx4," --TE ",inx6," --sortByPos --project ",
                   path, "count/", prex,"\n")
  cmd_13 <- paste0("STAR --runThreadN 10 --readFilesCommand zcat ",
                   "--outMultimapperOrder Random --outSAMtype BAM SortedByCoordinate ",
                   "--outFileNamePrefix ",path,"rRNA/",prex," ",
                   "--outSAMprimaryFlag AllBestScore --twopassMode Basic ",
                   "--genomeDir ",inxC," --readFilesIn ",lays,"\n")
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
                   "--genomeDir ",inxD," --readFilesIn ",lays,"\n")
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
{
  bw01 <- "/ChIP_seq_2/aaSHY/ghyuu/ripSEQ/Aef/single/bw/bcv/IP.bw"
  bw02 <- "/ChIP_seq_2/aaSHY/ghyuu/ripSEQ/Aef/single/bw/bcv/Input.bw"
  bw03 <- "/ChIP_seq_2/aaSHY/Npm1/ripSEQ/star/single/bw/bcv/Input-oth.bw"
  shell <- paste0(path,"src/run_f.sh")
  cat("#!/bin/bash\n", file = shell)
  cmd_01 <- paste0("computeMatrix scale-regions -p 10 -a 2000 -b 2000 -m 5000 ",
                   "--missingDataAsZero ", "-R ", c(inx7,inx8,inx9,inxA,inxB)," -S ",
                   paste(bw01,bw02,bw03, collapse = " ")," -o ", path,
                   "deeptools/", basename(c(inx7,inx8,inx9,inxA,inxB)), ".mat.gz", "\n")
  cmd_02 <- paste0("plotHeatmap --heatmapHeight 12 --colorList white,red -m ",path,
                   "deeptools/",basename(c(inx7,inx8,inx9,inxA,inxB)),
                   ".mat.gz -out ",path,"PLOT/",
                   basename(c(inx7,inx8,inx9,inxA,inxB)),".v2.pdf","\n")
  for (i in cmd_01) {cat(i, append = T, file = shell)}
  for (i in cmd_02) {cat(i, append = T, file = shell)}
  print(paste0("nohup bash ",shell," >",path,"Log/run_f.log"," 2>&1 &"))
}
################################################################################
#补充
{
  #shell <- paste0(path,"src/run_b.sh")
  #cat("#!/bin/bash\n", file = shell)
  #prex <- c("IP","Input")
  #lays <- paste0(path,"trim/",prex,"_trimmed.fq.gz")
  #cmd_21 <- paste0("STAR --runThreadN 10 --readFilesCommand zcat ",
  #                 "--outMultimapperOrder Random --outSAMtype BAM SortedByCoordinate ",
  #                 "--outFileNamePrefix ",path,"theAlus/",prex," ",
  #                 "--outSAMprimaryFlag AllBestScore --twopassMode Basic ",
  #                 "--genomeDir ",inxE," --readFilesIn ",lays,"\n")
  #cmd_22 <- paste0("mv ",path,"theAlus/",prex,"Aligned.sortedByCoord.out.bam ",path,
  #                 "theAlus/",prex,".bam","\n")
  #cmd_23 <- paste0("samtools index ",path,"theAlus/",prex,".bam","\n")
  #cmd_24 <- paste0("bamCoverage --normalizeUsing CPM -p 10 --minMappingQuality 1 ",
  #                 "--samFlagExclude 256 -of bigwig -bs 20 -b ",path,
  #                 "theAlus/",prex,".bam -o ",path,
  #                 "theAlus/",prex,".theAlus.bw\n")
  #for (i in cmd_21) {cat(i, append = T, file = shell)}
  #for (i in cmd_22) {cat(i, append = T, file = shell)}
  #for (i in cmd_23) {cat(i, append = T, file = shell)}
  #for (i in cmd_24) {cat(i, append = T, file = shell)}
  #print(paste0("nohup bash ",shell," >",path,"Log/run_b.log"," 2>&1 &"))
}
################################################################################
#04、计算CPM Excel
{
  cnt1 <- read.table(paste0(path, "count/Aef.cntTable"), header = T, sep = "\t", row.names = 1)
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
  openxlsx::write.xlsx(x = tmp1, file = paste0(path,"count/Aef.cpm.onlyTE.xlsx"))
  openxlsx::write.xlsx(x = tmp2, file = paste0(path,"count/Aef.cpm.onlyGene.xlsx"))
  openxlsx::write.xlsx(x = tmp3, file = paste0(path,"count/Aef.cpm.allGeneTE.xlsx"))
}
################################################################################
#05、CPM点图
{
  cnt1 <- openxlsx::read.xlsx(xlsxFile = paste0(path,"count/Aef.cpm.onlyTE.xlsx"))
  cnt1$sort <- (cnt1$IP+1)/(cnt1$Input+1)
  cnt1 <- cnt1[order(cnt1$sort, decreasing = T),]
  cnt1$mark <- "other"
  cnt1$mark[grep(x = cnt1$theName, pattern = "ERV1")] <- "ERV1"
  cnt1$mark[grep(x = cnt1$theName, pattern = "ERVL")] <- "ERVL"
  cnt1$mark[grep(x = cnt1$theName, pattern = "ERVK")] <- "ERVK"
  labe <- cnt1[which(cnt1$mark != "other"),]
  labe <- labe[which(labe$sort > 5),]
  p_01 <- ggplot() +
    geom_point(data = cnt1, aes(x = log2(Input +1), y = log2(IP +1), color = mark)) +
    scale_x_continuous(limits = c(0,21)) +
    scale_y_continuous(limits = c(0,21)) +
    scale_color_manual(values = c("blue","red","green","gray")) +
    geom_text_repel(data = labe, aes(x = log2(Input +1), y = log2(IP +1), label=theName)) +
    labs(x="Log2(CPM+1)(Input)",y="Log2(CPM+1)(IP)",title = "TE: IP vs Input") +
    theme_few() +
    theme(text = element_text(size = 12,family = "serif",color = "black"),
          plot.title = element_text(size = 16,family = "serif",color = "black"),
          legend.position = c(0.9,0.6),
          axis.title = element_text(size = 12,family = "serif",color = "black"),
          axis.text  = element_text(size = 12,family = "serif",color = "black")) +
    geom_abline(slope = 1,intercept = 0,lty="dashed")
  ggsave(plot = p_01,
         units = "cm", width = 20, height = 15,
         filename = paste0(path,"count/cpm.point.pdf"))
}
################################################################################
#06、peak Analysis (pvalue 0.001)
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
  pk01 <- suppressWarnings(ChIPseeker::annotatePeak(peak = paste0(path, "peak/AefN_peaks.narrowPeak"), 
                                                    tssRegion = c(-2000, 2000), TxDb = txb1))
  pdf(file = paste0(path, "PLOT/Aef-N-IP_chipseekerAnno_all.pdf"), width = 8, height = 8)
  ChIPseeker::plotAnnoPie(pk01, legend.position = "rightside")
  dev.off()
  write.csv(as.data.frame(pk01), paste0(path, "peak/Aef-N-IP_chipseekerAnno_all.csv"), quote = F, row.names = F)
  pk02 <- suppressWarnings(ChIPseeker::annotatePeak(peak = paste0(path, "peak/AefN_peaks.narrowPeak"), 
                                                    tssRegion = c(-2000, 2000), TxDb = txb2))
  pdf(file = paste0(path, "PLOT/Aef-N-IP_chipseekerAnno_twoC.pdf"), width = 8, height = 8)
  ChIPseeker::plotAnnoPie(pk02, legend.position = "rightside")
  dev.off()
  write.csv(as.data.frame(pk02), paste0(path, "peak/Aef-N-IP_chipseekerAnno_twoC.csv"), quote = F, row.names = F)
  #
  pk03 <- suppressWarnings(ChIPseeker::annotatePeak(peak = paste0(path, "peak/AefB_peaks.broadPeak"),
                                                    tssRegion=c(-2000, 2000), TxDb = txb1))
  pdf(file = paste0(path, "PLOT/Aef-B-IP_chipseekerAnno_all.pdf"), width = 8, height = 8)
  ChIPseeker::plotAnnoPie(pk03, legend.position = "rightside")
  dev.off()
  write.csv(as.data.frame(pk03), paste0(path, "peak/Aef-B-IP_chipseekerAnno_all.csv"), quote = F, row.names = F)
  pk04 <- suppressWarnings(ChIPseeker::annotatePeak(peak = paste0(path, "peak/AefB_peaks.broadPeak"),
                                                    tssRegion=c(-2000, 2000), TxDb = txb2))
  pdf(file = paste0(path, "PLOT/Aef-B-IP_chipseekerAnno_twoC.pdf"), width = 8, height = 8)
  ChIPseeker::plotAnnoPie(pk04, legend.position = "rightside")
  dev.off()
  write.csv(as.data.frame(pk04), paste0(path, "peak/Aef-B-IP_chipseekerAnno_twoC.csv"), quote = F, row.names = F)
}
################################################################################
#07、TE Binding situation [RIP-seq Specific]
{
  #Observed
  red1 <- read.table(file = paste0(path, "count/Aef.cntTable"), 
                     sep = "\t", stringsAsFactors = F, header = T, row.names = 1)
  colnames(red1) <- sapply(str_split(colnames(red1),"\\."),"[",9)
  red2 <- red1[grep(rownames(red1), pattern = ":", invert = F),]
  #Expected
  mmu1 <- read.table(file = "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/Aef/count/Aef.cntTable", sep = "\t", header = T, stringsAsFactors = F, row.names = 1)
  mmu1 <- mmu1[grep(rownames(mmu1), pattern = ":", invert = F),]
  mmu2 <- as.data.frame(apply(mmu1, 2, function(x){x/sum(x)*1000000}))[,c(3,4)]
  mmu2$cpm <- rowMeans(mmu2)
  mmu2 <- mmu2[which(mmu2$cpm>=1),]
  red2 <- red2[which(rownames(red2) %in% rownames(mmu2)),]
  mmu3 <- mmu2[which(rownames(mmu2) %in% rownames(red2)),]
  sumReads <- sum(rowMeans(mmu1[which(rownames(mmu1) %in% rownames(mmu3)),c(3,4)]))
  for (i in rownames(red2)){
    print(i)
    m = i
    red2[i,3] <- round(as.numeric(rowMeans(mmu1[m,c(3,4)])*sum(red2$IP))/sumReads,4)
    red2[i,4] <- round(as.numeric(rowMeans(mmu1[m,c(3,4)])*sum(red2$Input))/sumReads,4)
    red2[i,5] <- suppressWarnings(fisher.test(matrix(as.numeric(red2[i,1:4]),ncol = 2,byrow = F),alternative = "two.sided")$p.value)
    red2[i,6] <- suppressWarnings(fisher.test(matrix(as.numeric(red2[i,1:4]),ncol = 2,byrow = F),alternative = "greater")$p.value)
    red2[i,7] <- suppressWarnings(fisher.test(matrix(as.numeric(red2[i,1:4]),ncol = 2,byrow = F),alternative = "less")$p.value)
  }
  colnames(red2) <- c("AefIP","AefInput","expcIP","expcInput","pvalue","greate","lesser")
  red2$pMark[red2$pvalue<0.05] <- "arresting"
  red2$gMark[red2$greate<0.05] <- "arresting"
  red2$lMark[red2$lesser<0.05] <- "arresting"
  red2$geneName <- rownames(red2)
  red2 <- red2[order(red2$greate, decreasing = F),]
  openxlsx::write.xlsx(x = red2, asTable = T, file = paste0(path,"fisher/AefFisher_TE_v1.xlsx"))
}
################################################################################
#08、Gene Binding situation [RIP-seq Specific]
{
  #Observed
  red1 <- read.table(file = paste0(path, "count/Aef.cntTable"), sep = "\t", stringsAsFactors = F, header = T, row.names = 1)
  colnames(red1) <- sapply(str_split(colnames(red1),"\\."),"[",9)
  red2 <- red1[grep(rownames(red1), pattern = ":", invert = T),]
  #Expected
  mmu1 <- read.table(file = "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/Aef/count/Aef.cntTable", sep = "\t", header = T, stringsAsFactors = F, row.names = 1)
  mmu1 <- mmu1[grep(rownames(mmu1), pattern = ":", invert = T),]
  mmu2 <- as.data.frame(apply(mmu1, 2, function(x){x/sum(x)*1000000}))[,c(3,4)]
  mmu2$cpm <- rowMeans(mmu2)
  mmu2 <- mmu2[which(mmu2$cpm>=1),]
  red2 <- red2[which(rownames(red2) %in% rownames(mmu2)),]
  mmu3 <- mmu2[which(rownames(mmu2) %in% rownames(red2)),]
  sumReads <- sum(rowMeans(mmu1[which(rownames(mmu1) %in% rownames(mmu3)),c(3,4)]))
  for (i in rownames(red2)){
    print(i)
    m = i
    red2[i,3] <- round(as.numeric(rowMeans(mmu1[m,c(3,4)])*sum(red2$IP))/sumReads,4)
    red2[i,4] <- round(as.numeric(rowMeans(mmu1[m,c(3,4)])*sum(red2$Input))/sumReads,4)
    red2[i,5] <- suppressWarnings(fisher.test(matrix(as.numeric(red2[i,1:4]),ncol = 2,byrow = F),alternative = "two.sided")$p.value)
    red2[i,6] <- suppressWarnings(fisher.test(matrix(as.numeric(red2[i,1:4]),ncol = 2,byrow = F),alternative = "greater")$p.value)
    red2[i,7] <- suppressWarnings(fisher.test(matrix(as.numeric(red2[i,1:4]),ncol = 2,byrow = F),alternative = "less")$p.value)
  }
  colnames(red2) <- c("AefIP","AefInput","expcIP","expcInput","pvalue","greate","lesser")
  red2$pMark[red2$pvalue<0.05] <- "arresting"
  red2$gMark[red2$greate<0.05] <- "arresting"
  red2$lMark[red2$lesser<0.05] <- "arresting"
  red2$geneName <- rownames(red2)
  red2 <- red2[order(red2$greate, decreasing = F),]
  openxlsx::write.xlsx(x = red2, asTable = T, file = paste0(path,"fisher/AefFisher_Gene_v1.xlsx"))
}
#
#
#
#
#
#
#
################################################################################
#09[替换Input为ff23k Input]
{
  bam1 <- paste0(path,"bam/",prex[1],".bam")
  bam2 <- "/ChIP_seq_2/aaSHY/ghyuu/ripSEQ/ff23k/star/single/bam/Input-oth.bam"
  bw01 <- paste0(path,"bw/bcv/",prex,".bw")[1]
  bw02 <- "/ChIP_seq_2/aaSHY/ghyuu/ripSEQ/ff23k/star/single/bw/bcv/Input-oth.bw"
  shell <- paste0(path,"src/run_g.sh")
  cat("#!/bin/bash\n", file = shell)
  cmd_05 <- paste0("bamCompare -p 20 --scaleFactorsMethod SES --sampleLength 100000 ",
                   "-b1 ",bam1," -b2 ",bam2," -o ",path,
                   "bw/bcp/IP_Input-s.bw","\n")
  cmd_06 <- paste0("computeMatrix scale-regions -p 10 -a 2000 -b 2000 -m 5000 ",
                   "--missingDataAsZero ", "-R ", c(inx7,inx8,inx9,inxA,inxB)," -S ", path,
                   "bw/bcp/IP_Input-s.bw -o ", path,
                   "deeptools/", basename(c(inx7,inx8,inx9,inxA,inxB)), "-4.mat.gz", "\n")
  cmd_07 <- paste0("plotHeatmap --heatmapHeight 12 --colorList white,red -m ",path,
                   "deeptools/",basename(c(inx7,inx8,inx9,inxA,inxB)),
                   "-4.mat.gz -out ",path,"PLOT/",
                   basename(c(inx7,inx8,inx9,inxA,inxB)),"-bcp-s.pdf","\n")
  cmd_08 <- paste0("macs2 callpeak -f BAM -g mm -p 0.001 --keep-dup 1 -t ",
                   bam1," -c ",bam2," -n AefN-s --outdir ",
                   path,"peak/", "\n")
  cmd_09 <- paste0("macs2 callpeak -f BAM -g mm -p 0.001 --broad --broad-cutoff 0.001 ",
                   "--keep-dup 1 -t ",bam1," -c ",bam2," -n AefB-s ",
                   "--outdir ",path,"peak/", "\n")
  for (i in cmd_05) {cat(i, append = T, file = shell)}
  for (i in cmd_06) {cat(i, append = T, file = shell)}
  for (i in cmd_07) {cat(i, append = T, file = shell)}
  for (i in cmd_08) {cat(i, append = T, file = shell)}
  for (i in cmd_09) {cat(i, append = T, file = shell)}
  print(paste0("nohup bash ",shell," >",path,"Log/run_g.log"," 2>&1 &"))
}
################################################################################
#10、计算CPM Excel
{
  cnt1 <- read.table(paste0(path, "count/Aef-s.cntTable"), header = T, sep = "\t", row.names = 1)
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
  openxlsx::write.xlsx(x = tmp1, file = paste0(path,"count/Aef-s.cpm.onlyTE.xlsx"))
  openxlsx::write.xlsx(x = tmp2, file = paste0(path,"count/Aef-s.cpm.onlyGene.xlsx"))
  openxlsx::write.xlsx(x = tmp3, file = paste0(path,"count/Aef-s.cpm.allGeneTE.xlsx"))
}
################################################################################
#11、CPM点图
{
  cnt1 <- openxlsx::read.xlsx(xlsxFile = paste0(path,"count/Aef-s.cpm.onlyTE.xlsx"))
  cnt1$sort <- (cnt1$IP+1)/(cnt1$Input+1)
  cnt1 <- cnt1[order(cnt1$sort, decreasing = T),]
  cnt1$mark <- "other"
  cnt1$mark[grep(x = cnt1$theName, pattern = "ERV1")] <- "ERV1"
  cnt1$mark[grep(x = cnt1$theName, pattern = "ERVL")] <- "ERVL"
  cnt1$mark[grep(x = cnt1$theName, pattern = "ERVK")] <- "ERVK"
  labe <- cnt1[which(cnt1$mark != "other"),]
  labe <- labe[which(labe$sort > 5),]
  p_01 <- ggplot() +
    geom_point(data = cnt1, aes(x = log2(Input +1), y = log2(IP +1), color = mark)) +
    scale_x_continuous(limits = c(0,21)) +
    scale_y_continuous(limits = c(0,21)) +
    scale_color_manual(values = c("blue","red","green","gray")) +
    geom_text_repel(data = labe, aes(x = log2(Input +1), y = log2(IP +1), label=theName)) +
    labs(x="Log2(CPM+1)(Input)",y="Log2(CPM+1)(IP)",title = "TE: IP vs Input") +
    theme_few() +
    theme(text = element_text(size = 12,family = "serif",color = "black"),
          plot.title = element_text(size = 16,family = "serif",color = "black"),
          legend.position = c(0.9,0.6),
          axis.title = element_text(size = 12,family = "serif",color = "black"),
          axis.text  = element_text(size = 12,family = "serif",color = "black")) +
    geom_abline(slope = 1,intercept = 0,lty="dashed")
  ggsave(plot = p_01,
         units = "cm", width = 20, height = 15,
         filename = paste0(path,"count/cpm-s.point.pdf"))
}
################################################################################
#12、TE Binding situation [RIP-seq Specific]
{
  #Observed
  red1 <- read.table(file = paste0(path, "count/Aef-s.cntTable"), 
                     sep = "\t", stringsAsFactors = F, header = T, row.names = 1)
  colnames(red1) <- sapply(str_split(colnames(red1),"\\."),"[",9)
  red2 <- red1[grep(rownames(red1), pattern = ":", invert = F),]
  #Expected
  mmu1 <- read.table(file = "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/Aef/count/Aef.cntTable", sep = "\t", header = T, stringsAsFactors = F, row.names = 1)
  mmu1 <- mmu1[grep(rownames(mmu1), pattern = ":", invert = F),]
  mmu2 <- as.data.frame(apply(mmu1, 2, function(x){x/sum(x)*1000000}))[,c(3,4)]
  mmu2$cpm <- rowMeans(mmu2)
  mmu2 <- mmu2[which(mmu2$cpm>=1),]
  red2 <- red2[which(rownames(red2) %in% rownames(mmu2)),]
  mmu3 <- mmu2[which(rownames(mmu2) %in% rownames(red2)),]
  sumReads <- sum(rowMeans(mmu1[which(rownames(mmu1) %in% rownames(mmu3)),c(3,4)]))
  for (i in rownames(red2)){
    print(i)
    m = i
    red2[i,3] <- round(as.numeric(rowMeans(mmu1[m,c(3,4)])*sum(red2$IP))/sumReads,4)
    red2[i,4] <- round(as.numeric(rowMeans(mmu1[m,c(3,4)])*sum(red2$Input))/sumReads,4)
    red2[i,5] <- suppressWarnings(fisher.test(matrix(as.numeric(red2[i,1:4]),ncol = 2,byrow = F),alternative = "two.sided")$p.value)
    red2[i,6] <- suppressWarnings(fisher.test(matrix(as.numeric(red2[i,1:4]),ncol = 2,byrow = F),alternative = "greater")$p.value)
    red2[i,7] <- suppressWarnings(fisher.test(matrix(as.numeric(red2[i,1:4]),ncol = 2,byrow = F),alternative = "less")$p.value)
  }
  colnames(red2) <- c("AefIP","AefInput","expcIP","expcInput","pvalue","greate","lesser")
  red2$pMark[red2$pvalue<0.05] <- "arresting"
  red2$gMark[red2$greate<0.05] <- "arresting"
  red2$lMark[red2$lesser<0.05] <- "arresting"
  red2$geneName <- rownames(red2)
  red2 <- red2[order(red2$greate, decreasing = F),]
  openxlsx::write.xlsx(x = red2, asTable = T, file = paste0(path,"fisher/Aef-sFisher_TE_v1.xlsx"))
}
################################################################################
#13、Gene Binding situation [RIP-seq Specific]
{
  #Observed
  red1 <- read.table(file = paste0(path, "count/Aef-s.cntTable"), sep = "\t", stringsAsFactors = F, header = T, row.names = 1)
  colnames(red1) <- sapply(str_split(colnames(red1),"\\."),"[",9)
  red2 <- red1[grep(rownames(red1), pattern = ":", invert = T),]
  #Expected
  mmu1 <- read.table(file = "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/Aef/count/Aef.cntTable", sep = "\t", header = T, stringsAsFactors = F, row.names = 1)
  mmu1 <- mmu1[grep(rownames(mmu1), pattern = ":", invert = T),]
  mmu2 <- as.data.frame(apply(mmu1, 2, function(x){x/sum(x)*1000000}))[,c(3,4)]
  mmu2$cpm <- rowMeans(mmu2)
  mmu2 <- mmu2[which(mmu2$cpm>=1),]
  red2 <- red2[which(rownames(red2) %in% rownames(mmu2)),]
  mmu3 <- mmu2[which(rownames(mmu2) %in% rownames(red2)),]
  sumReads <- sum(rowMeans(mmu1[which(rownames(mmu1) %in% rownames(mmu3)),c(3,4)]))
  for (i in rownames(red2)){
    print(i)
    m = i
    red2[i,3] <- round(as.numeric(rowMeans(mmu1[m,c(3,4)])*sum(red2$IP))/sumReads,4)
    red2[i,4] <- round(as.numeric(rowMeans(mmu1[m,c(3,4)])*sum(red2$Input))/sumReads,4)
    red2[i,5] <- suppressWarnings(fisher.test(matrix(as.numeric(red2[i,1:4]),ncol = 2,byrow = F),alternative = "two.sided")$p.value)
    red2[i,6] <- suppressWarnings(fisher.test(matrix(as.numeric(red2[i,1:4]),ncol = 2,byrow = F),alternative = "greater")$p.value)
    red2[i,7] <- suppressWarnings(fisher.test(matrix(as.numeric(red2[i,1:4]),ncol = 2,byrow = F),alternative = "less")$p.value)
  }
  colnames(red2) <- c("AefIP","AefInput","expcIP","expcInput","pvalue","greate","lesser")
  red2$pMark[red2$pvalue<0.05] <- "arresting"
  red2$gMark[red2$greate<0.05] <- "arresting"
  red2$lMark[red2$lesser<0.05] <- "arresting"
  red2$geneName <- rownames(red2)
  red2 <- red2[order(red2$greate, decreasing = F),]
  openxlsx::write.xlsx(x = red2, asTable = T, file = paste0(path,"fisher/Aef-sFisher_Gene_v1.xlsx"))
}
################################################################################
#14、peak Analysis (pvalue 0.001)
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
  pk01 <- suppressWarnings(ChIPseeker::annotatePeak(peak = paste0(path, "peak/AefN-s_peaks.narrowPeak"), 
                                                    tssRegion = c(-2000, 2000), TxDb = txb1))
  pdf(file = paste0(path, "PLOT/AefN-s-IP_chipseekerAnno_all.pdf"), width = 8, height = 8)
  ChIPseeker::plotAnnoPie(pk01, legend.position = "rightside")
  dev.off()
  write.csv(as.data.frame(pk01), paste0(path, "peak/AefN-s-IP_chipseekerAnno_all.csv"), quote = F, row.names = F)
  pk02 <- suppressWarnings(ChIPseeker::annotatePeak(peak = paste0(path, "peak/AefN-s_peaks.narrowPeak"), 
                                                    tssRegion = c(-2000, 2000), TxDb = txb2))
  pdf(file = paste0(path, "PLOT/AefN-s-IP_chipseekerAnno_twoC.pdf"), width = 8, height = 8)
  ChIPseeker::plotAnnoPie(pk02, legend.position = "rightside")
  dev.off()
  write.csv(as.data.frame(pk02), paste0(path, "peak/AefN-s-IP_chipseekerAnno_twoC.csv"), quote = F, row.names = F)
  #
  pk03 <- suppressWarnings(ChIPseeker::annotatePeak(peak = paste0(path, "peak/AefB-s_peaks.broadPeak"),
                                                    tssRegion=c(-2000, 2000), TxDb = txb1))
  pdf(file = paste0(path, "PLOT/AefB-s-IP_chipseekerAnno_all.pdf"), width = 8, height = 8)
  ChIPseeker::plotAnnoPie(pk03, legend.position = "rightside")
  dev.off()
  write.csv(as.data.frame(pk03), paste0(path, "peak/AefB-s-IP_chipseekerAnno_all.csv"), quote = F, row.names = F)
  pk04 <- suppressWarnings(ChIPseeker::annotatePeak(peak = paste0(path, "peak/AefB-s_peaks.broadPeak"),
                                                    tssRegion=c(-2000, 2000), TxDb = txb2))
  pdf(file = paste0(path, "PLOT/AefB-s-IP_chipseekerAnno_twoC.pdf"), width = 8, height = 8)
  ChIPseeker::plotAnnoPie(pk04, legend.position = "rightside")
  dev.off()
  write.csv(as.data.frame(pk04), paste0(path, "peak/AefB-s-IP_chipseekerAnno_twoC.csv"), quote = F, row.names = F)
}
#
#
#
#
#
#
#
################################################################################
#15、ff23k和Aef在MM融合基因MM 起始的TSS上的富集
{
  tf01 <- import(inx5)
  aims <- tf01[grep("MM-int|MM22", tf01$gene_id)]
  bath <- "/ChIP_seq_1/aaSHY/rna2014Deng/myRun/"
  prex <- c("twoEarly","twoMid","twoLate")
  regs <- data.frame()
  for (i in prex) {
    print(i)
    asem <- import(paste0(bath,"stringtie/",i,".merge.gtf"))
    tran <- asem[which(asem$type=="transcript")]
    chm1 <- asem[which(asem$type=="exon")]
    chm1$exon_number <- as.numeric(chm1$exon_number)
    chm1 <- as.data.frame(chm1)
    chm1 <- chm1 %>%
      dplyr::group_by(transcript_id) %>%
      dplyr::mutate(nums = n())
    chm1 <- chm1[which(chm1$nums != 1),]
    chm1_p <- chm1 %>% 
      dplyr::filter(strand == "+")
    chm1_p <- chm1_p %>%
      dplyr::group_by(transcript_id) %>%
      slice_min(exon_number, n = 1)
    chm1_m <- chm1 %>% 
      dplyr::filter(strand == "-")
    chm1_m <- chm1_m %>%
      dplyr::group_by(transcript_id) %>%
      slice_max(exon_number, n = 1)
    chm2 <- rbind(chm1_p, chm1_m)
    chm2 <- makeGRangesFromDataFrame(df = chm2, keep.extra.columns = T)
    #
    ov_1 <- suppressWarnings(findOverlaps(chm2,aims,ignore.strand=T))
    reg1 <- unique(chm2$transcript_id[queryHits(ov_1)])
    reg1 <- tran[which(tran$transcript_id %in% reg1)]
    reg1 <- GenomicRanges::promoters(x = reg1, upstream = 500, downstream = 500)
    reg1 <- as.data.frame(reg1)[,c(1,2,3,11,4,5)]
    regs <- rbind(regs,reg1)
  }
  regs <- as.data.frame(distinct(regs))
  write.table(x = regs, file = "/ChIP_seq_2/aaSHY/ghyuu/ask/fuss_2/MM.tss.bed",
              quote = F, sep = "\t", col.names = F, row.names = F)
  #
  #
  #
  #
  #
  #
  cath <- ""
  bwss <- c("/ChIP_seq_2/aaSHY/ghyuu/ripSEQ/Aef/pairEND/bw/bcp/IP_Input-5.bw",
            "/ChIP_seq_2/aaSHY/ghyuu/ripSEQ/ff23k/star/pairEND/bw/bcp/IP_Input.bw")
  cmd_01 <- paste0("computeMatrix reference-point --referencePoint center -a 500 -b 500 ",
                   "--missingDataAsZero -p 20 -R ", 
                   "/ChIP_seq_2/aaSHY/ghyuu/ask/fuss_2/MM.tss.bed ",
                   "-S ",paste0(bwss, collapse = " ")," -o ",
                   "/ChIP_seq_2/aaSHY/ghyuu/ask/fuss_2/matrixxx.mat.gz","\n")
  cmd_02 <- paste0("plotHeatmap --heatmapHeight 12 --colorList white,red -m ",
                   "/ChIP_seq_2/aaSHY/ghyuu/ask/fuss_2/matrixxx.mat.gz ",
                   "--samplesLabel ",paste0(c("Aef","ff23k"), collapse = " ")," ",
                   "-out /ChIP_seq_2/aaSHY/ghyuu/ask/fuss_2/xx.pdf","\n")
}



















