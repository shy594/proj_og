#write by Sun Haiayng at 2024.08.29 [双端模式-调参]
#01、搭环境
{
  rm(list=ls())
  for (i in c("ReactomePA","ggplot2","stringr","ChIPseeker","DESeq2","clusterProfiler","topGO","Rgraphviz","pathview","org.Mm.eg.db","TxDb.Mmusculus.UCSC.mm10.knownGene")) {library(i, character.only = T)};rm(i)
}
#02、给索引
{
  path = "/ChIP_seq_2/aaSHY/ghyuu/ripSEQ/Aef/pairEND/"
  inx1 = "/Reference/aaSHY/Ensembl99/GRCm38/INDEX/bowtie2/mm10"
  inx2 = "/Reference/aaSHY/Ensembl99/GRCm38/INDEX/STARv279a"
  inx3 = "/Reference/aaSHY/Ensembl99/GRCm38/fa/GRCm38.dna.primary_assembly.fa"
  inx4 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/Mus_musculus.GRCm38.99.gtf"
  inx5 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf"
  inx6 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf.locInd"
  inx7 = "/Reference/aaSHY/Ensembl99/GRCm38/BED/TEsubFamily/MM::ERVL::LTR.bed"
  inx8 = "/Reference/aaSHY/Ensembl99/GRCm38/BED/special/mervl_zyw.bed"
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
  shell <- paste0(path,"src/run_a.sh")
  cat("#!/bin/bash\n", file = shell)
  prex <- c("IP","Input")
  lay1 <- paste0(path, "trim/", prex, "_1_val_1.fq.gz")
  lay2 <- paste0(path, "trim/", prex, "_2_val_2.fq.gz")
  cmd_01 <- paste0("#STAR --runThreadN 10 --readFilesCommand zcat ",
                   "--peOverlapNbasesMin 5 ",
                   "--outMultimapperOrder Random --outSAMtype BAM SortedByCoordinate ",
                   "--outFileNamePrefix ",path,"bam/",prex,"-5 ",
                   "--outSAMprimaryFlag AllBestScore --twopassMode Basic ",
                   "--genomeDir ",inx2," --readFilesIn ",lay1," ",lay2,"\n")
  cmd_02 <- paste0("#mv ",path,"bam/",prex,"-5Aligned.sortedByCoord.out.bam ",path,
                   "bam/",prex,".bam","\n")
  cmd_03 <- paste0("#samtools index ",path,"bam/",prex,"-5.bam","\n")
  cmd_04 <- paste0("#bamCoverage --normalizeUsing CPM -p 10 --minMappingQuality 1 ",
                   "--samFlagExclude 256 -of bigwig -bs 20 -b ",
                   path,"bam/",prex,"-5.bam -o ",
                   path,"bw/bcv/",prex,"-5.bw","\n")
  cmd_05 <- paste0("#bamCompare -p 10 --scaleFactorsMethod SES --sampleLength 100000 -b1 ",
                   path,"bam/",prex[1],"-5.bam -b2 ",
                   path,"bam/",prex[2],"-5.bam -o ",
                   path,"bw/bcp/IP_Input-5.bw","\n")
  cmd_06 <- paste0("computeMatrix scale-regions -p 10 -a 2000 -b 2000 -m 5000 ",
                   "--missingDataAsZero ", "-R ", c(inx7,inx8,inx9,inxA,inxB)," -S ", path,
                   "bw/bcp/IP_Input-5.bw -o ", path,
                   "deeptools/", basename(c(inx7,inx8,inx9,inxA,inxB)), "-1.mat.gz", "\n")
  cmd_07 <- paste0("plotHeatmap --heatmapHeight 12 --colorList white,red -m ",path,
                   "deeptools/",basename(c(inx7,inx8,inx9,inxA,inxB)),
                   "-1.mat.gz -out ",path,"PLOT/",
                   basename(c(inx7,inx8,inx9,inxA,inxB)),"-bcp-5.pdf","\n")
  cmd_08 <- paste0("#macs2 callpeak -f BAM -g mm -p 0.001 --keep-dup 1 -t ",
                   path,"bam/",prex[1],"-5.bam -c ",
                   path,"bam/",prex[2],"-5.bam -n AefN-5 --outdir ",
                   path,"peak/", "\n")
  cmd_09 <- paste0("#macs2 callpeak -f BAM -g mm -p 0.001 --broad --broad-cutoff 0.001 ",
                   "--keep-dup 1 -t ",path,"bam/", prex[1],
                   "-5.bam -c ", path,"bam/", prex[2],
                   "-5.bam -n AefB-5 --outdir ", path,"peak/", "\n")
  cmd_10 <- paste0("#TEtranscripts -t ",
                   path,"bam/",prex[1],"-5.bam"," -c ",
                   path,"bam/",prex[2],"-5.bam"," ",
                   "--GTF ",inx4," --TE ",inx5," --sortByPos --mode multi ",
                   "--project Aef-5 --outdir ",path,"count/","\n")
  cmd_11 <- paste0("#TElocal -b ",
                   path,"bam/",prex,"-5.bam ",
                   "--GTF ",inx4," --TE ",inx6," --sortByPos --project ",
                   path, "count/", prex,"-5","\n")
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
  print(paste0("nohup bash ",shell," >",path,"Log/run_a.log"," 2>&1 &"))
}
################################################################################
{
  bw01 <- paste0(path,"bw/bcv/",prex,"-5.bw")[1]
  bw02 <- paste0(path,"bw/bcv/",prex,"-5.bw")[2]
  bw03 <- "/ChIP_seq_2/aaSHY/ghyuu/ripSEQ/ff23k/star/pairEND/bw/bcv/Input-5s.bw"
  shell <- paste0(path,"src/run_f.sh")
  cat("#!/bin/bash\n", file = shell)
  cmd_0a <- paste0("computeMatrix scale-regions -p 10 -a 2000 -b 2000 -m 5000 ",
                   "--missingDataAsZero ", "-R ", 
                   c(inx7,inx8,inx9,inxA,inxB)," -S ",
                   paste(bw01,bw02,bw03, collapse = " ")," -o ", path,
                   "deeptools/", basename(c(inx7,inx8,inx9,inxA,inxB)), "-2.mat.gz", "\n")
  cmd_0b <- paste0("plotHeatmap --heatmapHeight 12 --colorList white,red -m ",path,
                   "deeptools/",basename(c(inx7,inx8,inx9,inxA,inxB)),
                   "-2.mat.gz -out ",path,"PLOT/",
                   basename(c(inx7,inx8,inx9,inxA,inxB)),"-bcv-5.pdf","\n")
  for (i in cmd_0a) {cat(i, append = T, file = shell)}
  for (i in cmd_0b) {cat(i, append = T, file = shell)}
  print(paste0("nohup bash ",shell," >",path,"Log/run_f.log"," 2>&1 &"))
}
################################################################################
#04、计算CPM Excel
{
  cnt1 <- read.table(paste0(path, "count/Aef-5.cntTable"), header = T, sep = "\t", row.names = 1)
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
  openxlsx::write.xlsx(x = tmp1, file = paste0(path,"count/Aef-5.cpm.onlyTE.xlsx"))
  openxlsx::write.xlsx(x = tmp2, file = paste0(path,"count/Aef-5.cpm.onlyGene.xlsx"))
  openxlsx::write.xlsx(x = tmp3, file = paste0(path,"count/Aef-5.cpm.allGeneTE.xlsx"))
}
################################################################################
#05、CPM点图
{
  cnt1 <- openxlsx::read.xlsx(xlsxFile = paste0(path,"count/Aef-5.cpm.onlyTE.xlsx"))
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
         filename = paste0(path,"count/cpm.point-5.pdf"))
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
  pk01 <- suppressWarnings(ChIPseeker::annotatePeak(peak = paste0(path, "peak/AefN-5_peaks.narrowPeak"), 
                                                    tssRegion = c(-2000, 2000), TxDb = txb1))
  pdf(file = paste0(path, "PLOT/AefN-5-IP_chipseekerAnno_all.pdf"), width = 8, height = 8)
  ChIPseeker::plotAnnoPie(pk01, legend.position = "rightside")
  dev.off()
  write.csv(as.data.frame(pk01), paste0(path, "peak/AefN-5-IP_chipseekerAnno_all.csv"), quote = F, row.names = F)
  pk02 <- suppressWarnings(ChIPseeker::annotatePeak(peak = paste0(path, "peak/AefN-5_peaks.narrowPeak"), 
                                                    tssRegion = c(-2000, 2000), TxDb = txb2))
  pdf(file = paste0(path, "PLOT/AefN-5-IP_chipseekerAnno_twoC.pdf"), width = 8, height = 8)
  ChIPseeker::plotAnnoPie(pk02, legend.position = "rightside")
  dev.off()
  write.csv(as.data.frame(pk02), paste0(path, "peak/AefN-5-IP_chipseekerAnno_twoC.csv"), quote = F, row.names = F)
  #
  pk03 <- suppressWarnings(ChIPseeker::annotatePeak(peak = paste0(path, "peak/AefB-5_peaks.broadPeak"),
                                                    tssRegion=c(-2000, 2000), TxDb = txb1))
  pdf(file = paste0(path, "PLOT/AefB-5-IP_chipseekerAnno_all.pdf"), width = 8, height = 8)
  ChIPseeker::plotAnnoPie(pk03, legend.position = "rightside")
  dev.off()
  write.csv(as.data.frame(pk03), paste0(path, "peak/AefB-5-IP_chipseekerAnno_all.csv"), quote = F, row.names = F)
  pk04 <- suppressWarnings(ChIPseeker::annotatePeak(peak = paste0(path, "peak/AefB-5_peaks.broadPeak"),
                                                    tssRegion=c(-2000, 2000), TxDb = txb2))
  pdf(file = paste0(path, "PLOT/AefB-5-IP_chipseekerAnno_twoC.pdf"), width = 8, height = 8)
  ChIPseeker::plotAnnoPie(pk04, legend.position = "rightside")
  dev.off()
  write.csv(as.data.frame(pk04), paste0(path, "peak/AefB-5-IP_chipseekerAnno_twoC.csv"), quote = F, row.names = F)
}
################################################################################
#07、TE Binding situation [RIP-seq Specific]
{
  #Observed
  red1 <- read.table(file = paste0(path, "count/Aef-5.cntTable"), 
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
  openxlsx::write.xlsx(x = red2, asTable = T, file = paste0(path,"fisher/Aef-5Fisher_TE_v1.xlsx"))
}
################################################################################
#08、Gene Binding situation [RIP-seq Specific]
{
  #Observed
  red1 <- read.table(file = paste0(path, "count/Aef-5.cntTable"), sep = "\t", stringsAsFactors = F, header = T, row.names = 1)
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
  openxlsx::write.xlsx(x = red2, asTable = T, file = paste0(path,"fisher/Aef-5Fisher_Gene_v1.xlsx"))
}
##
##
##
################################################################################
#01、[替换Input为ff23k Input]
{
  bam1 <- paste0(path,"bam/",prex[1],"-5.bam")
  bam2 <- "/ChIP_seq_2/aaSHY/ghyuu/ripSEQ/ff23k/star/pairEND/bam/Input-5.bam"
  bw01 <- paste0(path,"bw/bcv/",prex,"-5.bw")[1]
  bw02 <- "/ChIP_seq_2/aaSHY/ghyuu/ripSEQ/ff23k/star/pairEND/bw/bcv/Input-5s.bw"
  shell <- paste0(path,"src/run_g.sh")
  cat("#!/bin/bash\n", file = shell)
  cmd_05 <- paste0("bamCompare -p 10 --scaleFactorsMethod SES --sampleLength 100000 ",
                   "-b1 ",bam1," -b2 ",bam2," -o ",path,
                   "bw/bcp/IP_Input-5s.bw","\n")
  cmd_06 <- paste0("computeMatrix scale-regions -p 10 -a 2000 -b 2000 -m 5000 ",
                   "--missingDataAsZero ", "-R ", c(inx7,inx8,inx9,inxA,inxB)," -S ", path,
                   "bw/bcp/IP_Input-5s.bw -o ", path,
                   "deeptools/", basename(c(inx7,inx8,inx9,inxA,inxB)), "-3.mat.gz", "\n")
  cmd_07 <- paste0("plotHeatmap --heatmapHeight 12 --colorList white,red -m ",path,
                   "deeptools/",basename(c(inx7,inx8,inx9,inxA,inxB)),
                   "-3.mat.gz -out ",path,"PLOT/",
                   basename(c(inx7,inx8,inx9,inxA,inxB)),"-bcp-5s.pdf","\n")
  cmd_08 <- paste0("macs2 callpeak -f BAM -g mm -p 0.001 --keep-dup 1 -t ",
                   bam1," -c ",bam2," -n AefN-5s --outdir ",
                   path,"peak/", "\n")
  cmd_09 <- paste0("macs2 callpeak -f BAM -g mm -p 0.001 --broad --broad-cutoff 0.001 ",
                   "--keep-dup 1 -t ",bam1," -c ",bam2," -n AefB-5s ",
                   "--outdir ",path,"peak/", "\n")
  for (i in cmd_05) {cat(i, append = T, file = shell)}
  for (i in cmd_06) {cat(i, append = T, file = shell)}
  for (i in cmd_07) {cat(i, append = T, file = shell)}
  for (i in cmd_08) {cat(i, append = T, file = shell)}
  for (i in cmd_09) {cat(i, append = T, file = shell)}
  print(paste0("nohup bash ",shell," >",path,"Log/run_g.log"," 2>&1 &"))
}
################################################################################
#02、计算CPM Excel
{
  cnt1 <- read.table(paste0(path, "count/Aef-5s.cntTable"), header = T, sep = "\t", row.names = 1)
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
  openxlsx::write.xlsx(x = tmp1, file = paste0(path,"count/Aef-5s.cpm.onlyTE.xlsx"))
  openxlsx::write.xlsx(x = tmp2, file = paste0(path,"count/Aef-5s.cpm.onlyGene.xlsx"))
  openxlsx::write.xlsx(x = tmp3, file = paste0(path,"count/Aef-5s.cpm.allGeneTE.xlsx"))
}
################################################################################
#03、CPM点图
{
  cnt1 <- openxlsx::read.xlsx(xlsxFile = paste0(path,"count/Aef-5s.cpm.onlyTE.xlsx"))
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
         filename = paste0(path,"count/cpm.point-5s.pdf"))
}
################################################################################
#04、peak Analysis (pvalue 0.001)
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
  pk01 <- suppressWarnings(ChIPseeker::annotatePeak(peak = paste0(path, "peak/AefN-5s_peaks.narrowPeak"), 
                                                    tssRegion = c(-2000, 2000), TxDb = txb1))
  pdf(file = paste0(path, "PLOT/AefN-5s-IP_chipseekerAnno_all.pdf"), width = 8, height = 8)
  ChIPseeker::plotAnnoPie(pk01, legend.position = "rightside")
  dev.off()
  write.csv(as.data.frame(pk01), paste0(path, "peak/AefN-5s-IP_chipseekerAnno_all.csv"), quote = F, row.names = F)
  pk02 <- suppressWarnings(ChIPseeker::annotatePeak(peak = paste0(path, "peak/AefN-5s_peaks.narrowPeak"), 
                                                    tssRegion = c(-2000, 2000), TxDb = txb2))
  pdf(file = paste0(path, "PLOT/AefN-5s-IP_chipseekerAnno_twoC.pdf"), width = 8, height = 8)
  ChIPseeker::plotAnnoPie(pk02, legend.position = "rightside")
  dev.off()
  write.csv(as.data.frame(pk02), paste0(path, "peak/AefN-5s-IP_chipseekerAnno_twoC.csv"), quote = F, row.names = F)
  #
  pk03 <- suppressWarnings(ChIPseeker::annotatePeak(peak = paste0(path, "peak/AefB-5s_peaks.broadPeak"),
                                                    tssRegion=c(-2000, 2000), TxDb = txb1))
  pdf(file = paste0(path, "PLOT/AefB-5s-IP_chipseekerAnno_all.pdf"), width = 8, height = 8)
  ChIPseeker::plotAnnoPie(pk03, legend.position = "rightside")
  dev.off()
  write.csv(as.data.frame(pk03), paste0(path, "peak/AefB-5s-IP_chipseekerAnno_all.csv"), quote = F, row.names = F)
  pk04 <- suppressWarnings(ChIPseeker::annotatePeak(peak = paste0(path, "peak/AefB-5s_peaks.broadPeak"),
                                                    tssRegion=c(-2000, 2000), TxDb = txb2))
  pdf(file = paste0(path, "PLOT/AefB-5s-IP_chipseekerAnno_twoC.pdf"), width = 8, height = 8)
  ChIPseeker::plotAnnoPie(pk04, legend.position = "rightside")
  dev.off()
  write.csv(as.data.frame(pk04), paste0(path, "peak/AefB-5s-IP_chipseekerAnno_twoC.csv"), quote = F, row.names = F)
}
################################################################################
#05、TE Binding situation [RIP-seq Specific]
{
  #Observed
  red1 <- read.table(file = paste0(path, "count/Aef-5s.cntTable"), 
                     sep = "\t", stringsAsFactors = F, header = T, row.names = 1)
  colnames(red1) <- sapply(str_split(colnames(red1),"\\."),"[",9)
  red2 <- red1[grep(rownames(red1), pattern = ":", invert = F),]
  colnames(red2) <- prex
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
  openxlsx::write.xlsx(x = red2, asTable = T, file = paste0(path,"fisher/Aef-5sFisher_TE_v1.xlsx"))
}
################################################################################
#06、Gene Binding situation [RIP-seq Specific]
{
  #Observed
  red1 <- read.table(file = paste0(path, "count/Aef-5s.cntTable"), sep = "\t", stringsAsFactors = F, header = T, row.names = 1)
  colnames(red1) <- sapply(str_split(colnames(red1),"\\."),"[",9)
  red2 <- red1[grep(rownames(red1), pattern = ":", invert = T),]
  colnames(red2) <- prex
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
  openxlsx::write.xlsx(x = red2, asTable = T, file = paste0(path,"fisher/Aef-5sFisher_Gene_v1.xlsx"))
}
################################################################################
#
#
#
#
#
#
#
#
#
#
#
#01、画个venn图 [fisher pvalue < 0.05 + 2C gene + Aef KD上调基因]
{
  ff01 <- openxlsx::read.xlsx(xlsxFile = paste0(path,"fisher/Aef-5sFisher_Gene_v1.xlsx"))
  ff01 <- ff01[which(ff01$greate < 0.05),]
  ff02 <- read.table(file = "/Reference/aaSHY/BED/special/TWOC_gene.bed", header = F)
  ff03 <- read.csv(file = "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/Aef/DESeq2/Aef.geneOnly.csv")
  ff03 <- ff03[which(ff03$log2FoldChange > log2(1.5) & ff03$padj <0.05),]
  tmp1 <- list()
  tmp1[["a"]] <- unique(ff01$geneName)
  tmp1[["b"]] <- unique(ff02$V4)
  tmp1[["c"]] <- unique(ff03$X)
  p_01 <- VennDiagram::venn.diagram(x = tmp1, 
                                    filename = NULL, 
                                    main = "main~main~main",
                                    sub  = "sub~sub~sub~sub~sub~sub", 
                                    main.fontfamily = "serif", sub.fontfamily = "serif", 
                                    cex = 1.5, main.cex = 1, sub.cex = 1,
                                    print.mode = "raw",#c("percent", "raw"),
                                    category.names = c("fisherGreater","twoC","AefUp"),
                                    fill=c("#FFFFCC","#CCFFFF","#FFCCCC"),
                                    units = "cm", height = 12, width = 12)
  pdf(file = paste0(path, "PLOT/geneVenn.fisherG.twoC.upup.pdf"))
  grid.draw(p_01)
  dev.off()
  my01 <- intersect(intersect(tmp1[["a"]],tmp1[["b"]]),tmp1[["c"]])
  #要“交集基因”的表达
  my02 <- intersect(tmp1[["a"]],tmp1[["c"]])
  expr <- openxlsx::read.xlsx(xlsxFile = "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/Aef/count/Aef.cpm.onlyGene.xlsx")
  expr <- expr[which(expr$theName %in% my02),]
  openxlsx::write.xlsx(x = expr, file = paste0(path, "PLOT/geneVenn.fisherG.twoC.upup-expr.xlsx"))
}
################################################################################
#02、画个venn图 [fisher pvalue < 0.05 + 2C gene + Aef KD下调基因]
{
  ff01 <- openxlsx::read.xlsx(xlsxFile = paste0(path,"fisher/Aef-5sFisher_Gene_v1.xlsx"))
  ff01 <- ff01[which(ff01$greate < 0.05),]
  ff02 <- read.table(file = "/Reference/aaSHY/BED/special/TWOC_gene.bed", header = F)
  ff03 <- read.csv(file = "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/Aef/DESeq2/Aef.geneOnly.csv")
  ff03 <- ff03[which(ff03$log2FoldChange < -log2(1.5) & ff03$padj <0.05),]
  tmp1 <- list()
  tmp1[["a"]] <- unique(ff01$geneName)
  tmp1[["b"]] <- unique(ff02$V4)
  tmp1[["c"]] <- unique(ff03$X)
  p_01 <- VennDiagram::venn.diagram(x = tmp1, 
                                    filename = NULL, 
                                    main = "main~main~main",
                                    sub  = "sub~sub~sub~sub~sub~sub", 
                                    main.fontfamily = "serif", sub.fontfamily = "serif", 
                                    cex = 1.5, main.cex = 1, sub.cex = 1,
                                    print.mode = "raw",#c("percent", "raw"),
                                    category.names = c("fisherGreater","twoC","AefDown"),
                                    fill=c("#FFFFCC","#CCFFFF","#FFCCCC"),
                                    units = "cm", height = 12, width = 12)
  pdf(file = paste0(path, "PLOT/geneVenn.fisherG.twoC.down.pdf"))
  grid.draw(p_01)
  dev.off()
  my01 <- intersect(intersect(tmp1[["a"]],tmp1[["b"]]),tmp1[["c"]])
  #要“交集基因”的表达
  my02 <- intersect(tmp1[["a"]],tmp1[["c"]])
  expr <- openxlsx::read.xlsx(xlsxFile = "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/Aef/count/Aef.cpm.onlyGene.xlsx")
  expr <- expr[which(expr$theName %in% my02),]
  openxlsx::write.xlsx(x = expr, file = paste0(path, "PLOT/geneVenn.fisherG.twoC.down-expr.xlsx"))
}
################################################################################
#03、画个venn图 [fisher pvalue < 0.05 + Aef KD上调TE]
{
  ff01 <- openxlsx::read.xlsx(xlsxFile = paste0(path,"fisher/Aef-5sFisher_TE_v1.xlsx"))
  ff01 <- ff01[which(ff01$greate < 0.05),]
  ff03 <- read.csv(file = "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/Aef/DESeq2/Aef.teOnly.csv")
  ff03 <- ff03[which(ff03$log2FoldChange > log2(1.5) & ff03$padj <0.05),]
  tmp1 <- list()
  tmp1[["a"]] <- unique(ff01$geneName)
  tmp1[["b"]] <- unique(ff03$X)
  p_01 <- VennDiagram::venn.diagram(x = tmp1, 
                                    filename = NULL, 
                                    main = "main~main~main",
                                    sub  = "sub~sub~sub~sub~sub~sub", 
                                    main.fontfamily = "serif", sub.fontfamily = "serif", 
                                    cex = 1.5, main.cex = 1, sub.cex = 1,
                                    print.mode = "raw",#c("percent", "raw"),
                                    category.names = c("fisherGreater","AefUp"),
                                    fill=c("#FFFFCC","#CCFFFF"),
                                    units = "cm", height = 12, width = 12)
  pdf(file = paste0(path, "PLOT/teVenn.fisherG.upup.pdf"))
  grid.draw(p_01)
  dev.off()
  #要“交集基因”的表达
  my02 <- intersect(tmp1[["a"]],tmp1[["b"]])
  expr <- openxlsx::read.xlsx(xlsxFile = "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/Aef/count/Aef.cpm.onlyTE.xlsx")
  expr <- expr[which(expr$theName %in% my02),]
  openxlsx::write.xlsx(x = expr, file = paste0(path, "PLOT/teVenn.fisherG.upup-expr.xlsx"))
}
################################################################################
#04、画个venn图 [fisher pvalue < 0.05 + Aef KD下调TE]
{
  ff01 <- openxlsx::read.xlsx(xlsxFile = paste0(path,"fisher/Aef-5sFisher_TE_v1.xlsx"))
  ff01 <- ff01[which(ff01$greate < 0.05),]
  ff03 <- read.csv(file = "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/Aef/DESeq2/Aef.teOnly.csv")
  ff03 <- ff03[which(ff03$log2FoldChange < -log2(1.5) & ff03$padj <0.05),]
  tmp1 <- list()
  tmp1[["a"]] <- unique(ff01$geneName)
  tmp1[["b"]] <- unique(ff03$X)
  p_01 <- VennDiagram::venn.diagram(x = tmp1, 
                                    filename = NULL, 
                                    main = "main~main~main",
                                    sub  = "sub~sub~sub~sub~sub~sub", 
                                    main.fontfamily = "serif", sub.fontfamily = "serif", 
                                    cex = 1.5, main.cex = 1, sub.cex = 1,
                                    print.mode = "raw",#c("percent", "raw"),
                                    category.names = c("fisherGreater","AefDown"),
                                    fill=c("#FFFFCC","#CCFFFF"),
                                    units = "cm", height = 12, width = 12)
  pdf(file = paste0(path, "PLOT/teVenn.fisherG.down.pdf"))
  grid.draw(p_01)
  dev.off()
  #要“交集基因”的表达
  my02 <- intersect(tmp1[["a"]],tmp1[["b"]])
  expr <- openxlsx::read.xlsx(xlsxFile = "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/Aef/count/Aef.cpm.onlyTE.xlsx")
  expr <- expr[which(expr$theName %in% my02),]
  openxlsx::write.xlsx(x = expr, file = paste0(path, "PLOT/teVenn.fisherG.down-expr.xlsx"))
}
################################################################################
#05、画个venn图 [fisher pvalue < 0.05 + ghyuu KO上调 + Aef KD上调基因]
{
  ff01 <- openxlsx::read.xlsx(xlsxFile = paste0(path,"fisher/Aef-5sFisher_Gene_v1.xlsx"))
  ff01 <- ff01[which(ff01$greate < 0.05),]
  ff02 <- read.csv(file = "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuu/DESeq2/ghyuu.geneOnly.csv")
  ff02 <- ff02[which(ff02$log2FoldChange > log2(1.5) & ff02$padj <0.05),]
  ff03 <- read.csv(file = "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/Aef/DESeq2/Aef.geneOnly.csv")
  ff03 <- ff03[which(ff03$log2FoldChange > log2(1.5) & ff03$padj <0.05),]
  tmp1 <- list()
  tmp1[["a"]] <- unique(ff01$geneName)
  tmp1[["b"]] <- unique(ff02$X)
  tmp1[["c"]] <- unique(ff03$X)
  p_01 <- VennDiagram::venn.diagram(x = tmp1, 
                                    filename = NULL, 
                                    main = "main~main~main",
                                    sub  = "sub~sub~sub~sub~sub~sub", 
                                    main.fontfamily = "serif", sub.fontfamily = "serif", 
                                    cex = 1.5, main.cex = 1, sub.cex = 1,
                                    print.mode = "raw",#c("percent", "raw"),
                                    category.names = c("fisherGreater","ghyuuKO-Up","AefKD-Up"),
                                    fill=c("#FFFFCC","#CCFFFF","#FFCCCC"),
                                    units = "cm", height = 12, width = 12)
  pdf(file = paste0(path, "PLOT/geneVenn.fisherG.ghyuuKO.upup.pdf"))
  grid.draw(p_01)
  dev.off()
  my01 <- intersect(intersect(tmp1[["a"]],tmp1[["b"]]),tmp1[["c"]])
  #要“交集基因”的表达
  my02 <- intersect(tmp1[["a"]],tmp1[["c"]])
  expr <- openxlsx::read.xlsx(xlsxFile = "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/Aef/count/Aef.cpm.onlyGene.xlsx")
  expr <- expr[which(expr$theName %in% my02),]
  openxlsx::write.xlsx(x = expr, file = paste0(path, "PLOT/geneVenn.fisherG.ghyuuKO.upup-expr.xlsx"))
}
################################################################################
#06、画个venn图 [fisher pvalue < 0.05 + ghyuu KO下调 + Aef KD下调基因]
{
  ff01 <- openxlsx::read.xlsx(xlsxFile = paste0(path,"fisher/Aef-5sFisher_Gene_v1.xlsx"))
  ff01 <- ff01[which(ff01$greate < 0.05),]
  ff02 <- read.csv(file = "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuu/DESeq2/ghyuu.geneOnly.csv")
  ff02 <- ff02[which(ff02$log2FoldChange < -log2(1.5) & ff02$padj <0.05),]
  ff03 <- read.csv(file = "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/Aef/DESeq2/Aef.geneOnly.csv")
  ff03 <- ff03[which(ff03$log2FoldChange < -log2(1.5) & ff03$padj <0.05),]
  tmp1 <- list()
  tmp1[["a"]] <- unique(ff01$geneName)
  tmp1[["b"]] <- unique(ff02$X)
  tmp1[["c"]] <- unique(ff03$X)
  p_01 <- VennDiagram::venn.diagram(x = tmp1, 
                                    filename = NULL, 
                                    main = "main~main~main",
                                    sub  = "sub~sub~sub~sub~sub~sub", 
                                    main.fontfamily = "serif", sub.fontfamily = "serif", 
                                    cex = 1.5, main.cex = 1, sub.cex = 1,
                                    print.mode = "raw",#c("percent", "raw"),
                                    category.names = c("fisherGreater","ghyuuKO-Down","AefKD-Down"),
                                    fill=c("#FFFFCC","#CCFFFF","#FFCCCC"),
                                    units = "cm", height = 12, width = 12)
  pdf(file = paste0(path, "PLOT/geneVenn.fisherG.ghyuuKO.down.pdf"))
  grid.draw(p_01)
  dev.off()
  my01 <- intersect(intersect(tmp1[["a"]],tmp1[["b"]]),tmp1[["c"]])
  #要“交集基因”的表达
  my02 <- intersect(tmp1[["a"]],tmp1[["c"]])
  expr <- openxlsx::read.xlsx(xlsxFile = "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/Aef/count/Aef.cpm.onlyGene.xlsx")
  expr <- expr[which(expr$theName %in% my02),]
  openxlsx::write.xlsx(x = expr, file = paste0(path, "PLOT/geneVenn.fisherG.ghyuuKO.down-expr.xlsx"))
}




