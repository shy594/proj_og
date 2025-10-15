#write by Sun Haiayng at 2023.12.19
#01、搭环境
{
  rm(list=ls())
  for (i in c("homologene","data.table","ggridges","ggrepel","KEGG.db",
              "scales","stringr","Biostrings","tidyr","pheatmap","enrichplot",
              "DESeq2","patchwork","ggpubr","ggsci","ggtranscript","rtracklayer",
              "Mfuzz","BiocManager","dplyr","IRanges","GenomicRanges","ReactomePA","ggplot2","clusterProfiler","topGO","Rgraphviz","pathview","org.Mm.eg.db")){library(i, character.only = T,quietly = T)}; rm(i)
}
#02、给索引
{
  path = "/sc/aaSHY/ribo/riboITP/ribo/"
  inx1 = "/Reference/aaSHY/Ensembl99/GRCm38/fa/GRCm38.dna.primary_assembly.fa"
  inx2 = "/Reference/aaSHY/Ensembl99/GRCm38/INDEX/STARv279a"
  inx3 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/Mus_musculus.GRCm38.99.gtf"
  inx4 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf"
  inx5 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf.locInd"
  inx6 = "/Reference/aaSHY/zOther/INDEX/Bowtie2/GRCm38rRNA/rRNA"
  inx7 = "/Reference/aaSHY/Ensembl99/GRCm38/INDEX/bowtie2/mm10"
  inx8 = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/GenBank/rRNA.fa"
}
#03、跑命令
{
  for (i in c("src","fq","Log","stringtie","umi_tools","rRNAdust","bamCoverage","trim","dumpROOM","bam","count")) {
    if (dir.exists(paths = paste0(path, i))==F) {
      dir.create(path = paste0(path,i), recursive = T)
    }
  }
  shell <- paste0(path,"src/run_b.sh")
  cat("#!/bin/bash\n", file = shell)
  prex <- basename(list.files(path = paste0(path, "rawdata"), full.names = T))
  prex <- gsub(".sra","",prex)
  cmd_01 <- paste0("#parallel-fastq-dump -t 10 --split-3 --gzip -s ",path, 
                   "rawdata/", prex, ".sra -O ", path, "fq", "\n")
  cmd_02 <- paste0("#rename s/fastq/fq/ ",path,"fq/*.gz", "\n")
  cmd_03 <- paste0("#umi_tools extract -p ", '"^(?P<umi_1>.{12})(?P<discard_1>.{4}).+$" ',
                   "--extract-method=regex -I ", path, "fq/", prex, ".fq.gz",
                   " -S ", path, "umi_tools/", prex, ".fq.gz", "\n")
  cmd_04 <- paste0("#cutadapt -j 4 -a AAAAAAAAAACAAAAAAAAAA --overlap=4 ",
                   "--trimmed-only -m 8 -o ",path, 
                   "trim/", prex, ".fq ", path, "umi_tools/", prex, ".fq.gz", "\n")
  cmd_05 <- paste0("#rRNAdust -t 4 ", inx8, " ", path, "trim/", prex, ".fq ",
                   "-o ",  path, "rRNAdust/", prex, ".fq; ",
                   "gzip ",path, "rRNAdust/", prex, ".fq", "\n")
  lays <- paste0(path, "rRNAdust/", prex, ".fq.gz")
  cmd_06 <- paste0("#STAR --runThreadN 10 --readFilesCommand zcat ",
                   "--outMultimapperOrder Random --outSAMtype BAM SortedByCoordinate ",
                   "--outFileNamePrefix ",path,"bam/",prex," ",
                   "--outSAMprimaryFlag AllBestScore ",
                   "--genomeDir ",inx2," --readFilesIn ",lays,"\n")
  cmd_07 <- paste0("#mv ",path,"bam/",prex,"Aligned.sortedByCoord.out.bam ",path,
                   "bam/",prex,".bam","\n")
  cmd_08 <- paste0("#umi_tools dedup -I ",path,"bam/",prex,".bam -S ",path,
                   "bam/",prex,".dedup.bam","\n")
  cmd_09 <- paste0("#samtools index ",path,"bam/",prex,".dedup.bam","\n")
  for (i in cmd_01) {cat(i, append = T, file = shell)}
  for (i in cmd_02) {cat(i, append = T, file = shell)}
  for (i in cmd_03) {cat(i, append = T, file = shell)}
  for (i in cmd_04) {cat(i, append = T, file = shell)}
  for (i in cmd_05) {cat(i, append = T, file = shell)}
  for (i in cmd_06) {cat(i, append = T, file = shell)}
  for (i in cmd_07) {cat(i, append = T, file = shell)}
  for (i in cmd_08) {cat(i, append = T, file = shell)}
  for (i in cmd_09) {cat(i, append = T, file = shell)}
  #
  #
  #
  #
  aaaa <- unique(sapply(str_split(prex,"-"),"[",1))
  for (i in aaaa) {
    brex <- prex[grep(i, prex)]
    cmd_10 <- paste0("#samtools merge -@ 10 -f -o ", path,
                     "bam/", i, ".dedup.bam ",
                     paste0(path, "bam/", brex, ".dedup.bam", collapse = " "), "\n")
    cmd_11 <- paste0("#samtools index ", path, "bam/", i, ".dedup.bam", "\n")
    cmd_12 <- paste0("#bamCoverage --normalizeUsing CPM -p 20 --minMappingQuality 1 ",
                     "--samFlagExclude 256 -of bigwig -bs 20 -b ",path,
                     "bam/", i, ".dedup.bam -o ",path,
                     "bamCoverage/",i,".dedup.bw\n")
    for (i in cmd_10) {cat(i, append = T, file = shell)}
    for (i in cmd_11) {cat(i, append = T, file = shell)}
    for (i in cmd_12) {cat(i, append = T, file = shell)}
  }
  cmd_13 <- paste0("TEtranscripts -t ",
                   paste(paste0(path,"bam/",prex[1:22], ".dedup.bam"),collapse = " ")," -c ",
                   paste(paste0(path,"bam/",prex[23:43],".dedup.bam"),collapse = " ")," ",
                   "--GTF ",inx3," --TE ",inx4," --sortByPos --mode multi ",
                   "--project itpS_dedup --outdir ",path,"count", "\n")
  cmd_14 <- paste0("TEtranscripts -t ",
                   paste(paste0(path,"bam/",aaaa[1:3],".dedup.bam"), collapse = " ")," -c ",
                   paste(paste0(path,"bam/",aaaa[4:6],".dedup.bam"), collapse = " ")," ",
                   "--GTF ",inx3," --TE ",inx4," --sortByPos --mode multi ",
                   "--project itpM_dedup --outdir ", path, "count", "\n")
  for (i in cmd_13) {cat(i, append = T, file = shell)}
  for (i in cmd_14) {cat(i, append = T, file = shell)}
  print(paste0("nohup bash ",shell," >",path,"Log/run_b.log"," 2>&1 &"))
}
################################################################################
#04、DESeq2：Gene、TE、TE sites
{
  shy_diff <- function(mark, geneFile, VS){
    tmp1 <- read.table(geneFile, header = T, check.names = F, stringsAsFactors = F, row.names = 1)
    colnames(tmp1) <- sapply(str_split(basename(colnames(tmp1)),pattern = ".dedup.bam"), "[", 1)
    tmp1 <- tmp1[,grep(paste0(VS[2],"|",VS[3]), colnames(tmp1))]
    Class <- c(colnames(tmp1)[grep(VS[2],colnames(tmp1))],
               colnames(tmp1)[grep(VS[3],colnames(tmp1))])
    Class <- factor(sapply(str_split(Class,"-"),"[",1))
    dfclas <- data.frame(row.names = colnames(tmp1), Class, stringsAsFactors = T)
    geneAA <- tmp1[grep(rownames(tmp1), pattern = ":", invert = T),]
    teOnly <- tmp1[grep(rownames(tmp1), pattern = ":", invert = F),]
    teGene <- tmp1
    #
    geneAA <- DESeq2::DESeqDataSetFromMatrix(countData = geneAA, colData = dfclas, design = ~Class)
    #geneAA <- geneAA[rowSums(BiocGenerics::counts(geneAA) > 2) >= colLen, ]
    geneAA <- geneAA[rowMeans(BiocGenerics::counts(geneAA)) > 5,]
    geneAA <- DESeq2::DESeq(geneAA,)
    geneAA <- DESeq2::results(geneAA, contrast = VS)
    write.csv(geneAA, paste0(path,"DESeq2/",mark, ".geneOnly.csv"))
    upup <- subset(geneAA, log2FoldChange > log2(1.5) & padj < 0.05)
    down <- subset(geneAA, log2FoldChange < -log2(1.5)& padj < 0.05)
    write.csv(upup,paste0(path,"DESeq2/",mark, ".geneOnly.upup.csv"))
    write.csv(down,paste0(path,"DESeq2/",mark, ".geneOnly.down.csv"))
    #
    teOnly <- DESeq2::DESeqDataSetFromMatrix(countData = teOnly, colData = dfclas, design = ~Class)
    #teOnly <- teOnly[rowSums(BiocGenerics::counts(teOnly) > 2) >= colLen, ]#[rowMeans(BiocGenerics::counts(teOnly)) > 1,]
    teOnly <- teOnly[rowMeans(BiocGenerics::counts(teOnly)) > 5,]
    teOnly <- DESeq2::DESeq(teOnly,)
    teOnly <- DESeq2::results(teOnly, contrast = VS)
    write.csv(teOnly, paste0(path,"DESeq2/",mark, ".teOnly.csv"))
    upup <- subset(teOnly, log2FoldChange > log2(1.5) & padj < 0.05)
    down <- subset(teOnly, log2FoldChange < -log2(1.5)& padj < 0.05)
    write.csv(upup,paste0(path,"DESeq2/",mark, ".teOnly.upup.csv"))
    write.csv(down,paste0(path,"DESeq2/",mark, ".teOnly.down.csv"))
    #
    teGene <- DESeq2::DESeqDataSetFromMatrix(countData = teGene, colData = dfclas, design = ~Class)
    #teGene <- teGene[rowSums(BiocGenerics::counts(teGene) > 2) >= colLen, ]#[rowMeans(BiocGenerics::counts(teGene)) > 1,]
    teGene <- teGene[rowMeans(BiocGenerics::counts(teGene)) > 5,]
    teGene <- DESeq2::DESeq(teGene,)
    teGene <- DESeq2::results(teGene, contrast = VS)
    teGene <- teGene[grep(rownames(teGene), pattern = ":", invert = F),]
    write.csv(teGene, paste0(path,"DESeq2/",mark, ".te.GeneBackGround.csv"))
    upup <- subset(teGene, log2FoldChange > log2(1.5) & padj < 0.05)
    down <- subset(teGene, log2FoldChange < -log2(1.5)& padj < 0.05)
    write.csv(upup,paste0(path,"DESeq2/",mark, ".te.GeneBackGround.upup.csv"))
    write.csv(down,paste0(path,"DESeq2/",mark, ".te.GeneBackGround.down.csv"))
  }
  shy_diff(mark = "four_twos",
           geneFile = paste0(path, "count/itpS_dedup.cntTable"),
           VS = c("Class","four","twos"))
}
################################################################################
#05、TE的翻译效率
###mark <- "twos____zygote"
###red1 <- read.csv("/sc/aaSHY/ribo/riboITP/rnaSeq/DESeq2/twos_zygote.te.GeneBackGround.csv")
###red2 <- read.csv("/sc/aaSHY/ribo/riboITP/ribo/DESeq2/twos_zygote.te.GeneBackGround.csv")
{
  mark <- "four____twos"
  red1 <- read.csv("/sc/aaSHY/ribo/riboITP/rnaSeq/DESeq2/four_twos.te.GeneBackGround.csv")
  red2 <- read.csv("/sc/aaSHY/ribo/riboITP/ribo/DESeq2/four_twos.te.GeneBackGround.csv")
  colnames(red1)[3] <- "log2FC_RNA"
  colnames(red2)[3] <- "log2FC_Ribo"
  red3 <- left_join(red1,red2[,c(1,3)],by="X")
  #
  red4 <- na.omit(red3)
  red4$colo <- "heise"
  red4$colo[which(abs(red4$log2FC_RNA) < 1)] <- "huise"
  red4$colo[which(abs(red4$log2FC_Ribo) < 1)] <- "huise"
  labb <- red4[grep("MERVL-int|MT2_Mm",red4$X),]
  p_002 <- ggplot(data = red4) +
    geom_point(aes(x = log2FC_RNA, y = log2FC_Ribo, color = colo)) +
    geom_point(data = labb, aes(x = log2FC_RNA, y = log2FC_Ribo), color = "red") +
    scale_color_manual(values = c("black","gray"), labels = c("significant","not")) +
    geom_hline(aes(yintercept = -1), linetype="dashed", colour="royalblue") +
    geom_hline(aes(yintercept =  1), linetype="dashed", colour="royalblue") +
    geom_vline(aes(xintercept =  1), linetype="dashed", colour="royalblue") +
    geom_vline(aes(xintercept = -1), linetype="dashed", colour="royalblue") +
    ###scale_x_continuous(limits = c(floor(-max(red4[,c(3,8)])),ceiling(max(red4[,c(3,8)])))) +
    ###scale_y_continuous(limits = c(floor(-max(red4[,c(3,8)])),ceiling(max(red4[,c(3,8)])))) +
    scale_x_continuous(limits = c(floor(-max(abs(red4[,c(3,8)]))),ceiling(max(abs(red4[,c(3,8)]))))) +
    scale_y_continuous(limits = c(floor(-max(abs(red4[,c(3,8)]))),ceiling(max(abs(red4[,c(3,8)]))))) +
    geom_text_repel(data = labb, aes(x = log2FC_RNA, y = log2FC_Ribo, label = X)) +
    theme_few() +
    theme(plot.margin = margin(unit = "cm",1,1,1,1)) +
    labs(title = mark)
  p_002
  ggsave(plot = p_002,
         filename = "/sc/aaSHY/ribo/riboITP/ribo/PLOT/four__twos.TE_v2.pdf",
         units = "cm", height = 15, width = 20)
  ###write.csv(x = red4, quote = F,row.names = F,
  ###          file = "/sc/aaSHY/ribo/riboITP/ribo/DESeq2/twos_zygote.te.GeneBackGround.FC.csv")
  write.csv(x = red4, quote = F,row.names = F,
            file = "/sc/aaSHY/ribo/riboITP/ribo/DESeq2/four_twos.te.GeneBackGround.FC.csv")
}
################################################################################
#06、基因的翻译效率 [最后说, 那算了]
####与MERVL发生融合的基因
{
  red1 <- import(inx3); red1 <- red1[which(red1$type=="gene")]
  red2 <- import(inx4); red2 <- red2[grep("MERVL-int|MT2_Mm",red2$gene_id)]
  red3_1 <- import("/ChIP_seq_1/aaSHY/rna2014Deng/myRun/stringtie/twoEarly.merge.gtf")
  red3_2 <- import("/ChIP_seq_1/aaSHY/rna2014Deng/myRun/stringtie/twoMid.merge.gtf")
  red3_3 <- import("/ChIP_seq_1/aaSHY/rna2014Deng/myRun/stringtie/twoLate.merge.gtf")
  red3 <- c(red3_1,red3_2,red3_3)
  ff_exon <- red3[which(red3$type=="exon")]
  ff_tran <- red3[which(red3$type=="transcript")]
  #
  ov_1 <- suppressWarnings(findOverlaps(query = red2, subject = ff_exon, ignore.strand=T))
  ff_tran_ID <- ff_exon$transcript_id[unique(subjectHits(ov_1))]
  ff_tran_RR <- ff_tran[which(ff_tran$transcript_id %in% ff_tran_ID)]
  ov_2 <- suppressWarnings(findOverlaps(query = red1, subject = ff_tran_RR, ignore.strand=T))
  gggg <- unique(red1$gene_name[unique(queryHits(ov_2))])
}
#
{
  intersect(c("jjuuzz","Ddit4l"), red2$X)
  mark <- "four____twos"
  red1 <- read.csv("/sc/aaSHY/ribo/riboITP/rnaSeq/DESeq2/four_twos.geneOnly.csv")
  red2 <- read.csv("/sc/aaSHY/ribo/riboITP/ribo/DESeq2/four_twos.geneOnly.csv")
  colnames(red1)[3] <- "log2FC_RNA"
  colnames(red2)[3] <- "log2FC_Ribo"
  red3 <- left_join(red1,red2[,c(1,3)],by="X")
  #
  red4 <- na.omit(red3)
  red4$colo <- "heise"
  red4$colo[which(red4$X %in% gggg)] <- "huise"
  #red4$colo[which(abs(red4$log2FC_RNA) < 1)] <- "huise"
  #red4$colo[which(abs(red4$log2FC_Ribo) < 1)] <- "huise"
  labb <- red4[grep("jjuuzz|Ddit4l",red4$X),]
  p_002 <- ggplot(data = red4) +
    geom_point(aes(x = log2FC_RNA, y = log2FC_Ribo, color = colo)) +
    geom_point(data = labb, aes(x = log2FC_RNA, y = log2FC_Ribo), color = "red") +
    scale_color_manual(values = c("black","cyan"), labels = c("others","chimeric")) +
    geom_hline(aes(yintercept = -1), linetype="dashed", colour="royalblue") +
    geom_hline(aes(yintercept =  1), linetype="dashed", colour="royalblue") +
    geom_vline(aes(xintercept =  1), linetype="dashed", colour="royalblue") +
    geom_vline(aes(xintercept = -1), linetype="dashed", colour="royalblue") +
    ###scale_x_continuous(limits = c(floor(-max(red4[,c(3,8)])),ceiling(max(red4[,c(3,8)])))) +
    ###scale_y_continuous(limits = c(floor(-max(red4[,c(3,8)])),ceiling(max(red4[,c(3,8)])))) +
    scale_x_continuous(limits = c(floor(-max(abs(red4[,c(3,8)]))),ceiling(max(abs(red4[,c(3,8)]))))) +
    scale_y_continuous(limits = c(floor(-max(abs(red4[,c(3,8)]))),ceiling(max(abs(red4[,c(3,8)]))))) +
    geom_text_repel(data = labb, aes(x = log2FC_RNA, y = log2FC_Ribo, label = X)) +
    theme_few() +
    theme(plot.margin = margin(unit = "cm",1,1,1,1)) +
    labs(title = mark)
  p_002
  ggsave(plot = p_002,
         filename = "/sc/aaSHY/ribo/riboITP/ribo/PLOT/four__twos.TE_v2.pdf",
         units = "cm", height = 15, width = 20)
  write.csv(x = red4, quote = F,row.names = F,
            file = "/sc/aaSHY/ribo/riboITP/ribo/DESeq2/four_twos.te.GeneBackGround.FC.csv")
}
####为什么jjuuzz没有在Ribo-seq中被检测到? 因为它的reads在样本中非常少, 在DESeq2时被过滤掉了
{
  path = "/sc/aaSHY/ribo/riboITP/ribo/"
  path = "/sc/aaSHY/ribo/riboITP/rnaSeq/"
  tmp1 <- read.table(paste0(path, "count/itp_dedup.cntTable"), header = T, sep = "\t", row.names = 1)
  colnames(tmp1) <- paste0(sapply(str_split(colnames(tmp1),"\\."),"[",8),"-",
                           sapply(str_split(colnames(tmp1),"\\."),"[",9))
  tmp1["jjuuzz",grep("two|four",colnames(tmp1))]
}
#


