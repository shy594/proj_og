#write by Sun Haiayng at 2024.08.30
#01、搭环境
{
  rm(list=ls())
  for (i in c("tidyr","data.table","ggridges","ggrepel","scales","stringr","Biostrings","pheatmap","enrichplot","DESeq2","patchwork","ggpubr","ggsci","ggtranscript","rtracklayer","Mfuzz","BiocManager","dplyr","IRanges","GenomicRanges","ReactomePA","ggplot2","clusterProfiler","topGO","Rgraphviz","pathview","org.Mm.eg.db")){print(i);library(i, character.only = T,quietly = T)}; rm(i)
}
#02、给索引
{
  path = "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/thehuhuhu/"
  inx1 = "/Reference/aaSHY/Ensembl99/GRCm38/INDEX/STARv279a"
  inx3 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/Mus_musculus.GRCm38.99.gtf"
  inx4 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf"
  inx5 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf.locInd"
  inx6 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.saf"
}
#03、跑命令
{
  for (i in c("src","fq","Log","stringtie/count","bamCoverage","trim","DESeq2","dumpROOM","bam","count","PLOT")) {
    if (dir.exists(paths = paste0(path, i))==F) {
      dir.create(path = paste0(path,i), recursive = T)
    }
  }
  shell <- paste0(path,"src/run_a.sh")
  cat("#!/bin/bash\n", file = shell)
  ff01 <- list.files(path = paste0(path,"fq"),pattern = "_1.fq.gz")
  prex <- sapply(str_split(ff01,pattern = "_1.fq.gz"),"[",1)
  lay1 <- paste0(path,"fq/",prex,"_1.fq.gz")
  lay2 <- paste0(path,"fq/",prex,"_2.fq.gz")
  cmd_01 <- paste0("trim_galore --phred33 -j 2 --fastqc ",
                   "-o ", path, "trim --paired ",lay1," ",lay2,"\n")
  lay1 <- paste0(path,"trim/",prex,"_1_val_1.fq.gz")
  lay2 <- paste0(path,"trim/",prex,"_2_val_2.fq.gz")
  cmd_02 <- paste0("STAR --runThreadN 4 --readFilesCommand zcat ",
                   "--outMultimapperOrder Random --outSAMtype BAM SortedByCoordinate ",
                   "--outFileNamePrefix ",paste0(path,"bam/",prex)," ",
                   "--outSAMprimaryFlag AllBestScore --twopassMode Basic ",
                   "--genomeDir ",inx1," --readFilesIn ",lay1," ",lay2,"\n")
  cmd_03 <- paste0("mv ",paste0(path,"bam/",prex,"Aligned.sortedByCoord.out.bam "),paste0(path,"bam/",prex,".bam"),"\n")
  cmd_04 <- paste0("samtools index ",paste0(path,"bam/",prex,".bam"),"\n")
  cmd_05 <- paste0("bamCoverage --normalizeUsing CPM -p 10 --minMappingQuality 1 --samFlagExclude 256 -of bigwig -bs 20 -b ",
                   paste0(path,"bam/",prex,".bam")," -o ",paste0(path,"bamCoverage/",prex,".bw"),"\n")
  cmd_06 <- paste0("TEtranscripts -t ",
                   paste(paste0(path,"bam/",prex[1:10], ".bam"),collapse = " ")," -c ",
                   paste(paste0(path,"bam/",prex[11:20],".bam"),collapse = " "),
                   " --GTF ",inx3," --TE ",inx4," --sortByPos --mode multi ",
                   "--project huhuhu --outdir ",paste0(path,"count"),"\n")
  cmd_07 <- paste0("TElocal -b ",
                   paste0(path,"bam/",prex,".bam"),
                   " --GTF ",inx3," --TE ",inx5," --sortByPos --project ",
                   paste0(path, "count/",prex),"\n")
  cmd_08 <- paste0("stringtie ", paste0(path,"bam/",prex,".bam")," ",
                   "-p 4 -G ", inx3, " -o ", paste0(path,"stringtie/",prex,".gtf"),"\n")
  cmd_09 <- paste0("stringtie --merge -o ",
                   path,"stringtie/huhuhu.merge.gtf -G ",inx3," ",
                   paste(paste0(path,"stringtie/",prex,".gtf"), collapse = " "),"\n")
  cmd_10 <- paste0("cat ", paste0(path,"stringtie/huhuhu.merge.gtf "), 
                   "|sed 's/gene_name/gene_abcd/g; s/\\tgene_id/\\tgene_name/g' >",
                   paste0(path,"stringtie/huhuhu.merge.R.gtf"),"\n")
  cmd_11 <- paste0("TEtranscripts ",
                   "-t ",paste(paste0(path,"bam/",prex[1:10], ".bam"),collapse = " ")," ",
                   "-c ",paste(paste0(path,"bam/",prex[11:20],".bam"),collapse = " ")," ",
                   "--GTF ",path,"stringtie/huhuhu.merge.R.gtf ",
                   "--TE ", inx4," --sortByPos --mode multi --project huhuhu ",
                   "--outdir ",path,"stringtie/count","\n")
  for (i in cmd_01) {cat(i, file = shell, append = T)}
  for (i in cmd_02) {cat(i, file = shell, append = T)}
  for (i in cmd_03) {cat(i, file = shell, append = T)}
  for (i in cmd_04) {cat(i, file = shell, append = T)}
  for (i in cmd_05) {cat(i, file = shell, append = T)}
  for (i in cmd_06) {cat(i, file = shell, append = T)}
  for (i in cmd_07) {cat(i, file = shell, append = T)}
  for (i in cmd_08) {cat(i, file = shell, append = T)}
  for (i in cmd_09) {cat(i, file = shell, append = T)}
  for (i in cmd_10) {cat(i, file = shell, append = T)}
  for (i in cmd_11) {cat(i, file = shell, append = T)}
  #
  #two-pass mode是否有影响
  cmd_12 <- paste0("#STAR --runThreadN 10 --readFilesCommand zcat ",
                   "--outMultimapperOrder Random --outSAMtype BAM SortedByCoordinate ",
                   "--outFileNamePrefix ",path,"dumpROOM/",prex," ",
                   "--outSAMprimaryFlag AllBestScore ",
                   "--genomeDir ",inx1," --readFilesIn ",lay1," ",lay2,"\n")
  cmd_13 <- paste0("#mv ",
                   paste0(path,"dumpROOM/",prex,"Aligned.sortedByCoord.out.bam "),
                   paste0(path,"dumpROOM/",prex,".bam"),"\n")
  cmd_14 <- paste0("#samtools index ",paste0(path,"dumpROOM/",prex,".bam"),"\n")
  cmd_15 <- paste0("#TEtranscripts -t ",
                   paste(paste0(path,"dumpROOM/",prex[1:2],".bam"),collapse = " ")," -c ",
                   paste(paste0(path,"dumpROOM/",prex[3:4],".bam"),collapse = " "),
                   " --GTF ",inx3," --TE ",inx4," --sortByPos --mode multi ",
                   "--project huhuhu --outdir ",paste0(path,"dumpROOM"),"\n")
  for (i in cmd_12) {cat(i, file = shell, append = T)}
  for (i in cmd_13) {cat(i, file = shell, append = T)}
  for (i in cmd_14) {cat(i, file = shell, append = T)}
  for (i in cmd_15) {cat(i, file = shell, append = T)}
  print(paste0("nohup bash ",shell, " >",paste0(path,"Log/"),basename(shell),".log"," 2>&1 &"))
}
################################################################################
#05、CPM Excel
{
  cnt1 <- read.table(paste0(path,"count/huhuhu.cntTable"), header = T, sep = "\t", row.names = 1)
  colnames(cnt1) <- prex
  tmp1 <- cnt1[grep(x = rownames(cnt1), pattern=":", invert = F),]#TE
  tmp2 <- cnt1[grep(x = rownames(cnt1), pattern=":", invert = T),]#Gene
  tmp3 <- cnt1#All
  tmp1 <- as.data.frame(apply(tmp1, 2, function(x){x/sum(x)*1000000}))
  tmp2 <- as.data.frame(apply(tmp2, 2, function(x){x/sum(x)*1000000}))
  tmp3 <- as.data.frame(apply(tmp3, 2, function(x){x/sum(x)*1000000}))
  tmp1$theName <- rownames(tmp1)
  tmp2$theName <- rownames(tmp2)
  tmp3$theName <- rownames(tmp3)
  openxlsx::write.xlsx(x = tmp1, file = paste0(path, "count/huhuhu.cpm.onlyTE.xlsx"))
  openxlsx::write.xlsx(x = tmp2, file = paste0(path, "count/huhuhu.cpm.onlyGene.xlsx"))
  openxlsx::write.xlsx(x = tmp3, file = paste0(path, "count/huhuhu.cpm.allGeneTE.xlsx"))
}
################################################################################
#06、重命名文件列名
#CPM Excel
{
  tmps <- read.table("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/thehuhuhu/src/rename.txt",header = F, sep = "\t")
  tttt <- list.files(paste0(path,"count"),"xlsx",full.names = T)
  for (i in seq_along(tttt)) {
    print(i)
    tmp1 <- openxlsx::read.xlsx(tttt[i])
    tmp1 <- tmp1[,c("theName",tmps$V1)]
    colnames(tmp1) <- c("theName",tmps$V2)
    openxlsx::write.xlsx(x = tmp1,
                         file = gsub("xlsx","rename.xlsx",tttt[i]))
  }
}
#TEtranscripts
{
  tmps <- read.table("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/thehuhuhu/src/rename.txt",header = F, sep = "\t")
  tmp1 <- read.table(file = paste0(path,"count/huhuhu.cntTable"),
                     sep = "\t", header = T, row.names = 1)
  colnames(tmp1) <- gsub("X.ChIP_seq_2.aaSHY.ghyuu.rnaSEQ.thehuhuhu.bam.(.*).bam.*","\\1",colnames(tmp1))
  colnames(tmp1) <- gsub("\\.","-",colnames(tmp1))
  tmp1$genesName <- rownames(tmp1)
  tmp1 <- tmp1[,c("genesName",tmps$V1)]
  colnames(tmp1) <- c("genesName",tmps$V2)
  write.table(x = tmp1, file = paste0(path,"count/huhuhu.rename.cntTable"),
              sep = "\t", quote = F, col.names = T, row.names = F)
}
#TElocal
{
  paste0("paste ",paste0(path,"countUniq/",tmps$V1,".cntTable", collapse = " ")," |",
         "cut -f 1,",paste0(seq(2,40,2), collapse = ",")," >",path,
         "countUniq/","huhuhu.site.cntTable","\n")
  tmps <- read.table("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/thehuhuhu/src/rename.txt",header = F, sep = "\t")
  tmp1 <- read.table(file = paste0(path,"countUniq/huhuhu.site.cntTable"),
                     sep = "\t", header = T, row.names = 1)
  colnames(tmp1) <- gsub("X.ChIP_seq_2.aaSHY.ghyuu.rnaSEQ.thehuhuhu.bam.(.*).bam.*","\\1",colnames(tmp1))
  colnames(tmp1) <- gsub("\\.","-",colnames(tmp1))
  tmp1$genesName <- rownames(tmp1)
  tmp1 <- tmp1[,c("genesName",tmps$V1)]
  colnames(tmp1) <- c("genesName",tmps$V2)
  write.table(x = tmp1, file = paste0(path,"countUniq/huhuhu.site.rename.cntTable"),
              sep = "\t", quote = F, col.names = T, row.names = F)
}
#file rename
{
  namess <- read.table(paste0(path,"src/rename.txt"),header = F)
  aa <- paste0(path,"bam/",namess$V1,".bam.bai")
  bb <- paste0(path,"bam/",namess$V2,".bam.bai")
  file.rename(from = aa, to = bb)
  rm(namess,aa,bb)
}
################################################################################
#07、PCA
{
  #----------------------------------PCA-batch----------------------------------
  cnt1 <- read.table(paste0(path, "count/huhuhu.rename.cntTable"),
                     header = T, check.names = F, stringsAsFactors = F, row.names = 1)
  cnt1 <- cnt1[grep(rownames(cnt1), pattern = ":", invert = T),]
  #colnames(cnt1) <- sapply(str_split(basename(colnames(cnt1)), ".bam"),"[",1)
  #
  #base::load(file = "/Reference/aaSHY/DatabaseFiles/forPCA.RData")
  #myta <- cbind(cnt1, forPCA)
  #
  myta <- cnt1
  condition = str_sub(colnames(cnt1),1,6)
  btch = c(rep("ours",length(colnames(cnt1))))
  clas <- data.frame(names = colnames(myta), condition = condition, batch = btch)
  dds1 <- suppressWarnings(DESeq2::DESeqDataSetFromMatrix(myta, colData = clas, design = ~condition))#[rowSums(counts(vtax))>5]
  dds2 <- BiocGenerics::estimateSizeFactors(dds1)
  rldf <- DESeq2::rlog(dds2)
  pd_1 <- plotPCA(rldf, intgroup = c("condition"), returnData = T)
  rldf <- DESeq2::vst(dds2)
  pd_2 <- plotPCA(rldf, intgroup = c("condition"), returnData = T)
  #去除批次效应
  #SummarizedExperiment::assay(rldf) <- limma::removeBatchEffect(assay(rldf), batch=rldf$batch)
  #https://cloud.tencent.com/developer/article/1625223
  #豆包说“通常情况下，数据集小于30个样品时可以用rlog，数据集大于30个样品时用vst
  #vsdf <- vst(dds2)
  #rawx <- SummarizedExperiment(counts(dds2, normalized = F), colData=colData(dds2))
  #norx <- SummarizedExperiment(counts(dds2, normalized = T), colData=colData(dds2))
  #pd_3 <- plotPCA(DESeqTransform(rawx), intgroup = c("condition"), returnData = T)
  #pd_4 <- plotPCA(DESeqTransform(norx), intgroup = c("condition"), returnData = T)
  shy_linshi <- function(df, i){
    mark <- c("rlog", "vst")[i]
    perc <- round(100 * attr(df, "percentVar"))
    df$condition <- factor(as.character(df$condition), levels = as.character(df$condition),
                           labels = as.character(df$condition))
    df$repsnames <- rownames(df)
    ppp <<- ggplot(df, aes(PC1, PC2, color=condition)) +
      geom_point(size=5) +
      xlab(paste0("PC1: ",perc[1],"% variance")) +
      ylab(paste0("PC2: ",perc[2],"% variance")) +
      scale_color_d3() +
      theme_classic() +
      labs(color="Cell type",title = paste0("TITLE: plot of PCA ... ... ... ... ... ...",mark)) +
      theme(legend.text  = element_text(family = "serif", size = 12),
            legend.title = element_text(family = "serif", size = 12),
            legend.background = element_rect(fill = "gray95"),
            legend.position = c(0.8,0.50),
            axis.text    = element_text(family = "serif", size = 12),
            axis.title   = element_text(family = "serif", size = 12)) +
    geom_text_repel(aes(label = repsnames))
    print(ppp)
  }
  pdf(file = paste0(path,"PLOT/PCA.pdf"),width = 10, height = 8)
  shy_linshi(df = pd_1, i = 1)
  shy_linshi(df = pd_2, i = 2)
  dev.off()
}
################################################################################
#08、DESeq2
{
  shy_diff <- function(mark, grou,geneFile,siteFile, Class, VS){
    tmp1 <- read.table(geneFile, header = T, check.names = F, stringsAsFactors = F, row.names = 1)
    tmp1 <- tmp1[,grep(grou,colnames(tmp1))]
    dfclas <- data.frame(row.names = colnames(tmp1), Class, stringsAsFactors = T)
    colLen <- ceiling(length(colnames(tmp1))/2)
    geneAA <- tmp1[grep(rownames(tmp1), pattern = ":", invert = T),]
    teOnly <- tmp1[grep(rownames(tmp1), pattern = ":", invert = F),]
    teGene <- tmp1
    #
    geneAA <- DESeq2::DESeqDataSetFromMatrix(countData = geneAA, colData = dfclas, design = ~Class)
    #geneAA <- geneAA[rowSums(BiocGenerics::counts(geneAA) > 2) >= colLen, ]
    geneAA <- geneAA[rowMeans(BiocGenerics::counts(geneAA)) > 5,]
    geneAA <- DESeq2::DESeq(geneAA,)
    geneAA <- DESeq2::results(geneAA, contrast = VS)
    write.csv(geneAA, paste0(path,"DESeq2/",mark,".",grou,".geneOnly.csv"))
    upup <- subset(geneAA, log2FoldChange > log2(1.5) & padj < 0.05)
    down <- subset(geneAA, log2FoldChange < -log2(1.5)& padj < 0.05)
    write.csv(upup,paste0(path,"DESeq2/",mark,".",grou,".geneOnly.upup.csv"))
    write.csv(down,paste0(path,"DESeq2/",mark,".",grou,".geneOnly.down.csv"))
    #
    teOnly <- DESeq2::DESeqDataSetFromMatrix(countData = teOnly, colData = dfclas, design = ~Class)
    #teOnly <- teOnly[rowSums(BiocGenerics::counts(teOnly) > 2) >= colLen, ]#[rowMeans(BiocGenerics::counts(teOnly)) > 1,]
    teOnly <- teOnly[rowMeans(BiocGenerics::counts(teOnly)) > 5,]
    teOnly <- DESeq2::DESeq(teOnly,)
    teOnly <- DESeq2::results(teOnly, contrast = VS)
    write.csv(teOnly, paste0(path,"DESeq2/",mark,".",grou,".teOnly.csv"))
    upup <- subset(teOnly, log2FoldChange > log2(1.5) & padj < 0.05)
    down <- subset(teOnly, log2FoldChange < -log2(1.5)& padj < 0.05)
    write.csv(upup,paste0(path,"DESeq2/",mark,".",grou,".teOnly.upup.csv"))
    write.csv(down,paste0(path,"DESeq2/",mark,".",grou,".teOnly.down.csv"))
    #
    teGene <- DESeq2::DESeqDataSetFromMatrix(countData = teGene, colData = dfclas, design = ~Class)
    #teGene <- teGene[rowSums(BiocGenerics::counts(teGene) > 2) >= colLen, ]#[rowMeans(BiocGenerics::counts(teGene)) > 1,]
    teGene <- teGene[rowMeans(BiocGenerics::counts(teGene)) > 5,]
    teGene <- DESeq2::DESeq(teGene,)
    teGene <- DESeq2::results(teGene, contrast = VS)
    teGene <- teGene[grep(rownames(teGene), pattern = ":", invert = F),]
    write.csv(teGene, paste0(path,"DESeq2/",mark,".",grou,".te.GeneBackGround.csv"))
    upup <- subset(teGene, log2FoldChange > log2(1.5) & padj < 0.05)
    down <- subset(teGene, log2FoldChange < -log2(1.5)& padj < 0.05)
    write.csv(upup,paste0(path,"DESeq2/",mark,".",grou,".te.GeneBackGround.upup.csv"))
    write.csv(down,paste0(path,"DESeq2/",mark,".",grou,".te.GeneBackGround.down.csv"))
    #
    #
    tmp2 <- read.table(siteFile, header = T, check.names = F, stringsAsFactors = F, row.names = 1)
    tmp2 <- tmp2[,grep(grou,colnames(tmp2))]
    colLen <- ceiling(length(colnames(tmp2))/2)
    #
    teGene <- tmp2
    teGene <- DESeq2::DESeqDataSetFromMatrix(countData = teGene, colData = dfclas, design = ~Class)
    #teGene <- teGene[rowSums(BiocGenerics::counts(teGene) > 2) >= colLen, ]#[rowMeans(BiocGenerics::counts(teGene)) > 1,]
    teGene <- teGene[rowMeans(BiocGenerics::counts(teGene)) > 5,]
    teGene <- DESeq2::DESeq(teGene,)
    teGene <- DESeq2::results(teGene, contrast = VS)
    teGene <- teGene[grep(rownames(teGene), pattern = ":", invert = F),]
    write.csv(teGene, paste0(path,"DESeq2/",mark,".",grou,".te.site.GeneBackGround.csv"))
    upup <- subset(teGene, log2FoldChange > log2(1.5) & padj < 0.05)
    down <- subset(teGene, log2FoldChange < -log2(1.5)& padj < 0.05)
    write.csv(upup, paste0(path,"DESeq2/",mark,".",grou,".te.site.GeneBackGround.upup.csv"))
    write.csv(down, paste0(path,"DESeq2/",mark,".",grou,".te.site.GeneBackGround.down.csv"))
    #
    teGene <- tmp2[grep(rownames(tmp2), pattern = ":", invert = F),]
    teGene <- DESeq2::DESeqDataSetFromMatrix(countData = teGene, colData = dfclas, design = ~Class)
    #teGene <- teGene[rowSums(BiocGenerics::counts(teGene) > 2) >= colLen, ]#[rowMeans(BiocGenerics::counts(teGene)) > 1,]
    teGene <- teGene[rowMeans(BiocGenerics::counts(teGene)) > 5,]
    teGene <- DESeq2::DESeq(teGene,)
    teGene <- DESeq2::results(teGene, contrast = VS)
    write.csv(teGene, paste0(path,"DESeq2/",mark,".",grou,".te.site.teOnly.csv"))
    upup <- subset(teGene, log2FoldChange > log2(1.5) & padj < 0.05)
    down <- subset(teGene, log2FoldChange < -log2(1.5)& padj < 0.05)
    write.csv(upup, paste0(path,"DESeq2/",mark,".",grou,".te.site.teOnly.upup.csv"))
    write.csv(down, paste0(path,"DESeq2/",mark,".",grou,".te.site.teOnly.down.csv"))
  }
  sets <- paste0("group",seq(5))
  for (grou in sets) {
    print(grou)
    shy_diff(mark = "huhuhu", grou,
             geneFile = paste0(path, "count/huhuhu.rename.cntTable"),
             siteFile = paste0(path, "count/huhuhu.site.rename.cntTable"),
             Class = factor(c("OP","OP","WT","WT")), VS = c("Class","OP","WT"))
  }
}
################################################################################
#08、DESeq2, 给TE copies补上坐标
{
  ff01 <- list.files(paste0(path,"DESeq2"),pattern = "te.site", full.names = T)
  tf01 <- import(inx4)
  tf01 <- as.data.frame(tf01)
  tf01$positions <- paste0(tf01$seqnames,":",tf01$start,"-",tf01$end)
  #
  for (i in seq_along(ff01)) {
    print(i)
    tmp1 <- read.csv(ff01[i])
    tmp1$transcript_id <- sapply(str_split(tmp1$X,":"),"[",1)
    tmp2 <- dplyr::left_join(tmp1,tf01[,c("transcript_id","positions")],by="transcript_id")
    write.csv(tmp2,file = gsub("site","site.coord",ff01[i]),quote = F, row.names = F)
  }; rm(tmp1,ff01,i,sub1)
  tf01 <- import("/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf")
}
################################################################################
#09、pheatmap
{
  shy_pheatmap <- function(mark,grou,read_cntFile,gene_resFile,tete_resFile,savePath){
    print("reads counted by TEtranscripts")
    #
    data <- read.table(read_cntFile, header = T, row.names = 1, check.names = F, stringsAsFactors = F)
    data <- data[,grep(grou,colnames(data))]
    gene <- data[grep(rownames(data), pattern = ":", invert = T),]
    gene <- gene[rowSums(gene) > 5,]
    gene <- as.data.frame(apply(gene, 2, function(x){x/sum(x)*1000000}))
    tete <- data[grep(rownames(data), pattern = ":", invert = F),]
    tete <- tete[rowSums(tete) > 5,]
    tete <- as.data.frame(apply(tete, 2, function(x){x/sum(x)*1000000}))
    #
    resG <- as.data.frame(read.csv(gene_resFile, header = T, stringsAsFactors = F))
    resG <- na.omit(resG)
    sigG <- resG[which(abs(resG$log2FoldChange)>log2(1.5) & resG$padj<0.05), "X"]
    sigG_upup <- resG[which(resG$log2FoldChange>log2(1.5) & resG$padj<0.05), "X"]
    sigG_down <- resG[which(resG$log2FoldChange < -log2(1.5) & resG$padj<0.05), "X"]
    resG <- as.data.frame(read.csv(tete_resFile, header = T, stringsAsFactors = F))
    resG <- na.omit(resG)
    sigT <- resG[which(abs(resG$log2FoldChange)>log2(1.5) & resG$padj<0.05), "X"]
    sigT_upup <- resG[which(resG$log2FoldChange>log2(1.5) & resG$padj<0.05), "X"]
    sigT_down <- resG[which(resG$log2FoldChange < -log2(1.5) & resG$padj<0.05), "X"]
    ##开始画图
    pdf(paste0(savePath,mark,".",grou,".pheatmap.pdf"), width = 8, height = 10)
    ##gene表达量热图
    pheatmap(mat = gene,
             main="the heatmap of all Genes",
             scale = "row", cellwidth = 60, 
             show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = T,
             labels_col = colnames(gene), angle_col = "45",
             breaks = seq(-2,2,by=0.01), border = F,
             color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                        colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
    ##siginificant gene表达量热图
    pheatmap(mat = gene[which(rownames(gene) %in% sigG),],
             main="the heatmap of significant Genes",
             scale = "row", cellwidth = 60, 
             show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = T,
             labels_col = colnames(gene), angle_col = "45",
             breaks = seq(-2,2,by=0.01), border = F,
             color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                        colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
    ##Upregulated gene表达量热图
    pheatmap(mat = gene[which(rownames(gene) %in% sigG_upup),],
             main="the heatmap of Upregulated Genes",
             scale = "row", cellwidth = 60, 
             show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = T,
             labels_col = colnames(gene), angle_col = "45",
             breaks = seq(-2,2,by=0.01), border = F,
             color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                        colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
    ##Downregulated gene表达量热图
    pheatmap(mat = gene[which(rownames(gene) %in% sigG_down),],
             main="the heatmap of Downregulated Genes",
             scale = "row", cellwidth = 60, 
             show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = T,
             labels_col = colnames(gene), angle_col = "45",
             breaks = seq(-2,2,by=0.01), border = F,
             color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                        colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
    ##tete表达量热图
    pheatmap(mat = tete,
             main="the heatmap of all TEs",
             scale = "row", cellwidth = 60, 
             show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = T,
             labels_col = colnames(tete), angle_col = "45",
             breaks = seq(-2,2,by=0.01), border = F,
             color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                        colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
    ##siginificant tete表达量热图
    pheatmap(mat = tete[which(rownames(tete) %in% sigT),],
             main="the heatmap of significant TEs",
             scale = "row", cellwidth = 60, 
             show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = T,
             labels_col = colnames(tete), angle_col = "45",
             breaks = seq(-2,2,by=0.01), border = F,
             color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                        colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
    ##Upregulated tete表达量热图
    pheatmap(mat = tete[which(rownames(tete) %in% sigT_upup),],
             main="the heatmap of Upregulated TEs",
             scale = "row", cellwidth = 60, 
             show_rownames = T, show_colnames = T, cluster_row = T, cluster_cols = T,
             labels_col = colnames(tete), angle_col = "45",
             breaks = seq(-2,2,by=0.01), border = F,
             color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                        colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
    ##Downregulated tete表达量热图
    pheatmap(mat = tete[which(rownames(tete) %in% sigT_down),],
             main="the heatmap of Downregulated TEs",
             scale = "row", cellwidth = 60,
             show_rownames = T, show_colnames = T, cluster_row = T, cluster_cols = T,
             labels_col = colnames(tete), angle_col = "45",
             breaks = seq(-2,2,by=0.01), border = F,
             color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                        colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
    
    dev.off()
  }
  sets <- paste0("group",seq(5))
  for (grou in sets) {
    print(grou)
    shy_pheatmap(mark = "huhuhu", grou,
                 read_cntFile = paste0(path, "count/huhuhu.rename.cntTable"),
                 gene_resFile = paste0(path, "DESeq2/huhuhu.",grou,".geneOnly.csv"),
                 tete_resFile = paste0(path, "DESeq2/huhuhu.",grou,".teOnly.csv"),
                 savePath = paste0(path, "PLOT/pheatmap/"))
  }
}
################################################################################
#10、pheatmap of the sig genes across all groups
{
  shy_pheatmap <- function(mark,sets,read_cntFile,savePath){
    print("reads counted by TEtranscripts")
    #
    data <- read.table(read_cntFile,header = T, row.names = 1, check.names = F, stringsAsFactors = F)
    gene <- data[grep(rownames(data), pattern = ":", invert = T),]
    gene <- gene[rowSums(gene) > 5,]
    gene <- as.data.frame(apply(gene, 2, function(x){x/sum(x)*1000000}))
    tete <- data[grep(rownames(data), pattern = ":", invert = F),]
    tete <- tete[rowSums(tete) > 5,]
    tete <- as.data.frame(apply(tete, 2, function(x){x/sum(x)*1000000}))
    #
    resG <- data.frame()
    for (i in sets) {
      tmp1 <- as.data.frame(read.csv(paste0(path,"DESeq2/huhuhu.",i,".geneOnly.csv"), 
                                     header = T, stringsAsFactors = F))
      resG <- rbind(resG,tmp1)
    }
    resG <- na.omit(resG)
    sigG <- resG[which(abs(resG$log2FoldChange)>log2(1.5) & resG$padj<0.05), "X"]
    sigG_upup <- resG[which(resG$log2FoldChange>log2(1.5) & resG$padj<0.05), "X"]
    sigG_down <- resG[which(resG$log2FoldChange < -log2(1.5) & resG$padj<0.05), "X"]
    sigG <- unique(sigG)
    sigG_upup <- unique(sigG_upup)
    sigG_down <- unique(sigG_down)
    #
    resG <- data.frame()
    for (i in sets) {
      tmp1 <- as.data.frame(read.csv(paste0(path,"DESeq2/huhuhu.",i,".teOnly.csv"), 
                                     header = T, stringsAsFactors = F))
      resG <- rbind(resG,tmp1)
    }
    resG <- na.omit(resG)
    sigT <- resG[which(abs(resG$log2FoldChange)>log2(1.5) & resG$padj<0.05), "X"]
    sigT_upup <- resG[which(resG$log2FoldChange>log2(1.5) & resG$padj<0.05), "X"]
    sigT_down <- resG[which(resG$log2FoldChange < -log2(1.5) & resG$padj<0.05), "X"]
    sigT <- unique(sigT)
    sigT_upup <- unique(sigT_upup)
    sigT_down <- unique(sigT_down)
    ##开始画图
    pdf(paste0(savePath,mark,".groups.pheatmap.pdf"), width = 10, height = 10)
    ##gene表达量热图
    pheatmap(mat = gene,
             main="the heatmap of all Genes",
             scale = "row", cellwidth = 20, 
             show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = T,
             labels_col = colnames(gene), angle_col = "45",
             breaks = seq(-2,2,by=0.01), border = F,
             color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                        colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
    ##siginificant gene表达量热图
    pheatmap(mat = gene[which(rownames(gene) %in% sigG),],
             main="the heatmap of significant Genes",
             scale = "row", cellwidth = 20, 
             show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = T,
             labels_col = colnames(gene), angle_col = "45",
             breaks = seq(-2,2,by=0.01), border = F,
             color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                        colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
    ##Upregulated gene表达量热图
    pheatmap(mat = gene[which(rownames(gene) %in% sigG_upup),],
             main="the heatmap of Upregulated Genes",
             scale = "row", cellwidth = 20, 
             show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = T,
             labels_col = colnames(gene), angle_col = "45",
             breaks = seq(-2,2,by=0.01), border = F,
             color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                        colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
    ##Downregulated gene表达量热图
    pheatmap(mat = gene[which(rownames(gene) %in% sigG_down),],
             main="the heatmap of Downregulated Genes",
             scale = "row", cellwidth = 20, 
             show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = T,
             labels_col = colnames(gene), angle_col = "45",
             breaks = seq(-2,2,by=0.01), border = F,
             color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                        colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
    ##tete表达量热图
    pheatmap(mat = tete,
             main="the heatmap of all TEs",
             scale = "row", cellwidth = 20, 
             show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = T,
             labels_col = colnames(tete), angle_col = "45",
             breaks = seq(-2,2,by=0.01), border = F,
             color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                        colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
    ##siginificant tete表达量热图
    pheatmap(mat = tete[which(rownames(tete) %in% sigT),],
             main="the heatmap of significant TEs",
             scale = "row", cellwidth = 20, 
             show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = T,
             labels_col = colnames(tete), angle_col = "45",
             breaks = seq(-2,2,by=0.01), border = F,
             color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                        colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
    ##Upregulated tete表达量热图
    pheatmap(mat = tete[which(rownames(tete) %in% sigT_upup),],
             main="the heatmap of Upregulated TEs",
             scale = "row", cellwidth = 20, 
             show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = T,
             labels_col = colnames(tete), angle_col = "45",
             breaks = seq(-2,2,by=0.01), border = F,
             color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                        colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
    ##Downregulated tete表达量热图
    pheatmap(mat = tete[which(rownames(tete) %in% sigT_down),],
             main="the heatmap of Downregulated TEs",
             scale = "row", cellwidth = 20,
             show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = T,
             labels_col = colnames(tete), angle_col = "45",
             breaks = seq(-2,2,by=0.01), border = F,
             color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                        colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
    
    dev.off()
  }
  shy_pheatmap(mark = "huhuhu", sets = paste0("group",seq(5)),
               read_cntFile = paste0(path, "count/huhuhu.rename.cntTable"),
               savePath = paste0(path, "PLOT/pheatmap/"))
}
################################################################################
#11、按CPM画下基因表达的Boxplot看看高糖和低糖基因表达分布是否一样的
{
  tmp1 <- openxlsx::read.xlsx(paste0(path, "count/huhuhu.cpm.onlyGene.rename.xlsx"))
  tmp2 <- tmp1 %>% tidyr::pivot_longer(cols = -theName, names_to = "group", values_to = "cpm")
  tmp2$mark <- gsub(x = tmp2$group,pattern = "group.*-P[0-9]-(....).*",replacement = "\\1")
  ggplot(data = tmp2) +
    geom_boxplot(aes(x = group, y = log2(cpm), fill = mark)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.margin = margin(1,1,1,3,unit = "cm"))
  
  tmp1[1,]
  tmp3 <- as.data.frame(summary(tmp1[,-1]))
  tmp3 <- tmp3[,-1]
  tmp3 <- tmp3 %>% tidyr::separate(col = Freq, into = c("clas","vaue"), sep = ":")
  tmp4 <- tmp3 %>% tidyr::pivot_wider(names_from = clas, values_from = vaue)
  tmp4 <- as.data.frame(tmp4)
  View(tmp4)
}
################################################################################
#12、测试scale方法 [以group2为例]
{
  ff01 <- paste0(path,"testNormalization/huhuhu.group2.cntTable")
  mark = "huhuhu"
  clas = factor(c("OP","OP","WT","WT")); vsvs = c("clas","OP","WT")
  #
  #
  #
  #方案一, Quantile
  {
    #来自https://www.spsanderson.com/steveondata/posts/2024-03-28/index.html
    qnqn <- function(myta){
      data_sort <- apply(myta, 2, sort)
      rows_mean <- rowMeans(data_sort)
      data_sort <- matrix(rows_mean, 
                          nrow = nrow(data_sort), 
                          ncol = ncol(data_sort), 
                          byrow = F)
      inex_rank <- apply(myta, 2, order)
      norm_data <- matrix(nrow = nrow(myta), ncol = ncol(myta))
      for(i in 1:ncol(myta)){
        tmps <- data.frame(aaaa = data_sort[,i], bbbb = inex_rank[,i])
        tmps <- tmps[order(tmps$bbbb, decreasing = F),]
        norm_data[,i] <- tmps$aaaa
      }
      return(norm_data)
    }
    tempss <- read.table(ff01, header = T, check.names = F, stringsAsFactors = F, row.names = 1)
    tempss <- tempss[grep(rownames(tempss), pattern = ":", invert = F),]
    normed <- qnqn(myta = tempss); normed <- as.data.frame(normed)
    rownames(normed) <- rownames(tempss)
    colnames(normed) <- colnames(tempss)
    normed[grep("MERVL-int",rownames(normed)),]
    tempss[grep("MERVL-int",rownames(tempss)),]
    #
    #
    #
    normed <- normed[order(normed$`group2-E14-P0-0180+Pyru-rep2`, decreasing = T),]
    normed$rank <- seq_along(normed$`group2-E14-P0-0180+Pyru-rep1`)
    normed[grep("MERVL-int",rownames(normed)),]
    tempss <- read.table(ff01, header = T, check.names = F, stringsAsFactors = F, row.names = 1)
    tempss <- tempss[order(tempss$`group2-E14-P0-4500+Pyru-rep2`, decreasing = T),]
    tempss$rank <- seq_along(tempss$`group2-E14-P0-4500+Pyru-rep1`)
    tempss[grep("MERVL-int",rownames(tempss)),]
    #
    #
    #
    normed <- ceiling(normed); tempss <- normed
    dfclas <- data.frame(row.names = colnames(tempss), clas, stringsAsFactors = T)
    colLen <- ceiling(length(colnames(tempss))/2)
    geOnly <- tempss[grep(rownames(tempss), pattern = ":", invert = T),]
    teOnly <- tempss[grep(rownames(tempss), pattern = ":", invert = F),]
    teGene <- tempss
    #
    teGene <- DESeq2::DESeqDataSetFromMatrix(countData = teGene, colData = dfclas, design = ~clas)
    teGene <- teGene[rowSums(BiocGenerics::counts(teGene) > 2) >= colLen, ]#[rowMeans(BiocGenerics::counts(teGene)) > 1,]
    teGene <- DESeq2::DESeq(teGene)
    teGene <- DESeq2::results(teGene, contrast = vsvs)
    teGene <- teGene[grep(rownames(teGene), pattern = ":", invert = F),]
    write.csv(teGene, paste0(path,"testNormalization/",mark,".te.GeneBackGround.csv"))
    #
    geOnly <- DESeq2::DESeqDataSetFromMatrix(countData = geOnly, colData = dfclas, design = ~clas)
    geOnly <- geOnly[rowSums(BiocGenerics::counts(geOnly) > 2) >= colLen, ]#[rowMeans(BiocGenerics::counts(geOnly)) > 1,]
    geOnly <- DESeq2::DESeq(geOnly,)
    geOnly <- DESeq2::results(geOnly, contrast = vsvs)
    write.csv(geOnly, paste0(path,"testNormalization/",mark,".geneOnly.csv"))
    #
    teOnly <- DESeq2::DESeqDataSetFromMatrix(countData = teOnly, colData = dfclas, design = ~clas)
    teOnly <- teOnly[rowSums(BiocGenerics::counts(teOnly) > 2) >= colLen, ]#[rowMeans(BiocGenerics::counts(teOnly)) > 1,]
    teOnly <- DESeq2::DESeq(teOnly,)
    teOnly <- DESeq2::results(teOnly, contrast = vsvs)
    write.csv(teOnly, paste0(path,"testNormalization/",mark,".teOnly.csv"))
  }
  #方案二, Housekeeping gene, 8407个mouse的
  {
    tempss <- read.table(ff01, header = T, check.names = F, stringsAsFactors = F, row.names = 1)
    hkGene <- openxlsx::read.xlsx("/Reference/aaSHY/DatabaseFiles/houseKeep/set1/hk_mm10_15tissues.xlsx",sheet = "Sheet1")
    #
    normed <-tempss
    dfclas <- data.frame(row.names = colnames(tempss), clas, stringsAsFactors = T)
    colLen <- ceiling(length(colnames(tempss))/2)
    geOnly <- tempss[grep(rownames(tempss), pattern = ":", invert = T),]
    teOnly <- tempss[grep(rownames(tempss), pattern = ":", invert = F),]
    teGene <- tempss
    teGene <- DESeq2::DESeqDataSetFromMatrix(countData = teGene, colData = dfclas, design = ~clas)
    size_1 <- DESeq2::DESeq(teGene)
    size_1@colData
    #
    tttt <- teGene[which(row.names(teGene) %in% unique(hkGene$geneName)),]
    size_2 <- DESeq2::DESeq(tttt)
    size_2@colData
    
    
    
    
    
    
    
    
    teGene <- teGene[rowSums(BiocGenerics::counts(teGene) > 2) >= colLen, ]#[rowMeans(BiocGenerics::counts(teGene)) > 1,]
    teGene <- DESeq2::DESeq(teGene)
    teGene <- DESeq2::results(teGene, contrast = vsvs)
    teGene <- teGene[grep(rownames(teGene), pattern = ":", invert = F),]
    write.csv(teGene, paste0(path,"testNormalization/",mark,".te.GeneBackGround.csv"))
    
  }
  #方案三, Housekeeping gene, 3804个human的, 转换获得3460个小鼠同源基因
  {
    tempss <- read.table(ff01, header = T, check.names = F, stringsAsFactors = F, row.names = 1)
    dfclas <- data.frame(row.names = colnames(tempss), clas, stringsAsFactors = T)
    colLen <- ceiling(length(colnames(tempss))/2)
    teGene <- tempss
    teGene <- DESeq2::DESeqDataSetFromMatrix(countData = teGene, colData = dfclas, design = ~clas)
    teGene <- teGene[rowSums(BiocGenerics::counts(teGene) > 2) >= colLen, ]#[rowMeans(BiocGenerics::counts(teGene)) > 1,]
    size_1 <- DESeq2::DESeq(teGene); size_1@colData$sizeFactor
    mervls <- DESeq2::results(size_1, contrast = vsvs)
    mervls <- as.data.frame(mervls)
    mervls[grep("MERVL-int",rownames(mervls)),]
    #
    hkGene <- read.csv("/Reference/aaSHY/DatabaseFiles/houseKeep/set2/hk_human_3804.txt", header = F)
    library("babelgene")
    hkGene <- babelgene::orthologs(genes = unique(hkGene$V1), species = "mouse", human = T)
    hkGene <- unique(hkGene$symbol)
    tttt <- teGene[which(row.names(teGene) %in% hkGene),]
    #
    size_1 <- estimateSizeFactors(teGene)
    size_2 <- estimateSizeFactors(tttt)
    size_1$sizeFactor <- size_2$sizeFactor#更改size factor后, 继续跑
    size_1 <- estimateDispersions(size_1)
    size_1 <- nbinomWaldTest(size_1)
    #
    mervls <- DESeq2::results(size_1, contrast = vsvs)
    mervls <- as.data.frame(mervls)
    mervls[grep("MERVL-int",rownames(mervls)),]
  }
  #方案四, Housekeeping gene, 3804个human的, 转换获得3460个小鼠同源基因，并筛选高表达基因
  {
    tempss <- read.table(ff01, header = T, check.names = F, stringsAsFactors = F, row.names = 1)
    dfclas <- data.frame(row.names = colnames(tempss), clas, stringsAsFactors = T)
    colLen <- ceiling(length(colnames(tempss))/2)
    teGene <- tempss
    teGene <- DESeq2::DESeqDataSetFromMatrix(countData = teGene, colData = dfclas, design = ~clas)
    #teGene <- teGene[rowSums(BiocGenerics::counts(teGene) > 2) >= colLen, ]
    #teGene <- teGene[rowMeans(BiocGenerics::counts(teGene)) > 1,]
    size_1 <- DESeq2::DESeq(teGene); size_1@colData$sizeFactor
    mervls <- DESeq2::results(size_1, contrast = vsvs)
    mervls <- as.data.frame(mervls)
    mervls[grep("MERVL-int",rownames(mervls)),]
    #
    hkGene <- read.csv("/Reference/aaSHY/DatabaseFiles/houseKeep/set2/hk_human_3804.txt", header = F)
    library("babelgene")
    hkGene <- babelgene::orthologs(genes = unique(hkGene$V1), species = "mouse", human = T)
    hkGene <- unique(hkGene$symbol)
    tttt <- teGene[which(row.names(teGene) %in% hkGene),]
    #tttt <- tttt[rowSums(BiocGenerics::counts(tttt) > 10000) == length(colnames(tttt)), ]
    tttt <- tttt[rowMeans(BiocGenerics::counts(tttt)) > 50000,]
    row.names(tttt)
    #
    size_1 <- estimateSizeFactors(teGene)
    size_2 <- estimateSizeFactors(tttt)
    size_1$sizeFactor <- size_2$sizeFactor#更改size factor后, 继续跑
    size_1 <- estimateDispersions(size_1)
    size_1 <- nbinomWaldTest(size_1)
    #
    mervls <- DESeq2::results(size_1, contrast = vsvs)
    mervls <- as.data.frame(mervls)
    mervls[grep("MERVL-int",rownames(mervls)),]
  }
  #方案五, Housekeeping gene, 8407个mouse的，并筛选高表达基因
  {
    tempss <- read.table(ff01, header = T, check.names = F, stringsAsFactors = F, row.names = 1)
    dfclas <- data.frame(row.names = colnames(tempss), clas, stringsAsFactors = T)
    colLen <- ceiling(length(colnames(tempss))/2)
    teGene <- tempss
    teGene <- DESeq2::DESeqDataSetFromMatrix(countData = teGene, colData = dfclas, design = ~clas)
    #teGene <- teGene[rowSums(BiocGenerics::counts(teGene) > 2) >= colLen, ]
    #teGene <- teGene[rowMeans(BiocGenerics::counts(teGene)) > 1,]
    size_1 <- DESeq2::DESeq(teGene); size_1@colData$sizeFactor
    mervls <- DESeq2::results(size_1, contrast = vsvs)
    mervls <- as.data.frame(mervls)
    mervls[grep("MERVL-int",rownames(mervls)),]
    #
    hkGene <- openxlsx::read.xlsx("/Reference/aaSHY/DatabaseFiles/houseKeep/set1/hk_mm10_15tissues.xlsx",sheet = "Sheet1")
    tttt <- teGene[which(row.names(teGene) %in% hkGene$geneName),]
    tttt <- tttt[rowSums(BiocGenerics::counts(tttt) > 50000) == length(colnames(tttt)), ]
    #tttt <- tttt[rowMeans(BiocGenerics::counts(tttt)) > 50000,]
    row.names(tttt)
    #
    size_1 <- estimateSizeFactors(teGene)
    size_2 <- estimateSizeFactors(tttt)
    size_1$sizeFactor <- size_2$sizeFactor#更改size factor后, 继续跑
    size_1 <- estimateDispersions(size_1)
    size_1 <- nbinomWaldTest(size_1)
    #
    mervls <- DESeq2::results(size_1, contrast = vsvs)
    mervls <- as.data.frame(mervls)
    mervls[grep("MERVL-int",rownames(mervls)),]
  }
  #方案六, 按低糖组reads count取top 200个基因，做富集，找GO重要term
  {
    tempss <- read.table(ff01, header = T, check.names = F, stringsAsFactors = F, row.names = 1)
    geOnly <- tempss[grep(rownames(tempss), pattern = ":", invert = T),]
    xuan_1 <- geOnly[order(rowMeans(geOnly[,c(1,2)]), decreasing = T)[1:200],]
    #
    xuan_1 <- rownames(xuan_1)
    xuan_1 <- na.omit(mapIds(x=org.Mm.eg.db, keys=unique(xuan_1), keytype="SYMBOL", column="ENTREZID"))
    gogo_1 <- clusterProfiler::enrichGO(gene = xuan_1, OrgDb = org.Mm.eg.db, 
                                        keyType = "ENTREZID", ont = "MF",
                                        pvalueCutoff = 0.05,qvalueCutoff = 0.05, readable = T)
    gogo_1 <- gogo_1@result
    gogo_1 <- gogo_1[which(gogo_1$p.adjust <0.05),]
    #View(gogo_1), 得到：
    #GO:0002181, cytoplasmic translation
    #GO:0051169, nuclear transport
    #GO:0022613, ribonucleoprotein complex biogenesis
    #GO:0140694, non-membrane-bounded organelle assembly
    base::load(file = "/Reference/aaSHY/DatabaseFiles/GO_DATA.RData")
    gogo_s <- GO_DATA$PATHID2EXTID
    xuan_s <- c(gogo_s[["GO:0002181"]],gogo_s[["GO:0051169"]],gogo_s[["GO:0022613"]],gogo_s[["GO:0140694"]])
    xuan_s <- unique(xuan_s)
    #
    dfclas <- data.frame(row.names = colnames(tempss), clas, stringsAsFactors = T)
    colLen <- ceiling(length(colnames(tempss))/2)
    teGene <- tempss
    teGene <- DESeq2::DESeqDataSetFromMatrix(countData = teGene, colData = dfclas, design = ~clas)
    #teGene <- teGene[rowSums(BiocGenerics::counts(teGene) > 2) >= colLen, ]
    #teGene <- teGene[rowMeans(BiocGenerics::counts(teGene)) > 1,]
    size_1 <- DESeq2::DESeq(teGene); size_1@colData$sizeFactor
    mervls <- DESeq2::results(size_1, contrast = vsvs)
    mervls <- as.data.frame(mervls)
    mervls[grep("MERVL-int",rownames(mervls)),]
    #
    hkGene <- xuan_s
    tttt <- teGene[which(row.names(teGene) %in% hkGene),]
    tttt <- tttt[rowSums(BiocGenerics::counts(tttt) > 20000) == length(colnames(tttt)), ]
    #tttt <- tttt[rowMeans(BiocGenerics::counts(tttt)) > 50000,]
    row.names(tttt)
    #
    size_1 <- estimateSizeFactors(teGene)
    size_2 <- estimateSizeFactors(tttt)
    size_1$sizeFactor <- size_2$sizeFactor#更改size factor后, 继续跑
    size_1 <- estimateDispersions(size_1)
    size_1 <- nbinomWaldTest(size_1)
    #
    mervls <- DESeq2::results(size_1, contrast = vsvs)
    mervls <- as.data.frame(mervls)
    mervls[grep("MERVL-int",rownames(mervls)),]
  }
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
################################################################################
#13、MA plot, group2, 有、无基因背景
{
  shy_teMA <- function(difffile, savefile){
    rf01 <- read.csv(file = difffile, row.names = 1)
    rf01$subFam <- sapply(str_split(rownames(rf01),pattern = ":"),"[",1)
    rf01$family <- sapply(str_split(rownames(rf01),pattern = ":"),"[",2)
    rf01$classs <- sapply(str_split(rownames(rf01),pattern = ":"),"[",3)
    rf01$colo <- "other"
    for (i in c("ERV1","ERVK","ERVL","Gypsy")){
      rf01$colo[which(rf01$padj<0.05 & abs(rf01$log2FoldChange)>log2(1.5) & rf01$family==i)] <- i
    }
    rf02 <- rf01
    rf02$colo <- factor(rf02$colo,levels=c("ERV1","ERVK","ERVL","Gypsy","other"))
    myCo <- list(ERV1="sienna1",ERVK="steelblue3",ERVL="palevioletred2",Gypsy="green",Others="grey80")
    myCo <- unlist(myCo)
    #
    p_01 <- ggplot() +
      geom_point(data = rf02, 
                 aes(x=log2(baseMean+1), y=log2FoldChange, color=colo, shape=NULL),
                 size=1.5, alpha=1) +
      guides(color = guide_legend(title = "repFamliy")) +
      labs(x = "log2(Expression in baseMean)", y = "log2Foldchange") +
      geom_hline(aes(yintercept =  log2(1.5)), linetype="dashed", colour="grey") +
      geom_hline(aes(yintercept = -log2(1.5)), linetype="dashed", colour="grey") +
      geom_text_repel(data = rf02[which(rf02$colo!="other"),], 
                      aes(x=log2(baseMean+1), y=log2FoldChange, label=subFam), 
                      hjust=-0.125, size=4, colour="black", family="serif") +
      scale_color_manual(values = myCo)+
      scale_y_continuous(limits = c(floor(min(rf02$log2FoldChange)),ceiling(max(rf02$log2FoldChange))),
                         breaks = seq(floor(min(rf02$log2FoldChange)),ceiling(max(rf02$log2FoldChange)))) +
      theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.line = element_line(),
            axis.text = element_text(family = "serif"),
            panel.border = element_rect(colour = "black", linewidth = 0.8))+
      theme(legend.position.inside = c(0.91,0.25),legend.background = element_blank(),
            legend.text = element_text(family = "serif"))
    ggsave(plot = p_01,
           units = "cm", height = 12, width = 13.5577, filename = savefile)
  }
  shy_teMA(difffile = paste0(path,"DESeq2/huhuhu.group2.te.GeneBackGround.csv"),
           savefile = paste0(path,"PLOT/huhuhu.group2.te.GeneBackGround.MA.pdf"))
  shy_teMA(difffile = paste0(path,"DESeq2/huhuhu.group2.teOnly.csv"),
           savefile = paste0(path,"PLOT/huhuhu.group2.teOnly.MA.pdf"))
}
################################################################################
#14、2C gene GSEA, group2
{
  temp <- read.csv(paste0(path, "DESeq2/huhuhu.group2.geneOnly.csv"), header = T)
  temp <- temp[, c("X", "log2FoldChange")]
  temp <- temp[order(temp$log2FoldChange, decreasing = T),]
  write.table(x = temp, file = paste0(path, "DESeq2/huhuhu.group2.res.rnk"), quote = F, sep = "\t", col.names = F, row.names = F)
  cmd_12 <- paste0("gseapy prerank -g ","/Reference/aaSHY/GSEAgmt/TwoGenes.gmt -r ",
                   paste0(path, "DESeq2/huhuhu.group2.res.rnk"), " -p 10 -o ", 
                   paste0(path, "PLOT/enrichGSEA/"))#conda activate base gseapy 0.10.8
  print(cmd_12)
}
################################################################################
#15、TE copies CPM Excel, group2 【暂时不做了,因为要做unique map模式】
{
  cnt1 <- read.table("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/thehuhuhu/count/huhuhu.site.rename.cntTable", header = T, sep = "\t", row.names = 1)
  #tmp1 <- cnt1[grep(x = rownames(cnt1), pattern=":", invert = F),]#TE
  #tmp2 <- cnt1[grep(x = rownames(cnt1), pattern=":", invert = T),]#Gene
  tmp3 <- cnt1#All
  tmp1 <- as.data.frame(apply(tmp1, 2, function(x){x/sum(x)*1000000}))
  tmp2 <- as.data.frame(apply(tmp2, 2, function(x){x/sum(x)*1000000}))
  tmp3 <- as.data.frame(apply(tmp3, 2, function(x){x/sum(x)*1000000}))
  tmp1$theName <- rownames(tmp1)
  tmp2$theName <- rownames(tmp2)
  tmp3$theName <- rownames(tmp3)
  openxlsx::write.xlsx(x = tmp1, file = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202404/gata3OE/count/gata3OE.cpm.onlyTE.xlsx")
  openxlsx::write.xlsx(x = tmp2, file = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202404/gata3OE/count/gata3OE.cpm.onlyGene.xlsx")
  openxlsx::write.xlsx(x = tmp3, file = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202404/gata3OE/count/gata3OE.cpm.allGeneTE.xlsx")
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
################################################################################
#16、Unique map模式下跑count reads
{
  #跑TEtranscripts、TElocal的uniq模式、以及featureCounts的ignoreDup的模式
  for (i in c("stringtie/countUniq","DESeq2Uniq","countUniq","PLOT_uniq")) {
    if (dir.exists(paths = paste0(path, i))==F) {
      dir.create(path = paste0(path,i), recursive = T)
    }
  }
  shell <- paste0(path,"src/run_e.sh")
  cat("#!/bin/bash\n", file = shell)
  ff01 <- list.files(path = paste0(path,"fq"),pattern = "_1.fq.gz")
  prex <- sapply(str_split(ff01,pattern = "_1.fq.gz"),"[",1)
  cmd_01 <- paste0("TEtranscripts -t ",
                   paste(paste0(path,"bam/",prex[1:10], ".bam"),collapse=" "),
                   " -c ",
                   paste(paste0(path,"bam/",prex[11:20],".bam"),collapse=" "),
                   " --GTF ",inx3," --TE ",inx4," --sortByPos --mode uniq",
                   " --project huhuhu --outdir ",paste0(path,"countUniq"),"\n")
  cmd_02 <- paste0("TElocal -b ",
                   paste0(path,"bam/",prex,".bam"),
                   " --GTF ",inx3," --TE ",inx5," --sortByPos --mode uniq",
                   " --project ",paste0(path, "countUniq/",prex),"\n")
  cmd_03 <- paste0("TEtranscripts -t ",
                   paste(paste0(path,"bam/",prex[1:10], ".bam"),collapse = " "),
                   " -c ",
                   paste(paste0(path,"bam/",prex[11:20],".bam"),collapse = " "),
                   " --GTF ",path,"stringtie/huhuhu.merge.R.gtf",
                   " --TE ", inx4," --sortByPos --mode uniq --project huhuhu",
                   " --outdir ",path,"stringtie/countUniq","\n")
  cmd_04 <- paste0("featureCounts -T 16 -g 'gene_name' ",
                   "-F 'GTF' --ignoreDup -s 0 -p -P -d 0 -D 600 -B ",
                   "-a ",inx3," -o ", path,"countUniq/featureCounts_gene.txt ",
                   paste(paste0(path,"bam/",prex, ".bam"),collapse = " "),"\n")
  cmd_05 <- paste0("featureCounts -T 16 -g 'gene_id' ",
                   "-F 'GTF' --ignoreDup -s 0 -p -P -d 0 -D 600 -B ",
                   "-a ",inx4," -o ", path,"countUniq/featureCounts_tete.txt ",
                   paste(paste0(path,"bam/",prex, ".bam"),collapse = " "),"\n")
  cmd_06 <- paste0("featureCounts -T 16 -g 'transcript_id' ",
                   "-F 'GTF' --ignoreDup -s 0 -p -P -d 0 -D 600 -B ",
                   "-a ",inx4," -o ", path,"countUniq/featureCounts_tete_trans.txt ",
                   paste(paste0(path,"bam/",prex, ".bam"),collapse = " "),"\n")
  for (i in cmd_01) {cat(i, file = shell, append = T)}
  for (i in cmd_02) {cat(i, file = shell, append = T)}
  for (i in cmd_03) {cat(i, file = shell, append = T)}
  for (i in cmd_04) {cat(i, file = shell, append = T)}
  for (i in cmd_05) {cat(i, file = shell, append = T)}
  for (i in cmd_06) {cat(i, file = shell, append = T)}
  print(paste0("nohup bash ",shell, " >",paste0(path,"Log/"),basename(shell),".log"," 2>&1 &"))
}
################################################################################
#17、CPM, for featureCounts、TEtools
#####更改列名 (样本名)
{
  tmps <- read.table("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/thehuhuhu/src/rename.txt",header = F, sep = "\t")
  tmp1 <- read.table(file = paste0(path,"countUniq/huhuhu.cntTable"),
                     sep = "\t", header = T, row.names = 1)
  colnames(tmp1) <- gsub("X.ChIP_seq_2.aaSHY.ghyuu.rnaSEQ.thehuhuhu.bam.(.*).bam.*","\\1",colnames(tmp1))
  colnames(tmp1) <- gsub("\\.","-",colnames(tmp1))
  tmp1$genesName <- rownames(tmp1)
  tmp1 <- tmp1[,c("genesName",tmps$V1)]
  colnames(tmp1) <- c("genesName",tmps$V2)
  write.table(x = tmp1, file = paste0(path,"countUniq/huhuhu.rename.cntTable.txt"),
              sep = "\t", quote = F, col.names = T, row.names = F)
}
#####计算featureCounts的CPM
{
  #for featureCounts
  #step1: 更改TE名后与基因表合并, 使其变为与TEtranscripts一样的格式
  {
    tf01 <- import(inx4); tf01 <- as.data.frame(tf01); colnames(tf01)[10] <- "genesName"
    cnt1 <- read.table(paste0(path,"countUniq/featureCounts_gene.txt"), header = T, sep = "\t")
    cnt2 <- read.table(paste0(path,"countUniq/featureCounts_tete.txt"), header = T, sep = "\t")
    cnt3 <- read.table(paste0(path,"countUniq/featureCounts_tete_trans.txt"), header = T, sep = "\t")
    #
    cnt2 <- dplyr::left_join(cnt2, distinct(tf01[,c(10,12,13)]),by="genesName")
    cnt2$genesName <- paste0(cnt2$genesName,":",cnt2$family_id,":",cnt2$class_id)
    cnt2 <- cnt2[,-c(22,23)]
    wr01 <- rbind(cnt1,cnt2)
    #
    colnames(tf01)[c(10,11)] <- c("gene_id","genesName")
    cnt3 <- dplyr::left_join(cnt3, distinct(tf01[,c(10,11,12,13)]),by="genesName")
    cnt3$genesName <- paste0(cnt3$genesName,":",cnt3$gene_id,":",cnt3$family_id,":",cnt3$class_id)
    cnt3 <- cnt3[,-c(22,23,24)]
    wr02 <- rbind(cnt1,cnt3)
    #
    write.table(x = wr01, 
                file = paste0(path,"countUniq/featureCounts_genetete.txt"),
                col.names = T, row.names = F, quote = F, sep = "\t")
    write.table(x = wr02, 
                file = paste0(path,"countUniq/featureCounts_genetetetrans.txt"),
                col.names = T, row.names = F, quote = F, sep = "\t")
    rm(cnt1,cnt2,cnt3,wr01,wr02); tf01 <- import(inx4)
  }
  #step2: 计算CPM
  {
    cnt1 <- read.table(paste0(path,"countUniq/featureCounts_genetete.txt"), header = T, sep = "\t", row.names = 1)
    tmp1 <- cnt1[grep(x = rownames(cnt1), pattern=":", invert = F),]#TE
    tmp2 <- cnt1[grep(x = rownames(cnt1), pattern=":", invert = T),]#Gene
    tmp3 <- cnt1#All
    tmp1 <- as.data.frame(apply(tmp1, 2, function(x){x/sum(x)*1000000}))
    tmp2 <- as.data.frame(apply(tmp2, 2, function(x){x/sum(x)*1000000}))
    tmp3 <- as.data.frame(apply(tmp3, 2, function(x){x/sum(x)*1000000}))
    tmp1$theName <- rownames(tmp1)
    tmp2$theName <- rownames(tmp2)
    tmp3$theName <- rownames(tmp3)
    openxlsx::write.xlsx(x = tmp1, file = paste0(path, "countUniq/fC_huhuhu.cpm.onlyTE.xlsx"))
    openxlsx::write.xlsx(x = tmp2, file = paste0(path, "countUniq/fC_huhuhu.cpm.onlyGene.xlsx"))
    openxlsx::write.xlsx(x = tmp3, file = paste0(path, "countUniq/fC_huhuhu.cpm.allGeneTE.xlsx"))
  }
  #step3: 计算site的CPM
  {
    cnt1 <- read.table(paste0(path,"countUniq/featureCounts_tete_trans.txt"), header = T, sep = "\t", row.names = 1)
    colLen <- length(colnames(cnt1))
    tmp1 <- cnt1[rowSums(cnt1 == 0) < colLen, ]
    tmp1 <- as.data.frame(apply(tmp1, 2, function(x){x/sum(x)*1000000}))
    tmp1$theName <- rownames(tmp1)
    openxlsx::write.xlsx(x = tmp1, file = paste0(path, "countUniq/fC_huhuhu.cpm.onlyTE.site.xlsx"))
  }
}
#####计算TEtools的CPM
{
  #for TEtools
  cnt1 <- read.table(paste0(path,"countUniq/huhuhu.rename.cntTable"), header = T, sep = "\t", row.names = 1)
  tmp1 <- cnt1[grep(x = rownames(cnt1), pattern=":", invert = F),]#TE
  tmp2 <- cnt1[grep(x = rownames(cnt1), pattern=":", invert = T),]#Gene
  tmp3 <- cnt1#All
  tmp1 <- as.data.frame(apply(tmp1, 2, function(x){x/sum(x)*1000000}))
  tmp2 <- as.data.frame(apply(tmp2, 2, function(x){x/sum(x)*1000000}))
  tmp3 <- as.data.frame(apply(tmp3, 2, function(x){x/sum(x)*1000000}))
  tmp1$theName <- rownames(tmp1)
  tmp2$theName <- rownames(tmp2)
  tmp3$theName <- rownames(tmp3)
  openxlsx::write.xlsx(x = tmp1, file = paste0(path, "countUniq/TEtools_huhuhu.cpm.onlyTE.xlsx"))
  openxlsx::write.xlsx(x = tmp2, file = paste0(path, "countUniq/TEtools_huhuhu.cpm.onlyGene.xlsx"))
  openxlsx::write.xlsx(x = tmp3, file = paste0(path, "countUniq/TEtools_huhuhu.cpm.allGeneTE.xlsx"))
}
#####计算TEtools的site的CPM 【去除了在所有样本中reads count均为0的TE site】
{
  cnt1 <- read.table(paste0(path,"countUniq/TEtools/huhuhu.site.rename.cntTable"), header = T, sep = "\t", row.names = 1)
  colLen <- length(colnames(cnt1))
  cnt1 <- cnt1[rowSums(cnt1 == 0) < colLen, ]
  tmp1 <- cnt1[grep(x = rownames(cnt1), pattern=":", invert = F),]#TE
  tmp1 <- as.data.frame(apply(tmp1, 2, function(x){x/sum(x)*1000000}))
  tmp1$theName <- rownames(tmp1)
  openxlsx::write.xlsx(x = tmp1, file = paste0(path, "countUniq/TEtools_huhuhu.cpm.onlyTE.site.xlsx"))
}
################################################################################
#18、比较featureCount与TEtools是否有差异 [用group1的CPM的样本间均值]
{
  tmp1 <- read.table("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/thehuhuhu/countUniq/featureCounts_genetete.txt",header = T)
  tmp1 <- tmp1[grep(":",tmp1$genesName,invert = F),]
  tmp1$meanss <- rowMeans(tmp1[,2:5])
  tmp2 <- read.table("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/thehuhuhu/countUniq/huhuhu.rename.cntTable",header = T)
  tmp2 <- tmp2[grep(":",tmp2$genesName,invert = F),]
  tmp2$meanss <- rowMeans(tmp2[,2:5])
  cous <- dplyr::left_join(tmp1[,c(1,22)],tmp2[,c(1,22)],by="genesName")
  ggplot(cous) +
    geom_point(aes(x = meanss.x, y = meanss.y)) +
    labs(x = "featureCounts", y="TEtools") +
    geom_abline(intercept = 0) +
    theme(plot.margin = margin(unit = "cm",1,1,1,1))
}
################################################################################
#19、DESeq2 [featureCounts、TEtools]
{
  shy_diff <- function(mark,grou,geneFile,siteFile,clas,vsvs){
    tmp1 <- read.table(geneFile, header = T, check.names = F, stringsAsFactors = F, row.names = 1)
    tmp1 <- tmp1[,grep(grou,colnames(tmp1))]
    dfclas <- data.frame(row.names = colnames(tmp1), clas, stringsAsFactors = T)
    colLen <- ceiling(length(colnames(tmp1))/2)
    geneAA <- tmp1[grep(rownames(tmp1), pattern = ":", invert = T),]
    teOnly <- tmp1[grep(rownames(tmp1), pattern = ":", invert = F),]
    teGene <- tmp1
    #
    geneAA <- DESeq2::DESeqDataSetFromMatrix(countData = geneAA, colData = dfclas, design = ~clas)
    geneAA <- geneAA[rowSums(BiocGenerics::counts(geneAA) > 2) >= colLen, ]#[rowMeans(BiocGenerics::counts(geneAA)) > 1,]
    geneAA <- DESeq2::DESeq(geneAA,)
    geneAA <- DESeq2::results(geneAA, contrast = vsvs)
    write.csv(geneAA, paste0(path,"DESeq2Uniq/",mark,".",grou,".geneOnly.csv"))
    upup <- subset(geneAA, log2FoldChange > log2(1.5) & padj < 0.05)
    down <- subset(geneAA, log2FoldChange < -log2(1.5)& padj < 0.05)
    write.csv(upup,paste0(path,"DESeq2Uniq/",mark,".",grou,".geneOnly.upup.csv"))
    write.csv(down,paste0(path,"DESeq2Uniq/",mark,".",grou,".geneOnly.down.csv"))
    #
    teOnly <- DESeq2::DESeqDataSetFromMatrix(countData = teOnly, colData = dfclas, design = ~clas)
    teOnly <- teOnly[rowSums(BiocGenerics::counts(teOnly) > 2) >= colLen, ]#[rowMeans(BiocGenerics::counts(teOnly)) > 1,]
    teOnly <- DESeq2::DESeq(teOnly,)
    teOnly <- DESeq2::results(teOnly, contrast = vsvs)
    write.csv(teOnly, paste0(path,"DESeq2Uniq/",mark,".",grou,".teOnly.csv"))
    upup <- subset(teOnly, log2FoldChange > log2(1.5) & padj < 0.05)
    down <- subset(teOnly, log2FoldChange < -log2(1.5)& padj < 0.05)
    write.csv(upup,paste0(path,"DESeq2Uniq/",mark,".",grou,".teOnly.upup.csv"))
    write.csv(down,paste0(path,"DESeq2Uniq/",mark,".",grou,".teOnly.down.csv"))
    #
    teGene <- DESeq2::DESeqDataSetFromMatrix(countData = teGene, colData = dfclas, design = ~clas)
    teGene <- teGene[rowSums(BiocGenerics::counts(teGene) > 2) >= colLen, ]#[rowMeans(BiocGenerics::counts(teGene)) > 1,]
    teGene <- DESeq2::DESeq(teGene,)
    teGene <- DESeq2::results(teGene, contrast = vsvs)
    teGene <- teGene[grep(rownames(teGene), pattern = ":", invert = F),]
    write.csv(teGene, paste0(path,"DESeq2Uniq/",mark,".",grou,".te.GeneBackGround.csv"))
    upup <- subset(teGene, log2FoldChange > log2(1.5) & padj < 0.05)
    down <- subset(teGene, log2FoldChange < -log2(1.5)& padj < 0.05)
    write.csv(upup,paste0(path,"DESeq2Uniq/",mark,".",grou,".te.GeneBackGround.upup.csv"))
    write.csv(down,paste0(path,"DESeq2Uniq/",mark,".",grou,".te.GeneBackGround.down.csv"))
    #
    #
    tmp2 <- read.table(siteFile, header = T, check.names = F, stringsAsFactors = F, row.names = 1)
    tmp2 <- tmp2[,grep(grou,colnames(tmp2))]
    colLen <- ceiling(length(colnames(tmp2))/2)
    #
    teGene <- tmp2
    teGene <- DESeq2::DESeqDataSetFromMatrix(countData = teGene, colData = dfclas, design = ~clas)
    teGene <- teGene[rowSums(BiocGenerics::counts(teGene) > 2) >= colLen, ]#[rowMeans(BiocGenerics::counts(teGene)) > 1,]
    teGene <- DESeq2::DESeq(teGene,)
    teGene <- DESeq2::results(teGene, contrast = vsvs)
    teGene <- teGene[grep(rownames(teGene), pattern = ":", invert = F),]
    write.csv(teGene, paste0(path,"DESeq2Uniq/",mark,".",grou,".te.site.GeneBackGround.csv"))
    upup <- subset(teGene, log2FoldChange > log2(1.5) & padj < 0.05)
    down <- subset(teGene, log2FoldChange < -log2(1.5)& padj < 0.05)
    write.csv(upup, paste0(path,"DESeq2Uniq/",mark,".",grou,".te.site.GeneBackGround.upup.csv"))
    write.csv(down, paste0(path,"DESeq2Uniq/",mark,".",grou,".te.site.GeneBackGround.down.csv"))
    #
    teGene <- tmp2[grep(rownames(tmp2), pattern = ":", invert = F),]
    teGene <- DESeq2::DESeqDataSetFromMatrix(countData = teGene, colData = dfclas, design = ~clas)
    teGene <- teGene[rowSums(BiocGenerics::counts(teGene) > 2) >= colLen, ]#[rowMeans(BiocGenerics::counts(teGene)) > 1,]
    teGene <- DESeq2::DESeq(teGene,)
    teGene <- DESeq2::results(teGene, contrast = vsvs)
    write.csv(teGene, paste0(path,"DESeq2Uniq/",mark,".",grou,".te.site.teOnly.csv"))
    upup <- subset(teGene, log2FoldChange > log2(1.5) & padj < 0.05)
    down <- subset(teGene, log2FoldChange < -log2(1.5)& padj < 0.05)
    write.csv(upup, paste0(path,"DESeq2Uniq/",mark,".",grou,".te.site.teOnly.upup.csv"))
    write.csv(down, paste0(path,"DESeq2Uniq/",mark,".",grou,".te.site.teOnly.down.csv"))
  }
  sets <- paste0("group",seq(5))
  for (grou in sets) {
    print(grou)
    shy_diff(mark = "fC_huhuhu", grou,
             geneFile = paste0(path, "countUniq/featureCounts_genetete.txt"),
             siteFile = paste0(path, "countUniq/featureCounts_genetetetrans.txt"),
             clas = factor(c("OP","OP","WT","WT")), vsvs = c("clas","OP","WT"))
  }
  for (grou in sets) {
    print(grou)
    shy_diff(mark = "TEtools_huhuhu", grou,
             geneFile = paste0(path, "countUniq/huhuhu.rename.cntTable"),
             siteFile = paste0(path, "countUniq/huhuhu.site.rename.cntTable"),
             clas = factor(c("OP","OP","WT","WT")), vsvs = c("clas","OP","WT"))
  }
}
################################################################################
#20、DESeq2: 添加新的比对 [featureCounts、TEtools]
{
  shy_diff <- function(mark,grou,geneFile,siteFile,clas,vsvs){
    tmp1 <- read.table(geneFile, header = T, check.names = F, stringsAsFactors = F, row.names = 1)
    tmp1 <- cbind(tmp1[,grep(sets[[grou]][1],colnames(tmp1))],
                  tmp1[,grep(sets[[grou]][2],colnames(tmp1))])
    dfclas <- data.frame(row.names = colnames(tmp1), clas, stringsAsFactors = T)
    colLen <- ceiling(length(colnames(tmp1))/2)
    geneAA <- tmp1[grep(rownames(tmp1), pattern = ":", invert = T),]
    teOnly <- tmp1[grep(rownames(tmp1), pattern = ":", invert = F),]
    teGene <- tmp1
    #
    geneAA <- DESeq2::DESeqDataSetFromMatrix(countData = geneAA, colData = dfclas, design = ~clas)
    geneAA <- geneAA[rowSums(BiocGenerics::counts(geneAA) > 2) >= colLen, ]#[rowMeans(BiocGenerics::counts(geneAA)) > 1,]
    geneAA <- DESeq2::DESeq(geneAA,)
    geneAA <- DESeq2::results(geneAA, contrast = vsvs)
    write.csv(geneAA, paste0(path,"DESeq2Uniq/",mark,".",grou,".geneOnly.csv"))
    upup <- subset(geneAA, log2FoldChange > log2(1.5) & padj < 0.05)
    down <- subset(geneAA, log2FoldChange < -log2(1.5)& padj < 0.05)
    write.csv(upup,paste0(path,"DESeq2Uniq/",mark,".",grou,".geneOnly.upup.csv"))
    write.csv(down,paste0(path,"DESeq2Uniq/",mark,".",grou,".geneOnly.down.csv"))
    #
    teOnly <- DESeq2::DESeqDataSetFromMatrix(countData = teOnly, colData = dfclas, design = ~clas)
    teOnly <- teOnly[rowSums(BiocGenerics::counts(teOnly) > 2) >= colLen, ]#[rowMeans(BiocGenerics::counts(teOnly)) > 1,]
    teOnly <- DESeq2::DESeq(teOnly,)
    teOnly <- DESeq2::results(teOnly, contrast = vsvs)
    write.csv(teOnly, paste0(path,"DESeq2Uniq/",mark,".",grou,".teOnly.csv"))
    upup <- subset(teOnly, log2FoldChange > log2(1.5) & padj < 0.05)
    down <- subset(teOnly, log2FoldChange < -log2(1.5)& padj < 0.05)
    write.csv(upup,paste0(path,"DESeq2Uniq/",mark,".",grou,".teOnly.upup.csv"))
    write.csv(down,paste0(path,"DESeq2Uniq/",mark,".",grou,".teOnly.down.csv"))
    #
    teGene <- DESeq2::DESeqDataSetFromMatrix(countData = teGene, colData = dfclas, design = ~clas)
    teGene <- teGene[rowSums(BiocGenerics::counts(teGene) > 2) >= colLen, ]#[rowMeans(BiocGenerics::counts(teGene)) > 1,]
    teGene <- DESeq2::DESeq(teGene,)
    teGene <- DESeq2::results(teGene, contrast = vsvs)
    teGene <- teGene[grep(rownames(teGene), pattern = ":", invert = F),]
    write.csv(teGene, paste0(path,"DESeq2Uniq/",mark,".",grou,".te.GeneBackGround.csv"))
    upup <- subset(teGene, log2FoldChange > log2(1.5) & padj < 0.05)
    down <- subset(teGene, log2FoldChange < -log2(1.5)& padj < 0.05)
    write.csv(upup,paste0(path,"DESeq2Uniq/",mark,".",grou,".te.GeneBackGround.upup.csv"))
    write.csv(down,paste0(path,"DESeq2Uniq/",mark,".",grou,".te.GeneBackGround.down.csv"))
    #
    #
    tmp2 <- read.table(siteFile, header = T, check.names = F, stringsAsFactors = F, row.names = 1)
    tmp2 <- cbind(tmp2[,grep(sets[[grou]][1],colnames(tmp2))],
                  tmp2[,grep(sets[[grou]][2],colnames(tmp2))])
    colLen <- ceiling(length(colnames(tmp2))/2)
    #
    teGene <- tmp2
    teGene <- DESeq2::DESeqDataSetFromMatrix(countData = teGene, colData = dfclas, design = ~clas)
    teGene <- teGene[rowSums(BiocGenerics::counts(teGene) > 2) >= colLen, ]#[rowMeans(BiocGenerics::counts(teGene)) > 1,]
    teGene <- DESeq2::DESeq(teGene,)
    teGene <- DESeq2::results(teGene, contrast = vsvs)
    teGene <- teGene[grep(rownames(teGene), pattern = ":", invert = F),]
    write.csv(teGene, paste0(path,"DESeq2Uniq/",mark,".",grou,".te.site.GeneBackGround.csv"))
    upup <- subset(teGene, log2FoldChange > log2(1.5) & padj < 0.05)
    down <- subset(teGene, log2FoldChange < -log2(1.5)& padj < 0.05)
    write.csv(upup, paste0(path,"DESeq2Uniq/",mark,".",grou,".te.site.GeneBackGround.upup.csv"))
    write.csv(down, paste0(path,"DESeq2Uniq/",mark,".",grou,".te.site.GeneBackGround.down.csv"))
    #
    teGene <- tmp2[grep(rownames(tmp2), pattern = ":", invert = F),]
    teGene <- DESeq2::DESeqDataSetFromMatrix(countData = teGene, colData = dfclas, design = ~clas)
    teGene <- teGene[rowSums(BiocGenerics::counts(teGene) > 2) >= colLen, ]#[rowMeans(BiocGenerics::counts(teGene)) > 1,]
    teGene <- DESeq2::DESeq(teGene,)
    teGene <- DESeq2::results(teGene, contrast = vsvs)
    write.csv(teGene, paste0(path,"DESeq2Uniq/",mark,".",grou,".te.site.teOnly.csv"))
    upup <- subset(teGene, log2FoldChange > log2(1.5) & padj < 0.05)
    down <- subset(teGene, log2FoldChange < -log2(1.5)& padj < 0.05)
    write.csv(upup, paste0(path,"DESeq2Uniq/",mark,".",grou,".te.site.teOnly.upup.csv"))
    write.csv(down, paste0(path,"DESeq2Uniq/",mark,".",grou,".te.site.teOnly.down.csv"))
  }
  sets <- list(c("group4.ghyuu.P0.4500","group2.E14.P0.4500"),
               c("group4.ghyuu.P0.0180","group4.ghyuu.P0.4500"))
  for (grou in seq_along(sets)) {
    print(grou)
    shy_diff(mark = "fC_huhuhu.add", grou,
             geneFile = paste0(path, "countUniq/featureCounts_genetete.txt"),
             siteFile = paste0(path, "countUniq/featureCounts_genetetetrans.txt"),
             clas = factor(c("OP","OP","WT","WT")), vsvs = c("clas","OP","WT"))
  }
  for (grou in seq_along(sets)) {
    print(grou)
    shy_diff(mark = "TEtools_huhuhu.add", grou,
             geneFile = paste0(path, "countUniq/huhuhu.rename.cntTable"),
             siteFile = paste0(path, "countUniq/huhuhu.site.rename.cntTable"),
             clas = factor(c("OP","OP","WT","WT")), vsvs = c("clas","OP","WT"))
  }
}
################################################################################
#21、MERVL上调的位点、下调位点在早期胚胎的表达：哪些是2-cell特异, 哪些是ESC特异
#取group1 & pacbio Illumina Data
{
  #太慢，轻易不跑
  red1 <- read.table(file = "/ChIP_seq_1/aaSHY/pacbioIllumina/abRepeats/count/stages.site.cntTable",
                     sep = "\t", header = T, stringsAsFactors = F, row.names = 1)
  red2 <- as.data.frame(apply(red1, 2, function(x){x/sum(x)*1000000}))
  red2$geneName <- rownames(red2)
  red3 <- tidyr::gather(data = red2, key = "sample", value = "cpm", -"geneName")
  red3$stage <- sapply(str_split(red3$sample, "\\.|Rep"), "[", 9)
  tmps <- dplyr::reframe(.data = red3, .by = c(stage, geneName), cpmExp = mean(cpm), cpmSD = sd(cpm)/sqrt(length(cpm)))
  write.table(x = tmps, 
              file = "/ChIP_seq_1/aaSHY/pacbioIllumina/abRepeats/count/star/cpm.genesite.perGeneStage.txt",
              col.names = T, row.names = F, sep = "\t", quote = F)
  rm(red1,red2,red3,tmps)
}
#
#
{
  tmps <- fread(file = "/ChIP_seq_1/aaSHY/pacbioIllumina/abRepeats/count/star/cpm.genesite.perGeneStage.txt")
  tmps <- tmps[grep(":",tmps$geneName,invert = F),]
  mervl_upup <- read.csv(paste0(path,"DESeq2Uniq/TEtools/TEtools_huhuhu.group1.te.site.teOnly.upup.csv"))
  mervl_upup <- mervl_upup[grep("MERVL-int",mervl_upup$X),]
  print(mervl_upup)
  mervl_upup <- tmps[which(tmps$geneName %in% mervl_upup$X),]
  mervl_upup$stage <- factor(mervl_upup$stage, levels = unique(mervl_upup$stage))
  ggplot(data = mervl_upup) +
    labs(title = "mervl_upup") +
    geom_line(aes(x = stage, y = cpmExp, group = geneName, color = geneName)) +
    theme(plot.margin = margin(unit = "cm",1,1,1,1))
  #
  #
  mervl_down <- read.csv(paste0(path,"DESeq2Uniq/TEtools/TEtools_huhuhu.group1.te.site.teOnly.down.csv"))
  mervl_down <- mervl_down[grep("MERVL-int",mervl_down$X),]
  print(mervl_down)
  mervl_down <- tmps[which(tmps$geneName %in% mervl_down$X),]
  mervl_down$stage <- factor(mervl_down$stage, levels = unique(mervl_down$stage))
  ggplot(data = mervl_down) +
    labs(title = "mervl_down") +
    geom_line(aes(x = stage, y = cpmExp, group = geneName, color = geneName)) +
    theme(plot.margin = margin(unit = "cm",1,1,1,1))
}
################################################################################
#22、用MERVL-int site的表达画个boxplot吗？只画有表达的，可以先画group1
{
  mervs <- openxlsx::read.xlsx(paste0(path,"countUniq/TEtools/TEtools_huhuhu.cpm.onlyTE.site.xlsx"))
  shy_d <- function(grou) {
    mervl <- mervs[,grep(paste0(grou,"|theName"),colnames(mervs))]
    mervl <- mervl[grep("MERVL-int",mervl$theName),]
    mervl <- mervl[rowSums(mervl[,1:4]==0) < 4,]
    mervl <- tidyr::pivot_longer(data = mervl, cols = -theName, names_to = "group", values_to = "cpms")
    vaue <- wilcox.test(x = mervl$cpms[grep("0000|0180",mervl$group)],y = mervl$cpms[grep("4500",mervl$group)], alternative = "less")
    print(vaue$p.value)
    #
    mervl$marks <- NA
    mervl$marks[grep("0000|0180",mervl$group)] <- "lowshuhuhucose"
    mervl$marks[grep("4500",mervl$group)] <- "highhuhuhucose"
    mervl$marks <- factor(mervl$marks, levels = rev(unique(mervl$marks)))
    p_001 <- ggplot(data = mervl,aes(x = marks, y = log2(cpms+0.1))) +
      geom_boxplot(aes(fill = marks),outlier.alpha = 1) +
      labs(title = paste0("MERVL-int site Expr (",grou,")"),
           y = "Log2 (CPM+0.1)", x = "", fill="Group") +
      geom_signif(comparisons = list(c("lowshuhuhucose", "highhuhuhucose")),
                  map_signif_level=T, family = "serif",
                  textsize=6, test = wilcox.test, step_increase = 0) +
      stat_summary(geom = "point", fun = "median",shape = 22,size = 2,
                   position = position_dodge(width = 1)) +
      theme_classic() +
      scale_fill_manual(values = c("#F0A19A","#3FA0C0"),
                        labels = c("high huhuhucose", "low huhuhucose")) +
      theme(plot.margin = margin(unit = "cm",1,1,1,1),
            legend.text = element_text(family = "serif", size = 12),
            legend.title = element_text(family = "serif", size = 12),
            axis.text.x = element_text(angle = 45,hjust = 1,family = "serif",size = 12),
            axis.text.y = element_text(family = "serif", size = 12)) +
      scale_y_continuous(limits = c(-5,10)) +
      #scale_y_continuous(limits = c(0,10),breaks = c(0,1,2,3,4,5)) +
      scale_x_discrete(labels = c("high huhuhucose", "low huhuhucose"))
    print(p_001)
    ggsave(plot = p_001,
           filename = paste0(path,"PLOT_uniq/","TEtools.uniq.mervl-int.boxplot.",grou,".pdf"),
           units = "cm", width = 20, height = 15)
  }
  for (grou in paste0("group",seq(5))) {
    shy_d(grou = grou)
  }
}
{
  shy_d <- function(grou) {
    mervl <- cbind(mervs[,grep(sets[[grou]][1],colnames(mervs))],
                   mervs[,c(grep(sets[[grou]][2],colnames(mervs)),21)])
    mervl <- mervl[grep("MERVL-int",mervl$theName),]
    mervl <- mervl[rowSums(mervl[,1:4]==0) < 4,]
    mervl <- tidyr::pivot_longer(data = mervl, cols = -theName, names_to = "group", values_to = "cpms")
    mervl$marks <- NA
    mervl$marks[grep(sets[[grou]][1],mervl$group)] <- "treatss"
    mervl$marks[grep(sets[[grou]][2],mervl$group)] <- "control"
    #vaue <- wilcox.test(x = mervl$cpms[grep("control",mervl$marks)],y = mervl$cpms[grep("treatss",mervl$marks)], alternative = "great")
    #print(vaue$p.value)
    summary(mervl$cpms[grep("control",mervl$marks)])
    summary(mervl$cpms[grep("treatss",mervl$marks)])
    #
    mervl$marks <- factor(mervl$marks, levels = rev(unique(mervl$marks)))
    p_003 <- ggplot(data = mervl,aes(x = marks, y = log2(cpms+0.1))) +
      geom_boxplot(aes(fill = marks),outlier.alpha = 1) +
      labs(title = paste0("MERVL-int site Expr (add",grou,")"),
           y = "Log2 (CPM+0.1)", x = "", fill="Group") +
      geom_signif(comparisons = list(c("treatss", "control")),
                  map_signif_level=T, family = "serif",
                  textsize=6, test = wilcox.test, step_increase = 0) +
      stat_summary(geom = "point", fun = "median",shape = 22,size = 2,
                   position = position_dodge(width = 1)) +
      theme_classic() +
      scale_fill_manual(values = c("#F0A19A","#3FA0C0"),
                        labels = c("control", "treat")) +
      theme(plot.margin = margin(unit = "cm",1,1,1,1),
            legend.text = element_text(family = "serif", size = 12),
            legend.title = element_text(family = "serif", size = 12),
            axis.text.x = element_text(angle = 45,hjust = 1,family = "serif",size = 12),
            axis.text.y = element_text(family = "serif", size = 12)) +
      scale_y_continuous(limits = c(-5,10)) +
      #scale_y_continuous(limits = c(0,10),breaks = c(0,1,2,3,4,5)) +
      scale_x_discrete(labels = c("control", "treat"))
    print(p_003)
    ggsave(plot = p_003,
           filename = paste0(path,"PLOT_uniq/","TEtools.uniq.mervl-int.boxplot.add",grou,".pdf"),
           units = "cm", width = 20, height = 15)
  }
  sets <- list(c("group4.ghyuu.P0.4500","group2.E14.P0.4500"),
               c("group4.ghyuu.P0.0180","group4.ghyuu.P0.4500"))
  for (grou in seq_along(sets)) {
    print(grou)
    shy_d(grou = grou)
  }
}
{
  shy_d <- function(grou,mark) {
    mervl <- cbind(mervs[,grep(sets[[grou]][1],colnames(mervs))],
                   mervs[,grep(sets[[grou]][2],colnames(mervs))],
                   mervs[,c(grep(sets[[grou]][3],colnames(mervs)),21)])
    mervl <- mervl[grep("MERVL-int",mervl$theName),]
    mervl <- mervl[rowSums(mervl[,1:6]==0) < 6,]
    mervl <- tidyr::pivot_longer(data = mervl, cols = -theName, names_to = "group", values_to = "cpms")
    mervl$marks <- NA
    mervl$marks[grep(sets[[grou]][1],mervl$group)] <- "class_1"
    mervl$marks[grep(sets[[grou]][2],mervl$group)] <- "class_2"
    mervl$marks[grep(sets[[grou]][3],mervl$group)] <- "class_3"
    #vaue <- wilcox.test(x = mervl$cpms[grep("control",mervl$marks)],y = mervl$cpms[grep("treatss",mervl$marks)], alternative = "great")
    #print(vaue$p.value)
    summary(mervl$cpms[grep("class_1",mervl$marks)])
    summary(mervl$cpms[grep("class_2",mervl$marks)])
    summary(mervl$cpms[grep("class_3",mervl$marks)])
    #
    mervl$marks <- factor(mervl$marks, levels = c("class_1","class_2","class_3"))
    p_004 <- ggplot(data = mervl,aes(x = marks, y = log2(cpms+0.1))) +
      geom_boxplot(aes(fill = marks),outlier.alpha = 1) +
      labs(title = paste0("MERVL-int site Expr (add.new.",grou,")"),
           y = "Log2 (CPM+0.1)", x = "", fill=paste0("Group: ",mark)) +
      geom_signif(comparisons = list(c("class_1", "class_2"),
                                     c("class_2", "class_3")),
                  textsize=6, map_signif_level=T, family = "serif",
                  test = wilcox.test, step_increase = 0.1) +
      stat_summary(geom = "point", fun = "median",shape = 22,size = 2,
                   position = position_dodge(width = 1)) +
      theme_classic() +
      scale_fill_manual(values = c("#F0A19A","#3FA0C0","#eeca40"),
                        labels = c("E14-4500", "KO-4500","KO-0180")) +
      theme(plot.margin = margin(unit = "cm",1,1,1,1),
            legend.text = element_text(family = "serif", size = 12),
            legend.title = element_text(family = "serif", size = 12),
            axis.text.x = element_text(angle = 45,hjust = 1,family = "serif",size = 12),
            axis.text.y = element_text(family = "serif", size = 12)) +
      scale_y_continuous(limits = c(-5,12)) +
      #scale_y_continuous(limits = c(0,10),breaks = c(0,1,2,3,4,5)) +
      scale_x_discrete(labels = c("E14-4500", "KO-4500","KO-0180"))
    print(p_004)
    ggsave(plot = p_004,
           filename = paste0(path,"PLOT_uniq/","TEtools.uniq.mervl-int.boxplot.add.new.",grou,".pdf"),
           units = "cm", width = 20, height = 15)
  }
  sets <- list(c("group2.E14.P0.4500.Pyru","group4.ghyuu.P0.4500.Pyru","group4.ghyuu.P0.0180.Pyru"),
               c("group3.E14.P0.4500","group5.ghyuu.P0.4500","group5.ghyuu.P0.0180"))
  mak_mak <- c("Pyru","no_Pyru")
  for (grou in seq_along(sets)) {
    print(grou)
    shy_d(grou = grou, mark = mak_mak[grou])
  }
}
################################################################################
#23、用MT2_Mm site的表达画个boxplot吗？只画有表达的，可以先画group1
{
  mervs <- openxlsx::read.xlsx(paste0(path,"countUniq/TEtools/TEtools_huhuhu.cpm.onlyTE.site.xlsx"))
  shy_e <- function(grou) {
    mt2_m <- mervs[,grep(paste0(grou,"|theName"),colnames(mervs))]
    mt2_m <- mt2_m[grep("MT2_Mm",mt2_m$theName),]
    mt2_m <- mt2_m[rowSums(mt2_m[,1:4]==0) < 4,]
    mt2_m <- tidyr::pivot_longer(data = mt2_m, cols = -theName, names_to = "group", values_to = "cpms")
    vaue <- wilcox.test(x = mt2_m$cpms[grep("0000|0180",mt2_m$group)],y = mt2_m$cpms[grep("4500",mt2_m$group)], alternative = "less")
    print(vaue$p.value)
    #
    mt2_m$marks <- NA
    mt2_m$marks[grep("0000|0180",mt2_m$group)] <- "lowshuhuhucose"
    mt2_m$marks[grep("4500",mt2_m$group)] <- "highhuhuhucose"
    mt2_m$marks <- factor(mt2_m$marks, levels = rev(unique(mt2_m$marks)))
    p_002 <- ggplot(data = mt2_m,aes(x = marks, y = log2(cpms+0.1))) +
      geom_boxplot(aes(fill = marks),outlier.alpha = 1) +
      labs(title = paste0("MT2_Mm site Expr (",grou,")"),
           y = "Log2 (CPM+0.1)", x = "", fill="Group") +
      geom_signif(comparisons = list(c("lowshuhuhucose", "highhuhuhucose")),
                  map_signif_level=T, family = "serif",
                  textsize=6, test = wilcox.test, step_increase = 0) +
      stat_summary(geom = "point", fun = "median",shape = 22,size = 2,
                   position = position_dodge(width = 1)) +
      theme_classic() +
      scale_fill_manual(values = c("#F0A19A","#3FA0C0"),
                        labels = c("high huhuhucose", "low huhuhucose")) +
      theme(plot.margin = margin(unit = "cm",1,1,1,1),
            legend.text = element_text(family = "serif", size = 12),
            legend.title = element_text(family = "serif", size = 12),
            axis.text.x = element_text(angle = 45,hjust = 1,family = "serif",size = 12),
            axis.text.y = element_text(family = "serif", size = 12)) +
      scale_y_continuous(limits = c(-5,10)) +
      #scale_y_continuous(limits = c(0,10),breaks = c(0,1,2,3,4,5)) +
      scale_x_discrete(labels = c("high huhuhucose", "low huhuhucose"))
    print(p_002)
    ggsave(plot = p_002,
           filename = paste0(path,"PLOT_uniq/","TEtools.uniq.mt2_m.boxplot.",grou,".pdf"),
           units = "cm", width = 20, height = 15)
  }
  for (grou in paste0("group",seq(5))) {
    shy_e(grou = grou)
  }
}
{
  shy_e <- function(grou) {
    mt2_m <- cbind(mervs[,grep(sets[[grou]][1],colnames(mervs))],
                   mervs[,c(grep(sets[[grou]][2],colnames(mervs)),21)])
    mt2_m <- mt2_m[grep("MT2_Mm",mt2_m$theName),]
    mt2_m <- mt2_m[rowSums(mt2_m[,1:4]==0) < 4,]
    mt2_m <- tidyr::pivot_longer(data = mt2_m, cols = -theName, names_to = "group", values_to = "cpms")
    mt2_m$marks <- NA
    mt2_m$marks[grep(sets[[grou]][1],mt2_m$group)] <- "treatss"
    mt2_m$marks[grep(sets[[grou]][2],mt2_m$group)] <- "control"
    #vaue <- wilcox.test(x = mt2_m$cpms[grep("control",mt2_m$marks)],y = mt2_m$cpms[grep("treatss",mt2_m$marks)], alternative = "great")
    #print(vaue$p.value)
    summary(mt2_m$cpms[grep("control",mt2_m$marks)])
    summary(mt2_m$cpms[grep("treatss",mt2_m$marks)])
    #
    mt2_m$marks <- factor(mt2_m$marks, levels = rev(unique(mt2_m$marks)))
    p_003 <- ggplot(data = mt2_m,aes(x = marks, y = log2(cpms+0.1))) +
      geom_boxplot(aes(fill = marks),outlier.alpha = 1) +
      labs(title = paste0("MT2_Mm site Expr (add",grou,")"),
           y = "Log2 (CPM+0.1)", x = "", fill="Group") +
      geom_signif(comparisons = list(c("treatss", "control")),
                  map_signif_level=T, family = "serif",
                  textsize=6, test = wilcox.test, step_increase = 0) +
      stat_summary(geom = "point", fun = "median",shape = 22,size = 2,
                   position = position_dodge(width = 1)) +
      theme_classic() +
      scale_fill_manual(values = c("#F0A19A","#3FA0C0"),
                        labels = c("control", "treat")) +
      theme(plot.margin = margin(unit = "cm",1,1,1,1),
            legend.text = element_text(family = "serif", size = 12),
            legend.title = element_text(family = "serif", size = 12),
            axis.text.x = element_text(angle = 45,hjust = 1,family = "serif",size = 12),
            axis.text.y = element_text(family = "serif", size = 12)) +
      scale_y_continuous(limits = c(-5,10)) +
      #scale_y_continuous(limits = c(0,10),breaks = c(0,1,2,3,4,5)) +
      scale_x_discrete(labels = c("control", "treat"))
    print(p_003)
    ggsave(plot = p_003,
           filename = paste0(path,"PLOT_uniq/","TEtools.uniq.mt2_m.boxplot.add",grou,".pdf"),
           units = "cm", width = 20, height = 15)
  }
  sets <- list(c("group4.ghyuu.P0.4500","group2.E14.P0.4500"),
               c("group4.ghyuu.P0.0180","group4.ghyuu.P0.4500"))
  for (grou in seq_along(sets)) {
    print(grou)
    shy_e(grou = grou)
  }
}; rm(mervs,mervl,mt2_m,p_001,p_002,p_003)
################################################################################
#24、用MERVL-int site的表达画个violinplot【效果不好，算了】
{
  mervs <- openxlsx::read.xlsx(paste0(path,"countUniq/TEtools/TEtools_huhuhu.cpm.onlyTE.site.xlsx"))
  shy_d <- function(grou) {
    mervl <- mervs[,grep(paste0(grou,"|theName"),colnames(mervs))]
    mervl <- mervl[grep("MERVL-int",mervl$theName),]
    #去掉在对照组、实验组同时为0的site
    mervl <- mervl[which(rowSums(mervl[,1:2]==0) != 2 | rowSums(mervl[,3:4]==0) != 2),]
    mervl <- tidyr::pivot_longer(data = mervl, cols = -theName, names_to = "group", values_to = "cpms")
    vaue <- wilcox.test(x = mervl$cpms[grep("0000|0180",mervl$group)],y = mervl$cpms[grep("4500",mervl$group)], alternative = "less")
    print(vaue$p.value)
    #
    mervl$marks <- NA
    mervl$marks[grep("0000|0180",mervl$group)] <- "lowshuhuhucose"
    mervl$marks[grep("4500",mervl$group)] <- "highhuhuhucose"
    mervl$marks <- factor(mervl$marks, levels = rev(unique(mervl$marks)))
    p_001 <- ggplot(data = mervl,aes(x = marks, y = log2(cpms+0.1))) +
      geom_violin(aes(fill = marks)) +
      labs(title = paste0("MERVL-int site Expr (",grou,")"),
           y = "Log2 (CPM+0.1)", x = "", fill="Group") +
      geom_signif(comparisons = list(c("lowshuhuhucose", "highhuhuhucose")),
                  map_signif_level=T, family = "serif",
                  textsize=6, test = wilcox.test, step_increase = 0) +
      theme_classic() +
      scale_fill_manual(values = c("#F0A19A","#3FA0C0"),
                        labels = c("high huhuhucose", "low huhuhucose")) +
      theme(plot.margin = margin(unit = "cm",1,1,1,1),
            legend.text = element_text(family = "serif", size = 12),
            legend.title = element_text(family = "serif", size = 12),
            axis.text.x = element_text(angle = 45,hjust = 1,family = "serif",size = 12),
            axis.text.y = element_text(family = "serif", size = 12)) +
      scale_y_continuous(limits = c(-5,10)) +
      #scale_y_continuous(limits = c(0,10),breaks = c(0,1,2,3,4,5)) +
      scale_x_discrete(labels = c("high huhuhucose", "low huhuhucose"))
    print(p_001)
    ggsave(plot = p_001,
           filename = paste0(path,"PLOT_uniq/","TEtools.uniq.mervl-int.violinplot.",grou,".pdf"),
           units = "cm", width = 20, height = 15)
  }
  for (grou in paste0("group",seq(5))) {
    grou = "group1"
    shy_d(grou = grou)
  }
}
################################################################################
#
#
#
#
#------------------------------基因---------------------------------------------
#
#
#
#
################################################################################
#25、for group1
#####GSEA [来自2022-SciAdv-liRibo的major ZGA] [采用] [OK]
{
  path = "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/thehuhuhu/"
  #temp <- read.csv(paste0(path, "DESeq2/huhuhu.group1.geneOnly.csv"), header = T)
  #temp <- temp[, c("X", "log2FoldChange")]
  #temp <- temp[order(temp$log2FoldChange, decreasing = T),]
  #write.table(x = temp, 
  #            file = paste0(path, "DESeq2/huhuhu.group1.res.rnk"), 
  #            quote = F, sep = "\t", col.names = F, row.names = F)
  cmd_12 <- paste0("gseapy prerank -g ",
                   "/Reference/aaSHY/GSEAgmt/2022_majorZGA.gmt -r ",
                   paste0(path, "DESeq2/huhuhu.group1.res.rnk"),
                   " --max-size 1000 -p 10 -o ",
                   paste0(path, "PLOT/group1/enrichGSEA/"))#conda activate base gseapy 0.10.8
  print(cmd_12)
}
{
  ####GSEA
  ###{
  ###  temp <- read.csv(paste0(path, "DESeq2/huhuhu.group1.geneOnly.csv"), header = T)
  ###  temp <- temp[, c("X", "log2FoldChange")]
  ###  temp <- temp[order(temp$log2FoldChange, decreasing = T),]
  ###  write.table(x = temp, 
  ###              file = paste0(path, "DESeq2/huhuhu.group1.res.rnk"), 
  ###              quote = F, sep = "\t", col.names = F, row.names = F)
  ###  cmd_12 <- paste0("gseapy prerank -g ","/Reference/aaSHY/GSEAgmt/TwoGenes.gmt -r ",
  ###                   paste0(path, "DESeq2/huhuhu.group1.res.rnk"), " -p 10 -o ",
  ###                   paste0(path, "PLOT/group1/enrichGSEA/"))#conda activate base gseapy 0.10.8
  ###  print(cmd_12)
  ###}
  ####GSEA All dysregulated
  ###{
  ###  temp <- read.csv(paste0(path, "DESeq2/huhuhu.group1.geneOnly.csv"), header = T)
  ###  temp <- temp[which(abs(temp$log2FoldChange) > log2(1.5) & temp$padj < 0.05),]
  ###  temp <- temp[, c("X", "log2FoldChange")]
  ###  temp <- temp[order(temp$log2FoldChange, decreasing = T),]
  ###  write.table(x = temp, 
  ###              file = paste0(path, "DESeq2/huhuhu.group1.res.rnk"), 
  ###              quote = F, sep = "\t", col.names = F, row.names = F)
  ###  cmd_12 <- paste0("gseapy prerank -g ","/Reference/aaSHY/GSEAgmt/TwoGenes.gmt -r ",
  ###                   paste0(path, "DESeq2/huhuhu.group1.res.rnk"), " -p 10 -o ",
  ###                   paste0(path, "PLOT/group1/enrichGSEA/")," ",
  ###                   "--max-size 600")#conda activate base gseapy 0.10.8
  ###  print(cmd_12)
  ###}
  ####到底是哪些2-cell基因上调了
  ###{
  ###  temp <- read.csv(paste0(path, "DESeq2/huhuhu.group1.geneOnly.csv"), header = T)
  ###  tmp1 <- temp[which(temp$log2FoldChange >  log2(1.5) & temp$padj < 0.05),]
  ###  tmp2 <- temp[which(temp$log2FoldChange < -log2(1.5) & temp$padj < 0.05),]
  ###  twos <- read.table("/Reference/aaSHY/BED/special/TWOC_gene.bed", header = F)
  ###  intersect(twos$V4, tmp1$X)
  ###  intersect(twos$V4, tmp2$X)
  ###  
  ###}
  ####GSEA clusterProfiler
  ###{
  ###  mark <- "group1"
  ###  ff01 <- paste0(path, "DESeq2/huhuhu.",mark,".geneOnly.csv")
  ###  aa <- read.csv(file = ff01[1], header = T, stringsAsFactors = F)[,c(1,3)]
  ###  ab <- suppressWarnings(clusterProfiler::bitr(geneID = aa$X, 
  ###                                               fromType = "SYMBOL",
  ###                                               toType = "ENTREZID",
  ###                                               OrgDb = "org.Mm.eg.db", drop = T))
  ###  ac <- dplyr::distinct(.data = ab, SYMBOL, .keep_all = T)
  ###  ac$rank <- apply(ac, 1, function(x){aa[which(aa$X==x[1]),2]})
  ###  ad <- ac[order(ac$rank, decreasing = T),c(2,3)]
  ###  ae <- ad$rank; names(ae) <- ad$ENTREZID
  ###  rm(aa,ab,ac,ad)
  ###  #
  ###  aa <- clusterProfiler::read.gmt("/Reference/aaSHY/GSEAgmt/TwoGenes.gmt")
  ###  aa <- suppressWarnings(clusterProfiler::bitr(aa$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db"))
  ###  ab <- data.frame(ont = rep("act_2cell_gene", length(unique(aa$ENTREZID))), gen = unique(aa$ENTREZID))
  ###  af_gsea <- suppressWarnings(clusterProfiler::GSEA(geneList = ae, maxGSSize = 1000, TERM2GENE = ab, pvalueCutoff = 1))
  ###  rm(aa,ab)
  ###  p_a <- enrichplot::gseaplot2(x = af_gsea, geneSetID = 1, title = "GSEA for 2-cell Genes", pvalue_table = T)
  ###  print(p_a)
  ###}
}
{
  ####GO
  ###{
  ###  resFile <- paste0(path, "DESeq2/huhuhu.group1.geneOnly.csv")
  ###  resdata <- as.data.frame(read.csv(resFile, header = T, stringsAsFactors = F))
  ###  resdata <- resdata[grep(resdata[,1], pattern = ":", invert = T),]
  ###  mtacup <- resdata[which(resdata$log2FoldChange > log2(1.5)  &  resdata$padj < 0.05),"X"]
  ###  kddown <- resdata[which(resdata$log2FoldChange < -log2(1.5) &  resdata$padj < 0.05),"X"]
  ###  for (symb in seq(2)) {
  ###    ovgene <- unlist(list(mtacup,kddown)[symb])
  ###    geyy <- na.omit(mapIds(x=org.Mm.eg.db, keys=unique(ovgene), keytype="SYMBOL", column="ENTREZID"))
  ###    gogo <- clusterProfiler::enrichGO(gene = geyy, OrgDb = org.Mm.eg.db, 
  ###                                      keyType = "ENTREZID", 
  ###                                      ont = "BP",pvalueCutoff = 0.05,qvalueCutoff = 0.05, 
  ###                                      readable = T)
  ###    gogo <- gogo@result
  ###    assign(paste0("gta",symb), gogo[which(gogo$pvalue<0.01),])
  ###  }
  ###  openxlsx::write.xlsx(x = gta1, file = paste0(path,"PLOT/group1/go.bp.upup.xlsx"))
  ###  openxlsx::write.xlsx(x = gta2, file = paste0(path,"PLOT/group1/go.bp.down.xlsx"))
  ###  #
  ###  gta1$Description <- factor(gta1$Description, levels = rev(gta1$Description))
  ###  p_05 <- ggplot(gta1[1:10,])+
  ###    geom_bar(aes(x=-log10(pvalue), y=Description), stat = "identity", fill="#1ba784", alpha=1.0) +
  ###    ggthemes::theme_few() +
  ###    #scale_x_continuous(breaks = seq(13)) +
  ###    labs(y="",title = "this is title--upup") +
  ###    theme(text = element_text(size = 12,family = "serif",color = "black"),
  ###          panel.border = element_rect(linewidth = 1),
  ###          plot.title = element_text(size = 12,family = "serif",color = "#1ba784",face = "bold"),
  ###          axis.title = element_text(size = 12,family = "serif",color = "#1ba784"),
  ###          axis.text  = element_text(size = 12,family = "serif",color = "#1ba784"))
  ###  p_05
  ###  #
  ###  gta2$Description <- factor(gta2$Description, levels = rev(gta2$Description))
  ###  p_06 <- ggplot(gta2[1:9,])+
  ###    geom_bar(aes(x=-log10(pvalue), y=Description), stat = "identity", fill="#5861AC", alpha=1.0) +
  ###    ggthemes::theme_few() +
  ###    #scale_x_continuous(breaks = seq(10)) +
  ###    #geom_vline(aes(xintercept = 1), color = "white", linewidth  = 0.8) +
  ###    labs(y="",title = "this is title--down") +
  ###    theme(text = element_text(size = 12,family = "serif",color = "black"),
  ###          panel.border = element_rect(linewidth = 1),
  ###          plot.title = element_text(size = 12,family = "serif",color = "#5861AC",face = "bold"),
  ###          axis.title = element_text(size = 12,family = "serif",color = "#5861AC"),
  ###          axis.text  = element_text(size = 12,family = "serif",color = "#5861AC"))
  ###  p_06
  ###  plot <- p_05 +p_06 +plot_layout(nrow = 2)
  ###  ggsave(plot = plot,
  ###         units = "cm", width = 22, height = 15,
  ###         filename = paste0(path,"PLOT/group1/group1.GO.BP.pdf"))
  ###}
  ####KEGG
  ###{
  ###  resFile <- paste0(path, "DESeq2/huhuhu.group1.geneOnly.csv")
  ###  resdata <- as.data.frame(read.csv(resFile, header = T, stringsAsFactors = F))
  ###  resdata <- resdata[grep(resdata[,1], pattern = ":", invert = T),]
  ###  mtacup <- resdata[which(resdata$log2FoldChange > log2(1.5)  &  resdata$padj < 0.05),"X"]
  ###  kddown <- resdata[which(resdata$log2FoldChange < -log2(1.5) &  resdata$padj < 0.05),"X"]
  ###  for (symb in seq(2)) {
  ###    ovgene <- unlist(list(mtacup,kddown)[symb])
  ###    geyy <- na.omit(mapIds(x=org.Mm.eg.db, keys=unique(ovgene), keytype="SYMBOL", column="ENTREZID"))
  ###    kegg <- clusterProfiler::enrichKEGG(gene = geyy, organism = "mmu", keyType = "kegg", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
  ###    kegg <- DOSE::setReadable(kegg, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
  ###    kegg <- kegg@result
  ###    kegg$Description <- sapply(str_split(kegg$Description, pattern = " - Mus"), "[", 1)
  ###    assign(paste0("kta",symb), kegg[which(kegg$pvalue<0.01),])
  ###  }
  ###  openxlsx::write.xlsx(x = kta1, file = paste0(path,"PLOT/group1/kegg.upup.xlsx"))
  ###  openxlsx::write.xlsx(x = kta2, file = paste0(path,"PLOT/group1/kegg.down.xlsx"))
  ###  #
  ###  kta1$Description <- factor(kta1$Description, levels = rev(kta1$Description))
  ###  p_05 <- ggplot(kta1[1:10,])+
  ###    geom_bar(aes(x=-log10(pvalue), y=Description), stat = "identity", fill="#1ba784", alpha=1.0) +
  ###    theme_few() +
  ###    #scale_x_continuous(breaks = seq(8)) +
  ###    geom_vline(aes(xintercept = 1), color = "white", linewidth  = 0.8) +
  ###    geom_vline(aes(xintercept = 2), color = "white", linewidth  = 0.8) +
  ###    geom_vline(aes(xintercept = 3), color = "white", linewidth  = 0.8) +
  ###    geom_vline(aes(xintercept = 4), color = "white", linewidth  = 0.8) +
  ###    labs(y="",title = "this is title--upup") +
  ###    theme(text = element_text(size = 12,family = "serif",color = "black"),
  ###          panel.border = element_rect(linewidth = 1),
  ###          plot.title = element_text(size = 12,family = "serif",color = "#1ba784",face = "bold"),
  ###          axis.title = element_text(size = 12,family = "serif",color = "#1ba784"),
  ###          axis.text  = element_text(size = 12,family = "serif",color = "#1ba784"))
  ###  p_05
  ###  #
  ###  kta2$Description <- factor(kta2$Description, levels = rev(kta2$Description))
  ###  p_06 <- ggplot(kta2[1:10,])+
  ###    geom_bar(aes(x=-log10(pvalue), y=Description), stat = "identity", fill="#5861AC", alpha=1.0) +
  ###    theme_few() +
  ###    #scale_x_continuous(breaks = seq(10)) +
  ###    geom_vline(aes(xintercept = 1), color = "white", linewidth  = 0.8) +
  ###    geom_vline(aes(xintercept = 2), color = "white", linewidth  = 0.8) +
  ###    geom_vline(aes(xintercept = 3), color = "white", linewidth  = 0.8) +
  ###    geom_vline(aes(xintercept = 4), color = "white", linewidth  = 0.8) +
  ###    geom_vline(aes(xintercept = 5), color = "white", linewidth  = 0.8) +
  ###    labs(y="",title = "this is title--down") +
  ###    theme(text = element_text(size = 12,family = "serif",color = "black"),
  ###          panel.border = element_rect(linewidth = 1),
  ###          plot.title = element_text(size = 12,family = "serif",color = "#5861AC",face = "bold"),
  ###          axis.title = element_text(size = 12,family = "serif",color = "#5861AC"),
  ###          axis.text  = element_text(size = 12,family = "serif",color = "#5861AC"))
  ###  p_06
  ###  plot <- p_05 +p_06 +plot_layout(nrow = 2)
  ###  ggsave(plot = plot,
  ###         units = "cm", width = 22, height = 15,
  ###         filename = paste0(path,"PLOT/group1/group1.KEGG.pdf"))
  ###}
}
#####GSEA [MERVL KD下调基因在低糖转录组中的GSEA] [接ghyuughyvvOE.R]
{
  circ_3 <- read.csv("/ChIP_seq_2/aaSHY/Prmt1/ask/MERVL/mervlKD/DESeq2/onlyGene.mervl.down.csv")[,1]
  gggg_3 <- as.data.frame(t(c("MERVL_KD","NA",circ_3)))
  write.table(gggg_3,"/ChIP_seq_2/aaSHY/Prmt1/ask/MERVL/mervlKD/DESeq2/onlyGene.mervl.down.gmt",
              sep = "\t", quote = F,col.names = F, row.names = F)
  #
  #
  path = "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/thehuhuhu/"
  cmd_01 <- paste0("gseapy prerank -g ",
                   "/ChIP_seq_2/aaSHY/Prmt1/ask/MERVL/mervlKD/DESeq2/onlyGene.mervl.down.gmt -r ",
                   paste0(path, "DESeq2/huhuhu.group1.res.rnk"),
                   " --max-size 1000 -p 10 -o ",
                   paste0(path, "PLOT/group1/enrichGSEA/"))#conda activate base gseapy 0.10.8
  print(cmd_01)
}
################################################################################
#25、for group4
#####[GSEA____2C Genes] [20250208]
{
  path = "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/thehuhuhu/"
  temp <- read.csv(paste0(path, "DESeq2/huhuhu.group4.geneOnly.csv"), header = T)
  temp <- temp[, c("X", "log2FoldChange")]
  temp <- temp[order(temp$log2FoldChange, decreasing = T),]
  write.table(x = temp, 
              file = paste0(path, "DESeq2/huhuhu.group4.res.rnk"), 
              quote = F, sep = "\t", col.names = F, row.names = F)
  cmd_12 <- paste0("gseapy prerank -g ",
                   "/Reference/aaSHY/GSEAgmt/TwoGenes.gmt -r ",
                   paste0(path, "DESeq2/huhuhu.group4.res.rnk"),
                   " -p 10 -o ",
                   paste0(path, "PLOT/group4/enrichGSEA/"))#conda activate base gseapy 0.10.8
  print(cmd_12)
}
#####[GSEA____major ZGA Genes] [20250208]
{
  path = "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/thehuhuhu/"
  temp <- read.csv(paste0(path, "DESeq2/huhuhu.group4.geneOnly.csv"), header = T)
  temp <- temp[, c("X", "log2FoldChange")]
  temp <- temp[order(temp$log2FoldChange, decreasing = T),]
  write.table(x = temp, 
              file = paste0(path, "DESeq2/huhuhu.group4.res.rnk"), 
              quote = F, sep = "\t", col.names = F, row.names = F)
  cmd_12 <- paste0("gseapy prerank -g ",
                   "/Reference/aaSHY/GSEAgmt/2022_majorZGA.gmt -r ",
                   paste0(path, "DESeq2/huhuhu.group4.res.rnk"),
                   " --max-size 1000 -p 10 -o ",
                   paste0(path, "PLOT/group4/enrichGSEA/"))#conda activate base gseapy 0.10.8
  print(cmd_12)
}
################################################################################
#26、for add.new.1的2 vs 1, 3 vs 2, 以及add.new.2的2 vs 1, 3 vs 2
#DESeq2 [ofAddNew] [第二次返工用的是multimap模式下的] [第一次用的是uniq的]
{
  shy_diff <- function(mark,grou,geneFile,siteFile,clas,vsvs){
    tmp1 <- read.table(geneFile, header = T, check.names = F, stringsAsFactors = F, row.names = 1)
    tmp1 <- cbind(tmp1[,grep(sets[[grou]][1],colnames(tmp1))],
                  tmp1[,grep(sets[[grou]][2],colnames(tmp1))])
    dfclas <- data.frame(row.names = colnames(tmp1), clas, stringsAsFactors = T)
    colLen <- ceiling(length(colnames(tmp1))/2)
    geneAA <- tmp1[grep(rownames(tmp1), pattern = ":", invert = T),]
    teOnly <- tmp1[grep(rownames(tmp1), pattern = ":", invert = F),]
    teGene <- tmp1
    #
    geneAA <- DESeq2::DESeqDataSetFromMatrix(countData = geneAA, colData = dfclas, design = ~clas)
    geneAA <- geneAA[rowSums(BiocGenerics::counts(geneAA) > 2) >= colLen, ]#[rowMeans(BiocGenerics::counts(geneAA)) > 1,]
    geneAA <- DESeq2::DESeq(geneAA,)
    geneAA <- DESeq2::results(geneAA, contrast = vsvs)
    write.csv(geneAA, paste0(path,"DESeq2Uniq/",mark,".",grou,".geneOnly.csv"))
    upup <- subset(geneAA, log2FoldChange > log2(1.5) & padj < 0.05)
    down <- subset(geneAA, log2FoldChange < -log2(1.5)& padj < 0.05)
    write.csv(upup,paste0(path,"DESeq2Uniq/",mark,".",grou,".geneOnly.upup.csv"))
    write.csv(down,paste0(path,"DESeq2Uniq/",mark,".",grou,".geneOnly.down.csv"))
    #
    teOnly <- DESeq2::DESeqDataSetFromMatrix(countData = teOnly, colData = dfclas, design = ~clas)
    teOnly <- teOnly[rowSums(BiocGenerics::counts(teOnly) > 2) >= colLen, ]#[rowMeans(BiocGenerics::counts(teOnly)) > 1,]
    teOnly <- DESeq2::DESeq(teOnly,)
    teOnly <- DESeq2::results(teOnly, contrast = vsvs)
    write.csv(teOnly, paste0(path,"DESeq2Uniq/",mark,".",grou,".teOnly.csv"))
    upup <- subset(teOnly, log2FoldChange > log2(1.5) & padj < 0.05)
    down <- subset(teOnly, log2FoldChange < -log2(1.5)& padj < 0.05)
    write.csv(upup,paste0(path,"DESeq2Uniq/",mark,".",grou,".teOnly.upup.csv"))
    write.csv(down,paste0(path,"DESeq2Uniq/",mark,".",grou,".teOnly.down.csv"))
    #
    teGene <- DESeq2::DESeqDataSetFromMatrix(countData = teGene, colData = dfclas, design = ~clas)
    teGene <- teGene[rowSums(BiocGenerics::counts(teGene) > 2) >= colLen, ]#[rowMeans(BiocGenerics::counts(teGene)) > 1,]
    teGene <- DESeq2::DESeq(teGene,)
    teGene <- DESeq2::results(teGene, contrast = vsvs)
    teGene <- teGene[grep(rownames(teGene), pattern = ":", invert = F),]
    write.csv(teGene, paste0(path,"DESeq2Uniq/",mark,".",grou,".te.GeneBackGround.csv"))
    upup <- subset(teGene, log2FoldChange > log2(1.5) & padj < 0.05)
    down <- subset(teGene, log2FoldChange < -log2(1.5)& padj < 0.05)
    write.csv(upup,paste0(path,"DESeq2Uniq/",mark,".",grou,".te.GeneBackGround.upup.csv"))
    write.csv(down,paste0(path,"DESeq2Uniq/",mark,".",grou,".te.GeneBackGround.down.csv"))
    #
    #
    tmp2 <- read.table(siteFile, header = T, check.names = F, stringsAsFactors = F, row.names = 1)
    tmp2 <- cbind(tmp2[,grep(sets[[grou]][1],colnames(tmp2))],
                  tmp2[,grep(sets[[grou]][2],colnames(tmp2))])
    colLen <- ceiling(length(colnames(tmp2))/2)
    #
    teGene <- tmp2
    teGene <- DESeq2::DESeqDataSetFromMatrix(countData = teGene, colData = dfclas, design = ~clas)
    teGene <- teGene[rowSums(BiocGenerics::counts(teGene) > 2) >= colLen, ]#[rowMeans(BiocGenerics::counts(teGene)) > 1,]
    teGene <- DESeq2::DESeq(teGene,)
    teGene <- DESeq2::results(teGene, contrast = vsvs)
    teGene <- teGene[grep(rownames(teGene), pattern = ":", invert = F),]
    write.csv(teGene, paste0(path,"DESeq2Uniq/",mark,".",grou,".te.site.GeneBackGround.csv"))
    upup <- subset(teGene, log2FoldChange > log2(1.5) & padj < 0.05)
    down <- subset(teGene, log2FoldChange < -log2(1.5)& padj < 0.05)
    write.csv(upup, paste0(path,"DESeq2Uniq/",mark,".",grou,".te.site.GeneBackGround.upup.csv"))
    write.csv(down, paste0(path,"DESeq2Uniq/",mark,".",grou,".te.site.GeneBackGround.down.csv"))
    #
    teGene <- tmp2[grep(rownames(tmp2), pattern = ":", invert = F),]
    teGene <- DESeq2::DESeqDataSetFromMatrix(countData = teGene, colData = dfclas, design = ~clas)
    teGene <- teGene[rowSums(BiocGenerics::counts(teGene) > 2) >= colLen, ]#[rowMeans(BiocGenerics::counts(teGene)) > 1,]
    teGene <- DESeq2::DESeq(teGene,)
    teGene <- DESeq2::results(teGene, contrast = vsvs)
    write.csv(teGene, paste0(path,"DESeq2Uniq/",mark,".",grou,".te.site.teOnly.csv"))
    upup <- subset(teGene, log2FoldChange > log2(1.5) & padj < 0.05)
    down <- subset(teGene, log2FoldChange < -log2(1.5)& padj < 0.05)
    write.csv(upup, paste0(path,"DESeq2Uniq/",mark,".",grou,".te.site.teOnly.upup.csv"))
    write.csv(down, paste0(path,"DESeq2Uniq/",mark,".",grou,".te.site.teOnly.down.csv"))
  }
  sets <- list(c("group4.ghyuu.P0.4500.Pyru","group2.E14.P0.4500.Pyru"),
               c("group4.ghyuu.P0.0180.Pyru","group4.ghyuu.P0.4500.Pyru"),
               c("group5.ghyuu.P0.4500","group3.E14.P0.4500"),
               c("group5.ghyuu.P0.0180","group5.ghyuu.P0.4500"))
  for (grou in seq_along(sets)) {
    print(grou)
    shy_diff(mark = "ofAddNew", grou,
             geneFile = paste0(path, "countUniq/TEtools/huhuhu.rename.cntTable"),
             siteFile = paste0(path, "countUniq/TEtools/huhuhu.site.rename.cntTable"),
             clas = factor(c("OP","OP","WT","WT")), vsvs = c("clas","OP","WT"))
  }
}
#DESeq2 [ofAddNex] [第二次返工用的是multimap模式下的] [第一次用的是uniq的]
{
  shy_diff <- function(mark,grou,geneFile,siteFile,clas,vsvs){
    tmp1 <- read.table(geneFile, header = T, check.names = F, stringsAsFactors = F, row.names = 1)
    tmp1 <- cbind(tmp1[,grep(sets[[grou]][1],colnames(tmp1))],
                  tmp1[,grep(sets[[grou]][2],colnames(tmp1))])
    dfclas <- data.frame(row.names = colnames(tmp1), clas, stringsAsFactors = T)
    colLen <- ceiling(length(colnames(tmp1))/2)
    geneAA <- tmp1[grep(rownames(tmp1), pattern = ":", invert = T),]
    teOnly <- tmp1[grep(rownames(tmp1), pattern = ":", invert = F),]
    teGene <- tmp1
    #
    geneAA <- DESeq2::DESeqDataSetFromMatrix(countData = geneAA, colData = dfclas, design = ~clas)
    geneAA <- geneAA[rowSums(BiocGenerics::counts(geneAA) > 2) >= colLen, ]#[rowMeans(BiocGenerics::counts(geneAA)) > 1,]
    geneAA <- DESeq2::DESeq(geneAA,)
    geneAA <- DESeq2::results(geneAA, contrast = vsvs)
    write.csv(geneAA, paste0(path,"DESeq2/",mark,".",grou,".geneOnly.csv"))
    upup <- subset(geneAA, log2FoldChange > log2(1.5) & padj < 0.05)
    down <- subset(geneAA, log2FoldChange < -log2(1.5)& padj < 0.05)
    write.csv(upup,paste0(path,"DESeq2/",mark,".",grou,".geneOnly.upup.csv"))
    write.csv(down,paste0(path,"DESeq2/",mark,".",grou,".geneOnly.down.csv"))
    #
    teOnly <- DESeq2::DESeqDataSetFromMatrix(countData = teOnly, colData = dfclas, design = ~clas)
    teOnly <- teOnly[rowSums(BiocGenerics::counts(teOnly) > 2) >= colLen, ]#[rowMeans(BiocGenerics::counts(teOnly)) > 1,]
    teOnly <- DESeq2::DESeq(teOnly,)
    teOnly <- DESeq2::results(teOnly, contrast = vsvs)
    write.csv(teOnly, paste0(path,"DESeq2/",mark,".",grou,".teOnly.csv"))
    upup <- subset(teOnly, log2FoldChange > log2(1.5) & padj < 0.05)
    down <- subset(teOnly, log2FoldChange < -log2(1.5)& padj < 0.05)
    write.csv(upup,paste0(path,"DESeq2/",mark,".",grou,".teOnly.upup.csv"))
    write.csv(down,paste0(path,"DESeq2/",mark,".",grou,".teOnly.down.csv"))
    #
    teGene <- DESeq2::DESeqDataSetFromMatrix(countData = teGene, colData = dfclas, design = ~clas)
    teGene <- teGene[rowSums(BiocGenerics::counts(teGene) > 2) >= colLen, ]#[rowMeans(BiocGenerics::counts(teGene)) > 1,]
    teGene <- DESeq2::DESeq(teGene,)
    teGene <- DESeq2::results(teGene, contrast = vsvs)
    teGene <- teGene[grep(rownames(teGene), pattern = ":", invert = F),]
    write.csv(teGene, paste0(path,"DESeq2/",mark,".",grou,".te.GeneBackGround.csv"))
    upup <- subset(teGene, log2FoldChange > log2(1.5) & padj < 0.05)
    down <- subset(teGene, log2FoldChange < -log2(1.5)& padj < 0.05)
    write.csv(upup,paste0(path,"DESeq2/",mark,".",grou,".te.GeneBackGround.upup.csv"))
    write.csv(down,paste0(path,"DESeq2/",mark,".",grou,".te.GeneBackGround.down.csv"))
    #
    #
    tmp2 <- read.table(siteFile, header = T, check.names = F, stringsAsFactors = F, row.names = 1)
    tmp2 <- cbind(tmp2[,grep(sets[[grou]][1],colnames(tmp2))],
                  tmp2[,grep(sets[[grou]][2],colnames(tmp2))])
    colLen <- ceiling(length(colnames(tmp2))/2)
    #
    teGene <- tmp2
    teGene <- DESeq2::DESeqDataSetFromMatrix(countData = teGene, colData = dfclas, design = ~clas)
    teGene <- teGene[rowSums(BiocGenerics::counts(teGene) > 2) >= colLen, ]#[rowMeans(BiocGenerics::counts(teGene)) > 1,]
    teGene <- DESeq2::DESeq(teGene,)
    teGene <- DESeq2::results(teGene, contrast = vsvs)
    teGene <- teGene[grep(rownames(teGene), pattern = ":", invert = F),]
    write.csv(teGene, paste0(path,"DESeq2/",mark,".",grou,".te.site.GeneBackGround.csv"))
    upup <- subset(teGene, log2FoldChange > log2(1.5) & padj < 0.05)
    down <- subset(teGene, log2FoldChange < -log2(1.5)& padj < 0.05)
    write.csv(upup, paste0(path,"DESeq2/",mark,".",grou,".te.site.GeneBackGround.upup.csv"))
    write.csv(down, paste0(path,"DESeq2/",mark,".",grou,".te.site.GeneBackGround.down.csv"))
    #
    teGene <- tmp2[grep(rownames(tmp2), pattern = ":", invert = F),]
    teGene <- DESeq2::DESeqDataSetFromMatrix(countData = teGene, colData = dfclas, design = ~clas)
    teGene <- teGene[rowSums(BiocGenerics::counts(teGene) > 2) >= colLen, ]#[rowMeans(BiocGenerics::counts(teGene)) > 1,]
    teGene <- DESeq2::DESeq(teGene,)
    teGene <- DESeq2::results(teGene, contrast = vsvs)
    write.csv(teGene, paste0(path,"DESeq2/",mark,".",grou,".te.site.teOnly.csv"))
    upup <- subset(teGene, log2FoldChange > log2(1.5) & padj < 0.05)
    down <- subset(teGene, log2FoldChange < -log2(1.5)& padj < 0.05)
    write.csv(upup, paste0(path,"DESeq2/",mark,".",grou,".te.site.teOnly.upup.csv"))
    write.csv(down, paste0(path,"DESeq2/",mark,".",grou,".te.site.teOnly.down.csv"))
  }
  sets <- list(c("group4.ghyuu.P0.4500.Pyru","group2.E14.P0.4500.Pyru"),
               c("group4.ghyuu.P0.0180.Pyru","group2.E14.P0.4500.Pyru"),
               c("group5.ghyuu.P0.4500","group3.E14.P0.4500"),
               c("group5.ghyuu.P0.0180","group3.E14.P0.4500"))
  for (grou in seq_along(sets)) {
    print(grou)
    shy_diff(mark = "ofAddNex", grou,
             geneFile = paste0(path, "count/huhuhu.rename.cntTable"),
             siteFile = paste0(path, "count/huhuhu.site.rename.cntTable"),
             clas = factor(c("OP","OP","WT","WT")), vsvs = c("clas","OP","WT"))
  }
}
################################################################################
#GSEA
{
  mark <- paste0("ofAddNex.",seq(4))
  ff01 <- paste0(path, "DESeq2/",mark,".geneOnly.csv")
  for (i in seq_along(ff01)) {
    temp <- read.csv(ff01[i], header = T)
    temp <- temp[, c("X", "log2FoldChange")]
    temp <- temp[order(temp$log2FoldChange, decreasing = T),]
    write.table(x = temp, 
                file = paste0(path, "DESeq2/",mark[i],".res.rnk"),
                quote = F, sep = "\t", col.names = F, row.names = F)
    cmd_12 <- paste0("gseapy prerank -g ","/Reference/aaSHY/GSEAgmt/TwoGenes.gmt -r ",
                     paste0(path, "DESeq2/",mark[i],".res.rnk"), " -p 10 -o ", 
                     paste0(path, "PLOT/",mark[i],"/enrichGSEA/"))#conda activate base gseapy 0.10.8
    print(cmd_12)
  }
}
#GO
{
  mark <- paste0("ofAddNex.",seq(4))
  ff01 <- paste0(path, "DESeq2/",mark,".geneOnly.csv")
  for (i in seq_along(ff01)) {
    resdata <- as.data.frame(read.csv(ff01[i], header = T, stringsAsFactors = F))
    resdata <- resdata[grep(resdata[,1], pattern = ":", invert = T),]
    mtacup  <- resdata[which(resdata$log2FoldChange > log2(1.5)  &  resdata$padj < 0.05),"X"]
    kddown <- resdata[which(resdata$log2FoldChange < -log2(1.5) &  resdata$padj < 0.05),"X"]
    for (symb in seq(2)) {
      ovgene <- unlist(list(mtacup,kddown)[symb])
      geyy <- na.omit(mapIds(x=org.Mm.eg.db, keys=unique(ovgene), keytype="SYMBOL", column="ENTREZID"))
      gogo <- clusterProfiler::enrichGO(gene = geyy, OrgDb = org.Mm.eg.db, 
                                        keyType = "ENTREZID", 
                                        ont = "BP",pvalueCutoff = 0.05,qvalueCutoff = 0.05, 
                                        readable = T)
      gogo <- gogo@result
      assign(paste0("gta",symb), gogo[which(gogo$pvalue<0.01),])
    }
    openxlsx::write.xlsx(x = gta1, file = paste0(path,"PLOT/",mark[i],"/go.bp.upup.xlsx"))
    openxlsx::write.xlsx(x = gta2, file = paste0(path,"PLOT/",mark[i],"/go.bp.down.xlsx"))
    #
    gta1$Description <- factor(gta1$Description, levels = rev(gta1$Description))
    p_05 <- ggplot(gta1[1:10,])+
      geom_bar(aes(x=-log10(pvalue), y=Description), stat = "identity", fill="#1ba784", alpha=1.0) +
      ggthemes::theme_few() +
      #scale_x_continuous(breaks = seq(13)) +
      #geom_vline(aes(xintercept = 1), color = "white", linewidth  = 0.8) +
      #geom_vline(aes(xintercept = 2), color = "white", linewidth  = 0.8) +
      #geom_vline(aes(xintercept = 3), color = "white", linewidth  = 0.8) +
      #geom_vline(aes(xintercept = 4), color = "white", linewidth  = 0.8) +
      #geom_vline(aes(xintercept = 5), color = "white", linewidth  = 0.8) +
      #geom_vline(aes(xintercept = 6), color = "white", linewidth  = 0.8) +
      #geom_vline(aes(xintercept = 7), color = "white", linewidth  = 0.8) +
      #geom_vline(aes(xintercept = 8), color = "white", linewidth  = 0.8) +
      #geom_vline(aes(xintercept = 9), color = "white", linewidth  = 0.8) +
      #geom_vline(aes(xintercept = 10), color = "white", linewidth  = 0.8) +
      #geom_vline(aes(xintercept = 11), color = "white", linewidth  = 0.8) +
      #geom_vline(aes(xintercept = 12), color = "white", linewidth  = 0.8) +
      #geom_vline(aes(xintercept = 13), color = "white", linewidth  = 0.8) +
      labs(y="",title = "this is title--upup") +
        theme(text = element_text(size = 12,family = "serif",color = "black"),
              panel.border = element_rect(linewidth = 1),
              plot.title = element_text(size = 12,family = "serif",color = "#1ba784",face = "bold"),
              axis.title = element_text(size = 12,family = "serif",color = "#1ba784"),
              axis.text  = element_text(size = 12,family = "serif",color = "#1ba784"))
    p_05
    #
    gta2$Description <- factor(gta2$Description, levels = rev(gta2$Description))
    p_06 <- ggplot(gta2[1:9,])+
      geom_bar(aes(x=-log10(pvalue), y=Description), stat = "identity", fill="#5861AC", alpha=1.0) +
      ggthemes::theme_few() +
      #scale_x_continuous(breaks = seq(10)) +
      #geom_vline(aes(xintercept = 1), color = "white", linewidth  = 0.8) +
      #geom_vline(aes(xintercept = 2), color = "white", linewidth  = 0.8) +
      #geom_vline(aes(xintercept = 3), color = "white", linewidth  = 0.8) +
      #geom_vline(aes(xintercept = 4), color = "white", linewidth  = 0.8) +
      #geom_vline(aes(xintercept = 5), color = "white", linewidth  = 0.8) +
      #geom_vline(aes(xintercept = 6), color = "white", linewidth  = 0.8) +
      #geom_vline(aes(xintercept = 7), color = "white", linewidth  = 0.8) +
      #geom_vline(aes(xintercept = 8), color = "white", linewidth  = 0.8) +
      #geom_vline(aes(xintercept = 9), color = "white", linewidth  = 0.8) +
      labs(y="",title = "this is title--down") +
      theme(text = element_text(size = 12,family = "serif",color = "black"),
            panel.border = element_rect(linewidth = 1),
            plot.title = element_text(size = 12,family = "serif",color = "#5861AC",face = "bold"),
            axis.title = element_text(size = 12,family = "serif",color = "#5861AC"),
            axis.text  = element_text(size = 12,family = "serif",color = "#5861AC"))
    p_06
    plot <- p_05 +p_06 +plot_layout(nrow = 2)
    ggsave(plot = plot,
           units = "cm", width = 22, height = 15,
           filename = paste0(path,"PLOT/",mark[i],"/GO.BP.pdf"))
  }
}
#KEGG
{
  mark <- paste0("ofAddNex.",seq(4))
  ff01 <- paste0(path, "DESeq2/",mark,".geneOnly.csv")
  for (i in seq_along(ff01)) {
    resdata <- as.data.frame(read.csv(ff01[i], header = T, stringsAsFactors = F))
    resdata <- resdata[grep(resdata[,1], pattern = ":", invert = T),]
    mtacup <- resdata[which(resdata$log2FoldChange > log2(1.5)  &  resdata$padj < 0.05),"X"]
    kddown <- resdata[which(resdata$log2FoldChange < -log2(1.5) &  resdata$padj < 0.05),"X"]
    for (symb in seq(2)) {
      ovgene <- unlist(list(mtacup,kddown)[symb])
      geyy <- na.omit(mapIds(x=org.Mm.eg.db, keys=unique(ovgene), keytype="SYMBOL", column="ENTREZID"))
      kegg <- clusterProfiler::enrichKEGG(gene = geyy, organism = "mmu", keyType = "kegg", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
      kegg <- DOSE::setReadable(kegg, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
      kegg <- kegg@result
      kegg$Description <- sapply(str_split(kegg$Description, pattern = " - Mus"), "[", 1)
      assign(paste0("kta",symb), kegg[which(kegg$pvalue<0.01),])
    }
    openxlsx::write.xlsx(x = kta1, file = paste0(path,"PLOT/",mark[i],"/kegg.upup.xlsx"))
    openxlsx::write.xlsx(x = kta2, file = paste0(path,"PLOT/",mark[i],"/kegg.down.xlsx"))
    #
    kta1$Description <- factor(kta1$Description, levels = rev(kta1$Description))
    p_05 <- ggplot(kta1[1:10,])+
      geom_bar(aes(x=-log10(pvalue), y=Description), stat = "identity", fill="#1ba784", alpha=1.0) +
      theme_few() +
      #scale_x_continuous(breaks = seq(8)) +
      #geom_vline(aes(xintercept = 1), color = "white", linewidth  = 0.8) +
      #geom_vline(aes(xintercept = 2), color = "white", linewidth  = 0.8) +
      #geom_vline(aes(xintercept = 3), color = "white", linewidth  = 0.8) +
      #geom_vline(aes(xintercept = 4), color = "white", linewidth  = 0.8) +
      labs(y="",title = "this is title--upup") +
      theme(text = element_text(size = 12,family = "serif",color = "black"),
            panel.border = element_rect(linewidth = 1),
            plot.title = element_text(size = 12,family = "serif",color = "#1ba784",face = "bold"),
            axis.title = element_text(size = 12,family = "serif",color = "#1ba784"),
            axis.text  = element_text(size = 12,family = "serif",color = "#1ba784"))
    p_05
    #
    kta2$Description <- factor(kta2$Description, levels = rev(kta2$Description))
    p_06 <- ggplot(kta2[1:10,])+
      geom_bar(aes(x=-log10(pvalue), y=Description), stat = "identity", fill="#5861AC", alpha=1.0) +
      theme_few() +
      #scale_x_continuous(breaks = seq(10)) +
      #geom_vline(aes(xintercept = 1), color = "white", linewidth  = 0.8) +
      #geom_vline(aes(xintercept = 2), color = "white", linewidth  = 0.8) +
      #geom_vline(aes(xintercept = 3), color = "white", linewidth  = 0.8) +
      #geom_vline(aes(xintercept = 4), color = "white", linewidth  = 0.8) +
      #geom_vline(aes(xintercept = 5), color = "white", linewidth  = 0.8) +
      labs(y="",title = "this is title--down") +
      theme(text = element_text(size = 12,family = "serif",color = "black"),
            panel.border = element_rect(linewidth = 1),
            plot.title = element_text(size = 12,family = "serif",color = "#5861AC",face = "bold"),
            axis.title = element_text(size = 12,family = "serif",color = "#5861AC"),
            axis.text  = element_text(size = 12,family = "serif",color = "#5861AC"))
    p_06
    plot <- p_05 +p_06 +plot_layout(nrow = 2)
    ggsave(plot = plot,
           units = "cm", width = 22, height = 15,
           filename = paste0(path,"PLOT/",mark[i],"/KEGG.pdf"))
  }
}
################################################################################
#
#
#pheatmap
#
#
################################################################################
#26、for add.new.1的2 vs 1, 3 vs 2, 以及add.new.2的2 vs 1, 3 vs 2
#[在ofAddNew.1上调 + 在ofAddNew.2的下调]
{
  sets <- list(c("group2.E14.P0.4500.Pyru","group4.ghyuu.P0.4500.Pyru","group4.ghyuu.P0.0180.Pyru"),
               c("group3.E14.P0.4500","group5.ghyuu.P0.4500","group5.ghyuu.P0.0180"))
  gencpms <- openxlsx::read.xlsx(paste0(path,"countUniq/TEtools/TEtools_huhuhu.cpm.onlyGene.xlsx"))
  gencpm1 <- gencpms[,c(grep(sets[[1]][1],colnames(gencpms)),
                        grep(sets[[1]][2],colnames(gencpms)),
                        grep(sets[[1]][3],colnames(gencpms)))]
  rownames(gencpm1) <- gencpms$theName
  gencpm2 <- gencpms[,c(grep(sets[[2]][1],colnames(gencpms)),
                        grep(sets[[2]][2],colnames(gencpms)),
                        grep(sets[[2]][3],colnames(gencpms)))]
  rownames(gencpm2) <- gencpms$theName
  #
  mark <- paste0("ofAddNew.",seq(4))
  ff01 <- paste0(path, "DESeq2Uniq/TEtools/",mark,".geneOnly.csv")
  resdata <- as.data.frame(read.csv(ff01[1], header = T, stringsAsFactors = F))
  upup_01 <- resdata[which(resdata$log2FoldChange > log2(1.5)  &  resdata$padj < 0.05),"X"]
  resdata <- as.data.frame(read.csv(ff01[2], header = T, stringsAsFactors = F))
  down_02 <- resdata[which(resdata$log2FoldChange < -log2(1.5) &  resdata$padj < 0.05),"X"]
  toview1 <- unique(c(upup_01, down_02))
  toview1 <- gencpm1[which(rownames(gencpm1) %in% toview1),]
  resdata <- as.data.frame(read.csv(ff01[3], header = T, stringsAsFactors = F))
  down_03 <- resdata[which(resdata$log2FoldChange < -log2(1.5) &  resdata$padj < 0.05),"X"]
  resdata <- as.data.frame(read.csv(ff01[4], header = T, stringsAsFactors = F))
  upup_04 <- resdata[which(resdata$log2FoldChange > log2(1.5)  &  resdata$padj < 0.05),"X"]
  toview2 <- unique(c(down_03, upup_04))
  toview2 <- gencpm2[which(rownames(gencpm2) %in% toview2),]
  #
  pheatmap(mat = toview1,
           main="\nheatmap of genes: add.new.1",
           scale = "row", cellwidth = 20, 
           show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = F,
           labels_col = colnames(toview1), angle_col = "45",
           breaks = seq(-2,2,by=0.01), border = F,
           color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                      colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
  pheatmap(mat = toview2,
           main="\nheatmap of genes: add.new.2",
           scale = "row", cellwidth = 20, 
           show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = F,
           labels_col = colnames(toview2), angle_col = "45",
           breaks = seq(-2,2,by=0.01), border = F,
           color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                      colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
  #
  #
  #
  testss <- toview2[which(rownames(toview2) %in% down_03),]
  pheatmap(mat = testss,
           main="\nheatmap of genes: add.new.2",
           scale = "row", cellwidth = 20, 
           show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = F,
           labels_col = colnames(testss), angle_col = "45",
           breaks = seq(-2,2,by=0.01), border = F,
           color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                      colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
}
################################################################################
#27、for add.new.1的2 vs 1, 3 vs 1, 以及add.new.2的2 vs 1, 3 vs 1
#[在ofAddNex1上调 + 在ofAddNex2的下调]
{
  sets <- list(c("group2.E14.P0.4500.Pyru","group4.ghyuu.P0.4500.Pyru","group4.ghyuu.P0.0180.Pyru"),
               c("group3.E14.P0.4500","group5.ghyuu.P0.4500","group5.ghyuu.P0.0180"))
  mark <- paste0("ofAddNex.",seq(4))
  ff01 <- paste0(path, "DESeq2Uniq/TEtools/",mark,".geneOnly.csv")
  #
  res1 <- as.data.frame(read.csv(ff01[1], header = T, stringsAsFactors = F))
  upup_01 <- res1[which(res1$log2FoldChange > log2(1.5) & res1$padj < 0.05), "X"]
  res2 <- as.data.frame(read.csv(ff01[2], header = T, stringsAsFactors = F))
  down_02 <- res2[which(res2$log2FoldChange < -log2(1.5) & res2$padj < 0.05),"X"]
  res3 <- as.data.frame(read.csv(ff01[3], header = T, stringsAsFactors = F))
  down_03 <- res3[which(res3$log2FoldChange < -log2(1.5) & res3$padj < 0.05),"X"]
  res4 <- as.data.frame(read.csv(ff01[4], header = T, stringsAsFactors = F))
  upup_04 <- res4[which(res4$log2FoldChange > log2(1.5) & res4$padj < 0.05), "X"]
  #
  toview1 <- left_join(res1[which(res1$X %in% upup_01),c(1,3)],res2[which(res2$X %in% upup_01),c(1,3)],by="X")
  rownames(toview1) <- toview1$X; toview1 <- toview1[,-1]
  pheatmap(mat = toview1,
           main="\nheatmap of genes: add.new.1",
           scale = "none", cellwidth = 20, 
           show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = F,
           labels_col = colnames(toview1), angle_col = "45",
           breaks = seq(-2,2,by=0.01), border = F,
           color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                      colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
  #
  toview2 <- left_join(res3[which(res3$X %in% down_03),c(1,3)],res4[which(res4$X %in% down_03),c(1,3)],by="X")
  rownames(toview2) <- toview2$X; toview2 <- toview2[,-1]
  pheatmap(mat = toview2,
           main="\nheatmap of genes: add.new.2",
           scale = "none", cellwidth = 20, 
           show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = F,
           labels_col = colnames(toview2), angle_col = "45",
           breaks = seq(-2,2,by=0.01), border = F,
           color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                      colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
}
################################################################################
#28、“180里，2'C gene没有下调吗？” [说的是add.new.2的3 vs 1的GSEA]
{
  mark <- paste0("ofAddNex.",seq(4))
  ff01 <- paste0(path, "DESeq2/",mark,".geneOnly.csv")
  for (i in seq_along(ff01)) {
    temp <- read.csv(ff01[i], header = T)
    temp <- temp[, c("X", "log2FoldChange")]
    temp <- temp[order(temp$log2FoldChange, decreasing = T),]
    write.table(x = temp, 
                file = paste0(path, "DESeq2/",mark[i],".res.rnk"),
                quote = F, sep = "\t", col.names = F, row.names = F)
    cmd_12 <- paste0("gseapy prerank -g ","/Reference/aaSHY/GSEAgmt/TwoGenes.gmt -r ",
                     paste0(path, "DESeq2/",mark[i],".res.rnk"), " -p 10 -o ", 
                     paste0(path, "PLOT/",mark[i],"/enrichGSEA/"))#conda activate base gseapy 0.10.8
    print(cmd_12)
  }
}
{
  mark <- paste0("ofAddNew.",seq(4))
  ff01 <- paste0(path, "DESeq2/",mark,".geneOnly.csv")
  aa <- read.csv(file = ff01[4], header = T, stringsAsFactors = F)[,c(1,3)]
  ab <- suppressWarnings(clusterProfiler::bitr(geneID = aa$X, 
                                               fromType = "SYMBOL",
                                               toType = "ENTREZID",
                                               OrgDb = "org.Mm.eg.db", drop = T))
  ac <- dplyr::distinct(.data = ab, SYMBOL, .keep_all = T)
  ac$rank <- apply(ac, 1, function(x){aa[which(aa$X==x[1]),2]})
  ad <- ac[order(ac$rank, decreasing = T),c(2,3)]
  ae <- ad$rank; names(ae) <- ad$ENTREZID
  rm(aa,ab,ac,ad)
  #
  aa <- clusterProfiler::read.gmt("/Reference/aaSHY/GSEAgmt/TwoGenes.gmt")
  aa <- suppressWarnings(clusterProfiler::bitr(aa$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db"))
  ab <- data.frame(ont = rep("act_2cell_gene", length(unique(aa$ENTREZID))), gen = unique(aa$ENTREZID))
  af_gsea <- suppressWarnings(clusterProfiler::GSEA(geneList = ae, maxGSSize = 1000, TERM2GENE = ab, pvalueCutoff = 1))
  rm(aa,ab)
  p_a <- enrichplot::gseaplot2(x = af_gsea, geneSetID = 1, title = "GSEA for 2-cell Genes", pvalue_table = T)
  print(p_a)
}
################################################################################
#29、for add.new.1的2 vs 1, 3 vs 1, 以及add.new.2的2 vs 1, 3 vs 1 【对照组CPM设为1】
{
  sets <- list(c("group2.E14.P0.4500.Pyru","group4.ghyuu.P0.4500.Pyru","group4.ghyuu.P0.0180.Pyru"),
               c("group3.E14.P0.4500","group5.ghyuu.P0.4500","group5.ghyuu.P0.0180"))
  gencpms <- openxlsx::read.xlsx(paste0(path,"countUniq/TEtools/TEtools_huhuhu.cpm.onlyGene.xlsx"))
  gencpm1 <- gencpms[,c(grep(sets[[1]][1],colnames(gencpms)),
                        grep(sets[[1]][2],colnames(gencpms)),
                        grep(sets[[1]][3],colnames(gencpms)))]
  rownames(gencpm1) <- gencpms$theName
  gencpm2 <- gencpms[,c(grep(sets[[2]][1],colnames(gencpms)),
                        grep(sets[[2]][2],colnames(gencpms)),
                        grep(sets[[2]][3],colnames(gencpms)))]
  rownames(gencpm2) <- gencpms$theName
  #
  mark <- paste0("ofAddNew.",seq(4))
  ff01 <- paste0(path, "DESeq2Uniq/TEtools/",mark,".geneOnly.csv")
  resdata <- as.data.frame(read.csv(ff01[1], header = T, stringsAsFactors = F))
  upup_01 <- resdata[which(resdata$log2FoldChange > log2(1.5)  &  resdata$padj < 0.05),"X"]
  resdata <- as.data.frame(read.csv(ff01[2], header = T, stringsAsFactors = F))
  down_02 <- resdata[which(resdata$log2FoldChange < -log2(1.5) &  resdata$padj < 0.05),"X"]
  toview1 <- unique(c(upup_01, down_02))
  toview1 <- gencpm1[which(rownames(gencpm1) %in% toview1),]
  toview1 <- toview1 +0.001
  toview1 <- as.data.frame(apply(toview1,2,function(x){x/rowMeans(toview1[,c(1,2)])}))
  toview1 <- log2(toview1 +0.001)
  toview1 <- toview1[which(rowSums(is.na(toview1)) <1),]
  #
  resdata <- as.data.frame(read.csv(ff01[3], header = T, stringsAsFactors = F))
  down_03 <- resdata[which(resdata$log2FoldChange < -log2(1.5) &  resdata$padj < 0.05),"X"]
  resdata <- as.data.frame(read.csv(ff01[4], header = T, stringsAsFactors = F))
  upup_04 <- resdata[which(resdata$log2FoldChange > log2(1.5)  &  resdata$padj < 0.05),"X"]
  toview2 <- unique(c(down_03, upup_04))
  toview2 <- gencpm2[which(rownames(gencpm2) %in% toview2),]
  toview2 <- toview2 +0.001
  toview2 <- as.data.frame(apply(toview2,2,function(x){x/rowMeans(toview2[,c(1,2)])}))
  toview2 <- log2(toview2 +0.001)
  toview2 <- toview2[which(rowSums(is.na(toview2)) <1),]
  #
  pdf(file = paste0(path,"PLOT_uniq/heatmap.add.new.pdf"), height = 6, width = 6, family = "serif")
  pheatmap(mat = toview1,
           main="\nheatmap of genes: add.new.1",
           scale = "none", cellwidth = 20, 
           show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = F,
           labels_col = colnames(toview1), angle_col = "45",
           breaks = seq(-2,2,by=0.01), border = F,
           color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                      colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
  pheatmap(mat = toview2,
           main="\nheatmap of genes: add.new.2",
           scale = "none", cellwidth = 20, 
           show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = F,
           labels_col = colnames(toview2), angle_col = "45",
           breaks = seq(-2,2,by=0.01), border = F,
           color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                      colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
  dev.off()
}
################################################################################
#30、for add.new.1的2 vs 1, 3 vs 1, 以及add.new.2的2 vs 1, 3 vs 1 【对照组CPM设为1】【TE-无基因背景】
{
  sets <- list(c("group2.E14.P0.4500.Pyru","group4.ghyuu.P0.4500.Pyru","group4.ghyuu.P0.0180.Pyru"),
               c("group3.E14.P0.4500","group5.ghyuu.P0.4500","group5.ghyuu.P0.0180"))
  gencpms <- openxlsx::read.xlsx(paste0(path,"countUniq/TEtools/TEtools_huhuhu.cpm.onlyTE.xlsx"))
  gencpm1 <- gencpms[,c(grep(sets[[1]][1],colnames(gencpms)),
                        grep(sets[[1]][2],colnames(gencpms)),
                        grep(sets[[1]][3],colnames(gencpms)))]
  rownames(gencpm1) <- gencpms$theName
  gencpm2 <- gencpms[,c(grep(sets[[2]][1],colnames(gencpms)),
                        grep(sets[[2]][2],colnames(gencpms)),
                        grep(sets[[2]][3],colnames(gencpms)))]
  rownames(gencpm2) <- gencpms$theName
  #
  mark <- paste0("ofAddNew.",seq(4))
  ff01 <- paste0(path, "DESeq2Uniq/TEtools/",mark,".teOnly.csv")
  resdata <- as.data.frame(read.csv(ff01[1], header = T, stringsAsFactors = F))
  upup_01 <- resdata[which(resdata$log2FoldChange > log2(1.5)  &  resdata$padj < 0.05),"X"]
  resdata <- as.data.frame(read.csv(ff01[2], header = T, stringsAsFactors = F))
  down_02 <- resdata[which(resdata$log2FoldChange < -log2(1.5) &  resdata$padj < 0.05),"X"]
  toview1 <- unique(c(upup_01, down_02))
  toview1 <- gencpm1[which(rownames(gencpm1) %in% toview1),]
  toview1 <- toview1 +0.001
  toview1 <- as.data.frame(apply(toview1,2,function(x){x/rowMeans(toview1[,c(1,2)])}))
  toview1 <- log2(toview1 +0.001)
  toview1 <- toview1[which(rowSums(is.na(toview1)) <1),]
  #
  resdata <- as.data.frame(read.csv(ff01[3], header = T, stringsAsFactors = F))
  down_03 <- resdata[which(resdata$log2FoldChange < -log2(1.5) &  resdata$padj < 0.05),"X"]
  resdata <- as.data.frame(read.csv(ff01[4], header = T, stringsAsFactors = F))
  upup_04 <- resdata[which(resdata$log2FoldChange > log2(1.5)  &  resdata$padj < 0.05),"X"]
  toview2 <- unique(c(down_03, upup_04))
  toview2 <- gencpm2[which(rownames(gencpm2) %in% toview2),]
  toview2 <- toview2 +0.001
  toview2 <- as.data.frame(apply(toview2,2,function(x){x/rowMeans(toview2[,c(1,2)])}))
  toview2 <- log2(toview2 +0.001)
  toview2 <- toview2[which(rowSums(is.na(toview2)) <1),]
  #
  pdf(file = paste0(path,"PLOT_uniq/heatmap.add.new.TE.Only.pdf"), height = 6, width = 6, family = "serif")
  pheatmap(mat = toview1,
           main="\nheatmap of genes: add.new.1",
           scale = "none", cellwidth = 20, 
           show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = F,
           labels_col = colnames(toview1), angle_col = "45",
           breaks = seq(-2,2,by=0.01), border = F,
           color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                      colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
  pheatmap(mat = toview2,
           main="\nheatmap of genes: add.new.2",
           scale = "none", cellwidth = 20, 
           show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = F,
           labels_col = colnames(toview2), angle_col = "45",
           breaks = seq(-2,2,by=0.01), border = F,
           color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                      colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
  dev.off()
}
################################################################################
#31、for add.new.1的2 vs 1, 3 vs 1, 以及add.new.2的2 vs 1, 3 vs 1 【对照组CPM设为1】【TE-有基因背景】
{
  sets <- list(c("group2.E14.P0.4500.Pyru","group4.ghyuu.P0.4500.Pyru","group4.ghyuu.P0.0180.Pyru"),
               c("group3.E14.P0.4500","group5.ghyuu.P0.4500","group5.ghyuu.P0.0180"))
  gencpms <- openxlsx::read.xlsx(paste0(path,"countUniq/TEtools/TEtools_huhuhu.cpm.allGeneTE.xlsx"))
  gencpm1 <- gencpms[,c(grep(sets[[1]][1],colnames(gencpms)),
                        grep(sets[[1]][2],colnames(gencpms)),
                        grep(sets[[1]][3],colnames(gencpms)))]
  rownames(gencpm1) <- gencpms$theName
  gencpm2 <- gencpms[,c(grep(sets[[2]][1],colnames(gencpms)),
                        grep(sets[[2]][2],colnames(gencpms)),
                        grep(sets[[2]][3],colnames(gencpms)))]
  rownames(gencpm2) <- gencpms$theName
  #
  mark <- paste0("ofAddNew.",seq(4))
  ff01 <- paste0(path, "DESeq2Uniq/TEtools/",mark,".te.GeneBackGround.csv")
  resdata <- as.data.frame(read.csv(ff01[1], header = T, stringsAsFactors = F))
  upup_01 <- resdata[which(resdata$log2FoldChange > log2(1.5)  &  resdata$padj < 0.05),"X"]
  resdata <- as.data.frame(read.csv(ff01[2], header = T, stringsAsFactors = F))
  down_02 <- resdata[which(resdata$log2FoldChange < -log2(1.5) &  resdata$padj < 0.05),"X"]
  toview1 <- unique(c(upup_01, down_02))
  toview1 <- gencpm1[which(rownames(gencpm1) %in% toview1),]
  toview1 <- toview1 +0.001
  toview1 <- as.data.frame(apply(toview1,2,function(x){x/rowMeans(toview1[,c(1,2)])}))
  toview1 <- log2(toview1 +0.001)
  toview1 <- toview1[which(rowSums(is.na(toview1)) <1),]
  #
  resdata <- as.data.frame(read.csv(ff01[3], header = T, stringsAsFactors = F))
  down_03 <- resdata[which(resdata$log2FoldChange < -log2(1.5) &  resdata$padj < 0.05),"X"]
  resdata <- as.data.frame(read.csv(ff01[4], header = T, stringsAsFactors = F))
  upup_04 <- resdata[which(resdata$log2FoldChange > log2(1.5)  &  resdata$padj < 0.05),"X"]
  toview2 <- unique(c(down_03, upup_04))
  toview2 <- gencpm2[which(rownames(gencpm2) %in% toview2),]
  toview2 <- toview2 +0.001
  toview2 <- as.data.frame(apply(toview2,2,function(x){x/rowMeans(toview2[,c(1,2)])}))
  toview2 <- log2(toview2 +0.001)
  toview2 <- toview2[which(rowSums(is.na(toview2)) <1),]
  #
  pdf(file = paste0(path,"PLOT_uniq/heatmap.add.new.TE.GeneBackGround.pdf"), height = 6, width = 6, family = "serif")
  pheatmap(mat = toview1,
           main="\nheatmap of genes: add.new.1",
           scale = "none", cellwidth = 20, 
           show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = F,
           labels_col = colnames(toview1), angle_col = "45",
           breaks = seq(-2,2,by=0.01), border = F,
           color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                      colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
  pheatmap(mat = toview2,
           main="\nheatmap of genes: add.new.2",
           scale = "none", cellwidth = 20, 
           show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = F,
           labels_col = colnames(toview2), angle_col = "45",
           breaks = seq(-2,2,by=0.01), border = F,
           color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                      colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
  dev.off()
}
################################################################################
#
#
#
#multimap: 选所有“2 vs 1” dysregulated gene, 而非“2 vs 1 + 3 vs 2 或+ 3 vs 1”
#
#
#
################################################################################
#32、gene Only Pheatmap
{
  sets <- list(c("group2.E14.P0.4500.Pyru","group4.ghyuu.P0.4500.Pyru","group4.ghyuu.P0.0180.Pyru"),
               c("group3.E14.P0.4500","group5.ghyuu.P0.4500","group5.ghyuu.P0.0180"))
  gencpms <- openxlsx::read.xlsx(paste0(path,"count/huhuhu.cpm.onlyGene.rename.xlsx"))
  gencpm1 <- gencpms[,c(grep(sets[[1]][1],colnames(gencpms)),
                        grep(sets[[1]][2],colnames(gencpms)),
                        grep(sets[[1]][3],colnames(gencpms)))]
  rownames(gencpm1) <- gencpms$theName
  gencpm2 <- gencpms[,c(grep(sets[[2]][1],colnames(gencpms)),
                        grep(sets[[2]][2],colnames(gencpms)),
                        grep(sets[[2]][3],colnames(gencpms)))]
  rownames(gencpm2) <- gencpms$theName
  #
  mark <- paste0("ofAddNew.",seq(4))
  ff01 <- paste0(path, "DESeq2/",mark,".geneOnly.csv")
  resdata <- as.data.frame(read.csv(ff01[1], header = T, stringsAsFactors = F))
  upup_01 <- resdata[which(resdata$log2FoldChange > log2(1.5)  &  resdata$padj < 0.05),"X"]
  resdata <- as.data.frame(read.csv(ff01[1], header = T, stringsAsFactors = F))
  down_02 <- resdata[which(resdata$log2FoldChange < -log2(1.5) &  resdata$padj < 0.05),"X"]
  toview1 <- unique(c(upup_01, down_02))
  toview1 <- gencpm1[which(rownames(gencpm1) %in% toview1),]
  toview1 <- toview1 +0.001
  toview1 <- as.data.frame(apply(toview1,2,function(x){x/rowMeans(toview1[,c(1,2)])}))
  toview1 <- log2(toview1 +0.001)
  toview1 <- toview1[which(rowSums(is.na(toview1)) <1),]
  #
  resdata <- as.data.frame(read.csv(ff01[3], header = T, stringsAsFactors = F))
  down_03 <- resdata[which(resdata$log2FoldChange < -log2(1.5) &  resdata$padj < 0.05),"X"]
  resdata <- as.data.frame(read.csv(ff01[3], header = T, stringsAsFactors = F))
  upup_04 <- resdata[which(resdata$log2FoldChange > log2(1.5)  &  resdata$padj < 0.05),"X"]
  toview2 <- unique(c(down_03, upup_04))
  toview2 <- gencpm2[which(rownames(gencpm2) %in% toview2),]
  toview2 <- toview2 +0.001
  toview2 <- as.data.frame(apply(toview2,2,function(x){x/rowMeans(toview2[,c(1,2)])}))
  toview2 <- log2(toview2 +0.001)
  toview2 <- toview2[which(rowSums(is.na(toview2)) <1),]
  #
  pdf(file = paste0(path,"PLOT/heatmap.add.new.gene.Only.pdf"), height = 6, width = 6, family = "serif")
  pheatmap(mat = toview1,
           main="\nheatmap of genes: add.new.1",
           scale = "none", cellwidth = 20, 
           show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = F,
           labels_col = colnames(toview1), angle_col = "45",
           breaks = seq(-2,2,by=0.01), border = F,
           color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                      colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
  pheatmap(mat = toview2,
           main="\nheatmap of genes: add.new.2",
           scale = "none", cellwidth = 20, 
           show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = F,
           labels_col = colnames(toview2), angle_col = "45",
           breaks = seq(-2,2,by=0.01), border = F,
           color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                      colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
  dev.off()
}
################################################################################
#33、跑multimap模式下的add.new的TE的差异结果
#DESeq2
{
  shy_diff <- function(mark,grou,geneFile,siteFile,clas,vsvs){
    tmp1 <- read.table(geneFile, header = T, check.names = F, stringsAsFactors = F, row.names = 1)
    tmp1 <- cbind(tmp1[,grep(sets[[grou]][1],colnames(tmp1))],
                  tmp1[,grep(sets[[grou]][2],colnames(tmp1))])
    dfclas <- data.frame(row.names = colnames(tmp1), clas, stringsAsFactors = T)
    colLen <- ceiling(length(colnames(tmp1))/2)
    geneAA <- tmp1[grep(rownames(tmp1), pattern = ":", invert = T),]
    teOnly <- tmp1[grep(rownames(tmp1), pattern = ":", invert = F),]
    teGene <- tmp1
    #
    geneAA <- DESeq2::DESeqDataSetFromMatrix(countData = geneAA, colData = dfclas, design = ~clas)
    geneAA <- geneAA[rowSums(BiocGenerics::counts(geneAA) > 2) >= colLen, ]#[rowMeans(BiocGenerics::counts(geneAA)) > 1,]
    geneAA <- DESeq2::DESeq(geneAA,)
    geneAA <- DESeq2::results(geneAA, contrast = vsvs)
    write.csv(geneAA, paste0(path,"DESeq2/",mark,".",grou,".geneOnly.csv"))
    upup <- subset(geneAA, log2FoldChange > log2(1.5) & padj < 0.05)
    down <- subset(geneAA, log2FoldChange < -log2(1.5)& padj < 0.05)
    write.csv(upup,paste0(path,"DESeq2/",mark,".",grou,".geneOnly.upup.csv"))
    write.csv(down,paste0(path,"DESeq2/",mark,".",grou,".geneOnly.down.csv"))
    #
    teOnly <- DESeq2::DESeqDataSetFromMatrix(countData = teOnly, colData = dfclas, design = ~clas)
    teOnly <- teOnly[rowSums(BiocGenerics::counts(teOnly) > 2) >= colLen, ]#[rowMeans(BiocGenerics::counts(teOnly)) > 1,]
    teOnly <- DESeq2::DESeq(teOnly,)
    teOnly <- DESeq2::results(teOnly, contrast = vsvs)
    write.csv(teOnly, paste0(path,"DESeq2/",mark,".",grou,".teOnly.csv"))
    upup <- subset(teOnly, log2FoldChange > log2(1.5) & padj < 0.05)
    down <- subset(teOnly, log2FoldChange < -log2(1.5)& padj < 0.05)
    write.csv(upup,paste0(path,"DESeq2/",mark,".",grou,".teOnly.upup.csv"))
    write.csv(down,paste0(path,"DESeq2/",mark,".",grou,".teOnly.down.csv"))
    #
    teGene <- DESeq2::DESeqDataSetFromMatrix(countData = teGene, colData = dfclas, design = ~clas)
    teGene <- teGene[rowSums(BiocGenerics::counts(teGene) > 2) >= colLen, ]#[rowMeans(BiocGenerics::counts(teGene)) > 1,]
    teGene <- DESeq2::DESeq(teGene,)
    teGene <- DESeq2::results(teGene, contrast = vsvs)
    teGene <- teGene[grep(rownames(teGene), pattern = ":", invert = F),]
    write.csv(teGene, paste0(path,"DESeq2/",mark,".",grou,".te.GeneBackGround.csv"))
    upup <- subset(teGene, log2FoldChange > log2(1.5) & padj < 0.05)
    down <- subset(teGene, log2FoldChange < -log2(1.5)& padj < 0.05)
    write.csv(upup,paste0(path,"DESeq2/",mark,".",grou,".te.GeneBackGround.upup.csv"))
    write.csv(down,paste0(path,"DESeq2/",mark,".",grou,".te.GeneBackGround.down.csv"))
    #
    #
    tmp2 <- read.table(siteFile, header = T, check.names = F, stringsAsFactors = F, row.names = 1)
    tmp2 <- cbind(tmp2[,grep(sets[[grou]][1],colnames(tmp2))],
                  tmp2[,grep(sets[[grou]][2],colnames(tmp2))])
    colLen <- ceiling(length(colnames(tmp2))/2)
    #
    teGene <- tmp2
    teGene <- DESeq2::DESeqDataSetFromMatrix(countData = teGene, colData = dfclas, design = ~clas)
    teGene <- teGene[rowSums(BiocGenerics::counts(teGene) > 2) >= colLen, ]#[rowMeans(BiocGenerics::counts(teGene)) > 1,]
    teGene <- DESeq2::DESeq(teGene,)
    teGene <- DESeq2::results(teGene, contrast = vsvs)
    teGene <- teGene[grep(rownames(teGene), pattern = ":", invert = F),]
    write.csv(teGene, paste0(path,"DESeq2/",mark,".",grou,".te.site.GeneBackGround.csv"))
    upup <- subset(teGene, log2FoldChange > log2(1.5) & padj < 0.05)
    down <- subset(teGene, log2FoldChange < -log2(1.5)& padj < 0.05)
    write.csv(upup, paste0(path,"DESeq2/",mark,".",grou,".te.site.GeneBackGround.upup.csv"))
    write.csv(down, paste0(path,"DESeq2/",mark,".",grou,".te.site.GeneBackGround.down.csv"))
    #
    teGene <- tmp2[grep(rownames(tmp2), pattern = ":", invert = F),]
    teGene <- DESeq2::DESeqDataSetFromMatrix(countData = teGene, colData = dfclas, design = ~clas)
    teGene <- teGene[rowSums(BiocGenerics::counts(teGene) > 2) >= colLen, ]#[rowMeans(BiocGenerics::counts(teGene)) > 1,]
    teGene <- DESeq2::DESeq(teGene,)
    teGene <- DESeq2::results(teGene, contrast = vsvs)
    write.csv(teGene, paste0(path,"DESeq2/",mark,".",grou,".te.site.teOnly.csv"))
    upup <- subset(teGene, log2FoldChange > log2(1.5) & padj < 0.05)
    down <- subset(teGene, log2FoldChange < -log2(1.5)& padj < 0.05)
    write.csv(upup, paste0(path,"DESeq2/",mark,".",grou,".te.site.teOnly.upup.csv"))
    write.csv(down, paste0(path,"DESeq2/",mark,".",grou,".te.site.teOnly.down.csv"))
  }
  sets <- list(c("group4.ghyuu.P0.4500.Pyru","group2.E14.P0.4500.Pyru"),
               c("group4.ghyuu.P0.0180.Pyru","group4.ghyuu.P0.4500.Pyru"),
               c("group5.ghyuu.P0.4500","group3.E14.P0.4500"),
               c("group5.ghyuu.P0.0180","group5.ghyuu.P0.4500"))
  for (grou in seq_along(sets)) {
    print(grou)
    shy_diff(mark = "ofAddNew", grou,
             geneFile = paste0(path, "count/huhuhu.rename.cntTable"),
             siteFile = paste0(path, "count/huhuhu.site.rename.cntTable"),
             clas = factor(c("OP","OP","WT","WT")), vsvs = c("clas","OP","WT"))
  }
}
################################################################################
#33、te Pheatmap [无基因背景]
{
  sets <- list(c("group2.E14.P0.4500.Pyru","group4.ghyuu.P0.4500.Pyru","group4.ghyuu.P0.0180.Pyru"),
               c("group3.E14.P0.4500","group5.ghyuu.P0.4500","group5.ghyuu.P0.0180"))
  gencpms <- openxlsx::read.xlsx(paste0(path,"count/huhuhu.cpm.onlyTE.rename.xlsx"))
  gencpm1 <- gencpms[,c(grep(sets[[1]][1],colnames(gencpms)),
                        grep(sets[[1]][2],colnames(gencpms)),
                        grep(sets[[1]][3],colnames(gencpms)))]
  rownames(gencpm1) <- gencpms$theName
  gencpm2 <- gencpms[,c(grep(sets[[2]][1],colnames(gencpms)),
                        grep(sets[[2]][2],colnames(gencpms)),
                        grep(sets[[2]][3],colnames(gencpms)))]
  rownames(gencpm2) <- gencpms$theName
  #
  mark <- paste0("ofAddNew.",seq(4))
  ff01 <- paste0(path, "DESeq2/",mark,".teOnly.csv")
  resdata <- as.data.frame(read.csv(ff01[1], header = T, stringsAsFactors = F))
  upup_01 <- resdata[which(resdata$log2FoldChange > log2(1.5)  &  resdata$padj < 0.05),"X"]
  resdata <- as.data.frame(read.csv(ff01[1], header = T, stringsAsFactors = F))
  down_02 <- resdata[which(resdata$log2FoldChange < -log2(1.5) &  resdata$padj < 0.05),"X"]
  toview1 <- unique(c(upup_01, down_02))
  toview1 <- gencpm1[which(rownames(gencpm1) %in% toview1),]
  toview1 <- toview1 +0.001
  toview1 <- as.data.frame(apply(toview1,2,function(x){x/rowMeans(toview1[,c(1,2)])}))
  toview1 <- log2(toview1 +0.001)
  toview1 <- toview1[which(rowSums(is.na(toview1)) <1),]
  #
  resdata <- as.data.frame(read.csv(ff01[3], header = T, stringsAsFactors = F))
  down_03 <- resdata[which(resdata$log2FoldChange < -log2(1.5) &  resdata$padj < 0.05),"X"]
  resdata <- as.data.frame(read.csv(ff01[3], header = T, stringsAsFactors = F))
  upup_04 <- resdata[which(resdata$log2FoldChange > log2(1.5)  &  resdata$padj < 0.05),"X"]
  toview2 <- unique(c(down_03, upup_04))
  toview2 <- gencpm2[which(rownames(gencpm2) %in% toview2),]
  toview2 <- toview2 +0.001
  toview2 <- as.data.frame(apply(toview2,2,function(x){x/rowMeans(toview2[,c(1,2)])}))
  toview2 <- log2(toview2 +0.001)
  toview2 <- toview2[which(rowSums(is.na(toview2)) <1),]
  #
  pdf(file = paste0(path,"PLOT/heatmap.add.new.TE.Only.pdf"), height = 6, width = 6, family = "serif")
  pheatmap(mat = toview1,
           main="\nheatmap of TEs: add.new.1",
           scale = "none", cellwidth = 20, 
           show_rownames = T, show_colnames = T, cluster_row = T, cluster_cols = F,
           labels_col = colnames(toview1), angle_col = "45",
           breaks = seq(-4,4,by=0.01), border = F,
           color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-4,4,by=0.01))/2),
                      colorRampPalette(colors = c("white","red"))(length(seq(-4,4,by=0.01))/2)))
  pheatmap(mat = toview2,
           main="\nheatmap of TEs: add.new.2",
           scale = "none", cellwidth = 20, 
           show_rownames = T, show_colnames = T, cluster_row = T, cluster_cols = F,
           labels_col = colnames(toview2), angle_col = "45",
           breaks = seq(-4,4,by=0.01), border = F,
           color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-4,4,by=0.01))/2),
                      colorRampPalette(colors = c("white","red"))(length(seq(-4,4,by=0.01))/2)))
  dev.off()
}
################################################################################
#34、te Pheatmap [有基因背景] [采用] [OK]------------------------------------------------------------------
{
  sets <- list(c("group2.E14.P0.4500.Pyru","group4.ghyuu.P0.4500.Pyru","group4.ghyuu.P0.0180.Pyru"),
               c("group3.E14.P0.4500","group5.ghyuu.P0.4500","group5.ghyuu.P0.0180"))
  gencpms <- openxlsx::read.xlsx(paste0(path,"count/huhuhu.cpm.allGeneTE.rename.xlsx"))
  gencpm1 <- gencpms[,c(grep(sets[[1]][1],colnames(gencpms)),
                        grep(sets[[1]][2],colnames(gencpms)),
                        grep(sets[[1]][3],colnames(gencpms)))]
  rownames(gencpm1) <- gencpms$theName
  gencpm2 <- gencpms[,c(grep(sets[[2]][1],colnames(gencpms)),
                        grep(sets[[2]][2],colnames(gencpms)),
                        grep(sets[[2]][3],colnames(gencpms)))]
  rownames(gencpm2) <- gencpms$theName
  #
  mark <- paste0("ofAddNew.",seq(4))
  ff01 <- paste0(path, "DESeq2/",mark,".te.GeneBackGround.csv")
  resdata <- as.data.frame(read.csv(ff01[1], header = T, stringsAsFactors = F))
  upup_01 <- resdata[which(resdata$log2FoldChange > log2(1.5)  &  resdata$padj < 0.05),"X"]
  resdata <- as.data.frame(read.csv(ff01[1], header = T, stringsAsFactors = F))
  down_02 <- resdata[which(resdata$log2FoldChange < -log2(1.5) &  resdata$padj < 0.05),"X"]
  toview1 <- unique(c(upup_01, down_02))
  toview1 <- gencpm1[which(rownames(gencpm1) %in% toview1),]
  toview1 <- toview1 +0.001
  toview1 <- as.data.frame(apply(toview1,2,function(x){x/rowMeans(toview1[,c(1,2)])}))
  toview1 <- log2(toview1 +0.001)
  toview1 <- toview1[which(rowSums(is.na(toview1)) <1),]
  #
  resdata <- as.data.frame(read.csv(ff01[3], header = T, stringsAsFactors = F))
  down_03 <- resdata[which(resdata$log2FoldChange < -log2(1.5) &  resdata$padj < 0.05),"X"]
  resdata <- as.data.frame(read.csv(ff01[3], header = T, stringsAsFactors = F))
  upup_04 <- resdata[which(resdata$log2FoldChange > log2(1.5)  &  resdata$padj < 0.05),"X"]
  toview2 <- unique(c(down_03, upup_04))
  toview2 <- gencpm2[which(rownames(gencpm2) %in% toview2),]
  toview2 <- toview2 +0.001
  toview2 <- as.data.frame(apply(toview2,2,function(x){x/rowMeans(toview2[,c(1,2)])}))
  toview2 <- log2(toview2 +0.001)
  toview2 <- toview2[which(rowSums(is.na(toview2)) <1),]
  #
  pdf(file = paste0(path,"PLOT/heatmap.add.new.TE.GeneBackGround-color3.pdf"), height = 6, width = 6, family = "serif")
  pheatmap(mat = toview1,
           main="\nheatmap of TEs: add.new.1",
           scale = "none", cellwidth = 20, 
           show_rownames = T, show_colnames = T, cluster_row = T, cluster_cols = F,
           labels_col = colnames(toview1), angle_col = "45",
           breaks = seq(-3,3,by=0.01), border = F,
           color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-3,3,by=0.01))/2),
                      colorRampPalette(colors = c("white", "red"))(length(seq(-3,3,by=0.01))/2)))
  pheatmap(mat = toview2,
           main="\nheatmap of TEs: add.new.2",
           scale = "none", cellwidth = 20, 
           show_rownames = T, show_colnames = T, cluster_row = T, cluster_cols = F,
           labels_col = colnames(toview2), angle_col = "45",
           breaks = seq(-3,3,by=0.01), border = F,
           color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-3,3,by=0.01))/2),
                      colorRampPalette(colors = c("white", "red"))(length(seq(-3,3,by=0.01))/2)))
  dev.off()
}
################################################################################
#35、gene Only Pheatmap [把在ghyuu-4500里上调的2细胞基因摘出来画下] [采用] [OK]---------------------------------
{
  #
  twos <- read.table("/Reference/aaSHY/BED/special/TWOC_gene.bed")
  twos <- twos$V4
  #
  sets <- list(c("group2.E14.P0.4500.Pyru","group4.ghyuu.P0.4500.Pyru","group4.ghyuu.P0.0180.Pyru"),
               c("group3.E14.P0.4500","group5.ghyuu.P0.4500","group5.ghyuu.P0.0180"))
  gencpms <- openxlsx::read.xlsx(paste0(path,"count/huhuhu.cpm.onlyGene.rename.xlsx"))
  gencpm1 <- gencpms[,c(grep(sets[[1]][1],colnames(gencpms)),
                        grep(sets[[1]][2],colnames(gencpms)),
                        grep(sets[[1]][3],colnames(gencpms)))]
  rownames(gencpm1) <- gencpms$theName
  gencpm2 <- gencpms[,c(grep(sets[[2]][1],colnames(gencpms)),
                        grep(sets[[2]][2],colnames(gencpms)),
                        grep(sets[[2]][3],colnames(gencpms)))]
  rownames(gencpm2) <- gencpms$theName
  #
  mark <- paste0("ofAddNew.",seq(4))
  ff01 <- paste0(path, "DESeq2/",mark,".geneOnly.csv")
  resdata <- as.data.frame(read.csv(ff01[1], header = T, stringsAsFactors = F))
  upup_01 <- resdata[which(resdata$log2FoldChange > log2(1.5)  &  resdata$padj < 0.05),"X"]
  upup_01 <- intersect(upup_01, twos)
  toview1 <- unique(c(upup_01))
  toview1 <- gencpm1[which(rownames(gencpm1) %in% toview1),]
  toview1 <- toview1 +0.001
  toview1 <- as.data.frame(apply(toview1,2,function(x){x/rowMeans(toview1[,c(1,2)])}))
  toview1 <- log2(toview1 +0.001)
  toview1 <- toview1[which(rowSums(is.na(toview1)) <1),]
  #
  resdata <- as.data.frame(read.csv(ff01[3], header = T, stringsAsFactors = F))
  upup_04 <- resdata[which(resdata$log2FoldChange > log2(1.5)  &  resdata$padj < 0.05),"X"]
  upup_04 <- intersect(upup_04, twos)
  toview2 <- unique(c(upup_04))
  toview2 <- gencpm2[which(rownames(gencpm2) %in% toview2),]
  toview2 <- toview2 +0.001
  toview2 <- as.data.frame(apply(toview2,2,function(x){x/rowMeans(toview2[,c(1,2)])}))
  toview2 <- log2(toview2 +0.001)
  toview2 <- toview2[which(rowSums(is.na(toview2)) <1),]
  #
  pdf(file = paste0(path,"PLOT/heatmap.add.new.gene.Only-2cellgene-color8.pdf"), height = 6, width = 6, family = "serif")
  pheatmap(mat = toview1,
           main="\nheatmap of genes: add.new.1",
           scale = "none", cellwidth = 20, 
           show_rownames = T, show_colnames = T, cluster_row = T, cluster_cols = F,
           labels_col = colnames(toview1), angle_col = "45",
           breaks = seq(-4,8,by=0.01), border = F,
           color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-4,8,by=0.01))*0.33),
                      colorRampPalette(colors = c("white", "red"))(length(seq(-4,8,by=0.01))*0.67)))
  pheatmap(mat = toview2,
           main="\nheatmap of genes: add.new.2",
           scale = "none", cellwidth = 20, 
           show_rownames = T, show_colnames = T, cluster_row = T, cluster_cols = F,
           labels_col = colnames(toview2), angle_col = "45",
           breaks = seq(-4,8,by=0.01), border = F,
           color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-4,8,by=0.01))*0.33),
                      colorRampPalette(colors = c("white", "red"))(length(seq(-4,8,by=0.01))*0.67)))
  dev.off()
}
#
#
#
#
################################################################################
#36、“ghyuu-4500 vs E14-4500”所有dysregulated gene，分为上调和下调 [采用] [OK]-----------------------------------------------
#####ofAddNew.1、ofAddNew.3
{
  sets <- list(c("group2.E14.P0.4500.Pyru","group4.ghyuu.P0.4500.Pyru","group4.ghyuu.P0.0180.Pyru"),
               c("group3.E14.P0.4500","group5.ghyuu.P0.4500","group5.ghyuu.P0.0180"))
  gencpms <- openxlsx::read.xlsx(paste0(path,"count/huhuhu.cpm.onlyGene.rename.xlsx"))
  gencpm1 <- gencpms[,c(grep(sets[[1]][1],colnames(gencpms)),
                        grep(sets[[1]][2],colnames(gencpms)),
                        grep(sets[[1]][3],colnames(gencpms)))]
  rownames(gencpm1) <- gencpms$theName
  gencpm2 <- gencpms[,c(grep(sets[[2]][1],colnames(gencpms)),
                        grep(sets[[2]][2],colnames(gencpms)),
                        grep(sets[[2]][3],colnames(gencpms)))]
  rownames(gencpm2) <- gencpms$theName
  #
  mark <- paste0("ofAddNew.",seq(4))
  ff01 <- paste0(path, "DESeq2/",mark,".geneOnly.csv")
}
{
  resdata <- as.data.frame(read.csv(ff01[1], header = T, stringsAsFactors = F))
  upup_01 <- resdata[which(resdata$log2FoldChange > log2(1.5)  &  resdata$padj < 0.05),"X"]
  toview1 <- unique(c(upup_01))
  toview1 <- gencpm1[which(rownames(gencpm1) %in% toview1),]
  toview1 <- toview1 +0.001
  toview1 <- as.data.frame(apply(toview1,2,function(x){x/rowMeans(toview1[,c(1,2)])}))
  toview1 <- log2(toview1 +0.001)
  toview1 <- toview1[which(rowSums(is.na(toview1)) <1),]
  #
  resdata <- as.data.frame(read.csv(ff01[1], header = T, stringsAsFactors = F))
  down_01 <- resdata[which(resdata$log2FoldChange < -log2(1.5)  &  resdata$padj < 0.05),"X"]
  toview2 <- unique(c(down_01))
  toview2 <- gencpm2[which(rownames(gencpm2) %in% toview2),]
  toview2 <- toview2 +0.001
  toview2 <- as.data.frame(apply(toview2,2,function(x){x/rowMeans(toview2[,c(1,2)])}))
  toview2 <- log2(toview2 +0.001)
  toview2 <- toview2[which(rowSums(is.na(toview2)) <1),]
  #
  pdf(file = paste0(path,"PLOT/heatmap.add.new.1.dysreguGene-color8.pdf"), height = 6, width = 6, family = "serif")
  pheatmap(mat = toview1,
           main="\nadd.new.1 ghyuu-4500 vs E14-4500 (Pyru) [upup]",
           scale = "none", cellwidth = 20, 
           show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = F,
           labels_col = colnames(toview1), angle_col = "45",
           breaks = seq(-4,8,by=0.01), border = F,
           color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-4,8,by=0.01))*0.33),
                      colorRampPalette(colors = c("white", "red"))(length(seq(-4,8,by=0.01))*0.67)))
  pheatmap(mat = toview2,
           main="\nadd.new.1 ghyuu-4500 vs E14-4500 (Pyru) [down]",
           scale = "none", cellwidth = 20, 
           show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = F,
           labels_col = colnames(toview2), angle_col = "45",
           breaks = seq(-4,8,by=0.01), border = F,
           color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-4,8,by=0.01))*0.33),
                      colorRampPalette(colors = c("white", "red"))(length(seq(-4,8,by=0.01))*0.67)))
  dev.off()
}
{
  resdata <- as.data.frame(read.csv(ff01[3], header = T, stringsAsFactors = F))
  upup_01 <- resdata[which(resdata$log2FoldChange > log2(1.5)  &  resdata$padj < 0.05),"X"]
  toview1 <- unique(c(upup_01))
  toview1 <- gencpm1[which(rownames(gencpm1) %in% toview1),]
  toview1 <- toview1 +0.001
  toview1 <- as.data.frame(apply(toview1,2,function(x){x/rowMeans(toview1[,c(1,2)])}))
  toview1 <- log2(toview1 +0.001)
  toview1 <- toview1[which(rowSums(is.na(toview1)) <1),]
  #
  resdata <- as.data.frame(read.csv(ff01[3], header = T, stringsAsFactors = F))
  down_01 <- resdata[which(resdata$log2FoldChange < -log2(1.5)  &  resdata$padj < 0.05),"X"]
  toview2 <- unique(c(down_01))
  toview2 <- gencpm2[which(rownames(gencpm2) %in% toview2),]
  toview2 <- toview2 +0.001
  toview2 <- as.data.frame(apply(toview2,2,function(x){x/rowMeans(toview2[,c(1,2)])}))
  toview2 <- log2(toview2 +0.001)
  toview2 <- toview2[which(rowSums(is.na(toview2)) <1),]
  #
  pdf(file = paste0(path,"PLOT/heatmap.add.new.2.dysreguGene-color8.pdf"), height = 6, width = 6, family = "serif")
  pheatmap(mat = toview1,
           main="\nadd.new.2 ghyuu-4500 vs E14-4500 [upup]",
           scale = "none", cellwidth = 20, 
           show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = F,
           labels_col = colnames(toview1), angle_col = "45",
           breaks = seq(-4,8,by=0.01), border = F,
           color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-4,8,by=0.01))*0.33),
                      colorRampPalette(colors = c("white", "red"))(length(seq(-4,8,by=0.01))*0.67)))
  pheatmap(mat = toview2,
           main="\nadd.new.2 ghyuu-4500 vs E14-4500 [down]",
           scale = "none", cellwidth = 20, 
           show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = F,
           labels_col = colnames(toview2), angle_col = "45",
           breaks = seq(-4,8,by=0.01), border = F,
           color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-4,8,by=0.01))*0.33),
                      colorRampPalette(colors = c("white", "red"))(length(seq(-4,8,by=0.01))*0.67)))
  dev.off()
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
################################################################################
#37、GO/KEGG [11月9号] [KEGG画气泡图]
#####2 vs 1 上调，3 vs 2 下调【编号-1】
#####2 vs 1 下调，3 vs 2 上调【编号-2】
sets <- list(c("group4.ghyuu.P0.4500.Pyru","group2.E14.P0.4500.Pyru"),#对应ofAddNew.1, 即2 vs 1
             c("group4.ghyuu.P0.0180.Pyru","group4.ghyuu.P0.4500.Pyru"),#对应ofAddNew.2, 即3 vs 2
             c("group5.ghyuu.P0.4500","group3.E14.P0.4500"),#对应ofAddNew.3, 即2 vs 1
             c("group5.ghyuu.P0.0180","group5.ghyuu.P0.4500"))#对应ofAddNew.4, 即3 vs 2
mark <- paste0("ofAddNew.",seq(4))
ff01 <- paste0(path, "DESeq2/",mark,".geneOnly.csv")
#add.new.1
{
  #ghyuu 4500 vs E14 4500的上调基因与ghyuu 180 vs ghyuu 4500的下调基因的VENN
  resdata <- as.data.frame(read.csv(ff01[1], header = T, stringsAsFactors = F))
  upup_01 <- resdata[which(resdata$log2FoldChange >  log2(1.5)  &  resdata$padj < 0.05),"X"]
  upup_01 <- unique(upup_01)
  resdata <- as.data.frame(read.csv(ff01[2], header = T, stringsAsFactors = F))
  down_01 <- resdata[which(resdata$log2FoldChange < -log2(1.5)  &  resdata$padj < 0.05),"X"]
  down_01 <- unique(down_01)
  #
  #
  tmp1 <- list()
  tmp1[["a"]] <- upup_01
  tmp1[["b"]] <- down_01
  p_001 <- VennDiagram::venn.diagram(x = tmp1,
                                     filename = NULL,
                                     main = "\nadd.new.1 (Pyru)",
                                     main.fontfamily = "serif",
                                     sub  = "\n\nKO4500_WT4500-upup___KO180_KO4500-down\n",
                                     sub.fontfamily = "serif",
                                     cex = 1.5, main.cex = 1.5, sub.cex = 1,
                                     print.mode = "raw",#c("percent", "raw"),
                                     category.names = c("KO4500_WT4500-upup",
                                                        "\n\n\n\n\n\n\n\n\n\nKO180_KO4500-down"),
                                     cat.cex = 1.5,cat.dist = -0.05,
                                     fill=c("#eaafc8","#00b4db"),
                                     units = "cm", height = 12, width = 12)
  pdf(file = paste0(path,"PLOT/","venn.add.new.1.up_down",".pdf"), width = 8, height = 8)
  grid.draw(p_001)
  dev.off()
}
{
  #ghyuu 4500 vs E14 4500的下调基因与ghyuu 180 vs ghyuu 4500的上调基因的VENN
  resdata <- as.data.frame(read.csv(ff01[1], header = T, stringsAsFactors = F))
  down_02 <- resdata[which(resdata$log2FoldChange < -log2(1.5)  &  resdata$padj < 0.05),"X"]
  down_02 <- unique(down_02)
  resdata <- as.data.frame(read.csv(ff01[2], header = T, stringsAsFactors = F))
  upup_02 <- resdata[which(resdata$log2FoldChange >  log2(1.5)  &  resdata$padj < 0.05),"X"]
  upup_02 <- unique(upup_02)
  #
  #
  tmp2 <- list()
  tmp2[["a"]] <- down_02
  tmp2[["b"]] <- upup_02
  p_002 <- VennDiagram::venn.diagram(x = tmp2,
                                     filename = NULL,
                                     main = "\nadd.new.1 (Pyru)",
                                     main.fontfamily = "serif",
                                     sub  = "\n\nKO4500_WT4500-down___KO180_KO4500-upup\n",
                                     sub.fontfamily = "serif",
                                     cex = 1.5, main.cex = 1.5, sub.cex = 1,
                                     print.mode = "raw",#c("percent", "raw"),
                                     category.names = c("KO4500_WT4500-down",
                                                        "\n\n\n\n\n\n\n\n\n\nKO180_KO4500-upup"),
                                     cat.cex = 1.5,cat.dist = -0.05,
                                     fill=c("#eaafc8","#00b4db"),
                                     units = "cm", height = 12, width = 12)
  pdf(file = paste0(path,"PLOT/","venn.add.new.1.down_up",".pdf"), width = 8, height = 8)
  grid.draw(p_002)
  dev.off()
}
#add.new.2
{
  #ghyuu 4500 vs E14 4500的上调基因与ghyuu 180 vs ghyuu 4500的下调基因的VENN
  resdata <- as.data.frame(read.csv(ff01[3], header = T, stringsAsFactors = F))
  upup_01 <- resdata[which(resdata$log2FoldChange >  log2(1.5)  &  resdata$padj < 0.05),"X"]
  upup_01 <- unique(upup_01)
  resdata <- as.data.frame(read.csv(ff01[4], header = T, stringsAsFactors = F))
  down_01 <- resdata[which(resdata$log2FoldChange < -log2(1.5)  &  resdata$padj < 0.05),"X"]
  down_01 <- unique(down_01)
  #
  #
  tmp3 <- list()
  tmp3[["a"]] <- upup_01
  tmp3[["b"]] <- down_01
  p_001 <- VennDiagram::venn.diagram(x = tmp3,
                                     filename = NULL,
                                     main = "\nadd.new.2",
                                     main.fontfamily = "serif",
                                     sub  = "\n\nKO4500_WT4500-upup___KO180_KO4500-down\n",
                                     sub.fontfamily = "serif",
                                     cex = 1.5, main.cex = 1.5, sub.cex = 1,
                                     print.mode = "raw",#c("percent", "raw"),
                                     category.names = c("KO4500_WT4500-upup",
                                                        "\n\n\n\n\n\n\n\n\n\nKO180_KO4500-down"),
                                     cat.cex = 1.5,cat.dist = -0.05,
                                     fill=c("#eaafc8","#00b4db"),
                                     units = "cm", height = 12, width = 12)
  pdf(file = paste0(path,"PLOT/","venn.add.new.2.up_down",".pdf"), width = 8, height = 8)
  grid.draw(p_001)
  dev.off()
}
{
  #ghyuu 4500 vs E14 4500的下调基因与ghyuu 180 vs ghyuu 4500的上调基因的VENN
  resdata <- as.data.frame(read.csv(ff01[3], header = T, stringsAsFactors = F))
  down_02 <- resdata[which(resdata$log2FoldChange < -log2(1.5)  &  resdata$padj < 0.05),"X"]
  down_02 <- unique(down_02)
  resdata <- as.data.frame(read.csv(ff01[4], header = T, stringsAsFactors = F))
  upup_02 <- resdata[which(resdata$log2FoldChange >  log2(1.5)  &  resdata$padj < 0.05),"X"]
  upup_02 <- unique(upup_02)
  #
  #
  tmp4 <- list()
  tmp4[["a"]] <- down_02
  tmp4[["b"]] <- upup_02
  p_002 <- VennDiagram::venn.diagram(x = tmp4,
                                     filename = NULL,
                                     main = "\nadd.new.2",
                                     main.fontfamily = "serif",
                                     sub  = "\n\nKO4500_WT4500-down___KO180_KO4500-upup\n",
                                     sub.fontfamily = "serif",
                                     cex = 1.5, main.cex = 1.5, sub.cex = 1,
                                     print.mode = "raw",#c("percent", "raw"),
                                     category.names = c("KO4500_WT4500-down",
                                                        "\n\n\n\n\n\n\n\n\n\nKO180_KO4500-upup"),
                                     cat.cex = 1.5,cat.dist = -0.05,
                                     fill=c("#eaafc8","#00b4db"),
                                     units = "cm", height = 12, width = 12)
  pdf(file = paste0(path,"PLOT/","venn.add.new.2.down_up",".pdf"), width = 8, height = 8)
  grid.draw(p_002)
  dev.off()
}
#P值的计算
{
  shyP <- function(inter, bg, righ, left) {
    a=righ+inter
    b=left+inter
    pzhi <- stats::phyper(inter-1, a, bg-a, b, lower.tail = F)
    print(format(pzhi,scientific=T))
  }
  shyP(inter = 818, righ = 1060, left = 2354, bg = 20000)
  shyP(inter = 156, righ = 583,  left = 3004, bg = 20000)
  shyP(inter = 482, righ = 657,  left = 2816, bg = 20000)
  shyP(inter = 113, righ = 310,  left = 3618, bg = 20000)
}
#
#
#对上述#add.new.1、#add.new.2的【编号-1】的venn基因做富集分析
{
  for (mark in seq(2)) {
    gggg <- list(base::intersect(x = tmp1[["a"]], y = tmp1[["b"]]),
                 base::intersect(x = tmp3[["a"]], y = tmp3[["b"]]))
    ov01 <- gggg[[mark]]
    lens <- length(ov01)
    geyy <- na.omit(mapIds(x=org.Mm.eg.db, keys=unique(ov01), keytype="SYMBOL", column="ENTREZID"))
    gogo <- clusterProfiler::enrichGO(gene = geyy, OrgDb = org.Mm.eg.db, 
                                      keyType = "ENTREZID", 
                                      ont = "BP",pvalueCutoff = 0.05,qvalueCutoff = 0.05, 
                                      readable = T)
    gogo <- gogo@result
    openxlsx::write.xlsx(x = gogo, file = paste0(path,"PLOT/ofAddNew.",mark,"/GO.BP.Venn_",lens,".xlsx"))
    #
    gogo <- gogo[1:10,]
    ifelse(mark==1,gogo$Description[7] <- "negative regulation of \nnervous system development",gogo$Description[7])
    ifelse(mark==2,gogo$Description[8] <- "positive regulation of \nsmooth muscle cell proliferation",gogo$Description[8])
    gogo$Description <- factor(gogo$Description, levels = rev(gogo$Description))
    p_001 <- ggplot(gogo)+
      geom_bar(aes(x=-log10(pvalue), y=Description), stat = "identity", fill="#1ba784", alpha=1.0) +
      ggthemes::theme_few() +
      scale_x_continuous(breaks = c(0,seq(9))) +
      labs(y="",x="-Log10 (P-value)",
           title = paste0("this is title--venn_",lens)) +
      geom_vline(aes(xintercept = 1), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 2), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 3), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 4), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 5), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 6), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 7), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 8), color = "white", linewidth  = 0.8) +
      theme(text = element_text(size = 12,family = "serif",color = "black"),
            panel.border = element_rect(linewidth = 1),
            plot.title = element_text(size = 14,family = "serif",color = "#1ba784",face = "bold"),
            axis.title = element_text(size = 14,family = "serif",color = "#1ba784"),
            axis.text  = element_text(size = 12,family = "serif",color = "#1ba784"))
    p_001
    ggsave(plot = p_001,
           units = "cm", width = 20, height = 10,
           filename = paste0(path,"PLOT/ofAddNew.",mark,"/GO.BP.Venn_",lens,".pdf"))
    #
    #
    #KEGG
    geyy <- na.omit(mapIds(x=org.Mm.eg.db, keys=unique(ov01), keytype="SYMBOL", column="ENTREZID"))
    kegg <- clusterProfiler::enrichKEGG(gene = geyy, organism = "mmu", keyType = "kegg", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
    kegg <- DOSE::setReadable(kegg, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
    kegg <- kegg@result
    kegg$Description <- sapply(str_split(kegg$Description, pattern = " - Mus"), "[", 1)
    openxlsx::write.xlsx(x = kegg, file = paste0(path,"PLOT/ofAddNew.",mark,"/KEGG.Venn_",lens,".xlsx"))
    #
    kegg <- kegg[1:10,]
    ifelse(mark==2,kegg$Description[8] <- "Glycosaminoglycan biosynthesis \nkeratan sulfate",kegg$Description[8])
    ifelse(mark==2,kegg$Description[10] <- "Arrhythmogenic right ventricular \ncardiomyopathy",kegg$Description[10])
    kegg$Description <- factor(kegg$Description, levels = rev(kegg$Description))
    p_002 <- ggplot(kegg, aes(x=-log10(pvalue), y=Description))+
      geom_point(aes(color=-log10(pvalue),size=Count)) +
      theme_classic() +
      labs(title = paste0("this is title--venn_",lens)) +
      scale_color_gradient2(low = "white", high = "red3", midpoint = 1) +
      theme(text = element_text(size = 12,family = "serif",color = "black"),
            axis.title = element_text(size = 12,family = "serif",color = "black"),
            axis.text  = element_text(size = 12,family = "serif",color = "black"))
    p_002
    ggsave(plot = p_002,
           units = "cm", width = 16, height = 10,
           filename = paste0(path,"PLOT/ofAddNew.",mark,"/KEGG.Venn_",lens,".pdf"))
  }
}
################################################################################
#38、GO/KEGG [11月9号] [KEGG画气泡图] [for group1] [OK]
{
  #GO
  {
    resFile <- paste0(path, "DESeq2/huhuhu.group1.geneOnly.csv")
    resdata <- as.data.frame(read.csv(resFile, header = T, stringsAsFactors = F))
    resdata <- resdata[grep(resdata[,1], pattern = ":", invert = T),]
    mtacup <- resdata[which(resdata$log2FoldChange > log2(1.5)  &  resdata$padj < 0.05),"X"]
    kddown <- resdata[which(resdata$log2FoldChange < -log2(1.5) &  resdata$padj < 0.05),"X"]
    for (symb in seq(2)) {
      ovgene <- unlist(list(mtacup,kddown)[symb])
      geyy <- na.omit(mapIds(x=org.Mm.eg.db, keys=unique(ovgene), keytype="SYMBOL", column="ENTREZID"))
      gogo <- clusterProfiler::enrichGO(gene = geyy, OrgDb = org.Mm.eg.db, 
                                        keyType = "ENTREZID", 
                                        ont = "BP",pvalueCutoff = 0.05,qvalueCutoff = 0.05, 
                                        readable = T)
      gogo <- gogo@result
      assign(paste0("gta",symb), gogo[which(gogo$pvalue<0.05),])
    }
    openxlsx::write.xlsx(x = gta1, file = paste0(path,"PLOT/group1/GO.BP.upup.xlsx"))
    openxlsx::write.xlsx(x = gta2, file = paste0(path,"PLOT/group1/GO.BP.down.xlsx"))
    #
    gta1$Description <- factor(gta1$Description, levels = rev(gta1$Description))
    p_003 <- ggplot(gta1[1:10,])+
      geom_bar(aes(x=-log10(pvalue), y=Description), stat = "identity", fill="#1ba784", alpha=1.0) +
      ggthemes::theme_few() +
      scale_x_continuous(breaks = seq(0,75,by=10)) +
      labs(y="",title = "this is title--upup") +
      geom_vline(aes(xintercept = 5),  color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 10), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 15), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 20), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 25), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 30), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 35), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 40), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 45), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 50), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 55), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 60), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 65), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 70), color = "white", linewidth  = 0.8) +
      theme(text = element_text(size = 12,family = "serif",color = "black"),
            panel.border = element_rect(linewidth = 1),
            plot.title = element_text(size = 12,family = "serif",color = "black",face = "bold"),
            axis.title = element_text(size = 12,family = "serif",color = "black"),
            axis.text  = element_text(size = 12,family = "serif",color = "black"))
    p_003
    #
    gta2$Description <- factor(gta2$Description, levels = rev(gta2$Description))
    p_004 <- ggplot(gta2[1:10,])+
      geom_bar(aes(x=-log10(pvalue), y=Description), stat = "identity", fill="#5861AC", alpha=1.0) +
      ggthemes::theme_few() +
      scale_x_continuous(breaks = seq(0,16,by=3)) +
      geom_vline(aes(xintercept = 1),  color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 2),  color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 3),  color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 4),  color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 5),  color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 6),  color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 7),  color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 8),  color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 9),  color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 10), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 11), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 12), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 13), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 14), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 15), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 16), color = "white", linewidth  = 0.8) +
      labs(y="",title = "this is title--down") +
      theme(text = element_text(size = 12,family = "serif",color = "black"),
            panel.border = element_rect(linewidth = 1),
            plot.title = element_text(size = 12,family = "serif",color = "black",face = "bold"),
            axis.title = element_text(size = 12,family = "serif",color = "black"),
            axis.text  = element_text(size = 12,family = "serif",color = "black"))
    plot <- p_003 +p_004 +plot_layout(nrow = 2)
    ggsave(plot = plot,
           units = "cm", width = 18, height = 16,
           filename = paste0(path,"PLOT/group1/group1.GO.BP.pdf"))
  }
  #KEGG
  {
    resFile <- paste0(path, "DESeq2/huhuhu.group1.geneOnly.csv")
    resdata <- as.data.frame(read.csv(resFile, header = T, stringsAsFactors = F))
    resdata <- resdata[grep(resdata[,1], pattern = ":", invert = T),]
    mtacup <- resdata[which(resdata$log2FoldChange > log2(1.5)  &  resdata$padj < 0.05),"X"]
    kddown <- resdata[which(resdata$log2FoldChange < -log2(1.5) &  resdata$padj < 0.05),"X"]
    for (symb in seq(2)) {
      ovgene <- unlist(list(mtacup,kddown)[symb])
      geyy <- na.omit(mapIds(x=org.Mm.eg.db, keys=unique(ovgene), keytype="SYMBOL", column="ENTREZID"))
      kegg <- clusterProfiler::enrichKEGG(gene = geyy, organism = "mmu", keyType = "kegg", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
      kegg <- DOSE::setReadable(kegg, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
      kegg <- kegg@result
      kegg$Description <- sapply(str_split(kegg$Description, pattern = " - Mus"), "[", 1)
      assign(paste0("kta",symb), kegg[which(kegg$pvalue<0.05),])
    }
    openxlsx::write.xlsx(x = kta1, file = paste0(path,"PLOT/group1/kegg.upup.xlsx"))
    openxlsx::write.xlsx(x = kta2, file = paste0(path,"PLOT/group1/kegg.down.xlsx"))
    #
    kta1$Description <- factor(kta1$Description, levels = rev(kta1$Description))
    p_005 <- ggplot(kta1[1:10,], aes(x=-log10(pvalue), y=Description))+
      geom_point(aes(color=-log10(pvalue),size=Count)) +
      theme_classic() +
      labs(title = paste0("this is title--upup")) +
      scale_color_gradient2(low = "white", high = "red3", midpoint = 1) +
      theme(plot.margin = margin(unit = "cm", 1,1,1,1),
            text = element_text(size = 12,family = "serif",color = "black"),
            axis.title = element_text(size = 12,family = "serif",color = "black"),
            axis.text  = element_text(size = 12,family = "serif",color = "black"))
    #
    kta2$Description <- factor(kta2$Description, levels = rev(kta2$Description))
    p_006 <- ggplot(kta2[1:10,], aes(x=-log10(pvalue), y=Description))+
      geom_point(aes(color=-log10(pvalue),size=Count)) +
      theme_classic() +
      labs(title = paste0("this is title--down")) +
      scale_color_gradient2(low = "white", high = "red3", midpoint = 1) +
      theme(plot.margin = margin(unit = "cm", 1,1,1,1),
            text = element_text(size = 12,family = "serif",color = "black"),
            axis.title = element_text(size = 12,family = "serif",color = "black"),
            axis.text  = element_text(size = 12,family = "serif",color = "black"))
    plot <- p_005 +p_006 +plot_layout(nrow = 2)
    ggsave(plot = plot,
           units = "cm", width = 22, height = 20,
           filename = paste0(path,"PLOT/group1/group1.KEGG.pdf"))
  }
}
################################################################################
#39、在TSC基因上画GSEA [for group1]
{
  #gggg <- read.table("/disk5/zx/ChIP_seq/H3K9me12_E14_WT_ZX/TSCs.bed",
  #                   sep = "\t", header = F)
  #gggg <- as.data.frame(t(c("this is TSC genes","NA",gggg$V4)))
  #write.table(x = gggg,quote = F,sep = "\t",col.names = F, row.names = F,
  #            file = "/Reference/aaSHY/GSEAgmt/TSCs.gmt")
  #
  #
  #
  #
  path = "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/thehuhuhu/"
  #temp <- read.csv(paste0(path, "DESeq2/huhuhu.group1.geneOnly.csv"), header = T)
  #temp <- temp[, c("X", "log2FoldChange")]
  #temp <- temp[order(temp$log2FoldChange, decreasing = T),]
  #write.table(x = temp, 
  #            file = paste0(path, "DESeq2/huhuhu.group1.res.rnk"),
  #            quote = F, sep = "\t", col.names = F, row.names = F)
  cmd_12 <- paste0("gseapy prerank -g ",
                   "/Reference/aaSHY/GSEAgmt/TSCs.gmt -r ",
                   paste0(path, "DESeq2/huhuhu.group1.res.rnk"),
                   " --max-size 500 -p 10 -o ",
                   paste0(path, "PLOT/group1/enrichGSEA/"))#conda activate base gseapy 0.10.8
  print(cmd_12)
}
################################################################################
#40、火山图 [for group1] [OK]
{
  shy_volcano_gene <- function(resFile, saveFile){
    resdata <- as.data.frame(read.csv(resFile, header = T, stringsAsFactors = F))
    resdata <- na.omit(resdata)
    resdata$threshold[resdata$log2FoldChange > log2(1.5) & resdata$padj <0.05] = "Upregulated genes"
    resdata$threshold[resdata$log2FoldChange < -log2(1.5) & resdata$padj <0.05] = "Downregulated genes"
    resdata$threshold[!(resdata$log2FoldChange > log2(1.5) & resdata$padj <0.05) & !(resdata$log2FoldChange < -log2(1.5) & resdata$padj < 0.05)] = "Other genes"
    resdata$threshold <- factor(resdata$threshold, levels = c("Upregulated genes","Downregulated genes","Other genes"))
    p_1 <- ggplot(resdata, aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) +
      labs(col="", x="Log2FoldChange",y="-Log10(adjusted p-value)") + 
      geom_point(size=0.8) +
      scale_x_continuous(limits = c(-9, 9), breaks = seq(-9, 9, 2)) +
      geom_hline(aes(yintercept=-log10(0.05)), linetype="dashed", colour="royalblue") +
      geom_vline(aes(xintercept=-log2(1.5)), linetype="dashed", colour="royalblue") +
      geom_vline(aes(xintercept=log2(1.5)), linetype="dashed", colour="royalblue") +
      theme_bw() + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      scale_color_manual(values = c("red","blue","black"),
                         labels = c(paste0(names(table(resdata$threshold)[1])," (",unname(table(resdata$threshold)[1]),")"),
                                    paste0(names(table(resdata$threshold)[2])," (",unname(table(resdata$threshold)[2]),")"),
                                    paste0(names(table(resdata$threshold)[3])," (",unname(table(resdata$threshold)[3]),")"))) +
      theme(legend.position = c(0.25,0.9), 
            legend.background = element_rect(fill = NA),
            legend.spacing.x = unit(0,"mm"),
            legend.spacing.y = unit(0,"mm"),
            legend.key = element_blank(),
            legend.text = element_text(size = 13),
            axis.title = element_text(size = 13)) +
      guides(color = guide_legend(override.aes = list(size = 3)))
    ggsave(plot = p_1, saveFile, units = "cm", width = 18, height = 16)
  }
  shy_volcano_gene(resFile  = paste0(path, "DESeq2/huhuhu.group1.geneOnly.csv"),
                   saveFile = paste0(path, "PLOT/group1/group1.gene.volcano.pdf"))
}
################################################################################
#41、Boxplot [LINE1, SINE B1, SINE B2, ERV1, ERVK, ERVL, ERVL-MaLR]
{
  mervs <- openxlsx::read.xlsx(paste0(path,"countUniq/TEtools/TEtools_huhuhu.cpm.onlyTE.site.xlsx"))
  gggg <- data.frame(TE = sapply(str_split(mervs$theName,":"),"[",2),
                     family = sapply(str_split(mervs$theName,":"),"[",3),
                     class = sapply(str_split(mervs$theName,":"),"[",4))
  gggg <- as.data.frame(distinct(gggg))
  View(gggg)
  
  
  shy_d <- function(grou="group1",mark) {
    mervl <- mervs[,grep(paste0(grou,"|theName"),colnames(mervs))]
    mervl <- mervl[grep(mark,sapply(str_split(mervl$theName,":"),"[",4)),]
    mervl <- mervl[rowSums(mervl[,1:4]==0) < 4,]
    mervl <- tidyr::pivot_longer(data = mervl, cols = -theName, names_to = "group", values_to = "cpms")
    vaue <- wilcox.test(x = mervl$cpms[grep("0000|0180",mervl$group)],y = mervl$cpms[grep("4500",mervl$group)], alternative = "less")
    print(vaue$p.value)
    #
    mervl$marks <- NA
    mervl$marks[grep("0000|0180",mervl$group)] <- "lowshuhuhucose"
    mervl$marks[grep("4500",mervl$group)] <- "highhuhuhucose"
    mervl$marks <- factor(mervl$marks, levels = rev(unique(mervl$marks)))
    p_001 <- ggplot(data = mervl,aes(x = marks, y = log2(cpms+0.1))) +
      geom_boxplot(aes(fill = marks),outlier.alpha = 1) +
      labs(title = paste0(mark," site Expr (",grou,")"),
           y = "Log2 (CPM+0.1)", x = "", fill="Group") +
      geom_signif(comparisons = list(c("lowshuhuhucose", "highhuhuhucose")),
                  map_signif_level=T, family = "serif",
                  textsize=6, test = wilcox.test, step_increase = 0) +
      stat_summary(geom = "point", fun = "median",shape = 22,size = 2,
                   position = position_dodge(width = 1)) +
      theme_classic() +
      scale_fill_manual(values = c("#F0A19A","#3FA0C0"),
                        labels = c("high huhuhucose", "low huhuhucose")) +
      theme(plot.margin = margin(unit = "cm",1,1,1,1),
            legend.text = element_text(family = "serif", size = 12),
            legend.title = element_text(family = "serif", size = 12),
            axis.text.x = element_text(angle = 45,hjust = 1,family = "serif",size = 12),
            axis.text.y = element_text(family = "serif", size = 12)) +
      #scale_y_continuous(limits = c(-5,10)) +
      #scale_y_continuous(limits = c(0,10),breaks = c(0,1,2,3,4,5)) +
      scale_x_discrete(labels = c("high huhuhucose", "low huhuhucose"))
    print(p_001)
    ggsave(plot = p_001,
           filename = paste0(path,"PLOT_uniq/group1/",mark,".boxplot.pdf"),
           units = "cm", width = 20, height = 15)
  }
  for (mark in c("LINE","SINE","LTR")) {
    shy_d(mark = mark)
  }
}
################################################################################
#42、PCA with Dux [for group1]
{
  path = "/ChIP_seq_2/aaSHY/ghyuu/publicData/D1_2C_cells/"
  for (i in c("src","fq","Log","bamCoverage","trim","dumpROOM","bam","count")) {
    if (dir.exists(paths = paste0(path, i))==F) {
      dir.create(path = paste0(path,i), recursive = T)
    }
  }
  shell <- paste0(path,"src/run_pca.sh")
  cat("#!/bin/bash\n", file = shell)
  prex <- gsub(x = list.files(path = paste0(path,"rawdata"), pattern = ".sra"),
               pattern = ".sra", replacement = "")
  lay1 <- paste0(path,"fq/",prex,".fq.gz")
  lay2 <- paste0(path,"trim/",prex,"_trimmed.fq.gz")
  cmd_01 <- paste0("parallel-fastq-dump -t 10 --split-3 --gzip -s ",path,
                   "rawdata/",prex,".sra -O ", path, "fq", "\n")
  cmd_02 <- paste0("rename s/fastq/fq/ ", path, "fq/*gz", "\n")
  cmd_03 <- paste0("fastqc -t 4 -o ",path,"fq/ ", lay1,"\n")
  cmd_04 <- paste0("trim_galore --phred33 -j 2 ",
                   "--fastqc --clip_R1 9 --three_prime_clip_R1 1 -o ",path,
                   "trim ", lay1, "\n")
  cmd_05 <- paste0("STAR --runThreadN 4 --readFilesCommand zcat ",
                   "--outMultimapperOrder Random --outSAMtype BAM SortedByCoordinate ",
                   "--outFileNamePrefix ",path,"bam/",prex," ",
                   "--outSAMprimaryFlag AllBestScore --twopassMode Basic ",
                   "--genomeDir ",inx1," --readFilesIn ",lay2,"\n")
  cmd_06 <- paste0("mv ",path,"bam/",prex,"Aligned.sortedByCoord.out.bam ",path,
                   "bam/",prex,".bam","\n")
  cmd_07 <- paste0("samtools index ",path,"bam/",prex,".bam","\n")
  cmd_08 <- paste0("bamCoverage --normalizeUsing CPM -p 10 --minMappingQuality 1 ",
                   "--samFlagExclude 256 -of bigwig -bs 20 -b ",
                   paste0(path,"bam/",prex,".bam")," -o ",paste0(path,"bamCoverage/",prex,".bw"),"\n")
  cmd_09 <- paste0("TEtranscripts -t ",
                   paste(paste0(path,"bam/",prex[1:3],".bam"),collapse = " ")," -c ",
                   paste(paste0(path,"bam/",prex[4:6],".bam"),collapse = " "),
                   " --GTF ",inx3," --TE ",inx4," --sortByPos --mode multi ",
                   "--project duxs --outdir ",path,"count","\n")
  cmd_10 <- paste0("TElocal -b ",
                   path,"bam/",prex,".bam",
                   " --GTF ",inx3," --TE ",inx5," --sortByPos --project ",
                   paste0(path, "count/",prex),"\n")
  for (i in cmd_01) {cat(i, file = shell, append = T)}
  for (i in cmd_02) {cat(i, file = shell, append = T)}
  for (i in cmd_03) {cat(i, file = shell, append = T)}
  for (i in cmd_04) {cat(i, file = shell, append = T)}
  for (i in cmd_05) {cat(i, file = shell, append = T)}
  for (i in cmd_06) {cat(i, file = shell, append = T)}
  for (i in cmd_07) {cat(i, file = shell, append = T)}
  for (i in cmd_08) {cat(i, file = shell, append = T)}
  for (i in cmd_09) {cat(i, file = shell, append = T)}
  for (i in cmd_10) {cat(i, file = shell, append = T)}
  print(paste0("nohup bash ",shell, " >",paste0(path,"Log/"),basename(shell),".log"," 2>&1 &"))
  path = "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/thehuhuhu/"
}
################################################################################
#
#
#
################################################################################
#01、对于基因, 使用featureCounts定量
{
  shell <- paste0(path,"src/run_a.sh")
  cat("#!/bin/bash\n", file = shell)
  ff01 <- list.files(path = paste0(path,"fq"),pattern = "_1.fq.gz")
  prex <- sapply(str_split(ff01,pattern = "_1.fq.gz"),"[",1)
  cmd_01 <- paste0("#featureCounts -T 20 -g 'gene_name' ",
                   "-F 'GTF' -O -C -p -P --countReadPairs -d 50 -D 1000 ",
                   "-B -a ",inx3," -o ", path,"countFeature/genes-1.txt ",
                   paste(paste0(path,"bam/",prex, ".bam"),collapse = " "),"\n")
  cmd_01 <- paste0("featureCounts -T 20 -g 'gene_name' ",
                   "-F 'GTF' -p ",
                   "-a ",inx3," -o ", path,"countFeature/genes-1.txt ",
                   paste(paste0(path,"bam/",prex, ".bam"),collapse = " "),"\n")
  cmd_02 <- paste0("cut -f 1,7- ",path,"countFeature/genes-1.txt >",path,
                   "countFeature/genes-2.txt","\n")
  for (i in cmd_01) {cat(i, file = shell, append = T)}
  for (i in cmd_02) {cat(i, file = shell, append = T)}
  print(paste0("nohup bash ",shell, " >",paste0(path,"Log/"),basename(shell),".log"," 2>&1 &"))
}
################################################################################
#02、DESeq2
{
  shy_diff <- function(mark, geneFile, Class, VS){
    tmp1 <- read.table(geneFile, header = T, check.names = F, stringsAsFactors = F, row.names = 1)[,1:4]
    colnames(tmp1) <- sapply(str_split(basename(colnames(tmp1)),pattern = ".bam"), "[", 1)
    dfclas <- data.frame(row.names = colnames(tmp1), Class, stringsAsFactors = T)
    colLen <- ceiling(length(colnames(tmp1))/2)
    #
    tmp2 <- DESeq2::DESeqDataSetFromMatrix(countData = tmp1, colData = dfclas, design = ~Class)
    #tmp2 <- tmp2[rowSums(BiocGenerics::counts(tmp2) > 2) >= colLen, ]
    tmp2 <- tmp2[rowMeans(BiocGenerics::counts(tmp2)) > 5,]
    tmp2 <- DESeq2::DESeq(tmp2,)
    tmp2 <- DESeq2::results(tmp2, contrast = VS)
    write.csv(tmp2, paste0(path,"countFeature/DESeq2/",mark, ".geneOnly.csv"))
    upup <- subset(tmp2, log2FoldChange >  log2(1.5) & padj < 0.05)
    down <- subset(tmp2, log2FoldChange < -log2(1.5) & padj < 0.05)
    write.csv(upup,paste0(path,"countFeature/DESeq2/",mark, ".geneOnly.upup.csv"))
    write.csv(down,paste0(path,"countFeature/DESeq2/",mark, ".geneOnly.down.csv"))
  }
  shy_diff(mark = "ghyuu", 
           geneFile = paste0(path, "countFeature/genes-2.txt"),
           Class = factor(c("KO","KO","WT","WT")), VS = c("Class","KO","WT"))
}
################################################################################
#03、GSEA
{
  temp <- read.csv(paste0(path, "DESeq2/huhuhu.group1.geneOnly.csv"), header = T)
  temp <- temp[, c("X", "log2FoldChange")]
  temp <- temp[order(temp$log2FoldChange, decreasing = T),]
  write.table(x = temp, 
              file = paste0(path, "DESeq2/huhuhu.group1.res.rnk"), 
              quote = F, sep = "\t", col.names = F, row.names = F)
  cmd_12 <- paste0("gseapy prerank -g ","/Reference/aaSHY/GSEAgmt/TwoGenes.gmt -r ",
                   paste0(path, "DESeq2/huhuhu.group1.res.rnk"), " -p 10 -o ",
                   paste0(path, "PLOT/group1/enrichGSEA/"))#conda activate base gseapy 0.10.8
  print(cmd_12)
}
{
  cmd_12 <- paste0("gseapy prerank -g ",
                   "/Reference/aaSHY/GSEAgmt/2022_majorZGA.gmt -r ",
                   paste0(path, "DESeq2/huhuhu.group1.res.rnk"),
                   " --max-size 1000 -p 10 -o ",
                   paste0(path, "PLOT/group1/enrichGSEA/"))#conda activate base gseapy 0.10.8
  print(cmd_12)
}
{
  ###gggg <- read.csv("/ChIP_seq_2/aaSHY/Prmt1/ask/MERVL/mervlKD/DESeq2/onlyGene.mervl.down.csv")[,1]
  ###gggg <- as.data.frame(t(c("MERVL.KD.Down","NA",gggg)))
  ###write.table(x = gggg,quote = F,sep = "\t",col.names = F, row.names = F,
  ###            file = "/ChIP_seq_2/aaSHY/Prmt1/ask/MERVL/mervlKD/DESeq2/onlyGene.mervl.down.gmt")
  #
  #
  cmd_12 <- paste0("gseapy prerank -g ",
                   "/ChIP_seq_2/aaSHY/Prmt1/ask/MERVL/mervlKD/DESeq2/onlyGene.mervl.down.gmt -r ",
                   paste0(path, "DESeq2/huhuhu.group1.res.rnk"), 
                   " --max-size 1000 -p 10 -o ",
                   paste0(path, "PLOT/group1/enrichGSEA/"))#conda activate base gseapy 0.10.8
  print(cmd_12)
}
{
  ###gggg <- read.csv("/ChIP_seq_2/aaSHY/Prmt1/ask/MERVL/mervlKD/DESeq2/onlyGene.mt2.down.csv")[,1]
  ###gggg <- as.data.frame(t(c("MT2.KD.Down","NA",gggg)))
  ###write.table(x = gggg,quote = F,sep = "\t",col.names = F, row.names = F,
  ###            file = "/ChIP_seq_2/aaSHY/Prmt1/ask/MERVL/mervlKD/DESeq2/onlyGene.mt2.down.gmt")
  #
  #
  cmd_12 <- paste0("gseapy prerank -g ",
                   "/ChIP_seq_2/aaSHY/Prmt1/ask/MERVL/mervlKD/DESeq2/onlyGene.mt2.down.gmt -r ",
                   paste0(path, "DESeq2/huhuhu.group1.res.rnk"), " -p 10 -o ",
                   paste0(path, "PLOT/group1/enrichGSEA/"))#conda activate base gseapy 0.10.8
  print(cmd_12)
}
















