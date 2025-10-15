#write by Sun Haiayng at 2024.08.30
#01、搭环境
{
  rm(list=ls())
  for (i in c("data.table","ggridges","ggrepel","KEGG.db","scales","stringr","Biostrings","tidyr","pheatmap","enrichplot","DESeq2","patchwork","ggpubr","ggsci","ggtranscript","rtracklayer","Mfuzz","BiocManager","dplyr","IRanges","GenomicRanges","ReactomePA","ggplot2","clusterProfiler","topGO","Rgraphviz","pathview","org.Mm.eg.db")){library(i, character.only = T,quietly = T)}; rm(i)
}
#02、给索引
{
  path = "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/Aef/"
  inx1 = "/Reference/aaSHY/Ensembl99/GRCm38/INDEX/STARv279a"
  inx3 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/Mus_musculus.GRCm38.99.gtf"
  inx4 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf"
  inx5 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf.locInd"
}
#03、跑命令
{
  for (i in c("src","fq","Log","stringtie/count","bamCoverage","trim","DESeq2","dumpROOM","bam","count","PLOT")) {
    if (dir.exists(paths = paste0(path, i))==F) {
      dir.create(path = paste0(path,i), recursive = T)
    }
  }
  shell <- paste0(path,"src/run_d.sh")
  cat("#!/bin/bash\n", file = shell)
  prex <- c("KD-1","KD-2","WT-1","WT-2")
  lay1 <- paste0(path,"fq/",prex,"_1.fq.gz")
  lay2 <- paste0(path,"fq/",prex,"_2.fq.gz")
  cmd_01 <- paste0("#trim_galore --phred33 -j 4 ",
                   "-o ", path, "trim --paired ",lay1," ",lay2,"\n")
  lay1 <- paste0(path,"trim/",prex,"_1_val_1.fq.gz")
  lay2 <- paste0(path,"trim/",prex,"_2_val_2.fq.gz")
  cmd_02 <- paste0("STAR --runThreadN 10 --readFilesCommand zcat ",
                   "--outMultimapperOrder Random --outSAMtype BAM SortedByCoordinate ",
                   "--outFileNamePrefix ",paste0(path,"bam/",prex)," ",
                   "--outSAMprimaryFlag AllBestScore --twopassMode Basic ",
                   "--genomeDir ",inx1," --readFilesIn ",lay1," ",lay2,"\n")
  cmd_03 <- paste0("mv ",paste0(path,"bam/",prex,"Aligned.sortedByCoord.out.bam "),paste0(path,"bam/",prex,".bam"),"\n")
  cmd_04 <- paste0("samtools index ",paste0(path,"bam/",prex,".bam"),"\n")
  cmd_05 <- paste0("bamCoverage --normalizeUsing CPM -p 10 --minMappingQuality 1 --samFlagExclude 256 -of bigwig -bs 20 -b ",
                   paste0(path,"bam/",prex,".bam")," -o ",paste0(path,"bamCoverage/",prex,".bw"),"\n")
  cmd_06 <- paste0("TEtranscripts -t ",
                   paste(paste0(path,"bam/",prex[1:2],".bam"),collapse = " ")," -c ",
                   paste(paste0(path,"bam/",prex[3:4],".bam"),collapse = " "),
                   " --GTF ",inx3," --TE ",inx4," --sortByPos --mode multi ",
                   "--project Aef --outdir ",paste0(path,"count"),"\n")
  cmd_07 <- paste0("TElocal -b ",
                   paste0(path,"bam/",prex,".bam"),
                   " --GTF ",inx3," --TE ",inx5," --sortByPos --project ",
                   paste0(path, "count/",prex),"\n")
  cmd_aa <- paste0("TEtranscripts -t ",
                   paste(paste0(path,"bam/",prex[1:2],".bam"),collapse = " ")," -c ",
                   paste(paste0(path,"bam/",prex[3:4],".bam"),collapse = " "),
                   " --GTF ",inx3," --TE ",inx4," --sortByPos --mode uniq ",
                   "--project Aef --outdir ",paste0(path,"countUniq"),"\n")
  cmd_bb <- paste0("TElocal -b ",
                   paste0(path,"bam/",prex,".bam"),
                   " --GTF ",inx3," --TE ",inx5," --sortByPos --mode uniq ",
                   "--project ",paste0(path, "countUniq/",prex),"\n")
  cmd_08 <- paste0("stringtie ", paste0(path,"bam/",prex,".bam")," ",
                   "-p 20 -G ", inx3, " -o ", paste0(path,"stringtie/",prex,".gtf"),"\n")
  cmd_09 <- paste0("stringtie --merge -o ",
                   path,"stringtie/Aef.merge.gtf -G ",inx3," ",
                   paste(paste0(path,"stringtie/",prex,".gtf"), collapse = " "),"\n")
  cmd_10 <- paste0("cat ", paste0(path,"stringtie/Aef.merge.gtf "), 
                   "|sed 's/gene_name/gene_abcd/g; s/\\tgene_id/\\tgene_name/g' >",
                   paste0(path,"stringtie/Aef.merge.R.gtf"),"\n")
  cmd_11 <- paste0("TEtranscripts ",
                   "-t ",paste(paste0(path,"bam/",prex[1:2],".bam"),collapse = " ")," ",
                   "-c ",paste(paste0(path,"bam/",prex[3:4],".bam"),collapse = " ")," ",
                   "--GTF ",path,"stringtie/Aef.merge.R.gtf ",
                   "--TE ", inx4," --sortByPos --mode multi --project Aef ",
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
  for (i in cmd_aa) {cat(i, file = shell, append = T)}
  for (i in cmd_bb) {cat(i, file = shell, append = T)}
  #
  #two-pass mode是否有影响
  cmd_12 <- paste0("STAR --runThreadN 10 --readFilesCommand zcat ",
                   "--outMultimapperOrder Random --outSAMtype BAM SortedByCoordinate ",
                   "--outFileNamePrefix ",path,"dumpROOM/",prex," ",
                   "--outSAMprimaryFlag AllBestScore ",
                   "--genomeDir ",inx1," --readFilesIn ",lay1," ",lay2,"\n")
  cmd_13 <- paste0("mv ",
                   paste0(path,"dumpROOM/",prex,"Aligned.sortedByCoord.out.bam "),
                   paste0(path,"dumpROOM/",prex,".bam"),"\n")
  cmd_14 <- paste0("samtools index ",paste0(path,"dumpROOM/",prex,".bam"),"\n")
  cmd_15 <- paste0("TEtranscripts -t ",
                   paste(paste0(path,"dumpROOM/",prex[1:2],".bam"),collapse = " ")," -c ",
                   paste(paste0(path,"dumpROOM/",prex[3:4],".bam"),collapse = " "),
                   " --GTF ",inx3," --TE ",inx4," --sortByPos --mode multi ",
                   "--project Aef --outdir ",paste0(path,"dumpROOM"),"\n")
  for (i in cmd_12) {cat(i, file = shell, append = T)}
  for (i in cmd_13) {cat(i, file = shell, append = T)}
  for (i in cmd_14) {cat(i, file = shell, append = T)}
  for (i in cmd_15) {cat(i, file = shell, append = T)}
  print(paste0("nohup bash ",shell, " >",paste0(path,"Log/"),basename(shell),".log"," 2>&1 &"))
}
################################################################################
#0A、带和不带2-pass mode, 对差异倍数的影响
{
  #tmp1 <- paste0(path,"count/Aef_sigdiff_gene_TE.txt")
  #tmp1 <- read.table(tmp1)
  #tmp2 <- paste0(path,"dumpROOM/Aef_sigdiff_gene_TE.txt")
  #tmp2 <- read.table(tmp2)
  #tmp3 <- cbind(tmp1[,], tmp2[match(rownames(tmp1),rownames(tmp2)),])
  #tmp3 <- tmp3[,c(2,8)]
  #tmp3$mark <- tmp3$log2FoldChange/tmp3$log2FoldChange.1
  #tmp3 <- tmp3[order(tmp3$mark, decreasing = T),]
  #tmp3[1:10,]
  #ggplot(tmp3) +geom_point(aes(x = log2FoldChange, y = log2FoldChange.1)) +geom_abline(intercept = 0)
}
################################################################################
#04、DESeq2：Gene、TE、TE sites
{
  shy_diff <- function(mark, geneFile, siteFile, Class, VS){
    tmp1 <- read.table(geneFile, header = T, check.names = F, stringsAsFactors = F, row.names = 1)
    colnames(tmp1) <- sapply(str_split(basename(colnames(tmp1)),pattern = ".bam"), "[", 1)
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
    write.csv(geneAA, paste0(path,"DESeq2/",mark, ".geneOnly.csv"))
    upup <- subset(geneAA, log2FoldChange > log2(1.5) & padj < 0.05)
    down <- subset(geneAA, log2FoldChange < -log2(1.5)& padj < 0.05)
    write.csv(upup,paste0(path,"DESeq2/",mark, ".geneOnly.upup.csv"))
    write.csv(down,paste0(path,"DESeq2/",mark, ".geneOnly.down.csv"))
    #
    teOnly <- DESeq2::DESeqDataSetFromMatrix(countData = teOnly, colData = dfclas, design = ~Class)
    #teOnly <- teOnly[rowSums(BiocGenerics::counts(teOnly) > 2) >= colLen, ]
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
    #teGene <- teGene[rowSums(BiocGenerics::counts(teGene) > 2) >= colLen, ]
    teGene <- teGene[rowMeans(BiocGenerics::counts(teGene)) > 5,]
    teGene <- DESeq2::DESeq(teGene,)
    teGene <- DESeq2::results(teGene, contrast = VS)
    teGene <- teGene[grep(rownames(teGene), pattern = ":", invert = F),]
    write.csv(teGene, paste0(path,"DESeq2/",mark, ".te.GeneBackGround.csv"))
    upup <- subset(teGene, log2FoldChange > log2(1.5) & padj < 0.05)
    down <- subset(teGene, log2FoldChange < -log2(1.5)& padj < 0.05)
    write.csv(upup,paste0(path,"DESeq2/",mark, ".te.GeneBackGround.upup.csv"))
    write.csv(down,paste0(path,"DESeq2/",mark, ".te.GeneBackGround.down.csv"))
    #
    #
    tmp2 <- read.table(siteFile, header = T, check.names = F, stringsAsFactors = F, row.names = 1)
    colnames(tmp2) <- sapply(str_split(basename(colnames(tmp2)),pattern = ".bam"), "[", 1)
    colLen <- ceiling(length(colnames(tmp2))/2)
    #
    teGene <- tmp2
    teGene <- DESeq2::DESeqDataSetFromMatrix(countData = teGene, colData = dfclas, design = ~Class)
    #teGene <- teGene[rowSums(BiocGenerics::counts(teGene) > 2) >= colLen, ]
    teGene <- teGene[rowMeans(BiocGenerics::counts(teGene)) > 5,]
    teGene <- DESeq2::DESeq(teGene,)
    teGene <- DESeq2::results(teGene, contrast = VS)
    teGene <- teGene[grep(rownames(teGene), pattern = ":", invert = F),]
    write.csv(teGene, paste0(path,"DESeq2/",mark, ".te.site.GeneBackGround.csv"))
    upup <- subset(teGene, log2FoldChange > log2(1.5) & padj < 0.05)
    down <- subset(teGene, log2FoldChange < -log2(1.5)& padj < 0.05)
    write.csv(upup, paste0(path,"DESeq2/",mark, ".te.site.GeneBackGround.upup.csv"))
    write.csv(down, paste0(path,"DESeq2/",mark, ".te.site.GeneBackGround.down.csv"))
    #
    teGene <- tmp2[grep(rownames(tmp2), pattern = ":", invert = F),]
    teGene <- DESeq2::DESeqDataSetFromMatrix(countData = teGene, colData = dfclas, design = ~Class)
    #teGene <- teGene[rowSums(BiocGenerics::counts(teGene) > 2) >= colLen, ]
    teGene <- teGene[rowMeans(BiocGenerics::counts(teGene)) > 5,]
    teGene <- DESeq2::DESeq(teGene,)
    teGene <- DESeq2::results(teGene, contrast = VS)
    write.csv(teGene, paste0(path,"DESeq2/",mark, ".te.site.teOnly.csv"))
    upup <- subset(teGene, log2FoldChange > log2(1.5) & padj < 0.05)
    down <- subset(teGene, log2FoldChange < -log2(1.5)& padj < 0.05)
    write.csv(upup, paste0(path,"DESeq2/",mark, ".te.site.teOnly.upup.csv"))
    write.csv(down, paste0(path,"DESeq2/",mark, ".te.site.teOnly.down.csv"))
  }
  shy_diff(mark = "Aef", 
           geneFile = paste0(path, "count/Aef.cntTable"),
           siteFile = paste0(path, "count/Aef.site.cntTable"),
           Class = factor(c("KD","KD","WT","WT")), VS = c("Class","KD","WT"))
}
################################################################################
#05、CPM Excel
{
  prex <- c("KD-1","KD-2","WT-1","WT-2")
  cnt1 <- read.table(paste0(path,"count/Aef.cntTable"), header = T, sep = "\t", row.names = 1)
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
  openxlsx::write.xlsx(x = tmp1, file = paste0(path, "count/Aef.cpm.onlyTE.xlsx"))
  openxlsx::write.xlsx(x = tmp2, file = paste0(path, "count/Aef.cpm.onlyGene.xlsx"))
  openxlsx::write.xlsx(x = tmp3, file = paste0(path, "count/Aef.cpm.allGeneTE.xlsx"))
}
{
  prex <- c("KD-1","KD-2","WT-1","WT-2")
  cnt1 <- read.table(paste0(path,"countUniq/Aef.site.cntTable"), header = T, sep = "\t", row.names = 1)
  colnames(cnt1) <- prex
  colLen <- length(colnames(cnt1))
  cnt1 <- cnt1[rowSums(cnt1 == 0) < colLen, ]
  tmp1 <- cnt1[grep(x = rownames(cnt1), pattern=":", invert = F),]#TE
  tmp1 <- as.data.frame(apply(tmp1, 2, function(x){x/sum(x)*1000000}))
  tmp1$theName <- rownames(tmp1)
  openxlsx::write.xlsx(x = tmp1, file = paste0(path, "countUniq/Aef.cpm.onlyTE.site.xlsx"))
}
################################################################################
#06、CPM Point
{
  cnt1 <- read.table("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/Aef/count/Aef.cntTable", 
                     header = T, sep = "\t", row.names = 1)[,c(3,4,1,2)]
  colnames(cnt1) <- c("WT-1","WT-2","Aef-KD-1","Aef-KD-2")
  cnt1 <- cnt1[grep(x = rownames(cnt1), pattern=":", invert = F),]
  cnt1 <- as.data.frame(apply(cnt1, 2, function(x){x/sum(x)*1000000}))
  myta <- data.frame(WT = log2(rowMeans(cnt1[,c(1,2)]) +1),
                     KD = log2(rowMeans(cnt1[,c(3,4)]) +1))
  #
  #
  #
  #无基因背景
  {
    te01 <- read.csv(file = "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/Aef/DESeq2/Aef.teOnly.csv", row.names = 1)
    te01$mark <- NA
    te01$mark[which(te01$log2FoldChange > 0 & te01$padj < 0.05)] <- "upup"
    te01$mark[which(te01$log2FoldChange < 0 & te01$padj < 0.05)] <- "down"
    myta$mark <- "other"
    myta$mark[which(rownames(myta) %in% rownames(te01)[which(te01$mark=="upup")])] <- "upup"
    myta$mark[which(rownames(myta) %in% rownames(te01)[which(te01$mark=="down")])] <- "down"
    labe <- myta[grep(x = rownames(myta), pattern = "MERVL-int|MT2_Mm"),]
    myta$mark <- factor(x = myta$mark, levels = c("upup","down","other"))
    p_07 <- ggplot() +
      geom_point(data = myta[which(myta$mark =="other"),], aes(x = WT, y = KD), color = "gray") +
      geom_point(data = myta[which(myta$mark !="other"),], aes(x = WT, y = KD, color = mark)) +
      scale_color_manual(values = c("red","blue")) +
      geom_point(data = labe, aes(x = WT, y = KD),size=2,shape=1,color = "red") +
      geom_text_repel(data = labe, aes(x = WT, y = KD, label=rownames(labe))) +
      labs(x="Log2(CPM+1)(WT)",y="Log2(CPM+1)(KD)",
           title = "TE: Aef KD___WT (nono_geneBackground)") +
      guides(color = guide_legend(title = "Condition")) +
      theme_few() +
      theme(plot.margin = margin(unit = "cm",1,1,1,1),
            text = element_text(size = 12,family = "serif",color = "black"),
            plot.title = element_text(size = 16,family = "serif",color = "black"),
            legend.position = c(0.8,0.2),
            axis.title = element_text(size = 12,family = "serif",color = "black"),
            axis.text  = element_text(size = 12,family = "serif",color = "black")) +
      geom_abline(slope = 1,intercept = 0,lty="dashed")
    p_07
    ggsave(plot = p_07, 
           units = "cm", width = 12, height = 12,
           filename = paste0(path,"PLOT/","Aef.CPM.Point.TE.only.pdf"))
  }
  #
  #
  #
  #有基因背景
  {
    te02 <- read.csv("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/Aef/DESeq2/Aef.te.GeneBackGround.csv", 
                     row.names = 1)
    te02$mark <- NA
    te02$mark[which(te02$log2FoldChange > 0 & te02$padj < 0.05)] <- "upup"
    te02$mark[which(te02$log2FoldChange < 0 & te02$padj < 0.05)] <- "down"
    myta$mark <- "other"
    myta$mark[which(rownames(myta) %in% rownames(te02)[which(te02$mark=="upup")])] <- "upup"
    myta$mark[which(rownames(myta) %in% rownames(te02)[which(te02$mark=="down")])] <- "down"
    labe <- myta[grep(x = rownames(myta), pattern = "MERVL-int|MT2_Mm"),]
    myta$mark <- factor(x = myta$mark, levels = c("upup","down","other"))
    p_08 <- ggplot() +
      geom_point(data = myta[which(myta$mark =="other"),], aes(x = WT, y = KD), color = "gray") +
      geom_point(data = myta[which(myta$mark !="other"),], aes(x = WT, y = KD, color = mark)) +
      scale_color_manual(values = c("red","blue")) +
      geom_point(data = labe, aes(x = WT, y = KD),size=2,shape=1,color = "red") +
      geom_text_repel(data = labe, aes(x = WT, y = KD, label=rownames(labe))) +
      labs(x="Log2(CPM+1)(WT)",y="Log2(CPM+1)(KD)",
           title = "TE: Aef KD___WT (have_geneBackground)") +
      guides(color = guide_legend(title = "Condition")) +
      theme_few() +
      theme(plot.margin = margin(unit = "cm",1,1,1,1),
            text = element_text(size = 12,family = "serif",color = "black"),
            plot.title = element_text(size = 16,family = "serif",color = "black"),
            legend.position = c(0.8,0.2),
            axis.title = element_text(size = 12,family = "serif",color = "black"),
            axis.text  = element_text(size = 12,family = "serif",color = "black")) +
      geom_abline(slope = 1,intercept = 0,lty="dashed")
    p_08
    ggsave(plot = p_08, 
           units = "cm", width = 12, height = 12,
           filename = paste0(path,"PLOT/","Aef.CPM.Point.TE.geneBackground.pdf"))
  }
}
################################################################################
#06、Volcano for Genes
#----------------------------------出图了---------------------------------------
#----------------------------------出图了---------------------------------------
#----------------------------------出图了---------------------------------------
#----------------------------------出图了---------------------------------------
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
      scale_x_continuous(limits = c(-8.5,8.5), breaks = seq(-8,8,2)) +
      geom_hline(aes(yintercept=-log10(0.05)), linetype="dashed", colour="royalblue") +
      geom_vline(aes(xintercept=-log2(1.5)), linetype="dashed", colour="royalblue") +
      geom_vline(aes(xintercept=log2(1.5)), linetype="dashed", colour="royalblue") +
      theme_bw() + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      scale_color_manual(values = c("red","blue","black"),
                         labels = c(paste0(names(table(resdata$threshold)[1])," (",unname(table(resdata$threshold)[1]),")"),
                                    paste0(names(table(resdata$threshold)[2])," (",unname(table(resdata$threshold)[2]),")"),
                                    paste0(names(table(resdata$threshold)[3])," (",unname(table(resdata$threshold)[3]),")"))) +
      theme(legend.position = c(0.2,0.9), 
            legend.background = element_rect(fill = NA),
            legend.spacing.x = unit(0,"mm"),
            legend.spacing.y = unit(0,"mm"),
            legend.key = element_blank(),
            legend.text = element_text(size = 13),
            axis.title = element_text(size = 13)) +
      guides(color = guide_legend(override.aes = list(size = 3)))
    ggsave(plot = p_1, saveFile, units = "cm", width = 18, height = 16)
  }
  shy_volcano_gene(resFile  = paste0(path, "DESeq2/Aef.geneOnly.csv"),
                   saveFile = paste0(path, "PLOT/Aef.geneOnly.volcano.pdf"))
}
################################################################################
#07、Volcano for TEs
{
  shy_volcano_te <- function(resFile, saveFile){
    resdata <- as.data.frame(read.csv(resFile, header = T, stringsAsFactors = F))
    resdata <- na.omit(resdata)
    resdata$threshold[resdata$log2FoldChange > log2(1.5) &  resdata$padj <0.05] = "Upregulated TEs"
    resdata$threshold[resdata$log2FoldChange < -log2(1.5) & resdata$padj <0.05] = "Downregulated TEs"
    resdata$threshold[!(resdata$log2FoldChange > log2(1.5) & resdata$padj <0.05) & !(resdata$log2FoldChange < -log2(1.5) & resdata$padj < 0.05)] = "Other TEs"
    resdata$threshold <- factor(resdata$threshold, levels = c("Upregulated TEs","Downregulated TEs","Other TEs"))
    resdata[which(resdata$threshold=="Other TEs"),"X"] <- NA
    labb <- resdata
    if (unname(table(resdata$threshold)[1]) + unname(table(resdata$threshold)[2]) >20){
      labb <- labb[which(labb$padj<10^(-5)),]#or -10
    }
    p_2 <- ggplot(resdata, aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) +
      labs(col="", x="Log2FoldChange",y="-Log10(adjusted p-value)") + 
      geom_point(size=0.8) +
      scale_x_continuous(limits = c(-4,4)) +
      geom_hline(aes(yintercept=-log10(0.05)), linetype="dashed", colour="royalblue") +
      geom_vline(aes(xintercept=-log2(1.5)), linetype="dashed", colour="royalblue") +
      geom_vline(aes(xintercept=log2(1.5)), linetype="dashed", colour="royalblue") +#ylim(0,max(-log10(resdata$padj)[is.finite(-log10(resdata$padj))])+20) +
      theme_bw() + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      scale_color_manual(values = c("red","blue","black"),
                         labels = c(paste0(names(table(resdata$threshold)[1])," (",unname(table(resdata$threshold)[1]),")"),
                                    paste0(names(table(resdata$threshold)[2])," (",unname(table(resdata$threshold)[2]),")"),
                                    paste0(names(table(resdata$threshold)[3])," (",unname(table(resdata$threshold)[3]),")"))) +
      theme(legend.position = c(0.2,0.8), 
            legend.background = element_rect(fill = NA),
            legend.spacing.x = unit(0,"mm"),
            legend.spacing.y = unit(0,"mm"),
            legend.key = element_blank(),
            legend.text = element_text(size = 13),
            axis.title = element_text(size = 13)) +
      geom_text_repel(data = labb,
                      aes(x=log2FoldChange, y=-log10(padj), label=X), hjust=-0.125, size=4, colour="black", family="serif") +
      guides(color = guide_legend(override.aes = list(size = 3)))
    p_2
    ggsave(plot = p_2, saveFile, units = "cm", width = 18, height = 16)
  }
  shy_volcano_te(resFile  = paste0(path, "DESeq2/Aef.te.GeneBackGround.csv"),
                 saveFile = paste0(path, "PLOT/Aef.te.GeneBackGround.volcano.pdf"))
  shy_volcano_te(resFile  = paste0(path, "DESeq2/Aef.teOnly.csv"),
                 saveFile = paste0(path, "PLOT/Aef.teOnly.volcano.pdf"))
}
################################################################################
#08、MA
{
  shy_teMA <- function(difffile, savefile){
    rf01 <- read.csv(file = difffile, row.names = 1)
    rf01$subFam <- sapply(str_split(rownames(rf01),pattern = ":"),"[",1)
    rf01$family <- sapply(str_split(rownames(rf01),pattern = ":"),"[",2)
    rf01$classs <- sapply(str_split(rownames(rf01),pattern = ":"),"[",3)
    rf01$colo <- rf01$family
    rf01$colo[-(which(rf01$padj<0.05 & abs(rf01$log2FoldChange)>log2(1.5)))] <- "other"
    rf02 <- rf01
    rf02$colo <- factor(rf02$colo,
                        levels=c(unique(rf02$colo[grep(x = rf02$colo, pattern = "other", invert = T)]),"other"))
    #
    p_01 <- ggplot() +
      geom_point(data = rf02, 
                 aes(x=log2(baseMean+1), y=log2FoldChange, color=colo, shape=NULL),size=2.5) +
      guides(color = guide_legend(title = "repFamliy")) +
      labs(x = "log2(Expression in baseMean)", y = "log2Foldchange") +
      geom_hline(aes(yintercept =  log2(1.5)), linetype="dashed", colour="grey") +
      geom_hline(aes(yintercept = -log2(1.5)), linetype="dashed", colour="grey") +
      geom_text_repel(data = rf02[which(rf02$colo!="other"),], 
                      aes(x=log2(baseMean+1), y=log2FoldChange, label=subFam), 
                      hjust=-0.125, size=4, colour="black", family="serif") +
      scale_color_manual(values = c(ggsci::pal_d3(palette = "category20")(20)[1:c(length(unique(rf02$colo))-1)],"grey80")) +
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
  shy_teMA(difffile = paste0(path,"DESeq2/Aef.te.GeneBackGround.csv"),
           savefile = paste0(path,"PLOT/Aef.te.GeneBackGround.MA.pdf"))
  shy_teMA(difffile = paste0(path,"DESeq2/Aef.teOnly.csv"),
           savefile = paste0(path,"PLOT/Aef.teOnly.MA.pdf"))
}
################################################################################
#09、MA plot TE sites
{
  shy_repSite <- function(resdFile, saveFile){
    gf00 <- read.csv(file = resdFile, header = T, stringsAsFactors = F, row.names = 1)
    gf00 <- na.omit(gf00)
    gfaa <- gf00[grep(rownames(gf00), pattern = ".*MERVL-int.*"),]
    gfbb <- gf00[grep(rownames(gf00), pattern = ".*MT2_Mm.*"),]
    gfcc <- gf00[grep(rownames(gf00), pattern = ".*MERVL-int.*|.*MT2_Mm.*", invert = T),]
    #开始作图
    p_1 <- ggplot() +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", linewidth = 0.8),
            axis.line = element_line(), axis.text = element_text(family = "serif")) +
      ylab(label = expression(Log2Fold~Change)) +
      xlab(label = expression(Log2~(baseMean + 1))) +
      scale_color_manual(breaks = c("MERVL-int","MT2_Mm","Others"),values = c("red","purple","gray")) +
      geom_point(data=gfcc,aes(x=log2(baseMean+1),y=log2FoldChange,colour="Others"),size=1.2) +
      geom_point(data=gfbb,aes(x=log2(baseMean+1),y=log2FoldChange,colour="MT2_Mm"),size=1.2) +
      geom_point(data=gfaa,aes(x=log2(baseMean+1),y=log2FoldChange,colour="MERVL-int"),size=1.2) +
      theme(legend.position = c(0.7,0.9), 
            legend.background = element_rect(fill = NA),
            legend.spacing.x = unit(0,"mm"), legend.spacing.y = unit(0,"mm"),
            legend.key = element_blank(), legend.title = element_blank(),
            legend.text  = element_text(size = 13), axis.title   = element_text(size = 13)) +
      geom_hline(aes(yintercept=-log2(1.5)),   linetype="dashed",colour="royalblue") +
      geom_hline(aes(yintercept= log2(1.5)),   linetype="dashed",colour="royalblue") +
      guides(color = guide_legend(override.aes = list(size = 3),order=1))
    ggsave(plot = p_1, filename = saveFile, width = 15, height = 15, units = "cm")
  }
  shy_repSite(resdFile = paste0(path,"DESeq2/Aef.te.site.teOnly.csv"),
              saveFile = paste0(path,"PLOT/Aef.te.site.teOnly.mervl.pdf"))
}
################################################################################
#10、PCA
{
  #----------------------------------PCA-batch----------------------------------
  cnt1 <- read.table(paste0(path, "count/Aef.cntTable"),
                     header = T, check.names = F, stringsAsFactors = F, row.names = 1)
  cnt1 <- cnt1[grep(rownames(cnt1), pattern = ":", invert = T),]
  colnames(cnt1) <- sapply(str_split(basename(colnames(cnt1)), ".bam"),"[",1)
  #
  base::load(file = "/Reference/aaSHY/DatabaseFiles/forPCA.RData")
  myta <- cbind(cnt1, forPCA)
  #
  condition = sapply(str_split(colnames(myta),"-"), "[", 1)
  btch = c(rep("ours",length(colnames(cnt1))), condition[c(length(colnames(cnt1))+1):length(myta)])
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
            legend.position = c(0.8,0.75),
            axis.text    = element_text(family = "serif", size = 12),
            axis.title   = element_text(family = "serif", size = 12))
    print(ppp)
  }
  pdf(file = paste0(path,"PLOT/PCA.pdf"),width = 10, height = 8)
  shy_linshi(df = pd_1, i = 1)
  shy_linshi(df = pd_2, i = 2)
  dev.off()
}
################################################################################
#11、GSEA
#----------------------------------出图了---------------------------------------
#----------------------------------出图了---------------------------------------
#----------------------------------出图了---------------------------------------
#----------------------------------出图了---------------------------------------
{
  temp <- read.csv(paste0(path, "DESeq2/Aef.geneOnly.csv"), header = T)
  temp <- temp[, c("X", "log2FoldChange")]
  temp <- temp[order(temp$log2FoldChange, decreasing = T),]
  write.table(x = temp, 
              file = paste0(path, "DESeq2/Aef.res.rnk"), 
              quote = F, sep = "\t", col.names = F, row.names = F)
  cmd_01 <- paste0("gseapy prerank -g ",
                   "/Reference/aaSHY/GSEAgmt/TwoGenes.gmt -r ",
                   paste0(path, "DESeq2/Aef.res.rnk "),
                   "-p 10 -o ",
                   paste0(path, "PLOT/enrichGSEA/"))#conda activate base gseapy 0.10.8
  cmd_02 <- paste0("gseapy prerank -g ",
                   "/Reference/aaSHY/GSEAgmt/2022_majorZGA.gmt -r ",
                   paste0(path, "DESeq2/Aef.res.rnk "),
                   "--max-size 1000 -p 10 -o ",
                   paste0(path, "PLOT/enrichGSEA/"))#conda activate base gseapy 0.10.8
  print(cmd_01)
  print(cmd_02)
}
################################################################################
#12、Volcano [来自ghyuughyvvOE.R代码venn交集的59个基因]
{
  gggg
  tmp1 <- read.csv("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/Aef/DESeq2/Aef.geneOnly.csv")
  tmp1 <- tmp1[which(tmp1$X %in% gggg),]; #boxplot(tmp1$log2FoldChange)
  tmp1$threshold <- "Others"
  tmp1$threshold[tmp1$log2FoldChange > log2(1.5)  & tmp1$padj <0.05] = "Upregulated genes"
  tmp1$threshold[tmp1$log2FoldChange < -log2(1.5) & tmp1$padj <0.05] = "Downregulated genes"
  p_1 <- ggplot(tmp1, aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) +
    labs(col="", x="Log2FoldChange",y="-Log10(adjusted p-value)") + 
    geom_point(size=0.8) +
    scale_x_continuous(limits = c(-2, 4), breaks = seq(-2, 4, 1)) +
    geom_hline(aes(yintercept=-log10(0.05)), linetype="dashed", colour="royalblue") +
    geom_vline(aes(xintercept=-log2(1.5)), linetype="dashed", colour="royalblue") +
    geom_vline(aes(xintercept=log2(1.5)), linetype="dashed", colour="royalblue") +
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    scale_color_manual(values = c("red","blue","black"),
                       labels = c(paste0(names(table(tmp1$threshold)[1])," (",unname(table(tmp1$threshold)[1]),")"),
                                  paste0(names(table(tmp1$threshold)[2])," (",unname(table(tmp1$threshold)[2]),")"),
                                  paste0(names(table(tmp1$threshold)[3])," (",unname(table(tmp1$threshold)[3]),")"))) +
    theme(legend.position = c(0.2,0.9), 
          legend.background = element_rect(fill = NA),
          legend.spacing.x = unit(0,"mm"),
          legend.spacing.y = unit(0,"mm"),
          legend.key = element_blank(),
          legend.text = element_text(size = 13),
          axis.title = element_text(size = 13)) +
    guides(color = guide_legend(override.aes = list(size = 3)))
  p_1
}
################################################################################
#13、pheatmap
{
  shy_pheatmap <- function(readsFile, gene_resFile, tete_resFile, savePath){
    print("reads counted by TEtranscripts")
    #
    data <- read.table(readsFile, header = T, row.names = 1, check.names = F, stringsAsFactors = F)
    data <- data[,c(3,4,1,2)]
    colnames(data) <- sapply(str_split(basename(colnames(data)), pattern = ".bam"), "[", 1)
    gene <- data[grep(rownames(data), pattern = ":", invert = T),]
    gene <- gene[rowSums(gene) > 5,]
    gene <- as.data.frame(apply(gene, 2, function(x){x/sum(x)*1000000}))
    tete <- data[grep(rownames(data), pattern = ":", invert = F),]
    tete <- tete[rowSums(tete) > 5,]
    tete <- as.data.frame(apply(tete, 2, function(x){x/sum(x)*1000000}))
    #
    resG <- as.data.frame(read.csv(gene_resFile, header = T, stringsAsFactors = F))
    resG <- na.omit(resG)
    sigG <- resG[which(abs(resG$log2FoldChange) >  log2(1.5) & resG$padj<0.05), "X"]
    sigG_upup <- resG[which(resG$log2FoldChange >  log2(1.5) & resG$padj<0.05), "X"]
    sigG_down <- resG[which(resG$log2FoldChange < -log2(1.5) & resG$padj<0.05), "X"]
    resG <- as.data.frame(read.csv(tete_resFile, header = T, stringsAsFactors = F))
    resG <- na.omit(resG)
    sigT <- resG[which(abs(resG$log2FoldChange) >  log2(1.5) & resG$padj<0.05), "X"]
    sigT_upup <- resG[which(resG$log2FoldChange >  log2(1.5) & resG$padj<0.05), "X"]
    sigT_down <- resG[which(resG$log2FoldChange < -log2(1.5) & resG$padj<0.05), "X"]
    ##开始画图
    pdf(paste0(savePath,"Aef.pheatmap.v1.pdf"), width = 8, height = 10)
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
             show_rownames = T, show_colnames = T, cluster_row = T, cluster_cols = T,
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
  shy_pheatmap(readsFile = paste0(path, "count/Aef.cntTable"),
               gene_resFile = paste0(path, "DESeq2/Aef.geneOnly.csv"),
               tete_resFile = paste0(path, "DESeq2/Aef.teOnly.csv"),
               savePath = paste0(path, "PLOT/pheatmap/"))
}
################################################################################
#14、pheatmap
####fold change >0 && pvalue <0.05
{
  shy_pheatmap <- function(readsFile, tete_resFile, tete_resBack, savePath){
    source("~/shyScripts/rnaSEQ/annoPheatmap.R")
    data <- read.table(readsFile, header = T, row.names = 1, check.names = F, stringsAsFactors = F)
    data <- data[,c(3,4,1,2)]
    colnames(data) <- c("WT-1","WT-2","Aef-KD-1","Aef-KD-2")
    #
    #没背景
    te01 <- data[grep(rownames(data), pattern = ":", invert = F),]
    #te01 <- te01[rowSums(te01) > 5,]
    te01 <- as.data.frame(apply(te01, 2, function(x){x/sum(x)*1000000}))
    #有背景
    te02 <- data
    #te02 <- te02[rowSums(te02) > 5,]
    te02 <- as.data.frame(apply(te02, 2, function(x){x/sum(x)*1000000}))
    te02 <- te02[grep(rownames(te02), pattern = ":", invert = F),]
    #
    res1 <- as.data.frame(read.csv(tete_resFile, header = T, stringsAsFactors = F))
    res1 <- na.omit(res1)
    #sig1 <- res1[which(abs(res1$log2FoldChange) >  log2(1.5) & res1$padj < 0.05), "X"]
    sig1 <- res1[which(abs(res1$log2FoldChange) >  0 & res1$pvalue < 0.05), "X"]
    res2 <- as.data.frame(read.csv(tete_resBack, header = T, stringsAsFactors = F))
    res2 <- na.omit(res2)
    #sig2 <- res2[which(abs(res2$log2FoldChange) >  log2(1.5) & res2$padj < 0.05), "X"]
    sig2 <- res2[which(abs(res2$log2FoldChange) >  0 & res2$pvalue < 0.05), "X"]
    #
    pdf(paste0(savePath,"Aef.pheatmap.TE-v2.pdf"), width = 8, height = 10)
    #
    gggg <- te01[which(rownames(te01) %in% sig1),]
    rownames(gggg) <- sapply(str_split(rownames(gggg),":"),"[",1)
    pheatmap(mat = gggg,
             main="\n\nsignificant TEs (nono_GeneBackground)",
             scale = "row", cellwidth = 60, 
             show_rownames = T, show_colnames = T, cluster_row = T, cluster_cols = F,
             labels_col = colnames(gggg), angle_col = "45",
             breaks = seq(-2,2,by=0.01), border = F,
             color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                        colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
    gggg <- te02[which(rownames(te02) %in% sig2),]
    rownames(gggg) <- sapply(str_split(rownames(gggg),":"),"[",1)
    pheatmap(mat = gggg,
             main="\n\nsignificant TEs (have_GeneBackground)",
             scale = "row", cellwidth = 60, 
             show_rownames = T, show_colnames = T, cluster_row = T, cluster_cols = F,
             labels_col = colnames(gggg), angle_col = "45",
             breaks = seq(-2,2,by=0.01), border = F,
             color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                        colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
    
    dev.off()
  }
  shy_pheatmap(readsFile = paste0(path, "count/Aef.cntTable"),
               tete_resFile = paste0(path, "DESeq2/Aef.teOnly.csv"),
               tete_resBack = paste0(path, "DESeq2/Aef.te.GeneBackGround.csv"),
               savePath = paste0(path, "PLOT/pheatmap/"))
}
####fold change >0 && ERVL
{
  shy_pheatmap <- function(readsFile, tete_resFile, tete_resBack, savePath){
    source("~/shyScripts/rnaSEQ/annoPheatmap.R")
    data <- read.table(readsFile, header = T, row.names = 1, check.names = F, stringsAsFactors = F)
    data <- data[,c(3,4,1,2)]
    colnames(data) <- c("WT-1","WT-2","Aef-KD-1","Aef-KD-2")
    #
    #没背景
    te01 <- data[grep(rownames(data), pattern = ":", invert = F),]
    #te01 <- te01[rowSums(te01) > 5,]
    te01 <- as.data.frame(apply(te01, 2, function(x){x/sum(x)*1000000}))
    te01 <- te01[grep("ERVL",rownames(te01)),]
    #有背景
    te02 <- data
    #te02 <- te02[rowSums(te02) > 5,]
    te02 <- as.data.frame(apply(te02, 2, function(x){x/sum(x)*1000000}))
    te02 <- te02[grep(rownames(te02), pattern = ":", invert = F),]
    te02 <- te02[grep("ERVL",rownames(te02)),]
    #
    res1 <- as.data.frame(read.csv(tete_resFile, header = T, stringsAsFactors = F))
    res1 <- na.omit(res1)
    #sig1 <- res1[which(abs(res1$log2FoldChange) >  log2(1.5) & res1$padj < 0.05), "X"]
    sig1 <- res1[which(abs(res1$log2FoldChange) >  0), "X"]
    res2 <- as.data.frame(read.csv(tete_resBack, header = T, stringsAsFactors = F))
    res2 <- na.omit(res2)
    #sig2 <- res2[which(abs(res2$log2FoldChange) >  log2(1.5) & res2$padj < 0.05), "X"]
    sig2 <- res2[which(abs(res2$log2FoldChange) >  0), "X"]
    #
    pdf(paste0(savePath,"Aef.pheatmap.TE-v2.pdf"), width = 8, height = 10)
    #
    gggg <- te01[which(rownames(te01) %in% sig1),]
    rownames(gggg) <- sapply(str_split(rownames(gggg),":"),"[",1)
    inds_merv <- c("MERVL-int","MT2_Mm")
    p_01 <- pheatmap(mat = gggg,
                     main="\n\nsignificant TEs (nono_GeneBackground)",
                     scale = "row", cellwidth = 60,
                     show_rownames = T, show_colnames = T, cluster_row = T, cluster_cols = F,
                     right_annotation = labs,
                     labels_col = colnames(gggg), angle_col = "45",
                     breaks = seq(-2,2,by=0.01), border = F,
                     color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))*(2/4)),
                                colorRampPalette(colors = c("white","red"))(length(seq(-2, 2,by=0.01))*(2/4))))
    add.flag(p_01,
             kept.labels = inds_merv,
             repel.degree = 0.2)
    #
    #
    #
    gggg <- te02[which(rownames(te02) %in% sig2),]
    rownames(gggg) <- sapply(str_split(rownames(gggg),":"),"[",1)
    inds_merv <- c("MERVL-int","MT2_Mm")
    p_02 <- pheatmap(mat = gggg,
                     main="\n\nsignificant TEs (have_GeneBackground)",
                     scale = "row", cellwidth = 60, 
                     show_rownames = T, show_colnames = T, cluster_row = T, cluster_cols = F,
                     labels_col = colnames(gggg), angle_col = "45",
                     breaks = seq(-2,2,by=0.01), border = F,
                     color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                                colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
    add.flag(p_02,
             kept.labels = inds_merv,
             repel.degree = 0.2)
    
    dev.off()
  }
  shy_pheatmap(readsFile = paste0(path, "count/Aef.cntTable"),
               tete_resFile = paste0(path, "DESeq2/Aef.teOnly.csv"),
               tete_resBack = paste0(path, "DESeq2/Aef.te.GeneBackGround.csv"),
               savePath = paste0(path, "PLOT/pheatmap/"))
}
####fold change >0 && ERVL && set WT as 1
{
  shy_pheatmap <- function(readsFile, tete_resFile, tete_resBack, savePath){
    source("~/shyScripts/rnaSEQ/annoPheatmap.R")
    data <- read.table(readsFile, header = T, row.names = 1, check.names = F, stringsAsFactors = F)
    data <- data[,c(3,4,1,2)]
    colnames(data) <- c("WT-1","WT-2","Aef-KD-1","Aef-KD-2")
    #
    #没背景
    te01 <- data[grep(rownames(data), pattern = ":", invert = F),]
    #te01 <- te01[rowSums(te01) > 5,]
    te01 <- as.data.frame(apply(te01, 2, function(x){x/sum(x)*1000000}))
    te01 <- te01[grep("ERVL",rownames(te01)),]
    #有背景
    te02 <- data
    #te02 <- te02[rowSums(te02) > 5,]
    te02 <- as.data.frame(apply(te02, 2, function(x){x/sum(x)*1000000}))
    te02 <- te02[grep(rownames(te02), pattern = ":", invert = F),]
    te02 <- te02[grep("ERVL",rownames(te02)),]
    #
    res1 <- as.data.frame(read.csv(tete_resFile, header = T, stringsAsFactors = F))
    res1 <- na.omit(res1)
    #sig1 <- res1[which(abs(res1$log2FoldChange) >  log2(1.5) & res1$padj < 0.05), "X"]
    sig1 <- res1[which(abs(res1$log2FoldChange) >  0), "X"]
    res2 <- as.data.frame(read.csv(tete_resBack, header = T, stringsAsFactors = F))
    res2 <- na.omit(res2)
    #sig2 <- res2[which(abs(res2$log2FoldChange) >  log2(1.5) & res2$padj < 0.05), "X"]
    sig2 <- res2[which(abs(res2$log2FoldChange) >  0), "X"]
    #
    pdf(paste0(savePath,"Aef.pheatmap.TE-v3.pdf"), width = 8, height = 10)
    #
    gggg <- te01[which(rownames(te01) %in% sig1),]
    rownames(gggg) <- sapply(str_split(rownames(gggg),":"),"[",1)
    gggg <- gggg +0.001
    gggg <- as.data.frame(apply(gggg,2,function(x){x/rowMeans(gggg[,c(1,2)])}))
    gggg <- log2(gggg +0.001)
    gggg <- gggg[which(rowSums(is.na(gggg)) <1),]
    gggg <- gggg[grep("MT2|MERVL",rownames(gggg)),]
    inds_merv <- c("MERVL-int","MT2_Mm")
    p_01 <- pheatmap(mat = gggg,
                     main="\n\nERVL (foldChange >0 & nono_GeneBackground)",
                     scale = "none", cellwidth = 40,
                     show_rownames = T, show_colnames = T, cluster_row = T, cluster_cols = F,
                     right_annotation = labs,
                     labels_col = colnames(gggg), angle_col = "45",
                     breaks = seq(-1,1,by=0.01), border = F,
                     color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-1,1,by=0.01))*(2/4)),
                                colorRampPalette(colors = c("white","red"))(length(seq(-1, 1,by=0.01))*(2/4))))
    add.flag(p_01,
             kept.labels = inds_merv,
             repel.degree = 0.2)
    #
    #
    #
    gggg <- te02[which(rownames(te02) %in% sig2),]
    rownames(gggg) <- sapply(str_split(rownames(gggg),":"),"[",1)
    gggg <- gggg +0.001
    gggg <- as.data.frame(apply(gggg,2,function(x){x/rowMeans(gggg[,c(1,2)])}))
    gggg <- log2(gggg +0.001)
    gggg <- gggg[which(rowSums(is.na(gggg)) <1),]
    #gggg <- gggg[grep("MT2|MERVL",rownames(gggg)),]
    inds_merv <- c("MERVL-int","MT2_Mm")
    p_02 <- pheatmap(mat = gggg,
                     main="\n\nERVL (foldChange >0 & have_GeneBackground)",
                     scale = "none", cellwidth = 40, 
                     show_rownames = T, show_colnames = T, cluster_row = T, cluster_cols = F,
                     labels_col = colnames(gggg), angle_col = "45",
                     breaks = seq(-1,1,by=0.01), border = F,
                     color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-1,1,by=0.01))/2),
                                colorRampPalette(colors = c("white","red"))(length(seq(-1,1,by=0.01))/2)))
    add.flag(p_02,
             kept.labels = inds_merv,
             repel.degree = 0.2)
    dev.off()
  }
  shy_pheatmap(readsFile = paste0(path, "count/Aef.cntTable"),
               tete_resFile = paste0(path, "DESeq2/Aef.teOnly.csv"),
               tete_resBack = paste0(path, "DESeq2/Aef.te.GeneBackGround.csv"),
               savePath = paste0(path, "PLOT/pheatmap/"))
}
#
#
#
#
#
#
###ghyvv OE & ghyuu KO共同上调的122个基因在Aef的表达热图








################################################################################
#15、画MERVL-int sites的CPM的boxplot
#----------------------------------出图了---------------------------------------
#----------------------------------出图了---------------------------------------
#----------------------------------出图了---------------------------------------
#----------------------------------出图了---------------------------------------
{
  mervs <- openxlsx::read.xlsx(paste0(path,"countUniq/Aef.cpm.onlyTE.site.xlsx"))
  mervl <- mervs[grep("MERVL-int",mervs$theName),]
  mervl <- mervl[rowSums(mervl[,1:4]==0) < 4,]
  mervl <- tidyr::pivot_longer(data = mervl, cols = -theName, names_to = "group", values_to = "cpms")
  vaue <- wilcox.test(x = mervl$cpms[grep("WT",mervl$group)],
                      y = mervl$cpms[grep("KD",mervl$group)], alternative = "less")
  print(vaue$p.value)
  #
  mervl$marks <- sapply(str_split(mervl$group,"-"),"[",1)
  mervl$marks <- factor(mervl$marks, levels = rev(unique(mervl$marks)))
  p_001 <- ggplot(data = mervl,aes(x = marks, y = log2(cpms+0.1))) +
    geom_boxplot(aes(fill = marks)) +
    labs(title = paste0("MERVL-int site Expr"),
         y = "Log2 (CPM+0.1)", x = "", fill="Group") +
    geom_signif(comparisons = list(c("WT", "KD")),
                map_signif_level=T, family = "serif",
                textsize=6, test = wilcox.test, step_increase = 0) +
    stat_summary(geom = "point", fun = "median",shape = 22,size = 2,
                 position = position_dodge(width = 1)) +
    theme_classic() +
    scale_fill_manual(values = c("#3FA0C0","#F0A19A"),
                      labels = c("WT", "Aef-KD")) +
    theme(plot.margin = margin(unit = "cm",1,1,1,1),
          legend.text = element_text(family = "serif", size = 12),
          legend.title = element_text(family = "serif", size = 12),
          axis.text.x = element_text(angle = 45,hjust = 1,family = "serif",size = 12),
          axis.text.y = element_text(family = "serif", size = 12)) +
    scale_y_continuous(limits = c(-5,10)) +
    #scale_y_continuous(limits = c(0,10),breaks = c(0,1,2,3,4,5)) +
    scale_x_discrete(labels = c("WT", "Aef-KD"))
  print(p_001)
  ggsave(plot = p_001,
         filename = paste0(path,"PLOT/siteTE/","Aef.uniq.mervl-int.boxplot.pdf"),
         units = "cm", width = 20, height = 15)
}
#
#
#
#
#
#
#
################################################################################
#
#
#
#
#
#
#
#
{
  Class = factor(c("WT","WT","KD","KD")); VS = c("Class","KD","WT")
  #
  #data <- read.csv("/disk5/fx/RNA-seq/Aef_RNAseq/Count/batch/TE_Aef.csv", row.names = 1)
  #data <- read.table("/disk5/fx/RNA-seq/Aef_RNAseq/Count/hisat2/mm10_repname_M.saf.cut",
  #                   header = T, row.names = 1, skip = "#", stringsAsFactors = F,sep = "\t")
  #data <- read.table("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/Aef/countUniq/Aef.cntTable",
  #                   header = T, row.names = 1, skip = "#",stringsAsFactors = F,sep = "\t")
  data <- read.table("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/Aef/countTry/Aef.cntTable",
                     header = T, row.names = 1, skip = "#",stringsAsFactors = F,sep = "\t")
  tmp1 <- data[grep(":",rownames(data)),]
  #gggg <- data[grep(":",rownames(data),invert = T),]
  #tttt <- data[grep("MERVL-int|MT2_Mm",rownames(data)),]
  #tmp1 <- rbind(gggg,tttt)
  dfclas <- data.frame(row.names = colnames(tmp1), Class, stringsAsFactors = T)
  #
  teOnly <- tmp1
  teOnly <- DESeq2::DESeqDataSetFromMatrix(countData = teOnly, colData = dfclas, design = ~Class)
  teOnly <- teOnly[rowMeans(BiocGenerics::counts(teOnly)) > 5,]
  teOnly <- DESeq2::DESeq(teOnly,)
  teOnly <- DESeq2::results(teOnly, contrast = VS)
  teOnly <- as.data.frame(teOnly)
  teOnly[grep("MERVL",rownames(teOnly)),]
  teOnly[which(teOnly$log2FoldChange>0 & teOnly$pvalue<0.05),1]
}
################################################################################
#
#
#
#
#
#
#
#special genes在Aef DESeq2结果中的分布 [2、3有出图]-------------------------------------------
{
  red1 <- read.csv("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/Aef/DESeq2/Aef.geneOnly.csv")
  gg_1 <- read.csv(file = "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuu/DESeq2/ghyuu.geneOnly.upup.csv")
  gg_1 <- gg_1$X
  gg_2 <- read.csv(file = "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuu/DESeq2/ghyuu.geneOnly.down.csv")
  gg_2 <- gg_2$X
  red1$mark <- "upup"
  red2 <- red1[which(red1$X %in% c(gg_1,gg_2)),]
  red2$mark[which(red2$X %in% gg_2)] <- "down"
  #
  #
  red3 <- red2
  ggplot(data = red3) +
    geom_density(aes(x = log2FoldChange, color = mark))
  #
  #
  ggplot(data = red3) +
    stat_ecdf(aes(x = log2FoldChange, color = mark))
}
{
  red1 <- read.csv("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/Aef/DESeq2/Aef.geneOnly.csv")
  gg_1 <- read.csv(file = "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuughyvv/DESeq2/ghyvv.geneOnly.upup.csv")
  gg_1 <- gg_1$X
  gg_2 <- read.csv(file = "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuughyvv/DESeq2/ghyvv.geneOnly.down.csv")
  gg_2 <- gg_2$X
  red1$mark <- "upup"
  red2 <- red1[which(red1$X %in% c(gg_1,gg_2)),]
  red2$mark[which(red2$X %in% gg_2)] <- "down"
  #
  #
  red3 <- red2
  ###ggplot(data = red3) +
  ###  geom_density(aes(x = log2FoldChange, color = mark))
  #
  #
  red3$mark <- factor(red3$mark, levels = c("upup","down"))
  p_001 <- ggplot(data = red3) +
    stat_ecdf(aes(x = log2FoldChange, color = mark), size = 1) +
    theme_few() +
    labs(title = paste0("title--------------------")) +
    scale_x_continuous(limits = c(floor(min(red3$log2FoldChange)), ceiling(max(red3$log2FoldChange)))) +
    scale_color_manual(values = c("indianred3","steelblue3")) +
    theme(plot.margin = margin(unit = "cm", 1,1,1,1),
          text = element_text(size = 12,family = "serif",color = "black"),
          axis.title = element_blank(),
          axis.text  = element_text(size = 12,family = "serif",color = "black"))
  ggsave(plot = p_001, 
         units = "cm", width = 20, height = 16,
         filename = "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/Aef/PLOT/cumulative.ghyvvOE_v1.pdf")
}
{
  red1 <- read.csv("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/Aef/DESeq2/Aef.geneOnly.csv")
  gg_1 <- read.csv(file = "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/thesugas/DESeq2/sugas.group1.geneOnly.upup.csv")
  gg_1 <- gg_1$X
  gg_2 <- read.csv(file = "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/thesugas/DESeq2/sugas.group1.geneOnly.down.csv")
  gg_2 <- gg_2$X
  red1$mark <- "upup"
  red2 <- red1[which(red1$X %in% c(gg_1,gg_2)),]
  red2$mark[which(red2$X %in% gg_2)] <- "down"
  #
  #
  red3 <- red2
  ###ggplot(data = red3) +
  ###  geom_density(aes(x = log2FoldChange, color = mark))
  #
  #
  red3$mark <- factor(red3$mark, levels = c("upup","down"))
  p_001 <- ggplot(data = red3) +
    stat_ecdf(aes(x = log2FoldChange, color = mark), size = 1) +
  theme_few() +
    labs(title = paste0("title--------------------")) +
    scale_x_continuous(limits = c(floor(min(red3$log2FoldChange)), ceiling(max(red3$log2FoldChange)))) +
    scale_color_manual(values = c("indianred3","steelblue3")) +
    theme(plot.margin = margin(unit = "cm", 1,1,1,1),
          text = element_text(size = 12,family = "serif",color = "black"),
          axis.title = element_blank(),
          axis.text  = element_text(size = 12,family = "serif",color = "black"))
  ggsave(plot = p_001, 
         units = "cm", width = 20, height = 16,
         filename = "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/Aef/PLOT/cumulative.group1_v1.pdf")
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
#
#
################################################################################
{
  mervs <- openxlsx::read.xlsx(paste0(path,"countUniq/Aef.cpm.onlyTE.site.xlsx"))
  mervl <- mervs[grep("MERVL-int",mervs$theName),]
  ### mervl <- mervl[rowSums(mervl[,1:4]==0) < 4,]
  mervl <- mervl[rowMeans(mervl[,1:4]) > 5,]
  mervl <- tidyr::pivot_longer(data = mervl, cols = -theName, names_to = "group", values_to = "cpms")
  vaue <- wilcox.test(x = mervl$cpms[grep("WT",mervl$group)],
                      y = mervl$cpms[grep("KD",mervl$group)], alternative = "less")
  print(vaue$p.value)
  #
  mervl$marks <- sapply(str_split(mervl$group,"-"),"[",1)
  mervl$marks <- factor(mervl$marks, levels = rev(unique(mervl$marks)))
  p_001 <- ggplot(data = mervl,aes(x = marks, y = log2(cpms+0.1))) +
    geom_boxplot(aes(fill = marks)) +
    labs(title = paste0("MERVL-int site Expr"),
         y = "Log2 (CPM+0.1)", x = "", fill="Group") +
    geom_signif(comparisons = list(c("WT", "KD")),
                map_signif_level=T, family = "serif",
                textsize=6, test = wilcox.test, step_increase = 0) +
    stat_summary(geom = "point", fun = "median",shape = 22,size = 2,
                 position = position_dodge(width = 1)) +
    theme_classic() +
    scale_fill_manual(values = c("#3FA0C0","#F0A19A"),
                      labels = c("WT", "Aef-KD")) +
    theme(plot.margin = margin(unit = "cm",1,1,1,1),
          legend.text = element_text(family = "serif", size = 12),
          legend.title = element_text(family = "serif", size = 12),
          axis.text.x = element_text(angle = 45,hjust = 1,family = "serif",size = 12),
          axis.text.y = element_text(family = "serif", size = 12)) +
    scale_y_continuous(limits = c(-5,10)) +
    #scale_y_continuous(limits = c(0,10),breaks = c(0,1,2,3,4,5)) +
    scale_x_discrete(labels = c("WT", "Aef-KD"))
  print(p_001)
  ggsave(plot = p_001,
         filename = paste0(path,"PLOT/siteTE/","Aef.uniq.mervl-int.boxplot-revNsmb.pdf"),
         units = "cm", width = 15, height = 15)
  #######
  #
  #
  #
  #######
  mervl[1,]
  sed3 <- mervl %>%
    tidyr::pivot_wider(names_from = group, values_from = cpms, id_cols = theName)
  sed3 <- sed3[,c(1,4,5,2,3)]
  sed4 <- as.data.frame(sed3[,2:5])
  rownames(sed4) <- sed3$theName
  pdf(file = paste0(path,"PLOT/siteTE/","Aef.uniq.mervl-int.pheatmap-revNsmb.pdf"),
      height = 6, width = 6, family = "serif")
  pheatmap(mat = sed4,
           main="\n-------Aef KD---------",
           scale = "row", cellwidth = 20, 
           show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = F,
           labels_col = colnames(sed4), angle_col = "45",
           breaks = seq(-2,2,by=0.01), border = F,
           color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                      colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
  dev.off()
}










