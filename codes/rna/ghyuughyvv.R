#write by Sun Haiayng at 2024.08.08
#01、搭环境
{
  rm(list=ls())
  for (i in c("data.table","ggridges","ggrepel","scales","stringr",
              "Biostrings","tidyr","pheatmap","enrichplot","DESeq2","patchwork",
              "ggpubr","ggsci","ggtranscript","rtracklayer","Mfuzz","BiocManager",
              "dplyr","IRanges","GenomicRanges","ReactomePA","ggplot2","clusterProfiler",
              "topGO","Rgraphviz","pathview","org.Mm.eg.db","KEGG.db")){library(i, character.only = T,quietly = T)}; rm(i)
}
#02、给索引
{
  path = "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuughyvv/"
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
  shell <- paste0(path,"src/run_a.sh")
  cat("#!/bin/bash\n", file = shell)
  prex <- sapply(str_split(list.files((paste0(path,"fq/")),pattern = "*1.fq.gz$*"),"_1.fq.gz"), "[",1)
  lay1 <- paste0(path,"fq/",prex,"_1.fq.gz")
  lay2 <- paste0(path,"fq/",prex,"_2.fq.gz")
  cmd_01 <- paste0("trim_galore --phred33 -j 8 --clip_R1 10 --clip_R2 10 ",
                   "-a 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC' ",
                   "-a2 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT' ",
                   "--fastqc -o ", path, "trim --paired ",lay1," ",lay2,"\n")
  lay1 <- paste0(path,"trim/",prex,"_1_val_1.fq.gz")
  lay2 <- paste0(path,"trim/",prex,"_2_val_2.fq.gz")
  cmd_03 <- paste0("STAR --runThreadN 10 --readFilesCommand zcat ",
                   "--outMultimapperOrder Random --outSAMtype BAM SortedByCoordinate ",
                   "--outFileNamePrefix ",paste0(path,"bam/",prex)," ",
                   "--outSAMprimaryFlag AllBestScore --twopassMode Basic ",
                   "--genomeDir ",inx1," --readFilesIn ",lay1," ",lay2,"\n")
  cmd_04 <- paste0("mv ",paste0(path,"bam/",prex,"Aligned.sortedByCoord.out.bam "),paste0(path,"bam/",prex,".bam"),"\n")
  cmd_05 <- paste0("samtools index ",paste0(path,"bam/",prex,".bam"),"\n")
  cmd_06 <- paste0("bamCoverage --normalizeUsing CPM -p 10 --minMappingQuality 1 ",
                   "--samFlagExclude 256 -of bigwig -bs 20 -b ",path,
                   "bam/",prex,".bam -o ",path,"bamCoverage/",prex,".bw","\n")
  cmd_07 <- paste0("TEtranscripts -t ",
                   paste(paste0(path,"bam/",prex[1:4],".bam"),collapse = " ")," -c ",
                   paste(paste0(path,"bam/",prex[5:6],".bam"),collapse = " "),
                   " --GTF ",inx3," --TE ",inx4," --sortByPos --mode multi ",
                   "--project ghyuughyvv --outdir ",paste0(path,"count"),"\n")
  cmd_08 <- paste0("TElocal -b ",
                   paste0(path,"bam/",prex,".bam"),
                   " --GTF ",inx3," --TE ",inx5," --sortByPos --project ",
                   paste0(path, "count/",prex),"\n")
  cmd_09 <- paste0("stringtie ", paste0(path,"bam/",prex,".bam")," ",
                   "-p 10 -G ", inx3, " -o ", paste0(path,"stringtie/",prex,".gtf"),"\n")
  cmd_10 <- paste0("stringtie --merge -o ",
                   path,"stringtie/ghyuughyvv.merge.gtf -G ",inx3," ",
                   paste(paste0(path,"stringtie/",prex,".gtf"), collapse = " "),"\n")
  cmd_11 <- paste0("cat ", paste0(path,"stringtie/ghyuughyvv.merge.gtf "), 
                   "|sed 's/gene_name/gene_abcd/g; s/\\tgene_id/\\tgene_name/g' >",
                   paste0(path,"stringtie/ghyuughyvv.merge.R.gtf"),"\n")
  cmd_12 <- paste0("TEtranscripts ",
                   "-t ",paste(paste0(path,"bam/",prex[1:4],".bam"),collapse = " ")," ",
                   "-c ",paste(paste0(path,"bam/",prex[5:6],".bam"),collapse = " ")," ",
                   "--GTF ",path,"stringtie/ghyuughyvv.merge.R.gtf ",
                   "--TE ", inx4," --sortByPos --mode multi --project ghyuughyvv ",
                   "--outdir ",path,"stringtie/count","\n")
  for (i in cmd_01) {cat(i, file = shell, append = T)}
  for (i in cmd_03) {cat(i, file = shell, append = T)}
  for (i in cmd_04) {cat(i, file = shell, append = T)}
  for (i in cmd_05) {cat(i, file = shell, append = T)}
  for (i in cmd_06) {cat(i, file = shell, append = T)}
  for (i in cmd_07) {cat(i, file = shell, append = T)}
  for (i in cmd_08) {cat(i, file = shell, append = T)}
  for (i in cmd_09) {cat(i, file = shell, append = T)}
  for (i in cmd_10) {cat(i, file = shell, append = T)}
  for (i in cmd_11) {cat(i, file = shell, append = T)}
  for (i in cmd_12) {cat(i, file = shell, append = T)}
  #
  #two-pass mode是否有影响
  cmd_13 <- paste0("#STAR --runThreadN 10 --readFilesCommand zcat ",
                   "--outMultimapperOrder Random --outSAMtype BAM SortedByCoordinate ",
                   "--outFileNamePrefix ",path,"dumpROOM/",prex," ",
                   "--outSAMprimaryFlag AllBestScore ",
                   "--genomeDir ",inx1," --readFilesIn ",lay1," ",lay2,"\n")
  cmd_14 <- paste0("#mv ",
                   paste0(path,"dumpROOM/",prex,"Aligned.sortedByCoord.out.bam "),
                   paste0(path,"dumpROOM/",prex,".bam"),"\n")
  cmd_15 <- paste0("#samtools index ",paste0(path,"dumpROOM/",prex,".bam"),"\n")
  cmd_16 <- paste0("#TEtranscripts -t ",
                   paste(paste0(path,"dumpROOM/",prex[1:4],".bam"),collapse = " ")," -c ",
                   paste(paste0(path,"dumpROOM/",prex[5:6],".bam"),collapse = " "),
                   " --GTF ",inx3," --TE ",inx4," --sortByPos --mode multi ",
                   "--project ghyuughyvv --outdir ",paste0(path,"dumpROOM"),"\n")
  for (i in cmd_13) {cat(i, file = shell, append = T)}
  for (i in cmd_14) {cat(i, file = shell, append = T)}
  for (i in cmd_15) {cat(i, file = shell, append = T)}
  for (i in cmd_16) {cat(i, file = shell, append = T)}
  print(paste0("nohup bash ",shell, " >",paste0(path,"Log/"),basename(shell),".log"," 2>&1 &"))
}
################################################################################
#0A、带和不带2-pass mode, 对差异倍数的影响
{
  #tmp1 <- "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/202408707/theghyuughyvv/count/ghyuughyvv_sigdiff_gene_TE.txt"
  #tmp1 <- read.table(tmp1)
  #tmp2 <- "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/202408707/theghyuughyvv/dumpROOM/ghyuughyvv_sigdiff_gene_TE.txt"
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
    tmp1 <- tmp1[,c(grep(x = colnames(tmp1), pattern = paste0(mark,"-")),5,6)]
    colnames(tmp1) <- sapply(str_split(basename(colnames(tmp1)),pattern = ".bam"), "[", 1)
    dfclas <- data.frame(row.names = colnames(tmp1), Class, stringsAsFactors = T)
    colLen <- ceiling(length(colnames(tmp1))/2)
    geneAA <- tmp1[grep(rownames(tmp1), pattern = ":", invert = T),]
    teOnly <- tmp1[grep(rownames(tmp1), pattern = ":", invert = F),]
    teGene <- tmp1
    #
    geneAA <- DESeq2::DESeqDataSetFromMatrix(countData = geneAA, colData = dfclas, design = ~Class)
    #geneAA <- geneAA[rowSums(BiocGenerics::counts(geneAA) > 2) >= colLen, ]#[rowMeans(BiocGenerics::counts(geneAA)) > 1,]
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
    #
    #
    tmp2 <- read.table(siteFile, header = T, check.names = F, stringsAsFactors = F, row.names = 1)
    tmp2 <- tmp2[,c(grep(x = colnames(tmp2), pattern = paste0(mark,"-")),5,6)]
    colnames(tmp2) <- sapply(str_split(basename(colnames(tmp2)),pattern = ".bam"), "[", 1)
    colLen <- ceiling(length(colnames(tmp2))/2)
    #
    teGene <- tmp2
    teGene <- DESeq2::DESeqDataSetFromMatrix(countData = teGene, colData = dfclas, design = ~Class)
    #teGene <- teGene[rowSums(BiocGenerics::counts(teGene) > 2) >= colLen, ]#[rowMeans(BiocGenerics::counts(teGene)) > 1,]
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
    #teGene <- teGene[rowSums(BiocGenerics::counts(teGene) > 2) >= colLen, ]#[rowMeans(BiocGenerics::counts(teGene)) > 1,]
    teGene <- teGene[rowMeans(BiocGenerics::counts(teGene)) > 5,]
    teGene <- DESeq2::DESeq(teGene,)
    teGene <- DESeq2::results(teGene, contrast = VS)
    write.csv(teGene, paste0(path,"DESeq2/",mark, ".te.site.teOnly.csv"))
    upup <- subset(teGene, log2FoldChange > log2(1.5) & padj < 0.05)
    down <- subset(teGene, log2FoldChange < -log2(1.5)& padj < 0.05)
    write.csv(upup, paste0(path,"DESeq2/",mark, ".te.site.teOnly.upup.csv"))
    write.csv(down, paste0(path,"DESeq2/",mark, ".te.site.teOnly.down.csv"))
  }
  shy_diff(mark = "ghyuu", 
           geneFile = paste0(path, "count/ghyuughyvv.cntTable"),
           siteFile = paste0(path, "count/ghyuughyvv.site.cntTable"),
           Class = factor(c("OE","OE","WT","WT")), VS = c("Class","OE","WT"))
  shy_diff(mark = "ghyvv", 
           geneFile = paste0(path, "count/ghyuughyvv.cntTable"),
           siteFile = paste0(path, "count/ghyuughyvv.site.cntTable"),
           Class = factor(c("OE","OE","WT","WT")), VS = c("Class","OE","WT"))
}
################################################################################
#05、CPM Excel
{
  cnt1 <- read.table("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/202408707/theghyuughyvv/count/ghyuughyvv.cntTable", header = T, sep = "\t", row.names = 1)
  colnames(cnt1) <- c("ghyuuOE-1","ghyuuOE-2","ghyvvOE-1","ghyvvOE-2","ovctrl-1","ovctrl-2")
  tmp1 <- cnt1[grep(x = rownames(cnt1), pattern=":", invert = F),]#TE
  tmp2 <- cnt1[grep(x = rownames(cnt1), pattern=":", invert = T),]#Gene
  tmp3 <- cnt1#All
  tmp1 <- as.data.frame(apply(tmp1, 2, function(x){x/sum(x)*1000000}))
  tmp2 <- as.data.frame(apply(tmp2, 2, function(x){x/sum(x)*1000000}))
  tmp3 <- as.data.frame(apply(tmp3, 2, function(x){x/sum(x)*1000000}))
  tmp1$theName <- rownames(tmp1)
  tmp2$theName <- rownames(tmp2)
  tmp3$theName <- rownames(tmp3)
  openxlsx::write.xlsx(x = tmp1, file = "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/202408707/theghyuughyvv/count/ghyuughyvv.cpm.onlyTE.xlsx")
  openxlsx::write.xlsx(x = tmp2, file = "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/202408707/theghyuughyvv/count/ghyuughyvv.cpm.onlyGene.xlsx")
  openxlsx::write.xlsx(x = tmp3, file = "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/202408707/theghyuughyvv/count/ghyuughyvv.cpm.allGeneTE.xlsx")
}
################################################################################
#06、Volcano for Genes [OK]
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
      scale_x_continuous(limits = c(-8, 8), breaks = seq(-8, 8, 2)) +
      geom_hline(aes(yintercept=-log10(0.05)), linetype="dashed", colour="royalblue") +
      geom_vline(aes(xintercept=-log2(1.5)), linetype="dashed", colour="royalblue") +
      geom_vline(aes(xintercept=log2(1.5)), linetype="dashed", colour="royalblue") +
      theme_bw() + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      scale_color_manual(values = c("red","blue","black"),
                         labels = c(paste0(names(table(resdata$threshold)[1])," (",unname(table(resdata$threshold)[1]),")"),
                                    paste0(names(table(resdata$threshold)[2])," (",unname(table(resdata$threshold)[2]),")"),
                                    paste0(names(table(resdata$threshold)[3])," (",unname(table(resdata$threshold)[3]),")"))) +
      theme(legend.position = c(0.24,0.9), 
            legend.background = element_rect(fill = NA),
            legend.spacing.x = unit(0,"mm"),
            legend.spacing.y = unit(0,"mm"),
            legend.key = element_blank(),
            legend.text = element_text(size = 13),
            axis.title = element_text(size = 13)) +
      guides(color = guide_legend(override.aes = list(size = 3)))
    ggsave(plot = p_1, saveFile, units = "cm", width = 18, height = 16)
  }
  shy_volcano_gene(resFile  = paste0(path, "DESeq2/ghyuu.geneOnly.csv"),
                   saveFile = paste0(path, "PLOT/ghyuu.geneOnly.volcano.pdf"))
  shy_volcano_gene(resFile  = paste0(path, "DESeq2/ghyvv.geneOnly.csv"),
                   saveFile = paste0(path, "PLOT/ghyvv.geneOnly.volcano.pdf"))
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
  shy_volcano_te(resFile  = paste0(path, "DESeq2/ghyuu.te.GeneBackGround.csv"),
                 saveFile = paste0(path, "PLOT/ghyuu.te.GeneBackGround.volcano.pdf"))
  shy_volcano_te(resFile  = paste0(path, "DESeq2/ghyuu.teOnly.csv"),
                 saveFile = paste0(path, "PLOT/ghyuu.teOnly.volcano.pdf"))
  shy_volcano_te(resFile  = paste0(path, "DESeq2/ghyvv.te.GeneBackGround.csv"),
                 saveFile = paste0(path, "PLOT/ghyvv.te.GeneBackGround.volcano.pdf"))
  shy_volcano_te(resFile  = paste0(path, "DESeq2/ghyvv.teOnly.csv"),
                 saveFile = paste0(path, "PLOT/ghyvv.teOnly.volcano.pdf"))
}
################################################################################
#08、MA [OK]
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
  shy_teMA(difffile = paste0(path,"DESeq2/ghyuu.te.GeneBackGround.csv"),
           savefile = paste0(path,"PLOT/ghyuu.te.GeneBackGround.MA-v3.pdf"))
  shy_teMA(difffile = paste0(path,"DESeq2/ghyvv.te.GeneBackGround.csv"),
           savefile = paste0(path,"PLOT/ghyvv.te.GeneBackGround.MA-v3.pdf"))
  shy_teMA(difffile = paste0(path,"DESeq2/ghyuu.teOnly.csv"),
           savefile = paste0(path,"PLOT/ghyuu.teOnly.MA-v3.pdf"))
  shy_teMA(difffile = paste0(path,"DESeq2/ghyvv.teOnly.csv"),
           savefile = paste0(path,"PLOT/ghyvv.teOnly.MA-v3.pdf"))
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
  shy_repSite(resdFile = paste0(path,"DESeq2/ghyuu.te.site.teOnly.csv"),
              saveFile = paste0(path,"PLOT/ghyuu.te.site.teOnly.mervl.pdf"))
  shy_repSite(resdFile = paste0(path,"DESeq2/ghyvv.te.site.teOnly.csv"),
              saveFile = paste0(path,"PLOT/ghyvv.te.site.teOnly.mervl.pdf"))
}
################################################################################
#10、PCA [效果不好, 舍弃]
{
  #10、PCA
  {
    #----------------------------------PCA-batch----------------------------------
    cnt1 <- read.table(paste0(path, "count/ghyuughyvv.cntTable"),
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
  #10、PCA [revise] [FactoMineR]
  {
    #----------------------------------PCA-batch----------------------------------
    cnt1 <- read.table(paste0(path, "count/ghyuughyvv.cntTable"),
                       header = T, check.names = F, stringsAsFactors = F, row.names = 1)[,-c(1,2)]
    cnt1 <- cnt1[grep(rownames(cnt1), pattern = ":", invert = T),]
    colnames(cnt1) <- sapply(str_split(basename(colnames(cnt1)), ".bam"),"[",1)
    #
    cnt2 <- read.table("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuu/count/ghyuu.cntTable",
                       header = T, check.names = F, stringsAsFactors = F, row.names = 1)
    cnt2 <- cnt2[grep(rownames(cnt2), pattern = ":", invert = T),]
    colnames(cnt2) <- sapply(str_split(basename(colnames(cnt2)), ".bam"),"[",1)
    #
    cnt3 <- read.table("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/thehuhuhu/count/huhuhu.rename.cntTable",
                       header = T, check.names = F, stringsAsFactors = F, row.names = 1)[,1:4]
    cnt3 <- cnt3[grep(rownames(cnt3), pattern = ":", invert = T),]
    colnames(cnt3) <- sapply(str_split(basename(colnames(cnt3)), ".bam"),"[",1)
    #
    
  }
  ################################################################################
  #10、PCA [revise] [ghyuu_KO, ghyvv_OE, tdTomato, D1_2C]
  {
    #----------------------------------PCA-batch----------------------------------
    cnt1 <- read.table(paste0(path, "count/ghyuughyvv.cntTable"),
                       header = T, check.names = F, stringsAsFactors = F, row.names = 1)[,-c(1,2)]
    cnt1 <- cnt1[grep(rownames(cnt1), pattern = ":", invert = T),]
    colnames(cnt1) <- sapply(str_split(basename(colnames(cnt1)), ".bam"),"[",1)
    #
    cnt2 <- read.table("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuu/count/ghyuu.cntTable",
                       header = T, check.names = F, stringsAsFactors = F, row.names = 1)
    cnt2 <- cnt2[grep(rownames(cnt2), pattern = ":", invert = T),]
    colnames(cnt2) <- sapply(str_split(basename(colnames(cnt2)), ".bam"),"[",1)
    #
    cnt3 <- read.table("/ChIP_seq_2/aaSHY/ghyuu/publicData/D1_2C_cells/count/duxs.cntTable",
                       header = T, check.names = F, stringsAsFactors = F, row.names = 1)[,3:6]
    cnt3 <- cnt3[grep(rownames(cnt3), pattern = ":", invert = T),]
    colnames(cnt3) <- sapply(str_split(basename(colnames(cnt3)), ".bam"),"[",1)
    #
    base::load(file = "/Reference/aaSHY/DatabaseFiles/forPCA.RData")
    myta <- cbind(cnt1,cnt2,cnt3,forPCA)[1:13]
    #
    condition = sapply(str_split(colnames(myta),"-"), "[", 1)
    condition[c(3,4)] <- paste0("ghyvv_",condition[c(3,4)])
    condition[c(5:8)] <- c("ghyuu","ghyuu","ghyuu_WT","ghyuu_WT")
    btch = c(rep("ghyvv_OE",length(colnames(cnt1))),
             rep("ghyuu_KO",length(colnames(cnt2))),
             rep("D1__CC",length(colnames(cnt3))),
             condition[c(length(c(colnames(cnt1),colnames(cnt2),colnames(cnt3)))+1):length(myta)])
    clas <- data.frame(names = colnames(myta), condition = condition, batch = btch)
    dds1 <- suppressWarnings(DESeq2::DESeqDataSetFromMatrix(myta, colData = clas, design = ~condition))#[rowSums(counts(vtax))>5]
    dds2 <- BiocGenerics::estimateSizeFactors(dds1)
    rldf <- DESeq2::rlog(dds2)
    SummarizedExperiment::assay(rldf) <- limma::removeBatchEffect(assay(rldf), batch=rldf$batch)
    pd_1 <- plotPCA(rldf, intgroup = c("condition"), returnData = T)
    rldf <- DESeq2::vst(dds2)
    SummarizedExperiment::assay(rldf) <- limma::removeBatchEffect(assay(rldf), batch=rldf$batch)
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
              axis.title   = element_text(family = "serif", size = 12),
              plot.margin = margin(unit = "cm",1,1,1,1))
      print(ppp)
    }
    pdf(file = paste0(path,"PLOT/PCA-revise.pdf"),width = 10, height = 8)
    shy_linshi(df = pd_1, i = 1)
    shy_linshi(df = pd_2, i = 2)
    dev.off()
  }
  ################################################################################
  #10、PCA [revise] [ghyuu_KO, ghyvv_OE, tdTomato, group1]
  {
    #----------------------------------PCA-batch----------------------------------
    cnt1 <- read.table(paste0(path, "count/ghyuughyvv.cntTable"),
                       header = T, check.names = F, stringsAsFactors = F, row.names = 1)[,-c(1,2)]
    cnt1 <- cnt1[grep(rownames(cnt1), pattern = ":", invert = T),]
    colnames(cnt1) <- sapply(str_split(basename(colnames(cnt1)), ".bam"),"[",1)
    #
    cnt2 <- read.table("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuu/count/ghyuu.cntTable",
                       header = T, check.names = F, stringsAsFactors = F, row.names = 1)
    cnt2 <- cnt2[grep(rownames(cnt2), pattern = ":", invert = T),]
    colnames(cnt2) <- sapply(str_split(basename(colnames(cnt2)), ".bam"),"[",1)
    #
    cnt3 <- read.table("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/thehuhuhu/count/huhuhu.rename.cntTable",
                       header = T, check.names = F, stringsAsFactors = F, row.names = 1)[,1:4]
    cnt3 <- cnt3[grep(rownames(cnt3), pattern = ":", invert = T),]
    colnames(cnt3) <- c("huhuhuLows-rep1","huhuhuLows-rep2","huhuhuHigh-rep1","huhuhuHigh-rep2")
    #
    base::load(file = "/Reference/aaSHY/DatabaseFiles/forPCA.RData")
    myta <- cbind(cnt1,cnt2,cnt3,forPCA)[1:13]
    #
    condition = sapply(str_split(colnames(myta),"-"), "[", 1)
    condition[c(3,4)] <- paste0("ghyvv_",condition[c(3,4)])
    condition[c(5:8)] <- c("ghyuu","ghyuu","ghyuu_WT","ghyuu_WT")
    btch = c(rep("ghyvv_OE",length(colnames(cnt1))),
             rep("ghyuu_KO",length(colnames(cnt2))),
             rep("huhuhucose",length(colnames(cnt3))),
             condition[c(length(c(colnames(cnt1),colnames(cnt2),colnames(cnt3)))+1):length(myta)])
    clas <- data.frame(names = colnames(myta), condition = condition, batch = btch)
    dds1 <- suppressWarnings(DESeq2::DESeqDataSetFromMatrix(myta, colData = clas, design = ~condition))#[rowSums(counts(vtax))>5]
    dds2 <- BiocGenerics::estimateSizeFactors(dds1)
    rldf <- DESeq2::rlog(dds2)
    SummarizedExperiment::assay(rldf) <- limma::removeBatchEffect(assay(rldf), batch=rldf$batch)
    pd_1 <- plotPCA(rldf, intgroup = c("condition"), returnData = T)
    rldf <- DESeq2::vst(dds2)
    SummarizedExperiment::assay(rldf) <- limma::removeBatchEffect(assay(rldf), batch=rldf$batch)
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
              axis.title   = element_text(family = "serif", size = 12),
              plot.margin = margin(unit = "cm",1,1,1,1))
      print(ppp)
    }
    pdf(file = paste0(path,"PLOT/PCA-revise.pdf"),width = 10, height = 8)
    shy_linshi(df = pd_1, i = 1)
    shy_linshi(df = pd_2, i = 2)
    dev.off()
  }
  ################################################################################
  #
  #
  #
  #10、PCA [revise] [group1, D1_2C]
  {
    #----------------------------------PCA-batch----------------------------------
    cnt2 <- read.table("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/thehuhuhu/count/huhuhu.rename.cntTable",
                       header = T, check.names = F, stringsAsFactors = F, row.names = 1)[,1:4]
    cnt2 <- cnt2[grep(rownames(cnt2), pattern = ":", invert = T),]
    colnames(cnt2) <- c("huhuhuLows-rep1","huhuhuLows-rep2","huhuhuHigh-rep1","huhuhuHigh-rep2")
    #
    cnt3 <- read.table("/ChIP_seq_2/aaSHY/ghyuu/publicData/D1_2C_cells/count/duxs.cntTable",
                       header = T, check.names = F, stringsAsFactors = F, row.names = 1)[,3:6]
    cnt3 <- cnt3[grep(rownames(cnt3), pattern = ":", invert = T),]
    colnames(cnt3) <- sapply(str_split(basename(colnames(cnt3)), ".bam"),"[",1)
    #
    myta <- cbind(cnt2,cnt3)
    #
    condition = sapply(str_split(colnames(myta),"-"), "[", 1)
    btch = c(rep("huhuhucose",length(colnames(cnt2))),
             rep("D1_2C",length(colnames(cnt3))))
    clas <- data.frame(names = colnames(myta), condition = condition, batch = btch)
    dds1 <- suppressWarnings(DESeq2::DESeqDataSetFromMatrix(myta, colData = clas, design = ~condition))#[rowSums(counts(vtax))>5]
    dds2 <- BiocGenerics::estimateSizeFactors(dds1)
    rldf <- DESeq2::rlog(dds2)
    SummarizedExperiment::assay(rldf) <- limma::removeBatchEffect(assay(rldf), batch=rldf$batch)
    pd_1 <- plotPCA(rldf, intgroup = c("condition"), returnData = T)
    rldf <- DESeq2::vst(dds2)
    SummarizedExperiment::assay(rldf) <- limma::removeBatchEffect(assay(rldf), batch=rldf$batch)
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
              legend.position = c(0.8,0.25),
              axis.text    = element_text(family = "serif", size = 12),
              axis.title   = element_text(family = "serif", size = 12),
              plot.margin = margin(unit = "cm",1,1,1,1))
      print(ppp)
    }
    pdf(file = paste0(path,"PLOT/PCA-revise.pdf"),width = 10, height = 8)
    shy_linshi(df = pd_1, i = 1)
    shy_linshi(df = pd_2, i = 2)
    dev.off()
  }
  #
  #
  #
  #10、PCA [revise] [ghyuu_KO, ghyvv_OE, tdTomato, D1_2C]
  {
    #----------------------------------PCA-batch----------------------------------
    cnt1 <- read.table(paste0(path, "count/ghyuughyvv.cntTable"),
                       header = T, check.names = F, stringsAsFactors = F, row.names = 1)[,-c(1,2)]
    cnt1 <- cnt1[grep(rownames(cnt1), pattern = ":", invert = T),]
    colnames(cnt1) <- sapply(str_split(basename(colnames(cnt1)), ".bam"),"[",1)
    #
    cnt2 <- read.table("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuu/count/ghyuu.cntTable",
                       header = T, check.names = F, stringsAsFactors = F, row.names = 1)
    cnt2 <- cnt2[grep(rownames(cnt2), pattern = ":", invert = T),]
    colnames(cnt2) <- sapply(str_split(basename(colnames(cnt2)), ".bam"),"[",1)
    
    myta <- cbind(cnt1,cnt2)
    #
    condition = sapply(str_split(colnames(myta),"-"), "[", 1)
    condition[c(3,4)] <- paste0("ghyvv_",condition[c(3,4)])
    condition[c(5:8)] <- c("ghyuu","ghyuu","ghyuu_WT","ghyuu_WT")
    btch = c(rep("ghyvv_OE",length(colnames(cnt1))),
             rep("ghyuu_KO",length(colnames(cnt2))))
    clas <- data.frame(names = colnames(myta), condition = condition, batch = btch)
    dds1 <- suppressWarnings(DESeq2::DESeqDataSetFromMatrix(myta, colData = clas, design = ~condition))#[rowSums(counts(vtax))>5]
    dds2 <- BiocGenerics::estimateSizeFactors(dds1)
    rldf <- DESeq2::rlog(dds2)
    SummarizedExperiment::assay(rldf) <- limma::removeBatchEffect(assay(rldf), batch=rldf$batch)
    pd_1 <- plotPCA(rldf, intgroup = c("condition"), returnData = T)
    rldf <- DESeq2::vst(dds2)
    SummarizedExperiment::assay(rldf) <- limma::removeBatchEffect(assay(rldf), batch=rldf$batch)
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
              axis.title   = element_text(family = "serif", size = 12),
              plot.margin = margin(unit = "cm",1,1,1,1))
      print(ppp)
    }
    pdf(file = paste0(path,"PLOT/PCA-revise.pdf"),width = 10, height = 8)
    shy_linshi(df = pd_1, i = 1)
    shy_linshi(df = pd_2, i = 2)
    dev.off()
  }
  #
  #
  #
  #10、相关性热图
  {
    bath <- "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuu/"
    bw01 <- list.files("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuu/bamCoverage", 
                       full.names = T)
    bw02 <- list.files("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuughyvv/bamCoverage",
                       full.names = T)[c(5,6,1,2)]
    tags <- gsub(".bw","",basename(c(bw01,bw02)))
    cmd_01 <- paste0("multiBigwigSummary bins -p 20 -b ", 
                     paste0(c(bw01,bw02), collapse = " "),
                     " -l ",paste0(tags,collapse = " "),
                     " -o ",paste0(bath,"Corelation/pearsonCorre.npz","\n"))
    cmd_02 <- paste0("plotCorrelation -in ",bath,
                     "Corelation/pearsonCorre.npz ",
                     "-c pearson -p heatmap --plotFileFormat pdf --removeOutliers ",
                     "--colorMap bwr -o ", bath,"Corelation/pearsonCorre_v1.pdf ",
                     "--outFileCorMatrix ",bath,"Corelation/pearsonCorre_v1.index.txt","\n")
    #
    #
    #
    #
    bath <- "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuu/"
    bw01 <- list.files("/ChIP_seq_2/aaSHY/ghyuu/publicData/D1_2C_cells/bamCoverage", 
                       full.names = T)[3:6]
    bw02 <- list.files("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/thehuhuhu/bamCoverage",
                       full.names = T)[1:4]
    tags <- gsub(".bw","",basename(c(bw01,bw02)))
    cmd_03 <- paste0("multiBigwigSummary bins -p 20 -b ", 
                     paste0(c(bw01,bw02), collapse = " "),
                     " -l ",paste0(tags,collapse = " "),
                     " -o ",paste0(bath,"Corelation/pearsonCorre.npz","\n"))
    cmd_04 <- paste0("plotCorrelation -in ",bath,
                     "Corelation/pearsonCorre.npz --zMin 0.9 --zMax 1 ",
                     "-c pearson -p heatmap --plotFileFormat pdf --removeOutliers ",
                     "--colorMap bwr -o ", bath,"Corelation/pearsonCorre_v2.pdf ",
                     "--outFileCorMatrix ",bath,"Corelation/pearsonCorre_v2.index.txt","\n")
  }
  #
  #
  #
  ################################################################################
  #10、PCA [revise] [3D]
  {
    #----------------------------------PCA-batch----------------------------------
    cnt1 <- read.table(paste0(path, "count/ghyuughyvv.cntTable"),
                       header = T, check.names = F, stringsAsFactors = F, row.names = 1)[,-c(1,2)]
    cnt1 <- cnt1[grep(rownames(cnt1), pattern = ":", invert = T),]
    colnames(cnt1) <- sapply(str_split(basename(colnames(cnt1)), ".bam"),"[",1)
    #
    cnt2 <- read.table("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuu/count/ghyuu.cntTable",
                       header = T, check.names = F, stringsAsFactors = F, row.names = 1)
    cnt2 <- cnt2[grep(rownames(cnt2), pattern = ":", invert = T),]
    colnames(cnt2) <- sapply(str_split(basename(colnames(cnt2)), ".bam"),"[",1)
    #
    base::load(file = "/Reference/aaSHY/DatabaseFiles/forPCA.RData")
    myta <- cbind(cnt1,cnt2,forPCA)[1:9]
    #
    condition = sapply(str_split(colnames(myta),"-"), "[", 1)
    btch = c(rep("ghyvv_OE",length(colnames(cnt1))),
             rep("ghyuu_KO",length(colnames(cnt2))),
             condition[c(length(c(colnames(cnt1),colnames(cnt2)))+1):length(myta)])
    clas <- data.frame(names = colnames(myta), condition = condition, batch = btch)
    dds1 <- suppressWarnings(DESeq2::DESeqDataSetFromMatrix(myta, colData = clas, design = ~condition))#[rowSums(counts(vtax))>5]
    dds2 <- BiocGenerics::estimateSizeFactors(dds1)
    #
    rldf <- DESeq2::rlog(dds2)
    SummarizedExperiment::assay(rldf) <- limma::removeBatchEffect(assay(rldf), batch=rldf$batch)
    df01 <- as.data.frame(prcomp(SummarizedExperiment::assay(rldf))$rotation)
    df01 <- df01[,c(1,2,3)]; df01$cond <- clas$condition; df01$batc <- btch
    plotly::plot_ly(data = df01, type = "scatter3d", mode = "markers",
                    x = ~PC1, y = ~PC2, z = ~PC3, color = ~cond)
    #
    vsdf <- DESeq2::vst(dds2)
    df01 <- as.data.frame(prcomp(SummarizedExperiment::assay(vsdf))$rotation)
    df01 <- df01[,c(1,2,3)]; df01$cond <- clas$condition; df01$batc <- btch
    plotly::plot_ly(data = df01, type = "scatter3d", mode = "markers",
                    x = ~PC1, y = ~PC2, z = ~PC3, color = ~cond)
    #ggplot(df01, aes(x=PC1, y=PC2, z=PC3, color=cond)) + theme_void() + axes_3D() + stat_3D()
    
    
    
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
    pdf(file = paste0(path,"PLOT/PCA-revise.pdf"),width = 10, height = 8)
    shy_linshi(df = pd_1, i = 1)
    shy_linshi(df = pd_2, i = 2)
    dev.off()
  }
}
################################################################################
#11、pheatmap
{
  shy_pheatmap <- function(readsFile, gene_resFile, tete_resFile, savePath, mark){
    print("reads counted by TEtranscripts")
    #
    data <- read.table(readsFile, header = T, row.names = 1, check.names = F, stringsAsFactors = F)
    data <- data[,c(grep(x = colnames(data), pattern = paste0(mark,"-")),5,6)]
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
    sigG <- resG[which(abs(resG$log2FoldChange)>log2(1.5) & resG$padj<0.05), "X"]
    sigG_upup <- resG[which(resG$log2FoldChange>log2(1.5) & resG$padj<0.05), "X"]
    sigG_down <- resG[which(resG$log2FoldChange < -log2(1.5) & resG$padj<0.05), "X"]
    resG <- as.data.frame(read.csv(tete_resFile, header = T, stringsAsFactors = F))
    resG <- na.omit(resG)
    sigT <- resG[which(abs(resG$log2FoldChange)>log2(1.5) & resG$padj<0.05), "X"]
    sigT_upup <- resG[which(resG$log2FoldChange>log2(1.5) & resG$padj<0.05), "X"]
    sigT_down <- resG[which(resG$log2FoldChange < -log2(1.5) & resG$padj<0.05), "X"]
    ##开始画图
    pdf(paste0(savePath, str_split(basename(readsFile), pattern = "\\.")[[1]][1],
               ".by_",mark, "_pheatmap.v1.pdf"), width = 8, height = 10)
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
    ######pheatmap(mat = tete[which(rownames(tete) %in% sigT_down),],
    ######         main="the heatmap of Downregulated TEs",
    ######         scale = "row", cellwidth = 60, 
    ######         show_rownames = T, show_colnames = T, cluster_row = T, cluster_cols = T,
    ######         labels_col = colnames(tete), angle_col = "45",
    ######         breaks = seq(-2,2,by=0.01), border = F,
    ######         color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
    ######                    colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
    
    dev.off()
  }
  shy_pheatmap(mark = "ghyuu", readsFile = paste0(path, "count/ghyuughyvv.cntTable"),
               gene_resFile = paste0(path, "DESeq2/ghyuu.geneOnly.csv"),
               tete_resFile = paste0(path, "DESeq2/ghyuu.teOnly.csv"),
               savePath = paste0(path, "PLOT/pheatmap/"))
  shy_pheatmap(mark = "ghyvv", readsFile = paste0(path, "count/ghyuughyvv.cntTable"),
               gene_resFile = paste0(path, "DESeq2/ghyvv.geneOnly.csv"),
               tete_resFile = paste0(path, "DESeq2/ghyvv.teOnly.csv"),
               savePath = paste0(path, "PLOT/pheatmap/"))
}
################################################################################
#12、Up TEs loci number [OK]
{
  shy_lociNum <- function(resdFile, saveFile, mark){
    ress <- read.csv(file = resdFile, header = T, row.names = 1, stringsAsFactors = F)
    tenm <- as.data.frame(table(sapply(str_split(string = rownames(ress),pattern = ":|_dup"), "[",1)))
    rela <- paste(sapply(str_split(string = rownames(ress),pattern = ":"), "[",2),
                  sapply(str_split(string = rownames(ress),pattern = ":"), "[",3), sep = "__")
    tenm$clas <- sapply(str_split(as.character(as.data.frame(table(rela))[,1]),"__"),"[",2)
    tenm <- tenm[order(tenm$Freq,decreasing = T),]
    forP <- tenm[1:10,]
    forP <- forP[!(is.na(forP$Freq)),]
    forP$Var1 = factor(forP$Var1, levels = as.character(forP$Var1))
    p_1 <- ggplot(forP, aes(Var1, Freq, fill = clas)) + 
      geom_bar(stat = "identity",width = 0.6) +
      theme_classic() +
      theme(axis.title = element_text(family = "serif", color = "black", size = 14)) +
      theme(axis.text = element_text(size = 12, family = "serif", color = "black")) +
      theme(legend.title = element_blank(), 
            title = element_text(family = "serif", color = "black", size = 14)) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4)) +
      theme(legend.position = c(0.7,0.8)) +
      ylab("Number of Loci") + 
      xlab(NULL) +
      labs(title = expression(deRegulated~TEs~"in"~italic(geneName)^{"+/+"}~ESCs)) +
      scale_fill_manual(breaks = unique(forP$clas),values = ggsci::pal_aaas(palette = "default")(length(unique(forP$clas))))
    ggsave(plot = p_1, 
           filename = saveFile, units = "cm", width = 12, height = 12)
  }
  shy_lociNum(mark = "ghyuu",
              resdFile = paste0(path, "DESeq2/ghyuu.te.site.teOnly.upup.csv"),
              saveFile = paste0(path, "PLOT/TElociNum/ghyuu.lociNum.teOnly.upup.pdf"))
  shy_lociNum(mark = "ghyuu",
              resdFile = paste0(path, "DESeq2/ghyuu.te.site.teOnly.down.csv"),
              saveFile = paste0(path, "PLOT/TElociNum/ghyuu.lociNum.teOnly.down.pdf"))
  shy_lociNum(mark = "ghyvv",
              resdFile = paste0(path, "DESeq2/ghyvv.te.site.teOnly.upup.csv"),
              saveFile = paste0(path, "PLOT/TElociNum/ghyvv.lociNum.teOnly.upup.pdf"))
  shy_lociNum(mark = "ghyvv",
              resdFile = paste0(path, "DESeq2/ghyvv.te.site.teOnly.down.csv"),
              saveFile = paste0(path, "PLOT/TElociNum/ghyvv.lociNum.teOnly.down.pdf"))
  shy_lociNum(mark = "ghyvv",
              resdFile = paste0(path, "DESeq2/ghyvv.te.site.GeneBackGround.upup.csv"),
              saveFile = paste0(path, "PLOT/TElociNum/ghyvv.lociNum.GeneBackGround.upup.pdf"))
  shy_lociNum(mark = "ghyvv",
              resdFile = paste0(path, "DESeq2/ghyvv.te.site.GeneBackGround.down.csv"),
              saveFile = paste0(path, "PLOT/TElociNum/ghyvv.lociNum.GeneBackGround.down.pdf"))
}
################################################################################
#13、enrich GO BP
{
  resFile <- paste0(path,"DESeq2/ghyvv.geneOnly.csv")
  resdata <- as.data.frame(read.csv(resFile, header = T, stringsAsFactors = F))
  resdata <- resdata[grep(resdata[,1], pattern = ":", invert = T),]
  mtacup <- resdata[which(resdata$log2FoldChange > log2(1.5)  &  resdata$padj < 0.05),"X"]
  kddown <- resdata[which(resdata$log2FoldChange < -log2(1.5) &  resdata$padj < 0.05),"X"]
  for (symb in seq(2)) {
    ovgene <- unlist(list(mtacup,kddown)[symb])
    geyy <- na.omit(mapIds(x=org.Mm.eg.db,
                           keys=unique(ovgene), keytype="SYMBOL", column="ENTREZID"))
    gogo <- clusterProfiler::enrichGO(gene = geyy, 
                                      OrgDb = org.Mm.eg.db, keyType = "ENTREZID", 
                                      ont = "BP",pvalueCutoff = 0.05,qvalueCutoff = 0.05, readable = T)
    gogo <- gogo@result
    assign(paste0("gta",symb), gogo[which(gogo$pvalue<0.01),])
  }
  openxlsx::write.xlsx(x = gta1, file = paste0(path,"PLOT/ghyvv.GO.BP.upup.xlsx"))
  openxlsx::write.xlsx(x = gta2, file = paste0(path,"PLOT/ghyvv.GO.BP.down.xlsx"))
  #
  gta1 <- as.data.frame(gta1)
  gta1$Description[9] <- "siRNA-mediated retrotransposon silencing\n by heterochromatin formation"
  gta1$Description <- factor(gta1$Description, levels = rev(gta1$Description))
  p_05 <- ggplot(gta1[1:10,])+
    geom_bar(aes(x=-log10(pvalue), y=Description), stat = "identity", fill="#1ba784", alpha=1.0) +
    theme_few() +
    scale_x_continuous(breaks = seq(0,16,by=2)) +
    geom_vline(aes(xintercept = 1), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 2), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 3), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 4), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 5), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 6), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 7), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 8), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 9), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 10), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 11), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 12), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 13), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 14), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 15), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 16), color = "white", linewidth  = 0.8) +
    labs(y="",title = "this is title--upup") +
    theme(text = element_text(size = 12,family = "serif",color = "black"),
          panel.border = element_rect(linewidth = 1),
          plot.title = element_text(size = 12,family = "serif",color = "black",face = "bold"),
          axis.title = element_text(size = 12,family = "serif",color = "black"),
          axis.text  = element_text(size = 12,family = "serif",color = "black"))
  p_05
  #
  gta2$Description <- factor(gta2$Description, levels = rev(gta2$Description))
  p_06 <- ggplot(gta2[1:10,])+
    geom_bar(aes(x=-log10(pvalue), y=Description), stat = "identity", fill="#5861AC", alpha=1.0) +
    theme_few() +
    scale_x_continuous(breaks = seq(10)) +
    geom_vline(aes(xintercept = 1), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 2), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 3), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 4), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 5), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 6), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 7), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 8), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 9), color = "white", linewidth  = 0.8) +
    labs(y="",title = "this is title--down") +
    theme(text = element_text(size = 12,family = "serif",color = "black"),
          panel.border = element_rect(linewidth = 1),
          plot.title = element_text(size = 12,family = "serif",color = "black",face = "bold"),
          axis.title = element_text(size = 12,family = "serif",color = "black"),
          axis.text  = element_text(size = 12,family = "serif",color = "black"))
  p_06
  plot <- p_05 +p_06 +plot_layout(nrow = 2)
  ggsave(plot = plot,
         units = "cm", width = 22, height = 15,
         filename = paste0(path,"PLOT/","ghyvv.GOBP.pdf"))
}
################################################################################
#14、enrich KEGG [又不要条形图]
{
  resFile <- paste0(path,"DESeq2/ghyvv.geneOnly.csv")
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
    assign(paste0("kta",symb), kegg[which(kegg$pvalue<0.01),])
  }
  #
  kta1$Description <- factor(kta1$Description, levels = rev(kta1$Description))
  p_05 <- ggplot(kta1[1:7,])+
    geom_bar(aes(x=-log10(pvalue), y=Description), stat = "identity", fill="#1ba784", alpha=1.0) +
    theme_few() +
    scale_x_continuous(breaks = seq(8)) +
    geom_vline(aes(xintercept = 1), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 2), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 3), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 4), color = "white", linewidth  = 0.8) +
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
    scale_x_continuous(breaks = seq(10)) +
    geom_vline(aes(xintercept = 1), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 2), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 3), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 4), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 5), color = "white", linewidth  = 0.8) +
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
         filename = paste0(path,"PLOT/enrich/","ghyvv.KEGG.v1.pdf"))
}
################################################################################
#14、enrich KEGG [气泡图]
#KEGG
{
  resFile <- paste0(path, "DESeq2/ghyvv.geneOnly.csv")
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
  openxlsx::write.xlsx(x = kta1, file = paste0(path,"PLOT/ghyvv.kegg.upup.xlsx"))
  openxlsx::write.xlsx(x = kta2, file = paste0(path,"PLOT/ghyvv.kegg.down.xlsx"))
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
         filename = paste0(path,"PLOT/ghyvv.KEGG.pdf"))
}
################################################################################
#15、enrich GO MF
{
  resFile <- paste0(path,"DESeq2/ghyvv.geneOnly.csv")
  resdata <- as.data.frame(read.csv(resFile, header = T, stringsAsFactors = F))
  resdata <- resdata[grep(resdata[,1], pattern = ":", invert = T),]
  mtacup <- resdata[which(resdata$log2FoldChange > log2(1.5)  &  resdata$padj < 0.05),"X"]
  kddown <- resdata[which(resdata$log2FoldChange < -log2(1.5) &  resdata$padj < 0.05),"X"]
  for (symb in seq(2)) {
    ovgene <- unlist(list(mtacup,kddown)[symb])
    geyy <- na.omit(mapIds(x=org.Mm.eg.db, 
                           keys=unique(ovgene), keytype="SYMBOL", column="ENTREZID"))
    gogo <- clusterProfiler::enrichGO(gene = geyy, 
                                      OrgDb = org.Mm.eg.db, keyType = "ENTREZID", 
                                      ont = "MF",pvalueCutoff = 0.05,qvalueCutoff = 0.05, readable = T)
    gogo <- gogo@result
    assign(paste0("gta",symb), gogo[which(gogo$pvalue<0.01),])
  }
  #
  gta1$Description <- factor(gta1$Description, levels = rev(gta1$Description))
  p_05 <- ggplot(gta1[1:10,])+
    geom_bar(aes(x=-log10(pvalue), y=Description), stat = "identity", fill="#1ba784", alpha=1.0) +
    theme_few() +
    scale_x_continuous(breaks = seq(13)) +
    geom_vline(aes(xintercept = 1), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 2), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 3), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 4), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 5), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 6), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 7), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 8), color = "white", linewidth  = 0.8) +
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
  gta2$Description[8] <- "oxidoreductase activity, acting on the CH-NH group \nof donors, NAD or NADP as acceptor"
  gta2$Description <- factor(gta2$Description, levels = rev(gta2$Description))
  p_06 <- ggplot(gta2[1:10,])+
    geom_bar(aes(x=-log10(pvalue), y=Description), stat = "identity", fill="#5861AC", alpha=1.0) +
    theme_few() +
    scale_x_continuous(breaks = seq(10)) +
    geom_vline(aes(xintercept = 1), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 2), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 3), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 4), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 5), color = "white", linewidth  = 0.8) +
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
         filename = paste0(path,"PLOT/enrich/","ghyvv.GOMF.pdf"))
}
################################################################################
#16、GSEA
{
  shy_gsea <- function(resFile, savePath){
    aa <- read.csv(file = resFile, header = T, stringsAsFactors = F)[,c(1,3)]
    ab <- suppressWarnings(clusterProfiler::bitr(geneID = aa$X, fromType = "SYMBOL", 
                                                 toType = "ENTREZID", 
                                                 OrgDb = "org.Mm.eg.db", drop = T))
    ac <- dplyr::distinct(.data = ab, SYMBOL, .keep_all = T)
    ac$rank <- apply(ac, 1, function(x){aa[which(aa$X==x[1]),2]})
    ad <- ac[order(ac$rank, decreasing = T),c(2,3)]
    ae <- ad$rank; names(ae) <- ad$ENTREZID
    rm(aa,ab,ac,ad)
    af_gogo <- suppressWarnings(clusterProfiler::gseGO(geneList = ae, 
                                                       ont = "BP", seed = T, maxGSSize = 1000, verbose = T, keyType = "ENTREZID", pvalueCutoff = 1, OrgDb = "org.Mm.eg.db", by = "fgsea"))
    testxxx <- af_gogo@result
    testxxx[,10] <- apply(af_gogo@result, 1, function(x){chartr(old = ", ", new = "; ", x[10])})
    write.csv(testxxx, row.names = F, quote = F,
              file = paste0(savePath, str_split(basename(resFile),".csv")[[1]][1], ".result.GOBP.csv"))
    af_kegg <- suppressWarnings(clusterProfiler::gseKEGG(geneList = ae, maxGSSize = 1000, organism = "mmu", pvalueCutoff = 1, use_internal_data = F, seed = T, by = "fgsea"))
    temp <- af_kegg@result$Description
    temp <- sapply(str_split(temp, pattern = " - Mus"), "[", 1)
    af_kegg@result$Description <- temp
    testxxx <- af_kegg@result
    testxxx[,10] <- apply(af_kegg@result, 1, function(x){chartr(old = ", ", new = "; ", x[10])})
    write.csv(testxxx, row.names = F, quote = F,
              file = paste0(savePath, str_split(basename(resFile),".csv")[[1]][1], ".result.KEGG.csv"))
    ##2-cell Genes
    aa <- clusterProfiler::read.gmt("/Reference/aaSHY/GSEAgmt/TwoGenes.gmt")
    aa <- suppressWarnings(clusterProfiler::bitr(aa$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db"))
    ab <- data.frame(ont = rep("act_2cell_gene", length(unique(aa$ENTREZID))), gen = unique(aa$ENTREZID))
    af_gsea <- suppressWarnings(clusterProfiler::GSEA(geneList = ae, maxGSSize = 1000, TERM2GENE = ab, pvalueCutoff = 1))
    rm(aa,ab)
    #af_reac <- suppressWarnings(ReactomePA::gsePathway(geneList = ae, maxGSSize = 1000, eps = NA, organism = "mouse", pvalueCutoff = 1, verbose = T, seed = T, by = "fgsea"))
    ##开始画图
    #dev.new()
    p_a <- enrichplot::gseaplot2(x = af_gsea, geneSetID = 1, title = "GSEA for 2-cell Genes", pvalue_table = T)
    p_1 <- enrichplot::gseaplot2(x = af_gogo, 1:5, title = "GO_BP_GSEA: top 5", rel_heights = c(1.5, 0.3, 0.6), color = ggsci::pal_lancet()(5), pvalue_table = T)
    p_2 <- enrichplot::ridgeplot(x = af_gogo, showCategory = 50, fill = "pvalue", label_format = 100)
    p_3 <- enrichplot::cnetplot(setReadable(af_gogo, OrgDb = "org.Mm.eg.db", keyType = "ENTREZID"), color.params = list(foldChange = ae))
    p_4 <- enrichplot::gseaplot2(x = af_kegg, 1:5, title = "KEGG__GSEA: top 5", rel_heights = c(1.5, 0.3, 0.6), color = ggsci::pal_lancet()(5), pvalue_table = T)
    p_5 <- enrichplot::ridgeplot(x = af_kegg, showCategory = 50, fill = "pvalue", label_format = 100)
    p_6 <- enrichplot::cnetplot(setReadable(af_kegg, OrgDb = "org.Mm.eg.db", keyType = "ENTREZID"), color.params = list(foldChange = ae))
    pdf(file = paste0(savePath, str_split(basename(resFile),".csv")[[1]][1], ".GSEA.pdf"), width = 12, height = 10)
    print(p_a)
    print(p_1)
    print(p_2)
    print(p_3)
    print(p_4)
    print(p_5)
    print(p_6)
    dev.off()
    #path=paste0(savePath,"pathview/",str_split(basename(resFile), "\\.")[[1]][2])
    #dir.create(path, recursive = T); setwd(path)
    #for (i in seq(10)){
    #  suppressWarnings(pathview::pathview(gene.data = ae, pathway.id = af_kegg@result$ID[i], same.layer = F,
    #                                      pdf.size = c(12,12), kegg.native = F, species = "mmu", map.symbol = T, gene.annotpkg = "org.Mm.eg.db"))
    #}
  }
  shy_gsea(savePath = paste0(path, "PLOT/enrichGSEA/"),
           resFile = paste0(path, "DESeq2/ghyuu.geneOnly.csv"))
  shy_gsea(savePath = paste0(path, "PLOT/enrichGSEA/"),
           resFile = paste0(path, "DESeq2/ghyvv.geneOnly.csv"))
}
################################################################################
#17、python GSEA
{
  shy_gseapy <- function(mark) {
    temp <- read.csv(paste0(path, "DESeq2/",mark,".geneOnly.csv"), header = T)
    temp <- temp[, c("X", "log2FoldChange")]
    temp <- temp[order(temp$log2FoldChange, decreasing = T),]
    write.table(x = temp, 
                file = paste0(path, "DESeq2/",mark,".res.rnk"), quote = F, sep = "\t", col.names = F, row.names = F)
    cmd_12 <- paste0("gseapy prerank -g /Reference/aaSHY/GSEAgmt/TwoGenes.gmt -r ",
                     path,"DESeq2/",mark,".res.rnk -p 10 -o ", 
                     path,"PLOT/enrichGSEApy/")#conda activate base gseapy 0.10.8
    print(cmd_12)
  }
  shy_gseapy(mark = "ghyuu")
  shy_gseapy(mark = "ghyvv")
}
{
  shy_gseapy <- function(mark) {
    temp <- read.csv(paste0(path, "DESeq2/",mark,".geneOnly.csv"), header = T)
    temp <- temp[, c("X", "log2FoldChange")]
    temp <- temp[order(temp$log2FoldChange, decreasing = T),]
    write.table(x = temp, 
                file = paste0(path, "DESeq2/",mark,".res.rnk"), quote = F, sep = "\t", col.names = F, row.names = F)
    cmd_12 <- paste0("gseapy prerank -g /Reference/aaSHY/GSEAgmt/2022_majorZGA.gmt -r ",
                     path,"DESeq2/",mark,".res.rnk --max-size 1000 -p 10 -o ", 
                     path,"PLOT/enrichGSEApy/")#conda activate base gseapy 0.10.8
    print(cmd_12)
  }
  shy_gseapy(mark = "ghyuu")
  shy_gseapy(mark = "ghyvv")
}
################################################################################
#22、可变剪接
{
  shy_rMATS <- function(mark) {
    bam1 <- list.files(paste0(path,"bam"), pattern = paste0(mark,".*bam$"),full.names = T)
    bam2 <- list.files(paste0(path,"bam"), pattern = "ovc.*bam$",full.names = T)
    cat(file = paste0(path,"rMATS/",mark,"/opt.txt"), sep = ",", bam1)
    cat(file = paste0(path,"rMATS/",mark,"/ctr.txt"), sep = ",", bam2)
    paste0("rmats.py --nthread 10 -t paired --readLength 140 ",
           "--gtf ",inx3," --b1 ",path,"rMATS/",mark,"/opt.txt ",
           "--b2 ",path,"rMATS/",mark,"/ctr.txt ",
           "--od ",path,"rMATS/",mark," --tmp ",path,"rMATS/",mark,"/tmp","\n")
  }
  shy_rMATS(mark = "ghyuu")
  shy_rMATS(mark = "ghyvv")
}
################################################################################
#
#
#
#
#
#
#
################################################################################
#23、ghyvv TE sites不要MA的了, 要换成Boxplot的, 首先把count模式改成unique的
{
  prex <- sapply(str_split(list.files((paste0(path,"fq/")),pattern = "*1.fq.gz$*"),"_1.fq.gz"), "[",1)
  cmd_07 <- paste0("TEtranscripts ",
                   "-t ",paste(paste0(path,"bam/",prex[3:4],".bam"),collapse = " ")," ",
                   "-c ",paste(paste0(path,"bam/",prex[5:6],".bam"),collapse = " ")," ",
                   "--GTF ",inx3," --TE ",inx4," --sortByPos --mode uniq ",
                   "--project ghyvv --outdir ",paste0(path,"countUniq"),"\n")
  cmd_08 <- paste0("TElocal ",
                   "-b ",paste0(path,"bam/",prex[3:6],".bam")," ",
                   "--GTF ",inx3," --TE ",inx5," --sortByPos --mode uniq ",
                   "--project ",paste0(path, "countUniq/",prex[3:6]),"\n")
  shell <- paste0(path,"src/run_b.sh"); cat("#!/bin/bash\n", file = shell)
  for (i in cmd_07) {cat(i, file = shell, append = T)}
  print(paste0("nohup bash ",shell, " >",paste0(path,"Log/"),basename(shell),".log"," 2>&1 &"))
  shell <- paste0(path,"src/run_c.sh"); cat("#!/bin/bash\n", file = shell)
  for (i in cmd_08) {cat(i, file = shell, append = T)}
  print(paste0("nohup bash ",shell, " >",paste0(path,"Log/"),basename(shell),".log"," 2>&1 &"))
}
################################################################################
#24、计算CPM
{
  cnt1 <- read.table(paste0(path,"countUniq/ghyvv.site.cntTable"), header = T, sep = "\t", row.names = 1)
  colnames(cnt1) <- c("ghyvvOE-1","ghyvvOE-2","ovctrl-1","ovctrl-2")
  tmp1 <- cnt1[grep(x = rownames(cnt1), pattern=":", invert = F),]#TE
  tmp2 <- cnt1[grep(x = rownames(cnt1), pattern=":", invert = T),]#Gene
  tmp3 <- cnt1#All
  tmp1 <- as.data.frame(apply(tmp1, 2, function(x){x/sum(x)*1000000}))
  tmp2 <- as.data.frame(apply(tmp2, 2, function(x){x/sum(x)*1000000}))
  tmp3 <- as.data.frame(apply(tmp3, 2, function(x){x/sum(x)*1000000}))
  tmp1$theName <- rownames(tmp1)
  tmp2$theName <- rownames(tmp2)
  tmp3$theName <- rownames(tmp3)
  openxlsx::write.xlsx(x = tmp1, file = paste0(path,"countUniq/ghyvv.site.cpm.onlyTE.xlsx"))
  #openxlsx::write.xlsx(x = tmp2, file = paste0(path,"countUniq/ghyvv.site.cpm.onlyGene.xlsx"))
  #openxlsx::write.xlsx(x = tmp3, file = paste0(path,"countUniq/ghyvv.site.cpm.allGeneTE.xlsx"))
}
################################################################################
#25、画TE sites的boxplot图
{
  shy_d <- function(grou,mark) {
    mervs <- openxlsx::read.xlsx(paste0(path,"countUniq/ghyvv.site.cpm.onlyTE.xlsx"))
    mervl <- mervs[grep("MERVL-int",mervs$theName),]
    mervl <- mervl[rowSums(mervl[,1:4]==0) < 4,]
    mervl <- tidyr::pivot_longer(data = mervl, cols = -theName, names_to = "group", values_to = "cpms")
    mervl$marks <- NA
    mervl$marks[grep("ovctr",mervl$group)] <- "class_1"
    mervl$marks[grep("ghyvvOE",mervl$group)] <- "class_2"
    #vaue <- wilcox.test(x = mervl$cpms[grep("control",mervl$marks)],y = mervl$cpms[grep("treatss",mervl$marks)], alternative = "great")
    #print(vaue$p.value)
    summary(mervl$cpms[grep("class_1",mervl$marks)])
    summary(mervl$cpms[grep("class_2",mervl$marks)])
    #
    mervl$marks <- factor(mervl$marks, levels = c("class_1","class_2"))
    p_001 <- ggplot(data = mervl,aes(x = marks, y = log2(cpms+0.1))) +
      geom_boxplot(aes(fill = marks),outlier.alpha = 1) +
      labs(title = "MERVL-int site Expr (ghyvv OE)",
           y = "Log2 (CPM+0.1)", x = "", fill="Group") +
      geom_signif(comparisons = list(c("class_1", "class_2")),
                  textsize=6, map_signif_level=T, family = "serif",
                  test = wilcox.test, step_increase = 0.1) +
      stat_summary(geom = "point", fun = "median",shape = 22,size = 2,
                   position = position_dodge(width = 1)) +
      theme_classic() +
      scale_fill_manual(values = c("#F0A19A","#3FA0C0"),
                        labels = c("ghyvv -WT","ghyvv-OE")) +
      theme(plot.margin = margin(unit = "cm",1,1,1,1),
            legend.text = element_text(family = "serif", size = 12),
            legend.title = element_text(family = "serif", size = 12),
            axis.text.x = element_text(angle = 45,hjust = 1,family = "serif",size = 12),
            axis.text.y = element_text(family = "serif", size = 12)) +
      scale_y_continuous(limits = c(-5,12)) +
      #scale_y_continuous(limits = c(0,10),breaks = c(0,1,2,3,4,5)) +
      scale_x_discrete(labels = c("ghyvv -WT","ghyvv-OE"))
    print(p_001)
    ggsave(plot = p_001,
           filename = paste0(path,"countUniq/","TEtools.uniq.mervl-int.boxplot.ghyvvOE.pdf"),
           units = "cm", width = 20, height = 15)
  }
  shy_d(grou = grou)
}
################################################################################
#
#
#
#
#
#
#
################################################################################
#重合下ghyuu KO和ghyvv OE上调的基因和TE、
#重合下ghyuu KO和ghyvv OE下调的基因和TE、
#并制作相关的图
{
  #26、重合下ghyuu KO和ghyvv OE上调的基因
  {
    resdata <- as.data.frame(read.csv("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuu/DESeq2/ghyuu.geneOnly.csv", 
                                      header = T, stringsAsFactors = F))
    upup_01 <- resdata[which(resdata$log2FoldChange > log2(1.5) & resdata$padj < 0.05),"X"]
    upup_01 <- unique(upup_01)
    resdata <- as.data.frame(read.csv("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuughyvv/DESeq2/ghyvv.geneOnly.csv", 
                                      header = T, stringsAsFactors = F))
    upup_02 <- resdata[which(resdata$log2FoldChange > log2(1.5) & resdata$padj < 0.05),"X"]
    upup_02 <- unique(upup_02)
    #
    #
    tmp1 <- list()
    tmp1[["a"]] <- upup_01
    tmp1[["b"]] <- upup_02
    #gggg <- intersect(upup_01,upup_02)
    p_001 <- VennDiagram::venn.diagram(x = tmp1,
                                       filename = NULL,
                                       main = "\ngene",
                                       main.fontfamily = "serif",
                                       sub  = "\n\nghyuu_KO-upup___ghyvv_OE-upup\n",
                                       sub.fontfamily = "serif",
                                       cex = 1.5, main.cex = 1.5, sub.cex = 1,
                                       print.mode = "raw",#c("percent", "raw"),
                                       category.names = c("ghyuu_KO-upup",
                                                          "\n\n\n\n\nghyvv_OE-upup"),
                                       cat.cex = 1.5,cat.dist = -0.05,
                                       fill=c("#eaafc8","#00b4db"),
                                       units = "cm", height = 12, width = 12)
    pdf(file = paste0(path,"PLOT/","venn.gene.ghyuu.ghyvv.upup",".pdf"), width = 8, height = 8)
    grid.draw(p_001)
    dev.off()
  }
  ################################################################################
  #27、重合下ghyuu KO和ghyvv OE下调的基因
  {
    resdata <- as.data.frame(read.csv("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuu/DESeq2/ghyuu.geneOnly.csv", 
                                      header = T, stringsAsFactors = F))
    down_01 <- resdata[which(resdata$log2FoldChange < -log2(1.5) & resdata$padj < 0.05),"X"]
    down_01 <- unique(down_01)
    resdata <- as.data.frame(read.csv("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuughyvv/DESeq2/ghyvv.geneOnly.csv", 
                                      header = T, stringsAsFactors = F))
    down_02 <- resdata[which(resdata$log2FoldChange < -log2(1.5) & resdata$padj < 0.05),"X"]
    down_02 <- unique(down_02)
    #
    #
    tmp2 <- list()
    tmp2[["a"]] <- down_01
    tmp2[["b"]] <- down_02
    p_002 <- VennDiagram::venn.diagram(x = tmp2,
                                       filename = NULL,
                                       main = "\ngene",
                                       main.fontfamily = "serif",
                                       sub  = "\n\nghyuu_KO-down___ghyvv_OE-down\n",
                                       sub.fontfamily = "serif",
                                       cex = 1.5, main.cex = 1.5, sub.cex = 1,
                                       print.mode = "raw",#c("percent", "raw"),
                                       category.names = c("ghyuu_KO-down",
                                                          "\n\n\n\n\nghyvv_OE-down"),
                                       cat.cex = 1.5,cat.dist = -0.05,
                                       fill=c("#eaafc8","#00b4db"),
                                       units = "cm", height = 12, width = 12)
    pdf(file = paste0(path,"PLOT/","venn.gene.ghyuu.ghyvv.down",".pdf"), width = 8, height = 8)
    grid.draw(p_002)
    dev.off()
  }
  ################################################################################
  #28、重合下ghyuu KO和ghyvv OE上调的TE [基因背景]
  {
    resdata <- as.data.frame(read.csv("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuu/DESeq2/ghyuu.te.GeneBackGround.csv", 
                                      header = T, stringsAsFactors = F))
    upup_01 <- resdata[which(resdata$log2FoldChange > log2(1.5) & resdata$padj < 0.05),"X"]
    upup_01 <- unique(upup_01)
    resdata <- as.data.frame(read.csv("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuughyvv/DESeq2/ghyvv.te.GeneBackGround.csv", 
                                      header = T, stringsAsFactors = F))
    upup_02 <- resdata[which(resdata$log2FoldChange > log2(1.5) & resdata$padj < 0.05),"X"]
    upup_02 <- unique(upup_02)
    #
    #
    tmp1 <- list()
    tmp1[["a"]] <- upup_01
    tmp1[["b"]] <- upup_02
    p_001 <- VennDiagram::venn.diagram(x = tmp1,
                                       filename = NULL,
                                       main = "\nTE (gene background)",
                                       main.fontfamily = "serif",
                                       sub  = "\n\nghyuu_KO-upup___ghyvv_OE-upup\n",
                                       sub.fontfamily = "serif",
                                       cex = 1.5, main.cex = 1.5, sub.cex = 1,
                                       print.mode = "raw",#c("percent", "raw"),
                                       category.names = c("ghyuu_KO-upup",
                                                          "\n\n\n\n\nghyvv_OE-upup"),
                                       cat.cex = 1.5,cat.dist = -0.05,
                                       fill=c("#eaafc8","#00b4db"),
                                       units = "cm", height = 12, width = 12)
    pdf(file = paste0(path,"PLOT/","venn.TE.geneBG.ghyuu.ghyvv.upup",".pdf"), width = 8, height = 8)
    grid.draw(p_001)
    dev.off()
  }
  ################################################################################
  #29、重合下ghyuu KO和ghyvv OE下调的TE [基因背景]
  {
    resdata <- as.data.frame(read.csv("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuu/DESeq2/ghyuu.te.GeneBackGround.csv", 
                                      header = T, stringsAsFactors = F))
    down_01 <- resdata[which(resdata$log2FoldChange < -log2(1.5) & resdata$padj < 0.05),"X"]
    down_01 <- unique(down_01)
    resdata <- as.data.frame(read.csv("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuughyvv/DESeq2/ghyvv.te.GeneBackGround.csv", 
                                      header = T, stringsAsFactors = F))
    down_02 <- resdata[which(resdata$log2FoldChange < -log2(1.5) & resdata$padj < 0.05),"X"]
    down_02 <- unique(down_02)
    #
    #
    tmp2 <- list()
    tmp2[["a"]] <- down_01
    tmp2[["b"]] <- down_02
    p_002 <- VennDiagram::venn.diagram(x = tmp2,
                                       filename = NULL,
                                       main = "\nTE (gene background)",
                                       main.fontfamily = "serif",
                                       sub  = "\n\nghyuu_KO-down___ghyvv_OE-down\n",
                                       sub.fontfamily = "serif",
                                       cex = 1.5, main.cex = 1.5, sub.cex = 1,
                                       print.mode = "raw",#c("percent", "raw"),
                                       category.names = c("ghyuu_KO-down",
                                                          "\n\n\n\n\nghyvv_OE-down"),
                                       cat.cex = 1.5,cat.dist = -0.05,
                                       fill=c("#eaafc8","#00b4db"),
                                       units = "cm", height = 12, width = 12)
    pdf(file = paste0(path,"PLOT/","venn.TE.geneBG.ghyuu.ghyvv.down",".pdf"), width = 8, height = 8)
    grid.draw(p_002)
    dev.off()
  }
  ################################################################################
  #30、重合下ghyuu KO和ghyvv OE上调的TE [没有基因背景]
  {
    resdata <- as.data.frame(read.csv("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuu/DESeq2/ghyuu.teOnly.csv", 
                                      header = T, stringsAsFactors = F))
    upup_01 <- resdata[which(resdata$log2FoldChange > log2(1.5) & resdata$padj < 0.05),"X"]
    upup_01 <- unique(upup_01)
    resdata <- as.data.frame(read.csv("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuughyvv/DESeq2/ghyvv.teOnly.csv", 
                                      header = T, stringsAsFactors = F))
    upup_02 <- resdata[which(resdata$log2FoldChange > log2(1.5) & resdata$padj < 0.05),"X"]
    upup_02 <- unique(upup_02)
    #
    #
    tmp1 <- list()
    tmp1[["a"]] <- upup_01
    tmp1[["b"]] <- upup_02
    p_001 <- VennDiagram::venn.diagram(x = tmp1,
                                       filename = NULL,
                                       main = "\nTE",
                                       main.fontfamily = "serif",
                                       sub  = "\n\nghyuu_KO-upup___ghyvv_OE-upup\n",
                                       sub.fontfamily = "serif",
                                       cex = 1.5, main.cex = 1.5, sub.cex = 1,
                                       print.mode = "raw",#c("percent", "raw"),
                                       category.names = c("ghyuu_KO-upup",
                                                          "\n\n\n\n\nghyvv_OE-upup"),
                                       cat.cex = 1.5,cat.dist = -0.05,
                                       fill=c("#eaafc8","#00b4db"),
                                       units = "cm", height = 12, width = 12)
    pdf(file = paste0(path,"PLOT/","venn.TE.ghyuu.ghyvv.upup",".pdf"), width = 8, height = 8)
    grid.draw(p_001)
    dev.off()
  }
  ################################################################################
  #31、重合下ghyuu KO和ghyvv OE下调的TE [没有基因背景]
  {
    resdata <- as.data.frame(read.csv("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuu/DESeq2/ghyuu.teOnly.csv", 
                                      header = T, stringsAsFactors = F))
    down_01 <- resdata[which(resdata$log2FoldChange < -log2(1.5) & resdata$padj < 0.05),"X"]
    down_01 <- unique(down_01)
    resdata <- as.data.frame(read.csv("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuughyvv/DESeq2/ghyvv.teOnly.csv", 
                                      header = T, stringsAsFactors = F))
    down_02 <- resdata[which(resdata$log2FoldChange < -log2(1.5) & resdata$padj < 0.05),"X"]
    down_02 <- unique(down_02)
    #
    #
    tmp2 <- list()
    tmp2[["a"]] <- down_01
    tmp2[["b"]] <- down_02
    p_002 <- VennDiagram::venn.diagram(x = tmp2,
                                       filename = NULL,
                                       main = "\nTE",
                                       main.fontfamily = "serif",
                                       sub  = "\n\nghyuu_KO-down___ghyvv_OE-down\n",
                                       sub.fontfamily = "serif",
                                       cex = 1.5, main.cex = 1.5, sub.cex = 1,
                                       print.mode = "raw",#c("percent", "raw"),
                                       category.names = c("ghyuu_KO-down",
                                                          "\n\n\n\n\nghyvv_OE-down"),
                                       cat.cex = 1.5,cat.dist = -0.05,
                                       fill=c("#eaafc8","#00b4db"),
                                       units = "cm", height = 12, width = 12)
    pdf(file = paste0(path,"PLOT/","venn.TE.ghyuu.ghyvv.down",".pdf"), width = 8, height = 8)
    grid.draw(p_002)
    dev.off()
  }
}
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
#
#
#
#
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
  shy_diff <- function(mark, geneFile, Class, VS, cuts){
    tmp1 <- read.table(geneFile, header = T, check.names = F, stringsAsFactors = F, row.names = 1)[,cuts]
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
           cuts = c(1,2,5,6),
           geneFile = paste0(path, "countFeature/genes-2.txt"),
           Class = factor(c("OE","OE","WT","WT")), VS = c("Class","OE","WT"))
  shy_diff(mark = "ghyvv", 
           cuts = c(3,4,5,6),
           geneFile = paste0(path, "countFeature/genes-2.txt"),
           Class = factor(c("OE","OE","WT","WT")), VS = c("Class","OE","WT"))
}
################################################################################
#03、重合下ghyuu KO和ghyvv OE上调的基因和TE、重合下ghyuu KO和ghyvv OE下调的基因和TE
{
  #重合下ghyuu KO和ghyvv OE上调的基因, 做富集分析 [富集结果不好] [×]
  {
    #KEGG
    {
      ovgene <- intersect(upup_01,upup_02)
      geyy <- na.omit(mapIds(x=org.Mm.eg.db, keys=unique(ovgene), keytype="SYMBOL", column="ENTREZID"))
      kegg <- clusterProfiler::enrichKEGG(gene = geyy, organism = "mmu", keyType = "kegg", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
      kegg <- DOSE::setReadable(kegg, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
      kegg <- kegg@result
      kegg$Description <- sapply(str_split(kegg$Description, pattern = " - Mus"), "[", 1)
      kta1 <- kegg
      openxlsx::write.xlsx(x = kta1, file = paste0(path,"PLOT/venn/venn.ghyvvOEghyuuKO.upup_179.kegg.xlsx"))
      #
      kta1$Description <- factor(kta1$Description, levels = rev(kta1$Description))
      p_01 <- ggplot(kta1[1:10,])+
        geom_bar(aes(x=-log10(pvalue), y=Description), stat = "identity", fill="#1ba784", alpha=1.0) +
        theme_few() +
        scale_x_continuous(breaks = seq(8)) +
        geom_vline(aes(xintercept = 1), color = "white", linewidth  = 0.8) +
        geom_vline(aes(xintercept = 2), color = "white", linewidth  = 0.8) +
        geom_vline(aes(xintercept = 3), color = "white", linewidth  = 0.8) +
        geom_vline(aes(xintercept = 4), color = "white", linewidth  = 0.8) +
        labs(y="",title = "this is title--upup") +
        theme(text = element_text(size = 12,family = "serif",color = "black"),
              panel.border = element_rect(linewidth = 1),
              plot.title = element_text(size = 12,family = "serif",color = "#1ba784",face = "bold"),
              axis.title = element_text(size = 12,family = "serif",color = "#1ba784"),
              axis.text  = element_text(size = 12,family = "serif",color = "#1ba784"))
      p_01
    }
    #GO BP
    {
      ovgene <- intersect(upup_01,upup_02)
      geyy <- na.omit(mapIds(x=org.Mm.eg.db,
                             keys=unique(ovgene), keytype="SYMBOL", column="ENTREZID"))
      gogo <- clusterProfiler::enrichGO(gene = geyy, 
                                        OrgDb = org.Mm.eg.db, keyType = "ENTREZID", 
                                        ont = "BP",pvalueCutoff = 0.05,qvalueCutoff = 0.05, readable = T)
      gogo <- gogo@result
      gta1 <- gogo
      openxlsx::write.xlsx(x = gta1, file = paste0(path,"PLOT/venn/venn.ghyvvOEghyuuKO.upup_179.go.bp.xlsx"))
      #
      gta1$Description <- factor(gta1$Description, levels = rev(gta1$Description))
      p_02 <- ggplot(gta1[1:10,])+
        geom_bar(aes(x=-log10(pvalue), y=Description), stat = "identity", fill="#1ba784", alpha=1.0) +
        theme_few() +
        scale_x_continuous(breaks = seq(13)) +
        geom_vline(aes(xintercept = 1), color = "white", linewidth  = 0.8) +
        geom_vline(aes(xintercept = 2), color = "white", linewidth  = 0.8) +
        geom_vline(aes(xintercept = 3), color = "white", linewidth  = 0.8) +
        geom_vline(aes(xintercept = 4), color = "white", linewidth  = 0.8) +
        geom_vline(aes(xintercept = 5), color = "white", linewidth  = 0.8) +
        geom_vline(aes(xintercept = 6), color = "white", linewidth  = 0.8) +
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
      p_02
    }
  }
  #重合下ghyuu KO和ghyvv OE上调的基因
  {
    resdata <- as.data.frame(read.csv("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuu/DESeq2/ghyuu.geneOnly.csv", 
                                      header = T, stringsAsFactors = F))
    upup_01 <- resdata[which(resdata$log2FoldChange > log2(1.5) & resdata$padj < 0.05),"X"]
    upup_01 <- unique(upup_01)
    resdata <- as.data.frame(read.csv("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuughyvv/DESeq2/ghyvv.geneOnly.csv", 
                                      header = T, stringsAsFactors = F))
    upup_02 <- resdata[which(resdata$log2FoldChange > log2(1.5) & resdata$padj < 0.05),"X"]
    upup_02 <- unique(upup_02)
    #
    #
    tmp1 <- list()
    tmp1[["a"]] <- upup_01
    tmp1[["b"]] <- upup_02
    p_001 <- VennDiagram::venn.diagram(x = tmp1,
                                       filename = NULL,
                                       main = "\ngene",
                                       main.fontfamily = "serif",
                                       sub  = "\n\nghyuu_KO-upup___ghyvv_OE-upup\n",
                                       sub.fontfamily = "serif",
                                       cex = 1.5, main.cex = 1.5, sub.cex = 1,
                                       print.mode = "raw",#c("percent", "raw"),
                                       category.names = c("ghyuu_KO-upup",
                                                          "\n\n\n\n\nghyvv_OE-upup"),
                                       cat.cex = 1.5,cat.dist = -0.05,
                                       fill=c("#eaafc8","#00b4db"),
                                       units = "cm", height = 12, width = 12)
    pdf(file = paste0(path,"PLOT/","venn.ghyuuKO.ghyvvOE.coUp",".pdf"), width = 8, height = 8)
    grid.draw(p_001)
    dev.off()
    #
    #
    #
    circ_1 <- intersect(upup_01, upup_02)
    circ_2 <- read.csv("/ChIP_seq_2/aaSHY/Prmt1/ask/MERVL/mt2Activation/DESeq2/MT2.geneOnly.upup.csv")[,1]
    circ_3 <- read.csv("/ChIP_seq_2/aaSHY/Prmt1/ask/MERVL/mervlKD/DESeq2/onlyGene.mervl.down.csv")[,1]
    circ_4 <- read.table("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuughyvv/PLOT/v1/Ssrp1_TableS4_Geness.txt",header = F)[,1]
    circ_5 <- read.csv("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/Aef/DESeq2/Aef.geneOnly.upup.csv")[,1]
    circ_6 <- read.csv("/ChIP_seq_2/aaSHY/Prmt1/ask/MERVL/mervlKD/DESeq2/onlyGene.mt2.down.csv")[,1]
  }
  #重合下ghyuu KO和ghyvv OE上调的基因, pheatmap [基于embryogenesis数据]
  {
    ovgene <- intersect(upup_01,upup_02)
    #aims <- ovgene
    #aims <- c("Plagl1","Arntl2","Foxn4","Ryr3","Zcchc12","Ptpre","Slc25a31","Rhobtb1","Sh3kbp1","Wls","Peg10","AU022751","Cyp11a1","Nlrp4f","Trim75","Fam107b","Trpa1","Myh6","Lrrk2","Slc5a1","Taf7l","A530040E14Rik","AC168977.1","Sp110","Hormad1","Gm1995","Neto2","Gm13119","Ubtfl1","Gm8038","Gm8332","Nes","Gm13128","Lonrf3","Usp17lb","Gm2027","Gm5788","Arhgef26","Gm5989","Tcstv3","Gm29125","Gm6189","BB287469","Gm12800","Zscan4d","Gm18765","Gm21818","Zscan4c","Pramef25","Gm21761","Usp17lc","Gm16478","Gm20767","Ccser1","Usp26","Arg2","Stk26","Rhox13","Eif4a3l2","Usp17la","AC140365.1","Gm7957","Zscan4b","Gm4782","Usp17le","Zscan4e","Eif4a3l1","Gm12183","Gm4840","Gm5039","Zscan4f","Tmem92","Zscan4−ps2","Gm5662","Gm2035","Gm16368","Gm8300","Zfp352","Zscan4a","Gm4027","Gm50180","Gm4971","Zscan4−ps1","Gm12794","Zscan4−ps3","Gm2056","Gm21037","Gm2016","Gm6494","Cd34","Ccn1","Ccn2","Fkbp6","Pramel6","Ctsc","Pramel7","Pla2g4a","4930502E18Rik","Abcb5","Specc1","Plin2","Gask1b","Hal","Gm40977","Psd4","Maats1","Tnfrsf8","Clu","Hpcal4","Catip","Myh13","Myo7a","Bin1","Cacna1s","Aqp3","Calcoco2","Phlda1","Gprc5b","Rnf128","Tinagl1","Fam124a","Muc3")
    aims <- read.csv("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuughyvv/PLOT/heat.gene.ghyuu.ghyvv.upup_122_rank-xy.csv")
    aims <- aims$geneName
    base::load(file = "/Reference/aaSHY/DatabaseFiles/exprPlotHeatmap.RData")
    tmp2 <- exprHeat[which(rownames(exprHeat) %in% aims),c(1:11)]
    colnames(tmp2) <- gsub(x = colnames(tmp2), pattern = "G01_", replacement = "")
    tmp2 <- tmp2[match(aims,rownames(tmp2)),]
    tmp2 <- tmp2[!is.na(tmp2$oocyte),]
    #write.table(x = tmp2,
    #            paste0(path,"PLOT/","heat.gene.ghyuu.ghyvv.upup_122_rank",".csv"),
    #            sep = ",", quote = F, col.names = T, row.names = T)
    #
    #
    pdf(file = paste0(path,"PLOT/","heat.gene.ghyuu.ghyvv.upup_122_rank-2",".pdf"), width = 8, height = 12)
    pheatmap(mat = tmp2,
             main="\nheatmap of 122 genes (ghyvv OE & ghyuu KO co-up)",
             scale = "row", cellwidth = 20, 
             show_rownames = F, show_colnames = T, cluster_row = F, cluster_cols = F,
             row_names_justify = "left",
             labels_col = colnames(tmp2), angle_col = "45",
             breaks = seq(-2,2,by=0.01), border = F,
             color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                        colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
    dev.off()
  }
  #重合下ghyuu KO和ghyvv OE下调的基因
  {
    resdata <- as.data.frame(read.csv("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuu/DESeq2/ghyuu.geneOnly.csv", 
                                      header = T, stringsAsFactors = F))
    down_01 <- resdata[which(resdata$log2FoldChange < -log2(1.5) & resdata$padj < 0.05),"X"]
    down_01 <- unique(down_01)
    resdata <- as.data.frame(read.csv("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuughyvv/DESeq2/ghyvv.geneOnly.csv", 
                                      header = T, stringsAsFactors = F))
    down_02 <- resdata[which(resdata$log2FoldChange < -log2(1.5) & resdata$padj < 0.05),"X"]
    down_02 <- unique(down_02)
    #
    #
    tmp2 <- list()
    tmp2[["a"]] <- down_01
    tmp2[["b"]] <- down_02
    p_002 <- VennDiagram::venn.diagram(x = tmp2,
                                       filename = NULL,
                                       main = "\ngene",
                                       main.fontfamily = "serif",
                                       sub  = "\n\nghyuu_KO-down___ghyvv_OE-down\n",
                                       sub.fontfamily = "serif",
                                       cex = 1.5, main.cex = 1.5, sub.cex = 1,
                                       print.mode = "raw",#c("percent", "raw"),
                                       category.names = c("ghyuu_KO-down",
                                                          "\n\n\n\n\nghyvv_OE-down"),
                                       cat.cex = 1.5,cat.dist = -0.05,
                                       fill=c("#eaafc8","#00b4db"),
                                       units = "cm", height = 12, width = 12)
    pdf(file = paste0(path,"PLOT/","venn.gene.ghyuu.ghyvv.down",".pdf"), width = 8, height = 8)
    grid.draw(p_002)
    dev.off()
  }
  #重合下ghyuu KO和ghyvv OE下调的基因, pheatmap [基于embryogenesis数据]
  {
    ovgene <- intersect(down_01,down_02)
    aims <- ovgene
    #aims <- c("Plagl1","Arntl2","Foxn4","Ryr3","Zcchc12","Ptpre","Slc25a31","Rhobtb1","Sh3kbp1","Wls","Peg10","AU022751","Cyp11a1","Nlrp4f","Trim75","Fam107b","Trpa1","Myh6","Lrrk2","Slc5a1","Taf7l","A530040E14Rik","AC168977.1","Sp110","Hormad1","Gm1995","Neto2","Gm13119","Ubtfl1","Gm8038","Gm8332","Nes","Gm13128","Lonrf3","Usp17lb","Gm2027","Gm5788","Arhgef26","Gm5989","Tcstv3","Gm29125","Gm6189","BB287469","Gm12800","Zscan4d","Gm18765","Gm21818","Zscan4c","Pramef25","Gm21761","Usp17lc","Gm16478","Gm20767","Ccser1","Usp26","Arg2","Stk26","Rhox13","Eif4a3l2","Usp17la","AC140365.1","Gm7957","Zscan4b","Gm4782","Usp17le","Zscan4e","Eif4a3l1","Gm12183","Gm4840","Gm5039","Zscan4f","Tmem92","Zscan4−ps2","Gm5662","Gm2035","Gm16368","Gm8300","Zfp352","Zscan4a","Gm4027","Gm50180","Gm4971","Zscan4−ps1","Gm12794","Zscan4−ps3","Gm2056","Gm21037","Gm2016","Gm6494","Cd34","Ccn1","Ccn2","Fkbp6","Pramel6","Ctsc","Pramel7","Pla2g4a","4930502E18Rik","Abcb5","Specc1","Plin2","Gask1b","Hal","Gm40977","Psd4","Maats1","Tnfrsf8","Clu","Hpcal4","Catip","Myh13","Myo7a","Bin1","Cacna1s","Aqp3","Calcoco2","Phlda1","Gprc5b","Rnf128","Tinagl1","Fam124a","Muc3")
    #aims <- read.csv("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuughyvv/PLOT/heat.gene.ghyuu.ghyvv.upup_122_rank-xy.csv")
    #aims <- aims$geneName
    base::load(file = "/Reference/aaSHY/DatabaseFiles/exprPlotHeatmap.RData")
    tmp2 <- exprHeat[which(rownames(exprHeat) %in% aims),c(1:11)]
    colnames(tmp2) <- gsub(x = colnames(tmp2), pattern = "G01_", replacement = "")
    tmp2 <- tmp2[match(aims,rownames(tmp2)),]
    tmp2 <- tmp2[!is.na(tmp2$oocyte),]
    #write.table(x = tmp2,
    #            paste0(path,"PLOT/","heat.gene.ghyuu.ghyvv.upup_122_rank",".csv"),
    #            sep = ",", quote = F, col.names = T, row.names = T)
    #
    #
    pdf(file = paste0(path,"PLOT/","heat.gene.ghyuu.ghyvv.down.14-v1",".pdf"), width = 8, height = 12)
    pheatmap(mat = tmp2,
             main="\n------------14 genes (ghyvv OE & ghyuu KO co-down)",
             scale = "row", cellwidth = 20, 
             show_rownames = T, show_colnames = T, cluster_row = F, cluster_cols = F,
             row_names_justify = "left",
             labels_col = colnames(tmp2), angle_col = "45",
             breaks = seq(-2,2,by=0.01), border = F,
             color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                        colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
    dev.off()
  }
  #统计MERVL-int旁边的LTR组成
  {
    teff <- import(inx4)
    merv <- teff[grep("MERVL-int",teff$gene_id)]
    merv <- c(flank(merv,width = 10,start = T),flank(merv,width = 10,start = F))
    merv <- reduce(merv)
    ov01 <- suppressWarnings(findOverlaps(teff,merv,ignore.strand=T))
    tef1 <- teff[unique(queryHits(ov01))]
    tef2 <- as.data.frame(table(paste0(tef1$gene_id,":",tef1$family_id,":",tef1$class_id)))
    tef2 <- tef2[order(tef2$Freq,decreasing = T),]
    View(tef2)
  }
  #重合下ghyuu KO和ghyvv OE上调的基因, 融合基因 [基于embryogenesis数据] [MERVL|MT2C|ORR1A3-int] [也不好] [×]
  {
    gtff <- import(inx3)
    gtff <- gtff[which(gtff$type=="gene")]
    teff <- import(inx4)
    aims <- teff[grep("MT2_Mm|MERVL-int|MT2C_Mm|ORR1A3-int",teff$gene_id)]
    #
    #
    bath <- "/ChIP_seq_1/aaSHY/rna2014Deng/myRun/"
    prex <- c("twoEarly","twoMid","twoLate")
    df08 <- data.frame()
    for (i in prex) {
      print(i)
      asem <- import(paste0(bath,"stringtie/",i,".merge.gtf"))
      exon <- asem[which(asem$type=="exon")]
      chm1 <- asem[which(asem$type=="exon")]
      ov_1 <- suppressWarnings(findOverlaps(chm1,aims,ignore.strand=T))
      chm1 <- unique(chm1$transcript_id[unique(queryHits(ov_1))])
      exon <- exon[which(exon$transcript_id %in% chm1)]
      ov_2 <- suppressWarnings(findOverlaps(gtff,exon,ignore.strand=T))
      gggg <- unique(gtff$gene_name[unique(queryHits(ov_2))])
      tmps <- data.frame(mark = i,
                         gggg = gggg, tete = "MERVL")
      df08 <- rbind(df08,tmps)
    }; rm(asem,exon,chm1,ov_1,ov_2,gggg,tmps,i)
    gggg <- unique(df08$gggg)
    #
    #
    #
    #
    #
    tmp1 <- list()
    tmp1[["a"]] <- circ_1
    tmp1[["b"]] <- circ_3
    tmp1[["c"]] <- unique(df08$gggg)
    p_001 <- VennDiagram::venn.diagram(x = tmp1,
                                       filename = NULL,
                                       main = "\ngene",
                                       main.fontfamily = "serif",
                                       sub  = "\n\n---------------------\n\n\n",
                                       sub.fontfamily = "serif",
                                       cex = 1.5, main.cex = 1.5, sub.cex = 1,
                                       print.mode = "raw",#c("percent", "raw"),
                                       category.names = c("ghyuuKO.ghyvvOE.coUp",
                                                          "mervlKD.down","chimera.selectTE"),
                                       cat.cex = 1.5,cat.dist = -0.05,
                                       fill=c("#eaafc8","#00b4db","green4"),
                                       units = "cm", height = 12, width = 12)
    pdf(file = paste0(path,"PLOT/","venn.gene.ghyuu.ghyvv.upup",".pdf"), width = 8, height = 8)
    grid.draw(p_001)
    dev.off()
  }
  #重合下ghyuu KO和ghyvv OE上调的基因, 融合基因 [基于MT2 activatio数据] [MERVL] [也不好] [×]
  {
    gtff <- import(inx3)
    gtff <- gtff[which(gtff$type=="gene")]
    teff <- import(inx4)
    aims <- teff[grep("MT2_Mm|MERVL-int",teff$gene_id)]
    #
    #
    bath <- "/ChIP_seq_2/aaSHY/Prmt1/ask/MERVL/mt2Activation/"
    prex <- c("activation-1","activation-3","control-2","control-3")
    df09 <- data.frame()
    for (i in prex) {
      print(i)
      asem <- import(paste0(bath,"stringtie/",i,".gtf"))
      exon <- asem[which(asem$type=="exon")]
      chm1 <- asem[which(asem$type=="exon")]
      ov_1 <- suppressWarnings(findOverlaps(chm1,aims,ignore.strand=T))
      chm1 <- unique(chm1$transcript_id[unique(queryHits(ov_1))])
      exon <- exon[which(exon$transcript_id %in% chm1)]
      ov_2 <- suppressWarnings(findOverlaps(gtff,exon,ignore.strand=T))
      gggg <- unique(gtff$gene_name[unique(queryHits(ov_2))])
      tmps <- data.frame(mark = i,
                         gggg = gggg, tete = "MERVL")
      df09 <- rbind(df09,tmps)
    }; rm(asem,exon,chm1,ov_1,ov_2,gggg,tmps,i)
    gggg <- unique(df09$gggg)
    #
    #
    #
    #
    #
    tmp1 <- list()
    tmp1[["a"]] <- circ_1
    tmp1[["b"]] <- circ_3
    tmp1[["c"]] <- unique(df09$gggg)
    p_001 <- VennDiagram::venn.diagram(x = tmp1,
                                       filename = NULL,
                                       main = "\ngene",
                                       main.fontfamily = "serif",
                                       sub  = "\n\n---------------------\n\n\n",
                                       sub.fontfamily = "serif",
                                       cex = 1.5, main.cex = 1.5, sub.cex = 1,
                                       print.mode = "raw",#c("percent", "raw"),
                                       category.names = c("ghyuuKO.ghyvvOE.coUp",
                                                          "mervlKD.down","chimera.MERVL (MT2.act)"),
                                       cat.cex = 1.5,cat.dist = -0.05,
                                       fill=c("#eaafc8","#00b4db","grey75"),
                                       units = "cm", height = 12, width = 12)
    pdf(file = paste0(path,"PLOT/","venn.gene.ghyuu.ghyvv.upup",".pdf"), width = 8, height = 8)
    grid.draw(p_001)
    dev.off()
  }
  #重合下ghyuu KO和ghyvv OE上调的基因, 融合基因 [基于embryogenesis数据] [MERVL]
  {
    ###bath <- "/ChIP_seq_1/aaSHY/rna2014Deng/myRun/"
    ###shell <- paste0(bath,"src/run_c.sh")
    ###cat("#!/bin/bash\n", file = shell)
    ###prex <- c("oocyte","zygote","twoEarly","twoMid","twoLate","four","eight","sixteen","blastoEarly","blastoMid","blastoLate")
    ###for (mark in prex) {
    ###  gtfs <- list.files(paste0(bath,"stringtie"),mark,full.names = T)
    ###  cmd_01 <- paste0("stringtie --merge -o ",bath,
    ###                   "stringtie/",mark,".merge.gtf ",
    ###                   "-G ",inx3," ",paste(gtfs, collapse = " "),"\n")
    ###  for (i in cmd_01) {cat(i, file = shell, append = T)}
    ###}
    ###print(paste0("nohup bash ",shell, " >",paste0(bath,"Log/"),basename(shell),".log"," 2>&1 &"))
    gtff <- import(inx3)
    gtff <- gtff[which(gtff$type=="gene")]
    teff <- import(inx4)
    mt2s <- teff[grep("MT2_Mm|MERVL-int",teff$gene_id)]
    orra <- teff[grep("ORR1A",teff$gene_id)]
    orr1 <- orra[grep("int",orra$gene_id, invert = T)]
    orr2 <- orra[grep("int",orra$gene_id, invert = F)]
    #
    #
    bath <- "/ChIP_seq_1/aaSHY/rna2014Deng/myRun/"
    prex <- c("oocyte","zygote","twoEarly","twoMid","twoLate","four","eight","sixteen","blastoEarly","blastoMid","blastoLate")
    df01 <- data.frame()
    for (i in prex) {
      print(i)
      asem <- import(paste0(bath,"stringtie/",i,".merge.gtf"))
      exon <- asem[which(asem$type=="exon")]
      chm1 <- asem[which(asem$type=="exon")]
      ov_1 <- suppressWarnings(findOverlaps(chm1,mt2s,ignore.strand=T))
      chm1 <- unique(chm1$transcript_id[unique(queryHits(ov_1))])
      exon <- exon[which(exon$transcript_id %in% chm1)]
      ov_2 <- suppressWarnings(findOverlaps(gtff,exon,ignore.strand=T))
      gggg <- unique(gtff$gene_name[unique(queryHits(ov_2))])
      tmps <- data.frame(mark = i,
                         gggg = gggg, tete = "MERVL")
      df01 <- rbind(df01,tmps)
    }; rm(asem,exon,chm1,ov_1,ov_2,gggg,tmps,i)
    df02 <- data.frame()
    for (i in prex) {
      print(i)
      asem <- import(paste0(bath,"stringtie/",i,".merge.gtf"))
      exon <- asem[which(asem$type=="exon")]
      chm1 <- asem[which(asem$type=="exon")]
      ov_1 <- suppressWarnings(findOverlaps(chm1,orr1,ignore.strand=T))
      chm1 <- unique(chm1$transcript_id[unique(queryHits(ov_1))])
      exon <- exon[which(exon$transcript_id %in% chm1)]
      ov_2 <- suppressWarnings(findOverlaps(gtff,exon,ignore.strand=T))
      gggg <- unique(gtff$gene_name[unique(queryHits(ov_2))])
      tmps <- data.frame(mark = i,
                         gggg = gggg, tete = "ORR1A_LTR")
      df02 <- rbind(df02,tmps)
    }; rm(asem,exon,chm1,ov_1,ov_2,gggg,tmps,i)
    df03 <- data.frame()
    for (i in prex) {
      print(i)
      asem <- import(paste0(bath,"stringtie/",i,".merge.gtf"))
      exon <- asem[which(asem$type=="exon")]
      chm1 <- asem[which(asem$type=="exon")]
      ov_1 <- suppressWarnings(findOverlaps(chm1,orr2,ignore.strand=T))
      chm1 <- unique(chm1$transcript_id[unique(queryHits(ov_1))])
      exon <- exon[which(exon$transcript_id %in% chm1)]
      ov_2 <- suppressWarnings(findOverlaps(gtff,exon,ignore.strand=T))
      gggg <- unique(gtff$gene_name[unique(queryHits(ov_2))])
      tmps <- data.frame(mark = i,
                         gggg = gggg, tete = "ORR1A-int")
      df03 <- rbind(df03,tmps)
    }; rm(asem,exon,chm1,ov_1,ov_2,gggg,tmps,i)
    #
    #
    #
    #
    #
    #
    tmp1 <- list()
    tmp1[["a"]] <- upup_01
    tmp1[["b"]] <- upup_02
    tmp1[["c"]] <- unique(df01$gggg)
    #tmp1[["d"]] <- read.table("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuughyvv/PLOT/Ssrp1_TableS4_Geness.txt",header = F)[,1]
    p_001 <- VennDiagram::venn.diagram(x = tmp1,
                                       filename = NULL,
                                       main = "\ngene",
                                       main.fontfamily = "serif",
                                       sub  = "\n\n---------------------\n\n\n",
                                       sub.fontfamily = "serif",
                                       cex = 1.5, main.cex = 1.5, sub.cex = 1,
                                       print.mode = "raw",#c("percent", "raw"),
                                       category.names = c("ghyuuKO.Up",
                                                          "\n\n\n\n\nghyvvOE.up","chimera.MERVL"),
                                       cat.cex = 1.5,cat.dist = -0.05,
                                       fill=c("#eaafc8","#00b4db","green4"),
                                       units = "cm", height = 12, width = 12)
    pdf(file = paste0(path,"PLOT/","venn.gene.ghyuu.ghyvv.upup",".pdf"), width = 8, height = 8)
    grid.draw(p_001)
    dev.off()
  }
  #重合下ghyuu KO和ghyvv OE上调的基因, 融合基因 [基于Ssrp1_KO_notS数据] [MERVL]
  {
    ssrp_diff <- read.csv("/ChIP_seq_2/aaSHY/Ssrp1/rnaSeq/specificSNot/DESeq2/ssrp1KO.geneOnly.upup.csv")[,1]
    #
    #
    gtff <- import(inx3)
    gtff <- gtff[which(gtff$type=="gene")]
    teff <- import(inx4)
    mt2s <- teff[grep("MT2|MERVL",teff$gene_id)]
    orra <- teff[grep("ORR1A",teff$gene_id)]
    orr1 <- orra[grep("int",orra$gene_id, invert = T)]
    orr2 <- orra[grep("int",orra$gene_id, invert = F)]
    #
    #
    bath <- "/ChIP_seq_2/aaSHY/Ssrp1/rnaSeq/specificSNot/"
    prex <- c("ssrp1.merge","TE_fusion_transcripts_list")
    prex <- c("KO-1","KO-2","WT-1","WT-2")
    df04 <- data.frame()
    for (i in prex) {
      print(i)
      asem <- import(paste0(bath,"stringtie/",i,".gtf"))
      exon <- asem[which(asem$type=="exon")]
      #tran <- asem[which(asem$type=="transcript")]
      chm1 <- asem[which(asem$type=="exon")]
      ov_1 <- suppressWarnings(findOverlaps(chm1,mt2s,ignore.strand=T))
      chm1 <- unique(chm1$transcript_id[unique(queryHits(ov_1))])
      exon <- exon[which(exon$transcript_id %in% chm1)]
      #tran <- tran[which(tran$transcript_id %in% chm1)]
      ov_2 <- suppressWarnings(findOverlaps(gtff,exon,ignore.strand=T))
      #ov_2 <- suppressWarnings(findOverlaps(gtff,tran,ignore.strand=T))
      gggg <- unique(gtff$gene_name[unique(queryHits(ov_2))])
      tmps <- data.frame(mark = i,
                         gggg = gggg, tete = "MERVL")
      df04 <- rbind(df04,tmps)
    }; rm(asem,exon,chm1,ov_1,ov_2,gggg,tmps,i)
    df05 <- data.frame()
    for (i in prex) {
      print(i)
      asem <- import(paste0(bath,"stringtie/",i,".gtf"))
      exon <- asem[which(asem$type=="exon")]
      chm1 <- asem[which(asem$type=="exon")]
      ov_1 <- suppressWarnings(findOverlaps(chm1,orr1,ignore.strand=T))
      chm1 <- unique(chm1$transcript_id[unique(queryHits(ov_1))])
      exon <- exon[which(exon$transcript_id %in% chm1)]
      ov_2 <- suppressWarnings(findOverlaps(gtff,exon,ignore.strand=T))
      gggg <- unique(gtff$gene_name[unique(queryHits(ov_2))])
      tmps <- data.frame(mark = i,
                         gggg = gggg, tete = "ORR1A_LTR")
      df05 <- rbind(df05,tmps)
    }; rm(asem,exon,chm1,ov_1,ov_2,gggg,tmps,i)
    df06 <- data.frame()
    for (i in prex) {
      print(i)
      asem <- import(paste0(bath,"stringtie/",i,".gtf"))
      exon <- asem[which(asem$type=="exon")]
      chm1 <- asem[which(asem$type=="exon")]
      ov_1 <- suppressWarnings(findOverlaps(chm1,orr2,ignore.strand=T))
      chm1 <- unique(chm1$transcript_id[unique(queryHits(ov_1))])
      exon <- exon[which(exon$transcript_id %in% chm1)]
      ov_2 <- suppressWarnings(findOverlaps(gtff,exon,ignore.strand=T))
      gggg <- unique(gtff$gene_name[unique(queryHits(ov_2))])
      tmps <- data.frame(mark = i,
                         gggg = gggg, tete = "ORR1A-int")
      df06 <- rbind(df06,tmps)
    }; rm(asem,exon,chm1,ov_1,ov_2,gggg,tmps,i)
    ###up_chim <- unique(df01$gggg[which(df01$mark==prex[2])]) #first prex
    ###up_chim <- intersect(ssrp_diff,unique(df01$gggg[which(df01$mark==prex[1])]))
    #
    #
    #
    #
    #
    tmp1 <- list()
    tmp1[["a"]] <- circ_1
    tmp1[["b"]] <- circ_3
    tmp1[["c"]] <- unique(df04$gggg)
    p_001 <- VennDiagram::venn.diagram(x = tmp1,
                                       filename = NULL,
                                       main = "\ngene",
                                       main.fontfamily = "serif",
                                       sub  = "\n\n---------------------\n\n\n",
                                       sub.fontfamily = "serif",
                                       cex = 1.5, main.cex = 1.5, sub.cex = 1,
                                       print.mode = "raw",#c("percent", "raw"),
                                       category.names = c("ghyuuKO.ghyvvOE.coUp",
                                                          "mervlKD.down","chimera.MERVL"),
                                       cat.cex = 1.5,cat.dist = -0.05,
                                       fill=c("#eaafc8","#00b4db","green4"),
                                       units = "cm", height = 12, width = 12)
    pdf(file = paste0(path,"PLOT/","venn.gene.ghyuu.ghyvv.upup",".pdf"), width = 8, height = 8)
    grid.draw(p_001)
    dev.off()
  }
  #重合下ghyuu KO和ghyvv OE上调的基因, 融合基因 [基于embryogenesis数据] [All TE & Up in Ssrp1 KO]
  {
    ssrp_diff <- read.csv("/ChIP_seq_2/aaSHY/Ssrp1/rnaSeq/specificSNot/DESeq2/ssrp1KO.geneOnly.upup.csv")[,1]
    #
    #
    gtff <- import(inx3)
    gtff <- gtff[which(gtff$type=="gene")]
    gtff <- gtff[which(gtff$gene_name %in% ssrp_diff)]
    #gtff <- gtff[which(gtff$gene_name %in% circ_4)]
    teff <- import(inx4)
    ###mt2s <- teff[grep("MT2|MERVL",teff$gene_id)]
    ###orra <- teff[grep("ORR1A",teff$gene_id)]
    ###orr1 <- orra[grep("int",orra$gene_id, invert = T)]
    ###orr2 <- orra[grep("int",orra$gene_id, invert = F)]
    #
    #
    bath <- "/ChIP_seq_1/aaSHY/rna2014Deng/myRun/"
    prex <- c("oocyte","zygote","twoEarly","twoMid","twoLate","four","eight","sixteen","blastoEarly","blastoMid","blastoLate")
    df07 <- data.frame()
    for (i in prex) {
      print(i)
      asem <- import(paste0(bath,"stringtie/",i,".merge.gtf"))
      exon <- asem[which(asem$type=="exon")]
      chm1 <- asem[which(asem$type=="exon")]
      ov_1 <- suppressWarnings(findOverlaps(chm1,teff,ignore.strand=T))
      chm1 <- unique(chm1$transcript_id[unique(queryHits(ov_1))])
      exon <- exon[which(exon$transcript_id %in% chm1)]
      ov_2 <- suppressWarnings(findOverlaps(gtff,exon,ignore.strand=T))
      gggg <- unique(gtff$gene_name[unique(queryHits(ov_2))])
      tmps <- data.frame(mark = i,
                         gggg = gggg, tete = "TE")
      df07 <- rbind(df07,tmps)
    }; rm(asem,exon,chm1,ov_1,ov_2,gggg,tmps,i)
    #
    #
    #
    #
    #
    tmp1 <- list()
    tmp1[["a"]] <- circ_1
    tmp1[["b"]] <- circ_3
    tmp1[["c"]] <- unique(df07$gggg)
    p_001 <- VennDiagram::venn.diagram(x = tmp1,
                                       filename = NULL,
                                       main = "\ngene",
                                       main.fontfamily = "serif",
                                       sub  = "\n\n---------------------\n\n\n",
                                       sub.fontfamily = "serif",
                                       cex = 1.5, main.cex = 1.5, sub.cex = 1,
                                       print.mode = "raw",#c("percent", "raw"),
                                       category.names = c("ghyuuKO.ghyvvOE.coUp",
                                                          "mervlKD.down","chimera.TE"),
                                       cat.cex = 1.5,cat.dist = -0.05,
                                       fill=c("#eaafc8","#00b4db","green4"),
                                       units = "cm", height = 12, width = 12)
    pdf(file = paste0(path,"PLOT/","venn.gene.ghyuu.ghyvv.upup",".pdf"), width = 8, height = 8)
    grid.draw(p_001)
    dev.off()
  }
  #重合下ghyuu KO和ghyvv OE上调的基因, 上一步的融合基因 (unique(df07$gggg) ~840), 以及MT2 KD Down的基因
  {
    tmp1 <- list()
    tmp1[["a"]] <- circ_1
    tmp1[["b"]] <- unique(df07$gggg)
    tmp1[["c"]] <- circ_6
    p_001 <- VennDiagram::venn.diagram(x = tmp1,
                                       filename = NULL,
                                       main = "\ngene",
                                       main.fontfamily = "serif",
                                       sub  = "\n\n-------------------------\n",
                                       sub.fontfamily = "serif",
                                       cex = 1.5, main.cex = 1.5, sub.cex = 1,
                                       print.mode = "raw",#c("percent", "raw"),
                                       category.names = c("ghyuuKO.ghyvvOE.coUp",
                                                          "chimeras.TE",
                                                          "mt2KD.down"),
                                       cat.cex = 1.5,cat.dist = -0.05,
                                       fill=c("#eaafc8","#00b4db","grey75"),
                                       units = "cm", height = 12, width = 12)
    pdf(file = paste0(path,"PLOT/","venn.ghyuuKO.ghyvvOE.coUp",".pdf"), width = 8, height = 8)
    grid.draw(p_001)
    dev.off()
  }
  #重合下ghyuu KO和ghyvv OE上调的基因, 其中多少融合基因 [换一下检测fusion基因的方法：附近50 kb]
  {
    gtff <- import(inx3)
    gtff <- gtff[which(gtff$type=="gene")]
    teff <- import(inx4)
    mt2s <- teff[grep("MT2_Mm|MERVL-int",teff$gene_id)]
    orra <- teff[grep("ORR1A",teff$gene_id)]
    orr1 <- orra[grep("int",orra$gene_id, invert = T)]
    orr2 <- orra[grep("int",orra$gene_id, invert = F)]
    #
    #
    #
    mt2s <- mt2s +50000
    ov01 <- suppressWarnings(findOverlaps(gtff,mt2s,ignore.strand=T))
    gggg <- unique(gtff$gene_name[unique(queryHits(ov01))])
    #
    #
    #
    #
    #
    #
    #
    tmp1 <- list()
    tmp1[["a"]] <- circ_1
    tmp1[["b"]] <- circ_3
    tmp1[["c"]] <- gggg
    p_001 <- VennDiagram::venn.diagram(x = tmp1,
                                       filename = NULL,
                                       main = "\ngene",
                                       main.fontfamily = "serif",
                                       sub  = "\n\nghyuuKO_ghyvvOE_coUp.mervlKD-down.chimeras\n\n\n",
                                       sub.fontfamily = "serif",
                                       cex = 1.5, main.cex = 1.5, sub.cex = 1,
                                       print.mode = "raw",#c("percent", "raw"),
                                       category.names = c("ghyuuKO_ghyvvOE_coUp",
                                                          "mervlKD-down","chimeras"),
                                       cat.cex = 1.5,cat.dist = -0.05,
                                       fill=c("#eaafc8","#00b4db","green4"),
                                       units = "cm", height = 12, width = 12)
    pdf(file = paste0(path,"PLOT/","venn.gene.ghyuu.ghyvv.upup",".pdf"), width = 8, height = 8)
    grid.draw(p_001)
    dev.off()
  }
  #重合下ghyuu KO和ghyvv OE上调的基因, + MERVL KD 下调基因, + Ssrp1 TableS4基因
  {
    tmp1 <- list()
    tmp1[["a"]] <- circ_1
    tmp1[["b"]] <- circ_3
    tmp1[["c"]] <- circ_4
    tmp1[["d"]] <- circ_5
    p_001 <- VennDiagram::venn.diagram(x = tmp1,
                                       filename = NULL,
                                       main = "\ngene",
                                       main.fontfamily = "serif",
                                       sub  = "-------------------\n\n\n",
                                       sub.fontfamily = "serif",
                                       cex = 1.5, main.cex = 1.5, sub.cex = 1,
                                       print.mode = "raw",#c("percent", "raw"),
                                       category.names = c("ghyuuKO & ghyvvOE co-up",
                                                          "\n\n\n\n\nMERVL_KD down",
                                                          "Supplementary Table S4"),
                                       cat.cex = 1.5,cat.dist = -0.05,
                                       fill=c("#eaafc8","#00b4db","grey70"),#"#FD763F","#23BAC5","#B2DBB9"
                                       units = "cm", height = 12, width = 12)
    pdf(file = paste0(path,"PLOT/","venn.ghyuuKO_ghyvvOE_coUp.mervlKD_down.TableS4",".pdf"), width = 8, height = 8)
    grid.draw(p_001)
    dev.off()
  }
  #对上述交集的59个基因, pheatmap可视化其在Aef KD样本组的表达热图
  {
    gggg <- intersect(intersect(circ_1,circ_3),circ_4)
    df01 <- openxlsx::read.xlsx("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/Aef/count/Aef.cpm.onlyGene.xlsx")
    df01 <- df01[which(df01$theName %in% gggg),]
    rownames(df01) <- df01$theName; df01 <- df01[,-length(colnames(df01))]
    df01 <- df01[,c(3,4,1,2)]
    #
    #
    #
    pdf(file = paste0(path,"PLOT/","venn.ghyuuKO_ghyvvOE_coUp.mervlKD_down.TableS4-59Aef",".pdf"), width = 8, height = 12)
    pheatmap(mat = df01,
             main="\n\nheatmap of 59 genes\n (ghyvvOE.ghyuuKO.up__ssrp1TableS4__mervlKD.down)",
             scale = "row", cellwidth = 20, 
             show_rownames = T, show_colnames = T, cluster_row = T, cluster_cols = F,
             row_names_justify = "left",
             labels_col = colnames(df01), angle_col = "45",
             breaks = seq(-2,2,by=0.01), border = F,
             color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                        colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
    dev.off()
  }
  #关于Ssrp1 TableS4
  {
    circ_4 <- read.table("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuughyvv/PLOT/Ssrp1_TableS4_Geness.txt",header = F)[,1]
    tmp1 <- circ_4[1:728]
    tmp2 <- circ_4[729:895]
    ssrp <- read.table("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuughyvv/PLOT/Ssrp1_TableS4_Related.txt",
                       header = T, sep = "\t")
    tmp1 <- ssrp[which(ssrp$gene_name %in% tmp1),]
    tmp2 <- ssrp[which(ssrp$gene_name %in% tmp2),]
    #write.csv(x = tmp1, row.names = F, quote = F,
    #          "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuughyvv/PLOT/Ssrp1_TableS4_col1-728.csv")
    #write.csv(x = tmp2, row.names = F, quote = F,
    #          "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuughyvv/PLOT/Ssrp1_TableS4_col2-167.csv")
    freq_te <- as.data.frame(table(unlist(str_split(ssrp$promoter_down," "))))
    freq_te <- freq_te[order(freq_te$Freq,decreasing = T),]
    as.character(freq_te$Var1)[1:10]
  }
  #Aef: 如果不卡fold, 按p<0.05 (不是padj)算, 122, 65, 92里面分别有多少基因上调？
  set1 <- circ_1
  set2 <- intersect(circ_1,circ_3)
  set3 <- intersect(circ_1,circ_6)
  abce <- read.csv("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/Aef/DESeq2/Aef.geneOnly.csv")
  abce <- abce[which(abce$pvalue<0.05 & abce$log2FoldChange >0),]
  set1 <- intersect(set1,abce$X)
  set2 <- intersect(set2,abce$X)
  set3 <- intersect(set3,abce$X)
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
  #
  #
  #
  #
  #
  #
  #
  #
  #经过上面的试验, 决定出以下图
  {
    #figure-01 [venn-ghyuuKO.ghyvvOE.coUp-mervlKD.Down.pdf]
    {
      tmp1 <- list()
      tmp1[["a"]] <- circ_1
      tmp1[["b"]] <- circ_3
      p_001 <- VennDiagram::venn.diagram(x = tmp1,
                                         filename = NULL,
                                         main = "\ngene",
                                         main.fontfamily = "serif",
                                         sub  = "\n\n-------------------------\n",
                                         sub.fontfamily = "serif",
                                         cex = 1.5, main.cex = 1.5, sub.cex = 1,
                                         print.mode = "raw",#c("percent", "raw"),
                                         category.names = c("ghyuuKO.ghyvvOE.coUp",
                                                            "\nmervlKD.Down"),
                                         cat.cex = 1.5,cat.dist = -0.05,
                                         fill=c("#eaafc8","#00b4db"),
                                         units = "cm", height = 12, width = 12)
      pdf(file = paste0(path,"PLOT/","venn-ghyuuKO.ghyvvOE.coUp-mervlKD.Down",".pdf"), width = 8, height = 8)
      grid.draw(p_001)
      dev.off()
      #
      shyP <- function(inter, bg, righ, left) {
        a=righ+inter
        b=left+inter
        pzhi <- stats::phyper(inter-1, a, bg-a, b, lower.tail = F)
        print(format(pzhi,scientific=T))
      }
      shyP(inter = 92, righ = 30, left = 866, bg = 20000)
    }
    #
    #
    #
    #figure-02 [venn-ghyuuKO.ghyvvOE.coUp-mt2KD.Down.pdf]
    {
      tmp1 <- list()
      tmp1[["a"]] <- circ_1
      tmp1[["b"]] <- circ_6
      p_002 <- VennDiagram::venn.diagram(x = tmp1,
                                         filename = NULL,
                                         main = "\ngene",
                                         main.fontfamily = "serif",
                                         sub  = "\n\n-------------------------\n",
                                         sub.fontfamily = "serif",
                                         cex = 1.5, main.cex = 1.5, sub.cex = 1,
                                         print.mode = "raw",#c("percent", "raw"),
                                         category.names = c("ghyuuKO.ghyvvOE.coUp",
                                                            "\nmt2KD.Down"),
                                         cat.cex = 1.5,cat.dist = -0.05,
                                         fill=c("#eaafc8","#00b4db"),
                                         units = "cm", height = 12, width = 12)
      pdf(file = paste0(path,"PLOT/","venn-ghyuuKO.ghyvvOE.coUp-mt2KD.Down",".pdf"), width = 8, height = 8)
      grid.draw(p_002)
      dev.off()
      #
      shyP <- function(inter, bg, righ, left) {
        a=righ+inter
        b=left+inter
        pzhi <- stats::phyper(inter-1, a, bg-a, b, lower.tail = F)
        print(format(pzhi,scientific=T))
      }
      shyP(inter = 65, righ = 57, left = 134, bg = 20000)
    }
    #
    #
    #
    #figure-03 [在Aef KD样本的expression热图] [ghyuuKO.ghyvvOE.coUp__92__mervlKD.Down]
    {
      gggg <- intersect(circ_1,circ_3)
      hot1 <- openxlsx::read.xlsx("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/Aef/count/Aef.cpm.onlyGene.xlsx")
      hot1 <- hot1[which(hot1$theName %in% gggg),]
      rownames(hot1) <- hot1$theName; hot1 <- hot1[,-length(colnames(hot1))]
      hot1 <- hot1[,c(3,4,1,2)]
      colnames(hot1) <- c("WT-1","WT-2","Aef-KD-1","Aef-KD-2")
      #
      #
      #
      pdf(file = paste0(path,"PLOT/","venn-ghyuuKO.ghyvvOE.coUp-mervlKD.Down_92Aef",".pdf"), width = 8, height = 12)
      pheatmap(mat = hot1,
               main="\n\nheatmap of 92 genes\n (ghyvvOE.ghyuuKO.coUp___mervlKD.Down)",
               scale = "row", cellwidth = 20, 
               show_rownames = T, show_colnames = T, cluster_row = T, cluster_cols = F,
               row_names_justify = "left",
               labels_col = colnames(hot1), angle_col = "45",
               breaks = seq(-2,2,by=0.01), border = F,
               color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                          colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
      dev.off()
    }
    #
    #
    #
    #figure-04 [在Aef KD样本的expression热图] [ghyuuKO.ghyvvOE.coUp__65__mt2KD.Down]
    {
      gggg <- intersect(circ_1,circ_6)
      hot1 <- openxlsx::read.xlsx("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/Aef/count/Aef.cpm.onlyGene.xlsx")
      hot1 <- hot1[which(hot1$theName %in% gggg),]
      rownames(hot1) <- hot1$theName; hot1 <- hot1[,-length(colnames(hot1))]
      hot1 <- hot1[,c(3,4,1,2)]
      colnames(hot1) <- c("WT-1","WT-2","Aef-KD-1","Aef-KD-2")
      #
      #
      #
      pdf(file = paste0(path,"PLOT/","venn-ghyuuKO.ghyvvOE.coUp-mt2KD.Down_65Aef",".pdf"), width = 8, height = 12)
      pheatmap(mat = hot1,
               main="\n\nheatmap of 65 genes\n (ghyvvOE.ghyuuKO.coUp___mt2KD.Down)",
               scale = "row", cellwidth = 20, 
               show_rownames = T, show_colnames = T, cluster_row = T, cluster_cols = F,
               row_names_justify = "left",
               labels_col = colnames(hot1), angle_col = "45",
               breaks = seq(-2,2,by=0.01), border = F,
               color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                          colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
      dev.off()
    }
    #
    #
    #
    #figure-05 [在Aef KD样本的expression热图] [ghyuuKO.ghyvvOE.coUp__122__] [又不要了]
    {
      gggg <- circ_1
      hot1 <- openxlsx::read.xlsx("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/Aef/count/Aef.cpm.onlyGene.xlsx")
      hot1 <- hot1[which(hot1$theName %in% gggg),]
      rownames(hot1) <- hot1$theName; hot1 <- hot1[,-length(colnames(hot1))]
      hot1 <- hot1[,c(3,4,1,2)]
      colnames(hot1) <- c("WT-1","WT-2","Aef-KD-1","Aef-KD-2")
      #
      #
      #
      pdf(file = paste0(path,"PLOT/","venn-ghyuuKO.ghyvvOE.coUp_122Aef",".pdf"), width = 8, height = 12)
      pheatmap(mat = hot1,
               main="\n\nheatmap of 122 genes\n (ghyvvOE.ghyuuKO.coUp)",
               scale = "row", cellwidth = 20, 
               show_rownames = T, show_colnames = T, cluster_row = T, cluster_cols = F,
               row_names_justify = "left",
               labels_col = colnames(hot1), angle_col = "45",
               breaks = seq(-2,2,by=0.01), border = F,
               color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                          colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
      dev.off()
    }
  }
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
  #
  #
  #
  #
  #
  #
  #
  #
  ################################################################################
  #重合下ghyuu KO和ghyvv OE下调的基因
  {
    resdata <- as.data.frame(read.csv("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuu/DESeq2/ghyuu.geneOnly.csv", 
                                      header = T, stringsAsFactors = F))
    down_01 <- resdata[which(resdata$log2FoldChange < -log2(1.5) & resdata$padj < 0.05),"X"]
    down_01 <- unique(down_01)
    resdata <- as.data.frame(read.csv("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuughyvv/DESeq2/ghyvv.geneOnly.csv", 
                                      header = T, stringsAsFactors = F))
    down_02 <- resdata[which(resdata$log2FoldChange < -log2(1.5) & resdata$padj < 0.05),"X"]
    down_02 <- unique(down_02)
    #
    #
    tmp2 <- list()
    tmp2[["a"]] <- down_01
    tmp2[["b"]] <- down_02
    p_002 <- VennDiagram::venn.diagram(x = tmp2,
                                       filename = NULL,
                                       main = "\ngene",
                                       main.fontfamily = "serif",
                                       sub  = "\n\nghyuu_KO-down___ghyvv_OE-down\n",
                                       sub.fontfamily = "serif",
                                       cex = 1.5, main.cex = 1.5, sub.cex = 1,
                                       print.mode = "raw",#c("percent", "raw"),
                                       category.names = c("ghyuu_KO-down",
                                                          "\n\n\n\n\nghyvv_OE-down"),
                                       cat.cex = 1.5,cat.dist = -0.05,
                                       fill=c("#eaafc8","#00b4db"),
                                       units = "cm", height = 12, width = 12)
    pdf(file = paste0(path,"PLOT/","venn.gene.ghyuu.ghyvv.down",".pdf"), width = 8, height = 8)
    grid.draw(p_002)
    dev.off()
  }
  ################################################################################
  #重合下ghyuu KO和ghyvv OE上调的TE [没有基因背景]
  {
    resdata <- as.data.frame(read.csv("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuu/DESeq2/ghyuu.teOnly.csv", 
                                      header = T, stringsAsFactors = F))
    upup_01 <- resdata[which(resdata$log2FoldChange > log2(1.5) & resdata$padj < 0.05),"X"]
    upup_01 <- unique(upup_01)
    resdata <- as.data.frame(read.csv("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuughyvv/DESeq2/ghyvv.teOnly.csv", 
                                      header = T, stringsAsFactors = F))
    upup_02 <- resdata[which(resdata$log2FoldChange > log2(1.5) & resdata$padj < 0.05),"X"]
    upup_02 <- unique(upup_02)
    #
    #
    tmp1 <- list()
    tmp1[["a"]] <- upup_01
    tmp1[["b"]] <- upup_02
    p_001 <- VennDiagram::venn.diagram(x = tmp1,
                                       filename = NULL,
                                       main = "\nTE",
                                       main.fontfamily = "serif",
                                       sub  = "\n\nghyuu_KO-upup___ghyvv_OE-upup\n",
                                       sub.fontfamily = "serif",
                                       cex = 1.5, main.cex = 1.5, sub.cex = 1,
                                       print.mode = "raw",#c("percent", "raw"),
                                       category.names = c("ghyuu_KO-upup",
                                                          "\n\n\n\n\nghyvv_OE-upup"),
                                       cat.cex = 1.5,cat.dist = -0.05,
                                       fill=c("#eaafc8","#00b4db"),
                                       units = "cm", height = 12, width = 12)
    pdf(file = paste0(path,"PLOT/","venn.TE.ghyuu.ghyvv.upup",".pdf"), width = 8, height = 8)
    grid.draw(p_001)
    dev.off()
  }
  ################################################################################
  #重合下ghyuu KO和ghyvv OE下调的TE [没有基因背景]
  {
    resdata <- as.data.frame(read.csv("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuu/DESeq2/ghyuu.teOnly.csv", 
                                      header = T, stringsAsFactors = F))
    down_01 <- resdata[which(resdata$log2FoldChange < -log2(1.5) & resdata$padj < 0.05),"X"]
    down_01 <- unique(down_01)
    resdata <- as.data.frame(read.csv("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuughyvv/DESeq2/ghyvv.teOnly.csv", 
                                      header = T, stringsAsFactors = F))
    down_02 <- resdata[which(resdata$log2FoldChange < -log2(1.5) & resdata$padj < 0.05),"X"]
    down_02 <- unique(down_02)
    #
    #
    tmp2 <- list()
    tmp2[["a"]] <- down_01
    tmp2[["b"]] <- down_02
    p_002 <- VennDiagram::venn.diagram(x = tmp2,
                                       filename = NULL,
                                       main = "\nTE",
                                       main.fontfamily = "serif",
                                       sub  = "\n\nghyuu_KO-down___ghyvv_OE-down\n",
                                       sub.fontfamily = "serif",
                                       cex = 1.5, main.cex = 1.5, sub.cex = 1,
                                       print.mode = "raw",#c("percent", "raw"),
                                       category.names = c("ghyuu_KO-down",
                                                          "\n\n\n\n\nghyvv_OE-down"),
                                       cat.cex = 1.5,cat.dist = -0.05,
                                       fill=c("#eaafc8","#00b4db"),
                                       units = "cm", height = 12, width = 12)
    pdf(file = paste0(path,"PLOT/","venn.TE.ghyuu.ghyvv.down",".pdf"), width = 8, height = 8)
    grid.draw(p_002)
    dev.off()
  }
  ################################################################################
  #重合下ghyuu KO和ghyvv OE上调的TE [基因背景]
  {
    resdata <- as.data.frame(read.csv("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuu/DESeq2/ghyuu.te.GeneBackGround.csv", 
                                      header = T, stringsAsFactors = F))
    upup_01 <- resdata[which(resdata$log2FoldChange > log2(1.5) & resdata$padj < 0.05),"X"]
    upup_01 <- unique(upup_01)
    resdata <- as.data.frame(read.csv("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuughyvv/DESeq2/ghyvv.te.GeneBackGround.csv",
                                      header = T, stringsAsFactors = F))
    upup_02 <- resdata[which(resdata$log2FoldChange > log2(1.5) & resdata$padj < 0.05),"X"]
    upup_02 <- unique(upup_02)
    #
    #
    tmp1 <- list()
    tmp1[["a"]] <- upup_01
    tmp1[["b"]] <- upup_02
    p_001 <- VennDiagram::venn.diagram(x = tmp1,
                                       filename = NULL,
                                       main = "\nTE (gene Background)",
                                       main.fontfamily = "serif",
                                       sub  = "\n\nghyuu_KO-upup___ghyvv_OE-upup\n",
                                       sub.fontfamily = "serif",
                                       cex = 1.5, main.cex = 1.5, sub.cex = 1,
                                       print.mode = "raw",#c("percent", "raw"),
                                       category.names = c("ghyuu_KO-upup",
                                                          "\n\n\n\n\nghyvv_OE-upup"),
                                       cat.cex = 1.5,cat.dist = -0.05,
                                       fill=c("#eaafc8","#00b4db"),
                                       units = "cm", height = 12, width = 12)
    pdf(file = paste0(path,"PLOT/","venn.TE.ghyuu.ghyvv.upup",".pdf"), width = 8, height = 8)
    grid.draw(p_001)
    dev.off()
  }
  ################################################################################
  #重合下ghyuu KO和ghyvv OE下调的TE [基因背景]
  {
    resdata <- as.data.frame(read.csv("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuu/DESeq2/ghyuu.te.GeneBackGround.csv", 
                                      header = T, stringsAsFactors = F))
    down_01 <- resdata[which(resdata$log2FoldChange < -log2(1.5) & resdata$padj < 0.05),"X"]
    down_01 <- unique(down_01)
    resdata <- as.data.frame(read.csv("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuughyvv/DESeq2/ghyvv.te.GeneBackGround.csv",
                                      header = T, stringsAsFactors = F))
    down_02 <- resdata[which(resdata$log2FoldChange < -log2(1.5) & resdata$padj < 0.05),"X"]
    down_02 <- unique(down_02)
    #
    #
    tmp2 <- list()
    tmp2[["a"]] <- down_01
    tmp2[["b"]] <- down_02
    p_002 <- VennDiagram::venn.diagram(x = tmp2,
                                       filename = NULL,
                                       main = "\nTE (gene Background)",
                                       main.fontfamily = "serif",
                                       sub  = "\n\nghyuu_KO-down___ghyvv_OE-down\n",
                                       sub.fontfamily = "serif",
                                       cex = 1.5, main.cex = 1.5, sub.cex = 1,
                                       print.mode = "raw",#c("percent", "raw"),
                                       category.names = c("ghyuu_KO-down",
                                                          "\n\n\n\n\nghyvv_OE-down"),
                                       cat.cex = 1.5,cat.dist = -0.05,
                                       fill=c("#eaafc8","#00b4db"),
                                       units = "cm", height = 12, width = 12)
    pdf(file = paste0(path,"PLOT/","venn.TE.ghyuu.ghyvv.down",".pdf"), width = 8, height = 8)
    grid.draw(p_002)
    dev.off()
  }
}
#04、重合下ghyuu KO和ghyvv OE上调的基因、major ZGA基因
{
  resdata <- as.data.frame(read.csv("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuu/countFeature/DESeq2/ghyuu.geneOnly.csv", 
                                    header = T, stringsAsFactors = F))
  upup_01 <- resdata[which(resdata$log2FoldChange > log2(1.5) & resdata$padj < 0.05),"X"]
  upup_01 <- unique(upup_01)
  resdata <- as.data.frame(read.csv("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuughyvv/countFeature/DESeq2/ghyvv.geneOnly.csv", 
                                    header = T, stringsAsFactors = F))
  upup_02 <- resdata[which(resdata$log2FoldChange > log2(1.5) & resdata$padj < 0.05),"X"]
  upup_02 <- unique(upup_02)
  resdata <- read.table("/Reference/aaSHY/GSEAgmt/2022_majorZGA.gmt",sep = "\t",header = F)
  majoZGA <- as.data.frame(t(resdata))[,1][-c(1,2)]
  #
  #
  tmp1 <- list()
  tmp1[["a"]] <- upup_01
  tmp1[["b"]] <- upup_02
  tmp1[["c"]] <- majoZGA
  p_001 <- VennDiagram::venn.diagram(x = tmp1,
                                     filename = NULL,
                                     main = "\ngene",
                                     main.fontfamily = "serif",
                                     sub  = "\n\nghyuu_KO-upup___ghyvv_OE-upup___ZGA\n",
                                     sub.fontfamily = "serif",
                                     cex = 1.5, main.cex = 1.5, sub.cex = 1,
                                     print.mode = "raw",#c("percent", "raw"),
                                     category.names = c("ghyuu_KO-upup",
                                                        "\n\n\n\n\nghyvv_OE-upup","majorZGA"),
                                     cat.cex = 1.5,cat.dist = -0.05,
                                     fill=c("#eaafc8","#00b4db","green"),
                                     units = "cm", height = 12, width = 12)
  pdf(file = paste0(path,"PLOT/","venn.gene.ghyuu.ghyvv.upup",".pdf"), width = 8, height = 8)
  grid.draw(p_001)
  dev.off()
}
#04、重合下ghyuu KO和ghyvv OE上调的基因、major ZGA基因
{
  resdata <- as.data.frame(read.csv("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuu/DESeq2/ghyuu.geneOnly.csv", 
                                    header = T, stringsAsFactors = F))
  upup_01 <- resdata[which(resdata$log2FoldChange > log2(1.5) & resdata$padj < 0.05),"X"]
  upup_01 <- unique(upup_01)
  resdata <- as.data.frame(read.csv("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuughyvv/DESeq2/ghyvv.geneOnly.csv", 
                                    header = T, stringsAsFactors = F))
  upup_02 <- resdata[which(resdata$log2FoldChange > log2(1.5) & resdata$padj < 0.05),"X"]
  upup_02 <- unique(upup_02)
  resdata <- read.table("/Reference/aaSHY/GSEAgmt/2022_majorZGA.gmt",sep = "\t",header = F)
  majoZGA <- as.data.frame(t(resdata))[,1][-c(1,2)]
  #
  #
  tmp1 <- list()
  tmp1[["a"]] <- upup_01
  tmp1[["b"]] <- upup_02
  tmp1[["c"]] <- majoZGA
  p_001 <- VennDiagram::venn.diagram(x = tmp1,
                                     filename = NULL,
                                     main = "\ngene",
                                     main.fontfamily = "serif",
                                     sub  = "\n\nghyuu_KO-upup___ghyvv_OE-upup___ZGA\n",
                                     sub.fontfamily = "serif",
                                     cex = 1.5, main.cex = 1.5, sub.cex = 1,
                                     print.mode = "raw",#c("percent", "raw"),
                                     category.names = c("ghyuu_KO-upup",
                                                        "\n\n\n\n\nghyvv_OE-upup","majorZGA"),
                                     cat.cex = 1.5,cat.dist = -0.05,
                                     fill=c("#eaafc8","#00b4db","green"),
                                     units = "cm", height = 12, width = 12)
  pdf(file = paste0(path,"PLOT/","venn.gene.ghyuu.ghyvv.upup",".pdf"), width = 8, height = 8)
  grid.draw(p_001)
  dev.off()
}
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
#
#
#
#
################################################################################
#GO term Overlap之后, 挑了10个做了可视化 [出图了]-----------------------------------------------
{
  red1 <- openxlsx::read.xlsx("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuu/PLOT/GO_BP_Group1-KO-OE.xlsx",
                              sheet = "group1_down")
  red2 <- openxlsx::read.xlsx("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuu/PLOT/GO_BP_Group1-KO-OE.xlsx",
                              sheet = "ghyuuKO_up")
  red3 <- openxlsx::read.xlsx("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuu/PLOT/GO_BP_Group1-KO-OE.xlsx",
                              sheet = "ghyvvOE_up")
  ov_1 <- intersect(intersect(red1$ID,red2$ID),red3$ID)
  red1 <- red1[which(red1$ID %in% ov_1),]
  red2 <- red2[which(red2$ID %in% ov_1),]
  red3 <- red3[which(red3$ID %in% ov_1),]
  ###openxlsx::write.xlsx(x = red1, 
  ###                     file = "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuu/PLOT/GO_BP_Group1_down.xlsx")
  ###openxlsx::write.xlsx(x = red2, 
  ###                     file = "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuu/PLOT/GO_BP_ghyuuKO_up---.xlsx")
  ###openxlsx::write.xlsx(x = red3, 
  ###                     file = "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuu/PLOT/GO_BP_ghyvvOE_up---.xlsx")
  sele <- c("GO:0061982","GO:0007292","GO:0001890","GO:0007160","GO:0090130","GO:0033674","GO:0006979","GO:0051402","GO:0051090","GO:0072089")
  red1 <- red1[match(sele,red1$ID),c(1,2,8,12)]
  red2 <- red2[match(sele,red2$ID),c(1,2,8,12)]
  red3 <- red3[match(sele,red3$ID),c(1,2,8,12)]
  red1$clas <- "group1_down"
  red2$clas <- "ghyuuKO_up"
  red3$clas <- "ghyvvOE_up"
  tmps <- rbind(red1,red2,red3)
  tmps$pvalue[which(tmps$pvalue < 1e-10)] <- 1e-10
  #
  #
  tmps$Description <- factor(tmps$Description, levels = rev(unique(tmps$Description)))
  p_001 <- ggplot(tmps) +
    geom_point(aes(x = clas, 
                   y = Description, color = -log10(pvalue), size = Count)) +
    theme_classic() +
    labs(title = paste0("title--------------------")) +
    scale_size_continuous(range = c(1,5), breaks = c(10, 20, 30, 40, 50)) +
    scale_color_gradient2(low = "white", high = "red", midpoint = 1,
                          limits = c(0,10), breaks = c(0,2.5,5,7.5,10)) +
    theme(plot.margin = margin(unit = "cm", 1,1,1,1),
          text = element_text(size = 12,family = "serif",color = "black"),
          axis.title = element_blank(),
          axis.text  = element_text(size = 12,family = "serif",color = "black"),
          axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
  ggsave(plot = p_001, 
         units = "cm", width = 20, height = 16,
         filename = "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuu/PLOT/GO_BP_Group1-KO-OE____174.pdf")
}
####################
{
  red1 <- openxlsx::read.xlsx("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuu/PLOT/GO_BP_Group1-KO-OE.xlsx",
                              sheet = "group1_up")
  red2 <- openxlsx::read.xlsx("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuu/PLOT/GO_BP_Group1-KO-OE.xlsx",
                              sheet = "ghyuuKO_down")
  red3 <- openxlsx::read.xlsx("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuu/PLOT/GO_BP_Group1-KO-OE.xlsx",
                              sheet = "ghyvvOE_down")
  ov_1 <- intersect(intersect(red1$ID,red2$ID),red3$ID)
  red1 <- red1[which(red1$ID %in% ov_1),]
  red2 <- red2[which(red2$ID %in% ov_1),]
  red3 <- red3[which(red3$ID %in% ov_1),]
  openxlsx::write.xlsx(x = red1, 
                       file = "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuu/PLOT/GO_BP_Group1_up--.xlsx")
  openxlsx::write.xlsx(x = red2, 
                       file = "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuu/PLOT/GO_BP_ghyuuKO_down-.xlsx")
  openxlsx::write.xlsx(x = red3, 
                       file = "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuu/PLOT/GO_BP_ghyvvOE_down-.xlsx")
}


















