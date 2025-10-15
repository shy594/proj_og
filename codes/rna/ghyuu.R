#write by Sun Haiayng at 2024.08.30
#01、搭环境
{
  rm(list=ls())
  for (i in c("data.table","ggridges","ggrepel","KEGG.db","scales","stringr","Biostrings","tidyr","pheatmap","enrichplot","DESeq2","patchwork","ggpubr","ggsci","ggtranscript","rtracklayer","Mfuzz","BiocManager","dplyr","IRanges","GenomicRanges","ReactomePA","ggplot2","clusterProfiler","topGO","Rgraphviz","pathview","org.Mm.eg.db")){library(i, character.only = T,quietly = T)}; rm(i)
}
#02、给索引
{
  path = "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuu/"
  inx1 = "/Reference/aaSHY/Ensembl99/GRCm38/INDEX/STARv279a"
  inx2 <- "/Reference/aaSHY/Ensembl99/GRCm38/fa/GRCm38.dna.primary_assembly.fa"
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
  shell <- paste0(path,"src/run_e.sh")
  cat("#!/bin/bash\n", file = shell)
  prex <- c("KO-1","KO-2","WT-1","WT-2")
  lay1 <- paste0(path,"fq/",prex,"_1.fq.gz")
  lay2 <- paste0(path,"fq/",prex,"_2.fq.gz")
  cmd_01 <- paste0("trim_galore --phred33 -j 4 --fastqc ",
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
                   "--project ghyuu --outdir ",paste0(path,"count"),"\n")
  cmd_07 <- paste0("TElocal -b ",
                   paste0(path,"bam/",prex,".bam"),
                   " --GTF ",inx3," --TE ",inx5," --sortByPos --project ",
                   paste0(path, "count/",prex),"\n")
  cmd_aa <- paste0("TEtranscripts -t ",
                   paste(paste0(path,"bam/",prex[1:2],".bam"),collapse = " ")," -c ",
                   paste(paste0(path,"bam/",prex[3:4],".bam"),collapse = " "),
                   " --GTF ",inx3," --TE ",inx4," --sortByPos --mode uniq ",
                   "--project ghyuu --outdir ",paste0(path,"countUniq"),"\n")
  cmd_bb <- paste0("TElocal -b ",
                   paste0(path,"bam/",prex,".bam"),
                   " --GTF ",inx3," --TE ",inx5," --sortByPos --mode uniq ",
                   "--project ",paste0(path, "countUniq/",prex),"\n")
  cmd_08 <- paste0("stringtie ", paste0(path,"bam/",prex,".bam")," ",
                   "-p 20 -G ", inx3, " -o ", paste0(path,"stringtie/",prex,".gtf"),"\n")
  cmd_09 <- paste0("stringtie --merge -o ",
                   path,"stringtie/ghyuu.merge.gtf -G ",inx3," ",
                   paste(paste0(path,"stringtie/",prex,".gtf"), collapse = " "),"\n")
  cmd_10 <- paste0("cat ", paste0(path,"stringtie/ghyuu.merge.gtf "), 
                   "|sed 's/gene_name/gene_abcd/g; s/\\tgene_id/\\tgene_name/g' >",
                   paste0(path,"stringtie/ghyuu.merge.R.gtf"),"\n")
  cmd_11 <- paste0("TEtranscripts ",
                   "-t ",paste(paste0(path,"bam/",prex[1:2],".bam"),collapse = " ")," ",
                   "-c ",paste(paste0(path,"bam/",prex[3:4],".bam"),collapse = " ")," ",
                   "--GTF ",path,"stringtie/ghyuu.merge.R.gtf ",
                   "--TE ", inx4," --sortByPos --mode multi --project ghyuu ",
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
                   "--project ghyuu --outdir ",paste0(path,"dumpROOM"),"\n")
  for (i in cmd_12) {cat(i, file = shell, append = T)}
  for (i in cmd_13) {cat(i, file = shell, append = T)}
  for (i in cmd_14) {cat(i, file = shell, append = T)}
  for (i in cmd_15) {cat(i, file = shell, append = T)}
  print(paste0("nohup bash ",shell, " >",paste0(path,"Log/"),basename(shell),".log"," 2>&1 &"))
}
################################################################################
#0A、带和不带2-pass mode, 对差异倍数的影响
{
  #tmp1 <- paste0(path,"count/ghyuu_sigdiff_gene_TE.txt")
  #tmp1 <- read.table(tmp1)
  #tmp2 <- paste0(path,"dumpROOM/ghyuu_sigdiff_gene_TE.txt")
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
           geneFile = paste0(path, "count/ghyuu.cntTable"),
           siteFile = paste0(path, "count/ghyuu.site.cntTable"),
           Class = factor(c("KO","KO","WT","WT")), VS = c("Class","KO","WT"))
}
################################################################################
#05、CPM Excel
{
  cnt1 <- read.table(paste0(path,"count/ghyuu.cntTable"), header = T, sep = "\t", row.names = 1)
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
  openxlsx::write.xlsx(x = tmp1, file = paste0(path, "count/ghyuu.cpm.onlyTE.xlsx"))
  openxlsx::write.xlsx(x = tmp2, file = paste0(path, "count/ghyuu.cpm.onlyGene.xlsx"))
  openxlsx::write.xlsx(x = tmp3, file = paste0(path, "count/ghyuu.cpm.allGeneTE.xlsx"))
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
      theme(legend.position = c(0.22,0.9), 
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
}
################################################################################
#08、MA----------------[OK]
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
           savefile = paste0(path,"PLOT/ghyuu.te.GeneBackGround.MA.pdf"))
  shy_teMA(difffile = paste0(path,"DESeq2/ghyuu.teOnly.csv"),
           savefile = paste0(path,"PLOT/ghyuu.teOnly.MA.pdf"))
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
}
################################################################################
#10、PCA
{
  #----------------------------------PCA-batch----------------------------------
  cnt1 <- read.table(paste0(path, "count/ghyuu.cntTable"),
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
#11、GSEA, on 2-cell markers [OK]
{
  temp <- read.csv(paste0(path, "DESeq2/ghyuu.geneOnly.csv"), header = T)
  temp <- temp[, c("X", "log2FoldChange")]
  temp <- temp[order(temp$log2FoldChange, decreasing = T),]
  write.table(x = temp, file = paste0(path, "DESeq2/ghyuu.res.rnk"), quote = F, sep = "\t", col.names = F, row.names = F)
  cmd_01 <- paste0("gseapy prerank -g ","/Reference/aaSHY/GSEAgmt/TwoGenes.gmt -r ",
                   paste0(path, "DESeq2/ghyuu.res.rnk"), " -p 10 -o ", 
                   paste0(path, "PLOT/enrichGSEA/"))#conda activate base gseapy 0.10.8
  cmd_02 <- paste0("gseapy prerank -g ",
                   "/Reference/aaSHY/GSEAgmt/2022_majorZGA.gmt -r ",
                   paste0(path, "DESeq2/ghyuu.res.rnk "),
                   "--max-size 1000 -p 10 -o ",
                   paste0(path, "PLOT/enrichGSEA/"))#conda activate base gseapy 0.10.8
  print(cmd_01)
  print(cmd_02)
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
  shy_lociNum(mark = "ghyuu",
              resdFile = paste0(path, "DESeq2/ghyuu.te.site.GeneBackGround.upup.csv"),
              saveFile = paste0(path, "PLOT/TElociNum/ghyuu.lociNum.GeneBackGround.upup.pdf"))
  shy_lociNum(mark = "ghyuu",
              resdFile = paste0(path, "DESeq2/ghyuu.te.site.GeneBackGround.down.csv"),
              saveFile = paste0(path, "PLOT/TElociNum/ghyuu.lociNum.GeneBackGround.down.pdf"))
}
################################################################################
#13、GO/KEGG [KEGG画气泡图]
{
  #GO
  {
    resFile <- paste0(path, "DESeq2/ghyuu.geneOnly.csv")
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
    #View(gta1)
    openxlsx::write.xlsx(x = gta1, file = paste0(path,"PLOT/ghyuughyvv.122.GO.BP.xlsx"))
    openxlsx::write.xlsx(x = gta1, file = paste0(path,"PLOT/ghyuu.GO.BP.upup.xlsx"))
    openxlsx::write.xlsx(x = gta2, file = paste0(path,"PLOT/ghyuu.GO.BP.down.xlsx"))
    #
    gta1$Description <- factor(gta1$Description, levels = rev(gta1$Description))
    p_003 <- ggplot(gta1[1:10,])+
      geom_bar(aes(x=-log10(pvalue), y=Description), stat = "identity", fill="#1ba784", alpha=1.0) +
      ggthemes::theme_few() +
      scale_x_continuous(breaks = seq(0,11,by=1)) +
      labs(y="",title = "this is title--upup") +
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
      theme(text = element_text(size = 12,family = "serif",color = "black"),
            panel.border = element_rect(linewidth = 1),
            plot.title = element_text(size = 12,family = "serif",color = "#1ba784",face = "bold"),
            axis.title = element_text(size = 12,family = "serif",color = "#1ba784"),
            axis.text  = element_text(size = 12,family = "serif",color = "#1ba784"))
    #
    gta2$Description <- factor(gta2$Description, levels = rev(gta2$Description))
    p_004 <- ggplot(gta2[1:9,])+
      geom_bar(aes(x=-log10(pvalue), y=Description), stat = "identity", fill="#5861AC", alpha=1.0) +
      ggthemes::theme_few() +
      scale_x_continuous(breaks = seq(0,6,by=1)) +
      geom_vline(aes(xintercept = 1),  color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 2),  color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 3),  color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 4),  color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 5),  color = "white", linewidth  = 0.8) +
      labs(y="",title = "this is title--down") +
      theme(text = element_text(size = 12,family = "serif",color = "black"),
            panel.border = element_rect(linewidth = 1),
            plot.title = element_text(size = 12,family = "serif",color = "#5861AC",face = "bold"),
            axis.title = element_text(size = 12,family = "serif",color = "#5861AC"),
            axis.text  = element_text(size = 12,family = "serif",color = "#5861AC"))
    plot <- p_003 +p_004 +plot_layout(nrow = 2)
    ggsave(plot = plot,
           units = "cm", width = 18, height = 16,
           filename = paste0(path,"PLOT/ghyuu.GO.BP.pdf"))
  }
  #KEGG
  {
    resFile <- paste0(path, "DESeq2/ghyuu.geneOnly.csv")
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
    openxlsx::write.xlsx(x = kta1, file = paste0(path,"PLOT/ghyuu.kegg.upup.xlsx"))
    openxlsx::write.xlsx(x = kta2, file = paste0(path,"PLOT/ghyuu.kegg.down.xlsx"))
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
    p_006 <- ggplot(kta2[1:7,], aes(x=-log10(pvalue), y=Description))+
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
           filename = paste0(path,"PLOT/ghyuu.KEGG.pdf"))
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
#
#
#
#
################################################################################
#01、对于基因, 使用featureCounts定量
{
  shell <- paste0(path,"src/run_a.sh")
  cat("#!/bin/bash\n", file = shell)
  prex <- c("KO-1","KO-2","WT-1","WT-2")
  cmd_01 <- paste0("#featureCounts -T 20 -g 'gene_name' ",
                   "-F 'GTF' -O -C -p -P --countReadPairs -d 50 -D 1000 ",
                   "-B -a ",inx3," -o ", path,"countFeature/genes-1.txt ",
                   paste(paste0(path,"bam/",prex, ".bam"),collapse = " "),"\n")
  cmd_01 <- paste0("featureCounts -T 20 -g 'gene_name' ",
                   "-F 'GTF' ",
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
    tmp1 <- read.table(geneFile, header = T, check.names = F, stringsAsFactors = F, row.names = 1)
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
#给的序列是否在MERVL融合基因富集
#step 01：序列
{
  seq1 <- "GTCCAGGAATCAAGGGGTGGAAAACGGAATAGTTCCACTTACTATCACTCCTAGTGACCCTCTAGGAAAATTTTTGCTTCCTGTCCCCATAACTCTAGGTTCTGCTGGCCT"
  seq2 <- "AAGAGCCAAGACCTGCTGAGGTGCTGGCTGAAGGTGAAGGAAATACAGAATGGGTAGTAGAGGAAGGTAGTTATAAATACCAATTAAGGCCACGTAACCAGTTGCAGAAACGAGGA"
  quersss <- "/ChIP_seq_2/aaSHY/ghyuu/ask/fuss/query_1.fa"
  cat(">seq_1\n", file = quersss); cat(paste0(seq1,"\n"), file = quersss, append = T)
  cat(">seq_2\n", file = quersss, append = T); cat(paste0(seq2,"\n"), file = quersss, append = T)
}
#step 02：MERVL归类：不表达の + 融合的, 即multiple exon的isoform + solo isoform
{
  tf01 <- import(inx4)
  aims <- tf01[grep("MERVL-int", tf01$gene_id)]
  aims$mmak <- "Solo"
  bath <- "/ChIP_seq_1/aaSHY/rna2014Deng/myRun/"
  prex <- c("twoEarly","twoMid","twoLate")
  have <- c()
  mult <- c()
  for (i in prex) {
    print(i)
    asem <- import(paste0(bath,"stringtie/",i,".merge.gtf"))
    exon <- asem[which(asem$type=="exon")]
    chm1 <- asem[which(asem$type=="exon")]
    ov_1 <- suppressWarnings(findOverlaps(chm1,aims,ignore.strand=T))
    have <- c(have,subjectHits(ov_1))
    #
    chm2 <- chm1[which(chm1$transcript_id %in% unique(chm1$transcript_id[queryHits(ov_1)]))]
    chm2 <- as.data.frame(table(chm2$transcript_id))
    chm2 <- chm2$Var1[which(chm2$Freq > 1)]
    chm2 <- chm1[which(chm1$transcript_id %in% chm2)]
    ov_2 <- suppressWarnings(findOverlaps(chm2,aims,ignore.strand=T))
    mult <- c(mult,subjectHits(ov_2))
  }
  nono <- setdiff(seq_along(aims$source),unique(have))
  aims$mmak[nono] <- "Nono"
  aims$mmak[unique(mult)] <- "Mult"
  wr01 <- as.data.frame(aims)[,c(1,2,3,14,4,5)]
  write.table(x = wr01, file = "/ChIP_seq_2/aaSHY/ghyuu/ask/fuss/mervlInt.bed",
              quote = F, sep = "\t", col.names = F, row.names = F)
}
#step 03：blast建库
{
  cmd_01 <- paste0("bedtools getfasta -name -fi ",inx2," ",
                   "-bed /ChIP_seq_2/aaSHY/ghyuu/ask/fuss/mervlInt.bed -fo ",
                   "/ChIP_seq_2/aaSHY/ghyuu/ask/fuss/mervlInt.fa","\n")
  cmd_02 <- paste0("makeblastdb -dbtype nucl -in ",
                   "/ChIP_seq_2/aaSHY/ghyuu/ask/fuss/mervlInt.fa ",
                   "-input_type fasta -parse_seqids ",
                   "-out /ChIP_seq_2/aaSHY/ghyuu/ask/fuss/blastBase/shyDB","\n")
  cmd_03 <- paste0("blastn -num_threads 16 -task megablast -evalue 0.05 ",
                   "-query /ChIP_seq_2/aaSHY/ghyuu/ask/fuss/query_1.fa -db ",
                   "/ChIP_seq_2/aaSHY/ghyuu/ask/fuss/blastBase/shyDB ",
                   "-outfmt 6 -out ",
                   "/ChIP_seq_2/aaSHY/ghyuu/ask/fuss/query_1.out.tsv","\n")
  #
  red1 <- read.table("/ChIP_seq_2/aaSHY/ghyuu/ask/fuss/query_1.out.tsv")
  red1$clas <- sapply(str_split(red1$V2,"::"),"[",1)
  red2 <- red1 %>% 
    dplyr::group_by(V1,clas) %>%
    dplyr::summarise(cous = n(), .groups = "drop")
  ggplot(data = red2) +
    geom_bar(aes(x = clas, y = cous, fill = V1), stat = "identity", position = "dodge")
}
################################################################################
#画部分基因的表达
###画在Mouse SUPeR-seq表达折线图
{
  #用2015 SUPeR-seq的数据
  aims <- c("Aef","Exosc4","Pelo","Hbs1l","Dis3l","Dis3l2","Dis3","Skiv2l",
            "Ttc37","Smg1", "Smg5", "Smg7", "Upf1","Upf2", "Upf3a","Upf3b")
  aims <- c("ff23k", "Aef", "ghyuu", "ghyvv")
  aims <- c("Spz1","Zscan4d")
  #
  #
  #
  cn01 <- read.table("/disk5/aaSHY/nanoporeNGS/count/super_devs.cntTable", 
                     sep = "\t", header = T, 
                     row.names = 1, stringsAsFactors = F)
  prex <- paste0(sapply(str_split(colnames(cn01), pattern = "\\."), "[", 6),
                 "-",
                 sapply(str_split(colnames(cn01), pattern = "\\."), "[", 7))
  colnames(cn01) <- prex
  cn02 <- cn01[grep(rownames(cn01), pattern=":", invert = T),]#选基因？选TE？
  cn02 <- as.data.frame(apply(cn02, 2, function(x){x/sum(x)*1000000}))
  cn02$geneName <- rownames(cn02)
  #
  cn03 <- cn02
  df01 <- tidyr::gather(data = cn03, key = "sample", value = "cpm", -"geneName")
  df01$stage <- gsub("-[0-9]+","",df01$sample)
  df02 <- dplyr::reframe(.data = df01, .by = c(stage, geneName), 
                         cpmExp = mean(cpm), cpmSD = sd(cpm)/sqrt(length(cpm)))
  
  # View(df02)
  df03 <- data.frame()
  for (i in aims) {#, group = geneName
    temp <- df02[grep(df02$geneName, pattern = paste0("^", i,"$")),]
    df03 <- rbind(df03, temp)
  }
  df03$stage <- factor(x = df03$stage, 
                       levels = c("oocytes","zygotes","twos","four","eight","morula","blasto"))
  
  p_01 <- ggplot(data = df03) +
    geom_line(mapping = aes(x = stage, y = cpmExp, group=geneName), stat = "identity", color = "steelblue", linewidth = 1) +
    geom_smooth(aes(x = stage, y = cpmExp), se = T) +
    geom_point(mapping = aes(x = stage, y = cpmExp), color = "steelblue", size = 2) +
    geom_errorbar(mapping = aes(x = stage, ymin = cpmExp-cpmSD, ymax = cpmExp+cpmSD), 
                  width = 0.15, color = "black", linewidth = 0.8) +
    labs(x = "", y = "mean CPM", title = "2015 SUPeR-seq") +
    facet_wrap(.~ geneName, scales = "free",nrow = 2) +
    theme_classic() +
    theme(plot.margin = margin(2,2,2,2,unit = "cm"),
          plot.title = element_text(size = 16, color = "black", face = "bold"),
          axis.text.x = element_text(family = "serif", size = 12, angle = 60, vjust = 1, hjust = 1, colour = "black"),
          strip.text.x = element_text(size = 12, color = "red", face = "bold"),
          strip.background = element_blank(),
          panel.spacing.y = unit(0.8, "cm"))
  p_01
  ggsave(plot = p_01, 
         filename = "/ChIP_seq_2/aaSHY/ghyuu/ask/someGeneExpr/ghyuu_someGeneExpr-v2.pdf",
         units = "cm", height = 30, width = 30)
}
{
  #用2015 SUPeR-seq的数据
  aims <- c("Spz1","Zscan4d")
  #
  #
  #
  cn01 <- read.table(file = "/disk5/aaSHY/nanoporeNGS/count/super_devs.cntTable", sep = "\t", header = T, row.names = 1, stringsAsFactors = F)
  prex <- paste0(sapply(str_split(colnames(cn01), pattern = "\\."), "[", 6),
                 "-",
                 sapply(str_split(colnames(cn01), pattern = "\\."), "[", 7))
  colnames(cn01) <- prex
  cn01 <- cn01[,prex]
  cn02 <- cn01[grep(rownames(cn01), pattern=":", invert = T),]#选基因？选TE？
  cn02 <- as.data.frame(apply(cn02, 2, function(x){x/sum(x)*1000000}))
  cn02$geneName <- rownames(cn02)
  #
  cn03 <- cn02
  df01 <- tidyr::gather(data = cn03, key = "sample", value = "cpm", -"geneName")
  df01$stage <- gsub("-[0-9]+","",df01$sample)
  df02 <- dplyr::reframe(.data = df01, .by = c(stage, geneName), cpmExp = mean(cpm), cpmSD = sd(cpm)/sqrt(length(cpm)))
  df03 <- data.frame()
  for (i in aims) {#, group = geneName
    temp <- df02[grep(df02$geneName, pattern = paste0("^", i,"$")),]
    df03 <- rbind(df03, temp)
  }
  df03$stage <- factor(x = df03$stage, levels = unique(df03$stage)[c(5,7,6,3,2,4,1)])
  p_01 <- ggplot(data = df03) +
    geom_line(mapping = aes(x = stage, y = cpmExp, group=geneName,
                            color=geneName), stat = "identity", linewidth = 1) +
    geom_smooth(aes(x = stage, y = cpmExp), se = T) +
    scale_color_manual(values = c("indianred4","steelblue4")) +
    ###geom_point(mapping = aes(x = stage, y = cpmExp), color = "steelblue", size = 2) +
    ###geom_errorbar(mapping = aes(x = stage, ymin = cpmExp-cpmSD, ymax = cpmExp+cpmSD),
    ###              width = 0.15, color = "black", linewidth = 0.8) +
    labs(x = "", y = "mean CPM", title = "2015 SUPeR-seq") +
    #facet_wrap(.~ geneName, scales = "free",nrow = 4) +
    theme_classic() +
    theme(plot.margin = margin(1,1,1,1,unit = "cm"),
          plot.title = element_text(size = 16, color = "black", face = "bold"),
          axis.text.x = element_text(family = "serif", size = 12, angle = 60, vjust = 1, hjust = 1, colour = "black"),
          strip.text.x = element_text(size = 12, color = "red", face = "bold"),
          strip.background = element_blank(),
          panel.spacing.y = unit(0.8, "cm"))
  p_01
  ggsave(plot = p_01, 
         filename = "/ChIP_seq_2/aaSHY/ghyuu/ask/someGeneExpr/ghyuu_someGeneExpr-v3.pdf",
         units = "cm", height = 12, width = 14)
}
###画在PacBio Illumina-seq表达折线图
{
  #用PacBio Illumina-seq的数据
  aims <- c("Aef","Exosc4","Pelo","Hbs1l","Dis3l","Dis3l2","Dis3","Skiv2l",
            "Ttc37","Smg1", "Smg5", "Smg7", "Upf1","Upf2", "Upf3a","Upf3b")
  aims <- c("ff23k", "Aef", "ghyuu", "ghyvv")
  #
  #
  #
  cn01 <- read.table(file = "/ChIP_seq_1/aaSHY/pacbioIllumina/abRepeats/count/stage.cntTable", sep = "\t", header = T, row.names = 1, stringsAsFactors = F)
  pre1 <- c("oocyte","zygote","twos","four","eight","blasto")
  pre2 <- unique(paste0((sapply(str_split(colnames(cn01), pattern = "\\."), "[", 7)),"-",
                 sapply(str_split(colnames(cn01), pattern = "\\."), "[", 8)))
  prex <- c()
  for (i in pre1) {
    prex <- c(prex,pre2[grep(i,pre2)])
  }
  colnames(cn01) <- pre2
  cn01 <- cn01[,prex]
  cn02 <- cn01[grep(rownames(cn01), pattern=":", invert = T),]#选基因？选TE？
  cn02 <- as.data.frame(apply(cn02, 2, function(x){x/sum(x)*1000000}))
  cn02$geneName <- rownames(cn02)
  #
  cn03 <- cn02
  df01 <- tidyr::gather(data = cn03, key = "sample", value = "cpm", -"geneName")
  df01$stage <- gsub("-[0-9]+","",df01$sample)
  df02 <- dplyr::reframe(.data = df01, .by = c(stage, geneName), cpmExp = mean(cpm), cpmSD = sd(cpm)/sqrt(length(cpm)))
  df03 <- data.frame()
  for (i in aims) {#, group = geneName
    temp <- df02[grep(df02$geneName, pattern = paste0("^", i,"$")),]
    df03 <- rbind(df03, temp)
  }
  df03$stage <- factor(x = df03$stage, levels = pre1)
  p_01 <- ggplot(data = df03) +
    geom_line(mapping = aes(x = stage, y = cpmExp, group=geneName), stat = "identity", color = "steelblue", linewidth = 1) +
    geom_point(mapping = aes(x = stage, y = cpmExp), color = "steelblue", size = 1.5) +
    geom_errorbar(mapping = aes(x = stage, ymin = cpmExp-cpmSD, ymax = cpmExp+cpmSD), 
                  width = 0.15, color = "steelblue", linewidth = 0.8) +
    labs(x = "", y = "mean CPM", title = "2020__PacBio Illumina-seq") +
    facet_wrap(.~ geneName, scales = "free_y", nrow = 2) +
    theme_classic() +
    theme(plot.margin = margin(2,2,2,2,unit = "cm"),
          plot.title = element_text(size = 16, color = "black", face = "bold"),
          axis.text.x = element_text(family = "serif", size = 12, angle = 60, vjust = 1, hjust = 1, colour = "black"),
          strip.text.x = element_text(size = 12, color = "red", face = "bold"),
          strip.background = element_blank(),
          panel.spacing.y = unit(0.8, "cm"))
  p_01
  ggsave(plot = p_01, 
         filename = "/ChIP_seq_2/aaSHY/ghyuu/ask/someGeneExpr/xx.pdf",
         units = "cm", height = 15, width = 20)
}
###画在Mouse 2014 Deng et.al表达折线图 [采纳] [OK]----------------------------------------------------
{
  #用2014 Deng et.al的数据
  aims <- c("Aef","Exosc4","Pelo","Hbs1l","Dis3l","Dis3l2","Dis3","Skiv2l","Ttc37")
  aims <- c("Aef","Exosc4","Pelo","Hbs1l","Dis3l","Dis3l2","Dis3","Skiv2l",
            "Ttc37","Smg1", "Smg5", "Smg7", "Upf1","Upf2", "Upf3a","Upf3b")
  aims <- c("ff23k", "Aef", "ghyuu", "ghyvv")
  #
  #
  #
  cn01 <- read.table(file = "/ChIP_seq_1/aaSHY/rna2014Deng/myRun/count/stageALL.cntTable", sep = "\t", header = T, row.names = 1, stringsAsFactors = F)
  cn01 <- cn01[,colnames(cn01)[grep(x = colnames(cn01), pattern = "gene", invert = T)]]
  prex <- unique(sapply(str_split(colnames(cn01), pattern = "\\."), "[", 7))
  #prex <- prex[c(11,10,2,3,8,5,12,1,7,4,9,6)]
  prex <- prex[c(11,10,2,12,1,7,4)]
  prex[3] <- "twos"
  colnames(cn01) <- paste(sapply(str_split(colnames(cn01), pattern = "\\."), "[", 7),sapply(str_split(colnames(cn01), pattern = "\\."), "[", 8), sep = "-")
  colnames(cn01) <- gsub("two-","twos-",colnames(cn01))
  cn01 <- cn01[,grep(paste0(prex,collapse = "|"), colnames(cn01))]
  cn02 <- as.data.frame(apply(cn01, 2, function(x){x/sum(x)*1000000}))
  cn02$geneName <- rownames(cn02)
  cn03 <- cn02#[grep(rownames(cn02), pattern=":", invert = T),]#invert = F, 那么输出TE的CPM
  df01 <- tidyr::gather(data = cn03, key = "sample", value = "cpm", -"geneName")
  df01$stage <- sapply(str_split(df01$sample, "-"), "[", 1)
  df02 <- dplyr::reframe(.data = df01, .by = c(stage, geneName), cpmExp = mean(cpm), cpmSD = sd(cpm)/sqrt(length(cpm)))
  df03 <- data.frame()
  for (i in aims) {#, group = geneName
    temp <- df02[grep(df02$geneName, pattern = paste0("^", i,"$")),]
    df03 <- rbind(df03, temp)
  }
  df03$stage <- factor(x = df03$stage, levels = prex)
  p_01 <- ggplot(data = df03) +
    geom_line(mapping = aes(x = stage, y = cpmExp, group=geneName), stat = "identity", color = "steelblue", linewidth = 1) +
    geom_point(mapping = aes(x = stage, y = cpmExp), color = "steelblue", size = 1.5) +
    geom_errorbar(mapping = aes(x = stage, ymin = cpmExp-cpmSD, ymax = cpmExp+cpmSD), 
                  width = 0.15, color = "black", linewidth = 0.8) +
    labs(x = "", y = "mean CPM", title = "2014 Deng Early Mid Late") +
    facet_wrap(.~ geneName, scales = "free_y", nrow = 2) +
    theme_classic() +
    theme(plot.margin = margin(2,2,2,2,unit = "cm"),
          plot.title = element_text(size = 16, color = "black", face = "bold"),
          axis.text.x = element_text(family = "serif", size = 12, angle = 60, vjust = 1, hjust = 1, colour = "black"),
          strip.text.x = element_text(size = 12, color = "red", face = "bold"),
          strip.background = element_blank(),
          panel.spacing.y = unit(0.8, "cm"))
  p_01
  ggsave(plot = p_01,
         filename = "/ChIP_seq_2/aaSHY/ghyuu/ask/someGeneExpr/expr_some_stagesFromEarlyMidLate.pdf",
         units = "cm", height = 15, width = 20)
}
###画在2C+/-细胞中的表达的热图
{
  #用2C+/-细胞的数据
  aims <- c("Aef","Exosc4","Pelo","Hbs1l","Dis3l","Dis3l2","Dis3","Skiv2l","Ttc37","Dot1l")
  #
  #
  #
  red1 <- read.table("/ChIP_seq_2/aaSHY/ghyuu/publicData/D1_2C_cells/count/duxs.cntTable",
                     header = T, row.names = 1, check.names = F, stringsAsFactors = F)
  colnames(red1) <- sapply(str_split(basename(colnames(red1)),"\\."),"[",1)
  red1 <- red1[grep(rownames(red1), pattern = ":", invert = T),]
  red2 <- as.data.frame(apply(red1, 2, function(x){x/sum(x)*1000000}))
  red3 <- red2[aims,]
  pheatmap(mat = red3,
           main="\n---------2C+/-细胞------------",
           scale = "row", cellwidth = 20, 
           show_rownames = T, show_colnames = T, cluster_row = T, cluster_cols = F,
           row_names_justify = "left",
           labels_col = colnames(red3), angle_col = "45",
           breaks = seq(-2,2,by=0.01), border = F,
           color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                      colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
}















