###OE
#Write by Sun Haiayng at 2025.06.03
#01、搭环境
{
  rm(list=ls())
  for (i in c("data.table","ggridges","ggrepel","KEGG.db","scales","stringr","Biostrings","tidyr","pheatmap","enrichplot","DESeq2","patchwork","ggpubr","ggsci","ggtranscript","rtracklayer","Mfuzz","BiocManager","dplyr","IRanges","GenomicRanges","ReactomePA","ggplot2","clusterProfiler","topGO","Rgraphviz","pathview","org.Mm.eg.db")){library(i, character.only = T,quietly = T)}; rm(i)
}
#02、给索引
{
  #path = "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/thejjuuzzRY2/Chi/"
  path = "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/thejjuuzzRY2/ChiR/"
  inx1 = "/Reference/aaSHY/Ensembl99/GRCm38/INDEX/STARv279a"
  inx2 <- "/Reference/aaSHY/Ensembl99/GRCm38/fa/GRCm38.dna.primary_assembly.fa"
  inx3 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/Mus_musculus.GRCm38.99.gtf"
  inx4 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf"
  inx5 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf.locInd"
}
#03、跑命令
{
  for (i in c("src","fq","Log","stringtie/count",
              "countUniq","bamCoverage","trim","DESeq2","dumproom","bam","count","PLOT")) {
    if (dir.exists(paths = paste0(path, i))==F) {
      dir.create(path = paste0(path,i), recursive = T)
    }
  }
  shell <- paste0(path,"src/run_a.sh")
  cat("#!/bin/bash\n", file = shell)
  prex <- c("OE-1","OE-2","WT-1","WT-2")
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
  cmd_05 <- paste0("bamCoverage --normalizeUsing CPM -p 10 ",
                   "--minMappingQuality 1 --samFlagExclude 256 ",
                   "-of bigwig -bs 20 -b ",
                   paste0(path,"bam/",prex,".bam")," -o ",paste0(path,"bamCoverage/",prex,".bw"),"\n")
  cmd_06 <- paste0("TEtranscripts -t ",
                   paste(paste0(path,"bam/",prex[1:2],".bam"),collapse = " ")," -c ",
                   paste(paste0(path,"bam/",prex[3:4],".bam"),collapse = " "),
                   " --GTF ",inx3," --TE ",inx4," --sortByPos --mode multi ",
                   "--project jjuuzz --outdir ",paste0(path,"count"),"\n")
  cmd_07 <- paste0("TElocal -b ",
                   paste0(path,"bam/",prex,".bam"),
                   " --GTF ",inx3," --TE ",inx5," --sortByPos --project ",
                   paste0(path, "count/",prex),"\n")
  cmd_08 <- paste0("TEtranscripts -t ",
                   paste(paste0(path,"bam/",prex[1:2],".bam"),collapse = " ")," -c ",
                   paste(paste0(path,"bam/",prex[3:4],".bam"),collapse = " "),
                   " --GTF ",inx3," --TE ",inx4," --sortByPos --mode uniq ",
                   "--project jjuuzz --outdir ",paste0(path,"countUniq"),"\n")
  cmd_09 <- paste0("TElocal -b ",
                   paste0(path,"bam/",prex,".bam"),
                   " --GTF ",inx3," --TE ",inx5," --sortByPos --mode uniq ",
                   "--project ",paste0(path, "countUniq/",prex),"\n")
  cmd_10 <- paste0("stringtie ", paste0(path,"bam/",prex,".bam")," ",
                   "-p 20 -G ", inx3, " -o ", paste0(path,"stringtie/",prex,".gtf"),"\n")
  cmd_11 <- paste0("stringtie --merge -o ",
                   path,"stringtie/jjuuzz.merge.gtf -G ",inx3," ",
                   paste(paste0(path,"stringtie/",prex,".gtf"), collapse = " "),"\n")
  cmd_12 <- paste0("cat ", paste0(path,"stringtie/jjuuzz.merge.gtf "), 
                   "|sed 's/gene_name/gene_abcd/g; s/\\tgene_id/\\tgene_name/g' >",
                   paste0(path,"stringtie/jjuuzz.merge.R.gtf"),"\n")
  cmd_13 <- paste0("TEtranscripts ",
                   "-t ",paste(paste0(path,"bam/",prex[1:2],".bam"),collapse = " ")," ",
                   "-c ",paste(paste0(path,"bam/",prex[3:4],".bam"),collapse = " ")," ",
                   "--GTF ",path,"stringtie/jjuuzz.merge.R.gtf ",
                   "--TE ", inx4," --sortByPos --mode multi --project jjuuzz ",
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
  for (i in cmd_12) {cat(i, file = shell, append = T)}
  for (i in cmd_13) {cat(i, file = shell, append = T)}
  #
  #two-pass mode是否有影响
  cmd_14 <- paste0("#STAR --runThreadN 10 --readFilesCommand zcat ",
                   "--outMultimapperOrder Random --outSAMtype BAM SortedByCoordinate ",
                   "--outFileNamePrefix ",path,"dumproom/",prex," ",
                   "--outSAMprimaryFlag AllBestScore ",
                   "--genomeDir ",inx1," --readFilesIn ",lay1," ",lay2,"\n")
  cmd_15 <- paste0("#mv ",
                   paste0(path,"dumproom/",prex,"Aligned.sortedByCoord.out.bam "),
                   paste0(path,"dumproom/",prex,".bam"),"\n")
  cmd_16 <- paste0("#samtools index ",paste0(path,"dumproom/",prex,".bam"),"\n")
  cmd_17 <- paste0("#TEtranscripts -t ",
                   paste(paste0(path,"dumproom/",prex[1:2],".bam"),collapse = " ")," -c ",
                   paste(paste0(path,"dumproom/",prex[3:4],".bam"),collapse = " "),
                   " --GTF ",inx3," --TE ",inx4," --sortByPos --mode multi ",
                   "--project jjuuzz --outdir ",paste0(path,"dumproom"),"\n")
  cmd_18 <- paste0("#TElocal -b ",
                   paste0(path,"dumproom/",prex,".bam"),
                   " --GTF ",inx3," --TE ",inx5," --sortByPos --project ",
                   paste0(path, "dumproom/",prex),"\n")
  for (i in cmd_14) {cat(i, file = shell, append = T)}
  for (i in cmd_15) {cat(i, file = shell, append = T)}
  for (i in cmd_16) {cat(i, file = shell, append = T)}
  for (i in cmd_17) {cat(i, file = shell, append = T)}
  for (i in cmd_18) {cat(i, file = shell, append = T)}
  print(paste0("nohup bash ",shell, " >",paste0(path,"Log/"),basename(shell),".log"," 2>&1 &"))
}
################################################################################
#0A、带和不带2-pass mode, 对差异倍数的影响 [结论：几乎一模一样]
{
  #tmp1 <- paste0(path,"count/jjuuzz_sigdiff_gene_TE.txt")
  #tmp1 <- read.table(tmp1)
  #tmp2 <- paste0(path,"dumproom/jjuuzz_sigdiff_gene_TE.txt")
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
  shy_diff(mark = "jjuuzz",
           geneFile = paste0(path, "count/jjuuzz.cntTable"),
           siteFile = paste0(path, "count/jjuuzz.site.cntTable"),
           Class = factor(c("OE","OE","WT","WT")), VS = c("Class","OE","WT"))
}
################################################################################
#05、CPM Excel
{
  cnt1 <- read.table(paste0(path,"count/jjuuzz.cntTable"), header = T, sep = "\t", row.names = 1)
  prex <- c("OE-1","OE-2","WT-1","WT-2")
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
  openxlsx::write.xlsx(x = tmp1, file = paste0(path, "count/jjuuzz.cpm.onlyTE.xlsx"))
  openxlsx::write.xlsx(x = tmp2, file = paste0(path, "count/jjuuzz.cpm.onlyGene.xlsx"))
  openxlsx::write.xlsx(x = tmp3, file = paste0(path, "count/jjuuzz.cpm.allGeneTE.xlsx"))
}
################################################################################
#06、Volcano for Genes
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
      scale_x_continuous(limits = c(-5, 5), breaks = seq(-5, 5, 1)) +
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
  shy_volcano_gene(resFile  = paste0(path, "DESeq2/jjuuzz.geneOnly.csv"),
                   saveFile = paste0(path, "PLOT/jjuuzz.geneOnly.volcano.pdf"))
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
  shy_volcano_te(resFile  = paste0(path, "DESeq2/jjuuzz.te.GeneBackGround.csv"),
                 saveFile = paste0(path, "PLOT/jjuuzz.te.GeneBackGround.volcano.pdf"))
  shy_volcano_te(resFile  = paste0(path, "DESeq2/jjuuzz.teOnly.csv"),
                 saveFile = paste0(path, "PLOT/jjuuzz.teOnly.volcano.pdf"))
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
  shy_teMA(difffile = paste0(path,"DESeq2/jjuuzz.te.GeneBackGround.csv"),
           savefile = paste0(path,"PLOT/jjuuzz.te.GeneBackGround.MA.pdf"))
  shy_teMA(difffile = paste0(path,"DESeq2/jjuuzz.teOnly.csv"),
           savefile = paste0(path,"PLOT/jjuuzz.teOnly.MA.pdf"))
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
  shy_repSite(resdFile = paste0(path,"DESeq2/jjuuzz.te.site.teOnly.csv"),
              saveFile = paste0(path,"PLOT/jjuuzz.te.site.teOnly.mervl.pdf"))
}
################################################################################
#10、PCA
{
  #----------------------------------PCA-batch----------------------------------
  cnt1 <- read.table(paste0(path, "count/jjuuzz.cntTable"),
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
#11、GSEA, on 2-cell markers
{
  temp <- read.csv(paste0(path, "DESeq2/jjuuzz.geneOnly.csv"), header = T)
  temp <- temp[, c("X", "log2FoldChange")]
  temp <- temp[order(temp$log2FoldChange, decreasing = T),]
  write.table(x = temp, file = paste0(path, "DESeq2/jjuuzz.res.rnk"), quote = F, sep = "\t", col.names = F, row.names = F)
  cmd_01 <- paste0("gseapy prerank -g ","/Reference/aaSHY/GSEAgmt/TwoGenes.gmt -r ",
                   paste0(path, "DESeq2/jjuuzz.res.rnk"), " -p 10 -o ", 
                   paste0(path, "PLOT/enrichGSEA/"))#conda activate base gseapy 0.10.8
  cmd_02 <- paste0("gseapy prerank -g ",
                   "/Reference/aaSHY/GSEAgmt/2022_majorZGA.gmt -r ",
                   paste0(path, "DESeq2/jjuuzz.res.rnk "),
                   "--max-size 1000 -p 10 -o ",
                   paste0(path, "PLOT/enrichGSEA/"))#conda activate base gseapy 0.10.8
  print(cmd_01)
  print(cmd_02)
}
################################################################################
#12、Up TEs loci number
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
  shy_lociNum(mark = "jjuuzz",
              resdFile = paste0(path, "DESeq2/jjuuzz.te.site.teOnly.upup.csv"),
              saveFile = paste0(path, "PLOT/TElociNum/jjuuzz.lociNum.teOnly.upup.pdf"))
  shy_lociNum(mark = "jjuuzz",
              resdFile = paste0(path, "DESeq2/jjuuzz.te.site.teOnly.down.csv"),
              saveFile = paste0(path, "PLOT/TElociNum/jjuuzz.lociNum.teOnly.down.pdf"))
  shy_lociNum(mark = "jjuuzz",
              resdFile = paste0(path, "DESeq2/jjuuzz.te.site.GeneBackGround.upup.csv"),
              saveFile = paste0(path, "PLOT/TElociNum/jjuuzz.lociNum.GeneBackGround.upup.pdf"))
  shy_lociNum(mark = "jjuuzz",
              resdFile = paste0(path, "DESeq2/jjuuzz.te.site.GeneBackGround.down.csv"),
              saveFile = paste0(path, "PLOT/TElociNum/jjuuzz.lociNum.GeneBackGround.down.pdf"))
}
################################################################################
#13、GO/KEGG [KEGG画气泡图]
{
  #GO
  {
    resFile <- paste0(path, "DESeq2/jjuuzz.geneOnly.csv")
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
    openxlsx::write.xlsx(x = gta1, file = paste0(path,"PLOT/jjuuzz.GO.BP.upup.xlsx"))
    openxlsx::write.xlsx(x = gta2, file = paste0(path,"PLOT/jjuuzz.GO.BP.down.xlsx"))
    #
    gta1$Description[10] <- "siRNA-mediated retrotransposon\n silencing by heterochromatin formation"
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
      theme(text = element_text(size = 12,family = "serif",color = "black"),
            panel.border = element_rect(linewidth = 1),
            plot.title = element_text(size = 12,family = "serif",color = "#1ba784",face = "bold"),
            axis.title = element_text(size = 12,family = "serif",color = "#1ba784"),
            axis.text  = element_text(size = 12,family = "serif",color = "#1ba784"))
    #
    gta2$Description[1] <- "collagen-activated tyrosine kinase\n receptor signaling pathway"
    gta2$Description <- factor(gta2$Description, levels = rev(gta2$Description))
    p_004 <- ggplot(gta2[1:10,])+
      geom_bar(aes(x=-log10(pvalue), y=Description), stat = "identity", fill="#5861AC", alpha=1.0) +
      ggthemes::theme_few() +
      scale_x_continuous(breaks = seq(0,8,by=1)) +
      geom_vline(aes(xintercept = 1),  color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 2),  color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 3),  color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 4),  color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 5),  color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 6),  color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 7),  color = "white", linewidth  = 0.8) +
      labs(y="",title = "this is title--down") +
      theme(text = element_text(size = 12,family = "serif",color = "black"),
            panel.border = element_rect(linewidth = 1),
            plot.title = element_text(size = 12,family = "serif",color = "#5861AC",face = "bold"),
            axis.title = element_text(size = 12,family = "serif",color = "#5861AC"),
            axis.text  = element_text(size = 12,family = "serif",color = "#5861AC"))
    plot <- p_003 +p_004 +plot_layout(nrow = 2)
    ggsave(plot = plot,
           units = "cm", width = 18, height = 16,
           filename = paste0(path,"PLOT/jjuuzz.GO.BP.pdf"))
  }
  #KEGG
  {
    resFile <- paste0(path, "DESeq2/jjuuzz.geneOnly.csv")
    resdata <- as.data.frame(read.csv(resFile, header = T, stringsAsFactors = F))
    resdata <- resdata[grep(resdata[,1], pattern = ":", invert = T),]
    mtacup <- resdata[which(resdata$log2FoldChange > log2(1.5)  &  resdata$padj < 0.05),"X"]
    kddown <- resdata[which(resdata$log2FoldChange < -log2(1.5) &  resdata$padj < 0.05),"X"]
    for (symb in seq(2)) {
      ovgene <- unlist(list(mtacup,kddown)[symb])
      geyy <- na.omit(mapIds(x=org.Mm.eg.db, keys=unique(ovgene), keytype="SYMBOL", column="ENTREZID"))
      kegg <- clusterProfiler::enrichKEGG(gene = geyy, organism = "mmu", 
                                          keyType = "kegg", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
      kegg <- DOSE::setReadable(kegg, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
      kegg <- kegg@result
      kegg$Description <- sapply(str_split(kegg$Description, pattern = " - Mus"), "[", 1)
      assign(paste0("kta",symb), kegg[which(kegg$pvalue<0.05),])
    }
    openxlsx::write.xlsx(x = kta1, file = paste0(path,"PLOT/jjuuzz.kegg.upup.xlsx"))
    openxlsx::write.xlsx(x = kta2, file = paste0(path,"PLOT/jjuuzz.kegg.down.xlsx"))
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
           filename = paste0(path,"PLOT/jjuuzz.KEGG.pdf"))
  }
}
################################################################################
#14、计算TE site的CPM  [去除在所有样本中reads count均为0的TE site]
{
  cnt1 <- read.table(paste0(path,"countUniq/jjuuzz.siteUniq.cntTable"), header = T, sep = "\t", row.names = 1)
  colLen <- length(colnames(cnt1))
  cnt1 <- cnt1[rowSums(cnt1 == 0) < colLen, ]
  tmp1 <- cnt1[grep(x = rownames(cnt1), pattern=":", invert = F),]#TE
  tmp1 <- as.data.frame(apply(tmp1, 2, function(x){x/sum(x)*1000000}))
  tmp1$theName <- rownames(tmp1)
  openxlsx::write.xlsx(x = tmp1, file = paste0(path, "countUniq/jjuuzz.cpm.onlyTE.site.xlsx"))
}
################################################################################
#15、MERVL sites的Boxplot
{
  prex <- c("OE-1","OE-2","WT-1","WT-2")
  mervs <- openxlsx::read.xlsx(paste0(path,"countUniq/jjuuzz.cpm.onlyTE.site.xlsx"))
  colnames(mervs) <- c(prex,"theName")
  gggg <- data.frame(TE = sapply(str_split(mervs$theName,":"),"[",2),
                     family = sapply(str_split(mervs$theName,":"),"[",3),
                     class = sapply(str_split(mervs$theName,":"),"[",4))
  gggg <- as.data.frame(distinct(gggg))
  #View(gggg)
  
  
  shy_d <- function(grou="jjuuzz",mark) {
    mervl <- mervs
    mervl <- mervl[grep(mark,sapply(str_split(mervl$theName,":"),"[",2)),]
    mervl <- mervl[rowSums(mervl[,1:4]==0) < 4,]
    mervl <- tidyr::pivot_longer(data = mervl, cols = -theName, names_to = "group", values_to = "cpms")
    vaue <- wilcox.test(x = mervl$cpms[grep("OE",mervl$group)],
                        y = mervl$cpms[grep("WT",mervl$group)], alternative = "greater")
    print(vaue$p.value)
    #
    mervl$marks <- NA
    mervl$marks[grep("OE",mervl$group)] <- "Overexpression"
    mervl$marks[grep("WT",mervl$group)] <- "WildType"
    mervl$marks <- factor(mervl$marks, levels = c("WildType", "Overexpression"))
    p_001 <- ggplot(data = mervl,aes(x = marks, y = log2(cpms+0.1))) +
      geom_boxplot(aes(fill = marks),outlier.alpha = 1) +
      labs(title = paste0(mark," site Expr (",grou,") [",as.character(vaue$p.value),"]\n"),
           y = "Log2 (CPM+0.1)", x = "", fill="Group") +
      geom_signif(comparisons = list(c("Overexpression", "WildType")),
                  map_signif_level=T, family = "serif",
                  textsize=6, test = wilcox.test, step_increase = 0) +
      stat_summary(geom = "point", fun = "median",shape = 22,size = 2,
                   position = position_dodge(width = 1)) +
      theme_classic() +
      scale_fill_manual(values = c("#3FA0C0","#F0A19A"),
                        labels = c("WildType","Overexpression")) +
      theme(plot.margin = margin(unit = "cm",1,1,1,1),
            legend.text = element_text(family = "serif", size = 12),
            legend.title = element_text(family = "serif", size = 12),
            axis.text.x = element_text(angle = 45,hjust = 1,family = "serif",size = 12),
            axis.text.y = element_text(family = "serif", size = 12)) +
      #scale_y_continuous(limits = c(-5,10)) +
      #scale_y_continuous(limits = c(0,10),breaks = c(0,1,2,3,4,5)) +
      scale_x_discrete(labels = c("WildType", "Overexpression"))
    print(p_001)
    ggsave(plot = p_001,
           filename = paste0(path,"PLOT/",mark,".boxplot.pdf"),
           units = "cm", width = 20, height = 15)
  }
  for (mark in c("MERVL-int","MT2_Mm")) {
    shy_d(mark = mark)
  }
}
################################################################################
#
#
#
#
#
#------------------------------定制型需求
#
#
#
#
#
################################################################################
#01、失调基因交集
{
  set1_upup <- read.csv("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/thejjuuzzRY2/Chi/DESeq2/jjuuzz.geneOnly.upup.csv")
  set1_down <- read.csv("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/thejjuuzzRY2/Chi/DESeq2/jjuuzz.geneOnly.down.csv")
  set2_upup <- read.csv("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/thejjuuzzRY2/ChiR/DESeq2/jjuuzz.geneOnly.upup.csv")
  set2_down <- read.csv("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/thejjuuzzRY2/ChiR/DESeq2/jjuuzz.geneOnly.down.csv")
  sets <- list(set1_upup$X,set1_down$X,set2_upup$X,set2_down$X)
  tmp1_upup <- sets[c(1,3)]; intersect(set1_upup$X,set2_upup$X)
  tmp1_down <- sets[c(2,4)]; intersect(set1_down$X,set2_down$X)
  dev.off()
  ou_1 <- VennDiagram::venn.diagram(x = tmp1_down,
                                    filename = NULL,
                                    main = paste0("---------main--------"),
                                    main.fontfamily = "serif",
                                    sub  = paste0("----------sub-----------"), 
                                    sub.fontfamily = "serif",
                                    cex = 1.5, main.cex = 1.5, sub.cex = 1,
                                    print.mode = "raw",#c("percent", "raw"),
                                    category.names = c("Mut",
                                                       "Mut-Rescue"),
                                    cat.cex = 1.5,cat.dist = -0.05,
                                    fill=c("#eaafc8","#00b4db"),
                                    units = "cm", height = 12, width = 12)
  grid.draw(ou_1)
}
################################################################################
#02、样本重复性差
{
  red1 <- read.table("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/thejjuuzzRY2/Chi/count/jjuuzz.cntTable", header = T)
  red1 <- read.table("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/thejjuuzzRY2/ChiR/count/jjuuzz.cntTable", header = T)
  tmp1 <- data.frame(samp_OE_1 = red1[,2], samp_OE_2 = red1[,3])
  ggplot(tmp1) +
    geom_point(aes(x = samp_OE_1, y = samp_OE_2)) + scale_y_continuous(limits = c(0,2000)) +
    scale_x_continuous(limits = c(0,2000))
  tmp2 <- data.frame(samp_WT_1 = red1[,4], samp_WT_2 = red1[,5])
  ggplot(tmp2) +
    geom_point(aes(x = samp_WT_1, y = samp_WT_2)) + scale_y_continuous(limits = c(0,2000)) +
    scale_x_continuous(limits = c(0,2000))
  #
  #
  #
  #
  red1 <- read.xlsx("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/thejjuuzzRY2/ChiR/count/jjuuzz.cpm.allGeneTE.xlsx")
  colnames(red1) <- c("samp_OE_1","samp_OE_2","samp_WT_1","samp_WT_2","theName")
  maks <- red1[grep("MERVL-int|MT2_Mm",red1$theName),]
  ggplot() +
    labs(title = "By CPM value (experimental group)") +
    geom_point(data = red1, aes(x = samp_OE_1, y = samp_OE_2)) + 
    scale_y_continuous(limits = c(0,2000)) +
    scale_x_continuous(limits = c(0,2000)) +
    geom_point(data = maks, aes(x = samp_OE_1, y = samp_OE_2), color = "red") +
    geom_label_repel(data = maks, aes(x = samp_OE_1, y = samp_OE_2, label = theName)) +
    theme(plot.margin = margin(unit = "cm",1,1,1,1))
  ggplot() +
    labs(title = "By CPM value (control group)") +
    geom_point(data = red1, aes(x = samp_WT_1, y = samp_WT_2)) + 
    scale_y_continuous(limits = c(0,2000)) +
    scale_x_continuous(limits = c(0,2000)) +
    geom_point(data = maks, aes(x = samp_WT_1, y = samp_WT_2), color = "red") +
    geom_label_repel(data = maks, aes(x = samp_WT_1, y = samp_WT_2, label = theName)) +
    theme(plot.margin = margin(unit = "cm",1,1,1,1))
}
################################################################################
#03、共同富集条目:GOBP
{
  set1_upup <- read.xlsx("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/thejjuuzzRY2/Chi/PLOT/jjuuzz.GO.BP.upup.xlsx")
  set1_down <- read.xlsx("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/thejjuuzzRY2/Chi/PLOT/jjuuzz.GO.BP.down.xlsx")
  set2_upup <- read.xlsx("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/thejjuuzzRY2/ChiR/PLOT/jjuuzz.GO.BP.upup.xlsx")
  set2_down <- read.xlsx("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/thejjuuzzRY2/ChiR/PLOT/jjuuzz.GO.BP.down.xlsx")
  sets <- list(set1_upup$ID,set1_down$ID,set2_upup$ID,set2_down$ID)
  tmp1_upup <- sets[c(1,3)]; intersect(set1_upup$Description,set2_upup$Description)
  tmp1_down <- sets[c(2,4)]; intersect(set1_down$Description,set2_down$Description)
  dev.off()
  ou_1 <- VennDiagram::venn.diagram(x = tmp1_down,
                                    filename = NULL,
                                    main = paste0("---------main--------"),
                                    main.fontfamily = "serif",
                                    sub  = paste0("----------sub-----------"), 
                                    sub.fontfamily = "serif",
                                    cex = 1.5, main.cex = 1.5, sub.cex = 1,
                                    print.mode = "raw",#c("percent", "raw"),
                                    category.names = c("Mut",
                                                       "Mut-Rescue"),
                                    cat.cex = 1.5,cat.dist = -0.05,
                                    fill=c("#eaafc8","#00b4db"),
                                    units = "cm", height = 12, width = 12)
  grid.draw(ou_1)
}
################################################################################
#04、共同富集条目:KEGG
{
  set1_upup <- read.xlsx("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/thejjuuzzRY2/Chi/PLOT/jjuuzz.kegg.upup.xlsx")
  set1_down <- read.xlsx("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/thejjuuzzRY2/Chi/PLOT/jjuuzz.kegg.down.xlsx")
  set2_upup <- read.xlsx("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/thejjuuzzRY2/ChiR/PLOT/jjuuzz.kegg.upup.xlsx")
  set2_down <- read.xlsx("/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/thejjuuzzRY2/ChiR/PLOT/jjuuzz.kegg.down.xlsx")
  sets <- list(set1_upup$ID,set1_down$ID,set2_upup$ID,set2_down$ID)
  tmp1_upup <- sets[c(1,3)]; intersect(set1_upup$Description,set2_upup$Description)
  tmp1_down <- sets[c(2,4)]; intersect(set1_down$Description,set2_down$Description)
  dev.off()
  ou_1 <- VennDiagram::venn.diagram(x = tmp1_down,
                                    filename = NULL,
                                    main = paste0("---------main--------"),
                                    main.fontfamily = "serif",
                                    sub  = paste0("----------sub-----------"), 
                                    sub.fontfamily = "serif",
                                    cex = 1.5, main.cex = 1.5, sub.cex = 1,
                                    print.mode = "raw",#c("percent", "raw"),
                                    category.names = c("Mut",
                                                       "Mut-Rescue"),
                                    cat.cex = 1.5,cat.dist = -0.05,
                                    fill=c("#eaafc8","#00b4db"),
                                    units = "cm", height = 12, width = 12)
  grid.draw(ou_1)
}







