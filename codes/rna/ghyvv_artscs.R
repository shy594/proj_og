#write by Sun Haiayng at 2025.05.08
#01、搭环境
{
  rm(list=ls())
  for (i in c("data.table","stringr","EDASeq","dplyr","tidyr","RUVSeq","ggridges","ggrepel","KEGG.db","scales","stringr","Biostrings","tidyr","pheatmap","enrichplot","DESeq2","patchwork","ggpubr","ggsci","ggtranscript","rtracklayer","Mfuzz","BiocManager","dplyr","IRanges","GenomicRanges","ReactomePA","ggplot2","clusterProfiler","topGO","Rgraphviz","pathview","org.Mm.eg.db")){library(i, character.only = T,quietly = T)}; rm(i)
}
#02、给索引
{
  path = "/ChIP_seq_2/aaSHY/ghyuu/publicData/Nsmb/"
  inx1 = "/Reference/aaSHY/Ensembl99/GRCm38/INDEX/STARv279a"
  inx2 = "/Reference/aaSHY/zOther/Rn45s/INDEX/bowtie2/Rn45s"
  inx3 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/Mus_musculus.GRCm38.99.gtf"
  inx4 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf"
  inx5 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf.locInd"
  inx6 = "/Reference/aaSHY/Ensembl99/GRCm38/fa/GRCm38.dna.primary_assembly.fa"
}
################################################################################
#03、跑命令
{
  for (i in c("src","fq","Log","stringtie/count","Nsmb/bam","INDEX/STAR",
              "Nsmb/count","bamCoverage","trim","DESeq2","dumpROOM","bam","count","PLOT")) {
    if (dir.exists(paths = paste0(path, i))==F) {
      dir.create(path = paste0(path,i), recursive = T)
    }
  }
  shell <- paste0(path,"src/run_a.sh")
  cat("#!/bin/bash\n", file = shell)
  cmd_00 <- paste0("#STAR --runThreadN 10 --runMode genomeGenerate ",
                   "--genomeDir ",path,"INDEX/STAR ",
                   "--genomeFastaFiles ",inx6," --sjdbGTFfile ",inx3)
  inx1 <- "/ChIP_seq_2/aaSHY/ghyuu/publicData/Nsmb/INDEX/STAR"
  prex <- c("KO-1","KO-2","WT-1","WT-2")
  cmd_01 <- paste0("#parallel-fastq-dump -t 10 --split-3 --gzip -s ",
                   list.files(path = paste0(path, "rawdata"), pattern = "sra", recursive = T, full.names = T),
                   " -O ", path, "fq", "\n")
  cmd_02 <- paste0("#rename s/fastq/fq/ ", path, "fq/*gz", "\n")
  lays <- paste0(path,"fq/",prex,".fq.gz")
  cmd_03 <- paste0("#trim_galore --phred33 -j 4 -o ", path, "trim ",lays,"\n")
  lays <- paste0(path,"trim/",prex,"_trimmed.fq.gz")
  cmd_04 <- paste0("#STAR --runThreadN 10 --readFilesCommand zcat ",
                   "--outSAMtype BAM SortedByCoordinate --outFileNamePrefix ",path,
                   "Nsmb/bam/", prex, " --outFilterMultimapNmax 1 ",
                   "--sjdbOverhang 100 --genomeDir ",inx1," ",
                   "--readFilesIn ",lays,"\n")
  cmd_05 <- paste0("#mv ",path,"Nsmb/bam/",prex,"Aligned.sortedByCoord.out.bam ",path,
                   "Nsmb/bam/",prex,".bam","\n")
  cmd_06 <- paste0("#samtools index ",path,"Nsmb/bam/",prex,".bam","\n")
  cmd_07 <- paste0("#htseq-count -i gene_name -n 4 -s yes -r pos -a 10 ",path,
                   "Nsmb/bam/",prex,".bam ",inx3,"  >",path,
                   "Nsmb/count/",prex,".txt","\n")
  cmd_08 <- paste0("#TEtranscripts -t ",
                   paste(paste0(path,"Nsmb/bam/",prex[1:2],".bam"),collapse = " ")," -c ",
                   paste(paste0(path,"Nsmb/bam/",prex[3:4],".bam"),collapse = " "),
                   " --GTF ",inx3," --TE ",inx4," --sortByPos --mode multi ",
                   "--project ghyvvko --outdir ",path,"Nsmb/count","\n")
  cmd_09 <- paste0("#TElocal -b ",path,
                   "Nsmb/bam/",prex,".bam --mode multi ",
                   "--GTF ",inx3," --TE ",inx5," --sortByPos --project ",path,
                   "Nsmb/count/",prex,"\n")
  #
  #
  #
  cmd_10 <- paste0("#STAR --runThreadN 10 --readFilesCommand zcat ",
                   "--outMultimapperOrder Random --outSAMtype BAM SortedByCoordinate ",
                   "--outFileNamePrefix ",path,"bam/",prex," ",
                   "--outSAMprimaryFlag AllBestScore --twopassMode Basic ",
                   "--genomeDir ",inx1," --readFilesIn ",lays,"\n")
  cmd_11 <- paste0("#mv ",path,"bam/",prex,"Aligned.sortedByCoord.out.bam ",path,
                   "bam/",prex,".bam","\n")
  cmd_12 <- paste0("#samtools index ",path,"bam/",prex,".bam","\n")
  cmd_13 <- paste0("#bamCoverage --normalizeUsing CPM -p 10 ",
                   "--minMappingQuality 1 --samFlagExclude 256 -of ",
                   "bigwig -bs 20 -b ",path,"bam/",prex,".bam ",
                   "-o ",path,"bamCoverage/",prex,".bw","\n")
  cmd_14 <- paste0("#TEtranscripts -t ",
                   paste(paste0(path,"bam/",prex[1:2],".bam"),collapse = " ")," -c ",
                   paste(paste0(path,"bam/",prex[3:4],".bam"),collapse = " "),
                   " --GTF ",inx3," --TE ",inx4," --sortByPos --mode multi ",
                   "--project ghyvvko --outdir ",path,"count","\n")
  cmd_15 <- paste0("#TElocal -b ",path,"bam/",prex,".bam ",
                   "--GTF ",inx3," --TE ",inx5," --sortByPos --project ",path, 
                   "count/",prex,"\n")
  cmd_16 <- paste0("#stringtie ", path,"bam/",prex,".bam ",
                   "-p 10 -G ", inx3, " -o ",path,"stringtie/",prex,".gtf","\n")
  cmd_17 <- paste0("#stringtie --merge -o ",
                   path,"stringtie/ghyvvko.merge.gtf -G ",inx3," ",
                   paste(paste0(path,"stringtie/",prex,".gtf"), collapse = " "),"\n")
  cmd_18 <- paste0("#cat ", path,"stringtie/ghyvvko.merge.gtf ", 
                   "|sed 's/gene_name/gene_abcd/g; s/\\tgene_id/\\tgene_name/g' ",
                   ">",path,"stringtie/ghyvvko.merge.r.gtf","\n")
  cmd_19 <- paste0("#TEtranscripts ",
                   "-t ",paste(paste0(path,"bam/",prex[1:2],".bam"),collapse = " ")," ",
                   "-c ",paste(paste0(path,"bam/",prex[3:4],".bam"),collapse = " ")," ",
                   "--GTF ",path,"stringtie/ghyvvko.merge.r.gtf ",
                   "--TE ", inx4," --sortByPos --mode multi --project ghyvvko ",
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
  for (i in cmd_14) {cat(i, file = shell, append = T)}
  for (i in cmd_15) {cat(i, file = shell, append = T)}
  for (i in cmd_16) {cat(i, file = shell, append = T)}
  for (i in cmd_17) {cat(i, file = shell, append = T)}
  for (i in cmd_18) {cat(i, file = shell, append = T)}
  for (i in cmd_19) {cat(i, file = shell, append = T)}
  print(paste0("nohup bash ",shell, " >",paste0(path,"Log/"),basename(shell),".log"," 2>&1 &"))
}
################################################################################
#05、CPM Excel
{
  path = "/ChIP_seq_2/aaSHY/ghyuu/publicData/Nsmb/"
  cnt1 <- read.table(paste0(path,"Nsmb/htseq/ghyvvko.cntTable"), header = T, 
                     sep = "\t", row.names = 1)
  colnames(cnt1) <- c("KO-1","KO-2","WT-1","WT-2")
  tmp2 <- cnt1
  tmp2 <- as.data.frame(apply(tmp2, 2, function(x){x/sum(x)*1000000}))
  tmp2$theName <- rownames(tmp2)
  openxlsx::write.xlsx(x = tmp2, file = paste0(path,"Nsmb/htseq/ghyvvko.cpm.onlyGene.xlsx"))
}
################################################################################
#04、TE差异分析
{
  shy_diff <- function(prex,mark,htsqFile,geneFile, siteFile, Class, VS){
    tmp1 <- read.table(geneFile, header = T, check.names = F, stringsAsFactors = F, row.names = 1)
    colnames(tmp1) <- sapply(str_split(basename(colnames(tmp1)),pattern = ".bam"), "[", 1)
    dfclas <- data.frame(row.names = colnames(tmp1), Class, stringsAsFactors = T)
    geneAA <- tmp1[grep(rownames(tmp1), pattern = ":", invert = T),]
    teOnly <- tmp1[grep(rownames(tmp1), pattern = ":", invert = F),]
    teGene <- tmp1
    #
    geneAA <- DESeq2::DESeqDataSetFromMatrix(countData = geneAA, colData = dfclas, design = ~Class)
    #geneAA <- geneAA[rowMeans(BiocGenerics::counts(geneAA)) > 5,]
    geneAA <- geneAA[rowSums(BiocGenerics::counts(geneAA)) >= 10,]
    geneAA <- DESeq2::DESeq(geneAA, fitType = "parametric")
    geneAA <- DESeq2::results(geneAA, contrast = VS, alpha = 0.05, pAdjustMethod = "BH")
    write.csv(geneAA, paste0(path,"DESeq2/",mark, ".geneOnly.csv"))
    upup <- subset(geneAA, log2FoldChange >  log2(2) & padj < 0.05)
    down <- subset(geneAA, log2FoldChange < -log2(2) & padj < 0.05)
    write.csv(upup,paste0(path,"DESeq2/",mark, ".geneOnly.upup.csv"))
    write.csv(down,paste0(path,"DESeq2/",mark, ".geneOnly.down.csv"))
    #
    teOnly <- DESeq2::DESeqDataSetFromMatrix(countData = teOnly, colData = dfclas, design = ~Class)
    #teOnly <- teOnly[rowMeans(BiocGenerics::counts(teOnly)) > 5, ]
    teOnly <- teOnly[rowSums(BiocGenerics::counts(teOnly)) >= 10,]
    teOnly <- DESeq2::DESeq(teOnly, fitType = "parametric")
    teOnly <- DESeq2::results(teOnly, contrast = VS, alpha = 0.05, pAdjustMethod = "BH")
    write.csv(teOnly, paste0(path,"DESeq2/",mark, ".teOnly.csv"))
    upup <- subset(teOnly, log2FoldChange >  log2(2) & padj < 0.05)
    down <- subset(teOnly, log2FoldChange < -log2(2) & padj < 0.05)
    write.csv(upup,paste0(path,"DESeq2/",mark, ".teOnly.upup.csv"))
    write.csv(down,paste0(path,"DESeq2/",mark, ".teOnly.down.csv"))
    #
    teGene <- DESeq2::DESeqDataSetFromMatrix(countData = teGene, colData = dfclas, design = ~Class)
    #teGene <- teGene[rowMeans(BiocGenerics::counts(teGene)) > 5, ]
    teGene <- teGene[rowSums(BiocGenerics::counts(teGene)) >= 10,]
    teGene <- DESeq2::DESeq(teGene, fitType = "parametric")
    teGene <- DESeq2::results(teGene, contrast = VS, alpha = 0.05, pAdjustMethod = "BH")
    teGene <- teGene[grep(rownames(teGene), pattern = ":", invert = F),]
    write.csv(teGene, paste0(path,"DESeq2/",mark, ".te.GeneBackGround.csv"))
    upup <- subset(teGene, log2FoldChange >  log2(2) & padj < 0.05)
    down <- subset(teGene, log2FoldChange < -log2(2) & padj < 0.05)
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
    #teGene <- teGene[rowMeans(BiocGenerics::counts(teGene)) > 5, ]
    teGene <- teGene[rowSums(BiocGenerics::counts(teGene)) >= 10,]
    teGene <- DESeq2::DESeq(teGene, fitType = "parametric")
    teGene <- DESeq2::results(teGene, contrast = VS, alpha = 0.05, pAdjustMethod = "BH")
    teGene <- teGene[grep(rownames(teGene), pattern = ":", invert = F),]
    write.csv(teGene, paste0(path,"DESeq2/",mark, ".te.site.GeneBackGround.csv"))
    upup <- subset(teGene, log2FoldChange >  log2(2) & padj < 0.05)
    down <- subset(teGene, log2FoldChange < -log2(2) & padj < 0.05)
    write.csv(upup, paste0(path,"DESeq2/",mark, ".te.site.GeneBackGround.upup.csv"))
    write.csv(down, paste0(path,"DESeq2/",mark, ".te.site.GeneBackGround.down.csv"))
    #
    teGene <- tmp2[grep(rownames(tmp2), pattern = ":", invert = F),]
    teGene <- DESeq2::DESeqDataSetFromMatrix(countData = teGene, colData = dfclas, design = ~Class)
    #teGene <- teGene[rowMeans(BiocGenerics::counts(teGene)) > 5, ]
    teGene <- teGene[rowSums(BiocGenerics::counts(teGene)) >= 10,]
    teGene <- DESeq2::DESeq(teGene, fitType = "parametric")
    teGene <- DESeq2::results(teGene, contrast = VS, alpha = 0.05, pAdjustMethod = "BH")
    write.csv(teGene, paste0(path,"DESeq2/",mark, ".te.site.teOnly.csv"))
    upup <- subset(teGene, log2FoldChange >  log2(2) & padj < 0.05)
    down <- subset(teGene, log2FoldChange < -log2(2)& padj < 0.05)
    write.csv(upup, paste0(path,"DESeq2/",mark, ".te.site.teOnly.upup.csv"))
    write.csv(down, paste0(path,"DESeq2/",mark, ".te.site.teOnly.down.csv"))
    #
    #
    #
    #
    #补充HTseq的
    tmp3 <- read.table(htsqFile, header = F, 
                       check.names = F, stringsAsFactors = F, row.names = 1)
    colnames(tmp3) <- prex
    dfclas <- data.frame(row.names = colnames(tmp3), Class, stringsAsFactors = T)
    geneAA <- tmp3
    #
    geneAA <- DESeq2::DESeqDataSetFromMatrix(countData = geneAA, 
                                             colData = dfclas, design = ~Class)
    #geneAA <- geneAA[rowMeans(BiocGenerics::counts(geneAA)) > 5,]
    geneAA <- geneAA[rowSums(BiocGenerics::counts(geneAA)) >= 10,]
    
    geneAA <- DESeq2::DESeq(geneAA, fitType = "parametric")
    geneAA <- DESeq2::results(geneAA, contrast = VS, alpha = 0.05, pAdjustMethod = "BH")
    write.csv(geneAA, paste0(path,"DESeq2/htseq-",mark, ".geneOnly.csv"))
    upup <- subset(geneAA, log2FoldChange >  log2(2) & padj < 0.05)
    down <- subset(geneAA, log2FoldChange < -log2(2) & padj < 0.05)
    write.csv(upup,paste0(path,"DESeq2/htseq-",mark, ".geneOnly.upup.csv"))
    write.csv(down,paste0(path,"DESeq2/htseq-",mark, ".geneOnly.down.csv"))
  }
  path = "/ChIP_seq_2/aaSHY/ghyuu/publicData/Nsmb/Nsmb/"
  shy_diff(mark = "ghyvvko", 
           prex = c("KO-1","KO-2","WT-1","WT-2"),
           geneFile = paste0(path, "count/ghyvvko.cntTable"),
           htsqFile = paste0(path, "htseq/ghyvvko.cntTable"),
           siteFile = paste0(path, "count/ghyvvko.site.cntTable"),
           Class = factor(c("KO","KO","WT","WT")), VS = c("Class","KO","WT"))
}
################################################################################
#05、Volcano for TEs
{
  shy_volcano_te <- function(resFile, saveFile){
    resdata <- as.data.frame(read.csv(resFile, header = T, stringsAsFactors = F))
    resdata <- na.omit(resdata)
    resdata$threshold = "Other TEs"
    resdata$threshold[resdata$log2FoldChange > log2(2)  & resdata$padj <0.05] = "Upregulated TEs"
    resdata$threshold[resdata$log2FoldChange < -log2(2) & resdata$padj <0.05] = "Downregulated TEs"
    resdata$threshold <- factor(resdata$threshold, levels = c("Upregulated TEs","Downregulated TEs","Other TEs"))
    resdata[which(resdata$threshold=="Other TEs"),"X"] <- NA
    labb <- resdata#[which(resdata$threshold!="Other TEs"),]
    if (unname(table(resdata$threshold)[1]) + unname(table(resdata$threshold)[2]) >20){
      labb <- labb[which(labb$padj<10^(-5)),]#or -10
    } else {labb <- labb[!(is.na(labb)),]}
    p_2 <- ggplot(resdata, aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) +
      labs(col="", x="Log2FoldChange",y="-Log10(adjusted p-value)") + 
      geom_point(size=0.8) +
      scale_x_continuous(limits = c(-4,4)) +
      geom_hline(aes(yintercept=-log10(0.05)), linetype="dashed", colour="royalblue") +
      geom_vline(aes(xintercept=-log2(2)), linetype="dashed", colour="royalblue") +
      geom_vline(aes(xintercept=log2(2)), linetype="dashed", colour="royalblue") +#ylim(0,max(-log10(resdata$padj)[is.finite(-log10(resdata$padj))])+20) +
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
    #View(resdata)
    ggsave(plot = p_2, saveFile, units = "cm", width = 18, height = 16)
  }
  shy_volcano_te(resFile  = paste0(path, "DESeq2/ghyvvko.te.GeneBackGround.csv"),
                 saveFile = paste0(path, "PLOT/ghyvvko.te.GeneBackGround.volcano.pdf"))
  shy_volcano_te(resFile  = paste0(path, "DESeq2/ghyvvko.teOnly.csv"),
                 saveFile = paste0(path, "PLOT/ghyvvko.teOnly.volcano.pdf"))
}
################################################################################
#06、gene差异分析
{
  shy_volcano_gene <- function(resFile, saveFile){
    resdata <- as.data.frame(read.csv(resFile, header = T, stringsAsFactors = F))
    resdata <- na.omit(resdata)
    resdata$threshold = "Other genes"
    resdata$threshold[resdata$log2FoldChange >  log2(2.0) & resdata$padj <0.05] = "Upregulated genes"
    resdata$threshold[resdata$log2FoldChange < -log2(2.0) & resdata$padj <0.05] = "Downregulated genes"
    resdata$threshold <- factor(resdata$threshold, levels = c("Upregulated genes","Downregulated genes","Other genes"))
    p_1 <- ggplot(resdata, aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) +
      labs(col="", x="Log2FoldChange",y="-Log10(adjusted p-value)") + 
      geom_point(size=0.8) +
      scale_x_continuous(limits = c(-6, 9), breaks = seq(-6, 9, 2)) +
      geom_hline(aes(yintercept=-log10(0.05)), linetype="dashed", colour="royalblue") +
      geom_vline(aes(xintercept=-log2(2)), linetype="dashed", colour="royalblue") +
      geom_vline(aes(xintercept= log2(2)), linetype="dashed", colour="royalblue") +
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
  shy_volcano_gene(resFile  = paste0(path, "DESeq2/htseq-ghyvvko.geneOnly.csv"),
                   saveFile = paste0(path, "PLOT/htseq-ghyvvko.geneOnly.volcano.pdf"))
}
################################################################################
#07、富集分析 [KEGG]
{
  resFile <- paste0(path, "DESeq2/htseq-ghyvvko.geneOnly.csv")
  resdata <- as.data.frame(read.csv(resFile, header = T, stringsAsFactors = F))
  resdata <- resdata[grep(resdata[,1], pattern = ":", invert = T),]
  mtacup <- resdata[which(resdata$log2FoldChange >  log2(2.0) & resdata$padj < 0.05),"X"]
  kddown <- resdata[which(resdata$log2FoldChange < -log2(2.0) & resdata$padj < 0.05),"X"]
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
    geom_vline(aes(xintercept = 5), color = "white", linewidth  = 0.8) +
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
         filename = paste0(path,"PLOT/ghyvvko.kegg.pdf"))
}
################################################################################
#08、富集分析 [GO]
{
  resFile <- paste0(path, "DESeq2/htseq-ghyvvko.geneOnly.csv")
  resdata <- as.data.frame(read.csv(resFile, header = T, stringsAsFactors = F))
  resdata <- resdata[grep(resdata[,1], pattern = ":", invert = T),]
  mtacup <- resdata[which(resdata$log2FoldChange >  log2(2.0) & resdata$padj < 0.05),"X"]
  kddown <- resdata[which(resdata$log2FoldChange < -log2(2.0) & resdata$padj < 0.05),"X"]
  for (symb in seq(2)) {
    ovgene <- unlist(list(mtacup,kddown)[symb])
    geyy <- na.omit(mapIds(x=org.Mm.eg.db, keys=unique(ovgene), keytype="SYMBOL", column="ENTREZID"))
    gogo <- clusterProfiler::enrichGO(gene = geyy, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", 
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
         filename = paste0(path,"PLOT/ghyvvko.gomf.pdf"))
}
################################################################################
#09、指定term的GSEA
{
  resFile <- paste0(path, "DESeq2/htseq-ghyvvko.geneOnly.csv")
  aa <- read.csv(file = resFile, header = T, stringsAsFactors = F)[,c(1,3)]
  ab <- suppressWarnings(clusterProfiler::bitr(geneID = aa$X, fromType = "SYMBOL", 
                                               toType = "ENTREZID", 
                                               OrgDb = "org.Mm.eg.db", drop = T))
  ac <- dplyr::distinct(.data = ab, SYMBOL, .keep_all = T)
  ac$rank <- apply(ac, 1, function(x){aa[which(aa$X==x[1]),2]})
  ad <- ac[order(ac$rank, decreasing = T),c(2,3)]
  ae <- ad$rank; names(ae) <- ad$ENTREZID
  rm(aa,ab,ac,ad)
  af_gogo <- suppressWarnings(gseGO(geneList = ae, ont = "BP", seed = T,
                                    maxGSSize = 1000, verbose = T, 
                                    keyType = "ENTREZID", pvalueCutoff = 1, 
                                    OrgDb = "org.Mm.eg.db", by = "fgsea"))
  testxxx <- af_gogo@result
  testxxx[,10] <- apply(af_gogo@result, 1, function(x){chartr(old = ", ", new = "; ", x[10])})
  write.csv(testxxx, row.names = F, quote = F,
            file = paste0(path,"PLOT/ghyvvko.gseGO.csv"))
  #
  #
  aims <- af_gogo[grep("apoptotic",af_gogo@result$Description),]
  aims <- aims[order(aims$pvalue, decreasing = F),1]
  aims <- match(aims,af_gogo@result$ID)
  aim2 <- split(1:70, ceiling(seq_along(1:70)/7))
  for (i in seq_along(aim2)) {
    print(i)
    p_1 <- enrichplot::gseaplot2(x = af_gogo, 
                                 aims[aim2[[i]]], 
                                 title = paste0("gseGO: Apoptotic subgroup-",i), 
                                 rel_heights = c(1.5, 0.3, 0.6), 
                                 color = ggsci::pal_lancet()(7), pvalue_table = T)
    ggsave(plot = p_1,
           units = "cm", width = 35, height = 15,
           filename = paste0(path,"PLOT/ghyvvko.gseGO-",i,".pdf"))
  }
}
################################################################################
#10、凋亡相关基因的表达变化
{
  #GO_DATA <- get_GO_data("org.Mm.eg.db", "ALL", "SYMBOL"); base::save(list = c("GO_DATA"), file = "/Reference/aaSHY/DatabaseFiles/GO_DATA.RData")
  base::load(file = "/Reference/aaSHY/DatabaseFiles/GO_DATA.RData")
  #GO ID和GO ONT
  tmp1 <- data.frame(goID = base::names(GO_DATA$GO2ONT),
                     goONT = base::unname(GO_DATA$GO2ONT))
  #GO ID和GO Name
  tmp2 <- data.frame(goID = base::names(GO_DATA$PATHID2NAME),goName = base::unname(GO_DATA$PATHID2NAME))
  tmp2 <- tmp2[which(tmp2$goID %in% tmp1$goID[which(tmp1$goONT=="BP")]),]
  #GO ID和GO Genes Name
  tmp3 <- GO_DATA$PATHID2EXTID
  ##Gene和GO IDs
  tmp4 <- GO_DATA$EXTID2PATHID
  ##
  myid <- tmp2[grep(x = tmp2$goName, pattern = "apoptotic", ignore.case = T),1]
  set1 <- unique(unname(unlist(tmp3[myid])))
  ##
  red1 <- openxlsx::read.xlsx(xlsxFile = paste0(path,"Nsmb/htseq/ghyvvko.cpm.onlyGene.xlsx"))
  red1 <- red1[which(red1$theName %in% set1),]
  red2 <- red1[,1:4] %>%
    tidyr::pivot_longer(cols = 1:4, names_to = "type", values_to = "cpms")
  red3 <- red2
  red3$clas <- sapply(str_split(red3$type, "-"),"[",1)
  red3$type <- factor(red3$type, levels = prex[c(3,4,1,2)])
  ggplot(data = red3) +
    geom_boxplot(mapping = aes(x = type, y = log2(cpms + 0.1), fill = clas)) +
    theme_classic() +
    labs(x = "", 
         y = "Log2(CPM +0.1)", 
         title = "all genes appeared in apoptotic terms (retrived from GO database)",
         fill = "") +
    theme(axis.title = element_text(family = "serif", size=12),
          axis.text  = element_text(family = "serif", size=12),
          legend.text = element_text(family = "serif", size=12),
          plot.margin = margin(unit = "cm",1,1,1,1))
}
################################################################################
#11、凋亡相关基因的表达变化
{
  set1 <- c("Trp53","Bax","Bak1","Bcl2",
            "Bbc3", "Bad", "Bcl2l1", "Mcl1", "Birc2", "Birc3",
            "Casp3","Casp9","Fas","Apaf1","Bid","Xiap")
  ##
  red1 <- openxlsx::read.xlsx(xlsxFile = paste0(path,"Nsmb/htseq/ghyvvko.cpm.onlyGene.xlsx"))
  #red1 <- openxlsx::read.xlsx(xlsxFile = "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuu/count/ghyuu.cpm.onlyGene.xlsx")
  #red1 <- openxlsx::read.xlsx(xlsxFile = "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/theghyuughyvv/count/ghyuughyvv.cpm.onlyGene.xlsx")[,3:7]
  red1 <- red1[which(red1$theName %in% set1),]
  red2 <- red1 %>%
    tidyr::pivot_longer(cols = 1:4, names_to = "type", values_to = "cpms")
  red3 <- red2
  red3$clas <- sapply(str_split(red3$type, "-"),"[",1)
  red3$type <- factor(red3$type, levels = prex[c(3,4,1,2)])
  #red3$type <- factor(red3$type, levels = colnames(red1)[c(3,4,1,2)])
  ggplot(data = red3) +
    geom_bar(aes(x = type, y = cpms, fill = clas), stat = "identity") +
    facet_wrap(facets = ~theName, scales = "free_y", ncol = 3)
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
################################################################################
{
  ###path = "/ChIP_seq_2/aaSHY/ghyuu/publicData/Nsmb/Nsmb/"
  ###resFile <- paste0(path, "DESeq2/htseq-ghyvvko.geneOnly.csv")
  ###aa <- read.csv(file = resFile, header = T, stringsAsFactors = F)[,c(1,3)]
  ###ab <- suppressWarnings(clusterProfiler::bitr(geneID = aa$X, fromType = "SYMBOL", 
  ###                                             toType = "ENTREZID", 
  ###                                             OrgDb = "org.Mm.eg.db", drop = T))
  ###ac <- dplyr::distinct(.data = ab, SYMBOL, .keep_all = T)
  ###ac$rank <- apply(ac, 1, function(x){aa[which(aa$X==x[1]),2]})
  ###ad <- ac[order(ac$rank, decreasing = T),c(2,3)]
  ###ae <- ad$rank; names(ae) <- ad$ENTREZID
  ###rm(aa,ab,ac,ad)
  ###af_gogo <- suppressWarnings(gseGO(geneList = ae, ont = "BP", seed = T,
  ###                                  maxGSSize = 1000, verbose = T, 
  ###                                  keyType = "ENTREZID", pvalueCutoff = 1, 
  ###                                  OrgDb = "org.Mm.eg.db", by = "fgsea"))
  ###testxxx <- af_gogo@result
  ###testxxx[,10] <- apply(af_gogo@result, 1, function(x){chartr(old = ", ", new = "; ", x[10])})
  ###write.csv(testxxx, row.names = F, quote = F,
  ###          file = paste0(path,"PLOT/ghyvvko.gseGO.csv"))
  #
  #
  gseo <- read.csv(paste0(path,"PLOT/ghyvvko.gseGO.csv"))
  aims <- gseo[grep("apoptotic|apoptosis", gseo$Description),]
  aims <- aims[order(aims$pvalue, decreasing = F), 1:9]
  to_plot <- aims[1:7,]
  colnames(to_plot)
  to_plot$Description <- factor(to_plot$Description, levels = rev(to_plot$Description))
  P001 <- ggplot(data = to_plot) +
    geom_bar(aes(x = Description, y = NES,  fill = -log2(pvalue)),
             stat = "identity", width = 0.8) +
    scale_fill_gradient2() +
    labs(x = "gseGO terms", 
         y = "NES (Normalized Enrichment Score)",fill ="-log2 (P-value)",
         title = "gseGO that related to apoptotic") +
    geom_text(aes(label = Description, 
                  y = 0.05, x = Description), 
              hjust = 0, vjust = 0.5, size = 5, color = "white", family ="serif") +
    theme_classic() +
    theme(plot.margin = margin(unit = "cm", 1,1,1,1),
          plot.title  = element_text(family= "serif", size = 18, face = "bold"),
          legend.text = element_text(family= "serif", size = 14),
          legend.title= element_text(family= "serif", size = 14),
          axis.text   = element_text(family= "serif", size = 14),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title  = element_text(family= "serif", size = 16)) +
    coord_flip()
  P001
  ggsave(plot = P001, 
         filename = paste0(path,"PLOT_need/2025Nsmb-gseGO-apoptotic.pdf"),
         units = "cm", width = 24, height = 16)
  #
  #
  #
  #
  #
  #
  aims <- match(aims, gseo$ID)
  p_1 <- enrichplot::gseaplot2(x = gseo, 
                               aims[aim2[[i]]], 
                               title = paste0("gseGO: Apoptotic subgroup-",i), 
                               rel_heights = c(1.5, 0.3, 0.6), 
                               color = ggsci::pal_lancet()(7), pvalue_table = T)
  
  aim2 <- split(1:70, ceiling(seq_along(1:70)/7))
  for (i in seq_along(aim2)) {
    print(i)
    p_1 <- enrichplot::gseaplot2(x = gseo, 
                                 aims[aim2[[i]]], 
                                 title = paste0("gseGO: Apoptotic subgroup-",i), 
                                 rel_heights = c(1.5, 0.3, 0.6), 
                                 color = ggsci::pal_lancet()(7), pvalue_table = T)
    ggsave(plot = p_1,
           units = "cm", width = 35, height = 15,
           filename = paste0(path,"PLOT/ghyvvko.gseGO-",i,".pdf"))
  }
}





