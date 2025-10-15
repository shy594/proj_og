###OE
#Write by Sun Haiayng at 2025.06.03
#01、搭环境
{
  rm(list=ls())
  for (i in c("data.table","ggridges","ggrepel","KEGG.db","scales","stringr",
              "Biostrings","tidyr","pheatmap","enrichplot",
              "DESeq2","patchwork","ggpubr","ggsci","ggtranscript",
              "rtracklayer","Mfuzz","BiocManager","dplyr",
              "IRanges","GenomicRanges","ReactomePA","ggplot2",
              "clusterProfiler","topGO","Rgraphviz","pathview","org.Mm.eg.db")){library(i, character.only = T,quietly = T)}; rm(i)
}
#02、给索引
{
  path = "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/thejjuuzzRY2/Chi/"
  #path = "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/thejjuuzzRY2/ChiR/"
  inx1 = "/Reference/aaSHY/Ensembl99/GRCm38/INDEX/STARv279a"
  inx2 <- "/Reference/aaSHY/Ensembl99/GRCm38/fa/GRCm38.dna.primary_assembly.fa"
  inx3 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/Mus_musculus.GRCm38.99.gtf"
  inx4 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf"
  inx5 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf.locInd"
  #
  #
  for (i in c("testScale/Quantile","testScale/HK","testScale/DESeq2Self")) {
    if (dir.exists(paths = paste0(path, i))==F) {
      dir.create(path = paste0(path,i), recursive = T)
    }
  }
}
################################################################################
#03、尝试标准化
{
  #####预置变量
  {
    ff01 <- paste0(path,"count/jjuuzz.cntTable")
    ff02 <- paste0(path,"countUniq/jjuuzz.siteUniq.cntTable")
    mark = "jjuuzz"
    clas = factor(c("OE","OE","WT","WT")); vsvs = c("clas","OE","WT")
  }
  #####方案一, Quantile
  {
    shy_cpm <- function (read, nm_1) {
      cnt1 <- read
      tmp1 <- as.data.frame(apply(cnt1, 2, function(x){x/sum(x)*1000000}))
      tmp1$theName <- rownames(tmp1)
      openxlsx::write.xlsx(x = tmp1, 
                           paste0(path, "testScale/Quantile/",mark,".cpm.",nm_1,".xlsx"))
    }
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
      norm_data <- as.data.frame(norm_data)
      rownames(norm_data) <- rownames(myta)
      colnames(norm_data) <- colnames(myta)
      return(norm_data)
    }
    tempss <- read.table(ff01, header = T, check.names = F, stringsAsFactors = F, row.names = 1)
    colnames(tempss) <- gsub(".bam..","",basename(colnames(tempss)))
    teSite <- read.table(ff02, header = T, check.names = F, stringsAsFactors = F, row.names = 1)
    colnames(teSite) <- gsub(".bam..","",basename(colnames(teSite)))
    dfclas <- data.frame(row.names = colnames(tempss), clas, stringsAsFactors = T)
    #
    geOnly <- ceiling(qnqn(myta = tempss[grep(rownames(tempss), pattern = ":", invert = T),]))
    teOnly <- ceiling(qnqn(myta = tempss[grep(rownames(tempss), pattern = ":", invert = F),]))
    teGene <- ceiling(qnqn(myta = tempss))
    teSite <- ceiling(qnqn(myta = teSite[grep(rownames(teSite), pattern = ":", invert = F),]))
    shy_cpm(read = geOnly, nm_1 = "onlyGene")
    shy_cpm(read = teOnly, nm_1 = "onlyTE")
    shy_cpm(read = teGene, nm_1 = "allGeneTE")
    shy_cpm(read = teSite, nm_1 = "onlyTE.site")
    #
    teGene <- DESeq2::DESeqDataSetFromMatrix(countData = teGene, colData = dfclas, design = ~clas)
    teGene <- teGene[rowMeans(BiocGenerics::counts(teGene)) > 5,]
    teGene <- DESeq2::DESeq(teGene)
    teGene <- DESeq2::results(teGene, contrast = vsvs)
    teGene <- teGene[grep(rownames(teGene), pattern = ":", invert = F),]
    write.csv(teGene, paste0(path,"testScale/Quantile/",mark,".te.GeneBackGround.csv"))
    #
    geOnly <- DESeq2::DESeqDataSetFromMatrix(countData = geOnly, colData = dfclas, design = ~clas)
    geOnly <- geOnly[rowMeans(BiocGenerics::counts(geOnly)) > 5,]
    geOnly <- DESeq2::DESeq(geOnly,)
    geOnly <- DESeq2::results(geOnly, contrast = vsvs)
    write.csv(geOnly, paste0(path,"testScale/Quantile/",mark,".geneOnly.csv"))
    #
    teOnly <- DESeq2::DESeqDataSetFromMatrix(countData = teOnly, colData = dfclas, design = ~clas)
    teOnly <- teOnly[rowMeans(BiocGenerics::counts(teOnly)) > 5,]
    teOnly <- DESeq2::DESeq(teOnly,)
    teOnly <- DESeq2::results(teOnly, contrast = vsvs)
    write.csv(teOnly, paste0(path,"testScale/Quantile/",mark,".teOnly.csv"))
    #
    teSite <- DESeq2::DESeqDataSetFromMatrix(countData = teSite, colData = dfclas, design = ~clas)
    teSite <- teSite[rowMeans(BiocGenerics::counts(teSite)) > 5,]
    teSite <- DESeq2::DESeq(teSite,)
    teSite <- DESeq2::results(teSite, contrast = vsvs)
    write.csv(teOnly, paste0(path,"testScale/Quantile/",mark,".te.site.teOnly.csv"))
  }
  #####方案二, Housekeeping gene, 8407个mouse的 [太多了, 先不跑________________]
  {
    tempss <- read.table(ff01, header = T, check.names = F, stringsAsFactors = F, row.names = 1)
    hkGene <- openxlsx::read.xlsx("/Reference/aaSHY/DatabaseFiles/houseKeep/set1/hk_mm10_15tissues.xlsx",sheet = "Sheet1")
    #
    normed <- tempss
    dfclas <- data.frame(row.names = colnames(tempss), clas, stringsAsFactors = T)
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
    
    
    
    
    teGene <- teGene[rowMeans(BiocGenerics::counts(teGene)) > 5,]
    teGene <- DESeq2::DESeq(teGene)
    teGene <- DESeq2::results(teGene, contrast = vsvs)
    teGene <- teGene[grep(rownames(teGene), pattern = ":", invert = F),]
    write.csv(teGene, paste0(path,"testScale/",mark,".te.GeneBackGround.csv"))
    
  }
  #####方案三, Housekeeping gene, 3804个human的, 转换获得3460个小鼠同源基因
  {
    #管家基因与reads count文件读取
    hkGene <- read.csv("/Reference/aaSHY/DatabaseFiles/houseKeep/set2/hk_human_3804.txt", header = F)
    library("babelgene")
    hkGene <- babelgene::orthologs(genes = unique(hkGene$V1), species = "mouse", human = T)
    hkGene <- unique(hkGene$symbol)
    tempss <- read.table(ff01, header = T,
                         check.names  = F, stringsAsFactors = F, row.names = 1)
    colnames(tempss) <- gsub(".bam..","",basename(colnames(tempss)))
    dfclas <- data.frame(row.names = colnames(tempss), clas, stringsAsFactors = T)
    #
    #
    #获取HK size factors
    raws_1 <- tempss
    raws_1 <- DESeq2::DESeqDataSetFromMatrix(countData = raws_1,
                                             colData = dfclas, design = ~clas)
    raws_1 <- raws_1[rowMeans(BiocGenerics::counts(raws_1)) > 5,]
    hkGene <- raws_1[which(row.names(raws_1) %in% hkGene),]
    hkGene <- estimateSizeFactors(hkGene)
    hkGene <- hkGene$sizeFactor
    #
    #
    #
    #
    #基于Size factor计算CPM
    shx_cpm <- function (read, nm_1) {
      cnt1 <- read
      cnt2 <- data.frame(a = cnt1[,1]/hkGene[1],b = cnt1[,2]/hkGene[2],
                         c = cnt1[,3]/hkGene[3],d = cnt1[,4]/hkGene[4])
      colnames(cnt2) <- colnames(cnt1)
      tmp1 <- as.data.frame(apply(cnt2, 2, function(x){x/sum(x)*1000000}))
      tmp1$theName <- rownames(tmp1)
      openxlsx::write.xlsx(x = tmp1, 
                           paste0(path, "testScale/HK/",mark,".cpm.",nm_1,".xlsx"))
    }
    geOnly <- tempss[grep(rownames(tempss), pattern = ":", invert = T),]
    teOnly <- tempss[grep(rownames(tempss), pattern = ":", invert = F),]
    teGene <- tempss
    shx_cpm(read = geOnly, nm_1 = "onlyGene")
    shx_cpm(read = teOnly, nm_1 = "onlyTE")
    shx_cpm(read = teGene, nm_1 = "allGeneTE")
    #
    #
    #
    #
    #基于Size factor做差异分析
    shx_dif <- function(read, nm_1) {
      raws_1 <- read
      raws_1 <- DESeq2::DESeqDataSetFromMatrix(countData = raws_1,
                                               colData = dfclas, design = ~clas)
      raws_1 <- raws_1[rowMeans(BiocGenerics::counts(raws_1)) > 5,]
      raws_2 <- estimateSizeFactors(raws_1)
      raws_2$sizeFactor <- hkGene #更改size factor后, 继续跑
      raws_2 <- estimateDispersions(raws_2)
      raws_2 <- nbinomWaldTest(raws_2)
      deseqs <- results(raws_2, contrast = vsvs)
      write.csv(deseqs, paste0(path,"testScale/HK/",mark,".",nm_1,".csv"))
    }
    shx_dif(read = geOnly, nm_1 = "geneOnly")
    shx_dif(read = teOnly, nm_1 = "teOnly")
    shx_dif(read = teGene, nm_1 = "te.GeneBackGround")
  }
  #####方案四, DESeq2加上批次信息（仅限没有矫正突变组）
  #############iso2-2是先收的，iso-1是后收的；CR-1和CR-2是一起收的，OV-1和OV-2也是一起收的
  {
    shy_diff <- function(mark, geneFile, siteFile, Class, Batch, VS){
      tmp1 <- read.table(geneFile, header = T, check.names = F, stringsAsFactors = F, row.names = 1)
      colnames(tmp1) <- sapply(str_split(basename(colnames(tmp1)),pattern = ".bam"), "[", 1)
      dfclas <- data.frame(row.names = colnames(tmp1), Class, stringsAsFactors = T)
      colLen <- ceiling(length(colnames(tmp1))/2)
      geneAA <- tmp1[grep(rownames(tmp1), pattern = ":", invert = T),]
      teOnly <- tmp1[grep(rownames(tmp1), pattern = ":", invert = F),]
      teGene <- tmp1
      #
      geneAA <- DESeq2::DESeqDataSetFromMatrix(countData = geneAA, colData = dfclas, design = ~Class)
      geneAA <- geneAA[rowMeans(BiocGenerics::counts(geneAA)) > 5,]
      geneAA <- DESeq2::DESeq(geneAA,)
      geneAA <- DESeq2::results(geneAA, contrast = VS)
      write.csv(geneAA, paste0(bath,"DESeq2/",mark, ".geneOnly.csv"))
      upup <- subset(geneAA, log2FoldChange > log2(1.5) & padj < 0.05)
      down <- subset(geneAA, log2FoldChange < -log2(1.5)& padj < 0.05)
      write.csv(upup,paste0(bath,"DESeq2/",mark, ".geneOnly.upup.csv"))
      write.csv(down,paste0(bath,"DESeq2/",mark, ".geneOnly.down.csv"))
      #
      teOnly <- DESeq2::DESeqDataSetFromMatrix(countData = teOnly, colData = dfclas, design = ~Class)
      teOnly <- teOnly[rowMeans(BiocGenerics::counts(teOnly)) > 5,]
      teOnly <- DESeq2::DESeq(teOnly,)
      teOnly <- DESeq2::results(teOnly, contrast = VS)
      write.csv(teOnly, paste0(bath,"DESeq2/",mark, ".teOnly.csv"))
      upup <- subset(teOnly, log2FoldChange > log2(1.5) & padj < 0.05)
      down <- subset(teOnly, log2FoldChange < -log2(1.5)& padj < 0.05)
      write.csv(upup,paste0(bath,"DESeq2/",mark, ".teOnly.upup.csv"))
      write.csv(down,paste0(bath,"DESeq2/",mark, ".teOnly.down.csv"))
      #
      teGene <- DESeq2::DESeqDataSetFromMatrix(countData = teGene, colData = dfclas, design = ~Class)
      teGene <- teGene[rowMeans(BiocGenerics::counts(teGene)) > 5,]
      teGene <- DESeq2::DESeq(teGene,)
      teGene <- DESeq2::results(teGene, contrast = VS)
      teGene <- teGene[grep(rownames(teGene), pattern = ":", invert = F),]
      write.csv(teGene, paste0(bath,"DESeq2/",mark, ".te.GeneBackGround.csv"))
      upup <- subset(teGene, log2FoldChange > log2(1.5) & padj < 0.05)
      down <- subset(teGene, log2FoldChange < -log2(1.5)& padj < 0.05)
      write.csv(upup,paste0(bath,"DESeq2/",mark, ".te.GeneBackGround.upup.csv"))
      write.csv(down,paste0(bath,"DESeq2/",mark, ".te.GeneBackGround.down.csv"))
      #
      #
      tmp2 <- read.table(siteFile, header = T, check.names = F, stringsAsFactors = F, row.names = 1)
      colnames(tmp2) <- sapply(str_split(basename(colnames(tmp2)),pattern = ".bam"), "[", 1)
      colLen <- ceiling(length(colnames(tmp2))/2)
      #
      teGene <- tmp2
      teGene <- DESeq2::DESeqDataSetFromMatrix(countData = teGene, colData = dfclas, design = ~Class)
      teGene <- teGene[rowMeans(BiocGenerics::counts(teGene)) > 5,]
      teGene <- DESeq2::DESeq(teGene,)
      teGene <- DESeq2::results(teGene, contrast = VS)
      teGene <- teGene[grep(rownames(teGene), pattern = ":", invert = F),]
      write.csv(teGene, paste0(bath,"DESeq2/",mark, ".te.site.GeneBackGround.csv"))
      upup <- subset(teGene, log2FoldChange > log2(1.5) & padj < 0.05)
      down <- subset(teGene, log2FoldChange < -log2(1.5)& padj < 0.05)
      write.csv(upup, paste0(bath,"DESeq2/",mark, ".te.site.GeneBackGround.upup.csv"))
      write.csv(down, paste0(bath,"DESeq2/",mark, ".te.site.GeneBackGround.down.csv"))
      #
      teGene <- tmp2[grep(rownames(tmp2), pattern = ":", invert = F),]
      teGene <- DESeq2::DESeqDataSetFromMatrix(countData = teGene, colData = dfclas, design = ~Class)
      teGene <- teGene[rowMeans(BiocGenerics::counts(teGene)) > 5,]
      teGene <- DESeq2::DESeq(teGene,)
      teGene <- DESeq2::results(teGene, contrast = VS)
      write.csv(teGene, paste0(bath,"DESeq2/",mark, ".te.site.teOnly.csv"))
      upup <- subset(teGene, log2FoldChange > log2(1.5) & padj < 0.05)
      down <- subset(teGene, log2FoldChange < -log2(1.5)& padj < 0.05)
      write.csv(upup, paste0(bath,"DESeq2/",mark, ".te.site.teOnly.upup.csv"))
      write.csv(down, paste0(bath,"DESeq2/",mark, ".te.site.teOnly.down.csv"))
    }
    path = "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/thejjuuzzRY2/Chi/"
    bath = paste0(path,"testScale/DESeq2Self/")
    shy_diff(mark = "jjuuzz",
             geneFile = paste0(path, "count/jjuuzz.cntTable"),
             siteFile = paste0(path, "count/jjuuzz.site.cntTable"),
             Class = factor(c("OE","OE","WT","WT")), 
             Batch = factor(c("OE_1","OE_2","WT","WT")),
             VS = c("Class","OE","WT"))
  }
}
#04、观察火山图TE差异表达
{
  shy_volcano_te <- function(resFile, saveFile){
    resdata <- as.data.frame(read.csv(resFile, header = T, stringsAsFactors = F))
    resdata <- resdata[grep(":", resdata$X, invert = F),]
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
  shy_volcano_te(resFile  = paste0(bath, "DESeq2/jjuuzz.te.GeneBackGround.csv"),
                 saveFile = paste0(bath, "PLOT/jjuuzz.te.GeneBackGround.volcano.pdf"))
  shy_volcano_te(resFile  = paste0(bath, "DESeq2/jjuuzz.teOnly.csv"),
                 saveFile = paste0(bath, "PLOT/jjuuzz.teOnly.volcano.pdf"))
}
#
#
#
#
#
#
################################################################################
#老师说用Quantile有基因背景的去做图
{
  rm(list=ls())
  for (i in c("data.table","ggridges","ggrepel","KEGG.db","scales","stringr",
              "Biostrings","tidyr","pheatmap","enrichplot","DESeq2","patchwork","ggpubr","ggsci","ggtranscript","rtracklayer","Mfuzz","BiocManager","dplyr","IRanges","GenomicRanges","ReactomePA","ggplot2","clusterProfiler","topGO","Rgraphviz","pathview","org.Mm.eg.db")){library(i, character.only = T,quietly = T)}; rm(i)
  #path = "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/thejjuuzzRY2/Chi/testScale/Quantile/"
  path = "/ChIP_seq_2/aaSHY/ghyuu/rnaSEQ/thejjuuzzRY2/Chi/testScale/Quantile/"
}
#06、Volcano for Genes
################################################################################
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
      scale_x_continuous(breaks = seq(0,20,by=1)) +
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
      geom_vline(aes(xintercept = 10),  color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 11),  color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 12),  color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 13),  color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 14),  color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 15),  color = "white", linewidth  = 0.8) +
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
#15、MERVL sites的Boxplot
{
  prex <- c("OE-1","OE-2","WT-1","WT-2")
  mervs <- openxlsx::read.xlsx(paste0(path,"testScale/Quantile/jjuuzz.cpm.onlyTE.site.xlsx"))
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
           filename = paste0(path,"testScale/Quantile/PLOT/",mark,".boxplot.pdf"),
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
#
#
#
#
#
#
#
################################################################################
{
  prex <- c("OE-1","OE-2","WT-1","WT-2")
  mervs <- openxlsx::read.xlsx(paste0(path,"testScale/Quantile/jjuuzz.cpm.onlyTE.site.xlsx"))
  colnames(mervs) <- c(prex,"theName")
  gggg <- data.frame(TE = sapply(str_split(mervs$theName,":"),"[",2),
                     family = sapply(str_split(mervs$theName,":"),"[",3),
                     class = sapply(str_split(mervs$theName,":"),"[",4))
  gggg <- as.data.frame(distinct(gggg))
  #
  #
  #
  #
  mark = "MERVL-int"
  grou = "jjuuzz"
  mervl <- mervs
  mervl <- mervl[grep(mark,sapply(str_split(mervl$theName,":"),"[",2)),]
  ### mervl <- mervl[rowSums(mervl[,1:4]==0) < 4,]
  mervl <- mervl[rowMeans(mervl[,1:4]) > 5,]
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
    labs(title = paste0(mark," site Expr (MT2-",grou," OE) \n[",
                        "P-value is ",as.character(vaue$p.value),"]\n"),
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
         filename = paste0(path,"testScale/Quantile/PLOT/",mark,".boxplot-revNsmb.pdf"),
         units = "cm", width = 15, height = 15)
  #######
  #
  #
  #
  #######
  mervl[1,]
  sed1 <- mervl %>%
    tidyr::pivot_wider(names_from = group, values_from = cpms, id_cols = theName)
  sed1 <- sed1[,c(1,4,5,2,3)]
  sed2 <- as.data.frame(sed1[,2:5])
  rownames(sed2) <- sed1$theName
  pdf(file = paste0(path,"testScale/Quantile/PLOT/",mark,".pheatmap-revNsmb.pdf"),
      height = 6, width = 6, family = "serif")
  pheatmap(mat = sed2,
           main="\n-------MT2-jjuuzz OE---------",
           scale = "row", cellwidth = 20, 
           show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = F,
           labels_col = colnames(sed2), angle_col = "45",
           breaks = seq(-2,2,by=0.01), border = F,
           color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                      colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
  dev.off()
}
















