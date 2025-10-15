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
  for (i in c("src","fq","Log","stringtie/count","bamCoverage","trim","DESeq2","dumpROOM","bam","count","PLOT")) {
    if (dir.exists(paths = paste0(path, i))==F) {
      dir.create(path = paste0(path,i), recursive = T)
    }
  }
  ff01 <- list.files(path = paste0(path,"fq"),pattern = "_1.fq.gz")
  prex <- sapply(str_split(ff01,pattern = "_1.fq.gz"),"[",1)
}
################################################################################
#
#
#
#
#
#
################################################################################
{
  mervs <- openxlsx::read.xlsx(paste0(path,
                                      "countUniq/TEtools/TEtools_huhuhu.cpm.onlyTE.site.xlsx"))
  grou = "group1"
  mervl <- mervs[,grep(paste0(grou,"|theName"),colnames(mervs))]
  mervl <- mervl[grep("MERVL-int",mervl$theName),]
  ### mervl <- mervl[rowSums(mervl[,1:4]==0) < 4,]
  mervl <- mervl[rowMeans(mervl[,1:4]) > 5,]
  mervl <- tidyr::pivot_longer(data = mervl, cols = -theName, names_to = "group", values_to = "cpms")
  vaue <- wilcox.test(x = mervl$cpms[grep("0000|0180",mervl$group)], 
                      y = mervl$cpms[grep("4500",mervl$group)], alternative = "less")
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
    scale_fill_manual(values = c("#3FA0C0","#F0A19A"),
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
         paste0(path,"PLOT_uniq/","TEtools.uniq.mervl-int.boxplot.group1-revNsmb.pdf"),
         units = "cm", width = 15, height = 15)
  #######
  #
  #
  #
  #######
  red1 <- mervl %>%
    tidyr::pivot_wider(names_from = group, values_from = cpms, id_cols = theName)
  red1 <- red1[,c(1,4,5,2,3)]
  colnames(red1)[2:5] <- c("Highhuhuhucose-1","Highhuhuhucose-2",
                           "Nonehuhuhucose-1","Nonehuhuhucose-2")
  red2 <- as.data.frame(red1[,2:5])
  rownames(red2) <- red1$theName
  pdf(file = paste0(path,"PLOT_uniq/","TEtools.uniq.mervl-int.pheatmap.group1-revNsmb.pdf"), 
      height = 6, width = 6, family = "serif")
  pheatmap(mat = red2,
           main="\n-------Group1---------",
           scale = "row", cellwidth = 20, 
           show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = F,
           labels_col = colnames(red2), angle_col = "45",
           breaks = seq(-2,2,by=0.01), border = F,
           color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                      colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
  dev.off()
}
################################################################################
{
  sets <- list(c("group4.ghyuu.P0.4500","group2.E14.P0.4500"))
  grou = 1
  mervl <- cbind(mervs[,grep(sets[[grou]][1],colnames(mervs))],
                 mervs[,c(grep(sets[[grou]][2],colnames(mervs)),21)])
  mervl <- mervl[grep("MERVL-int",mervl$theName),]
  ###mervl <- mervl[rowSums(mervl[,1:4]==0) < 4,]
  mervl <- mervl[rowMeans(mervl[,1:4]) > 5,]
  mervl <- tidyr::pivot_longer(data = mervl, cols = -theName, 
                               names_to = "group", values_to = "cpms")
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
    scale_fill_manual(values = c("#3FA0C0","#F0A19A"),
                      labels = c("control", "treat")) +
    theme(plot.margin = margin(unit = "cm",1,1,1,1),
          legend.text = element_text(family = "serif", size = 12),
          legend.title = element_text(family = "serif", size = 12),
          axis.text.x = element_text(angle = 45,hjust = 1,family = "serif",size = 12),
          axis.text.y = element_text(family = "serif", size = 12)) +
    scale_y_continuous(limits = c(-5,10)) +
    #scale_y_continuous(limits = c(0,10),breaks = c(0,1,2,3,4,5)) +
    scale_x_discrete(labels = c("WT ESC+huhuhuc", "ghyuu-/- ESC+huhuhuc"))
  print(p_003)
  ggsave(plot = p_003,
         filename = paste0(path,"PLOT_uniq/",
                           "TEtools.uniq.mervl-int.boxplot.add",grou,"-revNsmb.pdf"),
         units = "cm", width = 15, height = 15)
  #######
  #
  #
  #
  #######
  mervl[1,]
  red3 <- mervl %>%
    tidyr::pivot_wider(names_from = group, values_from = cpms, id_cols = theName)
  red3 <- red3[,c(1,4,5,2,3)]
  colnames(red3)[2:5] <- c("WT-ESC+huhuhuc-1","WT-ESC+huhuhuc-2",
                           "ghyuu-KO+huhuhuc-1","ghyuu-KO+huhuhuc-2")
  red4 <- as.data.frame(red3[,2:5])
  rownames(red4) <- red3$theName
  pdf(file = paste0(path,"PLOT_uniq/","TEtools.uniq.mervl-int.pheatmap.add1-revNsmb.pdf"), 
      height = 6, width = 6, family = "serif")
  pheatmap(mat = red4,
           main="\n-------ghyuu KO---------",
           scale = "row", cellwidth = 20, 
           show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = F,
           labels_col = colnames(red4), angle_col = "45",
           breaks = seq(-2,2,by=0.01), border = F,
           color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                      colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
  dev.off()
}


































































