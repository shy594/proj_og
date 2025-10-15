#write by Sun Haiayng at 2024.12.30
#01、搭环境
{
  rm(list=ls())
  for (i in c("data.table","ggridges","ggrepel","KEGG.db","scales","stringr","Biostrings","tidyr","pheatmap","enrichplot","DESeq2","patchwork","ggpubr","ggsci","ggtranscript","rtracklayer","Mfuzz","BiocManager","dplyr","IRanges","GenomicRanges","ReactomePA","ggplot2","clusterProfiler","topGO","Rgraphviz","pathview","org.Mm.eg.db")){library(i, character.only = T,quietly = T)}; rm(i)
}
#02、给索引
{
  path = "/sc/aaSHY/ribo/riboITP/"
  inx1 = "/Reference/aaSHY/Ensembl99/GRCm38/INDEX/STARv279a"
  inx2 <- "/Reference/aaSHY/Ensembl99/GRCm38/fa/GRCm38.dna.primary_assembly.fa"
  inx3 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/Mus_musculus.GRCm38.99.gtf"
  inx4 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf"
  inx5 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf.locInd"
}
#03、做分析
{
  #MM----融合基因
  gtf1 <- import(inx4); gtf1 <- gtf1[grep("MM",gtf1$gene_id)]
  gtf2 <- import(inx3); gtf2_G <- gtf2[which(gtf2$type=="gene")]
  bath <- "/ChIP_seq_1/aaSHY/rna2014Deng/myRun/"
  prex <- c("twoEarly","twoMid","twoLate")
  gggg <- c()
  for (ii in prex) {
    print(ii)
    asem <- import(paste0(bath,"stringtie/",ii,".merge.gtf"))
    chm1 <- asem[which(asem$type=="exon")]
    ov_1 <- suppressWarnings(findOverlaps(chm1,gtf1,ignore.strand=T))
    idss <- unique(chm1$transcript_id[unique(queryHits(ov_1))])
    chm2 <- asem[which(asem$type=="transcript" & asem$transcript_id %in% idss)]
    ov_2 <- suppressWarnings(findOverlaps(chm2,gtf2_G,ignore.strand=T))
    gggg <- c(gggg,unique(gtf2_G$gene_name[subjectHits(ov_2)]))
  }
  aims <- unique(gggg)
  #
  #
  #
  shet <- c("2cell_to_1cell","4cell_to_2cell"); mark = shet[2]
  red1 <- openxlsx::read.xlsx("/sc/aaSHY/ribo/riboITP/suppl/tableS3.xlsx", sheet = mark)#"4cell_to_2cell"
  red1$transcript <- sapply(str_split(red1$transcript,"-"),"[",1)
  red2 <- red1 %>% 
    dplyr::group_by(transcript) %>%
    dplyr::summarise(theTEs = mean(TE_log2FoldChange),
                     theRNA = mean(RNA_log2FoldChange), .groups = "drop")
  red3 <- red2[rowSums(is.na(red2))==0,]
  red3$clas <- "theOthers"
  red3$clas[which(red3$transcript %in% aims)] <- "MMFusions"
  ggplot() +
    labs(title = mark, x = "RNA_Log2FC", y = "TE_Log2FC") +
    geom_point(data = red3[which(red3$clas=="theOthers"),], 
               aes(x = theRNA, y = theTEs), color = "grey") +
    geom_point(data = red3[which(red3$clas=="MMFusions"),], 
               aes(x = theRNA, y = theTEs), color = "red") +
    scale_x_continuous(limits = c(-12,12)) +
    scale_y_continuous(limits = c(-8,8)) +
    theme(plot.margin = margin(unit = "cm",1,1,1,1))
}















