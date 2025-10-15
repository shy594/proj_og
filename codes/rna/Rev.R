#write by Sun Haiayng at 2025.05.08
#01、搭环境
{
  rm(list=ls())
  for (i in c("data.table","stringr","EDASeq","dplyr","tidyr","RUVSeq","ggridges","ggrepel","KEGG.db","scales","stringr","Biostrings","tidyr","pheatmap","enrichplot","DESeq2","patchwork","ggpubr","ggsci","ggtranscript","rtracklayer","Mfuzz","BiocManager","dplyr","IRanges","GenomicRanges","ReactomePA","ggplot2","clusterProfiler","topGO","Rgraphviz","pathview","org.Mm.eg.db")){library(i, character.only = T,quietly = T)}; rm(i)
}
#02、给索引
{
  path = "/disk5/aaSHY/ghyuu/publicData/ppnn/"
  inx1 = "/Reference/aaSHY/Ensembl99/GRCm38/INDEX/STARv279a"
  inx2 = "/Reference/aaSHY/zOther/Rn45s/INDEX/bowtie2/Rn45s"
  inx3 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/Mus_musculus.GRCm38.99.gtf"
  inx4 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf"
  inx5 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf.locInd"
  inx6 = "/Reference/aaSHY/Ensembl99/GRCm38/BED/TEsubFamily/MERVL-int::ERVL::LTR.bed"
  inx7 = "/Reference/aaSHY/Ensembl99/GRCm38/BED/TEsubFamily/MT2_Mm::ERVL::LTR.bed"
  gtf1 <- rtracklayer::import(inx3)
  gtf2 <- as.data.frame(gtf1[which(gtf1$type == "gene")])
  gtf2 <- gtf2[,10:12]
}
################################################################################
#03、下数据 [先跑这个]
{
  for (i in c("src","fq","Log","bamCoverage","trim",
              "rDNA_Del","DESeq2","dumpROOM","bam","count","PLOT")) {
    if (dir.exists(paths = paste0(path, i))==F) {
      dir.create(path = paste0(path,i), recursive = T)
    }
  }
  shell <- paste0(path,"src/run_download.sh")
  cat("#!/bin/bash\n", file = shell)
  cmd_01 <- paste0("#wget -c -nc -i ", path,"src/ENA_Acc_List.txt ",
                   "-P ",path,"rawdata","\n")
  for (i in cmd_01) {cat(i, file = shell, append = T)}
  print(paste0("nohup bash ",shell, " >",path,"Log/",basename(shell),".log"," 2>&1 &"))
  #
  lup1 <- openxlsx::read.xlsx(paste0(path,"aa-README/data_access.xlsx"), sheet = "Selected")
  lup1 <- lup1[,1:4]
  colnames(lup1) <- c("samp","cond","stag","link")
  nam1 <- paste0(path,"rawdata/",gsub("fastq","fq",basename(lup1$link)))
  nam2 <- paste0(path,"rawdata/",
                 gsub("2-cell","twos",lup1$stag),"_",
                 sapply(str_split(lup1$cond,"_"),"[",1),"_",lup1$samp,"_",
                 rep(c("1","2"),length(nam1)/2),".fq.gz")
  #file.rename(from = nam1, to = nam2)
}
################################################################################
#04、根据作者提供的附件的基因水平的reads count做差异分析 [不去除批次效应]
{
  ### DESeq2的差异分析
  {
    lup2 <- lup1[,1:3]
    lup2$name <- paste0(gsub("2-cell","twos",lup1$stag),
                        "_",sapply(str_split(lup1$cond,"_"),"[",1))
    cup1 <- read.csv(paste0(path,"suppl/counts_table_all_stages.csv"), 
                     header = T, row.names = 1)
    cup2 <- cup1[,colnames(cup1)[which(colnames(cup1) %in% unique(lup2$samp))]]
    meta_01 <- data.frame(samp = colnames(cup2), clas = ".")
    meta_01 <- left_join(meta_01, distinct(lup2[,c(1,4)]), by = "samp")
    meta_01 <- na.omit(meta_01)
    meta_02 <- data.frame(row.names = meta_01$samp, Class = meta_01$name, stringsAsFactors = T)
    meta_02 <- tidyr::separate(meta_02, col = Class, sep = "_", into = c("stag","cond"))
    #gtf1 <- import(inx3)
    #gtf2 <- as.data.frame(gtf1[which(gtf1$type == "gene")])
    #gtf2 <- gtf2[,10:12]
    for (i in c("twos", "morula")) {
      print(i)
      for (j in c("dBtgh","none")) {
        print(j)
        vsvs = c("cond","Btgh",j)
        meta_03 <- meta_02[grep(i, x = meta_02$stag),]
        meta_03 <- meta_03[which(meta_03$cond %in% c("Btgh",j)),]
        meta_03$cond <- as.factor(meta_03$cond)
        cup3 <- cup2[,rownames(meta_03)]
        cup3 <- cup3 %>% 
          mutate(across(where(is.character) | where(is.factor), readr::parse_number))
        #
        diff_01 <- DESeq2::DESeqDataSetFromMatrix(countData = cup3, 
                                                  colData = meta_03, design = ~cond)
        diff_01 <- diff_01[rowSums(counts(diff_01)) >= 10,]
        ### diff_01 <- diff_01[rowMeans(BiocGenerics::counts(diff_01)) > 5,]
        diff_01 <- DESeq2::DESeq(diff_01,)
        diff_01 <- DESeq2::results(diff_01, contrast = vsvs)
        diff_02 <- as.data.frame(diff_01)
        diff_02$gene_id <- rownames(diff_02)
        diff_02 <- left_join(diff_02, gtf2[,c(1,3)], by = "gene_id")
        diff_02 <- diff_02[,c(7,8,1:6)]
        write.csv(diff_02, 
                  paste0(path,"suppl/DESeq2_",i,"_","Btgh-vs-",j, ".geneOnly.csv"),
                  quote = F, row.names = F)
        upup <- subset(diff_02, log2FoldChange > log2(1.5) & padj < 0.05)
        down <- subset(diff_02, log2FoldChange < -log2(1.5)& padj < 0.05)
        write.csv(upup, 
                  paste0(path,"suppl/DESeq2_",i,"_","Btgh-vs-",j, ".geneOnly.upup.csv"),
                  quote = F, row.names = F)
        write.csv(down, 
                  paste0(path,"suppl/DESeq2_",i,"_","Btgh-vs-",j, ".geneOnly.down.csv"),
                  quote = F, row.names = F)
      }
    }
  }
  ### GSEA
  {
    resd <- list.files(path = paste0(path, "suppl"),pattern = "csv", full.names = T)[c(8,14)]
    for (iu in resd) {
      temp <- read.csv(iu, header = T)
      temp <- temp[, c("gene_name", "log2FoldChange")]
      temp <- temp %>% 
        arrange(desc(log2FoldChange)) %>%
        distinct(gene_name, .keep_all = T)
      temp <- temp[order(temp$log2FoldChange, decreasing = T),]
      temp <- na.omit(temp)
      
      write.table(x = temp, 
                  file = paste0(path,"suppl/",gsub(".csv","",basename(iu)),".rnk"), 
                  quote = F, sep = "\t", col.names = F, row.names = F)
      cmd_01 <- paste0("gseapy prerank -g ",
                       "/Reference/aaSHY/GSEAgmt/TwoGenes.gmt -r ",
                       paste0(path,"suppl/",gsub(".csv","",basename(iu)),".rnk "),
                       "-p 10 -o ",
                       paste0(path, "suppl/enrichGSEA/"))#conda activate base gseapy 0.10.8
      cmd_02 <- paste0("gseapy prerank -g ",
                       "/Reference/aaSHY/GSEAgmt/2022_majorZGA.gmt -r ",
                       paste0(path,"suppl/",gsub(".csv","",basename(iu)),".rnk "),
                       "--max-size 1000 -p 10 -o ",
                       paste0(path, "suppl/enrichGSEA/"))#conda activate base gseapy 0.10.8\
      print(cmd_01)
      print(cmd_02)
    }
  }
}
################################################################################
#05、根据作者提供的附件的基因水平的reads count做差异分析 [RUVSeq去除批次效应]
{
  ### DESeq2的差异分析
  {
    des1 <- read.table(paste0(path,"suppl/twos_batch.txt"), sep = "\t", header = T)
    des1 <- des1[,c(2,3,5,6)]
    colnames(des1)[1] <- "samp"
    lup2 <- distinct(lup1[which(lup1$stag=="2-cell"),1:3])
    lup2$name <- paste0(gsub("2-cell","twos",lup2$stag),
                        "_",sapply(str_split(lup2$cond,"_"),"[",1))
    lup2 <- left_join(lup2, des1, by = "samp")
    #
    cup1 <- read.csv(paste0(path,"suppl/counts_table_all_stages.csv"), 
                     header = T, row.names = 1)
    cup2 <- cup1[,colnames(cup1)[which(colnames(cup1) %in% unique(lup2$samp))]]
    meta_01 <- lup2
    meta_01 <- na.omit(meta_01)
    meta_01 <- meta_01[which(meta_01$samp %in% colnames(cup2)),]
    rownames(meta_01) <- meta_01$samp; meta_01 <- meta_01[,-c(3,4,5)]
    
    meta_02 <- data.frame(row.names = colnames(meta_01), 
                          labelDescription = c("sampleID","conditions","batch","replicates"))
    meta_03 <- Biobase::AnnotatedDataFrame(data = meta_01, 
                                           varMetadata = meta_02)
    meta_03 <- meta_03[,c(2,3,4)]
    rep1 <- pData(meta_03)
    rep2 <- c()
    lens <- max(table(rep1$biol_rep))
    for (i in unique(rep1$biol_rep)) {
      tmp1 <- rownames(rep1)[which(rep1$biol_rep==i)]
      tmp2 <- if (length(tmp1) < lens) {
        c(tmp1, rep("-1",lens-length(tmp1)))
      } else {tmp1}
      rep2 <- c(rep2, tmp2)
    }
    rep3 <- matrix(rep2,byrow = T, ncol = lens)
    cup2 <- cup2[,meta_01$samp]
    cup3 <- cup2 %>% 
      mutate(across(where(is.character) | where(is.factor), readr::parse_number))
    cup3 <- as.matrix(cup3)
    cup4 <- EDASeq::newSeqExpressionSet(counts = cup3, phenoData = meta_03)
    ruvs <- RUVSeq::RUVs(cup4, cIdx = rownames(cup3), k = 3, scIdx = rep3)
    ruvs <- data.frame(W_1 = ruvs$W_1, W_2 = ruvs$W_2, W_3 = ruvs$W_3)
    meta_04 <- meta_01
    lens_02 <- length(colnames(meta_04))
    meta_04[,c(lens_02+1):c(lens_02 + 3)] <- ruvs
    meta_04$cond <- as.factor(meta_04$cond)
    meta_04$cond <- gsub("_injected","",meta_04$cond)
    #
    #
    #
    #
    #gtf1 <- rtracklayer::import(inx3)
    #gtf2 <- as.data.frame(gtf1[which(gtf1$type == "gene")])
    #gtf2 <- gtf2[,10:12]
    for (j in c("dBtgh","none")) {
      print(j)
      vsvs = c("cond","Btgh",j)
      meta_05 <- meta_04[which(meta_04$cond %in% c("Btgh",j)),]
      meta_05$cond <- as.factor(meta_05$cond)
      
      cup4 <- cup2[,rownames(meta_05)]
      cup4 <- cup4 %>% 
        mutate(across(where(is.character) | where(is.factor), readr::parse_number))
      #
      diff_01 <- DESeqDataSetFromMatrix(countData = cup4, 
                                        colData = meta_05, design = ~cond + W_1 + W_2 + W_3)
      diff_01 <- diff_01[rowSums(counts(diff_01)) >= 10,]
      ### diff_01 <- diff_01[rowMeans(BiocGenerics::counts(diff_01)) > 5,]
      diff_01 <- DESeq2::DESeq(diff_01,)
      diff_01 <- DESeq2::results(diff_01, contrast = vsvs)
      diff_02 <- as.data.frame(diff_01)
      diff_02$gene_id <- rownames(diff_02)
      diff_02 <- left_join(diff_02, gtf2[,c(1,3)], by = "gene_id")
      diff_02 <- diff_02[,c(7,8,1:6)]
      write.csv(diff_02, 
                paste0(path,"suppl/DESeq2_twos_","Btgh-vs-",j, ".geneOnly.rmBatch.csv"),
                quote = F, row.names = F)
      upup <- subset(diff_02, log2FoldChange > log2(1.5) & padj < 0.05)
      down <- subset(diff_02, log2FoldChange < -log2(1.5)& padj < 0.05)
      write.csv(upup, 
                paste0(path,"suppl/DESeq2_twos_","Btgh-vs-",j, ".geneOnly.upup.rmBatch.csv"),
                quote = F, row.names = F)
      write.csv(down, 
                paste0(path,"suppl/DESeq2_twos_","Btgh-vs-",j, ".geneOnly.down.rmBatch.csv"),
                quote = F, row.names = F)
    }
  }
  ### GSEA
  {
    resd <- list.files(path = paste0(path, "suppl"),pattern = "csv", full.names = T)[c(11,17)]
    for (iu in resd) {
      temp <- read.csv(iu, header = T)
      temp <- temp[, c("gene_name", "log2FoldChange")]
      temp <- temp %>% 
        arrange(desc(log2FoldChange)) %>%
        distinct(gene_name, .keep_all = T)
      temp <- temp[order(temp$log2FoldChange, decreasing = T),]
      temp <- na.omit(temp)
      
      write.table(x = temp, 
                  file = paste0(path,"suppl/",gsub(".csv","",basename(iu)),".rnk"), 
                  quote = F, sep = "\t", col.names = F, row.names = F)
      cmd_01 <- paste0("gseapy prerank -g ",
                       "/Reference/aaSHY/GSEAgmt/TwoGenes.gmt -r ",
                       paste0(path,"suppl/",gsub(".csv","",basename(iu)),".rnk "),
                       "-p 10 -o ",
                       paste0(path, "suppl/enrichGSEA/"))#conda activate base gseapy 0.10.8
      cmd_02 <- paste0("gseapy prerank -g ",
                       "/Reference/aaSHY/GSEAgmt/2022_majorZGA.gmt -r ",
                       paste0(path,"suppl/",gsub(".csv","",basename(iu)),".rnk "),
                       "--max-size 1000 -p 10 -o ",
                       paste0(path, "suppl/enrichGSEA/"))#conda activate base gseapy 0.10.8\
      print(cmd_01)
      print(cmd_02)
    }
  }
}
#
#
#
#
#
#
################################################################################
#06、跑命令 [for TE alignment]
{
  shell <- paste0(path,"src/run_a.sh")
  cat("#!/bin/bash\n", file = shell)
  prex <- unique(gsub("_[0-9].fq.gz","",basename(nam2)))
  lay1 <- paste0(path,"rawdata/",prex,"_1.fq.gz")
  lay2 <- paste0(path,"rawdata/",prex,"_2.fq.gz")
  cmd_01 <- paste0("#trim_galore --phred33 -j 4 ",
                   "-o ", path, "trim --paired ",lay1," ",lay2,"\n")
  for (i in cmd_01) {cat(i, file = shell, append = T)}
  #
  #
  #
  for (K in c("twos","morula")) {
    shell <- paste0(path,"src/run_",K,".sh")
    cat("#!/bin/bash\n", file = shell)
    pre1 <- prex[grep(K,prex)]
    lay1 <- paste0(path,"trim/",pre1,"_1_val_1.fq.gz")
    lay2 <- paste0(path,"trim/",pre1,"_2_val_2.fq.gz")
    cmd_02 <- paste0("#bowtie2 -p 20 --very-fast --no-mixed --no-discordant ",
                     "--un-conc-gz ",path,"rDNA_Del/",pre1,"_un_%.fq.gz ",
                     "-x ", inx2, " -1 ", lay1, " -2 ",lay2, " >/dev/null","\n")
    #for (i in cmd_02) {cat(i, file = shell, append = T)}
    lay1 <- paste0(path,"rDNA_Del/",pre1,"_un_1.fq.gz")
    lay2 <- paste0(path,"rDNA_Del/",pre1,"_un_2.fq.gz")
    cmd_03 <- paste0("STAR --runThreadN 10 --readFilesCommand zcat ",
                     "--genomeDir ",inx1," --readFilesIn ",lay1," ",lay2," ",
                     "--outFileNamePrefix ",path,"bam/",pre1," ",
                     "--runMode alignReads --alignEndsType EndToEnd ",
                     "--outSAMtype BAM Unsorted --outMultimapperOrder Random ",
                     "--genomeLoad LoadAndKeep ",
                     "--outFilterMultimapNmax 5000 --outSAMmultNmax 1 --alignIntronMax 1 ",
                     "--outFilterMismatchNmax 3 --alignMatesGapMax 350 ",
                     "--winAnchorMultimapNmax 5000 ",
                     "--seedSearchStartLmax 30 --alignTranscriptsPerReadNmax 30000 ",
                     "--alignWindowsPerReadNmax 30000 --alignTranscriptsPerWindowNmax 300 ",
                     "--seedPerReadNmax 3000 --seedPerWindowNmax 300 --seedNoneLociPerWindow 1000 ",
                     "--outReadsUnmapped None --outSAMunmapped None","\n")
    cmd_04 <- paste0("samtools sort -@ 10 -o ",path,"bam/",pre1,
                     ".bam ",path,"bam/",pre1,"Aligned.out.bam","\n")
    cmd_05 <- paste0("samtools index ", path,"bam/",pre1,".bam","\n")
    cmd_06 <- paste0("#rm ", path,"bam/", pre1,"Aligned.out.bam","\n")
    for (i in cmd_03) {cat(i, file = shell, append = T)}
    for (i in cmd_04) {cat(i, file = shell, append = T)}
    for (i in cmd_05) {cat(i, file = shell, append = T)}
    for (i in cmd_06) {cat(i, file = shell, append = T)}
    print(paste0("nohup bash ",shell, " >",path,"Log/",basename(shell),".log"," 2>&1 &"))
  }
}
################################################################################
#07、跑命令 [for TE counts、differential analysis]
{
  ### 制作符合要求的SAF
  {
    tef1 <- as.data.frame(import(inx4))
    tef2 <- tef1 %>%
      dplyr::filter(!(startsWith(gene_id, "L1")  & width <= 5000))
    tef3 <- tef2 %>%
      dplyr::filter(!(startsWith(gene_id, "IAP") & width <= 6000))
    tef4 <- tef3 %>%
      dplyr::filter(!(startsWith(gene_id, "MMERVK10C") & width <= 4500))
    tef5 <- tef4 %>%
      dplyr::filter(seqnames %in% c(seq(19),"X","Y")) %>%
      droplevels()
    tef6 <- makeGRangesFromDataFrame(df = tef5, keep.extra.columns = T)
    gtf1 <- import(inx3)
    ov_1 <- suppressWarnings(findOverlaps(tef6, gtf1, ignore.strand=T))
    tef7 <- tef6[setdiff(seq_along(tef6$source),unique(queryHits(ov_1)))]
    tef7 <- as.data.frame(tef7)
    wr01 <- tef7[,c(10,1:3,5)]
    colnames(wr01) <- c("GeneID","Chr","Start","End","Strand")
    write.table(x = wr01, 
                file = paste0(path,"count/TE-filtered.saf"),
                sep = "\t", quote = F, col.names = T, row.names = F)
    wr02 <- tef1[,c(10,1:3,5)]
    colnames(wr02) <- c("GeneID","Chr","Start","End","Strand")
    write.table(x = wr02, 
                file = paste0(path,"count/TE-full.saf"),
                sep = "\t", quote = F, col.names = T, row.names = F)
  }
  ### 跑featureCounts
  {
    inx_tmp_1 <- paste0(path,"count/TE-filtered.saf")
    inx_tmp_2 <- paste0(path,"count/TE-full.saf")
    shell <- paste0(path,"src/run_e.sh")
    cat("#!/bin/bash\n", file = shell)
    cmd_07 <- paste0("featureCounts -F SAF -p -B -s 0 --fracOverlap 1 ",
                     "-M -T 16 -a ",inx_tmp_1," -o ",path,
                     "count/TE-filtered.txt ",
                     paste0(path,"bam/",prex,".bam", collapse = " "),"\n")
    cmd_08 <- paste0("featureCounts -F SAF -p -B -s 0 --fracOverlap 1 ",
                     "-M -T 16 -a ",inx_tmp_2," -o ",path,
                     "count/TE-full.txt ",
                     paste0(path,"bam/",prex,".bam", collapse = " "),"\n")
    for (i in cmd_07) {cat(i, file = shell, append = T)}
    for (i in cmd_08) {cat(i, file = shell, append = T)}
    print(paste0("nohup bash ",shell, " >",path,"Log/",basename(shell),".log"," 2>&1 &"))
  }
  ### DESeq2 (2-cell)
  {
    #用meta_04 (在上面05的第一步)
    cte1 <- read.table(paste0(path,"count/TE-full-2.txt"), 
                       header = T, sep = "\t", row.names = 1)
    colnames(cte1) <- sapply(str_split(colnames(cte1),"\\.|_"),"[",10)
    cte2 <- cte1[,meta_04$samp]
    #
    meta_04$cond <- as.factor(meta_04$cond)
    diff_01 <- DESeqDataSetFromMatrix(countData = cte2, 
                                      colData = meta_04, design = ~cond + W_1 + W_2 + W_3)
    diff_01 <- DESeq2::DESeq(diff_01,)
    diff_02 <- counts(diff_01, normalized = T)
    diff_03 <- diff_02[grep("MERVL-int|MT2_Mm",rownames(diff_02)),]
    diff_03 <- as.data.frame(diff_03)
    diff_03$name <- rownames(diff_03)
    diff_04 <- diff_03 %>% 
      tidyr::pivot_longer(cols = 1:c(length(colnames(diff_03))-1), names_to = "samp", values_to = "norm")
    diff_04 <- left_join(diff_04, meta_04[,c(1,2)], by = "samp")
    diff_05 <- as.data.frame(diff_04)
    ggplot(data = diff_05) +
      geom_boxplot(aes(x = cond, y = norm, fill = cond)) +
      scale_fill_few() +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "", y = "normlized reads count",fill ="") +
      theme_classic() +
      theme(plot.margin = margin(unit = "cm", 1,1,1,1),
            legend.text = element_text(family= "serif", size = 14),
            axis.text = element_text(family= "serif", size = 14),
            axis.title = element_text(family= "serif", size = 16),
            strip.text = element_text(family= "serif", size = 16),
            strip.background = element_rect(fill = "gray85"))
    #
    for (j in c("dBtgh","none")) {
      print(j)
      vsvs = c("cond","Btgh",j)
      meta_05 <- meta_04[which(meta_04$cond %in% c("Btgh",j)),]
      meta_05$cond <- as.factor(meta_05$cond)
      
      cte3 <- cte2[,rownames(meta_05)]
      cte3 <- cte3 %>% 
        mutate(across(where(is.character) | where(is.factor), readr::parse_number))
      #
      diff_01 <- DESeqDataSetFromMatrix(countData = cte3, 
                                        colData = meta_05, design = ~cond + W_1 + W_2 + W_3)
      diff_01 <- diff_01[rowSums(counts(diff_01)) >= 10,]
      ### diff_01 <- diff_01[rowMeans(BiocGenerics::counts(diff_01)) > 5,]
      diff_01 <- DESeq2::DESeq(diff_01,)
      diff_01 <- DESeq2::results(diff_01, contrast = vsvs)
      diff_02 <- as.data.frame(diff_01)
      diff_02$gene_id <- rownames(diff_02)
      diff_02 <- left_join(diff_02, gtf2[,c(1,3)], by = "gene_id")
      diff_02 <- diff_02[,c(7,8,1:6)]
      write.csv(diff_02, 
                paste0(path,"suppl/DESeq2_twos_","Btgh-vs-",j, ".teOnly.rmBatch.csv"),
                quote = F, row.names = F)
      upup <- subset(diff_02, log2FoldChange > log2(1.5) & padj < 0.05)
      down <- subset(diff_02, log2FoldChange < -log2(1.5)& padj < 0.05)
      write.csv(upup, 
                paste0(path,"suppl/DESeq2_twos_","Btgh-vs-",j, ".teOnly.upup.rmBatch.csv"),
                quote = F, row.names = F)
      write.csv(down, 
                paste0(path,"suppl/DESeq2_twos_","Btgh-vs-",j, ".teOnly.down.rmBatch.csv"),
                quote = F, row.names = F)
    }
    
    
    

  }
  ### DESeq2 (morula)
  {
    lup2 <- lup1[,1:3]
    lup2$name <- paste0(gsub("2-cell","twos",lup1$stag),
                        "_",sapply(str_split(lup1$cond,"_"),"[",1))
    cup1 <- read.csv(paste0(path,"suppl/counts_table_all_stages.csv"), 
                     header = T, row.names = 1)
    cup2 <- cup1[,colnames(cup1)[which(colnames(cup1) %in% unique(lup2$samp))]]
    meta_01 <- data.frame(samp = colnames(cup2), clas = ".")
    meta_01 <- left_join(meta_01, distinct(lup2[,c(1,4)]), by = "samp")
    meta_01 <- na.omit(meta_01)
    meta_02 <- data.frame(row.names = meta_01$samp, Class = meta_01$name, stringsAsFactors = T)
    meta_02 <- tidyr::separate(meta_02, col = Class, sep = "_", into = c("stag","cond"))
    meta_03 <- meta_02[which(meta_02$stag=="morula"),]
    meta_03$samp <- rownames(meta_03)
    meta_04 <- meta_03
    #
    #
    cte1 <- read.table(paste0(path,"count/TE-full-2.txt"), 
                       header = T, sep = "\t", row.names = 1)
    colnames(cte1) <- sapply(str_split(colnames(cte1),"\\.|_"),"[",10)
    cte2 <- cte1[,meta_04$samp]
    #
    meta_04$cond <- as.factor(meta_04$cond)
    diff_01 <- DESeqDataSetFromMatrix(countData = cte2, 
                                      colData = meta_04, design = ~cond)
    diff_01 <- DESeq2::DESeq(diff_01,)
    diff_02 <- counts(diff_01, normalized = T)
    diff_03 <- diff_02[grep("MERVL-int|MT2_Mm",rownames(diff_02)),]
    diff_03 <- as.data.frame(diff_03)
    diff_03$name <- rownames(diff_03)
    diff_04 <- diff_03 %>% 
      tidyr::pivot_longer(cols = 1:c(length(colnames(diff_03))-1), names_to = "samp", values_to = "norm")
    diff_04 <- left_join(diff_04, meta_04[,c(3,2)], by = "samp")
    diff_05 <- as.data.frame(diff_04)
    ggplot(data = diff_05) +
      geom_boxplot(aes(x = cond, y = norm, fill = cond)) +
      scale_fill_few() +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "", y = "normlized reads count",fill ="") +
      theme_classic() +
      theme(plot.margin = margin(unit = "cm", 1,1,1,1),
            legend.text = element_text(family= "serif", size = 14),
            axis.text = element_text(family= "serif", size = 14),
            axis.title = element_text(family= "serif", size = 16),
            strip.text = element_text(family= "serif", size = 16),
            strip.background = element_rect(fill = "gray85"))
    #
    #
    #
    for (j in c("dBtgh","none")) {
      print(j)
      vsvs = c("cond","Btgh",j)
      meta_05 <- meta_04[which(meta_04$cond %in% c("Btgh",j)),]
      meta_05$cond <- as.factor(meta_05$cond)
      
      cte3 <- cte2[,rownames(meta_05)]
      cte3 <- cte3 %>% 
        mutate(across(where(is.character) | where(is.factor), readr::parse_number))
      #
      diff_01 <- DESeqDataSetFromMatrix(countData = cte3, 
                                        colData = meta_05, design = ~cond)
      diff_01 <- diff_01[rowSums(counts(diff_01)) >= 10,]
      ### diff_01 <- diff_01[rowMeans(BiocGenerics::counts(diff_01)) > 5,]
      diff_01 <- DESeq2::DESeq(diff_01,)
      diff_01 <- DESeq2::results(diff_01, contrast = vsvs)
      diff_02 <- as.data.frame(diff_01)
      diff_02$gene_id <- rownames(diff_02)
      diff_02 <- left_join(diff_02, gtf2[,c(1,3)], by = "gene_id")
      diff_02 <- diff_02[,c(7,8,1:6)]
      write.csv(diff_02, 
                paste0(path,"suppl/DESeq2_morula_","Btgh-vs-",j, ".teOnly.rmBatch.csv"),
                quote = F, row.names = F)
      upup <- subset(diff_02, log2FoldChange > log2(1.5) & padj < 0.05)
      down <- subset(diff_02, log2FoldChange < -log2(1.5)& padj < 0.05)
      write.csv(upup, 
                paste0(path,"suppl/DESeq2_morula_","Btgh-vs-",j, ".teOnly.upup.rmBatch.csv"),
                quote = F, row.names = F)
      write.csv(down, 
                paste0(path,"suppl/DESeq2_morula_","Btgh-vs-",j, ".teOnly.down.rmBatch.csv"),
                quote = F, row.names = F)
    }
  }
}
################################################################################
#
#
#
#
#
#
#08、跑命令 [按本课题组传统做法]
{
  for (i in c("ourLab/count","ourLab/bam")) {
    if (dir.exists(paths = paste0(path, i))==F) {
      dir.create(path = paste0(path,i), recursive = T)
    }
  }
  prex <- unique(gsub("_[0-9].fq.gz","",basename(nam2)))
  for (K in c("twos","morula")) {
    shell <- paste0(path,"src/run_shy_",K,".sh")
    cat("#!/bin/bash\n", file = shell)
    pre1 <- prex[grep(K,prex)]
    lay1 <- paste0(path,"trim/",pre1,"_1_val_1.fq.gz")
    lay2 <- paste0(path,"trim/",pre1,"_2_val_2.fq.gz")
    #cmd_01 <- paste0("bowtie2 -p 20 --very-fast --no-mixed --no-discordant ",
    #                 "--un-conc-gz ",path,"rDNA_Del/",pre1,"_un_%.fq.gz ",
    #                 "-x ", inx2, " -1 ", lay1, " -2 ",lay2, " >/dev/null","\n")
    #for (i in cmd_01) {cat(i, file = shell, append = T)}
    #lay1 <- paste0(path,"rDNA_Del/",pre1,"_un_1.fq.gz")
    #lay2 <- paste0(path,"rDNA_Del/",pre1,"_un_2.fq.gz")
    #cmd_02 <- paste0("STAR --runThreadN 10 --readFilesCommand zcat ",
    #                 "--outMultimapperOrder Random --outSAMtype BAM Unsorted ",
    #                 "--outFileNamePrefix ",paste0(path,"ourLab/bam/",pre1)," ",
    #                 "--outSAMprimaryFlag AllBestScore --twopassMode Basic ",
    #                 "--genomeDir ",inx1," --readFilesIn ",lay1," ",lay2,"\n")
    #cmd_03 <- paste0("samtools sort -@ 10 -o ",path,"ourLab/bam/",pre1,
    #                 ".bam ",path,"ourLab/bam/",pre1,"Aligned.out.bam","\n")
    #cmd_04 <- paste0("samtools index ", path,"ourLab/bam/",pre1,".bam","\n")
    cont <- pre1[1:9]
    trea <- pre1[10:length(pre1)]
    cmd_05 <- paste0("TEtranscripts -t ",
                     paste(paste0(path,"ourLab/bam/",trea,".bam"),collapse = " ")," -c ",
                     paste(paste0(path,"ourLab/bam/",cont,".bam"),collapse = " "),
                     " --GTF ",inx3," --TE ",inx4," --sortByPos --mode multi ",
                     "--project ",K," --outdir ",path,"ourLab/count","\n")
    #for (i in cmd_02) {cat(i, file = shell, append = T)}
    #for (i in cmd_03) {cat(i, file = shell, append = T)}
    #for (i in cmd_04) {cat(i, file = shell, append = T)}
    for (i in cmd_05) {cat(i, file = shell, append = T)}
    print(paste0("nohup bash ",shell, " >",path,"Log/",basename(shell),".log"," 2>&1 &"))
  }
}
#09、跑命令 [for TE differential analysis]
{
  ### DESeq2 (2-cell)
  {
    #用meta_04 (在上面05的第一步)
    meta_04 <- meta_04[which(meta_04$cond!="dBtgh"),]
    cte1 <- read.table(paste0(path,"ourLab/count/twos.cntTable"), 
                       header = T, sep = "\t", row.names = 1)
    cte1 <- cte1[grep(":",rownames(cte1)),]
    colnames(cte1) <- sapply(str_split(colnames(cte1),"\\.|_"),"[",11)
    cte2 <- cte1[,meta_04$samp]
    #
    meta_04$cond <- as.factor(meta_04$cond)
    diff_01 <- DESeqDataSetFromMatrix(countData = cte2, 
                                      colData = meta_04, design = ~cond + W_1 + W_2 + W_3)
    diff_01 <- DESeq2::DESeq(diff_01,)
    diff_02 <- counts(diff_01, normalized = T)
    diff_03 <- diff_02[grep("MERVL-int|MT2_Mm",rownames(diff_02)),]
    diff_03 <- as.data.frame(diff_03)
    diff_03$name <- rownames(diff_03)
    diff_04 <- diff_03 %>% 
      tidyr::pivot_longer(cols = 1:c(length(colnames(diff_03))-1), names_to = "samp", values_to = "norm")
    diff_04 <- left_join(diff_04, meta_04[,c(1,2)], by = "samp")
    diff_05 <- as.data.frame(diff_04)
    #
    #
    #
    #
    #
    aims <- "MERVL-int"
    aims <- "MT2_Mm"
    diff_05 <- diff_04[grep(aims,diff_04$name),]
    {
      vsvs = c("cond","Btgh","none")
      meta_05 <- meta_04[which(meta_04$cond %in% c("Btgh",j)),]
      meta_05$cond <- as.factor(meta_05$cond)
      
      cte3 <- cte2[,rownames(meta_05)]
      cte3 <- cte3 %>% 
        mutate(across(where(is.character) | where(is.factor), readr::parse_number))
      #
      ress_01 <- DESeqDataSetFromMatrix(countData = cte3, 
                                        colData = meta_05, design = ~cond + W_1 + W_2 + W_3)
      #ress_01 <- ress_01[rowSums(counts(ress_01)) >= 10,]
      ress_01 <- ress_01[rowMeans(BiocGenerics::counts(ress_01)) > 5,]
      ress_01 <- DESeq2::DESeq(ress_01,)
      ress_01 <- DESeq2::results(ress_01, contrast = vsvs)
      
      ress_02 <- as.data.frame(ress_01)
      ress_02$gene_id <- rownames(ress_02)
      ress_02 <- left_join(ress_02, gtf2[,c(1,3)], by = "gene_id")
      ress_02 <- ress_02[,c(7,8,1:6)]
      fcfc <- ress_02[grep(aims, ress_02$gene_id), 4]
      padj <- ress_02[grep(aims, ress_02$gene_id), 8]
    }
    #
    #
    #
    #
    #
    #
    P001 <- ggplot(data = diff_05) +
      geom_boxplot(aes(x = cond, y = norm, fill = cond)) +
      scale_fill_manual(values = c("#F0A19A","#3FA0C0")) +
      scale_x_discrete(labels = c("Inject ghyuu into Zygote","Control")) +
      guides(fill = "none") +
      labs(x = "", 
           y = "Normlized reads count",fill ="",
           title = paste0(aims," samples expression level (RNA-seq of 2-Cell stage)")) +
      theme_classic() +
      theme(plot.margin = margin(unit = "cm", 1,1,1,1),
            plot.title  = element_text(family= "serif", size = 18, face = "bold"),
            legend.text = element_text(family= "serif", size = 14),
            axis.text   = element_text(family= "serif", size = 14),
            axis.title  = element_text(family= "serif", size = 16),
            strip.text  = element_text(family= "serif", size = 16),
            strip.background = element_rect(fill = "gray85")) +
      annotate(geom = "text", 
               x = 1.5,
               y = max(diff_05$norm) * 1.05,
               label = paste0("Fold change is: ",round(fcfc,4),
                              "; P-adjust value is: ",round(padj,4)),
               color = "red", 
               size = 5, 
               family = "serif")
    P001
    ggsave(plot = P001, 
           filename = paste0(path,"PLOT_need/2025ppnn-samples.expr.",aims,".pdf"),
           units = "cm", width = 20, height = 18)
  }
  ### DESeq2 (morula)
  {
    #做出morula需要的meta data表
    {
      lup2 <- lup1[,1:3]
      lup2$name <- paste0(gsub("2-cell","twos",lup1$stag),
                          "_",sapply(str_split(lup1$cond,"_"),"[",1))
      cup1 <- read.csv(paste0(path,"suppl/counts_table_all_stages.csv"), 
                       header = T, row.names = 1)
      cup2 <- cup1[,colnames(cup1)[which(colnames(cup1) %in% unique(lup2$samp))]]
      meta_01 <- data.frame(samp = colnames(cup2), clas = ".")
      meta_01 <- left_join(meta_01, distinct(lup2[,c(1,4)]), by = "samp")
      meta_01 <- na.omit(meta_01)
      meta_02 <- data.frame(row.names = meta_01$samp, Class = meta_01$name, stringsAsFactors = T)
      meta_02 <- tidyr::separate(meta_02, col = Class, sep = "_", into = c("stag","cond"))
      meta_03 <- meta_02[which(meta_02$stag=="morula"),]
      meta_03$samp <- rownames(meta_03)
      meta_04 <- meta_03
    }
    #
    #
    cte1 <- read.table(paste0(path,"ourLab/count/morula.cntTable"), 
                       header = T, sep = "\t", row.names = 1)
    cte1 <- cte1[grep(":",rownames(cte1)),]
    colnames(cte1) <- sapply(str_split(colnames(cte1),"\\.|_"),"[",11)
    cte2 <- cte1[,meta_04$samp]
    #
    meta_04$cond <- as.factor(meta_04$cond)
    diff_01 <- DESeqDataSetFromMatrix(countData = cte2, 
                                      colData = meta_04, design = ~cond)
    diff_01 <- DESeq2::DESeq(diff_01,)
    diff_02 <- counts(diff_01, normalized = T)
    diff_03 <- diff_02[grep("MERVL-int|MT2_Mm",rownames(diff_02)),]
    diff_03 <- as.data.frame(diff_03)
    diff_03$name <- rownames(diff_03)
    diff_04 <- diff_03 %>% 
      tidyr::pivot_longer(cols = 1:c(length(colnames(diff_03))-1), names_to = "samp", values_to = "norm")
    diff_04 <- left_join(diff_04, meta_04[,c(3,2)], by = "samp")
    diff_05 <- as.data.frame(diff_04)
    ggplot(data = diff_05) +
      geom_boxplot(aes(x = cond, y = norm, fill = cond)) +
      scale_fill_few() +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "", y = "normlized reads count",fill ="") +
      theme_classic() +
      theme(plot.margin = margin(unit = "cm", 1,1,1,1),
            legend.text = element_text(family= "serif", size = 14),
            axis.text = element_text(family= "serif", size = 14),
            axis.title = element_text(family= "serif", size = 16),
            strip.text = element_text(family= "serif", size = 16),
            strip.background = element_rect(fill = "gray85"))
    #
    #
    #
    for (j in c("dBtgh","none")) {
      print(j)
      vsvs = c("cond","Btgh",j)
      meta_05 <- meta_04[which(meta_04$cond %in% c("Btgh",j)),]
      meta_05$cond <- as.factor(meta_05$cond)
      
      cte3 <- cte2[,rownames(meta_05)]
      cte3 <- cte3 %>% 
        mutate(across(where(is.character) | where(is.factor), readr::parse_number))
      #
      diff_01 <- DESeqDataSetFromMatrix(countData = cte3, 
                                        colData = meta_05, design = ~cond)
      diff_01 <- diff_01[rowSums(counts(diff_01)) >= 10,]
      ### diff_01 <- diff_01[rowMeans(BiocGenerics::counts(diff_01)) > 5,]
      diff_01 <- DESeq2::DESeq(diff_01,)
      diff_01 <- DESeq2::results(diff_01, contrast = vsvs)
      diff_02 <- as.data.frame(diff_01)
      diff_02$gene_id <- rownames(diff_02)
      diff_02 <- left_join(diff_02, gtf2[,c(1,3)], by = "gene_id")
      diff_02 <- diff_02[,c(7,8,1:6)]
      write.csv(diff_02, 
                paste0(path,"ourLab/DESeq2/DESeq2_morula_","Btgh-vs-",j, ".teOnly.rmBatch.csv"),
                quote = F, row.names = F)
      upup <- subset(diff_02, log2FoldChange > log2(1.5) & padj < 0.05)
      down <- subset(diff_02, log2FoldChange < -log2(1.5)& padj < 0.05)
      write.csv(upup, 
                paste0(path,"ourLab/DESeq2/DESeq2_morula_","Btgh-vs-",j, ".teOnly.upup.rmBatch.csv"),
                quote = F, row.names = F)
      write.csv(down, 
                paste0(path,"ourLab/DESeq2/DESeq2_morula_","Btgh-vs-",j, ".teOnly.down.rmBatch.csv"),
                quote = F, row.names = F)
    }
  }
}









################################################################################
#10、其他临时想法
{
  tmp1 <- read.csv("/disk5/aaSHY/ghyuu/publicData/ppnn/suppl/counts_table_all_stages.csv",
                   row.names = 1, header = T)
  tmp2 <- tmp1[,meta_01$samp]
  boxplot(as.numeric(tmp2["ENSMUSG00000034160",mmmm$Btgh_injected]),
          as.numeric(tmp2["ENSMUSG00000034160",mmmm$dBtgh_injected]),
          as.numeric(tmp2["ENSMUSG00000034160",mmmm$none]))
  
  grep("ENSMUSG00000034160", rownames(tmp2))
  tttt <- ruvs@assayData$normalizedCounts
  tttt <- tttt[9572,meta_01$samp]
  mmmm <- split(x = meta_01$samp, f = meta_01$cond)
  boxplot(tttt[mmmm$Btgh_injected],tttt[mmmm$dBtgh_injected],tttt[mmmm$none])
  
  gtf3[grep("^ghyvv|^ghyuu",ignore.case = T, x = gtf3$gene_nm),]
  mean(as.numeric(tmp2["ENSMUSG00000034160",meta_01$samp])) #ghyvv
  mean(as.numeric(tmp2["ENSMUSG00000025220",meta_01$samp])) #ghyuu
  
  
}
#####用CPM
#用meta_04 (在上面05的第一步)
{
  cte1 <- read.table(paste0(path,"ourLab/count/twos.cntTable"), 
                     header = T, sep = "\t", row.names = 1)
  #cte1 <- cte1[grep(":",rownames(cte1)),]
  colnames(cte1) <- sapply(str_split(colnames(cte1),"\\.|_"),"[",11)
  cte2 <- cte1[,meta_04$samp]
  cte3 <- as.data.frame(apply(cte2, 2, function(x){x/sum(x)*1000000}))
  
  #
  diff_03 <- cte3[grep("MERVL-int|MT2_Mm",rownames(cte3)),]
  diff_03 <- as.data.frame(diff_03)
  diff_03$name <- rownames(diff_03)
  diff_04 <- diff_03 %>% 
    tidyr::pivot_longer(cols = 1:c(length(colnames(diff_03))-1), names_to = "samp", values_to = "norm")
  diff_04 <- left_join(diff_04, meta_04[,c(1,2)], by = "samp")
  diff_05 <- as.data.frame(diff_04)
  ggplot(data = diff_05) +
    geom_boxplot(aes(x = cond, y = norm, fill = cond)) +
    scale_fill_few() +
    facet_wrap(~name, scales = "free_y") +
    labs(x = "", y = "normlized reads count",fill ="") +
    theme_classic() +
    theme(plot.margin = margin(unit = "cm", 1,1,1,1),
          legend.text = element_text(family= "serif", size = 14),
          axis.text   = element_text(family= "serif", size = 14),
          axis.title  = element_text(family= "serif", size = 16),
          strip.text  = element_text(family= "serif", size = 16),
          strip.background = element_rect(fill = "gray85"))
}
################################################################################
#11、其他临时想法
#####bw Correlations
{
  shell <- paste0(path,"src/run_shy_cor.sh")
  cat("#!/bin/bash\n", file = shell)
  shymed <- c("pearson","spearman")
  cmd_01 <- paste0("multiBigwigSummary bins -p 16 -b ",
                   paste0(path,"ourLab/bw/twos_",meta_04$cond,
                          "_",meta_04$samp,".bw ", collapse = ""),"-o ",path,
                   "ourLab/bw/heatmapCor-twos-geno.npz","\n")
  cmd_02 <- paste0("multiBigwigSummary BED-file -p 16 ",
                   "--BED ",inx6," ",inx7," ",
                   "-b ",paste0(path,"ourLab/bw/twos_",meta_04$cond,"_",
                                meta_04$samp,".bw ", collapse = ""),"-o ",path,
                   "ourLab/bw/heatmapCor-twos-beds.npz","\n")
  for (i in cmd_01) {cat(i, file = shell, append = T)}
  for (i in cmd_02) {cat(i, file = shell, append = T)}
  for (P in shymed) {
    cmd_03 <- paste0("plotCorrelation -in ",path,
                     "ourLab/bw/heatmapCor-twos-geno.npz ",
                     "--labels ",paste0(meta_04$cond,"-",meta_04$samp,collapse = " ")," ",
                     "--whatToPlot heatmap --corMethod ",P," -o ",path,
                     "ourLab/bw/heatmapCor-twos-geno-",P,".pdf","\n")
    cmd_04 <- paste0("plotCorrelation -in ",path,
                     "ourLab/bw/heatmapCor-twos-beds.npz ",
                     "--labels ",paste0(meta_04$cond,"-",meta_04$samp,collapse = " ")," ",
                     "--whatToPlot heatmap --corMethod ",P," -o ",path,
                     "ourLab/bw/heatmapCor-twos-beds-",P,".pdf","\n")
    for (i in cmd_03) {cat(i, file = shell, append = T)}
    for (i in cmd_04) {cat(i, file = shell, append = T)}
  }
  print(paste0("nohup bash ",shell, " >",path,"Log/",basename(shell),".log"," 2>&1 &"))
}


































