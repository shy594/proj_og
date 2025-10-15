#write by Sun Haiayng at 2025.01.17
#01、搭环境
{
  rm(list=ls())
  for (i in c("homologene","data.table","ggridges","ggrepel","KEGG.db","scales","stringr","Biostrings","tidyr","pheatmap","enrichplot","DESeq2","patchwork","ggpubr","ggsci","ggtranscript","rtracklayer","Mfuzz","BiocManager","dplyr","IRanges","GenomicRanges","ReactomePA","ggplot2","clusterProfiler","topGO","Rgraphviz","pathview","org.Mm.eg.db")){library(i, character.only = T,quietly = T)}; rm(i)
}
#02、给索引
{
  path = "/sc/aaSHY/ribo/riboITP/rnaSeq/"
  inx1 = "/Reference/aaSHY/Ensembl99/GRCm38/fa/GRCm38.dna.primary_assembly.fa"
  inx2 = "/Reference/aaSHY/Ensembl99/GRCm38/INDEX/STARv279a"
  inx3 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/Mus_musculus.GRCm38.99.gtf"
  inx4 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf"
  inx5 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf.locInd"
  inx6 = "/Reference/aaSHY/zOther/INDEX/Bowtie2/GRCm38rRNA/rRNA"
  inx7 = "/Reference/aaSHY/Ensembl99/GRCm38/INDEX/bowtie2/mm10"
  inx8 = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/GenBank/rRNA.fa"
}
################################################################################
#03、跑命令
{
  for (i in c("src","fq","Log","stringtie","umi_tools",
              "rRNAdust","bamCoverage","trim","dumpROOM","bam","count")) {
    if (dir.exists(paths = paste0(path, i))==F) {
      dir.create(path = paste0(path,i), recursive = T)
    }
  }
  shell <- paste0(path,"src/run_a.sh")
  cat("#!/bin/bash\n", file = shell)
  prex <- list.files(path = paste0(path, "rawdata"), full.names = T, recursive = T, pattern = "sra")
  prex <- gsub(x = basename(prex), pattern = ".sra", replacement = "", fixed = T)
  cmd_01 <- paste0("#parallel-fastq-dump -t 10 --split-3 --gzip -s ",
                   path, "rawdata/", prex, "/", prex, ".sra",
                   " -O ", path, "fq", "\n")
  cmd_02 <- paste0("#cutadapt -j 2 --trimmed-only -m 8 -g ATTGCGCAATG -o ",path,
                   "trim/", prex, "_1.fq.gz -p ",path,"trim/",prex,
                   "_2.fq.gz ", path, "fq/", prex, 
                   "_1.fq.gz ",path,"fq/",prex,"_2.fq.gz","\n")
  cmd_03 <- paste0("#umi_tools extract -p 'NNNNNNNN' ",
                   "-I ", path, "trim/", prex, "_1.fq.gz ",
                   "-S ", path, "umi_tools/", prex,"_1.fq.gz ", "--read2-in=",path,
                   "trim/",prex,"_2.fq.gz ", "--read2-out=",path,
                   "umi_tools/",prex,"_2.fq ","--log=",path,"Log/cmd_03.log", "\n")
  cmd_03 <- paste0("#rRNAdust -t 4 ", inx8, " ", path, "umi_tools/", prex, "_2.fq ", 
                   "-o ", path, "rRNAdust/", prex, "_2.fq", "\n")
  lays <- paste0(path, "rRNAdust/", prex, "_2.fq")
  cmd_04 <- paste0("#STAR --runThreadN 10 --readFilesCommand zcat ",
                   "--outMultimapperOrder Random --outSAMtype BAM SortedByCoordinate ",
                   "--outFileNamePrefix ",path,"bam/",prex," ",
                   "--outSAMprimaryFlag AllBestScore ",
                   "--genomeDir ",inx2," --readFilesIn ",lays,"\n")
  cmd_05 <- paste0("#mv ",path,"bam/",prex,"Aligned.sortedByCoord.out.bam ",path,
                   "bam/",prex,".bam","\n")
  cmd_06 <- paste0("umi_tools dedup -I ",path,"bam/",prex,".bam -S ",path,
                   "bam/",prex,".dedup.bam","\n")
  cmd_07 <- paste0("samtools index ",path,"bam/",prex,".dedup.bam","\n")
  cmd_08 <- paste0("TEtranscripts -t ",
                   paste(paste0(path,"bam/",prex[01:11],".dedup.bam"),collapse = " ")," -c ",
                   paste(paste0(path,"bam/",prex[12:22],".dedup.bam"),collapse = " ")," ",
                   "--GTF ",inx3," --TE ",inx4," --sortByPos --mode multi ",
                   "--project itp_dedup --outdir ",paste0(path,"count"),"\n")
  for (i in cmd_01) {cat(i, append = T, file = shell)}
  for (i in cmd_02) {cat(i, append = T, file = shell)}
  for (i in cmd_03) {cat(i, append = T, file = shell)}
  for (i in cmd_04) {cat(i, append = T, file = shell)}
  for (i in cmd_05) {cat(i, append = T, file = shell)}
  for (i in cmd_06) {cat(i, append = T, file = shell)}
  for (i in cmd_07) {cat(i, append = T, file = shell)}
  for (i in cmd_08) {cat(i, append = T, file = shell)}
  aaaa <- unique(sapply(str_split(prex,"-"),"[",1))
  for (i in aaaa) {
    brex <- prex[grep(i, prex)]
    cmd_09 <- paste0("samtools merge -@ 10 -f -o ", path,
                     "bam/", i, ".dedup.bam ",
                     paste0(path, "bam/", brex, ".dedup.bam", collapse = " "), "\n")
    cmd_10 <- paste0("samtools index ", path, "bam/", i, ".dedup.bam ", "\n")
    cmd_11 <- paste0("bamCoverage --normalizeUsing CPM -p 20 --minMappingQuality 1 ",
                     "--samFlagExclude 256 -of bigwig -bs 20 -b ",path,
                     "bam/", i, ".dedup.bam -o ",path,
                     "bamCoverage/",i,".dedup.bw\n")
    for (i in cmd_09) {cat(i, append = T, file = shell)}
    for (i in cmd_10) {cat(i, append = T, file = shell)}
    for (i in cmd_11) {cat(i, append = T, file = shell)}
  }
  print(paste0("nohup bash ",shell," >",path,"Log/run_a.log"," 2>&1 &"))
}
################################################################################
#04、DESeq2：Gene、TE、TE sites
{
  shy_diff <- function(mark, geneFile, VS){
    tmp1 <- read.table(geneFile, header = T, check.names = F, stringsAsFactors = F, row.names = 1)
    colnames(tmp1) <- sapply(str_split(basename(colnames(tmp1)),pattern = ".dedup.bam"), "[", 1)
    tmp1 <- tmp1[,grep(paste0(VS[2],"|",VS[3]), colnames(tmp1))]
    Class <- c(colnames(tmp1)[grep(VS[2],colnames(tmp1))],
               colnames(tmp1)[grep(VS[3],colnames(tmp1))])
    Class <- factor(sapply(str_split(Class,"-"),"[",1))
    dfclas <- data.frame(row.names = colnames(tmp1), Class, stringsAsFactors = T)
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
  }
  shy_diff(mark = "twos_zygote",
           geneFile = paste0(path, "count/itp_dedup.cntTable"),
           VS = c("Class","twos","zygote"))
}
