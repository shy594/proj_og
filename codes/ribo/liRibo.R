#write by Sun Haiayng at 2024.01.04
#01、搭环境
{
  rm(list=ls())
  for (i in c("homologene","data.table","ggridges","ggrepel","KEGG.db","scales","stringr","Biostrings","tidyr","pheatmap","enrichplot","DESeq2","patchwork","ggpubr","ggsci","ggtranscript","rtracklayer","Mfuzz","BiocManager","dplyr","IRanges","GenomicRanges","ReactomePA","ggplot2","clusterProfiler","topGO","Rgraphviz","pathview","org.Mm.eg.db")){library(i, character.only = T,quietly = T)}; rm(i)
}
#02、给索引
{
  path = "/ChIP_seq_1/aaSHY/ribo/liRibo/ribo/"
  inx1 = "/Reference/aaSHY/Ensembl99/GRCm38/fa/GRCm38.dna.primary_assembly.fa"
  inx2 = "/Reference/aaSHY/Ensembl99/GRCm38/INDEX/STARv279a"
  inx3 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/Mus_musculus.GRCm38.99.gtf"
  inx4 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf"
  inx5 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf.locInd"
  inx6 = "/ChIP_seq_1/aaSHY/ribo/liRibo/ribo/trRNAbt2Index/mm10"
  inx7 = "/Reference/aaSHY/Ensembl99/GRCm38/INDEX/bowtie2/mm10"
}
#03、跑命令
{
  for (i in c("src","fq","Log","stringtie","bamCoverage","trim","dumpROOM","bam","count")) {
    if (dir.exists(paths = paste0(path, i))==F) {
      dir.create(path = paste0(path,i), recursive = T)
    }
  }
  shell <- paste0(path,"src/run_a.sh")
  cat("#!/bin/bash\n", file = shell)
  prex <- list.files(path = paste0(path, "rawdata"), full.names = T, recursive = T, pattern = "SRR*")
  prex <- basename(prex)
  cmd_01 <- paste0("#parallel-fastq-dump -t 10 --split-3 --gzip -s ",
                   path, "rawdata/", prex, " -O ", path, "fq", "\n")
  cmd_02 <- paste0("#rename s/fastq/fq/ ", path, "fq/*gz", "\n")
  cmd_03 <- paste0("#cutadapt -m 20 -u 3 -a AAAAAAAA -j 10 -o ",
                   path, "trim/", prex, ".fq.gz ", path, "fq/", prex, ".fq.gz", "\n")
  cmd_04 <- paste0("#bowtie2 -x ", inx6, " -p 20 -U ", path, "trim/", prex, ".fq.gz ",
                   "|samtools view -h -f 4 -b >", path, "bam/", prex, "_v1.bam", "\n")
  cmd_05 <- paste0("STAR --genomeDir ",inx2, " --genomeLoad LoadAndExit","\n")
  cmd_06 <- paste0("STAR --runThreadN 10 --readFilesCommand samtools view ",
                   "--genomeLoad LoadAndKeep --limitBAMsortRAM 30000000000 ",
                   "--alignEndsType EndToEnd --readFilesType SAM SE ",
                   "--outMultimapperOrder Random --outSAMtype BAM SortedByCoordinate ",
                   "--outFilterMismatchNmax 2 --outFileNamePrefix ",path,"bam/",prex,"_v2 ",
                   "--outSAMprimaryFlag AllBestScore --outSAMmultNmax -1 ",
                   "--genomeDir ",inx2," --readFilesIn ", path, "bam/", prex, "_v1.bam","\n")
  cmd_07 <- paste0("mv ",path,"bam/",prex,"_v2Aligned.sortedByCoord.out.bam ",path,"bam/",prex,"_v2.bam","\n")
  cmd_08 <- paste0("samtools index ",path,"bam/",prex,"_v2.bam","\n")
  for (i in cmd_01) {cat(i, append = T, file = shell)}
  for (i in cmd_02) {cat(i, append = T, file = shell)}
  for (i in cmd_03) {cat(i, append = T, file = shell)}
  for (i in cmd_04) {cat(i, append = T, file = shell)}
  for (i in cmd_05) {cat(i, append = T, file = shell)}
  for (i in cmd_06) {cat(i, append = T, file = shell)}
  for (i in cmd_07) {cat(i, append = T, file = shell)}
  for (i in cmd_08) {cat(i, append = T, file = shell)}
  ff01 <- read.table(file = paste0(path, "rawdata/SraRunTable.txt"), sep = "\t", header = T, stringsAsFactors = F)
  for (i in unique(ff01$source_name)) {
    prex <- ff01$Run[which(ff01$source_name==i)]
    cmd_09 <- paste0("samtools merge -@ 10 -f -o ", path, "bam/", i, ".sorted.bam ",
                     paste0(path, "bam/", prex, "_v2.bam", collapse = " "), "\n")
    cmd_10 <- paste0("samtools index ", path, "bam/", i, ".sorted.bam ", "\n")
    cmd_11 <- paste0("bamCoverage --normalizeUsing CPM -p 20 --minMappingQuality 1 ",
                     "--samFlagExclude 256 -of bigwig -bs 20 -b ",path, 
                     "bam/", i, ".sorted.bam -o ",path,
                     "bamCoverage/",i,".bw\n")
    for (i in cmd_09) {cat(i, append = T, file = shell)}
    for (i in cmd_10) {cat(i, append = T, file = shell)}
    for (i in cmd_11) {cat(i, append = T, file = shell)}
  }
  print(paste0("nohup bash ",shell," >",path,"Log/run_a.log"," 2>&1 &"))
}






