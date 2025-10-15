#write by Sun Haiayng at 2023.12.13
#01、搭环境
{
  rm(list=ls())
  for (i in c("homologene","data.table","ggridges","ggrepel","KEGG.db","scales","stringr","Biostrings","tidyr","pheatmap","enrichplot","DESeq2","patchwork","ggpubr","ggsci","ggtranscript","rtracklayer","Mfuzz","BiocManager","dplyr","IRanges","GenomicRanges","ReactomePA","ggplot2","clusterProfiler","topGO","Rgraphviz","pathview","org.Mm.eg.db")){library(i, character.only = T,quietly = T)}; rm(i)
}
#02、给索引
{
  path = "/ChIP_seq_2/zx/jjuuzz/"
  inx1 = "/Reference/aaSHY/Ensembl99/GRCm38/fa/GRCm38.dna.primary_assembly.fa"
  inx2 = "/Reference/aaSHY/Ensembl99/GRCm38/INDEX/STARv279a"
  inx3 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/Mus_musculus.GRCm38.99.gtf"
  inx4 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf"
  inx5 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf.locInd"
}
#03、跑命令
prex <- c("jjuuzz-ChIP","jjuuzz-Input")
cmd_01 <- paste0("macs2 callpeak -f BAMPE -g mm -q 0.1 --keep-dup 1 -t ",
                 path,"bam/",prex[1],".bam -c ", path,"bam/",prex[2],".bam -n ",
                 "jjuuzz.macs2","-N --outdir ",
                 path,"Homer/", "\n")
cmd_02 <- paste0("findMotifsGenome.pl ", path, "Homer/byP/jjuuzz.macs2-N-p0.001_peaks.narrowPeak ",
                 "mm10 ", path, "Homer/byP -size 200 -len 8,10,12")
cmd_03 <- paste0("findMotifsGenome.pl ", path, "Homer/byQ/jjuuzz.macs2-N-q0.1_peaks.narrowPeak ",
                 "mm10 ", path, "Homer/byQ -size 200 -len 8,10,12")












