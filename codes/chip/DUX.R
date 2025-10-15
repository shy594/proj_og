#write by Sun Haiayng at 2024.02.01
#01、搭环境
{
  rm(list=ls())
  for (i in c("ReactomePA","ggplot2","stringr","ChIPseeker","DESeq2",
              "clusterProfiler","topGO","Rgraphviz","pathview","org.Mm.eg.db","TxDb.Mmusculus.UCSC.mm10.knownGene")) {library(i, character.only = T)};rm(i)
}
#02、给索引
{
  path = "/ChIP_seq_2/aaSHY/guuyvvMM/publicData/Dux__ChIP/"
  inx1 = "/Reference/aaSHY/Ensembl99/GRCm38/INDEX/bowtie2/mm10"
  inx2 = "/Reference/aaSHY/Ensembl99/GRCm38/INDEX/STARv279a"
  inx3 = "/Reference/aaSHY/Ensembl99/GRCm38/fa/GRCm38.dna.primary_assembly.fa"
  inx4 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/Mus_musculus.GRCm38.99.gtf"
  inx5 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf"
  inx6 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf.locInd"
  inx7 = "/Reference/aaSHY/Ensembl99/GRCm38/BED/TEsubFamily/guuyvvMM-int::ERVL::LTR.bed"
  inx8 = "/Reference/aaSHY/Ensembl99/GRCm38/BED/special/guuyvvMM_zyw.bed"
  inx9 = "/Reference/aaSHY/Ensembl99/GRCm38/BED/TEsubFamily/guuyvvMM::ERVL::LTR.bed"
  inxA = "/Reference/aaSHY/Ensembl99/GRCm38/BED/special/allGene.bed"
  inxB = "/Reference/aaSHY/Ensembl99/GRCm38/BED/special/twoC.bed"
  inxC = "/Reference/aaSHY/zOther/Rn45s/INDEX/bowtie2/Rn45s"
  inxD = "/Reference/aaSHY/zOther/TEconsensus/MuERVL/bt2/MuERVL"
}
################################################################################
#03、跑命令
{
  for (i in c("src","fq","Log","dumpROOM","bam/csem","trim","PLOT","bw/bcv","bw/bcp")) {
    if (dir.exists(paths = paste0(path, i))==F) {
      dir.create(path = paste0(path,i), recursive = T)
    }
  }
  shell <- paste0(path,"run_a.sh")
  cat("#!/bin/bash\n", file = shell)
  prex <- list.files(path = paste0(path, "rawdata"), pattern = "sra", recursive = T)
  prex <- unique(gsub(pattern = ".sra", replacement = "", fixed = F, x = prex))
  cmd_01 <- paste0("#parallel-fastq-dump -t 10 --split-3 --gzip -s ",path, 
                   "rawdata/",prex,".sra -O ", path, "fq", "\n")
  cmd_02 <- paste0("#rename s/fastq/fq/ ", path, "fq/*gz", "\n")
  lay1 <- paste0(path,"fq/",prex,"_1.fq.gz")
  lay2 <- paste0(path,"fq/",prex,"_2.fq.gz")
  cmd_03 <- paste0("#trim_galore --phred33 --illumina ",
                   "--paired ",lay1," ",lay2," -o ",path,"trim","\n")
  lay1 <- paste0(path,"trim/",prex,"_1_val_1.fq.gz")
  lay2 <- paste0(path,"trim/",prex,"_2_val_2.fq.gz")
  cmd_04 <- paste0("#bowtie2 -p 20 -q --sensitive -k 5 --reorder ",
                   "--no-mixed --no-discordant ",
                   "-x ",inx1," -1 ", lay1, " -2 ",lay2," ",
                   "|samtools view -F 4 -b > ", path, "bam/",prex,".bam","\n")
  cmd_05 <- paste0("#run-csem -p 20 --no-extending-reads --bam ",path,
                   "bam/",prex,".bam"," 120 ",path,
                   "bam/csem/",prex,"\n")
  cmd_06 <- paste0("bamCoverage --normalizeUsing CPM -p 20 --minMappingQuality 1 ",
                   "--samFlagExclude 256 -of bigwig -bs 20 -b ",
                   paste0(path,"bam/csem/",prex,".sorted.bam -o "),
                   paste0(path,"bw/bcv/",prex,".bw\n"))
  cmd_07 <- paste0("bamCompare -p 20 --scaleFactorsMethod SES ",
                   "--sampleLength 100000 -b1 ",
                   path,"bam/csem/",prex[2:5],".sorted.bam -b2 ",
                   path,"bam/csem/",prex[1],".sorted.bam -o ",
                   path,"bw/bcp/",prex[2:5],"_Input.bw","\n")
  for (i in cmd_01) {cat(i, append = T, file = shell)}
  for (i in cmd_02) {cat(i, append = T, file = shell)}
  for (i in cmd_03) {cat(i, append = T, file = shell)}
  for (i in cmd_04) {cat(i, append = T, file = shell)}
  for (i in cmd_05) {cat(i, append = T, file = shell)}
  for (i in cmd_06) {cat(i, append = T, file = shell)}
  for (i in cmd_07) {cat(i, append = T, file = shell)}
  print(paste0("nohup bash ",shell," >",path,"Log/run_a.log"," 2>&1 &"))
}
###~~~ 补充bamCompare
{

  
}

################################################################################
####文章里的参数
cmd_04 <- paste0("bowtie2 -p 20 ",
                 "-t -q -N 1 -L 25 -X 2000 --no-mixed --no-discordant ",
                 "-x ",inx1," -1 ", lay1, " -2 ",lay2," ",
                 "|samtools view -F 4 -b > ", path, "bam/",prex,".bam","\n")


