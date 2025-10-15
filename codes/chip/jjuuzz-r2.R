#write by Sun Haiayng at 2024.05.30
#01、搭环境
{
  rm(list=ls())
  for (i in c("ReactomePA","ggplot2","stringr","ChIPseeker",
              "DESeq2","clusterProfiler","topGO","Rgraphviz",
              "pathview","org.Mm.eg.db","TxDb.Mmusculus.UCSC.mm10.knownGene")) {library(i, character.only = T)};rm(i)
}
#02、给索引
{
  path = "/disk5/aaSHY/Oga/ChIPseq/jjuuzz/"
  inx1 = "/Reference/aaSHY/Ensembl99/GRCm38/INDEX/bowtie2/mm10"
  inx2 = "/Reference/aaSHY/Ensembl99/GRCm38/INDEX/STARv279a"
  inx3 = "/Reference/aaSHY/Ensembl99/GRCm38/fa/GRCm38.dna.primary_assembly.fa"
  inx4 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/Mus_musculus.GRCm38.99.gtf"
  inx5 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf"
  inx6 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensem0bl_rmsk_TE.gtf.locInd"
  inx7 = "/Reference/aaSHY/Ensembl99/GRCm38/BED/TEsubFamily/MM-int::ERVL::LTR.bed"
  inx8 = "/Reference/aaSHY/Ensembl99/GRCm38/BED/special/MM_zyw.bed"
  inx9 = "/Reference/aaSHY/Ensembl99/GRCm38/BED/TEsubFamily/MT2_Mm::ERVL::LTR.bed"
  inxA = "/Reference/aaSHY/Ensembl99/GRCm38/BED/special/allGene.bed"
  inxB = "/Reference/aaSHY/Ensembl99/GRCm38/BED/special/twoC.bed"
  inxC = "/Reference/aaSHY/zOther/Rn45s/INDEX/bowtie2/Rn45s"
  inxD = "/Reference/aaSHY/zOther/TEconsensus/INDEX/MuERVL/bt2/MuERVL"
}
################################################################################
#03、跑命令
{
  for (i in c("src","fq","Log","peak","dumpROOM","DESeq2","rRNA/csem","MM/csem","fisher","deeptools","bam/csem","trim","PLOT","bw/bcv","bw/bcp","count")) {
    if (dir.exists(paths = paste0(path, i))==F) {
      dir.create(path = paste0(path,i), recursive = T)
    }
  }
  shell <- paste0(path,"src/run_a.sh")
  cat("#!/bin/bash\n", file = shell)
  prex <- c("IP","Input")
  cmd_01 <- paste0("#parallel-fastq-dump -t 10 --split-3 --gzip -s ",
                   list.files(path = paste0(path, "rawdata"), pattern = "sra", recursive = T, full.names = T),
                   " -O ", path, "fq", "\n")
  cmd_02 <- paste0("#rename s/fastq/fq/ ", path, "fq/*gz", "\n")
  lay1 <- paste0(path,"fq/",prex,"_1.fq.gz")
  lay2 <- paste0(path,"fq/",prex,"_2.fq.gz")
  lay3 <- paste0(path,"trim/",prex,"_1_val_1.fq.gz")
  lay4 <- paste0(path,"trim/",prex,"_2_val_2.fq.gz")
  cmd_03 <- paste0("#trim_galore --phred33 -j 2 ",
                   "-a 'AGATCGGAAGAGC -a G{10}' -a2 'AGATCGGAAGAGC -a G{10}' --fastqc ",
                   "--clip_R1 8 --clip_R2 8 --three_prime_clip_R1 2 --three_prime_clip_R2 2 ",
                   "-o ", path, "trim --paired ",lay1," ",lay2, "\n")
  cmd_04 <- paste0("bowtie2 -q --sensitive -k 5 --no-mixed --no-discordant ",
                   "--reorder -p 20 -x ",inx1," -1 ",lay3," -2 ",lay4,
                   " |samtools view -F 4 -b > ", path,
                   "bam/",prex,".bam","\n")
  cmd_05 <- paste0("run-csem -p 20 --no-extending-reads --bam ",
                   paste0(path,"bam/",prex,".bam")," 140 ",
                   paste0(path,"bam/csem/",prex),"\n")
  cmd_06 <- paste0("bamCoverage --normalizeUsing CPM -p 20 --minMappingQuality 1 ",
                   "--samFlagExclude 256 -of bigwig -bs 20 -b ",
                   paste0(path,"bam/csem/",prex,".sorted.bam -o "),
                   paste0(path,"bw/bcv/",prex,".bw\n"))
  cmd_07 <- paste0("bamCompare -p 20 --scaleFactorsMethod SES ",
                   "--sampleLength 100000 -b1 ",
                   paste0(path,"bam/csem/",prex[1],".sorted.bam -b2 "),
                   paste0(path,"bam/csem/",prex[2],".sorted.bam -o "),
                   paste0(path,"bw/bcp/",prex[1],"_Input.bw"),"\n")
  cmd_08 <- paste0("computeMatrix scale-regions -p 20 -a 2000 -b 2000 -m 5000 ",
                   "--missingDataAsZero -R ", c(inx7,inx8,inx9,inxA,inxB),
                   " -S ", paste(path,"bw/bcp/",prex[1],"_Input.bw", sep = "", collapse = " "),
                   " -o ", path,"deeptools/", basename(c(inx7,inx8,inx9,inxA,inxB)), ".mat.gz","\n")
  cmd_09 <- paste0("plotHeatmap --heatmapHeight 12 --colorList white,red -m ",path,
                   "deeptools/",basename(c(inx7,inx8,inx9,inxA,inxB)),".mat.gz ",
                   "-out ",path,"PLOT/",
                   basename(c(inx7,inx8,inx9,inxA,inxB)),".pdf","\n")
  cmd_10 <- paste0("macs2 callpeak -f BAM -g mm -p 0.001 --keep-dup 1 -t ",
                   path,"bam/csem/",prex[1],".sorted.bam -c ",
                   path,"bam/csem/",prex[2],".sorted.bam -n ",
                   prex[1],"-N --outdir ",
                   path,"peak/", "\n")
  cmd_11 <- paste0("macs2 callpeak -f BAM -g mm -p 0.001 --broad --broad-cutoff 0.001 ",
                   "--keep-dup 1 -t ",path,"bam/csem/", prex[1],
                   ".sorted.bam -c ", path,"bam/csem/", prex[2],
                   ".sorted.bam -n ", prex[1],
                   "-B --outdir ", path,"peak/", "\n")
  cmd_12 <- paste0("#TEtranscripts -t ",
                   paste0(path,"bam/csem/",prex[1],".sorted.bam", collapse = " ")," -c ",
                   paste0(path,"bam/csem/",prex[2],".sorted.bam", collapse = " "),
                   " --GTF ",inx4," --TE ",inx5," --sortByPos --mode multi ",
                   "--project dot1 --outdir ",paste0(path,"count/"),"\n")
  cmd_13 <- paste0("#TElocal -b ",
                   path,"bam/csem/",prex,".sorted.bam ",
                   "--GTF ",inx4," --TE ",inx6," --sortByPos --project ",
                   path, "count/", prex,"\n")
  for (i in cmd_01) {cat(i, append = T, file = shell)}
  for (i in cmd_02) {cat(i, append = T, file = shell)}
  for (i in cmd_03) {cat(i, append = T, file = shell)}
  for (i in cmd_04) {cat(i, append = T, file = shell)}
  for (i in cmd_05) {cat(i, append = T, file = shell)}
  for (i in cmd_06) {cat(i, append = T, file = shell)}
  for (i in cmd_07) {cat(i, append = T, file = shell)}
  for (i in cmd_08) {cat(i, append = T, file = shell)}
  for (i in cmd_09) {cat(i, append = T, file = shell)}
  for (i in cmd_10) {cat(i, append = T, file = shell)}
  for (i in cmd_11) {cat(i, append = T, file = shell)}
  for (i in cmd_12) {cat(i, append = T, file = shell)}
  for (i in cmd_13) {cat(i, append = T, file = shell)}
  for (kk in c("MM","rRNA")) {
    if (kk == "MM") {
      inxM <- inxD
    }
    else {
      inxM <- inxC
    }
    cmd_14 <- paste0("#bowtie2 -q --sensitive -k 5 --no-mixed --no-discordant ",
                     "--reorder -p 20 -x ",inxM," -1 ", lay3," -2 ",lay4," ",
                     "|samtools view -F 4 -b > ", path, kk, "/",prex,".bam","\n")
    cmd_15 <- paste0("#run-csem -p 20 --no-extending-reads --bam ",
                     paste0(path,kk,"/",prex,".bam")," 140 ",
                     paste0(path,kk,"/csem/",prex),"\n")
    cmd_16 <- paste0("#bamCoverage --normalizeUsing CPM -p 20 --minMappingQuality 1 ",
                     "--samFlagExclude 256 -of bigwig -bs 20 -b ",
                     paste0(path,kk,"/csem/",prex,".sorted.bam -o "),
                     paste0(path,kk,"/",prex,".",kk,".bw\n"))
    for (i in cmd_14) {cat(i, append = T, file = shell)}
    for (i in cmd_15) {cat(i, append = T, file = shell)}
    for (i in cmd_16) {cat(i, append = T, file = shell)}
  }
  print(paste0("nohup bash ",shell," >",path,"Log/run_a.log"," 2>&1 &"))
}
################################################################################
#04、peak Analysis (pvalue 0.05)
{
  prex <- c("IP","Input")
  gtf1 <- rtracklayer::import(con = inx4, format = "gtf")
  gtf1 <- as.data.frame(gtf1)
  colnames(gtf1)[c(10,12)] <- c("gene_name","gene_id")
  gtf1 <- GenomicRanges::makeGRangesFromDataFrame(df = gtf1, keep.extra.columns = T)
  txb1 <- suppressWarnings(GenomicFeatures::makeTxDbFromGRanges(gr = gtf1, taxonomyId = 10090))
  twoC <- read.table(inxB)
  gtf2 <- gtf1[which(gtf1$gene_id %in% unique(twoC$V4))]
  txb2 <- suppressWarnings(GenomicFeatures::makeTxDbFromGRanges(gr = gtf2, taxonomyId = 10090))
  #
  for (i in (prex[1])) {
    pk01 <- suppressWarnings(ChIPseeker::annotatePeak(peak = paste0(path, "peak/",i,"-N_peaks.narrowPeak"), 
                                                      tssRegion = c(-2000, 2000), TxDb = txb1))
    pdf(file = paste0(path, "PLOT/",i,"-N_chipseekerAnno.pdf"), width = 8, height = 8)
    ChIPseeker::plotAnnoPie(pk01, legend.position = "rightside")
    dev.off()
    write.csv(as.data.frame(pk01), paste0(path, "peak/",i,"-N_chipseekerAnno.csv"), quote = F, row.names = F)
  }
}
################################################################################
# 补充为center画法
{
  cmd_01 <- paste0("computeMatrix reference-point --referencePoint center ",
                   "-p 20 -a 3000 -b 3000 ", 
                   "--skipZeros --missingDataAsZero -bs 10 ",
                   "-R ",inx8, " -S ",paste0(path,"bw/bcp/",prex[1],"_Input.bw")," ",
                   "-o ",path,"deeptools/", basename(c(inx8)), ".mat.gz", "\n")
  cmd_02 <- paste0("plotHeatmap --heatmapHeight 12 --colorList white,red -m ",path,
                   "deeptools/",basename(c(inx8)),".mat.gz ",
                   "-out ",path,"PLOT/",
                   basename(c(inx8)),"-center.pdf","\n")
  print(cmd_01)
  print(cmd_02)
}























































