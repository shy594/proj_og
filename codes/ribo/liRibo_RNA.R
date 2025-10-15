#write by Sun Haiayng at 2024.01.28
#01、搭环境
{
  rm(list=ls())
  for (i in c("homologene","data.table","ggridges","ggrepel","KEGG.db","scales","stringr","Biostrings","tidyr","pheatmap","enrichplot","DESeq2","patchwork","ggpubr","ggsci","ggtranscript","rtracklayer","Mfuzz","BiocManager","dplyr","IRanges","GenomicRanges","ReactomePA","ggplot2","clusterProfiler","topGO","Rgraphviz","pathview","org.Mm.eg.db")){library(i, character.only = T,quietly = T)}; rm(i)
}
#02、给索引
{
  path = "/sc/aaSHY/ribo/liRibo/totalRNA/"
  inx1 = "/Reference/aaSHY/Ensembl99/GRCm38/fa/GRCm38.dna.primary_assembly.fa"
  inx2 = "/Reference/aaSHY/Ensembl99/GRCm38/INDEX/STARv279a"
  inx3 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/Mus_musculus.GRCm38.99.gtf"
  inx4 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf"
  inx5 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf.locInd"
}
################################################################################
#03、跑命令
{
  for (i in c("src","fq","Log","stringtie","bamCoverage","trim","dumpROOM","bam","count")) {
    if (dir.exists(paths = paste0(path, i))==F) {
      dir.create(path = paste0(path,i), recursive = T)
    }
  }
  shell <- paste0(path,"src/run_a.sh")
  cat("#!/bin/bash\n", file = shell)
  prex <- list.files(paste0(path, "rawdata"), recursive = T, full.names = T, pattern = "sra$")
  prex <- unique(sapply(str_split(basename(prex), pattern = "-"), "[", 1))
  cmd_01 <- paste0("#parallel-fastq-dump -t 10 --split-3 --gzip -s ",
                   list.files(paste0(path, "rawdata"), recursive = T, full.names = T, pattern = "sra$"),
                   " -O ", path, "fq", "\n")
  cmd_02 <- paste0("#rename s/fastq/fq/ ", path, "fq/*gz", "\n")
  arex <- paste0(prex, "-",c(rep(1,6), rep(2,6)))
  lay1 <- paste0(path, "fq/", arex, rep("_1",12), ".fq.gz")
  lay2 <- paste0(path, "fq/", arex, rep("_2",12), ".fq.gz")
  cmd_03 <- paste0("#trim_galore --phred33 -j 8 --paired -o ", path, "trim ", 
                   lay1, " ", lay2, "\n")
  lay1 <- paste0(path, "trim/", arex, rep("_1_val_1",12), ".fq.gz")
  lay2 <- paste0(path, "trim/", arex, rep("_2_val_2",12), ".fq.gz")
  cmd_04 <- paste0("STAR --runThreadN 20 --genomeDir ", inx2," ",
                   "--readFilesCommand zcat --outSAMtype BAM SortedByCoordinate ",
                   "--outSAMstrandField intronMotif --outFilterType BySJout ",
                   "--outSAMprimaryFlag AllBestScore --twopassMode Basic ",
                   "--outSAMattributes All --outFileNamePrefix ",path,
                   "bam/",arex," --readFilesIn ",lay1," ",lay2,"\n")
  cmd_05 <- paste0("mv ",path,"bam/",arex,
                   "Aligned.sortedByCoord.out.bam ",path,"bam/",arex,".bam","\n")
  cmd_06 <- paste0("samtools index ",path,"bam/",arex,".bam","\n")
  cmd_07 <- paste0("#TEtranscripts -t ",
                   paste(paste0(path,"bam/",arex[1:6],".bam"),collapse = " ")," -c ",
                   paste(paste0(path,"bam/",arex[7:12],".bam"),collapse = " "),
                   " --GTF ",inx3," --TE ",inx4," --sortByPos --mode multi ",
                   "--project stage",i," --outdir ",path,"count","\n")
  cmd_08 <- paste0("#TElocal -b ",
                   path,"bam/",arex,".bam",
                   " --GTF ",inx3," --TE ",inx5," --sortByPos --project ",
                   paste0(path, "count/",arex),"\n")
  cmd_09 <- paste0("#paste ", paste0(path, "count/",
                                    arex[c(1,7,2,8,3,9,4,10,5,11,6,12)], ".cntTable", 
                                    collapse = " "), " |",
                   "cut -f ", paste0(c(1,seq(2,length(arex)*2, by=2)), collapse = ","), " >", path, 
                   "count/stageALL.cntTable","\n")
  cmd_10 <- paste0("#bamCoverage --normalizeUsing CPM -p 20 ",
                   "--minMappingQuality 1 --samFlagExclude 256 -of bigwig -bs 20 -b ",
                   path,"bam/",arex,".bam -o ",path,"bamCoverage/",arex,".bw","\n")
  cmd_11 <- paste0("stringtie ", path,"bam/",arex,".bam ",
                   "-p 20 -G ", inx3, " -o ", path,"stringtie/",arex,".gtf","\n")
  cmd_12 <- paste0("stringtie --merge -G ",
                   inx3," -o ",path,"stringtie/stage.merge.gtf ",
                   paste0(path,"stringtie/",arex,".gtf", collapse = " "),"\n")
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
  for (i in seq(6)) {
    cmd_13 <- paste0("#stringtie --merge -G ", inx3, " -o ", path, "stringtie/", prex[i], 
                     ".gtf ", paste0(path, "stringtie/",c(arex[i],arex[i+6]),".gtf", collapse = " "), "\n")
    for (i in cmd_13) {cat(i, append = T, file = shell)}
  }; rm(i)
  print(paste0("nohup bash ",shell," >",path,"Log/run_a.log"," 2>&1 &"))
}
################################################################################
#04、表达量
{
  cn01 <- read.table(file = paste0(path, "count/stageALL.cntTable"), sep = "\t", header = T, row.names = 1, stringsAsFactors = F)
  cn01 <- cn01[,colnames(cn01)[grep(x = colnames(cn01), pattern = "gene", invert = T)]]
  prex <- unique(sapply(str_split(colnames(cn01), pattern = "\\."), "[", 7))
  prex <- prex[c(11,10,2,3,8,5,12,1,7,4,9,6)]
  colnames(cn01) <- paste(sapply(str_split(colnames(cn01), pattern = "\\."), "[", 7),sapply(str_split(colnames(cn01), pattern = "\\."), "[", 8), sep = "-")
  cn02 <- as.data.frame(apply(cn01, 2, function(x){x/sum(x)*1000000}))
  cn02$geneName <- rownames(cn02)
  cn03 <- cn02[grep(rownames(cn02), pattern=":", invert = T),]#invert = F, 那么输出TE的CPM
  df01 <- tidyr::gather(data = cn03, key = "sample", value = "cpm", -"geneName")
  df01$stage <- sapply(str_split(df01$sample, "-"), "[", 1)
  df02 <- dplyr::reframe(.data = df01, .by = c(stage, geneName), cpmExp = mean(cpm), cpmSD = sd(cpm)/sqrt(length(cpm)))
  wr01 <- df02
  wr01$cpmAndSD <- paste0(round(wr01$cpmExp, 3), " (", round(wr01$cpmSD, 3), ")")
  wr01 <- tidyr::spread(data = wr01[,-c(3,4)], key = "stage", value = "cpmAndSD")
  openxlsx::write.xlsx(x = wr01, file = paste0(path, "count/gene.cpmAndSD.xlsx"))
}
################################################################################
#05、可视化
{
  df03 <- df02[which(df02$geneName %in% c("Fcc28","Rau55")),]
  df03 <- df02[which(df02$geneName %in% c("Hox4","Dot1")),]
  df03 <- df02[which(df02$geneName %in% unique(ms01$Gene)[17:32]),]
  #
  df03$stage <- factor(x = df03$stage, levels = prex)
  p_01 <- ggplot(data = df03) +
    geom_line(mapping = aes(x = stage, y = cpmExp, group = geneName), color = "steelblue", linewidth = 1) +
    geom_point(mapping = aes(x = stage, y = cpmExp), color = "steelblue", size = 1.5) +
    geom_errorbar(mapping = aes(x = stage, ymin = cpmExp-cpmSD, ymax = cpmExp+cpmSD), 
                  width = 0.15, color = "steelblue", linewidth = 0.8) +
    labs(x = "", y = "mean CPM") +
    facet_wrap(.~ geneName, scales = "free_y", nrow = 4) +
    theme_classic() +
    theme(plot.margin = margin(1,1,1,1,unit = "cm"),
          axis.text.x = element_text(family = "serif", size = 12, angle = 60, vjust = 1, hjust = 1, colour = "black"),
          strip.text.x = element_text(size = 12, color = "red", face = "bold"),
          strip.background = element_blank(),
          panel.spacing.y = unit(0.8, "cm"))
  p_01
  ggsave(plot = p_01, 
         filename = paste0(path, "result02/geneExpr.pdf"), units = "cm", width = 12, height = 20)
}
################################################################################
#06、stringtie组装
{
  paste0("stringtie --merge -G ",inx3," -o ",path,"stringtie/stage.merge.gtf ",
         paste0(path,"stringtie/",arex,".gtf", collapse = " "),"\n")
}













