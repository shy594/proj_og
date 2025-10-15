#write by Sun Haiayng at 2023.10.27
#1、搭环境
{
  rm(list=ls())
  for (i in c("ReactomePA","ggplot2","stringr","ChIPseeker","DESeq2","clusterProfiler","topGO","Rgraphviz","pathview","org.Mm.eg.db","TxDb.Mmusculus.UCSC.mm10.knownGene")) {library(i, character.only = T)};rm(i)
}
#2、给索引
{
  path = "/ChIP_seq_2/aaSHY/ghyuu/ripseq/20231025/"
  inx1 = "/Reference/aaSHY/Ensembl99/GRCm38/INDEX/bowtie2/mm10"
  inx2 = "/Reference/aaSHY/Ensembl99/GRCm38/INDEX/STARv279a"
  inx3 = "/Reference/aaSHY/Ensembl99/GRCm38/fa/GRCm38.dna.primary_assembly.fa"
  inx4 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/Mus_musculus.GRCm38.99.gtf"
  inx5 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf"
  inx6 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf.locInd"
  inx7 = "/Reference/aaSHY/Ensembl99/GRCm38/BED/TEsubFamily/MM::ERVL::LTR.bed"
  inx8 = "/Reference/aaSHY/Ensembl99/GRCm38/BED/special/MM_zyw.bed"
  inx9 = "/Reference/aaSHY/Ensembl99/GRCm38/BED/TEsubFamily/MM22::ERVL::LTR.bed"
  inxA = "/Reference/aaSHY/Ensembl99/GRCm38/BED/special/allGene.bed"
  inxB = "/Reference/aaSHY/Ensembl99/GRCm38/BED/special/twoC.bed"
  inxC = "/Reference/aaSHY/zOther/Rn45s/INDEX/bowtie2/Rn45s"
  inxD = "/Reference/aaSHY/zOther/TEconsensus/MuERVL/bt2/MuERVL"
}
################################################################################
#3、跑命令
{
  for (i in c("src","fq","Log","peak","rRNA","rRNA/csem","fisher","deeptools","bam","bam/csem","trim","PLOT","bw/bcv","bw/bcp","count")) {
    if (dir.exists(paths = paste0(path, i))==F) {
      dir.create(path = paste0(path,i))
    }
  }
  shell <- paste0(path,"src/run_a.sh")
  cat("#!/bin/bash\n", file = shell)
  cmd_01 <- paste0("#cp ", path,"rawdata/RawFq/*.fq.gz ", path, "fq/", "\n")
  cmd_02 <- paste0("#rename s/E14_// ", path, "fq/*gz", "\n")
  prex <- gsub(pattern = "E14_", replacement = "", 
               x = sapply(str_split(list.files((paste0(path,"rawdata/RawFq")), pattern = "*_1.fq.gz$*"),"_1.fq.gz"), "[",1), fixed = T)
  cmd_03 <- paste0("#trim_galore --phred33 --fastqc --illumina ",
                   "--clip_R1 10 --clip_R2 10 --three_prime_clip_R1 3 --three_prime_clip_R2 3 ",
                   "--paired ", path,"fq/",prex,"_1.fq.gz ",path, 
                   "fq/", prex, "_2.fq.gz -o ", path, "trim", "\n")
  lay1 <- paste0(path,"trim/",prex,"_1_val_1.fq.gz")
  lay2 <- paste0(path,"trim/",prex,"_2_val_2.fq.gz")
  cmd_04 <- paste0("#bowtie2 -q --very-sensitive -k 5 --reorder -p 20 -x ",
                   inx1," --no-mixed --no-discordant -1 ",lay1," -2 ",lay2,
                   " |samtools view -F 4 -b > ", path,
                   "bam/",prex,".bam","\n")
  cmd_05 <- paste0("#run-csem -p 20 --no-extending-reads --bam ",
                  paste0(path,"bam/",prex,".bam")," 140 ",
                  paste0(path,"bam/csem/",prex),"\n")
  cmd_06 <- paste0("#bamCoverage --normalizeUsing CPM -p 20 --minMappingQuality 1 ",
                   "--samFlagExclude 256 -of bigwig -bs 20 -b ",
                   paste0(path,"bam/csem/",prex,".sorted.bam -o "),
                   paste0(path,"bw/bcv/",prex,".bw\n"))
  cmd_07 <- paste0("#bamCompare -p 20 --scaleFactorsMethod SES --sampleLength 100000 -b1 ",
                   paste0(path,"bam/csem/",prex[1],".sorted.bam -b2 "),
                   paste0(path,"bam/csem/",prex[2],".sorted.bam -o "),
                   paste0(path,"bw/bcp/ipVSinput.bw"),"\n")
  cmd_08 <- paste0("#computeMatrix scale-regions -p 20 -a 2000 -b 2000 -m 5000 ",
                   "--missingDataAsZero -p 20 ", "-R ", c(inx7,inx8), " -S ", path,
                   "bw/bcp/ipVSinput.bw -o ", path,"deeptools/", 
                   c(basename(inx7),basename(inx8)), ".mat.gz","\n")
  cmd_09 <- paste0("#computeMatrix scale-regions -p 20 -a 1000 -b 1000 -m 500 ",
                   "--missingDataAsZero -p 20 -R ", inx9, " -S ",path,
                   "bw/bcp/ipVSinput.bw -o ",path,"deeptools/",
                   basename(inx9),".mat.gz","\n")
  cmd_10 <- paste0("#computeMatrix scale-regions -p 20 -a 2000 -b 2000 -m 5000 ",
                   "--missingDataAsZero -p 20 -R ", inxA, " -S ",path,
                   "bw/bcp/ipVSinput.bw -o ",path,"deeptools/",
                   basename(inxA),".mat.gz","\n")
  cmd_11 <- paste0("#computeMatrix scale-regions -p 20 -a 5000 -b 5000 -m 30000 ",
                   "--missingDataAsZero -p 20 -R ", inxB, " -S ",path,
                   "bw/bcp/ipVSinput.bw -o ",path,"deeptools/",
                   basename(inxB),".mat.gz","\n")
  cmd_12 <- paste0("#plotHeatmap --heatmapHeight 12 --colorList white,red -m ",path,
                   "deeptools/",c(basename(inx7),basename(inx8),basename(inx9),
                                  basename(inxA),basename(inxB)),
                   ".mat.gz -out ",path,"PLOT/",
                   c(basename(inx7),basename(inx8),
                     basename(inx9),basename(inxA),basename(inxB)),".pdf","\n")
  cmd_13 <- paste0("#macs2 callpeak -f BAMPE -g mm -p 0.05 --keep-dup 1 -t ",
                   path,"bam/csem/",prex[1],".sorted.bam -c ",
                   path,"bam/csem/",prex[2],".sorted.bam -n ghyuuN --outdir ",
                   path,"peak/", "\n")
  cmd_14 <- paste0("#macs2 callpeak -f BAMPE -g mm -p 0.05 --broad --broad-cutoff 0.05 ",
                   "--keep-dup 1 -t ",path,"bam/csem/", prex[1],
                   ".sorted.bam -c ", path,"bam/csem/", prex[2],
                   ".sorted.bam -n ghyuuB --outdir ", path,"peak/", "\n")
  cmd_15 <- paste0("#TEtranscripts -t ",
                   paste0(path,"bam/csem/",prex[1],".sorted.bam")," -c ",
                   paste0(path,"bam/csem/",prex[2],".sorted.bam"),
                   " --GTF ",inx4," --TE ",inx5," --sortByPos --mode multi ",
                   "--project ccghyuu --outdir ",paste0(path,"count/"),"\n")
  cmd_16 <- paste0("#TElocal -b ",
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
  for (i in cmd_14) {cat(i, append = T, file = shell)}
  for (i in cmd_15) {cat(i, append = T, file = shell)}
  for (i in cmd_16) {cat(i, append = T, file = shell)}
  print(paste0("nohup bash ",shell," >",path,"Log/run_a.log"," 2>&1 &"))
  #
  shell <- paste0(path,"src/run_c.sh")
  cat("#!/bin/bash\n", file = shell)
  cmd_17 <- paste0("bowtie2 -q --very-sensitive -k 5 --reorder -p 20 -x ",
                   inxD," --no-mixed --no-discordant -1 ",lay1," -2 ",lay2,
                   " |samtools view -F 4 -b > ", path,
                   "MM/",prex,".bam","\n")
  cmd_18 <- paste0("run-csem -p 20 --no-extending-reads --bam ",
                   paste0(path,"MM/",prex,".bam")," 140 ",
                   paste0(path,"MM/csem/",prex),"\n")
  cmd_19 <- paste0("bamCoverage --normalizeUsing CPM -p 20 --minMappingQuality 1 ",
                   "--samFlagExclude 256 -of bigwig -bs 20 -b ",
                   paste0(path,"MM/csem/",prex,".sorted.bam -o "),
                   paste0(path,"MM/",prex,".MM.bw\n"))
  for (i in cmd_17) {cat(i, append = T, file = shell)}
  for (i in cmd_18) {cat(i, append = T, file = shell)}
  for (i in cmd_19) {cat(i, append = T, file = shell)}
  print(paste0("nohup bash ",shell," >",path,"Log/run_c.log"," 2>&1 &"))
}
################################################################################
#4、peak Analysis (pvalue 0.05)
{
  gtf1 <- rtracklayer::import(con = inx4, format = "gtf")
  gtf1 <- as.data.frame(gtf1)
  colnames(gtf1)[c(10,12)] <- c("gene_name","gene_id")
  gtf1 <- GenomicRanges::makeGRangesFromDataFrame(df = gtf1, keep.extra.columns = T)
  txb1 <- suppressWarnings(GenomicFeatures::makeTxDbFromGRanges(gr = gtf1, taxonomyId = 10090))
  twoC <- read.table(inxB)
  gtf2 <- gtf1[which(gtf1$gene_id %in% unique(twoC$V4))]
  txb2 <- suppressWarnings(GenomicFeatures::makeTxDbFromGRanges(gr = gtf2, taxonomyId = 10090))
  #
  pk01 <- suppressWarnings(ChIPseeker::annotatePeak(peak = paste0(path, "peak/ghyuuN_peaks.narrowPeak"), 
                                                    tssRegion = c(-2000, 2000), TxDb = txb1))
  pdf(file = paste0(path, "PLOT/ghyuuN_chipseekerAnno_all.pdf"), width = 8, height = 8)
  ChIPseeker::plotAnnoPie(pk01, legend.position = "other")
  dev.off()
  write.csv(as.data.frame(pk01), paste0(path, "peak/ghyuuN_chipseekerAnno_all.csv"), quote = F, row.names = F)
  pk02 <- suppressWarnings(ChIPseeker::annotatePeak(peak = paste0(path, "peak/ghyuuN_peaks.narrowPeak"), 
                                                    tssRegion = c(-2000, 2000), TxDb = txb2))
  pdf(file = paste0(path, "PLOT/ghyuuN_chipseekerAnno_twoC.pdf"), width = 8, height = 8)
  ChIPseeker::plotAnnoPie(pk02, legend.position = "rightside")
  dev.off()
  write.csv(as.data.frame(pk02), paste0(path, "peak/ghyuuN_chipseekerAnno_twoC.csv"), quote = F, row.names = F)
  #
  pk03 <- suppressWarnings(ChIPseeker::annotatePeak(peak = paste0(path, "peak/ghyuuB_peaks.broadPeak"),
                                                    tssRegion=c(-2000, 2000), TxDb = txb1))
  pdf(file = paste0(path, "PLOT/ghyuuB_chipseekerAnno_all.pdf"), width = 8, height = 8)
  ChIPseeker::plotAnnoPie(pk03, legend.position = "other")
  dev.off()
  write.csv(as.data.frame(pk03), paste0(path, "peak/ghyuuB_chipseekerAnno_all.csv"), quote = F, row.names = F)
  pk04 <- suppressWarnings(ChIPseeker::annotatePeak(peak = paste0(path, "peak/ghyuuB_peaks.broadPeak"),
                                                    tssRegion=c(-2000, 2000), TxDb = txb2))
  pdf(file = paste0(path, "PLOT/ghyuuB_chipseekerAnno_twoC.pdf"), width = 8, height = 8)
  ChIPseeker::plotAnnoPie(pk04, legend.position = "rightside")
  dev.off()
  write.csv(as.data.frame(pk04), paste0(path, "peak/ghyuuB_chipseekerAnno_twoC.csv"), quote = F, row.names = F)
}
################################################################################
#5、TE Binding situation [适用于ChIP-seq]
{
  #Observed
  red1 <- read.table(file = paste0(path,"count/ccghyuu.cntTable"), sep = "\t", stringsAsFactors = F, header = T, row.names = 1)
  colnames(red1) <- sapply(str_split(colnames(red1),"\\."),"[",9)
  red2 <- red1[grep(rownames(red1), pattern = ":", invert = F),]
  #Expected
  mmu1 <- read.table(file = paste0(inx3,".fai"))
  refL <- sum(mmu1$V2)
  mmu2 <- read.table(file = inx5, sep = "\t", stringsAsFactors = F, header = F)
  mmu3 <- mmu2
  mmu3[,9] <- gsub("; transcript_id.*family_id |; class_id ",":",gsub("gene_id |;$","",mmu3[,9]))
  mmu3[,8] <- mmu3[,5]-mmu3[,4]+1
  for (i in seq_along(rownames(red2))){
    print(i)
    red2[i,3] <- round(sum(mmu3[which(mmu3[,9] == rownames(red2)[i]),8])/refL*sum(red1[,1]),4)
    red2[i,4] <- round(sum(mmu3[which(mmu3[,9] == rownames(red2)[i]),8])/refL*sum(red1[,2]),4)
    red2[i,5] <- suppressWarnings(fisher.test(matrix(as.numeric(red2[i,1:4]),ncol = 2,byrow = F),alternative = "two.sided")$p.value)
    red2[i,6] <- suppressWarnings(fisher.test(matrix(as.numeric(red2[i,1:4]),ncol = 2,byrow = F),alternative = "greater")$p.value)
    red2[i,7] <- suppressWarnings(fisher.test(matrix(as.numeric(red2[i,1:4]),ncol = 2,byrow = F),alternative = "less")$p.value)
  }
  colnames(red2) <- c("ghyuuIP","ghyuuInput","expcIP","expcInput","pvalue","greate","lesser")
  #
  red2$pMark[red2$pvalue<0.05] <- "arresting"
  red2$gMark[red2$greate<0.05] <- "arresting"
  red2$lMark[red2$lesser<0.05] <- "arresting"
  red2$repeatName <- rownames(red2)
  openxlsx::write.xlsx(x = red2, asTable = T, file = paste0(path,"fisher/ghyuuFisher.xlsx"))
}
################################################################################
#6、TE Binding situation [RIP-seq Specific]
{
  #Observed
  red1 <- read.table(file = paste0(path,"count/ccghyuu.cntTable"), sep = "\t", stringsAsFactors = F, header = T, row.names = 1)
  colnames(red1) <- sapply(str_split(colnames(red1),"\\."),"[",9)
  red2 <- red1[grep(rownames(red1), pattern = ":", invert = F),]
  #Expected
  mmu1 <- read.table(file = "/rna_seq_1/zx/ghyuu_RNA-LXM/Count/mm10_repname.saf_cut", sep = "\t", header = T, stringsAsFactors = F, row.names = 1)
  mmu2 <- as.data.frame(apply(mmu1, 2, function(x){x/sum(x)*1000000}))[,c(1,2)]
  mmu2$cpm <- rowMeans(mmu2)
  mmu2 <- mmu2[which(mmu2$cpm>=1),]
  red2 <- red2[which(sapply(str_split(rownames(red2), pattern = ":"), "[", 1) %in% rownames(mmu2)),]
  mmu3 <- mmu2[which(rownames(mmu2) %in% sapply(str_split(rownames(red2), pattern = ":"), "[", 1)),]
  sumReads <- sum(rowMeans(mmu1[which(rownames(mmu1) %in% rownames(mmu3)),c(1,2)]))
  for (i in rownames(red2)){
    print(i)
    m = str_split(i, pattern = ":")[[1]][1]
    red2[i,3] <- round(as.numeric(rowMeans(mmu1[m,c(1,2)])*sum(red2$ff23k_IP))/sumReads,4)
    red2[i,4] <- round(as.numeric(rowMeans(mmu1[m,c(1,2)])*sum(red2$Input))/sumReads,4)
    red2[i,5] <- suppressWarnings(fisher.test(matrix(as.numeric(red2[i,1:4]),ncol = 2,byrow = F),alternative = "two.sided")$p.value)
    red2[i,6] <- suppressWarnings(fisher.test(matrix(as.numeric(red2[i,1:4]),ncol = 2,byrow = F),alternative = "greater")$p.value)
    red2[i,7] <- suppressWarnings(fisher.test(matrix(as.numeric(red2[i,1:4]),ncol = 2,byrow = F),alternative = "less")$p.value)
  }
  colnames(red2) <- c("ghyuuIP","ghyuuInput","expcIP","expcInput","pvalue","greate","lesser")
  red2$pMark[red2$pvalue<0.05] <- "arresting"
  red2$gMark[red2$greate<0.05] <- "arresting"
  red2$lMark[red2$lesser<0.05] <- "arresting"
  red2$geneName <- rownames(red2)
  red2 <- red2[order(red2$greate, decreasing = F),]
  openxlsx::write.xlsx(x = red2, asTable = T, file = paste0(path,"fisher/ghyuuFisher_TE_v2.xlsx"))
}
################################################################################
#7、基因的Binding situation [适用于ChIP-seq]
{
  #Observed
  red1 <- read.table(file = paste0(path,"count/ccghyuu.cntTable"), sep = "\t", stringsAsFactors = F, header = T, row.names = 1)
  colnames(red1) <- sapply(str_split(colnames(red1),"\\."),"[",9)
  red2 <- red1[grep(rownames(red1), pattern = ":", invert = T),]
  #Expected
  mmu1 <- read.table(file = paste0(inx3,".fai"))
  refL <- sum(mmu1$V2)
  mmu2 <- rtracklayer::import(con = inx4,format = "gtf")
  mmu3 <- as.data.frame(mmu2[which(mmu2$type=="exon")])
  mmu3[1,]
  for (i in seq_along(rownames(red2))){
    print(i)
    red2[i,3] <- round(sum(mmu3[which(mmu3[,12] == rownames(red2)[i]),4])/refL*sum(red1[,1]),4)
    red2[i,4] <- round(sum(mmu3[which(mmu3[,12] == rownames(red2)[i]),4])/refL*sum(red1[,2]),4)
    red2[i,5] <- suppressWarnings(fisher.test(matrix(as.numeric(red2[i,1:4]),ncol = 2,byrow = F),alternative = "two.sided")$p.value)
    red2[i,6] <- suppressWarnings(fisher.test(matrix(as.numeric(red2[i,1:4]),ncol = 2,byrow = F),alternative = "greater")$p.value)
    red2[i,7] <- suppressWarnings(fisher.test(matrix(as.numeric(red2[i,1:4]),ncol = 2,byrow = F),alternative = "less")$p.value)
  }
  colnames(red2) <- c("ghyuuIP","ghyuuInput","expcIP","expcInput","pvalue","greate","lesser")
  red2$pMark[red2$pvalue<0.05] <- "arresting"
  red2$gMark[red2$greate<0.05] <- "arresting"
  red2$lMark[red2$lesser<0.05] <- "arresting"
  red2$geneName <- rownames(red2)
  red2 <- red2[order(red2$greate, decreasing = F),]
  openxlsx::write.xlsx(x = red2, asTable = T, file = paste0(path,"fisher/ghyuuFisher_gene.xlsx"))
  ##
  gen1 <- na.omit(mapIds(x=org.Mm.eg.db, keys=unique(red2$geneName[1:1000]), keytype="SYMBOL", column="ENTREZID"))
  pat1 <- clusterProfiler::enrichKEGG(gene = gen1, organism = "mmu", 
                                      keyType = "kegg", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
  kegg <- pat1@result[1:15,]
  kegg$Description <- gsub(pattern = " - Mus musculus (house mouse)", replacement = "", x = kegg$Description, fixed = T)
  pat2 <- ReactomePA::enrichPathway(gene = gen1, 
                                    organism = "mouse", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
  reac <- pat2@result[1:15,]
  gogo <- clusterProfiler::enrichGO(gene = gen1, OrgDb = org.Mm.eg.db, 
                                    keyType = "ENTREZID", ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
  term <- gogo@result[1:15,]
  kegg$Description <- factor(kegg$Description, levels = rev(rbind(kegg$Description)))
  reac$Description <- factor(reac$Description, levels = rev(rbind(reac$Description)))
  term$Description <- factor(term$Description, levels = rev(rbind(term$Description)))
  p_01 <- ggplot(kegg, aes(x=-log10(pvalue), y=Description))+
    geom_bar(stat = "identity", fill="royalblue", width = 0.8) +
    theme_classic() +
    labs(title = "fisher gene top 1000: KEGG Pathway (top 15)", y="") +
    theme(text = element_text(size = 12,family = "serif",color = "black"),
          axis.title = element_text(size = 12,family = "serif",color = "black"),
          axis.text  = element_text(size = 12,family = "serif",color = "black"))
  p_02 <- ggplot(reac, aes(x=-log10(pvalue), y=Description))+
    geom_bar(stat = "identity", fill="royalblue", width = 0.8) +
    theme_classic() +
    labs(title = "fisher gene top 1000: Reactome Pathway (top 15)", y="") +
    theme(text = element_text(size = 12,family = "serif",color = "black"),
          axis.title = element_text(size = 12,family = "serif",color = "black"),
          axis.text  = element_text(size = 12,family = "serif",color = "black"))
  p_03 <- ggplot(term, aes(x=-log10(pvalue), y=Description))+
    geom_bar(stat = "identity", fill="royalblue", width = 0.8) +
    theme_classic() +
    labs(title = "fisher gene top 1000: GO BP term (top 15)", y="") +
    theme(text = element_text(size = 12,family = "serif",color = "black"),
          axis.title = element_text(size = 12,family = "serif",color = "black"),
          axis.text  = element_text(size = 12,family = "serif",color = "black"))
  ggsave(plot = p_01, units = "cm", width = 20, height = 16, filename = paste0(path, "fisher/kegg.pdf"))
  ggsave(plot = p_02, units = "cm", width = 20, height = 16, filename = paste0(path, "fisher/reac.pdf"))
  ggsave(plot = p_03, units = "cm", width = 20, height = 16, filename = paste0(path, "fisher/gobp.pdf"))
}
################################################################################
#8、基因的Binding situation [RIP-seq Specific]
{
  #Observed
  red1 <- read.table(file = paste0(path,"count/ccghyuu.cntTable"), sep = "\t", stringsAsFactors = F, header = T, row.names = 1)
  colnames(red1) <- sapply(str_split(colnames(red1),"\\."),"[",9)
  red2 <- red1[grep(rownames(red1), pattern = ":", invert = T),]
  #Expected
  mmu1 <- read.table(file = "/rna_seq_1/zx/ghyuu_RNA-LXM/Count/gencode.vM21.primary_assembly.annotation.gtf_cut", sep = "\t", header = T, stringsAsFactors = F, row.names = 1)
  mmu2 <- as.data.frame(apply(mmu1, 2, function(x){x/sum(x)*1000000}))[,c(1,2)]
  mmu2$cpm <- rowMeans(mmu2)
  mmu2 <- mmu2[which(mmu2$cpm>=1),]
  red2 <- red2[which(rownames(red2) %in% rownames(mmu2)),]
  mmu3 <- mmu2[which(rownames(mmu2) %in% rownames(red2)),]
  sumReads <- sum(rowMeans(mmu1[which(rownames(mmu1) %in% rownames(mmu3)),c(1,2)]))
  for (i in rownames(red2)){
    print(i)
    red2[i,3] <- round(as.numeric(rowMeans(mmu1[i,c(1,2)])*sum(red2$ff23k_IP))/sumReads,4)
    red2[i,4] <- round(as.numeric(rowMeans(mmu1[i,c(1,2)])*sum(red2$Input))/sumReads,4)
    red2[i,5] <- suppressWarnings(fisher.test(matrix(as.numeric(red2[i,1:4]),ncol = 2,byrow = F),alternative = "two.sided")$p.value)
    red2[i,6] <- suppressWarnings(fisher.test(matrix(as.numeric(red2[i,1:4]),ncol = 2,byrow = F),alternative = "greater")$p.value)
    red2[i,7] <- suppressWarnings(fisher.test(matrix(as.numeric(red2[i,1:4]),ncol = 2,byrow = F),alternative = "less")$p.value)
  }
  colnames(red2) <- c("ghyuuIP","ghyuuInput","expcIP","expcInput","pvalue","greate","lesser")
  red2$pMark[red2$pvalue<0.05] <- "arresting"
  red2$gMark[red2$greate<0.05] <- "arresting"
  red2$lMark[red2$lesser<0.05] <- "arresting"
  red2$geneName <- rownames(red2)
  red2 <- red2[order(red2$greate, decreasing = F),]
  openxlsx::write.xlsx(x = red2, asTable = T, file = paste0(path,"fisher/ghyuuFisher_gene_v4.xlsx"))
}
################################################################################
#9、RIP、RNA-seq基因的交集
{
  rips <- openxlsx::read.xlsx(xlsxFile = paste0(path, "fisher/ghyuuFisher_TE_v2.xlsx"))
  rips <- rips[which(rips$greate <0.001),]
  ghyuur <- read.csv(file = "/rna_seq_1/zx/ghyuu_RNA-LXM/TE/ghyuu_KO_TE_G.csv")
  ghyuu1 <- ghyuur[which(ghyuur$log2FoldChange >=0 & ghyuur$pvalue <0.05),]
  ghyuu2 <- ghyuur[which(ghyuur$log2FoldChange < 0 & ghyuur$pvalue <0.05),]
  abce <- read.csv(file = "/disk5/fx/RNA-seq/Aef_RNAseq/Count/repeat_analysis/Nogene/repNogeneres_Aef_h_M.csv")
  abc1 <- abce[which(abce$log2FoldChange >=0 & abce$pvalue <0.05),]
  abc2 <- abce[which(abce$log2FoldChange < 0 & abce$pvalue <0.05),]
  #
  temp <- list()
  temp[["x"]] = sapply(str_split(rips$geneName, pattern = ":"), "[", 1)
  temp[["y"]] = ghyuu2$X
  temp[["z"]] = abc2$X
  p_01 <- VennDiagram::venn.diagram(x = temp, filename = NULL, 
                                    main = "VENN DIAGRAM (TE)",
                                    sub = "This is subtitle...(ghyuuRIP & ghyuuUP & AefUP)", 
                                    main.fontfamily = "serif", sub.fontfamily = "serif", 
                                    cex = 1.5, main.cex = 2, sub.cex = 2,
                                    print.mode = "raw",#c("percent", "raw"),
                                    category.names = c("RIP-seqFisher", "ghyuuUP", "AefUP"), fill=c("#FFFFCC","#CCFFFF","#FFCCCC"),
                                    units = "cm", height = 12, width = 12)#c('#FFFFCC','#CCFFFF',"#FFCCCC")#c('#FFFFCC','#CCFFFF',"#FFCCCC","#CCCCFF", "#CCFFCC")
  pdf(file = paste0(path, "fisher/VENN/ghyuuRIP_ghyuuUP_AefUP_TE_V2.pdf"))
  grid.draw(p_01)
  dev.off()
  #
  intersect(temp[["x"]], temp[["y"]])
  intersect(temp[["x"]], temp[["z"]])
  intersect(temp[["y"]], temp[["z"]])
  intersect(intersect(temp[["x"]], temp[["y"]]),temp[["z"]])
}










