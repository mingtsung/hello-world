#!/usr/bin/env Rscript

rm(list = ls()) # remove all variables in workspace

# ---------------------------------------------------------------------------------------------------------------------------------------------------
#current_folder = "~/Project/DrugPairing"
current_folder = "/Users/mingtsung/Desktop/Project/DrugPairing"
setwd(current_folder)

#install.packages("mgsub") # mgsub: Safe, Multiple, Simultaneous String Substitution
#install.packages("RColorBrewer") # RColorBrewer: ColorBrewer Palettes
#install.packages("corrplot") # corrplot: Visualization of a Correlation Matrix
#install.packages("openxlsx") # openxlsx: Read, Write and Edit xlsx Files
#install.packages("gtools")  # gtools: Various R Programming Tools
#install.packages("stringr") # stringr: Simple, Consistent Wrappers for Common String Operations

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("ComplexHeatmap")

#suppressMessages(library(mgsub)) # mgsub()
suppressMessages(library(RColorBrewer))
suppressMessages(library(corrplot))
suppressMessages(library(openxlsx)) # read.xlsx()
suppressMessages(library(gtools)) # mixedsort(), mixedorder()
suppressMessages(library(stringr)) # str_to_title(), str_to_upper()
suppressMessages(library(ComplexHeatmap)) # Heatmap()

mainDir <- current_folder
subDir <- "Result"
if (! file.exists(file.path(mainDir, subDir))) {
    dir.create(file.path(mainDir, subDir))
}


# --- Source function -------------------------------------------------------------------------------------------------------------------------------
source('./Script/function/QC/qc_data.R')
source('./Script/function/Monocle3/Monocle3_analysis.R')
source('./Script/function/Seurat/Seurat_analysis.R')
source('./Script/function/Seurat/StackedVlnPlot.R')
source('./Script/function/plot/volcano_plot.R')


# --- QC criterion (gene cutoff) --------------------------------------------------------------------------------------------------------------------
gene_cutoff <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
gene_cutoff.default = gene_cutoff[2] # min_gene_frac=0.1


# --- Data source -----------------------------------------------------------------------------------------------------------------------------------
Species <- c('Mouse', 'Human')


# --- Day information -------------------------------------------------------------------------------------------------------------------------------
sample <- c("control", 
            "ch02", # Day 5
            
            # /// Day 3 ///
            "d3_p1_1", # Day 3, repeat 1-1
            "20190919_d3_p1_1", # Day 3, repeat 1-2
            "poor_d3_p1_1", # Day 3, repeat 1, poor
            
            "d3_p1_2", # Day 3, repeat 2
            
            # /// Day 5 ///
            "d5_p1", # Day 5, subculture 1
            "poor_d5_p1", # Day 5, subculture 1, poor
            
            "d5_p2", # Day 5, subculture 2
            
            # /// Day 7 ///
            "d7_p1", # Day 7, subculture 1
            "poor_d7_p1", # Day 7, subculture 1, poor
            
            "d7_p2_1", # Day 7, subculture 2, repeat 1
            "d7_p2_2") # Day 7, subculture 2, repeat 2

# /// Day 3 ///
d3_p1_1.good <- c("d3_p1_1", # Day 3, repeat 1-1
                  "20190919_d3_p1_1") # Day 3, repeat 1-2
d3_p1_1.poor <- c("poor_d3_p1_1") # Day 3, repeat 1, poor
#d3_p1_1 <- c(d3_p1_1.good)
d3_p1_1 <- c(d3_p1_1.good, d3_p1_1.poor)
d3_p1 <- c(d3_p1_1, 
           "d3_p1_2") # Day 3, repeat 2
d3 <- d3_p1

# /// Day 5 ///
d5_p1.good <- c("d5_p1") # Day 5, subculture 1
d5_p1.poor <- c("poor_d5_p1") # Day 5, subculture 1, poor
#d5_p1 <- c(d5_p1.good)
d5_p1 <- c(d5_p1.good, d5_p1.poor)
d5 <- c(d5_p1, 
        "d5_p2") # Day 5, subculture 2

# /// Day 7 ///
d7_p1.good <- c("d7_p1") # Day 7, subculture 1
d7_p1.poor <- c("poor_d7_p1") # Day 7, subculture 1, poor
#d7_p1 <- c(d7_p1.good)
d7_p1 <- c(d7_p1.good, d7_p1.poor)
d7_p2 <- c("d7_p2_1", # Day 7, subculture 2, repeat 1
           "d7_p2_2") # Day 7, subculture 2, repeat 2
d7 <- c(d7_p1, 
        d7_p2)

# [R] How to insert one element into a vector?
# https://stat.ethz.ch/pipermail/r-help/2004-November/061337.html
insert <- function(v,e,pos){
    if (pos <= length(v)) {
        new.v = c(v[1:(pos-1)],e,v[(pos):length(v)])
    } else {
        new.v = c(v[1:(pos-1)],e)
    }
    return(new.v)
}
sample.merged <- c("d3_p1_1_merged", # Day 3, no "d3_p1_2"
                   "d3", # Day 3
                   
                   "d5_p1_merged", # Day 5, no "d5_p2"
                   "d5", # Day 5
                   
                   "d7_p1_merged", # Day 7, no "d7_p2_1" and "d7_p2_2"
                   "d7_p2", # Day 7, no "d7_p1"
                   "d7") # Day 7
sample.reorder <- insert(sample, sample.merged[which(sample.merged %in% c("d3_p1_1_merged"))], ifelse(length(which(sample %in% c(d3_p1_1.good, d3_p1_1.poor))) > 0, max(which(sample %in% c(d3_p1_1.good, d3_p1_1.poor))) + 1, length(sample) + 1))
sample.reorder <- insert(sample.reorder, sample.merged[which(sample.merged %in% c("d3"))], ifelse(length(which(sample.reorder %in% c(d3, sample.merged[which(sample.merged %in% c("d3_p1_1_merged"))]))) > 0, max(which(sample.reorder %in% c(d3, sample.merged[which(sample.merged %in% c("d3_p1_1_merged"))]))) + 1, length(sample.reorder) + 1))
sample.reorder <- insert(sample.reorder, sample.merged[which(sample.merged %in% c("d5_p1_merged"))], ifelse(length(which(sample.reorder %in% c(d5_p1.good, d5_p1.poor))) > 0, max(which(sample.reorder %in% c(d5_p1.good, d5_p1.poor))) + 1, length(sample.reorder) + 1))
sample.reorder <- insert(sample.reorder, sample.merged[which(sample.merged %in% c("d5"))], ifelse(length(which(sample.reorder %in% c(d5, sample.merged[which(sample.merged %in% c("d5_p1_merged"))]))) > 0, max(which(sample.reorder %in% c(d5, sample.merged[which(sample.merged %in% c("d5_p1_merged"))]))) + 1, length(sample.reorder) + 1))
sample.reorder <- insert(sample.reorder, sample.merged[which(sample.merged %in% c("d7_p1_merged"))], ifelse(length(which(sample.reorder %in% c(d7_p1.good, d7_p1.poor))) > 0, max(which(sample.reorder %in% c(d7_p1.good, d7_p1.poor))) + 1, length(sample.reorder) + 1))
sample.reorder <- insert(sample.reorder, sample.merged[which(sample.merged %in% c("d7_p2"))], ifelse(length(which(sample.reorder %in% d7_p2)) > 0, max(which(sample.reorder %in% d7_p2)) + 1, length(sample.reorder) + 1))
sample.reorder <- insert(sample.reorder, sample.merged[which(sample.merged %in% c("d7"))], ifelse(length(which(sample.reorder %in% c(d7, sample.merged[which(sample.merged %in% c("d7_p1_merged","d7_p2"))]))) > 0, max(which(sample.reorder %in% c(d7, sample.merged[which(sample.merged %in% c("d7_p1_merged","d7_p2"))]))) + 1, length(sample.reorder) + 1))


# --- cell type -------------------------------------------------------------------------------------------------------------------------------------
# healthy (acinar, ductal, beta)
# cancer (double-mutant: PCC8, PCC9)(triple-mutant: MIN6)

# /// counts_csv ///
acinar.counts <- c("acinar")
ductal.counts <- c("counts_kcppc_new")
cancer_tripleMut.counts <- c("cancer")
beta.counts <- c("min6")
cancer_doubleMut.counts <- c("pcc8","pcc9")
healthy.counts <- c(acinar.counts, ductal.counts, beta.counts)
cancer.counts <- c(cancer_doubleMut.counts, cancer_tripleMut.counts)
control.counts <- c(healthy.counts, cancer.counts)

# /// old ///
acinar.old <- c("1_acinar")
ductal.old <- c("1_ductal") # old ductal cell (old processing pipeline)
cancer_tripleMut.old <- c("1_cancer")
beta.old <- c("16726_MIN6","16728_MIN6")
cancer_doubleMut.old <- c("16728_PCC8","16726_PCC9")
healthy.old <- c(acinar.old, ductal.old, beta.old)
cancer.old <- c(cancer_doubleMut.old, cancer_tripleMut.old)
control.old <- c(healthy.old, cancer.old)

# /// new ///
acinar.new <- c("1_acinar_new")
ductal.new <- c("kcppc_bam_counts_new") # old ductal cell (new processing pipeline, with gene_names)
cancer_tripleMut.new <- c("1_cancer_new")
beta.new <- c("16726_min6","16728_min6")
cancer_doubleMut.new <- c("16728_pcc8","16726_pcc9")
healthy.new <- c(acinar.new, ductal.new, beta.new)
cancer.new <- c(cancer_doubleMut.new, cancer_tripleMut.new)
control.new <- c(healthy.new, cancer.new)

# /// paper ///
paper.fibroblast <- c("SRR8872328","SRR8872329")    # Fibroblast-enriched fraction (DAPI-CD45-CD31-EPCAM-)
paper.viable <- c("SRR8872326","SRR8872327")        # Viable cell fraction (DAPI-)
paper.cell <- c(paper.fibroblast, paper.viable)

# /// control ///
#control <- control.counts
#control <- control.old
#control <- control.new
control <- c(control.counts, paper.cell)
#control <- c(control.old, paper.cell)
#control <- c(control.new, paper.cell)


# --- sc26 time-course sgRNA count matrices ---------------------------------------------------------------------------------------------------------
sgRNA.barcode <- c("1_TCGCGACACAAAAGTGTA", 
                   "2_TTAAGAAGAGAACGGTGT", 
                   "3_GTCGTACCGAAGGGGCGC", 
                   "4_AAACGCTCAATCGTCGGG", 
                   "5_CATACCCGTTGCGAGAGA", 
                   "6_CCCGGCTAGATAACTTGG", 
                   "7_AGCGCTGACGCCGCATC", 
                   "8_AACATGATATGCCGCGCG", 
                   "9_GCATTGCTAGCGGGCCTT", 
                   "10_ATGACATCAAGTCGCGTA", 
                   "11_CCCAGGGAGTTATTCGAT", # AB2
                   "12_ATCCAAGTACCACCCCGA", # AB2
                   "13_TAAAATACAATCCCTCGG", # toxic sgRNA (strong toxic)
                   "14_CCACAAGAATACATGCGT", # toxic sgRNA (strong toxic)
                   "15_GTGCACACTAACAATGAG", # negative control sgRNA 1
                   "16_GGCTAGAATGTGTACCAT") # negative control sgRNA 2

sgRNA.barcode.full <- c("1_CACCCCTACACTCCTCCG", "2_CATGGGAAACGGTGTAAG", "3_CGTCGATCGGGTCCAAG", "4_CTAAGCTGCCGATCCGTG", "5_CTCCAGCACTCTACCTGA", 
                        "6_CCGGTTACAGTATTACTG", "7_CATCCGGCCAACCTCGCG", "8_AGGTAATTGGCAGATCCT", "9_AAATGGCTAAGGGTGGTG", "10_ATTACGCCGAAGAAAGAG", 
                        "11_TGATCGAATGTGGCAAGG", "12_GGACCTGATTCATAGGAG", "13_AGAATACGGAGTCCGACA", "14_GCCAGACTTATTTCGGTC", "15_TCGCGACACAAAAGTGTA", 
                        "16_TTAAGAAGAGAACGGTGT", "17_TCTAAAACCATGGCCTC", "18_AAGAACAGTGGCATCCGG", "19_AGGCTTACACCATCTGTG", "20_TATTTGCAACAGCGCCAG", 
                        "21_TCGATACCAGTTACAGTC", "22_TCGCGTAGGACATGTGCG", "23_GTCGTACCGAAGGGGCGC", "24_AAACGCTCAATCGTCGGG", "25_GGGTCCCCCCTCCCTGCG", 
                        "26_TCCTATATTGCGCTACGA", "27_TCGTAGCTCATATGAATG", "28_TCCCTCCCCCTGCGGTCT", "29_ACCAAGATACGTTCCGAC", "30_TGCGTCACTAAACTTATG", 
                        "31_AGCGCTGACGCCGCATC", "32_AACATGATATGCCGCGCG", "33_CCCGCACCGTGAGGCCTT", "34_CATGTATGTGGCTTCCGG", "35_CATACCCGTTGCGAGAGA", 
                        "36_CCCGGCTAGATAACTTGG", "37_GCATTGCTAGCGGGCCTT", "38_ATGACATCAAGTCGCGTA", 
                        "39_CGGATGCTAATGCACGTG", # toxic sgRNA
                        "40_AGACGTCGTATTTGTTGT", 
                        "41_ATAGCGATGATAAAGACG", "42_TAAGTGACGGCTTTGATG", "43_CTGCACTTCGGGGTCGGT", "44_AATAACCGTAAATTGAGT", "45_CCGCCGAAGCATGAAAAA", 
                        "46_CGGACCGTCAATGGTGG", "47_AGCGTCAACCGACTTCTG", "48_CACGCTCGTGGCACCGTC", "49_CCCACACACACAGCGTGT", "50_CCTGAAGCGAGTGTGAGG", 
                        "51_TGTGTACCTCGGCGCAGA", "52_TTCCGCTATTGTCCGATG", "53_TAGAAAGCGATTATACGC", "54_ACGCGAAACCGGGCGTTA", "55_GTCACAGTAACTCGGCAT", 
                        "56_TAAGACGGCTATTGACAA", "57_ACGTGGATCATAGCACGC", "58_TTAGTAGACTGCGAACTG", "59_GCGCTACGACCGGAATCG", "60_GCCCTGCAACCTGTTTTA", 
                        "61_CGAATGCTGAGAATAGTC", "62_CACGCGTCGGCCTATCAC", "63_CCCTCGATACGCAATTCA", 
                        "64_CCCAGGGAGTTATTCGAT", # AB2
                        "65_ATCCAAGTACCACCCCGA", # AB2
                        "66_TAAAATACAATCCCTCGG", # toxic sgRNA (strong toxic)
                        "67_GATGTATCGGAGTAGTTG", # toxic sgRNA (strong toxic)
                        "68_CCACAAGAATACATGCGT", # toxic sgRNA (strong toxic)
                        "69_GGAGAATCGAGATGGTGG", # toxic sgRNA (strong toxic)
                        "70_TCGCTACTACCTACTAGG", 
                        "71_GTGCACACTAACAATGAG", # negative control sgRNA 1
                        "72_GGCTAGAATGTGTACCAT") # negative control sgRNA 2

#sgRNA.barcode = sgRNA.barcode
sgRNA.barcode = sgRNA.barcode.full

# https://www.r-bloggers.com/how-to-expand-color-palette-with-ggplot-and-rcolorbrewer/
colourCount = length(sgRNA.barcode)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
sgRNA.barcode.color <- getPalette(colourCount)


# --- sgRNA cell barcodes ---------------------------------------------------------------------------------------------------------------------------
sample_folder = "crispr"
filename = paste('ch02 sg cell barcodes','_','12May2021', sep = '')
sg_cell_barcodes <- read.xlsx(paste(current_folder,'/Data','/','Mouse','/',sample_folder,'/',filename,'.xlsx', sep = ''), sheet = filename, startRow = 1, check.names = FALSE, sep.names = " ")

barcode.num.table <- c()

sg_cell_barcode.all <- list()
for (i in 1:nrow(sg_cell_barcodes)) {
    cell_barcode.idx <- which(colnames(sg_cell_barcodes) == "cell barcode")
    sg_cell_barcode <- sg_cell_barcodes[i,c(cell_barcode.idx:ncol(sg_cell_barcodes))]
    sg_cell_barcode <- sg_cell_barcode[which(!is.na(sg_cell_barcode))]
    sg_cell_barcode <- as.character(sg_cell_barcode)
    sg_cell_barcode <- trimws(sg_cell_barcode)
    
    barcode.num = length(sg_cell_barcode)
    names(barcode.num) <- paste(sg_cell_barcodes[i,"X1"],'_',trimws(sg_cell_barcodes[i,"sgrna barcode"]), sep = '')
    barcode.num.table <- c(barcode.num.table, barcode.num)
    
    sg_cell_barcode <- list(sg_cell_barcode)
    names(sg_cell_barcode) = paste(sg_cell_barcodes[i,"X1"],'_',trimws(sg_cell_barcodes[i,"sgrna barcode"]), sep = '')
    
    sg_cell_barcode.all <- c(sg_cell_barcode.all, sg_cell_barcode)
}

cat('[','ch02','] Cell barcode number in "',filename,'.xlsx','": ',sum(barcode.num.table),'\n', sep = '')
print(barcode.num.table)
cat('\n')


# --- cell_metadata ---------------------------------------------------------------------------------------------------------------------------------
cell_metadata.all <- c()
for (i in 1:length(sg_cell_barcode.all)) {
    cell_metadata <- data.frame("sgRNA" = rep(names(sg_cell_barcode.all[i]), length(sg_cell_barcode.all[[i]])), 
                                "sgRNA_id" = rep(unlist(strsplit(names(sg_cell_barcode.all[i]), '_'))[1], length(sg_cell_barcode.all[[i]])), 
                                "sgRNA_barcode" = rep(unlist(strsplit(names(sg_cell_barcode.all[i]), '_'))[2], length(sg_cell_barcode.all[[i]])), 
                                "cell_barcode" = sg_cell_barcode.all[[i]], 
                                stringsAsFactors = FALSE)
    
    cell_metadata[,"sgRNA"] <- factor(cell_metadata[,"sgRNA"], levels = unique(cell_metadata[,"sgRNA"]))
    cell_metadata[,"sgRNA_id"] <- factor(cell_metadata[,"sgRNA_id"], levels = unique(cell_metadata[,"sgRNA_id"]))
    cell_metadata[,"sgRNA_barcode"] <- factor(cell_metadata[,"sgRNA_barcode"], levels = unique(cell_metadata[,"sgRNA_barcode"]))
    
    rownames(cell_metadata) <- cell_metadata[,"cell_barcode"]
    
    cell_metadata.all <- rbind(cell_metadata.all, cell_metadata)
}


# --- crispr ----------------------------------------------------------------------------------------------------------------------------------------
sample_folder = "crispr"
sample_file = "ch02_counts"
filename = paste(sample_file,'_new', sep = '')
sgRNA <- read.csv(paste(current_folder,'/Data','/','Mouse','/',sample_folder,'/',filename,'.csv', sep = ''), sep = ',', header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

total_cell_barcodes <- colnames(sgRNA)[which(! colnames(sgRNA) %in% c("ensembl","genes","gene_names"))]

mainDir <- paste(current_folder,'/Result', sep = '')
subDir <- 'Mouse'
if (! file.exists(file.path(mainDir, subDir))) {
    dir.create(file.path(mainDir, subDir))
}

mainDir <- paste(current_folder,'/Result','/','Mouse', sep = '')
subDir <- sample_folder
if (! file.exists(file.path(mainDir, subDir))) {
    dir.create(file.path(mainDir, subDir))
}

mainDir <- paste(current_folder,'/Result','/','Mouse','/',sample_folder, sep = '')
subDir <- "ch02"
if (! file.exists(file.path(mainDir, subDir))) {
    dir.create(file.path(mainDir, subDir))
}

barcode.num.new.table <- c()

for (i in 1:length(sg_cell_barcode.all)) {
    #if (length(sg_cell_barcode.all[[i]]) > 0) {
    #    sgRNA.colidx <- sapply(sg_cell_barcode.all[[i]], FUN = function(X) which(colnames(sgRNA) %in% X))
    #    sg_cell_barcode.idx <- which(sapply(sg_cell_barcode.all[[i]], FUN = function(X) X %in% colnames(sgRNA)))
    #}
    
    sgRNA.single <- sgRNA[,c("genes","gene_names",sg_cell_barcode.all[[i]])]
    
    sgRNA.single[,"ensembl"] = rownames(sgRNA.single)
    sgRNA.single <- sgRNA.single[,c("ensembl","genes","gene_names",sg_cell_barcode.all[[i]])]
    
    barcode.num.new = ncol(sgRNA.single) - 3
    names(barcode.num.new) <- names(sg_cell_barcode.all[i])
    barcode.num.new.table <- c(barcode.num.new.table, barcode.num.new)
    
    write.csv(sgRNA.single, file = paste(current_folder,'/Result','/','Mouse','/',sample_folder,'/','ch02','/',names(sg_cell_barcode.all[i]),'.csv', sep = ''), row.names = FALSE, quote = TRUE)
}

cat('[','ch02','] Cell barcode number in "',filename,'.csv','": ',length(total_cell_barcodes),'\n', sep = '')
cat('[','ch02','] Sum of cell barcode number in each sgRNA cell barcode',': ',sum(barcode.num.new.table),'\n', sep = '')
print(barcode.num.new.table)
cat('\n')

absent_cell_barcode.idx <- which(sapply(total_cell_barcodes, FUN = function(X) ! X %in% unlist(sg_cell_barcode.all)))
absent_cell_barcode <- total_cell_barcodes[absent_cell_barcode.idx]

if (length(absent_cell_barcode) > 0) {
    if (length(absent_cell_barcode) == 1) {
        write.table(absent_cell_barcode, file = paste(current_folder,'/Result','/','Mouse','/',sample_folder,'/','absent_cell_barcode','.txt', sep = ''), sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    } else if (length(absent_cell_barcode) > 1) {
        write.table(absent_cell_barcode, file = paste(current_folder,'/Result','/','Mouse','/',sample_folder,'/','absent_cell_barcodes','.txt', sep = ''), sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    }
}

sgRNA.all <- sgRNA[,c("genes","gene_names",unlist(sg_cell_barcode.all))]
sgRNA.all[,"ensembl"] = rownames(sgRNA.all)
sgRNA.all <- sgRNA.all[,c("ensembl","genes","gene_names",unlist(sg_cell_barcode.all))]

write.csv(sgRNA.all, file = paste(current_folder,'/Result','/','Mouse','/',sample_folder,'/',filename,'.','subset','.csv', sep = ''), row.names = FALSE, quote = TRUE)

# /// rename ///
sgRNA.all <- sgRNA.all[,which(colnames(sgRNA.all) != "ensembl")]

genes.idx <- which(colnames(sgRNA.all) == "genes")
if (length(genes.idx) == 0) {
    genes.NA = 1
    sgRNA.all[,"genes"] <- rownames(sgRNA.all)
} else {
    genes.NA = 0
    NA.idx <- sort(unique(c(which(sgRNA.all[,"genes"] == ""), grep('Row', sgRNA.all[,"genes"], value = FALSE))))
    sgRNA.all[NA.idx,"genes"] <- rownames(sgRNA.all)[NA.idx]
}
genes.idx <- which(colnames(sgRNA.all) == "genes")
if (genes.idx != 1) {
    col.idx <- c(genes.idx, 1:(ncol(sgRNA.all)-1))
    sgRNA.all <- sgRNA.all[, col.idx]
}
genes.idx <- which(colnames(sgRNA.all) == "genes")

gene_names.idx <- which(colnames(sgRNA.all) == "gene_names")
if (length(gene_names.idx) == 0) {
    gene_names.NA = 1
    sgRNA.all[,"gene_names"] <- rownames(sgRNA.all)
} else {
    gene_names.NA = 0
    NA.idx <- sort(unique(c(which(sgRNA.all[,"gene_names"] == ""), grep('Row', sgRNA.all[,"gene_names"], value = FALSE))))
    sgRNA.all[NA.idx,"gene_names"] <- rownames(sgRNA.all)[NA.idx]
}
gene_names.idx <- which(colnames(sgRNA.all) == "gene_names")
if (gene_names.idx != 2) {
    col.idx <- c(1, gene_names.idx, 2:(ncol(sgRNA.all)-1))
    sgRNA.all <- sgRNA.all[, col.idx]
}
gene_names.idx <- which(colnames(sgRNA.all) == "gene_names")

# /// gene: replace empty gene_names with ID and append ID to duplicated gene_names ///
gene.dup.idx <- which(duplicated(sgRNA.all[,"gene_names"]))
if (length(gene.dup.idx) > 0) {
    for (d in 1:length(gene.dup.idx)) {
        gene_name <- sgRNA.all[which(sgRNA.all[,"gene_names"] == sgRNA.all[gene.dup.idx[d],"gene_names"]),"gene_names"]
        gene_name <- unique(gene_name)
        
        if (length(grep(')$', sgRNA.all[which(sgRNA.all[,"gene_names"] == gene_name),"gene_names"], value = TRUE)) == 0) {
            sgRNA.all[which(sgRNA.all[,"gene_names"] == gene_name),"gene_names"] <- paste(sgRNA.all[which(sgRNA.all[,"gene_names"] == gene_name),"gene_names"], "(", sgRNA.all[which(sgRNA.all[,"gene_names"] == gene_name),"genes"], ")", sep = '')
        }
    }
}

rownames(sgRNA.all) <- sgRNA.all[,"genes"]
sgRNA.all[,"ensembl"] = rownames(sgRNA.all)
sgRNA.all <- sgRNA.all[,c("ensembl","genes","gene_names",unlist(sg_cell_barcode.all))]

write.csv(sgRNA.all, file = paste(current_folder,'/Result','/','Mouse','/',sample_folder,'/',filename,'.','subset','.','rename','.csv', sep = ''), row.names = FALSE, quote = TRUE)


# --- Read in "gene(s) of interest" and "signature gene(s)" -----------------------------------------------------------------------------------------
# /// query gene(s) ///
# *** KRAS ***
genesets_of_interest <- list()

geneset1_of_interest <- c("Kras","St6galnac6")
geneset1_of_interest <- list(unique(geneset1_of_interest))
names(geneset1_of_interest) = "geneset1_of_interest"
genesets_of_interest <- c(genesets_of_interest, geneset1_of_interest)

geneset2_of_interest <- c("Brca1","Brca2","Atm")
geneset2_of_interest <- list(unique(geneset2_of_interest))
names(geneset2_of_interest) = "geneset2_of_interest"
genesets_of_interest <- c(genesets_of_interest, geneset2_of_interest)

genes_of_interest <- genesets_of_interest
KRAS.genes_of_interest = genes_of_interest

# /// signature gene(s) ///
# *** KRAS ***
KRAS.signature <- c("SINGH_KRAS_DEPENDENCY_SIGNATURE","FRIDMAN_SENESCENCE_UP","HALLMARK_KRAS_SIGNALING_UP")

KRAS.human.list <- list()
KRAS.mouse.list <- list()
for (sg in 1:length(KRAS.signature)) {
    
    KRAS <- read.table(paste(current_folder,'/Data/Signature/','KRAS','/',KRAS.signature[sg],'.txt', sep = ''), sep = '\t', header = FALSE, skip = 2, check.names = FALSE, stringsAsFactors = FALSE)
    KRAS <- trimws(KRAS[,1])
    
    KRAS.human <- str_to_upper(KRAS)
    KRAS.mouse <- str_to_title(KRAS)
    
    KRAS.human.list <- c(KRAS.human.list, list(KRAS.human))
    KRAS.mouse.list <- c(KRAS.mouse.list, list(KRAS.mouse))
}
names(KRAS.human.list) = KRAS.signature
names(KRAS.mouse.list) = KRAS.signature

signature_genes = KRAS.mouse.list
KRAS.signature_genes = signature_genes

# *** PML ***
# PML: premalignant lesion
# Beane, Jennifer E., et al. "Molecular subtyping reveals immune alterations associated with progression of bronchial premalignant lesions." Nature communications 10.1 (2019): 1-13.
PML <- read.xlsx(paste(current_folder,'/Data/Signature/','PML','/','41467_2019_9834_MOESM3_ESM','.xlsx', sep = ''), startRow = 2, check.names = FALSE, sep.names = " ")

PML.ModuleNumber <- PML[,"Module Number"]
PML.NumberOfUniqueGeneSymbols <- PML[,"Number of Unique Gene Symbols"]
PML.CorrelationWithModuleEigengene <- trimws(PML[,"Correlation with Module Eigengene"])
PML.Database <- trimws(PML[,"Database"])
PML.Category <- trimws(PML[,"Category"])
PML.Pvalue <- PML[,"P-value"]
PML.FDR <- PML[,"FDR"]
PML.Genes <- trimws(PML[,"Genes"])

PML.ModuleGenes.pos.human.list <- list()
PML.ModuleGenes.neg.human.list <- list()
PML.ModuleGenes.human.list <- list()
PML.ModuleGenes.pos.mouse.list <- list()
PML.ModuleGenes.neg.mouse.list <- list()
PML.ModuleGenes.mouse.list <- list()

PML.ModuleID <- unique(PML.ModuleNumber)
for (M.No in 1:length(PML.ModuleID)) {
    PML.ModuleNumber.pos_idx <- which(PML.ModuleNumber == PML.ModuleID[M.No] & PML.CorrelationWithModuleEigengene == 'Positive')
    PML.ModuleNumber.neg_idx <- which(PML.ModuleNumber == PML.ModuleID[M.No] & PML.CorrelationWithModuleEigengene == 'Negative')
    
    PML.ModuleGenes.pos <- unique(trimws(unlist(strsplit(paste(PML.Genes[PML.ModuleNumber.pos_idx], collapse = ","), ","))))
    PML.ModuleGenes.pos <- mixedsort(PML.ModuleGenes.pos)
    
    PML.ModuleGenes.neg <- unique(trimws(unlist(strsplit(paste(PML.Genes[PML.ModuleNumber.neg_idx], collapse = ","), ","))))
    PML.ModuleGenes.neg <- mixedsort(PML.ModuleGenes.neg)
    
    # /// human ///
    PML.ModuleGenes.pos.human <- str_to_upper(PML.ModuleGenes.pos)
    PML.ModuleGenes.pos.human = list(PML.ModuleGenes.pos.human)
    names(PML.ModuleGenes.pos.human) = paste('PML','.','Module',PML.ModuleID[M.No],'.','Positive', sep = '')
    PML.ModuleGenes.pos.human.list <- c(PML.ModuleGenes.pos.human.list, PML.ModuleGenes.pos.human)
    
    PML.ModuleGenes.neg.human <- str_to_upper(PML.ModuleGenes.neg)
    PML.ModuleGenes.neg.human = list(PML.ModuleGenes.neg.human)
    names(PML.ModuleGenes.neg.human) = paste('PML','.','Module',PML.ModuleID[M.No],'.','Negative', sep = '')
    PML.ModuleGenes.neg.human.list <- c(PML.ModuleGenes.neg.human.list, PML.ModuleGenes.neg.human)
    
    PML.ModuleGenes.human <- unique(c(unlist(PML.ModuleGenes.pos.human), unlist(PML.ModuleGenes.neg.human)))
    PML.ModuleGenes.human <- mixedsort(PML.ModuleGenes.human)
    PML.ModuleGenes.human = list(PML.ModuleGenes.human)
    names(PML.ModuleGenes.human) = paste('PML','.','Module',PML.ModuleID[M.No], sep = '')
    PML.ModuleGenes.human.list <- c(PML.ModuleGenes.human.list, PML.ModuleGenes.human)
    
    # /// mouse ///
    PML.ModuleGenes.pos.mouse <- str_to_title(PML.ModuleGenes.pos)
    PML.ModuleGenes.pos.mouse = list(PML.ModuleGenes.pos.mouse)
    names(PML.ModuleGenes.pos.mouse) = paste('PML','.','Module',PML.ModuleID[M.No],'.','Positive', sep = '')
    PML.ModuleGenes.pos.mouse.list <- c(PML.ModuleGenes.pos.mouse.list, PML.ModuleGenes.pos.mouse)
    
    PML.ModuleGenes.neg.mouse <- str_to_title(PML.ModuleGenes.neg)
    PML.ModuleGenes.neg.mouse = list(PML.ModuleGenes.neg.mouse)
    names(PML.ModuleGenes.neg.mouse) = paste('PML','.','Module',PML.ModuleID[M.No],'.','Negative', sep = '')
    PML.ModuleGenes.neg.mouse.list <- c(PML.ModuleGenes.neg.mouse.list, PML.ModuleGenes.neg.mouse)
    
    PML.ModuleGenes.mouse <- unique(c(unlist(PML.ModuleGenes.pos.mouse), unlist(PML.ModuleGenes.neg.mouse)))
    PML.ModuleGenes.mouse <- mixedsort(PML.ModuleGenes.mouse)
    PML.ModuleGenes.mouse = list(PML.ModuleGenes.mouse)
    names(PML.ModuleGenes.mouse) = paste('PML','.','Module',PML.ModuleID[M.No], sep = '')
    PML.ModuleGenes.mouse.list <- c(PML.ModuleGenes.mouse.list, PML.ModuleGenes.mouse)
}

PML.ModuleGenes.pos.human.list.num <- sapply(PML.ModuleGenes.pos.human.list, FUN = function(X) length(X))
PML.ModuleGenes.neg.human.list.num <- sapply(PML.ModuleGenes.neg.human.list, FUN = function(X) length(X))
PML.ModuleGenes.human.list.num <- sapply(PML.ModuleGenes.human.list, FUN = function(X) length(X))
PML.ModuleGenes.pos.mouse.list.num <- sapply(PML.ModuleGenes.pos.mouse.list, FUN = function(X) length(X))
PML.ModuleGenes.neg.mouse.list.num <- sapply(PML.ModuleGenes.neg.mouse.list, FUN = function(X) length(X))
PML.ModuleGenes.mouse.list.num <- sapply(PML.ModuleGenes.mouse.list, FUN = function(X) length(X))

signature_genes = PML.ModuleGenes.mouse.list
signature_genes.idx <- unlist(sapply(c(paste('PML','.','Module',3, sep = ''),paste('PML','.','Module',4, sep = '')), FUN = function(X) which(names(signature_genes) == X)))

signature_genes = signature_genes[signature_genes.idx]
signature_genes.idx <- unlist(sapply(c(paste('PML','.','Module',3, sep = ''),paste('PML','.','Module',4, sep = '')), FUN = function(X) which(names(signature_genes) == X)))

names(signature_genes)[which(names(signature_genes) == paste('PML','.','Module',3, sep = ''))] = "Proliferative"
names(signature_genes)[which(names(signature_genes) == paste('PML','.','Module',4, sep = ''))] = "Inflammatory"

PML.signature_genes = signature_genes


# --- Analysis --------------------------------------------------------------------------------------------------------------------------------------
mainDir <- paste(current_folder,'/Result', sep = '')
subDir <- 'Mouse'
if (! file.exists(file.path(mainDir, subDir))) {
    dir.create(file.path(mainDir, subDir))
}

exp.name.all <- c()
cell_metadata.all <- c()

#for (s in 1:length(sample)) {
for (s in c(1,2)) {
    
    gc() # Garbage Collection
    
    mainDir <- paste(current_folder,'/Result','/','Mouse', sep = '')
    subDir <- 'library'
    if (! file.exists(file.path(mainDir, subDir))) {
        dir.create(file.path(mainDir, subDir))
    }
    
    mainDir <- paste(current_folder,'/Result','/','Mouse','/','library', sep = '')
    subDir <- sample[s]
    if (! file.exists(file.path(mainDir, subDir))) {
        dir.create(file.path(mainDir, subDir))
    }
    
    mainDir <- paste(current_folder,'/Result','/','Mouse','/','library','/',sample[s], sep = '')
    subDir <- 'raw'
    if (! file.exists(file.path(mainDir, subDir))) {
        dir.create(file.path(mainDir, subDir))
    }
    
    if (sample[s] == "control") {
        sample_folder = sample[s]
        file_num = length(control)
    } else if (sample[s] == "ch02") {
        #sample_folder = "20200910 Old dataset ch02"
        #sample_folder = "20200910 Old dataset full ch02"
        #sample_folder = "20201116_ch02_counts_new"
        sample_folder = "crispr"
        file_num = length(sgRNA.barcode)
    } else if (sample[s] == "d5_p1" | sample[s] == "d5_p2") {
        sample_folder = paste(sample[s],' reversed',' sgRNA count matrices', sep = '')
        file_num = length(sgRNA.barcode)
    } else {
        sample_folder = paste(sample[s],' sgRNA count matrices', sep = '')
        file_num = length(sgRNA.barcode)
    }
    
    cell.num.table <- c()
    cell.num.new.table <- c()
    
    gene.num.table <- c()
    gene.num.new.table <- c()
    
    cell.barcode.all <- c()
    cell.barcode.list <- list()
    
    df.all <- c()
    
    in.idx = 0
    
    for (i in 1:file_num) {
        
        gc() # Garbage Collection
        
        if (sample[s] == "control") {
            sample_file = control[i]
        } else {
            sample_file = sgRNA.barcode[i]
        }
        
        control.idx <- which(control %in% sample_file)
        sgRNA.barcode.idx <- which(sgRNA.barcode %in% sample_file)
        if (length(control.idx) == 1) {
            if (control.idx > 1) {
                presample_file <- control[control.idx - 1]
            }
        } else if (length(sgRNA.barcode.idx) == 1) {
            if (sgRNA.barcode.idx > 1) {
                presample_file <- sgRNA.barcode[sgRNA.barcode.idx - 1]
            }
        } else {
            presample_file <- NULL
        }
        
        if (sample_file %in% control.counts) {
            #cat('[control.counts] ',sample_file,'\n', sep = '')
            
            if (sample_file %in% c("acinar")) {
                sgRNA <- read.csv(paste('./Data','/','Mouse','/',sample_folder,'/','counts_csv','/','counts_',sample_file,'_final','.csv', sep = ''), sep = ',', header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
                
                # "counts_acinar_final.csv" loses "ensembl" column.
                cancer.sgRNA <- read.csv(paste('./Data','/','Mouse','/',sample_folder,'/','counts_csv','/','counts_','cancer','_final','.csv', sep = ''), sep = ',', header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
                
                NA.idx <- which(sgRNA[,"genes"] == "")
                sgRNA[NA.idx,"genes"] <- paste('Row',rownames(sgRNA)[NA.idx], sep = '')
                NA.idx <- which(sgRNA[,"gene_names"] == "")
                sgRNA[NA.idx,"gene_names"] <- paste('Row',rownames(sgRNA)[NA.idx], sep = '')
                rownames(sgRNA) <- sgRNA[,"genes"]
                
                row.idx <- grep('Row', rownames(sgRNA), value = FALSE)
                rownames(sgRNA)[row.idx] = rownames(cancer.sgRNA)[row.idx]
            } else if (sample_file %in% c("counts_kcppc_new")) {
                sgRNA <- read.csv(paste('./Data','/','Mouse','/',sample_folder,'/','counts_csv','/',sample_file,'.csv', sep = ''), sep = ',', header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
            } else if (sample_file %in% c(cancer_tripleMut.counts)) {
                sgRNA <- read.csv(paste('./Data','/','Mouse','/',sample_folder,'/','counts_csv','/','counts_',sample_file,'_final','.csv', sep = ''), sep = ',', header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
            } else if (sample_file %in% c(beta.counts, cancer_doubleMut.counts)) {
                sgRNA <- read.csv(paste('./Data','/','Mouse','/',sample_folder,'/','counts_csv','/','counts_',sample_file,'_final','.csv', sep = ''), sep = ',', header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
            }
        } else if (sample_file %in% control.old) {
            #cat('[control.old] ',sample_file,'\n', sep = '')
            
            if (sample_file %in% c(acinar.old, ductal.old, cancer_tripleMut.old)) {
                sgRNA <- read.csv(paste('./Data','/','Mouse','/',sample_folder,'/','old','/',sample_file,'.csv', sep = ''), sep = ',', header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
            } else if (sample_file %in% c(beta.old, cancer_doubleMut.old)) {
                sgRNA <- read.csv(paste('./Data','/','Mouse','/',sample_folder,'/','old','/','Chun-Hao 2021Q1','/','data_Chun-Hao 2021Q1_newest_analysis_results_',sample_file,'_filtered','.csv', sep = ''), sep = ',', header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
            }
        } else if (sample_file %in% control.new) {
            #cat('[control.new] ',sample_file,'\n', sep = '')
            
            if (sample_file %in% c(acinar.new, ductal.new, cancer_tripleMut.new)) {
                sgRNA <- read.csv(paste('./Data','/','Mouse','/',sample_folder,'/','new','/',sample_file,'.csv', sep = ''), sep = ',', header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
            } else if (sample_file %in% c(beta.new, cancer_doubleMut.new)) {
                sgRNA <- read.csv(paste('./Data','/','Mouse','/',sample_folder,'/','new','/','counts_new','/',sample_file,'_new','.csv', sep = ''), sep = ',', header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
            }
        } else if (sample_file %in% paper.cell) {
            #cat('[paper.cell] ',sample_file,'\n', sep = '')
            
            if (sample_file %in% c(paper.fibroblast)) {
                sgRNA <- read.csv(paste('./Data','/','NCBI_data','/','Tuveson_SRP191615','/','fibroblast','/',sample_file,'_counts_new','.csv', sep = ''), sep = ',', header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
            } else if (sample_file %in% c(paper.viable)) {
                sgRNA <- read.csv(paste('./Data','/','NCBI_data','/','Tuveson_SRP191615','/','viable','/',sample_file,'_counts_new','.csv', sep = ''), sep = ',', header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
            }
        } else {
            if (sample[s] == "ch02") {
                if (sample_folder == "20200910 Old dataset ch02" | sample_folder == "20200910 Old dataset full ch02") {
                    sgRNA <- read.csv(paste('./Data','/','Mouse','/',sample_folder,'/',sample_file,'.csv', sep = ''), sep = ',', header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
                } else if (sample_folder == "20201116_ch02_counts_new") {
                    sgRNA <- read.csv(paste('./Result','/','Mouse','/',sample_folder,'/','ch02_20201116','/',sample_file,'.csv', sep = ''), sep = ',', header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
                } else if (sample_folder == "crispr") {
                    sgRNA <- read.csv(paste('./Result','/','Mouse','/',sample_folder,'/','ch02','/',sample_file,'.csv', sep = ''), sep = ',', header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
                }
            } else {
                sgRNA <- read.csv(paste('./Data','/','Mouse','/',sample_folder,'/',sample_file,'.csv', sep = ''), sep = ',', header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
            }
        }
        #print(sample_file); print(dim(sgRNA))
        
        genes.idx <- which(colnames(sgRNA) == "genes")
        if (length(genes.idx) == 0) {
            genes.NA = 1
            sgRNA[,"genes"] <- rownames(sgRNA)
        } else {
            genes.NA = 0
            NA.idx <- sort(unique(c(which(sgRNA[,"genes"] == ""), grep('Row', sgRNA[,"genes"], value = FALSE))))
            sgRNA[NA.idx,"genes"] <- rownames(sgRNA)[NA.idx]
        }
        genes.idx <- which(colnames(sgRNA) == "genes")
        if (genes.idx != 1) {
            col.idx <- c(genes.idx, 1:(ncol(sgRNA)-1))
            sgRNA <- sgRNA[, col.idx]
        }
        genes.idx <- which(colnames(sgRNA) == "genes")
        
        gene_names.idx <- which(colnames(sgRNA) == "gene_names")
        if (length(gene_names.idx) == 0) {
            gene_names.NA = 1
            sgRNA[,"gene_names"] <- rownames(sgRNA)
        } else {
            gene_names.NA = 0
            NA.idx <- sort(unique(c(which(sgRNA[,"gene_names"] == ""), grep('Row', sgRNA[,"gene_names"], value = FALSE))))
            sgRNA[NA.idx,"gene_names"] <- rownames(sgRNA)[NA.idx]
        }
        gene_names.idx <- which(colnames(sgRNA) == "gene_names")
        if (gene_names.idx != 2) {
            col.idx <- c(1, gene_names.idx, 2:(ncol(sgRNA)-1))
            sgRNA <- sgRNA[, col.idx]
        }
        gene_names.idx <- which(colnames(sgRNA) == "gene_names")
        
        # /// duplicated cell barcode in the same file ///
        dup.idx <- which(duplicated(colnames(sgRNA)))
        dup.name <- unique(colnames(sgRNA)[dup.idx])
        if (length(dup.name) > 0) {
            if (in.idx == 0) {
                cat('[',sample[s],'] Duplicated cell barcode in "',sample[s],'" dataset: ','\n', sep = '')
                in.idx = 1
            }
            cat('Cell barcodes are duplicated in "',paste(sample_file,'.csv', sep = ''),'". ', sep = '')
            if (length(dup.name) == 1) {
                cat('Duplicated cell barcode: ',paste(dup.name, collapse = ', '),'\n', sep = '')
            } else if (length(dup.name) > 1) {
                cat('Duplicated cell barcodes: ',paste(dup.name, collapse = ', '),'\n', sep = '')
            }
            
            for (d in 1:length(dup.name)) {
                col.idx <- which(colnames(sgRNA) == dup.name[d])
                sgRNA[, col.idx] = rowSums(sgRNA[, col.idx], na.rm = TRUE)
            }
            sgRNA <- sgRNA[, -dup.idx, drop = FALSE]
        }
        
        cell.num = ncol(sgRNA) - 2
        names(cell.num) <- sample_file
        cell.num.table <- c(cell.num.table, cell.num)
        
        gene.num = nrow(sgRNA)
        names(gene.num) <- sample_file
        gene.num.table <- c(gene.num.table, gene.num)
        
        if (cell.num > 0) {
            
            cell.barcode <- trimws(colnames(sgRNA)[3:ncol(sgRNA)])
            
            if (length(grep('-1$', cell.barcode, value = TRUE)) > 0) { # "kcppc_bam_counts" # old ductal cell (new processing pipeline, without gene_names)
                col.idx <- grep('-1$', cell.barcode, value = FALSE)
                cell.barcode[col.idx] <- gsub('-1$', '', cell.barcode[col.idx])
            }
            
            cell.barcode.all.name <- sapply(cell.barcode.all, FUN = function(X) unlist(strsplit(X, '\\.'))[[1]])
            
            
            # --- QC 1. Check duplicated cell barcode -----------------------------------------------------------------------------------------------------------
            cell.barcode.all <- c(cell.barcode.all, cell.barcode)
            
            cell.barcode.list[[i]] <- cell.barcode
            names(cell.barcode.list)[i] <- sample_file
            
            
            # --- Data curation ---------------------------------------------------------------------------------------------------------------------------------
            # /// cell: append sample name as suffix to duplicated cells ///
            cell.barcode_dup.idx <- which(cell.barcode %in% cell.barcode.all.name)
            cell.barcode.all_dup.idx <- which(cell.barcode.all.name %in% cell.barcode)
            
            if (length(cell.barcode.all_dup.idx) > 0) {
                for (d in 1:length(cell.barcode.all_dup.idx)) {
                    
                    #print(d)
                    #print(cell.barcode.all_dup.idx[d])
                    #print(cell.barcode.all[cell.barcode.all_dup.idx[d]])
                    
                    if (length(grep('\\.', cell.barcode.all[cell.barcode.all_dup.idx[d]], value = TRUE)) == 0) {
                        #cell.barcode.all[cell.barcode.all_dup.idx[d]] <- paste(cell.barcode.all[cell.barcode.all_dup.idx[d]], '.', presample_file, sep = '')
                        
                        if (length(control.idx) == 1) {
                            if (control.idx > 1) {
                                presample <- control
                                presample.idx <- control.idx - 1
                            }
                        } else if (length(sgRNA.barcode.idx) == 1) {
                            if (sgRNA.barcode.idx > 1) {
                                presample <- sgRNA.barcode
                                presample.idx <- sgRNA.barcode.idx - 1
                            }
                        } else {
                            presample <- NULL
                            presample.idx <- NULL
                        }
                        
                        col.num = 0
                        for (k in 1:presample.idx) {
                            presample_file <- presample[k]
                            
                            presample_file.df <- eval(parse(text=paste('sample','.',sample[s],'.',presample_file,'.','df', sep = '')))
                            presample_file.df.colname <- sapply(colnames(presample_file.df), FUN = function(X) unlist(strsplit(X, '\\.'))[[1]])
                            col.idx <- which(presample_file.df.colname == cell.barcode.all[cell.barcode.all_dup.idx[d]])
                            if (length(col.idx) > 0) {
                                if (length(grep('\\.', colnames(presample_file.df)[col.idx], value = TRUE)) == 0) {
                                    colnames(presample_file.df)[col.idx] <- paste(colnames(presample_file.df)[col.idx], '.', presample_file, sep = '')
                                }
                            }
                            assign(paste('sample','.',sample[s],'.',presample_file,'.','df', sep = ''), presample_file.df)
                            
                            presample_file.df_new <- eval(parse(text=paste('sample','.',sample[s],'.',presample_file,'.','df_new', sep = '')))
                            presample_file.df_new.colname <- sapply(colnames(presample_file.df_new), FUN = function(X) unlist(strsplit(X, '\\.'))[[1]])
                            col.idx <- which(presample_file.df_new.colname == cell.barcode.all[cell.barcode.all_dup.idx[d]])
                            if (length(col.idx) > 0) {
                                if (length(grep('\\.', colnames(presample_file.df_new)[col.idx], value = TRUE)) == 0) {
                                    colnames(presample_file.df_new)[col.idx] <- paste(colnames(presample_file.df_new)[col.idx], '.', presample_file, sep = '')
                                }
                            }
                            assign(paste('sample','.',sample[s],'.',presample_file,'.','df_new', sep = ''), presample_file.df_new)
                            col.num = col.num + ncol(presample_file.df_new)
                            
                            #presample_file.df_new_normalized <- eval(parse(text=paste('sample','.',sample[s],'.',presample_file,'.','df_new_normalized', sep = '')))
                            #presample_file.df_new_normalized.colname <- sapply(colnames(presample_file.df_new_normalized), FUN = function(X) unlist(strsplit(X, '\\.'))[[1]])
                            #col.idx <- which(presample_file.df_new_normalized.colname == cell.barcode.all[cell.barcode.all_dup.idx[d]])
                            #if (length(col.idx) > 0) {
                            #    if (length(grep('\\.', colnames(presample_file.df_new_normalized)[col.idx], value = TRUE)) == 0) {
                            #        colnames(presample_file.df_new_normalized)[col.idx] <- paste(colnames(presample_file.df_new_normalized)[col.idx], '.', presample_file, sep = '')
                            #    }
                            #}
                            #assign(paste('sample','.',sample[s],'.',presample_file,'.','df_new_normalized', sep = ''), presample_file.df_new_normalized)
                            
                            presample_file.df_all <- eval(parse(text=paste('sample','.',sample[s],'.','df_all', sep = '')))
                            presample_file.df_all.colname <- sapply(colnames(presample_file.df_all), FUN = function(X) unlist(strsplit(X, '\\.'))[[1]])
                            col.idx <- which(presample_file.df_all.colname == cell.barcode.all[cell.barcode.all_dup.idx[d]])
                            if (length(col.idx) > 0) {
                                if (col.idx <= col.num) {
                                    if (length(grep('\\.', colnames(presample_file.df_all)[col.idx], value = TRUE)) == 0) {
                                        #print(colnames(presample_file.df_all)[col.idx])
                                        colnames(presample_file.df_all)[col.idx] <- paste(colnames(presample_file.df_all)[col.idx], '.', presample_file, sep = '')
                                        #print(colnames(presample_file.df_all)[col.idx])
                                    }
                                    
                                    if (length(cell.barcode_dup.idx) > 0) {
                                        cell.barcode2cell.barcode.all_idx <- which(cell.barcode[cell.barcode_dup.idx] %in% cell.barcode.all[cell.barcode.all_dup.idx[d]])
                                        
                                        if (length(cell.barcode2cell.barcode.all_idx) > 0) {
                                            for (dd in 1:length(cell.barcode2cell.barcode.all_idx)) {
                                                if (length(grep('\\.', cell.barcode[cell.barcode_dup.idx[cell.barcode2cell.barcode.all_idx[dd]]], value = TRUE)) == 0) {
                                                    #print(cell.barcode[cell.barcode_dup.idx[cell.barcode2cell.barcode.all_idx[dd]]])
                                                    cell.barcode[cell.barcode_dup.idx[cell.barcode2cell.barcode.all_idx[dd]]] <- paste(cell.barcode[cell.barcode_dup.idx[cell.barcode2cell.barcode.all_idx[dd]]], '.', sample_file, sep = '')
                                                    #print(cell.barcode[cell.barcode_dup.idx[cell.barcode2cell.barcode.all_idx[dd]]])
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                            assign(paste('sample','.',sample[s],'.','df_all', sep = ''), presample_file.df_all)
                            assign("df.all", presample_file.df_all)
                        }
                    }
                }
            }
            
            # /// gene: replace empty gene_names with ID and append ID to duplicated gene_names ///
            df <- sgRNA[, 3:ncol(sgRNA), drop = FALSE]
            
            colnames(df) <- cell.barcode
            
            gene.dup.idx <- which(duplicated(sgRNA[,"gene_names"]))
            if (length(gene.dup.idx) > 0) {
                for (d in 1:length(gene.dup.idx)) {
                    gene_name <- sgRNA[which(sgRNA[,"gene_names"] == sgRNA[gene.dup.idx[d],"gene_names"]),"gene_names"]
                    gene_name <- unique(gene_name)
                    
                    if (length(grep(')$', sgRNA[which(sgRNA[,"gene_names"] == gene_name),"gene_names"], value = TRUE)) == 0) {
                        sgRNA[which(sgRNA[,"gene_names"] == gene_name),"gene_names"] <- paste(sgRNA[which(sgRNA[,"gene_names"] == gene_name),"gene_names"], "(", sgRNA[which(sgRNA[,"gene_names"] == gene_name),"genes"], ")", sep = '')
                    }
                }
            }
            
            all_gene_id <- sgRNA[,"genes"]
            
            rownames(df) <- sgRNA[,"gene_names"]
            all_gene_name <- rownames(df)
            
            if (sample[s] == "control") {
                control_gene_id <- all_gene_id
                control_gene_name <- all_gene_name
            }
            
            write.csv(df, file = paste(current_folder,'/Result','/','Mouse','/','library','/',sample[s],'/','raw','/',sample[s],'.',sample_file,'.csv', sep = ''), row.names = TRUE, quote = TRUE)
            
            
            # --- Data filtering --------------------------------------------------------------------------------------------------------------------------------
            for (g in 1:length(gene_cutoff)) {
                df_qc <- qc_data(df, normalize = FALSE, min_gene_frac = gene_cutoff[g])
                df_qc.normalized <- qc_data(df, normalize = TRUE, min_gene_frac = gene_cutoff[g])
                
                df.new <- df_qc
                df.new.normalized <- df_qc.normalized
                
                if (! is.null(df.new)) {
                    cell.num.new = ncol(df.new)
                    gene.num.new = nrow(df.new)
                } else {
                    cell.num.new = 0
                    gene.num.new = 0
                }
                
                assign(paste('sample','.',sample[s],'.',sample_file,'.',gene_cutoff[g],'.','cell_num_new', sep = ''), cell.num.new)
                assign(paste('sample','.',sample[s],'.',sample_file,'.',gene_cutoff[g],'.','gene_num_new', sep = ''), gene.num.new)
            }
            
            df_qc <- qc_data(df, normalize = FALSE, min_gene_frac = gene_cutoff.default)
            df_qc.normalized <- qc_data(df, normalize = TRUE, min_gene_frac = gene_cutoff.default)
            
            df.new <- df_qc
            df.new.normalized <- df_qc.normalized
            
            assign(paste('sample','.',sample[s],'.',sample_file,'.','df', sep = ''), df)
            assign(paste('sample','.',sample[s],'.',sample_file,'.','df_new', sep = ''), df.new)
            #assign(paste('sample','.',sample[s],'.',sample_file,'.','df_new_normalized', sep = ''), df.new.normalized)
            
            if (! is.null(df.new)) {
                cell.num.new = ncol(df.new)
                gene.num.new = nrow(df.new)
                
                df.all <- merge(df.all, df.new, by = "row.names", all = TRUE, sort = TRUE)
                rownames(df.all) <- df.all$Row.names
                df.all$Row.names <- NULL
                
                assign(paste('sample','.',sample[s],'.','df_all', sep = ''), df.all)
            } else {
                cell.num.new = 0
                gene.num.new = 0
            }
        } else {
            
            cell.barcode.list[[i]] <- NA
            names(cell.barcode.list)[i] <- sample_file
            
            assign(paste('sample','.',sample[s],'.',sample_file,'.','df', sep = ''), NULL)
            assign(paste('sample','.',sample[s],'.',sample_file,'.','df_new', sep = ''), NULL)
            #assign(paste('sample','.',sample[s],'.',sample_file,'.','df_new_normalized', sep = ''), NULL)
            
            for (g in 1:length(gene_cutoff)) {
                cell.num.new = 0
                gene.num.new = 0
                
                assign(paste('sample','.',sample[s],'.',sample_file,'.',gene_cutoff[g],'.','cell_num_new', sep = ''), cell.num.new)
                assign(paste('sample','.',sample[s],'.',sample_file,'.',gene_cutoff[g],'.','gene_num_new', sep = ''), gene.num.new)
            }
            
            cell.num.new = 0
            gene.num.new = 0
            
        }
        
        names(cell.num.new) <- sample_file
        cell.num.new.table <- c(cell.num.new.table, cell.num.new)
        
        names(gene.num.new) <- sample_file
        gene.num.new.table <- c(gene.num.new.table, gene.num.new)
        
    }
    cat('\n')
    
    #df.all[is.na(df.all)] <- 0
    
    row.idx <- unlist(sapply(all_gene_name, FUN = function(X) which(rownames(df.all) %in% X)))
    df.all <- df.all[row.idx, , drop = FALSE]
    
    assign(paste('sample','.',sample[s],'.','df_all', sep = ''), df.all)
    
    
    RData.name = unlist(strsplit(current_folder, '/'))
    RData.name = RData.name[length(RData.name)]
    save.image(file = paste('./../',RData.name,'_','part1','.RData', sep = ''))
    
    
    # --- QC 1. Check duplicated cell barcode -----------------------------------------------------------------------------------------------------------
    #print(length(cell.barcode.all))
    #print(cell.barcode.list)
    
    # /// duplicated cell barcode in different files ///
    dup.idx <- which(duplicated(cell.barcode.all))
    if (length(dup.idx) > 0) {
        cell.barcode.dup <- cell.barcode.all[dup.idx]
        
        cell.barcode.dup.idx <- which(duplicated(cell.barcode.dup))
        if (length(cell.barcode.dup.idx) > 0) {
            rm.idx.all <- c()
            for (d in 1:length(cell.barcode.dup.idx)) {
                rm.idx <- which(cell.barcode.dup %in% cell.barcode.dup[cell.barcode.dup.idx[d]])
                rm.idx <- rm.idx[-which.max(rm.idx)]
                rm.idx.all <- c(rm.idx.all, rm.idx)
            }
            rm.idx.all <- sort(unique(rm.idx.all))
            dup.idx <- dup.idx[-rm.idx.all]
        }
        
        if (length(dup.idx) > 0) {
            if (in.idx == 0) {
                cat('[',sample[s],'] Duplicated cell barcode in "',sample[s],'" dataset: ','\n', sep = '')
                in.idx = 1
            }
            for (d in 1:length(dup.idx)) {
                name.idx <- which(sapply(cell.barcode.list, FUN = function(X) cell.barcode.all[dup.idx[d]] %in% X))
                cat('"',cell.barcode.all[dup.idx[d]],'" was found in "',paste(paste(unique(names(name.idx)),'.csv', sep = ''), collapse = '";"'),'".','\n', sep = '')
            }
            cat('\n')
        }
    }
    
    assign(paste('sample','.',sample[s],'.','cell_num_table', sep = ''), cell.num.table)
    assign(paste('sample','.',sample[s],'.','cell_num_new_table', sep = ''), cell.num.new.table)
    
    cat('[',sample[s],'] Cell number before QC: ',sum(cell.num.table),'\n', sep = '')
    print(cell.num.table)
    cat('\n')
    cat('[',sample[s],'] Cell number after QC: ',sum(cell.num.new.table),'\n', sep = '')
    print(cell.num.new.table)
    cat('\n')
    #cat('\n')
    
    assign(paste('sample','.',sample[s],'.','gene_num_table', sep = ''), gene.num.table)
    assign(paste('sample','.',sample[s],'.','gene_num_new_table', sep = ''), gene.num.new.table)
    
    cat('[',sample[s],'] Gene number before QC: ',sum(gene.num.table),'\n', sep = '')
    print(gene.num.table)
    cat('\n')
    cat('[',sample[s],'] Gene number after QC: ',sum(gene.num.new.table),'\n', sep = '')
    print(gene.num.new.table)
    cat('\n')
    cat('\n')
    
    
    # --- cell_metadata ---------------------------------------------------------------------------------------------------------------------------------
    if (sample[s] == "control") {
        exp.name <- control
    } else {
        exp.name <- sgRNA.barcode
    }
    exp.name.all <- c(exp.name.all, exp.name)
    exp.name.all <- unique(exp.name.all)
    
    exp.idx <- lapply(exp.name, FUN = function(X) which(colnames(df.all) %in% colnames(eval(parse(text=paste('sample','.',sample[s],'.',X,'.','df_new', sep = ''))))))
    names(exp.idx) <- exp.name
    
    exp.num <- sapply(exp.idx, FUN = function(X) length(X))
    
    cell_metadata <- data.frame("cell_barcode" = colnames(df.all)[unlist(exp.idx)], 
                                "sample" = rep(sample[s], ncol(df.all)), 
                                "sgRNA" = rep(names(exp.num), exp.num), 
                                stringsAsFactors = FALSE)
    
    cell_metadata[,"sample"] <- factor(cell_metadata[,"sample"], levels = unique(cell_metadata[,"sample"]))
    cell_metadata[,"sgRNA"] <- factor(cell_metadata[,"sgRNA"], levels = unique(cell_metadata[,"sgRNA"]))
    
    cell_metadata[,"Condition"] <- paste(cell_metadata[,"sample"], cell_metadata[,"sgRNA"], sep = '.')
    cell_metadata[,"Condition"] <- factor(cell_metadata[,"Condition"], levels = unique(cell_metadata[,"Condition"]))
    
    cell_metadata[,"ID"] <- cell_metadata[,"cell_barcode"]
    rownames(cell_metadata) <- cell_metadata[,"cell_barcode"]
    
    cell_metadata[,"cell_barcode"] <- sapply(cell_metadata[,"cell_barcode"], FUN = function(X) unlist(strsplit(X, '\\.'))[[1]])
    
    #assign(paste('sample','.',sample[s],'.','cell_metadata', sep = ''), cell_metadata)
    
    cell_metadata.all <- rbind(cell_metadata.all, cell_metadata)
    
    dup.idx <- which(duplicated(cell_metadata.all[,"ID"]))
    if (length(dup.idx) > 0) {
        for (d in 1:length(dup.idx)) {
            row.idx <- which(cell_metadata.all[,"ID"] == cell_metadata.all[dup.idx[d],"ID"])
            
            if (length(grep('\\.', cell_metadata.all[row.idx,"ID"], value = TRUE)) == 0) {
                cell_metadata.all[row.idx,"ID"] <- paste(cell_metadata.all[row.idx,"ID"], '.', cell_metadata.all[row.idx,"sample"], sep = '')
            }
        }
    }
    
    rownames(cell_metadata.all) <- cell_metadata.all[,"ID"]
    cell_metadata.all <- unique(cell_metadata.all)
    
    
    # --- Data correction and combination ---------------------------------------------------------------------------------------------------------------
    # /// replace NA with true value ///
    mainDir <- paste(current_folder,'/Result','/','Mouse','/','library','/',sample[s], sep = '')
    subDir <- 'QC'
    if (! file.exists(file.path(mainDir, subDir))) {
        dir.create(file.path(mainDir, subDir))
    }
    
    initial.idx = 1
    initial2end.idx <- c()
    df.subset.all <- c()
    for (i in 1:length(exp.name)) {
        df <- eval(parse(text=paste('sample','.',sample[s],'.',exp.name[i],'.','df', sep = '')))
        df.origin <- df
        
        df.all.subset <- df.all[, exp.idx[[i]], drop = FALSE]
        
        if (! is.null(df)) {
            colnames(df) = sapply(colnames(df), FUN = function(X) unlist(strsplit(X, '\\.'))[[1]])
            colnames(df.all.subset) = sapply(colnames(df.all.subset), FUN = function(X) unlist(strsplit(X, '\\.'))[[1]])
            
            row.idx <- unlist(sapply(rownames(df.all.subset), FUN = function(X) which(rownames(df) %in% X)))
            col.idx <- unlist(sapply(colnames(df.all.subset), FUN = function(X) which(colnames(df) %in% X)))
            
            df.subset <- df.origin[row.idx, col.idx, drop = FALSE]
            
            #assign(paste('sample','.',sample[s],'.',exp.name[i],'.','df_subset', sep = ''), df.subset)
            #write.csv(df.subset, file = paste(current_folder,'/Result','/','Mouse','/','library','/',sample[s],'/','QC','/',sample[s],'.QC','.',exp.name[i],'.csv', sep = ''), row.names = TRUE, quote = TRUE)
            
            initial2end.idx <- c(initial2end.idx, list(seq(initial.idx, length.out = ncol(df.subset))))
            names(initial2end.idx)[length(initial2end.idx)] <- exp.name[i]
            initial.idx = initial.idx + ncol(df.subset)
            
            if (! is.null(df.subset.all)) {
                df.subset.all.colnames = sapply(colnames(df.subset.all), FUN = function(X) unlist(strsplit(X, '\\.'))[[1]])
                df.subset.colnames = sapply(colnames(df.subset), FUN = function(X) unlist(strsplit(X, '\\.'))[[1]])
                
                df.subset.all.dup.idx <- which(sapply(df.subset.all.colnames, FUN = function(X) X %in% df.subset.colnames))
                df.subset.dup.idx <- which(sapply(df.subset.colnames, FUN = function(X) X %in% df.subset.all.colnames))
                if (length(df.subset.all.dup.idx) > 0) {
                    for (d in 1:length(df.subset.all.dup.idx)) {
                        if (length(grep('\\.', colnames(df.subset.all)[df.subset.all.dup.idx[d]], value = TRUE)) == 0) {
                            sample.name <- names(which(sapply(initial2end.idx, FUN = function(X) df.subset.all.dup.idx[d] %in% unlist(X))))
                            colnames(df.subset.all)[df.subset.all.dup.idx[d]] <- paste(colnames(df.subset.all)[df.subset.all.dup.idx[d]], '.', sample.name, sep = '')
                        }
                    }
                }
                if (length(df.subset.dup.idx) > 0) {
                    for (d in 1:length(df.subset.dup.idx)) {
                        if (length(grep('\\.', colnames(df.subset)[df.subset.dup.idx[d]], value = TRUE)) == 0) {
                            sample.name <- exp.name[i]
                            colnames(df.subset)[df.subset.dup.idx[d]] <- paste(colnames(df.subset)[df.subset.dup.idx[d]], '.', sample.name, sep = '')
                        }
                    }
                }
            }
            
            df.subset.all <- merge(df.subset.all, df.subset, by = "row.names", all = TRUE, sort = TRUE)
            rownames(df.subset.all) <- df.subset.all$Row.names
            df.subset.all$Row.names <- NULL
        } else {
            initial2end.idx <- c(initial2end.idx, NA)
            names(initial2end.idx)[length(initial2end.idx)] <- exp.name[i]
        }
    }
    
    if (! is.null(df.subset.all)) {
        df.subset.all[is.na(df.subset.all)] <- 0
        
        row.idx <- unlist(sapply(all_gene_name, FUN = function(X) which(rownames(df.subset.all) %in% X)))
        df.subset.all <- df.subset.all[row.idx, , drop = FALSE]
        
        assign(paste('sample','.',sample[s],'.','df_subset_all', sep = ''), df.subset.all)
        write.csv(df.subset.all, file = paste(current_folder,'/Result','/','Mouse','/','library','/',sample[s],'/',sample[s],'.QC','.csv', sep = ''), row.names = TRUE, quote = TRUE)
    }
    
    
    # --- gene_metadata ---------------------------------------------------------------------------------------------------------------------------------
    expression_matrix = as.matrix(df.subset.all)
    
    row.idx <- unlist(sapply(rownames(df.subset.all), FUN = function(X) which(all_gene_name %in% X)))
    all_gene_id.subset.all <- all_gene_id[row.idx]
    
    gene_metadata <- data.frame("gene_short_name" = rownames(df.subset.all), 
                                "gene_id" = all_gene_id.subset.all, 
                                stringsAsFactors = FALSE)
    rownames(gene_metadata) <- gene_metadata[,"gene_short_name"]
    
    
    RData.name = unlist(strsplit(current_folder, '/'))
    RData.name = RData.name[length(RData.name)]
    save.image(file = paste('./../',RData.name,'_','part2','.RData', sep = ''))
    
    
    # --- Monocle3 & Seurat -----------------------------------------------------------------------------------------------------------------------------
    if (sample[s] == "control" | sample[s] == "ch02") {
        
        subdata <- list()
        subdata.name <- c()
        for (comb in 1:2) {
            if (comb == 1) {
                subdata.counts <- c(acinar.counts, ductal.counts, cancer_tripleMut.counts)
                subdata.old <- c(acinar.old, ductal.old, cancer_tripleMut.old)
                subdata.new <- c(acinar.new, ductal.new, cancer_tripleMut.new)
            } else if (comb == 2) {
                subdata.counts <- c(ductal.counts, cancer_tripleMut.counts)
                subdata.old <- c(ductal.old, cancer_tripleMut.old)
                subdata.new <- c(ductal.new, cancer_tripleMut.new)
            }
            subdata <- c(subdata, list(subdata.counts))
            #subdata <- c(subdata, list(subdata.old))
            #subdata <- c(subdata, list(subdata.new))
            
            subdata.name <- c(subdata.name, paste('subset',comb, sep = ''))
        }
        names(subdata) = subdata.name
        
        #for (sb in -1:length(subdata)) {
        for (sb in -1:-1) {
            
            expression_matrix.new = expression_matrix
            gene_metadata.new = gene_metadata
            cell_metadata.new = cell_metadata
            
            project.name = unlist(strsplit(current_folder, '/'))
            project.name = project.name[length(project.name)]
            
            if (sb == 0) {         # all
                cells.input = control
                project.input = paste(project.name,'.','all', sep = '')
            } else if (sb == -1) { # all_combine_paper
                cells.input = control
                project.input = paste(project.name,'.','all_combine_paper', sep = '')
            } else if (sb > 0) {   # subset
                cells.input = sapply(unlist(subdata[sb]), FUN = function(X) control[which(control == X)])
                project.input = paste(project.name,'.',names(subdata[sb]), sep = '')
                
                col.idx <- sapply(unlist(subdata[sb]), FUN = function(X) which(cell_metadata.new[,"sgRNA"] == X))
                names(col.idx) = unlist(subdata[sb])
                col.idx <- unlist(col.idx)
                
                expression_matrix.new = expression_matrix.new[, col.idx, drop = FALSE]
                cell_metadata.new = cell_metadata.new[col.idx, , drop = FALSE]
                for (k in 1:ncol(cell_metadata.new)) {
                    if (is.factor(cell_metadata.new[,k])) {
                        cell_metadata.new[,k] <- factor(cell_metadata.new[,k], levels = unique(cell_metadata.new[,k]))
                    }
                }
                
                Zero.idx <- which(rowSums(expression_matrix.new) == 0)
                if (length(Zero.idx) > 0) {
                    row.idx <- which(rowSums(expression_matrix.new) != 0)
                    
                    expression_matrix.new = expression_matrix.new[row.idx, , drop = FALSE]
                    gene_metadata.new = gene_metadata.new[row.idx, , drop = FALSE]
                }
            }
            
            acinar.idx = which(cell_metadata.new[,"sgRNA"] %in% c("acinar","1_acinar","1_acinar_new"))
            ductal.idx = which(cell_metadata.new[,"sgRNA"] %in% c("counts_kcppc_new","1_ductal","kcppc_bam_counts_new"))
            cancer.idx = which(cell_metadata.new[,"sgRNA"] %in% c("cancer","1_cancer","1_cancer_new"))
            
            beta.idx = c()
            min6.idx = which(cell_metadata.new[,"sgRNA"] %in% "min6")
            if (length(min6.idx) > 0) {
                beta.idx = min6.idx
            }
            min6_1.idx = which(cell_metadata.new[,"sgRNA"] %in% c("16726_MIN6","16726_min6"))
            min6_2.idx = which(cell_metadata.new[,"sgRNA"] %in% c("16728_MIN6","16728_min6"))
            if (length(min6_1.idx) > 0) {
                beta.idx = min6_1.idx
            }
            if (length(min6_2.idx) > 0) {
                beta.idx = min6_2.idx
            }
            if (length(min6_1.idx) > 0 & length(min6_2.idx) > 0) {
                beta.idx = unique(c(min6_1.idx, min6_2.idx))
            }
            
            PCC8.idx = which(cell_metadata.new[,"sgRNA"] %in% "16728_PCC8")
            pcc8.idx = which(cell_metadata.new[,"sgRNA"] %in% c("pcc8","16728_pcc8"))
            
            PCC9.idx = which(cell_metadata.new[,"sgRNA"] %in% "16726_PCC9")
            pcc9.idx = which(cell_metadata.new[,"sgRNA"] %in% c("pcc9","16726_pcc9"))
            
            SRR8872326.idx = which(cell_metadata.new[,"sgRNA"] %in% "SRR8872326")
            SRR8872327.idx = which(cell_metadata.new[,"sgRNA"] %in% "SRR8872327")
            SRR8872328.idx = which(cell_metadata.new[,"sgRNA"] %in% "SRR8872328")
            SRR8872329.idx = which(cell_metadata.new[,"sgRNA"] %in% "SRR8872329")
            paper.fibroblast.idx = unlist(sapply(paper.fibroblast, FUN = function(X) which(cell_metadata.new[,"sgRNA"] %in% X)))
            paper.viable.idx = unlist(sapply(paper.viable, FUN = function(X) which(cell_metadata.new[,"sgRNA"] %in% X)))
            paper.cell.idx = unlist(sapply(paper.cell, FUN = function(X) which(cell_metadata.new[,"sgRNA"] %in% X)))
            
            sgRNA_all.idx = which(cell_metadata.new[,"sgRNA"] %in% sgRNA.barcode)
            for (b in 1:length(sgRNA.barcode)) {
                assign(paste('sgRNA','.',sgRNA.barcode[b],'.','idx', sep = ''), which(cell_metadata.new[,"sgRNA"] %in% sgRNA.barcode[b]))
            }
            
            if (all(control.counts %in% control)) {
                cell_metadata.new[acinar.idx,"batch"] = 1
                cell_metadata.new[ductal.idx,"batch"] = 1
                cell_metadata.new[beta.idx,"batch"] = 1
                cell_metadata.new[PCC8.idx,"batch"] = 1
                cell_metadata.new[pcc8.idx,"batch"] = 1
                cell_metadata.new[PCC9.idx,"batch"] = 1
                cell_metadata.new[pcc9.idx,"batch"] = 1
                cell_metadata.new[cancer.idx,"batch"] = 1
                cell_metadata.new[paper.fibroblast.idx,"batch"] = 3
                cell_metadata.new[paper.viable.idx,"batch"] = 2
                cell_metadata.new[sgRNA_all.idx,"batch"] = 4
            } else if (all(control.old %in% control)) {
                cell_metadata.new[acinar.idx,"batch"] = 1
                cell_metadata.new[ductal.idx,"batch"] = 1
                cell_metadata.new[beta.idx,"batch"] = 2
                cell_metadata.new[PCC8.idx,"batch"] = 2
                cell_metadata.new[pcc8.idx,"batch"] = 2
                cell_metadata.new[PCC9.idx,"batch"] = 2
                cell_metadata.new[pcc9.idx,"batch"] = 2
                cell_metadata.new[cancer.idx,"batch"] = 1
                cell_metadata.new[paper.fibroblast.idx,"batch"] = 4
                cell_metadata.new[paper.viable.idx,"batch"] = 3
                cell_metadata.new[sgRNA_all.idx,"batch"] = 5
            } else if (all(control.new %in% control)) {
                cell_metadata.new[acinar.idx,"batch"] = 1
                cell_metadata.new[ductal.idx,"batch"] = 1
                cell_metadata.new[beta.idx,"batch"] = 2
                cell_metadata.new[PCC8.idx,"batch"] = 2
                cell_metadata.new[pcc8.idx,"batch"] = 2
                cell_metadata.new[PCC9.idx,"batch"] = 2
                cell_metadata.new[pcc9.idx,"batch"] = 2
                cell_metadata.new[cancer.idx,"batch"] = 1
                cell_metadata.new[paper.fibroblast.idx,"batch"] = 4
                cell_metadata.new[paper.viable.idx,"batch"] = 3
                cell_metadata.new[sgRNA_all.idx,"batch"] = 5
            }
            cell_metadata.new[,"batch"] <- factor(cell_metadata.new[,"batch"], levels = sort(unique(cell_metadata.new[,"batch"])))
            
            cell_metadata.new[acinar.idx,"cell.type"] = "acinar"
            cell_metadata.new[ductal.idx,"cell.type"] = "ductal"
            cell_metadata.new[beta.idx,"cell.type"] = "beta"
            cell_metadata.new[PCC8.idx,"cell.type"] = "PCC8"
            cell_metadata.new[pcc8.idx,"cell.type"] = "pcc8"
            cell_metadata.new[PCC9.idx,"cell.type"] = "PCC9"
            cell_metadata.new[pcc9.idx,"cell.type"] = "pcc9"
            cell_metadata.new[cancer.idx,"cell.type"] = "cancer"
            cell_metadata.new[paper.fibroblast.idx,"cell.type"] = "fibroblast"  # Fibroblast-enriched fraction (DAPI-CD45-CD31-EPCAM-)
            cell_metadata.new[paper.viable.idx,"cell.type"] = "viable"          # Viable cell fraction (DAPI-)
            for (b in 1:length(sgRNA.barcode)) {
                sgRNA_single.idx <- eval(parse(text=paste('sgRNA','.',sgRNA.barcode[b],'.','idx', sep = '')))
                cell_metadata.new[sgRNA_single.idx,"cell.type"] = paste('sgRNA','_',b, sep = '')
            }
            cell_metadata.new[,"cell.type"] <- factor(cell_metadata.new[,"cell.type"], levels = unique(cell_metadata.new[,"cell.type"]))
            
            cell_metadata.new[acinar.idx,"cell.type.color"] = "#E77D72"
            cell_metadata.new[ductal.idx,"cell.type.color"] = "#B39F33"
            cell_metadata.new[beta.idx,"cell.type.color"] = "#53B74C"
            cell_metadata.new[PCC8.idx,"cell.type.color"] = "#00FDFF"
            cell_metadata.new[pcc8.idx,"cell.type.color"] = "#00FDFF"
            cell_metadata.new[PCC9.idx,"cell.type.color"] = "#6E9AF8"
            cell_metadata.new[pcc9.idx,"cell.type.color"] = "#6E9AF8"
            cell_metadata.new[cancer.idx,"cell.type.color"] = "#E46DDD"
            cell_metadata.new[paper.fibroblast.idx,"cell.type.color"] = "#53B49C"
            cell_metadata.new[paper.viable.idx,"cell.type.color"] = "#D73A35"
            cell_metadata.new[sgRNA_all.idx,"cell.type.color"] = "#A6A6A6"
            for (b in 1:length(sgRNA.barcode)) {
                sgRNA_single.idx <- eval(parse(text=paste('sgRNA','.',sgRNA.barcode[b],'.','idx', sep = '')))
                if (sgRNA.barcode[b] == "39_CGGATGCTAATGCACGTG") { # toxic sgRNA
                    cell_metadata.new[sgRNA_single.idx,"cell.type.color"] = "#D8A584"
                } else if (sgRNA.barcode[b] == "64_CCCAGGGAGTTATTCGAT") { # AB2
                    cell_metadata.new[sgRNA_single.idx,"cell.type.color"] = "#FFC000"
                } else if (sgRNA.barcode[b] == "65_ATCCAAGTACCACCCCGA") { # AB2
                    cell_metadata.new[sgRNA_single.idx,"cell.type.color"] = "#FFFF00"
                } else if (sgRNA.barcode[b] == "66_TAAAATACAATCCCTCGG") { # toxic sgRNA (strong toxic)
                    cell_metadata.new[sgRNA_single.idx,"cell.type.color"] = "#AF5AFF"
                } else if (sgRNA.barcode[b] == "67_GATGTATCGGAGTAGTTG") { # toxic sgRNA (strong toxic)
                    cell_metadata.new[sgRNA_single.idx,"cell.type.color"] = "#863CBD"
                } else if (sgRNA.barcode[b] == "68_CCACAAGAATACATGCGT") { # toxic sgRNA (strong toxic)
                    cell_metadata.new[sgRNA_single.idx,"cell.type.color"] = "#5B2A82"
                } else if (sgRNA.barcode[b] == "69_GGAGAATCGAGATGGTGG") { # toxic sgRNA (strong toxic)
                    cell_metadata.new[sgRNA_single.idx,"cell.type.color"] = "#35184C"
                } else if (sgRNA.barcode[b] == "71_GTGCACACTAACAATGAG") { # negative control sgRNA 1
                    cell_metadata.new[sgRNA_single.idx,"cell.type.color"] = "#54667E"
                } else if (sgRNA.barcode[b] == "72_GGCTAGAATGTGTACCAT") { # negative control sgRNA 2
                    cell_metadata.new[sgRNA_single.idx,"cell.type.color"] = "#6D9385"
                }
            }
            cell_metadata.new[,"cell.type.color"] <- factor(cell_metadata.new[,"cell.type.color"], levels = unique(cell_metadata.new[,"cell.type.color"]))
            
            if (all(paper.cell %in% control)) {
                cell_metadata.new[acinar.idx,"stage"] = 1
                cell_metadata.new[ductal.idx,"stage"] = 2
                cell_metadata.new[beta.idx,"stage"] = 1
                
                cell_metadata.new[PCC8.idx,"stage"] = 4
                cell_metadata.new[pcc8.idx,"stage"] = 4
                cell_metadata.new[PCC9.idx,"stage"] = 4
                cell_metadata.new[pcc9.idx,"stage"] = 4
                
                cell_metadata.new[cancer.idx,"stage"] = 7
                
                cell_metadata.new[paper.fibroblast.idx,"stage"] = 3
                cell_metadata.new[paper.viable.idx,"stage"] = 5
                
                cell_metadata.new[sgRNA_all.idx,"stage"] = 6
            } else {
                cell_metadata.new[acinar.idx,"stage"] = 1
                cell_metadata.new[ductal.idx,"stage"] = 2
                cell_metadata.new[beta.idx,"stage"] = 1
                
                cell_metadata.new[PCC8.idx,"stage"] = 3
                cell_metadata.new[pcc8.idx,"stage"] = 3
                cell_metadata.new[PCC9.idx,"stage"] = 3
                cell_metadata.new[pcc9.idx,"stage"] = 3
                
                cell_metadata.new[cancer.idx,"stage"] = 5
                
                cell_metadata.new[sgRNA_all.idx,"stage"] = 4
            }
            cell_metadata.new[,"stage"] <- factor(cell_metadata.new[,"stage"], levels = sort(unique(cell_metadata.new[,"stage"])))
            
            # /// Monocle3 ///
            # *** no genes of interest and signature genes ***
            Monocle3_analysis(expression_matrix = expression_matrix.new, 
                              cell_metadata = cell_metadata.new, 
                              gene_metadata = gene_metadata.new, 
                              batch = "batch", 
                              cell_type = "cell.type", 
                              time = "stage", 
                              cells = cells.input, 
                              project = project.input, 
                              analysis_type = "trajectories", 
                              redDim_method = "UMAP", 
                              genes_of_interest = NULL, 
                              signature_genes = NULL)
            
            # *** KRAS ***
            Monocle3_analysis(expression_matrix = expression_matrix.new, 
                              cell_metadata = cell_metadata.new, 
                              gene_metadata = gene_metadata.new, 
                              batch = "batch", 
                              cell_type = "cell.type", 
                              time = "stage", 
                              cells = cells.input, 
                              project = project.input, 
                              analysis_type = "trajectories", 
                              redDim_method = "UMAP", 
                              genes_of_interest = KRAS.genes_of_interest, 
                              signature_genes = KRAS.signature_genes)
            
            # *** PML ***
            Monocle3_analysis(expression_matrix = expression_matrix.new, 
                              cell_metadata = cell_metadata.new, 
                              gene_metadata = gene_metadata.new, 
                              batch = "batch", 
                              cell_type = "cell.type", 
                              time = "stage", 
                              cells = cells.input, 
                              project = project.input, 
                              analysis_type = "trajectories", 
                              redDim_method = "UMAP", 
                              genes_of_interest = NULL, 
                              signature_genes = PML.signature_genes)
            
            # /// Seurat ///
            min.cells = 0 # min.cells: Include features detected in at least this many cells. Will subset the counts matrix as well. To reintroduce excluded features, create a new object with a lower cutoff.
            # To filter out low-quality cells, we first removed cells for which less than 500 genes were detected.
            min.features = 500 # min.features: Include cells where at least this many features are detected.
            
            nFeature_RNA.lower_cutoff = min.features
            nFeature_RNA.higher_cutoff = 2500
            # To filter out low-quality cells, we removed cells for which over 10% genes derived from mitochondrial genome.
            percent.mt.cutoff = 10
            
            Seurat.obj <- Seurat_analysis(expression_matrix = expression_matrix.new, 
                                          cell_metadata = cell_metadata.new, 
                                          gene_metadata = gene_metadata.new, 
                                          project = project.input, 
                                          min.cells = 0, 
                                          min.features = 0, 
                                          nFeature_RNA.lower_cutoff = 0, 
                                          nFeature_RNA.higher_cutoff = nFeature_RNA.higher_cutoff, 
                                          percent.mt.cutoff = 100 + 0.1)
            
        }
    }
    
    for (i in 1:length(exp.name)) {
        df.new <- eval(parse(text=paste('sample','.',sample[s],'.',exp.name[i],'.','df_new', sep = '')))
        
        if (! is.null(df.new)) {
            row.idx <- which(rownames(df.subset.all) %in% rownames(df.new))
            col.idx <- which(colnames(df.subset.all) %in% colnames(df.new))
            
            df.subset <- df.subset.all[row.idx, col.idx, drop = FALSE]
            
            write.csv(df.subset, file = paste(current_folder,'/Result','/','Mouse','/','library','/',sample[s],'/','QC','/',sample[s],'.QC','.',exp.name[i],'.csv', sep = ''), row.names = TRUE, quote = TRUE)
        }
    }
    
    if (sample[s] != "control") {
        
        # /// combine sgRNA and control ///
        mainDir <- paste(current_folder,'/Result','/','Mouse','/','library','/',sample[s], sep = '')
        subDir <- 'control'
        if (! file.exists(file.path(mainDir, subDir))) {
            dir.create(file.path(mainDir, subDir))
        }
        
        df.subset.all.control <- eval(parse(text=paste('sample','.','control','.','df_subset_all', sep = '')))
        
        df.subset.all.new <- df.subset.all
        if (! is.null(df.subset.all.new)) {
            rownames(df.subset.all.new)[which(rownames(df.subset.all.new) == "Tead2")] = "Tead2(ENSMUSG00000030796)"
            rownames(df.subset.all.new)[which(rownames(df.subset.all.new) == "l7Rn6")] = "l7Rn6(ENSMUSG00000062797)"
            
            row.idx <- unlist(sapply(control_gene_name, FUN = function(X) which(rownames(df.subset.all.new) %in% X)))
            df.subset.all.new <- df.subset.all.new[row.idx, , drop = FALSE]
            
            df.subset.all.new.dup.idx <- which(sapply(colnames(df.subset.all.new), FUN = function(X) X %in% colnames(df.subset.all.control)))
            df.subset.all.control.dup.idx <- which(sapply(colnames(df.subset.all.control), FUN = function(X) X %in% colnames(df.subset.all.new)))
            if (length(df.subset.all.new.dup.idx) > 0) {
                colnames(df.subset.all.new)[df.subset.all.new.dup.idx] <- paste(colnames(df.subset.all.new)[df.subset.all.new.dup.idx], '.', sample[s], sep = '')
            }
            if (length(df.subset.all.control.dup.idx) > 0) {
                colnames(df.subset.all.control)[df.subset.all.control.dup.idx] <- paste(colnames(df.subset.all.control)[df.subset.all.control.dup.idx], '.', 'control', sep = '')
            }
        }
        
        df2control <- merge(df.subset.all.new, df.subset.all.control, by = "row.names", all = TRUE, sort = TRUE)
        #df2control <- merge(df.subset.all.control, df.subset.all.new, by = "row.names", all = TRUE, sort = TRUE)
        rownames(df2control) <- df2control$Row.names
        df2control$Row.names <- NULL
        
        all_gene_name.new <- all_gene_name
        all_gene_name.new[which(all_gene_name.new == "Tead2")] = "Tead2(ENSMUSG00000030796)"
        all_gene_name.new[which(all_gene_name.new == "l7Rn6")] = "l7Rn6(ENSMUSG00000062797)"
        
        row.idx <- unlist(sapply(unique(c(all_gene_name.new, control_gene_name)), FUN = function(X) which(rownames(df2control) %in% X)))
        df2control <- df2control[row.idx, , drop = FALSE]
        
        # /// replace NA with true value ///
        row.idx1 <- unlist(sapply(rownames(df2control), FUN = function(X) which(control_gene_name %in% X)))
        row.idx2 <- unlist(sapply(rownames(df2control), FUN = function(X) which(all_gene_name.new %in% X)))
        
        preSum <- 0
        cell.num.new.table.new.idx.all <- c()
        for (run in 1:2) {
            if (run == 1) {
                cell.num.new.table.new <- eval(parse(text=paste('sample','.',sample[s],'.','cell_num_new_table', sep = '')))
            } else if (run == 2) {
                cell.num.new.table.new <- eval(parse(text=paste('sample','.','control','.','cell_num_new_table', sep = '')))
            }
            
            Zero.idx <- which(cell.num.new.table.new == 0)
            if (length(Zero.idx) > 0) {
                cell.num.new.table.new <- cell.num.new.table.new[-Zero.idx]
            }
            
            if (length(cell.num.new.table.new) > 0) {
                cell.num.new.table.new.idx <- sapply(c(1:length(cell.num.new.table.new)), FUN = function(X) c(1:cell.num.new.table.new[X]) + cumsum(cell.num.new.table.new)[X] - cell.num.new.table.new[X] + preSum)
                names(cell.num.new.table.new.idx) <- names(cell.num.new.table.new)
                
                cell.num.new.table.new.idx.all <- c(cell.num.new.table.new.idx.all, cell.num.new.table.new.idx)
            }
            
            preSum <- preSum + sum(cell.num.new.table.new)
        }
        
        sgRNA2cell <- sapply(cell.num.new.table.new.idx.all, FUN = function(X) colnames(df2control)[X])
        cell2sgRNA <- sapply(colnames(df2control), FUN = function(X) names(which(sapply(sgRNA2cell, FUN = function(Y) X %in% unlist(Y)))))
        col.name <- cell2sgRNA
        
        for (i in 1:length(exp.name.all)) {
            col.idx <- which(col.name == exp.name.all[i])
            
            match.idx <- grep('\\.', names(col.idx), value = FALSE)
            match.name <- grep('\\.', names(col.idx), value = TRUE)
            
            if (length(match.name) > 0) {
                match.idx <- match.idx[which(sapply(match.name, FUN = function(X) unlist(strsplit(X, '\\.'))[[2]] %in% sample))]
                if (length(match.idx) > 0) {
                    names(col.idx)[match.idx] = sapply(names(col.idx)[match.idx], FUN = function(X) unlist(strsplit(X, '\\.'))[[1]])
                }
            }
            
            if (exp.name.all[i] %in% control) {
                sample.name <- "control"
                sample_file <- exp.name.all[i]
                row.idx <- row.idx1
            } else {
                sample.name <- sample[s]
                sample_file <- exp.name.all[i]
                row.idx <- row.idx2
            }
            
            df <- eval(parse(text=paste('sample','.',sample.name,'.',sample_file,'.','df', sep = '')))
            
            if (! is.null(df)) {
                rownames(df)[which(rownames(df) == "Tead2")] = "Tead2(ENSMUSG00000030796)"
                rownames(df)[which(rownames(df) == "l7Rn6")] = "l7Rn6(ENSMUSG00000062797)"
                
                NA.idx <- which(! rownames(df2control) %in% rownames(df))
                df2control[NA.idx, col.idx][is.na(df2control[NA.idx, col.idx])] <- 0
                
                colnames(df) = sapply(colnames(df), FUN = function(X) unlist(strsplit(X, '\\.'))[[1]])
                names(col.idx) = sapply(names(col.idx), FUN = function(X) unlist(strsplit(X, '\\.'))[[1]])
                
                df2control[names(row.idx), col.idx] = df[row.idx, names(col.idx)]
            }
        }
        
        #df2control[is.na(df2control)] <- 0
        
        row.idx <- which(rownames(df2control) %in% intersect(all_gene_name.new, control_gene_name))
        df2control <- df2control[row.idx, , drop = FALSE]
        
        #assign(paste('sample','.',sample[s],'.','df2control', sep = ''), df2control)
        #write.csv(df2control, file = paste(current_folder,'/Result','/','Mouse','/','library','/',sample[s],'/',sample[s],'.control','.csv', sep = ''), row.names = TRUE, quote = TRUE)
        
        # /// split into individual sgRNA and control ///
        for (i in 1:length(exp.name.all)) {
            col.idx <- which(col.name == exp.name.all[i])
            
            if (length(col.idx) > 0) {
                df2control.subset <- df2control[, col.idx, drop = FALSE]
                #assign(paste('sample','.',sample[s],'.',exp.name.all[i],'.','df2control_subset', sep = ''), df2control.subset)
                write.csv(df2control.subset, file = paste(current_folder,'/Result','/','Mouse','/','library','/',sample[s],'/','control','/',sample[s],'.control','.',exp.name.all[i],'.csv', sep = ''), row.names = TRUE, quote = TRUE)
            }
        }
        
        
        # --- combine control and ch02 ----------------------------------------------------------------------------------------------------------------------
        if (sample[s] == "ch02") {
            
            sgRNAs_of_interest <- c(64, 65,             # AB2
                                    39,                 # toxic sgRNA
                                    66, 67, 68, 69,     # toxic sgRNA (strong toxic)
                                    71, 72)             # negative control sgRNA
            
            exp.name.all.control_idx <- which(exp.name.all %in% control)
            exp.name.all.control <- exp.name.all[exp.name.all.control_idx]
            
            exp.name.all.sgRNA_idx <- which(exp.name.all %in% sgRNA.barcode)
            exp.name.all.sgRNA <- exp.name.all[exp.name.all.sgRNA_idx]
            
            exp.name.all.sgRNAs_of_interest <- sapply(sgRNAs_of_interest, FUN = function(X) exp.name.all.sgRNA[which(sapply(exp.name.all.sgRNA, FUN = function(Y) unlist(strsplit(Y, '_'))[[1]]) %in% X)])
            
            col.idx <- sapply(c(exp.name.all.control, exp.name.all.sgRNAs_of_interest), FUN = function(X) which(col.name %in% X))
            col.idx <- unlist(col.idx)
            
            if (length(col.idx) > 0) {
                df2control.subset <- df2control[, col.idx, drop = FALSE]
                #assign(paste('sample','.',sample[s],'.','sgRNAs_of_interest','.','df2control_subset', sep = ''), df2control.subset)
                write.csv(df2control.subset, file = paste(current_folder,'/Result','/','Mouse','/','library','/',sample[s],'/','control','/',sample[s],'.control','.','sgRNAs_of_interest','.csv', sep = ''), row.names = TRUE, quote = TRUE)
            }
            
            # --- gene_metadata ---------------------------------------------------------------------------------------------------------------------------------
            expression_matrix = as.matrix(df2control.subset)
            
            row.idx <- unlist(sapply(rownames(df2control.subset), FUN = function(X) which(all_gene_name %in% X)))
            all_gene_id.subset.all <- all_gene_id[row.idx]
            
            gene_metadata <- data.frame("gene_short_name" = rownames(df2control.subset), 
                                        "gene_id" = all_gene_id.subset.all, 
                                        stringsAsFactors = FALSE)
            rownames(gene_metadata) <- gene_metadata[,"gene_short_name"]
            
            # --- cell_metadata ---------------------------------------------------------------------------------------------------------------------------------
            row.idx <- sapply(colnames(df2control.subset), FUN = function(X) which(rownames(cell_metadata.all) %in% X))
            cell_metadata <- cell_metadata.all[row.idx, , drop = FALSE]
            
            
            # --- Monocle3 & Seurat -----------------------------------------------------------------------------------------------------------------------------
            subdata <- list()
            subdata.name <- c()
            for (comb in 1:2) {
                if (comb == 1) {
                    subdata.counts <- c(acinar.counts, ductal.counts, cancer_tripleMut.counts)
                    subdata.old <- c(acinar.old, ductal.old, cancer_tripleMut.old)
                    subdata.new <- c(acinar.new, ductal.new, cancer_tripleMut.new)
                } else if (comb == 2) {
                    subdata.counts <- c(ductal.counts, cancer_tripleMut.counts)
                    subdata.old <- c(ductal.old, cancer_tripleMut.old)
                    subdata.new <- c(ductal.new, cancer_tripleMut.new)
                }
                subdata <- c(subdata, list(subdata.counts))
                #subdata <- c(subdata, list(subdata.old))
                #subdata <- c(subdata, list(subdata.new))
                
                subdata.name <- c(subdata.name, paste('subset',comb, sep = ''))
            }
            names(subdata) = subdata.name
            
            #for (sb in -1:length(subdata)) {
            for (sb in -1:-1) {
                
                expression_matrix.new = expression_matrix
                gene_metadata.new = gene_metadata
                cell_metadata.new = cell_metadata
                
                project.name = unlist(strsplit(current_folder, '/'))
                project.name = project.name[length(project.name)]
                
                if (sb == 0) {         # all
                    cells.input = control
                    project.input = paste(project.name,'.','all', sep = '')
                } else if (sb == -1) { # all_combine_paper
                    cells.input = control
                    project.input = paste(project.name,'.','all_combine_paper', sep = '')
                } else if (sb > 0) {   # subset
                    cells.input = sapply(unlist(subdata[sb]), FUN = function(X) control[which(control == X)])
                    project.input = paste(project.name,'.',names(subdata[sb]), sep = '')
                    
                    col.idx <- sapply(unlist(subdata[sb]), FUN = function(X) which(cell_metadata.new[,"sgRNA"] == X))
                    names(col.idx) = unlist(subdata[sb])
                    col.idx <- unlist(col.idx)
                    
                    expression_matrix.new = expression_matrix.new[, col.idx, drop = FALSE]
                    cell_metadata.new = cell_metadata.new[col.idx, , drop = FALSE]
                    for (k in 1:ncol(cell_metadata.new)) {
                        if (is.factor(cell_metadata.new[,k])) {
                            cell_metadata.new[,k] <- factor(cell_metadata.new[,k], levels = unique(cell_metadata.new[,k]))
                        }
                    }
                    
                    Zero.idx <- which(rowSums(expression_matrix.new) == 0)
                    if (length(Zero.idx) > 0) {
                        row.idx <- which(rowSums(expression_matrix.new) != 0)
                        
                        expression_matrix.new = expression_matrix.new[row.idx, , drop = FALSE]
                        gene_metadata.new = gene_metadata.new[row.idx, , drop = FALSE]
                    }
                }
                
                acinar.idx = which(cell_metadata.new[,"sgRNA"] %in% c("acinar","1_acinar","1_acinar_new"))
                ductal.idx = which(cell_metadata.new[,"sgRNA"] %in% c("counts_kcppc_new","1_ductal","kcppc_bam_counts_new"))
                cancer.idx = which(cell_metadata.new[,"sgRNA"] %in% c("cancer","1_cancer","1_cancer_new"))
                
                beta.idx = c()
                min6.idx = which(cell_metadata.new[,"sgRNA"] %in% "min6")
                if (length(min6.idx) > 0) {
                    beta.idx = min6.idx
                }
                min6_1.idx = which(cell_metadata.new[,"sgRNA"] %in% c("16726_MIN6","16726_min6"))
                min6_2.idx = which(cell_metadata.new[,"sgRNA"] %in% c("16728_MIN6","16728_min6"))
                if (length(min6_1.idx) > 0) {
                    beta.idx = min6_1.idx
                }
                if (length(min6_2.idx) > 0) {
                    beta.idx = min6_2.idx
                }
                if (length(min6_1.idx) > 0 & length(min6_2.idx) > 0) {
                    beta.idx = unique(c(min6_1.idx, min6_2.idx))
                }
                
                PCC8.idx = which(cell_metadata.new[,"sgRNA"] %in% "16728_PCC8")
                pcc8.idx = which(cell_metadata.new[,"sgRNA"] %in% c("pcc8","16728_pcc8"))
                
                PCC9.idx = which(cell_metadata.new[,"sgRNA"] %in% "16726_PCC9")
                pcc9.idx = which(cell_metadata.new[,"sgRNA"] %in% c("pcc9","16726_pcc9"))
                
                SRR8872326.idx = which(cell_metadata.new[,"sgRNA"] %in% "SRR8872326")
                SRR8872327.idx = which(cell_metadata.new[,"sgRNA"] %in% "SRR8872327")
                SRR8872328.idx = which(cell_metadata.new[,"sgRNA"] %in% "SRR8872328")
                SRR8872329.idx = which(cell_metadata.new[,"sgRNA"] %in% "SRR8872329")
                paper.fibroblast.idx = unlist(sapply(paper.fibroblast, FUN = function(X) which(cell_metadata.new[,"sgRNA"] %in% X)))
                paper.viable.idx = unlist(sapply(paper.viable, FUN = function(X) which(cell_metadata.new[,"sgRNA"] %in% X)))
                paper.cell.idx = unlist(sapply(paper.cell, FUN = function(X) which(cell_metadata.new[,"sgRNA"] %in% X)))
                
                sgRNA_all.idx = which(cell_metadata.new[,"sgRNA"] %in% sgRNA.barcode)
                for (b in 1:length(sgRNA.barcode)) {
                    assign(paste('sgRNA','.',sgRNA.barcode[b],'.','idx', sep = ''), which(cell_metadata.new[,"sgRNA"] %in% sgRNA.barcode[b]))
                }
                
                if (all(control.counts %in% control)) {
                    cell_metadata.new[acinar.idx,"batch"] = 1
                    cell_metadata.new[ductal.idx,"batch"] = 1
                    cell_metadata.new[beta.idx,"batch"] = 1
                    cell_metadata.new[PCC8.idx,"batch"] = 1
                    cell_metadata.new[pcc8.idx,"batch"] = 1
                    cell_metadata.new[PCC9.idx,"batch"] = 1
                    cell_metadata.new[pcc9.idx,"batch"] = 1
                    cell_metadata.new[cancer.idx,"batch"] = 1
                    cell_metadata.new[paper.fibroblast.idx,"batch"] = 3
                    cell_metadata.new[paper.viable.idx,"batch"] = 2
                    cell_metadata.new[sgRNA_all.idx,"batch"] = 4
                } else if (all(control.old %in% control)) {
                    cell_metadata.new[acinar.idx,"batch"] = 1
                    cell_metadata.new[ductal.idx,"batch"] = 1
                    cell_metadata.new[beta.idx,"batch"] = 2
                    cell_metadata.new[PCC8.idx,"batch"] = 2
                    cell_metadata.new[pcc8.idx,"batch"] = 2
                    cell_metadata.new[PCC9.idx,"batch"] = 2
                    cell_metadata.new[pcc9.idx,"batch"] = 2
                    cell_metadata.new[cancer.idx,"batch"] = 1
                    cell_metadata.new[paper.fibroblast.idx,"batch"] = 4
                    cell_metadata.new[paper.viable.idx,"batch"] = 3
                    cell_metadata.new[sgRNA_all.idx,"batch"] = 5
                } else if (all(control.new %in% control)) {
                    cell_metadata.new[acinar.idx,"batch"] = 1
                    cell_metadata.new[ductal.idx,"batch"] = 1
                    cell_metadata.new[beta.idx,"batch"] = 2
                    cell_metadata.new[PCC8.idx,"batch"] = 2
                    cell_metadata.new[pcc8.idx,"batch"] = 2
                    cell_metadata.new[PCC9.idx,"batch"] = 2
                    cell_metadata.new[pcc9.idx,"batch"] = 2
                    cell_metadata.new[cancer.idx,"batch"] = 1
                    cell_metadata.new[paper.fibroblast.idx,"batch"] = 4
                    cell_metadata.new[paper.viable.idx,"batch"] = 3
                    cell_metadata.new[sgRNA_all.idx,"batch"] = 5
                }
                cell_metadata.new[,"batch"] <- factor(cell_metadata.new[,"batch"], levels = sort(unique(cell_metadata.new[,"batch"])))
                
                cell_metadata.new[acinar.idx,"cell.type"] = "acinar"
                cell_metadata.new[ductal.idx,"cell.type"] = "ductal"
                cell_metadata.new[beta.idx,"cell.type"] = "beta"
                cell_metadata.new[PCC8.idx,"cell.type"] = "PCC8"
                cell_metadata.new[pcc8.idx,"cell.type"] = "pcc8"
                cell_metadata.new[PCC9.idx,"cell.type"] = "PCC9"
                cell_metadata.new[pcc9.idx,"cell.type"] = "pcc9"
                cell_metadata.new[cancer.idx,"cell.type"] = "cancer"
                cell_metadata.new[paper.fibroblast.idx,"cell.type"] = "fibroblast"  # Fibroblast-enriched fraction (DAPI-CD45-CD31-EPCAM-)
                cell_metadata.new[paper.viable.idx,"cell.type"] = "viable"          # Viable cell fraction (DAPI-)
                for (b in 1:length(sgRNA.barcode)) {
                    sgRNA_single.idx <- eval(parse(text=paste('sgRNA','.',sgRNA.barcode[b],'.','idx', sep = '')))
                    cell_metadata.new[sgRNA_single.idx,"cell.type"] = paste('sgRNA','_',b, sep = '')
                }
                cell_metadata.new[,"cell.type"] <- factor(cell_metadata.new[,"cell.type"], levels = unique(cell_metadata.new[,"cell.type"]))
                
                cell_metadata.new[acinar.idx,"cell.type.color"] = "#E77D72"
                cell_metadata.new[ductal.idx,"cell.type.color"] = "#B39F33"
                cell_metadata.new[beta.idx,"cell.type.color"] = "#53B74C"
                cell_metadata.new[PCC8.idx,"cell.type.color"] = "#00FDFF"
                cell_metadata.new[pcc8.idx,"cell.type.color"] = "#00FDFF"
                cell_metadata.new[PCC9.idx,"cell.type.color"] = "#6E9AF8"
                cell_metadata.new[pcc9.idx,"cell.type.color"] = "#6E9AF8"
                cell_metadata.new[cancer.idx,"cell.type.color"] = "#E46DDD"
                cell_metadata.new[paper.fibroblast.idx,"cell.type.color"] = "#53B49C"
                cell_metadata.new[paper.viable.idx,"cell.type.color"] = "#D73A35"
                cell_metadata.new[sgRNA_all.idx,"cell.type.color"] = "#A6A6A6"
                for (b in 1:length(sgRNA.barcode)) {
                    sgRNA_single.idx <- eval(parse(text=paste('sgRNA','.',sgRNA.barcode[b],'.','idx', sep = '')))
                    if (sgRNA.barcode[b] == "39_CGGATGCTAATGCACGTG") { # toxic sgRNA
                        cell_metadata.new[sgRNA_single.idx,"cell.type.color"] = "#D8A584"
                    } else if (sgRNA.barcode[b] == "64_CCCAGGGAGTTATTCGAT") { # AB2
                        cell_metadata.new[sgRNA_single.idx,"cell.type.color"] = "#FFC000"
                    } else if (sgRNA.barcode[b] == "65_ATCCAAGTACCACCCCGA") { # AB2
                        cell_metadata.new[sgRNA_single.idx,"cell.type.color"] = "#FFFF00"
                    } else if (sgRNA.barcode[b] == "66_TAAAATACAATCCCTCGG") { # toxic sgRNA (strong toxic)
                        cell_metadata.new[sgRNA_single.idx,"cell.type.color"] = "#AF5AFF"
                    } else if (sgRNA.barcode[b] == "67_GATGTATCGGAGTAGTTG") { # toxic sgRNA (strong toxic)
                        cell_metadata.new[sgRNA_single.idx,"cell.type.color"] = "#863CBD"
                    } else if (sgRNA.barcode[b] == "68_CCACAAGAATACATGCGT") { # toxic sgRNA (strong toxic)
                        cell_metadata.new[sgRNA_single.idx,"cell.type.color"] = "#5B2A82"
                    } else if (sgRNA.barcode[b] == "69_GGAGAATCGAGATGGTGG") { # toxic sgRNA (strong toxic)
                        cell_metadata.new[sgRNA_single.idx,"cell.type.color"] = "#35184C"
                    } else if (sgRNA.barcode[b] == "71_GTGCACACTAACAATGAG") { # negative control sgRNA 1
                        cell_metadata.new[sgRNA_single.idx,"cell.type.color"] = "#54667E"
                    } else if (sgRNA.barcode[b] == "72_GGCTAGAATGTGTACCAT") { # negative control sgRNA 2
                        cell_metadata.new[sgRNA_single.idx,"cell.type.color"] = "#6D9385"
                    }
                }
                cell_metadata.new[,"cell.type.color"] <- factor(cell_metadata.new[,"cell.type.color"], levels = unique(cell_metadata.new[,"cell.type.color"]))
                
                if (all(paper.cell %in% control)) {
                    cell_metadata.new[acinar.idx,"stage"] = 1
                    cell_metadata.new[ductal.idx,"stage"] = 2
                    cell_metadata.new[beta.idx,"stage"] = 1
                    
                    cell_metadata.new[PCC8.idx,"stage"] = 4
                    cell_metadata.new[pcc8.idx,"stage"] = 4
                    cell_metadata.new[PCC9.idx,"stage"] = 4
                    cell_metadata.new[pcc9.idx,"stage"] = 4
                    
                    cell_metadata.new[cancer.idx,"stage"] = 7
                    
                    cell_metadata.new[paper.fibroblast.idx,"stage"] = 3
                    cell_metadata.new[paper.viable.idx,"stage"] = 5
                    
                    cell_metadata.new[sgRNA_all.idx,"stage"] = 6
                } else {
                    cell_metadata.new[acinar.idx,"stage"] = 1
                    cell_metadata.new[ductal.idx,"stage"] = 2
                    cell_metadata.new[beta.idx,"stage"] = 1
                    
                    cell_metadata.new[PCC8.idx,"stage"] = 3
                    cell_metadata.new[pcc8.idx,"stage"] = 3
                    cell_metadata.new[PCC9.idx,"stage"] = 3
                    cell_metadata.new[pcc9.idx,"stage"] = 3
                    
                    cell_metadata.new[cancer.idx,"stage"] = 5
                    
                    cell_metadata.new[sgRNA_all.idx,"stage"] = 4
                }
                cell_metadata.new[,"stage"] <- factor(cell_metadata.new[,"stage"], levels = sort(unique(cell_metadata.new[,"stage"])))
                
                # /// Monocle3 ///
                # *** no genes of interest and signature genes ***
                Monocle3_analysis(expression_matrix = expression_matrix.new, 
                                  cell_metadata = cell_metadata.new, 
                                  gene_metadata = gene_metadata.new, 
                                  batch = "batch", 
                                  cell_type = "cell.type", 
                                  time = "stage", 
                                  cells = cells.input, 
                                  project = project.input, 
                                  analysis_type = "trajectories", 
                                  redDim_method = "UMAP", 
                                  genes_of_interest = NULL, 
                                  signature_genes = NULL)
                
                # /// Seurat ///
                min.cells = 0 # min.cells: Include features detected in at least this many cells. Will subset the counts matrix as well. To reintroduce excluded features, create a new object with a lower cutoff.
                # To filter out low-quality cells, we first removed cells for which less than 500 genes were detected.
                min.features = 500 # min.features: Include cells where at least this many features are detected.
                
                nFeature_RNA.lower_cutoff = min.features
                nFeature_RNA.higher_cutoff = 2500
                # To filter out low-quality cells, we removed cells for which over 10% genes derived from mitochondrial genome.
                percent.mt.cutoff = 10
                
                Seurat.obj <- Seurat_analysis(expression_matrix = expression_matrix.new, 
                                              cell_metadata = cell_metadata.new, 
                                              gene_metadata = gene_metadata.new, 
                                              project = project.input, 
                                              min.cells = 0, 
                                              min.features = 0, 
                                              nFeature_RNA.lower_cutoff = 0, 
                                              nFeature_RNA.higher_cutoff = nFeature_RNA.higher_cutoff, 
                                              percent.mt.cutoff = 100 + 0.1)
                
            }
        }
    }
}


RData.name = unlist(strsplit(current_folder, '/'))
RData.name = RData.name[length(RData.name)]
save.image(file = paste('./../',RData.name,'_','part3','.RData', sep = ''))


# --- Merge data with the same day ------------------------------------------------------------------------------------------------------------------
for (b in 1:length(sgRNA.barcode)) {
    
    gc() # Garbage Collection
    
    for (m in 1:length(sample.merged)) {
        rm.var_name <- c()
        if (exists(paste('sample','.',sample.merged[m],'.',sgRNA.barcode[b],'.','df', sep = ''))) {
            rm.var_name <- paste('sample','.',sample.merged[m],'.',sgRNA.barcode[b],'.','df', sep = '')
        }
        rm(list = eval(rm.var_name))
    }
    
    for (s in 1:length(sample)) {
        
        if (! sample[s] %in% c("control", "ch02")) {
            
            df <- eval(parse(text=paste('sample','.',sample[s],'.',sgRNA.barcode[b],'.','df', sep = '')))
            
            for (m in 1:length(sample.merged)) {
                
                df.merged <- c()
                
                run = 0
                
                if (sample.merged[m] == "d3_p1_1_merged") {
                    if (sample[s] %in% d3_p1_1) {
                        run = 1
                    }
                } else if (sample.merged[m] == "d3") {
                    if (sample[s] %in% d3) {
                        run = 1
                    }
                } else if (sample.merged[m] == "d5_p1_merged") {
                    if (sample[s] %in% d5_p1) {
                        run = 1
                    }
                } else if (sample.merged[m] == "d5") {
                    if (sample[s] %in% d5) {
                        run = 1
                    }
                } else if (sample.merged[m] == "d7_p1_merged") {
                    if (sample[s] %in% d7_p1) {
                        run = 1
                    }
                } else if (sample.merged[m] == "d7_p2") {
                    if (sample[s] %in% d7_p2) {
                        run = 1
                    }
                } else if (sample.merged[m] == "d7") {
                    if (sample[s] %in% d7) {
                        run = 1
                    }
                }
                
                if (run == 1) {
                    if (exists(paste('sample','.',sample.merged[m],'.',sgRNA.barcode[b],'.','df', sep = ''))) {
                        df.merged <- eval(parse(text=paste('sample','.',sample.merged[m],'.',sgRNA.barcode[b],'.','df', sep = '')))
                    } else {
                        df.merged <- c()
                    }
                    df.merged <- merge(df.merged, df, by = "row.names", all = TRUE, sort = TRUE)
                    rownames(df.merged) <- df.merged$Row.names
                    df.merged$Row.names <- NULL
                    
                    assign(paste('sample','.',sample.merged[m],'.',sgRNA.barcode[b],'.','df', sep = ''), df.merged)
                }
            }
        }
    }
    
    mainDir <- paste(current_folder,'/Result','/','Mouse','/','library', sep = '')
    subDir <- 'merged'
    if (! file.exists(file.path(mainDir, subDir))) {
        dir.create(file.path(mainDir, subDir))
    }
    
    for (m in 1:length(sample.merged)) {
        
        mainDir <- paste(current_folder,'/Result','/','Mouse','/','library','/','merged', sep = '')
        subDir <- sample.merged[m]
        if (! file.exists(file.path(mainDir, subDir))) {
            dir.create(file.path(mainDir, subDir))
        }
        
        mainDir <- paste(current_folder,'/Result','/','Mouse','/','library','/','merged','/',sample.merged[m], sep = '')
        subDir <- 'raw'
        if (! file.exists(file.path(mainDir, subDir))) {
            dir.create(file.path(mainDir, subDir))
        }
        
        #mainDir <- paste(current_folder,'/Result','/','Mouse','/','library','/','merged','/',sample.merged[m], sep = '')
        #subDir <- 'QC'
        #if (! file.exists(file.path(mainDir, subDir))) {
        #    dir.create(file.path(mainDir, subDir))
        #}
        
        if (exists(paste('sample','.',sample.merged[m],'.',sgRNA.barcode[b],'.','df', sep = ''))) {
            df.merged <- eval(parse(text=paste('sample','.',sample.merged[m],'.',sgRNA.barcode[b],'.','df', sep = '')))
            
            colnames(df.merged) <- sapply(colnames(df.merged), FUN = function(X) unlist(strsplit(X, '\\.'))[[1]])
            
            dup.idx <- which(duplicated(colnames(df.merged)))
            if (length(dup.idx) > 0) {
                rm.col.idx.all <- c()
                for (d in 1:length(dup.idx)) {
                    col.idx <- which(colnames(df.merged) == colnames(df.merged)[dup.idx[d]])
                    
                    keep.col.idx = col.idx[1]
                    rm.col.idx = col.idx[which(! col.idx %in% keep.col.idx)]
                    rm.col.idx.all <- c(rm.col.idx.all, rm.col.idx)
                    
                    df.merged[, keep.col.idx] = apply(df.merged[,col.idx], 1, sum) # for a matrix 1 indicates rows, 2 indicates columns
                }
                df.merged <- df.merged[, -rm.col.idx.all, drop = FALSE]
            }
            
            if (! is.null(df.merged)) {
                cell.num = ncol(df.merged)
                gene.num = nrow(df.merged)
                
                if (nrow(df.merged) > 0) {
                    df.merged[is.na(df.merged)] <- 0
                    
                    row.idx <- unlist(sapply(all_gene_name, FUN = function(X) which(rownames(df.merged) %in% X)))
                    df.merged <- df.merged[row.idx, , drop = FALSE]
                    assign(paste('sample','.',sample.merged[m],'.',sgRNA.barcode[b],'.','df', sep = ''), df.merged)
                    write.csv(df.merged, file = paste(current_folder,'/Result','/','Mouse','/','library','/','merged','/',sample.merged[m],'/','raw','/',sample.merged[m],'.',sgRNA.barcode[b],'.csv', sep = ''), row.names = TRUE, quote = TRUE)
                    
                    
                    # --- Data filtering --------------------------------------------------------------------------------------------------------------------------------
                    for (g in 1:length(gene_cutoff)) {
                        df.merged_qc <- qc_data(df.merged, normalize = FALSE, min_gene_frac = gene_cutoff[g])
                        df.merged_qc.normalized <- qc_data(df.merged, normalize = TRUE, min_gene_frac = gene_cutoff[g])
                        
                        df.merged.new <- df.merged_qc
                        df.merged.new.normalized <- df.merged_qc.normalized
                        
                        if (! is.null(df.merged.new)) {
                            cell.num.new = ncol(df.merged.new)
                            gene.num.new = nrow(df.merged.new)
                        } else {
                            cell.num.new = 0
                            gene.num.new = 0
                        }
                        
                        assign(paste('sample','.',sample.merged[m],'.',sgRNA.barcode[b],'.',gene_cutoff[g],'.','cell_num_new', sep = ''), cell.num.new)
                        assign(paste('sample','.',sample.merged[m],'.',sgRNA.barcode[b],'.',gene_cutoff[g],'.','gene_num_new', sep = ''), gene.num.new)
                    }
                    
                    df.merged_qc <- qc_data(df.merged, normalize = FALSE, min_gene_frac = gene_cutoff.default)
                    df.merged_qc.normalized <- qc_data(df.merged, normalize = TRUE, min_gene_frac = gene_cutoff.default)
                    
                    df.merged.new <- df.merged_qc
                    df.merged.new.normalized <- df.merged_qc.normalized
                    
                    assign(paste('sample','.',sample.merged[m],'.',sgRNA.barcode[b],'.','df_new', sep = ''), df.merged.new)
                    #assign(paste('sample','.',sample.merged[m],'.',sgRNA.barcode[b],'.','df_new_normalized', sep = ''), df.merged.new.normalized)
                    
                    if (! is.null(df.merged.new)) {
                        cell.num.new = ncol(df.merged.new)
                        gene.num.new = nrow(df.merged.new)
                        
                        #write.csv(df.merged.new, file = paste(current_folder,'/Result','/','Mouse','/','library','/','merged','/',sample.merged[m],'/','QC','/',sample.merged[m],'.QC','.',sgRNA.barcode[b],'.csv', sep = ''), row.names = TRUE, quote = TRUE)
                    } else {
                        cell.num.new = 0
                        gene.num.new = 0
                    }
                } else {
                    assign(paste('sample','.',sample.merged[m],'.',sgRNA.barcode[b],'.','df', sep = ''), NULL)
                    assign(paste('sample','.',sample.merged[m],'.',sgRNA.barcode[b],'.','df_new', sep = ''), NULL)
                    #assign(paste('sample','.',sample.merged[m],'.',sgRNA.barcode[b],'.','df_new_normalized', sep = ''), NULL)
                    
                    for (g in 1:length(gene_cutoff)) {
                        cell.num.new = 0
                        gene.num.new = 0
                        
                        assign(paste('sample','.',sample.merged[m],'.',sgRNA.barcode[b],'.',gene_cutoff[g],'.','cell_num_new', sep = ''), cell.num.new)
                        assign(paste('sample','.',sample.merged[m],'.',sgRNA.barcode[b],'.',gene_cutoff[g],'.','gene_num_new', sep = ''), gene.num.new)
                    }
                    
                    cell.num.new = 0
                    gene.num.new = 0
                }
            } else {
                cell.num = 0
                gene.num = 0
                
                assign(paste('sample','.',sample.merged[m],'.',sgRNA.barcode[b],'.','df', sep = ''), NULL)
                assign(paste('sample','.',sample.merged[m],'.',sgRNA.barcode[b],'.','df_new', sep = ''), NULL)
                #assign(paste('sample','.',sample.merged[m],'.',sgRNA.barcode[b],'.','df_new_normalized', sep = ''), NULL)
                
                for (g in 1:length(gene_cutoff)) {
                    cell.num.new = 0
                    gene.num.new = 0
                    
                    assign(paste('sample','.',sample.merged[m],'.',sgRNA.barcode[b],'.',gene_cutoff[g],'.','cell_num_new', sep = ''), cell.num.new)
                    assign(paste('sample','.',sample.merged[m],'.',sgRNA.barcode[b],'.',gene_cutoff[g],'.','gene_num_new', sep = ''), gene.num.new)
                }
                
                cell.num.new = 0
                gene.num.new = 0
            }
        } else {
            cell.num = 0
            gene.num = 0
            
            assign(paste('sample','.',sample.merged[m],'.',sgRNA.barcode[b],'.','df', sep = ''), NULL)
            assign(paste('sample','.',sample.merged[m],'.',sgRNA.barcode[b],'.','df_new', sep = ''), NULL)
            #assign(paste('sample','.',sample.merged[m],'.',sgRNA.barcode[b],'.','df_new_normalized', sep = ''), NULL)
            
            for (g in 1:length(gene_cutoff)) {
                cell.num.new = 0
                gene.num.new = 0
                
                assign(paste('sample','.',sample.merged[m],'.',sgRNA.barcode[b],'.',gene_cutoff[g],'.','cell_num_new', sep = ''), cell.num.new)
                assign(paste('sample','.',sample.merged[m],'.',sgRNA.barcode[b],'.',gene_cutoff[g],'.','gene_num_new', sep = ''), gene.num.new)
            }
            
            cell.num.new = 0
            gene.num.new = 0
        }
        
        names(cell.num) <- sgRNA.barcode[b]
        assign(paste('sample','.',sample.merged[m],'.',sgRNA.barcode[b],'.','cell_num', sep = ''), cell.num)
        
        names(gene.num) <- sgRNA.barcode[b]
        assign(paste('sample','.',sample.merged[m],'.',sgRNA.barcode[b],'.','gene_num', sep = ''), gene.num)
        
        names(cell.num.new) <- sgRNA.barcode[b]
        assign(paste('sample','.',sample.merged[m],'.',sgRNA.barcode[b],'.','cell_num_new', sep = ''), cell.num.new)
        
        names(gene.num.new) <- sgRNA.barcode[b]
        assign(paste('sample','.',sample.merged[m],'.',sgRNA.barcode[b],'.','gene_num_new', sep = ''), gene.num.new)
        
    }
}

for (m in 1:length(sample.merged)) {
    
    cell.num.table <- c()
    cell.num.new.table <- c()
    
    gene.num.table <- c()
    gene.num.new.table <- c()
    
    for (b in 1:length(sgRNA.barcode)) {
        
        cell.num <- eval(parse(text=paste('sample','.',sample.merged[m],'.',sgRNA.barcode[b],'.','cell_num', sep = '')))
        gene.num <- eval(parse(text=paste('sample','.',sample.merged[m],'.',sgRNA.barcode[b],'.','gene_num', sep = '')))
        
        cell.num.new <- eval(parse(text=paste('sample','.',sample.merged[m],'.',sgRNA.barcode[b],'.','cell_num_new', sep = '')))
        gene.num.new <- eval(parse(text=paste('sample','.',sample.merged[m],'.',sgRNA.barcode[b],'.','gene_num_new', sep = '')))
        
        cell.num.table <- c(cell.num.table, cell.num)
        gene.num.table <- c(gene.num.table, gene.num)
        
        cell.num.new.table <- c(cell.num.new.table, cell.num.new)
        gene.num.new.table <- c(gene.num.new.table, gene.num.new)
        
    }
    
    assign(paste('sample','.',sample.merged[m],'.','cell_num_table', sep = ''), cell.num.table)
    assign(paste('sample','.',sample.merged[m],'.','cell_num_new_table', sep = ''), cell.num.new.table)
    
    cat('[',sample.merged[m],'] Merged cell number before QC: ',sum(cell.num.table),'\n', sep = '')
    print(cell.num.table)
    cat('\n')
    cat('[',sample.merged[m],'] Merged cell number after QC: ',sum(cell.num.new.table),'\n', sep = '')
    print(cell.num.new.table)
    cat('\n')
    #cat('\n')
    
    assign(paste('sample','.',sample.merged[m],'.','gene_num_table', sep = ''), gene.num.table)
    assign(paste('sample','.',sample.merged[m],'.','gene_num_new_table', sep = ''), gene.num.new.table)
    
    cat('[',sample.merged[m],'] Merged gene number before QC: ',sum(gene.num.table),'\n', sep = '')
    print(gene.num.table)
    cat('\n')
    cat('[',sample.merged[m],'] Merged gene number after QC: ',sum(gene.num.new.table),'\n', sep = '')
    print(gene.num.new.table)
    cat('\n')
    cat('\n')
    
}

exp.name <- sgRNA.barcode

exp.name.all <- c(control, sgRNA.barcode)
exp.name.all <- unique(exp.name.all)

for (m in 1:length(sample.merged)) {
    
    # --- Data correction and combination ---------------------------------------------------------------------------------------------------------------
    # /// replace NA with true value ///
    mainDir <- paste(current_folder,'/Result','/','Mouse','/','library','/','merged','/',sample.merged[m], sep = '')
    subDir <- 'QC'
    if (! file.exists(file.path(mainDir, subDir))) {
        dir.create(file.path(mainDir, subDir))
    }
    
    initial.idx = 1
    initial2end.idx <- c()
    df.subset.all <- c()
    for (i in 1:length(exp.name)) {
        df <- eval(parse(text=paste('sample','.',sample.merged[m],'.',exp.name[i],'.','df', sep = '')))
        df.origin <- df
        
        df.new <- eval(parse(text=paste('sample','.',sample.merged[m],'.',exp.name[i],'.','df_new', sep = '')))
        
        if (! is.null(df)) {
            colnames(df) = sapply(colnames(df), FUN = function(X) unlist(strsplit(X, '\\.'))[[1]])
            colnames(df.new) = sapply(colnames(df.new), FUN = function(X) unlist(strsplit(X, '\\.'))[[1]])
            
            row.idx <- unlist(sapply(rownames(df.new), FUN = function(X) which(rownames(df) %in% X)))
            col.idx <- unlist(sapply(colnames(df.new), FUN = function(X) which(colnames(df) %in% X)))
            
            df.subset <- df.origin[row.idx, col.idx, drop = FALSE]
            
            #assign(paste('sample','.',sample.merged[m],'.',exp.name[i],'.','df_subset', sep = ''), df.subset)
            #write.csv(df.subset, file = paste(current_folder,'/Result','/','Mouse','/','library','/','merged','/',sample.merged[m],'/','QC','/',sample.merged[m],'.QC','.',exp.name[i],'.csv', sep = ''), row.names = TRUE, quote = TRUE)
            
            initial2end.idx <- c(initial2end.idx, list(seq(initial.idx, length.out = ncol(df.subset))))
            names(initial2end.idx)[length(initial2end.idx)] <- exp.name[i]
            initial.idx = initial.idx + ncol(df.subset)
            
            if (! is.null(df.subset.all)) {
                df.subset.all.colnames = sapply(colnames(df.subset.all), FUN = function(X) unlist(strsplit(X, '\\.'))[[1]])
                df.subset.colnames = sapply(colnames(df.subset), FUN = function(X) unlist(strsplit(X, '\\.'))[[1]])
                
                df.subset.all.dup.idx <- which(sapply(df.subset.all.colnames, FUN = function(X) X %in% df.subset.colnames))
                df.subset.dup.idx <- which(sapply(df.subset.colnames, FUN = function(X) X %in% df.subset.all.colnames))
                if (length(df.subset.all.dup.idx) > 0) {
                    for (d in 1:length(df.subset.all.dup.idx)) {
                        if (length(grep('\\.', colnames(df.subset.all)[df.subset.all.dup.idx[d]], value = TRUE)) == 0) {
                            sample.name <- names(which(sapply(initial2end.idx, FUN = function(X) df.subset.all.dup.idx[d] %in% unlist(X))))
                            colnames(df.subset.all)[df.subset.all.dup.idx[d]] <- paste(colnames(df.subset.all)[df.subset.all.dup.idx[d]], '.', sample.name, sep = '')
                        }
                    }
                }
                if (length(df.subset.dup.idx) > 0) {
                    for (d in 1:length(df.subset.dup.idx)) {
                        if (length(grep('\\.', colnames(df.subset)[df.subset.dup.idx[d]], value = TRUE)) == 0) {
                            sample.name <- exp.name[i]
                            colnames(df.subset)[df.subset.dup.idx[d]] <- paste(colnames(df.subset)[df.subset.dup.idx[d]], '.', sample.name, sep = '')
                        }
                    }
                }
            }
            
            df.subset.all <- merge(df.subset.all, df.subset, by = "row.names", all = TRUE, sort = TRUE)
            rownames(df.subset.all) <- df.subset.all$Row.names
            df.subset.all$Row.names <- NULL
        } else {
            initial2end.idx <- c(initial2end.idx, NA)
            names(initial2end.idx)[length(initial2end.idx)] <- exp.name[i]
        }
    }
    
    if (! is.null(df.subset.all)) {
        df.subset.all[is.na(df.subset.all)] <- 0
        
        row.idx <- unlist(sapply(all_gene_name, FUN = function(X) which(rownames(df.subset.all) %in% X)))
        df.subset.all <- df.subset.all[row.idx, , drop = FALSE]
        
        assign(paste('sample','.',sample.merged[m],'.','df_subset_all', sep = ''), df.subset.all)
        #write.csv(df.subset.all, file = paste(current_folder,'/Result','/','Mouse','/','library','/','merged','/',sample.merged[m],'/',sample.merged[m],'.QC','.csv', sep = ''), row.names = TRUE, quote = TRUE)
        
        for (i in 1:length(exp.name)) {
            df.new <- eval(parse(text=paste('sample','.',sample.merged[m],'.',exp.name[i],'.','df_new', sep = '')))
            
            if (! is.null(df.new)) {
                row.idx <- which(rownames(df.subset.all) %in% rownames(df.new))
                col.idx <- which(colnames(df.subset.all) %in% colnames(df.new))
                
                df.subset <- df.subset.all[row.idx, col.idx, drop = FALSE]
                
                write.csv(df.subset, file = paste(current_folder,'/Result','/','Mouse','/','library','/','merged','/',sample.merged[m],'/','QC','/',sample.merged[m],'.QC','.',exp.name[i],'.csv', sep = ''), row.names = TRUE, quote = TRUE)
            }
        }
    }
    
    if (sample.merged[m] != "control") {
        
        # /// combine sgRNA and control ///
        mainDir <- paste(current_folder,'/Result','/','Mouse','/','library','/','merged','/',sample.merged[m], sep = '')
        subDir <- 'control'
        if (! file.exists(file.path(mainDir, subDir))) {
            dir.create(file.path(mainDir, subDir))
        }
        
        df.subset.all.control <- eval(parse(text=paste('sample','.','control','.','df_subset_all', sep = '')))
        
        df.subset.all.new <- df.subset.all
        if (! is.null(df.subset.all.new)) {
            rownames(df.subset.all.new)[which(rownames(df.subset.all.new) == "Tead2")] = "Tead2(ENSMUSG00000030796)"
            rownames(df.subset.all.new)[which(rownames(df.subset.all.new) == "l7Rn6")] = "l7Rn6(ENSMUSG00000062797)"
            
            row.idx <- unlist(sapply(control_gene_name, FUN = function(X) which(rownames(df.subset.all.new) %in% X)))
            df.subset.all.new <- df.subset.all.new[row.idx, , drop = FALSE]
            
            df.subset.all.new.dup.idx <- which(sapply(colnames(df.subset.all.new), FUN = function(X) X %in% colnames(df.subset.all.control)))
            df.subset.all.control.dup.idx <- which(sapply(colnames(df.subset.all.control), FUN = function(X) X %in% colnames(df.subset.all.new)))
            if (length(df.subset.all.new.dup.idx) > 0) {
                colnames(df.subset.all.new)[df.subset.all.new.dup.idx] <- paste(colnames(df.subset.all.new)[df.subset.all.new.dup.idx], '.', sample.merged[m], sep = '')
            }
            if (length(df.subset.all.control.dup.idx) > 0) {
                colnames(df.subset.all.control)[df.subset.all.control.dup.idx] <- paste(colnames(df.subset.all.control)[df.subset.all.control.dup.idx], '.', 'control', sep = '')
            }
        }
        
        df2control <- merge(df.subset.all.new, df.subset.all.control, by = "row.names", all = TRUE, sort = TRUE)
        #df2control <- merge(df.subset.all.control, df.subset.all.new, by = "row.names", all = TRUE, sort = TRUE)
        rownames(df2control) <- df2control$Row.names
        df2control$Row.names <- NULL
        
        all_gene_name.new <- all_gene_name
        all_gene_name.new[which(all_gene_name.new == "Tead2")] = "Tead2(ENSMUSG00000030796)"
        all_gene_name.new[which(all_gene_name.new == "l7Rn6")] = "l7Rn6(ENSMUSG00000062797)"
        
        row.idx <- unlist(sapply(unique(c(all_gene_name.new, control_gene_name)), FUN = function(X) which(rownames(df2control) %in% X)))
        df2control <- df2control[row.idx, , drop = FALSE]
        
        # /// replace NA with true value ///
        row.idx1 <- unlist(sapply(rownames(df2control), FUN = function(X) which(control_gene_name %in% X)))
        row.idx2 <- unlist(sapply(rownames(df2control), FUN = function(X) which(all_gene_name.new %in% X)))
        
        preSum <- 0
        cell.num.new.table.new.idx.all <- c()
        for (run in 1:2) {
            if (run == 1) {
                cell.num.new.table.new <- eval(parse(text=paste('sample','.',sample.merged[m],'.','cell_num_new_table', sep = '')))
            } else if (run == 2) {
                cell.num.new.table.new <- eval(parse(text=paste('sample','.','control','.','cell_num_new_table', sep = '')))
            }
            
            Zero.idx <- which(cell.num.new.table.new == 0)
            if (length(Zero.idx) > 0) {
                cell.num.new.table.new <- cell.num.new.table.new[-Zero.idx]
            }
            
            if (length(cell.num.new.table.new) > 0) {
                cell.num.new.table.new.idx <- sapply(c(1:length(cell.num.new.table.new)), FUN = function(X) c(1:cell.num.new.table.new[X]) + cumsum(cell.num.new.table.new)[X] - cell.num.new.table.new[X] + preSum)
                names(cell.num.new.table.new.idx) <- names(cell.num.new.table.new)
                
                cell.num.new.table.new.idx.all <- c(cell.num.new.table.new.idx.all, cell.num.new.table.new.idx)
            }
            
            preSum <- preSum + sum(cell.num.new.table.new)
        }
        
        sgRNA2cell <- sapply(cell.num.new.table.new.idx.all, FUN = function(X) colnames(df2control)[X])
        cell2sgRNA <- sapply(colnames(df2control), FUN = function(X) names(which(sapply(sgRNA2cell, FUN = function(Y) X %in% unlist(Y)))))
        col.name <- cell2sgRNA
        
        for (i in 1:length(exp.name.all)) {
            col.idx <- which(col.name == exp.name.all[i])
            
            match.idx <- grep('\\.', names(col.idx), value = FALSE)
            match.name <- grep('\\.', names(col.idx), value = TRUE)
            
            if (length(match.name) > 0) {
                match.idx <- match.idx[which(sapply(match.name, FUN = function(X) unlist(strsplit(X, '\\.'))[[2]] %in% sample.reorder))]
                if (length(match.idx) > 0) {
                    names(col.idx)[match.idx] = sapply(names(col.idx)[match.idx], FUN = function(X) unlist(strsplit(X, '\\.'))[[1]])
                }
            }
            
            if (exp.name.all[i] %in% control) {
                sample.name <- "control"
                sample_file <- exp.name.all[i]
                row.idx <- row.idx1
            } else {
                sample.name <- sample.merged[m]
                sample_file <- exp.name.all[i]
                row.idx <- row.idx2
            }
            
            df <- eval(parse(text=paste('sample','.',sample.name,'.',sample_file,'.','df', sep = '')))
            
            if (! is.null(df)) {
                rownames(df)[which(rownames(df) == "Tead2")] = "Tead2(ENSMUSG00000030796)"
                rownames(df)[which(rownames(df) == "l7Rn6")] = "l7Rn6(ENSMUSG00000062797)"
                
                NA.idx <- which(! rownames(df2control) %in% rownames(df))
                df2control[NA.idx, col.idx][is.na(df2control[NA.idx, col.idx])] <- 0
                
                colnames(df) = sapply(colnames(df), FUN = function(X) unlist(strsplit(X, '\\.'))[[1]])
                names(col.idx) = sapply(names(col.idx), FUN = function(X) unlist(strsplit(X, '\\.'))[[1]])
                
                df2control[names(row.idx), col.idx] = df[row.idx, names(col.idx)]
            }
        }
        
        #df2control[is.na(df2control)] <- 0
        
        row.idx <- which(rownames(df2control) %in% intersect(all_gene_name.new, control_gene_name))
        df2control <- df2control[row.idx, , drop = FALSE]
        
        #assign(paste('sample','.',sample.merged[m],'.','df2control', sep = ''), df2control)
        #write.csv(df2control, file = paste(current_folder,'/Result','/','Mouse','/','library','/','merged','/',sample.merged[m],'/',sample.merged[m],'.control','.csv', sep = ''), row.names = TRUE, quote = TRUE)
        
        # /// split into individual sgRNA and control ///
        for (i in 1:length(exp.name.all)) {
            col.idx <- which(col.name == exp.name.all[i])
            
            if (length(col.idx) > 0) {
                df2control.subset <- df2control[, col.idx, drop = FALSE]
                #assign(paste('sample','.',sample.merged[m],'.',exp.name.all[i],'.','df2control_subset', sep = ''), df2control.subset)
                write.csv(df2control.subset, file = paste(current_folder,'/Result','/','Mouse','/','library','/','merged','/',sample.merged[m],'/','control','/',sample.merged[m],'.control','.',exp.name.all[i],'.csv', sep = ''), row.names = TRUE, quote = TRUE)
            }
        }
    }
}


# --- Merge data with different days ----------------------------------------------------------------------------------------------------------------
# /// sgRNA across different days (e.g., D3, D5, D7, ...) ///
for (b in 1:length(sgRNA.barcode)) {
    
    gc() # Garbage Collection
    
    mainDir <- paste(current_folder,'/Result','/','Mouse', sep = '')
    subDir <- 'sgRNA'
    if (! file.exists(file.path(mainDir, subDir))) {
        dir.create(file.path(mainDir, subDir))
    }
    
    initial.idx = 1
    initial2end.idx <- c()
    df.all <- c()
    for (s in 1:length(sample.reorder)) {
        
        if (! sample.reorder[s] %in% c("control", "ch02", "poor_d3_p1_1", "poor_d5_p1", "poor_d7_p1")) {
            
            df.new <- eval(parse(text=paste('sample','.',sample.reorder[s],'.',sgRNA.barcode[b],'.','df_new', sep = '')))
            
            if (! is.null(df.new)) {
                initial2end.idx <- c(initial2end.idx, list(seq(initial.idx, length.out = ncol(df.new))))
                names(initial2end.idx)[length(initial2end.idx)] <- sample.reorder[s]
                initial.idx = initial.idx + ncol(df.new)
                
                if (! is.null(df.all)) {
                    df.all.colnames = sapply(colnames(df.all), FUN = function(X) unlist(strsplit(X, '\\.'))[[1]])
                    df.new.colnames = sapply(colnames(df.new), FUN = function(X) unlist(strsplit(X, '\\.'))[[1]])
                    
                    df.all.dup.idx <- which(sapply(df.all.colnames, FUN = function(X) X %in% df.new.colnames))
                    df.new.dup.idx <- which(sapply(df.new.colnames, FUN = function(X) X %in% df.all.colnames))
                    if (length(df.all.dup.idx) > 0) {
                        for (d in 1:length(df.all.dup.idx)) {
                            if (length(grep('\\.', colnames(df.all)[df.all.dup.idx[d]], value = TRUE)) == 0) {
                                sample.name <- names(which(sapply(initial2end.idx, FUN = function(X) df.all.dup.idx[d] %in% unlist(X))))
                                colnames(df.all)[df.all.dup.idx[d]] <- paste(colnames(df.all)[df.all.dup.idx[d]], '.', sample.name, sep = '')
                            }
                        }
                    }
                    if (length(df.new.dup.idx) > 0) {
                        colnames(df.new)[df.new.dup.idx] <- paste(colnames(df.new)[df.new.dup.idx], '.', sample.reorder[s], sep = '')
                    }
                }
                
                df.all <- merge(df.all, df.new, by = "row.names", all = TRUE, sort = TRUE)
                rownames(df.all) <- df.all$Row.names
                df.all$Row.names <- NULL
            }
        }
    }
    
    if (! is.null(df.all)) {
        
        if (nrow(df.all) > 0) {
            
            mainDir <- paste(current_folder,'/Result','/','Mouse','/','sgRNA', sep = '')
            subDir <- sgRNA.barcode[b]
            if (! file.exists(file.path(mainDir, subDir))) {
                dir.create(file.path(mainDir, subDir))
            }
            
            #df.all[is.na(df.all)] <- 0
            
            row.idx <- unlist(sapply(all_gene_name, FUN = function(X) which(rownames(df.all) %in% X)))
            df.all <- df.all[row.idx, , drop = FALSE]
            
            assign(paste('sgRNA','.',sgRNA.barcode[b],'.','df_all', sep = ''), df.all)
            
            
            # --- Data correction and combination ---------------------------------------------------------------------------------------------------------------
            exp.name <- sample.reorder[which(! sample.reorder %in% c("control", "ch02", "poor_d3_p1_1", "poor_d5_p1", "poor_d7_p1"))]
            
            exp.idx <- lapply(exp.name, FUN = function(X) which(colnames(df.all) %in% paste(colnames(eval(parse(text=paste('sample','.',X,'.',sgRNA.barcode[b],'.','df_new', sep = '')))),'.',X, sep = '')))
            names(exp.idx) <- exp.name
            
            exp.num <- sapply(exp.idx, FUN = function(X) length(X))
            
            # /// replace NA with true value ///
            mainDir <- paste(current_folder,'/Result','/','Mouse','/','sgRNA','/',sgRNA.barcode[b], sep = '')
            subDir <- 'QC'
            if (! file.exists(file.path(mainDir, subDir))) {
                dir.create(file.path(mainDir, subDir))
            }
            
            initial.idx = 1
            initial2end.idx <- c()
            df.subset.all <- c()
            for (i in 1:length(exp.name)) {
                df <- eval(parse(text=paste('sample','.',exp.name[i],'.',sgRNA.barcode[b],'.','df', sep = '')))
                df.origin <- df
                
                df.all.subset <- df.all[, exp.idx[[i]], drop = FALSE]
                
                if (! is.null(df)) {
                    colnames(df) = sapply(colnames(df), FUN = function(X) unlist(strsplit(X, '\\.'))[[1]])
                    colnames(df.all.subset) = sapply(colnames(df.all.subset), FUN = function(X) unlist(strsplit(X, '\\.'))[[1]])
                    
                    row.idx <- unlist(sapply(rownames(df.all.subset), FUN = function(X) which(rownames(df) %in% X)))
                    col.idx <- unlist(sapply(colnames(df.all.subset), FUN = function(X) which(colnames(df) %in% X)))
                    
                    df.subset <- df.origin[row.idx, col.idx, drop = FALSE]
                    
                    #assign(paste('sgRNA','.',sgRNA.barcode[b],'.',exp.name[i],'.','df_subset', sep = ''), df.subset)
                    #write.csv(df.subset, file = paste(current_folder,'/Result','/','Mouse','/','sgRNA','/',sgRNA.barcode[b],'/','QC','/',sgRNA.barcode[b],'.QC','.',exp.name[i],'.csv', sep = ''), row.names = TRUE, quote = TRUE)
                    
                    initial2end.idx <- c(initial2end.idx, list(seq(initial.idx, length.out = ncol(df.subset))))
                    names(initial2end.idx)[length(initial2end.idx)] <- exp.name[i]
                    initial.idx = initial.idx + ncol(df.subset)
                    
                    if (! is.null(df.subset.all)) {
                        df.subset.all.colnames = sapply(colnames(df.subset.all), FUN = function(X) unlist(strsplit(X, '\\.'))[[1]])
                        df.subset.colnames = sapply(colnames(df.subset), FUN = function(X) unlist(strsplit(X, '\\.'))[[1]])
                        
                        df.subset.all.dup.idx <- which(sapply(df.subset.all.colnames, FUN = function(X) X %in% df.subset.colnames))
                        df.subset.dup.idx <- which(sapply(df.subset.colnames, FUN = function(X) X %in% df.subset.all.colnames))
                        if (length(df.subset.all.dup.idx) > 0) {
                            for (d in 1:length(df.subset.all.dup.idx)) {
                                if (length(grep('\\.', colnames(df.subset.all)[df.subset.all.dup.idx[d]], value = TRUE)) == 0) {
                                    sample.name <- names(which(sapply(initial2end.idx, FUN = function(X) df.subset.all.dup.idx[d] %in% unlist(X))))
                                    colnames(df.subset.all)[df.subset.all.dup.idx[d]] <- paste(colnames(df.subset.all)[df.subset.all.dup.idx[d]], '.', sample.name, sep = '')
                                }
                            }
                        }
                        if (length(df.subset.dup.idx) > 0) {
                            for (d in 1:length(df.subset.dup.idx)) {
                                if (length(grep('\\.', colnames(df.subset)[df.subset.dup.idx[d]], value = TRUE)) == 0) {
                                    sample.name <- exp.name[i]
                                    colnames(df.subset)[df.subset.dup.idx[d]] <- paste(colnames(df.subset)[df.subset.dup.idx[d]], '.', sample.name, sep = '')
                                }
                            }
                        }
                    }
                    
                    df.subset.all <- merge(df.subset.all, df.subset, by = "row.names", all = TRUE, sort = TRUE)
                    rownames(df.subset.all) <- df.subset.all$Row.names
                    df.subset.all$Row.names <- NULL
                } else {
                    initial2end.idx <- c(initial2end.idx, NA)
                    names(initial2end.idx)[length(initial2end.idx)] <- exp.name[i]
                }
            }
            
            if (! is.null(df.subset.all)) {
                df.subset.all[is.na(df.subset.all)] <- 0
                
                row.idx <- unlist(sapply(all_gene_name, FUN = function(X) which(rownames(df.subset.all) %in% X)))
                df.subset.all <- df.subset.all[row.idx, , drop = FALSE]
                
                assign(paste('sgRNA','.',sgRNA.barcode[b],'.','df_subset_all', sep = ''), df.subset.all)
                #write.csv(df.subset.all, file = paste(current_folder,'/Result','/','Mouse','/','sgRNA','/',sgRNA.barcode[b],'/',sgRNA.barcode[b],'.QC','.csv', sep = ''), row.names = TRUE, quote = TRUE)
                
                for (i in 1:length(exp.name)) {
                    df.new <- eval(parse(text=paste('sample','.',exp.name[i],'.',sgRNA.barcode[b],'.','df_new', sep = '')))
                    
                    if (! is.null(df.new)) {
                        row.idx <- which(rownames(df.subset.all) %in% rownames(df.new))
                        col.idx <- which(colnames(df.subset.all) %in% paste(colnames(df.new),'.',exp.name[i], sep = ''))
                        
                        df.subset <- df.subset.all[row.idx, col.idx, drop = FALSE]
                        
                        write.csv(df.subset, file = paste(current_folder,'/Result','/','Mouse','/','sgRNA','/',sgRNA.barcode[b],'/','QC','/',sgRNA.barcode[b],'.QC','.',exp.name[i],'.csv', sep = ''), row.names = TRUE, quote = TRUE)
                    }
                }
            }
            
            # /// combine library and control ///
            mainDir <- paste(current_folder,'/Result','/','Mouse','/','sgRNA','/',sgRNA.barcode[b], sep = '')
            subDir <- 'control'
            if (! file.exists(file.path(mainDir, subDir))) {
                dir.create(file.path(mainDir, subDir))
            }
            
            df.subset.all.control <- eval(parse(text=paste('sample','.','control','.','df_subset_all', sep = '')))
            
            df.subset.all.new <- df.subset.all
            if (! is.null(df.subset.all.new)) {
                rownames(df.subset.all.new)[which(rownames(df.subset.all.new) == "Tead2")] = "Tead2(ENSMUSG00000030796)"
                rownames(df.subset.all.new)[which(rownames(df.subset.all.new) == "l7Rn6")] = "l7Rn6(ENSMUSG00000062797)"
                
                row.idx <- unlist(sapply(control_gene_name, FUN = function(X) which(rownames(df.subset.all.new) %in% X)))
                df.subset.all.new <- df.subset.all.new[row.idx, , drop = FALSE]
                
                df.subset.all.new.dup.idx <- which(sapply(colnames(df.subset.all.new), FUN = function(X) X %in% colnames(df.subset.all.control)))
                df.subset.all.control.dup.idx <- which(sapply(colnames(df.subset.all.control), FUN = function(X) X %in% colnames(df.subset.all.new)))
                if (length(df.subset.all.new.dup.idx) > 0) {
                    for (d in 1:length(df.subset.all.new.dup.idx)) {
                        if (length(grep('\\.', colnames(df.subset.all.new)[df.subset.all.new.dup.idx[d]], value = TRUE)) == 0) {
                            sample.name <- names(which(sapply(exp.idx, FUN = function(X) df.subset.all.new.dup.idx[d] %in% unlist(X))))
                            colnames(df.subset.all.new)[df.subset.all.new.dup.idx[d]] <- paste(colnames(df.subset.all.new)[df.subset.all.new.dup.idx[d]], '.', sample.name, sep = '')
                        }
                    }
                }
                if (length(df.subset.all.control.dup.idx) > 0) {
                    colnames(df.subset.all.control)[df.subset.all.control.dup.idx] <- paste(colnames(df.subset.all.control)[df.subset.all.control.dup.idx], '.', 'control', sep = '')
                }
            }
            
            df2control <- merge(df.subset.all.new, df.subset.all.control, by = "row.names", all = TRUE, sort = TRUE)
            #df2control <- merge(df.subset.all.control, df.subset.all.new, by = "row.names", all = TRUE, sort = TRUE)
            rownames(df2control) <- df2control$Row.names
            df2control$Row.names <- NULL
            
            all_gene_name.new <- all_gene_name
            all_gene_name.new[which(all_gene_name.new == "Tead2")] = "Tead2(ENSMUSG00000030796)"
            all_gene_name.new[which(all_gene_name.new == "l7Rn6")] = "l7Rn6(ENSMUSG00000062797)"
            
            row.idx <- unlist(sapply(unique(c(all_gene_name.new, control_gene_name)), FUN = function(X) which(rownames(df2control) %in% X)))
            df2control <- df2control[row.idx, , drop = FALSE]
            
            # /// replace NA with true value ///
            row.idx1 <- unlist(sapply(rownames(df2control), FUN = function(X) which(control_gene_name %in% X)))
            row.idx2 <- unlist(sapply(rownames(df2control), FUN = function(X) which(all_gene_name.new %in% X)))
            
            exp.idx1 <- lapply(control, FUN = function(X) which(colnames(df2control) %in% colnames(eval(parse(text=paste('sample','.','control','.',X,'.','df_new', sep = ''))))))
            names(exp.idx1) <- control
            
            exp.idx1.NA <- lapply(control, FUN = function(X) colnames(eval(parse(text=paste('sample','.','control','.',X,'.','df_new', sep = ''))))[which(! colnames(eval(parse(text=paste('sample','.','control','.',X,'.','df_new', sep = '')))) %in% colnames(df2control))])
            exp.idx1.NA <- lapply(exp.idx1.NA, FUN = function(Y) sapply(Y, FUN = function(X) unlist(strsplit(X, '\\.'))[[1]]))
            exp.idx1.NA <- lapply(1:length(control), FUN = function(X) which(colnames(df2control) %in% exp.idx1.NA[[X]]))
            names(exp.idx1.NA) <- control
            
            exp.idx1 <- lapply(1:length(control), FUN = function(X) sort(unique(c(exp.idx1[[X]], exp.idx1.NA[[X]]))))
            names(exp.idx1) <- control
            
            exp.idx2 <- lapply(exp.name, FUN = function(X) which(colnames(df2control) %in% paste(colnames(eval(parse(text=paste('sample','.',X,'.',sgRNA.barcode[b],'.','df_new', sep = '')))),'.',X, sep = '')))
            names(exp.idx2) <- paste(exp.name, '.', sgRNA.barcode[b], sep = '')
            
            exp.idx12 <- c(exp.idx1, exp.idx2)
            
            for (k in 1:length(exp.idx12)) {
                col.idx <- unlist(exp.idx12[k])
                names(col.idx) <- colnames(df2control)[col.idx]
                
                match.idx <- grep('\\.', names(col.idx), value = FALSE)
                match.name <- grep('\\.', names(col.idx), value = TRUE)
                
                if (length(match.name) > 0) {
                    match.idx <- match.idx[which(sapply(match.name, FUN = function(X) unlist(strsplit(X, '\\.'))[[2]] %in% sample.reorder))]
                    if (length(match.idx) > 0) {
                        names(col.idx)[match.idx] = sapply(names(col.idx)[match.idx], FUN = function(X) unlist(strsplit(X, '\\.'))[[1]])
                    }
                }
                
                if (names(exp.idx12[k]) %in% control) {
                    sample.name <- "control"
                    sample_file <- names(exp.idx12[k])
                    row.idx <- row.idx1
                } else {
                    sample.name <- unlist(strsplit(names(exp.idx12[k]), '\\.'))[1]
                    sample_file <- unlist(strsplit(names(exp.idx12[k]), '\\.'))[2]
                    row.idx <- row.idx2
                }
                
                df <- eval(parse(text=paste('sample','.',sample.name,'.',sample_file,'.','df', sep = '')))
                
                if (! is.null(df)) {
                    rownames(df)[which(rownames(df) == "Tead2")] = "Tead2(ENSMUSG00000030796)"
                    rownames(df)[which(rownames(df) == "l7Rn6")] = "l7Rn6(ENSMUSG00000062797)"
                    
                    NA.idx <- which(! rownames(df2control) %in% rownames(df))
                    df2control[NA.idx, col.idx][is.na(df2control[NA.idx, col.idx])] <- 0
                    
                    colnames(df) = sapply(colnames(df), FUN = function(X) unlist(strsplit(X, '\\.'))[[1]])
                    names(col.idx) = sapply(names(col.idx), FUN = function(X) unlist(strsplit(X, '\\.'))[[1]])
                    
                    df2control[names(row.idx), col.idx] = df[row.idx, names(col.idx)]
                }
            }
            
            #df2control[is.na(df2control)] <- 0
            
            row.idx <- which(rownames(df2control) %in% intersect(all_gene_name.new, control_gene_name))
            df2control <- df2control[row.idx, , drop = FALSE]
            
            #assign(paste('sgRNA','.',sgRNA.barcode[b],'.','df2control', sep = ''), df2control)
            #write.csv(df2control, file = paste(current_folder,'/Result','/','Mouse','/','sgRNA','/',sgRNA.barcode[b],'/',sgRNA.barcode[b],'.control','.csv', sep = ''), row.names = TRUE, quote = TRUE)
            
            # /// split into individual library and control ///
            for (k in 1:length(exp.idx12)) {
                col.idx <- unlist(exp.idx12[k])
                names(col.idx) <- colnames(df2control)[col.idx]
                
                sample.name <- unlist(strsplit(names(exp.idx12[k]), '\\.'))[1]
                
                if (length(col.idx) > 0) {
                    df2control.subset <- df2control[, col.idx, drop = FALSE]
                    #assign(paste('sgRNA','.',sgRNA.barcode[b],'.',sample.name,'.','df2control_subset', sep = ''), df2control.subset)
                    write.csv(df2control.subset, file = paste(current_folder,'/Result','/','Mouse','/','sgRNA','/',sgRNA.barcode[b],'/','control','/',sgRNA.barcode[b],'.control','.',sample.name,'.csv', sep = ''), row.names = TRUE, quote = TRUE)
                }
            }
        }
    }
}


# --- sample2number (sample2cellnumber & sample2genenumber) -----------------------------------------------------------------------------------------
mainDir <- paste(current_folder,'/Result','/','Mouse','/','library', sep = '')
subDir <- 'sample2number'
if (! file.exists(file.path(mainDir, subDir))) {
    dir.create(file.path(mainDir, subDir))
}

sample.type <- c("sample", "sample.merged", "sample.reorder", "d3", "d5", "d7")

for (n in 1:length(sample.type)) {
    
    mainDir <- paste(current_folder,'/Result','/','Mouse','/','library','/','sample2number', sep = '')
    subDir <- sample.type[n]
    if (! file.exists(file.path(mainDir, subDir))) {
        dir.create(file.path(mainDir, subDir))
    }
    
    if (sample.type[n] == "sample") {
        sample.list = sample
        prename = "sample_single"
        width.size = length(sample.list) + 7
    } else if (sample.type[n] == "sample.merged") {
        sample.list = sample.merged
        prename = "sample_merged"
        width.size = length(sample.list) + 7
    } else if (sample.type[n] == "sample.reorder") {
        sample.list = sample.reorder
        prename = "sample_all"
        width.size = length(sample.list) + 11
    } else if (sample.type[n] == "d3") {
        sample.list = sample.reorder[which(sample.reorder %in% c(d3, sample.merged[which(sample.merged %in% c("d3_p1_1_merged","d3"))]))]
        prename = "sample_d3"
        width.size = length(sample.list) + 7
    } else if (sample.type[n] == "d5") {
        sample.list = sample.reorder[which(sample.reorder %in% c(d5, sample.merged[which(sample.merged %in% c("d5_p1_merged","d5"))]))]
        prename = "sample_d5"
        width.size = length(sample.list) + 7
    } else if (sample.type[n] == "d7") {
        sample.list = sample.reorder[which(sample.reorder %in% c(d7, sample.merged[which(sample.merged %in% c("d7_p1_merged","d7_p2","d7"))]))]
        prename = "sample_d7"
        width.size = length(sample.list) + 7
    }
    
    # --- different gene cutoff -------------------------------------------------------------------------------------------------------------------------
    for (s in 1:length(sample.list)) {
        
        if (sample.list[s] == "control") {
            file_num = length(control)
        } else {
            file_num = length(sgRNA.barcode)
        }
        
        xy.all <- c()
        
        for (i in 1:file_num) {
            
            if (sample.list[s] == "control") {
                sample_file = control[i]
            } else {
                sample_file = sgRNA.barcode[i]
            }
            
            cell.num.new.all <- c()
            gene.num.new.all <- c()
            
            for (g in 1:length(gene_cutoff)) {
                cell.num.new <- eval(parse(text=paste('sample','.',sample.list[s],'.',sample_file,'.',gene_cutoff[g],'.','cell_num_new', sep = '')))
                gene.num.new <- eval(parse(text=paste('sample','.',sample.list[s],'.',sample_file,'.',gene_cutoff[g],'.','gene_num_new', sep = '')))
                
                if (g == 1) {
                    cell.num.new.all <- c(cell.num.new.all, cell.num.new)
                } else {
                    if (cell.num.new > cell.num.new.all[length(cell.num.new.all)]) {
                        cell.num.new.all <- c(cell.num.new.all, 0)
                    } else {
                        cell.num.new.all <- c(cell.num.new.all, cell.num.new)
                    }
                }
                
                if (g == 1) {
                    gene.num.new.all <- c(gene.num.new.all, gene.num.new)
                } else {
                    if (gene.num.new > gene.num.new.all[length(gene.num.new.all)]) {
                        gene.num.new.all <- c(gene.num.new.all, 0)
                    } else {
                        gene.num.new.all <- c(gene.num.new.all, gene.num.new)
                    }
                }
            }
            
            x = gene_cutoff
            y = gene.num.new.all
            
            xy <- data.frame(gene_cutoff = x, gene_num = y, group = sample_file)
            
            xy.all <- rbind(xy.all, xy)
            
        }
        
        
        # https://dannagifford.com/2020/04/14/r-and-ggplot2-make-zero-print-as-0/
        # Make zeros print as "0" always
        library(stringr)
        prettyZero <- function(l){
            max.decimals = max(nchar(str_extract(l, "\\.[0-9]+")), na.rm = T)-1
            lnew = formatC(l, replace.zero = T, zero.print = "0",
                           digits = max.decimals, format = "f", preserve.width=T)
            return(lnew)
        }
        
        library(ggplot2)
        
        xbreaks = c(0, gene_cutoff)
        
        ybreaks = seq(0, max(xy.all[,"gene_num"]), 1000)
        ybreaks = unique(c(ybreaks, ybreaks[length(ybreaks)]+1000))
        
        p <- ggplot(data = xy.all, aes(x = gene_cutoff, y = gene_num, group = group)) + 
            geom_line(aes(color = group)) + 
            geom_point(aes(color = group))
        p <- p + 
            theme(panel.border = element_blank(), axis.line = element_line()) + 
            theme(panel.background = element_rect(fill = "white", colour = "white")) + 
            theme(panel.grid = element_line(colour = "grey80"), panel.grid.minor = element_blank())
        p <- p + 
            ggtitle(sample.list[s]) + theme(plot.title = element_text(hjust = 0.5, size = 24, family = "Helvetica")) + 
            xlab("gene cutoff") + theme(axis.title.x = element_text(size = 20, family = "Helvetica"), axis.text.x = element_text(size = 16, family = "Helvetica")) + 
            ylab("gene number") + theme(axis.title.y = element_text(size = 20, family = "Helvetica"), axis.text.y = element_text(size = 16, family = "Helvetica")) + 
            labs(col = "sgRNA") + theme(legend.title = element_text(size = 20, family = "Helvetica"), legend.text = element_text(size = 12, family = "Courier"))
        p <- p + 
            scale_x_continuous(limits = c(0, max(xbreaks)), expand = c(0, 0), breaks = xbreaks, labels = prettyZero) + 
            scale_y_continuous(limits = c(0, max(ybreaks)), expand = c(0, 0), breaks = ybreaks)
        
        if (sample.list[s] %in% sample.merged) {
            png(file = paste(current_folder,'/Result','/','Mouse','/','library','/','merged','/',sample.list[s],'/',sample.list[s],'.','gene_cutoff2gene_number','.png', sep = ''), width = 10, height = 6, units = 'in', res = 600)
        } else {
            png(file = paste(current_folder,'/Result','/','Mouse','/','library','/',sample.list[s],'/',sample.list[s],'.','gene_cutoff2gene_number','.png', sep = ''), width = 10, height = 6, units = 'in', res = 600)
        }
        print(p)
        dev.off()
        
    }
    
    # --- dynamics --------------------------------------------------------------------------------------------------------------------------------------
    for (g in 1:length(gene_cutoff)) {
        
        file_num = length(sgRNA.barcode)
        
        xy.all <- c()
        xz.all <- c()
        
        for (i in 1:file_num) {
            
            sample_file = sgRNA.barcode[i]
            
            cell.num.new.all <- c()
            gene.num.new.all <- c()
            
            for (s in 1:length(sample.list)) {
                
                if (sample.list[s] != "control") {
                    
                    cell.num.new <- eval(parse(text=paste('sample','.',sample.list[s],'.',sample_file,'.',gene_cutoff[g],'.','cell_num_new', sep = '')))
                    gene.num.new <- eval(parse(text=paste('sample','.',sample.list[s],'.',sample_file,'.',gene_cutoff[g],'.','gene_num_new', sep = '')))
                    
                    cell.num.new.all <- c(cell.num.new.all, cell.num.new)
                    gene.num.new.all <- c(gene.num.new.all, gene.num.new)
                    
                }
            }
            
            x = sample.list[sample.list != "control"]
            rename.idx <- which(x == "20190919_d3_p1_1")
            if (length(rename.idx) > 0) x[rename.idx] = "19919_d3p11"
            rename.idx <- which(x == "poor_d3_p1_1")
            if (length(rename.idx) > 0) x[rename.idx] = "poor_d3p11"
            rename.idx <- which(x == "poor_d5_p1")
            if (length(rename.idx) > 0) x[rename.idx] = "poor_d5p1"
            rename.idx <- which(x == "poor_d7_p1")
            if (length(rename.idx) > 0) x[rename.idx] = "poor_d7p1"
            rename.idx <- which(x == "d3_p1_1_merged")
            if (length(rename.idx) > 0) x[rename.idx] = "d3p11_merged"
            rename.idx <- which(x == "d5_p1_merged")
            if (length(rename.idx) > 0) x[rename.idx] = "d5p1_merged"
            rename.idx <- which(x == "d7_p1_merged")
            if (length(rename.idx) > 0) x[rename.idx] = "d7p1_merged"
            
            y = gene.num.new.all
            z = cell.num.new.all
            
            xy <- data.frame(library = x, gene_num = y, group = sample_file)
            xz <- data.frame(library = x, cell_num = z, group = sample_file)
            
            xy.all <- rbind(xy.all, xy)
            xz.all <- rbind(xz.all, xz)
            
        }
        
        xy.all[,1] <- factor(xy.all[,1], levels = unique(xy.all[,1]))
        xz.all[,1] <- factor(xz.all[,1], levels = unique(xz.all[,1]))
        
        
        library(ggplot2)
        
        # /// gene dynamics ///
        ybreaks = seq(0, max(xy.all[,"gene_num"]), 1000)
        ybreaks = unique(c(ybreaks, ybreaks[length(ybreaks)]+1000))
        
        p <- ggplot(data = xy.all, aes(x = library, y = gene_num, group = group)) + 
            geom_line(aes(color = group)) + 
            geom_point(aes(color = group))
        p <- p + 
            theme(panel.border = element_blank(), axis.line = element_line()) + 
            theme(panel.background = element_rect(fill = "white", colour = "white")) + 
            theme(panel.grid = element_line(colour = "grey80"), panel.grid.minor = element_blank())
        p <- p + 
            ggtitle(paste("sc26"," (","gene cutoff"," = ",gene_cutoff[g],")", sep = '')) + theme(plot.title = element_text(hjust = 0.5, size = 24, family = "Helvetica")) + 
            xlab("sample library") + theme(axis.title.x = element_text(size = 20, family = "Helvetica"), axis.text.x = element_text(size = 16, family = "Helvetica")) + 
            ylab("gene number") + theme(axis.title.y = element_text(size = 20, family = "Helvetica"), axis.text.y = element_text(size = 16, family = "Helvetica")) + 
            labs(col = "sgRNA") + theme(legend.title = element_text(size = 20, family = "Helvetica"), legend.text = element_text(size = 12, family = "Courier"))
        p <- p + 
            scale_y_continuous(limits = c(0, max(ybreaks)), expand = c(0, 0), breaks = ybreaks)
        
        png(file = paste(current_folder,'/Result','/','Mouse','/','library','/','sample2number','/',sample.type[n],'/',prename,'2','gene_number','.','gene_cutoff','=',gene_cutoff[g],'.png', sep = ''), width = width.size, height = 6, units = 'in', res = 600)
        print(p)
        dev.off()
        
        # /// cell dynamics ///
        if (gene_cutoff[g] == 0.1) {
            
            ybreaks = seq(0, max(xz.all[,"cell_num"]), 1000)
            ybreaks = unique(c(ybreaks, ybreaks[length(ybreaks)]+1000))
            
            p <- ggplot(data = xz.all, aes(x = library, y = cell_num, group = group)) + 
                geom_line(aes(color = group)) + 
                geom_point(aes(color = group))
            p <- p + 
                theme(panel.border = element_blank(), axis.line = element_line()) + 
                theme(panel.background = element_rect(fill = "white", colour = "white")) + 
                theme(panel.grid = element_line(colour = "grey80"), panel.grid.minor = element_blank())
            p <- p + 
                ggtitle("sc26") + theme(plot.title = element_text(hjust = 0.5, size = 24, family = "Helvetica")) + 
                xlab("sample library") + theme(axis.title.x = element_text(size = 20, family = "Helvetica"), axis.text.x = element_text(size = 16, family = "Helvetica")) + 
                ylab("cell number") + theme(axis.title.y = element_text(size = 20, family = "Helvetica"), axis.text.y = element_text(size = 16, family = "Helvetica")) + 
                labs(col = "sgRNA") + theme(legend.title = element_text(size = 20, family = "Helvetica"), legend.text = element_text(size = 12, family = "Courier"))
            p <- p + 
                scale_y_continuous(limits = c(0, max(ybreaks)), expand = c(0, 0), breaks = ybreaks)
            
            png(file = paste(current_folder,'/Result','/','Mouse','/','library','/','sample2number','/',sample.type[n],'/',prename,'2','cell_number','.png', sep = ''), width = width.size, height = 6, units = 'in', res = 600)
            print(p)
            dev.off()
            
        }
    }
}

