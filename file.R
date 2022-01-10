#!/usr/bin/env Rscript

rm(list = ls()) # remove all variables in workspace

# ---------------------------------------------------------------------------------------------------------------------------------------------------
#current_folder = "~/Project/Eczema"
current_folder = "/Users/mingtsung/Desktop/Project/Eczema"
setwd(current_folder)

project.name = unlist(strsplit(current_folder, '/'))
project.name = project.name[length(project.name)]

#install.packages("mgsub") # mgsub: Safe, Multiple, Simultaneous String Substitution
#install.packages("flextable") # flextable: Functions for Tabular Reporting
#install.packages("gtools") # gtools: Various R Programming Tools
#install.packages("stringr") # stringr: Simple, Consistent Wrappers for Common String Operations
#install.packages("dendextend") # dendextend: Extending 'dendrogram' Functionality in R
#install.packages("gplots") # gplots: Various R Programming Tools for Plotting Data
#install.packages("pheatmap") # pheatmap: Pretty Heatmaps
#install.packages("Rtsne") # Rtsne: T-Distributed Stochastic Neighbor Embedding using a Barnes-Hut Implementation
#install.packages("coin") # coin: Conditional Inference Procedures in a Permutation Test Framework
#install.packages("openxlsx") # openxlsx: Read, Write and Edit xlsx Files
#install.packages("readxl") # readxl: Read Excel Files
#install.packages("Seurat") # Seurat: Tools for Single Cell Genomics (https://satijalab.org/seurat/vignettes.html)

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("ComplexHeatmap")
#BiocManager::install("edgeR")
#BiocManager::install("DESeq2")

#suppressMessages(library(mgsub)) # mgsub()
#suppressMessages(library(flextable)) # italic()
suppressMessages(library(gtools)) # mixedsort(), mixedorder()
suppressMessages(library(stringr)) # str_to_title(), str_to_upper()
suppressMessages(library(dendextend)) # Beautiful dendrogram visualizations in R: http://www.sthda.com/english/wiki/beautiful-dendrogram-visualizations-in-r-5-must-known-methods-unsupervised-machine-learning
suppressMessages(library(gplots)) # heatmap.2()
suppressMessages(library(pheatmap)) # pheatmap()
suppressMessages(library(ComplexHeatmap)) # Heatmap()
suppressMessages(library(Rtsne)) # Rtsne()
suppressMessages(library(coin)) # independence_test()
suppressMessages(library(openxlsx)) # read.xlsx()
suppressMessages(library(readxl)) # read_excel()
suppressMessages(library(Seurat)) # Read10X()
#suppressMessages(library(edgeR))
#suppressMessages(library(DESeq2))

mainDir <- current_folder
subDir <- "Result"
if (! file.exists(file.path(mainDir, subDir))) {
    dir.create(file.path(mainDir, subDir))
}


# --- Source function -------------------------------------------------------------------------------------------------------------------------------
source('./Script/function/Monocle3/Monocle3_analysis.R')
source('./Script/function/Seurat/Seurat_analysis.R')
source('./Script/function/Seurat/StackedVlnPlot.R')
source('./Script/function/plot/volcano_plot.R')


# --- Read in "gene(s) of interest" and "signature gene(s)" -----------------------------------------------------------------------------------------
# /// query gene(s) ///
geneSets_of_interest.human <- list()

# *** geneSet1_of_interest ***
geneSet1_of_interest <- c("KRT5")

geneSet1_of_interest.human <- str_to_upper(geneSet1_of_interest)

geneSet1_of_interest.human <- list(unique(geneSet1_of_interest.human))
names(geneSet1_of_interest.human) = "geneSet1_of_interest"
geneSets_of_interest.human <- c(geneSets_of_interest.human, geneSet1_of_interest.human)

genes_of_interest <- geneSets_of_interest.human
KRT5.genes_of_interest = genes_of_interest


# --- Read in GSA -----------------------------------------------------------------------------------------------------------------------------------
Proj_ID <- c("NG-27918")
Dataset <- c(Proj_ID)

# /// cutoff ///
pval_cutoff = 10^-3
padj_cutoff = 10^-3
FDR_cutoff = padj_cutoff
qval_cutoff = padj_cutoff

FC_cutoff = 2
log2FC_cutoff = log(FC_cutoff, 2)
logFC_cutoff = log(FC_cutoff, 2)

#for (i in 1:length(Proj_ID)) {
for (i in 1) {
    
    mainDir <- paste(current_folder,'/Result', sep = '')
    subDir <- Proj_ID[i]
    if (! file.exists(file.path(mainDir, subDir))) {
        dir.create(file.path(mainDir, subDir))
    }
    
    if (Proj_ID[i] == "NG-27918") {
        sample.metadata <- read.xlsx(paste(current_folder,'/Data','/','LeoPharma','/','bulk_RNAseq','/','Bulk-RNAseq metadata_2Nov2021','.xlsx', sep = ''), startRow = 3, check.names = FALSE, sep.names = " ")
        
        # +++ swap error label to correct label (2021.11.30) ++++++++++++++++++++++++++++++++++++++++++++++
        # *** swap unstimulated and stimulated labels of T12 **********************************************
        sample.metadata[,6] = gsub("Time 12h, unstimulated","Time 12h, treat with unstimulated",sample.metadata[,6])
        sample.metadata[,6] = gsub("Time 12h, stimulated","Time 12h, treat with stimulated",sample.metadata[,6])
        sample.metadata[,6] = gsub("Time 12h, treat with unstimulated","Time 12h, stimulated",sample.metadata[,6])
        sample.metadata[,6] = gsub("Time 12h, treat with stimulated","Time 12h, unstimulated",sample.metadata[,6])
        
        T12.unstimulated.idx <- grep("Time 12h, unstimulated", sample.metadata[,6], value = FALSE)
        T12.stimulated.idx <- grep("Time 12h, stimulated", sample.metadata[,6], value = FALSE)
        
        row.idx <- c(1:(min(T12.unstimulated.idx, T12.stimulated.idx) - 1), 
                     T12.unstimulated.idx, 
                     T12.stimulated.idx, 
                     (max(T12.unstimulated.idx, T12.stimulated.idx) + 1):nrow(sample.metadata))
        
        sample.metadata <- sample.metadata[row.idx, , drop = FALSE]
        # *************************************************************************************************
        
        sample.id <- sample.metadata[,4]
        sample.name <- sample.metadata[,5]
        sample.annotation <- sample.metadata[,6]
        
        sample.name <- sapply(sample.name, FUN = function(X) gsub("\\.","-",gsub("-","_",unlist(strsplit(X, '_'))[[1]])))
        names(sample.id) = sample.name
        names(sample.name) = sample.name
        names(sample.annotation) = sample.name
        
        sample.time <- sapply(sample.name, FUN = function(X) unlist(strsplit(X, '_'))[[1]])
        
        sample.treatment <- sapply(sapply(sample.name, FUN = function(X) unlist(strsplit(X, '_'))[[2]]), FUN = function(Y) unlist(strsplit(Y, '-'))[[1]])
        sample.treatment[which(sample.treatment == 1)] = "unstimulated"
        sample.treatment[which(sample.treatment == 2)] = "stimulated"
        
        # *** swap unstimulated and stimulated labels of T12 **********************************************
        T12_1.idx <- grep("T12_1", names(sample.treatment), value = FALSE)
        T12_2.idx <- grep("T12_2", names(sample.treatment), value = FALSE)
        sample.treatment[T12_2.idx] = "unstimulated"
        sample.treatment[T12_1.idx] = "stimulated"
        # *************************************************************************************************
        
        sample.bio_replicate <- sapply(sapply(sample.name, FUN = function(X) unlist(strsplit(X, '_'))[[2]]), FUN = function(Y) unlist(strsplit(Y, '-'))[[2]])
        sample.bio_replicate[which(sample.bio_replicate == 1)] = "Sample 1"
        sample.bio_replicate[which(sample.bio_replicate == 2)] = "Sample 2"
        sample.bio_replicate[which(sample.bio_replicate == 3)] = "Sample 3"
        
        sample.id.all <- sample.id
        sample.name.all <- sample.name
        sample.annotation.all <- sample.annotation
        sample.time.all <- sample.time
        sample.treatment.all <- sample.treatment
        sample.bio_replicate.all <- sample.bio_replicate
    }
    
    sample <- sample.name.all
    
    sample.id.new <- c()
    sample.name.new <- c()
    sample.annotation.new <- c()
    sample.time.new <- c()
    sample.treatment.new <- c()
    sample.bio_replicate.new <- c()
    
    cell.num.table <- c()
    cell.num.new.table <- c()
    
    gene.num.table <- c()
    gene.num.new.table <- c()
    
    df.all <- c()
    
    for (j in 1:length(sample.name.all)) {
        
        gc() # Garbage Collection
        
        print(j)
        
        sample.id = sample.id.all[j]
        sample.name = sample.name.all[j]
        sample.annotation = sample.annotation.all[j]
        sample.time = sample.time.all[j]
        sample.treatment = sample.treatment.all[j]
        sample.bio_replicate = sample.bio_replicate.all[j]
        
        sample.name.origin <- sample.name
        print(sample.name)
        
        project.name = 'LeoPharma'
        project.subname = paste(Proj_ID[i],'.',gsub("\\.| ","_",sample.name), sep = '')
        project.input = paste(project.name,'.',project.subname, sep = '')
        
        # --- read from count matrix ------------------------------------------------------------------------------------------------------------------------
        if (sample.name == "T2_2-1") {
            sample.run = 2
        } else if (sample.name == "T12_1-1") {
            sample.run = 2
        } else if (sample.name == "T12_1-3") {
            sample.run = 2
        } else {
            sample.run = 1
        }
        
        for (sR in 1:sample.run) {
            sample.name <- sample.name.origin
            
            if (sample.name == "T2_2-1") {
                if (sR == 1) {
                    sample.name = paste(sample.name,"-7541", sep = '')
                } else if (sR == 2) {
                    sample.name = paste(sample.name,"-7550", sep = '')
                }
            } else if (sample.name == "T12_1-1") {
                if (sR == 1) {
                    sample.name = paste(sample.name,"-7541", sep = '')
                } else if (sR == 2) {
                    sample.name = paste(sample.name,"-7550", sep = '')
                }
            } else if (sample.name == "T12_1-3") {
                if (sR == 1) {
                    sample.name = paste(sample.name,"-7541", sep = '')
                } else if (sR == 2) {
                    sample.name = paste(sample.name,"-7550", sep = '')
                }
            }
            names(sample.id) = sample.name
            names(sample.name) = sample.name
            names(sample.annotation) = sample.name
            names(sample.time) = sample.name
            names(sample.treatment) = sample.name
            names(sample.bio_replicate) = sample.name
            
            sample.id.new <- c(sample.id.new, sample.id)
            sample.name.new <- c(sample.name.new, sample.name)
            sample.annotation.new <- c(sample.annotation.new, sample.annotation)
            sample.time.new <- c(sample.time.new, sample.time)
            sample.treatment.new <- c(sample.treatment.new, sample.treatment)
            sample.bio_replicate.new <- c(sample.bio_replicate.new, sample.bio_replicate)
            
            
            count_matrix.name = paste(sample.name,'_','reverse','_','counts_new', sep = '')
            matrix.data <- read.csv(paste(current_folder,'/Data','/','LeoPharma','/','bulk_2021','/',paste(Proj_ID[i],'_',sample.time, sep = ''),'/',paste(sample.time,'_','count', sep = ''),'/',count_matrix.name,'.csv', sep = ''), sep = ',', header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
            
            genes.idx <- which(colnames(matrix.data) == "genes")
            if (length(genes.idx) == 0) {
                genes.NA = 1
                matrix.data[,"genes"] <- rownames(matrix.data)
            } else {
                genes.NA = 0
                NA.idx <- sort(unique(c(which(matrix.data[,"genes"] == ""), grep('Row', matrix.data[,"genes"], value = FALSE))))
                matrix.data[NA.idx,"genes"] <- rownames(matrix.data)[NA.idx]
            }
            genes.idx <- which(colnames(matrix.data) == "genes")
            if (genes.idx != 1) {
                col.idx <- c(genes.idx, 1:(ncol(matrix.data)-1))
                matrix.data <- matrix.data[, col.idx]
            }
            genes.idx <- which(colnames(matrix.data) == "genes")
            
            gene_names.idx <- which(colnames(matrix.data) == "gene_names")
            if (length(gene_names.idx) == 0) {
                gene_names.NA = 1
                matrix.data[,"gene_names"] <- rownames(matrix.data)
            } else {
                gene_names.NA = 0
                NA.idx <- sort(unique(c(which(matrix.data[,"gene_names"] == ""), grep('Row', matrix.data[,"gene_names"], value = FALSE))))
                matrix.data[NA.idx,"gene_names"] <- rownames(matrix.data)[NA.idx]
            }
            gene_names.idx <- which(colnames(matrix.data) == "gene_names")
            if (gene_names.idx != 2) {
                col.idx <- c(1, gene_names.idx, 2:(ncol(matrix.data)-1))
                matrix.data <- matrix.data[, col.idx]
            }
            gene_names.idx <- which(colnames(matrix.data) == "gene_names")
            
            # /// gene: replace empty gene_names with ID and append ID to duplicated gene_names ///
            gene.dup.idx <- which(duplicated(matrix.data[,"gene_names"]))
            if (length(gene.dup.idx) > 0) {
                for (d in 1:length(gene.dup.idx)) {
                    gene_name <- matrix.data[which(matrix.data[,"gene_names"] == matrix.data[gene.dup.idx[d],"gene_names"]),"gene_names"]
                    gene_name <- unique(gene_name)
                    
                    if (length(grep(')$', matrix.data[which(matrix.data[,"gene_names"] == gene_name),"gene_names"], value = TRUE)) == 0) {
                        matrix.data[which(matrix.data[,"gene_names"] == gene_name),"gene_names"] <- paste(matrix.data[which(matrix.data[,"gene_names"] == gene_name),"gene_names"], "(", matrix.data[which(matrix.data[,"gene_names"] == gene_name),"genes"], ")", sep = '')
                    }
                }
            }
            
            features.id <- matrix.data[,"genes"]
            features.name <- matrix.data[,"gene_names"]
            
            # /// rownames: gene_name ///
            matrix.data.row_name <- matrix.data
            rownames(matrix.data.row_name) <- matrix.data[,"gene_names"]
            # /// rownames: gene_id ///
            matrix.data.row_id <- matrix.data
            rownames(matrix.data.row_id) <- matrix.data[,"genes"]
            
            
            id2id <- features.id
            id2name <- features.name
            
            id.table <- data.frame(id2id, id2name)
            colnames(id.table) <- c("gene_id","gene_name")
            
            # /// gene: replace empty Symbol with ID and append ID to duplicated Symbol ///
            gene.dup.idx <- which(duplicated(id2name))
            if (length(gene.dup.idx) > 0) {
                for (d in 1:length(gene.dup.idx)) {
                    gene_name <- id2name[gene.dup.idx[d]]
                    gene_name <- unique(gene_name)
                    
                    if (length(grep(')$', id2name[which(id2name == gene_name)], value = TRUE)) == 0) {
                        id2name[which(id2name == gene_name)] <- paste(id2name[which(id2name == gene_name)], "(", id2id[which(id2name == gene_name)], ")", sep = '')
                    }
                }
            }
            
            id.table.new <- data.frame(id2id, id2name)
            colnames(id.table.new) <- c("gene_id","gene_name")
            
            matrix.data.origin <- matrix.data.row_name
            
            
            all_gene_id <- rownames(matrix.data.row_id)
            all_gene_name <- rownames(matrix.data.row_name)
            
            DATA <- matrix.data.row_name # read from count matrix
            
            genes.idx <- which(colnames(DATA) == "genes")
            gene_names.idx <- which(colnames(DATA) == "gene_names")
            
            DATA <- DATA[, -c(genes.idx, gene_names.idx), drop = FALSE]
            
            cell.barcode <- colnames(DATA)
            
            # /// duplicated cell barcode in the same file ///
            dup.idx <- which(duplicated(cell.barcode))
            if (length(dup.idx) > 0) {
                cell.barcode.dup <- cell.barcode[dup.idx]
                
                if (length(cell.barcode.dup) > 0) {
                    mainDir <- paste(current_folder,'/Result','/',Proj_ID[i], sep = '')
                    subDir <- 'cell_barcodes'
                    if (! file.exists(file.path(mainDir, subDir))) {
                        dir.create(file.path(mainDir, subDir))
                    }
                    
                    if (length(cell.barcode.dup) == 1) {
                        cat(length(cell.barcode.dup),' duplicated cell barcode in "',sample.name,'" dataset: ','\n', sep = '')
                        write.table(cell.barcode.dup, file = paste(current_folder,'/Result','/',Proj_ID[i],'/','cell_barcodes','/','duplicated_barcode (',sample.name,').csv', sep = ''), sep = ',', row.names = FALSE, col.names = FALSE, quote = FALSE)
                    } else if (length(cell.barcode.dup) > 1) {
                        cat(length(cell.barcode.dup),' duplicated cell barcodes in "',sample.name,'" dataset: ','\n', sep = '')
                        write.table(cell.barcode.dup, file = paste(current_folder,'/Result','/',Proj_ID[i],'/','cell_barcodes','/','duplicated_barcodes (',sample.name,').csv', sep = ''), sep = ',', row.names = FALSE, col.names = FALSE, quote = FALSE)
                    }
                    for (d in 1:length(cell.barcode.dup)) {
                        cat('"',cell.barcode.dup[d],'" was found in "',sample.name,'".','\n', sep = '')
                    }
                    cat('\n')
                    
                    for (d in 1:length(cell.barcode.dup)) {
                        col.idx <- which(cell.barcode == cell.barcode.dup[d])
                        DATA[, col.idx] = rowSums(DATA[, col.idx], na.rm = TRUE)
                    }
                    DATA <- DATA[, -dup.idx, drop = FALSE]
                } else {
                    cat('\n')
                }
            } else {
                #cat('\n')
            }
            
            DATA.origin <- DATA
            
            if (! Proj_ID[i] == "GSE147424") {
                colnames(DATA) <- paste(sample.name,'_',colnames(DATA), sep = '')
            }
            
            cell.num = ncol(DATA)
            names(cell.num) <- sample.name
            cell.num.table <- c(cell.num.table, cell.num)
            
            gene.num = nrow(DATA)
            names(gene.num) <- sample.name
            gene.num.table <- c(gene.num.table, gene.num)
            
            df <- DATA
            #write.csv(df, file = paste(current_folder,'/Result','/',Proj_ID[i],'/',sample.name,'.','count_matrix','.csv', sep = ''), row.names = TRUE, quote = TRUE)
            
            
            # --- Data filtering --------------------------------------------------------------------------------------------------------------------------------
            df_qc <- df
            df_qc.normalized <- df
            
            df.new <- df_qc
            df.new.normalized <- df_qc.normalized
            
            assign(paste('sample','.',gsub("\\.| |-","_",sample.name),'.','df', sep = ''), df)
            assign(paste('sample','.',gsub("\\.| |-","_",sample.name),'.','df_new', sep = ''), df.new)
            #assign(paste('sample','.',gsub("\\.| |-","_",sample.name),'.','df_new_normalized', sep = ''), df.new.normalized)
            
            if (! is.null(df.new)) {
                cell.num.new = ncol(df.new)
                gene.num.new = nrow(df.new)
                
                df.all <- merge(df.all, df.new, by = "row.names", all = TRUE, sort = TRUE)
                rownames(df.all) <- df.all$Row.names
                df.all$Row.names <- NULL
                
                assign(paste('sample','.',gsub("\\.| |-","_",sample.name),'.','df_all', sep = ''), df.all)
            } else {
                cell.num.new = 0
                gene.num.new = 0
            }
            
            names(cell.num.new) <- sample.name
            cell.num.new.table <- c(cell.num.new.table, cell.num.new)
            
            names(gene.num.new) <- sample.name
            gene.num.new.table <- c(gene.num.new.table, gene.num.new)
            
            
            # --- cell_metadata ---------------------------------------------------------------------------------------------------------------------------------
            DATA.sample <- sapply(colnames(DATA), FUN = function(X) paste(unlist(strsplit(X, '_'))[[1]],'_',unlist(strsplit(X, '_'))[[2]], sep = ''))
            DATA.barcode <- sapply(colnames(DATA), FUN = function(X) unlist(strsplit(X, '_'))[[3]])
            
            if (DATA.sample == "T2_2-1-7541" | DATA.sample == "T2_2-1-7550") {
                DATA.sample <- "T2_2-1"
            } else if (DATA.sample == "T12_1-1-7541" | DATA.sample == "T12_1-1-7550") {
                DATA.sample <- "T12_1-1"
            } else if (DATA.sample == "T12_1-3-7541" | DATA.sample == "T12_1-3-7550") {
                DATA.sample <- "T12_1-3"
            }
            
            exp.name <- sample
            
            exp.idx <- lapply(exp.name, FUN = function(X) which(DATA.sample %in% X))
            names(exp.idx) <- exp.name
            
            exp.num <- sapply(exp.idx, FUN = function(X) length(X))
            
            cell_metadata <- data.frame("sample.cell_barcode" = colnames(DATA)[unlist(exp.idx)], 
                                        "sample" = DATA.sample[unlist(exp.idx)], 
                                        "cell_barcode" = DATA.barcode[unlist(exp.idx)], 
                                        stringsAsFactors = FALSE)
            
            cell_metadata[,"sample"] <- factor(cell_metadata[,"sample"], levels = unique(cell_metadata[,"sample"]))
            
            cell_metadata[,"sample_id"] <- sapply(cell_metadata[,"sample"], FUN = function(X) sample.id.all[which(names(sample.id.all) %in% X)])
            cell_metadata[,"sample_annotation"] <- sapply(cell_metadata[,"sample"], FUN = function(X) sample.annotation.all[which(names(sample.annotation.all) %in% X)])
            cell_metadata[,"sample_time"] <- sapply(cell_metadata[,"sample"], FUN = function(X) sample.time.all[which(names(sample.time.all) %in% X)])
            cell_metadata[,"sample_treatment"] <- sapply(cell_metadata[,"sample"], FUN = function(X) sample.treatment.all[which(names(sample.treatment.all) %in% X)])
            cell_metadata[,"sample_time_treatment"] <- sapply(cell_metadata[,"sample"], FUN = function(X) paste(sample.time.all[which(names(sample.time.all) %in% X)],'_',sample.treatment.all[which(names(sample.treatment.all) %in% X)], sep = ''))
            cell_metadata[,"sample_bio_replicate"] <- sapply(cell_metadata[,"sample"], FUN = function(X) sample.bio_replicate.all[which(names(sample.bio_replicate.all) %in% X)])
            
            cell_metadata[,"sample_id"] <- factor(cell_metadata[,"sample_id"], levels = unique(cell_metadata[,"sample_id"]))
            cell_metadata[,"sample_annotation"] <- factor(cell_metadata[,"sample_annotation"], levels = unique(cell_metadata[,"sample_annotation"]))
            cell_metadata[,"sample_time"] <- factor(cell_metadata[,"sample_time"], levels = unique(cell_metadata[,"sample_time"]))
            cell_metadata[,"sample_treatment"] <- factor(cell_metadata[,"sample_treatment"], levels = unique(cell_metadata[,"sample_treatment"]))
            cell_metadata[,"sample_time_treatment"] <- factor(cell_metadata[,"sample_time_treatment"], levels = unique(cell_metadata[,"sample_time_treatment"]))
            cell_metadata[,"sample_bio_replicate"] <- factor(cell_metadata[,"sample_bio_replicate"], levels = unique(cell_metadata[,"sample_bio_replicate"]))
            
            if (Proj_ID[i] == "NG-27918") {
                old.LV = levels(cell_metadata[,"sample_time_treatment"])
                old.LV.name1 <- unique(sapply(old.LV, FUN = function(X) unlist(strsplit(X, '_'))[[1]]))
                old.LV.name2 <- unique(sapply(old.LV, FUN = function(X) unlist(strsplit(X, '_'))[[2]]))
                
                #new.LV = apply(expand.grid(old.LV.name1, old.LV.name2), 1, paste, collapse = '_')
                new.LV = as.vector(outer(old.LV.name1, old.LV.name2, paste, sep = '_'))
                
                cell_metadata[,"sample_time_treatment"] <- factor(cell_metadata[,"sample_time_treatment"], levels = new.LV)
            }
            
            cell_metadata[,"ID"] <- cell_metadata[,"sample.cell_barcode"]
            rownames(cell_metadata) <- cell_metadata[,"ID"]
            
            assign(paste('cell_metadata','.',gsub("\\.| |-","_",sample.name), sep = ''), cell_metadata)
            
            
            # --- expression_matrix (count matrix) --------------------------------------------------------------------------------------------------------------
            #DATA.mtx <- as.matrix(DATA)
            #DATA.df <- as.data.frame(DATA, stringsAsFactors = FALSE)
            
            expression_matrix <- DATA
            expression_matrix <- expression_matrix[, unlist(exp.idx), drop = FALSE]
            
            assign(paste('expression_matrix','.',gsub("\\.| |-","_",sample.name), sep = ''), expression_matrix)
            
            
            # --- gene_metadata ---------------------------------------------------------------------------------------------------------------------------------
            gene_metadata <- id.table.new
            gene_metadata <- gene_metadata[,c("gene_name","gene_id")]
            
            colnames(gene_metadata)[which(colnames(gene_metadata) == "gene_name")] = "gene_short_name"
            rownames(gene_metadata) <- gene_metadata[,"gene_short_name"]
            
            row.idx <- sapply(rownames(expression_matrix), FUN = function(X) which(rownames(gene_metadata) %in% X))
            row.idx <- unlist(row.idx)
            gene_metadata <- gene_metadata[row.idx, , drop = FALSE]
            
            assign(paste('gene_metadata','.',gsub("\\.| |-","_",sample.name), sep = ''), gene_metadata)
            
            sample.name <- sample.name.origin
        }
    }
    cat('\n')
    
    #df.all[is.na(df.all)] <- 0
    
    row.idx <- sapply(all_gene_name, FUN = function(X) which(rownames(df.all) %in% X))
    row.idx <- unlist(row.idx)
    df.all <- df.all[row.idx, , drop = FALSE]
    
    assign(paste('sample','.',gsub("\\.| |-","_",Proj_ID[i]),'.','df_all', sep = ''), df.all)
    
    assign(paste('sample','.',gsub("\\.| |-","_",Proj_ID[i]),'.','cell_num_table', sep = ''), cell.num.table)
    assign(paste('sample','.',gsub("\\.| |-","_",Proj_ID[i]),'.','cell_num_new_table', sep = ''), cell.num.new.table)
    
    cat('[',Proj_ID[i],'] Cell number: ',sum(cell.num.table),'\n', sep = '')
    print(cell.num.table)
    cat('\n')
    #cat('\n')
    
    assign(paste('sample','.',gsub("\\.| |-","_",Proj_ID[i]),'.','gene_num_table', sep = ''), gene.num.table)
    assign(paste('sample','.',gsub("\\.| |-","_",Proj_ID[i]),'.','gene_num_new_table', sep = ''), gene.num.new.table)
    
    cat('[',Proj_ID[i],'] Gene number: ',sum(gene.num.table),'\n', sep = '')
    print(gene.num.table)
    cat('\n')
    cat('\n')
    
    
    RData.name = unlist(strsplit(current_folder, '/'))
    RData.name = RData.name[length(RData.name)]
    save.image(file = paste('./../',RData.name,'_','LeoPharma','_','part1','.RData', sep = ''))
    
    
    # --- cell_metadata ---------------------------------------------------------------------------------------------------------------------------------
    exp.name <- sample.name.new
    
    exp.idx <- lapply(c(1:length(sample.name.new)), FUN = function(X) which(colnames(df.all) %in% colnames(eval(parse(text=paste('sample','.',gsub("\\.| |-","_",sample.name.new[X]),'.','df_new', sep = ''))))))
    names(exp.idx) <- exp.name
    
    exp.num <- sapply(exp.idx, FUN = function(X) length(X))
    
    cell_metadata.all <- data.frame("sample.cell_barcode" = colnames(df.all)[unlist(exp.idx)], 
                                    stringsAsFactors = FALSE)
    
    cell_metadata.all[,"sample"] <- sapply(cell_metadata.all[,"sample.cell_barcode"], FUN = function(X) paste(unlist(strsplit(X, '_'))[[1]],'_',unlist(strsplit(X, '_'))[[2]], sep = ''))
    cell_metadata.all[,"sample"] <- factor(cell_metadata.all[,"sample"], levels = unique(cell_metadata.all[,"sample"]))
    
    cell_metadata.all[,"cell_barcode"] <- sapply(cell_metadata.all[,"sample.cell_barcode"], FUN = function(X) unlist(strsplit(X, '_'))[[3]])
    
    cell_metadata.all[,"sample_id"] <- sapply(cell_metadata.all[,"sample"], FUN = function(X) sample.id.new[which(names(sample.id.new) %in% X)])
    cell_metadata.all[,"sample_annotation"] <- sapply(cell_metadata.all[,"sample"], FUN = function(X) sample.annotation.new[which(names(sample.annotation.new) %in% X)])
    cell_metadata.all[,"sample_time"] <- sapply(cell_metadata.all[,"sample"], FUN = function(X) sample.time.new[which(names(sample.time.new) %in% X)])
    cell_metadata.all[,"sample_treatment"] <- sapply(cell_metadata.all[,"sample"], FUN = function(X) sample.treatment.new[which(names(sample.treatment.new) %in% X)])
    cell_metadata.all[,"sample_time_treatment"] <- sapply(cell_metadata.all[,"sample"], FUN = function(X) paste(sample.time.new[which(names(sample.time.new) %in% X)],'_',sample.treatment.new[which(names(sample.treatment.new) %in% X)], sep = ''))
    cell_metadata.all[,"sample_bio_replicate"] <- sapply(cell_metadata.all[,"sample"], FUN = function(X) sample.bio_replicate.new[which(names(sample.bio_replicate.new) %in% X)])
    
    cell_metadata.all[,"sample_id"] <- factor(cell_metadata.all[,"sample_id"], levels = unique(cell_metadata.all[,"sample_id"]))
    cell_metadata.all[,"sample_annotation"] <- factor(cell_metadata.all[,"sample_annotation"], levels = unique(cell_metadata.all[,"sample_annotation"]))
    cell_metadata.all[,"sample_time"] <- factor(cell_metadata.all[,"sample_time"], levels = unique(cell_metadata.all[,"sample_time"]))
    cell_metadata.all[,"sample_treatment"] <- factor(cell_metadata.all[,"sample_treatment"], levels = unique(cell_metadata.all[,"sample_treatment"]))
    cell_metadata.all[,"sample_time_treatment"] <- factor(cell_metadata.all[,"sample_time_treatment"], levels = unique(cell_metadata.all[,"sample_time_treatment"]))
    cell_metadata.all[,"sample_bio_replicate"] <- factor(cell_metadata.all[,"sample_bio_replicate"], levels = unique(cell_metadata.all[,"sample_bio_replicate"]))
    
    if (Proj_ID[i] == "NG-27918") {
        old.LV = levels(cell_metadata.all[,"sample_time_treatment"])
        old.LV.name1 <- unique(sapply(old.LV, FUN = function(X) unlist(strsplit(X, '_'))[[1]]))
        old.LV.name2 <- unique(sapply(old.LV, FUN = function(X) unlist(strsplit(X, '_'))[[2]]))
        
        #new.LV = apply(expand.grid(old.LV.name1, old.LV.name2), 1, paste, collapse = '_')
        new.LV = as.vector(outer(old.LV.name1, old.LV.name2, paste, sep = '_'))
        
        cell_metadata.all[,"sample_time_treatment"] <- factor(cell_metadata.all[,"sample_time_treatment"], levels = new.LV)
    }
    
    cell_metadata.all[,"ID"] <- cell_metadata.all[,"sample.cell_barcode"]
    rownames(cell_metadata.all) <- cell_metadata.all[,"ID"]
    
    assign(paste('sample','.',gsub("\\.| |-","_",Proj_ID[i]),'.','cell_metadata','.','all', sep = ''), cell_metadata.all)
    
    
    # --- expression_matrix -----------------------------------------------------------------------------------------------------------------------------
    expression_matrix.all = as.matrix(df.all)
    
    assign(paste('sample','.',gsub("\\.| |-","_",Proj_ID[i]),'.','expression_matrix','.','all', sep = ''), expression_matrix.all)
    
    
    # --- gene_metadata ---------------------------------------------------------------------------------------------------------------------------------
    row.idx <- sapply(rownames(df.all), FUN = function(X) which(all_gene_name %in% X))
    row.idx <- unlist(row.idx)
    all_gene_id.all <- all_gene_id[row.idx]
    
    gene_metadata.all <- data.frame("gene_short_name" = rownames(df.all), 
                                    "gene_id" = all_gene_id.all, 
                                    stringsAsFactors = FALSE)
    
    rownames(gene_metadata.all) <- gene_metadata.all[,"gene_short_name"]
    
    assign(paste('sample','.',gsub("\\.| |-","_",Proj_ID[i]),'.','gene_metadata','.','all', sep = ''), gene_metadata.all)
    
    
    # --- Monocle3 & Seurat -----------------------------------------------------------------------------------------------------------------------------
    ## /// all ///
    #project.name = 'LeoPharma'
    #project.subname = Proj_ID[i]
    #project.input = paste(project.name,'.',project.subname, sep = '')
    #
    ## /// Monocle3 ///
    ## *** no genes of interest and signature genes ***
    #Monocle3_analysis(expression_matrix = expression_matrix.all, 
    #                  cell_metadata = cell_metadata.all, 
    #                  gene_metadata = gene_metadata.all, 
    #                  #batch = "sample_bio_replicate", 
    #                  cell_type = c("sample", "sample_time", "sample_treatment", "sample_time_treatment", "sample_bio_replicate"), 
    #                  time = "sample_time", 
    #                  project = project.input, 
    #                  analysis_type = "trajectories", 
    #                  redDim_method = "UMAP", 
    #                  redDim_minDim = 10, 
    #                  #redDim_maxDim = 10, 
    #                  top_gene_num = 10)
    
    
    # /// sample set(s) of interest ///
    sampleSets_of_interest <- list()
    
    if (Proj_ID[i] == "NG-27918") {
        # *** sampleSet1_of_interest ***
        sub_sample.name.new <- sample.name.new[grep("-7550", sample.name.new, value = FALSE, invert = TRUE)]
        sampleSet1_of_interest <- sub_sample.name.new
        
        sampleSet1_of_interest <- list(unique(sampleSet1_of_interest))
        names(sampleSet1_of_interest) = "sampleSet1_of_interest"
        sampleSets_of_interest <- c(sampleSets_of_interest, sampleSet1_of_interest)
        
        # *** sampleSet2_of_interest ***
        sub_sample.name.new <- sample.name.new[grep("-7541", sample.name.new, value = FALSE, invert = TRUE)]
        sampleSet2_of_interest <- sub_sample.name.new
        
        sampleSet2_of_interest <- list(unique(sampleSet2_of_interest))
        names(sampleSet2_of_interest) = "sampleSet2_of_interest"
        sampleSets_of_interest <- c(sampleSets_of_interest, sampleSet2_of_interest)
        
        # *** sampleSet3_of_interest ***
        sub_sample.idx <- which(sample.treatment.new == "unstimulated")
        sub_sample.idx <- sub_sample.idx[grep("-7550", names(sub_sample.idx), value = FALSE, invert = TRUE)]
        sub_sample.idx <- c(sub_sample.idx, which(sample.time.new == "T0")) # T0 with both unstimulated and stimulated
        rm.idx <- which(duplicated(sub_sample.idx))
        if (length(rm.idx) > 0) sub_sample.idx <- sub_sample.idx[-rm.idx]
        sub_sample.idx <- sort(sub_sample.idx)
        
        sub_sample.name.new <- sample.name.new[sub_sample.idx]
        sampleSet3_of_interest <- sub_sample.name.new
        
        sampleSet3_of_interest <- list(unique(sampleSet3_of_interest))
        names(sampleSet3_of_interest) = "unstimulated_with_bothT0"
        sampleSets_of_interest <- c(sampleSets_of_interest, sampleSet3_of_interest)
        
        # *** sampleSet4_of_interest ***
        sub_sample.idx <- which(sample.treatment.new == "unstimulated")
        sub_sample.idx <- sub_sample.idx[grep("-7550", names(sub_sample.idx), value = FALSE, invert = TRUE)]
        #sub_sample.idx <- c(sub_sample.idx, which(sample.time.new == "T0")) # T0 with both unstimulated and stimulated
        rm.idx <- which(duplicated(sub_sample.idx))
        if (length(rm.idx) > 0) sub_sample.idx <- sub_sample.idx[-rm.idx]
        sub_sample.idx <- sort(sub_sample.idx)
        
        sub_sample.name.new <- sample.name.new[sub_sample.idx]
        sampleSet4_of_interest <- sub_sample.name.new
        
        sampleSet4_of_interest <- list(unique(sampleSet4_of_interest))
        names(sampleSet4_of_interest) = "unstimulated"
        sampleSets_of_interest <- c(sampleSets_of_interest, sampleSet4_of_interest)
        
        # *** sampleSet5_of_interest ***
        sub_sample.idx <- which(sample.treatment.new == "stimulated")
        sub_sample.idx <- sub_sample.idx[grep("-7550", names(sub_sample.idx), value = FALSE, invert = TRUE)]
        sub_sample.idx <- c(sub_sample.idx, which(sample.time.new == "T0")) # T0 with both unstimulated and stimulated
        rm.idx <- which(duplicated(sub_sample.idx))
        if (length(rm.idx) > 0) sub_sample.idx <- sub_sample.idx[-rm.idx]
        sub_sample.idx <- sort(sub_sample.idx)
        
        sub_sample.name.new <- sample.name.new[sub_sample.idx]
        sampleSet5_of_interest <- sub_sample.name.new
        
        sampleSet5_of_interest <- list(unique(sampleSet5_of_interest))
        names(sampleSet5_of_interest) = "stimulated_with_bothT0"
        sampleSets_of_interest <- c(sampleSets_of_interest, sampleSet5_of_interest)
        
        # *** sampleSet6_of_interest ***
        sub_sample.idx <- which(sample.treatment.new == "stimulated")
        sub_sample.idx <- sub_sample.idx[grep("-7550", names(sub_sample.idx), value = FALSE, invert = TRUE)]
        #sub_sample.idx <- c(sub_sample.idx, which(sample.time.new == "T0")) # T0 with both unstimulated and stimulated
        rm.idx <- which(duplicated(sub_sample.idx))
        if (length(rm.idx) > 0) sub_sample.idx <- sub_sample.idx[-rm.idx]
        sub_sample.idx <- sort(sub_sample.idx)
        
        sub_sample.name.new <- sample.name.new[sub_sample.idx]
        sampleSet6_of_interest <- sub_sample.name.new
        
        sampleSet6_of_interest <- list(unique(sampleSet6_of_interest))
        names(sampleSet6_of_interest) = "stimulated"
        sampleSets_of_interest <- c(sampleSets_of_interest, sampleSet6_of_interest)
    }
    
    #for (sX in 1:length(sampleSets_of_interest)) {
    #for (sX in 0:length(sampleSets_of_interest)) {
    for (sX in c(1,3:6)) {
        
        if (sX == 0) {
            cell_metadata.samples <- cell_metadata.all
            expression_matrix.samples <- expression_matrix.all
            gene_metadata.samples <- gene_metadata.all
        } else {
            #print(names(sampleSets_of_interest[sX]))
            
            suffix.name = paste(gsub("\\.| ","_",names(sampleSets_of_interest[sX])), sep = '')
            suffix.name = paste(gsub("_of_interest","",suffix.name), sep = '')
            
            samples <- sampleSets_of_interest[[sX]]
            samples.idx <- which(sapply(samples, FUN = function(X) X %in% levels(cell_metadata.all[,"sample"])))
            if (length(samples.idx) > 0) {
                samples = samples[samples.idx]
                
                
                # --- cell_metadata ---------------------------------------------------------------------------------------------------------------------------------
                row.idx <- sapply(samples, FUN = function(X) which(cell_metadata.all[,"sample"] %in% X))
                row.idx <- unlist(row.idx)
                cell_metadata.samples <- cell_metadata.all[row.idx, , drop = FALSE]
                
                cell_metadata.samples[,"sample"] <- sapply(cell_metadata.samples[,"sample.cell_barcode"], FUN = function(X) paste(unlist(strsplit(X, '_'))[[1]],'_',unlist(strsplit(X, '_'))[[2]], sep = ''))
                cell_metadata.samples[,"sample"] <- factor(cell_metadata.samples[,"sample"], levels = unique(cell_metadata.samples[,"sample"]))
                
                cell_metadata.samples[,"cell_barcode"] <- sapply(cell_metadata.samples[,"sample.cell_barcode"], FUN = function(X) unlist(strsplit(X, '_'))[[3]])
                
                cell_metadata.samples[,"sample_id"] <- sapply(cell_metadata.samples[,"sample"], FUN = function(X) sample.id.new[which(names(sample.id.new) %in% X)])
                cell_metadata.samples[,"sample_annotation"] <- sapply(cell_metadata.samples[,"sample"], FUN = function(X) sample.annotation.new[which(names(sample.annotation.new) %in% X)])
                cell_metadata.samples[,"sample_time"] <- sapply(cell_metadata.samples[,"sample"], FUN = function(X) sample.time.new[which(names(sample.time.new) %in% X)])
                cell_metadata.samples[,"sample_treatment"] <- sapply(cell_metadata.samples[,"sample"], FUN = function(X) sample.treatment.new[which(names(sample.treatment.new) %in% X)])
                cell_metadata.samples[,"sample_time_treatment"] <- sapply(cell_metadata.samples[,"sample"], FUN = function(X) paste(sample.time.new[which(names(sample.time.new) %in% X)],'_',sample.treatment.new[which(names(sample.treatment.new) %in% X)], sep = ''))
                cell_metadata.samples[,"sample_bio_replicate"] <- sapply(cell_metadata.samples[,"sample"], FUN = function(X) sample.bio_replicate.new[which(names(sample.bio_replicate.new) %in% X)])
                
                cell_metadata.samples[,"sample_id"] <- factor(cell_metadata.samples[,"sample_id"], levels = unique(cell_metadata.samples[,"sample_id"]))
                cell_metadata.samples[,"sample_annotation"] <- factor(cell_metadata.samples[,"sample_annotation"], levels = unique(cell_metadata.samples[,"sample_annotation"]))
                cell_metadata.samples[,"sample_time"] <- factor(cell_metadata.samples[,"sample_time"], levels = unique(cell_metadata.samples[,"sample_time"]))
                cell_metadata.samples[,"sample_treatment"] <- factor(cell_metadata.samples[,"sample_treatment"], levels = unique(cell_metadata.samples[,"sample_treatment"]))
                cell_metadata.samples[,"sample_time_treatment"] <- factor(cell_metadata.samples[,"sample_time_treatment"], levels = unique(cell_metadata.samples[,"sample_time_treatment"]))
                cell_metadata.samples[,"sample_bio_replicate"] <- factor(cell_metadata.samples[,"sample_bio_replicate"], levels = unique(cell_metadata.samples[,"sample_bio_replicate"]))
                
                if (Proj_ID[i] == "NG-27918") {
                    if (suffix.name == "sampleSet1") {
                        levels(cell_metadata.samples[,"sample"]) = gsub("-7541","",levels(cell_metadata.samples[,"sample"]))
                    } else if (suffix.name == "sampleSet2") {
                        levels(cell_metadata.samples[,"sample"]) = gsub("-7550","",levels(cell_metadata.samples[,"sample"]))
                    }
                }
                
                if (Proj_ID[i] == "NG-27918") {
                    old.LV = levels(cell_metadata.samples[,"sample_time_treatment"])
                    old.LV.name1 <- unique(sapply(old.LV, FUN = function(X) unlist(strsplit(X, '_'))[[1]]))
                    old.LV.name2 <- unique(sapply(old.LV, FUN = function(X) unlist(strsplit(X, '_'))[[2]]))
                    
                    #new.LV = apply(expand.grid(old.LV.name1, old.LV.name2), 1, paste, collapse = '_')
                    new.LV = as.vector(outer(old.LV.name1, old.LV.name2, paste, sep = '_'))
                    
                    cell_metadata.samples[,"sample_time_treatment"] <- factor(cell_metadata.samples[,"sample_time_treatment"], levels = new.LV)
                }
                
                cell_metadata.samples[,"ID"] <- cell_metadata.samples[,"sample.cell_barcode"]
                rownames(cell_metadata.samples) <- cell_metadata.samples[,"ID"]
                
                assign(paste('sample','.',gsub("\\.| |-","_",Proj_ID[i]),'.','cell_metadata','.',gsub("\\.| |-","_",suffix.name), sep = ''), cell_metadata.samples)
                
                
                # --- expression_matrix -----------------------------------------------------------------------------------------------------------------------------
                col.idx <- sapply(cell_metadata.samples[,"ID"], FUN = function(X) which(colnames(expression_matrix.all) %in% X))
                col.idx <- unlist(col.idx)
                expression_matrix.samples <- expression_matrix.all[, col.idx, drop = FALSE]
                
                assign(paste('sample','.',gsub("\\.| |-","_",Proj_ID[i]),'.','expression_matrix','.',gsub("\\.| |-","_",suffix.name), sep = ''), expression_matrix.samples)
                
                
                # --- gene_metadata ---------------------------------------------------------------------------------------------------------------------------------
                row.idx <- sapply(rownames(expression_matrix.samples), FUN = function(X) which(rownames(gene_metadata.all) %in% X))
                row.idx <- unlist(row.idx)
                gene_metadata.samples <- gene_metadata.all[row.idx, , drop = FALSE]
                
                assign(paste('sample','.',gsub("\\.| |-","_",Proj_ID[i]),'.','gene_metadata','.',gsub("\\.| |-","_",suffix.name), sep = ''), gene_metadata.samples)
            } else {
                next
            }
        }
        
        if (sX == 0) {
            expression_matrix.new <- expression_matrix.all
            cell_metadata.new <- cell_metadata.all
            gene_metadata.new <- gene_metadata.all
        } else {
            expression_matrix.new <- expression_matrix.samples
            cell_metadata.new <- cell_metadata.samples
            gene_metadata.new <- gene_metadata.samples
        }
        
        # /// all ///
        project.name = 'LeoPharma'
        if (sX == 0) {
            project.subname = Proj_ID[i]
        } else {
            project.subname = paste(Proj_ID[i],'.',suffix.name, sep = '')
        }
        project.input = paste(project.name,'.',project.subname, sep = '')
        
        # /// Monocle3 ///
        # *** no genes of interest and signature genes ***
        Monocle3_analysis(expression_matrix = expression_matrix.new, 
                          cell_metadata = cell_metadata.new, 
                          gene_metadata = gene_metadata.new, 
                          #batch = "sample_bio_replicate", 
                          cell_type = c("sample", "sample_time", "sample_treatment", "sample_time_treatment", "sample_bio_replicate"), 
                          time = "sample_time", 
                          project = project.input, 
                          analysis_type = "trajectories", 
                          redDim_method = "UMAP", 
                          redDim_minDim = 10, 
                          #redDim_maxDim = 10, 
                          top_gene_num = 10)
    }
    
    
    #RData.name = unlist(strsplit(current_folder, '/'))
    #RData.name = RData.name[length(RData.name)]
    #save.image(file = paste('./../',RData.name,'_','LeoPharma','_','part2','.RData', sep = ''))
    
    
    mainDir <- paste(current_folder,'/Result','/',Proj_ID[i], sep = '')
    if (file.exists(file.path(mainDir))) {
        if (length(list.files(path = mainDir, full.names = TRUE, recursive = FALSE)) == 0) {
            unlink(mainDir, recursive = TRUE)
        }
    }
}

