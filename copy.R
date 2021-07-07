#!/usr/bin/env Rscript

rm(list = ls()) # remove all variables in workspace

# ---------------------------------------------------------------------------------------------------------------------------------------------------
#current_folder = "~/Project/Pancreas"
current_folder = "/Users/mingtsung/Desktop/Project/Pancreas"
setwd(current_folder)

#install.packages("mgsub") # mgsub: Safe, Multiple, Simultaneous String Substitution
#install.packages("RColorBrewer") # RColorBrewer: ColorBrewer Palettes
#install.packages("corrplot") # corrplot: Visualization of a Correlation Matrix

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("ComplexHeatmap")

#suppressMessages(library(mgsub)) # mgsub()
suppressMessages(library(RColorBrewer))
suppressMessages(library(corrplot))
suppressMessages(library(ComplexHeatmap)) # Heatmap()

mainDir <- current_folder
subDir <- "Result"
if (! file.exists(file.path(mainDir, subDir))) {
    dir.create(file.path(mainDir, subDir))
}


# --- Source function -------------------------------------------------------------------------------------------------------------------------------
source('./Script/function/QC/qc_data.R')

gene_cutoff <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
gene_cutoff.default = gene_cutoff[2] # min_gene_frac=0.1


# --- Read in DATA ----------------------------------------------------------------------------------------------------------------------------------
sample <- c("HPNE", # human, pancreas, duct
            "MIA") # human, pancreas, carcinoma

# /// control ///
control <- c("16726_hpne_6h_dmso_counts_new") # human, pancreas, duct

# /// sc26 time-course DATA count matrices ///
disease <- c("mia_6h_dmso_counts_new") # human, pancreas, carcinoma

# https://www.r-bloggers.com/how-to-expand-color-palette-with-ggplot-and-rcolorbrewer/
colourCount = length(disease)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
disease.color <- getPalette(colourCount)

exp.name.all <- c()
cell_metadata.all <- c()

for (s in 1:length(sample)) {
    
    gc() # Garbage Collection
    
    mainDir <- paste(current_folder,'/Result', sep = '')
    subDir <- 'library'
    if (! file.exists(file.path(mainDir, subDir))) {
        dir.create(file.path(mainDir, subDir))
    }
    
    mainDir <- paste(current_folder,'/Result','/','library', sep = '')
    subDir <- sample[s]
    if (! file.exists(file.path(mainDir, subDir))) {
        dir.create(file.path(mainDir, subDir))
    }
    
    mainDir <- paste(current_folder,'/Result','/','library','/',sample[s], sep = '')
    subDir <- 'raw'
    if (! file.exists(file.path(mainDir, subDir))) {
        dir.create(file.path(mainDir, subDir))
    }
    
    if (sample[s] == "HPNE") {
        file_num = length(control)
        sample_folder = 'MIA_HPME_6h_dmso_counts'
    } else if (sample[s] == "MIA") {
        file_num = length(disease)
        sample_folder = 'MIA_HPME_6h_dmso_counts'
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
        
        if (sample[s] == "HPNE") {
            sample_file = control[i]
        } else {
            sample_file = disease[i]
        }
        
        control.idx <- which(control %in% sample_file)
        disease.idx <- which(disease %in% sample_file)
        if (length(control.idx) == 1) {
            if (control.idx > 1) {
                presample_file <- control[control.idx - 1]
            }
        } else if (length(disease.idx) == 1) {
            if (disease.idx > 1) {
                presample_file <- disease[disease.idx - 1]
            }
        } else {
            presample_file <- NULL
        }
        
        DATA <- read.csv(paste('./Data/',sample_folder,'/',sample_file,'.csv', sep = ''), sep = ',', header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
        #print(sample_file); print(dim(DATA))
        
        genes.idx <- which(colnames(DATA) == "genes")
        if (length(genes.idx) == 0) {
            genes.NA = 1
            DATA[,"genes"] <- rownames(DATA)
        } else {
            genes.NA = 0
            NA.idx <- which(DATA[,"genes"] == "")
            DATA[NA.idx,"genes"] <- rownames(DATA)[NA.idx]
        }
        genes.idx <- which(colnames(DATA) == "genes")
        if (genes.idx != 1) {
            col.idx <- c(genes.idx, 1:(ncol(DATA)-1))
            DATA <- DATA[, col.idx]
        }
        
        gene_names.idx <- which(colnames(DATA) == "gene_names")
        if (length(gene_names.idx) == 0) {
            gene_names.NA = 1
            DATA[,"gene_names"] <- rownames(DATA)
        } else {
            gene_names.NA = 0
            NA.idx <- which(DATA[,"gene_names"] == "")
            DATA[NA.idx,"gene_names"] <- rownames(DATA)[NA.idx]
        }
        gene_names.idx <- which(colnames(DATA) == "gene_names")
        if (gene_names.idx != 2) {
            col.idx <- c(1, gene_names.idx, 2:(ncol(DATA)-1))
            DATA <- DATA[, col.idx]
        }
        
        cell.num = ncol(DATA) - 2
        names(cell.num) <- sample_file
        cell.num.table <- c(cell.num.table, cell.num)
        
        gene.num = nrow(DATA)
        names(gene.num) <- sample_file
        gene.num.table <- c(gene.num.table, gene.num)
        
        if (cell.num > 0) {
            
            cell.barcode <- trimws(colnames(DATA)[3:ncol(DATA)])
            
            if (length(grep('-1$', cell.barcode, value = TRUE)) > 0) { # "kcppc_bam_counts" # old ductal cell (new processing pipeline, without gene_names)
                col.idx <- grep('-1$', cell.barcode, value = FALSE)
                cell.barcode[col.idx] <- gsub('-1$', '', cell.barcode[col.idx])
            }
            
            cell.barcode.all.name <- sapply(cell.barcode.all, FUN = function(X) unlist(strsplit(X, '\\.'))[[1]])
            
            
            # --- QC 1. Check duplicated cell barcode -----------------------------------------------------------------------------------------------------------
            cell.barcode.all <- c(cell.barcode.all, cell.barcode)
            
            cell.barcode.list[[i]] <- cell.barcode
            names(cell.barcode.list)[i] <- sample_file
            
            # /// duplicated cell barcode in the same file ///
            dup.idx <- which(duplicated(cell.barcode))
            if (length(dup.idx) > 0) {
                if (in.idx == 0) {
                    cat('[',sample[s],'] Duplicated cell barcode in "',sample[s],'" dataset: ','\n', sep = '')
                    in.idx = 1
                }
                cat('Cell barcodes are duplicated in "',paste(sample_file,'.csv', sep = ''),'". ', sep = '')
                if (length(dup.idx) == 1) {
                    cat('Duplicated cell barcode: ',paste(cell.barcode[dup.idx], collapse = ', '),'\n', sep = '')
                } else if (length(dup.idx) > 1) {
                    cat('Duplicated cell barcodes: ',paste(cell.barcode[dup.idx], collapse = ', '),'\n', sep = '')
                }
            }
            
            
            # --- Data curation ---------------------------------------------------------------------------------------------------------------------------------
            # /// cell: append sample name as suffix to duplicated cells ///
            cell.barcode_dup.idx <- which(cell.barcode %in% cell.barcode.all.name)
            cell.barcode.all_dup.idx <- which(cell.barcode.all.name %in% cell.barcode)
            
            if (length(cell.barcode_dup.idx) > 0) {
                for (j in 1:length(cell.barcode_dup.idx)) {
                    if (length(grep('\\.', cell.barcode[cell.barcode_dup.idx[j]], value = TRUE)) == 0) {
                        cell.barcode[cell.barcode_dup.idx[j]] <- paste(cell.barcode[cell.barcode_dup.idx[j]], '.', sample_file, sep = '')
                    }
                }
            }
            if (length(cell.barcode.all_dup.idx) > 0) {
                for (j in 1:length(cell.barcode.all_dup.idx)) {
                    if (length(grep('\\.', cell.barcode.all[cell.barcode.all_dup.idx[j]], value = TRUE)) == 0) {
                        #cell.barcode.all[cell.barcode.all_dup.idx[j]] <- paste(cell.barcode.all[cell.barcode.all_dup.idx[j]], '.', presample_file, sep = '')
                        
                        if (length(control.idx) == 1) {
                            if (control.idx > 1) {
                                presample <- control
                                presample.idx <- control.idx - 1
                            }
                        } else if (length(disease.idx) == 1) {
                            if (disease.idx > 1) {
                                presample <- disease
                                presample.idx <- disease.idx - 1
                            }
                        } else {
                            presample <- NULL
                            presample.idx <- NULL
                        }
                        
                        col.num = 0
                        for (k in 1:presample.idx) {
                            presample_file <- presample[k]
                            
                            presample_file.df <- eval(parse(text=paste(sample[s],'.',presample_file,'.','df', sep = '')))
                            col.idx <- which(colnames(presample_file.df) == cell.barcode.all[cell.barcode.all_dup.idx[j]])
                            colnames(presample_file.df)[col.idx] <- paste(colnames(presample_file.df)[col.idx], '.', presample_file, sep = '')
                            assign(paste(sample[s],'.',presample_file,'.','df', sep = ''), presample_file.df)
                            
                            presample_file.df_new <- eval(parse(text=paste(sample[s],'.',presample_file,'.','df_new', sep = '')))
                            col.idx <- which(colnames(presample_file.df_new) == cell.barcode.all[cell.barcode.all_dup.idx[j]])
                            colnames(presample_file.df_new)[col.idx] <- paste(colnames(presample_file.df_new)[col.idx], '.', presample_file, sep = '')
                            assign(paste(sample[s],'.',presample_file,'.','df_new', sep = ''), presample_file.df_new)
                            col.num = col.num + ncol(presample_file.df_new)
                            
                            #presample_file.df_new_normalized <- eval(parse(text=paste(sample[s],'.',presample_file,'.','df_new_normalized', sep = '')))
                            #col.idx <- which(colnames(presample_file.df_new_normalized) == cell.barcode.all[cell.barcode.all_dup.idx[j]])
                            #colnames(presample_file.df_new_normalized)[col.idx] <- paste(colnames(presample_file.df_new_normalized)[col.idx], '.', presample_file, sep = '')
                            #assign(paste(sample[s],'.',presample_file,'.','df_new_normalized', sep = ''), presample_file.df_new_normalized)
                            
                            presample_file.df_all <- eval(parse(text=paste(sample[s],'.','df_all', sep = '')))
                            col.idx <- which(colnames(presample_file.df_all) == cell.barcode.all[cell.barcode.all_dup.idx[j]])
                            if (length(col.idx) > 0) {
                                if (col.idx <= col.num) {
                                    colnames(presample_file.df_all)[col.idx] <- paste(colnames(presample_file.df_all)[col.idx], '.', presample_file, sep = '')
                                }
                            }
                            assign(paste(sample[s],'.','df_all', sep = ''), presample_file.df_all)
                            assign("df.all", presample_file.df_all)
                        }
                    }
                }
            }
            
            # /// gene: replace empty gene_names with ID and append ID to duplicated gene_names ///
            df <- DATA[, 3:ncol(DATA), drop = FALSE]
            
            colnames(df) <- cell.barcode
            
            gene.dup.idx <- which(duplicated(DATA[,"gene_names"]))
            if (length(gene.dup.idx) > 0) {
                for (j in 1:length(gene.dup.idx)) {
                    gene_name <- DATA[which(DATA[,"gene_names"] == DATA[gene.dup.idx[j],"gene_names"]),"gene_names"]
                    gene_name <- unique(gene_name)
                    
                    if (length(grep(')$', DATA[which(DATA[,"gene_names"] == gene_name),"gene_names"], value = TRUE)) == 0) {
                        DATA[which(DATA[,"gene_names"] == gene_name),"gene_names"] <- paste(DATA[which(DATA[,"gene_names"] == gene_name),"gene_names"], '(', DATA[which(DATA[,"gene_names"] == gene_name),"genes"], ')', sep = '')
                    }
                }
            }
            
            rownames(df) <- DATA[,"gene_names"]
            all_gene_name <- rownames(df)
            
            if (sample[s] == "HPNE") {
                control_gene_name <- all_gene_name
            }
            
            write.csv(df, file = paste(current_folder,'/Result','/','library','/',sample[s],'/','raw','/',sample[s],'.',sample_file,'.csv', sep = ''), row.names = TRUE, quote = TRUE)
            
            
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
                
                assign(paste(sample[s],'.',sample_file,'.',gene_cutoff[g],'.','cell_num_new', sep = ''), cell.num.new)
                assign(paste(sample[s],'.',sample_file,'.',gene_cutoff[g],'.','gene_num_new', sep = ''), gene.num.new)
            }
            
            df_qc <- qc_data(df, normalize = FALSE, min_gene_frac = gene_cutoff.default)
            df_qc.normalized <- qc_data(df, normalize = TRUE, min_gene_frac = gene_cutoff.default)
            
            df.new <- df_qc
            df.new.normalized <- df_qc.normalized
            
            assign(paste(sample[s],'.',sample_file,'.','df', sep = ''), df)
            assign(paste(sample[s],'.',sample_file,'.','df_new', sep = ''), df.new)
            #assign(paste(sample[s],'.',sample_file,'.','df_new_normalized', sep = ''), df.new.normalized)
            
            if (! is.null(df.new)) {
                cell.num.new = ncol(df.new)
                gene.num.new = nrow(df.new)
                
                df.all <- merge(df.all, df.new, by = "row.names", all = TRUE, sort = TRUE)
                rownames(df.all) <- df.all$Row.names
                df.all$Row.names <- NULL
                
                assign(paste(sample[s],'.','df_all', sep = ''), df.all)
            } else {
                cell.num.new = 0
                gene.num.new = 0
            }
        } else {
            
            cell.barcode.list[[i]] <- NA
            names(cell.barcode.list)[i] <- sample_file
            
            assign(paste(sample[s],'.',sample_file,'.','df', sep = ''), NULL)
            assign(paste(sample[s],'.',sample_file,'.','df_new', sep = ''), NULL)
            #assign(paste(sample[s],'.',sample_file,'.','df_new_normalized', sep = ''), NULL)
            
            for (g in 1:length(gene_cutoff)) {
                cell.num.new = 0
                gene.num.new = 0
                
                assign(paste(sample[s],'.',sample_file,'.',gene_cutoff[g],'.','cell_num_new', sep = ''), cell.num.new)
                assign(paste(sample[s],'.',sample_file,'.',gene_cutoff[g],'.','gene_num_new', sep = ''), gene.num.new)
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
    
    assign(paste(sample[s],'.','df_all', sep = ''), df.all)
    
    
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
    
    assign(paste(sample[s],'.','cell_num_table', sep = ''), cell.num.table)
    assign(paste(sample[s],'.','cell_num_new_table', sep = ''), cell.num.new.table)
    
    cat('[',sample[s],'] Cell number before QC: ',sum(cell.num.table),'\n', sep = '')
    print(cell.num.table)
    cat('\n')
    cat('[',sample[s],'] Cell number after QC: ',sum(cell.num.new.table),'\n', sep = '')
    print(cell.num.new.table)
    cat('\n')
    #cat('\n')
    
    assign(paste(sample[s],'.','gene_num_table', sep = ''), gene.num.table)
    assign(paste(sample[s],'.','gene_num_new_table', sep = ''), gene.num.new.table)
    
    cat('[',sample[s],'] Gene number before QC: ',sum(gene.num.table),'\n', sep = '')
    print(gene.num.table)
    cat('\n')
    cat('[',sample[s],'] Gene number after QC: ',sum(gene.num.new.table),'\n', sep = '')
    print(gene.num.new.table)
    cat('\n')
    cat('\n')
    
    
    # --- cell_metadata ---------------------------------------------------------------------------------------------------------------------------------
    if (sample[s] == "HPNE") {
        exp.name <- control
    } else {
        exp.name <- disease
    }
    exp.name.all <- c(exp.name.all, exp.name)
    exp.name.all <- unique(exp.name.all)
    
    exp.idx <- lapply(exp.name, FUN = function(X) which(colnames(df.all) %in% colnames(eval(parse(text=paste(sample[s],'.',X,'.','df_new', sep = ''))))))
    names(exp.idx) <- exp.name
    
    exp.num <- sapply(exp.idx, FUN = function(X) length(X))
    
    cell_metadata <- data.frame("cell_barcode" = colnames(df.all)[unlist(exp.idx)], 
                                "sample" = rep(sample[s], ncol(df.all)), 
                                "DATA" = rep(names(exp.num), exp.num), 
                                stringsAsFactors = FALSE)
    rownames(cell_metadata) <- cell_metadata[,"cell_barcode"]
    
    cell_metadata[,"ID"] <- cell_metadata[,"cell_barcode"]
    cell_metadata[,"Condition"] <- paste(cell_metadata[,"sample"], cell_metadata[,"DATA"], sep = '.')
    
    cell_metadata[,"cell_barcode"] <- sapply(cell_metadata[,"cell_barcode"], FUN = function(X) unlist(strsplit(X, '\\.'))[[1]])
    
    #assign(paste(sample[s],'.','cell_metadata', sep = ''), cell_metadata)
    
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
    mainDir <- paste(current_folder,'/Result','/','library','/',sample[s], sep = '')
    subDir <- 'QC'
    if (! file.exists(file.path(mainDir, subDir))) {
        dir.create(file.path(mainDir, subDir))
    }
    
    initial.idx = 1
    initial2end.idx <- c()
    df.subset.all <- c()
    for (i in 1:length(exp.name)) {
        df <- eval(parse(text=paste(sample[s],'.',exp.name[i],'.','df', sep = '')))
        
        df.all.subset <- df.all[, exp.idx[[i]], drop = FALSE]
        
        if (! is.null(df)) {
            
            
            
            row.idx <- unlist(sapply(rownames(df.all.subset), FUN = function(X) which(rownames(df) %in% X)))
            col.idx <- unlist(sapply(colnames(df.all.subset), FUN = function(X) which(colnames(df) %in% X)))
            
            df.subset <- df[row.idx, col.idx, drop = FALSE]
            
            #assign(paste(sample[s],'.',exp.name[i],'.','df_subset', sep = ''), df.subset)
            #write.csv(df.subset, file = paste(current_folder,'/Result','/','library','/',sample[s],'/','QC','/',sample[s],'.QC','.',exp.name[i],'.csv', sep = ''), row.names = TRUE, quote = TRUE)
            
            df.subset.all <- merge(df.subset.all, df.subset, by = "row.names", all = TRUE, sort = TRUE)
            rownames(df.subset.all) <- df.subset.all$Row.names
            df.subset.all$Row.names <- NULL
        }
    }
    row.idx <- unlist(sapply(all_gene_name, FUN = function(X) which(rownames(df.subset.all) %in% X)))
    df.subset.all <- df.subset.all[row.idx, , drop = FALSE]
    assign(paste(sample[s],'.','df_subset_all', sep = ''), df.subset.all)
    #write.csv(df.subset.all, file = paste(current_folder,'/Result','/','library','/',sample[s],'/',sample[s],'.QC','.csv', sep = ''), row.names = TRUE, quote = TRUE)
    
    for (i in 1:length(exp.name)) {
        df.new <- eval(parse(text=paste(sample[s],'.',exp.name[i],'.','df_new', sep = '')))
        
        if (! is.null(df.new)) {
            row.idx <- which(rownames(df.subset.all) %in% rownames(df.new))
            col.idx <- which(colnames(df.subset.all) %in% colnames(df.new))
            
            df.subset <- df.subset.all[row.idx, col.idx, drop = FALSE]
            
            write.csv(df.subset, file = paste(current_folder,'/Result','/','library','/',sample[s],'/','QC','/',sample[s],'.QC','.',exp.name[i],'.csv', sep = ''), row.names = TRUE, quote = TRUE)
        }
    }
    
    if (sample[s] != "HPNE") {
        
        # /// combine DATA and control ///
        mainDir <- paste(current_folder,'/Result','/','library','/',sample[s], sep = '')
        subDir <- 'control'
        if (! file.exists(file.path(mainDir, subDir))) {
            dir.create(file.path(mainDir, subDir))
        }
        
        df.subset.all.control <- eval(parse(text=paste('HPNE','.','df_subset_all', sep = '')))
        
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
                colnames(df.subset.all.control)[df.subset.all.control.dup.idx] <- paste(colnames(df.subset.all.control)[df.subset.all.control.dup.idx], '.', 'HPNE', sep = '')
            }
        }
        
        df2control <- merge(df.subset.all.new, df.subset.all.control, by = "row.names", all = TRUE, sort = TRUE)
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
                cell.num.new.table.new <- eval(parse(text=paste(sample[s],'.','cell_num_new_table', sep = '')))
            } else if (run == 2) {
                cell.num.new.table.new <- eval(parse(text=paste('HPNE','.','cell_num_new_table', sep = '')))
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
        
        DATA2cell <- sapply(cell.num.new.table.new.idx.all, FUN = function(X) colnames(df2control)[X])
        cell2DATA <- sapply(colnames(df2control), FUN = function(X) names(which(sapply(DATA2cell, FUN = function(Y) X %in% unlist(Y)))))
        col.name <- cell2DATA
        
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
                sample.name <- "HPNE"
                sample_file <- exp.name.all[i]
                row.idx <- row.idx1
            } else {
                sample.name <- sample[s]
                sample_file <- exp.name.all[i]
                row.idx <- row.idx2
            }
            
            df <- eval(parse(text=paste(sample.name,'.',sample_file,'.','df', sep = '')))
            
            if (! is.null(df)) {
                rownames(df)[which(rownames(df) == "Tead2")] = "Tead2(ENSMUSG00000030796)"
                rownames(df)[which(rownames(df) == "l7Rn6")] = "l7Rn6(ENSMUSG00000062797)"
                
                NA.idx <- which(! rownames(df2control) %in% rownames(df))
                df2control[NA.idx, col.idx][is.na(df2control[NA.idx, col.idx])] <- 0
                
                
                
                
                df2control[names(row.idx), col.idx] = df[row.idx, names(col.idx)]
            }
        }
        
        #df2control[is.na(df2control)] <- 0
        
        row.idx <- which(rownames(df2control) %in% intersect(all_gene_name.new, control_gene_name))
        df2control <- df2control[row.idx, , drop = FALSE]
        
        #assign(paste(sample[s],'.','df2control', sep = ''), df2control)
        #write.csv(df2control, file = paste(current_folder,'/Result','/','library','/',sample[s],'/',sample[s],'.control','.csv', sep = ''), row.names = TRUE, quote = TRUE)
        
        # /// split into individual DATA and control ///
        for (i in 1:length(exp.name.all)) {
            col.idx <- which(col.name == exp.name.all[i])
            
            if (length(col.idx) > 0) {
                df2control.subset <- df2control[, col.idx, drop = FALSE]
                #assign(paste(sample[s],'.',exp.name.all[i],'.','df2control_subset', sep = ''), df2control.subset)
                write.csv(df2control.subset, file = paste(current_folder,'/Result','/','library','/',sample[s],'/','control','/',sample[s],'.control','.',exp.name.all[i],'.csv', sep = ''), row.names = TRUE, quote = TRUE)
            }
        }
    }
}


    # --- different gene cutoff -------------------------------------------------------------------------------------------------------------------------
    for (s in 1:length(sample)) {
        
        if (sample[s] == "HPNE") {
            file_num = length(control)
        } else {
            file_num = length(disease)
        }
        
        xy.all <- c()
        
        for (i in 1:file_num) {
            
            if (sample[s] == "HPNE") {
                sample_file = control[i]
            } else {
                sample_file = disease[i]
            }
            
            cell.num.new.all <- c()
            gene.num.new.all <- c()
            
            for (g in 1:length(gene_cutoff)) {
                cell.num.new <- eval(parse(text=paste(sample[s],'.',sample_file,'.',gene_cutoff[g],'.','cell_num_new', sep = '')))
                gene.num.new <- eval(parse(text=paste(sample[s],'.',sample_file,'.',gene_cutoff[g],'.','gene_num_new', sep = '')))
                
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
        
        p <- ggplot(data = xy.all, aes(x = gene_cutoff, y = gene_num, col = group, group = group)) + 
            geom_line(aes(color = group)) + 
            geom_point(aes(color = group))
        p <- p + 
            theme(panel.border = element_blank(), axis.line = element_line()) + 
            theme(panel.background = element_rect(fill = "white", colour = "white")) + 
            theme(panel.grid = element_line(colour = "grey80"), panel.grid.minor = element_blank())
        p <- p + 
            ggtitle(sample[s]) + theme(plot.title = element_text(hjust = 0.5, size = 24, family = "Helvetica")) + 
            xlab("gene cutoff") + theme(axis.title.x = element_text(size = 20, family = "Helvetica"), axis.text.x = element_text(size = 16, family = "Helvetica")) + 
            ylab("gene number") + theme(axis.title.y = element_text(size = 20, family = "Helvetica"), axis.text.y = element_text(size = 16, family = "Helvetica")) + 
            labs(col = "DATA") + theme(legend.title = element_text(size = 20, family = "Helvetica"), legend.text = element_text(size = 12, family = "Courier"))
        p <- p + 
            scale_x_continuous(limits = c(0, max(xbreaks)), expand = c(0, 0), breaks = xbreaks, labels = prettyZero) + 
            scale_y_continuous(limits = c(0, max(ybreaks)), expand = c(0, 0), breaks = ybreaks)
        
        png(file = paste(current_folder,'/Result','/','library','/',sample[s],'/',sample[s],'.','gene_cutoff2gene_number','.png', sep = ''), width = 10, height = 6, units = 'in', res = 600)
        print(p)
        dev.off()
        
    }
    
    # --- dynamics --------------------------------------------------------------------------------------------------------------------------------------
    for (g in 1:length(gene_cutoff)) {
        
        file_num = length(disease)
        
        xy.all <- c()
        xz.all <- c()
        
        for (i in 1:file_num) {
            
            sample_file = disease[i]
            
            cell.num.new.all <- c()
            gene.num.new.all <- c()
            
            for (s in 1:length(sample)) {
                
                if (sample[s] != "HPNE") {
                    
                    cell.num.new <- eval(parse(text=paste(sample[s],'.',sample_file,'.',gene_cutoff[g],'.','cell_num_new', sep = '')))
                    gene.num.new <- eval(parse(text=paste(sample[s],'.',sample_file,'.',gene_cutoff[g],'.','gene_num_new', sep = '')))
                    
                    cell.num.new.all <- c(cell.num.new.all, cell.num.new)
                    gene.num.new.all <- c(gene.num.new.all, gene.num.new)
                    
                }
            }
            
            x = sample[sample != "HPNE"]
            
            
            
            
            
            
            
            y = gene.num.new.all
            z = cell.num.new.all
            
            xy <- data.frame(library = x, gene_num = y, group = sample_file)
            xz <- data.frame(library = x, cell_num = z, group = sample_file)
            
            xy.all <- rbind(xy.all, xy)
            xz.all <- rbind(xz.all, xz)
            
        }
        
        
        
        
        
        library(ggplot2)
        
        # /// gene dynamics ///
        ybreaks = seq(0, max(xy.all[,"gene_num"]), 1000)
        ybreaks = unique(c(ybreaks, ybreaks[length(ybreaks)]+1000))
        
        p <- ggplot(data = xy.all, aes(x = factor(library), y = gene_num, col = group, group = group)) + 
            geom_line(aes(color = group)) + 
            geom_point(aes(color = group))
        p <- p + 
            theme(panel.border = element_blank(), axis.line = element_line()) + 
            theme(panel.background = element_rect(fill = "white", colour = "white")) + 
            theme(panel.grid = element_line(colour = "grey80"), panel.grid.minor = element_blank())
        p <- p + 
            ggtitle(paste('sc26',' (','gene cutoff',' = ',gene_cutoff[g],')', sep = '')) + theme(plot.title = element_text(hjust = 0.5, size = 24, family = "Helvetica")) + 
            xlab("sample library") + theme(axis.title.x = element_text(size = 20, family = "Helvetica"), axis.text.x = element_text(size = 16, family = "Helvetica")) + 
            ylab("gene number") + theme(axis.title.y = element_text(size = 20, family = "Helvetica"), axis.text.y = element_text(size = 16, family = "Helvetica")) + 
            labs(col = "DATA") + theme(legend.title = element_text(size = 20, family = "Helvetica"), legend.text = element_text(size = 12, family = "Courier"))
        p <- p + 
            scale_y_continuous(limits = c(0, max(ybreaks)), expand = c(0, 0), breaks = ybreaks)
        
        png(file = paste(current_folder,'/Result','/','library','/','sample_library2gene_number','.','gene_cutoff','=',gene_cutoff[g],'.png', sep = ''), width = 12, height = 6, units = 'in', res = 600)
        print(p)
        dev.off()
        
        # /// cell dynamics ///
        if (gene_cutoff[g] == 0.1) {
            
            ybreaks = seq(0, max(xz.all[,"cell_num"]), 1000)
            ybreaks = unique(c(ybreaks, ybreaks[length(ybreaks)]+1000))
            
            p <- ggplot(data = xz.all, aes(x = factor(library), y = cell_num, col = group, group = group)) + 
                geom_line(aes(color = group)) + 
                geom_point(aes(color = group))
            p <- p + 
                theme(panel.border = element_blank(), axis.line = element_line()) + 
                theme(panel.background = element_rect(fill = "white", colour = "white")) + 
                theme(panel.grid = element_line(colour = "grey80"), panel.grid.minor = element_blank())
            p <- p + 
                ggtitle('sc26') + theme(plot.title = element_text(hjust = 0.5, size = 24, family = "Helvetica")) + 
                xlab("sample library") + theme(axis.title.x = element_text(size = 20, family = "Helvetica"), axis.text.x = element_text(size = 16, family = "Helvetica")) + 
                ylab("cell number") + theme(axis.title.y = element_text(size = 20, family = "Helvetica"), axis.text.y = element_text(size = 16, family = "Helvetica")) + 
                labs(col = "DATA") + theme(legend.title = element_text(size = 20, family = "Helvetica"), legend.text = element_text(size = 12, family = "Courier"))
            p <- p + 
                scale_y_continuous(limits = c(0, max(ybreaks)), expand = c(0, 0), breaks = ybreaks)
            
            png(file = paste(current_folder,'/Result','/','library','/','sample_library2cell_number','.png', sep = ''), width = 12, height = 6, units = 'in', res = 600)
            print(p)
            dev.off()
            
        }
    }


# --- DATA across different days (e.g., D3, D5, D7, ...) --------------------------------------------------------------------------------------------
# Error: vector memory exhausted (limit reached?)

# https://stackoverflow.com/questions/51295402/r-on-macos-error-vector-memory-exhausted-limit-reached
#Step 1: Open terminal
#Step 2: 
#    cd ~
#    touch .Renviron
#    open .Renviron
#Step 3: Save the following as the first line of .Renviron
#    R_MAX_VSIZE=100Gb

# https://stackoverflow.com/questions/11624885/remove-multiple-objects-with-rm

for (i in 1:length(disease)) {
    
    gc() # Garbage Collection
    
    mainDir <- paste(current_folder,'/Result', sep = '')
    subDir <- 'DATA'
    if (! file.exists(file.path(mainDir, subDir))) {
        dir.create(file.path(mainDir, subDir))
    }
    
    initial.idx = 1
    initial2end.idx <- c()
    df.all <- c()
    for (s in 1:length(sample)) {
        
        if (sample[s] != "HPNE" & sample[s] != "ch02") {
            
            df.new <- eval(parse(text=paste(sample[s],'.',disease[i],'.','df_new', sep = '')))
            
                df.all <- merge(df.all, df.new, by = "row.names", all = TRUE, sort = TRUE)
                rownames(df.all) <- df.all$Row.names
                df.all$Row.names <- NULL
            
        }
    }
    
        if (ncol(df.all) > 0) {
            
            mainDir <- paste(current_folder,'/Result','/','DATA', sep = '')
            subDir <- disease[i]
            if (! file.exists(file.path(mainDir, subDir))) {
                dir.create(file.path(mainDir, subDir))
            }
            
            #df.all[is.na(df.all)] <- 0
            
            row.idx <- unlist(sapply(all_gene_name, FUN = function(X) which(rownames(df.all) %in% X)))
            df.all <- df.all[row.idx, , drop = FALSE]
            
            assign(paste('DATA','.',disease[i],'.','df_all', sep = ''), df.all)
            
            
            # --- Data correction and combination ---------------------------------------------------------------------------------------------------------------
            exp.name <- sample[which(sample != "HPNE")]
            
            exp.idx <- lapply(exp.name, FUN = function(X) which(colnames(df.all) %in% colnames(eval(parse(text=paste(X,'.',disease[i],'.','df_new', sep = ''))))))
            names(exp.idx) <- exp.name
            
            exp.num <- sapply(exp.idx, FUN = function(X) length(X))
            
            # /// replace NA with true value ///
            mainDir <- paste(current_folder,'/Result','/','DATA','/',disease[i], sep = '')
            subDir <- 'QC'
            if (! file.exists(file.path(mainDir, subDir))) {
                dir.create(file.path(mainDir, subDir))
            }
            
            initial.idx = 1
            initial2end.idx <- c()
            df.subset.all <- c()
            for (s in 1:length(exp.name)) {
                df <- eval(parse(text=paste(exp.name[s],'.',disease[i],'.','df', sep = '')))
                
                df.all.subset <- df.all[, exp.idx[[s]], drop = FALSE]
                
                if (! is.null(df)) {
                    
                    
                    
                    row.idx <- unlist(sapply(rownames(df.all.subset), FUN = function(X) which(rownames(df) %in% X)))
                    col.idx <- unlist(sapply(colnames(df.all.subset), FUN = function(X) which(colnames(df) %in% X)))
                    
                    df.subset <- df[row.idx, col.idx, drop = FALSE]
                    
                    #assign(paste('DATA','.',disease[i],'.',exp.name[s],'.','df_subset', sep = ''), df.subset)
                    #write.csv(df.subset, file = paste(current_folder,'/Result','/','DATA','/',disease[i],'/','QC','/',disease[i],'.QC','.',exp.name[s],'.csv', sep = ''), row.names = TRUE, quote = TRUE)
                    
                    df.subset.all <- merge(df.subset.all, df.subset, by = "row.names", all = TRUE, sort = TRUE)
                    rownames(df.subset.all) <- df.subset.all$Row.names
                    df.subset.all$Row.names <- NULL
                }
            }
            row.idx <- unlist(sapply(all_gene_name, FUN = function(X) which(rownames(df.subset.all) %in% X)))
            df.subset.all <- df.subset.all[row.idx, , drop = FALSE]
            assign(paste('DATA','.',disease[i],'.','df_subset_all', sep = ''), df.subset.all)
            #write.csv(df.subset.all, file = paste(current_folder,'/Result','/','DATA','/',disease[i],'/',disease[i],'.QC','.csv', sep = ''), row.names = TRUE, quote = TRUE)
            
            for (s in 1:length(exp.name)) {
                df.new <- eval(parse(text=paste(exp.name[s],'.',disease[i],'.','df_new', sep = '')))
                
                if (! is.null(df.new)) {
                    row.idx <- which(rownames(df.subset.all) %in% rownames(df.new))
                    col.idx <- which(colnames(df.subset.all) %in% colnames(df.new))
                    
                    df.subset <- df.subset.all[row.idx, col.idx, drop = FALSE]
                    
                    write.csv(df.subset, file = paste(current_folder,'/Result','/','DATA','/',disease[i],'/','QC','/',disease[i],'.QC','.',exp.name[s],'.csv', sep = ''), row.names = TRUE, quote = TRUE)
                }
            }
            
            # /// combine library and control ///
            mainDir <- paste(current_folder,'/Result','/','DATA','/',disease[i], sep = '')
            subDir <- 'control'
            if (! file.exists(file.path(mainDir, subDir))) {
                dir.create(file.path(mainDir, subDir))
            }
            
            df.subset.all.control <- eval(parse(text=paste('HPNE','.','df_subset_all', sep = '')))
            
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
                            
                            sample.name <- names(which(sapply(exp.idx, FUN = function(X) df.subset.all.new.dup.idx[d] %in% unlist(X))))
                            colnames(df.subset.all.new)[df.subset.all.new.dup.idx[d]] <- paste(colnames(df.subset.all.new)[df.subset.all.new.dup.idx[d]], '.', sample.name, sep = '')
                            
                    }
                }
                if (length(df.subset.all.control.dup.idx) > 0) {
                    colnames(df.subset.all.control)[df.subset.all.control.dup.idx] <- paste(colnames(df.subset.all.control)[df.subset.all.control.dup.idx], '.', 'HPNE', sep = '')
                }
            }
            
            df2control <- merge(df.subset.all.new, df.subset.all.control, by = "row.names", all = TRUE, sort = TRUE)
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
            
            exp.idx1 <- lapply(control, FUN = function(X) which(colnames(df2control) %in% colnames(eval(parse(text=paste('HPNE','.',X,'.','df_new', sep = ''))))))
            names(exp.idx1) <- control
            
            
            
            
            
            
            
            
            
            exp.idx2 <- lapply(exp.name, FUN = function(X) which(colnames(df2control) %in% colnames(eval(parse(text=paste(X,'.',disease[i],'.','df_new', sep = ''))))))
            names(exp.idx2) <- paste(exp.name, '.', disease[i], sep = '')
            
            exp.idx12 <- c(exp.idx1, exp.idx2)
            
            for (k in 1:length(exp.idx12)) {
                col.idx <- unlist(exp.idx12[k])
                names(col.idx) <- colnames(df2control)[col.idx]
                
                match.idx <- grep('\\.', names(col.idx), value = FALSE)
                match.name <- grep('\\.', names(col.idx), value = TRUE)
                
                if (length(match.name) > 0) {
                    match.idx <- match.idx[which(sapply(match.name, FUN = function(X) unlist(strsplit(X, '\\.'))[[2]] %in% sample))]
                    if (length(match.idx) > 0) {
                        names(col.idx)[match.idx] = sapply(names(col.idx)[match.idx], FUN = function(X) unlist(strsplit(X, '\\.'))[[1]])
                    }
                }
                
                if (names(exp.idx12[k]) %in% control) {
                    sample.name <- "HPNE"
                    sample_file <- names(exp.idx12[k])
                    row.idx <- row.idx1
                } else {
                    sample.name <- unlist(strsplit(names(exp.idx12[k]), '\\.'))[1]
                    sample_file <- unlist(strsplit(names(exp.idx12[k]), '\\.'))[2]
                    row.idx <- row.idx2
                }
                
                df <- eval(parse(text=paste(sample.name,'.',sample_file,'.','df', sep = '')))
                
                if (! is.null(df)) {
                    rownames(df)[which(rownames(df) == "Tead2")] = "Tead2(ENSMUSG00000030796)"
                    rownames(df)[which(rownames(df) == "l7Rn6")] = "l7Rn6(ENSMUSG00000062797)"
                    
                    NA.idx <- which(! rownames(df2control) %in% rownames(df))
                    df2control[NA.idx, col.idx][is.na(df2control[NA.idx, col.idx])] <- 0
                    
                    
                    
                    
                    df2control[names(row.idx), col.idx] = df[row.idx, names(col.idx)]
                }
            }
            
            #df2control[is.na(df2control)] <- 0
            
            row.idx <- which(rownames(df2control) %in% intersect(all_gene_name.new, control_gene_name))
            df2control <- df2control[row.idx, , drop = FALSE]
            
            #assign(paste('DATA','.',disease[i],'.','df2control', sep = ''), df2control)
            #write.csv(df2control, file = paste(current_folder,'/Result','/','DATA','/',disease[i],'/',disease[i],'.control','.csv', sep = ''), row.names = TRUE, quote = TRUE)
            
            # /// split into individual library and control ///
            for (k in 1:length(exp.idx12)) {
                col.idx <- unlist(exp.idx12[k])
                names(col.idx) <- colnames(df2control)[col.idx]
                
                sample.name <- unlist(strsplit(names(exp.idx12[k]), '\\.'))[1]
                
                if (length(col.idx) > 0) {
                    df2control.subset <- df2control[, col.idx, drop = FALSE]
                    #assign(paste('DATA','.',disease[i],'.',sample.name,'.','df2control_subset', sep = ''), df2control.subset)
                    write.csv(df2control.subset, file = paste(current_folder,'/Result','/','DATA','/',disease[i],'/','control','/',disease[i],'.control','.',sample.name,'.csv', sep = ''), row.names = TRUE, quote = TRUE)
                }
            }
        }
}

