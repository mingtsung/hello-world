Monocle3_analysis <- function(expression_matrix, cell_metadata, gene_metadata, cds = NULL, cell_subset = NULL, gene_subset = NULL, cell_choose = FALSE, batch = NULL, cell_type = NULL, time = NULL, cells = NULL, project = "Monocle3_project", analysis_type = c("clustering", "trajectories"), redDim_method = c("UMAP", "tSNE"), redDim_minDim = 2, redDim_maxDim = 10, top_gene_num = 3, genes_of_interest = NULL, signature_genes = NULL) {
    
    current_folder = "."
    setwd(current_folder)
    
    
    #project = "tutorial"
    
    project.name = unlist(strsplit(project, '\\.'))[1]
    if (length(unlist(strsplit(project, '\\.'))) > 1) {
        project.subname <- c()
        for (pn in 2:length(unlist(strsplit(project, '\\.')))) {
            project.subname <- c(project.subname, unlist(strsplit(project, '\\.'))[pn])
        }
        project.subname <- paste(project.subname, collapse = '.')
    }
    
    if (project.name %in% "GSM4150378") {
        if (exists("project.subname")) {
            if (project.subname %in% c("A549", "K562", "MCF7")) {
                graph.run = 1
            } else {
                graph.run = 0
            }
        } else {
            graph.run = 0
        }
    } else {
        graph.run = 1
    }
    
    if (project.name %in% c("GSM4150377", "GSM4150378")) {
        if (project.name %in% "GSM4150377") {
            dose_character.color <- c("0" = "gray", 
                                      "0.1" = "#1D1147FF", 
                                      "0.5" = "#51127CFF", 
                                      "1" = "#822681FF", 
                                      "5" = "#B63679FF", 
                                      "10" = "#E65164FF", 
                                      "50" = "#FB8861FF", 
                                      "100" = "#FEC287FF")
            
            dose_character.unit <- paste("Dose"," ","[µM]", sep = '')
        } else if (project.name %in% "GSM4150378") {
            dose_character.color <- c("0" = "gray", 
                                      "10" = "#1D1147FF", 
                                      "100" = "#822681FF", 
                                      "1000" = "#E65164FF", 
                                      "10000" = "#FEC287FF")
            
            dose_character.unit <- paste("Dose"," ","[nM]", sep = '')
        }
    } else {
        dose_character.color <- ""
        dose_character.unit <- ""
    }
    
    
    
    # === Monocle =======================================================================================================================================
    # Monocle 3: An analysis toolkit for single-cell RNA-seq. https://cole-trapnell-lab.github.io/monocle3/
    # Monocle 3 is designed for use with absolute transcript counts (e.g. from UMI experiments).
    # Monocle 3 works "out-of-the-box" with the transcript count matrices produced by Cell Ranger, the software pipeline for analyzing experiments from the 10X Genomics Chromium instrument.
    # Monocle 3 also works well with data from other RNA-Seq workflows such as sci-RNA-Seq (single-cell combinatorial indexing RNA sequencing) and instruments like the Biorad ddSEQ.
    
    
    
    # === Major updates in Monocle 3 ====================================================================================================================
    # https://cole-trapnell-lab.github.io/monocle3/docs/updates/
    
    # Most of the algorithmic details in Monocle 3 are described in Cao & Spielmann et al.
    # Cao, Junyue, et al. "The single-cell transcriptional landscape of mammalian organogenesis." Nature 566.7745 (2019): 496-502.
    # https://www.nature.com/articles/s41586-019-0969-x
    # https://cole-trapnell-lab.github.io/pdfs/papers/cao-spielmann-mouse-emb.pdf
    # Mouse Organogenesis Cell Atlas (MOCA): https://oncoscape.v3.sttrcancer.org/atlas.gs.washington.edu.mouse.rna/landing
    
    
    
    # === Citations and Acknowledgements ================================================================================================================
    # https://cole-trapnell-lab.github.io/monocle3/docs/citations/
    
    
    
    # === Installing Monocle 3 ==========================================================================================================================
    # https://cole-trapnell-lab.github.io/monocle3/docs/installation/
    
    #if (!requireNamespace("BiocManager", quietly = TRUE))
    #    install.packages("BiocManager")
    #BiocManager::install()
    
    #BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
    #                       'limma', 'S4Vectors', 'SingleCellExperiment',
    #                       'SummarizedExperiment', 'batchelor', 'Matrix.utils'))
    
    #install.packages("devtools")
    #devtools::install_github('cole-trapnell-lab/leidenbase')
    #devtools::install_github('cole-trapnell-lab/monocle3')
    
    #devtools::install_github('cole-trapnell-lab/monocle3', ref="develop") # the develop branch of monocle3
    
    #install.packages("ggplot2")        # ggplot2: Create Elegant Data Visualisations Using the Grammar of Graphics
    #install.packages("dplyr")          # dplyr: A Grammar of Data Manipulation
    #install.packages("Matrix")         # Matrix: Sparse and Dense Matrix Classes and Methods
    #install.packages("gtools")         # gtools: Various R Programming Tools
    
    suppressPackageStartupMessages({
        #library(monocle3)
        library(ggplot2)               # ggplot()
        library(dplyr)                 # mutate() # imported for some downstream data manipulation; tbl_df()  # `tbl_df()` was deprecated in dplyr 1.0.0. Please use `tibble::as_tibble()` instead.
        library(Matrix)                # To work with your data in a sparse format, simply provide it to Monocle 3 as a sparse matrix from the Matrix package
        library(gtools)                # mixedsort(), mixedorder()
    })
    
    
    
    # === source monocle3 ===============================================================================================================================
    #for (f in 1:length(list.files(path = "./Script/function/Monocle3/monocle3/monocle3-master/R", pattern = ".R", full.names = TRUE, recursive = FALSE))) {
    #    source(list.files(path = "./Script/function/Monocle3/monocle3/monocle3-master/R", pattern = ".R", full.names = TRUE, recursive = FALSE)[f])
    #}
    #source("./Script/function/Monocle3/monocle3/monocle3-master/R/plotting_modified.R")
    
    #install.packages('./Script/function/Monocle3/leidenbase/leidenbase-master', repos = NULL, type = "source")
    #install.packages('./Script/function/Monocle3/monocle3/monocle3-master', repos = NULL, type = "source")
    
    #install.packages("devtools")       # devtools: Tools to Make Developing R Packages Easier
    
    suppressPackageStartupMessages({
        library(devtools)              # load_all()
    })
    
    load_all('./Script/function/Monocle3/leidenbase/leidenbase-master', quiet = TRUE)
    load_all('./Script/function/Monocle3/monocle3/monocle3-master', quiet = TRUE)
    
    
    
    # === User defined function =========================================================================================================================
    append_umap_coordinates <- function(cds) {
        #num_dim = dim(cds@reducedDims[["UMAP"]])[2]
        num_dim = dim(reducedDims(cds)[["UMAP"]])[2]
        for (i in seq(1, num_dim)) {
            new_col_name = paste("umap", i, sep = "")
            #colData(cds)[[new_col_name]] = cds@reducedDims[["UMAP"]][, i]
            colData(cds)[[new_col_name]] = reducedDims(cds)[["UMAP"]][, i]
        }
        return(cds)
    }
    
    
    
    # === Install packages ==============================================================================================================================
    #install.packages("plotly")         # plotly: Create Interactive Web Graphics via 'plotly.js'
    #install.packages("htmlwidgets")    # htmlwidgets: HTML Widgets for R
    #install.packages("webshot")        # webshot: Take Screenshots of Web Pages
    #install.packages("stringr")        # stringr: Simple, Consistent Wrappers for Common String Operations
    #install.packages("reshape2")       # reshape2: Flexibly Reshape Data: A Reboot of the Reshape Package
    #install.packages("viridis")        # viridis: Colorblind-Friendly Color Maps for R
    #install.packages("hrbrthemes")     # hrbrthemes: Additional Themes and Theme Components for 'ggplot2'
    
    suppressPackageStartupMessages({
        library(plotly)
        library(htmlwidgets)
        library(webshot)
        library(stringr)               # str_to_title(), str_to_upper()
        library(reshape2)              # melt()
        library(viridis)               # scale_color_viridis()
        library(hrbrthemes)            # theme_ipsum()
    })
    
    
    
    # === plot ==========================================================================================================================================
    # /// 2d ///
    # https://ggrepel.slowkow.com/articles/examples.html
    op <- options()
    options(ggrepel.max.overlaps = Inf)
    
    # /// 3d ///
    # https://plotly.com/r/axes/
    axis <- list(title = "", zeroline = FALSE, showline = FALSE, showticklabels = FALSE, showgrid = FALSE) # Hide Axes Title, Lines, Ticks, and Labels
    scene <- list(xaxis = axis, yaxis = axis, zaxis = axis) # Modifying Axes for 3D Plots
    
    
    
    # === input argument ================================================================================================================================
    cds.input = cds
    
    # /// subset ///
    cell_subset.input = cell_subset
    gene_subset.input = gene_subset
    
    if (! is.null(cell_subset)) {
        cell_subset.idx <- which(sapply(cell_subset, FUN = function(X) X %in% rownames(cell_metadata)))
        if (length(cell_subset.idx) > 0) {
            cell_subset <- cell_subset[cell_subset.idx]
        } else {
            cell_subset <- NULL
        }
    }
    
    if (length(cell_subset) == 0 & cell_choose == FALSE) {
        SUB.run = 0
    } else if (length(cell_subset) > 0 & cell_choose == FALSE) {
        SUB.run = c(0, 1)
    } else if (length(cell_subset) == 0 & cell_choose == TRUE) {
        SUB.run = c(0, 2)
    } else if (length(cell_subset) > 0 & cell_choose == TRUE) {
        SUB.run = c(0, 1, 2)
    }
    
    batch.input = batch
    cell_type.input = cell_type
    
    time.input = time
    stage.input = time
    
    cells.input = cells
    
    analysis_type.input = analysis_type
    
    redDim_method.input = redDim_method
    redDim.method <- c("UMAP", "tSNE")
    #print(redDim.method)
    redDim.method.idx <- which(sapply(redDim.method, FUN = function(X) X %in% redDim_method.input))
    redDim.method = redDim.method[redDim.method.idx]
    #print(redDim.method.idx)
    #print(redDim.method)
    
    if (redDim_maxDim < 3) {
        redDim_maxDim.input = 3
    } else {
        redDim_maxDim.input = redDim_maxDim
    }
    
    if (! is.null(redDim_minDim)) {
        if (! is.na(redDim_minDim)) {
            if (redDim_minDim < redDim_maxDim.input | redDim_minDim > redDim_maxDim.input) {
                redDim_minDim = 2
            }
        } else {
            redDim_minDim = redDim_maxDim.input
        }
    } else {
        redDim_minDim = redDim_maxDim.input
    }
    redDim_minDim.input = redDim_minDim
    
    if (top_gene_num < 1) {
        top_gene_num = 1
    } else {
        top_gene_num = ceiling(top_gene_num)
    }
    
    # *** query gene(s) ***
    genes_of_interest.input = genes_of_interest
    genes_of_interest.input.title = names(genes_of_interest)
    genes_of_interest.titles = c(genes_of_interest.input.title, "genes_of_interest")
    
    # *** signature gene(s) ***
    signature_genes.input = signature_genes
    signature_genes.input.title = names(signature_genes)
    signature_genes.titles = c(signature_genes.input.title, "signature_genes")
    
    
    
    # === input data ====================================================================================================================================
    mainDir <- current_folder
    subDir <- "Result"
    if (! file.exists(file.path(mainDir, subDir))) {
        dir.create(file.path(mainDir, subDir))
    }
    mainDir <- paste(current_folder,'/Result', sep = '')
    subDir <- "Monocle3"
    if (! file.exists(file.path(mainDir, subDir))) {
        dir.create(file.path(mainDir, subDir))
    }
    mainDir <- paste(current_folder,'/Result','/','Monocle3', sep = '')
    subDir <- project.name
    if (! file.exists(file.path(mainDir, subDir))) {
        dir.create(file.path(mainDir, subDir))
    }
    
    if (project.name == "tutorial") {
        cell_type = "cao_cell_type"
    } else {
        cell_type = cell_type.input
    }
    
    if (! is.null(cell_type)) {
        cell_type.colidx <- which(sapply(paste(unique(c(cell_type, unlist(strsplit(cell_type, '\\.|_'))[1])),'.','color', sep = ''), FUN = function(X) X %in% colnames(cell_metadata)))
        if (length(cell_type.colidx) > 0) {
            color_set = names(cell_type.colidx)
            color_palette <- list()
            for (c.set in 1:length(color_set)) {
                color_palette.set = list(levels(cell_metadata[,color_set[c.set]])[which(sapply(levels(cell_metadata[,color_set[c.set]]), FUN = function(X) X %in% sort(unique(cell_metadata[,color_set[c.set]]))))])
                names(color_palette.set) = color_set[c.set]
                color_palette <- c(color_palette, color_palette.set)
            }
        } else {
            color_palette <- NULL
        }
    }
    
    
    subdata.type = c("KRAS", # KRAS: KRAS-related pathway genes
                     "PML") # PML: premalignant lesion (Beane, Jennifer E., et al. "Molecular subtyping reveals immune alterations associated with progression of bronchial premalignant lesions." Nature communications 10.1 (2019): 1-13.)
    
    
    #for (subdata.idx in 0:length(subdata.type)) {
    for (subdata.idx in 0:0) {
        
        if (subdata.idx == 0) {
            subdata.folder <- 'ALL'
        } else {
            subdata.folder <- subdata.type[subdata.idx]
        }
        
        
        
        # === Getting started with Monocle 3 ================================================================================================================
        # https://cole-trapnell-lab.github.io/monocle3/docs/starting/
        
        
        if ("clustering" %in% analysis_type.input) {
            
            # === Clustering and classifying your cells =========================================================================================================
            # https://cole-trapnell-lab.github.io/monocle3/docs/clustering/
            
            mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name, sep = '')
            subDir <- 'clustering'
            if (! file.exists(file.path(mainDir, subDir))) {
                dir.create(file.path(mainDir, subDir))
            }
            
            mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering', sep = '')
            subDir <- subdata.folder
            if (! file.exists(file.path(mainDir, subDir))) {
                dir.create(file.path(mainDir, subDir))
            }
            
            mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder, sep = '')
            subDir <- project.subname
            if (! file.exists(file.path(mainDir, subDir))) {
                dir.create(file.path(mainDir, subDir))
            }
            
            if (project.name == "tutorial") {
                batch = "plate"
                cell_type = "cao_cell_type"
            } else {
                batch = batch.input
                cell_type = cell_type.input
            }
            
            if (! is.null(batch)) {
                batch.name = paste('.',gsub("\\.","_",batch),'_','align','.', sep = '')
            } else {
                batch.name = '.'
            }
            
            
            # Load the data
            # expression_matrix, a numeric matrix of expression values, where rows are genes, and columns are cells
            # cell_metadata, a data frame, where rows are cells, and columns are cell attributes (such as cell type, culture condition, day captured, etc.)
            # gene_metadata, an data frame, where rows are features (e.g. genes), and columns are gene attributes, such as biotype, gc content, etc.
            # Additionally: one of the columns of the gene_metadata should be named "gene_short_name", which represents the gene symbol or simple name (generally used for plotting) for each gene.
            
            if (project.name == "tutorial") {
                #expression_matrix <- readRDS(url("http://staff.washington.edu/hpliner/data/cao_l2_expression.rds"))
                #cell_metadata <- readRDS(url("http://staff.washington.edu/hpliner/data/cao_l2_colData.rds"))
                #gene_metadata <- readRDS(url("http://staff.washington.edu/hpliner/data/cao_l2_rowData.rds"))
                expression_matrix <- readRDS(paste(current_folder,'/Script/function','/','Monocle3','/','tutorial','/',"cao_l2_expression.rds", sep = ''))
                cell_metadata <- readRDS(paste(current_folder,'/Script/function','/','Monocle3','/','tutorial','/',"cao_l2_colData.rds", sep = ''))
                gene_metadata <- readRDS(paste(current_folder,'/Script/function','/','Monocle3','/','tutorial','/',"cao_l2_rowData.rds", sep = ''))
            } else {
                if (is.matrix(expression_matrix) | is.data.frame(expression_matrix)) {
                    expression_matrix <- Matrix(as.matrix(expression_matrix), sparse = TRUE)
                } else {
                    expression_matrix <- expression_matrix
                }
                cell_metadata <- cell_metadata
                gene_metadata <- gene_metadata
            }
            expression_matrix[is.na(expression_matrix)] <- 0
            #dim(expression_matrix)
            #dim(cell_metadata)
            #dim(gene_metadata)
            
            #colnames(cell_metadata)
            #colnames(gene_metadata)
            
            # --- max_components = maximum dimension ------------------------------------------------------------------------------------------------------------
            if (redDim_maxDim.input == 2 | redDim_maxDim.input == 3) {
                dim.run.idx1 = 1
                dim.run.idx2 = 2
            } else {
                if (redDim_minDim.input == redDim_maxDim.input) {
                    dim.run.idx1 = 3
                    dim.run.idx2 = 3
                } else {
                    dim.run.idx1 = 1
                    dim.run.idx2 = 3
                }
            }
            
            for (dim.run in dim.run.idx1:dim.run.idx2) {
                
                if (dim.run == 1) {
                    #reduce_dimension() # max_components: the dimensionality of the reduced space. Default is 2.
                    assign(paste('UMAP','.','dim_max', sep = ''), 2)
                    assign(paste('tSNE','.','dim_max', sep = ''), 2)
                } else if (dim.run == 2) {
                    assign(paste('UMAP','.','dim_max', sep = ''), 3)
                    assign(paste('tSNE','.','dim_max', sep = ''), 3)
                } else if (dim.run == 3) {
                    assign(paste('UMAP','.','dim_max', sep = ''), redDim_maxDim.input)
                    assign(paste('tSNE','.','dim_max', sep = ''), 3)
                }
                
                # Store data in a cell_data_set object: https://cole-trapnell-lab.github.io/monocle3/docs/starting/#cell_data_set
                # Make the CDS object
                cds.origin <- new_cell_data_set(expression_matrix,
                                                cell_metadata = cell_metadata,
                                                gene_metadata = gene_metadata)
                
                # Working with large data sets
                # To work with your data in a sparse format, simply provide it to Monocle 3 as a sparse matrix from the Matrix package
                # The output from a number of RNA-Seq pipelines, including Cell Ranger, is already in a sparseMatrix format (e.g. MTX).
                #cds.origin <- new_cell_data_set(as(umi_matrix, "sparseMatrix"),
                #                         cell_metadata = cell_metadata,
                #                         gene_metadata = gene_metadata)
                
                if (subdata.idx == 0) {
                    cds <- cds.origin
                } else {
                    cds <- cds.origin
                    
                    #rowidx <- 
                    #cds <- cds.origin[rowidx, , drop = FALSE]
                }
                
                if (ncol(cds) > 0) {
                    
                    # Combining CDS objects
                    # make a fake second cds object for demonstration
                    #cds2 <- cds[1:100,]
                    #big_cds <- combine_cds(list(cds, cds2))
                    
                    
                    if (ncol(cds) <= 100) {
                        cell_size_2d = 1.5
                        cell_size_3d = 250
                    } else {
                        cell_size_2d = 0.35 # default for plot_cells()
                        cell_size_3d = 25 # default for plot_cells_3d()
                    }
                    
                    for (redDim in 1:length(redDim.method)) {
                        mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',project.subname, sep = '')
                        subDir <- redDim.method[redDim]
                        if (! file.exists(file.path(mainDir, subDir))) {
                            dir.create(file.path(mainDir, subDir))
                        }
                        
                        dim.folder = paste(eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))),'D', sep = '')
                        
                        mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim], sep = '')
                        subDir <- dim.folder
                        if (! file.exists(file.path(mainDir, subDir))) {
                            dir.create(file.path(mainDir, subDir))
                        }
                    }
                    
                    
                    ## Step 1: Normalize and pre-process the data *****************************************************
                    # Pre-process the data
                    if (project.name %in% "GSM4150377") {
                        # Identify over-dispersed genes to be used as feature selection
                        dispersion_test = estimateDispersionsForCellDataSet(sciPlex_cds, min_cells_detected = 100, removeOutliers = TRUE)
                        
                        disp_table = dispersionTable(dispersion_test)
                        disp_table = disp_table %>% 
                            mutate(excess_disp = (dispersion_empirical - dispersion_fit) / dispersion_fit) %>% 
                            dplyr::arrange(plyr::desc(excess_disp))
                        
                        unsup_clustering_genes = as.character(head(disp_table, 1000)$gene_id)
                        ordering_genes = unsup_clustering_genes
                        
                        cds <- preprocess_cds(cds, method = "PCA", num_dim = 35, norm_method = "log", residualModelFormulaStr = "~log(n.umi)", use_genes = ordering_genes, verbose = FALSE)
                    } else if (project.name %in% "GSM4150378") {
                        cds <- preprocess_cds(cds, method = "PCA", num_dim = 25, norm_method = "log", verbose = FALSE)
                    } else {
                        if (ncol(cds) / 2 <= 50) {
                            # Warning in (function (A, nv = 5, nu = nv, maxit = 1000, work = nv + 7, reorth = TRUE,  :
                            # You're computing too large a percentage of total singular values, use a standard svd instead.
                            if (floor(ncol(cds) / 2) == ncol(cds) / 2) {
                                cds <- preprocess_cds(cds, method = "PCA", num_dim = floor(ncol(cds) / 2) - 1)
                            } else {
                                cds <- preprocess_cds(cds, method = "PCA", num_dim = floor(ncol(cds) / 2))
                            }
                        } else {
                            cds <- preprocess_cds(cds, method = "PCA", num_dim = 50) # method: Default is "PCA". method = c("PCA", "LSI") # "PCA": Principal Components Analysis (the standard for RNA-seq), "LSI": Latent Semantic Indexing (common in ATAC-seq)
                        }
                    }
                    
                    png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',project.subname,'/',project,'.','clustering','.','plot_pc_variance_explained','.png', sep = ''), width = 6, height = 4, units = 'in', res = 600)
                    p <- plot_pc_variance_explained(cds)
                    print(p)
                    dev.off()
                    
                    
                    ## Step 2: Reduce the dimensions using UMAP *******************************************************
                    # Reduce dimensionality and visualize the cells
                    for (redDim in 1:length(redDim.method)) {
                        #reduce_dimension() # max_components: the dimensionality of the reduced space. Default is 2.
                        #umap() # n_components: The dimension of the space to embed into. This defaults to 2 to provide easy visualization, but can reasonably be set to any integer value in the range 2 to 100.
                        # UMAP: Monocle 3 uses UMAP by default, as we feel that it is both faster and better suited for clustering and trajectory analysis in RNA-seq.
                        # tSNE: Error in .check_tsne_params(nrow(X), dims = dims, perplexity = perplexity,  : dims should be either 1, 2 or 3
                        cds <- reduce_dimension(cds, reduction_method = redDim.method[redDim], max_components = eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))))
                        
                        #cds@int_colData@listData$reducedDims@listData[[redDim.method[redDim]]]
                        
                        #if (redDim.method[redDim] == "UMAP") {
                        #    # Faster clustering with UMAP: 
                        #    # If you have a relatively large dataset (with >10,000 cells or more), you may want to take advantage of options that can accelerate UMAP.
                        #    # However, invoking reduce_dimension() with either of these options will make it produce slighly different output each time you run it.
                        #    cds <- reduce_dimension(cds, reduction_method = redDim.method[redDim], max_components = eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))), umap.fast_sgd = TRUE, cores = 2)
                        #}
                        
                        dim.folder = paste(eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))),'D', sep = '')
                        dim.finalX = eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) - (2-1)
                        
                        if (dim.finalX >= 1) {
                            
                            for (dim.comb in 1:dim.finalX) {
                                dim.x <- dim.comb
                                dim.y <- dim.x + 1
                                
                                if (dim.finalX > 1) {
                                    if (dim.x == 1) {
                                        comb.run = 2
                                    } else {
                                        comb.run = 1
                                    }
                                } else {
                                    comb.run = 1
                                }
                                
                                for (dim.comb.run in 1:comb.run) {
                                    if (dim.comb.run > 1) {
                                        dim.y <- dim.y + 1
                                    }
                                    
                                    if (eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) == 2) {
                                        dim.xy <- paste('D',dim.x,'-','D',dim.y, sep = '')
                                    } else {
                                        dim.xy <- paste('D',dim.x,'-','D',dim.y, sep = '')
                                    }
                                    #print(dim.xy)
                                    
                                    mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder, sep = '')
                                    subDir <- dim.xy
                                    if (! file.exists(file.path(mainDir, subDir))) {
                                        dir.create(file.path(mainDir, subDir))
                                    }
                                    
                                    # /// no color_cells ///
                                    png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',project,'.','clustering','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,'.png', sep = ''), width = 4.075, height = 4, units = 'in', res = 600)
                                    p <- plot_cells(cds, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d)
                                    print(p)
                                    dev.off()
                                    
                                    # /// color_cells_by = cell_type ///
                                    if (length(cell_type) > 0) {
                                        for (color.plt in 1:length(cell_type)) {
                                            legend.max_nchar = max(nchar(levels(colData(cds)[, cell_type[color.plt]])))
                                            legend.multiply_factor = 0.5 + legend.max_nchar * 0.05
                                            legend.num = length(levels(colData(cds)[, cell_type[color.plt]]))
                                            legend.ncol = ceiling(legend.num / 10)
                                            width.size = 4.075 + (legend.multiply_factor * legend.ncol)
                                            
                                            png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',project,'.','clustering','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,'.',gsub("\\.","_",cell_type[color.plt]),'.png', sep = ''), width = 4.075, height = 4, units = 'in', res = 600)
                                            p <- plot_cells(cds, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = cell_type[color.plt], label_groups_by_cluster = FALSE) # Finally, you can disable this labeling policy, making plot_cells() behave like it did before we called cluster_cells()
                                            if (! is.null(color_palette)) {
                                                cell_type.colidx <- which(paste(unique(c(cell_type[color.plt], unlist(strsplit(cell_type[color.plt], '\\.|_'))[1])),'.','color', sep = '') %in% names(color_palette))
                                                if (length(cell_type.colidx) == 1) {
                                                    p <- p + scale_color_manual(values = color_palette[[which(names(color_palette) %in% paste(unique(c(cell_type[color.plt], unlist(strsplit(cell_type[color.plt], '\\.|_'))[1])),'.','color', sep = ''))]])
                                                }
                                            }
                                            print(p)
                                            dev.off()
                                            
                                            png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',project,'.','clustering','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,'.',gsub("\\.","_",cell_type[color.plt]),'_','legend','.png', sep = ''), width = width.size, height = 4, units = 'in', res = 600)
                                            p <- plot_cells(cds, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = cell_type[color.plt], label_groups_by_cluster = FALSE, label_cell_groups = FALSE) # Finally, you can disable this labeling policy, making plot_cells() behave like it did before we called cluster_cells()
                                            if (! is.null(color_palette)) {
                                                cell_type.colidx <- which(paste(unique(c(cell_type[color.plt], unlist(strsplit(cell_type[color.plt], '\\.|_'))[1])),'.','color', sep = '') %in% names(color_palette))
                                                if (length(cell_type.colidx) == 1) {
                                                    p <- p + scale_color_manual(values = color_palette[[which(names(color_palette) %in% paste(unique(c(cell_type[color.plt], unlist(strsplit(cell_type[color.plt], '\\.|_'))[1])),'.','color', sep = ''))]])
                                                }
                                            }
                                            p <- p + theme(legend.position = "right")
                                            p$guides$colour$ncol = legend.ncol
                                            if (project.name %in% c("GSE139944", c("GSM4150376", "GSM4150377", "GSM4150378", "GSM4150379"))) {
                                                if (cell_type[color.plt] == "pathway_level_1") {
                                                    p$guides$colour$title = "Pathway"
                                                } else if (cell_type[color.plt] == "product_name") {
                                                    p$guides$colour$title = "Treatment"
                                                } else if (cell_type[color.plt] == "dose_character") {
                                                    p$guides$colour$title = dose_character.unit
                                                } else {
                                                    p$guides$colour$title = "Treatment"
                                                }
                                            }
                                            print(p)
                                            dev.off()
                                        }
                                    }
                                }
                            }
                        }
                    }
                    
                    if (project.name == "tutorial") {
                        # /// marker.genes ///
                        gene_set.list <- list()
                        
                        example_genes <- list(unique(c("cpna-2", "egl-21", "ram-2", "inos-1")))
                        if (length(unlist(example_genes)) == 1) {
                            names(example_genes) = "example_gene"
                        } else if (length(unlist(example_genes)) > 1) {
                            names(example_genes) = "example_genes"
                        }
                        gene_set.list <- c(gene_set.list, example_genes)
                        
                        if (length(gene_set.list) > 0) {
                            for (N in 1:length(gene_set.list)) {
                                #print(names(gene_set.list[N]))
                                
                                marker.genes <- gene_set.list[[N]]
                                marker.genes.idx <- which(sapply(marker.genes, FUN = function(X) X %in% rowData(cds)[,"gene_short_name"]))
                                if (length(marker.genes.idx) > 0) {
                                    marker.genes = marker.genes[marker.genes.idx]
                                    
                                    panel.num = length(unlist(unique(marker.genes)))
                                    if (panel.num == 1) {
                                        width.size = 7
                                        height.size = 6
                                    } else if (panel.num == 2) {
                                        width.size = 7 + (5 * 1)
                                        height.size = 6
                                    } else if (panel.num == 3) {
                                        width.size = 7 + (5 * 2)
                                        height.size = 6
                                    } else if (panel.num > 3) {
                                        width.size = 7 + (5 * 3)
                                        height.size = width.size - 1
                                    }
                                    
                                    for (redDim in 1:length(redDim.method)) {
                                        
                                        # /// 2d ///
                                        dim.folder = paste(eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))),'D', sep = '')
                                        dim.finalX = eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) - (2-1)
                                        
                                        if (dim.finalX >= 1) {
                                            
                                            for (dim.comb in 1:dim.finalX) {
                                                dim.x <- dim.comb
                                                dim.y <- dim.x + 1
                                                
                                                if (dim.finalX > 1) {
                                                    if (dim.x == 1) {
                                                        comb.run = 2
                                                    } else {
                                                        comb.run = 1
                                                    }
                                                } else {
                                                    comb.run = 1
                                                }
                                                
                                                for (dim.comb.run in 1:comb.run) {
                                                    if (dim.comb.run > 1) {
                                                        dim.y <- dim.y + 1
                                                    }
                                                    
                                                    if (eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) == 2) {
                                                        dim.xy <- paste('D',dim.x,'-','D',dim.y, sep = '')
                                                    } else {
                                                        dim.xy <- paste('D',dim.x,'-','D',dim.y, sep = '')
                                                    }
                                                    #print(dim.xy)
                                                    
                                                    png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',project,'.','clustering','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,'.',gsub("\\.","_",names(gene_set.list[N])),'.png', sep = ''), width = width.size, height = height.size, units = 'in', res = 600)
                                                    p <- plot_cells(cds, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, genes = marker.genes, norm_method = "log")
                                                    print(p)
                                                    dev.off()
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    
                    
                    ## Step 3: Remove batch effects with cell alignment ***********************************************
                    # (Optional) Remove batch effects: https://cole-trapnell-lab.github.io/monocle3/docs/clustering/#batch-effects
                    
                    # Check for and remove batch effects
                    # You should always check for batch effects when you perform dimensionality reduction.
                    if (! is.null(batch)) {
                        if (project.name == "tutorial") {
                            width.size = 4.9
                            height.size = 4
                        } else {
                            legend.max_nchar = max(nchar(levels(colData(cds)[, batch])))
                            legend.multiply_factor = 0.5 + legend.max_nchar * 0.05
                            legend.num = length(levels(colData(cds)[, batch]))
                            legend.ncol = ceiling(legend.num / 10)
                            width.size = 4.075 + (legend.multiply_factor * legend.ncol)
                            height.size = 4
                        }
                        
                        for (redDim in 1:length(redDim.method)) {
                            dim.folder = paste(eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))),'D', sep = '')
                            dim.finalX = eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) - (2-1)
                            
                            if (dim.finalX >= 1) {
                                
                                for (dim.comb in 1:dim.finalX) {
                                    dim.x <- dim.comb
                                    dim.y <- dim.x + 1
                                    
                                    if (dim.finalX > 1) {
                                        if (dim.x == 1) {
                                            comb.run = 2
                                        } else {
                                            comb.run = 1
                                        }
                                    } else {
                                        comb.run = 1
                                    }
                                    
                                    for (dim.comb.run in 1:comb.run) {
                                        if (dim.comb.run > 1) {
                                            dim.y <- dim.y + 1
                                        }
                                        
                                        if (eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) == 2) {
                                            dim.xy <- paste('D',dim.x,'-','D',dim.y, sep = '')
                                        } else {
                                            dim.xy <- paste('D',dim.x,'-','D',dim.y, sep = '')
                                        }
                                        #print(dim.xy)
                                        
                                        png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',project,'.','clustering','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,'.',gsub("\\.","_",batch),'.png', sep = ''), width = width.size, height = height.size, units = 'in', res = 600)
                                        p <- plot_cells(cds, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = batch, label_cell_groups = FALSE)
                                        p <- p + theme(legend.position = "right")
                                        p$guides$colour$ncol = legend.ncol
                                        print(p)
                                        dev.off()
                                    }
                                }
                            }
                        }
                        
                        # When run with the alignment_group argument, align_cds() tries to remove batch effects using mutual nearest neighbor alignment, a technique introduced by John Marioni's lab. 
                        # Monocle 3 does so by calling Aaron Lun's excellent package batchelor.
                        # batchelor: Single-Cell Batch Correction Methods (https://bioconductor.org/packages/release/bioc/html/batchelor.html)
                        cds <- suppressMessages(
                            align_cds(cds, alignment_group = batch, num_dim = 50) # to using the alignment_group argument to align_cds(), which aligns groups of cells (i.e. batches)
                        )
                        
                        for (redDim in 1:length(redDim.method)) {
                            dim.folder = paste(eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))),'D', sep = '')
                            dim.finalX = eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) - (2-1)
                            
                            if (dim.finalX >= 1) {
                                
                                for (dim.comb in 1:dim.finalX) {
                                    dim.x <- dim.comb
                                    dim.y <- dim.x + 1
                                    
                                    if (dim.finalX > 1) {
                                        if (dim.x == 1) {
                                            comb.run = 2
                                        } else {
                                            comb.run = 1
                                        }
                                    } else {
                                        comb.run = 1
                                    }
                                    
                                    for (dim.comb.run in 1:comb.run) {
                                        if (dim.comb.run > 1) {
                                            dim.y <- dim.y + 1
                                        }
                                        
                                        if (eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) == 2) {
                                            dim.xy <- paste('D',dim.x,'-','D',dim.y, sep = '')
                                        } else {
                                            dim.xy <- paste('D',dim.x,'-','D',dim.y, sep = '')
                                        }
                                        #print(dim.xy)
                                        
                                        png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',project,'.','clustering','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'png', sep = ''), width = width.size, height = height.size, units = 'in', res = 600)
                                        p <- plot_cells(cds, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = batch, label_cell_groups = FALSE)
                                        p <- p + theme(legend.position = "right")
                                        p$guides$colour$ncol = legend.ncol
                                        print(p)
                                        dev.off()
                                    }
                                }
                            }
                        }
                    }
                    
                    
                    ## Step 4: Cluster the cells **********************************************************************
                    # Cluster your cells: https://cole-trapnell-lab.github.io/monocle3/docs/clustering/
                    
                    cds <- cluster_cells(cds, reduction_method = "PCA")
                    colData(cds)$Cluster = clusters(cds, reduction_method = "PCA")
                    
                    for (redDim in 1:length(redDim.method)) {
                        # Group cells into clusters
                        # Monocle uses a technique called community detection to group cells. 
                        # This approach was introduced by Levine et al as part of the phenoGraph algorithm. (https://github.com/jacoblevine/PhenoGraph)
                        cds <- cluster_cells(cds, reduction_method = redDim.method[redDim]) # resolution: Parameter that controls the resolution of clustering. If NULL (Default), the parameter is determined automatically.
                        
                        if (redDim.method[redDim] == "UMAP") {
                            colData(cds)$louvain_component = cds@clusters[[redDim.method[redDim]]]$partitions
                        }
                        
                        # /// 2d ///
                        dim.folder = paste(eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))),'D', sep = '')
                        dim.finalX = eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) - (2-1)
                        
                        if (dim.finalX >= 1) {
                            
                            for (dim.comb in 1:dim.finalX) {
                                dim.x <- dim.comb
                                dim.y <- dim.x + 1
                                
                                if (dim.finalX > 1) {
                                    if (dim.x == 1) {
                                        comb.run = 2
                                    } else {
                                        comb.run = 1
                                    }
                                } else {
                                    comb.run = 1
                                }
                                
                                for (dim.comb.run in 1:comb.run) {
                                    if (dim.comb.run > 1) {
                                        dim.y <- dim.y + 1
                                    }
                                    
                                    if (eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) == 2) {
                                        dim.xy <- paste('D',dim.x,'-','D',dim.y, sep = '')
                                    } else {
                                        dim.xy <- paste('D',dim.x,'-','D',dim.y, sep = '')
                                    }
                                    #print(dim.xy)
                                    
                                    mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder, sep = '')
                                    subDir <- dim.xy
                                    if (! file.exists(file.path(mainDir, subDir))) {
                                        dir.create(file.path(mainDir, subDir))
                                    }
                                    
                                    legend.max_nchar = max(nchar(levels(cds@clusters[[redDim.method[redDim]]]$clusters)))
                                    legend.multiply_factor = 0.5 + legend.max_nchar * 0.05
                                    legend.num = length(levels(cds@clusters[[redDim.method[redDim]]]$clusters))
                                    legend.ncol = ceiling(legend.num / 10)
                                    width.size = 4.075 + (legend.multiply_factor * legend.ncol)
                                    
                                    png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',project,'.','clustering','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'cluster','.png', sep = ''), width = 4.075, height = 4, units = 'in', res = 600)
                                    p <- plot_cells(cds, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d) # call plot_cells() with no arguments, it colors the cells by cluster according to default.
                                    print(p)
                                    dev.off()
                                    
                                    png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',project,'.','clustering','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'cluster','_','legend','.png', sep = ''), width = width.size, height = 4, units = 'in', res = 600)
                                    p <- plot_cells(cds, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d) # call plot_cells() with no arguments, it colors the cells by cluster according to default.
                                    p <- p + theme(legend.position = "right")
                                    p$guides$colour$ncol = legend.ncol
                                    print(p)
                                    dev.off()
                                    
                                    if (length(cell_type) > 0) {
                                        for (color.plt in 1:length(cell_type)) {
                                            legend.max_nchar = max(nchar(levels(colData(cds)[, cell_type[color.plt]])))
                                            legend.multiply_factor = 0.5 + legend.max_nchar * 0.05
                                            legend.num = length(levels(colData(cds)[, cell_type[color.plt]]))
                                            legend.ncol = ceiling(legend.num / 10)
                                            width.size = 4.075 + (legend.multiply_factor * legend.ncol)
                                            
                                            png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',project,'.','clustering','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'cluster','.',gsub("\\.","_",cell_type[color.plt]),'.png', sep = ''), width = 4.075, height = 4, units = 'in', res = 600)
                                            p <- plot_cells(cds, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = cell_type[color.plt]) # each cluster is labeled according the most common annotation within it
                                            if (! is.null(color_palette)) {
                                                cell_type.colidx <- which(paste(unique(c(cell_type[color.plt], unlist(strsplit(cell_type[color.plt], '\\.|_'))[1])),'.','color', sep = '') %in% names(color_palette))
                                                if (length(cell_type.colidx) == 1) {
                                                    p <- p + scale_color_manual(values = color_palette[[which(names(color_palette) %in% paste(unique(c(cell_type[color.plt], unlist(strsplit(cell_type[color.plt], '\\.|_'))[1])),'.','color', sep = ''))]])
                                                }
                                            }
                                            print(p)
                                            dev.off()
                                            
                                            png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',project,'.','clustering','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'cluster','.',gsub("\\.","_",cell_type[color.plt]),'_','legend','.png', sep = ''), width = width.size, height = 4, units = 'in', res = 600)
                                            p <- plot_cells(cds, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = cell_type[color.plt], label_cell_groups = FALSE) # each cluster is labeled according the most common annotation within it
                                            if (! is.null(color_palette)) {
                                                cell_type.colidx <- which(paste(unique(c(cell_type[color.plt], unlist(strsplit(cell_type[color.plt], '\\.|_'))[1])),'.','color', sep = '') %in% names(color_palette))
                                                if (length(cell_type.colidx) == 1) {
                                                    p <- p + scale_color_manual(values = color_palette[[which(names(color_palette) %in% paste(unique(c(cell_type[color.plt], unlist(strsplit(cell_type[color.plt], '\\.|_'))[1])),'.','color', sep = ''))]])
                                                }
                                            }
                                            p <- p + theme(legend.position = "right")
                                            p$guides$colour$ncol = legend.ncol
                                            if (project.name %in% c("GSE139944", c("GSM4150376", "GSM4150377", "GSM4150378", "GSM4150379"))) {
                                                if (cell_type[color.plt] == "pathway_level_1") {
                                                    p$guides$colour$title = "Pathway"
                                                } else if (cell_type[color.plt] == "product_name") {
                                                    p$guides$colour$title = "Treatment"
                                                } else if (cell_type[color.plt] == "dose_character") {
                                                    p$guides$colour$title = dose_character.unit
                                                } else {
                                                    p$guides$colour$title = "Treatment"
                                                }
                                            }
                                            print(p)
                                            dev.off()
                                            
                                            png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',project,'.','clustering','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'cluster','.',gsub("\\.","_",cell_type[color.plt]),'.','disableLabel','.png', sep = ''), width = 4.075, height = 4, units = 'in', res = 600)
                                            p <- plot_cells(cds, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = cell_type[color.plt], label_groups_by_cluster = FALSE) # Finally, you can disable this labeling policy, making plot_cells() behave like it did before we called cluster_cells()
                                            if (! is.null(color_palette)) {
                                                cell_type.colidx <- which(paste(unique(c(cell_type[color.plt], unlist(strsplit(cell_type[color.plt], '\\.|_'))[1])),'.','color', sep = '') %in% names(color_palette))
                                                if (length(cell_type.colidx) == 1) {
                                                    p <- p + scale_color_manual(values = color_palette[[which(names(color_palette) %in% paste(unique(c(cell_type[color.plt], unlist(strsplit(cell_type[color.plt], '\\.|_'))[1])),'.','color', sep = ''))]])
                                                }
                                            }
                                            print(p)
                                            dev.off()
                                            
                                            png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',project,'.','clustering','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'cluster','.',gsub("\\.","_",cell_type[color.plt]),'_','legend','.','disableLabel','.png', sep = ''), width = width.size, height = 4, units = 'in', res = 600)
                                            p <- plot_cells(cds, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = cell_type[color.plt], label_groups_by_cluster = FALSE, label_cell_groups = FALSE) # Finally, you can disable this labeling policy, making plot_cells() behave like it did before we called cluster_cells()
                                            if (! is.null(color_palette)) {
                                                cell_type.colidx <- which(paste(unique(c(cell_type[color.plt], unlist(strsplit(cell_type[color.plt], '\\.|_'))[1])),'.','color', sep = '') %in% names(color_palette))
                                                if (length(cell_type.colidx) == 1) {
                                                    p <- p + scale_color_manual(values = color_palette[[which(names(color_palette) %in% paste(unique(c(cell_type[color.plt], unlist(strsplit(cell_type[color.plt], '\\.|_'))[1])),'.','color', sep = ''))]])
                                                }
                                            }
                                            p <- p + theme(legend.position = "right")
                                            p$guides$colour$ncol = legend.ncol
                                            if (project.name %in% c("GSE139944", c("GSM4150376", "GSM4150377", "GSM4150378", "GSM4150379"))) {
                                                if (cell_type[color.plt] == "pathway_level_1") {
                                                    p$guides$colour$title = "Pathway"
                                                } else if (cell_type[color.plt] == "product_name") {
                                                    p$guides$colour$title = "Treatment"
                                                } else if (cell_type[color.plt] == "dose_character") {
                                                    p$guides$colour$title = dose_character.unit
                                                } else {
                                                    p$guides$colour$title = "Treatment"
                                                }
                                            }
                                            print(p)
                                            dev.off()
                                        }
                                    }
                                    
                                    # The cluster_cells() also divides the cells into larger, more well separated groups called partitions, 
                                    # using a statistical test from Alex Wolf et al, introduced as part of their PAGA algorithm. (https://github.com/theislab/paga)
                                    png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',project,'.','clustering','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'partition','.png', sep = ''), width = 4.075, height = 4, units = 'in', res = 600)
                                    p <- plot_cells(cds, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = "partition")
                                    print(p)
                                    dev.off()
                                    
                                    png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',project,'.','clustering','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'partition','.','group_cells','.png', sep = ''), width = 4.075, height = 4, units = 'in', res = 600)
                                    p <- plot_cells(cds, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = "partition", group_cells_by = "partition")
                                    print(p)
                                    dev.off()
                                }
                            }
                        }
                        
                        
                        # Find marker genes expressed by each cluster
                        # Once cells have been clustered, we can ask what genes makes them different from one another.
                        # We could group the cells according to cluster, partition, or any categorical variable in colData(cds).
                        mapping.idx = sapply(rownames(colData(cds)), FUN = function(X) which(names(cds@clusters[[redDim.method[redDim]]]$clusters) %in% X))
                        cds@clusters[[redDim.method[redDim]]]$clusters = cds@clusters[[redDim.method[redDim]]]$clusters[mapping.idx]
                        cds@clusters[[redDim.method[redDim]]]$clusters = factor(cds@clusters[[redDim.method[redDim]]]$clusters, levels = mixedsort(as.character(unique(cds@clusters[[redDim.method[redDim]]]$clusters))))
                        colData(cds)$cluster = cds@clusters[[redDim.method[redDim]]]$clusters
                        
                        mapping.idx = sapply(rownames(colData(cds)), FUN = function(X) which(names(cds@clusters[[redDim.method[redDim]]]$partitions) %in% X))
                        cds@clusters[[redDim.method[redDim]]]$partitions = cds@clusters[[redDim.method[redDim]]]$partitions[mapping.idx]
                        cds@clusters[[redDim.method[redDim]]]$partitions = factor(cds@clusters[[redDim.method[redDim]]]$partitions, levels = mixedsort(as.character(unique(cds@clusters[[redDim.method[redDim]]]$partitions))))
                        colData(cds)$partition = cds@clusters[[redDim.method[redDim]]]$partitions
                        
                        mapping.idx = sapply(rownames(colData(cds)), FUN = function(X) which(rownames(cds@principal_graph_aux[[redDim.method[redDim]]]$pr_graph_cell_proj_closest_vertex) %in% X))
                        cds@principal_graph_aux[[redDim.method[redDim]]]$pr_graph_cell_proj_closest_vertex = cds@principal_graph_aux[[redDim.method[redDim]]]$pr_graph_cell_proj_closest_vertex[mapping.idx, , drop = FALSE]
                        
                        xt.all <- c("cluster", "partition", batch)
                        
                        for (x in 1:length(xt.all)) {
                            
                            xt = xt.all[x]
                            
                            # /// "violin" plot ///
                            if (xt == "cluster") {
                                run = xt
                            } else if (xt == "partition") {
                                run = xt
                            } else if (xt == batch) {
                                run = xt
                            }
                            
                            for (r in 1:length(run)) {
                                
                                group_var = run[r]
                                
                                if (is.factor(colData(cds)[, group_var])) {
                                    
                                    if (length(levels(colData(cds)[, group_var])) > 1) {
                                        
                                        # The data frame marker_test_res contains a number of metrics for how specifically expressed each gene is in each partition.
                                        marker_test_res <- top_markers(cds, reduction_method = redDim.method[redDim], group_cells_by = group_var)
                                        
                                        row.orderidx = order(marker_test_res[,"pseudo_R2"], 
                                                             rev(marker_test_res[,"marker_test_q_value"]), 
                                                             rev(marker_test_res[,"marker_test_p_value"]), 
                                                             marker_test_res[,"marker_score"], 
                                                             marker_test_res[,"mean_expression"], 
                                                             decreasing = TRUE)
                                        marker_test_res <- marker_test_res[row.orderidx, , drop = FALSE]
                                        
                                        #top_specific_markers <- marker_test_res %>% filter(fraction_expressing >= 0.10) %>% group_by(cell_group) %>% top_n(1, pseudo_R2)
                                        top_specific_markers <- marker_test_res %>% filter(fraction_expressing >= 0.10) %>% group_by(cell_group)
                                        
                                        #top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))
                                        #top_specific_marker_names <- unique(top_specific_markers %>% pull(gene_short_name))
                                        
                                        wt = "pseudo_R2"
                                        top_No <- c(1, 3, 5, 10, 20, 30, 50)
                                        top_No.idx <- which(top_No <= top_gene_num)
                                        if (length(top_No.idx) > 0) {
                                            top_No <- top_No[top_No.idx]
                                        } else {
                                            top_No <- c(1)
                                        }
                                        for (No in 1:length(top_No)) {
                                            top_specific_markers.top <- c()
                                            cell_group.id <- mixedsort(unlist(unique(top_specific_markers[,"cell_group"])))
                                            for (LV in 1:length(cell_group.id)) {
                                                LV.idx = which(top_specific_markers[,"cell_group"] == cell_group.id[LV])
                                                top_specific_markers.LV <- top_specific_markers[LV.idx, , drop = FALSE]
                                                
                                                top_specific_markers.LV <- top_specific_markers.LV[order(top_specific_markers.LV[,wt], decreasing = TRUE), ]
                                                top_specific_markers.LV.top <- head(top_specific_markers.LV, top_No[No])
                                                
                                                top_specific_markers.top <- rbind(top_specific_markers.top, top_specific_markers.LV.top)
                                            }
                                            assign(paste('top_specific_markers.top',top_No[No], sep = ''), top_specific_markers.top)
                                            
                                            if (top_No[No] == 1) {
                                                
                                                top_specific_marker_names <- unique(eval(parse(text=paste('top_specific_markers.top',top_No[No], sep = '')))[,"gene_short_name"])
                                                top_specific_marker_names <- unique(as.character(unlist(top_specific_marker_names)))
                                                
                                                # /// marker.genes ///
                                                gene_set.list <- list()
                                                
                                                # *** Select serial list ***
                                                serial_num <- c(1:top_gene_num)
                                                serial_num.idx <- which(serial_num <= length(top_specific_marker_names))
                                                if (length(serial_num.idx) > 0) {
                                                    serial_num <- serial_num[serial_num.idx]
                                                } else {
                                                    serial_num <- c()
                                                }
                                                if (length(serial_num) > 0) {
                                                    for (serial in 1:length(serial_num)) {
                                                        serial_genes <- list(top_specific_marker_names[serial_num[serial]])
                                                        names(serial_genes) = paste('gene',serial_num[serial], sep = '')
                                                        
                                                        gene_set.list <- c(gene_set.list, serial_genes)
                                                    }
                                                }
                                                
                                                # *** Select top list ***
                                                top_num <- c(1, 3, 5, 10, 20, 30, 50)
                                                top_num.idx <- which(top_num <= top_gene_num)
                                                if (length(top_num.idx) > 0) {
                                                    top_num <- top_num[top_num.idx]
                                                } else {
                                                    top_num <- c(1)
                                                }
                                                top_num.idx <- which(top_num <= length(top_specific_marker_names))
                                                if (length(top_num.idx) > 0) {
                                                    top_num <- top_num[top_num.idx]
                                                } else {
                                                    top_num <- c()
                                                }
                                                if (length(top_num) > 0) {
                                                    for (top in 1:length(top_num)) {
                                                        top_genes <- list(unique(head(top_specific_marker_names, top_num[top])))
                                                        if (top_num[top] == 1) {
                                                            names(top_genes) = paste('top',top_num[top],'_','gene', sep = '')
                                                        } else if (top_num[top] > 1) {
                                                            names(top_genes) = paste('top',top_num[top],'_','genes', sep = '')
                                                        }
                                                        
                                                        gene_set.list <- c(gene_set.list, top_genes)
                                                    }
                                                }
                                                
                                                # *** query gene(s) ***
                                                if (! is.null(genes_of_interest.input)) {
                                                    for (gX in 1:length(genes_of_interest.input)) {
                                                        for (gY in 1:length(genes_of_interest.input[[gX]])) {
                                                            genes_of_interest.gX.gY <- list(genes_of_interest.input[[gX]][gY])
                                                            names(genes_of_interest.gX.gY) = paste(names(genes_of_interest.input[gX]),'.',genes_of_interest.input[[gX]][gY], sep = '')
                                                            gene_set.list <- c(gene_set.list, genes_of_interest.gX.gY)
                                                        }
                                                    }
                                                    genes_of_interest <- genes_of_interest.input
                                                    gene_set.list <- c(gene_set.list, genes_of_interest)
                                                }
                                                
                                                # *** signature gene(s) ***
                                                if (! is.null(signature_genes.input)) {
                                                    for (gX in 1:length(signature_genes.input)) {
                                                        for (gY in 1:length(signature_genes.input[[gX]])) {
                                                            signature_genes.gX.gY <- list(signature_genes.input[[gX]][gY])
                                                            names(signature_genes.gX.gY) = paste(names(signature_genes.input[gX]),'.',signature_genes.input[[gX]][gY], sep = '')
                                                            gene_set.list <- c(gene_set.list, signature_genes.gX.gY)
                                                        }
                                                    }
                                                    signature_genes <- signature_genes.input
                                                    gene_set.list <- c(gene_set.list, signature_genes)
                                                }
                                                
                                                if (length(gene_set.list) > 0) {
                                                    for (N in 1:length(gene_set.list)) {
                                                        #print(names(gene_set.list[N]))
                                                        
                                                        marker.genes <- gene_set.list[[N]]
                                                        marker.genes.idx <- which(sapply(marker.genes, FUN = function(X) X %in% rowData(cds)[,"gene_short_name"]))
                                                        if (length(marker.genes.idx) > 0) {
                                                            marker.genes = marker.genes[marker.genes.idx]
                                                            
                                                            rowidx <- sapply(marker.genes, FUN = function(X) which(rowData(cds)[,"gene_short_name"] %in% X))
                                                            rowidx <- unlist(rowidx)
                                                            
                                                            cds_subset <- cds[rowidx, , drop = FALSE]
                                                            
                                                            if (nrow(cds_subset) > 100) { # Error: cds_subset has more than 100 genes - pass only the subset of the CDS to be plotted.
                                                                cds_subset <- head(cds_subset, 100)
                                                            }
                                                            
                                                            #rowData(cds_subset)[,"gene_short_name"] = factor(rowData(cds_subset)[,"gene_short_name"], levels = unique(rowData(cds_subset)[,"gene_short_name"]))
                                                            
                                                            if (group_var == "cluster") {
                                                                cds_subset.level = levels(cds_subset@clusters[[redDim.method[redDim]]]$clusters)
                                                            } else if (group_var == "partition") {
                                                                cds_subset.level = levels(cds_subset@clusters[[redDim.method[redDim]]]$partitions)
                                                            } else if (group_var == batch) {
                                                                #cds_subset.level = levels(colData(cds_subset)[, group_var])
                                                                cds_subset.level = names(which(sapply(levels(colData(cds_subset)[, group_var]), FUN = function(X) X %in% unique(colData(cds_subset)[, group_var]))))
                                                            }
                                                            
                                                            if (length(cds_subset.level) > 20) {
                                                                ncol.num = 1
                                                                width.size = length(cds_subset.level) / 4 + 0.75
                                                                height.size = ceiling(nrow(cds_subset) / ncol.num) * 2 + 0.75
                                                            } else {
                                                                ncol.num = 2
                                                                width.size = length(cds_subset.level) / 1.5 + 0.75
                                                                height.size = ceiling(nrow(cds_subset) / ncol.num) * 1.25 + 0.75
                                                            }
                                                            
                                                            if (max(nchar(cds_subset.level)) > 20) {
                                                                angle.num = 60
                                                            } else {
                                                                angle.num = 45
                                                            }
                                                            
                                                            dim.folder = paste(eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))),'D', sep = '')
                                                            
                                                            mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim], sep = '')
                                                            subDir <- dim.folder
                                                            if (! file.exists(file.path(mainDir, subDir))) {
                                                                dir.create(file.path(mainDir, subDir))
                                                            }
                                                            
                                                            mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder, sep = '')
                                                            subDir <- 'differential'
                                                            if (! file.exists(file.path(mainDir, subDir))) {
                                                                dir.create(file.path(mainDir, subDir))
                                                            }
                                                            mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/','differential', sep = '')
                                                            subDir <- 'gene'
                                                            if (! file.exists(file.path(mainDir, subDir))) {
                                                                dir.create(file.path(mainDir, subDir))
                                                            }
                                                            
                                                            if (length(grep(paste(genes_of_interest.titles, collapse = '|'), names(gene_set.list[N]), value = TRUE, invert = FALSE)) > 0) {
                                                                mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/','differential', sep = '')
                                                                subDir <- 'genes_of_interest'
                                                                if (! file.exists(file.path(mainDir, subDir))) {
                                                                    dir.create(file.path(mainDir, subDir))
                                                                }
                                                                mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/','differential','/','genes_of_interest', sep = '')
                                                                subDir <- unlist(strsplit(names(gene_set.list[N]), '\\.'))[1]
                                                                if (! file.exists(file.path(mainDir, subDir))) {
                                                                    dir.create(file.path(mainDir, subDir))
                                                                }
                                                                mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/','differential','/','genes_of_interest','/',unlist(strsplit(names(gene_set.list[N]), '\\.'))[1], sep = '')
                                                                subDir <- 'gene'
                                                                if (! file.exists(file.path(mainDir, subDir))) {
                                                                    dir.create(file.path(mainDir, subDir))
                                                                }
                                                                
                                                                file.folder = paste('genes_of_interest','/',unlist(strsplit(names(gene_set.list[N]), '\\.'))[1],'/','gene', sep = '')
                                                                if (names(gene_set.list[N]) %in% genes_of_interest.titles) {
                                                                    file.name.suffix = names(gene_set.list[N])
                                                                } else {
                                                                    file.name.suffix = gsub("\\.","_",unlist(strsplit(names(gene_set.list[N]), '\\.'))[2])
                                                                }
                                                            } else if (length(grep(paste(signature_genes.titles, collapse = '|'), names(gene_set.list[N]), value = TRUE, invert = FALSE)) > 0) {
                                                                mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/','differential', sep = '')
                                                                subDir <- 'signature_genes'
                                                                if (! file.exists(file.path(mainDir, subDir))) {
                                                                    dir.create(file.path(mainDir, subDir))
                                                                }
                                                                mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/','differential','/','signature_genes', sep = '')
                                                                subDir <- unlist(strsplit(names(gene_set.list[N]), '\\.'))[1]
                                                                if (! file.exists(file.path(mainDir, subDir))) {
                                                                    dir.create(file.path(mainDir, subDir))
                                                                }
                                                                mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/','differential','/','signature_genes','/',unlist(strsplit(names(gene_set.list[N]), '\\.'))[1], sep = '')
                                                                subDir <- 'gene'
                                                                if (! file.exists(file.path(mainDir, subDir))) {
                                                                    dir.create(file.path(mainDir, subDir))
                                                                }
                                                                
                                                                file.folder = paste('signature_genes','/',unlist(strsplit(names(gene_set.list[N]), '\\.'))[1],'/','gene', sep = '')
                                                                if (names(gene_set.list[N]) %in% signature_genes.titles) {
                                                                    file.name.suffix = names(gene_set.list[N])
                                                                } else {
                                                                    file.name.suffix = gsub("\\.","_",unlist(strsplit(names(gene_set.list[N]), '\\.'))[2])
                                                                }
                                                            } else {
                                                                file.folder = 'gene'
                                                                file.name.suffix = gsub("\\.","_",names(gene_set.list[N]))
                                                            }
                                                            
                                                            mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/','differential','/',file.folder, sep = '')
                                                            subDir <- gsub("\\.","_",group_var)
                                                            if (! file.exists(file.path(mainDir, subDir))) {
                                                                dir.create(file.path(mainDir, subDir))
                                                            }
                                                            mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/','differential','/',file.folder,'/',gsub("\\.","_",group_var), sep = '')
                                                            subDir <- 'violin'
                                                            if (! file.exists(file.path(mainDir, subDir))) {
                                                                dir.create(file.path(mainDir, subDir))
                                                            }
                                                            
                                                            if (is.null(file.name.suffix)) {
                                                                png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/','differential','/',file.folder,'/',gsub("\\.","_",group_var),'/','violin','/',project,'.',gsub("\\.","_",subdata.folder),'.',gsub("\\.","_",project.subname),'.',redDim.method[redDim],'.','clustering','.','differential','.','plot_genes_violin','.',gsub("\\.","_",group_var),'.png', sep = ''), width = width.size, height = height.size, units = 'in', res = 600)
                                                            } else {
                                                                if (is.na(file.name.suffix) | file.name.suffix == "") {
                                                                    png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/','differential','/',file.folder,'/',gsub("\\.","_",group_var),'/','violin','/',project,'.',gsub("\\.","_",subdata.folder),'.',gsub("\\.","_",project.subname),'.',redDim.method[redDim],'.','clustering','.','differential','.','plot_genes_violin','.',gsub("\\.","_",group_var),'.png', sep = ''), width = width.size, height = height.size, units = 'in', res = 600)
                                                                } else {
                                                                    png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/','differential','/',file.folder,'/',gsub("\\.","_",group_var),'/','violin','/',project,'.',gsub("\\.","_",subdata.folder),'.',gsub("\\.","_",project.subname),'.',redDim.method[redDim],'.','clustering','.','differential','.','plot_genes_violin','.',file.name.suffix,'.',gsub("\\.","_",group_var),'.png', sep = ''), width = width.size, height = height.size, units = 'in', res = 600)
                                                                }
                                                            }
                                                            
                                                            p <- plot_genes_violin(cds_subset, group_cells_by = group_var, ncol = ncol.num) + theme(axis.text.x = element_text(angle=angle.num, hjust=1))
                                                            suppressWarnings(print(p))
                                                            dev.off()
                                                            
                                                            # /// plot_genes_by_group ///
                                                            if (length(marker.genes) > 1) {
                                                                
                                                                mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/','differential','/',file.folder,'/',gsub("\\.","_",group_var), sep = '')
                                                                subDir <- 'group'
                                                                if (! file.exists(file.path(mainDir, subDir))) {
                                                                    dir.create(file.path(mainDir, subDir))
                                                                }
                                                                
                                                                for (order.tp in 1:3) {
                                                                    if (order.tp == 1) {
                                                                        ordering.type = "cluster_row_col"
                                                                    } else if (order.tp == 2) {
                                                                        ordering.type = "maximal_on_diag"
                                                                    } else if (order.tp == 3) {
                                                                        ordering.type = "none"
                                                                    }
                                                                    
                                                                    for (xy.axis in 1:2) {
                                                                        if (xy.axis == 1) {
                                                                            axis_order.type = "group_marker"
                                                                            width.size = length(levels(colData(cds)[, group_var])) + 1
                                                                            height.size = length(marker.genes) / 1.5 + 2.5
                                                                        } else if (xy.axis == 2) {
                                                                            axis_order.type = "marker_group"
                                                                            width.size = length(marker.genes) / 1.5 + 1
                                                                            height.size = length(levels(colData(cds)[, group_var])) + 2.5
                                                                        }
                                                                        
                                                                        if (is.null(file.name.suffix)) {
                                                                            png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/','differential','/',file.folder,'/',gsub("\\.","_",group_var),'/','group','/',project,'.',gsub("\\.","_",subdata.folder),'.',gsub("\\.","_",project.subname),'.',redDim.method[redDim],'.','clustering','.','differential','.','plot_genes_by_group','.',axis_order.type,'.',ordering.type,'.',gsub("\\.","_",group_var),'.png', sep = ''), width = width.size, height = height.size, units = 'in', res = 600)
                                                                        } else {
                                                                            if (is.na(file.name.suffix) | file.name.suffix == "") {
                                                                                png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/','differential','/',file.folder,'/',gsub("\\.","_",group_var),'/','group','/',project,'.',gsub("\\.","_",subdata.folder),'.',gsub("\\.","_",project.subname),'.',redDim.method[redDim],'.','clustering','.','differential','.','plot_genes_by_group','.',axis_order.type,'.',ordering.type,'.',gsub("\\.","_",group_var),'.png', sep = ''), width = width.size, height = height.size, units = 'in', res = 600)
                                                                            } else {
                                                                                png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/','differential','/',file.folder,'/',gsub("\\.","_",group_var),'/','group','/',project,'.',gsub("\\.","_",subdata.folder),'.',gsub("\\.","_",project.subname),'.',redDim.method[redDim],'.','clustering','.','differential','.','plot_genes_by_group','.',axis_order.type,'.',ordering.type,'.',file.name.suffix,'.',gsub("\\.","_",group_var),'.png', sep = ''), width = width.size, height = height.size, units = 'in', res = 600)
                                                                            }
                                                                        }
                                                                        
                                                                        p <- plot_genes_by_group(cds, 
                                                                                                 reduction_method = redDim.method[redDim], 
                                                                                                 markers = marker.genes, 
                                                                                                 norm_method = "log", 
                                                                                                 group_cells_by = group_var, 
                                                                                                 ordering_type = ordering.type, 
                                                                                                 axis_order = axis_order.type, 
                                                                                                 max.size = 12)
                                                                        print(p)
                                                                        dev.off()
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    
                    
                    # --- subset ----------------------------------------------------------------------------------------------------------------------------------------
                    for (redDim in 1:length(redDim.method)) {
                        dim.folder = paste(eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))),'D', sep = '')
                        dim.finalX = eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) - (2-1)
                        
                        if (dim.finalX >= 1) {
                            
                            for (dim.comb in 1:dim.finalX) {
                                dim.x <- dim.comb
                                dim.y <- dim.x + 1
                                
                                if (dim.finalX > 1) {
                                    if (dim.x == 1) {
                                        comb.run = 2
                                    } else {
                                        comb.run = 1
                                    }
                                } else {
                                    comb.run = 1
                                }
                                
                                for (dim.comb.run in 1:comb.run) {
                                    if (dim.comb.run > 1) {
                                        dim.y <- dim.y + 1
                                    }
                                    
                                    if (eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) == 2) {
                                        dim.xy <- paste('D',dim.x,'-','D',dim.y, sep = '')
                                    } else {
                                        dim.xy <- paste('D',dim.x,'-','D',dim.y, sep = '')
                                    }
                                    #print(dim.xy)
                                    
                                    
                                    # /// subset ///
                                    for (SUB in SUB.run) {
                                        if (SUB == 0) {
                                            next # next halts the processing of the current iteration and advances the looping index.
                                        } else if (SUB == 1) {
                                            SUB.prefolder = 'cell_subset'
                                            SUB.prename = 'cell_subset'
                                        } else if (SUB == 2) {
                                            if (length(cell_subset) == 0 & cell_choose == TRUE) {
                                                SUB.prefolder = 'choose'
                                                SUB.prename = 'choose'
                                            } else if (length(cell_subset) > 0 & cell_choose == TRUE) {
                                                SUB.prefolder = paste('cell_subset','/','choose', sep = '')
                                                SUB.prename = paste('cell_subset','.','choose', sep = '')
                                            }
                                        }
                                        
                                        mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy, sep = '')
                                        subDir <- SUB.prefolder
                                        if (! file.exists(file.path(mainDir, subDir))) {
                                            dir.create(file.path(mainDir, subDir), recursive = TRUE)
                                        }
                                        
                                        if (dim.comb == 1 & dim.comb.run == 1) {
                                            if (SUB == 1) {
                                                cds_choose.cell <- cell_subset
                                                
                                                cds.colidx <- which(sapply(colnames(cds), FUN = function(X) X %in% cds_choose.cell))
                                                if (length(cds.colidx) > 0) {
                                                    cds_choose <- cds[, cds.colidx]
                                                } else {
                                                    cds_choose <- cds
                                                }
                                            } else if (SUB == 2) {
                                                if (length(cell_subset) == 0 & cell_choose == TRUE) {
                                                    cds4choose <- cds
                                                } else if (length(cell_subset) > 0 & cell_choose == TRUE) {
                                                    cds4choose <- cds_choose
                                                }
                                                
                                                # Choose cells interactively to subset a cds
                                                # clear_cds: Logical, clear CDS slots before returning. After clearing the cds, re-run processing from preprocess_cds(), ... Default is FALSE.
                                                # return_list: Logical, return a list of cells instead of a subsetted CDS object.
                                                #cds_choose.before <- choose_cells(cds4choose, reduction_method = redDim.method[redDim], clear_cds = TRUE, return_list = FALSE)
                                                #cds_choose.after <- choose_cells(cds4choose, reduction_method = redDim.method[redDim], clear_cds = FALSE, return_list = FALSE)
                                                cds_choose.cell <- choose_cells(cds4choose, reduction_method = redDim.method[redDim], clear_cds = FALSE, return_list = TRUE)
                                                
                                                if (length(cell_subset) > 0) {
                                                    cds_choose.cell <- intersect(cell_subset, cds_choose.cell)
                                                }
                                                
                                                cds4choose.colidx <- which(sapply(colnames(cds4choose), FUN = function(X) X %in% cds_choose.cell))
                                                if (length(cds4choose.colidx) > 0) {
                                                    cds_choose <- cds4choose[, cds4choose.colidx]
                                                } else {
                                                    cds_choose <- cds4choose
                                                }
                                            }
                                            
                                            mapping.idx = sapply(rownames(colData(cds_choose)), FUN = function(X) which(names(cds_choose@clusters[[redDim.method[redDim]]]$clusters) %in% X))
                                            cds_choose@clusters[[redDim.method[redDim]]]$clusters = cds_choose@clusters[[redDim.method[redDim]]]$clusters[mapping.idx]
                                            cds_choose@clusters[[redDim.method[redDim]]]$clusters = factor(cds_choose@clusters[[redDim.method[redDim]]]$clusters, levels = mixedsort(as.character(unique(cds_choose@clusters[[redDim.method[redDim]]]$clusters))))
                                            colData(cds_choose)$cluster = cds_choose@clusters[[redDim.method[redDim]]]$clusters
                                            
                                            mapping.idx = sapply(rownames(colData(cds_choose)), FUN = function(X) which(names(cds_choose@clusters[[redDim.method[redDim]]]$partitions) %in% X))
                                            cds_choose@clusters[[redDim.method[redDim]]]$partitions = cds_choose@clusters[[redDim.method[redDim]]]$partitions[mapping.idx]
                                            cds_choose@clusters[[redDim.method[redDim]]]$partitions = factor(cds_choose@clusters[[redDim.method[redDim]]]$partitions, levels = mixedsort(as.character(unique(cds_choose@clusters[[redDim.method[redDim]]]$partitions))))
                                            colData(cds_choose)$partition = cds_choose@clusters[[redDim.method[redDim]]]$partitions
                                            
                                            mapping.idx = sapply(rownames(colData(cds_choose)), FUN = function(X) which(rownames(cds_choose@principal_graph_aux[[redDim.method[redDim]]]$pr_graph_cell_proj_closest_vertex) %in% X))
                                            cds_choose@principal_graph_aux[[redDim.method[redDim]]]$pr_graph_cell_proj_closest_vertex = cds_choose@principal_graph_aux[[redDim.method[redDim]]]$pr_graph_cell_proj_closest_vertex[mapping.idx, , drop = FALSE]
                                            
                                            saveRDS(cds_choose, paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',project,'.',SUB.prename,'.','clustering','.','cds_choose','.rds', sep = ''))
                                            #cds_choose <- readRDS(paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',project,'.',SUB.prename,'.','clustering','.','cds_choose','.rds', sep = ''))
                                            
                                            if (SUB == 1) {
                                                cds_choose.cell_subset <- cds_choose
                                            } else if (SUB == 2) {
                                                cds_choose.cell_choose <- cds_choose
                                            }
                                        } else {
                                            if (SUB == 1) {
                                                cds_choose <- cds_choose.cell_subset
                                            } else if (SUB == 2) {
                                                cds_choose <- cds_choose.cell_choose
                                            }
                                        }
                                        
                                        
                                        ## Step 2: Reduce the dimensions using UMAP *******************************************************
                                        # /// no color_cells ///
                                        png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',SUB.prefolder,'/',project,'.',SUB.prename,'.','clustering','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,'.png', sep = ''), width = 4.075, height = 4, units = 'in', res = 600)
                                        p <- plot_cells(cds_choose, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d)
                                        print(p)
                                        dev.off()
                                        
                                        # /// color_cells_by = cell_type ///
                                        if (length(cell_type) > 0) {
                                            for (color.plt in 1:length(cell_type)) {
                                                legend.max_nchar = max(nchar(levels(colData(cds_choose)[, cell_type[color.plt]])))
                                                legend.multiply_factor = 0.5 + legend.max_nchar * 0.05
                                                legend.num = length(levels(colData(cds_choose)[, cell_type[color.plt]]))
                                                legend.ncol = ceiling(legend.num / 10)
                                                width.size = 4.075 + (legend.multiply_factor * legend.ncol)
                                                
                                                png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',SUB.prefolder,'/',project,'.',SUB.prename,'.','clustering','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,'.',gsub("\\.","_",cell_type[color.plt]),'.png', sep = ''), width = 4.075, height = 4, units = 'in', res = 600)
                                                p <- plot_cells(cds_choose, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = cell_type[color.plt], label_groups_by_cluster = FALSE) # Finally, you can disable this labeling policy, making plot_cells() behave like it did before we called cluster_cells()
                                                if (! is.null(color_palette)) {
                                                    cell_type.colidx <- which(paste(unique(c(cell_type[color.plt], unlist(strsplit(cell_type[color.plt], '\\.|_'))[1])),'.','color', sep = '') %in% names(color_palette))
                                                    if (length(cell_type.colidx) == 1) {
                                                        p <- p + scale_color_manual(values = color_palette[[which(names(color_palette) %in% paste(unique(c(cell_type[color.plt], unlist(strsplit(cell_type[color.plt], '\\.|_'))[1])),'.','color', sep = ''))]])
                                                    }
                                                }
                                                print(p)
                                                dev.off()
                                                
                                                png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',SUB.prefolder,'/',project,'.',SUB.prename,'.','clustering','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,'.',gsub("\\.","_",cell_type[color.plt]),'_','legend','.png', sep = ''), width = width.size, height = 4, units = 'in', res = 600)
                                                p <- plot_cells(cds_choose, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = cell_type[color.plt], label_groups_by_cluster = FALSE, label_cell_groups = FALSE) # Finally, you can disable this labeling policy, making plot_cells() behave like it did before we called cluster_cells()
                                                if (! is.null(color_palette)) {
                                                    cell_type.colidx <- which(paste(unique(c(cell_type[color.plt], unlist(strsplit(cell_type[color.plt], '\\.|_'))[1])),'.','color', sep = '') %in% names(color_palette))
                                                    if (length(cell_type.colidx) == 1) {
                                                        p <- p + scale_color_manual(values = color_palette[[which(names(color_palette) %in% paste(unique(c(cell_type[color.plt], unlist(strsplit(cell_type[color.plt], '\\.|_'))[1])),'.','color', sep = ''))]])
                                                    }
                                                }
                                                p <- p + theme(legend.position = "right")
                                                p$guides$colour$ncol = legend.ncol
                                                if (project.name %in% c("GSE139944", c("GSM4150376", "GSM4150377", "GSM4150378", "GSM4150379"))) {
                                                    if (cell_type[color.plt] == "pathway_level_1") {
                                                        p$guides$colour$title = "Pathway"
                                                    } else if (cell_type[color.plt] == "product_name") {
                                                        p$guides$colour$title = "Treatment"
                                                    } else if (cell_type[color.plt] == "dose_character") {
                                                        p$guides$colour$title = dose_character.unit
                                                    } else {
                                                        p$guides$colour$title = "Treatment"
                                                    }
                                                }
                                                print(p)
                                                dev.off()
                                            }
                                        }
                                        
                                        
                                        ## Step 3: Remove batch effects with cell alignment ***********************************************
                                        if (! is.null(batch)) {
                                            if (project.name == "tutorial") {
                                                width.size = 4.9
                                                height.size = 4
                                            } else {
                                                legend.max_nchar = max(nchar(levels(colData(cds_choose)[, batch])))
                                                legend.multiply_factor = 0.5 + legend.max_nchar * 0.05
                                                legend.num = length(levels(colData(cds_choose)[, batch]))
                                                legend.ncol = ceiling(legend.num / 10)
                                                width.size = 4.075 + (legend.multiply_factor * legend.ncol)
                                                height.size = 4
                                            }
                                            
                                            png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',SUB.prefolder,'/',project,'.',SUB.prename,'.','clustering','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'png', sep = ''), width = width.size, height = height.size, units = 'in', res = 600)
                                            p <- plot_cells(cds_choose, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = batch, label_cell_groups = FALSE)
                                            p <- p + theme(legend.position = "right")
                                            p$guides$colour$ncol = legend.ncol
                                            print(p)
                                            dev.off()
                                        }
                                        
                                        
                                        ## Step 4: Cluster the cells **********************************************************************
                                        legend.max_nchar = max(nchar(levels(cds_choose@clusters[[redDim.method[redDim]]]$clusters)))
                                        legend.multiply_factor = 0.5 + legend.max_nchar * 0.05
                                        legend.num = length(levels(cds_choose@clusters[[redDim.method[redDim]]]$clusters))
                                        legend.ncol = ceiling(legend.num / 10)
                                        width.size = 4.075 + (legend.multiply_factor * legend.ncol)
                                        
                                        png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',SUB.prefolder,'/',project,'.',SUB.prename,'.','clustering','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'cluster','.png', sep = ''), width = 4.075, height = 4, units = 'in', res = 600)
                                        p <- plot_cells(cds_choose, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d) # call plot_cells() with no arguments, it colors the cells by cluster according to default.
                                        print(p)
                                        dev.off()
                                        
                                        png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',SUB.prefolder,'/',project,'.',SUB.prename,'.','clustering','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'cluster','_','legend','.png', sep = ''), width = width.size, height = 4, units = 'in', res = 600)
                                        p <- plot_cells(cds_choose, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d) # call plot_cells() with no arguments, it colors the cells by cluster according to default.
                                        p <- p + theme(legend.position = "right")
                                        p$guides$colour$ncol = legend.ncol
                                        print(p)
                                        dev.off()
                                        
                                        if (length(cell_type) > 0) {
                                            for (color.plt in 1:length(cell_type)) {
                                                legend.max_nchar = max(nchar(levels(colData(cds_choose)[, cell_type[color.plt]])))
                                                legend.multiply_factor = 0.5 + legend.max_nchar * 0.05
                                                legend.num = length(levels(colData(cds_choose)[, cell_type[color.plt]]))
                                                legend.ncol = ceiling(legend.num / 10)
                                                width.size = 4.075 + (legend.multiply_factor * legend.ncol)
                                                
                                                png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',SUB.prefolder,'/',project,'.',SUB.prename,'.','clustering','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'cluster','.',gsub("\\.","_",cell_type[color.plt]),'.png', sep = ''), width = 4.075, height = 4, units = 'in', res = 600)
                                                p <- plot_cells(cds_choose, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = cell_type[color.plt]) # each cluster is labeled according the most common annotation within it
                                                if (! is.null(color_palette)) {
                                                    cell_type.colidx <- which(paste(unique(c(cell_type[color.plt], unlist(strsplit(cell_type[color.plt], '\\.|_'))[1])),'.','color', sep = '') %in% names(color_palette))
                                                    if (length(cell_type.colidx) == 1) {
                                                        p <- p + scale_color_manual(values = color_palette[[which(names(color_palette) %in% paste(unique(c(cell_type[color.plt], unlist(strsplit(cell_type[color.plt], '\\.|_'))[1])),'.','color', sep = ''))]])
                                                    }
                                                }
                                                print(p)
                                                dev.off()
                                                
                                                png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',SUB.prefolder,'/',project,'.',SUB.prename,'.','clustering','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'cluster','.',gsub("\\.","_",cell_type[color.plt]),'_','legend','.png', sep = ''), width = width.size, height = 4, units = 'in', res = 600)
                                                p <- plot_cells(cds_choose, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = cell_type[color.plt], label_cell_groups = FALSE) # each cluster is labeled according the most common annotation within it
                                                if (! is.null(color_palette)) {
                                                    cell_type.colidx <- which(paste(unique(c(cell_type[color.plt], unlist(strsplit(cell_type[color.plt], '\\.|_'))[1])),'.','color', sep = '') %in% names(color_palette))
                                                    if (length(cell_type.colidx) == 1) {
                                                        p <- p + scale_color_manual(values = color_palette[[which(names(color_palette) %in% paste(unique(c(cell_type[color.plt], unlist(strsplit(cell_type[color.plt], '\\.|_'))[1])),'.','color', sep = ''))]])
                                                    }
                                                }
                                                p <- p + theme(legend.position = "right")
                                                p$guides$colour$ncol = legend.ncol
                                                if (project.name %in% c("GSE139944", c("GSM4150376", "GSM4150377", "GSM4150378", "GSM4150379"))) {
                                                    if (cell_type[color.plt] == "pathway_level_1") {
                                                        p$guides$colour$title = "Pathway"
                                                    } else if (cell_type[color.plt] == "product_name") {
                                                        p$guides$colour$title = "Treatment"
                                                    } else if (cell_type[color.plt] == "dose_character") {
                                                        p$guides$colour$title = dose_character.unit
                                                    } else {
                                                        p$guides$colour$title = "Treatment"
                                                    }
                                                }
                                                print(p)
                                                dev.off()
                                                
                                                png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',SUB.prefolder,'/',project,'.',SUB.prename,'.','clustering','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'cluster','.',gsub("\\.","_",cell_type[color.plt]),'.','disableLabel','.png', sep = ''), width = 4.075, height = 4, units = 'in', res = 600)
                                                p <- plot_cells(cds_choose, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = cell_type[color.plt], label_groups_by_cluster = FALSE) # Finally, you can disable this labeling policy, making plot_cells() behave like it did before we called cluster_cells()
                                                if (! is.null(color_palette)) {
                                                    cell_type.colidx <- which(paste(unique(c(cell_type[color.plt], unlist(strsplit(cell_type[color.plt], '\\.|_'))[1])),'.','color', sep = '') %in% names(color_palette))
                                                    if (length(cell_type.colidx) == 1) {
                                                        p <- p + scale_color_manual(values = color_palette[[which(names(color_palette) %in% paste(unique(c(cell_type[color.plt], unlist(strsplit(cell_type[color.plt], '\\.|_'))[1])),'.','color', sep = ''))]])
                                                    }
                                                }
                                                print(p)
                                                dev.off()
                                                
                                                png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',SUB.prefolder,'/',project,'.',SUB.prename,'.','clustering','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'cluster','.',gsub("\\.","_",cell_type[color.plt]),'_','legend','.','disableLabel','.png', sep = ''), width = width.size, height = 4, units = 'in', res = 600)
                                                p <- plot_cells(cds_choose, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = cell_type[color.plt], label_groups_by_cluster = FALSE, label_cell_groups = FALSE) # Finally, you can disable this labeling policy, making plot_cells() behave like it did before we called cluster_cells()
                                                if (! is.null(color_palette)) {
                                                    cell_type.colidx <- which(paste(unique(c(cell_type[color.plt], unlist(strsplit(cell_type[color.plt], '\\.|_'))[1])),'.','color', sep = '') %in% names(color_palette))
                                                    if (length(cell_type.colidx) == 1) {
                                                        p <- p + scale_color_manual(values = color_palette[[which(names(color_palette) %in% paste(unique(c(cell_type[color.plt], unlist(strsplit(cell_type[color.plt], '\\.|_'))[1])),'.','color', sep = ''))]])
                                                    }
                                                }
                                                p <- p + theme(legend.position = "right")
                                                p$guides$colour$ncol = legend.ncol
                                                if (project.name %in% c("GSE139944", c("GSM4150376", "GSM4150377", "GSM4150378", "GSM4150379"))) {
                                                    if (cell_type[color.plt] == "pathway_level_1") {
                                                        p$guides$colour$title = "Pathway"
                                                    } else if (cell_type[color.plt] == "product_name") {
                                                        p$guides$colour$title = "Treatment"
                                                    } else if (cell_type[color.plt] == "dose_character") {
                                                        p$guides$colour$title = dose_character.unit
                                                    } else {
                                                        p$guides$colour$title = "Treatment"
                                                    }
                                                }
                                                print(p)
                                                dev.off()
                                            }
                                        }
                                        
                                        # The cluster_cells() also divides the cells into larger, more well separated groups called partitions, 
                                        # using a statistical test from Alex Wolf et al, introduced as part of their PAGA algorithm. (https://github.com/theislab/paga)
                                        png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',SUB.prefolder,'/',project,'.',SUB.prename,'.','clustering','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'partition','.png', sep = ''), width = 4.075, height = 4, units = 'in', res = 600)
                                        p <- plot_cells(cds_choose, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = "partition")
                                        print(p)
                                        dev.off()
                                        
                                        png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',SUB.prefolder,'/',project,'.',SUB.prename,'.','clustering','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'partition','.','group_cells','.png', sep = ''), width = 4.075, height = 4, units = 'in', res = 600)
                                        p <- plot_cells(cds_choose, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = "partition", group_cells_by = "partition")
                                        print(p)
                                        dev.off()
                                    }
                                }
                            }
                        }
                    }
                    
                    
                    # --- Graph-autocorrelation analysis for comparing clusters ---------------------------------------
                    # With graph autocorrelation:
                    
                    if (project.name == "tutorial") {
                        # In the L2 worm data, we identified a number of clusters that were very distinct as neurons:
                        colData(cds)$assigned_cell_type <- as.character(partitions(cds, reduction_method = "UMAP")) # (can also use clusters(cds) depending on your dataset)
                        
                        # Subset just the neurons:
                        # There are many subtypes of neurons, so perhaps the different neuron clusters correspond to different subtypes.
                        colData(cds)$assigned_cell_type <- dplyr::recode(colData(cds)$assigned_cell_type,
                                                                         "1"="Germline",
                                                                         "2"="Body wall muscle",
                                                                         "3"="Unclassified neurons",
                                                                         "4"="Vulval precursors",
                                                                         "5"="Failed QC",
                                                                         "6"="Seam cells",
                                                                         "7"="Pharyngeal epithelia",
                                                                         "8"="Coelomocytes",
                                                                         "9"="Am/PH sheath cells",
                                                                         "10"="Failed QC",
                                                                         "11"="Touch receptor neurons",
                                                                         "12"="Intestinal/rectal muscle",
                                                                         "13"="Pharyngeal neurons",
                                                                         "14"="NA",
                                                                         "15"="flp-1(+) interneurons",
                                                                         "16"="Canal associated neurons",
                                                                         "17"="Ciliated sensory neurons",
                                                                         "18"="Other interneurons",
                                                                         "19"="Pharyngeal gland",
                                                                         "20"="Failed QC",
                                                                         "21"="Ciliated sensory neurons",
                                                                         "22"="Oxygen sensory neurons",
                                                                         "23"="Ciliated sensory neurons",
                                                                         "24"="Ciliated sensory neurons",
                                                                         "25"="Ciliated sensory neurons",
                                                                         "26"="Ciliated sensory neurons",
                                                                         "27"="Oxygen sensory neurons",
                                                                         "28"="Ciliated sensory neurons",
                                                                         "29"="Unclassified neurons",
                                                                         "30"="Socket cells",
                                                                         "31"="Failed QC",
                                                                         "32"="Pharyngeal gland",
                                                                         "33"="Ciliated sensory neurons",
                                                                         "34"="Ciliated sensory neurons",
                                                                         "35"="Ciliated sensory neurons",
                                                                         "36"="Failed QC",
                                                                         "37"="Ciliated sensory neurons",
                                                                         "38"="Pharyngeal muscle")
                        subpopulation.type = c("neuron", "pharyngeal")
                    } else {
                        colData(cds)$assigned_cell_type <- NULL
                        subpopulation.type = c()
                    }
                    
                    
                    for (subpopulation.idx in 0:length(subpopulation.type)) {
                        
                        if (subpopulation.idx == 0) {
                            cds_subpopulation <- cds
                        } else {
                            cds_subpopulation <- cds[, grepl(subpopulation.type[subpopulation.idx], colData(cds)$assigned_cell_type, ignore.case=TRUE)]
                        }
                        
                        cds_subpopulation <- cds_subpopulation
                        cds_subpopulation.origin.processed = cds_subpopulation
                        assign(paste('cds_subpopulation','.',subpopulation.idx,'.origin.processed', sep = ''), cds_subpopulation.origin.processed)
                        
                        
                        if (ncol(cds_subpopulation) > 0) {
                            
                            if (length(subpopulation.type[subpopulation.idx]) > 0) {
                                subpopulation.folder <- subpopulation.type[subpopulation.idx]
                            } else {
                                subpopulation.folder <- project.subname
                            }
                            
                            mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder, sep = '')
                            subDir <- subpopulation.folder
                            if (! file.exists(file.path(mainDir, subDir))) {
                                dir.create(file.path(mainDir, subDir))
                            }
                            
                            for (redDim in 1:length(redDim.method)) {
                                
                                dim.folder = paste(eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))),'D', sep = '')
                                
                                mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder, sep = '')
                                subDir <- redDim.method[redDim]
                                if (! file.exists(file.path(mainDir, subDir))) {
                                    dir.create(file.path(mainDir, subDir))
                                }
                                
                                mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim], sep = '')
                                subDir <- dim.folder
                                if (! file.exists(file.path(mainDir, subDir))) {
                                    dir.create(file.path(mainDir, subDir))
                                }
                                
                                
                                # /// subset ///
                                for (SUB in SUB.run) {
                                    if (SUB == 0) {
                                        next # next halts the processing of the current iteration and advances the looping index.
                                    } else if (SUB == 1) {
                                        SUB.prefolder = 'cell_subset'
                                        SUB.prename = 'cell_subset'
                                    } else if (SUB == 2) {
                                        if (length(cell_subset) == 0 & cell_choose == TRUE) {
                                            SUB.prefolder = 'choose'
                                            SUB.prename = 'choose'
                                        } else if (length(cell_subset) > 0 & cell_choose == TRUE) {
                                            SUB.prefolder = paste('cell_subset','/','choose', sep = '')
                                            SUB.prename = paste('cell_subset','.','choose', sep = '')
                                        }
                                    }
                                    
                                    if (SUB == 1) {
                                        cds_subpopulation_choose.cell <- cell_subset
                                        
                                        cds_subpopulation.colidx <- which(sapply(colnames(cds_subpopulation), FUN = function(X) X %in% cds_subpopulation_choose.cell))
                                        if (length(cds_subpopulation.colidx) > 0) {
                                            cds_subpopulation_choose <- cds_subpopulation[, cds_subpopulation.colidx]
                                        } else {
                                            cds_subpopulation_choose <- cds_subpopulation
                                        }
                                    } else if (SUB == 2) {
                                        if (length(cell_subset) == 0 & cell_choose == TRUE) {
                                            cds_subpopulation4choose <- cds_subpopulation
                                        } else if (length(cell_subset) > 0 & cell_choose == TRUE) {
                                            cds_subpopulation4choose <- cds_subpopulation_choose
                                        }
                                        
                                        #cds_subpopulation_choose.before <- choose_cells(cds_subpopulation4choose, reduction_method = redDim.method[redDim], clear_cds = TRUE, return_list = FALSE)
                                        #cds_subpopulation_choose.after <- choose_cells(cds_subpopulation4choose, reduction_method = redDim.method[redDim], clear_cds = FALSE, return_list = FALSE)
                                        cds_subpopulation_choose.cell <- choose_cells(cds_subpopulation4choose, reduction_method = redDim.method[redDim], clear_cds = FALSE, return_list = TRUE)
                                        
                                        if (length(cell_subset) > 0) {
                                            cds_subpopulation_choose.cell <- intersect(cell_subset, cds_subpopulation_choose.cell)
                                        }
                                        
                                        cds_subpopulation4choose.colidx <- which(sapply(colnames(cds_subpopulation4choose), FUN = function(X) X %in% cds_subpopulation_choose.cell))
                                        if (length(cds_subpopulation4choose.colidx) > 0) {
                                            cds_subpopulation_choose <- cds_subpopulation4choose[, cds_subpopulation4choose.colidx]
                                        } else {
                                            cds_subpopulation_choose <- cds_subpopulation4choose
                                        }
                                    }
                                    
                                    mapping.idx = sapply(rownames(colData(cds_subpopulation_choose)), FUN = function(X) which(names(cds_subpopulation_choose@clusters[[redDim.method[redDim]]]$clusters) %in% X))
                                    cds_subpopulation_choose@clusters[[redDim.method[redDim]]]$clusters = cds_subpopulation_choose@clusters[[redDim.method[redDim]]]$clusters[mapping.idx]
                                    cds_subpopulation_choose@clusters[[redDim.method[redDim]]]$clusters = factor(cds_subpopulation_choose@clusters[[redDim.method[redDim]]]$clusters, levels = mixedsort(as.character(unique(cds_subpopulation_choose@clusters[[redDim.method[redDim]]]$clusters))))
                                    colData(cds_subpopulation_choose)$cluster = cds_subpopulation_choose@clusters[[redDim.method[redDim]]]$clusters
                                    
                                    mapping.idx = sapply(rownames(colData(cds_subpopulation_choose)), FUN = function(X) which(names(cds_subpopulation_choose@clusters[[redDim.method[redDim]]]$partitions) %in% X))
                                    cds_subpopulation_choose@clusters[[redDim.method[redDim]]]$partitions = cds_subpopulation_choose@clusters[[redDim.method[redDim]]]$partitions[mapping.idx]
                                    cds_subpopulation_choose@clusters[[redDim.method[redDim]]]$partitions = factor(cds_subpopulation_choose@clusters[[redDim.method[redDim]]]$partitions, levels = mixedsort(as.character(unique(cds_subpopulation_choose@clusters[[redDim.method[redDim]]]$partitions))))
                                    colData(cds_subpopulation_choose)$partition = cds_subpopulation_choose@clusters[[redDim.method[redDim]]]$partitions
                                    
                                    mapping.idx = sapply(rownames(colData(cds_subpopulation_choose)), FUN = function(X) which(rownames(cds_subpopulation_choose@principal_graph_aux[[redDim.method[redDim]]]$pr_graph_cell_proj_closest_vertex) %in% X))
                                    cds_subpopulation_choose@principal_graph_aux[[redDim.method[redDim]]]$pr_graph_cell_proj_closest_vertex = cds_subpopulation_choose@principal_graph_aux[[redDim.method[redDim]]]$pr_graph_cell_proj_closest_vertex[mapping.idx, , drop = FALSE]
                                    
                                    saveRDS(cds_subpopulation_choose, paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',project,'.',SUB.prename,'.','clustering','.','cds_subpopulation_choose','.rds', sep = ''))
                                    #cds_subpopulation_choose <- readRDS(paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',project,'.',SUB.prename,'.','clustering','.','cds_subpopulation_choose','.rds', sep = ''))
                                    
                                    if (SUB == 1) {
                                        cds_subpopulation_choose.cell_subset <- cds_subpopulation_choose
                                    } else if (SUB == 2) {
                                        cds_subpopulation_choose.cell_choose <- cds_subpopulation_choose
                                    }
                                    
                                    if (SUB == 1) {
                                        cds_subpopulation_choose.cell_subset <- cds_subpopulation_choose.cell_subset
                                        assign(paste('cds_subpopulation_choose','.',subpopulation.idx,'.cell_subset', sep = ''), cds_subpopulation_choose.cell_subset)
                                    } else if (SUB == 2) {
                                        cds_subpopulation_choose.cell_choose <- cds_subpopulation_choose.cell_choose
                                        assign(paste('cds_subpopulation_choose','.',subpopulation.idx,'.cell_choose', sep = ''), cds_subpopulation_choose.cell_choose)
                                    }
                                }
                                
                                
                                # /// 2d ///
                                dim.folder = paste(eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))),'D', sep = '')
                                dim.finalX = eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) - (2-1)
                                
                                if (dim.finalX >= 1) {
                                    
                                    for (dim.comb in 1:dim.finalX) {
                                        dim.x <- dim.comb
                                        dim.y <- dim.x + 1
                                        
                                        if (dim.finalX > 1) {
                                            if (dim.x == 1) {
                                                comb.run = 2
                                            } else {
                                                comb.run = 1
                                            }
                                        } else {
                                            comb.run = 1
                                        }
                                        
                                        for (dim.comb.run in 1:comb.run) {
                                            if (dim.comb.run > 1) {
                                                dim.y <- dim.y + 1
                                            }
                                            
                                            if (eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) == 2) {
                                                dim.xy <- paste('D',dim.x,'-','D',dim.y, sep = '')
                                            } else {
                                                dim.xy <- paste('D',dim.x,'-','D',dim.y, sep = '')
                                            }
                                            #print(dim.xy)
                                            
                                            mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder, sep = '')
                                            subDir <- dim.xy
                                            if (! file.exists(file.path(mainDir, subDir))) {
                                                dir.create(file.path(mainDir, subDir))
                                            }
                                            
                                            png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',project,'.',gsub("\\.","_",subdata.folder),'.',gsub("\\.","_",subpopulation.folder),'.','clustering','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,'.','partition','.','group_cells','.png', sep = ''), width = 4.075, height = 4, units = 'in', res = 600)
                                            p <- plot_cells(cds_subpopulation, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = "partition", group_cells_by = "partition")
                                            print(p)
                                            dev.off()
                                        }
                                    }
                                }
                                
                                
                                # /// 3d ///
                                # Working with 3D plots
                                if (eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) == 2) {
                                    dim.folder = paste(eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))),'D', sep = '')
                                    #dim.finalX = eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) - (3-1)
                                    #dim.folder = paste(3,'D', sep = '')
                                    dim.finalX = 3 - (3-1)
                                    
                                    mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim], sep = '')
                                    subDir <- dim.folder
                                    if (! file.exists(file.path(mainDir, subDir))) {
                                        dir.create(file.path(mainDir, subDir))
                                    }
                                    
                                    if (redDim == 1) {
                                        cds_3d <- new_cell_data_set(expression_matrix,
                                                                    cell_metadata = cell_metadata,
                                                                    gene_metadata = gene_metadata)
                                        
                                        if (ncol(cds_3d) / 2 <= 50) {
                                            # Warning in (function (A, nv = 5, nu = nv, maxit = 1000, work = nv + 7, reorth = TRUE,  :
                                            # You're computing too large a percentage of total singular values, use a standard svd instead.
                                            if (floor(ncol(cds_3d) / 2) == ncol(cds_3d) / 2) {
                                                cds_3d <- preprocess_cds(cds_3d, method = "PCA", num_dim = floor(ncol(cds_3d) / 2) - 1)
                                            } else {
                                                cds_3d <- preprocess_cds(cds_3d, method = "PCA", num_dim = floor(ncol(cds_3d) / 2))
                                            }
                                        } else {
                                            cds_3d <- preprocess_cds(cds_3d, method = "PCA", num_dim = 50) # method: Default is "PCA". method = c("PCA", "LSI") # "PCA": Principal Components Analysis (the standard for RNA-seq), "LSI": Latent Semantic Indexing (common in ATAC-seq)
                                        }
                                    }
                                    
                                    cds_3d <- reduce_dimension(cds_3d, reduction_method = redDim.method[redDim], max_components = 3)
                                    cds_3d <- cluster_cells(cds_3d, reduction_method = redDim.method[redDim])
                                    
                                    if (subpopulation.idx == 0) {
                                        cds_3d_subpopulation <- cds_3d
                                    } else {
                                        cds_3d_subpopulation <- cds_3d[, grepl(subpopulation.type[subpopulation.idx], colData(cds)$assigned_cell_type, ignore.case=TRUE)]
                                    }
                                    
                                    cds_subpopulation_3d <- cds_3d_subpopulation
                                    cds_subpopulation_3d.origin.processed = cds_subpopulation_3d
                                    assign(paste('cds_subpopulation_3d','.',subpopulation.idx,'.origin.processed', sep = ''), cds_subpopulation_3d.origin.processed)
                                    
                                    
                                    # /// subset ///
                                    for (SUB in SUB.run) {
                                        if (SUB == 0) {
                                            next # next halts the processing of the current iteration and advances the looping index.
                                        } else if (SUB == 1) {
                                            SUB.prefolder = 'cell_subset'
                                            SUB.prename = 'cell_subset'
                                        } else if (SUB == 2) {
                                            if (length(cell_subset) == 0 & cell_choose == TRUE) {
                                                SUB.prefolder = 'choose'
                                                SUB.prename = 'choose'
                                            } else if (length(cell_subset) > 0 & cell_choose == TRUE) {
                                                SUB.prefolder = paste('cell_subset','/','choose', sep = '')
                                                SUB.prename = paste('cell_subset','.','choose', sep = '')
                                            }
                                        }
                                        
                                        if (SUB == 1) {
                                            cds_3d_subpopulation_choose.cell <- cell_subset
                                            
                                            cds_3d_subpopulation.colidx <- which(sapply(colnames(cds_3d_subpopulation), FUN = function(X) X %in% cds_3d_subpopulation_choose.cell))
                                            if (length(cds_3d_subpopulation.colidx) > 0) {
                                                cds_3d_subpopulation_choose <- cds_3d_subpopulation[, cds_3d_subpopulation.colidx]
                                            } else {
                                                cds_3d_subpopulation_choose <- cds_3d_subpopulation
                                            }
                                        } else if (SUB == 2) {
                                            if (length(cell_subset) == 0 & cell_choose == TRUE) {
                                                cds_3d_subpopulation4choose <- cds_3d_subpopulation
                                            } else if (length(cell_subset) > 0 & cell_choose == TRUE) {
                                                cds_3d_subpopulation4choose <- cds_3d_subpopulation_choose
                                            }
                                            
                                            #cds_3d_subpopulation_choose.before <- choose_cells(cds_3d_subpopulation4choose, reduction_method = redDim.method[redDim], clear_cds = TRUE, return_list = FALSE)
                                            #cds_3d_subpopulation_choose.after <- choose_cells(cds_3d_subpopulation4choose, reduction_method = redDim.method[redDim], clear_cds = FALSE, return_list = FALSE)
                                            cds_3d_subpopulation_choose.cell <- choose_cells(cds_3d_subpopulation4choose, reduction_method = redDim.method[redDim], clear_cds = FALSE, return_list = TRUE)
                                            
                                            if (length(cell_subset) > 0) {
                                                cds_3d_subpopulation_choose.cell <- intersect(cell_subset, cds_3d_subpopulation_choose.cell)
                                            }
                                            
                                            cds_3d_subpopulation4choose.colidx <- which(sapply(colnames(cds_3d_subpopulation4choose), FUN = function(X) X %in% cds_3d_subpopulation_choose.cell))
                                            if (length(cds_3d_subpopulation4choose.colidx) > 0) {
                                                cds_3d_subpopulation_choose <- cds_3d_subpopulation4choose[, cds_3d_subpopulation4choose.colidx]
                                            } else {
                                                cds_3d_subpopulation_choose <- cds_3d_subpopulation4choose
                                            }
                                        }
                                        
                                        mapping.idx = sapply(rownames(colData(cds_3d_subpopulation_choose)), FUN = function(X) which(names(cds_3d_subpopulation_choose@clusters[[redDim.method[redDim]]]$clusters) %in% X))
                                        cds_3d_subpopulation_choose@clusters[[redDim.method[redDim]]]$clusters = cds_3d_subpopulation_choose@clusters[[redDim.method[redDim]]]$clusters[mapping.idx]
                                        cds_3d_subpopulation_choose@clusters[[redDim.method[redDim]]]$clusters = factor(cds_3d_subpopulation_choose@clusters[[redDim.method[redDim]]]$clusters, levels = mixedsort(as.character(unique(cds_3d_subpopulation_choose@clusters[[redDim.method[redDim]]]$clusters))))
                                        colData(cds_3d_subpopulation_choose)$cluster = cds_3d_subpopulation_choose@clusters[[redDim.method[redDim]]]$clusters
                                        
                                        mapping.idx = sapply(rownames(colData(cds_3d_subpopulation_choose)), FUN = function(X) which(names(cds_3d_subpopulation_choose@clusters[[redDim.method[redDim]]]$partitions) %in% X))
                                        cds_3d_subpopulation_choose@clusters[[redDim.method[redDim]]]$partitions = cds_3d_subpopulation_choose@clusters[[redDim.method[redDim]]]$partitions[mapping.idx]
                                        cds_3d_subpopulation_choose@clusters[[redDim.method[redDim]]]$partitions = factor(cds_3d_subpopulation_choose@clusters[[redDim.method[redDim]]]$partitions, levels = mixedsort(as.character(unique(cds_3d_subpopulation_choose@clusters[[redDim.method[redDim]]]$partitions))))
                                        colData(cds_3d_subpopulation_choose)$partition = cds_3d_subpopulation_choose@clusters[[redDim.method[redDim]]]$partitions
                                        
                                        mapping.idx = sapply(rownames(colData(cds_3d_subpopulation_choose)), FUN = function(X) which(rownames(cds_3d_subpopulation_choose@principal_graph_aux[[redDim.method[redDim]]]$pr_graph_cell_proj_closest_vertex) %in% X))
                                        cds_3d_subpopulation_choose@principal_graph_aux[[redDim.method[redDim]]]$pr_graph_cell_proj_closest_vertex = cds_3d_subpopulation_choose@principal_graph_aux[[redDim.method[redDim]]]$pr_graph_cell_proj_closest_vertex[mapping.idx, , drop = FALSE]
                                        
                                        saveRDS(cds_3d_subpopulation_choose, paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',project,'.',SUB.prename,'.','clustering','.','cds_3d_subpopulation_choose','.rds', sep = ''))
                                        #cds_3d_subpopulation_choose <- readRDS(paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',project,'.',SUB.prename,'.','clustering','.','cds_3d_subpopulation_choose','.rds', sep = ''))
                                        
                                        if (SUB == 1) {
                                            cds_3d_subpopulation_choose.cell_subset <- cds_3d_subpopulation_choose
                                        } else if (SUB == 2) {
                                            cds_3d_subpopulation_choose.cell_choose <- cds_3d_subpopulation_choose
                                        }
                                        
                                        if (SUB == 1) {
                                            cds_subpopulation_choose_3d.cell_subset <- cds_3d_subpopulation_choose.cell_subset
                                            assign(paste('cds_subpopulation_choose_3d','.',subpopulation.idx,'.cell_subset', sep = ''), cds_subpopulation_choose_3d.cell_subset)
                                        } else if (SUB == 2) {
                                            cds_subpopulation_choose_3d.cell_choose <- cds_3d_subpopulation_choose.cell_choose
                                            assign(paste('cds_subpopulation_choose_3d','.',subpopulation.idx,'.cell_choose', sep = ''), cds_subpopulation_choose_3d.cell_choose)
                                        }
                                    }
                                } else {
                                    dim.folder = paste(eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))),'D', sep = '')
                                    dim.finalX = eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) - (3-1)
                                    
                                    cds_subpopulation_3d <- cds_subpopulation
                                    cds_subpopulation_3d.origin.processed = cds_subpopulation_3d
                                    assign(paste('cds_subpopulation_3d','.',subpopulation.idx,'.origin.processed', sep = ''), cds_subpopulation_3d.origin.processed)
                                }
                                
                                
                                # /// subset ///
                                for (SUB in SUB.run) {
                                    if (SUB == 0) {
                                        #next # next halts the processing of the current iteration and advances the looping index.
                                    } else if (SUB == 1) {
                                        SUB.prefolder = 'cell_subset'
                                        SUB.prename = 'cell_subset'
                                    } else if (SUB == 2) {
                                        if (length(cell_subset) == 0 & cell_choose == TRUE) {
                                            SUB.prefolder = 'choose'
                                            SUB.prename = 'choose'
                                        } else if (length(cell_subset) > 0 & cell_choose == TRUE) {
                                            SUB.prefolder = paste('cell_subset','/','choose', sep = '')
                                            SUB.prename = paste('cell_subset','.','choose', sep = '')
                                        }
                                    }
                                    
                                    if (eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) == 2) {
                                        if (SUB == 1) {
                                            cds_subpopulation_choose_3d <- cds_subpopulation_choose_3d.cell_subset
                                        } else if (SUB == 2) {
                                            cds_subpopulation_choose_3d <- cds_subpopulation_choose_3d.cell_choose
                                        }
                                    } else {
                                        if (SUB == 1) {
                                            cds_subpopulation_choose_3d <- cds_subpopulation_choose.cell_subset
                                        } else if (SUB == 2) {
                                            cds_subpopulation_choose_3d <- cds_subpopulation_choose.cell_choose
                                        }
                                    }
                                    
                                    if (SUB == 0) {
                                        cds_subpopulation_3d = cds_subpopulation_3d.origin.processed
                                        file.prename = 'clustering'
                                    } else if (SUB == 1) {
                                        cds_subpopulation_3d = cds_subpopulation_choose_3d
                                        file.prename = paste(SUB.prename,'.','clustering', sep = '')
                                    } else if (SUB == 2) {
                                        cds_subpopulation_3d = cds_subpopulation_choose_3d
                                        file.prename = paste(SUB.prename,'.','clustering', sep = '')
                                    }
                                    
                                    if (dim.finalX >= 1) {
                                        
                                        for (dim.comb in 1:dim.finalX) {
                                            dim.x <- dim.comb
                                            dim.y <- dim.x + 1
                                            dim.z <- dim.y + 1
                                            
                                            if (eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) == 2) {
                                                dim.xyz <- paste('D',dim.x,'-','D',dim.y,'-','D',dim.z, sep = '')
                                            } else {
                                                dim.xyz <- paste('D',dim.x,'-','D',dim.y,'-','D',dim.z, sep = '')
                                            }
                                            #print(dim.xyz)
                                            
                                            mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder, sep = '')
                                            subDir <- dim.xyz
                                            if (! file.exists(file.path(mainDir, subDir))) {
                                                dir.create(file.path(mainDir, subDir))
                                            }
                                            
                                            if (SUB == 0) {
                                                file.prefolder = dim.xyz
                                            } else {
                                                mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xyz, sep = '')
                                                subDir <- SUB.prefolder
                                                if (! file.exists(file.path(mainDir, subDir))) {
                                                    dir.create(file.path(mainDir, subDir))
                                                }
                                                file.prefolder = paste(dim.xyz,'/',SUB.prefolder, sep = '')
                                            }
                                            
                                            plt.3d <- unique(c(cell_type.input, stage.input))
                                            
                                            if (redDim.method[redDim] == "UMAP") {
                                                if (! is.null(time)) {
                                                    plt.3d.num = length(plt.3d)
                                                } else { # No pseudotime calculated for reduction_method = UMAP. Please first run order_cells with reduction_method = UMAP.
                                                    plt.3d.num = length(plt.3d)
                                                }
                                            } else if (redDim.method[redDim] == "tSNE") {
                                                plt.3d.num = length(plt.3d)
                                            }
                                            
                                            if (plt.3d.num > 0) {
                                                for (color.plt in 1:plt.3d.num) {
                                                    if (color.plt <= length(plt.3d)) {
                                                        if (project.name %in% c("Pancreas", "DrugPairing", "CRISPRINT")) {
                                                            color_cells_by.attr = plt.3d[color.plt]
                                                        } else if (project.name %in% c("LeoPharma", "Eczema")) {
                                                            color_cells_by.attr = plt.3d[color.plt]
                                                        } else if (project.name %in% c("Paper")) {
                                                            color_cells_by.attr = plt.3d[color.plt]
                                                        } else if (project.name %in% c("GSE139944", c("GSM4150376", "GSM4150377", "GSM4150378", "GSM4150379"))) {
                                                            color_cells_by.attr = plt.3d[color.plt]
                                                        } else {
                                                            color_cells_by.attr = "partition"
                                                        }
                                                    } else if (color.plt > length(plt.3d)) {
                                                        #color_cells_by.attr = "pseudotime"
                                                    }
                                                    
                                                    #if (redDim.method[redDim] == "UMAP") {
                                                    #    trj.plt = 1
                                                    #} else if (redDim.method[redDim] == "tSNE") {
                                                    trj.plt = 0
                                                    #}
                                                    if (graph.run == 0) {
                                                        trj.plt = 0
                                                    }
                                                    
                                                    for (trj in 0:trj.plt) {
                                                        if (trj == 0) {
                                                            trj.prt.endidx = 1
                                                        } else if (trj == 1) {
                                                            trj.prt.endidx = 0
                                                        }
                                                        
                                                        #for (trj.prt in 1:1) {
                                                        for (trj.prt in 1:trj.prt.endidx) {
                                                            if (trj.prt == 1) {
                                                                trj.prt.name = 'with_partition'
                                                                cds_subpopulation_3d.plt = cds_subpopulation_3d
                                                            } else if (trj.prt == 0) {
                                                                #trj.prt.name = 'without_partition'
                                                                #cds_subpopulation_3d.plt = cds_subpopulation_3d.no_partition
                                                            }
                                                            
                                                            if (! is.null(cds_subpopulation_3d.plt)) {
                                                                
                                                                for (grid.plt in 1:0) {
                                                                    if (grid.plt == 1) {
                                                                        grid.name = 'with_grid'
                                                                    } else if (grid.plt == 0) {
                                                                        grid.name = 'without_grid'
                                                                    }
                                                                    
                                                                    if (trj == 0) {
                                                                        show_trajectory = FALSE
                                                                        file.name = paste(project,'.',file.prename,'.','plot_cells_3d','.',redDim.method[redDim],'_',dim.xyz,'.',gsub("\\.","_",color_cells_by.attr),'.','without_trajectory', sep = '')
                                                                        file.name = paste(file.name,'.',grid.name, sep = '')
                                                                    } else if (trj == 1) {
                                                                        show_trajectory = TRUE
                                                                        file.name = paste(project,'.',file.prename,'.','plot_cells_3d','.',redDim.method[redDim],'_',dim.xyz,'.',gsub("\\.","_",color_cells_by.attr),'.','with_trajectory', sep = '')
                                                                        #file.name = paste(file.name,'.',trj.prt.name, sep = '')
                                                                        file.name = paste(file.name,'.',grid.name, sep = '')
                                                                    }
                                                                    
                                                                    file.html = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/',file.name,'.html', sep = '')
                                                                    file.png = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/',file.name,'.png', sep = '')
                                                                    
                                                                    if (! is.null(color_palette)) {
                                                                        cell_type.colidx <- which(paste(unique(c(color_cells_by.attr, unlist(strsplit(color_cells_by.attr, '\\.|_'))[1])),'.','color', sep = '') %in% names(color_palette))
                                                                        if (length(cell_type.colidx) == 1) {
                                                                            cds_subpopulation_3d_plot_obj <- plot_cells_3d(cds_subpopulation_3d.plt, reduction_method = redDim.method[redDim], dims = c(dim.x, dim.y, dim.z), cell_size = cell_size_3d, show_trajectory_graph = show_trajectory, color_cells_by = color_cells_by.attr, color_palette = color_palette[[which(names(color_palette) %in% paste(unique(c(color_cells_by.attr, unlist(strsplit(color_cells_by.attr, '\\.|_'))[1])),'.','color', sep = ''))]])
                                                                        } else {
                                                                            cds_subpopulation_3d_plot_obj <- plot_cells_3d(cds_subpopulation_3d.plt, reduction_method = redDim.method[redDim], dims = c(dim.x, dim.y, dim.z), cell_size = cell_size_3d, show_trajectory_graph = show_trajectory, color_cells_by = color_cells_by.attr)
                                                                        }
                                                                    } else {
                                                                        cds_subpopulation_3d_plot_obj <- plot_cells_3d(cds_subpopulation_3d.plt, reduction_method = redDim.method[redDim], dims = c(dim.x, dim.y, dim.z), cell_size = cell_size_3d, show_trajectory_graph = show_trajectory, color_cells_by = color_cells_by.attr)
                                                                    }
                                                                    
                                                                    if (grid.plt == 0) {
                                                                        cds_subpopulation_3d_plot_obj <- cds_subpopulation_3d_plot_obj %>% plotly::layout(scene = scene)
                                                                    }
                                                                    
                                                                    # https://community.rstudio.com/t/save-viewer-object-rendered-in-rstudio-as-image/32796
                                                                    saveWidget(cds_subpopulation_3d_plot_obj, file.html)
                                                                    #webshot(file.html, file.png)
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    
                                    cds_subpopulation_3d = cds_subpopulation_3d.origin.processed
                                }
                            }
                        }
                    }
                    
                    
                    
                    # === Differential expression analysis ==============================================================================================================
                    # https://cole-trapnell-lab.github.io/monocle3/docs/differential/
                    
                    # There are two approaches for differential analysis in Monocle:
                    # Regression analysis: using fit_models(), you can evaluate whether each gene depends on variables such as time, treatments, etc.
                    # Graph-autocorrelation analysis: using graph_test(), you can find genes that vary over a trajectory or between clusters.
                    
                    # Monocle also comes with specialized functions for finding co-regulated modules of differentially expressed genes. 
                    # Monocle also allows you to interactively interrogate specific clusters or regions of a trajectory (e.g. branch points) for genes that vary within them.
                    
                    # (Optional) Perform differential expression analysis: https://cole-trapnell-lab.github.io/monocle3/docs/differential/
                    
                    # *** query gene(s) ***
                    if (! is.null(genes_of_interest.input)) {
                        # /// query.gene_group_df ///
                        genes_of_interest.new.all <- list()
                        for (gX in 1:length(genes_of_interest.input)) {
                            marker.genes <- genes_of_interest.input[[gX]]
                            marker.genes.idx <- which(sapply(marker.genes, FUN = function(X) X %in% rowData(cds)[,"gene_short_name"]))
                            if (length(marker.genes.idx) > 0) {
                                marker.genes = marker.genes[marker.genes.idx]
                                
                                genes_of_interest.new <- list(marker.genes)
                                names(genes_of_interest.new) = names(genes_of_interest.input[gX])
                                
                                genes_of_interest.new.all <- c(genes_of_interest.new.all, genes_of_interest.new)
                            }
                        }
                        genes_of_interest.new = genes_of_interest.new.all
                        
                        query.gene_group_df.all <- c()
                        for (gX in 1:length(genes_of_interest.new)) {
                            marker.genes <- genes_of_interest.new[[gX]]
                            marker.genes.idx <- which(sapply(marker.genes, FUN = function(X) X %in% rowData(cds)[,"gene_short_name"]))
                            if (length(marker.genes.idx) > 0) {
                                marker.genes = marker.genes[marker.genes.idx]
                                
                                query.gene_group_df <- data.frame("id" = marker.genes, 
                                                                  "module" = names(genes_of_interest.new[gX]), 
                                                                  "supermodule" = "genes_of_interest", 
                                                                  stringsAsFactors = FALSE)
                                for (dim_ in 1:eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = '')))) {
                                    query.gene_group_df <- data.frame(query.gene_group_df, 
                                                                      as.numeric(NA), 
                                                                      stringsAsFactors = FALSE)
                                    colnames(query.gene_group_df)[ncol(query.gene_group_df)] = paste('dim_',dim_, sep = '')
                                }
                                query.gene_group_df <- data.frame(query.gene_group_df, 
                                                                  "gene_short_name" = marker.genes, 
                                                                  stringsAsFactors = FALSE)
                                
                                query.gene_group_df[,"module"] <- factor(query.gene_group_df[,"module"], levels = unique(query.gene_group_df[,"module"]))
                                query.gene_group_df[,"supermodule"] <- factor(query.gene_group_df[,"supermodule"], levels = unique(query.gene_group_df[,"supermodule"]))
                                
                                query.gene_group_df.all <- rbind(query.gene_group_df.all, query.gene_group_df)
                            }
                        }
                        #query.gene_group_df.all <- tbl_df(query.gene_group_df.all)
                        query.gene_group_df.all <- tibble::as_tibble(query.gene_group_df.all)
                    }
                    
                    # *** signature gene(s) ***
                    if (! is.null(signature_genes.input)) {
                        # /// signature.gene_group_df ///
                        signature_genes.new.all <- list()
                        for (gX in 1:length(signature_genes.input)) {
                            marker.genes <- signature_genes.input[[gX]]
                            marker.genes.idx <- which(sapply(marker.genes, FUN = function(X) X %in% rowData(cds)[,"gene_short_name"]))
                            if (length(marker.genes.idx) > 0) {
                                marker.genes = marker.genes[marker.genes.idx]
                                
                                signature_genes.new <- list(marker.genes)
                                names(signature_genes.new) = names(signature_genes.input[gX])
                                
                                signature_genes.new.all <- c(signature_genes.new.all, signature_genes.new)
                            }
                        }
                        signature_genes.new = signature_genes.new.all
                        
                        signature.gene_group_df.all <- c()
                        for (gX in 1:length(signature_genes.new)) {
                            marker.genes <- signature_genes.new[[gX]]
                            marker.genes.idx <- which(sapply(marker.genes, FUN = function(X) X %in% rowData(cds)[,"gene_short_name"]))
                            if (length(marker.genes.idx) > 0) {
                                marker.genes = marker.genes[marker.genes.idx]
                                
                                signature.gene_group_df <- data.frame("id" = marker.genes, 
                                                                      "module" = names(signature_genes.new[gX]), 
                                                                      "supermodule" = "signature_genes", 
                                                                      stringsAsFactors = FALSE)
                                for (dim_ in 1:eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = '')))) {
                                    signature.gene_group_df <- data.frame(signature.gene_group_df, 
                                                                          as.numeric(NA), 
                                                                          stringsAsFactors = FALSE)
                                    colnames(signature.gene_group_df)[ncol(signature.gene_group_df)] = paste('dim_',dim_, sep = '')
                                }
                                signature.gene_group_df <- data.frame(signature.gene_group_df, 
                                                                      "gene_short_name" = marker.genes, 
                                                                      stringsAsFactors = FALSE)
                                
                                signature.gene_group_df[,"module"] <- factor(signature.gene_group_df[,"module"], levels = unique(signature.gene_group_df[,"module"]))
                                signature.gene_group_df[,"supermodule"] <- factor(signature.gene_group_df[,"supermodule"], levels = unique(signature.gene_group_df[,"supermodule"]))
                                
                                signature.gene_group_df.all <- rbind(signature.gene_group_df.all, signature.gene_group_df)
                            }
                        }
                        #signature.gene_group_df.all <- tbl_df(signature.gene_group_df.all)
                        signature.gene_group_df.all <- tibble::as_tibble(signature.gene_group_df.all)
                    }
                    
                    
                    for (subpopulation.idx in 0:length(subpopulation.type)) {
                        
                        cds_subpopulation.origin.processed = eval(parse(text=paste('cds_subpopulation','.',subpopulation.idx,'.origin.processed', sep = '')))
                        cds_subpopulation_3d.origin.processed = eval(parse(text=paste('cds_subpopulation_3d','.',subpopulation.idx,'.origin.processed', sep = '')))
                        
                        for (redDim in 1:length(redDim.method)) {
                            
                            # /// subset ///
                            for (SUB in SUB.run) {
                                if (SUB == 0) {
                                    #next # next halts the processing of the current iteration and advances the looping index.
                                } else if (SUB == 1) {
                                    SUB.prefolder = 'cell_subset'
                                    SUB.prename = 'cell_subset'
                                } else if (SUB == 2) {
                                    if (length(cell_subset) == 0 & cell_choose == TRUE) {
                                        SUB.prefolder = 'choose'
                                        SUB.prename = 'choose'
                                    } else if (length(cell_subset) > 0 & cell_choose == TRUE) {
                                        SUB.prefolder = paste('cell_subset','/','choose', sep = '')
                                        SUB.prename = paste('cell_subset','.','choose', sep = '')
                                    }
                                }
                                
                                if (SUB == 1) {
                                    cds_subpopulation_choose_3d.cell_subset = eval(parse(text=paste('cds_subpopulation_choose_3d','.',subpopulation.idx,'.cell_subset', sep = '')))
                                    cds_subpopulation_choose.cell_subset = eval(parse(text=paste('cds_subpopulation_choose','.',subpopulation.idx,'.cell_subset', sep = '')))
                                } else if (SUB == 2) {
                                    cds_subpopulation_choose_3d.cell_choose = eval(parse(text=paste('cds_subpopulation_choose_3d','.',subpopulation.idx,'.cell_choose', sep = '')))
                                    cds_subpopulation_choose.cell_choose = eval(parse(text=paste('cds_subpopulation_choose','.',subpopulation.idx,'.cell_choose', sep = '')))
                                }
                                
                                if (eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) == 2) {
                                    if (SUB == 1) {
                                        cds_subpopulation_choose_3d <- cds_subpopulation_choose_3d.cell_subset
                                    } else if (SUB == 2) {
                                        cds_subpopulation_choose_3d <- cds_subpopulation_choose_3d.cell_choose
                                    }
                                } else {
                                    if (SUB == 1) {
                                        cds_subpopulation_choose_3d <- cds_subpopulation_choose.cell_subset
                                    } else if (SUB == 2) {
                                        cds_subpopulation_choose_3d <- cds_subpopulation_choose.cell_choose
                                    }
                                }
                                
                                if (SUB == 0) {
                                    cds_subpopulation = cds_subpopulation.origin.processed
                                    cds_subpopulation_3d = cds_subpopulation_3d.origin.processed
                                    file.prename = 'clustering'
                                } else if (SUB == 1) {
                                    cds_subpopulation = cds_subpopulation_choose.cell_subset
                                    cds_subpopulation_3d = cds_subpopulation_choose_3d
                                    file.prename = paste(SUB.prename,'.','clustering', sep = '')
                                } else if (SUB == 2) {
                                    cds_subpopulation = cds_subpopulation_choose.cell_choose
                                    cds_subpopulation_3d = cds_subpopulation_choose_3d
                                    file.prename = paste(SUB.prename,'.','clustering', sep = '')
                                }
                                
                                
                                if (ncol(cds_subpopulation) > 0) {
                                    
                                    if (length(subpopulation.type[subpopulation.idx]) > 0) {
                                        subpopulation.folder <- subpopulation.type[subpopulation.idx]
                                    } else {
                                        subpopulation.folder <- project.subname
                                    }
                                    
                                    dim.folder = paste(eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))),'D', sep = '')
                                    
                                    mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim], sep = '')
                                    subDir <- dim.folder
                                    if (! file.exists(file.path(mainDir, subDir))) {
                                        dir.create(file.path(mainDir, subDir))
                                    }
                                    
                                    mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder, sep = '')
                                    if (eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) == 2) {
                                        #differential.folder <- paste('D12','/','differential', sep = '')
                                        differential.folder <- 'differential'
                                    } else {
                                        differential.folder <- 'differential'
                                    }
                                    subDir <- differential.folder
                                    if (! file.exists(file.path(mainDir, subDir))) {
                                        dir.create(file.path(mainDir, subDir))
                                    }
                                    
                                    if (SUB == 0) {
                                        differential.prefolder = differential.folder
                                    } else {
                                        mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',differential.folder, sep = '')
                                        subDir <- SUB.prefolder
                                        if (! file.exists(file.path(mainDir, subDir))) {
                                            dir.create(file.path(mainDir, subDir), recursive = TRUE)
                                        }
                                        differential.prefolder = paste(differential.folder,'/',SUB.prefolder, sep = '')
                                    }
                                    
                                    mapping.idx = sapply(rownames(colData(cds_subpopulation)), FUN = function(X) which(names(cds_subpopulation@clusters[[redDim.method[redDim]]]$clusters) %in% X))
                                    cds_subpopulation@clusters[[redDim.method[redDim]]]$clusters = cds_subpopulation@clusters[[redDim.method[redDim]]]$clusters[mapping.idx]
                                    cds_subpopulation@clusters[[redDim.method[redDim]]]$clusters = factor(cds_subpopulation@clusters[[redDim.method[redDim]]]$clusters, levels = mixedsort(as.character(unique(cds_subpopulation@clusters[[redDim.method[redDim]]]$clusters))))
                                    colData(cds_subpopulation)$cluster = cds_subpopulation@clusters[[redDim.method[redDim]]]$clusters
                                    
                                    mapping.idx = sapply(rownames(colData(cds_subpopulation)), FUN = function(X) which(names(cds_subpopulation@clusters[[redDim.method[redDim]]]$partitions) %in% X))
                                    cds_subpopulation@clusters[[redDim.method[redDim]]]$partitions = cds_subpopulation@clusters[[redDim.method[redDim]]]$partitions[mapping.idx]
                                    cds_subpopulation@clusters[[redDim.method[redDim]]]$partitions = factor(cds_subpopulation@clusters[[redDim.method[redDim]]]$partitions, levels = mixedsort(as.character(unique(cds_subpopulation@clusters[[redDim.method[redDim]]]$partitions))))
                                    colData(cds_subpopulation)$partition = cds_subpopulation@clusters[[redDim.method[redDim]]]$partitions
                                    
                                    mapping.idx = sapply(rownames(colData(cds_subpopulation)), FUN = function(X) which(rownames(cds_subpopulation@principal_graph_aux[[redDim.method[redDim]]]$pr_graph_cell_proj_closest_vertex) %in% X))
                                    cds_subpopulation@principal_graph_aux[[redDim.method[redDim]]]$pr_graph_cell_proj_closest_vertex = cds_subpopulation@principal_graph_aux[[redDim.method[redDim]]]$pr_graph_cell_proj_closest_vertex[mapping.idx, , drop = FALSE]
                                    
                                    
                                    # To investigate which genes are expressed differentially across the clusters, we could use the regression analysis tools discussed above. 
                                    # However, Monocle provides an alternative way of finding genes that vary between groups of cells in UMAP or t-SNE space.
                                    
                                    if (! project.name %in% c("GSM4150377", "GSM4150378")) {
                                        
                                        # The function graph_test() uses a statistic from spatial autocorrelation analysis called Moran's I, which Cao & Spielmann et al showed to be effective in finding genes that vary in single-cell RNA-seq datasets.
                                        # Moran's I: https://en.wikipedia.org/wiki/Moran%27s_I
                                        # The data frame pr_graph_test_res has the Moran's I test results for each gene in the cell_data_set. 
                                        # Significant values much less than zero are generally rare.
                                        
                                        if (redDim.method[redDim] == "UMAP") {
                                            # /// neighbor_graph = "knn" ///
                                            # reduction_method: character, the method used to reduce dimension. Currently only supported for "UMAP".
                                            pr_graph_test_res.knn <- graph_test(cds_subpopulation, reduction_method = "UMAP", neighbor_graph = "knn", cores = 8) # neighbor_graph: String indicating what neighbor graph to use. "principal_graph" and "knn" are supported. Default is "knn"
                                            
                                            # If you'd like to rank the genes by effect size, sort this table by the morans_I column, which ranges from -1 to +1. 
                                            morans_I.orderidx = order(pr_graph_test_res.knn[,"morans_I"], decreasing = TRUE)
                                            pr_graph_test_res.knn <- pr_graph_test_res.knn[morans_I.orderidx,]
                                            
                                            # A value of 0 indicates no effect, while +1 indicates perfect positive autocorrelation and suggests that nearby cells have very similar values of a gene's expression. 
                                            # Positive values indicate a gene is expressed in a focal region of the UMAP space (e.g. specific to one or more clusters).
                                            pr_graph_test_res.knn.pos <- pr_graph_test_res.knn[which(pr_graph_test_res.knn[,"morans_I"] > 0),]
                                            pr_graph_test_res.knn.zero <- pr_graph_test_res.knn[which(pr_graph_test_res.knn[,"morans_I"] == 0),]
                                            pr_graph_test_res.knn.neg <- pr_graph_test_res.knn[which(pr_graph_test_res.knn[,"morans_I"] < 0),]
                                            pr_graph_test_res.knn.NA <- pr_graph_test_res.knn[which(is.na(pr_graph_test_res.knn[,"morans_I"])),]
                                            
                                            pr_deg.knn <- subset(pr_graph_test_res.knn, q_value < 0.05)
                                            pr_deg.knn.pos <- subset(pr_graph_test_res.knn.pos, q_value < 0.05)
                                            pr_deg.knn.zero <- subset(pr_graph_test_res.knn.zero, q_value < 0.05)
                                            pr_deg.knn.neg <- subset(pr_graph_test_res.knn.neg, q_value < 0.05)
                                            pr_deg.knn.NA <- subset(pr_graph_test_res.knn.NA, q_value < 0.05)
                                            
                                            pr_deg_ids.knn <- rownames(pr_deg.knn)
                                            pr_deg_ids.knn.pos <- rownames(pr_deg.knn.pos)
                                            pr_deg_ids.knn.zero <- rownames(pr_deg.knn.zero)
                                            pr_deg_ids.knn.neg <- rownames(pr_deg.knn.neg)
                                            pr_deg_ids.knn.NA <- rownames(pr_deg.knn.NA)
                                            
                                            pr_deg_names.knn <- as.character(pr_deg.knn[,"gene_short_name"])
                                            pr_deg_names.knn.pos <- as.character(pr_deg.knn.pos[,"gene_short_name"])
                                            pr_deg_names.knn.zero <- as.character(pr_deg.knn.zero[,"gene_short_name"])
                                            pr_deg_names.knn.neg <- as.character(pr_deg.knn.neg[,"gene_short_name"])
                                            pr_deg_names.knn.NA <- as.character(pr_deg.knn.NA[,"gene_short_name"])
                                        }
                                    }
                                    
                                    
                                    # Here are a couple of interesting genes that score as highly significant according to graph_test():
                                    # /// marker.genes ///
                                    gene_set.list <- list()
                                    
                                    if (exists("pr_deg_names.knn")) {
                                        # *** Select serial list ***
                                        serial_num <- c(1:top_gene_num)
                                        serial_num.idx <- which(serial_num <= length(pr_deg_names.knn))
                                        if (length(serial_num.idx) > 0) {
                                            serial_num <- serial_num[serial_num.idx]
                                        } else {
                                            serial_num <- c()
                                        }
                                        if (length(serial_num) > 0) {
                                            for (serial in 1:length(serial_num)) {
                                                serial_genes <- list(pr_deg_names.knn[serial_num[serial]])
                                                names(serial_genes) = paste('gene',serial_num[serial], sep = '')
                                                
                                                gene_set.list <- c(gene_set.list, serial_genes)
                                            }
                                        }
                                        
                                        # *** Select top list ***
                                        top_num <- c(1, 3, 5, 10, 20, 30, 50)
                                        top_num.idx <- which(top_num <= top_gene_num)
                                        if (length(top_num.idx) > 0) {
                                            top_num <- top_num[top_num.idx]
                                        } else {
                                            top_num <- c(1)
                                        }
                                        top_num.idx <- which(top_num <= length(pr_deg_names.knn))
                                        if (length(top_num.idx) > 0) {
                                            top_num <- top_num[top_num.idx]
                                        } else {
                                            top_num <- c()
                                        }
                                        if (length(top_num) > 0) {
                                            for (top in 1:length(top_num)) {
                                                top_genes <- list(unique(head(pr_deg_names.knn, top_num[top])))
                                                if (top_num[top] == 1) {
                                                    names(top_genes) = paste('top',top_num[top],'_','gene', sep = '')
                                                } else if (top_num[top] > 1) {
                                                    names(top_genes) = paste('top',top_num[top],'_','genes', sep = '')
                                                }
                                                
                                                gene_set.list <- c(gene_set.list, top_genes)
                                            }
                                        }
                                    }
                                    
                                    # *** query gene(s) ***
                                    if (! is.null(genes_of_interest.input)) {
                                        for (gX in 1:length(genes_of_interest.input)) {
                                            for (gY in 1:length(genes_of_interest.input[[gX]])) {
                                                genes_of_interest.gX.gY <- list(genes_of_interest.input[[gX]][gY])
                                                names(genes_of_interest.gX.gY) = paste(names(genes_of_interest.input[gX]),'.',genes_of_interest.input[[gX]][gY], sep = '')
                                                gene_set.list <- c(gene_set.list, genes_of_interest.gX.gY)
                                            }
                                        }
                                        genes_of_interest <- genes_of_interest.input
                                        gene_set.list <- c(gene_set.list, genes_of_interest)
                                    }
                                    
                                    # *** signature gene(s) ***
                                    if (! is.null(signature_genes.input)) {
                                        for (gX in 1:length(signature_genes.input)) {
                                            for (gY in 1:length(signature_genes.input[[gX]])) {
                                                signature_genes.gX.gY <- list(signature_genes.input[[gX]][gY])
                                                names(signature_genes.gX.gY) = paste(names(signature_genes.input[gX]),'.',signature_genes.input[[gX]][gY], sep = '')
                                                gene_set.list <- c(gene_set.list, signature_genes.gX.gY)
                                            }
                                        }
                                        signature_genes <- signature_genes.input
                                        gene_set.list <- c(gene_set.list, signature_genes)
                                    }
                                    
                                    if (length(gene_set.list) > 0) {
                                        for (N in 1:length(gene_set.list)) {
                                            #print(names(gene_set.list[N]))
                                            
                                            marker.genes <- gene_set.list[[N]]
                                            marker.genes.idx <- which(sapply(marker.genes, FUN = function(X) X %in% rowData(cds_subpopulation)[,"gene_short_name"]))
                                            if (length(marker.genes.idx) > 0) {
                                                marker.genes = marker.genes[marker.genes.idx]
                                                
                                                panel.num = length(unlist(unique(marker.genes)))
                                                if (panel.num == 1) {
                                                    width.size = 7
                                                    height.size = 6
                                                } else if (panel.num == 2) {
                                                    width.size = 7 + (5 * 1)
                                                    height.size = 6
                                                } else if (panel.num == 3) {
                                                    width.size = 7 + (5 * 2)
                                                    height.size = 6
                                                } else if (panel.num > 3) {
                                                    width.size = 7 + (5 * 3)
                                                    height.size = width.size - 1
                                                }
                                                
                                                #for (redDim in 1:length(redDim.method)) {
                                                
                                                # /// 2d ///
                                                dim.folder = paste(eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))),'D', sep = '')
                                                dim.finalX = eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) - (2-1)
                                                
                                                if (dim.finalX >= 1) {
                                                    
                                                    for (dim.comb in 1:dim.finalX) {
                                                        dim.x <- dim.comb
                                                        dim.y <- dim.x + 1
                                                        
                                                        if (dim.finalX > 1) {
                                                            if (dim.x == 1) {
                                                                comb.run = 2
                                                            } else {
                                                                comb.run = 1
                                                            }
                                                        } else {
                                                            comb.run = 1
                                                        }
                                                        
                                                        for (dim.comb.run in 1:comb.run) {
                                                            if (dim.comb.run > 1) {
                                                                dim.y <- dim.y + 1
                                                            }
                                                            
                                                            if (eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) == 2) {
                                                                dim.xy <- paste('D',dim.x,'-','D',dim.y, sep = '')
                                                            } else {
                                                                dim.xy <- paste('D',dim.x,'-','D',dim.y, sep = '')
                                                            }
                                                            #print(dim.xy)
                                                            
                                                            if (SUB == 0) {
                                                                file.prefolder = dim.xy
                                                            } else {
                                                                file.prefolder = paste(dim.xy,'/',SUB.prefolder, sep = '')
                                                            }
                                                            
                                                            mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder, sep = '')
                                                            subDir <- 'differential'
                                                            if (! file.exists(file.path(mainDir, subDir))) {
                                                                dir.create(file.path(mainDir, subDir))
                                                            }
                                                            mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential', sep = '')
                                                            subDir <- 'gene'
                                                            if (! file.exists(file.path(mainDir, subDir))) {
                                                                dir.create(file.path(mainDir, subDir))
                                                            }
                                                            
                                                            if (length(grep(paste(genes_of_interest.titles, collapse = '|'), names(gene_set.list[N]), value = TRUE, invert = FALSE)) > 0) {
                                                                mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential', sep = '')
                                                                subDir <- 'genes_of_interest'
                                                                if (! file.exists(file.path(mainDir, subDir))) {
                                                                    dir.create(file.path(mainDir, subDir))
                                                                }
                                                                mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential','/','genes_of_interest', sep = '')
                                                                subDir <- unlist(strsplit(names(gene_set.list[N]), '\\.'))[1]
                                                                if (! file.exists(file.path(mainDir, subDir))) {
                                                                    dir.create(file.path(mainDir, subDir))
                                                                }
                                                                mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential','/','genes_of_interest','/',unlist(strsplit(names(gene_set.list[N]), '\\.'))[1], sep = '')
                                                                subDir <- 'gene'
                                                                if (! file.exists(file.path(mainDir, subDir))) {
                                                                    dir.create(file.path(mainDir, subDir))
                                                                }
                                                                
                                                                file.folder = paste('genes_of_interest','/',unlist(strsplit(names(gene_set.list[N]), '\\.'))[1],'/','gene', sep = '')
                                                                if (names(gene_set.list[N]) %in% genes_of_interest.titles) {
                                                                    file.name.suffix = names(gene_set.list[N])
                                                                } else {
                                                                    file.name.suffix = gsub("\\.","_",unlist(strsplit(names(gene_set.list[N]), '\\.'))[2])
                                                                }
                                                            } else if (length(grep(paste(signature_genes.titles, collapse = '|'), names(gene_set.list[N]), value = TRUE, invert = FALSE)) > 0) {
                                                                mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential', sep = '')
                                                                subDir <- 'signature_genes'
                                                                if (! file.exists(file.path(mainDir, subDir))) {
                                                                    dir.create(file.path(mainDir, subDir))
                                                                }
                                                                mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential','/','signature_genes', sep = '')
                                                                subDir <- unlist(strsplit(names(gene_set.list[N]), '\\.'))[1]
                                                                if (! file.exists(file.path(mainDir, subDir))) {
                                                                    dir.create(file.path(mainDir, subDir))
                                                                }
                                                                mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential','/','signature_genes','/',unlist(strsplit(names(gene_set.list[N]), '\\.'))[1], sep = '')
                                                                subDir <- 'gene'
                                                                if (! file.exists(file.path(mainDir, subDir))) {
                                                                    dir.create(file.path(mainDir, subDir))
                                                                }
                                                                
                                                                file.folder = paste('signature_genes','/',unlist(strsplit(names(gene_set.list[N]), '\\.'))[1],'/','gene', sep = '')
                                                                if (names(gene_set.list[N]) %in% signature_genes.titles) {
                                                                    file.name.suffix = names(gene_set.list[N])
                                                                } else {
                                                                    file.name.suffix = gsub("\\.","_",unlist(strsplit(names(gene_set.list[N]), '\\.'))[2])
                                                                }
                                                            } else {
                                                                file.folder = 'gene'
                                                                file.name.suffix = gsub("\\.","_",names(gene_set.list[N]))
                                                            }
                                                            
                                                            if (is.null(file.name.suffix)) {
                                                                file.name = paste(project,'.',file.prename,'.','plot_cells','.',redDim.method[redDim],'_',dim.xy,'.','differential', sep = '')
                                                            } else {
                                                                if (is.na(file.name.suffix) | file.name.suffix == "") {
                                                                    file.name = paste(project,'.',file.prename,'.','plot_cells','.',redDim.method[redDim],'_',dim.xy,'.','differential', sep = '')
                                                                } else {
                                                                    file.name = paste(project,'.',file.prename,'.','plot_cells','.',redDim.method[redDim],'_',dim.xy,'.','differential','.',file.name.suffix, sep = '')
                                                                }
                                                            }
                                                            
                                                            png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential','/',file.folder,'/',file.name,'.png', sep = ''), width = width.size, height = height.size, units = 'in', res = 600)
                                                            p <- plot_cells(cds_subpopulation, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, genes = marker.genes, norm_method = "log", label_cell_groups = FALSE, show_trajectory_graph = FALSE, label_leaves = FALSE)
                                                            print(p)
                                                            dev.off()
                                                        }
                                                    }
                                                }
                                                
                                                # /// 3d ///
                                                if (length(grep('top', names(gene_set.list[N]), value = TRUE, invert = TRUE)) > 0) {
                                                    
                                                    if (length(marker.genes) == 1) {
                                                        
                                                        if (eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) == 2) {
                                                            dim.folder = paste(eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))),'D', sep = '')
                                                            #dim.finalX = eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) - (3-1)
                                                            #dim.folder = paste(3,'D', sep = '')
                                                            dim.finalX = 3 - (3-1)
                                                        } else {
                                                            dim.folder = paste(eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))),'D', sep = '')
                                                            dim.finalX = eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) - (3-1)
                                                        }
                                                        
                                                        if (dim.finalX >= 1) {
                                                            
                                                            for (dim.comb in 1:dim.finalX) {
                                                                dim.x <- dim.comb
                                                                dim.y <- dim.x + 1
                                                                dim.z <- dim.y + 1
                                                                
                                                                if (eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) == 2) {
                                                                    dim.xyz <- paste('D',dim.x,'-','D',dim.y,'-','D',dim.z, sep = '')
                                                                } else {
                                                                    dim.xyz <- paste('D',dim.x,'-','D',dim.y,'-','D',dim.z, sep = '')
                                                                }
                                                                #print(dim.xyz)
                                                                
                                                                if (SUB == 0) {
                                                                    file.prefolder = dim.xyz
                                                                } else {
                                                                    file.prefolder = paste(dim.xyz,'/',SUB.prefolder, sep = '')
                                                                }
                                                                
                                                                mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder, sep = '')
                                                                subDir <- 'differential'
                                                                if (! file.exists(file.path(mainDir, subDir))) {
                                                                    dir.create(file.path(mainDir, subDir))
                                                                }
                                                                mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential', sep = '')
                                                                subDir <- 'gene'
                                                                if (! file.exists(file.path(mainDir, subDir))) {
                                                                    dir.create(file.path(mainDir, subDir))
                                                                }
                                                                
                                                                if (length(grep(paste(genes_of_interest.titles, collapse = '|'), names(gene_set.list[N]), value = TRUE, invert = FALSE)) > 0) {
                                                                    mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential', sep = '')
                                                                    subDir <- 'genes_of_interest'
                                                                    if (! file.exists(file.path(mainDir, subDir))) {
                                                                        dir.create(file.path(mainDir, subDir))
                                                                    }
                                                                    mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential','/','genes_of_interest', sep = '')
                                                                    subDir <- unlist(strsplit(names(gene_set.list[N]), '\\.'))[1]
                                                                    if (! file.exists(file.path(mainDir, subDir))) {
                                                                        dir.create(file.path(mainDir, subDir))
                                                                    }
                                                                    mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential','/','genes_of_interest','/',unlist(strsplit(names(gene_set.list[N]), '\\.'))[1], sep = '')
                                                                    subDir <- 'gene'
                                                                    if (! file.exists(file.path(mainDir, subDir))) {
                                                                        dir.create(file.path(mainDir, subDir))
                                                                    }
                                                                    
                                                                    file.folder = paste('genes_of_interest','/',unlist(strsplit(names(gene_set.list[N]), '\\.'))[1],'/','gene', sep = '')
                                                                    if (names(gene_set.list[N]) %in% genes_of_interest.titles) {
                                                                        file.name.suffix = names(gene_set.list[N])
                                                                    } else {
                                                                        file.name.suffix = gsub("\\.","_",unlist(strsplit(names(gene_set.list[N]), '\\.'))[2])
                                                                    }
                                                                } else if (length(grep(paste(signature_genes.titles, collapse = '|'), names(gene_set.list[N]), value = TRUE, invert = FALSE)) > 0) {
                                                                    mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential', sep = '')
                                                                    subDir <- 'signature_genes'
                                                                    if (! file.exists(file.path(mainDir, subDir))) {
                                                                        dir.create(file.path(mainDir, subDir))
                                                                    }
                                                                    mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential','/','signature_genes', sep = '')
                                                                    subDir <- unlist(strsplit(names(gene_set.list[N]), '\\.'))[1]
                                                                    if (! file.exists(file.path(mainDir, subDir))) {
                                                                        dir.create(file.path(mainDir, subDir))
                                                                    }
                                                                    mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential','/','signature_genes','/',unlist(strsplit(names(gene_set.list[N]), '\\.'))[1], sep = '')
                                                                    subDir <- 'gene'
                                                                    if (! file.exists(file.path(mainDir, subDir))) {
                                                                        dir.create(file.path(mainDir, subDir))
                                                                    }
                                                                    
                                                                    file.folder = paste('signature_genes','/',unlist(strsplit(names(gene_set.list[N]), '\\.'))[1],'/','gene', sep = '')
                                                                    if (names(gene_set.list[N]) %in% signature_genes.titles) {
                                                                        file.name.suffix = names(gene_set.list[N])
                                                                    } else {
                                                                        file.name.suffix = gsub("\\.","_",unlist(strsplit(names(gene_set.list[N]), '\\.'))[2])
                                                                    }
                                                                } else {
                                                                    file.folder = 'gene'
                                                                    file.name.suffix = gsub("\\.","_",names(gene_set.list[N]))
                                                                }
                                                                
                                                                #if (redDim.method[redDim] == "UMAP") {
                                                                #    trj.plt = 1
                                                                #} else if (redDim.method[redDim] == "tSNE") {
                                                                trj.plt = 0
                                                                #}
                                                                if (graph.run == 0) {
                                                                    trj.plt = 0
                                                                }
                                                                
                                                                for (trj in 0:trj.plt) {
                                                                    if (trj == 0) {
                                                                        trj.prt.endidx = 1
                                                                    } else if (trj == 1) {
                                                                        trj.prt.endidx = 0
                                                                    }
                                                                    
                                                                    #for (trj.prt in 1:1) {
                                                                    for (trj.prt in 1:trj.prt.endidx) {
                                                                        if (trj.prt == 1) {
                                                                            trj.prt.name = 'with_partition'
                                                                            cds_subpopulation_3d.plt = cds_subpopulation_3d
                                                                        } else if (trj.prt == 0) {
                                                                            #trj.prt.name = 'without_partition'
                                                                            #cds_subpopulation_3d.plt = cds_subpopulation_3d.no_partition
                                                                        }
                                                                        
                                                                        if (! is.null(cds_subpopulation_3d.plt)) {
                                                                            
                                                                            for (grid.plt in 1:0) {
                                                                                if (grid.plt == 1) {
                                                                                    grid.name = 'with_grid'
                                                                                } else if (grid.plt == 0) {
                                                                                    grid.name = 'without_grid'
                                                                                }
                                                                                
                                                                                if (trj == 0) {
                                                                                    show_trajectory = FALSE
                                                                                    if (length(marker.genes) == 1) {
                                                                                        if (is.null(file.name.suffix)) {
                                                                                            file.name = paste(project,'.',file.prename,'.','plot_cells_3d','.',redDim.method[redDim],'_',dim.xyz,'.','differential','.',marker.genes,'.','without_trajectory', sep = '')
                                                                                        } else {
                                                                                            if (is.na(file.name.suffix) | file.name.suffix == "") {
                                                                                                file.name = paste(project,'.',file.prename,'.','plot_cells_3d','.',redDim.method[redDim],'_',dim.xyz,'.','differential','.',marker.genes,'.','without_trajectory', sep = '')
                                                                                            } else {
                                                                                                file.name = paste(project,'.',file.prename,'.','plot_cells_3d','.',redDim.method[redDim],'_',dim.xyz,'.','differential','.',paste(unique(c(file.name.suffix, marker.genes)), collapse = '_'),'.','without_trajectory', sep = '')
                                                                                            }
                                                                                        }
                                                                                    } else {
                                                                                        if (is.null(file.name.suffix)) {
                                                                                            file.name = paste(project,'.',file.prename,'.','plot_cells_3d','.',redDim.method[redDim],'_',dim.xyz,'.','differential','.','without_trajectory', sep = '')
                                                                                        } else {
                                                                                            if (is.na(file.name.suffix) | file.name.suffix == "") {
                                                                                                file.name = paste(project,'.',file.prename,'.','plot_cells_3d','.',redDim.method[redDim],'_',dim.xyz,'.','differential','.','without_trajectory', sep = '')
                                                                                            } else {
                                                                                                file.name = paste(project,'.',file.prename,'.','plot_cells_3d','.',redDim.method[redDim],'_',dim.xyz,'.','differential','.',file.name.suffix,'.','without_trajectory', sep = '')
                                                                                            }
                                                                                        }
                                                                                    }
                                                                                    file.name = paste(file.name,'.',grid.name, sep = '')
                                                                                } else if (trj == 1) {
                                                                                    show_trajectory = TRUE
                                                                                    if (length(marker.genes) == 1) {
                                                                                        if (is.null(file.name.suffix)) {
                                                                                            file.name = paste(project,'.',file.prename,'.','plot_cells_3d','.',redDim.method[redDim],'_',dim.xyz,'.','differential','.',marker.genes,'.','with_trajectory', sep = '')
                                                                                        } else {
                                                                                            if (is.na(file.name.suffix) | file.name.suffix == "") {
                                                                                                file.name = paste(project,'.',file.prename,'.','plot_cells_3d','.',redDim.method[redDim],'_',dim.xyz,'.','differential','.',marker.genes,'.','with_trajectory', sep = '')
                                                                                            } else {
                                                                                                file.name = paste(project,'.',file.prename,'.','plot_cells_3d','.',redDim.method[redDim],'_',dim.xyz,'.','differential','.',paste(unique(c(file.name.suffix, marker.genes)), collapse = '_'),'.','with_trajectory', sep = '')
                                                                                            }
                                                                                        }
                                                                                    } else {
                                                                                        if (is.null(file.name.suffix)) {
                                                                                            file.name = paste(project,'.',file.prename,'.','plot_cells_3d','.',redDim.method[redDim],'_',dim.xyz,'.','differential','.','with_trajectory', sep = '')
                                                                                        } else {
                                                                                            if (is.na(file.name.suffix) | file.name.suffix == "") {
                                                                                                file.name = paste(project,'.',file.prename,'.','plot_cells_3d','.',redDim.method[redDim],'_',dim.xyz,'.','differential','.','with_trajectory', sep = '')
                                                                                            } else {
                                                                                                file.name = paste(project,'.',file.prename,'.','plot_cells_3d','.',redDim.method[redDim],'_',dim.xyz,'.','differential','.',file.name.suffix,'.','with_trajectory', sep = '')
                                                                                            }
                                                                                        }
                                                                                    }
                                                                                    #file.name = paste(file.name,'.',trj.prt.name, sep = '')
                                                                                    file.name = paste(file.name,'.',grid.name, sep = '')
                                                                                }
                                                                                
                                                                                file.html = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential','/',file.folder,'/',file.name,'.html', sep = '')
                                                                                file.png = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential','/',file.folder,'/',file.name,'.png', sep = '')
                                                                                
                                                                                cds_subpopulation_3d_plot_obj <- plot_cells_3d(cds_subpopulation_3d.plt, reduction_method = redDim.method[redDim], dims = c(dim.x, dim.y, dim.z), cell_size = cell_size_3d, genes = marker.genes, norm_method = "log", show_trajectory_graph = show_trajectory)
                                                                                
                                                                                if (grid.plt == 0) {
                                                                                    cds_subpopulation_3d_plot_obj <- cds_subpopulation_3d_plot_obj %>% plotly::layout(scene = scene)
                                                                                }
                                                                                
                                                                                # https://community.rstudio.com/t/save-viewer-object-rendered-in-rstudio-as-image/32796
                                                                                if (dim.comb == 1) { # only save D1-D2-D3 (no D2-D3-D4, D3-D4-D5, ...)
                                                                                    saveWidget(cds_subpopulation_3d_plot_obj, file.html)
                                                                                    #webshot(file.html, file.png)
                                                                                } else { # remove D2-D3-D4, D3-D4-D5, ...
                                                                                    mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential', sep = '')
                                                                                    if (file.exists(file.path(mainDir))) {
                                                                                        unlink(mainDir, recursive = TRUE)
                                                                                    }
                                                                                }
                                                                            }
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                                #}
                                            }
                                        }
                                    }
                                    
                                    
                                    # --- Finding modules of co-regulated genes -------------------------------------------------------
                                    # How do we associate genes with clusters? 
                                    # This section explains how to collect genes into modules that have similar patterns of expression and associate them with clusters.
                                    
                                    # Once you have a set of genes that vary in some interesting way across the clusters, Monocle provides a means of grouping them into modules.
                                    #for (redDim in 1:length(redDim.method)) {
                                    
                                    dim.folder = paste(eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))),'D', sep = '')
                                    
                                    mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim], sep = '')
                                    subDir <- dim.folder
                                    if (! file.exists(file.path(mainDir, subDir))) {
                                        dir.create(file.path(mainDir, subDir))
                                    }
                                    
                                    mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder, sep = '')
                                    if (eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) == 2) {
                                        #differential.folder <- paste('D12','/','differential', sep = '')
                                        differential.folder <- 'differential'
                                    } else {
                                        differential.folder <- 'differential'
                                    }
                                    subDir <- differential.folder
                                    if (! file.exists(file.path(mainDir, subDir))) {
                                        dir.create(file.path(mainDir, subDir))
                                    }
                                    
                                    if (SUB == 0) {
                                        differential.prefolder = differential.folder
                                    } else {
                                        mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',differential.folder, sep = '')
                                        subDir <- SUB.prefolder
                                        if (! file.exists(file.path(mainDir, subDir))) {
                                            dir.create(file.path(mainDir, subDir))
                                        }
                                        differential.prefolder = paste(differential.folder,'/',SUB.prefolder, sep = '')
                                    }
                                    
                                    if (redDim.method[redDim] == "UMAP") {
                                        if (exists("pr_deg_ids.knn") & exists("pr_deg_names.knn")) {
                                            if (length(pr_deg_ids.knn) > 0) {
                                                differential.run = 1
                                            } else {
                                                differential.run = 0
                                            }
                                        } else {
                                            differential.run = 0
                                        }
                                    } else if (redDim.method[redDim] == "tSNE") {
                                        differential.run = 0
                                    }
                                    
                                    if (differential.run == 1) {
                                        # You can call find_gene_modules(), which essentially runs UMAP on the genes (as opposed to the cells) and then groups them into modules using Louvain community analysis:
                                        # reduction_method: The dimensionality reduction method used to generate the lower dimensional space in which genes will be clustered. Currently only UMAP is supported.
                                        # max_components: The number of dimensions in which to cluster genes into modules.
                                        # resolution: Resolution parameter passed to Louvain. Can be a list. If so, this method will evaluate modularity at each resolution and use the one with the highest value.
                                        gene_group_df <- find_gene_modules(cds_subpopulation[pr_deg_ids.knn,], reduction_method = "UMAP", max_components = eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))), resolution = c(10^seq(-6,-1))) # The data frame gene_group_df contains a row for each gene and identifies the module it belongs to.
                                        
                                        pr_deg_ids.knn.idx <- sapply(unlist(gene_group_df[,"id"]), FUN = function(X) which(pr_deg_ids.knn %in% X))
                                        gene_group_df[,"gene_short_name"] = pr_deg_names.knn[pr_deg_ids.knn.idx]
                                        
                                        row.orderidx = order(as.numeric(unlist(gene_group_df[,"supermodule"])), 
                                                             as.numeric(unlist(gene_group_df[,"module"])), 
                                                             decreasing = FALSE)
                                        gene_group_df <- gene_group_df[row.orderidx, , drop = FALSE]
                                        gene_group_df.origin <- gene_group_df
                                        
                                        if (is.null(genes_of_interest.input) & is.null(signature_genes.input)) {
                                            gene_group_df <- gene_group_df.origin
                                        } else if (! is.null(genes_of_interest.input) & is.null(signature_genes.input)) {
                                            gene_group_df <- rbind(gene_group_df.origin, query.gene_group_df.all)
                                        } else if (is.null(genes_of_interest.input) & ! is.null(signature_genes.input)) {
                                            gene_group_df <- rbind(gene_group_df.origin, signature.gene_group_df.all)
                                        } else if (! is.null(genes_of_interest.input) & ! is.null(signature_genes.input)) {
                                            gene_group_df <- rbind(gene_group_df.origin, query.gene_group_df.all, signature.gene_group_df.all)
                                        }
                                        
                                        if (! is.null(gene_group_df)) {
                                            write.csv(gene_group_df, file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',differential.prefolder,'/',project,'.',gsub("\\.","_",subdata.folder),'.',gsub("\\.","_",subpopulation.folder),'.',redDim.method[redDim],'.',file.prename,'.','differential','.','geneModule','.csv', sep = ''), row.names = TRUE, quote = FALSE)
                                        }
                                        
                                        
                                        # To see which modules are expressed in which clusters or partitions you can use two different approaches for visualization. 
                                        
                                        # 1. The first is just to make a simple table that shows the aggregate expression of all genes in each module across all the clusters. 
                                        
                                        # Some modules are highly specific to certain partitions of cells, while others are shared across multiple partitions. 
                                        # Note that aggregate_gene_expression can work with arbitrary groupings of cells and genes. You're not limited to looking at modules from find_gene_modules(), clusters(), and partitions().
                                        if (length(cell_type) > 0) {
                                            for (color.plt in 1:length(cell_type)) {
                                                mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',differential.prefolder, sep = '')
                                                subDir <- gsub("\\.","_",cell_type[color.plt])
                                                if (! file.exists(file.path(mainDir, subDir))) {
                                                    dir.create(file.path(mainDir, subDir))
                                                }
                                                
                                                if (project.name == "tutorial") {
                                                    cell_group_df <- tibble::tibble(cell = rownames(colData(cds_subpopulation)), cell_group = partitions(cds_subpopulation, reduction_method = redDim.method[redDim])[colnames(cds_subpopulation)])
                                                } else {
                                                    cell_group_df <- tibble::tibble(cell = rownames(colData(cds_subpopulation)), cell_group = colData(cds_subpopulation)[, cell_type[color.plt]])
                                                }
                                                
                                                agg_mat <- aggregate_gene_expression(cds_subpopulation, gene_group_df, cell_group_df)
                                                
                                                a = grep(paste(c(genes_of_interest.titles, signature_genes.titles), collapse = '|'), rownames(agg_mat), value = TRUE, invert = TRUE)
                                                a.idx = grep(paste(c(genes_of_interest.titles, signature_genes.titles), collapse = '|'), rownames(agg_mat), value = FALSE, invert = TRUE)
                                                if (length(grep(paste(genes_of_interest.titles, collapse = '|'), rownames(agg_mat), value = TRUE, invert = FALSE)) > 0) {
                                                    b.genes_of_interest = grep(paste(genes_of_interest.titles, collapse = '|'), rownames(agg_mat), value = TRUE, invert = FALSE)
                                                    b.genes_of_interest.idx = grep(paste(genes_of_interest.titles, collapse = '|'), rownames(agg_mat), value = FALSE, invert = FALSE)
                                                } else {
                                                    b.genes_of_interest = NULL
                                                    b.genes_of_interest.idx = NULL
                                                }
                                                if (length(grep(paste(signature_genes.titles, collapse = '|'), rownames(agg_mat), value = TRUE, invert = FALSE)) > 0) {
                                                    b.signature_genes = grep(paste(signature_genes.titles, collapse = '|'), rownames(agg_mat), value = TRUE, invert = FALSE)
                                                    b.signature_genes.idx = grep(paste(signature_genes.titles, collapse = '|'), rownames(agg_mat), value = FALSE, invert = FALSE)
                                                } else {
                                                    b.signature_genes = NULL
                                                    b.signature_genes.idx = NULL
                                                }
                                                row.ordername = c(a[mixedorder(a, decreasing = FALSE)], 
                                                                  names(which(sapply(genes_of_interest.titles, FUN = function(X) X %in% b.genes_of_interest))), 
                                                                  names(which(sapply(signature_genes.titles, FUN = function(X) X %in% b.signature_genes))))
                                                row.orderidx = sapply(row.ordername, FUN = function(X) which(rownames(agg_mat) %in% X))
                                                
                                                if (project.name %in% c("Pancreas", "DrugPairing", "CRISPRINT")) {
                                                    cell_type.name <- c()
                                                    for (stg in 1:length(levels(colData(cds_subpopulation)[, stage.input]))) {
                                                        stage.name = as.character(unique(colData(cds_subpopulation)[which(colData(cds_subpopulation)[, stage.input] == levels(colData(cds_subpopulation)[, stage.input])[stg]), cell_type[color.plt]]))
                                                        cell_type.name <- c(cell_type.name, stage.name)
                                                    }
                                                    cell_type.name <- unique(cell_type.name)
                                                } else if (project.name %in% c("Paper")) {
                                                    cell_type.name <- mixedsort(as.character(unique(colData(cds_subpopulation)[, cell_type[color.plt]])))
                                                    cell_type.name <- cell_type.name[which(!is.na(cell_type.name))]
                                                    
                                                    if (exists("project.subname")) {
                                                        project.subname.prename = unlist(strsplit(project.subname, '\\.'))[1]
                                                        
                                                        if (project.subname.prename == "GSE147424" | project.subname.prename == "SRP253773") {
                                                            if (project.subname.prename == "GSE147424") {
                                                                sample_type.LV <- c("H", "NL", "LS")
                                                            } else if (project.subname.prename == "SRP253773") {
                                                                sample_type.LV <- c("Healthy", "Non-lesional", "Lesional")
                                                            }
                                                            sample_type.LV.idx <- which(sapply(sample_type.LV, FUN = function(X) X %in% cell_type.name))
                                                            
                                                            if (length(sample_type.LV.idx) > 0) {
                                                                sample_type.LV = sample_type.LV[sample_type.LV.idx]
                                                                cell_type.name <- sample_type.LV
                                                            }
                                                        }
                                                    }
                                                } else {
                                                    cell_type.name <- as.character(unique(colData(cds_subpopulation)[, cell_type[color.plt]]))
                                                    cell_type.name <- cell_type.name[which(!is.na(cell_type.name))]
                                                }
                                                if (project.name == "tutorial") {
                                                    col.orderidx = mixedorder(colnames(agg_mat), decreasing = FALSE)
                                                } else {
                                                    col.orderidx = sapply(cell_type.name, FUN = function(X) which(colnames(agg_mat) %in% X))
                                                }
                                                
                                                agg_mat <- agg_mat[row.orderidx, col.orderidx, drop = FALSE]
                                                
                                                rownames(agg_mat)[a.idx] <- stringr::str_c("Module ", rownames(agg_mat)[a.idx])
                                                if (! is.null(genes_of_interest.input)) {
                                                    for (gX in 1:length(genes_of_interest.input)) {
                                                        rownames(agg_mat)[which(rownames(agg_mat) == paste("geneSet",gX,"_of_interest", sep = ''))] = paste("geneSet",gX, sep = '')
                                                    }
                                                }
                                                rownames(agg_mat)[which(rownames(agg_mat) == "SINGH_KRAS_DEPENDENCY_SIGNATURE")] = "SINGH"
                                                rownames(agg_mat)[which(rownames(agg_mat) == "FRIDMAN_SENESCENCE_UP")] = "FRIDMAN"
                                                rownames(agg_mat)[which(rownames(agg_mat) == "HALLMARK_KRAS_SIGNALING_UP")] = "HALLMARK"
                                                
                                                colnames(agg_mat)[which(colnames(agg_mat) == "Acinar cell")] = "Acinar"
                                                colnames(agg_mat)[which(colnames(agg_mat) == "Ductal cell type 1")] = "Ductal 1"
                                                colnames(agg_mat)[which(colnames(agg_mat) == "Ductal cell type 2")] = "Ductal 2"
                                                colnames(agg_mat)[which(colnames(agg_mat) == "Endocrine cell")] = "Endocrine"
                                                colnames(agg_mat)[which(colnames(agg_mat) == "Endothelial cell")] = "Endothelial"
                                                colnames(agg_mat)[which(colnames(agg_mat) == "Fibroblast cell")] = "Fibroblast"
                                                colnames(agg_mat)[which(colnames(agg_mat) == "Macrophage cell")] = "Macrophage"
                                                colnames(agg_mat)[which(colnames(agg_mat) == "Stellate cell")] = "Stellate"
                                                colnames(agg_mat)[which(colnames(agg_mat) == "B cell")] = "B"
                                                colnames(agg_mat)[which(colnames(agg_mat) == "T cell")] = "T"
                                                
                                                if (exists("project.subname")) {
                                                    project.subname.prename = unlist(strsplit(project.subname, '\\.'))[1]
                                                    
                                                    if (project.subname.prename == "GSE147424" | project.subname.prename == "SRP253773") {
                                                        agg_mat.new <- agg_mat
                                                        
                                                        colnames(agg_mat.new) = gsub("sample|sample |Sample|Sample ","S",colnames(agg_mat.new))
                                                        
                                                        colnames(agg_mat.new)[which(colnames(agg_mat.new) == "Healthy")] = "H"
                                                        colnames(agg_mat.new)[which(colnames(agg_mat.new) == "Non-lesional")] = "NL"
                                                        colnames(agg_mat.new)[which(colnames(agg_mat.new) == "Lesional")] = "LS"
                                                        
                                                        agg_mat <- agg_mat.new
                                                    } else if (project.subname.prename == "NG-27918") {
                                                        agg_mat.new <- agg_mat
                                                        
                                                        #colnames(agg_mat.new) = gsub("sample|sample |Sample|Sample ","S",colnames(agg_mat.new))
                                                        
                                                        colnames(agg_mat.new) = gsub("unstimulated","unsti",colnames(agg_mat.new))
                                                        colnames(agg_mat.new) = gsub("stimulated","sti",colnames(agg_mat.new))
                                                        
                                                        agg_mat <- agg_mat.new
                                                    }
                                                }
                                                
                                                if (project.name == "tutorial") {
                                                    colnames(agg_mat) <- stringr::str_c("Partition ", colnames(agg_mat))
                                                } else {
                                                    colnames(agg_mat) <- stringr::str_c(colnames(agg_mat))
                                                }
                                                
                                                agg_mat <- as.data.frame(agg_mat)
                                                
                                                if (length(grep("CRA001160", project, value = TRUE)) > 0) {
                                                    if (cell_type[color.plt] == "sample") {
                                                        diagnoses2sample <- sapply(levels(cell_metadata[,"sample_pathologic_diagnoses"]), FUN = function(X) mixedsort(unique(cell_metadata[which(cell_metadata[,"sample_pathologic_diagnoses"] %in% X),"sample"])))
                                                        diagnoses2sample <- unlist(diagnoses2sample)
                                                        diagnoses2sample <- factor(diagnoses2sample, levels = unique(diagnoses2sample))
                                                        
                                                        agg_mat <- agg_mat[, as.character(diagnoses2sample)]
                                                    }
                                                }
                                                
                                                write.csv(agg_mat, file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',differential.prefolder,'/',gsub("\\.","_",cell_type[color.plt]),'/',project,'.',gsub("\\.","_",subdata.folder),'.',gsub("\\.","_",subpopulation.folder),'.',redDim.method[redDim],'.',file.prename,'.','differential','.','geneModule','.',gsub("\\.","_",cell_type[color.plt]),'.csv', sep = ''), row.names = TRUE, quote = FALSE)
                                                
                                                
                                                # /// pheatmap ///
                                                for (clu in 1:4) {
                                                    
                                                    if (clu == 1) {
                                                        row.clustering = FALSE
                                                        col.clustering = FALSE
                                                    } else if (clu == 2) {
                                                        row.clustering = FALSE
                                                        col.clustering = TRUE
                                                    } else if (clu == 3) {
                                                        row.clustering = TRUE
                                                        col.clustering = FALSE
                                                    } else if (clu == 4) {
                                                        row.clustering = TRUE
                                                        col.clustering = TRUE
                                                    }
                                                    
                                                    plt = 0
                                                    if (col.clustering == FALSE) {
                                                        plt = 1
                                                    } else if (col.clustering == TRUE) {
                                                        if (ncol(agg_mat) >= 2) {
                                                            plt = 1
                                                        }
                                                    }
                                                    
                                                    if (plt == 1) {
                                                        if (row.clustering == FALSE & col.clustering == FALSE) {
                                                            suffix.name = paste('geneModule','2',gsub("\\.","_",cell_type[color.plt]), sep = '')
                                                            width.size = (ncol(agg_mat) / 4)
                                                            height.size = (nrow(agg_mat) / 10)
                                                        } else if (row.clustering == FALSE & col.clustering == TRUE) {
                                                            suffix.name = paste('geneModule','2',gsub("\\.","_",cell_type[color.plt]),'_','clustering', sep = '')
                                                            width.size = (ncol(agg_mat) / 4)
                                                            height.size = (nrow(agg_mat) / 10) + 0.5
                                                        } else if (row.clustering == TRUE & col.clustering == FALSE) {
                                                            suffix.name = paste('geneModule','_','clustering','2',gsub("\\.","_",cell_type[color.plt]), sep = '')
                                                            width.size = (ncol(agg_mat) / 4) + 0.625
                                                            height.size = (nrow(agg_mat) / 10)
                                                        } else if (row.clustering == TRUE & col.clustering == TRUE) {
                                                            suffix.name = paste('geneModule','_','clustering','2',gsub("\\.","_",cell_type[color.plt]),'_','clustering', sep = '')
                                                            width.size = (ncol(agg_mat) / 4) + 0.625
                                                            height.size = (nrow(agg_mat) / 10) + 0.5
                                                        }
                                                        
                                                        if ((ncol(agg_mat) / 4) <= 1.5) {
                                                            width.size.fold = 1.5 / (ncol(agg_mat) / 4)
                                                            width.size = width.size * width.size.fold
                                                        }
                                                        
                                                        fontsize = 6
                                                        
                                                        png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',differential.prefolder,'/',gsub("\\.","_",cell_type[color.plt]),'/',project,'.',gsub("\\.","_",subdata.folder),'.',gsub("\\.","_",subpopulation.folder),'.',redDim.method[redDim],'.',file.prename,'.','differential','.','pheatmap','.',suffix.name,'.png', sep = ''), width = width.size, height = height.size, units = 'in', res = 600)
                                                        p <- pheatmap::pheatmap(agg_mat, cluster_rows = row.clustering, cluster_cols = col.clustering, 
                                                                                scale = "column", # scale: character indicating if the values should be centered and scaled in either the row direction or the column direction, or none. Corresponding values are "row", "column" and "none"
                                                                                clustering_method = "ward.D2", # hclust(): method: the agglomeration method to be used. This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
                                                                                fontsize = fontsize)
                                                        print(p)
                                                        dev.off()
                                                        
                                                        if (row.clustering == TRUE & col.clustering == FALSE) {
                                                            p.tree_row.order <- p$tree_row$order
                                                        }
                                                    }
                                                }
                                                
                                                
                                                if (ncol(agg_mat) > 1) {
                                                    # --- spaghetti -------------------------------------------------------------------------------------------------------------------------------------
                                                    # THE SPAGHETTI PLOT: https://www.data-to-viz.com/caveat/spaghetti.html
                                                    for (od in 0:1) {
                                                        if (od == 0) {
                                                            agg_mat.differential <- agg_mat
                                                        } else if (od == 1) {
                                                            agg_mat.differential <- agg_mat[p.tree_row.order, , drop = FALSE]
                                                        }
                                                        
                                                        agg_mat.differential <- data.frame("module" = rownames(agg_mat.differential), agg_mat.differential, check.names = FALSE)
                                                        colnames(agg_mat.differential)[1] = "module"
                                                        
                                                        if (length(grep("CRA001160", project, value = TRUE)) > 0) {
                                                            if (cell_type[color.plt] == "sample") {
                                                                diagnoses2sample <- sapply(levels(cell_metadata[,"sample_pathologic_diagnoses"]), FUN = function(X) mixedsort(unique(cell_metadata[which(cell_metadata[,"sample_pathologic_diagnoses"] %in% X),"sample"])))
                                                                diagnoses2sample <- unlist(diagnoses2sample)
                                                                diagnoses2sample <- factor(diagnoses2sample, levels = unique(diagnoses2sample))
                                                                
                                                                agg_mat.differential <- agg_mat.differential[, c("module", as.character(diagnoses2sample))]
                                                            }
                                                        }
                                                        
                                                        agg_mat.differential.origin <- agg_mat.differential
                                                        
                                                        agg_mat.differential.mtx <- agg_mat.differential
                                                        rownames(agg_mat.differential.mtx) = agg_mat.differential.mtx[,"module"]
                                                        module.idx <- which(colnames(agg_mat.differential.mtx) %in% "module")
                                                        agg_mat.differential.mtx <- agg_mat.differential.mtx[, -module.idx, drop = FALSE]
                                                        agg_mat.differential.mtx.origin <- agg_mat.differential.mtx
                                                        
                                                        
                                                        # --- trend -----------------------------------------------------------------------------------------------------------------------------------------
                                                        agg_mat.differential <- agg_mat.differential.origin
                                                        agg_mat.differential.mtx <- agg_mat.differential.mtx.origin
                                                        
                                                        if (length(grep("CRA001160", project, value = TRUE)) > 0) {
                                                            if (cell_type[color.plt] == "sample") {
                                                                sample_type2sample <- sapply(levels(cell_metadata[,"sample_type"]), FUN = function(X) mixedsort(unique(cell_metadata[which(cell_metadata[,"sample_type"] %in% X),"sample"])))
                                                                
                                                                col.idx1 <- unlist(sapply(sample_type2sample, FUN = function(X) which(colnames(agg_mat.differential.mtx) %in% X))[1])
                                                                col.idx2 <- unlist(sapply(sample_type2sample, FUN = function(X) which(colnames(agg_mat.differential.mtx) %in% X))[2])
                                                                
                                                                col.idx.order <- col.idx2
                                                            } else {
                                                                col.idx1 <- NULL
                                                                col.idx2 <- NULL
                                                                
                                                                col.idx.order <- c(1:ncol(agg_mat.differential.mtx))
                                                            }
                                                        } else {
                                                            col.idx1 <- NULL
                                                            col.idx2 <- NULL
                                                            
                                                            col.idx.order <- c(1:ncol(agg_mat.differential.mtx))
                                                        }
                                                        
                                                        low2high.idx <- c()
                                                        increasing.idx <- c()
                                                        low2high_increasing.idx <- c()
                                                        
                                                        high2low.idx <- c()
                                                        decreasing.idx <- c()
                                                        high2low_decreasing.idx <- c()
                                                        
                                                        for (m in 1:nrow(agg_mat.differential.mtx)) {
                                                            col.idx.inc <- col.idx.order[order(agg_mat.differential.mtx[m, col.idx.order], decreasing = FALSE)]
                                                            col.idx.dec <- col.idx.order[order(agg_mat.differential.mtx[m, col.idx.order], decreasing = TRUE)]
                                                            
                                                            # /// increasing ///
                                                            if (! is.null(col.idx1) & ! is.null(col.idx2)) {
                                                                if (max(agg_mat.differential.mtx[m, col.idx1], na.rm = TRUE) <= min(agg_mat.differential.mtx[m, col.idx2], na.rm = TRUE)) {
                                                                    #cat(rownames(agg_mat.differential.mtx[m, ]), sep = "\n")
                                                                    
                                                                    #cat("Group 1 <= Group 2\n")
                                                                    low2high.idx <- c(low2high.idx, m)
                                                                    
                                                                    if (all(col.idx.inc == col.idx.order)) {
                                                                        #cat("increasing\n")
                                                                        increasing.idx <- c(increasing.idx, m)
                                                                        
                                                                        low2high_increasing.idx <- c(low2high_increasing.idx, m)
                                                                    }
                                                                    
                                                                    #print(agg_mat.differential.mtx[m, ])
                                                                    #cat("\n")
                                                                    #cat("\n")
                                                                } else {
                                                                    if (all(col.idx.inc == col.idx.order)) {
                                                                        #cat(rownames(agg_mat.differential.mtx[m, ]), sep = "\n")
                                                                        
                                                                        #cat("increasing\n")
                                                                        increasing.idx <- c(increasing.idx, m)
                                                                        
                                                                        #print(agg_mat.differential.mtx[m, col.idx.order])
                                                                        #cat("\n")
                                                                        #cat("\n")
                                                                    }
                                                                }
                                                            } else {
                                                                if (all(col.idx.inc == col.idx.order)) {
                                                                    #cat(rownames(agg_mat.differential.mtx[m, ]), sep = "\n")
                                                                    
                                                                    #cat("increasing\n")
                                                                    increasing.idx <- c(increasing.idx, m)
                                                                    
                                                                    #print(agg_mat.differential.mtx[m, col.idx.order])
                                                                    #cat("\n")
                                                                    #cat("\n")
                                                                }
                                                            }
                                                            
                                                            # /// decreasing ///
                                                            if (! is.null(col.idx1) & ! is.null(col.idx2)) {
                                                                if (min(agg_mat.differential.mtx[m, col.idx1], na.rm = TRUE) >= max(agg_mat.differential.mtx[m, col.idx2], na.rm = TRUE)) {
                                                                    #cat(rownames(agg_mat.differential.mtx[m, ]), sep = "\n")
                                                                    
                                                                    #cat("Group 1 >= Group 2\n")
                                                                    high2low.idx <- c(high2low.idx, m)
                                                                    
                                                                    if (all(col.idx.dec == col.idx.order)) {
                                                                        #cat("decreasing\n")
                                                                        decreasing.idx <- c(decreasing.idx, m)
                                                                        
                                                                        high2low_decreasing.idx <- c(high2low_decreasing.idx, m)
                                                                    }
                                                                    
                                                                    #print(agg_mat.differential.mtx[m, ])
                                                                    #cat("\n")
                                                                    #cat("\n")
                                                                } else {
                                                                    if (all(col.idx.dec == col.idx.order)) {
                                                                        #cat(rownames(agg_mat.differential.mtx[m, ]), sep = "\n")
                                                                        
                                                                        #cat("decreasing\n")
                                                                        decreasing.idx <- c(decreasing.idx, m)
                                                                        
                                                                        #print(agg_mat.differential.mtx[m, col.idx.order])
                                                                        #cat("\n")
                                                                        #cat("\n")
                                                                    }
                                                                }
                                                            } else {
                                                                if (all(col.idx.dec == col.idx.order)) {
                                                                    #cat(rownames(agg_mat.differential.mtx[m, ]), sep = "\n")
                                                                    
                                                                    #cat("decreasing\n")
                                                                    decreasing.idx <- c(decreasing.idx, m)
                                                                    
                                                                    #print(agg_mat.differential.mtx[m, col.idx.order])
                                                                    #cat("\n")
                                                                    #cat("\n")
                                                                }
                                                            }
                                                        }
                                                        
                                                        if (od == 0) {
                                                            #cat("Group 1 <= Group 2\n")
                                                            #print(low2high.idx)
                                                            #cat("increasing\n")
                                                            #print(increasing.idx)
                                                            #cat("Group 1 <= Group 2 & increasing\n")
                                                            #print(low2high_increasing.idx)
                                                            #cat("\n")
                                                            
                                                            #cat("Group 1 >= Group 2\n")
                                                            #print(high2low.idx)
                                                            #cat("decreasing\n")
                                                            #print(decreasing.idx)
                                                            #cat("Group 1 >= Group 2 & decreasing\n")
                                                            #print(high2low_decreasing.idx)
                                                            #cat("\n")
                                                            
                                                            if (! is.null(low2high.idx)) {
                                                                agg_mat.differential.mtx.low2high <- agg_mat.differential.mtx[low2high.idx, , drop = FALSE]
                                                                write.csv(agg_mat.differential.mtx.low2high, file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',differential.prefolder,'/',gsub("\\.","_",cell_type[color.plt]),'/',project,'.',gsub("\\.","_",subdata.folder),'.',gsub("\\.","_",subpopulation.folder),'.',redDim.method[redDim],'.',file.prename,'.','differential','.','geneModule','.',gsub("\\.","_",cell_type[color.plt]),'.','low2high','.csv', sep = ''), row.names = TRUE, quote = FALSE)
                                                            }
                                                            if (! is.null(increasing.idx)) {
                                                                agg_mat.differential.mtx.inc <- agg_mat.differential.mtx[increasing.idx, , drop = FALSE]
                                                                write.csv(agg_mat.differential.mtx.inc, file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',differential.prefolder,'/',gsub("\\.","_",cell_type[color.plt]),'/',project,'.',gsub("\\.","_",subdata.folder),'.',gsub("\\.","_",subpopulation.folder),'.',redDim.method[redDim],'.',file.prename,'.','differential','.','geneModule','.',gsub("\\.","_",cell_type[color.plt]),'.','increasing','.csv', sep = ''), row.names = TRUE, quote = FALSE)
                                                            }
                                                            if (! is.null(low2high_increasing.idx)) {
                                                                agg_mat.differential.mtx.low2high_inc <- agg_mat.differential.mtx[low2high_increasing.idx, , drop = FALSE]
                                                                write.csv(agg_mat.differential.mtx.low2high_inc, file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',differential.prefolder,'/',gsub("\\.","_",cell_type[color.plt]),'/',project,'.',gsub("\\.","_",subdata.folder),'.',gsub("\\.","_",subpopulation.folder),'.',redDim.method[redDim],'.',file.prename,'.','differential','.','geneModule','.',gsub("\\.","_",cell_type[color.plt]),'.','low2high_increasing','.csv', sep = ''), row.names = TRUE, quote = FALSE)
                                                            }
                                                            
                                                            if (! is.null(high2low.idx)) {
                                                                agg_mat.differential.mtx.high2low <- agg_mat.differential.mtx[high2low.idx, , drop = FALSE]
                                                                write.csv(agg_mat.differential.mtx.high2low, file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',differential.prefolder,'/',gsub("\\.","_",cell_type[color.plt]),'/',project,'.',gsub("\\.","_",subdata.folder),'.',gsub("\\.","_",subpopulation.folder),'.',redDim.method[redDim],'.',file.prename,'.','differential','.','geneModule','.',gsub("\\.","_",cell_type[color.plt]),'.','high2low','.csv', sep = ''), row.names = TRUE, quote = FALSE)
                                                            }
                                                            if (! is.null(decreasing.idx)) {
                                                                agg_mat.differential.mtx.dec <- agg_mat.differential.mtx[decreasing.idx, , drop = FALSE]
                                                                write.csv(agg_mat.differential.mtx.dec, file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',differential.prefolder,'/',gsub("\\.","_",cell_type[color.plt]),'/',project,'.',gsub("\\.","_",subdata.folder),'.',gsub("\\.","_",subpopulation.folder),'.',redDim.method[redDim],'.',file.prename,'.','differential','.','geneModule','.',gsub("\\.","_",cell_type[color.plt]),'.','decreasing','.csv', sep = ''), row.names = TRUE, quote = FALSE)
                                                            }
                                                            if (! is.null(high2low_decreasing.idx)) {
                                                                agg_mat.differential.mtx.high2low_dec <- agg_mat.differential.mtx[high2low_decreasing.idx, , drop = FALSE]
                                                                write.csv(agg_mat.differential.mtx.high2low_dec, file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',differential.prefolder,'/',gsub("\\.","_",cell_type[color.plt]),'/',project,'.',gsub("\\.","_",subdata.folder),'.',gsub("\\.","_",subpopulation.folder),'.',redDim.method[redDim],'.',file.prename,'.','differential','.','geneModule','.',gsub("\\.","_",cell_type[color.plt]),'.','high2low_decreasing','.csv', sep = ''), row.names = TRUE, quote = FALSE)
                                                            }
                                                        }
                                                        
                                                        
                                                        # /// spaghetti plot ///
                                                        if (subpopulation.idx == 0) {
                                                            if (SUB == 0) {
                                                                main.title = paste(project.name, gsub("\\.","_",project.subname), sep = '.')
                                                            } else {
                                                                main.title = paste(paste(project.name, gsub("\\.","_",project.subname), sep = '.'),'.',SUB.prename, sep = '')
                                                            }
                                                        } else {
                                                            if (SUB == 0) {
                                                                main.title = paste(project.name, gsub("\\.","_",project.subname), gsub("\\.","_",subpopulation.type[subpopulation.idx]), sep = '.')
                                                            } else {
                                                                main.title = paste(paste(project.name, gsub("\\.","_",project.subname), gsub("\\.","_",subpopulation.type[subpopulation.idx]), sep = '.'),'.',SUB.prename, sep = '')
                                                            }
                                                        }
                                                        
                                                        xlab.name = unlist(strsplit(cell_type[color.plt], '_'))[1]
                                                        xlab.subname <- c()
                                                        if (length(unlist(strsplit(cell_type[color.plt], '_'))) > 1) {
                                                            xlab.subname <- c()
                                                            for (xn in 2:length(unlist(strsplit(cell_type[color.plt], '_')))) {
                                                                xlab.subname <- c(xlab.subname, unlist(strsplit(cell_type[color.plt], '_'))[xn])
                                                            }
                                                            xlab.subname <- paste(xlab.subname, collapse = '_')
                                                            xlab.name = paste(xlab.name, xlab.subname, sep = ' ')
                                                        }
                                                        xlab.title = str_to_title(xlab.name)
                                                        
                                                        for (trend in 0:6) {
                                                            agg_mat.differential <- agg_mat.differential.origin
                                                            agg_mat.differential.mtx <- agg_mat.differential.mtx.origin
                                                            
                                                            if (trend == 0) {
                                                                agg_mat.differential <- agg_mat.differential.origin
                                                                agg_mat.differential.mtx <- agg_mat.differential.mtx.origin
                                                                suffix.name = gsub("\\.","_",cell_type[color.plt])
                                                            } else if (trend == 1) {
                                                                if (! is.null(low2high.idx)) {
                                                                    agg_mat.differential <- agg_mat.differential.origin[low2high.idx, , drop = FALSE]
                                                                    agg_mat.differential.mtx <- agg_mat.differential.mtx.origin[low2high.idx, , drop = FALSE]
                                                                    suffix.name = paste(gsub("\\.","_",cell_type[color.plt]),'.','low2high', sep = '')
                                                                } else {
                                                                    next # next halts the processing of the current iteration and advances the looping index.
                                                                }
                                                            } else if (trend == 2) {
                                                                if (! is.null(increasing.idx)) {
                                                                    agg_mat.differential <- agg_mat.differential.origin[increasing.idx, , drop = FALSE]
                                                                    agg_mat.differential.mtx <- agg_mat.differential.mtx.origin[increasing.idx, , drop = FALSE]
                                                                    suffix.name = paste(gsub("\\.","_",cell_type[color.plt]),'.','increasing', sep = '')
                                                                } else {
                                                                    next # next halts the processing of the current iteration and advances the looping index.
                                                                }
                                                            } else if (trend == 3) {
                                                                if (! is.null(low2high_increasing.idx)) {
                                                                    agg_mat.differential <- agg_mat.differential.origin[low2high_increasing.idx, , drop = FALSE]
                                                                    agg_mat.differential.mtx <- agg_mat.differential.mtx.origin[low2high_increasing.idx, , drop = FALSE]
                                                                    suffix.name = paste(gsub("\\.","_",cell_type[color.plt]),'.','low2high_increasing', sep = '')
                                                                } else {
                                                                    next # next halts the processing of the current iteration and advances the looping index.
                                                                }
                                                            } else if (trend == 4) {
                                                                if (! is.null(high2low.idx)) {
                                                                    agg_mat.differential <- agg_mat.differential.origin[high2low.idx, , drop = FALSE]
                                                                    agg_mat.differential.mtx <- agg_mat.differential.mtx.origin[high2low.idx, , drop = FALSE]
                                                                    suffix.name = paste(gsub("\\.","_",cell_type[color.plt]),'.','high2low', sep = '')
                                                                } else {
                                                                    next # next halts the processing of the current iteration and advances the looping index.
                                                                }
                                                            } else if (trend == 5) {
                                                                if (! is.null(decreasing.idx)) {
                                                                    agg_mat.differential <- agg_mat.differential.origin[decreasing.idx, , drop = FALSE]
                                                                    agg_mat.differential.mtx <- agg_mat.differential.mtx.origin[decreasing.idx, , drop = FALSE]
                                                                    suffix.name = paste(gsub("\\.","_",cell_type[color.plt]),'.','decreasing', sep = '')
                                                                } else {
                                                                    next # next halts the processing of the current iteration and advances the looping index.
                                                                }
                                                            } else if (trend == 6) {
                                                                if (! is.null(high2low_decreasing.idx)) {
                                                                    agg_mat.differential <- agg_mat.differential.origin[high2low_decreasing.idx, , drop = FALSE]
                                                                    agg_mat.differential.mtx <- agg_mat.differential.mtx.origin[high2low_decreasing.idx, , drop = FALSE]
                                                                    suffix.name = paste(gsub("\\.","_",cell_type[color.plt]),'.','high2low_decreasing', sep = '')
                                                                } else {
                                                                    next # next halts the processing of the current iteration and advances the looping index.
                                                                }
                                                            }
                                                            
                                                            agg_mat.differential <- suppressMessages(melt(agg_mat.differential))
                                                            
                                                            agg_mat.differential[,"module"] <- factor(agg_mat.differential[,"module"], levels = unique(agg_mat.differential[,"module"]))
                                                            agg_mat.differential[,"variable"] <- factor(agg_mat.differential[,"variable"], levels = unique(agg_mat.differential[,"variable"]))
                                                            
                                                            colnames(agg_mat.differential) = str_to_title(colnames(agg_mat.differential))
                                                            
                                                            if (od == 0) {
                                                                legend.num = nrow(agg_mat.differential.mtx)
                                                                legend.ncol = ceiling(legend.num / 18)
                                                                legend.nrow = ceiling(legend.num / legend.ncol)
                                                                width.size = 10 + (0.5 * legend.ncol)
                                                                if (legend.nrow > 15) {
                                                                    height.size = 5 + ((legend.nrow - 15) * 0.2)
                                                                } else {
                                                                    height.size = 5
                                                                }
                                                                
                                                                png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',differential.prefolder,'/',gsub("\\.","_",cell_type[color.plt]),'/',project,'.',gsub("\\.","_",subdata.folder),'.',gsub("\\.","_",subpopulation.folder),'.',redDim.method[redDim],'.',file.prename,'.','differential','.','geneModule','.',suffix.name,'.png', sep = ''), width = width.size, height = height.size, units = 'in', res = 600)
                                                                p <- ggplot(agg_mat.differential, aes(x=Variable, y=Value, group=Module, color=Module)) + geom_line()
                                                                p <- p + theme_light()
                                                                #p <- p + theme_ipsum()
                                                                p <- p + ggtitle(main.title)
                                                                p <- p + xlab(xlab.title)
                                                                print(p)
                                                                dev.off()
                                                            }
                                                            
                                                            legend.num = nrow(agg_mat.differential.mtx)
                                                            legend.ncol = ceiling(sqrt(legend.num))
                                                            legend.nrow = ceiling(legend.num / legend.ncol)
                                                            Variable.num <- ncol(agg_mat.differential.mtx)
                                                            Value.range <- ceiling(max(agg_mat.differential.mtx, na.rm = TRUE)) - floor(min(agg_mat.differential.mtx, na.rm = TRUE))
                                                            width.size = legend.ncol * 2 + ceiling(Variable.num / 4) * 2 + 2
                                                            height.size = legend.nrow * 2 + ceiling(Value.range / 4) * 2 - 2
                                                            if (height.size < 5) {
                                                                if (legend.nrow <= 1) {
                                                                    height.size = 4
                                                                } else {
                                                                    height.size = 5
                                                                }
                                                            }
                                                            
                                                            if (od == 0) {
                                                                png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',differential.prefolder,'/',gsub("\\.","_",cell_type[color.plt]),'/',project,'.',gsub("\\.","_",subdata.folder),'.',gsub("\\.","_",subpopulation.folder),'.',redDim.method[redDim],'.',file.prename,'.','differential','.','geneModule','.',suffix.name,'.spaghetti','.png', sep = ''), width = width.size, height = height.size, units = 'in', res = 600)
                                                            } else if (od == 1) {
                                                                png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',differential.prefolder,'/',gsub("\\.","_",cell_type[color.plt]),'/',project,'.',gsub("\\.","_",subdata.folder),'.',gsub("\\.","_",subpopulation.folder),'.',redDim.method[redDim],'.',file.prename,'.','differential','.','geneModule','.',suffix.name,'.spaghetti','.','geneModule','_','clustering','.png', sep = ''), width = width.size, height = height.size, units = 'in', res = 600)
                                                            }
                                                            data <- agg_mat.differential
                                                            p <- data %>% mutate(Module2=Module)
                                                            p <- p %>% 
                                                                ggplot(aes(x=Variable, y=Value, group=Module)) + 
                                                                geom_line(data=p %>% dplyr::select(-Module), aes(group=Module2), color="grey", size=0.5, alpha=0.5) + 
                                                                geom_line(aes(color=Module), color="red", size=1.2) + 
                                                                scale_color_viridis(discrete = TRUE) + 
                                                                theme_ipsum() + 
                                                                theme(
                                                                    legend.position="none", 
                                                                    plot.title = element_text(size=22), 
                                                                    panel.grid = element_blank()
                                                                ) + 
                                                                ggtitle(main.title) + 
                                                                facet_wrap(~Module)
                                                            p <- p + xlab(xlab.title)
                                                            p <- p + 
                                                                theme(
                                                                    strip.text = element_text(size=20), 
                                                                    axis.title.x = element_text(size=18), 
                                                                    axis.title.y = element_text(size=18), 
                                                                    axis.text.x = element_text(size=14), 
                                                                    axis.text.y = element_text(size=14)
                                                                )
                                                            print(p)
                                                            dev.off()
                                                            
                                                            agg_mat.differential <- agg_mat.differential.origin
                                                            agg_mat.differential.mtx <- agg_mat.differential.mtx.origin
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    } else if (differential.run == 0) {
                                        if (is.null(genes_of_interest.input) & is.null(signature_genes.input)) {
                                            gene_group_df <- NULL
                                        } else if (! is.null(genes_of_interest.input) & is.null(signature_genes.input)) {
                                            gene_group_df <- query.gene_group_df.all
                                        } else if (is.null(genes_of_interest.input) & ! is.null(signature_genes.input)) {
                                            gene_group_df <- signature.gene_group_df.all
                                        } else if (! is.null(genes_of_interest.input) & ! is.null(signature_genes.input)) {
                                            gene_group_df <- rbind(query.gene_group_df.all, signature.gene_group_df.all)
                                        }
                                        
                                        if (! is.null(gene_group_df)) {
                                            write.csv(gene_group_df, file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',differential.prefolder,'/',project,'.',gsub("\\.","_",subdata.folder),'.',gsub("\\.","_",subpopulation.folder),'.',redDim.method[redDim],'.',file.prename,'.','differential','.','geneModule','.csv', sep = ''), row.names = TRUE, quote = FALSE)
                                        }
                                    }
                                    
                                    mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',differential.prefolder, sep = '')
                                    if (file.exists(file.path(mainDir))) {
                                        if (length(list.files(path = mainDir, full.names = TRUE, recursive = FALSE)) == 0) {
                                            unlink(mainDir, recursive = TRUE)
                                        }
                                    }
                                    
                                    mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',differential.folder, sep = '')
                                    if (file.exists(file.path(mainDir))) {
                                        if (length(list.files(path = mainDir, full.names = TRUE, recursive = FALSE)) == 0) {
                                            unlink(mainDir, recursive = TRUE)
                                        }
                                    }
                                    
                                    mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder, sep = '')
                                    if (file.exists(file.path(mainDir))) {
                                        if (length(list.files(path = mainDir, full.names = TRUE, recursive = FALSE)) == 0) {
                                            unlink(mainDir, recursive = TRUE)
                                        }
                                    }
                                    #}
                                    
                                    
                                    # 2. The second way of looking at modules and their expression is to pass gene_group_df directly to plot_cells(). 
                                    # If there are many modules, it can be hard to see where each one is expressed, so we'll just look at a subset of them:
                                    if (! is.null(gene_group_df)) {
                                        
                                        #module.type = c("module", "supermodule")
                                        module.type = c("module")
                                        
                                        for (M in 1:length(module.type)) {
                                            
                                            # /// marker.genes ///
                                            module.levels <- levels(unlist(gene_group_df[,module.type[M]]))
                                            
                                            module_set.list <- list()
                                            
                                            for (LV in 1:length(module.levels)) {
                                                module.levels.LV <- list(module.levels[LV])
                                                if (module.levels[LV] %in% c(genes_of_interest.titles, signature_genes.titles)) {
                                                    names(module.levels.LV) = module.levels[LV]
                                                } else {
                                                    names(module.levels.LV) = paste(module.type[M],module.levels[LV], sep = '')
                                                }
                                                module_set.list <- c(module_set.list, module.levels.LV)
                                            }
                                            
                                            # *** query gene(s) ***
                                            if (! is.null(genes_of_interest.input)) {
                                                for (gX in 1:length(genes_of_interest.input)) {
                                                    genes_of_interest.genes <- genes_of_interest.input[[gX]]
                                                    genes_of_interest.genes.idx <- which(sapply(genes_of_interest.genes, FUN = function(X) X %in% unlist(gene_group_df[,"gene_short_name"])))
                                                    if (length(genes_of_interest.genes.idx) > 0) {
                                                        genes_of_interest.genes = genes_of_interest.genes[genes_of_interest.genes.idx]
                                                        
                                                        genes_of_interest.genes.idx <- lapply(genes_of_interest.genes, FUN = function(X) which(unlist(gene_group_df[,"gene_short_name"]) %in% X))
                                                        names(genes_of_interest.genes.idx) = genes_of_interest.genes
                                                        
                                                        marker.genes = gene_group_df[unlist(genes_of_interest.genes.idx),]
                                                        
                                                        marker.modules <- marker.genes[,module.type[M]]
                                                        genes_of_interest.modules <- as.character(unlist(unique(marker.modules)))
                                                        genes_of_interest.modules <- genes_of_interest.modules[which(! genes_of_interest.modules %in% genes_of_interest.input.title[which(! genes_of_interest.input.title %in% names(genes_of_interest.input[gX]))])]
                                                        if (! is.null(genes_of_interest.modules)) {
                                                            a = grep(paste(genes_of_interest.titles, collapse = '|'), genes_of_interest.modules, value = TRUE, invert = TRUE)
                                                            a.idx = grep(paste(genes_of_interest.titles, collapse = '|'), genes_of_interest.modules, value = FALSE, invert = TRUE)
                                                            if (length(grep(paste(genes_of_interest.titles, collapse = '|'), genes_of_interest.modules, value = TRUE, invert = FALSE)) > 0) {
                                                                b.genes_of_interest = grep(paste(genes_of_interest.titles, collapse = '|'), genes_of_interest.modules, value = TRUE, invert = FALSE)
                                                                b.genes_of_interest.idx = grep(paste(genes_of_interest.titles, collapse = '|'), genes_of_interest.modules, value = FALSE, invert = FALSE)
                                                            } else {
                                                                b.genes_of_interest = NULL
                                                                b.genes_of_interest.idx = NULL
                                                            }
                                                            genes_of_interest.modules = c(a[mixedorder(a, decreasing = FALSE)], 
                                                                                          names(which(sapply(genes_of_interest.titles, FUN = function(X) X %in% b.genes_of_interest))))
                                                            
                                                            for (gY in 1:length(genes_of_interest.modules)) {
                                                                genes_of_interest.modules.gX.gY <- list(genes_of_interest.modules[gY])
                                                                if (genes_of_interest.modules[gY] %in% genes_of_interest.titles) {
                                                                    names(genes_of_interest.modules.gX.gY) = paste(names(genes_of_interest.input[gX]),'.',genes_of_interest.modules[gY], sep = '')
                                                                } else {
                                                                    names(genes_of_interest.modules.gX.gY) = paste(names(genes_of_interest.input[gX]),'.',module.type[M],genes_of_interest.modules[gY], sep = '')
                                                                }
                                                                module_set.list <- c(module_set.list, genes_of_interest.modules.gX.gY)
                                                            }
                                                            genes_of_interest.modules <- list(unique(genes_of_interest.modules))
                                                            #if (length(unlist(genes_of_interest.modules)) == 1) {
                                                            #    names(genes_of_interest.modules) = paste("query_",module.type[M],"_of_interest", sep = '')
                                                            #} else if (length(unlist(genes_of_interest.modules)) > 1) {
                                                            #    names(genes_of_interest.modules) = paste("query_",module.type[M],"s","_of_interest", sep = '')
                                                            #}
                                                            names(genes_of_interest.modules) = paste(names(genes_of_interest.input[gX]),'.',"ALL", sep = '')
                                                            module_set.list <- c(module_set.list, genes_of_interest.modules)
                                                        }
                                                    }
                                                }
                                            }
                                            
                                            # *** signature gene(s) ***
                                            if (! is.null(signature_genes.input)) {
                                                for (gX in 1:length(signature_genes.input)) {
                                                    signature_genes.genes <- signature_genes.input[[gX]]
                                                    signature_genes.genes.idx <- which(sapply(signature_genes.genes, FUN = function(X) X %in% unlist(gene_group_df[,"gene_short_name"])))
                                                    if (length(signature_genes.genes.idx) > 0) {
                                                        signature_genes.genes = signature_genes.genes[signature_genes.genes.idx]
                                                        
                                                        signature_genes.genes.idx <- lapply(signature_genes.genes, FUN = function(X) which(unlist(gene_group_df[,"gene_short_name"]) %in% X))
                                                        names(signature_genes.genes.idx) = signature_genes.genes
                                                        
                                                        marker.genes = gene_group_df[unlist(signature_genes.genes.idx),]
                                                        
                                                        marker.modules <- marker.genes[,module.type[M]]
                                                        signature_genes.modules <- as.character(unlist(unique(marker.modules)))
                                                        signature_genes.modules <- signature_genes.modules[which(! signature_genes.modules %in% signature_genes.input.title[which(! signature_genes.input.title %in% names(signature_genes.input[gX]))])]
                                                        if (! is.null(signature_genes.modules)) {
                                                            a = grep(paste(signature_genes.titles, collapse = '|'), signature_genes.modules, value = TRUE, invert = TRUE)
                                                            a.idx = grep(paste(signature_genes.titles, collapse = '|'), signature_genes.modules, value = FALSE, invert = TRUE)
                                                            if (length(grep(paste(signature_genes.titles, collapse = '|'), signature_genes.modules, value = TRUE, invert = FALSE)) > 0) {
                                                                b.signature_genes = grep(paste(signature_genes.titles, collapse = '|'), signature_genes.modules, value = TRUE, invert = FALSE)
                                                                b.signature_genes.idx = grep(paste(signature_genes.titles, collapse = '|'), signature_genes.modules, value = FALSE, invert = FALSE)
                                                            } else {
                                                                b.signature_genes = NULL
                                                                b.signature_genes.idx = NULL
                                                            }
                                                            signature_genes.modules = c(a[mixedorder(a, decreasing = FALSE)], 
                                                                                        names(which(sapply(signature_genes.titles, FUN = function(X) X %in% b.signature_genes))))
                                                            
                                                            for (gY in 1:length(signature_genes.modules)) {
                                                                signature_genes.modules.gX.gY <- list(signature_genes.modules[gY])
                                                                if (signature_genes.modules[gY] %in% signature_genes.titles) {
                                                                    names(signature_genes.modules.gX.gY) = paste(names(signature_genes.input[gX]),'.',signature_genes.modules[gY], sep = '')
                                                                } else {
                                                                    names(signature_genes.modules.gX.gY) = paste(names(signature_genes.input[gX]),'.',module.type[M],signature_genes.modules[gY], sep = '')
                                                                }
                                                                module_set.list <- c(module_set.list, signature_genes.modules.gX.gY)
                                                            }
                                                            signature_genes.modules <- list(unique(signature_genes.modules))
                                                            #if (length(unlist(signature_genes.modules)) == 1) {
                                                            #    names(signature_genes.modules) = paste("signature_",module.type[M],"_of_interest", sep = '')
                                                            #} else if (length(unlist(signature_genes.modules)) > 1) {
                                                            #    names(signature_genes.modules) = paste("signature_",module.type[M],"s","_of_interest", sep = '')
                                                            #}
                                                            names(signature_genes.modules) = paste(names(signature_genes.input[gX]),'.',"ALL", sep = '')
                                                            module_set.list <- c(module_set.list, signature_genes.modules)
                                                        }
                                                    }
                                                }
                                            }
                                            
                                            if (module.type[M] == "module") {
                                                example_modules <- list(unique(c(8, 28, 33, 37)))
                                                if (length(unlist(example_modules)) == 1) {
                                                    names(example_modules) = "example_module"
                                                } else if (length(unlist(example_modules)) > 1) {
                                                    names(example_modules) = "example_modules"
                                                }
                                                module_set.list <- c(module_set.list, example_modules)
                                            }
                                            
                                            if (length(module_set.list) > 0) {
                                                for (N in 1:length(module_set.list)) {
                                                    #print(names(module_set.list[N]))
                                                    
                                                    marker.modules <- module_set.list[[N]]
                                                    marker.modules.idx <- which(sapply(marker.modules, FUN = function(X) X %in% unlist(gene_group_df[,module.type[M]])))
                                                    if (length(marker.modules.idx) > 0) {
                                                        marker.modules = marker.modules[marker.modules.idx]
                                                        
                                                        marker.modules.idx <- lapply(marker.modules, FUN = function(X) which(unlist(gene_group_df[,module.type[M]]) %in% X))
                                                        names(marker.modules.idx) = marker.modules
                                                        
                                                        marker.genes = gene_group_df[unlist(marker.modules.idx),]
                                                        
                                                        gene_id.dup.idx <- which(duplicated(marker.genes[,"id"]))
                                                        if (length(gene_id.dup.idx) > 0) {
                                                            rm.idx.all <- c()
                                                            for (d in 1:length(gene_id.dup.idx)) {
                                                                gene_id.idx <- which(marker.genes[,"id"] == unlist(marker.genes[gene_id.dup.idx[d],"id"]))
                                                                
                                                                gene_id <- marker.genes[gene_id.idx,"id"]
                                                                gene_module <- marker.genes[gene_id.idx,"module"]
                                                                #gene_supermodule <- marker.genes[gene_id.idx,"supermodule"]
                                                                #gene_name <- marker.genes[gene_id.idx,"gene_short_name"]
                                                                
                                                                gene_id <- unlist(gene_id)
                                                                gene_module <- unlist(gene_module)
                                                                #gene_supermodule <- unlist(gene_supermodule)
                                                                #gene_name <- unlist(gene_name)
                                                                
                                                                gene_module.idx <- which(gene_module %in% c(genes_of_interest.titles, signature_genes.titles))
                                                                #gene_supermodule.idx <- which(gene_supermodule %in% c(genes_of_interest.titles, signature_genes.titles))
                                                                
                                                                rm.idx <- gene_id.idx[gene_module.idx]
                                                                
                                                                rm.idx.all <- c(rm.idx.all, rm.idx)
                                                            }
                                                            rm.idx.all <- unique(rm.idx.all)
                                                            #print(rm.idx.all)
                                                            
                                                            if (length(rm.idx.all) > 0) {
                                                                marker.genes <- marker.genes[-rm.idx.all, , drop = FALSE]
                                                            }
                                                        }
                                                        
                                                        if (module.type[M] == "module") {
                                                            marker.modules <- marker.genes[,"module"]
                                                        } else if (module.type[M] == "supermodule") {
                                                            marker.modules <- marker.genes[,"module"]
                                                        }
                                                        marker.modules.elements <- as.character(unlist(marker.modules))
                                                        if (any(module_set.list[[N]] %in% c("genes_of_interest", "signature_genes"))) {
                                                            marker.modules.elements <- marker.modules.elements
                                                        } else {
                                                            marker.modules.elements <- marker.modules.elements[which(! marker.modules.elements %in% c(genes_of_interest.input.title, signature_genes.input.title)[which(! c(genes_of_interest.input.title, signature_genes.input.title) %in% unlist(strsplit(names(module_set.list[N]), '\\.'))[1])])]
                                                        }
                                                        if (! is.null(marker.modules.elements)) {
                                                            if (! is.null(genes_of_interest.input)) {
                                                                for (gX in 1:length(genes_of_interest.input)) {
                                                                    marker.modules.elements[which(marker.modules.elements == paste("geneSet",gX,"_of_interest", sep = ''))] = paste("geneSet",gX, sep = '')
                                                                }
                                                            }
                                                            marker.modules.elements[which(marker.modules.elements == "SINGH_KRAS_DEPENDENCY_SIGNATURE")] = "SINGH"
                                                            marker.modules.elements[which(marker.modules.elements == "FRIDMAN_SENESCENCE_UP")] = "FRIDMAN"
                                                            marker.modules.elements[which(marker.modules.elements == "HALLMARK_KRAS_SIGNALING_UP")] = "HALLMARK"
                                                            
                                                            if (module.type[M] == "module") {
                                                                marker.genes[,"module"] = factor(marker.modules.elements, levels = unique(marker.modules.elements))
                                                            } else if (module.type[M] == "supermodule") {
                                                                marker.genes[,"module"] = factor(marker.modules.elements, levels = unique(marker.modules.elements))
                                                            }
                                                        }
                                                        
                                                        if (module.type[M] == "module") {
                                                            panel.num = length(unlist(unique(marker.genes[,"module"])))
                                                        } else if (module.type[M] == "supermodule") {
                                                            panel.num = length(unlist(unique(marker.genes[,"module"])))
                                                        }
                                                        if (panel.num == 1) {
                                                            width.size = 7
                                                            height.size = 6
                                                        } else if (panel.num == 2) {
                                                            width.size = 7 + (5 * 1)
                                                            height.size = 6
                                                        } else if (panel.num == 3) {
                                                            width.size = 7 + (5 * 2)
                                                            height.size = 6
                                                        } else if (panel.num > 3) {
                                                            width.size = 7 + (5 * 3)
                                                            height.size = width.size - 1
                                                        }
                                                        
                                                        #for (redDim in 1:length(redDim.method)) {
                                                        
                                                        # /// 2d ///
                                                        dim.folder = paste(eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))),'D', sep = '')
                                                        dim.finalX = eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) - (2-1)
                                                        
                                                        if (dim.finalX >= 1) {
                                                            
                                                            for (dim.comb in 1:dim.finalX) {
                                                                dim.x <- dim.comb
                                                                dim.y <- dim.x + 1
                                                                
                                                                if (dim.finalX > 1) {
                                                                    if (dim.x == 1) {
                                                                        comb.run = 2
                                                                    } else {
                                                                        comb.run = 1
                                                                    }
                                                                } else {
                                                                    comb.run = 1
                                                                }
                                                                
                                                                for (dim.comb.run in 1:comb.run) {
                                                                    if (dim.comb.run > 1) {
                                                                        dim.y <- dim.y + 1
                                                                    }
                                                                    
                                                                    if (eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) == 2) {
                                                                        dim.xy <- paste('D',dim.x,'-','D',dim.y, sep = '')
                                                                    } else {
                                                                        dim.xy <- paste('D',dim.x,'-','D',dim.y, sep = '')
                                                                    }
                                                                    #print(dim.xy)
                                                                    
                                                                    if (SUB == 0) {
                                                                        file.prefolder = dim.xy
                                                                    } else {
                                                                        file.prefolder = paste(dim.xy,'/',SUB.prefolder, sep = '')
                                                                    }
                                                                    
                                                                    mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder, sep = '')
                                                                    subDir <- 'differential'
                                                                    if (! file.exists(file.path(mainDir, subDir))) {
                                                                        dir.create(file.path(mainDir, subDir))
                                                                    }
                                                                    mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential', sep = '')
                                                                    subDir <- module.type[M]
                                                                    if (! file.exists(file.path(mainDir, subDir))) {
                                                                        dir.create(file.path(mainDir, subDir))
                                                                    }
                                                                    
                                                                    if (length(grep(paste(genes_of_interest.titles, collapse = '|'), names(module_set.list[N]), value = TRUE, invert = FALSE)) > 0) {
                                                                        mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential', sep = '')
                                                                        subDir <- 'genes_of_interest'
                                                                        if (! file.exists(file.path(mainDir, subDir))) {
                                                                            dir.create(file.path(mainDir, subDir))
                                                                        }
                                                                        mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential','/','genes_of_interest', sep = '')
                                                                        subDir <- unlist(strsplit(names(module_set.list[N]), '\\.'))[1]
                                                                        if (! file.exists(file.path(mainDir, subDir))) {
                                                                            dir.create(file.path(mainDir, subDir))
                                                                        }
                                                                        mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential','/','genes_of_interest','/',unlist(strsplit(names(module_set.list[N]), '\\.'))[1], sep = '')
                                                                        subDir <- module.type[M]
                                                                        if (! file.exists(file.path(mainDir, subDir))) {
                                                                            dir.create(file.path(mainDir, subDir))
                                                                        }
                                                                        
                                                                        file.folder = paste('genes_of_interest','/',unlist(strsplit(names(module_set.list[N]), '\\.'))[1],'/',module.type[M], sep = '')
                                                                        if (names(module_set.list[N]) %in% genes_of_interest.titles) {
                                                                            file.name.suffix = names(module_set.list[N])
                                                                        } else {
                                                                            file.name.suffix = gsub("\\.","_",unlist(strsplit(names(module_set.list[N]), '\\.'))[2])
                                                                        }
                                                                    } else if (length(grep(paste(signature_genes.titles, collapse = '|'), names(module_set.list[N]), value = TRUE, invert = FALSE)) > 0) {
                                                                        mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential', sep = '')
                                                                        subDir <- 'signature_genes'
                                                                        if (! file.exists(file.path(mainDir, subDir))) {
                                                                            dir.create(file.path(mainDir, subDir))
                                                                        }
                                                                        mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential','/','signature_genes', sep = '')
                                                                        subDir <- unlist(strsplit(names(module_set.list[N]), '\\.'))[1]
                                                                        if (! file.exists(file.path(mainDir, subDir))) {
                                                                            dir.create(file.path(mainDir, subDir))
                                                                        }
                                                                        mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential','/','signature_genes','/',unlist(strsplit(names(module_set.list[N]), '\\.'))[1], sep = '')
                                                                        subDir <- module.type[M]
                                                                        if (! file.exists(file.path(mainDir, subDir))) {
                                                                            dir.create(file.path(mainDir, subDir))
                                                                        }
                                                                        
                                                                        file.folder = paste('signature_genes','/',unlist(strsplit(names(module_set.list[N]), '\\.'))[1],'/',module.type[M], sep = '')
                                                                        if (names(module_set.list[N]) %in% signature_genes.titles) {
                                                                            file.name.suffix = names(module_set.list[N])
                                                                        } else {
                                                                            file.name.suffix = gsub("\\.","_",unlist(strsplit(names(module_set.list[N]), '\\.'))[2])
                                                                        }
                                                                    } else {
                                                                        file.folder = module.type[M]
                                                                        file.name.suffix = gsub("\\.","_",names(module_set.list[N]))
                                                                    }
                                                                    
                                                                    if (is.null(file.name.suffix)) {
                                                                        file.name = paste(project,'.',file.prename,'.','plot_cells','.',redDim.method[redDim],'_',dim.xy,'.','differential', sep = '')
                                                                    } else {
                                                                        if (is.na(file.name.suffix) | file.name.suffix == "") {
                                                                            file.name = paste(project,'.',file.prename,'.','plot_cells','.',redDim.method[redDim],'_',dim.xy,'.','differential', sep = '')
                                                                        } else {
                                                                            file.name = paste(project,'.',file.prename,'.','plot_cells','.',redDim.method[redDim],'_',dim.xy,'.','differential','.',file.name.suffix, sep = '')
                                                                        }
                                                                    }
                                                                    
                                                                    png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential','/',file.folder,'/',file.name,'.png', sep = ''), width = width.size, height = height.size, units = 'in', res = 600)
                                                                    p <- plot_cells(cds_subpopulation, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, 
                                                                                    genes = marker.genes, norm_method = "log", 
                                                                                    color_cells_by = "partition", group_cells_by = "partition", 
                                                                                    label_cell_groups = FALSE, 
                                                                                    show_trajectory_graph = FALSE)
                                                                    print(p)
                                                                    dev.off()
                                                                }
                                                            }
                                                        }
                                                        
                                                        # /// 3d ///
                                                        if (module.type[M] == "module") {
                                                            
                                                            genes_of_interest.input.title.idx <- which(module.levels %in% genes_of_interest.titles)
                                                            signature_genes.input.title.idx <- which(module.levels %in% signature_genes.titles)
                                                            
                                                            if ((N <= 3) | (N %in% c(genes_of_interest.input.title.idx, signature_genes.input.title.idx)) | (N > length(module.levels))) { # keep at most 3 modules (not save more than 3 modules) for saving disk space, please modify or un-comment this condition if needed
                                                                
                                                                if (eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) == 2) {
                                                                    dim.folder = paste(eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))),'D', sep = '')
                                                                    #dim.finalX = eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) - (3-1)
                                                                    #dim.folder = paste(3,'D', sep = '')
                                                                    dim.finalX = 3 - (3-1)
                                                                } else {
                                                                    dim.folder = paste(eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))),'D', sep = '')
                                                                    dim.finalX = eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) - (3-1)
                                                                }
                                                                
                                                                if (dim.finalX >= 1) {
                                                                    
                                                                    for (dim.comb in 1:dim.finalX) {
                                                                        dim.x <- dim.comb
                                                                        dim.y <- dim.x + 1
                                                                        dim.z <- dim.y + 1
                                                                        
                                                                        if (eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) == 2) {
                                                                            dim.xyz <- paste('D',dim.x,'-','D',dim.y,'-','D',dim.z, sep = '')
                                                                        } else {
                                                                            dim.xyz <- paste('D',dim.x,'-','D',dim.y,'-','D',dim.z, sep = '')
                                                                        }
                                                                        #print(dim.xyz)
                                                                        
                                                                        if (SUB == 0) {
                                                                            file.prefolder = dim.xyz
                                                                        } else {
                                                                            file.prefolder = paste(dim.xyz,'/',SUB.prefolder, sep = '')
                                                                        }
                                                                        
                                                                        mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder, sep = '')
                                                                        subDir <- 'differential'
                                                                        if (! file.exists(file.path(mainDir, subDir))) {
                                                                            dir.create(file.path(mainDir, subDir))
                                                                        }
                                                                        mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential', sep = '')
                                                                        subDir <- module.type[M]
                                                                        if (! file.exists(file.path(mainDir, subDir))) {
                                                                            dir.create(file.path(mainDir, subDir))
                                                                        }
                                                                        
                                                                        if (length(grep(paste(genes_of_interest.titles, collapse = '|'), names(module_set.list[N]), value = TRUE, invert = FALSE)) > 0) {
                                                                            mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential', sep = '')
                                                                            subDir <- 'genes_of_interest'
                                                                            if (! file.exists(file.path(mainDir, subDir))) {
                                                                                dir.create(file.path(mainDir, subDir))
                                                                            }
                                                                            mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential','/','genes_of_interest', sep = '')
                                                                            subDir <- unlist(strsplit(names(module_set.list[N]), '\\.'))[1]
                                                                            if (! file.exists(file.path(mainDir, subDir))) {
                                                                                dir.create(file.path(mainDir, subDir))
                                                                            }
                                                                            mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential','/','genes_of_interest','/',unlist(strsplit(names(module_set.list[N]), '\\.'))[1], sep = '')
                                                                            subDir <- module.type[M]
                                                                            if (! file.exists(file.path(mainDir, subDir))) {
                                                                                dir.create(file.path(mainDir, subDir))
                                                                            }
                                                                            
                                                                            file.folder = paste('genes_of_interest','/',unlist(strsplit(names(module_set.list[N]), '\\.'))[1],'/',module.type[M], sep = '')
                                                                            if (names(module_set.list[N]) %in% genes_of_interest.titles) {
                                                                                file.name.suffix = names(module_set.list[N])
                                                                            } else {
                                                                                file.name.suffix = gsub("\\.","_",unlist(strsplit(names(module_set.list[N]), '\\.'))[2])
                                                                            }
                                                                        } else if (length(grep(paste(signature_genes.titles, collapse = '|'), names(module_set.list[N]), value = TRUE, invert = FALSE)) > 0) {
                                                                            mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential', sep = '')
                                                                            subDir <- 'signature_genes'
                                                                            if (! file.exists(file.path(mainDir, subDir))) {
                                                                                dir.create(file.path(mainDir, subDir))
                                                                            }
                                                                            mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential','/','signature_genes', sep = '')
                                                                            subDir <- unlist(strsplit(names(module_set.list[N]), '\\.'))[1]
                                                                            if (! file.exists(file.path(mainDir, subDir))) {
                                                                                dir.create(file.path(mainDir, subDir))
                                                                            }
                                                                            mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential','/','signature_genes','/',unlist(strsplit(names(module_set.list[N]), '\\.'))[1], sep = '')
                                                                            subDir <- module.type[M]
                                                                            if (! file.exists(file.path(mainDir, subDir))) {
                                                                                dir.create(file.path(mainDir, subDir))
                                                                            }
                                                                            
                                                                            file.folder = paste('signature_genes','/',unlist(strsplit(names(module_set.list[N]), '\\.'))[1],'/',module.type[M], sep = '')
                                                                            if (names(module_set.list[N]) %in% signature_genes.titles) {
                                                                                file.name.suffix = names(module_set.list[N])
                                                                            } else {
                                                                                file.name.suffix = gsub("\\.","_",unlist(strsplit(names(module_set.list[N]), '\\.'))[2])
                                                                            }
                                                                        } else {
                                                                            file.folder = module.type[M]
                                                                            file.name.suffix = gsub("\\.","_",names(module_set.list[N]))
                                                                        }
                                                                        
                                                                        #if (redDim.method[redDim] == "UMAP") {
                                                                        #    trj.plt = 1
                                                                        #} else if (redDim.method[redDim] == "tSNE") {
                                                                        trj.plt = 0
                                                                        #}
                                                                        if (graph.run == 0) {
                                                                            trj.plt = 0
                                                                        }
                                                                        
                                                                        for (trj in 0:trj.plt) {
                                                                            if (trj == 0) {
                                                                                trj.prt.endidx = 1
                                                                            } else if (trj == 1) {
                                                                                trj.prt.endidx = 0
                                                                            }
                                                                            
                                                                            #for (trj.prt in 1:1) {
                                                                            for (trj.prt in 1:trj.prt.endidx) {
                                                                                if (trj.prt == 1) {
                                                                                    trj.prt.name = 'with_partition'
                                                                                    cds_subpopulation_3d.plt = cds_subpopulation_3d
                                                                                } else if (trj.prt == 0) {
                                                                                    #trj.prt.name = 'without_partition'
                                                                                    #cds_subpopulation_3d.plt = cds_subpopulation_3d.no_partition
                                                                                }
                                                                                
                                                                                if (! is.null(cds_subpopulation_3d.plt)) {
                                                                                    
                                                                                    for (grid.plt in 1:0) {
                                                                                        if (grid.plt == 1) {
                                                                                            grid.name = 'with_grid'
                                                                                        } else if (grid.plt == 0) {
                                                                                            grid.name = 'without_grid'
                                                                                        }
                                                                                        
                                                                                        if (trj == 0) {
                                                                                            show_trajectory = FALSE
                                                                                            if (is.null(file.name.suffix)) {
                                                                                                file.name = paste(project,'.',file.prename,'.','plot_cells_3d','.',redDim.method[redDim],'_',dim.xyz,'.','differential','.','without_trajectory', sep = '')
                                                                                            } else {
                                                                                                if (is.na(file.name.suffix) | file.name.suffix == "") {
                                                                                                    file.name = paste(project,'.',file.prename,'.','plot_cells_3d','.',redDim.method[redDim],'_',dim.xyz,'.','differential','.','without_trajectory', sep = '')
                                                                                                } else {
                                                                                                    file.name = paste(project,'.',file.prename,'.','plot_cells_3d','.',redDim.method[redDim],'_',dim.xyz,'.','differential','.',file.name.suffix,'.','without_trajectory', sep = '')
                                                                                                }
                                                                                            }
                                                                                            file.name = paste(file.name,'.',grid.name, sep = '')
                                                                                        } else if (trj == 1) {
                                                                                            show_trajectory = TRUE
                                                                                            if (is.null(file.name.suffix)) {
                                                                                                file.name = paste(project,'.',file.prename,'.','plot_cells_3d','.',redDim.method[redDim],'_',dim.xyz,'.','differential','.','with_trajectory', sep = '')
                                                                                            } else {
                                                                                                if (is.na(file.name.suffix) | file.name.suffix == "") {
                                                                                                    file.name = paste(project,'.',file.prename,'.','plot_cells_3d','.',redDim.method[redDim],'_',dim.xyz,'.','differential','.','with_trajectory', sep = '')
                                                                                                } else {
                                                                                                    file.name = paste(project,'.',file.prename,'.','plot_cells_3d','.',redDim.method[redDim],'_',dim.xyz,'.','differential','.',file.name.suffix,'.','with_trajectory', sep = '')
                                                                                                }
                                                                                            }
                                                                                            #file.name = paste(file.name,'.',trj.prt.name, sep = '')
                                                                                            file.name = paste(file.name,'.',grid.name, sep = '')
                                                                                        }
                                                                                        
                                                                                        file.html = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential','/',file.folder,'/',file.name,'.html', sep = '')
                                                                                        file.png = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential','/',file.folder,'/',file.name,'.png', sep = '')
                                                                                        
                                                                                        cds_subpopulation_3d_plot_obj <- plot_cells_3d(cds_subpopulation_3d.plt, reduction_method = redDim.method[redDim], dims = c(dim.x, dim.y, dim.z), cell_size = cell_size_3d, genes = marker.genes, norm_method = "log", show_trajectory_graph = show_trajectory)
                                                                                        
                                                                                        if (grid.plt == 0) {
                                                                                            cds_subpopulation_3d_plot_obj <- cds_subpopulation_3d_plot_obj %>% plotly::layout(scene = scene)
                                                                                        }
                                                                                        
                                                                                        # https://community.rstudio.com/t/save-viewer-object-rendered-in-rstudio-as-image/32796
                                                                                        if (dim.comb == 1) { # only save D1-D2-D3 (no D2-D3-D4, D3-D4-D5, ...)
                                                                                            if (length(unlist(strsplit(names(module_set.list[N]), '\\.'))) == 1) { # keep part of modules (not save single module) for saving disk space, please modify or un-comment this condition if needed
                                                                                                saveWidget(cds_subpopulation_3d_plot_obj, file.html)
                                                                                                #webshot(file.html, file.png)
                                                                                            } else {
                                                                                                # *** query gene(s) ***
                                                                                                if (unlist(strsplit(names(module_set.list[N]), '\\.'))[1] %in% genes_of_interest.titles) {
                                                                                                    saveWidget(cds_subpopulation_3d_plot_obj, file.html)
                                                                                                    #webshot(file.html, file.png)
                                                                                                }
                                                                                                if (unlist(strsplit(names(module_set.list[N]), '\\.'))[2] %in% c(genes_of_interest.titles, "ALL")) {
                                                                                                    saveWidget(cds_subpopulation_3d_plot_obj, file.html)
                                                                                                    #webshot(file.html, file.png)
                                                                                                }
                                                                                                # *** signature gene(s) ***
                                                                                                if (unlist(strsplit(names(module_set.list[N]), '\\.'))[2] %in% c(signature_genes.titles, "ALL")) {
                                                                                                    saveWidget(cds_subpopulation_3d_plot_obj, file.html)
                                                                                                    #webshot(file.html, file.png)
                                                                                                }
                                                                                            }
                                                                                        } else { # remove D2-D3-D4, D3-D4-D5, ...
                                                                                            mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','clustering','/',subdata.folder,'/',subpopulation.folder,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential', sep = '')
                                                                                            if (file.exists(file.path(mainDir))) {
                                                                                                unlink(mainDir, recursive = TRUE)
                                                                                            }
                                                                                        }
                                                                                    }
                                                                                }
                                                                            }
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                        #}
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                                
                                cds_subpopulation = cds_subpopulation.origin.processed
                                cds_subpopulation_3d = cds_subpopulation_3d.origin.processed
                            }
                        }
                    }
                }
                
                cds_subpopulation.clustering <- cds_subpopulation.origin.processed
                #get_citations(cds_subpopulation.clustering)
                
                cds_subpopulation_3d.clustering <- cds_subpopulation_3d.origin.processed
                #get_citations(cds_subpopulation_3d.clustering)
                
            }
        }
        
        
        if ("trajectories" %in% analysis_type.input) {
            
            # === Constructing single-cell trajectories =========================================================================================================
            # https://cole-trapnell-lab.github.io/monocle3/docs/trajectories/
            
            # Analysis of single cell RNA-seq data: 11. Trajectory inference
            # https://biocellgen-public.svi.edu.au/mig_2019_scrnaseq-workshop/public/trajectory-inference.html#monocle-3
            
            mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name, sep = '')
            subDir <- 'trajectories'
            if (! file.exists(file.path(mainDir, subDir))) {
                dir.create(file.path(mainDir, subDir))
            }
            
            mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories', sep = '')
            subDir <- subdata.folder
            if (! file.exists(file.path(mainDir, subDir))) {
                dir.create(file.path(mainDir, subDir))
            }
            
            mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder, sep = '')
            subDir <- project.subname
            if (! file.exists(file.path(mainDir, subDir))) {
                dir.create(file.path(mainDir, subDir))
            }
            
            if (project.name == "tutorial") {
                batch = "batch"
                cell_type = "cell.type"
            } else {
                batch = batch.input
                cell_type = cell_type.input
            }
            
            if (! is.null(batch)) {
                batch.name = paste('.',gsub("\\.","_",batch),'_','align','.', sep = '')
            } else {
                batch.name = '.'
            }
            
            
            if (project.name == "tutorial") {
                #expression_matrix <- readRDS(url("http://staff.washington.edu/hpliner/data/packer_embryo_expression.rds"))
                #cell_metadata <- readRDS(url("http://staff.washington.edu/hpliner/data/packer_embryo_colData.rds"))
                #gene_metadata <- readRDS(url("http://staff.washington.edu/hpliner/data/packer_embryo_rowData.rds"))
                expression_matrix <- readRDS(paste(current_folder,'/Script/function','/','Monocle3','/','tutorial','/',"packer_embryo_expression.rds", sep = ''))
                cell_metadata <- readRDS(paste(current_folder,'/Script/function','/','Monocle3','/','tutorial','/',"packer_embryo_colData.rds", sep = ''))
                gene_metadata <- readRDS(paste(current_folder,'/Script/function','/','Monocle3','/','tutorial','/',"packer_embryo_rowData.rds", sep = ''))
            } else {
                if (is.matrix(expression_matrix) | is.data.frame(expression_matrix)) {
                    expression_matrix <- Matrix(as.matrix(expression_matrix), sparse = TRUE)
                } else {
                    expression_matrix <- expression_matrix
                }
                cell_metadata <- cell_metadata
                gene_metadata <- gene_metadata
            }
            expression_matrix[is.na(expression_matrix)] <- 0
            #dim(expression_matrix)
            #dim(cell_metadata)
            #dim(gene_metadata)
            
            #colnames(cell_metadata)
            #colnames(gene_metadata)
            
            # --- max_components = maximum dimension ------------------------------------------------------------------------------------------------------------
            if (redDim_maxDim.input == 2 | redDim_maxDim.input == 3) {
                dim.run.idx1 = 1
                dim.run.idx2 = 2
            } else {
                if (redDim_minDim.input == redDim_maxDim.input) {
                    dim.run.idx1 = 3
                    dim.run.idx2 = 3
                } else {
                    dim.run.idx1 = 1
                    dim.run.idx2 = 3
                }
            }
            
            for (dim.run in dim.run.idx1:dim.run.idx2) {
                
                if (dim.run == 1) {
                    #reduce_dimension() # max_components: the dimensionality of the reduced space. Default is 2.
                    assign(paste('UMAP','.','dim_max', sep = ''), 2)
                    assign(paste('tSNE','.','dim_max', sep = ''), 2)
                } else if (dim.run == 2) {
                    assign(paste('UMAP','.','dim_max', sep = ''), 3)
                    assign(paste('tSNE','.','dim_max', sep = ''), 3)
                } else if (dim.run == 3) {
                    assign(paste('UMAP','.','dim_max', sep = ''), redDim_maxDim.input)
                    assign(paste('tSNE','.','dim_max', sep = ''), 3)
                }
                
                # Make the CDS object
                cds.origin <- new_cell_data_set(expression_matrix,
                                                cell_metadata = cell_metadata,
                                                gene_metadata = gene_metadata)
                
                if (subdata.idx == 0) {
                    cds <- cds.origin
                } else {
                    cds <- cds.origin
                    
                    #rowidx <- 
                    #cds <- cds.origin[rowidx, , drop = FALSE]
                }
                
                if (ncol(cds) > 0) {
                    
                    if (ncol(cds) <= 100) {
                        cell_size_2d = 1.5
                        cell_size_3d = 250
                    } else {
                        cell_size_2d = 0.35 # default for plot_cells()
                        cell_size_3d = 25 # default for plot_cells_3d()
                    }
                    
                    for (redDim in 1:length(redDim.method)) {
                        mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname, sep = '')
                        subDir <- redDim.method[redDim]
                        if (! file.exists(file.path(mainDir, subDir))) {
                            dir.create(file.path(mainDir, subDir))
                        }
                        
                        dim.folder = paste(eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))),'D', sep = '')
                        
                        mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim], sep = '')
                        subDir <- dim.folder
                        if (! file.exists(file.path(mainDir, subDir))) {
                            dir.create(file.path(mainDir, subDir))
                        }
                    }
                    
                    
                    ## Step 1: Normalize and pre-process the data *****************************************************
                    # Pre-process the data
                    # This time, we will use a different strategy for batch correction, which includes what Packer & Zhu et al did in their original analysis
                    # Note: Your data will not have the loading batch information demonstrated here, you will correct batch using your own batch information.
                    if (project.name %in% "GSM4150377") {
                        # Identify over-dispersed genes to be used as feature selection
                        dispersion_test = estimateDispersionsForCellDataSet(sciPlex_cds, min_cells_detected = 100, removeOutliers = TRUE)
                        
                        disp_table = dispersionTable(dispersion_test)
                        disp_table = disp_table %>% 
                            mutate(excess_disp = (dispersion_empirical - dispersion_fit) / dispersion_fit) %>% 
                            dplyr::arrange(plyr::desc(excess_disp))
                        
                        unsup_clustering_genes = as.character(head(disp_table, 1000)$gene_id)
                        ordering_genes = unsup_clustering_genes
                        
                        cds <- preprocess_cds(cds, method = "PCA", num_dim = 35, norm_method = "log", residualModelFormulaStr = "~log(n.umi)", use_genes = ordering_genes, verbose = FALSE)
                    } else if (project.name %in% "GSM4150378") {
                        cds <- preprocess_cds(cds, method = "PCA", num_dim = 25, norm_method = "log", verbose = FALSE)
                    } else {
                        if (ncol(cds) / 2 <= 50) {
                            # Warning in (function (A, nv = 5, nu = nv, maxit = 1000, work = nv + 7, reorth = TRUE,  :
                            # You're computing too large a percentage of total singular values, use a standard svd instead.
                            if (floor(ncol(cds) / 2) == ncol(cds) / 2) {
                                cds <- preprocess_cds(cds, method = "PCA", num_dim = floor(ncol(cds) / 2) - 1)
                            } else {
                                cds <- preprocess_cds(cds, method = "PCA", num_dim = floor(ncol(cds) / 2))
                            }
                        } else {
                            cds <- preprocess_cds(cds, method = "PCA", num_dim = 50) # method: Default is "PCA". method = c("PCA", "LSI") # "PCA": Principal Components Analysis (the standard for RNA-seq), "LSI": Latent Semantic Indexing (common in ATAC-seq)
                        }
                    }
                    
                    png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',project,'.','trajectories','.','plot_pc_variance_explained','.png', sep = ''), width = 6, height = 4, units = 'in', res = 600)
                    p <- plot_pc_variance_explained(cds)
                    print(p)
                    dev.off()
                    
                    
                    ## Step 2: Reduce the dimensions using UMAP *******************************************************
                    # Reduce dimensionality and visualize the results
                    # unlike clustering, which works well with both UMAP and t-SNE, here we strongly urge you to use UMAP, the default method
                    for (redDim in 1:length(redDim.method)) {
                        #reduce_dimension() # max_components: the dimensionality of the reduced space. Default is 2.
                        #umap() # n_components: The dimension of the space to embed into. This defaults to 2 to provide easy visualization, but can reasonably be set to any integer value in the range 2 to 100.
                        # UMAP: Monocle 3 uses UMAP by default, as we feel that it is both faster and better suited for clustering and trajectory analysis in RNA-seq.
                        # tSNE: Error in .check_tsne_params(nrow(X), dims = dims, perplexity = perplexity,  : dims should be either 1, 2 or 3
                        if (project.name %in% "GSM4150377") {
                            cds <- reduce_dimension(cds, reduction_method = redDim.method[redDim], max_components = eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))), umap.metric = "cosine", umap.n_neighbors = 50, umap.min_dist = 0.5, umap.fast_sgd = FALSE, cores = 1, verbose = FALSE)
                        } else if (project.name %in% "GSM4150378") {
                            cds <- reduce_dimension(cds, reduction_method = redDim.method[redDim], max_components = eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))), umap.metric = "cosine", umap.n_neighbors = 10, umap.min_dist = 0.1, verbose = FALSE)
                        } else {
                            cds <- reduce_dimension(cds, reduction_method = redDim.method[redDim], max_components = eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))))
                        }
                        
                        #cds@int_colData@listData$reducedDims@listData[[redDim.method[redDim]]]
                        
                        dim.folder = paste(eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))),'D', sep = '')
                        dim.finalX = eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) - (2-1)
                        
                        if (dim.finalX >= 1) {
                            
                            for (dim.comb in 1:dim.finalX) {
                                dim.x <- dim.comb
                                dim.y <- dim.x + 1
                                
                                if (dim.finalX > 1) {
                                    if (dim.x == 1) {
                                        comb.run = 2
                                    } else {
                                        comb.run = 1
                                    }
                                } else {
                                    comb.run = 1
                                }
                                
                                for (dim.comb.run in 1:comb.run) {
                                    if (dim.comb.run > 1) {
                                        dim.y <- dim.y + 1
                                    }
                                    
                                    if (eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) == 2) {
                                        dim.xy <- paste('D',dim.x,'-','D',dim.y, sep = '')
                                    } else {
                                        dim.xy <- paste('D',dim.x,'-','D',dim.y, sep = '')
                                    }
                                    #print(dim.xy)
                                    
                                    mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder, sep = '')
                                    subDir <- dim.xy
                                    if (! file.exists(file.path(mainDir, subDir))) {
                                        dir.create(file.path(mainDir, subDir))
                                    }
                                    
                                    # /// no color_cells ///
                                    png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',project,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,'.png', sep = ''), width = 4.075, height = 4, units = 'in', res = 600)
                                    p <- plot_cells(cds, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d)
                                    print(p)
                                    dev.off()
                                    
                                    # /// color_cells_by = cell_type ///
                                    if (length(cell_type) > 0) {
                                        for (color.plt in 1:length(cell_type)) {
                                            legend.max_nchar = max(nchar(levels(colData(cds)[, cell_type[color.plt]])))
                                            legend.multiply_factor = 0.5 + legend.max_nchar * 0.05
                                            legend.num = length(levels(colData(cds)[, cell_type[color.plt]]))
                                            legend.ncol = ceiling(legend.num / 10)
                                            width.size = 4.075 + (legend.multiply_factor * legend.ncol)
                                            
                                            png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',project,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,'.',gsub("\\.","_",cell_type[color.plt]),'.png', sep = ''), width = 4.075, height = 4, units = 'in', res = 600)
                                            p <- plot_cells(cds, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = cell_type[color.plt], label_groups_by_cluster = FALSE) # Finally, you can disable this labeling policy, making plot_cells() behave like it did before we called cluster_cells()
                                            if (! is.null(color_palette)) {
                                                cell_type.colidx <- which(paste(unique(c(cell_type[color.plt], unlist(strsplit(cell_type[color.plt], '\\.|_'))[1])),'.','color', sep = '') %in% names(color_palette))
                                                if (length(cell_type.colidx) == 1) {
                                                    p <- p + scale_color_manual(values = color_palette[[which(names(color_palette) %in% paste(unique(c(cell_type[color.plt], unlist(strsplit(cell_type[color.plt], '\\.|_'))[1])),'.','color', sep = ''))]])
                                                }
                                            }
                                            print(p)
                                            dev.off()
                                            
                                            png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',project,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,'.',gsub("\\.","_",cell_type[color.plt]),'_','legend','.png', sep = ''), width = width.size, height = 4, units = 'in', res = 600)
                                            p <- plot_cells(cds, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = cell_type[color.plt], label_groups_by_cluster = FALSE, label_cell_groups = FALSE) # Finally, you can disable this labeling policy, making plot_cells() behave like it did before we called cluster_cells()
                                            if (! is.null(color_palette)) {
                                                cell_type.colidx <- which(paste(unique(c(cell_type[color.plt], unlist(strsplit(cell_type[color.plt], '\\.|_'))[1])),'.','color', sep = '') %in% names(color_palette))
                                                if (length(cell_type.colidx) == 1) {
                                                    p <- p + scale_color_manual(values = color_palette[[which(names(color_palette) %in% paste(unique(c(cell_type[color.plt], unlist(strsplit(cell_type[color.plt], '\\.|_'))[1])),'.','color', sep = ''))]])
                                                }
                                            }
                                            p <- p + theme(legend.position = "right")
                                            p$guides$colour$ncol = legend.ncol
                                            if (project.name %in% c("GSE139944", c("GSM4150376", "GSM4150377", "GSM4150378", "GSM4150379"))) {
                                                if (cell_type[color.plt] == "pathway_level_1") {
                                                    p$guides$colour$title = "Pathway"
                                                } else if (cell_type[color.plt] == "product_name") {
                                                    p$guides$colour$title = "Treatment"
                                                } else if (cell_type[color.plt] == "dose_character") {
                                                    p$guides$colour$title = dose_character.unit
                                                } else {
                                                    p$guides$colour$title = "Treatment"
                                                }
                                            }
                                            print(p)
                                            dev.off()
                                        }
                                    }
                                }
                            }
                        }
                    }
                    
                    if (project.name == "tutorial") {
                        # /// marker.genes ///
                        gene_set.list <- list()
                        
                        ciliated_genes <- list(unique(c("che-1", "hlh-17", "nhr-6", "dmd-6", "ceh-36", "ham-1")))
                        if (length(unlist(ciliated_genes)) == 1) {
                            names(ciliated_genes) = "ciliated_gene"
                        } else if (length(unlist(ciliated_genes)) > 1) {
                            names(ciliated_genes) = "ciliated_genes"
                        }
                        gene_set.list <- c(gene_set.list, ciliated_genes)
                        
                        if (length(gene_set.list) > 0) {
                            for (N in 1:length(gene_set.list)) {
                                #print(names(gene_set.list[N]))
                                
                                marker.genes <- gene_set.list[[N]]
                                marker.genes.idx <- which(sapply(marker.genes, FUN = function(X) X %in% rowData(cds)[,"gene_short_name"]))
                                if (length(marker.genes.idx) > 0) {
                                    marker.genes = marker.genes[marker.genes.idx]
                                    
                                    panel.num = length(unlist(unique(marker.genes)))
                                    if (panel.num == 1) {
                                        width.size = 7
                                        height.size = 6
                                    } else if (panel.num == 2) {
                                        width.size = 7 + (5 * 1)
                                        height.size = 6
                                    } else if (panel.num == 3) {
                                        width.size = 7 + (5 * 2)
                                        height.size = 6
                                    } else if (panel.num > 3) {
                                        width.size = 7 + (5 * 3)
                                        height.size = width.size - 1
                                    }
                                    
                                    for (redDim in 1:length(redDim.method)) {
                                        
                                        # /// 2d ///
                                        dim.folder = paste(eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))),'D', sep = '')
                                        dim.finalX = eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) - (2-1)
                                        
                                        if (dim.finalX >= 1) {
                                            
                                            for (dim.comb in 1:dim.finalX) {
                                                dim.x <- dim.comb
                                                dim.y <- dim.x + 1
                                                
                                                if (dim.finalX > 1) {
                                                    if (dim.x == 1) {
                                                        comb.run = 2
                                                    } else {
                                                        comb.run = 1
                                                    }
                                                } else {
                                                    comb.run = 1
                                                }
                                                
                                                for (dim.comb.run in 1:comb.run) {
                                                    if (dim.comb.run > 1) {
                                                        dim.y <- dim.y + 1
                                                    }
                                                    
                                                    if (eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) == 2) {
                                                        dim.xy <- paste('D',dim.x,'-','D',dim.y, sep = '')
                                                    } else {
                                                        dim.xy <- paste('D',dim.x,'-','D',dim.y, sep = '')
                                                    }
                                                    #print(dim.xy)
                                                    
                                                    png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',project,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,'.',gsub("\\.","_",names(gene_set.list[N])),'.png', sep = ''), width = width.size, height = height.size, units = 'in', res = 600)
                                                    p <- plot_cells(cds, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, genes = marker.genes, norm_method = "log", label_cell_groups = FALSE, show_trajectory_graph = FALSE, label_leaves = FALSE)
                                                    print(p)
                                                    dev.off()
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    
                    
                    ## Step 3: Remove batch effects with cell alignment ***********************************************
                    if (! is.null(batch)) {
                        if (project.name == "tutorial") {
                            width.size = 6.52
                            height.size = 4
                        } else {
                            legend.max_nchar = max(nchar(levels(colData(cds)[, batch])))
                            legend.multiply_factor = 0.5 + legend.max_nchar * 0.05
                            legend.num = length(levels(colData(cds)[, batch]))
                            legend.ncol = ceiling(legend.num / 10)
                            width.size = 4.075 + (legend.multiply_factor * legend.ncol)
                            height.size = 4
                        }
                        
                        for (redDim in 1:length(redDim.method)) {
                            dim.folder = paste(eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))),'D', sep = '')
                            dim.finalX = eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) - (2-1)
                            
                            if (dim.finalX >= 1) {
                                
                                for (dim.comb in 1:dim.finalX) {
                                    dim.x <- dim.comb
                                    dim.y <- dim.x + 1
                                    
                                    if (dim.finalX > 1) {
                                        if (dim.x == 1) {
                                            comb.run = 2
                                        } else {
                                            comb.run = 1
                                        }
                                    } else {
                                        comb.run = 1
                                    }
                                    
                                    for (dim.comb.run in 1:comb.run) {
                                        if (dim.comb.run > 1) {
                                            dim.y <- dim.y + 1
                                        }
                                        
                                        if (eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) == 2) {
                                            dim.xy <- paste('D',dim.x,'-','D',dim.y, sep = '')
                                        } else {
                                            dim.xy <- paste('D',dim.x,'-','D',dim.y, sep = '')
                                        }
                                        #print(dim.xy)
                                        
                                        png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',project,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,'.',gsub("\\.","_",batch),'.png', sep = ''), width = width.size, height = height.size, units = 'in', res = 600)
                                        p <- plot_cells(cds, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = batch, label_cell_groups = FALSE)
                                        p <- p + theme(legend.position = "right")
                                        p$guides$colour$ncol = legend.ncol
                                        print(p)
                                        dev.off()
                                    }
                                }
                            }
                        }
                        
                        # residual_model_formula_str: This argument is for subtracting continuous effects. You can use this to control for things like the fraction of mitochondrial reads in each cell, which is sometimes used as a QC metric for each cell.
                        # Passing these colums as terms in the residual_model_formula_str tells align_cds() to subtract these signals prior to dimensionality reduction, clustering, and trajectory inference.
                        # Note that you can call align_cds() with alignment_group, residual_model_formula, or both.
                        if (project.name == "tutorial") {
                            cds <- suppressMessages(
                                align_cds(cds, alignment_group = batch, residual_model_formula_str = "~ bg.300.loading + bg.400.loading + bg.500.1.loading + bg.500.2.loading + bg.r17.loading + bg.b01.loading + bg.b02.loading") # Each of the columns bg.300.loading, bg.400.loading, corresponds to a background signal that a cell might be contaminated with.
                            )
                        } else {
                            cds <- suppressMessages(
                                align_cds(cds, alignment_group = batch, num_dim = 50) # to using the alignment_group argument to align_cds(), which aligns groups of cells (i.e. batches)
                            )
                        }
                        
                        for (redDim in 1:length(redDim.method)) {
                            dim.folder = paste(eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))),'D', sep = '')
                            dim.finalX = eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) - (2-1)
                            
                            if (dim.finalX >= 1) {
                                
                                for (dim.comb in 1:dim.finalX) {
                                    dim.x <- dim.comb
                                    dim.y <- dim.x + 1
                                    
                                    if (dim.finalX > 1) {
                                        if (dim.x == 1) {
                                            comb.run = 2
                                        } else {
                                            comb.run = 1
                                        }
                                    } else {
                                        comb.run = 1
                                    }
                                    
                                    for (dim.comb.run in 1:comb.run) {
                                        if (dim.comb.run > 1) {
                                            dim.y <- dim.y + 1
                                        }
                                        
                                        if (eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) == 2) {
                                            dim.xy <- paste('D',dim.x,'-','D',dim.y, sep = '')
                                        } else {
                                            dim.xy <- paste('D',dim.x,'-','D',dim.y, sep = '')
                                        }
                                        #print(dim.xy)
                                        
                                        png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',project,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'png', sep = ''), width = width.size, height = height.size, units = 'in', res = 600)
                                        p <- plot_cells(cds, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = batch, label_cell_groups = FALSE)
                                        p <- p + theme(legend.position = "right")
                                        p$guides$colour$ncol = legend.ncol
                                        print(p)
                                        dev.off()
                                    }
                                }
                            }
                        }
                    }
                    
                    
                    ## Step 4: Cluster the cells **********************************************************************
                    cds <- cluster_cells(cds, reduction_method = "PCA")
                    colData(cds)$Cluster = clusters(cds, reduction_method = "PCA")
                    
                    for (redDim in 1:length(redDim.method)) {
                        # Cluster your cells
                        # Recall that we run cluster_cells(), each cell is assigned not only to a cluster but also to a partition. 
                        # When you are learning trajectories, each partition will eventually become a separate trajectory.
                        cds <- cluster_cells(cds, reduction_method = redDim.method[redDim]) # resolution: Parameter that controls the resolution of clustering. If NULL (Default), the parameter is determined automatically.
                        
                        if (redDim.method[redDim] == "UMAP") {
                            colData(cds)$louvain_component = cds@clusters[[redDim.method[redDim]]]$partitions
                        }
                        
                        # /// 2d ///
                        dim.folder = paste(eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))),'D', sep = '')
                        dim.finalX = eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) - (2-1)
                        
                        if (dim.finalX >= 1) {
                            
                            for (dim.comb in 1:dim.finalX) {
                                dim.x <- dim.comb
                                dim.y <- dim.x + 1
                                
                                if (dim.finalX > 1) {
                                    if (dim.x == 1) {
                                        comb.run = 2
                                    } else {
                                        comb.run = 1
                                    }
                                } else {
                                    comb.run = 1
                                }
                                
                                for (dim.comb.run in 1:comb.run) {
                                    if (dim.comb.run > 1) {
                                        dim.y <- dim.y + 1
                                    }
                                    
                                    if (eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) == 2) {
                                        dim.xy <- paste('D',dim.x,'-','D',dim.y, sep = '')
                                    } else {
                                        dim.xy <- paste('D',dim.x,'-','D',dim.y, sep = '')
                                    }
                                    #print(dim.xy)
                                    
                                    mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder, sep = '')
                                    subDir <- dim.xy
                                    if (! file.exists(file.path(mainDir, subDir))) {
                                        dir.create(file.path(mainDir, subDir))
                                    }
                                    
                                    legend.max_nchar = max(nchar(levels(cds@clusters[[redDim.method[redDim]]]$clusters)))
                                    legend.multiply_factor = 0.5 + legend.max_nchar * 0.05
                                    legend.num = length(levels(cds@clusters[[redDim.method[redDim]]]$clusters))
                                    legend.ncol = ceiling(legend.num / 10)
                                    width.size = 4.075 + (legend.multiply_factor * legend.ncol)
                                    
                                    png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',project,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'cluster','.png', sep = ''), width = 4.075, height = 4, units = 'in', res = 600)
                                    p <- plot_cells(cds, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d) # call plot_cells() with no arguments, it colors the cells by cluster according to default.
                                    print(p)
                                    dev.off()
                                    
                                    png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',project,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'cluster','_','legend','.png', sep = ''), width = width.size, height = 4, units = 'in', res = 600)
                                    p <- plot_cells(cds, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d) # call plot_cells() with no arguments, it colors the cells by cluster according to default.
                                    p <- p + theme(legend.position = "right")
                                    p$guides$colour$ncol = legend.ncol
                                    print(p)
                                    dev.off()
                                    
                                    if (length(cell_type) > 0) {
                                        for (color.plt in 1:length(cell_type)) {
                                            legend.max_nchar = max(nchar(levels(colData(cds)[, cell_type[color.plt]])))
                                            legend.multiply_factor = 0.5 + legend.max_nchar * 0.05
                                            legend.num = length(levels(colData(cds)[, cell_type[color.plt]]))
                                            legend.ncol = ceiling(legend.num / 10)
                                            width.size = 4.075 + (legend.multiply_factor * legend.ncol)
                                            
                                            png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',project,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'cluster','.',gsub("\\.","_",cell_type[color.plt]),'.png', sep = ''), width = 4.075, height = 4, units = 'in', res = 600)
                                            p <- plot_cells(cds, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = cell_type[color.plt]) # each cluster is labeled according the most common annotation within it
                                            if (! is.null(color_palette)) {
                                                cell_type.colidx <- which(paste(unique(c(cell_type[color.plt], unlist(strsplit(cell_type[color.plt], '\\.|_'))[1])),'.','color', sep = '') %in% names(color_palette))
                                                if (length(cell_type.colidx) == 1) {
                                                    p <- p + scale_color_manual(values = color_palette[[which(names(color_palette) %in% paste(unique(c(cell_type[color.plt], unlist(strsplit(cell_type[color.plt], '\\.|_'))[1])),'.','color', sep = ''))]])
                                                }
                                            }
                                            print(p)
                                            dev.off()
                                            
                                            png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',project,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'cluster','.',gsub("\\.","_",cell_type[color.plt]),'_','legend','.png', sep = ''), width = width.size, height = 4, units = 'in', res = 600)
                                            p <- plot_cells(cds, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = cell_type[color.plt], label_cell_groups = FALSE) # each cluster is labeled according the most common annotation within it
                                            if (! is.null(color_palette)) {
                                                cell_type.colidx <- which(paste(unique(c(cell_type[color.plt], unlist(strsplit(cell_type[color.plt], '\\.|_'))[1])),'.','color', sep = '') %in% names(color_palette))
                                                if (length(cell_type.colidx) == 1) {
                                                    p <- p + scale_color_manual(values = color_palette[[which(names(color_palette) %in% paste(unique(c(cell_type[color.plt], unlist(strsplit(cell_type[color.plt], '\\.|_'))[1])),'.','color', sep = ''))]])
                                                }
                                            }
                                            p <- p + theme(legend.position = "right")
                                            p$guides$colour$ncol = legend.ncol
                                            if (project.name %in% c("GSE139944", c("GSM4150376", "GSM4150377", "GSM4150378", "GSM4150379"))) {
                                                if (cell_type[color.plt] == "pathway_level_1") {
                                                    p$guides$colour$title = "Pathway"
                                                } else if (cell_type[color.plt] == "product_name") {
                                                    p$guides$colour$title = "Treatment"
                                                } else if (cell_type[color.plt] == "dose_character") {
                                                    p$guides$colour$title = dose_character.unit
                                                } else {
                                                    p$guides$colour$title = "Treatment"
                                                }
                                            }
                                            print(p)
                                            dev.off()
                                            
                                            png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',project,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'cluster','.',gsub("\\.","_",cell_type[color.plt]),'.','disableLabel','.png', sep = ''), width = 4.075, height = 4, units = 'in', res = 600)
                                            p <- plot_cells(cds, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = cell_type[color.plt], label_groups_by_cluster = FALSE) # Finally, you can disable this labeling policy, making plot_cells() behave like it did before we called cluster_cells()
                                            if (! is.null(color_palette)) {
                                                cell_type.colidx <- which(paste(unique(c(cell_type[color.plt], unlist(strsplit(cell_type[color.plt], '\\.|_'))[1])),'.','color', sep = '') %in% names(color_palette))
                                                if (length(cell_type.colidx) == 1) {
                                                    p <- p + scale_color_manual(values = color_palette[[which(names(color_palette) %in% paste(unique(c(cell_type[color.plt], unlist(strsplit(cell_type[color.plt], '\\.|_'))[1])),'.','color', sep = ''))]])
                                                }
                                            }
                                            print(p)
                                            dev.off()
                                            
                                            png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',project,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'cluster','.',gsub("\\.","_",cell_type[color.plt]),'_','legend','.','disableLabel','.png', sep = ''), width = width.size, height = 4, units = 'in', res = 600)
                                            p <- plot_cells(cds, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = cell_type[color.plt], label_groups_by_cluster = FALSE, label_cell_groups = FALSE) # Finally, you can disable this labeling policy, making plot_cells() behave like it did before we called cluster_cells()
                                            if (! is.null(color_palette)) {
                                                cell_type.colidx <- which(paste(unique(c(cell_type[color.plt], unlist(strsplit(cell_type[color.plt], '\\.|_'))[1])),'.','color', sep = '') %in% names(color_palette))
                                                if (length(cell_type.colidx) == 1) {
                                                    p <- p + scale_color_manual(values = color_palette[[which(names(color_palette) %in% paste(unique(c(cell_type[color.plt], unlist(strsplit(cell_type[color.plt], '\\.|_'))[1])),'.','color', sep = ''))]])
                                                }
                                            }
                                            p <- p + theme(legend.position = "right")
                                            p$guides$colour$ncol = legend.ncol
                                            if (project.name %in% c("GSE139944", c("GSM4150376", "GSM4150377", "GSM4150378", "GSM4150379"))) {
                                                if (cell_type[color.plt] == "pathway_level_1") {
                                                    p$guides$colour$title = "Pathway"
                                                } else if (cell_type[color.plt] == "product_name") {
                                                    p$guides$colour$title = "Treatment"
                                                } else if (cell_type[color.plt] == "dose_character") {
                                                    p$guides$colour$title = dose_character.unit
                                                } else {
                                                    p$guides$colour$title = "Treatment"
                                                }
                                            }
                                            print(p)
                                            dev.off()
                                        }
                                    }
                                    
                                    # The cluster_cells() also divides the cells into larger, more well separated groups called partitions, 
                                    # using a statistical test from Alex Wolf et al, introduced as part of their PAGA algorithm. (https://github.com/theislab/paga)
                                    png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',project,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'partition','.png', sep = ''), width = 4.075, height = 4, units = 'in', res = 600)
                                    p <- plot_cells(cds, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = "partition")
                                    print(p)
                                    dev.off()
                                    
                                    png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',project,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'partition','.','group_cells','.png', sep = ''), width = 4.075, height = 4, units = 'in', res = 600)
                                    p <- plot_cells(cds, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = "partition", group_cells_by = "partition")
                                    print(p)
                                    dev.off()
                                }
                            }
                        }
                    }
                    
                    
                    # --- subset ----------------------------------------------------------------------------------------------------------------------------------------
                    for (redDim in 1:length(redDim.method)) {
                        dim.folder = paste(eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))),'D', sep = '')
                        dim.finalX = eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) - (2-1)
                        
                        if (dim.finalX >= 1) {
                            
                            for (dim.comb in 1:dim.finalX) {
                                dim.x <- dim.comb
                                dim.y <- dim.x + 1
                                
                                if (dim.finalX > 1) {
                                    if (dim.x == 1) {
                                        comb.run = 2
                                    } else {
                                        comb.run = 1
                                    }
                                } else {
                                    comb.run = 1
                                }
                                
                                for (dim.comb.run in 1:comb.run) {
                                    if (dim.comb.run > 1) {
                                        dim.y <- dim.y + 1
                                    }
                                    
                                    if (eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) == 2) {
                                        dim.xy <- paste('D',dim.x,'-','D',dim.y, sep = '')
                                    } else {
                                        dim.xy <- paste('D',dim.x,'-','D',dim.y, sep = '')
                                    }
                                    #print(dim.xy)
                                    
                                    
                                    # /// subset ///
                                    for (SUB in SUB.run) {
                                        if (SUB == 0) {
                                            next # next halts the processing of the current iteration and advances the looping index.
                                        } else if (SUB == 1) {
                                            SUB.prefolder = 'cell_subset'
                                            SUB.prename = 'cell_subset'
                                        } else if (SUB == 2) {
                                            if (length(cell_subset) == 0 & cell_choose == TRUE) {
                                                SUB.prefolder = 'choose'
                                                SUB.prename = 'choose'
                                            } else if (length(cell_subset) > 0 & cell_choose == TRUE) {
                                                SUB.prefolder = paste('cell_subset','/','choose', sep = '')
                                                SUB.prename = paste('cell_subset','.','choose', sep = '')
                                            }
                                        }
                                        
                                        mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy, sep = '')
                                        subDir <- SUB.prefolder
                                        if (! file.exists(file.path(mainDir, subDir))) {
                                            dir.create(file.path(mainDir, subDir), recursive = TRUE)
                                        }
                                        
                                        if (dim.comb == 1 & dim.comb.run == 1) {
                                            if (SUB == 1) {
                                                cds_choose.cell <- cell_subset
                                                
                                                cds.colidx <- which(sapply(colnames(cds), FUN = function(X) X %in% cds_choose.cell))
                                                if (length(cds.colidx) > 0) {
                                                    cds_choose <- cds[, cds.colidx]
                                                } else {
                                                    cds_choose <- cds
                                                }
                                            } else if (SUB == 2) {
                                                if (length(cell_subset) == 0 & cell_choose == TRUE) {
                                                    cds4choose <- cds
                                                } else if (length(cell_subset) > 0 & cell_choose == TRUE) {
                                                    cds4choose <- cds_choose
                                                }
                                                
                                                # Choose cells interactively to subset a cds
                                                # clear_cds: Logical, clear CDS slots before returning. After clearing the cds, re-run processing from preprocess_cds(), ... Default is FALSE.
                                                # return_list: Logical, return a list of cells instead of a subsetted CDS object.
                                                #cds_choose.before <- choose_cells(cds4choose, reduction_method = redDim.method[redDim], clear_cds = TRUE, return_list = FALSE)
                                                #cds_choose.after <- choose_cells(cds4choose, reduction_method = redDim.method[redDim], clear_cds = FALSE, return_list = FALSE)
                                                cds_choose.cell <- choose_cells(cds4choose, reduction_method = redDim.method[redDim], clear_cds = FALSE, return_list = TRUE)
                                                
                                                if (length(cell_subset) > 0) {
                                                    cds_choose.cell <- intersect(cell_subset, cds_choose.cell)
                                                }
                                                
                                                cds4choose.colidx <- which(sapply(colnames(cds4choose), FUN = function(X) X %in% cds_choose.cell))
                                                if (length(cds4choose.colidx) > 0) {
                                                    cds_choose <- cds4choose[, cds4choose.colidx]
                                                } else {
                                                    cds_choose <- cds4choose
                                                }
                                            }
                                            
                                            mapping.idx = sapply(rownames(colData(cds_choose)), FUN = function(X) which(names(cds_choose@clusters[[redDim.method[redDim]]]$clusters) %in% X))
                                            cds_choose@clusters[[redDim.method[redDim]]]$clusters = cds_choose@clusters[[redDim.method[redDim]]]$clusters[mapping.idx]
                                            cds_choose@clusters[[redDim.method[redDim]]]$clusters = factor(cds_choose@clusters[[redDim.method[redDim]]]$clusters, levels = mixedsort(as.character(unique(cds_choose@clusters[[redDim.method[redDim]]]$clusters))))
                                            colData(cds_choose)$cluster = cds_choose@clusters[[redDim.method[redDim]]]$clusters
                                            
                                            mapping.idx = sapply(rownames(colData(cds_choose)), FUN = function(X) which(names(cds_choose@clusters[[redDim.method[redDim]]]$partitions) %in% X))
                                            cds_choose@clusters[[redDim.method[redDim]]]$partitions = cds_choose@clusters[[redDim.method[redDim]]]$partitions[mapping.idx]
                                            cds_choose@clusters[[redDim.method[redDim]]]$partitions = factor(cds_choose@clusters[[redDim.method[redDim]]]$partitions, levels = mixedsort(as.character(unique(cds_choose@clusters[[redDim.method[redDim]]]$partitions))))
                                            colData(cds_choose)$partition = cds_choose@clusters[[redDim.method[redDim]]]$partitions
                                            
                                            mapping.idx = sapply(rownames(colData(cds_choose)), FUN = function(X) which(rownames(cds_choose@principal_graph_aux[[redDim.method[redDim]]]$pr_graph_cell_proj_closest_vertex) %in% X))
                                            cds_choose@principal_graph_aux[[redDim.method[redDim]]]$pr_graph_cell_proj_closest_vertex = cds_choose@principal_graph_aux[[redDim.method[redDim]]]$pr_graph_cell_proj_closest_vertex[mapping.idx, , drop = FALSE]
                                            
                                            saveRDS(cds_choose, paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',project,'.',SUB.prename,'.','trajectories','.','cds_choose','.rds', sep = ''))
                                            #cds_choose <- readRDS(paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',project,'.',SUB.prename,'.','trajectories','.','cds_choose','.rds', sep = ''))
                                            
                                            if (SUB == 1) {
                                                cds_choose.cell_subset <- cds_choose
                                            } else if (SUB == 2) {
                                                cds_choose.cell_choose <- cds_choose
                                            }
                                        } else {
                                            if (SUB == 1) {
                                                cds_choose <- cds_choose.cell_subset
                                            } else if (SUB == 2) {
                                                cds_choose <- cds_choose.cell_choose
                                            }
                                        }
                                        
                                        
                                        ## Step 2: Reduce the dimensions using UMAP *******************************************************
                                        # /// no color_cells ///
                                        png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',SUB.prefolder,'/',project,'.',SUB.prename,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,'.png', sep = ''), width = 4.075, height = 4, units = 'in', res = 600)
                                        p <- plot_cells(cds_choose, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d)
                                        print(p)
                                        dev.off()
                                        
                                        # /// color_cells_by = cell_type ///
                                        if (length(cell_type) > 0) {
                                            for (color.plt in 1:length(cell_type)) {
                                                legend.max_nchar = max(nchar(levels(colData(cds_choose)[, cell_type[color.plt]])))
                                                legend.multiply_factor = 0.5 + legend.max_nchar * 0.05
                                                legend.num = length(levels(colData(cds_choose)[, cell_type[color.plt]]))
                                                legend.ncol = ceiling(legend.num / 10)
                                                width.size = 4.075 + (legend.multiply_factor * legend.ncol)
                                                
                                                png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',SUB.prefolder,'/',project,'.',SUB.prename,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,'.',gsub("\\.","_",cell_type[color.plt]),'.png', sep = ''), width = 4.075, height = 4, units = 'in', res = 600)
                                                p <- plot_cells(cds_choose, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = cell_type[color.plt], label_groups_by_cluster = FALSE) # Finally, you can disable this labeling policy, making plot_cells() behave like it did before we called cluster_cells()
                                                if (! is.null(color_palette)) {
                                                    cell_type.colidx <- which(paste(unique(c(cell_type[color.plt], unlist(strsplit(cell_type[color.plt], '\\.|_'))[1])),'.','color', sep = '') %in% names(color_palette))
                                                    if (length(cell_type.colidx) == 1) {
                                                        p <- p + scale_color_manual(values = color_palette[[which(names(color_palette) %in% paste(unique(c(cell_type[color.plt], unlist(strsplit(cell_type[color.plt], '\\.|_'))[1])),'.','color', sep = ''))]])
                                                    }
                                                }
                                                print(p)
                                                dev.off()
                                                
                                                png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',SUB.prefolder,'/',project,'.',SUB.prename,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,'.',gsub("\\.","_",cell_type[color.plt]),'_','legend','.png', sep = ''), width = width.size, height = 4, units = 'in', res = 600)
                                                p <- plot_cells(cds_choose, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = cell_type[color.plt], label_groups_by_cluster = FALSE, label_cell_groups = FALSE) # Finally, you can disable this labeling policy, making plot_cells() behave like it did before we called cluster_cells()
                                                if (! is.null(color_palette)) {
                                                    cell_type.colidx <- which(paste(unique(c(cell_type[color.plt], unlist(strsplit(cell_type[color.plt], '\\.|_'))[1])),'.','color', sep = '') %in% names(color_palette))
                                                    if (length(cell_type.colidx) == 1) {
                                                        p <- p + scale_color_manual(values = color_palette[[which(names(color_palette) %in% paste(unique(c(cell_type[color.plt], unlist(strsplit(cell_type[color.plt], '\\.|_'))[1])),'.','color', sep = ''))]])
                                                    }
                                                }
                                                p <- p + theme(legend.position = "right")
                                                p$guides$colour$ncol = legend.ncol
                                                if (project.name %in% c("GSE139944", c("GSM4150376", "GSM4150377", "GSM4150378", "GSM4150379"))) {
                                                    if (cell_type[color.plt] == "pathway_level_1") {
                                                        p$guides$colour$title = "Pathway"
                                                    } else if (cell_type[color.plt] == "product_name") {
                                                        p$guides$colour$title = "Treatment"
                                                    } else if (cell_type[color.plt] == "dose_character") {
                                                        p$guides$colour$title = dose_character.unit
                                                    } else {
                                                        p$guides$colour$title = "Treatment"
                                                    }
                                                }
                                                print(p)
                                                dev.off()
                                            }
                                        }
                                        
                                        
                                        ## Step 3: Remove batch effects with cell alignment ***********************************************
                                        if (! is.null(batch)) {
                                            if (project.name == "tutorial") {
                                                width.size = 6.52
                                                height.size = 4
                                            } else {
                                                legend.max_nchar = max(nchar(levels(colData(cds_choose)[, batch])))
                                                legend.multiply_factor = 0.5 + legend.max_nchar * 0.05
                                                legend.num = length(levels(colData(cds_choose)[, batch]))
                                                legend.ncol = ceiling(legend.num / 10)
                                                width.size = 4.075 + (legend.multiply_factor * legend.ncol)
                                                height.size = 4
                                            }
                                            
                                            png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',SUB.prefolder,'/',project,'.',SUB.prename,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'png', sep = ''), width = width.size, height = height.size, units = 'in', res = 600)
                                            p <- plot_cells(cds_choose, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = batch, label_cell_groups = FALSE)
                                            p <- p + theme(legend.position = "right")
                                            p$guides$colour$ncol = legend.ncol
                                            print(p)
                                            dev.off()
                                        }
                                        
                                        
                                        ## Step 4: Cluster the cells **********************************************************************
                                        legend.max_nchar = max(nchar(levels(cds_choose@clusters[[redDim.method[redDim]]]$clusters)))
                                        legend.multiply_factor = 0.5 + legend.max_nchar * 0.05
                                        legend.num = length(levels(cds_choose@clusters[[redDim.method[redDim]]]$clusters))
                                        legend.ncol = ceiling(legend.num / 10)
                                        width.size = 4.075 + (legend.multiply_factor * legend.ncol)
                                        
                                        png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',SUB.prefolder,'/',project,'.',SUB.prename,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'cluster','.png', sep = ''), width = 4.075, height = 4, units = 'in', res = 600)
                                        p <- plot_cells(cds_choose, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d) # call plot_cells() with no arguments, it colors the cells by cluster according to default.
                                        print(p)
                                        dev.off()
                                        
                                        png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',SUB.prefolder,'/',project,'.',SUB.prename,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'cluster','_','legend','.png', sep = ''), width = width.size, height = 4, units = 'in', res = 600)
                                        p <- plot_cells(cds_choose, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d) # call plot_cells() with no arguments, it colors the cells by cluster according to default.
                                        p <- p + theme(legend.position = "right")
                                        p$guides$colour$ncol = legend.ncol
                                        print(p)
                                        dev.off()
                                        
                                        if (length(cell_type) > 0) {
                                            for (color.plt in 1:length(cell_type)) {
                                                legend.max_nchar = max(nchar(levels(colData(cds_choose)[, cell_type[color.plt]])))
                                                legend.multiply_factor = 0.5 + legend.max_nchar * 0.05
                                                legend.num = length(levels(colData(cds_choose)[, cell_type[color.plt]]))
                                                legend.ncol = ceiling(legend.num / 10)
                                                width.size = 4.075 + (legend.multiply_factor * legend.ncol)
                                                
                                                png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',SUB.prefolder,'/',project,'.',SUB.prename,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'cluster','.',gsub("\\.","_",cell_type[color.plt]),'.png', sep = ''), width = 4.075, height = 4, units = 'in', res = 600)
                                                p <- plot_cells(cds_choose, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = cell_type[color.plt]) # each cluster is labeled according the most common annotation within it
                                                if (! is.null(color_palette)) {
                                                    cell_type.colidx <- which(paste(unique(c(cell_type[color.plt], unlist(strsplit(cell_type[color.plt], '\\.|_'))[1])),'.','color', sep = '') %in% names(color_palette))
                                                    if (length(cell_type.colidx) == 1) {
                                                        p <- p + scale_color_manual(values = color_palette[[which(names(color_palette) %in% paste(unique(c(cell_type[color.plt], unlist(strsplit(cell_type[color.plt], '\\.|_'))[1])),'.','color', sep = ''))]])
                                                    }
                                                }
                                                print(p)
                                                dev.off()
                                                
                                                png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',SUB.prefolder,'/',project,'.',SUB.prename,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'cluster','.',gsub("\\.","_",cell_type[color.plt]),'_','legend','.png', sep = ''), width = width.size, height = 4, units = 'in', res = 600)
                                                p <- plot_cells(cds_choose, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = cell_type[color.plt], label_cell_groups = FALSE) # each cluster is labeled according the most common annotation within it
                                                if (! is.null(color_palette)) {
                                                    cell_type.colidx <- which(paste(unique(c(cell_type[color.plt], unlist(strsplit(cell_type[color.plt], '\\.|_'))[1])),'.','color', sep = '') %in% names(color_palette))
                                                    if (length(cell_type.colidx) == 1) {
                                                        p <- p + scale_color_manual(values = color_palette[[which(names(color_palette) %in% paste(unique(c(cell_type[color.plt], unlist(strsplit(cell_type[color.plt], '\\.|_'))[1])),'.','color', sep = ''))]])
                                                    }
                                                }
                                                p <- p + theme(legend.position = "right")
                                                p$guides$colour$ncol = legend.ncol
                                                if (project.name %in% c("GSE139944", c("GSM4150376", "GSM4150377", "GSM4150378", "GSM4150379"))) {
                                                    if (cell_type[color.plt] == "pathway_level_1") {
                                                        p$guides$colour$title = "Pathway"
                                                    } else if (cell_type[color.plt] == "product_name") {
                                                        p$guides$colour$title = "Treatment"
                                                    } else if (cell_type[color.plt] == "dose_character") {
                                                        p$guides$colour$title = dose_character.unit
                                                    } else {
                                                        p$guides$colour$title = "Treatment"
                                                    }
                                                }
                                                print(p)
                                                dev.off()
                                                
                                                png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',SUB.prefolder,'/',project,'.',SUB.prename,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'cluster','.',gsub("\\.","_",cell_type[color.plt]),'.','disableLabel','.png', sep = ''), width = 4.075, height = 4, units = 'in', res = 600)
                                                p <- plot_cells(cds_choose, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = cell_type[color.plt], label_groups_by_cluster = FALSE) # Finally, you can disable this labeling policy, making plot_cells() behave like it did before we called cluster_cells()
                                                if (! is.null(color_palette)) {
                                                    cell_type.colidx <- which(paste(unique(c(cell_type[color.plt], unlist(strsplit(cell_type[color.plt], '\\.|_'))[1])),'.','color', sep = '') %in% names(color_palette))
                                                    if (length(cell_type.colidx) == 1) {
                                                        p <- p + scale_color_manual(values = color_palette[[which(names(color_palette) %in% paste(unique(c(cell_type[color.plt], unlist(strsplit(cell_type[color.plt], '\\.|_'))[1])),'.','color', sep = ''))]])
                                                    }
                                                }
                                                print(p)
                                                dev.off()
                                                
                                                png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',SUB.prefolder,'/',project,'.',SUB.prename,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'cluster','.',gsub("\\.","_",cell_type[color.plt]),'_','legend','.','disableLabel','.png', sep = ''), width = width.size, height = 4, units = 'in', res = 600)
                                                p <- plot_cells(cds_choose, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = cell_type[color.plt], label_groups_by_cluster = FALSE, label_cell_groups = FALSE) # Finally, you can disable this labeling policy, making plot_cells() behave like it did before we called cluster_cells()
                                                if (! is.null(color_palette)) {
                                                    cell_type.colidx <- which(paste(unique(c(cell_type[color.plt], unlist(strsplit(cell_type[color.plt], '\\.|_'))[1])),'.','color', sep = '') %in% names(color_palette))
                                                    if (length(cell_type.colidx) == 1) {
                                                        p <- p + scale_color_manual(values = color_palette[[which(names(color_palette) %in% paste(unique(c(cell_type[color.plt], unlist(strsplit(cell_type[color.plt], '\\.|_'))[1])),'.','color', sep = ''))]])
                                                    }
                                                }
                                                p <- p + theme(legend.position = "right")
                                                p$guides$colour$ncol = legend.ncol
                                                if (project.name %in% c("GSE139944", c("GSM4150376", "GSM4150377", "GSM4150378", "GSM4150379"))) {
                                                    if (cell_type[color.plt] == "pathway_level_1") {
                                                        p$guides$colour$title = "Pathway"
                                                    } else if (cell_type[color.plt] == "product_name") {
                                                        p$guides$colour$title = "Treatment"
                                                    } else if (cell_type[color.plt] == "dose_character") {
                                                        p$guides$colour$title = dose_character.unit
                                                    } else {
                                                        p$guides$colour$title = "Treatment"
                                                    }
                                                }
                                                print(p)
                                                dev.off()
                                            }
                                        }
                                        
                                        # The cluster_cells() also divides the cells into larger, more well separated groups called partitions, 
                                        # using a statistical test from Alex Wolf et al, introduced as part of their PAGA algorithm. (https://github.com/theislab/paga)
                                        png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',SUB.prefolder,'/',project,'.',SUB.prename,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'partition','.png', sep = ''), width = 4.075, height = 4, units = 'in', res = 600)
                                        p <- plot_cells(cds_choose, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = "partition")
                                        print(p)
                                        dev.off()
                                        
                                        png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',SUB.prefolder,'/',project,'.',SUB.prename,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'partition','.','group_cells','.png', sep = ''), width = 4.075, height = 4, units = 'in', res = 600)
                                        p <- plot_cells(cds_choose, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = "partition", group_cells_by = "partition")
                                        print(p)
                                        dev.off()
                                    }
                                }
                            }
                        }
                    }
                    
                    
                    if (graph.run == 1) {
                        
                        ## Step 5: Learn a graph **************************************************************************
                        # (Optional) Order cells in pseudotime along a trajectory: https://cole-trapnell-lab.github.io/monocle3/docs/trajectories/
                        
                        # Learn the trajectory graph
                        # Next, we will fit a principal graph within each parition using the learn_graph() function
                        # This graph will be used in many downstream steps, such as branch analysis and differential expression.
                        
                        # Use custom dimension reduction for learn_graph: https://github.com/cole-trapnell-lab/monocle3/issues/241
                        
                        for (redDim in 1:length(redDim.method)) {
                            
                            if (redDim.method[redDim] == "UMAP") {
                                
                            } else if (redDim.method[redDim] == "tSNE") {
                                next # next halts the processing of the current iteration and advances the looping index.
                            }
                            
                            cds.before_graph <- cds
                            
                            if (project.name %in% c("GSM4150377", "GSM4150378")) {
                                graph_parameters <- list()
                                graph_parameters[["minimal_branch_len"]] = 10
                                graph_parameters[["ncenter"]] = 750
                                
                                cds <- learn_graph(cds.before_graph, use_partition = TRUE, close_loop = FALSE, learn_graph_control = graph_parameters) #Please run reduce_dimension with reduction_method = UMAP and cluster_cells before running learn_graph.
                                cds.no_partition <- learn_graph(cds.before_graph, use_partition = FALSE, close_loop = FALSE, learn_graph_control = graph_parameters) #Please run reduce_dimension with reduction_method = UMAP and cluster_cells before running learn_graph.
                                
                                cds <- append_umap_coordinates(cds)
                            } else {
                                cds <- learn_graph(cds.before_graph, use_partition = TRUE) #Please run reduce_dimension with reduction_method = UMAP and cluster_cells before running learn_graph.
                                cds.no_partition <- learn_graph(cds.before_graph, use_partition = FALSE) #Please run reduce_dimension with reduction_method = UMAP and cluster_cells before running learn_graph.
                            }
                            
                            # refer to select_cells.R
                            sel <- c()
                            sel$nodes <- igraph::V(principal_graph(cds)[[redDim.method[redDim]]])$name
                            sel$cells <- colnames(cds)
                            principal_graph(cds)[[redDim.method[redDim]]] <- igraph::induced_subgraph(principal_graph(cds)[[redDim.method[redDim]]], sel$nodes)
                            
                            sel <- c()
                            sel$nodes <- igraph::V(principal_graph(cds.no_partition)[[redDim.method[redDim]]])$name
                            sel$cells <- colnames(cds.no_partition)
                            principal_graph(cds.no_partition)[[redDim.method[redDim]]] <- igraph::induced_subgraph(principal_graph(cds.no_partition)[[redDim.method[redDim]]], sel$nodes)
                            
                            cds.origin.processed = cds
                            cds.no_partition.origin.processed = cds.no_partition
                            
                            
                            dim.folder = paste(eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))),'D', sep = '')
                            dim.finalX = eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) - (2-1)
                            
                            if (dim.finalX >= 1) {
                                
                                for (dim.comb in 1:dim.finalX) {
                                    dim.x <- dim.comb
                                    dim.y <- dim.x + 1
                                    
                                    if (dim.finalX > 1) {
                                        if (dim.x == 1) {
                                            comb.run = 2
                                        } else {
                                            comb.run = 1
                                        }
                                    } else {
                                        comb.run = 1
                                    }
                                    
                                    for (dim.comb.run in 1:comb.run) {
                                        if (dim.comb.run > 1) {
                                            dim.y <- dim.y + 1
                                        }
                                        
                                        if (eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) == 2) {
                                            dim.xy <- paste('D',dim.x,'-','D',dim.y, sep = '')
                                        } else {
                                            dim.xy <- paste('D',dim.x,'-','D',dim.y, sep = '')
                                        }
                                        #print(dim.xy)
                                        
                                        mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder, sep = '')
                                        subDir <- dim.xy
                                        if (! file.exists(file.path(mainDir, subDir))) {
                                            dir.create(file.path(mainDir, subDir))
                                        }
                                        
                                        #for (trj.prt in 1:1) {
                                        for (trj.prt in 1:0) {
                                            if (trj.prt == 1) {
                                                trj.prt.name = 'with_partition'
                                                cds.plt = cds
                                            } else if (trj.prt == 0) {
                                                trj.prt.name = 'without_partition'
                                                cds.plt = cds.no_partition
                                            }
                                            
                                            if (! is.null(cds.plt)) {
                                                
                                                if (length(cell_type) > 0) {
                                                    for (color.plt in 1:length(cell_type)) {
                                                        legend.max_nchar = max(nchar(levels(colData(cds.plt)[, cell_type[color.plt]])))
                                                        legend.multiply_factor = 0.5 + legend.max_nchar * 0.05
                                                        legend.num = length(levels(colData(cds.plt)[, cell_type[color.plt]]))
                                                        legend.ncol = ceiling(legend.num / 10)
                                                        width.size = 4.075 + (legend.multiply_factor * legend.ncol)
                                                        
                                                        png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',project,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'graph','_',trj.prt.name,'.',gsub("\\.","_",cell_type[color.plt]),'.','disableLabel','.png', sep = ''), width = 4.075, height = 4, units = 'in', res = 600)
                                                        p <- plot_cells(cds.plt, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = cell_type[color.plt], label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)
                                                        if (! is.null(color_palette)) {
                                                            cell_type.colidx <- which(paste(unique(c(cell_type[color.plt], unlist(strsplit(cell_type[color.plt], '\\.|_'))[1])),'.','color', sep = '') %in% names(color_palette))
                                                            if (length(cell_type.colidx) == 1) {
                                                                p <- p + scale_color_manual(values = color_palette[[which(names(color_palette) %in% paste(unique(c(cell_type[color.plt], unlist(strsplit(cell_type[color.plt], '\\.|_'))[1])),'.','color', sep = ''))]])
                                                            }
                                                        }
                                                        print(p)
                                                        dev.off()
                                                        
                                                        png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',project,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'graph','_',trj.prt.name,'.',gsub("\\.","_",cell_type[color.plt]),'_','legend','.','disableLabel','.png', sep = ''), width = width.size, height = 4, units = 'in', res = 600)
                                                        p <- plot_cells(cds.plt, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = cell_type[color.plt], label_groups_by_cluster = FALSE, label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE)
                                                        if (! is.null(color_palette)) {
                                                            cell_type.colidx <- which(paste(unique(c(cell_type[color.plt], unlist(strsplit(cell_type[color.plt], '\\.|_'))[1])),'.','color', sep = '') %in% names(color_palette))
                                                            if (length(cell_type.colidx) == 1) {
                                                                p <- p + scale_color_manual(values = color_palette[[which(names(color_palette) %in% paste(unique(c(cell_type[color.plt], unlist(strsplit(cell_type[color.plt], '\\.|_'))[1])),'.','color', sep = ''))]])
                                                            }
                                                        }
                                                        p <- p + theme(legend.position = "right")
                                                        p$guides$colour$ncol = legend.ncol
                                                        if (project.name %in% c("GSE139944", c("GSM4150376", "GSM4150377", "GSM4150378", "GSM4150379"))) {
                                                            if (cell_type[color.plt] == "pathway_level_1") {
                                                                p$guides$colour$title = "Pathway"
                                                            } else if (cell_type[color.plt] == "product_name") {
                                                                p$guides$colour$title = "Treatment"
                                                            } else if (cell_type[color.plt] == "dose_character") {
                                                                p$guides$colour$title = dose_character.unit
                                                            } else {
                                                                p$guides$colour$title = "Treatment"
                                                            }
                                                        }
                                                        print(p)
                                                        dev.off()
                                                    }
                                                }
                                            }
                                        }
                                        
                                        
                                        # /// subset ///
                                        for (SUB in SUB.run) {
                                            if (SUB == 0) {
                                                next # next halts the processing of the current iteration and advances the looping index.
                                            } else if (SUB == 1) {
                                                SUB.prefolder = 'cell_subset'
                                                SUB.prename = 'cell_subset'
                                            } else if (SUB == 2) {
                                                if (length(cell_subset) == 0 & cell_choose == TRUE) {
                                                    SUB.prefolder = 'choose'
                                                    SUB.prename = 'choose'
                                                } else if (length(cell_subset) > 0 & cell_choose == TRUE) {
                                                    SUB.prefolder = paste('cell_subset','/','choose', sep = '')
                                                    SUB.prename = paste('cell_subset','.','choose', sep = '')
                                                }
                                            }
                                            
                                            ## Step 5: Learn a graph **************************************************************************
                                            if (dim.comb == 1 & dim.comb.run == 1) {
                                                if (SUB == 1) {
                                                    cds_choose.before_graph <- cds_choose.cell_subset
                                                } else if (SUB == 2) {
                                                    cds_choose.before_graph <- cds_choose.cell_choose
                                                }
                                                
                                                if (project.name %in% c("GSM4150377", "GSM4150378")) {
                                                    graph_parameters <- list()
                                                    graph_parameters[["minimal_branch_len"]] = 10
                                                    graph_parameters[["ncenter"]] = 750
                                                    
                                                    cds_choose <- learn_graph(cds_choose.before_graph, use_partition = TRUE, close_loop = FALSE, learn_graph_control = graph_parameters) #Please run reduce_dimension with reduction_method = UMAP and cluster_cells before running learn_graph.
                                                    cds_choose.no_partition <- learn_graph(cds_choose.before_graph, use_partition = FALSE, close_loop = FALSE, learn_graph_control = graph_parameters) #Please run reduce_dimension with reduction_method = UMAP and cluster_cells before running learn_graph.
                                                    
                                                    cds_choose <- append_umap_coordinates(cds_choose)
                                                } else {
                                                    cds_choose <- learn_graph(cds_choose.before_graph, use_partition = TRUE) #Please run reduce_dimension with reduction_method = UMAP and cluster_cells before running learn_graph.
                                                    cds_choose.no_partition <- learn_graph(cds_choose.before_graph, use_partition = FALSE) #Please run reduce_dimension with reduction_method = UMAP and cluster_cells before running learn_graph.
                                                }
                                                
                                                # refer to select_cells.R
                                                sel <- c()
                                                sel$nodes <- igraph::V(principal_graph(cds_choose)[[redDim.method[redDim]]])$name
                                                sel$cells <- colnames(cds_choose)
                                                principal_graph(cds_choose)[[redDim.method[redDim]]] <- igraph::induced_subgraph(principal_graph(cds_choose)[[redDim.method[redDim]]], sel$nodes)
                                                
                                                sel <- c()
                                                sel$nodes <- igraph::V(principal_graph(cds_choose.no_partition)[[redDim.method[redDim]]])$name
                                                sel$cells <- colnames(cds_choose.no_partition)
                                                principal_graph(cds_choose.no_partition)[[redDim.method[redDim]]] <- igraph::induced_subgraph(principal_graph(cds_choose.no_partition)[[redDim.method[redDim]]], sel$nodes)
                                                
                                                if (SUB == 1) {
                                                    cds_choose.cell_subset <- cds_choose
                                                } else if (SUB == 2) {
                                                    cds_choose.cell_choose <- cds_choose
                                                }
                                                
                                                if (SUB == 1) {
                                                    cds_choose.no_partition.cell_subset <- cds_choose.no_partition
                                                } else if (SUB == 2) {
                                                    cds_choose.no_partition.cell_choose <- cds_choose.no_partition
                                                }
                                                
                                                
                                                # Subset cells by branch
                                                # It is often useful to subset cells based on their branch in the trajectory.
                                                # The function choose_graph_segments allows you to do so interactively.
                                                # reduction_method: The reduction method to plot while choosing cells. Currently only "UMAP" is supported.
                                                
                                                # Choose cells along the path of a principal graph
                                                # clear_cds: Logical, clear CDS slots before returning. After clearing the cds, re-run processing from preprocess_cds(), ... Default is TRUE.
                                                # return_list: Logical, return a list of cells instead of a subsetted CDS object.
                                                
                                                # /// use_partition = TRUE ///
                                                if (SUB == 1) {
                                                    cds_choose.graph <- cds_choose
                                                }
                                                
                                                if (SUB == 2) {
                                                    if (length(cell_subset) == 0 & cell_choose == TRUE) {
                                                        cds4choose <- cds
                                                    } else if (length(cell_subset) > 0 & cell_choose == TRUE) {
                                                        #cds4choose <- cds_choose
                                                        cds4choose <- cds_choose.cell_subset
                                                    }
                                                    
                                                    #cds_choose.graph <- choose_graph_segments(cds4choose, reduction_method = "UMAP", clear_cds = TRUE, return_list = FALSE)
                                                    cds_choose.graph <- choose_graph_segments(cds4choose, reduction_method = "UMAP", clear_cds = FALSE, return_list = FALSE)
                                                    
                                                    while (ncol(cds_choose.graph) > length(cds_choose.cell)) {
                                                        cat('[Warning] Please re-choose cells for trajectory graph. The cells you chose are more than the chosen cells before trajectory plot.','\n', sep = '')
                                                        
                                                        #cds_choose.graph <- choose_graph_segments(cds4choose, reduction_method = "UMAP", clear_cds = TRUE, return_list = FALSE)
                                                        cds_choose.graph <- choose_graph_segments(cds4choose, reduction_method = "UMAP", clear_cds = FALSE, return_list = FALSE)
                                                    }
                                                    while (! all(colnames(cds_choose.graph) %in% cds_choose.cell)) {
                                                        cat('[Warning] Please re-choose cells for trajectory graph. The cells you chose are not the same or the subset of the chosen cells before trajectory plot.','\n', sep = '')
                                                        
                                                        #cds_choose.graph <- choose_graph_segments(cds4choose, reduction_method = "UMAP", clear_cds = TRUE, return_list = FALSE)
                                                        cds_choose.graph <- choose_graph_segments(cds4choose, reduction_method = "UMAP", clear_cds = FALSE, return_list = FALSE)
                                                    }
                                                }
                                                
                                                mapping.idx = sapply(rownames(colData(cds_choose.graph)), FUN = function(X) which(names(cds_choose.graph@clusters[[redDim.method[redDim]]]$clusters) %in% X))
                                                cds_choose.graph@clusters[[redDim.method[redDim]]]$clusters = cds_choose.graph@clusters[[redDim.method[redDim]]]$clusters[mapping.idx]
                                                cds_choose.graph@clusters[[redDim.method[redDim]]]$clusters = factor(cds_choose.graph@clusters[[redDim.method[redDim]]]$clusters, levels = mixedsort(as.character(unique(cds_choose.graph@clusters[[redDim.method[redDim]]]$clusters))))
                                                colData(cds_choose.graph)$cluster = cds_choose.graph@clusters[[redDim.method[redDim]]]$clusters
                                                
                                                mapping.idx = sapply(rownames(colData(cds_choose.graph)), FUN = function(X) which(names(cds_choose.graph@clusters[[redDim.method[redDim]]]$partitions) %in% X))
                                                cds_choose.graph@clusters[[redDim.method[redDim]]]$partitions = cds_choose.graph@clusters[[redDim.method[redDim]]]$partitions[mapping.idx]
                                                cds_choose.graph@clusters[[redDim.method[redDim]]]$partitions = factor(cds_choose.graph@clusters[[redDim.method[redDim]]]$partitions, levels = mixedsort(as.character(unique(cds_choose.graph@clusters[[redDim.method[redDim]]]$partitions))))
                                                colData(cds_choose.graph)$partition = cds_choose.graph@clusters[[redDim.method[redDim]]]$partitions
                                                
                                                mapping.idx = sapply(rownames(colData(cds_choose.graph)), FUN = function(X) which(rownames(cds_choose.graph@principal_graph_aux[[redDim.method[redDim]]]$pr_graph_cell_proj_closest_vertex) %in% X))
                                                cds_choose.graph@principal_graph_aux[[redDim.method[redDim]]]$pr_graph_cell_proj_closest_vertex = cds_choose.graph@principal_graph_aux[[redDim.method[redDim]]]$pr_graph_cell_proj_closest_vertex[mapping.idx, , drop = FALSE]
                                                
                                                saveRDS(cds_choose.graph, paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',project,'.',SUB.prename,'.','trajectories','.','cds_choose.graph','.rds', sep = ''))
                                                #cds_choose.graph <- readRDS(paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',project,'.',SUB.prename,'.','trajectories','.','cds_choose.graph','.rds', sep = ''))
                                                
                                                if (SUB == 1) {
                                                    cds_choose.cell_subset.graph <- cds_choose.graph
                                                } else if (SUB == 2) {
                                                    cds_choose.cell_choose.graph <- cds_choose.graph
                                                }
                                                
                                                # /// use_partition = FALSE ///
                                                if (SUB == 1) {
                                                    cds_choose.no_partition.graph <- cds_choose.no_partition
                                                }
                                                
                                                if (SUB == 2) {
                                                    if (length(cell_subset) == 0 & cell_choose == TRUE) {
                                                        cds4choose <- cds.no_partition
                                                    } else if (length(cell_subset) > 0 & cell_choose == TRUE) {
                                                        #cds4choose <- cds_choose.no_partition
                                                        cds4choose <- cds_choose.no_partition.cell_subset
                                                    }
                                                    
                                                    #cds_choose.no_partition.graph <- choose_graph_segments(cds4choose, reduction_method = "UMAP", clear_cds = TRUE, return_list = FALSE)
                                                    cds_choose.no_partition.graph <- choose_graph_segments(cds4choose, reduction_method = "UMAP", clear_cds = FALSE, return_list = FALSE)
                                                    
                                                    while (ncol(cds_choose.no_partition.graph) > length(cds_choose.cell)) {
                                                        cat('[Warning] Please re-choose cells for trajectory graph. The cells you chose are more than the chosen cells before trajectory plot.','\n', sep = '')
                                                        
                                                        #cds_choose.no_partition.graph <- choose_graph_segments(cds4choose, reduction_method = "UMAP", clear_cds = TRUE, return_list = FALSE)
                                                        cds_choose.no_partition.graph <- choose_graph_segments(cds4choose, reduction_method = "UMAP", clear_cds = FALSE, return_list = FALSE)
                                                    }
                                                    while (! all(colnames(cds_choose.no_partition.graph) %in% cds_choose.cell)) {
                                                        cat('[Warning] Please re-choose cells for trajectory graph. The cells you chose are not the same or the subset of the chosen cells before trajectory plot.','\n', sep = '')
                                                        
                                                        #cds_choose.no_partition.graph <- choose_graph_segments(cds4choose, reduction_method = "UMAP", clear_cds = TRUE, return_list = FALSE)
                                                        cds_choose.no_partition.graph <- choose_graph_segments(cds4choose, reduction_method = "UMAP", clear_cds = FALSE, return_list = FALSE)
                                                    }
                                                }
                                                
                                                mapping.idx = sapply(rownames(colData(cds_choose.no_partition.graph)), FUN = function(X) which(names(cds_choose.no_partition.graph@clusters[[redDim.method[redDim]]]$clusters) %in% X))
                                                cds_choose.no_partition.graph@clusters[[redDim.method[redDim]]]$clusters = cds_choose.no_partition.graph@clusters[[redDim.method[redDim]]]$clusters[mapping.idx]
                                                cds_choose.no_partition.graph@clusters[[redDim.method[redDim]]]$clusters = factor(cds_choose.no_partition.graph@clusters[[redDim.method[redDim]]]$clusters, levels = mixedsort(as.character(unique(cds_choose.no_partition.graph@clusters[[redDim.method[redDim]]]$clusters))))
                                                colData(cds_choose.no_partition.graph)$cluster = cds_choose.no_partition.graph@clusters[[redDim.method[redDim]]]$clusters
                                                
                                                mapping.idx = sapply(rownames(colData(cds_choose.no_partition.graph)), FUN = function(X) which(names(cds_choose.no_partition.graph@clusters[[redDim.method[redDim]]]$partitions) %in% X))
                                                cds_choose.no_partition.graph@clusters[[redDim.method[redDim]]]$partitions = cds_choose.no_partition.graph@clusters[[redDim.method[redDim]]]$partitions[mapping.idx]
                                                cds_choose.no_partition.graph@clusters[[redDim.method[redDim]]]$partitions = factor(cds_choose.no_partition.graph@clusters[[redDim.method[redDim]]]$partitions, levels = mixedsort(as.character(unique(cds_choose.no_partition.graph@clusters[[redDim.method[redDim]]]$partitions))))
                                                colData(cds_choose.no_partition.graph)$partition = cds_choose.no_partition.graph@clusters[[redDim.method[redDim]]]$partitions
                                                
                                                mapping.idx = sapply(rownames(colData(cds_choose.no_partition.graph)), FUN = function(X) which(rownames(cds_choose.no_partition.graph@principal_graph_aux[[redDim.method[redDim]]]$pr_graph_cell_proj_closest_vertex) %in% X))
                                                cds_choose.no_partition.graph@principal_graph_aux[[redDim.method[redDim]]]$pr_graph_cell_proj_closest_vertex = cds_choose.no_partition.graph@principal_graph_aux[[redDim.method[redDim]]]$pr_graph_cell_proj_closest_vertex[mapping.idx, , drop = FALSE]
                                                
                                                saveRDS(cds_choose.no_partition.graph, paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',project,'.',SUB.prename,'.','trajectories','.','cds_choose.no_partition.graph','.rds', sep = ''))
                                                #cds_choose.no_partition.graph <- readRDS(paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',project,'.',SUB.prename,'.','trajectories','.','cds_choose.no_partition.graph','.rds', sep = ''))
                                                
                                                if (SUB == 1) {
                                                    cds_choose.no_partition.cell_subset.graph <- cds_choose.no_partition.graph
                                                } else if (SUB == 2) {
                                                    cds_choose.no_partition.cell_choose.graph <- cds_choose.no_partition.graph
                                                }
                                            } else {
                                                # /// use_partition = TRUE ///
                                                if (SUB == 1) {
                                                    cds_choose.graph <- cds_choose.cell_subset.graph
                                                } else if (SUB == 2) {
                                                    cds_choose.graph <- cds_choose.cell_choose.graph
                                                }
                                                
                                                # /// use_partition = FALSE ///
                                                if (SUB == 1) {
                                                    cds_choose.no_partition.graph <- cds_choose.no_partition.cell_subset.graph
                                                } else if (SUB == 2) {
                                                    cds_choose.no_partition.graph <- cds_choose.no_partition.cell_choose.graph
                                                }
                                            }
                                            
                                            #for (trj.prt in 1:1) {
                                            for (trj.prt in 1:0) {
                                                if (trj.prt == 1) {
                                                    trj.prt.name = 'with_partition'
                                                    cds_choose.plt = cds_choose.graph
                                                } else if (trj.prt == 0) {
                                                    trj.prt.name = 'without_partition'
                                                    cds_choose.plt = cds_choose.no_partition.graph
                                                }
                                                
                                                if (! is.null(cds_choose.plt)) {
                                                    
                                                    if (length(cell_type) > 0) {
                                                        for (color.plt in 1:length(cell_type)) {
                                                            legend.max_nchar = max(nchar(levels(colData(cds_choose.plt)[, cell_type[color.plt]])))
                                                            legend.multiply_factor = 0.5 + legend.max_nchar * 0.05
                                                            legend.num = length(levels(colData(cds_choose.plt)[, cell_type[color.plt]]))
                                                            legend.ncol = ceiling(legend.num / 10)
                                                            width.size = 4.075 + (legend.multiply_factor * legend.ncol)
                                                            
                                                            png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',SUB.prefolder,'/',project,'.',SUB.prename,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'graph','_',trj.prt.name,'.',gsub("\\.","_",cell_type[color.plt]),'.','disableLabel','.png', sep = ''), width = 4.075, height = 4, units = 'in', res = 600)
                                                            p <- plot_cells(cds_choose.plt, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = cell_type[color.plt], label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)
                                                            if (! is.null(color_palette)) {
                                                                cell_type.colidx <- which(paste(unique(c(cell_type[color.plt], unlist(strsplit(cell_type[color.plt], '\\.|_'))[1])),'.','color', sep = '') %in% names(color_palette))
                                                                if (length(cell_type.colidx) == 1) {
                                                                    p <- p + scale_color_manual(values = color_palette[[which(names(color_palette) %in% paste(unique(c(cell_type[color.plt], unlist(strsplit(cell_type[color.plt], '\\.|_'))[1])),'.','color', sep = ''))]])
                                                                }
                                                            }
                                                            print(p)
                                                            dev.off()
                                                            
                                                            png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',SUB.prefolder,'/',project,'.',SUB.prename,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'graph','_',trj.prt.name,'.',gsub("\\.","_",cell_type[color.plt]),'_','legend','.','disableLabel','.png', sep = ''), width = width.size, height = 4, units = 'in', res = 600)
                                                            p <- plot_cells(cds_choose.plt, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = cell_type[color.plt], label_groups_by_cluster = FALSE, label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE)
                                                            if (! is.null(color_palette)) {
                                                                cell_type.colidx <- which(paste(unique(c(cell_type[color.plt], unlist(strsplit(cell_type[color.plt], '\\.|_'))[1])),'.','color', sep = '') %in% names(color_palette))
                                                                if (length(cell_type.colidx) == 1) {
                                                                    p <- p + scale_color_manual(values = color_palette[[which(names(color_palette) %in% paste(unique(c(cell_type[color.plt], unlist(strsplit(cell_type[color.plt], '\\.|_'))[1])),'.','color', sep = ''))]])
                                                                }
                                                            }
                                                            p <- p + theme(legend.position = "right")
                                                            p$guides$colour$ncol = legend.ncol
                                                            if (project.name %in% c("GSE139944", c("GSM4150376", "GSM4150377", "GSM4150378", "GSM4150379"))) {
                                                                if (cell_type[color.plt] == "pathway_level_1") {
                                                                    p$guides$colour$title = "Pathway"
                                                                } else if (cell_type[color.plt] == "product_name") {
                                                                    p$guides$colour$title = "Treatment"
                                                                } else if (cell_type[color.plt] == "dose_character") {
                                                                    p$guides$colour$title = dose_character.unit
                                                                } else {
                                                                    p$guides$colour$title = "Treatment"
                                                                }
                                                            }
                                                            print(p)
                                                            dev.off()
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        
                        
                        ## Step 6: Order cells ****************************************************************************
                        # Order the cells in pseudotime
                        # Once we've learned a graph, we are ready to order the cells according to their progress through the developmental program. 
                        # Monocle measures this progress in pseudotime.
                        
                        # What is pseudotime?
                        # Pseudotime is a measure of how much progress an individual cell has made through a process such as cell differentiation.
                        # By ordering each cell according to its progress along a learned trajectory, Monocle alleviates the problems that arise due to asynchrony.
                        # Instead of tracking changes in expression as a function of time, Monocle tracks changes as a function of progress along the trajectory, which we term "pseudotime".
                        # Pseudotime is an abstract unit of progress: it's simply the distance between a cell and the start of the trajectory, measured along the shortest path.
                        # The trajectory's total length is defined in terms of the total amount of transcriptional change that a cell undergoes as it moves from the starting state to the end state.
                        
                        # In order to place the cells in order, we need to tell Monocle where the "beginning" of the biological process is.
                        # We do so by choosing regions of the graph that we mark as "roots" of the trajectory.
                        # In time series experiments, this can usually be accomplished by finding spots in the UMAP space that are occupied by cells from early time points
                        
                        # The black lines show the structure of the graph. 
                        # Note that the graph is not fully connected: cells in different partitions are in distinct components of the graph. 
                        # The circles with numbers in them denote special points within the graph. 
                        # Each leaf, denoted by light gray circles, corresponds to a different outcome (i.e. cell fate) of the trajectory. 
                        # Black circles indicate branch nodes, in which cells can travel to one of several outcomes.
                        
                        # Now that we have a sense of where the early cells fall, we can call order_cells(), which will calculate where each cell falls in pseudotime. 
                        # In order to do so order_cells() needs you to specify the root nodes of the trajectory graph.
                        
                        #order_cells() # reduction_method: a string specifying the reduced dimension method to use when ordering cells. Currently only "UMAP" is supported.
                        
                        # Note that some of the cells are gray. This means they have infinite pseudotime, because they were not reachable from the root nodes that were picked. 
                        # In general, any cell on a parition that lacks a root node will be assigned an infinite pseudotime. 
                        # In general, you should choose at least one root per partition.
                        
                        # /// manually picking it ///
                        # If you don't provide them as an argument, it will launch a graphical user interface for selecting one or more root nodes.
                        #cds <- order_cells(cds)
                        
                        # /// to specify the root of the trajectory programmatically ///
                        # The function below does so by first grouping the cells according to which trajectory graph node they are nearest to. 
                        # Then, it calculates what fraction of the cells at each node come from the earliest time point. 
                        # Then it picks the node that is most heavily occupied by early cells and returns that as the root.
                        
                        # a helper function to identify the root principal points:
                        get_earliest_principal_node <- function(cds, cds.graph_node.idx = NULL, time = NULL, time_bin = NULL, root.num = 1) {
                            closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
                            closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
                            
                            if (! is.null(time_bin)) {
                                cell_ids <- which(colData(cds)[, time] %in% time_bin)
                                if (length(cell_ids) > 0) {
                                    #root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
                                    
                                    node.idx <- as.numeric(names(
                                        head(sort( # Then it picks the node that is most heavily occupied by early cells and returns that as the root.
                                            table( # Then, it calculates what fraction of the cells at each node come from the earliest time point. 
                                                closest_vertex[cell_ids,] # The function below does so by first grouping the cells according to which trajectory graph node they are nearest to. 
                                            ), 
                                            decreasing = TRUE), root.num)
                                    ))
                                    
                                    if (! is.null(cds.graph_node.idx)) {
                                        node.idx <- sapply(node.idx, FUN = function(X) which(cds.graph_node.idx %in% X))
                                    }
                                    
                                    root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[node.idx]
                                } else {
                                    root_pr_nodes <- NULL
                                }
                            } else {
                                root_pr_nodes <- NULL
                            }
                            
                            root_pr_nodes
                        }
                        
                        if (project.name == "tutorial") {
                            time = "embryo.time.bin"
                            level.num = 2
                        } else {
                            time = time.input
                            level.num = 1
                        }
                        
                        if (! is.null(time)) {
                            
                            for (redDim in 1:length(redDim.method)) {
                                
                                if (redDim.method[redDim] == "UMAP") {
                                    
                                } else if (redDim.method[redDim] == "tSNE") {
                                    next # next halts the processing of the current iteration and advances the looping index.
                                }
                                
                                cds = cds.origin.processed
                                cds.no_partition = cds.no_partition.origin.processed
                                
                                
                                dim.folder = paste(eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))),'D', sep = '')
                                dim.finalX = eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) - (2-1)
                                
                                if (dim.finalX >= 1) {
                                    
                                    for (dim.comb in 1:dim.finalX) {
                                        dim.x <- dim.comb
                                        dim.y <- dim.x + 1
                                        
                                        if (dim.finalX > 1) {
                                            if (dim.x == 1) {
                                                comb.run = 2
                                            } else {
                                                comb.run = 1
                                            }
                                        } else {
                                            comb.run = 1
                                        }
                                        
                                        for (dim.comb.run in 1:comb.run) {
                                            if (dim.comb.run > 1) {
                                                dim.y <- dim.y + 1
                                            }
                                            
                                            if (eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) == 2) {
                                                dim.xy <- paste('D',dim.x,'-','D',dim.y, sep = '')
                                            } else {
                                                dim.xy <- paste('D',dim.x,'-','D',dim.y, sep = '')
                                            }
                                            #print(dim.xy)
                                            
                                            mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder, sep = '')
                                            subDir <- dim.xy
                                            if (! file.exists(file.path(mainDir, subDir))) {
                                                dir.create(file.path(mainDir, subDir))
                                            }
                                            
                                            #for (trj.prt in 1:1) {
                                            for (trj.prt in 1:0) {
                                                if (trj.prt == 1) {
                                                    trj.prt.name = 'with_partition'
                                                    cds.plt = cds
                                                } else if (trj.prt == 0) {
                                                    trj.prt.name = 'without_partition'
                                                    cds.plt = cds.no_partition
                                                }
                                                
                                                if (! is.null(cds.plt)) {
                                                    
                                                    if (project.name == "tutorial") {
                                                        width.size = 5.4625
                                                        height.size = 4
                                                    } else {
                                                        legend.max_nchar = max(nchar(levels(colData(cds.plt)[, time])))
                                                        legend.multiply_factor = 0.5 + legend.max_nchar * 0.05
                                                        legend.num = length(levels(colData(cds.plt)[, time]))
                                                        legend.ncol = ceiling(legend.num / 10)
                                                        width.size = 4.075 + (legend.multiply_factor * legend.ncol)
                                                        height.size = 4
                                                    }
                                                    
                                                    png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',project,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'graph','_',trj.prt.name,'.',gsub("\\.","_",time),'.','disableRootLabel','.','with_trajectory','.png', sep = ''), width = width.size, height = height.size, units = 'in', res = 600)
                                                    p <- plot_cells(cds.plt, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = time, label_cell_groups = FALSE, group_label_size = 3, label_leaves = FALSE, label_branch_points = FALSE, label_roots = FALSE, show_trajectory_graph = TRUE, graph_label_size = 3)
                                                    if (! is.null(color_palette)) {
                                                        time.colidx <- which(paste(unique(c(gsub("\\.","_",time), unlist(strsplit(gsub("\\.","_",time), '\\.|_'))[2])),'.','color', sep = '') %in% names(color_palette))
                                                        if (length(time.colidx) == 1) {
                                                            p <- p + scale_color_manual(values = color_palette[[which(names(color_palette) %in% paste(unique(c(gsub("\\.","_",time), unlist(strsplit(gsub("\\.","_",time), '\\.|_'))[2])),'.','color', sep = ''))]])
                                                        }
                                                    }
                                                    p$guides$colour$ncol = legend.ncol
                                                    if (project.name %in% c("GSE139944", c("GSM4150376", "GSM4150377", "GSM4150378", "GSM4150379"))) {
                                                        p$guides$colour$title = dose_character.unit
                                                    }
                                                    print(p)
                                                    dev.off()
                                                    
                                                    png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',project,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'graph','_',trj.prt.name,'.',gsub("\\.","_",time),'.','disableRootLabel','.','without_trajectory','.png', sep = ''), width = width.size, height = height.size, units = 'in', res = 600)
                                                    p <- plot_cells(cds.plt, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = time, label_cell_groups = FALSE, group_label_size = 3, label_leaves = FALSE, label_branch_points = FALSE, label_roots = FALSE, show_trajectory_graph = FALSE, graph_label_size = 3)
                                                    if (! is.null(color_palette)) {
                                                        time.colidx <- which(paste(unique(c(gsub("\\.","_",time), unlist(strsplit(gsub("\\.","_",time), '\\.|_'))[2])),'.','color', sep = '') %in% names(color_palette))
                                                        if (length(time.colidx) == 1) {
                                                            p <- p + scale_color_manual(values = color_palette[[which(names(color_palette) %in% paste(unique(c(gsub("\\.","_",time), unlist(strsplit(gsub("\\.","_",time), '\\.|_'))[2])),'.','color', sep = ''))]])
                                                        }
                                                    }
                                                    p$guides$colour$ncol = legend.ncol
                                                    if (project.name %in% c("GSE139944", c("GSM4150376", "GSM4150377", "GSM4150378", "GSM4150379"))) {
                                                        p$guides$colour$title = dose_character.unit
                                                    }
                                                    print(p)
                                                    dev.off()
                                                    
                                                    png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',project,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'graph','_',trj.prt.name,'.',gsub("\\.","_",time),'.','branch','.png', sep = ''), width = width.size, height = height.size, units = 'in', res = 600)
                                                    p <- plot_cells(cds.plt, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = time, label_cell_groups = FALSE, group_label_size = 3, label_leaves = FALSE, label_branch_points = TRUE, label_roots = TRUE, show_trajectory_graph = TRUE, graph_label_size = 3)
                                                    if (! is.null(color_palette)) {
                                                        time.colidx <- which(paste(unique(c(gsub("\\.","_",time), unlist(strsplit(gsub("\\.","_",time), '\\.|_'))[2])),'.','color', sep = '') %in% names(color_palette))
                                                        if (length(time.colidx) == 1) {
                                                            p <- p + scale_color_manual(values = color_palette[[which(names(color_palette) %in% paste(unique(c(gsub("\\.","_",time), unlist(strsplit(gsub("\\.","_",time), '\\.|_'))[2])),'.','color', sep = ''))]])
                                                        }
                                                    }
                                                    p$guides$colour$ncol = legend.ncol
                                                    if (project.name %in% c("GSE139944", c("GSM4150376", "GSM4150377", "GSM4150378", "GSM4150379"))) {
                                                        p$guides$colour$title = dose_character.unit
                                                    }
                                                    print(p)
                                                    dev.off()
                                                    
                                                    png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',project,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'graph','_',trj.prt.name,'.',gsub("\\.","_",time),'.','leaves','.png', sep = ''), width = width.size, height = height.size, units = 'in', res = 600)
                                                    p <- plot_cells(cds.plt, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = time, label_cell_groups = FALSE, group_label_size = 3, label_leaves = TRUE, label_branch_points = FALSE, label_roots = TRUE, show_trajectory_graph = TRUE, graph_label_size = 3)
                                                    if (! is.null(color_palette)) {
                                                        time.colidx <- which(paste(unique(c(gsub("\\.","_",time), unlist(strsplit(gsub("\\.","_",time), '\\.|_'))[2])),'.','color', sep = '') %in% names(color_palette))
                                                        if (length(time.colidx) == 1) {
                                                            p <- p + scale_color_manual(values = color_palette[[which(names(color_palette) %in% paste(unique(c(gsub("\\.","_",time), unlist(strsplit(gsub("\\.","_",time), '\\.|_'))[2])),'.','color', sep = ''))]])
                                                        }
                                                    }
                                                    p$guides$colour$ncol = legend.ncol
                                                    if (project.name %in% c("GSE139944", c("GSM4150376", "GSM4150377", "GSM4150378", "GSM4150379"))) {
                                                        p$guides$colour$title = dose_character.unit
                                                    }
                                                    print(p)
                                                    dev.off()
                                                    
                                                    png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',project,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'graph','_',trj.prt.name,'.',gsub("\\.","_",time),'.','nodes','.png', sep = ''), width = width.size, height = height.size, units = 'in', res = 600)
                                                    p <- plot_cells(cds.plt, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = time, label_cell_groups = FALSE, group_label_size = 3, label_leaves = TRUE, label_branch_points = TRUE, label_roots = TRUE, show_trajectory_graph = TRUE, graph_label_size = 3)
                                                    if (! is.null(color_palette)) {
                                                        time.colidx <- which(paste(unique(c(gsub("\\.","_",time), unlist(strsplit(gsub("\\.","_",time), '\\.|_'))[2])),'.','color', sep = '') %in% names(color_palette))
                                                        if (length(time.colidx) == 1) {
                                                            p <- p + scale_color_manual(values = color_palette[[which(names(color_palette) %in% paste(unique(c(gsub("\\.","_",time), unlist(strsplit(gsub("\\.","_",time), '\\.|_'))[2])),'.','color', sep = ''))]])
                                                        }
                                                    }
                                                    p$guides$colour$ncol = legend.ncol
                                                    if (project.name %in% c("GSE139944", c("GSM4150376", "GSM4150377", "GSM4150378", "GSM4150379"))) {
                                                        p$guides$colour$title = dose_character.unit
                                                    }
                                                    print(p)
                                                    dev.off()
                                                    
                                                    png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',project,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'graph','_',trj.prt.name,'.',gsub("\\.","_",time),'.','label_cell_groups','.png', sep = ''), width = 4.075, height = 4, units = 'in', res = 600)
                                                    p <- plot_cells(cds.plt, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = time, label_cell_groups = TRUE, group_label_size = 3, label_leaves = FALSE, label_branch_points = FALSE, label_roots = FALSE, show_trajectory_graph = TRUE, graph_label_size = 3)
                                                    if (! is.null(color_palette)) {
                                                        time.colidx <- which(paste(unique(c(gsub("\\.","_",time), unlist(strsplit(gsub("\\.","_",time), '\\.|_'))[2])),'.','color', sep = '') %in% names(color_palette))
                                                        if (length(time.colidx) == 1) {
                                                            p <- p + scale_color_manual(values = color_palette[[which(names(color_palette) %in% paste(unique(c(gsub("\\.","_",time), unlist(strsplit(gsub("\\.","_",time), '\\.|_'))[2])),'.','color', sep = ''))]])
                                                        }
                                                    }
                                                    p$guides$colour$ncol = legend.ncol
                                                    if (project.name %in% c("GSE139944", c("GSM4150376", "GSM4150377", "GSM4150378", "GSM4150379"))) {
                                                        p$guides$colour$title = dose_character.unit
                                                    }
                                                    print(p)
                                                    dev.off()
                                                }
                                            }
                                            
                                            if (project.name %in% c("GSE139944", c("GSM4150376", "GSM4150377", "GSM4150378", "GSM4150379"))) {
                                                pseudo.name = "pseudodose"
                                            } else {
                                                pseudo.name = "pseudotime"
                                            }
                                            
                                            # /// all ///
                                            # /// use_partition = TRUE ///
                                            cds.level = names(which(sapply(levels(colData(cds)[, time]), FUN = function(X) X %in% unique(colData(cds)[, time]))))
                                            if (length(cds.level) > level.num) {
                                                time_bin = head(cds.level, level.num)
                                            } else {
                                                time_bin = head(cds.level, length(cds.level))
                                            }
                                            
                                            root.num = 2
                                            
                                            cds <- order_cells(cds, reduction_method = "UMAP", root_pr_nodes = get_earliest_principal_node(cds, time = time, time_bin = time_bin, root.num = root.num))
                                            
                                            if (project.name %in% c("GSM4150377", "GSM4150378")) {
                                                pData(cds)$Pseudotime = cds@principal_graph_aux[["UMAP"]]$pseudotime
                                                pData(cds)$Pseudotime_bin = cut(pData(cds)$Pseudotime, breaks = unique(quantile(pData(cds)$Pseudotime, seq(0, 1, 0.1))), labels = FALSE)
                                            }
                                            
                                            cds.origin.processed = cds
                                            
                                            # /// use_partition = FALSE ///
                                            cds.no_partition.level = names(which(sapply(levels(colData(cds.no_partition)[, time]), FUN = function(X) X %in% unique(colData(cds.no_partition)[, time]))))
                                            if (length(cds.no_partition.level) > level.num) {
                                                time_bin = head(cds.no_partition.level, level.num)
                                            } else {
                                                time_bin = head(cds.no_partition.level, length(cds.no_partition.level))
                                            }
                                            
                                            root.num = 2
                                            
                                            cds.no_partition <- order_cells(cds.no_partition, reduction_method = "UMAP", root_pr_nodes = get_earliest_principal_node(cds.no_partition, time = time, time_bin = time_bin, root.num = root.num))
                                            
                                            cds.no_partition.origin.processed = cds.no_partition
                                            
                                            
                                            #for (trj.prt in 1:1) {
                                            for (trj.prt in 1:0) {
                                                if (trj.prt == 1) {
                                                    trj.prt.name = 'with_partition'
                                                    cds.plt = cds
                                                } else if (trj.prt == 0) {
                                                    trj.prt.name = 'without_partition'
                                                    cds.plt = cds.no_partition
                                                }
                                                
                                                if (! is.null(cds.plt)) {
                                                    
                                                    png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',project,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'graph','_',trj.prt.name,'.','order','.',pseudo.name,'.png', sep = ''), width = 5.165, height = 4, units = 'in', res = 600)
                                                    p <- plot_cells(cds.plt, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = "pseudotime", label_cell_groups = FALSE, group_label_size = 3, label_leaves = FALSE, label_branch_points = FALSE, label_roots = TRUE, show_trajectory_graph = TRUE, graph_label_size = 3)
                                                    print(p)
                                                    dev.off()
                                                    
                                                    png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',project,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'graph','_',trj.prt.name,'.','order','.',pseudo.name,'.','disableRootLabel','.','with_trajectory','.png', sep = ''), width = 5.165, height = 4, units = 'in', res = 600)
                                                    p <- plot_cells(cds.plt, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = "pseudotime", label_cell_groups = FALSE, group_label_size = 3, label_leaves = FALSE, label_branch_points = FALSE, label_roots = FALSE, show_trajectory_graph = TRUE, graph_label_size = 3)
                                                    print(p)
                                                    dev.off()
                                                    
                                                    png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',project,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'graph','_',trj.prt.name,'.','order','.',pseudo.name,'.','disableRootLabel','.','without_trajectory','.png', sep = ''), width = 5.165, height = 4, units = 'in', res = 600)
                                                    p <- plot_cells(cds.plt, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = "pseudotime", label_cell_groups = FALSE, group_label_size = 3, label_leaves = FALSE, label_branch_points = FALSE, label_roots = FALSE, show_trajectory_graph = FALSE, graph_label_size = 3)
                                                    print(p)
                                                    dev.off()
                                                    
                                                    # /// partitions ///
                                                    cds_pr <- partitions(cds.plt, reduction_method = redDim.method[redDim]) # Note that we could easily do this on a per-partition basis by first grouping the cells by partition using the partitions() function. This would result in all cells being assigned a finite pseudotime.
                                                    
                                                    root_pr_nodes.all <- c()
                                                    for (pr in 1:length(levels(cds_pr))) {
                                                        
                                                        sub.idx = which(cds_pr == levels(cds_pr)[pr])
                                                        
                                                        if (length(sub.idx) > 0) {
                                                            
                                                            cds_pr_sub = cds.plt[, sub.idx]
                                                            
                                                            cds_pr_sub.level = names(which(sapply(levels(colData(cds_pr_sub)[, time]), FUN = function(X) X %in% unique(colData(cds_pr_sub)[, time]))))
                                                            if (length(cds_pr_sub.level) > level.num) {
                                                                time_bin = head(cds_pr_sub.level, level.num)
                                                            } else {
                                                                time_bin = head(cds_pr_sub.level, length(cds_pr_sub.level))
                                                            }
                                                            
                                                            root.num = 1
                                                            if (length(sub.idx) < root.num) {
                                                                root.num = length(sub.idx)
                                                            }
                                                            
                                                            root_pr_nodes = get_earliest_principal_node(cds_pr_sub, time = time, time_bin = time_bin, root.num = root.num)
                                                            
                                                            root_pr_nodes.all <- c(root_pr_nodes.all, root_pr_nodes)
                                                        }
                                                    }
                                                    root_pr_nodes.all <- unique(root_pr_nodes.all)
                                                    cds_pr <- order_cells(cds.plt, reduction_method = "UMAP", root_pr_nodes = root_pr_nodes.all)
                                                    
                                                    png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',project,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'graph','_',trj.prt.name,'.','order','.',pseudo.name,'.','partition','.png', sep = ''), width = 5.165, height = 4, units = 'in', res = 600)
                                                    p <- plot_cells(cds_pr, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = "pseudotime", label_cell_groups = FALSE, group_label_size = 3, label_leaves = FALSE, label_branch_points = FALSE, label_roots = TRUE, show_trajectory_graph = TRUE, graph_label_size = 3)
                                                    print(p)
                                                    dev.off()
                                                    
                                                    png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',project,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'graph','_',trj.prt.name,'.','order','.',pseudo.name,'.','partition','.','disableRootLabel','.','with_trajectory','.png', sep = ''), width = 5.165, height = 4, units = 'in', res = 600)
                                                    p <- plot_cells(cds_pr, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = "pseudotime", label_cell_groups = FALSE, group_label_size = 3, label_leaves = FALSE, label_branch_points = FALSE, label_roots = FALSE, show_trajectory_graph = TRUE, graph_label_size = 3)
                                                    print(p)
                                                    dev.off()
                                                    
                                                    png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',project,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'graph','_',trj.prt.name,'.','order','.',pseudo.name,'.','partition','.','disableRootLabel','.','without_trajectory','.png', sep = ''), width = 5.165, height = 4, units = 'in', res = 600)
                                                    p <- plot_cells(cds_pr, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = "pseudotime", label_cell_groups = FALSE, group_label_size = 3, label_leaves = FALSE, label_branch_points = FALSE, label_roots = FALSE, show_trajectory_graph = FALSE, graph_label_size = 3)
                                                    print(p)
                                                    dev.off()
                                                }
                                            }
                                            
                                            
                                            # /// subset ///
                                            for (SUB in SUB.run) {
                                                if (SUB == 0) {
                                                    next # next halts the processing of the current iteration and advances the looping index.
                                                } else if (SUB == 1) {
                                                    SUB.prefolder = 'cell_subset'
                                                    SUB.prename = 'cell_subset'
                                                } else if (SUB == 2) {
                                                    if (length(cell_subset) == 0 & cell_choose == TRUE) {
                                                        SUB.prefolder = 'choose'
                                                        SUB.prename = 'choose'
                                                    } else if (length(cell_subset) > 0 & cell_choose == TRUE) {
                                                        SUB.prefolder = paste('cell_subset','/','choose', sep = '')
                                                        SUB.prename = paste('cell_subset','.','choose', sep = '')
                                                    }
                                                }
                                                
                                                ## Step 6: Order cells ****************************************************************************
                                                # /// use_partition = TRUE ///
                                                if (SUB == 1) {
                                                    cds_choose <- cds_choose.cell_subset
                                                } else if (SUB == 2) {
                                                    cds_choose <- cds_choose.cell_choose
                                                }
                                                
                                                # /// use_partition = FALSE ///
                                                if (SUB == 1) {
                                                    cds_choose.no_partition <- cds_choose.no_partition.cell_subset
                                                } else if (SUB == 2) {
                                                    cds_choose.no_partition <- cds_choose.no_partition.cell_choose
                                                }
                                                
                                                #for (trj.prt in 1:1) {
                                                for (trj.prt in 1:0) {
                                                    if (trj.prt == 1) {
                                                        trj.prt.name = 'with_partition'
                                                        cds_choose.plt = cds_choose
                                                    } else if (trj.prt == 0) {
                                                        trj.prt.name = 'without_partition'
                                                        cds_choose.plt = cds_choose.no_partition
                                                    }
                                                    
                                                    if (! is.null(cds_choose.plt)) {
                                                        
                                                        if (project.name == "tutorial") {
                                                            width.size = 5.4625
                                                            height.size = 4
                                                        } else {
                                                            legend.max_nchar = max(nchar(levels(colData(cds_choose.plt)[, time])))
                                                            legend.multiply_factor = 0.5 + legend.max_nchar * 0.05
                                                            legend.num = length(levels(colData(cds_choose.plt)[, time]))
                                                            legend.ncol = ceiling(legend.num / 10)
                                                            width.size = 4.075 + (legend.multiply_factor * legend.ncol)
                                                            height.size = 4
                                                        }
                                                        
                                                        png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',SUB.prefolder,'/',project,'.',SUB.prename,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'graph','_',trj.prt.name,'.',gsub("\\.","_",time),'.','disableRootLabel','.','with_trajectory','.png', sep = ''), width = width.size, height = height.size, units = 'in', res = 600)
                                                        p <- plot_cells(cds_choose.plt, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = time, label_cell_groups = FALSE, group_label_size = 3, label_leaves = FALSE, label_branch_points = FALSE, label_roots = FALSE, show_trajectory_graph = TRUE, graph_label_size = 3)
                                                        if (! is.null(color_palette)) {
                                                            time.colidx <- which(paste(unique(c(gsub("\\.","_",time), unlist(strsplit(gsub("\\.","_",time), '\\.|_'))[2])),'.','color', sep = '') %in% names(color_palette))
                                                            if (length(time.colidx) == 1) {
                                                                p <- p + scale_color_manual(values = color_palette[[which(names(color_palette) %in% paste(unique(c(gsub("\\.","_",time), unlist(strsplit(gsub("\\.","_",time), '\\.|_'))[2])),'.','color', sep = ''))]])
                                                            }
                                                        }
                                                        p$guides$colour$ncol = legend.ncol
                                                        if (project.name %in% c("GSE139944", c("GSM4150376", "GSM4150377", "GSM4150378", "GSM4150379"))) {
                                                            p$guides$colour$title = dose_character.unit
                                                        }
                                                        print(p)
                                                        dev.off()
                                                        
                                                        png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',SUB.prefolder,'/',project,'.',SUB.prename,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'graph','_',trj.prt.name,'.',gsub("\\.","_",time),'.','disableRootLabel','.','without_trajectory','.png', sep = ''), width = width.size, height = height.size, units = 'in', res = 600)
                                                        p <- plot_cells(cds_choose.plt, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = time, label_cell_groups = FALSE, group_label_size = 3, label_leaves = FALSE, label_branch_points = FALSE, label_roots = FALSE, show_trajectory_graph = FALSE, graph_label_size = 3)
                                                        if (! is.null(color_palette)) {
                                                            time.colidx <- which(paste(unique(c(gsub("\\.","_",time), unlist(strsplit(gsub("\\.","_",time), '\\.|_'))[2])),'.','color', sep = '') %in% names(color_palette))
                                                            if (length(time.colidx) == 1) {
                                                                p <- p + scale_color_manual(values = color_palette[[which(names(color_palette) %in% paste(unique(c(gsub("\\.","_",time), unlist(strsplit(gsub("\\.","_",time), '\\.|_'))[2])),'.','color', sep = ''))]])
                                                            }
                                                        }
                                                        p$guides$colour$ncol = legend.ncol
                                                        if (project.name %in% c("GSE139944", c("GSM4150376", "GSM4150377", "GSM4150378", "GSM4150379"))) {
                                                            p$guides$colour$title = dose_character.unit
                                                        }
                                                        print(p)
                                                        dev.off()
                                                        
                                                        png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',SUB.prefolder,'/',project,'.',SUB.prename,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'graph','_',trj.prt.name,'.',gsub("\\.","_",time),'.','branch','.png', sep = ''), width = width.size, height = height.size, units = 'in', res = 600)
                                                        p <- plot_cells(cds_choose.plt, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = time, label_cell_groups = FALSE, group_label_size = 3, label_leaves = FALSE, label_branch_points = TRUE, label_roots = TRUE, show_trajectory_graph = TRUE, graph_label_size = 3)
                                                        if (! is.null(color_palette)) {
                                                            time.colidx <- which(paste(unique(c(gsub("\\.","_",time), unlist(strsplit(gsub("\\.","_",time), '\\.|_'))[2])),'.','color', sep = '') %in% names(color_palette))
                                                            if (length(time.colidx) == 1) {
                                                                p <- p + scale_color_manual(values = color_palette[[which(names(color_palette) %in% paste(unique(c(gsub("\\.","_",time), unlist(strsplit(gsub("\\.","_",time), '\\.|_'))[2])),'.','color', sep = ''))]])
                                                            }
                                                        }
                                                        p$guides$colour$ncol = legend.ncol
                                                        if (project.name %in% c("GSE139944", c("GSM4150376", "GSM4150377", "GSM4150378", "GSM4150379"))) {
                                                            p$guides$colour$title = dose_character.unit
                                                        }
                                                        print(p)
                                                        dev.off()
                                                        
                                                        png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',SUB.prefolder,'/',project,'.',SUB.prename,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'graph','_',trj.prt.name,'.',gsub("\\.","_",time),'.','leaves','.png', sep = ''), width = width.size, height = height.size, units = 'in', res = 600)
                                                        p <- plot_cells(cds_choose.plt, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = time, label_cell_groups = FALSE, group_label_size = 3, label_leaves = TRUE, label_branch_points = FALSE, label_roots = TRUE, show_trajectory_graph = TRUE, graph_label_size = 3)
                                                        if (! is.null(color_palette)) {
                                                            time.colidx <- which(paste(unique(c(gsub("\\.","_",time), unlist(strsplit(gsub("\\.","_",time), '\\.|_'))[2])),'.','color', sep = '') %in% names(color_palette))
                                                            if (length(time.colidx) == 1) {
                                                                p <- p + scale_color_manual(values = color_palette[[which(names(color_palette) %in% paste(unique(c(gsub("\\.","_",time), unlist(strsplit(gsub("\\.","_",time), '\\.|_'))[2])),'.','color', sep = ''))]])
                                                            }
                                                        }
                                                        p$guides$colour$ncol = legend.ncol
                                                        if (project.name %in% c("GSE139944", c("GSM4150376", "GSM4150377", "GSM4150378", "GSM4150379"))) {
                                                            p$guides$colour$title = dose_character.unit
                                                        }
                                                        print(p)
                                                        dev.off()
                                                        
                                                        png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',SUB.prefolder,'/',project,'.',SUB.prename,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'graph','_',trj.prt.name,'.',gsub("\\.","_",time),'.','nodes','.png', sep = ''), width = width.size, height = height.size, units = 'in', res = 600)
                                                        p <- plot_cells(cds_choose.plt, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = time, label_cell_groups = FALSE, group_label_size = 3, label_leaves = TRUE, label_branch_points = TRUE, label_roots = TRUE, show_trajectory_graph = TRUE, graph_label_size = 3)
                                                        if (! is.null(color_palette)) {
                                                            time.colidx <- which(paste(unique(c(gsub("\\.","_",time), unlist(strsplit(gsub("\\.","_",time), '\\.|_'))[2])),'.','color', sep = '') %in% names(color_palette))
                                                            if (length(time.colidx) == 1) {
                                                                p <- p + scale_color_manual(values = color_palette[[which(names(color_palette) %in% paste(unique(c(gsub("\\.","_",time), unlist(strsplit(gsub("\\.","_",time), '\\.|_'))[2])),'.','color', sep = ''))]])
                                                            }
                                                        }
                                                        p$guides$colour$ncol = legend.ncol
                                                        if (project.name %in% c("GSE139944", c("GSM4150376", "GSM4150377", "GSM4150378", "GSM4150379"))) {
                                                            p$guides$colour$title = dose_character.unit
                                                        }
                                                        print(p)
                                                        dev.off()
                                                        
                                                        png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',SUB.prefolder,'/',project,'.',SUB.prename,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'graph','_',trj.prt.name,'.',gsub("\\.","_",time),'.','label_cell_groups','.png', sep = ''), width = 4.075, height = 4, units = 'in', res = 600)
                                                        p <- plot_cells(cds_choose.plt, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = time, label_cell_groups = TRUE, group_label_size = 3, label_leaves = FALSE, label_branch_points = FALSE, label_roots = FALSE, show_trajectory_graph = TRUE, graph_label_size = 3)
                                                        if (! is.null(color_palette)) {
                                                            time.colidx <- which(paste(unique(c(gsub("\\.","_",time), unlist(strsplit(gsub("\\.","_",time), '\\.|_'))[2])),'.','color', sep = '') %in% names(color_palette))
                                                            if (length(time.colidx) == 1) {
                                                                p <- p + scale_color_manual(values = color_palette[[which(names(color_palette) %in% paste(unique(c(gsub("\\.","_",time), unlist(strsplit(gsub("\\.","_",time), '\\.|_'))[2])),'.','color', sep = ''))]])
                                                            }
                                                        }
                                                        p$guides$colour$ncol = legend.ncol
                                                        if (project.name %in% c("GSE139944", c("GSM4150376", "GSM4150377", "GSM4150378", "GSM4150379"))) {
                                                            p$guides$colour$title = dose_character.unit
                                                        }
                                                        print(p)
                                                        dev.off()
                                                    }
                                                }
                                                
                                                if (project.name %in% c("GSE139944", c("GSM4150376", "GSM4150377", "GSM4150378", "GSM4150379"))) {
                                                    pseudo.name = "pseudodose"
                                                } else {
                                                    pseudo.name = "pseudotime"
                                                }
                                                
                                                # /// all ///
                                                # /// use_partition = TRUE ///
                                                cds_choose.level = names(which(sapply(levels(colData(cds_choose)[, time]), FUN = function(X) X %in% unique(colData(cds_choose)[, time]))))
                                                if (length(cds_choose.level) > level.num) {
                                                    time_bin = head(cds_choose.level, level.num)
                                                } else {
                                                    time_bin = head(cds_choose.level, length(cds_choose.level))
                                                }
                                                
                                                root.num = 2
                                                
                                                cds.graph_node <- igraph::V(principal_graph(cds)[["UMAP"]])$name
                                                cds_choose.graph_node <- igraph::V(principal_graph(cds_choose)[["UMAP"]])$name
                                                cds.graph_node.idx <- sapply(cds_choose.graph_node, FUN = function(X) which(cds.graph_node %in% X))
                                                
                                                cds_choose <- order_cells(cds_choose, reduction_method = "UMAP", root_pr_nodes = get_earliest_principal_node(cds_choose, cds.graph_node.idx = cds.graph_node.idx, time = time, time_bin = time_bin, root.num = root.num))
                                                
                                                if (project.name %in% c("GSM4150377", "GSM4150378")) {
                                                    pData(cds_choose)$Pseudotime = cds_choose@principal_graph_aux[["UMAP"]]$pseudotime
                                                    pData(cds_choose)$Pseudotime_bin = cut(pData(cds_choose)$Pseudotime, breaks = unique(quantile(pData(cds_choose)$Pseudotime, seq(0, 1, 0.1))), labels = FALSE)
                                                }
                                                
                                                if (SUB == 1) {
                                                    cds_choose.cell_subset <- cds_choose
                                                } else if (SUB == 2) {
                                                    cds_choose.cell_choose <- cds_choose
                                                }
                                                
                                                # /// use_partition = FALSE ///
                                                cds_choose.no_partition.level = names(which(sapply(levels(colData(cds_choose.no_partition)[, time]), FUN = function(X) X %in% unique(colData(cds_choose.no_partition)[, time]))))
                                                if (length(cds_choose.no_partition.level) > level.num) {
                                                    time_bin = head(cds_choose.no_partition.level, level.num)
                                                } else {
                                                    time_bin = head(cds_choose.no_partition.level, length(cds_choose.no_partition.level))
                                                }
                                                
                                                root.num = 2
                                                
                                                cds.no_partition.graph_node <- igraph::V(principal_graph(cds.no_partition)[["UMAP"]])$name
                                                cds_choose.no_partition.graph_node <- igraph::V(principal_graph(cds_choose.no_partition)[["UMAP"]])$name
                                                cds.no_partition.graph_node.idx <- sapply(cds_choose.no_partition.graph_node, FUN = function(X) which(cds.no_partition.graph_node %in% X))
                                                
                                                cds_choose.no_partition <- order_cells(cds_choose.no_partition, reduction_method = "UMAP", root_pr_nodes = get_earliest_principal_node(cds_choose.no_partition, cds.graph_node.idx = cds.no_partition.graph_node.idx, time = time, time_bin = time_bin, root.num = root.num))
                                                
                                                if (SUB == 1) {
                                                    cds_choose.no_partition.cell_subset <- cds_choose.no_partition
                                                } else if (SUB == 2) {
                                                    cds_choose.no_partition.cell_choose <- cds_choose.no_partition
                                                }
                                                
                                                
                                                #for (trj.prt in 1:1) {
                                                for (trj.prt in 1:0) {
                                                    if (trj.prt == 1) {
                                                        trj.prt.name = 'with_partition'
                                                        cds.plt = cds
                                                        cds_choose.plt = cds_choose
                                                    } else if (trj.prt == 0) {
                                                        trj.prt.name = 'without_partition'
                                                        cds.plt = cds.no_partition
                                                        cds_choose.plt = cds_choose.no_partition
                                                    }
                                                    
                                                    if (! is.null(cds_choose.plt)) {
                                                        
                                                        png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',SUB.prefolder,'/',project,'.',SUB.prename,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'graph','_',trj.prt.name,'.','order','.',pseudo.name,'.png', sep = ''), width = 5.165, height = 4, units = 'in', res = 600)
                                                        p <- plot_cells(cds_choose.plt, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = "pseudotime", label_cell_groups = FALSE, group_label_size = 3, label_leaves = FALSE, label_branch_points = FALSE, label_roots = TRUE, show_trajectory_graph = TRUE, graph_label_size = 3)
                                                        print(p)
                                                        dev.off()
                                                        
                                                        png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',SUB.prefolder,'/',project,'.',SUB.prename,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'graph','_',trj.prt.name,'.','order','.',pseudo.name,'.','disableRootLabel','.','with_trajectory','.png', sep = ''), width = 5.165, height = 4, units = 'in', res = 600)
                                                        p <- plot_cells(cds_choose.plt, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = "pseudotime", label_cell_groups = FALSE, group_label_size = 3, label_leaves = FALSE, label_branch_points = FALSE, label_roots = FALSE, show_trajectory_graph = TRUE, graph_label_size = 3)
                                                        print(p)
                                                        dev.off()
                                                        
                                                        png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',SUB.prefolder,'/',project,'.',SUB.prename,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'graph','_',trj.prt.name,'.','order','.',pseudo.name,'.','disableRootLabel','.','without_trajectory','.png', sep = ''), width = 5.165, height = 4, units = 'in', res = 600)
                                                        p <- plot_cells(cds_choose.plt, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = "pseudotime", label_cell_groups = FALSE, group_label_size = 3, label_leaves = FALSE, label_branch_points = FALSE, label_roots = FALSE, show_trajectory_graph = FALSE, graph_label_size = 3)
                                                        print(p)
                                                        dev.off()
                                                        
                                                        # /// partitions ///
                                                        cds_choose_pr <- partitions(cds_choose.plt, reduction_method = redDim.method[redDim]) # Note that we could easily do this on a per-partition basis by first grouping the cells by partition using the partitions() function. This would result in all cells being assigned a finite pseudotime.
                                                        
                                                        root_pr_nodes.all <- c()
                                                        for (pr in 1:length(levels(cds_choose_pr))) {
                                                            
                                                            sub.idx = which(cds_choose_pr == levels(cds_choose_pr)[pr])
                                                            
                                                            if (length(sub.idx) > 0) {
                                                                
                                                                cds_choose_pr_sub = cds_choose.plt[, sub.idx]
                                                                
                                                                cds_choose_pr_sub.level = names(which(sapply(levels(colData(cds_choose_pr_sub)[, time]), FUN = function(X) X %in% unique(colData(cds_choose_pr_sub)[, time]))))
                                                                if (length(cds_choose_pr_sub.level) > level.num) {
                                                                    time_bin = head(cds_choose_pr_sub.level, level.num)
                                                                } else {
                                                                    time_bin = head(cds_choose_pr_sub.level, length(cds_choose_pr_sub.level))
                                                                }
                                                                
                                                                root.num = 1
                                                                if (length(sub.idx) < root.num) {
                                                                    root.num = length(sub.idx)
                                                                }
                                                                
                                                                cds.plt.graph_node <- igraph::V(principal_graph(cds.plt)[["UMAP"]])$name
                                                                cds_choose_pr_sub.graph_node <- igraph::V(principal_graph(cds_choose_pr_sub)[["UMAP"]])$name
                                                                cds.plt.graph_node.idx <- sapply(cds_choose_pr_sub.graph_node, FUN = function(X) which(cds.plt.graph_node %in% X))
                                                                
                                                                root_pr_nodes = get_earliest_principal_node(cds_choose_pr_sub, cds.graph_node.idx = cds.plt.graph_node.idx, time = time, time_bin = time_bin, root.num = root.num)
                                                                
                                                                root_pr_nodes.all <- c(root_pr_nodes.all, root_pr_nodes)
                                                            }
                                                        }
                                                        root_pr_nodes.all <- unique(root_pr_nodes.all)
                                                        cds_choose_pr <- order_cells(cds_choose.plt, reduction_method = "UMAP", root_pr_nodes = root_pr_nodes.all)
                                                        
                                                        png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',SUB.prefolder,'/',project,'.',SUB.prename,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'graph','_',trj.prt.name,'.','order','.',pseudo.name,'.','partition','.png', sep = ''), width = 5.165, height = 4, units = 'in', res = 600)
                                                        p <- plot_cells(cds_choose_pr, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = "pseudotime", label_cell_groups = FALSE, group_label_size = 3, label_leaves = FALSE, label_branch_points = FALSE, label_roots = TRUE, show_trajectory_graph = TRUE, graph_label_size = 3)
                                                        print(p)
                                                        dev.off()
                                                        
                                                        png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',SUB.prefolder,'/',project,'.',SUB.prename,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'graph','_',trj.prt.name,'.','order','.',pseudo.name,'.','partition','.','disableRootLabel','.','with_trajectory','.png', sep = ''), width = 5.165, height = 4, units = 'in', res = 600)
                                                        p <- plot_cells(cds_choose_pr, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = "pseudotime", label_cell_groups = FALSE, group_label_size = 3, label_leaves = FALSE, label_branch_points = FALSE, label_roots = FALSE, show_trajectory_graph = TRUE, graph_label_size = 3)
                                                        print(p)
                                                        dev.off()
                                                        
                                                        png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',SUB.prefolder,'/',project,'.',SUB.prename,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,batch.name,'graph','_',trj.prt.name,'.','order','.',pseudo.name,'.','partition','.','disableRootLabel','.','without_trajectory','.png', sep = ''), width = 5.165, height = 4, units = 'in', res = 600)
                                                        p <- plot_cells(cds_choose_pr, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, color_cells_by = "pseudotime", label_cell_groups = FALSE, group_label_size = 3, label_leaves = FALSE, label_branch_points = FALSE, label_roots = FALSE, show_trajectory_graph = FALSE, graph_label_size = 3)
                                                        print(p)
                                                        dev.off()
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    
                    
                    if (project.name %in% c("GSM4150377", "GSM4150378")) {
                        if (project.name %in% "GSM4150377") {
                            colData(cds)$product_name_alt = 
                                plyr::revalue((pData(cds)$product_name), 
                                              c("BMS345541" = "BMS345541", 
                                                "Dex" = "Dex", 
                                                "Nutlin3A" = "Nutlin3A", 
                                                "SAHA" = "SAHA"))
                        } else if (project.name %in% "GSM4150378") {
                            colData(cds)$product_name_alt = 
                                plyr::revalue((pData(cds)$product_name), 
                                              c("PD98059" = "PD98059", 
                                                "SL-327" = "SL-327", 
                                                "Sorafenib Tosylate" = "Sorafenib Tosylate", 
                                                "Trametinib (GSK1120212)" = "Trametinib (GSK1120212)", 
                                                "Vehicle" = "Vehicle"))
                        }
                        
                        for (redDim in 1:length(redDim.method)) {
                            dim.folder = paste(eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))),'D', sep = '')
                            dim.finalX = eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) - (2-1)
                            
                            if (dim.finalX >= 1) {
                                
                                for (dim.comb in 1:dim.finalX) {
                                    dim.x <- dim.comb
                                    dim.y <- dim.x + 1
                                    
                                    if (dim.finalX > 1) {
                                        if (dim.x == 1) {
                                            comb.run = 2
                                        } else {
                                            comb.run = 1
                                        }
                                    } else {
                                        comb.run = 1
                                    }
                                    
                                    for (dim.comb.run in 1:comb.run) {
                                        if (dim.comb.run > 1) {
                                            dim.y <- dim.y + 1
                                        }
                                        
                                        if (eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) == 2) {
                                            dim.xy <- paste('D',dim.x,'-','D',dim.y, sep = '')
                                        } else {
                                            dim.xy <- paste('D',dim.x,'-','D',dim.y, sep = '')
                                        }
                                        #print(dim.xy)
                                        
                                        cell_size = 0.35
                                        cell_stroke = I(cell_size / 2)
                                        
                                        # /// Dose ///
                                        png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',project,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,'.','dose','.png', sep = ''), width = 4.075, height = 4, units = 'in', res = 600)
                                        p <- colData(cds) %>% as.data.frame() %>% 
                                            ggplot() + 
                                            geom_point(aes(x = umap1, y = umap2, color = dose_character), size = I(cell_size), stroke = I(cell_stroke), na.rm = TRUE, alpha = 1) + 
                                            monocle:::monocle_theme_opts() + 
                                            scale_color_manual(dose_character.unit, values = dose_character.color) + 
                                            theme(legend.position = "none", text = element_text(size = 11)) + 
                                            xlab("UMAP 1") + 
                                            ylab("UMAP 2")
                                        print(p)
                                        dev.off()
                                        
                                        png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',project,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,'.','dose','_','legend','.png', sep = ''), width = 4.075 + 1, height = 4, units = 'in', res = 600)
                                        p <- colData(cds) %>% as.data.frame() %>% 
                                            ggplot() + 
                                            geom_point(aes(x = umap1, y = umap2, color = dose_character), size = I(cell_size), stroke = I(cell_stroke), na.rm = TRUE, alpha = 1) + 
                                            monocle:::monocle_theme_opts() + 
                                            scale_color_manual(dose_character.unit, values = dose_character.color) + 
                                            theme(legend.position = "right", text = element_text(size = 11)) + 
                                            xlab("UMAP 1") + 
                                            ylab("UMAP 2")
                                        p <- p + guides(color = guide_legend(override.aes = list(size = 4)))
                                        print(p)
                                        dev.off()
                                        
                                        # /// split into one by one ///
                                        for (prod.idx in 0:length(sort(unique(colData(cds)[,"product_name"])))) {
                                            if (prod.idx == 0) {
                                                prod.name = NULL
                                                cds_dataset = colData(cds) %>% as.data.frame()
                                                
                                                width.size = 4.075 + 1
                                                height.size = 4 + 2
                                                
                                                if (exists("project.subname")) {
                                                    filename = paste(project,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,'.','dose','_','legend','.',project.subname, sep = '')
                                                } else {
                                                    filename = paste(project,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,'.','dose','_','legend', sep = '')
                                                }
                                            } else {
                                                prod.name = sort(unique(colData(cds)[,"product_name"]))[prod.idx]
                                                prod.name = as.character(prod.name)
                                                cds_dataset = colData(cds) %>% as.data.frame() %>% filter(product_name_alt == prod.name)
                                                
                                                width.size = 4.075 + 1
                                                height.size = 4
                                                
                                                if (exists("project.subname")) {
                                                    filename = paste(project,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,'.','dose','_','legend','.',project.subname,'_',prod.name, sep = '')
                                                } else {
                                                    filename = paste(project,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,'.','dose','_','legend','.',prod.name, sep = '')
                                                }
                                            }
                                            
                                            if (is.null(prod.name)) {
                                                dose.plt = 1
                                            } else {
                                                if (prod.name != "Vehicle") {
                                                    dose.plt = 1
                                                } else {
                                                    #dose.plt = 0
                                                    dose.plt = 1
                                                }
                                            }
                                            
                                            if (dose.plt == 1) {
                                                png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',filename,'.png', sep = ''), width = width.size, height = height.size, units = 'in', res = 600)
                                                p <- cds_dataset %>% 
                                                    ggplot() + 
                                                    geom_point(aes(x = umap1, y = umap2, color = dose_character), size = I(cell_size), stroke = I(cell_stroke), na.rm = TRUE, alpha = 1) + 
                                                    monocle:::monocle_theme_opts() + 
                                                    scale_color_manual(dose_character.unit, values = dose_character.color) + 
                                                    facet_wrap(~product_name, ncol = 2) + 
                                                    theme(legend.position = "right", text = element_text(size = 11)) + 
                                                    xlab("UMAP 1") + 
                                                    ylab("UMAP 2")
                                                p <- p + guides(color = guide_legend(override.aes = list(size = 4)))
                                                print(p)
                                                dev.off()
                                            }
                                        }
                                        
                                        # /// Density (Cell Type) ///
                                        if (exists("project.subname")) {
                                            filename = paste(project,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,'.','cell_type_density','.',project.subname, sep = '')
                                        } else {
                                            filename = paste(project,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,'.','cell_type_density', sep = '')
                                        }
                                        png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',filename,'.png', sep = ''), width = 4.075, height = 4, units = 'in', res = 600)
                                        p <- colData(cds) %>% as.data.frame() %>% 
                                            ggplot() + 
                                            geom_point(data = pData(cds) %>% as.data.frame() %>% dplyr::select(-cell_type), 
                                                       aes(x = umap1, y = umap2), 
                                                       color = "grey80", 
                                                       size = I(cell_size), stroke = I(cell_stroke), na.rm = TRUE, alpha = 1) + 
                                            geom_density_2d(aes(x = umap1, y = umap2, color = cell_type), 
                                                            size = I(cell_size), na.rm = TRUE, alpha = 1) + 
                                            monocle:::monocle_theme_opts() + 
                                            scale_color_manual("Cell Line", values = c("A549" = "#438FCD", 
                                                                                       "K562" = "#D6C939", 
                                                                                       "MCF7" = "#DB3C6A")) + 
                                            #facet_wrap(~cell_type) + 
                                            theme(legend.position = "none", text = element_text(size = 11), strip.text.x = element_text(size = 11)) + 
                                            xlab("UMAP 1") + 
                                            ylab("UMAP 2")
                                        print(p)
                                        dev.off()
                                        
                                        # Figure 4B
                                        # product name in MAPK pathway
                                        pseudotime_bin_summary = pData(cds) %>% as.data.frame() %>% group_by(Pseudotime_bin) %>% summarise(num_in_bin = n())
                                        pseudotime_bin_summary = left_join(as.data.frame(pData(cds)), pseudotime_bin_summary, by = "Pseudotime_bin")
                                        
                                        pseudotime_bin_summary = 
                                            pseudotime_bin_summary %>% 
                                            group_by(Pseudotime_bin, cell_type) %>% 
                                            add_tally() %>% 
                                            mutate(fraction_cell_type = n/num_in_bin) %>% 
                                            dplyr::select(-n) %>% 
                                            ungroup()
                                        
                                        pseudotime_bin_summary = 
                                            pseudotime_bin_summary %>% 
                                            group_by(Pseudotime_bin, product_name_alt) %>% 
                                            add_tally() %>% 
                                            mutate(fraction_product_name = n/num_in_bin) %>% 
                                            dplyr::select(-n) %>% 
                                            ungroup()
                                        
                                        pseudotime_bin_summary = 
                                            pseudotime_bin_summary %>% 
                                            group_by(Pseudotime_bin, dose_character) %>% 
                                            add_tally() %>% 
                                            mutate(fraction_dose = n/num_in_bin) %>% 
                                            dplyr::select(-n) %>% 
                                            ungroup()
                                        
                                        # /// Ridge ///
                                        unique.keys.df = 
                                            pData(cds) %>% 
                                            as.data.frame() %>% 
                                            dplyr::select(cell_type, product_name_alt) %>% 
                                            unique()
                                        
                                        vehicle.df = 
                                            pData(cds) %>% 
                                            as.data.frame() %>% 
                                            dplyr::select(Pseudotime, dose_character, product_name_alt, cell_type) %>% 
                                            filter(product_name_alt == "Vehicle") %>% 
                                            dplyr::select(-product_name_alt)
                                        
                                        non.vehicle.df = 
                                            pData(cds) %>% 
                                            as.data.frame() %>% 
                                            dplyr::select(Pseudotime, dose_character, product_name_alt, cell_type) %>% 
                                            filter(product_name_alt != "Vehicle")
                                        
                                        duplicated.vehicle.df = inner_join(unique.keys.df, vehicle.df, by = "cell_type")
                                        
                                        duplicated.vehicle.df = 
                                            duplicated.vehicle.df %>% 
                                            filter(product_name_alt != "Vehicle") %>% 
                                            dplyr::select(cell_type, product_name_alt, Pseudotime, dose_character)
                                        
                                        non.vehicle.df = 
                                            non.vehicle.df %>% 
                                            dplyr::select(cell_type, product_name_alt, Pseudotime, dose_character)
                                        
                                        combined.df = rbind(duplicated.vehicle.df, non.vehicle.df)
                                        
                                        if (exists("project.subname")) {
                                            filename = paste(project,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,'.','ridge','.',project.subname, sep = '')
                                        } else {
                                            filename = paste(project,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,'.','ridge', sep = '')
                                        }
                                        png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',filename,'.png', sep = ''), width = 8, height = 4, units = 'in', res = 600)
                                        p <- combined.df %>% 
                                            filter(is.finite(Pseudotime)) %>% 
                                            ggplot() + 
                                            geom_density_ridges(aes(x = Pseudotime, y = dose_character, fill = dose_character), size = 0.25) + 
                                            monocle:::monocle_theme_opts() + 
                                            scale_fill_manual(dose_character.unit, values = dose_character.color) + 
                                            facet_wrap(product_name_alt~cell_type, ncol = 6) + 
                                            theme(legend.position = "none", 
                                                  text = element_text(size = 11), 
                                                  axis.title.x = element_blank(), 
                                                  axis.text.x = element_blank(), 
                                                  axis.ticks.x = element_blank()) + 
                                            xlab("") + 
                                            ylab("")
                                        print(p)
                                        dev.off()
                                        
                                        for (prod.idx in 1:length(sort(unique(colData(cds)[,"product_name"])))) {
                                            prod.name = sort(unique(colData(cds)[,"product_name"]))[prod.idx]
                                            prod.name = as.character(prod.name)
                                            
                                            if (prod.name != "Vehicle") {
                                                if (exists("project.subname")) {
                                                    filename = paste(project,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,'.','ridge','.',project.subname,'_',prod.name, sep = '')
                                                } else {
                                                    filename = paste(project,'.','trajectories','.','plot_cells','.',redDim.method[redDim],'_',dim.xy,'.','ridge','.',prod.name, sep = '')
                                                }
                                                png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xy,'/',filename,'.png', sep = ''), width = 2.4, height = 3.2, units = 'in', res = 600)
                                                p <- combined.df %>% 
                                                    filter(is.finite(Pseudotime)) %>% 
                                                    filter(product_name_alt == prod.name) %>% 
                                                    ggplot() + 
                                                    geom_density_ridges(aes(x = Pseudotime, y = dose_character, fill = dose_character), size = 0.25) + 
                                                    monocle:::monocle_theme_opts() + 
                                                    scale_fill_manual(dose_character.unit, values = dose_character.color) + 
                                                    facet_wrap(~cell_type, ncol = 3) + 
                                                    theme(legend.position = "none", 
                                                          text = element_text(size = 11), 
                                                          axis.title.x = element_blank(), 
                                                          axis.text.x = element_blank(), 
                                                          axis.ticks.x = element_blank(), 
                                                          strip.text.x = element_blank()) + 
                                                    xlab("") + 
                                                    ylab("")
                                                print(p)
                                                dev.off()
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    
                    
                    # /// 3d ///
                    # Working with 3D trajectories
                    for (redDim in 1:length(redDim.method)) {
                        
                        if (eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) == 2) {
                            dim.folder = paste(eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))),'D', sep = '')
                            #dim.finalX = eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) - (3-1)
                            #dim.folder = paste(3,'D', sep = '')
                            dim.finalX = 3 - (3-1)
                            
                            mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim], sep = '')
                            subDir <- dim.folder
                            if (! file.exists(file.path(mainDir, subDir))) {
                                dir.create(file.path(mainDir, subDir))
                            }
                            
                            if (redDim == 1) {
                                cds_3d <- new_cell_data_set(expression_matrix,
                                                            cell_metadata = cell_metadata,
                                                            gene_metadata = gene_metadata)
                                
                                if (ncol(cds_3d) / 2 <= 50) {
                                    # Warning in (function (A, nv = 5, nu = nv, maxit = 1000, work = nv + 7, reorth = TRUE,  :
                                    # You're computing too large a percentage of total singular values, use a standard svd instead.
                                    if (floor(ncol(cds_3d) / 2) == ncol(cds_3d) / 2) {
                                        cds_3d <- preprocess_cds(cds_3d, method = "PCA", num_dim = floor(ncol(cds_3d) / 2) - 1)
                                    } else {
                                        cds_3d <- preprocess_cds(cds_3d, method = "PCA", num_dim = floor(ncol(cds_3d) / 2))
                                    }
                                } else {
                                    cds_3d <- preprocess_cds(cds_3d, method = "PCA", num_dim = 50) # method: Default is "PCA". method = c("PCA", "LSI") # "PCA": Principal Components Analysis (the standard for RNA-seq), "LSI": Latent Semantic Indexing (common in ATAC-seq)
                                }
                            }
                            
                            cds_3d <- reduce_dimension(cds_3d, reduction_method = redDim.method[redDim], max_components = 3)
                            cds_3d <- cluster_cells(cds_3d, reduction_method = redDim.method[redDim])
                            
                            cds_3d <- cds_3d
                            cds_3d.origin.processed = cds_3d
                            
                            
                            # /// subset ///
                            for (SUB in SUB.run) {
                                if (SUB == 0) {
                                    next # next halts the processing of the current iteration and advances the looping index.
                                } else if (SUB == 1) {
                                    SUB.prefolder = 'cell_subset'
                                    SUB.prename = 'cell_subset'
                                } else if (SUB == 2) {
                                    if (length(cell_subset) == 0 & cell_choose == TRUE) {
                                        SUB.prefolder = 'choose'
                                        SUB.prename = 'choose'
                                    } else if (length(cell_subset) > 0 & cell_choose == TRUE) {
                                        SUB.prefolder = paste('cell_subset','/','choose', sep = '')
                                        SUB.prename = paste('cell_subset','.','choose', sep = '')
                                    }
                                }
                                
                                if (SUB == 1) {
                                    cds_3d_choose.cell <- cell_subset
                                    
                                    cds_3d.colidx <- which(sapply(colnames(cds_3d), FUN = function(X) X %in% cds_3d_choose.cell))
                                    if (length(cds_3d.colidx) > 0) {
                                        cds_3d_choose <- cds_3d[, cds_3d.colidx]
                                    } else {
                                        cds_3d_choose <- cds_3d
                                    }
                                } else if (SUB == 2) {
                                    if (length(cell_subset) == 0 & cell_choose == TRUE) {
                                        cds_3d4choose <- cds_3d
                                    } else if (length(cell_subset) > 0 & cell_choose == TRUE) {
                                        cds_3d4choose <- cds_3d_choose
                                    }
                                    
                                    #cds_3d_choose.before <- choose_cells(cds_3d4choose, reduction_method = redDim.method[redDim], clear_cds = TRUE, return_list = FALSE)
                                    #cds_3d_choose.after <- choose_cells(cds_3d4choose, reduction_method = redDim.method[redDim], clear_cds = FALSE, return_list = FALSE)
                                    cds_3d_choose.cell <- choose_cells(cds_3d4choose, reduction_method = redDim.method[redDim], clear_cds = FALSE, return_list = TRUE)
                                    
                                    if (length(cell_subset) > 0) {
                                        cds_3d_choose.cell <- intersect(cell_subset, cds_3d_choose.cell)
                                    }
                                    
                                    cds_3d4choose.colidx <- which(sapply(colnames(cds_3d4choose), FUN = function(X) X %in% cds_3d_choose.cell))
                                    if (length(cds_3d4choose.colidx) > 0) {
                                        cds_3d_choose <- cds_3d4choose[, cds_3d4choose.colidx]
                                    } else {
                                        cds_3d_choose <- cds_3d4choose
                                    }
                                }
                                
                                mapping.idx = sapply(rownames(colData(cds_3d_choose)), FUN = function(X) which(names(cds_3d_choose@clusters[[redDim.method[redDim]]]$clusters) %in% X))
                                cds_3d_choose@clusters[[redDim.method[redDim]]]$clusters = cds_3d_choose@clusters[[redDim.method[redDim]]]$clusters[mapping.idx]
                                cds_3d_choose@clusters[[redDim.method[redDim]]]$clusters = factor(cds_3d_choose@clusters[[redDim.method[redDim]]]$clusters, levels = mixedsort(as.character(unique(cds_3d_choose@clusters[[redDim.method[redDim]]]$clusters))))
                                colData(cds_3d_choose)$cluster = cds_3d_choose@clusters[[redDim.method[redDim]]]$clusters
                                
                                mapping.idx = sapply(rownames(colData(cds_3d_choose)), FUN = function(X) which(names(cds_3d_choose@clusters[[redDim.method[redDim]]]$partitions) %in% X))
                                cds_3d_choose@clusters[[redDim.method[redDim]]]$partitions = cds_3d_choose@clusters[[redDim.method[redDim]]]$partitions[mapping.idx]
                                cds_3d_choose@clusters[[redDim.method[redDim]]]$partitions = factor(cds_3d_choose@clusters[[redDim.method[redDim]]]$partitions, levels = mixedsort(as.character(unique(cds_3d_choose@clusters[[redDim.method[redDim]]]$partitions))))
                                colData(cds_3d_choose)$partition = cds_3d_choose@clusters[[redDim.method[redDim]]]$partitions
                                
                                mapping.idx = sapply(rownames(colData(cds_3d_choose)), FUN = function(X) which(rownames(cds_3d_choose@principal_graph_aux[[redDim.method[redDim]]]$pr_graph_cell_proj_closest_vertex) %in% X))
                                cds_3d_choose@principal_graph_aux[[redDim.method[redDim]]]$pr_graph_cell_proj_closest_vertex = cds_3d_choose@principal_graph_aux[[redDim.method[redDim]]]$pr_graph_cell_proj_closest_vertex[mapping.idx, , drop = FALSE]
                                
                                saveRDS(cds_3d_choose, paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',project,'.',SUB.prename,'.','trajectories','.','cds_3d_choose','.rds', sep = ''))
                                #cds_3d_choose <- readRDS(paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',project,'.',SUB.prename,'.','trajectories','.','cds_3d_choose','.rds', sep = ''))
                                
                                if (SUB == 1) {
                                    cds_3d_choose.cell_subset <- cds_3d_choose
                                } else if (SUB == 2) {
                                    cds_3d_choose.cell_choose <- cds_3d_choose
                                }
                            }
                            
                            
                            if (graph.run == 1) {
                                
                                if (redDim.method[redDim] == "UMAP") {
                                    
                                    ## Step 5: Learn a graph **************************************************************************
                                    cds_3d.before_graph <- cds_3d
                                    
                                    if (project.name %in% c("GSM4150377", "GSM4150378")) {
                                        graph_parameters <- list()
                                        graph_parameters[["minimal_branch_len"]] = 10
                                        graph_parameters[["ncenter"]] = 750
                                        
                                        cds_3d <- learn_graph(cds_3d.before_graph, use_partition = TRUE, close_loop = FALSE, learn_graph_control = graph_parameters) #Please run reduce_dimension with reduction_method = UMAP and cluster_cells before running learn_graph.
                                        cds_3d.no_partition <- learn_graph(cds_3d.before_graph, use_partition = FALSE, close_loop = FALSE, learn_graph_control = graph_parameters) #Please run reduce_dimension with reduction_method = UMAP and cluster_cells before running learn_graph.
                                        
                                        cds_3d <- append_umap_coordinates(cds_3d)
                                    } else {
                                        cds_3d <- learn_graph(cds_3d.before_graph, use_partition = TRUE) #Please run reduce_dimension with reduction_method = UMAP and cluster_cells before running learn_graph.
                                        cds_3d.no_partition <- learn_graph(cds_3d.before_graph, use_partition = FALSE) #Please run reduce_dimension with reduction_method = UMAP and cluster_cells before running learn_graph.
                                    }
                                    
                                    # refer to select_cells.R
                                    sel <- c()
                                    sel$nodes <- igraph::V(principal_graph(cds_3d)[[redDim.method[redDim]]])$name
                                    sel$cells <- colnames(cds_3d)
                                    principal_graph(cds_3d)[[redDim.method[redDim]]] <- igraph::induced_subgraph(principal_graph(cds_3d)[[redDim.method[redDim]]], sel$nodes)
                                    
                                    sel <- c()
                                    sel$nodes <- igraph::V(principal_graph(cds_3d.no_partition)[[redDim.method[redDim]]])$name
                                    sel$cells <- colnames(cds_3d.no_partition)
                                    principal_graph(cds_3d.no_partition)[[redDim.method[redDim]]] <- igraph::induced_subgraph(principal_graph(cds_3d.no_partition)[[redDim.method[redDim]]], sel$nodes)
                                    
                                    cds_3d.origin.processed = cds_3d
                                    cds_3d.no_partition.origin.processed = cds_3d.no_partition
                                    
                                    
                                    # /// subset ///
                                    for (SUB in SUB.run) {
                                        if (SUB == 0) {
                                            next # next halts the processing of the current iteration and advances the looping index.
                                        } else if (SUB == 1) {
                                            SUB.prefolder = 'cell_subset'
                                            SUB.prename = 'cell_subset'
                                        } else if (SUB == 2) {
                                            if (length(cell_subset) == 0 & cell_choose == TRUE) {
                                                SUB.prefolder = 'choose'
                                                SUB.prename = 'choose'
                                            } else if (length(cell_subset) > 0 & cell_choose == TRUE) {
                                                SUB.prefolder = paste('cell_subset','/','choose', sep = '')
                                                SUB.prename = paste('cell_subset','.','choose', sep = '')
                                            }
                                        }
                                        
                                        ## Step 5: Learn a graph **************************************************************************
                                        if (SUB == 1) {
                                            cds_3d_choose.before_graph <- cds_3d_choose.cell_subset
                                        } else if (SUB == 2) {
                                            cds_3d_choose.before_graph <- cds_3d_choose.cell_choose
                                        }
                                        
                                        if (project.name %in% c("GSM4150377", "GSM4150378")) {
                                            graph_parameters <- list()
                                            graph_parameters[["minimal_branch_len"]] = 10
                                            graph_parameters[["ncenter"]] = 750
                                            
                                            cds_3d_choose <- learn_graph(cds_3d_choose.before_graph, use_partition = TRUE, close_loop = FALSE, learn_graph_control = graph_parameters) #Please run reduce_dimension with reduction_method = UMAP and cluster_cells before running learn_graph.
                                            cds_3d_choose.no_partition <- learn_graph(cds_3d_choose.before_graph, use_partition = FALSE, close_loop = FALSE, learn_graph_control = graph_parameters) #Please run reduce_dimension with reduction_method = UMAP and cluster_cells before running learn_graph.
                                            
                                            cds_3d_choose <- append_umap_coordinates(cds_3d_choose)
                                        } else {
                                            cds_3d_choose <- learn_graph(cds_3d_choose.before_graph, use_partition = TRUE) #Please run reduce_dimension with reduction_method = UMAP and cluster_cells before running learn_graph.
                                            cds_3d_choose.no_partition <- learn_graph(cds_3d_choose.before_graph, use_partition = FALSE) #Please run reduce_dimension with reduction_method = UMAP and cluster_cells before running learn_graph.
                                        }
                                        
                                        # refer to select_cells.R
                                        sel <- c()
                                        sel$nodes <- igraph::V(principal_graph(cds_3d_choose)[[redDim.method[redDim]]])$name
                                        sel$cells <- colnames(cds_3d_choose)
                                        principal_graph(cds_3d_choose)[[redDim.method[redDim]]] <- igraph::induced_subgraph(principal_graph(cds_3d_choose)[[redDim.method[redDim]]], sel$nodes)
                                        
                                        sel <- c()
                                        sel$nodes <- igraph::V(principal_graph(cds_3d_choose.no_partition)[[redDim.method[redDim]]])$name
                                        sel$cells <- colnames(cds_3d_choose.no_partition)
                                        principal_graph(cds_3d_choose.no_partition)[[redDim.method[redDim]]] <- igraph::induced_subgraph(principal_graph(cds_3d_choose.no_partition)[[redDim.method[redDim]]], sel$nodes)
                                        
                                        if (SUB == 1) {
                                            cds_3d_choose.cell_subset <- cds_3d_choose
                                        } else if (SUB == 2) {
                                            cds_3d_choose.cell_choose <- cds_3d_choose
                                        }
                                        
                                        if (SUB == 1) {
                                            cds_3d_choose.no_partition.cell_subset <- cds_3d_choose.no_partition
                                        } else if (SUB == 2) {
                                            cds_3d_choose.no_partition.cell_choose <- cds_3d_choose.no_partition
                                        }
                                        
                                        
                                        # Subset cells by branch
                                        # It is often useful to subset cells based on their branch in the trajectory.
                                        # The function choose_graph_segments allows you to do so interactively.
                                        # reduction_method: The reduction method to plot while choosing cells. Currently only "UMAP" is supported.
                                        
                                        # Choose cells along the path of a principal graph
                                        # clear_cds: Logical, clear CDS slots before returning. After clearing the cds, re-run processing from preprocess_cds(), ... Default is TRUE.
                                        # return_list: Logical, return a list of cells instead of a subsetted CDS object.
                                        
                                        # /// use_partition = TRUE ///
                                        if (SUB == 1) {
                                            cds_3d_choose.graph <- cds_3d_choose
                                        }
                                        
                                        if (SUB == 2) {
                                            if (length(cell_subset) == 0 & cell_choose == TRUE) {
                                                cds_3d4choose <- cds_3d
                                            } else if (length(cell_subset) > 0 & cell_choose == TRUE) {
                                                #cds_3d4choose <- cds_3d_choose
                                                cds_3d4choose <- cds_3d_choose.cell_subset
                                            }
                                            
                                            #cds_3d_choose.graph <- choose_graph_segments(cds_3d4choose, reduction_method = "UMAP", clear_cds = TRUE, return_list = FALSE)
                                            cds_3d_choose.graph <- choose_graph_segments(cds_3d4choose, reduction_method = "UMAP", clear_cds = FALSE, return_list = FALSE)
                                            
                                            while (ncol(cds_3d_choose.graph) > length(cds_3d_choose.cell)) {
                                                cat('[Warning] Please re-choose cells for trajectory graph. The cells you chose are more than the chosen cells before trajectory plot.','\n', sep = '')
                                                
                                                #cds_3d_choose.graph <- choose_graph_segments(cds_3d4choose, reduction_method = "UMAP", clear_cds = TRUE, return_list = FALSE)
                                                cds_3d_choose.graph <- choose_graph_segments(cds_3d4choose, reduction_method = "UMAP", clear_cds = FALSE, return_list = FALSE)
                                            }
                                            while (! all(colnames(cds_3d_choose.graph) %in% cds_3d_choose.cell)) {
                                                cat('[Warning] Please re-choose cells for trajectory graph. The cells you chose are not the same or the subset of the chosen cells before trajectory plot.','\n', sep = '')
                                                
                                                #cds_3d_choose.graph <- choose_graph_segments(cds_3d4choose, reduction_method = "UMAP", clear_cds = TRUE, return_list = FALSE)
                                                cds_3d_choose.graph <- choose_graph_segments(cds_3d4choose, reduction_method = "UMAP", clear_cds = FALSE, return_list = FALSE)
                                            }
                                        }
                                        
                                        mapping.idx = sapply(rownames(colData(cds_3d_choose.graph)), FUN = function(X) which(names(cds_3d_choose.graph@clusters[[redDim.method[redDim]]]$clusters) %in% X))
                                        cds_3d_choose.graph@clusters[[redDim.method[redDim]]]$clusters = cds_3d_choose.graph@clusters[[redDim.method[redDim]]]$clusters[mapping.idx]
                                        cds_3d_choose.graph@clusters[[redDim.method[redDim]]]$clusters = factor(cds_3d_choose.graph@clusters[[redDim.method[redDim]]]$clusters, levels = mixedsort(as.character(unique(cds_3d_choose.graph@clusters[[redDim.method[redDim]]]$clusters))))
                                        colData(cds_3d_choose.graph)$cluster = cds_3d_choose.graph@clusters[[redDim.method[redDim]]]$clusters
                                        
                                        mapping.idx = sapply(rownames(colData(cds_3d_choose.graph)), FUN = function(X) which(names(cds_3d_choose.graph@clusters[[redDim.method[redDim]]]$partitions) %in% X))
                                        cds_3d_choose.graph@clusters[[redDim.method[redDim]]]$partitions = cds_3d_choose.graph@clusters[[redDim.method[redDim]]]$partitions[mapping.idx]
                                        cds_3d_choose.graph@clusters[[redDim.method[redDim]]]$partitions = factor(cds_3d_choose.graph@clusters[[redDim.method[redDim]]]$partitions, levels = mixedsort(as.character(unique(cds_3d_choose.graph@clusters[[redDim.method[redDim]]]$partitions))))
                                        colData(cds_3d_choose.graph)$partition = cds_3d_choose.graph@clusters[[redDim.method[redDim]]]$partitions
                                        
                                        mapping.idx = sapply(rownames(colData(cds_3d_choose.graph)), FUN = function(X) which(rownames(cds_3d_choose.graph@principal_graph_aux[[redDim.method[redDim]]]$pr_graph_cell_proj_closest_vertex) %in% X))
                                        cds_3d_choose.graph@principal_graph_aux[[redDim.method[redDim]]]$pr_graph_cell_proj_closest_vertex = cds_3d_choose.graph@principal_graph_aux[[redDim.method[redDim]]]$pr_graph_cell_proj_closest_vertex[mapping.idx, , drop = FALSE]
                                        
                                        saveRDS(cds_3d_choose.graph, paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',project,'.',SUB.prename,'.','trajectories','.','cds_3d_choose.graph','.rds', sep = ''))
                                        #cds_3d_choose.graph <- readRDS(paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',project,'.',SUB.prename,'.','trajectories','.','cds_3d_choose.graph','.rds', sep = ''))
                                        
                                        if (SUB == 1) {
                                            cds_3d_choose.cell_subset.graph <- cds_3d_choose.graph
                                        } else if (SUB == 2) {
                                            cds_3d_choose.cell_choose.graph <- cds_3d_choose.graph
                                        }
                                        
                                        # /// use_partition = FALSE ///
                                        if (SUB == 1) {
                                            cds_3d_choose.no_partition.graph <- cds_3d_choose.no_partition
                                        }
                                        
                                        if (SUB == 2) {
                                            if (length(cell_subset) == 0 & cell_choose == TRUE) {
                                                cds_3d4choose <- cds_3d.no_partition
                                            } else if (length(cell_subset) > 0 & cell_choose == TRUE) {
                                                #cds_3d4choose <- cds_3d_choose.no_partition
                                                cds_3d4choose <- cds_3d_choose.no_partition.cell_subset
                                            }
                                            
                                            #cds_3d_choose.no_partition.graph <- choose_graph_segments(cds_3d4choose, reduction_method = "UMAP", clear_cds = TRUE, return_list = FALSE)
                                            cds_3d_choose.no_partition.graph <- choose_graph_segments(cds_3d4choose, reduction_method = "UMAP", clear_cds = FALSE, return_list = FALSE)
                                            
                                            while (ncol(cds_3d_choose.no_partition.graph) > length(cds_3d_choose.cell)) {
                                                cat('[Warning] Please re-choose cells for trajectory graph. The cells you chose are more than the chosen cells before trajectory plot.','\n', sep = '')
                                                
                                                #cds_3d_choose.no_partition.graph <- choose_graph_segments(cds_3d4choose, reduction_method = "UMAP", clear_cds = TRUE, return_list = FALSE)
                                                cds_3d_choose.no_partition.graph <- choose_graph_segments(cds_3d4choose, reduction_method = "UMAP", clear_cds = FALSE, return_list = FALSE)
                                            }
                                            while (! all(colnames(cds_3d_choose.no_partition.graph) %in% cds_3d_choose.cell)) {
                                                cat('[Warning] Please re-choose cells for trajectory graph. The cells you chose are not the same or the subset of the chosen cells before trajectory plot.','\n', sep = '')
                                                
                                                #cds_3d_choose.no_partition.graph <- choose_graph_segments(cds_3d4choose, reduction_method = "UMAP", clear_cds = TRUE, return_list = FALSE)
                                                cds_3d_choose.no_partition.graph <- choose_graph_segments(cds_3d4choose, reduction_method = "UMAP", clear_cds = FALSE, return_list = FALSE)
                                            }
                                        }
                                        
                                        mapping.idx = sapply(rownames(colData(cds_3d_choose.no_partition.graph)), FUN = function(X) which(names(cds_3d_choose.no_partition.graph@clusters[[redDim.method[redDim]]]$clusters) %in% X))
                                        cds_3d_choose.no_partition.graph@clusters[[redDim.method[redDim]]]$clusters = cds_3d_choose.no_partition.graph@clusters[[redDim.method[redDim]]]$clusters[mapping.idx]
                                        cds_3d_choose.no_partition.graph@clusters[[redDim.method[redDim]]]$clusters = factor(cds_3d_choose.no_partition.graph@clusters[[redDim.method[redDim]]]$clusters, levels = mixedsort(as.character(unique(cds_3d_choose.no_partition.graph@clusters[[redDim.method[redDim]]]$clusters))))
                                        colData(cds_3d_choose.no_partition.graph)$cluster = cds_3d_choose.no_partition.graph@clusters[[redDim.method[redDim]]]$clusters
                                        
                                        mapping.idx = sapply(rownames(colData(cds_3d_choose.no_partition.graph)), FUN = function(X) which(names(cds_3d_choose.no_partition.graph@clusters[[redDim.method[redDim]]]$partitions) %in% X))
                                        cds_3d_choose.no_partition.graph@clusters[[redDim.method[redDim]]]$partitions = cds_3d_choose.no_partition.graph@clusters[[redDim.method[redDim]]]$partitions[mapping.idx]
                                        cds_3d_choose.no_partition.graph@clusters[[redDim.method[redDim]]]$partitions = factor(cds_3d_choose.no_partition.graph@clusters[[redDim.method[redDim]]]$partitions, levels = mixedsort(as.character(unique(cds_3d_choose.no_partition.graph@clusters[[redDim.method[redDim]]]$partitions))))
                                        colData(cds_3d_choose.no_partition.graph)$partition = cds_3d_choose.no_partition.graph@clusters[[redDim.method[redDim]]]$partitions
                                        
                                        mapping.idx = sapply(rownames(colData(cds_3d_choose.no_partition.graph)), FUN = function(X) which(rownames(cds_3d_choose.no_partition.graph@principal_graph_aux[[redDim.method[redDim]]]$pr_graph_cell_proj_closest_vertex) %in% X))
                                        cds_3d_choose.no_partition.graph@principal_graph_aux[[redDim.method[redDim]]]$pr_graph_cell_proj_closest_vertex = cds_3d_choose.no_partition.graph@principal_graph_aux[[redDim.method[redDim]]]$pr_graph_cell_proj_closest_vertex[mapping.idx, , drop = FALSE]
                                        
                                        saveRDS(cds_3d_choose.no_partition.graph, paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',project,'.',SUB.prename,'.','trajectories','.','cds_3d_choose.no_partition.graph','.rds', sep = ''))
                                        #cds_3d_choose.no_partition.graph <- readRDS(paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',project,'.',SUB.prename,'.','trajectories','.','cds_3d_choose.no_partition.graph','.rds', sep = ''))
                                        
                                        if (SUB == 1) {
                                            cds_3d_choose.no_partition.cell_subset.graph <- cds_3d_choose.no_partition.graph
                                        } else if (SUB == 2) {
                                            cds_3d_choose.no_partition.cell_choose.graph <- cds_3d_choose.no_partition.graph
                                        }
                                    }
                                    
                                    
                                    ## Step 6: Order cells ****************************************************************************
                                    # /// all ///
                                    if (project.name == "tutorial") {
                                        time = "embryo.time.bin"
                                        level.num = 2
                                    } else {
                                        time = time.input
                                        level.num = 1
                                    }
                                    
                                    if (! is.null(time)) {
                                        
                                        # /// use_partition = TRUE ///
                                        cds_3d.level = names(which(sapply(levels(colData(cds_3d)[, time]), FUN = function(X) X %in% unique(colData(cds_3d)[, time]))))
                                        if (length(cds_3d.level) > level.num) {
                                            time_bin = head(cds_3d.level, level.num)
                                        } else {
                                            time_bin = head(cds_3d.level, length(cds_3d.level))
                                        }
                                        
                                        root.num = 2
                                        
                                        cds_3d <- order_cells(cds_3d, reduction_method = "UMAP", root_pr_nodes = get_earliest_principal_node(cds_3d, time = time, time_bin = time_bin, root.num = root.num))
                                        
                                        if (project.name %in% c("GSM4150377", "GSM4150378")) {
                                            pData(cds_3d)$Pseudotime = cds_3d@principal_graph_aux[["UMAP"]]$pseudotime
                                            pData(cds_3d)$Pseudotime_bin = cut(pData(cds_3d)$Pseudotime, breaks = unique(quantile(pData(cds_3d)$Pseudotime, seq(0, 1, 0.1))), labels = FALSE)
                                        }
                                        
                                        cds_3d.origin.processed = cds_3d
                                        
                                        # /// use_partition = FALSE ///
                                        cds_3d.no_partition.level = names(which(sapply(levels(colData(cds_3d.no_partition)[, time]), FUN = function(X) X %in% unique(colData(cds_3d.no_partition)[, time]))))
                                        if (length(cds_3d.no_partition.level) > level.num) {
                                            time_bin = head(cds_3d.no_partition.level, level.num)
                                        } else {
                                            time_bin = head(cds_3d.no_partition.level, length(cds_3d.no_partition.level))
                                        }
                                        
                                        root.num = 2
                                        
                                        cds_3d.no_partition <- order_cells(cds_3d.no_partition, reduction_method = "UMAP", root_pr_nodes = get_earliest_principal_node(cds_3d.no_partition, time = time, time_bin = time_bin, root.num = root.num))
                                        
                                        cds_3d.no_partition.origin.processed = cds_3d.no_partition
                                        
                                        
                                        # /// subset ///
                                        for (SUB in SUB.run) {
                                            if (SUB == 0) {
                                                next # next halts the processing of the current iteration and advances the looping index.
                                            } else if (SUB == 1) {
                                                SUB.prefolder = 'cell_subset'
                                                SUB.prename = 'cell_subset'
                                            } else if (SUB == 2) {
                                                if (length(cell_subset) == 0 & cell_choose == TRUE) {
                                                    SUB.prefolder = 'choose'
                                                    SUB.prename = 'choose'
                                                } else if (length(cell_subset) > 0 & cell_choose == TRUE) {
                                                    SUB.prefolder = paste('cell_subset','/','choose', sep = '')
                                                    SUB.prename = paste('cell_subset','.','choose', sep = '')
                                                }
                                            }
                                            
                                            ## Step 6: Order cells ****************************************************************************
                                            # /// use_partition = TRUE ///
                                            if (SUB == 1) {
                                                cds_3d_choose <- cds_3d_choose.cell_subset
                                            } else if (SUB == 2) {
                                                cds_3d_choose <- cds_3d_choose.cell_choose
                                            }
                                            
                                            # /// use_partition = FALSE ///
                                            if (SUB == 1) {
                                                cds_3d_choose.no_partition <- cds_3d_choose.no_partition.cell_subset
                                            } else if (SUB == 2) {
                                                cds_3d_choose.no_partition <- cds_3d_choose.no_partition.cell_choose
                                            }
                                            
                                            # /// all ///
                                            # /// use_partition = TRUE ///
                                            cds_3d_choose.level = names(which(sapply(levels(colData(cds_3d_choose)[, time]), FUN = function(X) X %in% unique(colData(cds_3d_choose)[, time]))))
                                            if (length(cds_3d_choose.level) > level.num) {
                                                time_bin = head(cds_3d_choose.level, level.num)
                                            } else {
                                                time_bin = head(cds_3d_choose.level, length(cds_3d_choose.level))
                                            }
                                            
                                            root.num = 2
                                            
                                            cds_3d.graph_node <- igraph::V(principal_graph(cds_3d)[["UMAP"]])$name
                                            cds_3d_choose.graph_node <- igraph::V(principal_graph(cds_3d_choose)[["UMAP"]])$name
                                            cds_3d.graph_node.idx <- sapply(cds_3d_choose.graph_node, FUN = function(X) which(cds_3d.graph_node %in% X))
                                            
                                            cds_3d_choose <- order_cells(cds_3d_choose, reduction_method = "UMAP", root_pr_nodes = get_earliest_principal_node(cds_3d_choose, cds.graph_node.idx = cds_3d.graph_node.idx, time = time, time_bin = time_bin, root.num = root.num))
                                            
                                            if (project.name %in% c("GSM4150377", "GSM4150378")) {
                                                pData(cds_3d_choose)$Pseudotime = cds_3d_choose@principal_graph_aux[["UMAP"]]$pseudotime
                                                pData(cds_3d_choose)$Pseudotime_bin = cut(pData(cds_3d_choose)$Pseudotime, breaks = unique(quantile(pData(cds_3d_choose)$Pseudotime, seq(0, 1, 0.1))), labels = FALSE)
                                            }
                                            
                                            if (SUB == 1) {
                                                cds_3d_choose.cell_subset <- cds_3d_choose
                                            } else if (SUB == 2) {
                                                cds_3d_choose.cell_choose <- cds_3d_choose
                                            }
                                            
                                            # /// use_partition = FALSE ///
                                            cds_3d_choose.no_partition.level = names(which(sapply(levels(colData(cds_3d_choose.no_partition)[, time]), FUN = function(X) X %in% unique(colData(cds_3d_choose.no_partition)[, time]))))
                                            if (length(cds_3d_choose.no_partition.level) > level.num) {
                                                time_bin = head(cds_3d_choose.no_partition.level, level.num)
                                            } else {
                                                time_bin = head(cds_3d_choose.no_partition.level, length(cds_3d_choose.no_partition.level))
                                            }
                                            
                                            root.num = 2
                                            
                                            cds_3d.no_partition.graph_node <- igraph::V(principal_graph(cds_3d.no_partition)[["UMAP"]])$name
                                            cds_3d_choose.no_partition.graph_node <- igraph::V(principal_graph(cds_3d_choose.no_partition)[["UMAP"]])$name
                                            cds_3d.no_partition.graph_node.idx <- sapply(cds_3d_choose.no_partition.graph_node, FUN = function(X) which(cds_3d.no_partition.graph_node %in% X))
                                            
                                            cds_3d_choose.no_partition <- order_cells(cds_3d_choose.no_partition, reduction_method = "UMAP", root_pr_nodes = get_earliest_principal_node(cds_3d_choose.no_partition, cds.graph_node.idx = cds_3d.no_partition.graph_node.idx, time = time, time_bin = time_bin, root.num = root.num))
                                            
                                            if (SUB == 1) {
                                                cds_3d_choose.no_partition.cell_subset <- cds_3d_choose.no_partition
                                            } else if (SUB == 2) {
                                                cds_3d_choose.no_partition.cell_choose <- cds_3d_choose.no_partition
                                            }
                                        }
                                    }
                                }
                            }
                        } else {
                            dim.folder = paste(eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))),'D', sep = '')
                            dim.finalX = eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) - (3-1)
                            
                            cds_3d <- cds
                            cds_3d.origin.processed = cds_3d
                            
                            cds_3d.no_partition <- cds.no_partition
                            cds_3d.no_partition.origin.processed = cds_3d.no_partition
                        }
                        
                        
                        # /// subset ///
                        for (SUB in SUB.run) {
                            if (SUB == 0) {
                                #next # next halts the processing of the current iteration and advances the looping index.
                            } else if (SUB == 1) {
                                SUB.prefolder = 'cell_subset'
                                SUB.prename = 'cell_subset'
                            } else if (SUB == 2) {
                                if (length(cell_subset) == 0 & cell_choose == TRUE) {
                                    SUB.prefolder = 'choose'
                                    SUB.prename = 'choose'
                                } else if (length(cell_subset) > 0 & cell_choose == TRUE) {
                                    SUB.prefolder = paste('cell_subset','/','choose', sep = '')
                                    SUB.prename = paste('cell_subset','.','choose', sep = '')
                                }
                            }
                            
                            if (eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) == 2) {
                                if (SUB == 1) {
                                    cds_choose_3d <- cds_3d_choose.cell_subset
                                    cds_choose_3d.no_partition <- cds_3d_choose.no_partition.cell_subset
                                } else if (SUB == 2) {
                                    cds_choose_3d <- cds_3d_choose.cell_choose
                                    cds_choose_3d.no_partition <- cds_3d_choose.no_partition.cell_choose
                                }
                            } else {
                                if (SUB == 1) {
                                    cds_choose_3d <- cds_choose.cell_subset
                                    cds_choose_3d.no_partition <- cds_choose.no_partition.cell_subset
                                } else if (SUB == 2) {
                                    cds_choose_3d <- cds_choose.cell_choose
                                    cds_choose_3d.no_partition <- cds_choose.no_partition.cell_choose
                                }
                            }
                            
                            if (SUB == 0) {
                                cds_3d = cds_3d.origin.processed
                                cds_3d.no_partition = cds_3d.no_partition.origin.processed
                                file.prename = 'trajectories'
                            } else if (SUB == 1) {
                                cds_3d = cds_choose_3d
                                cds_3d.no_partition = cds_choose_3d.no_partition
                                file.prename = paste(SUB.prename,'.','trajectories', sep = '')
                            } else if (SUB == 2) {
                                cds_3d = cds_choose_3d
                                cds_3d.no_partition = cds_choose_3d.no_partition
                                file.prename = paste(SUB.prename,'.','trajectories', sep = '')
                            }
                            
                            if (dim.finalX >= 1) {
                                
                                for (dim.comb in 1:dim.finalX) {
                                    dim.x <- dim.comb
                                    dim.y <- dim.x + 1
                                    dim.z <- dim.y + 1
                                    
                                    if (eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) == 2) {
                                        dim.xyz <- paste('D',dim.x,'-','D',dim.y,'-','D',dim.z, sep = '')
                                    } else {
                                        dim.xyz <- paste('D',dim.x,'-','D',dim.y,'-','D',dim.z, sep = '')
                                    }
                                    #print(dim.xyz)
                                    
                                    mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder, sep = '')
                                    subDir <- dim.xyz
                                    if (! file.exists(file.path(mainDir, subDir))) {
                                        dir.create(file.path(mainDir, subDir))
                                    }
                                    
                                    if (SUB == 0) {
                                        file.prefolder = dim.xyz
                                    } else {
                                        mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',dim.xyz, sep = '')
                                        subDir <- SUB.prefolder
                                        if (! file.exists(file.path(mainDir, subDir))) {
                                            dir.create(file.path(mainDir, subDir))
                                        }
                                        file.prefolder = paste(dim.xyz,'/',SUB.prefolder, sep = '')
                                    }
                                    
                                    plt.3d <- unique(c(cell_type.input, stage.input))
                                    
                                    if (redDim.method[redDim] == "UMAP") {
                                        if (! is.null(time)) {
                                            plt.3d.num = length(plt.3d) + 1
                                        } else { # No pseudotime calculated for reduction_method = UMAP. Please first run order_cells with reduction_method = UMAP.
                                            plt.3d.num = length(plt.3d)
                                        }
                                    } else if (redDim.method[redDim] == "tSNE") {
                                        plt.3d.num = length(plt.3d)
                                    }
                                    
                                    if (plt.3d.num > 0) {
                                        for (color.plt in 1:plt.3d.num) {
                                            if (color.plt <= length(plt.3d)) {
                                                if (project.name %in% c("Pancreas", "DrugPairing", "CRISPRINT")) {
                                                    color_cells_by.attr = plt.3d[color.plt]
                                                } else if (project.name %in% c("LeoPharma", "Eczema")) {
                                                    color_cells_by.attr = plt.3d[color.plt]
                                                } else if (project.name %in% c("Paper")) {
                                                    color_cells_by.attr = plt.3d[color.plt]
                                                } else if (project.name %in% c("GSE139944", c("GSM4150376", "GSM4150377", "GSM4150378", "GSM4150379"))) {
                                                    color_cells_by.attr = plt.3d[color.plt]
                                                } else {
                                                    color_cells_by.attr = "partition"
                                                }
                                            } else if (color.plt > length(plt.3d)) {
                                                color_cells_by.attr = "pseudotime"
                                            }
                                            
                                            if (redDim.method[redDim] == "UMAP") {
                                                trj.plt = 1
                                            } else if (redDim.method[redDim] == "tSNE") {
                                                trj.plt = 0
                                            }
                                            if (graph.run == 0) {
                                                trj.plt = 0
                                            }
                                            
                                            for (trj in 0:trj.plt) {
                                                if (trj == 0) {
                                                    trj.prt.endidx = 1
                                                } else if (trj == 1) {
                                                    trj.prt.endidx = 0
                                                }
                                                
                                                #for (trj.prt in 1:1) {
                                                for (trj.prt in 1:trj.prt.endidx) {
                                                    if (trj.prt == 1) {
                                                        trj.prt.name = 'with_partition'
                                                        cds_3d.plt = cds_3d
                                                    } else if (trj.prt == 0) {
                                                        trj.prt.name = 'without_partition'
                                                        cds_3d.plt = cds_3d.no_partition
                                                    }
                                                    
                                                    if (! is.null(cds_3d.plt)) {
                                                        
                                                        for (grid.plt in 1:0) {
                                                            if (grid.plt == 1) {
                                                                grid.name = 'with_grid'
                                                            } else if (grid.plt == 0) {
                                                                grid.name = 'without_grid'
                                                            }
                                                            
                                                            if (trj == 0) {
                                                                show_trajectory = FALSE
                                                                file.name = paste(project,'.',file.prename,'.','plot_cells_3d','.',redDim.method[redDim],'_',dim.xyz,'.',gsub("\\.","_",color_cells_by.attr),'.','without_trajectory', sep = '')
                                                                file.name = paste(file.name,'.',grid.name, sep = '')
                                                            } else if (trj == 1) {
                                                                show_trajectory = TRUE
                                                                file.name = paste(project,'.',file.prename,'.','plot_cells_3d','.',redDim.method[redDim],'_',dim.xyz,'.',gsub("\\.","_",color_cells_by.attr),'.','with_trajectory', sep = '')
                                                                file.name = paste(file.name,'.',trj.prt.name, sep = '')
                                                                file.name = paste(file.name,'.',grid.name, sep = '')
                                                            }
                                                            
                                                            file.html = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/',file.name,'.html', sep = '')
                                                            file.png = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/',file.name,'.png', sep = '')
                                                            
                                                            if (! is.null(color_palette)) {
                                                                cell_type.colidx <- which(paste(unique(c(color_cells_by.attr, unlist(strsplit(color_cells_by.attr, '\\.|_'))[1])),'.','color', sep = '') %in% names(color_palette))
                                                                if (length(cell_type.colidx) == 1) {
                                                                    cds_3d_plot_obj <- plot_cells_3d(cds_3d.plt, reduction_method = redDim.method[redDim], dims = c(dim.x, dim.y, dim.z), cell_size = cell_size_3d, show_trajectory_graph = show_trajectory, color_cells_by = color_cells_by.attr, color_palette = color_palette[[which(names(color_palette) %in% paste(unique(c(color_cells_by.attr, unlist(strsplit(color_cells_by.attr, '\\.|_'))[1])),'.','color', sep = ''))]])
                                                                } else {
                                                                    cds_3d_plot_obj <- plot_cells_3d(cds_3d.plt, reduction_method = redDim.method[redDim], dims = c(dim.x, dim.y, dim.z), cell_size = cell_size_3d, show_trajectory_graph = show_trajectory, color_cells_by = color_cells_by.attr)
                                                                }
                                                            } else {
                                                                cds_3d_plot_obj <- plot_cells_3d(cds_3d.plt, reduction_method = redDim.method[redDim], dims = c(dim.x, dim.y, dim.z), cell_size = cell_size_3d, show_trajectory_graph = show_trajectory, color_cells_by = color_cells_by.attr)
                                                            }
                                                            
                                                            if (grid.plt == 0) {
                                                                cds_3d_plot_obj <- cds_3d_plot_obj %>% plotly::layout(scene = scene)
                                                            }
                                                            
                                                            # https://community.rstudio.com/t/save-viewer-object-rendered-in-rstudio-as-image/32796
                                                            saveWidget(cds_3d_plot_obj, file.html)
                                                            #webshot(file.html, file.png)
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                            
                            cds_3d = cds_3d.origin.processed
                            cds_3d.no_partition = cds_3d.no_partition.origin.processed
                        }
                    }
                    
                    
                    
                    # === Differential expression analysis ==============================================================================================================
                    # https://cole-trapnell-lab.github.io/monocle3/docs/differential/
                    
                    # There are two approaches for differential analysis in Monocle:
                    # Regression analysis: using fit_models(), you can evaluate whether each gene depends on variables such as time, treatments, etc.
                    # Graph-autocorrelation analysis: using graph_test(), you can find genes that vary over a trajectory or between clusters.
                    
                    # Monocle also comes with specialized functions for finding co-regulated modules of differentially expressed genes. 
                    # Monocle also allows you to interactively interrogate specific clusters or regions of a trajectory (e.g. branch points) for genes that vary within them.
                    
                    # (Optional) Perform differential expression analysis: https://cole-trapnell-lab.github.io/monocle3/docs/differential/
                    
                    # *** query gene(s) ***
                    if (! is.null(genes_of_interest.input)) {
                        # /// query.gene_group_df ///
                        genes_of_interest.new.all <- list()
                        for (gX in 1:length(genes_of_interest.input)) {
                            marker.genes <- genes_of_interest.input[[gX]]
                            marker.genes.idx <- which(sapply(marker.genes, FUN = function(X) X %in% rowData(cds)[,"gene_short_name"]))
                            if (length(marker.genes.idx) > 0) {
                                marker.genes = marker.genes[marker.genes.idx]
                                
                                genes_of_interest.new <- list(marker.genes)
                                names(genes_of_interest.new) = names(genes_of_interest.input[gX])
                                
                                genes_of_interest.new.all <- c(genes_of_interest.new.all, genes_of_interest.new)
                            }
                        }
                        genes_of_interest.new = genes_of_interest.new.all
                        
                        query.gene_group_df.all <- c()
                        for (gX in 1:length(genes_of_interest.new)) {
                            marker.genes <- genes_of_interest.new[[gX]]
                            marker.genes.idx <- which(sapply(marker.genes, FUN = function(X) X %in% rowData(cds)[,"gene_short_name"]))
                            if (length(marker.genes.idx) > 0) {
                                marker.genes = marker.genes[marker.genes.idx]
                                
                                query.gene_group_df <- data.frame("id" = marker.genes, 
                                                                  "module" = names(genes_of_interest.new[gX]), 
                                                                  "supermodule" = "genes_of_interest", 
                                                                  stringsAsFactors = FALSE)
                                for (dim_ in 1:eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = '')))) {
                                    query.gene_group_df <- data.frame(query.gene_group_df, 
                                                                      as.numeric(NA), 
                                                                      stringsAsFactors = FALSE)
                                    colnames(query.gene_group_df)[ncol(query.gene_group_df)] = paste('dim_',dim_, sep = '')
                                }
                                query.gene_group_df <- data.frame(query.gene_group_df, 
                                                                  "gene_short_name" = marker.genes, 
                                                                  stringsAsFactors = FALSE)
                                
                                query.gene_group_df[,"module"] <- factor(query.gene_group_df[,"module"], levels = unique(query.gene_group_df[,"module"]))
                                query.gene_group_df[,"supermodule"] <- factor(query.gene_group_df[,"supermodule"], levels = unique(query.gene_group_df[,"supermodule"]))
                                
                                query.gene_group_df.all <- rbind(query.gene_group_df.all, query.gene_group_df)
                            }
                        }
                        #query.gene_group_df.all <- tbl_df(query.gene_group_df.all)
                        query.gene_group_df.all <- tibble::as_tibble(query.gene_group_df.all)
                    }
                    
                    # *** signature gene(s) ***
                    if (! is.null(signature_genes.input)) {
                        # /// signature.gene_group_df ///
                        signature_genes.new.all <- list()
                        for (gX in 1:length(signature_genes.input)) {
                            marker.genes <- signature_genes.input[[gX]]
                            marker.genes.idx <- which(sapply(marker.genes, FUN = function(X) X %in% rowData(cds)[,"gene_short_name"]))
                            if (length(marker.genes.idx) > 0) {
                                marker.genes = marker.genes[marker.genes.idx]
                                
                                signature_genes.new <- list(marker.genes)
                                names(signature_genes.new) = names(signature_genes.input[gX])
                                
                                signature_genes.new.all <- c(signature_genes.new.all, signature_genes.new)
                            }
                        }
                        signature_genes.new = signature_genes.new.all
                        
                        signature.gene_group_df.all <- c()
                        for (gX in 1:length(signature_genes.new)) {
                            marker.genes <- signature_genes.new[[gX]]
                            marker.genes.idx <- which(sapply(marker.genes, FUN = function(X) X %in% rowData(cds)[,"gene_short_name"]))
                            if (length(marker.genes.idx) > 0) {
                                marker.genes = marker.genes[marker.genes.idx]
                                
                                signature.gene_group_df <- data.frame("id" = marker.genes, 
                                                                      "module" = names(signature_genes.new[gX]), 
                                                                      "supermodule" = "signature_genes", 
                                                                      stringsAsFactors = FALSE)
                                for (dim_ in 1:eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = '')))) {
                                    signature.gene_group_df <- data.frame(signature.gene_group_df, 
                                                                          as.numeric(NA), 
                                                                          stringsAsFactors = FALSE)
                                    colnames(signature.gene_group_df)[ncol(signature.gene_group_df)] = paste('dim_',dim_, sep = '')
                                }
                                signature.gene_group_df <- data.frame(signature.gene_group_df, 
                                                                      "gene_short_name" = marker.genes, 
                                                                      stringsAsFactors = FALSE)
                                
                                signature.gene_group_df[,"module"] <- factor(signature.gene_group_df[,"module"], levels = unique(signature.gene_group_df[,"module"]))
                                signature.gene_group_df[,"supermodule"] <- factor(signature.gene_group_df[,"supermodule"], levels = unique(signature.gene_group_df[,"supermodule"]))
                                
                                signature.gene_group_df.all <- rbind(signature.gene_group_df.all, signature.gene_group_df)
                            }
                        }
                        #signature.gene_group_df.all <- tbl_df(signature.gene_group_df.all)
                        signature.gene_group_df.all <- tibble::as_tibble(signature.gene_group_df.all)
                    }
                    
                    
                    for (redDim in 1:length(redDim.method)) {
                        
                        # /// subset ///
                        for (SUB in SUB.run) {
                            if (SUB == 0) {
                                #next # next halts the processing of the current iteration and advances the looping index.
                            } else if (SUB == 1) {
                                SUB.prefolder = 'cell_subset'
                                SUB.prename = 'cell_subset'
                            } else if (SUB == 2) {
                                if (length(cell_subset) == 0 & cell_choose == TRUE) {
                                    SUB.prefolder = 'choose'
                                    SUB.prename = 'choose'
                                } else if (length(cell_subset) > 0 & cell_choose == TRUE) {
                                    SUB.prefolder = paste('cell_subset','/','choose', sep = '')
                                    SUB.prename = paste('cell_subset','.','choose', sep = '')
                                }
                            }
                            
                            if (eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) == 2) {
                                if (SUB == 1) {
                                    cds_choose_3d <- cds_3d_choose.cell_subset
                                    cds_choose_3d.no_partition <- cds_3d_choose.no_partition.cell_subset
                                } else if (SUB == 2) {
                                    cds_choose_3d <- cds_3d_choose.cell_choose
                                    cds_choose_3d.no_partition <- cds_3d_choose.no_partition.cell_choose
                                }
                            } else {
                                if (SUB == 1) {
                                    cds_choose_3d <- cds_choose.cell_subset
                                    cds_choose_3d.no_partition <- cds_choose.no_partition.cell_subset
                                } else if (SUB == 2) {
                                    cds_choose_3d <- cds_choose.cell_choose
                                    cds_choose_3d.no_partition <- cds_choose.no_partition.cell_choose
                                }
                            }
                            
                            if (SUB == 0) {
                                cds = cds.origin.processed
                                cds_3d = cds_3d.origin.processed
                                cds_3d.no_partition = cds_3d.no_partition.origin.processed
                                file.prename = 'trajectories'
                            } else if (SUB == 1) {
                                cds = cds_choose.cell_subset
                                cds_3d = cds_choose_3d
                                cds_3d.no_partition = cds_choose_3d.no_partition
                                file.prename = paste(SUB.prename,'.','trajectories', sep = '')
                            } else if (SUB == 2) {
                                cds = cds_choose.cell_choose
                                cds_3d = cds_choose_3d
                                cds_3d.no_partition = cds_choose_3d.no_partition
                                file.prename = paste(SUB.prename,'.','trajectories', sep = '')
                            }
                            
                            dim.folder = paste(eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))),'D', sep = '')
                            
                            mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim], sep = '')
                            subDir <- dim.folder
                            if (! file.exists(file.path(mainDir, subDir))) {
                                dir.create(file.path(mainDir, subDir))
                            }
                            
                            mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder, sep = '')
                            if (eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) == 2) {
                                #differential.folder <- paste('D12','/','differential', sep = '')
                                differential.folder <- 'differential'
                            } else {
                                differential.folder <- 'differential'
                            }
                            subDir <- differential.folder
                            if (! file.exists(file.path(mainDir, subDir))) {
                                dir.create(file.path(mainDir, subDir))
                            }
                            
                            if (SUB == 0) {
                                differential.prefolder = differential.folder
                            } else {
                                mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',differential.folder, sep = '')
                                subDir <- SUB.prefolder
                                if (! file.exists(file.path(mainDir, subDir))) {
                                    dir.create(file.path(mainDir, subDir), recursive = TRUE)
                                }
                                differential.prefolder = paste(differential.folder,'/',SUB.prefolder, sep = '')
                            }
                            
                            mapping.idx = sapply(rownames(colData(cds)), FUN = function(X) which(names(cds@clusters[[redDim.method[redDim]]]$clusters) %in% X))
                            cds@clusters[[redDim.method[redDim]]]$clusters = cds@clusters[[redDim.method[redDim]]]$clusters[mapping.idx]
                            cds@clusters[[redDim.method[redDim]]]$clusters = factor(cds@clusters[[redDim.method[redDim]]]$clusters, levels = mixedsort(as.character(unique(cds@clusters[[redDim.method[redDim]]]$clusters))))
                            colData(cds)$cluster = cds@clusters[[redDim.method[redDim]]]$clusters
                            
                            mapping.idx = sapply(rownames(colData(cds)), FUN = function(X) which(names(cds@clusters[[redDim.method[redDim]]]$partitions) %in% X))
                            cds@clusters[[redDim.method[redDim]]]$partitions = cds@clusters[[redDim.method[redDim]]]$partitions[mapping.idx]
                            cds@clusters[[redDim.method[redDim]]]$partitions = factor(cds@clusters[[redDim.method[redDim]]]$partitions, levels = mixedsort(as.character(unique(cds@clusters[[redDim.method[redDim]]]$partitions))))
                            colData(cds)$partition = cds@clusters[[redDim.method[redDim]]]$partitions
                            
                            mapping.idx = sapply(rownames(colData(cds)), FUN = function(X) which(rownames(cds@principal_graph_aux[[redDim.method[redDim]]]$pr_graph_cell_proj_closest_vertex) %in% X))
                            cds@principal_graph_aux[[redDim.method[redDim]]]$pr_graph_cell_proj_closest_vertex = cds@principal_graph_aux[[redDim.method[redDim]]]$pr_graph_cell_proj_closest_vertex[mapping.idx, , drop = FALSE]
                            
                            
                            # --- Regression analysis -------------------------------------------------------------------------
                            # With regression:
                            # Monocle can analyze several thousands of genes even in large experiments, making it useful for discovering dynamically regulated genes during the biological process you're studying.
                            
                            # Monocle works by fitting a regression model to each gene. 
                            # You can specify this model to account for various factors in your experiment (time, treatment, and so on).
                            
                            # /// Choosing a distribution for modeling gene expression ///
                            # Monocle uses generalized linear models to capture how a gene's expression depends on each variable in the experiment. 
                            # These models require you to specify a distribution that describes gene expression values. 
                            # Most studies that use this approach to analyze their gene expression data use the negative binomial distribution, which is often appropriate for sequencing read or UMI count data. 
                            # The negative binomial is at the core of many packages for RNA-seq analysis, such as DESeq2.
                            
                            # Monocle's fit_models() supports the negative binomial distribution and several others listed in the table below. 
                            # The default is the "quasipoisson", which is very similar to the negative binomial. 
                            # Quasipoisson is a a bit less accurate than the negative binomial but much faster to fit, making it well suited to datasets with thousands of cells.
                            
                            # There are several allowed values for expression_family:
                            # fit_models(..., expression_family, ...) # expression_family: Specifies the family function used for expression responses. Can be one of "quasipoisson", "negbinomial", "poisson", "binomial", "gaussian", "zipoisson", or "zinegbinomial". Default is "quasipoisson".
                            # expression_family = "quasipoisson" # Quasi-poisson Distribution: https://en.wikipedia.org/wiki/Poisson_regression   # Default for fit_models(). Recommended for most users.
                            # expression_family = "negbinomial" # Negative binomial Distribution: https://en.wikipedia.org/wiki/Negative_binomial_distribution
                            # expression_family = "poisson" # Poisson Distribution: https://en.wikipedia.org/wiki/Poisson_regression
                            # expression_family = "binomial" # Binomial Distribution: https://en.wikipedia.org/wiki/Logistic_regression
                            
                            # Likelihood based analysis and quasipoisson
                            # The quasi-poisson distribution doesn't have a real likelihood function, so some of Monocle's methods won't work with it. 
                            # Several of the columns in results tables from evaluate_fits() and compare_models() will be NA.
                            
                            
                            # Let's begin with a small set of genes that we know are important in ciliated neurons to demonstrate Monocle's capabilities:
                            if (project.name == "tutorial") {
                                # /// marker.genes ///
                                gene_set.list <- list()
                                
                                ciliated_genes <- list(unique(c("che-1", "hlh-17", "nhr-6", "dmd-6", "ceh-36", "ham-1")))
                                if (length(unlist(ciliated_genes)) == 1) {
                                    names(ciliated_genes) = "ciliated_gene"
                                } else if (length(unlist(ciliated_genes)) > 1) {
                                    names(ciliated_genes) = "ciliated_genes"
                                }
                                gene_set.list <- c(gene_set.list, ciliated_genes)
                                
                                if (length(gene_set.list) > 0) {
                                    for (N in 1:length(gene_set.list)) {
                                        #print(names(gene_set.list[N]))
                                        
                                        marker.genes <- gene_set.list[[N]]
                                        marker.genes.idx <- which(sapply(marker.genes, FUN = function(X) X %in% rowData(cds)[,"gene_short_name"]))
                                        if (length(marker.genes.idx) > 0) {
                                            marker.genes = marker.genes[marker.genes.idx]
                                            
                                            rowidx <- sapply(marker.genes, FUN = function(X) which(rowData(cds)[,"gene_short_name"] %in% X))
                                            rowidx <- unlist(rowidx)
                                            
                                            cds_subset <- cds[rowidx, , drop = FALSE]
                                            
                                            #rowData(cds_subset)[,"gene_short_name"] = factor(rowData(cds_subset)[,"gene_short_name"], levels = unique(rowData(cds_subset)[,"gene_short_name"]))
                                            
                                            # if you use cluster_cells, you can test for genes that differ between clusters and partitions by using ~cluster or ~partition (respectively) as your model formula. 
                                            if (project.name == "tutorial") {
                                                xt.all <- c("embryo.time", "cluster", "partition", batch)
                                            } else {
                                                xt.all <- c("cluster", "partition", batch)
                                            }
                                            
                                            for (x in 1:length(xt.all)) {
                                                
                                                xt = xt.all[x]
                                                
                                                # a generalized linear model: log(yi) = β0 + βt xt
                                                # where 
                                                # yi is a random variable corresponding to the expression values of gene i
                                                # xt is the time each cell was collected (in minutes)
                                                # βt capture the effect of time on expression
                                                # β0 is an intercept term
                                                # testing whether its βt is significantly different from zero
                                                gene_fits <- fit_models(cds_subset, model_formula_str = paste('~',xt, sep = ''))
                                                
                                                # coefficient_table() tests whether each coefficient differs significantly from zero under the Wald test. 
                                                # By default, coefficient_table() adjusts these p-values for multiple hypothesis testing using the method of Benjamini and Hochberg. 
                                                # These adjusted values can be found in the q_value column. 
                                                fit_coefs <- coefficient_table(gene_fits)
                                                #print(fit_coefs[,"term"])
                                                #print(unique(fit_coefs[,"term"]))
                                                
                                                #emb_time_terms <- fit_coefs %>% filter(term == xt)
                                                #print(emb_time_terms)
                                                emb_time_terms <- fit_coefs %>% filter(term != "(Intercept)") # We generally don't care about the intercept term β0
                                                #print(emb_time_terms)
                                                
                                                emb_time_terms <- emb_time_terms %>% mutate(q_value = p.adjust(p_value))
                                                #print(emb_time_terms)
                                                
                                                # We can filter the results and control the false discovery rate as follows:
                                                emb_time_terms.subset <- emb_time_terms %>% filter (q_value < 0.05) %>% select(gene_short_name, term, q_value, estimate)
                                                #print(emb_time_terms.subset)
                                                
                                                sig_genes <- emb_time_terms %>% filter (q_value < 0.05) %>% pull(gene_short_name)
                                                sig_genes <- unique(sig_genes)
                                                #print(sig_genes)
                                                
                                                # /// "violin" plot ///
                                                if (xt == "embryo.time") {
                                                    if (project.name == "tutorial") {
                                                        time = "embryo.time.bin"
                                                    } else {
                                                        time = time.input
                                                    }
                                                    run = c(xt, time)
                                                } else if (xt == "cluster") {
                                                    run = xt
                                                } else if (xt == "partition") {
                                                    run = xt
                                                } else if (xt == batch) {
                                                    run = xt
                                                }
                                                
                                                for (r in 1:length(run)) {
                                                    
                                                    group_var = run[r]
                                                    
                                                    if (group_var == "embryo.time") {
                                                        cds_subset.level = sort(unique(colData(cds_subset)[, group_var]), decreasing = FALSE)
                                                    } else if (group_var == "embryo.time.bin") {
                                                        #cds_subset.level = levels(colData(cds_subset)[, group_var])
                                                        cds_subset.level = names(which(sapply(levels(colData(cds_subset)[, group_var]), FUN = function(X) X %in% unique(colData(cds_subset)[, group_var]))))
                                                    } else if (group_var == "cluster") {
                                                        cds_subset.level = levels(cds_subset@clusters[[redDim.method[redDim]]]$clusters)
                                                    } else if (group_var == "partition") {
                                                        cds_subset.level = levels(cds_subset@clusters[[redDim.method[redDim]]]$partitions)
                                                    } else if (group_var == batch) {
                                                        #cds_subset.level = levels(colData(cds_subset)[, group_var])
                                                        cds_subset.level = names(which(sapply(levels(colData(cds_subset)[, group_var]), FUN = function(X) X %in% unique(colData(cds_subset)[, group_var]))))
                                                    }
                                                    
                                                    if (length(cds_subset.level) > 20) {
                                                        ncol.num = 1
                                                        width.size = length(cds_subset.level) / 4 + 0.75
                                                        height.size = nrow(cds_subset) * 1.5 + 0.75
                                                    } else {
                                                        ncol.num = 2
                                                        width.size = length(cds_subset.level) / 1.5 + 0.75
                                                        height.size = nrow(cds_subset) * 0.75 + 0.75
                                                    }
                                                    
                                                    if (max(nchar(cds_subset.level)) > 20) {
                                                        angle.num = 60
                                                    } else {
                                                        angle.num = 45
                                                    }
                                                    
                                                    png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',differential.prefolder,'/',project,'.',gsub("\\.","_",subdata.folder),'.',gsub("\\.","_",project.subname),'.',redDim.method[redDim],'.',file.prename,'.','differential','.','plot_genes_violin','.',gsub("\\.","_",names(gene_set.list[N])),'.',gsub("\\.","_",group_var),'.png', sep = ''), width = width.size, height = height.size, units = 'in', res = 600)
                                                    p <- plot_genes_violin(cds_subset, group_cells_by = group_var, ncol = ncol.num) + theme(axis.text.x = element_text(angle=angle.num, hjust=1))
                                                    suppressWarnings(print(p))
                                                    dev.off()
                                                }
                                            }
                                            
                                            
                                            if (project.name == "tutorial") {
                                                time = "embryo.time"
                                                batch = "batch"
                                            } else {
                                                time = time.input
                                                batch = batch.input
                                            }
                                            
                                            if (! is.null(batch)) {
                                                batch.name = paste('.',gsub("\\.","_",batch),'_','align','.', sep = '')
                                            } else {
                                                batch.name = '.'
                                            }
                                            
                                            
                                            # --- Controlling for batch effects and other factors ---------------------------------------------
                                            # You can also include multiple variables, for example ~embryo.time + batch, which can be very helpful for subtracting unwanted effects.
                                            gene_fits <- fit_models(cds_subset, model_formula_str = paste('~',time,' + ',batch, sep = ''))
                                            fit_coefs <- coefficient_table(gene_fits)
                                            
                                            emb_time_terms <- fit_coefs %>% filter(term != "(Intercept)")
                                            #print(emb_time_terms)
                                            
                                            emb_time_terms.subset <- emb_time_terms %>% select(gene_short_name, term, q_value, estimate)
                                            #print(emb_time_terms.subset)
                                            
                                            
                                            # --- Evaluating models of gene expression --------------------------------------------------------
                                            # How good are these models at "explaining" gene expression? 
                                            # We can evaluate the fits of each model using the evaluate_fits() function:
                                            fit_evaluate <- evaluate_fits(gene_fits)
                                            
                                            # Should we include the batch term in our model of gene expression or not? 
                                            # Monocle provides a function compare_models() that can help you decide. 
                                            # Compare models takes two models and returns the result of a likelihood ratio test between them. 
                                            # Any time you add terms to a model, it will improve the fit. 
                                            # But we should always to use the simplest model we can to explain our data. 
                                            # The likelihood ratio test helps us decide whether the improvement in fit is large enough to justify the complexity our extra terms introduce.
                                            
                                            # full model (model_tbl_full): This model is essentially a way of predicting the expression value of each gene in a given cell knowing both what time it was collected and which batch of cells it came from.
                                            time_batch_models <- fit_models(cds_subset, model_formula_str = paste('~',time,' + ',batch, sep = ''), expression_family = "negbinomial")
                                            
                                            # reduced model (model_tbl_reduced): it only knows about the time each cell was collected.
                                            time_models <- fit_models(cds_subset, model_formula_str = paste('~',time, sep = ''), expression_family = "negbinomial")
                                            
                                            # Because the full model has more information about each cell, it will do a better job of predicting the expression of the gene in each cell.
                                            # The question Monocle must answer for each gene is how much better the full model's prediction is than the reduced model's.
                                            # The greater the improvement that comes from knowing the batch of each cell, the more significant the result of the likelihood ratio test.
                                            # As we can see, all of the genes' likelihood ratio tests are significant, indicating that there are substantial batch effects in the data. We are therefore justified in adding the batch term to our model.
                                            model_comp <- compare_models(time_batch_models, time_models) %>% select(gene_short_name, q_value)
                                        }
                                    }
                                }
                            }
                            
                            
                            # --- Finding genes that change as a function of pseudotime ---------------------------------------
                            # Identifying the genes that change as cells progress along a trajectory is a core objective of this type of analysis. 
                            # Knowing the order in which genes go on and off can inform new models of development.
                            
                            # Let's return to the embryo data:
                            # How do we find the genes that are differentially expressed on the different paths through the trajectory? 
                            # How do we find the ones that are restricted to the beginning of the trajectory? Or excluded from it?
                            
                            if (! project.name %in% c("GSM4150377", "GSM4150378")) {
                                
                                # Once again, we turn to graph_test(), this time passing it neighbor_graph="principal_graph", which tells it to test whether cells at similar positions on the trajectory have correlated expression:
                                
                                # The function graph_test() uses a statistic from spatial autocorrelation analysis called Moran's I, which Cao & Spielmann et al showed to be effective in finding genes that vary in single-cell RNA-seq datasets.
                                # Moran's I: https://en.wikipedia.org/wiki/Moran%27s_I
                                # The data frame pr_graph_test_res has the Moran's I test results for each gene in the cell_data_set. 
                                # Significant values much less than zero are generally rare.
                                
                                if (redDim.method[redDim] == "UMAP") {
                                    # /// neighbor_graph = "principal_graph" ///
                                    # reduction_method: character, the method used to reduce dimension. Currently only supported for "UMAP".
                                    pr_graph_test_res.principal_graph <- graph_test(cds, reduction_method = "UMAP", neighbor_graph = "principal_graph", cores = 8) # neighbor_graph: String indicating what neighbor graph to use. "principal_graph" and "knn" are supported. "principal_graph" is recommended for trajectory analysis
                                    
                                    # If you'd like to rank the genes by effect size, sort this table by the morans_I column, which ranges from -1 to +1. 
                                    morans_I.orderidx = order(pr_graph_test_res.principal_graph[,"morans_I"], decreasing = TRUE)
                                    pr_graph_test_res.principal_graph <- pr_graph_test_res.principal_graph[morans_I.orderidx,]
                                    
                                    # A value of 0 indicates no effect, while +1 indicates perfect positive autocorrelation and suggests that nearby cells have very similar values of a gene's expression. 
                                    # Positive values indicate a gene is expressed in a focal region of the UMAP space (e.g. specific to one or more clusters).
                                    pr_graph_test_res.principal_graph.pos <- pr_graph_test_res.principal_graph[which(pr_graph_test_res.principal_graph[,"morans_I"] > 0),]
                                    pr_graph_test_res.principal_graph.zero <- pr_graph_test_res.principal_graph[which(pr_graph_test_res.principal_graph[,"morans_I"] == 0),]
                                    pr_graph_test_res.principal_graph.neg <- pr_graph_test_res.principal_graph[which(pr_graph_test_res.principal_graph[,"morans_I"] < 0),]
                                    pr_graph_test_res.principal_graph.NA <- pr_graph_test_res.principal_graph[which(is.na(pr_graph_test_res.principal_graph[,"morans_I"])),]
                                    
                                    pr_deg.principal_graph <- subset(pr_graph_test_res.principal_graph, q_value < 0.05)
                                    pr_deg.principal_graph.pos <- subset(pr_graph_test_res.principal_graph.pos, q_value < 0.05)
                                    pr_deg.principal_graph.zero <- subset(pr_graph_test_res.principal_graph.zero, q_value < 0.05)
                                    pr_deg.principal_graph.neg <- subset(pr_graph_test_res.principal_graph.neg, q_value < 0.05)
                                    pr_deg.principal_graph.NA <- subset(pr_graph_test_res.principal_graph.NA, q_value < 0.05)
                                    
                                    pr_deg_ids.principal_graph <- rownames(pr_deg.principal_graph)
                                    pr_deg_ids.principal_graph.pos <- rownames(pr_deg.principal_graph.pos)
                                    pr_deg_ids.principal_graph.zero <- rownames(pr_deg.principal_graph.zero)
                                    pr_deg_ids.principal_graph.neg <- rownames(pr_deg.principal_graph.neg)
                                    pr_deg_ids.principal_graph.NA <- rownames(pr_deg.principal_graph.NA)
                                    
                                    pr_deg_names.principal_graph <- as.character(pr_deg.principal_graph[,"gene_short_name"])
                                    pr_deg_names.principal_graph.pos <- as.character(pr_deg.principal_graph.pos[,"gene_short_name"])
                                    pr_deg_names.principal_graph.zero <- as.character(pr_deg.principal_graph.zero[,"gene_short_name"])
                                    pr_deg_names.principal_graph.neg <- as.character(pr_deg.principal_graph.neg[,"gene_short_name"])
                                    pr_deg_names.principal_graph.NA <- as.character(pr_deg.principal_graph.NA[,"gene_short_name"])
                                }
                            }
                            
                            
                            # Here are a couple of interesting genes that score as highly significant according to graph_test():
                            # /// marker.genes ///
                            gene_set.list <- list()
                            
                            if (project.name == "tutorial") {
                                example_genes <- list(unique(c("hlh-4", "gcy-8", "dac-1", "oig-8")))
                                if (length(unlist(example_genes)) == 1) {
                                    names(example_genes) = "example_gene"
                                } else if (length(unlist(example_genes)) > 1) {
                                    names(example_genes) = "example_genes"
                                }
                                gene_set.list <- c(gene_set.list, example_genes)
                            }
                            
                            if (exists("pr_deg_names.principal_graph")) {
                                # *** Select serial list ***
                                serial_num <- c(1:top_gene_num)
                                serial_num.idx <- which(serial_num <= length(pr_deg_names.principal_graph))
                                if (length(serial_num.idx) > 0) {
                                    serial_num <- serial_num[serial_num.idx]
                                } else {
                                    serial_num <- c()
                                }
                                if (length(serial_num) > 0) {
                                    for (serial in 1:length(serial_num)) {
                                        serial_genes <- list(pr_deg_names.principal_graph[serial_num[serial]])
                                        names(serial_genes) = paste('gene',serial_num[serial], sep = '')
                                        
                                        gene_set.list <- c(gene_set.list, serial_genes)
                                    }
                                }
                                
                                # *** Select top list ***
                                top_num <- c(1, 3, 5, 10, 20, 30, 50)
                                top_num.idx <- which(top_num <= top_gene_num)
                                if (length(top_num.idx) > 0) {
                                    top_num <- top_num[top_num.idx]
                                } else {
                                    top_num <- c(1)
                                }
                                top_num.idx <- which(top_num <= length(pr_deg_names.principal_graph))
                                if (length(top_num.idx) > 0) {
                                    top_num <- top_num[top_num.idx]
                                } else {
                                    top_num <- c()
                                }
                                if (length(top_num) > 0) {
                                    for (top in 1:length(top_num)) {
                                        top_genes <- list(unique(head(pr_deg_names.principal_graph, top_num[top])))
                                        if (top_num[top] == 1) {
                                            names(top_genes) = paste('top',top_num[top],'_','gene', sep = '')
                                        } else if (top_num[top] > 1) {
                                            names(top_genes) = paste('top',top_num[top],'_','genes', sep = '')
                                        }
                                        
                                        gene_set.list <- c(gene_set.list, top_genes)
                                    }
                                }
                            }
                            
                            # *** query gene(s) ***
                            if (! is.null(genes_of_interest.input)) {
                                for (gX in 1:length(genes_of_interest.input)) {
                                    for (gY in 1:length(genes_of_interest.input[[gX]])) {
                                        genes_of_interest.gX.gY <- list(genes_of_interest.input[[gX]][gY])
                                        names(genes_of_interest.gX.gY) = paste(names(genes_of_interest.input[gX]),'.',genes_of_interest.input[[gX]][gY], sep = '')
                                        gene_set.list <- c(gene_set.list, genes_of_interest.gX.gY)
                                    }
                                }
                                genes_of_interest <- genes_of_interest.input
                                gene_set.list <- c(gene_set.list, genes_of_interest)
                            }
                            
                            # *** signature gene(s) ***
                            if (! is.null(signature_genes.input)) {
                                for (gX in 1:length(signature_genes.input)) {
                                    for (gY in 1:length(signature_genes.input[[gX]])) {
                                        signature_genes.gX.gY <- list(signature_genes.input[[gX]][gY])
                                        names(signature_genes.gX.gY) = paste(names(signature_genes.input[gX]),'.',signature_genes.input[[gX]][gY], sep = '')
                                        gene_set.list <- c(gene_set.list, signature_genes.gX.gY)
                                    }
                                }
                                signature_genes <- signature_genes.input
                                gene_set.list <- c(gene_set.list, signature_genes)
                            }
                            
                            if (length(gene_set.list) > 0) {
                                for (N in 1:length(gene_set.list)) {
                                    #print(names(gene_set.list[N]))
                                    
                                    marker.genes <- gene_set.list[[N]]
                                    marker.genes.idx <- which(sapply(marker.genes, FUN = function(X) X %in% rowData(cds)[,"gene_short_name"]))
                                    if (length(marker.genes.idx) > 0) {
                                        marker.genes = marker.genes[marker.genes.idx]
                                        
                                        panel.num = length(unlist(unique(marker.genes)))
                                        if (panel.num == 1) {
                                            width.size = 7
                                            height.size = 6
                                        } else if (panel.num == 2) {
                                            width.size = 7 + (5 * 1)
                                            height.size = 6
                                        } else if (panel.num == 3) {
                                            width.size = 7 + (5 * 2)
                                            height.size = 6
                                        } else if (panel.num > 3) {
                                            width.size = 7 + (5 * 3)
                                            height.size = width.size - 1
                                        }
                                        
                                        #for (redDim in 1:length(redDim.method)) {
                                        
                                        # /// 2d ///
                                        dim.folder = paste(eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))),'D', sep = '')
                                        dim.finalX = eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) - (2-1)
                                        
                                        if (dim.finalX >= 1) {
                                            
                                            for (dim.comb in 1:dim.finalX) {
                                                dim.x <- dim.comb
                                                dim.y <- dim.x + 1
                                                
                                                if (dim.finalX > 1) {
                                                    if (dim.x == 1) {
                                                        comb.run = 2
                                                    } else {
                                                        comb.run = 1
                                                    }
                                                } else {
                                                    comb.run = 1
                                                }
                                                
                                                for (dim.comb.run in 1:comb.run) {
                                                    if (dim.comb.run > 1) {
                                                        dim.y <- dim.y + 1
                                                    }
                                                    
                                                    if (eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) == 2) {
                                                        dim.xy <- paste('D',dim.x,'-','D',dim.y, sep = '')
                                                    } else {
                                                        dim.xy <- paste('D',dim.x,'-','D',dim.y, sep = '')
                                                    }
                                                    #print(dim.xy)
                                                    
                                                    if (SUB == 0) {
                                                        file.prefolder = dim.xy
                                                    } else {
                                                        file.prefolder = paste(dim.xy,'/',SUB.prefolder, sep = '')
                                                    }
                                                    
                                                    mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder, sep = '')
                                                    subDir <- 'differential'
                                                    if (! file.exists(file.path(mainDir, subDir))) {
                                                        dir.create(file.path(mainDir, subDir))
                                                    }
                                                    mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential', sep = '')
                                                    subDir <- 'gene'
                                                    if (! file.exists(file.path(mainDir, subDir))) {
                                                        dir.create(file.path(mainDir, subDir))
                                                    }
                                                    
                                                    if (length(grep(paste(genes_of_interest.titles, collapse = '|'), names(gene_set.list[N]), value = TRUE, invert = FALSE)) > 0) {
                                                        mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential', sep = '')
                                                        subDir <- 'genes_of_interest'
                                                        if (! file.exists(file.path(mainDir, subDir))) {
                                                            dir.create(file.path(mainDir, subDir))
                                                        }
                                                        mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential','/','genes_of_interest', sep = '')
                                                        subDir <- unlist(strsplit(names(gene_set.list[N]), '\\.'))[1]
                                                        if (! file.exists(file.path(mainDir, subDir))) {
                                                            dir.create(file.path(mainDir, subDir))
                                                        }
                                                        mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential','/','genes_of_interest','/',unlist(strsplit(names(gene_set.list[N]), '\\.'))[1], sep = '')
                                                        subDir <- 'gene'
                                                        if (! file.exists(file.path(mainDir, subDir))) {
                                                            dir.create(file.path(mainDir, subDir))
                                                        }
                                                        
                                                        file.folder = paste('genes_of_interest','/',unlist(strsplit(names(gene_set.list[N]), '\\.'))[1],'/','gene', sep = '')
                                                        if (names(gene_set.list[N]) %in% genes_of_interest.titles) {
                                                            file.name.suffix = names(gene_set.list[N])
                                                        } else {
                                                            file.name.suffix = gsub("\\.","_",unlist(strsplit(names(gene_set.list[N]), '\\.'))[2])
                                                        }
                                                    } else if (length(grep(paste(signature_genes.titles, collapse = '|'), names(gene_set.list[N]), value = TRUE, invert = FALSE)) > 0) {
                                                        mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential', sep = '')
                                                        subDir <- 'signature_genes'
                                                        if (! file.exists(file.path(mainDir, subDir))) {
                                                            dir.create(file.path(mainDir, subDir))
                                                        }
                                                        mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential','/','signature_genes', sep = '')
                                                        subDir <- unlist(strsplit(names(gene_set.list[N]), '\\.'))[1]
                                                        if (! file.exists(file.path(mainDir, subDir))) {
                                                            dir.create(file.path(mainDir, subDir))
                                                        }
                                                        mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential','/','signature_genes','/',unlist(strsplit(names(gene_set.list[N]), '\\.'))[1], sep = '')
                                                        subDir <- 'gene'
                                                        if (! file.exists(file.path(mainDir, subDir))) {
                                                            dir.create(file.path(mainDir, subDir))
                                                        }
                                                        
                                                        file.folder = paste('signature_genes','/',unlist(strsplit(names(gene_set.list[N]), '\\.'))[1],'/','gene', sep = '')
                                                        if (names(gene_set.list[N]) %in% signature_genes.titles) {
                                                            file.name.suffix = names(gene_set.list[N])
                                                        } else {
                                                            file.name.suffix = gsub("\\.","_",unlist(strsplit(names(gene_set.list[N]), '\\.'))[2])
                                                        }
                                                    } else {
                                                        file.folder = 'gene'
                                                        file.name.suffix = gsub("\\.","_",names(gene_set.list[N]))
                                                    }
                                                    
                                                    if (is.null(file.name.suffix)) {
                                                        file.name = paste(project,'.',file.prename,'.','plot_cells','.',redDim.method[redDim],'_',dim.xy,'.','differential', sep = '')
                                                    } else {
                                                        if (is.na(file.name.suffix) | file.name.suffix == "") {
                                                            file.name = paste(project,'.',file.prename,'.','plot_cells','.',redDim.method[redDim],'_',dim.xy,'.','differential', sep = '')
                                                        } else {
                                                            file.name = paste(project,'.',file.prename,'.','plot_cells','.',redDim.method[redDim],'_',dim.xy,'.','differential','.',file.name.suffix, sep = '')
                                                        }
                                                    }
                                                    
                                                    png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential','/',file.folder,'/',file.name,'.png', sep = ''), width = width.size, height = height.size, units = 'in', res = 600)
                                                    p <- plot_cells(cds, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, genes = marker.genes, norm_method = "log", label_cell_groups = FALSE, show_trajectory_graph = FALSE, label_leaves = FALSE)
                                                    print(p)
                                                    dev.off()
                                                }
                                            }
                                        }
                                        
                                        # /// 3d ///
                                        if (length(grep('top', names(gene_set.list[N]), value = TRUE, invert = TRUE)) > 0) {
                                            
                                            if (length(marker.genes) == 1) {
                                                
                                                if (eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) == 2) {
                                                    dim.folder = paste(eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))),'D', sep = '')
                                                    #dim.finalX = eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) - (3-1)
                                                    #dim.folder = paste(3,'D', sep = '')
                                                    dim.finalX = 3 - (3-1)
                                                } else {
                                                    dim.folder = paste(eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))),'D', sep = '')
                                                    dim.finalX = eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) - (3-1)
                                                }
                                                
                                                if (dim.finalX >= 1) {
                                                    
                                                    for (dim.comb in 1:dim.finalX) {
                                                        dim.x <- dim.comb
                                                        dim.y <- dim.x + 1
                                                        dim.z <- dim.y + 1
                                                        
                                                        if (eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) == 2) {
                                                            dim.xyz <- paste('D',dim.x,'-','D',dim.y,'-','D',dim.z, sep = '')
                                                        } else {
                                                            dim.xyz <- paste('D',dim.x,'-','D',dim.y,'-','D',dim.z, sep = '')
                                                        }
                                                        #print(dim.xyz)
                                                        
                                                        if (SUB == 0) {
                                                            file.prefolder = dim.xyz
                                                        } else {
                                                            file.prefolder = paste(dim.xyz,'/',SUB.prefolder, sep = '')
                                                        }
                                                        
                                                        mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder, sep = '')
                                                        subDir <- 'differential'
                                                        if (! file.exists(file.path(mainDir, subDir))) {
                                                            dir.create(file.path(mainDir, subDir))
                                                        }
                                                        mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential', sep = '')
                                                        subDir <- 'gene'
                                                        if (! file.exists(file.path(mainDir, subDir))) {
                                                            dir.create(file.path(mainDir, subDir))
                                                        }
                                                        
                                                        if (length(grep(paste(genes_of_interest.titles, collapse = '|'), names(gene_set.list[N]), value = TRUE, invert = FALSE)) > 0) {
                                                            mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential', sep = '')
                                                            subDir <- 'genes_of_interest'
                                                            if (! file.exists(file.path(mainDir, subDir))) {
                                                                dir.create(file.path(mainDir, subDir))
                                                            }
                                                            mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential','/','genes_of_interest', sep = '')
                                                            subDir <- unlist(strsplit(names(gene_set.list[N]), '\\.'))[1]
                                                            if (! file.exists(file.path(mainDir, subDir))) {
                                                                dir.create(file.path(mainDir, subDir))
                                                            }
                                                            mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential','/','genes_of_interest','/',unlist(strsplit(names(gene_set.list[N]), '\\.'))[1], sep = '')
                                                            subDir <- 'gene'
                                                            if (! file.exists(file.path(mainDir, subDir))) {
                                                                dir.create(file.path(mainDir, subDir))
                                                            }
                                                            
                                                            file.folder = paste('genes_of_interest','/',unlist(strsplit(names(gene_set.list[N]), '\\.'))[1],'/','gene', sep = '')
                                                            if (names(gene_set.list[N]) %in% genes_of_interest.titles) {
                                                                file.name.suffix = names(gene_set.list[N])
                                                            } else {
                                                                file.name.suffix = gsub("\\.","_",unlist(strsplit(names(gene_set.list[N]), '\\.'))[2])
                                                            }
                                                        } else if (length(grep(paste(signature_genes.titles, collapse = '|'), names(gene_set.list[N]), value = TRUE, invert = FALSE)) > 0) {
                                                            mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential', sep = '')
                                                            subDir <- 'signature_genes'
                                                            if (! file.exists(file.path(mainDir, subDir))) {
                                                                dir.create(file.path(mainDir, subDir))
                                                            }
                                                            mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential','/','signature_genes', sep = '')
                                                            subDir <- unlist(strsplit(names(gene_set.list[N]), '\\.'))[1]
                                                            if (! file.exists(file.path(mainDir, subDir))) {
                                                                dir.create(file.path(mainDir, subDir))
                                                            }
                                                            mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential','/','signature_genes','/',unlist(strsplit(names(gene_set.list[N]), '\\.'))[1], sep = '')
                                                            subDir <- 'gene'
                                                            if (! file.exists(file.path(mainDir, subDir))) {
                                                                dir.create(file.path(mainDir, subDir))
                                                            }
                                                            
                                                            file.folder = paste('signature_genes','/',unlist(strsplit(names(gene_set.list[N]), '\\.'))[1],'/','gene', sep = '')
                                                            if (names(gene_set.list[N]) %in% signature_genes.titles) {
                                                                file.name.suffix = names(gene_set.list[N])
                                                            } else {
                                                                file.name.suffix = gsub("\\.","_",unlist(strsplit(names(gene_set.list[N]), '\\.'))[2])
                                                            }
                                                        } else {
                                                            file.folder = 'gene'
                                                            file.name.suffix = gsub("\\.","_",names(gene_set.list[N]))
                                                        }
                                                        
                                                        if (redDim.method[redDim] == "UMAP") {
                                                            trj.plt = 1
                                                        } else if (redDim.method[redDim] == "tSNE") {
                                                            trj.plt = 0
                                                        }
                                                        if (graph.run == 0) {
                                                            trj.plt = 0
                                                        }
                                                        
                                                        for (trj in 0:trj.plt) {
                                                            if (trj == 0) {
                                                                trj.prt.endidx = 1
                                                            } else if (trj == 1) {
                                                                trj.prt.endidx = 0
                                                            }
                                                            
                                                            #for (trj.prt in 1:1) {
                                                            for (trj.prt in 1:trj.prt.endidx) {
                                                                if (trj.prt == 1) {
                                                                    trj.prt.name = 'with_partition'
                                                                    cds_3d.plt = cds_3d
                                                                } else if (trj.prt == 0) {
                                                                    trj.prt.name = 'without_partition'
                                                                    cds_3d.plt = cds_3d.no_partition
                                                                }
                                                                
                                                                if (! is.null(cds_3d.plt)) {
                                                                    
                                                                    for (grid.plt in 1:0) {
                                                                        if (grid.plt == 1) {
                                                                            grid.name = 'with_grid'
                                                                        } else if (grid.plt == 0) {
                                                                            grid.name = 'without_grid'
                                                                        }
                                                                        
                                                                        if (trj == 0) {
                                                                            show_trajectory = FALSE
                                                                            if (length(marker.genes) == 1) {
                                                                                if (is.null(file.name.suffix)) {
                                                                                    file.name = paste(project,'.',file.prename,'.','plot_cells_3d','.',redDim.method[redDim],'_',dim.xyz,'.','differential','.',marker.genes,'.','without_trajectory', sep = '')
                                                                                } else {
                                                                                    if (is.na(file.name.suffix) | file.name.suffix == "") {
                                                                                        file.name = paste(project,'.',file.prename,'.','plot_cells_3d','.',redDim.method[redDim],'_',dim.xyz,'.','differential','.',marker.genes,'.','without_trajectory', sep = '')
                                                                                    } else {
                                                                                        file.name = paste(project,'.',file.prename,'.','plot_cells_3d','.',redDim.method[redDim],'_',dim.xyz,'.','differential','.',paste(unique(c(file.name.suffix, marker.genes)), collapse = '_'),'.','without_trajectory', sep = '')
                                                                                    }
                                                                                }
                                                                            } else {
                                                                                if (is.null(file.name.suffix)) {
                                                                                    file.name = paste(project,'.',file.prename,'.','plot_cells_3d','.',redDim.method[redDim],'_',dim.xyz,'.','differential','.','without_trajectory', sep = '')
                                                                                } else {
                                                                                    if (is.na(file.name.suffix) | file.name.suffix == "") {
                                                                                        file.name = paste(project,'.',file.prename,'.','plot_cells_3d','.',redDim.method[redDim],'_',dim.xyz,'.','differential','.','without_trajectory', sep = '')
                                                                                    } else {
                                                                                        file.name = paste(project,'.',file.prename,'.','plot_cells_3d','.',redDim.method[redDim],'_',dim.xyz,'.','differential','.',file.name.suffix,'.','without_trajectory', sep = '')
                                                                                    }
                                                                                }
                                                                            }
                                                                            file.name = paste(file.name,'.',grid.name, sep = '')
                                                                        } else if (trj == 1) {
                                                                            show_trajectory = TRUE
                                                                            if (length(marker.genes) == 1) {
                                                                                if (is.null(file.name.suffix)) {
                                                                                    file.name = paste(project,'.',file.prename,'.','plot_cells_3d','.',redDim.method[redDim],'_',dim.xyz,'.','differential','.',marker.genes,'.','with_trajectory', sep = '')
                                                                                } else {
                                                                                    if (is.na(file.name.suffix) | file.name.suffix == "") {
                                                                                        file.name = paste(project,'.',file.prename,'.','plot_cells_3d','.',redDim.method[redDim],'_',dim.xyz,'.','differential','.',marker.genes,'.','with_trajectory', sep = '')
                                                                                    } else {
                                                                                        file.name = paste(project,'.',file.prename,'.','plot_cells_3d','.',redDim.method[redDim],'_',dim.xyz,'.','differential','.',paste(unique(c(file.name.suffix, marker.genes)), collapse = '_'),'.','with_trajectory', sep = '')
                                                                                    }
                                                                                }
                                                                            } else {
                                                                                if (is.null(file.name.suffix)) {
                                                                                    file.name = paste(project,'.',file.prename,'.','plot_cells_3d','.',redDim.method[redDim],'_',dim.xyz,'.','differential','.','with_trajectory', sep = '')
                                                                                } else {
                                                                                    if (is.na(file.name.suffix) | file.name.suffix == "") {
                                                                                        file.name = paste(project,'.',file.prename,'.','plot_cells_3d','.',redDim.method[redDim],'_',dim.xyz,'.','differential','.','with_trajectory', sep = '')
                                                                                    } else {
                                                                                        file.name = paste(project,'.',file.prename,'.','plot_cells_3d','.',redDim.method[redDim],'_',dim.xyz,'.','differential','.',file.name.suffix,'.','with_trajectory', sep = '')
                                                                                    }
                                                                                }
                                                                            }
                                                                            file.name = paste(file.name,'.',trj.prt.name, sep = '')
                                                                            file.name = paste(file.name,'.',grid.name, sep = '')
                                                                        }
                                                                        
                                                                        file.html = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential','/',file.folder,'/',file.name,'.html', sep = '')
                                                                        file.png = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential','/',file.folder,'/',file.name,'.png', sep = '')
                                                                        
                                                                        cds_3d_plot_obj <- plot_cells_3d(cds_3d.plt, reduction_method = redDim.method[redDim], dims = c(dim.x, dim.y, dim.z), cell_size = cell_size_3d, genes = marker.genes, norm_method = "log", show_trajectory_graph = show_trajectory)
                                                                        
                                                                        if (grid.plt == 0) {
                                                                            cds_3d_plot_obj <- cds_3d_plot_obj %>% plotly::layout(scene = scene)
                                                                        }
                                                                        
                                                                        # https://community.rstudio.com/t/save-viewer-object-rendered-in-rstudio-as-image/32796
                                                                        if (dim.comb == 1) { # only save D1-D2-D3 (no D2-D3-D4, D3-D4-D5, ...)
                                                                            saveWidget(cds_3d_plot_obj, file.html)
                                                                            #webshot(file.html, file.png)
                                                                        } else { # remove D2-D3-D4, D3-D4-D5, ...
                                                                            mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential', sep = '')
                                                                            if (file.exists(file.path(mainDir))) {
                                                                                unlink(mainDir, recursive = TRUE)
                                                                            }
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                        #}
                                    }
                                }
                            }
                            
                            
                            # As before, we can collect the trajectory-variable genes into modules:
                            #for (redDim in 1:length(redDim.method)) {
                            
                            dim.folder = paste(eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))),'D', sep = '')
                            
                            mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim], sep = '')
                            subDir <- dim.folder
                            if (! file.exists(file.path(mainDir, subDir))) {
                                dir.create(file.path(mainDir, subDir))
                            }
                            
                            mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder, sep = '')
                            if (eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) == 2) {
                                #differential.folder <- paste('D12','/','differential', sep = '')
                                differential.folder <- 'differential'
                            } else {
                                differential.folder <- 'differential'
                            }
                            subDir <- differential.folder
                            if (! file.exists(file.path(mainDir, subDir))) {
                                dir.create(file.path(mainDir, subDir))
                            }
                            
                            if (SUB == 0) {
                                differential.prefolder = differential.folder
                            } else {
                                mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',differential.folder, sep = '')
                                subDir <- SUB.prefolder
                                if (! file.exists(file.path(mainDir, subDir))) {
                                    dir.create(file.path(mainDir, subDir))
                                }
                                differential.prefolder = paste(differential.folder,'/',SUB.prefolder, sep = '')
                            }
                            
                            if (redDim.method[redDim] == "UMAP") {
                                if (exists("pr_deg_ids.principal_graph") & exists("pr_deg_names.principal_graph")) {
                                    if (length(pr_deg_ids.principal_graph) > 0) {
                                        differential.run = 1
                                    } else {
                                        differential.run = 0
                                    }
                                } else {
                                    differential.run = 0
                                }
                            } else if (redDim.method[redDim] == "tSNE") {
                                differential.run = 0
                            }
                            
                            if (differential.run == 1) {
                                # You can call find_gene_modules(), which essentially runs UMAP on the genes (as opposed to the cells) and then groups them into modules using Louvain community analysis:
                                # reduction_method: The dimensionality reduction method used to generate the lower dimensional space in which genes will be clustered. Currently only UMAP is supported.
                                # max_components: The number of dimensions in which to cluster genes into modules.
                                # resolution: Resolution parameter passed to Louvain. Can be a list. If so, this method will evaluate modularity at each resolution and use the one with the highest value.
                                gene_group_df <- find_gene_modules(cds[pr_deg_ids.principal_graph,], reduction_method = "UMAP", max_components = eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))), resolution = c(10^seq(-6,-1))) # The data frame gene_group_df contains a row for each gene and identifies the module it belongs to.
                                
                                pr_deg_ids.principal_graph.idx <- sapply(unlist(gene_group_df[,"id"]), FUN = function(X) which(pr_deg_ids.principal_graph %in% X))
                                gene_group_df[,"gene_short_name"] = pr_deg_names.principal_graph[pr_deg_ids.principal_graph.idx]
                                
                                row.orderidx = order(as.numeric(unlist(gene_group_df[,"supermodule"])), 
                                                     as.numeric(unlist(gene_group_df[,"module"])), 
                                                     decreasing = FALSE)
                                gene_group_df <- gene_group_df[row.orderidx, , drop = FALSE]
                                gene_group_df.origin <- gene_group_df
                                
                                if (is.null(genes_of_interest.input) & is.null(signature_genes.input)) {
                                    gene_group_df <- gene_group_df.origin
                                } else if (! is.null(genes_of_interest.input) & is.null(signature_genes.input)) {
                                    gene_group_df <- rbind(gene_group_df.origin, query.gene_group_df.all)
                                } else if (is.null(genes_of_interest.input) & ! is.null(signature_genes.input)) {
                                    gene_group_df <- rbind(gene_group_df.origin, signature.gene_group_df.all)
                                } else if (! is.null(genes_of_interest.input) & ! is.null(signature_genes.input)) {
                                    gene_group_df <- rbind(gene_group_df.origin, query.gene_group_df.all, signature.gene_group_df.all)
                                }
                                
                                if (! is.null(gene_group_df)) {
                                    write.csv(gene_group_df, file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',differential.prefolder,'/',project,'.',gsub("\\.","_",subdata.folder),'.',gsub("\\.","_",project.subname),'.',redDim.method[redDim],'.',file.prename,'.','differential','.','geneModule','.csv', sep = ''), row.names = TRUE, quote = FALSE)
                                }
                                
                                
                                # Here we plot the aggregate module scores within each group of cell types as annotated by Packer & Zhu et al:
                                if (length(cell_type) > 0) {
                                    for (color.plt in 1:length(cell_type)) {
                                        mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',differential.prefolder, sep = '')
                                        subDir <- gsub("\\.","_",cell_type[color.plt])
                                        if (! file.exists(file.path(mainDir, subDir))) {
                                            dir.create(file.path(mainDir, subDir))
                                        }
                                        
                                        cell_group_df <- tibble::tibble(cell = rownames(colData(cds)), cell_group = colData(cds)[, cell_type[color.plt]])
                                        
                                        agg_mat <- aggregate_gene_expression(cds, gene_group_df, cell_group_df)
                                        
                                        a = grep(paste(c(genes_of_interest.titles, signature_genes.titles), collapse = '|'), rownames(agg_mat), value = TRUE, invert = TRUE)
                                        a.idx = grep(paste(c(genes_of_interest.titles, signature_genes.titles), collapse = '|'), rownames(agg_mat), value = FALSE, invert = TRUE)
                                        if (length(grep(paste(genes_of_interest.titles, collapse = '|'), rownames(agg_mat), value = TRUE, invert = FALSE)) > 0) {
                                            b.genes_of_interest = grep(paste(genes_of_interest.titles, collapse = '|'), rownames(agg_mat), value = TRUE, invert = FALSE)
                                            b.genes_of_interest.idx = grep(paste(genes_of_interest.titles, collapse = '|'), rownames(agg_mat), value = FALSE, invert = FALSE)
                                        } else {
                                            b.genes_of_interest = NULL
                                            b.genes_of_interest.idx = NULL
                                        }
                                        if (length(grep(paste(signature_genes.titles, collapse = '|'), rownames(agg_mat), value = TRUE, invert = FALSE)) > 0) {
                                            b.signature_genes = grep(paste(signature_genes.titles, collapse = '|'), rownames(agg_mat), value = TRUE, invert = FALSE)
                                            b.signature_genes.idx = grep(paste(signature_genes.titles, collapse = '|'), rownames(agg_mat), value = FALSE, invert = FALSE)
                                        } else {
                                            b.signature_genes = NULL
                                            b.signature_genes.idx = NULL
                                        }
                                        row.ordername = c(a[mixedorder(a, decreasing = FALSE)], 
                                                          names(which(sapply(genes_of_interest.titles, FUN = function(X) X %in% b.genes_of_interest))), 
                                                          names(which(sapply(signature_genes.titles, FUN = function(X) X %in% b.signature_genes))))
                                        row.orderidx = sapply(row.ordername, FUN = function(X) which(rownames(agg_mat) %in% X))
                                        
                                        if (project.name %in% c("Pancreas", "DrugPairing", "CRISPRINT")) {
                                            cell_type.name <- c()
                                            for (stg in 1:length(levels(colData(cds)[, stage.input]))) {
                                                stage.name = as.character(unique(colData(cds)[which(colData(cds)[, stage.input] == levels(colData(cds)[, stage.input])[stg]), cell_type[color.plt]]))
                                                cell_type.name <- c(cell_type.name, stage.name)
                                            }
                                            cell_type.name <- unique(cell_type.name)
                                        } else if (project.name %in% c("Paper")) {
                                            cell_type.name <- mixedsort(as.character(unique(colData(cds)[, cell_type[color.plt]])))
                                            cell_type.name <- cell_type.name[which(!is.na(cell_type.name))]
                                            
                                            if (exists("project.subname")) {
                                                project.subname.prename = unlist(strsplit(project.subname, '\\.'))[1]
                                                
                                                if (project.subname.prename == "GSE147424" | project.subname.prename == "SRP253773") {
                                                    if (project.subname.prename == "GSE147424") {
                                                        sample_type.LV <- c("H", "NL", "LS")
                                                    } else if (project.subname.prename == "SRP253773") {
                                                        sample_type.LV <- c("Healthy", "Non-lesional", "Lesional")
                                                    }
                                                    sample_type.LV.idx <- which(sapply(sample_type.LV, FUN = function(X) X %in% cell_type.name))
                                                    
                                                    if (length(sample_type.LV.idx) > 0) {
                                                        sample_type.LV = sample_type.LV[sample_type.LV.idx]
                                                        cell_type.name <- sample_type.LV
                                                    }
                                                }
                                            }
                                        } else {
                                            cell_type.name <- as.character(unique(colData(cds)[, cell_type[color.plt]]))
                                            cell_type.name <- cell_type.name[which(!is.na(cell_type.name))]
                                        }
                                        if (project.name == "tutorial") {
                                            col.orderidx = mixedorder(colnames(agg_mat), decreasing = FALSE)
                                        } else {
                                            col.orderidx = sapply(cell_type.name, FUN = function(X) which(colnames(agg_mat) %in% X))
                                        }
                                        
                                        agg_mat <- agg_mat[row.orderidx, col.orderidx, drop = FALSE]
                                        
                                        rownames(agg_mat)[a.idx] <- stringr::str_c("Module ", rownames(agg_mat)[a.idx])
                                        if (! is.null(genes_of_interest.input)) {
                                            for (gX in 1:length(genes_of_interest.input)) {
                                                rownames(agg_mat)[which(rownames(agg_mat) == paste("geneSet",gX,"_of_interest", sep = ''))] = paste("geneSet",gX, sep = '')
                                            }
                                        }
                                        rownames(agg_mat)[which(rownames(agg_mat) == "SINGH_KRAS_DEPENDENCY_SIGNATURE")] = "SINGH"
                                        rownames(agg_mat)[which(rownames(agg_mat) == "FRIDMAN_SENESCENCE_UP")] = "FRIDMAN"
                                        rownames(agg_mat)[which(rownames(agg_mat) == "HALLMARK_KRAS_SIGNALING_UP")] = "HALLMARK"
                                        
                                        colnames(agg_mat)[which(colnames(agg_mat) == "Acinar cell")] = "Acinar"
                                        colnames(agg_mat)[which(colnames(agg_mat) == "Ductal cell type 1")] = "Ductal 1"
                                        colnames(agg_mat)[which(colnames(agg_mat) == "Ductal cell type 2")] = "Ductal 2"
                                        colnames(agg_mat)[which(colnames(agg_mat) == "Endocrine cell")] = "Endocrine"
                                        colnames(agg_mat)[which(colnames(agg_mat) == "Endothelial cell")] = "Endothelial"
                                        colnames(agg_mat)[which(colnames(agg_mat) == "Fibroblast cell")] = "Fibroblast"
                                        colnames(agg_mat)[which(colnames(agg_mat) == "Macrophage cell")] = "Macrophage"
                                        colnames(agg_mat)[which(colnames(agg_mat) == "Stellate cell")] = "Stellate"
                                        colnames(agg_mat)[which(colnames(agg_mat) == "B cell")] = "B"
                                        colnames(agg_mat)[which(colnames(agg_mat) == "T cell")] = "T"
                                        
                                        if (exists("project.subname")) {
                                            project.subname.prename = unlist(strsplit(project.subname, '\\.'))[1]
                                            
                                            if (project.subname.prename == "GSE147424" | project.subname.prename == "SRP253773") {
                                                agg_mat.new <- agg_mat
                                                
                                                colnames(agg_mat.new) = gsub("sample|sample |Sample|Sample ","S",colnames(agg_mat.new))
                                                
                                                colnames(agg_mat.new)[which(colnames(agg_mat.new) == "Healthy")] = "H"
                                                colnames(agg_mat.new)[which(colnames(agg_mat.new) == "Non-lesional")] = "NL"
                                                colnames(agg_mat.new)[which(colnames(agg_mat.new) == "Lesional")] = "LS"
                                                
                                                agg_mat <- agg_mat.new
                                            } else if (project.subname.prename == "NG-27918") {
                                                agg_mat.new <- agg_mat
                                                
                                                #colnames(agg_mat.new) = gsub("sample|sample |Sample|Sample ","S",colnames(agg_mat.new))
                                                
                                                colnames(agg_mat.new) = gsub("unstimulated","unsti",colnames(agg_mat.new))
                                                colnames(agg_mat.new) = gsub("stimulated","sti",colnames(agg_mat.new))
                                                
                                                agg_mat <- agg_mat.new
                                            }
                                        }
                                        
                                        colnames(agg_mat) <- stringr::str_c(colnames(agg_mat))
                                        
                                        agg_mat <- as.data.frame(agg_mat)
                                        
                                        if (length(grep("CRA001160", project, value = TRUE)) > 0) {
                                            if (cell_type[color.plt] == "sample") {
                                                diagnoses2sample <- sapply(levels(cell_metadata[,"sample_pathologic_diagnoses"]), FUN = function(X) mixedsort(unique(cell_metadata[which(cell_metadata[,"sample_pathologic_diagnoses"] %in% X),"sample"])))
                                                diagnoses2sample <- unlist(diagnoses2sample)
                                                diagnoses2sample <- factor(diagnoses2sample, levels = unique(diagnoses2sample))
                                                
                                                agg_mat <- agg_mat[, as.character(diagnoses2sample)]
                                            }
                                        }
                                        
                                        write.csv(agg_mat, file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',differential.prefolder,'/',gsub("\\.","_",cell_type[color.plt]),'/',project,'.',gsub("\\.","_",subdata.folder),'.',gsub("\\.","_",project.subname),'.',redDim.method[redDim],'.',file.prename,'.','differential','.','geneModule','.',gsub("\\.","_",cell_type[color.plt]),'.csv', sep = ''), row.names = TRUE, quote = FALSE)
                                        
                                        
                                        # /// pheatmap ///
                                        for (clu in 1:4) {
                                            
                                            if (clu == 1) {
                                                row.clustering = FALSE
                                                col.clustering = FALSE
                                            } else if (clu == 2) {
                                                row.clustering = FALSE
                                                col.clustering = TRUE
                                            } else if (clu == 3) {
                                                row.clustering = TRUE
                                                col.clustering = FALSE
                                            } else if (clu == 4) {
                                                row.clustering = TRUE
                                                col.clustering = TRUE
                                            }
                                            
                                            plt = 0
                                            if (col.clustering == FALSE) {
                                                plt = 1
                                            } else if (col.clustering == TRUE) {
                                                if (ncol(agg_mat) >= 2) {
                                                    plt = 1
                                                }
                                            }
                                            
                                            if (plt == 1) {
                                                if (row.clustering == FALSE & col.clustering == FALSE) {
                                                    suffix.name = paste('geneModule','2',gsub("\\.","_",cell_type[color.plt]), sep = '')
                                                    width.size = (ncol(agg_mat) / 4)
                                                    height.size = (nrow(agg_mat) / 10)
                                                } else if (row.clustering == FALSE & col.clustering == TRUE) {
                                                    suffix.name = paste('geneModule','2',gsub("\\.","_",cell_type[color.plt]),'_','clustering', sep = '')
                                                    width.size = (ncol(agg_mat) / 4)
                                                    height.size = (nrow(agg_mat) / 10) + 0.5
                                                } else if (row.clustering == TRUE & col.clustering == FALSE) {
                                                    suffix.name = paste('geneModule','_','clustering','2',gsub("\\.","_",cell_type[color.plt]), sep = '')
                                                    width.size = (ncol(agg_mat) / 4) + 0.625
                                                    height.size = (nrow(agg_mat) / 10)
                                                } else if (row.clustering == TRUE & col.clustering == TRUE) {
                                                    suffix.name = paste('geneModule','_','clustering','2',gsub("\\.","_",cell_type[color.plt]),'_','clustering', sep = '')
                                                    width.size = (ncol(agg_mat) / 4) + 0.625
                                                    height.size = (nrow(agg_mat) / 10) + 0.5
                                                }
                                                
                                                if ((ncol(agg_mat) / 4) <= 1.5) {
                                                    width.size.fold = 1.5 / (ncol(agg_mat) / 4)
                                                    width.size = width.size * width.size.fold
                                                }
                                                
                                                fontsize = 6
                                                
                                                png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',differential.prefolder,'/',gsub("\\.","_",cell_type[color.plt]),'/',project,'.',gsub("\\.","_",subdata.folder),'.',gsub("\\.","_",project.subname),'.',redDim.method[redDim],'.',file.prename,'.','differential','.','pheatmap','.',suffix.name,'.png', sep = ''), width = width.size, height = height.size, units = 'in', res = 600)
                                                p <- pheatmap::pheatmap(agg_mat, cluster_rows = row.clustering, cluster_cols = col.clustering, 
                                                                        scale = "column", # scale: character indicating if the values should be centered and scaled in either the row direction or the column direction, or none. Corresponding values are "row", "column" and "none"
                                                                        clustering_method = "ward.D2", # hclust(): method: the agglomeration method to be used. This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
                                                                        fontsize = fontsize)
                                                print(p)
                                                dev.off()
                                                
                                                if (row.clustering == TRUE & col.clustering == FALSE) {
                                                    p.tree_row.order <- p$tree_row$order
                                                }
                                            }
                                        }
                                        
                                        
                                        if (ncol(agg_mat) > 1) {
                                            # --- spaghetti -------------------------------------------------------------------------------------------------------------------------------------
                                            # THE SPAGHETTI PLOT: https://www.data-to-viz.com/caveat/spaghetti.html
                                            for (od in 0:1) {
                                                if (od == 0) {
                                                    agg_mat.differential <- agg_mat
                                                } else if (od == 1) {
                                                    agg_mat.differential <- agg_mat[p.tree_row.order, , drop = FALSE]
                                                }
                                                
                                                agg_mat.differential <- data.frame("module" = rownames(agg_mat.differential), agg_mat.differential, check.names = FALSE)
                                                colnames(agg_mat.differential)[1] = "module"
                                                
                                                if (length(grep("CRA001160", project, value = TRUE)) > 0) {
                                                    if (cell_type[color.plt] == "sample") {
                                                        diagnoses2sample <- sapply(levels(cell_metadata[,"sample_pathologic_diagnoses"]), FUN = function(X) mixedsort(unique(cell_metadata[which(cell_metadata[,"sample_pathologic_diagnoses"] %in% X),"sample"])))
                                                        diagnoses2sample <- unlist(diagnoses2sample)
                                                        diagnoses2sample <- factor(diagnoses2sample, levels = unique(diagnoses2sample))
                                                        
                                                        agg_mat.differential <- agg_mat.differential[, c("module", as.character(diagnoses2sample))]
                                                    }
                                                }
                                                
                                                agg_mat.differential.origin <- agg_mat.differential
                                                
                                                agg_mat.differential.mtx <- agg_mat.differential
                                                rownames(agg_mat.differential.mtx) = agg_mat.differential.mtx[,"module"]
                                                module.idx <- which(colnames(agg_mat.differential.mtx) %in% "module")
                                                agg_mat.differential.mtx <- agg_mat.differential.mtx[, -module.idx, drop = FALSE]
                                                agg_mat.differential.mtx.origin <- agg_mat.differential.mtx
                                                
                                                
                                                # --- trend -----------------------------------------------------------------------------------------------------------------------------------------
                                                agg_mat.differential <- agg_mat.differential.origin
                                                agg_mat.differential.mtx <- agg_mat.differential.mtx.origin
                                                
                                                if (length(grep("CRA001160", project, value = TRUE)) > 0) {
                                                    if (cell_type[color.plt] == "sample") {
                                                        sample_type2sample <- sapply(levels(cell_metadata[,"sample_type"]), FUN = function(X) mixedsort(unique(cell_metadata[which(cell_metadata[,"sample_type"] %in% X),"sample"])))
                                                        
                                                        col.idx1 <- unlist(sapply(sample_type2sample, FUN = function(X) which(colnames(agg_mat.differential.mtx) %in% X))[1])
                                                        col.idx2 <- unlist(sapply(sample_type2sample, FUN = function(X) which(colnames(agg_mat.differential.mtx) %in% X))[2])
                                                        
                                                        col.idx.order <- col.idx2
                                                    } else {
                                                        col.idx1 <- NULL
                                                        col.idx2 <- NULL
                                                        
                                                        col.idx.order <- c(1:ncol(agg_mat.differential.mtx))
                                                    }
                                                } else {
                                                    col.idx1 <- NULL
                                                    col.idx2 <- NULL
                                                    
                                                    col.idx.order <- c(1:ncol(agg_mat.differential.mtx))
                                                }
                                                
                                                low2high.idx <- c()
                                                increasing.idx <- c()
                                                low2high_increasing.idx <- c()
                                                
                                                high2low.idx <- c()
                                                decreasing.idx <- c()
                                                high2low_decreasing.idx <- c()
                                                
                                                for (m in 1:nrow(agg_mat.differential.mtx)) {
                                                    col.idx.inc <- col.idx.order[order(agg_mat.differential.mtx[m, col.idx.order], decreasing = FALSE)]
                                                    col.idx.dec <- col.idx.order[order(agg_mat.differential.mtx[m, col.idx.order], decreasing = TRUE)]
                                                    
                                                    # /// increasing ///
                                                    if (! is.null(col.idx1) & ! is.null(col.idx2)) {
                                                        if (max(agg_mat.differential.mtx[m, col.idx1], na.rm = TRUE) <= min(agg_mat.differential.mtx[m, col.idx2], na.rm = TRUE)) {
                                                            #cat(rownames(agg_mat.differential.mtx[m, ]), sep = "\n")
                                                            
                                                            #cat("Group 1 <= Group 2\n")
                                                            low2high.idx <- c(low2high.idx, m)
                                                            
                                                            if (all(col.idx.inc == col.idx.order)) {
                                                                #cat("increasing\n")
                                                                increasing.idx <- c(increasing.idx, m)
                                                                
                                                                low2high_increasing.idx <- c(low2high_increasing.idx, m)
                                                            }
                                                            
                                                            #print(agg_mat.differential.mtx[m, ])
                                                            #cat("\n")
                                                            #cat("\n")
                                                        } else {
                                                            if (all(col.idx.inc == col.idx.order)) {
                                                                #cat(rownames(agg_mat.differential.mtx[m, ]), sep = "\n")
                                                                
                                                                #cat("increasing\n")
                                                                increasing.idx <- c(increasing.idx, m)
                                                                
                                                                #print(agg_mat.differential.mtx[m, col.idx.order])
                                                                #cat("\n")
                                                                #cat("\n")
                                                            }
                                                        }
                                                    } else {
                                                        if (all(col.idx.inc == col.idx.order)) {
                                                            #cat(rownames(agg_mat.differential.mtx[m, ]), sep = "\n")
                                                            
                                                            #cat("increasing\n")
                                                            increasing.idx <- c(increasing.idx, m)
                                                            
                                                            #print(agg_mat.differential.mtx[m, col.idx.order])
                                                            #cat("\n")
                                                            #cat("\n")
                                                        }
                                                    }
                                                    
                                                    # /// decreasing ///
                                                    if (! is.null(col.idx1) & ! is.null(col.idx2)) {
                                                        if (min(agg_mat.differential.mtx[m, col.idx1], na.rm = TRUE) >= max(agg_mat.differential.mtx[m, col.idx2], na.rm = TRUE)) {
                                                            #cat(rownames(agg_mat.differential.mtx[m, ]), sep = "\n")
                                                            
                                                            #cat("Group 1 >= Group 2\n")
                                                            high2low.idx <- c(high2low.idx, m)
                                                            
                                                            if (all(col.idx.dec == col.idx.order)) {
                                                                #cat("decreasing\n")
                                                                decreasing.idx <- c(decreasing.idx, m)
                                                                
                                                                high2low_decreasing.idx <- c(high2low_decreasing.idx, m)
                                                            }
                                                            
                                                            #print(agg_mat.differential.mtx[m, ])
                                                            #cat("\n")
                                                            #cat("\n")
                                                        } else {
                                                            if (all(col.idx.dec == col.idx.order)) {
                                                                #cat(rownames(agg_mat.differential.mtx[m, ]), sep = "\n")
                                                                
                                                                #cat("decreasing\n")
                                                                decreasing.idx <- c(decreasing.idx, m)
                                                                
                                                                #print(agg_mat.differential.mtx[m, col.idx.order])
                                                                #cat("\n")
                                                                #cat("\n")
                                                            }
                                                        }
                                                    } else {
                                                        if (all(col.idx.dec == col.idx.order)) {
                                                            #cat(rownames(agg_mat.differential.mtx[m, ]), sep = "\n")
                                                            
                                                            #cat("decreasing\n")
                                                            decreasing.idx <- c(decreasing.idx, m)
                                                            
                                                            #print(agg_mat.differential.mtx[m, col.idx.order])
                                                            #cat("\n")
                                                            #cat("\n")
                                                        }
                                                    }
                                                }
                                                
                                                if (od == 0) {
                                                    #cat("Group 1 <= Group 2\n")
                                                    #print(low2high.idx)
                                                    #cat("increasing\n")
                                                    #print(increasing.idx)
                                                    #cat("Group 1 <= Group 2 & increasing\n")
                                                    #print(low2high_increasing.idx)
                                                    #cat("\n")
                                                    
                                                    #cat("Group 1 >= Group 2\n")
                                                    #print(high2low.idx)
                                                    #cat("decreasing\n")
                                                    #print(decreasing.idx)
                                                    #cat("Group 1 >= Group 2 & decreasing\n")
                                                    #print(high2low_decreasing.idx)
                                                    #cat("\n")
                                                    
                                                    if (! is.null(low2high.idx)) {
                                                        agg_mat.differential.mtx.low2high <- agg_mat.differential.mtx[low2high.idx, , drop = FALSE]
                                                        write.csv(agg_mat.differential.mtx.low2high, file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',differential.prefolder,'/',gsub("\\.","_",cell_type[color.plt]),'/',project,'.',gsub("\\.","_",subdata.folder),'.',gsub("\\.","_",project.subname),'.',redDim.method[redDim],'.',file.prename,'.','differential','.','geneModule','.',gsub("\\.","_",cell_type[color.plt]),'.','low2high','.csv', sep = ''), row.names = TRUE, quote = FALSE)
                                                    }
                                                    if (! is.null(increasing.idx)) {
                                                        agg_mat.differential.mtx.inc <- agg_mat.differential.mtx[increasing.idx, , drop = FALSE]
                                                        write.csv(agg_mat.differential.mtx.inc, file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',differential.prefolder,'/',gsub("\\.","_",cell_type[color.plt]),'/',project,'.',gsub("\\.","_",subdata.folder),'.',gsub("\\.","_",project.subname),'.',redDim.method[redDim],'.',file.prename,'.','differential','.','geneModule','.',gsub("\\.","_",cell_type[color.plt]),'.','increasing','.csv', sep = ''), row.names = TRUE, quote = FALSE)
                                                    }
                                                    if (! is.null(low2high_increasing.idx)) {
                                                        agg_mat.differential.mtx.low2high_inc <- agg_mat.differential.mtx[low2high_increasing.idx, , drop = FALSE]
                                                        write.csv(agg_mat.differential.mtx.low2high_inc, file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',differential.prefolder,'/',gsub("\\.","_",cell_type[color.plt]),'/',project,'.',gsub("\\.","_",subdata.folder),'.',gsub("\\.","_",project.subname),'.',redDim.method[redDim],'.',file.prename,'.','differential','.','geneModule','.',gsub("\\.","_",cell_type[color.plt]),'.','low2high_increasing','.csv', sep = ''), row.names = TRUE, quote = FALSE)
                                                    }
                                                    
                                                    if (! is.null(high2low.idx)) {
                                                        agg_mat.differential.mtx.high2low <- agg_mat.differential.mtx[high2low.idx, , drop = FALSE]
                                                        write.csv(agg_mat.differential.mtx.high2low, file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',differential.prefolder,'/',gsub("\\.","_",cell_type[color.plt]),'/',project,'.',gsub("\\.","_",subdata.folder),'.',gsub("\\.","_",project.subname),'.',redDim.method[redDim],'.',file.prename,'.','differential','.','geneModule','.',gsub("\\.","_",cell_type[color.plt]),'.','high2low','.csv', sep = ''), row.names = TRUE, quote = FALSE)
                                                    }
                                                    if (! is.null(decreasing.idx)) {
                                                        agg_mat.differential.mtx.dec <- agg_mat.differential.mtx[decreasing.idx, , drop = FALSE]
                                                        write.csv(agg_mat.differential.mtx.dec, file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',differential.prefolder,'/',gsub("\\.","_",cell_type[color.plt]),'/',project,'.',gsub("\\.","_",subdata.folder),'.',gsub("\\.","_",project.subname),'.',redDim.method[redDim],'.',file.prename,'.','differential','.','geneModule','.',gsub("\\.","_",cell_type[color.plt]),'.','decreasing','.csv', sep = ''), row.names = TRUE, quote = FALSE)
                                                    }
                                                    if (! is.null(high2low_decreasing.idx)) {
                                                        agg_mat.differential.mtx.high2low_dec <- agg_mat.differential.mtx[high2low_decreasing.idx, , drop = FALSE]
                                                        write.csv(agg_mat.differential.mtx.high2low_dec, file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',differential.prefolder,'/',gsub("\\.","_",cell_type[color.plt]),'/',project,'.',gsub("\\.","_",subdata.folder),'.',gsub("\\.","_",project.subname),'.',redDim.method[redDim],'.',file.prename,'.','differential','.','geneModule','.',gsub("\\.","_",cell_type[color.plt]),'.','high2low_decreasing','.csv', sep = ''), row.names = TRUE, quote = FALSE)
                                                    }
                                                }
                                                
                                                
                                                # /// spaghetti plot ///
                                                if (SUB == 0) {
                                                    main.title = paste(project.name, gsub("\\.","_",project.subname), sep = '.')
                                                } else {
                                                    main.title = paste(paste(project.name, gsub("\\.","_",project.subname), sep = '.'),'.',SUB.prename, sep = '')
                                                }
                                                
                                                xlab.name = unlist(strsplit(cell_type[color.plt], '_'))[1]
                                                xlab.subname <- c()
                                                if (length(unlist(strsplit(cell_type[color.plt], '_'))) > 1) {
                                                    xlab.subname <- c()
                                                    for (xn in 2:length(unlist(strsplit(cell_type[color.plt], '_')))) {
                                                        xlab.subname <- c(xlab.subname, unlist(strsplit(cell_type[color.plt], '_'))[xn])
                                                    }
                                                    xlab.subname <- paste(xlab.subname, collapse = '_')
                                                    xlab.name = paste(xlab.name, xlab.subname, sep = ' ')
                                                }
                                                xlab.title = str_to_title(xlab.name)
                                                
                                                for (trend in 0:6) {
                                                    agg_mat.differential <- agg_mat.differential.origin
                                                    agg_mat.differential.mtx <- agg_mat.differential.mtx.origin
                                                    
                                                    if (trend == 0) {
                                                        agg_mat.differential <- agg_mat.differential.origin
                                                        agg_mat.differential.mtx <- agg_mat.differential.mtx.origin
                                                        suffix.name = gsub("\\.","_",cell_type[color.plt])
                                                    } else if (trend == 1) {
                                                        if (! is.null(low2high.idx)) {
                                                            agg_mat.differential <- agg_mat.differential.origin[low2high.idx, , drop = FALSE]
                                                            agg_mat.differential.mtx <- agg_mat.differential.mtx.origin[low2high.idx, , drop = FALSE]
                                                            suffix.name = paste(gsub("\\.","_",cell_type[color.plt]),'.','low2high', sep = '')
                                                        } else {
                                                            next # next halts the processing of the current iteration and advances the looping index.
                                                        }
                                                    } else if (trend == 2) {
                                                        if (! is.null(increasing.idx)) {
                                                            agg_mat.differential <- agg_mat.differential.origin[increasing.idx, , drop = FALSE]
                                                            agg_mat.differential.mtx <- agg_mat.differential.mtx.origin[increasing.idx, , drop = FALSE]
                                                            suffix.name = paste(gsub("\\.","_",cell_type[color.plt]),'.','increasing', sep = '')
                                                        } else {
                                                            next # next halts the processing of the current iteration and advances the looping index.
                                                        }
                                                    } else if (trend == 3) {
                                                        if (! is.null(low2high_increasing.idx)) {
                                                            agg_mat.differential <- agg_mat.differential.origin[low2high_increasing.idx, , drop = FALSE]
                                                            agg_mat.differential.mtx <- agg_mat.differential.mtx.origin[low2high_increasing.idx, , drop = FALSE]
                                                            suffix.name = paste(gsub("\\.","_",cell_type[color.plt]),'.','low2high_increasing', sep = '')
                                                        } else {
                                                            next # next halts the processing of the current iteration and advances the looping index.
                                                        }
                                                    } else if (trend == 4) {
                                                        if (! is.null(high2low.idx)) {
                                                            agg_mat.differential <- agg_mat.differential.origin[high2low.idx, , drop = FALSE]
                                                            agg_mat.differential.mtx <- agg_mat.differential.mtx.origin[high2low.idx, , drop = FALSE]
                                                            suffix.name = paste(gsub("\\.","_",cell_type[color.plt]),'.','high2low', sep = '')
                                                        } else {
                                                            next # next halts the processing of the current iteration and advances the looping index.
                                                        }
                                                    } else if (trend == 5) {
                                                        if (! is.null(decreasing.idx)) {
                                                            agg_mat.differential <- agg_mat.differential.origin[decreasing.idx, , drop = FALSE]
                                                            agg_mat.differential.mtx <- agg_mat.differential.mtx.origin[decreasing.idx, , drop = FALSE]
                                                            suffix.name = paste(gsub("\\.","_",cell_type[color.plt]),'.','decreasing', sep = '')
                                                        } else {
                                                            next # next halts the processing of the current iteration and advances the looping index.
                                                        }
                                                    } else if (trend == 6) {
                                                        if (! is.null(high2low_decreasing.idx)) {
                                                            agg_mat.differential <- agg_mat.differential.origin[high2low_decreasing.idx, , drop = FALSE]
                                                            agg_mat.differential.mtx <- agg_mat.differential.mtx.origin[high2low_decreasing.idx, , drop = FALSE]
                                                            suffix.name = paste(gsub("\\.","_",cell_type[color.plt]),'.','high2low_decreasing', sep = '')
                                                        } else {
                                                            next # next halts the processing of the current iteration and advances the looping index.
                                                        }
                                                    }
                                                    
                                                    agg_mat.differential <- suppressMessages(melt(agg_mat.differential))
                                                    
                                                    agg_mat.differential[,"module"] <- factor(agg_mat.differential[,"module"], levels = unique(agg_mat.differential[,"module"]))
                                                    agg_mat.differential[,"variable"] <- factor(agg_mat.differential[,"variable"], levels = unique(agg_mat.differential[,"variable"]))
                                                    
                                                    colnames(agg_mat.differential) = str_to_title(colnames(agg_mat.differential))
                                                    
                                                    if (od == 0) {
                                                        legend.num = nrow(agg_mat.differential.mtx)
                                                        legend.ncol = ceiling(legend.num / 18)
                                                        legend.nrow = ceiling(legend.num / legend.ncol)
                                                        width.size = 10 + (0.5 * legend.ncol)
                                                        if (legend.nrow > 15) {
                                                            height.size = 5 + ((legend.nrow - 15) * 0.2)
                                                        } else {
                                                            height.size = 5
                                                        }
                                                        
                                                        png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',differential.prefolder,'/',gsub("\\.","_",cell_type[color.plt]),'/',project,'.',gsub("\\.","_",subdata.folder),'.',gsub("\\.","_",project.subname),'.',redDim.method[redDim],'.',file.prename,'.','differential','.','geneModule','.',suffix.name,'.png', sep = ''), width = width.size, height = height.size, units = 'in', res = 600)
                                                        p <- ggplot(agg_mat.differential, aes(x=Variable, y=Value, group=Module, color=Module)) + geom_line()
                                                        p <- p + theme_light()
                                                        #p <- p + theme_ipsum()
                                                        p <- p + ggtitle(main.title)
                                                        p <- p + xlab(xlab.title)
                                                        print(p)
                                                        dev.off()
                                                    }
                                                    
                                                    legend.num = nrow(agg_mat.differential.mtx)
                                                    legend.ncol = ceiling(sqrt(legend.num))
                                                    legend.nrow = ceiling(legend.num / legend.ncol)
                                                    Variable.num <- ncol(agg_mat.differential.mtx)
                                                    Value.range <- ceiling(max(agg_mat.differential.mtx, na.rm = TRUE)) - floor(min(agg_mat.differential.mtx, na.rm = TRUE))
                                                    width.size = legend.ncol * 2 + ceiling(Variable.num / 4) * 2 + 2
                                                    height.size = legend.nrow * 2 + ceiling(Value.range / 4) * 2 - 2
                                                    if (height.size < 5) {
                                                        if (legend.nrow <= 1) {
                                                            height.size = 4
                                                        } else {
                                                            height.size = 5
                                                        }
                                                    }
                                                    
                                                    if (od == 0) {
                                                        png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',differential.prefolder,'/',gsub("\\.","_",cell_type[color.plt]),'/',project,'.',gsub("\\.","_",subdata.folder),'.',gsub("\\.","_",project.subname),'.',redDim.method[redDim],'.',file.prename,'.','differential','.','geneModule','.',suffix.name,'.spaghetti','.png', sep = ''), width = width.size, height = height.size, units = 'in', res = 600)
                                                    } else if (od == 1) {
                                                        png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',differential.prefolder,'/',gsub("\\.","_",cell_type[color.plt]),'/',project,'.',gsub("\\.","_",subdata.folder),'.',gsub("\\.","_",project.subname),'.',redDim.method[redDim],'.',file.prename,'.','differential','.','geneModule','.',suffix.name,'.spaghetti','.','geneModule','_','clustering','.png', sep = ''), width = width.size, height = height.size, units = 'in', res = 600)
                                                    }
                                                    data <- agg_mat.differential
                                                    p <- data %>% mutate(Module2=Module)
                                                    p <- p %>% 
                                                        ggplot(aes(x=Variable, y=Value, group=Module)) + 
                                                        geom_line(data=p %>% dplyr::select(-Module), aes(group=Module2), color="grey", size=0.5, alpha=0.5) + 
                                                        geom_line(aes(color=Module), color="red", size=1.2) + 
                                                        scale_color_viridis(discrete = TRUE) + 
                                                        theme_ipsum() + 
                                                        theme(
                                                            legend.position="none", 
                                                            plot.title = element_text(size=22), 
                                                            panel.grid = element_blank()
                                                        ) + 
                                                        ggtitle(main.title) + 
                                                        facet_wrap(~Module)
                                                    p <- p + xlab(xlab.title)
                                                    p <- p + 
                                                        theme(
                                                            strip.text = element_text(size=20), 
                                                            axis.title.x = element_text(size=18), 
                                                            axis.title.y = element_text(size=18), 
                                                            axis.text.x = element_text(size=14), 
                                                            axis.text.y = element_text(size=14)
                                                        )
                                                    print(p)
                                                    dev.off()
                                                    
                                                    agg_mat.differential <- agg_mat.differential.origin
                                                    agg_mat.differential.mtx <- agg_mat.differential.mtx.origin
                                                }
                                            }
                                        }
                                    }
                                }
                            } else if (differential.run == 0) {
                                if (is.null(genes_of_interest.input) & is.null(signature_genes.input)) {
                                    gene_group_df <- NULL
                                } else if (! is.null(genes_of_interest.input) & is.null(signature_genes.input)) {
                                    gene_group_df <- query.gene_group_df.all
                                } else if (is.null(genes_of_interest.input) & ! is.null(signature_genes.input)) {
                                    gene_group_df <- signature.gene_group_df.all
                                } else if (! is.null(genes_of_interest.input) & ! is.null(signature_genes.input)) {
                                    gene_group_df <- rbind(query.gene_group_df.all, signature.gene_group_df.all)
                                }
                                
                                if (! is.null(gene_group_df)) {
                                    write.csv(gene_group_df, file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',differential.prefolder,'/',project,'.',gsub("\\.","_",subdata.folder),'.',gsub("\\.","_",project.subname),'.',redDim.method[redDim],'.',file.prename,'.','differential','.','geneModule','.csv', sep = ''), row.names = TRUE, quote = FALSE)
                                }
                            }
                            
                            mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',differential.prefolder, sep = '')
                            if (file.exists(file.path(mainDir))) {
                                if (length(list.files(path = mainDir, full.names = TRUE, recursive = FALSE)) == 0) {
                                    unlink(mainDir, recursive = TRUE)
                                }
                            }
                            
                            mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',differential.folder, sep = '')
                            if (file.exists(file.path(mainDir))) {
                                if (length(list.files(path = mainDir, full.names = TRUE, recursive = FALSE)) == 0) {
                                    unlink(mainDir, recursive = TRUE)
                                }
                            }
                            
                            mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder, sep = '')
                            if (file.exists(file.path(mainDir))) {
                                if (length(list.files(path = mainDir, full.names = TRUE, recursive = FALSE)) == 0) {
                                    unlink(mainDir, recursive = TRUE)
                                }
                            }
                            #}
                            
                            
                            # We can also pass gene_group_df to plot_cells() as we did when we compared clusters in the L2 data above.
                            if (! is.null(gene_group_df)) {
                                
                                #module.type = c("module", "supermodule")
                                module.type = c("module")
                                
                                for (M in 1:length(module.type)) {
                                    
                                    # /// marker.genes ///
                                    module.levels <- levels(unlist(gene_group_df[,module.type[M]]))
                                    
                                    module_set.list <- list()
                                    
                                    for (LV in 1:length(module.levels)) {
                                        module.levels.LV <- list(module.levels[LV])
                                        if (module.levels[LV] %in% c(genes_of_interest.titles, signature_genes.titles)) {
                                            names(module.levels.LV) = module.levels[LV]
                                        } else {
                                            names(module.levels.LV) = paste(module.type[M],module.levels[LV], sep = '')
                                        }
                                        module_set.list <- c(module_set.list, module.levels.LV)
                                    }
                                    
                                    # *** query gene(s) ***
                                    if (! is.null(genes_of_interest.input)) {
                                        for (gX in 1:length(genes_of_interest.input)) {
                                            genes_of_interest.genes <- genes_of_interest.input[[gX]]
                                            genes_of_interest.genes.idx <- which(sapply(genes_of_interest.genes, FUN = function(X) X %in% unlist(gene_group_df[,"gene_short_name"])))
                                            if (length(genes_of_interest.genes.idx) > 0) {
                                                genes_of_interest.genes = genes_of_interest.genes[genes_of_interest.genes.idx]
                                                
                                                genes_of_interest.genes.idx <- lapply(genes_of_interest.genes, FUN = function(X) which(unlist(gene_group_df[,"gene_short_name"]) %in% X))
                                                names(genes_of_interest.genes.idx) = genes_of_interest.genes
                                                
                                                marker.genes = gene_group_df[unlist(genes_of_interest.genes.idx),]
                                                
                                                marker.modules <- marker.genes[,module.type[M]]
                                                genes_of_interest.modules <- as.character(unlist(unique(marker.modules)))
                                                genes_of_interest.modules <- genes_of_interest.modules[which(! genes_of_interest.modules %in% genes_of_interest.input.title[which(! genes_of_interest.input.title %in% names(genes_of_interest.input[gX]))])]
                                                if (! is.null(genes_of_interest.modules)) {
                                                    a = grep(paste(genes_of_interest.titles, collapse = '|'), genes_of_interest.modules, value = TRUE, invert = TRUE)
                                                    a.idx = grep(paste(genes_of_interest.titles, collapse = '|'), genes_of_interest.modules, value = FALSE, invert = TRUE)
                                                    if (length(grep(paste(genes_of_interest.titles, collapse = '|'), genes_of_interest.modules, value = TRUE, invert = FALSE)) > 0) {
                                                        b.genes_of_interest = grep(paste(genes_of_interest.titles, collapse = '|'), genes_of_interest.modules, value = TRUE, invert = FALSE)
                                                        b.genes_of_interest.idx = grep(paste(genes_of_interest.titles, collapse = '|'), genes_of_interest.modules, value = FALSE, invert = FALSE)
                                                    } else {
                                                        b.genes_of_interest = NULL
                                                        b.genes_of_interest.idx = NULL
                                                    }
                                                    genes_of_interest.modules = c(a[mixedorder(a, decreasing = FALSE)], 
                                                                                  names(which(sapply(genes_of_interest.titles, FUN = function(X) X %in% b.genes_of_interest))))
                                                    
                                                    for (gY in 1:length(genes_of_interest.modules)) {
                                                        genes_of_interest.modules.gX.gY <- list(genes_of_interest.modules[gY])
                                                        if (genes_of_interest.modules[gY] %in% genes_of_interest.titles) {
                                                            names(genes_of_interest.modules.gX.gY) = paste(names(genes_of_interest.input[gX]),'.',genes_of_interest.modules[gY], sep = '')
                                                        } else {
                                                            names(genes_of_interest.modules.gX.gY) = paste(names(genes_of_interest.input[gX]),'.',module.type[M],genes_of_interest.modules[gY], sep = '')
                                                        }
                                                        module_set.list <- c(module_set.list, genes_of_interest.modules.gX.gY)
                                                    }
                                                    genes_of_interest.modules <- list(unique(genes_of_interest.modules))
                                                    #if (length(unlist(genes_of_interest.modules)) == 1) {
                                                    #    names(genes_of_interest.modules) = paste("query_",module.type[M],"_of_interest", sep = '')
                                                    #} else if (length(unlist(genes_of_interest.modules)) > 1) {
                                                    #    names(genes_of_interest.modules) = paste("query_",module.type[M],"s","_of_interest", sep = '')
                                                    #}
                                                    names(genes_of_interest.modules) = paste(names(genes_of_interest.input[gX]),'.',"ALL", sep = '')
                                                    module_set.list <- c(module_set.list, genes_of_interest.modules)
                                                }
                                            }
                                        }
                                    }
                                    
                                    # *** signature gene(s) ***
                                    if (! is.null(signature_genes.input)) {
                                        for (gX in 1:length(signature_genes.input)) {
                                            signature_genes.genes <- signature_genes.input[[gX]]
                                            signature_genes.genes.idx <- which(sapply(signature_genes.genes, FUN = function(X) X %in% unlist(gene_group_df[,"gene_short_name"])))
                                            if (length(signature_genes.genes.idx) > 0) {
                                                signature_genes.genes = signature_genes.genes[signature_genes.genes.idx]
                                                
                                                signature_genes.genes.idx <- lapply(signature_genes.genes, FUN = function(X) which(unlist(gene_group_df[,"gene_short_name"]) %in% X))
                                                names(signature_genes.genes.idx) = signature_genes.genes
                                                
                                                marker.genes = gene_group_df[unlist(signature_genes.genes.idx),]
                                                
                                                marker.modules <- marker.genes[,module.type[M]]
                                                signature_genes.modules <- as.character(unlist(unique(marker.modules)))
                                                signature_genes.modules <- signature_genes.modules[which(! signature_genes.modules %in% signature_genes.input.title[which(! signature_genes.input.title %in% names(signature_genes.input[gX]))])]
                                                if (! is.null(signature_genes.modules)) {
                                                    a = grep(paste(signature_genes.titles, collapse = '|'), signature_genes.modules, value = TRUE, invert = TRUE)
                                                    a.idx = grep(paste(signature_genes.titles, collapse = '|'), signature_genes.modules, value = FALSE, invert = TRUE)
                                                    if (length(grep(paste(signature_genes.titles, collapse = '|'), signature_genes.modules, value = TRUE, invert = FALSE)) > 0) {
                                                        b.signature_genes = grep(paste(signature_genes.titles, collapse = '|'), signature_genes.modules, value = TRUE, invert = FALSE)
                                                        b.signature_genes.idx = grep(paste(signature_genes.titles, collapse = '|'), signature_genes.modules, value = FALSE, invert = FALSE)
                                                    } else {
                                                        b.signature_genes = NULL
                                                        b.signature_genes.idx = NULL
                                                    }
                                                    signature_genes.modules = c(a[mixedorder(a, decreasing = FALSE)], 
                                                                                names(which(sapply(signature_genes.titles, FUN = function(X) X %in% b.signature_genes))))
                                                    
                                                    for (gY in 1:length(signature_genes.modules)) {
                                                        signature_genes.modules.gX.gY <- list(signature_genes.modules[gY])
                                                        if (signature_genes.modules[gY] %in% signature_genes.titles) {
                                                            names(signature_genes.modules.gX.gY) = paste(names(signature_genes.input[gX]),'.',signature_genes.modules[gY], sep = '')
                                                        } else {
                                                            names(signature_genes.modules.gX.gY) = paste(names(signature_genes.input[gX]),'.',module.type[M],signature_genes.modules[gY], sep = '')
                                                        }
                                                        module_set.list <- c(module_set.list, signature_genes.modules.gX.gY)
                                                    }
                                                    signature_genes.modules <- list(unique(signature_genes.modules))
                                                    #if (length(unlist(signature_genes.modules)) == 1) {
                                                    #    names(signature_genes.modules) = paste("signature_",module.type[M],"_of_interest", sep = '')
                                                    #} else if (length(unlist(signature_genes.modules)) > 1) {
                                                    #    names(signature_genes.modules) = paste("signature_",module.type[M],"s","_of_interest", sep = '')
                                                    #}
                                                    names(signature_genes.modules) = paste(names(signature_genes.input[gX]),'.',"ALL", sep = '')
                                                    module_set.list <- c(module_set.list, signature_genes.modules)
                                                }
                                            }
                                        }
                                    }
                                    
                                    if (module.type[M] == "module") {
                                        example_modules <- list(unique(c(27, 10, 7, 30)))
                                        if (length(unlist(example_modules)) == 1) {
                                            names(example_modules) = "example_module"
                                        } else if (length(unlist(example_modules)) > 1) {
                                            names(example_modules) = "example_modules"
                                        }
                                        module_set.list <- c(module_set.list, example_modules)
                                    }
                                    
                                    if (length(module_set.list) > 0) {
                                        for (N in 1:length(module_set.list)) {
                                            #print(names(module_set.list[N]))
                                            
                                            marker.modules <- module_set.list[[N]]
                                            marker.modules.idx <- which(sapply(marker.modules, FUN = function(X) X %in% unlist(gene_group_df[,module.type[M]])))
                                            if (length(marker.modules.idx) > 0) {
                                                marker.modules = marker.modules[marker.modules.idx]
                                                
                                                marker.modules.idx <- lapply(marker.modules, FUN = function(X) which(unlist(gene_group_df[,module.type[M]]) %in% X))
                                                names(marker.modules.idx) = marker.modules
                                                
                                                marker.genes = gene_group_df[unlist(marker.modules.idx),]
                                                
                                                gene_id.dup.idx <- which(duplicated(marker.genes[,"id"]))
                                                if (length(gene_id.dup.idx) > 0) {
                                                    rm.idx.all <- c()
                                                    for (d in 1:length(gene_id.dup.idx)) {
                                                        gene_id.idx <- which(marker.genes[,"id"] == unlist(marker.genes[gene_id.dup.idx[d],"id"]))
                                                        
                                                        gene_id <- marker.genes[gene_id.idx,"id"]
                                                        gene_module <- marker.genes[gene_id.idx,"module"]
                                                        #gene_supermodule <- marker.genes[gene_id.idx,"supermodule"]
                                                        #gene_name <- marker.genes[gene_id.idx,"gene_short_name"]
                                                        
                                                        gene_id <- unlist(gene_id)
                                                        gene_module <- unlist(gene_module)
                                                        #gene_supermodule <- unlist(gene_supermodule)
                                                        #gene_name <- unlist(gene_name)
                                                        
                                                        gene_module.idx <- which(gene_module %in% c(genes_of_interest.titles, signature_genes.titles))
                                                        #gene_supermodule.idx <- which(gene_supermodule %in% c(genes_of_interest.titles, signature_genes.titles))
                                                        
                                                        rm.idx <- gene_id.idx[gene_module.idx]
                                                        
                                                        rm.idx.all <- c(rm.idx.all, rm.idx)
                                                    }
                                                    rm.idx.all <- unique(rm.idx.all)
                                                    #print(rm.idx.all)
                                                    
                                                    if (length(rm.idx.all) > 0) {
                                                        marker.genes <- marker.genes[-rm.idx.all, , drop = FALSE]
                                                    }
                                                }
                                                
                                                if (module.type[M] == "module") {
                                                    marker.modules <- marker.genes[,"module"]
                                                } else if (module.type[M] == "supermodule") {
                                                    marker.modules <- marker.genes[,"module"]
                                                }
                                                marker.modules.elements <- as.character(unlist(marker.modules))
                                                if (any(module_set.list[[N]] %in% c("genes_of_interest", "signature_genes"))) {
                                                    marker.modules.elements <- marker.modules.elements
                                                } else {
                                                    marker.modules.elements <- marker.modules.elements[which(! marker.modules.elements %in% c(genes_of_interest.input.title, signature_genes.input.title)[which(! c(genes_of_interest.input.title, signature_genes.input.title) %in% unlist(strsplit(names(module_set.list[N]), '\\.'))[1])])]
                                                }
                                                if (! is.null(marker.modules.elements)) {
                                                    if (! is.null(genes_of_interest.input)) {
                                                        for (gX in 1:length(genes_of_interest.input)) {
                                                            marker.modules.elements[which(marker.modules.elements == paste("geneSet",gX,"_of_interest", sep = ''))] = paste("geneSet",gX, sep = '')
                                                        }
                                                    }
                                                    marker.modules.elements[which(marker.modules.elements == "SINGH_KRAS_DEPENDENCY_SIGNATURE")] = "SINGH"
                                                    marker.modules.elements[which(marker.modules.elements == "FRIDMAN_SENESCENCE_UP")] = "FRIDMAN"
                                                    marker.modules.elements[which(marker.modules.elements == "HALLMARK_KRAS_SIGNALING_UP")] = "HALLMARK"
                                                    
                                                    if (module.type[M] == "module") {
                                                        marker.genes[,"module"] = factor(marker.modules.elements, levels = unique(marker.modules.elements))
                                                    } else if (module.type[M] == "supermodule") {
                                                        marker.genes[,"module"] = factor(marker.modules.elements, levels = unique(marker.modules.elements))
                                                    }
                                                }
                                                
                                                if (module.type[M] == "module") {
                                                    panel.num = length(unlist(unique(marker.genes[,"module"])))
                                                } else if (module.type[M] == "supermodule") {
                                                    panel.num = length(unlist(unique(marker.genes[,"module"])))
                                                }
                                                if (panel.num == 1) {
                                                    width.size = 7
                                                    height.size = 6
                                                } else if (panel.num == 2) {
                                                    width.size = 7 + (5 * 1)
                                                    height.size = 6
                                                } else if (panel.num == 3) {
                                                    width.size = 7 + (5 * 2)
                                                    height.size = 6
                                                } else if (panel.num > 3) {
                                                    width.size = 7 + (5 * 3)
                                                    height.size = width.size - 1
                                                }
                                                
                                                #for (redDim in 1:length(redDim.method)) {
                                                
                                                # /// 2d ///
                                                dim.folder = paste(eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))),'D', sep = '')
                                                dim.finalX = eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) - (2-1)
                                                
                                                if (dim.finalX >= 1) {
                                                    
                                                    for (dim.comb in 1:dim.finalX) {
                                                        dim.x <- dim.comb
                                                        dim.y <- dim.x + 1
                                                        
                                                        if (dim.finalX > 1) {
                                                            if (dim.x == 1) {
                                                                comb.run = 2
                                                            } else {
                                                                comb.run = 1
                                                            }
                                                        } else {
                                                            comb.run = 1
                                                        }
                                                        
                                                        for (dim.comb.run in 1:comb.run) {
                                                            if (dim.comb.run > 1) {
                                                                dim.y <- dim.y + 1
                                                            }
                                                            
                                                            if (eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) == 2) {
                                                                dim.xy <- paste('D',dim.x,'-','D',dim.y, sep = '')
                                                            } else {
                                                                dim.xy <- paste('D',dim.x,'-','D',dim.y, sep = '')
                                                            }
                                                            #print(dim.xy)
                                                            
                                                            if (SUB == 0) {
                                                                file.prefolder = dim.xy
                                                            } else {
                                                                file.prefolder = paste(dim.xy,'/',SUB.prefolder, sep = '')
                                                            }
                                                            
                                                            mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder, sep = '')
                                                            subDir <- 'differential'
                                                            if (! file.exists(file.path(mainDir, subDir))) {
                                                                dir.create(file.path(mainDir, subDir))
                                                            }
                                                            mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential', sep = '')
                                                            subDir <- module.type[M]
                                                            if (! file.exists(file.path(mainDir, subDir))) {
                                                                dir.create(file.path(mainDir, subDir))
                                                            }
                                                            
                                                            if (length(grep(paste(genes_of_interest.titles, collapse = '|'), names(module_set.list[N]), value = TRUE, invert = FALSE)) > 0) {
                                                                mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential', sep = '')
                                                                subDir <- 'genes_of_interest'
                                                                if (! file.exists(file.path(mainDir, subDir))) {
                                                                    dir.create(file.path(mainDir, subDir))
                                                                }
                                                                mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential','/','genes_of_interest', sep = '')
                                                                subDir <- unlist(strsplit(names(module_set.list[N]), '\\.'))[1]
                                                                if (! file.exists(file.path(mainDir, subDir))) {
                                                                    dir.create(file.path(mainDir, subDir))
                                                                }
                                                                mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential','/','genes_of_interest','/',unlist(strsplit(names(module_set.list[N]), '\\.'))[1], sep = '')
                                                                subDir <- module.type[M]
                                                                if (! file.exists(file.path(mainDir, subDir))) {
                                                                    dir.create(file.path(mainDir, subDir))
                                                                }
                                                                
                                                                file.folder = paste('genes_of_interest','/',unlist(strsplit(names(module_set.list[N]), '\\.'))[1],'/',module.type[M], sep = '')
                                                                if (names(module_set.list[N]) %in% genes_of_interest.titles) {
                                                                    file.name.suffix = names(module_set.list[N])
                                                                } else {
                                                                    file.name.suffix = gsub("\\.","_",unlist(strsplit(names(module_set.list[N]), '\\.'))[2])
                                                                }
                                                            } else if (length(grep(paste(signature_genes.titles, collapse = '|'), names(module_set.list[N]), value = TRUE, invert = FALSE)) > 0) {
                                                                mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential', sep = '')
                                                                subDir <- 'signature_genes'
                                                                if (! file.exists(file.path(mainDir, subDir))) {
                                                                    dir.create(file.path(mainDir, subDir))
                                                                }
                                                                mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential','/','signature_genes', sep = '')
                                                                subDir <- unlist(strsplit(names(module_set.list[N]), '\\.'))[1]
                                                                if (! file.exists(file.path(mainDir, subDir))) {
                                                                    dir.create(file.path(mainDir, subDir))
                                                                }
                                                                mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential','/','signature_genes','/',unlist(strsplit(names(module_set.list[N]), '\\.'))[1], sep = '')
                                                                subDir <- module.type[M]
                                                                if (! file.exists(file.path(mainDir, subDir))) {
                                                                    dir.create(file.path(mainDir, subDir))
                                                                }
                                                                
                                                                file.folder = paste('signature_genes','/',unlist(strsplit(names(module_set.list[N]), '\\.'))[1],'/',module.type[M], sep = '')
                                                                if (names(module_set.list[N]) %in% signature_genes.titles) {
                                                                    file.name.suffix = names(module_set.list[N])
                                                                } else {
                                                                    file.name.suffix = gsub("\\.","_",unlist(strsplit(names(module_set.list[N]), '\\.'))[2])
                                                                }
                                                            } else {
                                                                file.folder = module.type[M]
                                                                file.name.suffix = gsub("\\.","_",names(module_set.list[N]))
                                                            }
                                                            
                                                            if (is.null(file.name.suffix)) {
                                                                file.name = paste(project,'.',file.prename,'.','plot_cells','.',redDim.method[redDim],'_',dim.xy,'.','differential', sep = '')
                                                            } else {
                                                                if (is.na(file.name.suffix) | file.name.suffix == "") {
                                                                    file.name = paste(project,'.',file.prename,'.','plot_cells','.',redDim.method[redDim],'_',dim.xy,'.','differential', sep = '')
                                                                } else {
                                                                    file.name = paste(project,'.',file.prename,'.','plot_cells','.',redDim.method[redDim],'_',dim.xy,'.','differential','.',file.name.suffix, sep = '')
                                                                }
                                                            }
                                                            
                                                            png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential','/',file.folder,'/',file.name,'.png', sep = ''), width = width.size, height = height.size, units = 'in', res = 600)
                                                            p <- plot_cells(cds, reduction_method = redDim.method[redDim], x = dim.x, y = dim.y, cell_size = cell_size_2d, 
                                                                            genes = marker.genes, norm_method = "log", 
                                                                            #color_cells_by = "partition", group_cells_by = "partition", 
                                                                            label_cell_groups = FALSE, 
                                                                            show_trajectory_graph = FALSE)
                                                            print(p)
                                                            dev.off()
                                                        }
                                                    }
                                                }
                                                
                                                # /// 3d ///
                                                if (module.type[M] == "module") {
                                                    
                                                    genes_of_interest.input.title.idx <- which(module.levels %in% genes_of_interest.titles)
                                                    signature_genes.input.title.idx <- which(module.levels %in% signature_genes.titles)
                                                    
                                                    if ((N <= 3) | (N %in% c(genes_of_interest.input.title.idx, signature_genes.input.title.idx)) | (N > length(module.levels))) { # keep at most 3 modules (not save more than 3 modules) for saving disk space, please modify or un-comment this condition if needed
                                                        
                                                        if (eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) == 2) {
                                                            dim.folder = paste(eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))),'D', sep = '')
                                                            #dim.finalX = eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) - (3-1)
                                                            #dim.folder = paste(3,'D', sep = '')
                                                            dim.finalX = 3 - (3-1)
                                                        } else {
                                                            dim.folder = paste(eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))),'D', sep = '')
                                                            dim.finalX = eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) - (3-1)
                                                        }
                                                        
                                                        if (dim.finalX >= 1) {
                                                            
                                                            for (dim.comb in 1:dim.finalX) {
                                                                dim.x <- dim.comb
                                                                dim.y <- dim.x + 1
                                                                dim.z <- dim.y + 1
                                                                
                                                                if (eval(parse(text=paste(redDim.method[redDim],'.','dim_max', sep = ''))) == 2) {
                                                                    dim.xyz <- paste('D',dim.x,'-','D',dim.y,'-','D',dim.z, sep = '')
                                                                } else {
                                                                    dim.xyz <- paste('D',dim.x,'-','D',dim.y,'-','D',dim.z, sep = '')
                                                                }
                                                                #print(dim.xyz)
                                                                
                                                                if (SUB == 0) {
                                                                    file.prefolder = dim.xyz
                                                                } else {
                                                                    file.prefolder = paste(dim.xyz,'/',SUB.prefolder, sep = '')
                                                                }
                                                                
                                                                mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder, sep = '')
                                                                subDir <- 'differential'
                                                                if (! file.exists(file.path(mainDir, subDir))) {
                                                                    dir.create(file.path(mainDir, subDir))
                                                                }
                                                                mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential', sep = '')
                                                                subDir <- module.type[M]
                                                                if (! file.exists(file.path(mainDir, subDir))) {
                                                                    dir.create(file.path(mainDir, subDir))
                                                                }
                                                                
                                                                if (length(grep(paste(genes_of_interest.titles, collapse = '|'), names(module_set.list[N]), value = TRUE, invert = FALSE)) > 0) {
                                                                    mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential', sep = '')
                                                                    subDir <- 'genes_of_interest'
                                                                    if (! file.exists(file.path(mainDir, subDir))) {
                                                                        dir.create(file.path(mainDir, subDir))
                                                                    }
                                                                    mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential','/','genes_of_interest', sep = '')
                                                                    subDir <- unlist(strsplit(names(module_set.list[N]), '\\.'))[1]
                                                                    if (! file.exists(file.path(mainDir, subDir))) {
                                                                        dir.create(file.path(mainDir, subDir))
                                                                    }
                                                                    mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential','/','genes_of_interest','/',unlist(strsplit(names(module_set.list[N]), '\\.'))[1], sep = '')
                                                                    subDir <- module.type[M]
                                                                    if (! file.exists(file.path(mainDir, subDir))) {
                                                                        dir.create(file.path(mainDir, subDir))
                                                                    }
                                                                    
                                                                    file.folder = paste('genes_of_interest','/',unlist(strsplit(names(module_set.list[N]), '\\.'))[1],'/',module.type[M], sep = '')
                                                                    if (names(module_set.list[N]) %in% genes_of_interest.titles) {
                                                                        file.name.suffix = names(module_set.list[N])
                                                                    } else {
                                                                        file.name.suffix = gsub("\\.","_",unlist(strsplit(names(module_set.list[N]), '\\.'))[2])
                                                                    }
                                                                } else if (length(grep(paste(signature_genes.titles, collapse = '|'), names(module_set.list[N]), value = TRUE, invert = FALSE)) > 0) {
                                                                    mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential', sep = '')
                                                                    subDir <- 'signature_genes'
                                                                    if (! file.exists(file.path(mainDir, subDir))) {
                                                                        dir.create(file.path(mainDir, subDir))
                                                                    }
                                                                    mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential','/','signature_genes', sep = '')
                                                                    subDir <- unlist(strsplit(names(module_set.list[N]), '\\.'))[1]
                                                                    if (! file.exists(file.path(mainDir, subDir))) {
                                                                        dir.create(file.path(mainDir, subDir))
                                                                    }
                                                                    mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential','/','signature_genes','/',unlist(strsplit(names(module_set.list[N]), '\\.'))[1], sep = '')
                                                                    subDir <- module.type[M]
                                                                    if (! file.exists(file.path(mainDir, subDir))) {
                                                                        dir.create(file.path(mainDir, subDir))
                                                                    }
                                                                    
                                                                    file.folder = paste('signature_genes','/',unlist(strsplit(names(module_set.list[N]), '\\.'))[1],'/',module.type[M], sep = '')
                                                                    if (names(module_set.list[N]) %in% signature_genes.titles) {
                                                                        file.name.suffix = names(module_set.list[N])
                                                                    } else {
                                                                        file.name.suffix = gsub("\\.","_",unlist(strsplit(names(module_set.list[N]), '\\.'))[2])
                                                                    }
                                                                } else {
                                                                    file.folder = module.type[M]
                                                                    file.name.suffix = gsub("\\.","_",names(module_set.list[N]))
                                                                }
                                                                
                                                                if (redDim.method[redDim] == "UMAP") {
                                                                    trj.plt = 1
                                                                } else if (redDim.method[redDim] == "tSNE") {
                                                                    trj.plt = 0
                                                                }
                                                                if (graph.run == 0) {
                                                                    trj.plt = 0
                                                                }
                                                                
                                                                for (trj in 0:trj.plt) {
                                                                    if (trj == 0) {
                                                                        trj.prt.endidx = 1
                                                                    } else if (trj == 1) {
                                                                        trj.prt.endidx = 0
                                                                    }
                                                                    
                                                                    #for (trj.prt in 1:1) {
                                                                    for (trj.prt in 1:trj.prt.endidx) {
                                                                        if (trj.prt == 1) {
                                                                            trj.prt.name = 'with_partition'
                                                                            cds_3d.plt = cds_3d
                                                                        } else if (trj.prt == 0) {
                                                                            trj.prt.name = 'without_partition'
                                                                            cds_3d.plt = cds_3d.no_partition
                                                                        }
                                                                        
                                                                        if (! is.null(cds_3d.plt)) {
                                                                            
                                                                            for (grid.plt in 1:0) {
                                                                                if (grid.plt == 1) {
                                                                                    grid.name = 'with_grid'
                                                                                } else if (grid.plt == 0) {
                                                                                    grid.name = 'without_grid'
                                                                                }
                                                                                
                                                                                if (trj == 0) {
                                                                                    show_trajectory = FALSE
                                                                                    if (is.null(file.name.suffix)) {
                                                                                        file.name = paste(project,'.',file.prename,'.','plot_cells_3d','.',redDim.method[redDim],'_',dim.xyz,'.','differential','.','without_trajectory', sep = '')
                                                                                    } else {
                                                                                        if (is.na(file.name.suffix) | file.name.suffix == "") {
                                                                                            file.name = paste(project,'.',file.prename,'.','plot_cells_3d','.',redDim.method[redDim],'_',dim.xyz,'.','differential','.','without_trajectory', sep = '')
                                                                                        } else {
                                                                                            file.name = paste(project,'.',file.prename,'.','plot_cells_3d','.',redDim.method[redDim],'_',dim.xyz,'.','differential','.',file.name.suffix,'.','without_trajectory', sep = '')
                                                                                        }
                                                                                    }
                                                                                    file.name = paste(file.name,'.',grid.name, sep = '')
                                                                                } else if (trj == 1) {
                                                                                    show_trajectory = TRUE
                                                                                    if (is.null(file.name.suffix)) {
                                                                                        file.name = paste(project,'.',file.prename,'.','plot_cells_3d','.',redDim.method[redDim],'_',dim.xyz,'.','differential','.','with_trajectory', sep = '')
                                                                                    } else {
                                                                                        if (is.na(file.name.suffix) | file.name.suffix == "") {
                                                                                            file.name = paste(project,'.',file.prename,'.','plot_cells_3d','.',redDim.method[redDim],'_',dim.xyz,'.','differential','.','with_trajectory', sep = '')
                                                                                        } else {
                                                                                            file.name = paste(project,'.',file.prename,'.','plot_cells_3d','.',redDim.method[redDim],'_',dim.xyz,'.','differential','.',file.name.suffix,'.','with_trajectory', sep = '')
                                                                                        }
                                                                                    }
                                                                                    file.name = paste(file.name,'.',trj.prt.name, sep = '')
                                                                                    file.name = paste(file.name,'.',grid.name, sep = '')
                                                                                }
                                                                                
                                                                                file.html = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential','/',file.folder,'/',file.name,'.html', sep = '')
                                                                                file.png = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential','/',file.folder,'/',file.name,'.png', sep = '')
                                                                                
                                                                                cds_3d_plot_obj <- plot_cells_3d(cds_3d.plt, reduction_method = redDim.method[redDim], dims = c(dim.x, dim.y, dim.z), cell_size = cell_size_3d, genes = marker.genes, norm_method = "log", show_trajectory_graph = show_trajectory)
                                                                                
                                                                                if (grid.plt == 0) {
                                                                                    cds_3d_plot_obj <- cds_3d_plot_obj %>% plotly::layout(scene = scene)
                                                                                }
                                                                                
                                                                                # https://community.rstudio.com/t/save-viewer-object-rendered-in-rstudio-as-image/32796
                                                                                if (dim.comb == 1) { # only save D1-D2-D3 (no D2-D3-D4, D3-D4-D5, ...)
                                                                                    if (length(unlist(strsplit(names(module_set.list[N]), '\\.'))) == 1) { # keep part of modules (not save single module) for saving disk space, please modify or un-comment this condition if needed
                                                                                        saveWidget(cds_3d_plot_obj, file.html)
                                                                                        #webshot(file.html, file.png)
                                                                                    } else {
                                                                                        # *** query gene(s) ***
                                                                                        if (unlist(strsplit(names(module_set.list[N]), '\\.'))[1] %in% genes_of_interest.titles) {
                                                                                            saveWidget(cds_3d_plot_obj, file.html)
                                                                                            #webshot(file.html, file.png)
                                                                                        }
                                                                                        if (unlist(strsplit(names(module_set.list[N]), '\\.'))[2] %in% c(genes_of_interest.titles, "ALL")) {
                                                                                            saveWidget(cds_3d_plot_obj, file.html)
                                                                                            #webshot(file.html, file.png)
                                                                                        }
                                                                                        # *** signature gene(s) ***
                                                                                        if (unlist(strsplit(names(module_set.list[N]), '\\.'))[2] %in% c(signature_genes.titles, "ALL")) {
                                                                                            saveWidget(cds_3d_plot_obj, file.html)
                                                                                            #webshot(file.html, file.png)
                                                                                        }
                                                                                    }
                                                                                } else { # remove D2-D3-D4, D3-D4-D5, ...
                                                                                    mainDir <- paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/',redDim.method[redDim],'/',dim.folder,'/',file.prefolder,'/','differential', sep = '')
                                                                                    if (file.exists(file.path(mainDir))) {
                                                                                        unlink(mainDir, recursive = TRUE)
                                                                                    }
                                                                                }
                                                                            }
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                                #}
                                            }
                                        }
                                    }
                                }
                            }
                            
                            
                            # Monocle offers another plotting function that can sometimes give a clearer view of a gene's dynamics along a single path.
                            # You can select a path with choose_cells() or by subsetting the cell data set by cluster, cell type, or other annotation that's restricted to the path. 
                            # Let's pick one such path, the AFD cells:
                            # /// marker.genes ///
                            if (project.name == "tutorial") {
                                cells <- c("AFD")
                            } else {
                                cells <- cells.input
                            }
                            
                            if (length(cells) > 0) {
                                
                                gene_set.list <- list()
                                
                                for (N in 1:length(cells)) {
                                    if (cells[N] == "AFD") {
                                        cell_genes <- list(unique(c("gcy-8", "dac-1", "oig-8"))) # You can see that dac-1 is activated before the other two genes.
                                    } else {
                                        cell_genes <- list(unique(c("")))
                                    }
                                    
                                    if (length(unlist(cell_genes)) == 1) {
                                        names(cell_genes) = paste(cells[N],'_','gene', sep = '')
                                        assign(paste(cells[N],'_','gene', sep = ''), cell_genes)
                                    } else if (length(unlist(cell_genes)) > 1) {
                                        names(cell_genes) = paste(cells[N],'_','genes', sep = '')
                                        assign(paste(cells[N],'_','genes', sep = ''), cell_genes)
                                    }
                                    
                                    gene_set.list <- c(gene_set.list, cell_genes)
                                }
                                
                                if (length(gene_set.list) > 0) {
                                    for (N in 1:length(gene_set.list)) {
                                        #print(names(gene_set.list[N]))
                                        
                                        marker.genes <- gene_set.list[[N]]
                                        marker.genes.idx <- which(sapply(marker.genes, FUN = function(X) X %in% rowData(cds)[,"gene_short_name"]))
                                        if (length(marker.genes.idx) > 0) {
                                            marker.genes = marker.genes[marker.genes.idx]
                                            
                                            rowidx <- sapply(marker.genes, FUN = function(X) which(rowData(cds)[,"gene_short_name"] %in% X))
                                            rowidx <- unlist(rowidx)
                                            
                                            colidx <- lapply(cells[N], FUN = function(X) which(colData(cds)$cell.type %in% X))
                                            names(colidx) = cells[N]
                                            colidx <- unlist(colidx)
                                            
                                            cds_subset <- cds[rowidx, colidx, drop = FALSE]
                                            
                                            #rowData(cds_subset)[,"gene_short_name"] = factor(rowData(cds_subset)[,"gene_short_name"], levels = unique(rowData(cds_subset)[,"gene_short_name"]))
                                            
                                            # The function plot_genes_in_pseudotime() takes a small set of genes and shows you their dynamics as a function of pseudotime:
                                            if (project.name == "tutorial") {
                                                time = "embryo.time.bin"
                                            } else {
                                                time = time.input
                                            }
                                            
                                            if (nrow(cds_subset) <= 3) {
                                                ncol.num = 5/5
                                                fontsize = 12
                                                height.size = 3
                                            } else {
                                                if (nrow(cds_subset) > 3 & nrow(cds_subset) <= 5) {
                                                    ncol.num = 5/5
                                                    fontsize = 12
                                                } else if (nrow(cds_subset) > 5 & nrow(cds_subset) <= 10) {
                                                    ncol.num = 10/5
                                                    fontsize = 12 + (4 * 1)
                                                } else if (nrow(cds_subset) > 10 & nrow(cds_subset) <= 15) {
                                                    ncol.num = 15/5
                                                    fontsize = 12 + (4 * 2)
                                                } else if (nrow(cds_subset) > 15 & nrow(cds_subset) <= 20) {
                                                    ncol.num = 20/5
                                                    fontsize = 12 + (4 * 3)
                                                } else if (nrow(cds_subset) > 20) {
                                                    ncol.num = 25/5
                                                    fontsize = 12 + (4 * 4)
                                                }
                                                height.size = nrow(cds_subset) * 0.9
                                            }
                                            width.size = ((5 * ncol.num) + 0.4625) * 0.9
                                            
                                            png(file = paste(current_folder,'/Result','/','Monocle3','/',project.name,'/','trajectories','/',subdata.folder,'/',project.subname,'/','UMAP','/',dim.folder,'/',differential.prefolder,'/',project,'.',gsub("\\.","_",subdata.folder),'.',gsub("\\.","_",project.subname),'.','UMAP','.',file.prename,'.','differential','.','plot_genes_in_pseudotime','.',gsub("\\.","_",time),'.',gsub("\\.","_",names(gene_set.list[N])),'.png', sep = ''), width = width.size, height = height.size, units = 'in', res = 600)
                                            p <- plot_genes_in_pseudotime(cds_subset, color_cells_by = time, min_expr = 0.5, ncol = ncol.num)
                                            p <- p + theme(text = element_text(size = fontsize))
                                            print(p)
                                            dev.off()
                                        }
                                    }
                                }
                            }
                            
                            
                            # --- Analyzing branches in single-cell trajectories ----------------------------------------------
                            
                            
                            
                            
                            
                            
                            
                            cds = cds.origin.processed
                            cds_3d = cds_3d.origin.processed
                            cds_3d.no_partition = cds_3d.no_partition.origin.processed
                        }
                    }
                }
                
                cds.trajectories <- cds.origin.processed
                #get_citations(cds.trajectories)
                
                cds_3d.trajectories <- cds_3d.origin.processed
                #get_citations(cds_3d.trajectories)
                
                cds_3d.no_partition.trajectories <- cds_3d.no_partition.origin.processed
                #get_citations(cds_3d.no_partition.trajectories)
                
            }
        }
    }
    
    options(ggrepel.max.overlaps = NULL)
    #options(op) # reset (all) initial options
    
}
