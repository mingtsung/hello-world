rm(list=ls())	# remove all variables in workspace

setwd("D:/PD/Analysis/miRNA/")

mainDir <- '.'
subDir <- 'Data'
if(file.exists(file.path(mainDir, subDir))) {
} else {
	dir.create(file.path(mainDir, subDir))
}
mainDir <- './Data'
subDir <- '02_significant'
if(file.exists(file.path(mainDir, subDir))) {
} else {
	dir.create(file.path(mainDir, subDir))
}
mainDir <- './Data/02_significant'
subDir <- 'piRNA'
if(file.exists(file.path(mainDir, subDir))) {
} else {
	dir.create(file.path(mainDir, subDir))
}

#install.packages("openxlsx")	# openxlsx: Read, Write and Edit XLSX Files
library(openxlsx)

# --- GroupInfo -----------------------------------------------------------------------------------
diseaseStage.ALL <- read.xlsx("./Data/01_datatable/piRNA_datatable.xlsx", sheet="sample information", colNames=TRUE, rowNames=FALSE, check.names=FALSE)
colnames(diseaseStage.ALL) = gsub('\\.', ' ', colnames(diseaseStage.ALL))
colnames(diseaseStage.ALL)[which(colnames(diseaseStage.ALL)=='p s')] = 'p.s'
colnames(diseaseStage.ALL)[which(colnames(diseaseStage.ALL)=='No ')] = 'No.'

# HC: Healthy control; PD: Parkinson's disease (ND: Non-dementia, MCI: Mild Cognitive Impairment, D: Dementia); MSA: Multiple System Atrophy
diseaseStage <- c('HC','PDND','PD-MCI','PDD','MSA')
diseaseGroup <- c('HC','PDND','PDMCI','PDD','PD','MSA')
diseaseStage.HC.idx = which(diseaseStage.ALL[,'Symp Group']=='HC')
diseaseStage.PDND.idx = which(diseaseStage.ALL[,'Symp Group']=='PDND')
diseaseStage.PDMCI.idx = which(diseaseStage.ALL[,'Symp Group']=='PD-MCI')
diseaseStage.PDD.idx = which(diseaseStage.ALL[,'Symp Group']=='PDD')
diseaseStage.PD.idx = c(diseaseStage.PDND.idx, diseaseStage.PDMCI.idx, diseaseStage.PDD.idx)
diseaseStage.MSA.idx = which(diseaseStage.ALL[,'Symp Group']=='MSA')

diseaseStage.HC = diseaseStage.ALL[diseaseStage.HC.idx, , drop=FALSE]
diseaseStage.PDND = diseaseStage.ALL[diseaseStage.PDND.idx, , drop=FALSE]
diseaseStage.PDMCI = diseaseStage.ALL[diseaseStage.PDMCI.idx, , drop=FALSE]
diseaseStage.PDD = diseaseStage.ALL[diseaseStage.PDD.idx, , drop=FALSE]
diseaseStage.PD = diseaseStage.ALL[diseaseStage.PD.idx, , drop=FALSE]
diseaseStage.MSA = diseaseStage.ALL[diseaseStage.MSA.idx, , drop=FALSE]

#diseaseStage.color <- c('limegreen','orangered1','orangered3','orangered4','gray40')
diseaseStage.color <- c('darkgreen','indianred','red','darkred','black')
names(diseaseStage.color) = diseaseStage

# --- piRNA datatable -----------------------------------------------------------------------------
piRNA.READsData <- read.xlsx("./Data/01_datatable/piRNA_datatable.xlsx", sheet="piRNA(READs)", colNames=TRUE, rowNames=FALSE, check.names=FALSE)
piRNA.UMIsData <- read.xlsx("./Data/01_datatable/piRNA_datatable.xlsx", sheet="piRNA(UMIs)", colNames=TRUE, rowNames=FALSE, check.names=FALSE)
piRNA.rawData <- read.xlsx("./Data/01_datatable/piRNA_datatable.xlsx", sheet="piRNA(raw)", colNames=TRUE, rowNames=FALSE, check.names=FALSE)
piRNA.normalizedData <- read.xlsx("./Data/01_datatable/piRNA_datatable.xlsx", sheet="piRNA(normalized)", colNames=TRUE, rowNames=FALSE, check.names=FALSE)
piRNA.normalizedData.HC <- read.xlsx("./Data/01_datatable/piRNA_datatable.xlsx", sheet="piRNA(normalized)(HC)", colNames=TRUE, rowNames=FALSE, check.names=FALSE)
piRNA.normalizedData.PDND <- read.xlsx("./Data/01_datatable/piRNA_datatable.xlsx", sheet="piRNA(normalized)(PDND)", colNames=TRUE, rowNames=FALSE, check.names=FALSE)
piRNA.normalizedData.PDMCI <- read.xlsx("./Data/01_datatable/piRNA_datatable.xlsx", sheet="piRNA(normalized)(PD-MCI)", colNames=TRUE, rowNames=FALSE, check.names=FALSE)
piRNA.normalizedData.PDD <- read.xlsx("./Data/01_datatable/piRNA_datatable.xlsx", sheet="piRNA(normalized)(PDD)", colNames=TRUE, rowNames=FALSE, check.names=FALSE)
piRNA.normalizedData.PD <- read.xlsx("./Data/01_datatable/piRNA_datatable.xlsx", sheet="piRNA(normalized)(PD)", colNames=TRUE, rowNames=FALSE, check.names=FALSE)
piRNA.normalizedData.MSA <- read.xlsx("./Data/01_datatable/piRNA_datatable.xlsx", sheet="piRNA(normalized)(MSA)", colNames=TRUE, rowNames=FALSE, check.names=FALSE)

# comparison
combn.idx <- combn(length(diseaseGroup), 2)
for(m in 1:ncol(combn.idx))
{
	gc()	# Garbage Collection (for release memory)
	
	Group <- c(diseaseGroup[combn.idx[1,m]], diseaseGroup[combn.idx[2,m]])
	combnName <- paste(diseaseGroup[combn.idx[1,m]], diseaseGroup[combn.idx[2,m]], sep='_')
	
	mainDir <- './Data/02_significant/piRNA'
	subDir <- paste(Group[1],' vs ',Group[2],sep='')
	if(file.exists(file.path(mainDir, subDir))) {
	} else {
		dir.create(file.path(mainDir, subDir))
	}
	
	diseaseStage.Group1 <- eval(parse(text=paste('diseaseStage.', Group[1], sep='')))
	diseaseStage.Group2 <- eval(parse(text=paste('diseaseStage.', Group[2], sep='')))
	diseaseStage.Group1_Group2 <- unique(rbind(diseaseStage.Group1, diseaseStage.Group2))
	
	piRNA.normalizedData.Group1 <- eval(parse(text=paste('piRNA.normalizedData.', Group[1], sep='')))
	piRNA.normalizedData.Group2 <- eval(parse(text=paste('piRNA.normalizedData.', Group[2], sep='')))
	
	# piRNA.normalizedData.HC_PDND, piRNA.normalizedData.HC_PDMCI, piRNA.normalizedData.HC_PDD, piRNA.normalizedData.HC_PD, piRNA.normalizedData.HC_MSA
	# piRNA.normalizedData.PDND_PDMCI, piRNA.normalizedData.PDND_PDD, piRNA.normalizedData.PDND_PD, piRNA.normalizedData.PDND_MSA
	# piRNA.normalizedData.PDMCI_PDD, piRNA.normalizedData.PDMCI_PD, piRNA.normalizedData.PDMCI_MSA
	# piRNA.normalizedData.PDD_PD, piRNA.normalizedData.PDD_MSA
	# piRNA.normalizedData.PD_MSA
	assign(paste('piRNA.normalizedData.', combnName, sep=''), merge(piRNA.normalizedData.Group1, piRNA.normalizedData.Group2, sort=FALSE))
	piRNA.normalizedData.Comp <- eval(parse(text=paste('piRNA.normalizedData.', combnName, sep='')))
	rownames(piRNA.normalizedData.Comp) = unlist(sapply(strsplit(piRNA.normalizedData.Comp[,'Name'], '/'), FUN=function(X) X[1]))
	
	# col.idx
	piRNA.annotation = c('type','Name','Accession_ID')
	piRNA.annotation.idx = c(
		which(colnames(piRNA.normalizedData.Comp)==piRNA.annotation[1]), 
		which(colnames(piRNA.normalizedData.Comp)==piRNA.annotation[2]), 
		which(colnames(piRNA.normalizedData.Comp)==piRNA.annotation[3]))
	
	samples_serial_num = sapply(strsplit(colnames(piRNA.normalizedData.Comp)[(length(piRNA.annotation)+1):ncol(piRNA.normalizedData.Comp)], '_'), FUN=function(X) X[1])
	samples_id = sapply(strsplit(colnames(piRNA.normalizedData.Comp)[(length(piRNA.annotation)+1):ncol(piRNA.normalizedData.Comp)], '_'), FUN=function(X) X[2])
	
	Group1.colidx = unlist(sapply(diseaseStage.Group1[,'Sample ID'], FUN=function(X) which(samples_id %in% X) + length(piRNA.annotation)))
	Group2.colidx = unlist(sapply(diseaseStage.Group2[,'Sample ID'], FUN=function(X) which(samples_id %in% X) + length(piRNA.annotation)))
	
	if(Group[1] == 'PD' | Group[2] == 'PD') {
		if(Group[1] == 'PD') {
			PDND.colidx = unlist(sapply(diseaseStage.Group1[which(diseaseStage.Group1[,'Symp Group']=='PDND'),'Sample ID'], FUN=function(X) which(samples_id %in% X) + length(piRNA.annotation)))
			PDMCI.colidx = unlist(sapply(diseaseStage.Group1[which(diseaseStage.Group1[,'Symp Group']=='PD-MCI'),'Sample ID'], FUN=function(X) which(samples_id %in% X) + length(piRNA.annotation)))
			PDD.colidx = unlist(sapply(diseaseStage.Group1[which(diseaseStage.Group1[,'Symp Group']=='PDD'),'Sample ID'], FUN=function(X) which(samples_id %in% X) + length(piRNA.annotation)))
		}
		if(Group[2] == 'PD') {
			PDND.colidx = unlist(sapply(diseaseStage.Group2[which(diseaseStage.Group2[,'Symp Group']=='PDND'),'Sample ID'], FUN=function(X) which(samples_id %in% X) + length(piRNA.annotation)))
			PDMCI.colidx = unlist(sapply(diseaseStage.Group2[which(diseaseStage.Group2[,'Symp Group']=='PD-MCI'),'Sample ID'], FUN=function(X) which(samples_id %in% X) + length(piRNA.annotation)))
			PDD.colidx = unlist(sapply(diseaseStage.Group2[which(diseaseStage.Group2[,'Symp Group']=='PDD'),'Sample ID'], FUN=function(X) which(samples_id %in% X) + length(piRNA.annotation)))
		}
	} else {
		PDND.colidx = NULL
		PDMCI.colidx = NULL
		PDD.colidx = NULL
	}
	
	# --- Boxplot -------------------------------------------------------------------------------------
	gc()	# Garbage Collection (for release memory)
	
	mainDir <- paste('./Data/02_significant/piRNA/',Group[1],' vs ',Group[2],sep='')
	subDir <- 'Boxplot'
	if(file.exists(file.path(mainDir, subDir))) {
	} else {
		dir.create(file.path(mainDir, subDir))
	}
	
	#install.packages("ggplot2")		# ggplot2: Create Elegant Data Visualisations Using the Grammar of Graphics
	suppressPackageStartupMessages(library(ggplot2))
	
	# subset
	Group1.piRNA.rowidx <- c()
	Group2.piRNA.rowidx <- c()
	Group1_Group2.piRNA.rowidx <- c()
	NA.piRNA.rowidx <- c()
	common_piRNA.ttest.pvalue <- c()
	#common_piRNA.ANOVA.pvalue <- c()
	common_piRNA.log2_FC <- c()
	common_piRNA.FC <- c()
	for(i in 1:nrow(piRNA.normalizedData.Comp))
	{
		Group1.value = as.numeric(piRNA.normalizedData.Comp[i, Group1.colidx])
		Group2.value = as.numeric(piRNA.normalizedData.Comp[i, Group2.colidx])
		Group1.value[which(is.na(Group1.value))] = 0
		Group2.value[which(is.na(Group2.value))] = 0
		
		# All idx
		#Group1.value.valididx = c(1:length(Group1.value))
		#Group2.value.valididx = c(1:length(Group2.value))
		# Non-zero idx
		Group1.value.valididx = which(Group1.value!=0)
		Group2.value.valididx = which(Group2.value!=0)
		
		main_title = unlist(strsplit(piRNA.normalizedData.Comp[i,'Name'], '/'))[1]
		
		sampleValue = piRNA.normalizedData.Comp[i, , drop=FALSE]
		sampleGroup = c(rep(NA, length(piRNA.annotation)), rep(Group[1], length(Group1.colidx)), rep(Group[2], length(Group2.colidx)))
		
		sampleValueGroup = as.data.frame(t(rbind(sampleValue, sampleGroup)[, -piRNA.annotation.idx, drop=FALSE]), stringsAsFactors=FALSE)
		colnames(sampleValueGroup) = c('Value','Group')
		sampleValueGroup[,'Value'] = as.numeric(sampleValueGroup[,'Value'])
		sampleValueGroup[,'Group'] = as.factor(sampleValueGroup[,'Group'])
		sampleValueGroup[,'Group'] = factor(sampleValueGroup[,'Group'], levels(sampleValueGroup[,'Group'])[c(which(levels(sampleValueGroup[,'Group'])==Group[1]), 
																											 which(levels(sampleValueGroup[,'Group'])==Group[2]))])
		
		ggplot.plotidx = 0
		if(sum(Group1.value, na.rm=TRUE)!=0 & sum(Group2.value, na.rm=TRUE)==0){
			Group1.piRNA.rowidx <- c(Group1.piRNA.rowidx, i)
			
			mainDir <- paste('./Data/02_significant/piRNA/',Group[1],' vs ',Group[2],'/','Boxplot',sep='')
			subDir <- Group[1]
			if(file.exists(file.path(mainDir, subDir))) {
			} else {
				dir.create(file.path(mainDir, subDir))
			}
			
			# /// boxplot {graphics} ///
			#png(paste('./Data/02_significant/piRNA/',Group[1],' vs ',Group[2],'/','Boxplot','/',Group[1],'/',gsub('/', '.', piRNA.normalizedData.Comp[i,'Name']),'.png',sep=''), width=4, height=4, units="in", res=600)
			#boxplot(Value~Group, data=sampleValueGroup, main=main_title, xlab="Symp Group", ylab="piRNA(normalized)")
			#dev.off()
			
			# /// ggplot {ggplot2} ///
			png(paste('./Data/02_significant/piRNA/',Group[1],' vs ',Group[2],'/','Boxplot','/',Group[1],'/',gsub('/', '.', piRNA.normalizedData.Comp[i,'Name']),'.png',sep=''), width=4, height=4, units="in", res=600)
			ggplot.plotidx = 1
		} else if(sum(Group1.value, na.rm=TRUE)==0 & sum(Group2.value, na.rm=TRUE)!=0){
			Group2.piRNA.rowidx <- c(Group2.piRNA.rowidx, i)
			
			mainDir <- paste('./Data/02_significant/piRNA/',Group[1],' vs ',Group[2],'/','Boxplot',sep='')
			subDir <- Group[2]
			if(file.exists(file.path(mainDir, subDir))) {
			} else {
				dir.create(file.path(mainDir, subDir))
			}
			
			# /// boxplot {graphics} ///
			#png(paste('./Data/02_significant/piRNA/',Group[1],' vs ',Group[2],'/','Boxplot','/',Group[2],'/',gsub('/', '.', piRNA.normalizedData.Comp[i,'Name']),'.png',sep=''), width=4, height=4, units="in", res=600)
			#boxplot(Value~Group, data=sampleValueGroup, main=main_title, xlab="Symp Group", ylab="piRNA(normalized)")
			#dev.off()
			
			# /// ggplot {ggplot2} ///
			png(paste('./Data/02_significant/piRNA/',Group[1],' vs ',Group[2],'/','Boxplot','/',Group[2],'/',gsub('/', '.', piRNA.normalizedData.Comp[i,'Name']),'.png',sep=''), width=4, height=4, units="in", res=600)
			ggplot.plotidx = 1
		} else if(sum(Group1.value, na.rm=TRUE)!=0 & sum(Group2.value, na.rm=TRUE)!=0){
			Group1_Group2.piRNA.rowidx <- c(Group1_Group2.piRNA.rowidx, i)
			
			# t-test
			piRNA.ttest = t.test(x=Group1.value, y=Group2.value, alternative="two.sided")
			piRNA.ttest.pvalue = piRNA.ttest$p.value
			piRNA.ttest$method		# "Welch Two Sample t-test"
			common_piRNA.ttest.pvalue <- c(common_piRNA.ttest.pvalue, piRNA.ttest.pvalue)
			
			# one-way ANOVA test
			#piRNA.ANOVA = aov(formula=Value~Group, data=sampleValueGroup)
			#piRNA.ANOVA.Summary = summary(piRNA.ANOVA)
			#piRNA.ANOVA.pvalue = piRNA.ANOVA.Summary[[1]]$'Pr(>F)'[1]
			#common_piRNA.ANOVA.pvalue <- c(common_piRNA.ANOVA.pvalue, piRNA.ANOVA.pvalue)
			
			# Fold Change
			Group1.value.Mean = mean(Group1.value[Group1.value.valididx], na.rm=TRUE)
			Group2.value.Mean = mean(Group2.value[Group2.value.valididx], na.rm=TRUE)
			
			log2_FC = Group2.value.Mean - Group1.value.Mean
			if(log2_FC >= 0) {			FC = 2^abs(log2_FC)
			} else if(log2_FC < 0) {	FC = -2^abs(log2_FC)
			}
			
			common_piRNA.log2_FC <- c(common_piRNA.log2_FC, log2_FC)
			common_piRNA.FC <- c(common_piRNA.FC, FC)
		} else if(sum(Group1.value, na.rm=TRUE)==0 & sum(Group2.value, na.rm=TRUE)==0){
			NA.piRNA.rowidx <- c(NA.piRNA.rowidx, i)
		}
		
		if(ggplot.plotidx == 1) {
			pic <- ggplot(sampleValueGroup, aes(x = Group, y = Value)) + 
				geom_boxplot(linetype = "dashed", outlier.shape = NA, color = "black", na.rm = TRUE) + 
				stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..), outlier.shape = NA, color = "black", na.rm = TRUE) + 
				stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width = 0.3, color = "black", na.rm = TRUE) + 
				stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..), width = 0.3, color = "black", na.rm = TRUE)
			
			#pic <- pic + geom_dotplot(binaxis = 'y', stackdir = 'center', stackratio = 1.25, dotsize = 0.15, binwidth = 0.6, alpha = 0.5, na.rm = TRUE)
			pic <- pic + geom_jitter(width = 0.15, height = 0, alpha = 0.5, na.rm = TRUE)
			
			pic <- pic + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
			pic <- pic + 
				labs(title = main_title, x = "Symp Group", y = "piRNA(normalized)", fill = "Symp Group", color = "Symp Group") + 
				theme(axis.text.x = element_text(color = "black", size = 12, margin = margin(3,0,0,0)), 
					  axis.text.y = element_text(color = "black", size = 12, margin = margin(0,3,0,0)), 
					  axis.title.x = element_text(color = "black", size = 14, margin = margin(10,0,0,0)), 
					  axis.title.y = element_text(color = "black", size = 14, margin = margin(0,10,0,0)), 
					  plot.title = element_text(color = "black", size = 16, margin = margin(0,0,15,0), hjust = 0.5, face = "bold"), 
					  legend.text = element_text(color = "black", size = 8), 
					  legend.title = element_text(color = "black", size = 8))
			
			pic <- pic + stat_summary(fun.y = mean, geom = "point", shape = 4, size = 1.5, stroke = 1.0, color = "red", na.rm = TRUE)
			
			print(pic)
			
			dev.off()
		}
	}
	
	# --- zero value distribution ---------------------------------------------------------------------
	gc()	# Garbage Collection (for release memory)
	
	#install.packages("reshape2")		# reshape2: Flexibly Reshape Data: A Reboot of the Reshape Package
	suppressPackageStartupMessages(library(reshape2))
	
	if(Group[1] == 'PDMCI' | Group[2] == 'PDMCI') {
		if(Group[1] == 'PDMCI') {
			Group1.colname = 'PD-MCI'
			Group2.colname = Group[2]
		}
		if(Group[2] == 'PDMCI') {
			Group1.colname = Group[1]
			Group2.colname = 'PD-MCI'
		}
	} else {
		Group1.colname = Group[1]
		Group2.colname = Group[2]
	}
	
	# all.ZeroCount
	Group1.all.ZeroCount = table(sapply(c(1:nrow(piRNA.normalizedData.Comp)), FUN=function(X) length(which(piRNA.normalizedData.Comp[X, Group1.colidx] == 0 | is.na(piRNA.normalizedData.Comp[X, Group1.colidx])))))
	Group2.all.ZeroCount = table(sapply(c(1:nrow(piRNA.normalizedData.Comp)), FUN=function(X) length(which(piRNA.normalizedData.Comp[X, Group2.colidx] == 0 | is.na(piRNA.normalizedData.Comp[X, Group2.colidx])))))
	Group1_Group2.all.ZeroCount = merge(as.data.frame(Group1.all.ZeroCount, stringsAsFactors=FALSE), as.data.frame(Group2.all.ZeroCount, stringsAsFactors=FALSE), by = "Var1", all = TRUE, sort = FALSE)
	colnames(Group1_Group2.all.ZeroCount) = c('all.ZeroCount',Group1.colname,Group2.colname)
	Group1_Group2.all.ZeroCount[,'all.ZeroCount'] = as.integer(Group1_Group2.all.ZeroCount[,'all.ZeroCount'])
	Group1_Group2.all.ZeroCount = Group1_Group2.all.ZeroCount[order(Group1_Group2.all.ZeroCount[,'all.ZeroCount']), ]
	rownames(Group1_Group2.all.ZeroCount) = c(1:nrow(Group1_Group2.all.ZeroCount))
	
	if(Group[1] == 'PD' | Group[2] == 'PD') {
		PDND.all.ZeroCount = table(sapply(c(1:nrow(piRNA.normalizedData.Comp)), FUN=function(X) length(which(piRNA.normalizedData.Comp[X, PDND.colidx] == 0 | is.na(piRNA.normalizedData.Comp[X, PDND.colidx])))))
		PDMCI.all.ZeroCount = table(sapply(c(1:nrow(piRNA.normalizedData.Comp)), FUN=function(X) length(which(piRNA.normalizedData.Comp[X, PDMCI.colidx] == 0 | is.na(piRNA.normalizedData.Comp[X, PDMCI.colidx])))))
		PDD.all.ZeroCount = table(sapply(c(1:nrow(piRNA.normalizedData.Comp)), FUN=function(X) length(which(piRNA.normalizedData.Comp[X, PDD.colidx] == 0 | is.na(piRNA.normalizedData.Comp[X, PDD.colidx])))))
		if(Group[1] == 'PD') {
			X = merge(as.data.frame(PDND.all.ZeroCount, stringsAsFactors=FALSE), as.data.frame(PDMCI.all.ZeroCount, stringsAsFactors=FALSE), by = "Var1", all = TRUE, sort = FALSE)
			colnames(X)[which(colnames(X)=='Freq.x')] = 'PDND'
			colnames(X)[which(colnames(X)=='Freq.y')] = 'PD-MCI'
			Y = merge(X, as.data.frame(PDD.all.ZeroCount, stringsAsFactors=FALSE), by = "Var1", all = TRUE, sort = FALSE)
			colnames(Y)[which(colnames(Y)=='Freq')] = 'PDD'
			Z = merge(Y, as.data.frame(Group2.all.ZeroCount, stringsAsFactors=FALSE), by = "Var1", all = TRUE, sort = FALSE)
			colnames(Z)[which(colnames(Z)=='Freq')] = Group2.colname
		}
		if(Group[2] == 'PD') {
			X = merge(as.data.frame(Group1.all.ZeroCount, stringsAsFactors=FALSE), as.data.frame(PDND.all.ZeroCount, stringsAsFactors=FALSE), by = "Var1", all = TRUE, sort = FALSE)
			colnames(X)[which(colnames(X)=='Freq.x')] = Group1.colname
			colnames(X)[which(colnames(X)=='Freq.y')] = 'PDND'
			Y = merge(X, as.data.frame(PDMCI.all.ZeroCount, stringsAsFactors=FALSE), by = "Var1", all = TRUE, sort = FALSE)
			colnames(Y)[which(colnames(Y)=='Freq')] = 'PD-MCI'
			Z = merge(Y, as.data.frame(PDD.all.ZeroCount, stringsAsFactors=FALSE), by = "Var1", all = TRUE, sort = FALSE)
			colnames(Z)[which(colnames(Z)=='Freq')] = 'PDD'
		}
		colnames(Z)[which(colnames(Z)=='Var1')] = 'all.ZeroCount'
		Z[,'all.ZeroCount'] = as.integer(Z[,'all.ZeroCount'])
		Z = Z[order(Z[,'all.ZeroCount']), ]
		rownames(Z) = c(1:nrow(Z))
		Group1_Group2.all.ZeroCount = Z
	}
	
	# common.ZeroCount
	Group1.common.ZeroCount = table(sapply(c(1:length(Group1_Group2.piRNA.rowidx)), FUN=function(X) length(which(piRNA.normalizedData.Comp[Group1_Group2.piRNA.rowidx[X], Group1.colidx] == 0 | is.na(piRNA.normalizedData.Comp[Group1_Group2.piRNA.rowidx[X], Group1.colidx])))))
	Group2.common.ZeroCount = table(sapply(c(1:length(Group1_Group2.piRNA.rowidx)), FUN=function(X) length(which(piRNA.normalizedData.Comp[Group1_Group2.piRNA.rowidx[X], Group2.colidx] == 0 | is.na(piRNA.normalizedData.Comp[Group1_Group2.piRNA.rowidx[X], Group2.colidx])))))
	Group1_Group2.common.ZeroCount = merge(as.data.frame(Group1.common.ZeroCount, stringsAsFactors=FALSE), as.data.frame(Group2.common.ZeroCount, stringsAsFactors=FALSE), by = "Var1", all = TRUE, sort = FALSE)
	colnames(Group1_Group2.common.ZeroCount) = c('common.ZeroCount',Group1.colname,Group2.colname)
	Group1_Group2.common.ZeroCount[,'common.ZeroCount'] = as.integer(Group1_Group2.common.ZeroCount[,'common.ZeroCount'])
	Group1_Group2.common.ZeroCount = Group1_Group2.common.ZeroCount[order(Group1_Group2.common.ZeroCount[,'common.ZeroCount']), ]
	rownames(Group1_Group2.common.ZeroCount) = c(1:nrow(Group1_Group2.common.ZeroCount))
	
	if(Group[1] == 'PD' | Group[2] == 'PD') {
		PDND.common.ZeroCount = table(sapply(c(1:length(Group1_Group2.piRNA.rowidx)), FUN=function(X) length(which(piRNA.normalizedData.Comp[Group1_Group2.piRNA.rowidx[X], PDND.colidx] == 0 | is.na(piRNA.normalizedData.Comp[Group1_Group2.piRNA.rowidx[X], PDND.colidx])))))
		PDMCI.common.ZeroCount = table(sapply(c(1:length(Group1_Group2.piRNA.rowidx)), FUN=function(X) length(which(piRNA.normalizedData.Comp[Group1_Group2.piRNA.rowidx[X], PDMCI.colidx] == 0 | is.na(piRNA.normalizedData.Comp[Group1_Group2.piRNA.rowidx[X], PDMCI.colidx])))))
		PDD.common.ZeroCount = table(sapply(c(1:length(Group1_Group2.piRNA.rowidx)), FUN=function(X) length(which(piRNA.normalizedData.Comp[Group1_Group2.piRNA.rowidx[X], PDD.colidx] == 0 | is.na(piRNA.normalizedData.Comp[Group1_Group2.piRNA.rowidx[X], PDD.colidx])))))
		if(Group[1] == 'PD') {
			X = merge(as.data.frame(PDND.common.ZeroCount, stringsAsFactors=FALSE), as.data.frame(PDMCI.common.ZeroCount, stringsAsFactors=FALSE), by = "Var1", all = TRUE, sort = FALSE)
			colnames(X)[which(colnames(X)=='Freq.x')] = 'PDND'
			colnames(X)[which(colnames(X)=='Freq.y')] = 'PD-MCI'
			Y = merge(X, as.data.frame(PDD.common.ZeroCount, stringsAsFactors=FALSE), by = "Var1", all = TRUE, sort = FALSE)
			colnames(Y)[which(colnames(Y)=='Freq')] = 'PDD'
			Z = merge(Y, as.data.frame(Group2.common.ZeroCount, stringsAsFactors=FALSE), by = "Var1", all = TRUE, sort = FALSE)
			colnames(Z)[which(colnames(Z)=='Freq')] = Group2.colname
		}
		if(Group[2] == 'PD') {
			X = merge(as.data.frame(Group1.common.ZeroCount, stringsAsFactors=FALSE), as.data.frame(PDND.common.ZeroCount, stringsAsFactors=FALSE), by = "Var1", all = TRUE, sort = FALSE)
			colnames(X)[which(colnames(X)=='Freq.x')] = Group1.colname
			colnames(X)[which(colnames(X)=='Freq.y')] = 'PDND'
			Y = merge(X, as.data.frame(PDMCI.common.ZeroCount, stringsAsFactors=FALSE), by = "Var1", all = TRUE, sort = FALSE)
			colnames(Y)[which(colnames(Y)=='Freq')] = 'PD-MCI'
			Z = merge(Y, as.data.frame(PDD.common.ZeroCount, stringsAsFactors=FALSE), by = "Var1", all = TRUE, sort = FALSE)
			colnames(Z)[which(colnames(Z)=='Freq')] = 'PDD'
		}
		colnames(Z)[which(colnames(Z)=='Var1')] = 'common.ZeroCount'
		Z[,'common.ZeroCount'] = as.integer(Z[,'common.ZeroCount'])
		Z = Z[order(Z[,'common.ZeroCount']), ]
		rownames(Z) = c(1:nrow(Z))
		Group1_Group2.common.ZeroCount = Z
	}
	
	for(z in 1:2)
	{
		if(z == 1) {
			Group1_Group2.ZeroCount = Group1_Group2.all.ZeroCount
			colnames(Group1_Group2.ZeroCount)[which(colnames(Group1_Group2.ZeroCount)=='all.ZeroCount')] = 'ZeroCount'
			png(paste('./Data/02_significant/piRNA/',Group[1],' vs ',Group[2],'/','piRNA_ZeroValueDistribution ','(',Group[1],' vs ',Group[2],' (All)',')','.png',sep=''), width=12, height=6, units="in", res=600)
		} else if(z == 2) {
			Group1_Group2.ZeroCount = Group1_Group2.common.ZeroCount
			colnames(Group1_Group2.ZeroCount)[which(colnames(Group1_Group2.ZeroCount)=='common.ZeroCount')] = 'ZeroCount'
			png(paste('./Data/02_significant/piRNA/',Group[1],' vs ',Group[2],'/','piRNA_ZeroValueDistribution ','(',Group[1],' vs ',Group[2],' (Common)',')','.png',sep=''), width=12, height=6, units="in", res=600)
		}
		
		if(Group[1] == 'PD' | Group[2] == 'PD') {
			if(Group[1] == 'PD') {
				samples_num = c(length(PDND.colidx),length(PDMCI.colidx),length(PDD.colidx),length(Group2.colidx))
				NA.idx = setdiff(c(0:max(samples_num)), Group1_Group2.ZeroCount[,'ZeroCount'])
				NA.df = data.frame(NA.idx, rep(NA,length(NA.idx)), rep(NA,length(NA.idx)), rep(NA,length(NA.idx)), rep(NA,length(NA.idx)), stringsAsFactors=FALSE)
				colnames(NA.df) = c('ZeroCount',c('PDND','PD-MCI','PDD'),Group2.colname)
				NA.df[,'ZeroCount'] = as.integer(NA.df[,'ZeroCount'])
			}
			if(Group[2] == 'PD') {
				samples_num = c(length(Group1.colidx),length(PDND.colidx),length(PDMCI.colidx),length(PDD.colidx))
				NA.idx = setdiff(c(0:max(samples_num)), Group1_Group2.ZeroCount[,'ZeroCount'])
				NA.df = data.frame(NA.idx, rep(NA,length(NA.idx)), rep(NA,length(NA.idx)), rep(NA,length(NA.idx)), rep(NA,length(NA.idx)), stringsAsFactors=FALSE)
				colnames(NA.df) = c('ZeroCount',Group1.colname,c('PDND','PD-MCI','PDD'))
				NA.df[,'ZeroCount'] = as.integer(NA.df[,'ZeroCount'])
			}
		} else {
			samples_num = c(length(Group1.colidx),length(Group2.colidx))
			NA.idx = setdiff(c(0:max(samples_num)), Group1_Group2.ZeroCount[,'ZeroCount'])
			NA.df = data.frame(NA.idx, rep(NA,length(NA.idx)), rep(NA,length(NA.idx)), stringsAsFactors=FALSE)
			colnames(NA.df) = c('ZeroCount',Group1.colname,Group2.colname)
			NA.df[,'ZeroCount'] = as.integer(NA.df[,'ZeroCount'])
		}
		
		Group1_Group2.ZeroCount = rbind(Group1_Group2.ZeroCount, NA.df)
		Group1_Group2.ZeroCount = Group1_Group2.ZeroCount[order(Group1_Group2.ZeroCount[,'ZeroCount']), ]
		rownames(Group1_Group2.ZeroCount) = c(1:nrow(Group1_Group2.ZeroCount))
		Group1_Group2.ZeroCount[,'ZeroCount'] = as.factor(Group1_Group2.ZeroCount[,'ZeroCount'])
		
		Group1_Group2.ZeroCount = melt(Group1_Group2.ZeroCount, id.vars = "ZeroCount", na.rm = FALSE)
		colnames(Group1_Group2.ZeroCount) = c('ZeroCount','Group','Value')
		Group1_Group2.ZeroCount = Group1_Group2.ZeroCount[, c('ZeroCount','Value','Group')]
		
		pic <- ggplot(Group1_Group2.ZeroCount, aes(x = ZeroCount, y = Value, fill = Group)) + geom_bar(stat = "identity", position = position_dodge(preserve = 'single'), na.rm = TRUE)
		pic <- pic + scale_fill_manual("Symp Group", values = diseaseStage.color)
		
		pic <- pic + theme_bw() #+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
		pic <- pic + 
			labs(title = "Zero value distribution", x = "Number of sample with zero value", y = "Number of piRNA") + 
			theme(axis.text.x = element_text(color = "black", size = 16, margin = margin(3,0,0,0)), 
				  axis.text.y = element_text(color = "black", size = 16, margin = margin(0,3,0,0)), 
				  axis.title.x = element_text(color = "black", size = 18, margin = margin(10,0,0,0)), 
				  axis.title.y = element_text(color = "black", size = 18, margin = margin(0,10,0,0)), 
				  plot.title = element_text(color = "black", size = 20, margin = margin(0,0,15,0), hjust = 0.5, face = "bold"), 
				  legend.text = element_text(color = "black", size = 12), 
				  legend.title = element_text(color = "black", size = 12))
		
		print(pic)
		
		dev.off()
	}
	
	# --- significance --------------------------------------------------------------------------------
	pvalue.cutoff = 0.01
	pvalue_FDR.cutoff = 0.01
	
	FC.cutoff = 1.5
	
	common_piRNA.ttest.pvalue_adjust = p.adjust(common_piRNA.ttest.pvalue, method="fdr")
	#common_piRNA.ANOVA.pvalue_adjust = p.adjust(common_piRNA.ANOVA.pvalue, method="fdr")
	
	piRNA.normalizedData.Comp.stat = cbind(piRNA.normalizedData.Comp[Group1_Group2.piRNA.rowidx, , drop=FALSE], common_piRNA.ttest.pvalue, common_piRNA.ttest.pvalue_adjust, common_piRNA.log2_FC, common_piRNA.FC, stringsAsFactors=FALSE)
	colnames(piRNA.normalizedData.Comp.stat)[(ncol(piRNA.normalizedData.Comp.stat)-3):ncol(piRNA.normalizedData.Comp.stat)] = c('pvalue','pvalue_FDR','log2_FC','FC')
	piRNA.normalizedData.Comp.stat = piRNA.normalizedData.Comp.stat[, c(piRNA.annotation.idx, (ncol(piRNA.normalizedData.Comp.stat)-3):ncol(piRNA.normalizedData.Comp.stat), (length(piRNA.annotation.idx)+1):(ncol(piRNA.normalizedData.Comp.stat)-4)), drop=FALSE]
	
	# significant index
	significant_piRNA.idx = Group1_Group2.piRNA.rowidx[which(piRNA.normalizedData.Comp.stat[,'pvalue'] < pvalue.cutoff)]
	nonsignificant_piRNA.idx = Group1_Group2.piRNA.rowidx[which(piRNA.normalizedData.Comp.stat[,'pvalue'] >= pvalue.cutoff)]
	
	significant_piRNA_FCup.idx = Group1_Group2.piRNA.rowidx[which(piRNA.normalizedData.Comp.stat[,'pvalue'] < pvalue.cutoff & piRNA.normalizedData.Comp.stat[,'FC'] > 0)]
	significant_piRNA_FCdown.idx = Group1_Group2.piRNA.rowidx[which(piRNA.normalizedData.Comp.stat[,'pvalue'] < pvalue.cutoff & piRNA.normalizedData.Comp.stat[,'FC'] < 0)]
	significant_piRNA_noFC.idx = Group1_Group2.piRNA.rowidx[which(piRNA.normalizedData.Comp.stat[,'pvalue'] < pvalue.cutoff & piRNA.normalizedData.Comp.stat[,'FC'] == 0)]
	
	if(length(significant_piRNA.idx) > 0) {
		mainDir <- paste('./Data/02_significant/piRNA/',Group[1],' vs ',Group[2],'/','Boxplot',sep='')
		subDir <- paste(Group[1],' vs ',Group[2],sep='')
		if(file.exists(file.path(mainDir, subDir))) {
		} else {
			dir.create(file.path(mainDir, subDir))
		}
		mainDir <- paste('./Data/02_significant/piRNA/',Group[1],' vs ',Group[2],'/','Boxplot','/',Group[1],' vs ',Group[2],sep='')
		subDir <- 'significant'
		if(file.exists(file.path(mainDir, subDir))) {
		} else {
			dir.create(file.path(mainDir, subDir))
		}
		
		for(i in 1:length(significant_piRNA.idx))
		{
			main_title = unlist(strsplit(piRNA.normalizedData.Comp[significant_piRNA.idx[i],'Name'], '/'))[1]
			
			sampleValue = piRNA.normalizedData.Comp[significant_piRNA.idx[i], , drop=FALSE]
			sampleGroup = c(rep(NA, length(piRNA.annotation)), rep(Group[1], length(Group1.colidx)), rep(Group[2], length(Group2.colidx)))
			
			sampleValueGroup = as.data.frame(t(rbind(sampleValue, sampleGroup)[, -piRNA.annotation.idx, drop=FALSE]), stringsAsFactors=FALSE)
			colnames(sampleValueGroup) = c('Value','Group')
			sampleValueGroup[,'Value'] = as.numeric(sampleValueGroup[,'Value'])
			sampleValueGroup[,'Group'] = as.factor(sampleValueGroup[,'Group'])
			sampleValueGroup[,'Group'] = factor(sampleValueGroup[,'Group'], levels(sampleValueGroup[,'Group'])[c(which(levels(sampleValueGroup[,'Group'])==Group[1]), 
																												 which(levels(sampleValueGroup[,'Group'])==Group[2]))])
			
			ggplot.plotidx = 0
			if(significant_piRNA.idx[i] %in% significant_piRNA_FCup.idx) {
				mainDir <- paste('./Data/02_significant/piRNA/',Group[1],' vs ',Group[2],'/','Boxplot','/',Group[1],' vs ',Group[2],'/significant',sep='')
				subDir <- 'FCup'
				if(file.exists(file.path(mainDir, subDir))) {
				} else {
					dir.create(file.path(mainDir, subDir))
				}
				
				# /// boxplot {graphics} ///
				#png(paste('./Data/02_significant/piRNA/',Group[1],' vs ',Group[2],'/','Boxplot','/',Group[1],' vs ',Group[2],'/significant/FCup/',gsub('/', '.', piRNA.normalizedData.Comp[significant_piRNA.idx[i],'Name']),'.png',sep=''), width=4, height=4, units="in", res=600)
				#boxplot(Value~Group, data=sampleValueGroup, main=main_title, xlab="Symp Group", ylab="piRNA(normalized)")
				#dev.off()
				
				# /// ggplot {ggplot2} ///
				png(paste('./Data/02_significant/piRNA/',Group[1],' vs ',Group[2],'/','Boxplot','/',Group[1],' vs ',Group[2],'/significant/FCup/',gsub('/', '.', piRNA.normalizedData.Comp[significant_piRNA.idx[i],'Name']),'.png',sep=''), width=4, height=4, units="in", res=600)
				ggplot.plotidx = 1
			} else if(significant_piRNA.idx[i] %in% significant_piRNA_FCdown.idx) {
				mainDir <- paste('./Data/02_significant/piRNA/',Group[1],' vs ',Group[2],'/','Boxplot','/',Group[1],' vs ',Group[2],'/significant',sep='')
				subDir <- 'FCdown'
				if(file.exists(file.path(mainDir, subDir))) {
				} else {
					dir.create(file.path(mainDir, subDir))
				}
				
				# /// boxplot {graphics} ///
				#png(paste('./Data/02_significant/piRNA/',Group[1],' vs ',Group[2],'/','Boxplot','/',Group[1],' vs ',Group[2],'/significant/FCdown/',gsub('/', '.', piRNA.normalizedData.Comp[significant_piRNA.idx[i],'Name']),'.png',sep=''), width=4, height=4, units="in", res=600)
				#boxplot(Value~Group, data=sampleValueGroup, main=main_title, xlab="Symp Group", ylab="piRNA(normalized)")
				#dev.off()
				
				# /// ggplot {ggplot2} ///
				png(paste('./Data/02_significant/piRNA/',Group[1],' vs ',Group[2],'/','Boxplot','/',Group[1],' vs ',Group[2],'/significant/FCdown/',gsub('/', '.', piRNA.normalizedData.Comp[significant_piRNA.idx[i],'Name']),'.png',sep=''), width=4, height=4, units="in", res=600)
				ggplot.plotidx = 1
			} else if(significant_piRNA.idx[i] %in% significant_piRNA_noFC.idx) {
				mainDir <- paste('./Data/02_significant/piRNA/',Group[1],' vs ',Group[2],'/','Boxplot','/',Group[1],' vs ',Group[2],'/significant',sep='')
				subDir <- 'noFC'
				if(file.exists(file.path(mainDir, subDir))) {
				} else {
					dir.create(file.path(mainDir, subDir))
				}
				
				# /// boxplot {graphics} ///
				#png(paste('./Data/02_significant/piRNA/',Group[1],' vs ',Group[2],'/','Boxplot','/',Group[1],' vs ',Group[2],'/significant/noFC/',gsub('/', '.', piRNA.normalizedData.Comp[significant_piRNA.idx[i],'Name']),'.png',sep=''), width=4, height=4, units="in", res=600)
				#boxplot(Value~Group, data=sampleValueGroup, main=main_title, xlab="Symp Group", ylab="piRNA(normalized)")
				#dev.off()
				
				# /// ggplot {ggplot2} ///
				png(paste('./Data/02_significant/piRNA/',Group[1],' vs ',Group[2],'/','Boxplot','/',Group[1],' vs ',Group[2],'/significant/noFC/',gsub('/', '.', piRNA.normalizedData.Comp[significant_piRNA.idx[i],'Name']),'.png',sep=''), width=4, height=4, units="in", res=600)
				ggplot.plotidx = 1
			}
			
			if(ggplot.plotidx == 1) {
				pic <- ggplot(sampleValueGroup, aes(x = Group, y = Value)) + 
					geom_boxplot(linetype = "dashed", outlier.shape = NA, color = "black", na.rm = TRUE) + 
					stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..), outlier.shape = NA, color = "black", na.rm = TRUE) + 
					stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width = 0.3, color = "black", na.rm = TRUE) + 
					stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..), width = 0.3, color = "black", na.rm = TRUE)
				
				#pic <- pic + geom_dotplot(binaxis = 'y', stackdir = 'center', stackratio = 1.25, dotsize = 0.15, binwidth = 0.6, alpha = 0.5, na.rm = TRUE)
				pic <- pic + geom_jitter(width = 0.15, height = 0, alpha = 0.5, na.rm = TRUE)
				
				pic <- pic + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
				pic <- pic + 
					labs(title = main_title, x = "Symp Group", y = "piRNA(normalized)", fill = "Symp Group", color = "Symp Group") + 
					theme(axis.text.x = element_text(color = "black", size = 12, margin = margin(3,0,0,0)), 
						  axis.text.y = element_text(color = "black", size = 12, margin = margin(0,3,0,0)), 
						  axis.title.x = element_text(color = "black", size = 14, margin = margin(10,0,0,0)), 
						  axis.title.y = element_text(color = "black", size = 14, margin = margin(0,10,0,0)), 
						  plot.title = element_text(color = "black", size = 16, margin = margin(0,0,15,0), hjust = 0.5, face = "bold"), 
						  legend.text = element_text(color = "black", size = 8), 
						  legend.title = element_text(color = "black", size = 8))
				
				pic <- pic + stat_summary(fun.y = mean, geom = "point", shape = 4, size = 1.5, stroke = 1.0, color = "red", na.rm = TRUE)
				
				print(pic)
				
				dev.off()
			}
		}
	}
	
	# --- Volcano plot --------------------------------------------------------------------------------
	gc()	# Garbage Collection (for release memory)
	
	mainDir <- paste('./Data/02_significant/piRNA/',Group[1],' vs ',Group[2],sep='')
	subDir <- 'Volcano'
	if(file.exists(file.path(mainDir, subDir))) {
	} else {
		dir.create(file.path(mainDir, subDir))
	}
	
	sampleStat = piRNA.normalizedData.Comp.stat[, c('Name','pvalue','pvalue_FDR','log2_FC','FC'), drop=FALSE]
	sampleStat[,'Name'] = as.factor(sampleStat[,'Name'])
	sampleStat = cbind(sampleStat, neg_log10_pvalue = -log10(sampleStat[,'pvalue']), neg_log10_pvalue_FDR = -log10(sampleStat[,'pvalue_FDR']), stringsAsFactors=FALSE)
	
	pvalue.Region.idx = which(sampleStat[,'pvalue'] < pvalue.cutoff)
	FC.Region.idx = which(sampleStat[,'FC'] > FC.cutoff | sampleStat[,'FC'] < -FC.cutoff)
	pvalue_FC.Region.idx = sort(intersect(pvalue.Region.idx, FC.Region.idx), decreasing = FALSE)
	sampleStat.pvalue.Region = sampleStat[pvalue.Region.idx, , drop=FALSE]
	sampleStat.FC.Region = sampleStat[FC.Region.idx, , drop=FALSE]
	sampleStat.pvalue_FC.Region = sampleStat[pvalue_FC.Region.idx, , drop=FALSE]
	
	# https://digibio.blogspot.com/2018/05/volcano-plot-with-ggplot2-and-basic.html
	# https://www.bioconductor.org/packages/release/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html
	# https://cran.r-project.org/web/packages/ggrepel/vignettes/ggrepel.html
	
	#install.packages("dplyr")			# dplyr: A Grammar of Data Manipulation
	#install.packages("ggplot2")		# ggplot2: Create Elegant Data Visualisations Using the Grammar of Graphics
	#install.packages("ggrepel")		# ggrepel: Automatically Position Non-Overlapping Text Labels with 'ggplot2'
	suppressPackageStartupMessages(library(dplyr))
	suppressPackageStartupMessages(library(ggplot2))
	suppressPackageStartupMessages(library(ggrepel))
	
	png(paste('./Data/02_significant/piRNA/',Group[1],' vs ',Group[2],'/','Volcano','/','piRNA_volcano ','(',Group[1],' vs ',Group[2],')','.png',sep=''), width=6, height=6, units="in", res=600)
	
	pic <- ggplot(sampleStat) + geom_point(data = sampleStat, aes(x = log2_FC, y = neg_log10_pvalue), color = "gray", alpha = 0.5, cex = 2)
	
	pic <- pic + geom_point(data = sampleStat.pvalue.Region, aes(x = log2_FC, y = neg_log10_pvalue), color = "blue", alpha = 0.5, cex = 2)
	pic <- pic + geom_point(data = sampleStat.FC.Region, aes(x = log2_FC, y = neg_log10_pvalue), color = "green", alpha = 0.5, cex = 2)
	pic <- pic + geom_point(data = sampleStat.pvalue_FC.Region, aes(x = log2_FC, y = neg_log10_pvalue), color = "red", alpha = 0.5, cex = 2)
	
	pic <- pic + 
		geom_vline(xintercept = log2(FC.cutoff), col = "red", linetype = "dotted", size = 1) + 
		geom_vline(xintercept = -log2(FC.cutoff), col = "red", linetype = "dotted", size = 1) + 
		geom_hline(yintercept = -log10(pvalue.cutoff), col = "red", linetype = "dashed", size = 1)
	
	label_count = max(c(nrow(subset(sampleStat.pvalue_FC.Region, FC > FC.cutoff)), nrow(subset(sampleStat.pvalue_FC.Region, FC < -FC.cutoff))))
	if(label_count <= 30) {								cex_num = 3
	} else if(label_count > 30 & label_count <= 40) {	cex_num = 2.5
	} else if(label_count > 40 & label_count <= 50) {	cex_num = 2
	} else if(label_count > 50) {						cex_num = 1.5
	}
	pic <- pic + 
		geom_text_repel(
			data          = subset(sampleStat.pvalue_FC.Region, FC > 3.0), 
			aes(x = log2_FC, y = neg_log10_pvalue, label = rownames(subset(sampleStat.pvalue_FC.Region, FC > 3.0))), 
			nudge_x       = 6.5 - subset(sampleStat.pvalue_FC.Region, FC > 3.0)$log2_FC, 
			segment.size  = 0.2, 
			segment.color = "grey50", 
			direction     = "y", 
			hjust         = 0, 
			cex           = cex_num
		) + 
		geom_text_repel(
			data          = subset(sampleStat.pvalue_FC.Region, FC > 2.0 & FC <= 3.0), 
			aes(x = log2_FC, y = neg_log10_pvalue, label = rownames(subset(sampleStat.pvalue_FC.Region, FC > 2.0 & FC <= 3.0))), 
			nudge_x       = 4.5 - subset(sampleStat.pvalue_FC.Region, FC > 2.0 & FC <= 3.0)$log2_FC, 
			segment.size  = 0.2, 
			segment.color = "grey50", 
			direction     = "y", 
			hjust         = 0, 
			cex           = cex_num
		) + 
		geom_text_repel(
			data          = subset(sampleStat.pvalue_FC.Region, FC > 1.5 & FC <= 2.0), 
			aes(x = log2_FC, y = neg_log10_pvalue, label = rownames(subset(sampleStat.pvalue_FC.Region, FC > 1.5 & FC <= 2.0))), 
			nudge_x       = 2.5 - subset(sampleStat.pvalue_FC.Region, FC > 1.5 & FC <= 2.0)$log2_FC, 
			segment.size  = 0.2, 
			segment.color = "grey50", 
			direction     = "y", 
			hjust         = 0, 
			cex           = cex_num
		) + 
		geom_text_repel(
			data          = subset(sampleStat.pvalue_FC.Region, FC < -1.5 & FC >= -2.0), 
			aes(x = log2_FC, y = neg_log10_pvalue, label = rownames(subset(sampleStat.pvalue_FC.Region, FC < -1.5 & FC >= -2.0))), 
			nudge_x       = -2.5 - subset(sampleStat.pvalue_FC.Region, FC < -1.5 & FC >= -2.0)$log2_FC, 
			segment.size  = 0.2, 
			segment.color = "grey50", 
			direction     = "y", 
			hjust         = 1, 
			cex           = cex_num
		) + 
		geom_text_repel(
			data          = subset(sampleStat.pvalue_FC.Region, FC < -2.0 & FC >= -3.0), 
			aes(x = log2_FC, y = neg_log10_pvalue, label = rownames(subset(sampleStat.pvalue_FC.Region, FC < -2.0 & FC >= -3.0))), 
			nudge_x       = -4.5 - subset(sampleStat.pvalue_FC.Region, FC < -2.0 & FC >= -3.0)$log2_FC, 
			segment.size  = 0.2, 
			segment.color = "grey50", 
			direction     = "y", 
			hjust         = 1, 
			cex           = cex_num
		) + 
		geom_text_repel(
			data          = subset(sampleStat.pvalue_FC.Region, FC < -3.0), 
			aes(x = log2_FC, y = neg_log10_pvalue, label = rownames(subset(sampleStat.pvalue_FC.Region, FC < -3.0))), 
			nudge_x       = -6.5 - subset(sampleStat.pvalue_FC.Region, FC < -3.0)$log2_FC, 
			segment.size  = 0.2, 
			segment.color = "grey50", 
			direction     = "y", 
			hjust         = 1, 
			cex           = cex_num
		)
	
	xlimit_num = (ceiling(max(c(abs(max(sampleStat[,'log2_FC'])), abs(min(sampleStat[,'log2_FC'])))) / 2)) * 2
	if(xlimit_num == 4)	xlimit_num = xlimit_num + 2
	ylimit_num = (ceiling(max(c(abs(max(sampleStat[,'neg_log10_pvalue'])), abs(min(sampleStat[,'neg_log10_pvalue'])))) / 2)) * 2
	if(ylimit_num < 6)	ylimit_num = 6
	pic <- pic + 
		scale_x_continuous(
			breaks = seq(-xlimit_num, xlimit_num), 
			limits = c(-xlimit_num, xlimit_num)
		) + 
		scale_y_continuous(
			breaks = seq(0, ylimit_num), 
			limits = c(0, ylimit_num)
		)
	
	pic <- pic + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	pic <- pic + 
		ggtitle(paste(Group[1],' vs ',Group[2],sep='')) + 
		xlab(expression(log[2]~'fold change')) + ylab(expression(-log[10]~'pvalue')) + 
		theme(axis.text.x = element_text(color = "black", size = 16, margin = margin(3,0,0,0)), 
			  axis.text.y = element_text(color = "black", size = 16, margin = margin(0,3,0,0)), 
			  axis.title.x = element_text(color = "black", size = 18, margin = margin(10,0,0,0)), 
			  axis.title.y = element_text(color = "black", size = 18, margin = margin(0,10,0,0)), 
			  plot.title = element_text(color = "black", size = 20, margin = margin(0,0,15,0), hjust = 0.5, face = "bold"), 
			  legend.text = element_text(color = "black", size = 12), 
			  legend.title = element_text(color = "black", size = 12))
	
	print(pic)
	
	dev.off()
	
	# combine specific piRNA and significant common piRNA
	# Group1
	piRNA.normalizedData.Comp.Group1 <- c()
	if(length(Group1.piRNA.rowidx) > 0) {
		# specific piRNA ordered by "non-zero sample number count" and "mean of non-zero sample value"
		Group1.nonZeroCount = sapply(c(1:length(Group1.piRNA.rowidx)), FUN=function(X) length(which(piRNA.normalizedData.Comp[Group1.piRNA.rowidx[X], Group1.colidx] != 0)))
		Group1.nonZeroMean = sapply(c(1:length(Group1.piRNA.rowidx)), FUN=function(X) mean(as.numeric(piRNA.normalizedData.Comp[Group1.piRNA.rowidx[X], Group1.colidx][which(piRNA.normalizedData.Comp[Group1.piRNA.rowidx[X], Group1.colidx] != 0)]), na.rm=TRUE))
		Group1.piRNA.rowidx = Group1.piRNA.rowidx[order(Group1.nonZeroCount, Group1.nonZeroMean, decreasing=TRUE)]
		
		# specific piRNA
		piRNA.normalizedData.Comp.Group1 = cbind(piRNA.normalizedData.Comp[Group1.piRNA.rowidx, , drop=FALSE], pvalue=rep(NA, length(Group1.piRNA.rowidx)), pvalue_FDR=rep(NA, length(Group1.piRNA.rowidx)), log2_FC=rep(NA, length(Group1.piRNA.rowidx)), FC=rep(NA, length(Group1.piRNA.rowidx)), stringsAsFactors=FALSE)
		piRNA.normalizedData.Comp.Group1 = piRNA.normalizedData.Comp.Group1[, c(piRNA.annotation.idx, (ncol(piRNA.normalizedData.Comp.Group1)-3):ncol(piRNA.normalizedData.Comp.Group1), (length(piRNA.annotation.idx)+1):(ncol(piRNA.normalizedData.Comp.Group1)-4)), drop=FALSE]
	}
	# Group2
	piRNA.normalizedData.Comp.Group2 <- c()
	if(length(Group2.piRNA.rowidx) > 0) {
		# specific piRNA ordered by "non-zero sample number count" and "mean of non-zero sample value"
		Group2.nonZeroCount = sapply(c(1:length(Group2.piRNA.rowidx)), FUN=function(X) length(which(piRNA.normalizedData.Comp[Group2.piRNA.rowidx[X], Group2.colidx] != 0)))
		Group2.nonZeroMean = sapply(c(1:length(Group2.piRNA.rowidx)), FUN=function(X) mean(as.numeric(piRNA.normalizedData.Comp[Group2.piRNA.rowidx[X], Group2.colidx][which(piRNA.normalizedData.Comp[Group2.piRNA.rowidx[X], Group2.colidx] != 0)]), na.rm=TRUE))
		Group2.piRNA.rowidx = Group2.piRNA.rowidx[order(Group2.nonZeroCount, Group2.nonZeroMean, decreasing=TRUE)]
		
		# specific piRNA
		piRNA.normalizedData.Comp.Group2 = cbind(piRNA.normalizedData.Comp[Group2.piRNA.rowidx, , drop=FALSE], pvalue=rep(NA, length(Group2.piRNA.rowidx)), pvalue_FDR=rep(NA, length(Group2.piRNA.rowidx)), log2_FC=rep(NA, length(Group2.piRNA.rowidx)), FC=rep(NA, length(Group2.piRNA.rowidx)), stringsAsFactors=FALSE)
		piRNA.normalizedData.Comp.Group2 = piRNA.normalizedData.Comp.Group2[, c(piRNA.annotation.idx, (ncol(piRNA.normalizedData.Comp.Group2)-3):ncol(piRNA.normalizedData.Comp.Group2), (length(piRNA.annotation.idx)+1):(ncol(piRNA.normalizedData.Comp.Group2)-4)), drop=FALSE]
	}
	
	# significant common piRNA
	piRNA.normalizedData.Comp.stat.FCup = piRNA.normalizedData.Comp.stat[which(piRNA.normalizedData.Comp.stat[,'pvalue'] < pvalue.cutoff & piRNA.normalizedData.Comp.stat[,'FC'] > 0), , drop=FALSE]
	piRNA.normalizedData.Comp.stat.FCup = piRNA.normalizedData.Comp.stat.FCup[order(piRNA.normalizedData.Comp.stat.FCup[,'pvalue']), , drop=FALSE]
	piRNA.normalizedData.Comp.stat.FCdown = piRNA.normalizedData.Comp.stat[which(piRNA.normalizedData.Comp.stat[,'pvalue'] < pvalue.cutoff & piRNA.normalizedData.Comp.stat[,'FC'] < 0), , drop=FALSE]
	piRNA.normalizedData.Comp.stat.FCdown = piRNA.normalizedData.Comp.stat.FCdown[order(piRNA.normalizedData.Comp.stat.FCdown[,'pvalue']), , drop=FALSE]
	
	# Group1 > Group2
	if(!is.null(piRNA.normalizedData.Comp.Group1)) {
		df1 = piRNA.normalizedData.Comp.Group1
		df2 = piRNA.normalizedData.Comp.stat.FCdown
		df3 = rbind(
			data.frame(c(df1, sapply(setdiff(names(df2), names(df1)), function(x) NA)), check.names = FALSE, stringsAsFactors = FALSE), 
			data.frame(c(df2, sapply(setdiff(names(df1), names(df2)), function(x) NA)), check.names = FALSE, stringsAsFactors = FALSE)
		)
		rownames(df3) = unlist(sapply(strsplit(df3[,'Name'], '/'), FUN=function(X) X[1]))
		if(which.max(c(ncol(df1), ncol(df2))) == 1) {
			df3 = df3[, colnames(df1), drop=FALSE]
		} else if(which.max(c(ncol(df1), ncol(df2))) == 2) {
			df3 = df3[, colnames(df2), drop=FALSE]
		}
		piRNA.normalizedData.Comp.UPinGroup1 = df3
	} else {
		piRNA.normalizedData.Comp.UPinGroup1 = piRNA.normalizedData.Comp.stat.FCdown
	}
	# Group1 < Group2
	if(!is.null(piRNA.normalizedData.Comp.Group2)) {
		df1 = piRNA.normalizedData.Comp.Group2
		df2 = piRNA.normalizedData.Comp.stat.FCup
		df3 = rbind(
			data.frame(c(df1, sapply(setdiff(names(df2), names(df1)), function(x) NA)), check.names = FALSE, stringsAsFactors = FALSE), 
			data.frame(c(df2, sapply(setdiff(names(df1), names(df2)), function(x) NA)), check.names = FALSE, stringsAsFactors = FALSE)
		)
		rownames(df3) = unlist(sapply(strsplit(df3[,'Name'], '/'), FUN=function(X) X[1]))
		if(which.max(c(ncol(df1), ncol(df2))) == 1) {
			df3 = df3[, colnames(df1), drop=FALSE]
		} else if(which.max(c(ncol(df1), ncol(df2))) == 2) {
			df3 = df3[, colnames(df2), drop=FALSE]
		}
		piRNA.normalizedData.Comp.UPinGroup2 = df3
	} else {
		piRNA.normalizedData.Comp.UPinGroup2 = piRNA.normalizedData.Comp.stat.FCup
	}
	
	# --- Analysis ------------------------------------------------------------------------------------
	for(n in 1:6)
	{
		gc()	# Garbage Collection (for release memory)
		
		Group1.colidx.new = Group1.colidx + length(c('pvalue','pvalue_FDR','log2_FC','FC'))
		Group2.colidx.new = Group2.colidx + length(c('pvalue','pvalue_FDR','log2_FC','FC'))
		
		if(n == 1) {
			mat1 = as.matrix(piRNA.normalizedData.Comp.UPinGroup1[, Group1.colidx.new, drop=FALSE])
			mat2 = as.matrix(piRNA.normalizedData.Comp.UPinGroup1[, setdiff(Group2.colidx.new,Group1.colidx.new), drop=FALSE])
			mat3 = as.matrix(piRNA.normalizedData.Comp.UPinGroup1[, unique(c(Group1.colidx.new,Group2.colidx.new)), drop=FALSE])
			
			heatmap.sep.plot <- paste('piRNA_heatmap_by_','High',Group[1],' ','(',Group[1],' + ',Group[2],' (Specific+Common)',')',sep='')
			heatmap.mix.plot <- paste('piRNA_heatmap_by_','High',Group[1],' ','(',Group[1],' vs ',Group[2],' (Specific+Common)',')',sep='')
			dendrogram.plot <- paste('piRNA_dendrogram_by_','High',Group[1],' ','(',Group[1],' vs ',Group[2],' (Specific+Common)',')',sep='')
			pca2D.plot <- paste('piRNA_pca2D_by_','High',Group[1],' ','(',Group[1],' vs ',Group[2],' (Specific+Common)',')',sep='')
			pca3D.plot <- paste('piRNA_pca3D_by_','High',Group[1],' ','(',Group[1],' vs ',Group[2],' (Specific+Common)',')',sep='')
			plsda2D.plot <- paste('piRNA_plsda2D_by_','High',Group[1],' ','(',Group[1],' vs ',Group[2],' (Specific+Common)',')',sep='')
			plsda3D.plot <- paste('piRNA_plsda3D_by_','High',Group[1],' ','(',Group[1],' vs ',Group[2],' (Specific+Common)',')',sep='')
		} else if(n == 2) {
			mat1 = as.matrix(piRNA.normalizedData.Comp.UPinGroup2[, Group1.colidx.new, drop=FALSE])
			mat2 = as.matrix(piRNA.normalizedData.Comp.UPinGroup2[, setdiff(Group2.colidx.new,Group1.colidx.new), drop=FALSE])
			mat3 = as.matrix(piRNA.normalizedData.Comp.UPinGroup2[, unique(c(Group1.colidx.new,Group2.colidx.new)), drop=FALSE])
			
			heatmap.sep.plot <- paste('piRNA_heatmap_by_','High',Group[2],' ','(',Group[1],' + ',Group[2],' (Specific+Common)',')',sep='')
			heatmap.mix.plot <- paste('piRNA_heatmap_by_','High',Group[2],' ','(',Group[1],' vs ',Group[2],' (Specific+Common)',')',sep='')
			dendrogram.plot <- paste('piRNA_dendrogram_by_','High',Group[2],' ','(',Group[1],' vs ',Group[2],' (Specific+Common)',')',sep='')
			pca2D.plot <- paste('piRNA_pca2D_by_','High',Group[2],' ','(',Group[1],' vs ',Group[2],' (Specific+Common)',')',sep='')
			pca3D.plot <- paste('piRNA_pca3D_by_','High',Group[2],' ','(',Group[1],' vs ',Group[2],' (Specific+Common)',')',sep='')
			plsda2D.plot <- paste('piRNA_plsda2D_by_','High',Group[2],' ','(',Group[1],' vs ',Group[2],' (Specific+Common)',')',sep='')
			plsda3D.plot <- paste('piRNA_plsda3D_by_','High',Group[2],' ','(',Group[1],' vs ',Group[2],' (Specific+Common)',')',sep='')
		} else if(n == 3) {
			mat1 = as.matrix(rbind(piRNA.normalizedData.Comp.UPinGroup1, piRNA.normalizedData.Comp.UPinGroup2)[, Group1.colidx.new, drop=FALSE])
			mat2 = as.matrix(rbind(piRNA.normalizedData.Comp.UPinGroup1, piRNA.normalizedData.Comp.UPinGroup2)[, setdiff(Group2.colidx.new,Group1.colidx.new), drop=FALSE])
			mat3 = as.matrix(rbind(piRNA.normalizedData.Comp.UPinGroup1, piRNA.normalizedData.Comp.UPinGroup2)[, unique(c(Group1.colidx.new,Group2.colidx.new)), drop=FALSE])
			
			heatmap.sep.plot <- paste('piRNA_heatmap_by_','High',Group[1],'+','High',Group[2],' ','(',Group[1],' + ',Group[2],' (Specific+Common)',')',sep='')
			heatmap.mix.plot <- paste('piRNA_heatmap_by_','High',Group[1],'+','High',Group[2],' ','(',Group[1],' vs ',Group[2],' (Specific+Common)',')',sep='')
			dendrogram.plot <- paste('piRNA_dendrogram_by_','High',Group[1],'+','High',Group[2],' ','(',Group[1],' vs ',Group[2],' (Specific+Common)',')',sep='')
			pca2D.plot <- paste('piRNA_pca2D_by_','High',Group[1],'+','High',Group[2],' ','(',Group[1],' vs ',Group[2],' (Specific+Common)',')',sep='')
			pca3D.plot <- paste('piRNA_pca3D_by_','High',Group[1],'+','High',Group[2],' ','(',Group[1],' vs ',Group[2],' (Specific+Common)',')',sep='')
			plsda2D.plot <- paste('piRNA_plsda2D_by_','High',Group[1],'+','High',Group[2],' ','(',Group[1],' vs ',Group[2],' (Specific+Common)',')',sep='')
			plsda3D.plot <- paste('piRNA_plsda3D_by_','High',Group[1],'+','High',Group[2],' ','(',Group[1],' vs ',Group[2],' (Specific+Common)',')',sep='')
		} else if(n == 4) {
			mat1 = as.matrix(piRNA.normalizedData.Comp.stat.FCdown[, Group1.colidx.new, drop=FALSE])
			mat2 = as.matrix(piRNA.normalizedData.Comp.stat.FCdown[, setdiff(Group2.colidx.new,Group1.colidx.new), drop=FALSE])
			mat3 = as.matrix(piRNA.normalizedData.Comp.stat.FCdown[, unique(c(Group1.colidx.new,Group2.colidx.new)), drop=FALSE])
			
			heatmap.sep.plot <- paste('piRNA_heatmap_by_','High',Group[1],' ','(',Group[1],' + ',Group[2],' (Common)',')',sep='')
			heatmap.mix.plot <- paste('piRNA_heatmap_by_','High',Group[1],' ','(',Group[1],' vs ',Group[2],' (Common)',')',sep='')
			dendrogram.plot <- paste('piRNA_dendrogram_by_','High',Group[1],' ','(',Group[1],' vs ',Group[2],' (Common)',')',sep='')
			pca2D.plot <- paste('piRNA_pca2D_by_','High',Group[1],' ','(',Group[1],' vs ',Group[2],' (Common)',')',sep='')
			pca3D.plot <- paste('piRNA_pca3D_by_','High',Group[1],' ','(',Group[1],' vs ',Group[2],' (Common)',')',sep='')
			plsda2D.plot <- paste('piRNA_plsda2D_by_','High',Group[1],' ','(',Group[1],' vs ',Group[2],' (Common)',')',sep='')
			plsda3D.plot <- paste('piRNA_plsda3D_by_','High',Group[1],' ','(',Group[1],' vs ',Group[2],' (Common)',')',sep='')
		} else if(n == 5) {
			mat1 = as.matrix(piRNA.normalizedData.Comp.stat.FCup[, Group1.colidx.new, drop=FALSE])
			mat2 = as.matrix(piRNA.normalizedData.Comp.stat.FCup[, setdiff(Group2.colidx.new,Group1.colidx.new), drop=FALSE])
			mat3 = as.matrix(piRNA.normalizedData.Comp.stat.FCup[, unique(c(Group1.colidx.new,Group2.colidx.new)), drop=FALSE])
			
			heatmap.sep.plot <- paste('piRNA_heatmap_by_','High',Group[2],' ','(',Group[1],' + ',Group[2],' (Common)',')',sep='')
			heatmap.mix.plot <- paste('piRNA_heatmap_by_','High',Group[2],' ','(',Group[1],' vs ',Group[2],' (Common)',')',sep='')
			dendrogram.plot <- paste('piRNA_dendrogram_by_','High',Group[2],' ','(',Group[1],' vs ',Group[2],' (Common)',')',sep='')
			pca2D.plot <- paste('piRNA_pca2D_by_','High',Group[2],' ','(',Group[1],' vs ',Group[2],' (Common)',')',sep='')
			pca3D.plot <- paste('piRNA_pca3D_by_','High',Group[2],' ','(',Group[1],' vs ',Group[2],' (Common)',')',sep='')
			plsda2D.plot <- paste('piRNA_plsda2D_by_','High',Group[2],' ','(',Group[1],' vs ',Group[2],' (Common)',')',sep='')
			plsda3D.plot <- paste('piRNA_plsda3D_by_','High',Group[2],' ','(',Group[1],' vs ',Group[2],' (Common)',')',sep='')
		} else if(n == 6) {
			mat1 = as.matrix(rbind(piRNA.normalizedData.Comp.stat.FCdown, piRNA.normalizedData.Comp.stat.FCup)[, Group1.colidx.new, drop=FALSE])
			mat2 = as.matrix(rbind(piRNA.normalizedData.Comp.stat.FCdown, piRNA.normalizedData.Comp.stat.FCup)[, setdiff(Group2.colidx.new,Group1.colidx.new), drop=FALSE])
			mat3 = as.matrix(rbind(piRNA.normalizedData.Comp.stat.FCdown, piRNA.normalizedData.Comp.stat.FCup)[, unique(c(Group1.colidx.new,Group2.colidx.new)), drop=FALSE])
			
			heatmap.sep.plot <- paste('piRNA_heatmap_by_','High',Group[1],'+','High',Group[2],' ','(',Group[1],' + ',Group[2],' (Common)',')',sep='')
			heatmap.mix.plot <- paste('piRNA_heatmap_by_','High',Group[1],'+','High',Group[2],' ','(',Group[1],' vs ',Group[2],' (Common)',')',sep='')
			dendrogram.plot <- paste('piRNA_dendrogram_by_','High',Group[1],'+','High',Group[2],' ','(',Group[1],' vs ',Group[2],' (Common)',')',sep='')
			pca2D.plot <- paste('piRNA_pca2D_by_','High',Group[1],'+','High',Group[2],' ','(',Group[1],' vs ',Group[2],' (Common)',')',sep='')
			pca3D.plot <- paste('piRNA_pca3D_by_','High',Group[1],'+','High',Group[2],' ','(',Group[1],' vs ',Group[2],' (Common)',')',sep='')
			plsda2D.plot <- paste('piRNA_plsda2D_by_','High',Group[1],'+','High',Group[2],' ','(',Group[1],' vs ',Group[2],' (Common)',')',sep='')
			plsda3D.plot <- paste('piRNA_plsda3D_by_','High',Group[1],'+','High',Group[2],' ','(',Group[1],' vs ',Group[2],' (Common)',')',sep='')
		}
		
		na.idx = which(is.na(mat1), arr.ind=TRUE)
		if(nrow(na.idx) > 0)
		{
			for(ri in 1:nrow(na.idx))
			{
				mat1[na.idx[ri,'row'],na.idx[ri,'col']] = 0
			}
		}
		na.idx = which(is.na(mat2), arr.ind=TRUE)
		if(nrow(na.idx) > 0)
		{
			for(ri in 1:nrow(na.idx))
			{
				mat2[na.idx[ri,'row'],na.idx[ri,'col']] = 0
			}
		}
		na.idx = which(is.na(mat3), arr.ind=TRUE)
		if(nrow(na.idx) > 0)
		{
			for(ri in 1:nrow(na.idx))
			{
				mat3[na.idx[ri,'row'],na.idx[ri,'col']] = 0
			}
		}
		
		# --- Heatmap -------------------------------------------------------------------------------------
		gc()	# Garbage Collection (for release memory)
		
		mainDir <- paste('./Data/02_significant/piRNA/',Group[1],' vs ',Group[2],sep='')
		subDir <- 'Heatmap'
		if(file.exists(file.path(mainDir, subDir))) {
		} else {
			dir.create(file.path(mainDir, subDir))
		}
		
		# https://jokergoo.github.io/ComplexHeatmap-reference/book/
		
		#install_github("jokergoo/ComplexHeatmap")	# ComplexHeatmap: Make Complex Heatmaps
		suppressPackageStartupMessages(library(ComplexHeatmap))
		
		Group1.colorname = Group[1]
		if(Group[1] == 'PDMCI') {
			Group1.colorname = 'PD-MCI'
		} else if(Group[1] == 'PD') {
			if(Group[2] == 'PDMCI') {
				Group1.colorname = setdiff(c('PDND','PDMCI','PDD'), Group[2])
			} else {
				Group1.colorname = setdiff(c('PDND','PD-MCI','PDD'), Group[2])
			}
		}
		Group1.coloridx = sapply(Group1.colorname, FUN=function(X) which(names(diseaseStage.color)==X))
		
		Group2.colorname = Group[2]
		if(Group[2] == 'PDMCI') {
			Group2.colorname = 'PD-MCI'
		} else if(Group[2] == 'PD') {
			if(Group[1] == 'PDMCI') {
				Group2.colorname = setdiff(c('PDND','PDMCI','PDD'), Group[1])
			} else {
				Group2.colorname = setdiff(c('PDND','PD-MCI','PDD'), Group[1])
			}
		}
		Group2.coloridx = sapply(Group2.colorname, FUN=function(X) which(names(diseaseStage.color)==X))
		
		# /// separate heatmap ///
		if(nrow(mat1) > 0) {
			column_ha = HeatmapAnnotation("Symp Group" = diseaseStage.Group1[which((diseaseStage.Group1[,'Symp Group'] %in% Group1.colorname) & (diseaseStage.Group1[,'No.']!='')),'Symp Group'], 
										  Gender = diseaseStage.Group1[which((diseaseStage.Group1[,'Symp Group'] %in% Group1.colorname) & (diseaseStage.Group1[,'No.']!='')),'Gender'],
										  show_annotation_name = FALSE, 
										  annotation_legend_param = list(
										  	"Symp Group" = list(
										  		title = "Symp Group", 
										  		at = Group1.colorname, 
										  		labels = Group1.colorname
										  	), 
										  	Gender = list(
										  		title = "Gender", 
										  		at = c("M", "F"), 
										  		labels = c("Male", "Female")
										  	)
										  ), 
										  col = list(
										  	"Symp Group" = diseaseStage.color[Group1.coloridx], 
										  	Gender = c("M" = "royalblue1", "F" = "palevioletred1")
										  ))
			ht1 = Heatmap(mat1, 
						  name = Group[1], 
						  column_title = Group[1], 
						  column_title_gp = gpar(fontsize = 20, fontface = "bold"), 
						  
						  cluster_rows = TRUE, 
						  clustering_distance_rows = "euclidean", 
						  clustering_method_rows = "ward.D2", 
						  show_row_dend = TRUE, 
						  row_dend_reorder = TRUE, 
						  
						  cluster_columns = TRUE, 
						  clustering_distance_columns = "euclidean", 
						  clustering_method_columns = "ward.D2", 
						  show_column_dend = TRUE, 
						  column_dend_reorder = TRUE, 
						  
						  top_annotation = column_ha)
			column_ha = HeatmapAnnotation("Symp Group" = diseaseStage.Group2[which((diseaseStage.Group2[,'Symp Group'] %in% Group2.colorname) & (diseaseStage.Group2[,'No.']!='')),'Symp Group'], 
										  Gender = diseaseStage.Group2[which((diseaseStage.Group2[,'Symp Group'] %in% Group2.colorname) & (diseaseStage.Group2[,'No.']!='')),'Gender'],
										  annotation_legend_param = list(
										  	"Symp Group" = list(
										  		title = "Symp Group", 
										  		at = Group2.colorname, 
										  		labels = Group2.colorname
										  	), 
										  	Gender = list(
										  		title = "Gender", 
										  		at = c("M", "F"), 
										  		labels = c("Male", "Female")
										  	)
										  ), 
										  col = list(
										  	"Symp Group" = diseaseStage.color[Group2.coloridx], 
										  	Gender = c("M" = "royalblue1", "F" = "palevioletred1")
										  ))
			ht2 = Heatmap(mat2, 
						  name = Group[2], 
						  column_title = Group[2], 
						  column_title_gp = gpar(fontsize = 20, fontface = "bold"), 
						  
						  cluster_rows = TRUE, 
						  clustering_distance_rows = "euclidean", 
						  clustering_method_rows = "ward.D2", 
						  show_row_dend = TRUE, 
						  row_dend_reorder = TRUE, 
						  
						  cluster_columns = TRUE, 
						  clustering_distance_columns = "euclidean", 
						  clustering_method_columns = "ward.D2", 
						  show_column_dend = TRUE, 
						  column_dend_reorder = TRUE, 
						  
						  top_annotation = column_ha)
			if(ncol(mat1)+ncol(mat2) > 0 & ncol(mat1)+ncol(mat2) <= 10) {			width.inch = 9
			} else if(ncol(mat1)+ncol(mat2) > 10 & ncol(mat1)+ncol(mat2) <= 20) {	width.inch = 10
			} else if(ncol(mat1)+ncol(mat2) > 20 & ncol(mat1)+ncol(mat2) <= 30) {	width.inch = 11
			} else if(ncol(mat1)+ncol(mat2) > 30 & ncol(mat1)+ncol(mat2) <= 40) {	width.inch = 12
			} else if(ncol(mat1)+ncol(mat2) > 40 & ncol(mat1)+ncol(mat2) <= 50) {	width.inch = 13
			} else if(ncol(mat1)+ncol(mat2) > 50 & ncol(mat1)+ncol(mat2) <= 60) {	width.inch = 14
			} else if(ncol(mat1)+ncol(mat2) > 60) {									width.inch = 15
			}
			if(nrow(mat1) > 0 & nrow(mat1) <= 3) {				height.inch = 3
			} else if(nrow(mat1) > 3 & nrow(mat1) <= 6) {		height.inch = 4
			} else if(nrow(mat1) > 6 & nrow(mat1) <= 12) {		height.inch = 5
			} else if(nrow(mat1) > 12 & nrow(mat1) <= 15) {		height.inch = 6
			} else if(nrow(mat1) > 15 & nrow(mat1) <= 18) {		height.inch = 7
			} else if(nrow(mat1) > 18 & nrow(mat1) <= 24) {		height.inch = 8
			} else if(nrow(mat1) > 24 & nrow(mat1) <= 28) {		height.inch = 9
			} else if(nrow(mat1) > 28 & nrow(mat1) <= 32) {		height.inch = 10
			} else if(nrow(mat1) > 32 & nrow(mat1) <= 36) {		height.inch = 11
			} else if(nrow(mat1) > 36 & nrow(mat1) <= 40) {		height.inch = 12
			} else if(nrow(mat1) > 40 & nrow(mat1) <= 50) {		height.inch = 13
			} else if(nrow(mat1) > 50 & nrow(mat1) <= 60) {		height.inch = 14
			} else if(nrow(mat1) > 60 & nrow(mat1) <= 70) {		height.inch = 15
			} else if(nrow(mat1) > 70 & nrow(mat1) <= 80) {		height.inch = 16
			} else if(nrow(mat1) > 80 & nrow(mat1) <= 90) {		height.inch = 18
			} else if(nrow(mat1) > 90 & nrow(mat1) <= 100) {	height.inch = 20
			} else if(nrow(mat1) > 100 & nrow(mat1) <= 110) {	height.inch = 22
			} else if(nrow(mat1) > 110 & nrow(mat1) <= 120) {	height.inch = 24
			} else if(nrow(mat1) > 120) {						height.inch = 25
			}
			png(paste('./Data/02_significant/piRNA/',Group[1],' vs ',Group[2],'/','Heatmap','/',heatmap.sep.plot,'.png',sep=''), width=width.inch, height=height.inch, units="in", res=600)
			ht_list = ht1 + ht2
			draw(ht_list, heatmap_legend_side = "left", annotation_legend_side = "right")
			dev.off()
		}
		
		# /// mix heatmap ///
		if(nrow(mat3) > 0) {
			column_ha = HeatmapAnnotation("Symp Group" = diseaseStage.Group1_Group2[which((diseaseStage.Group1_Group2[,'Symp Group'] %in% c(Group1.colorname, Group2.colorname)) & (diseaseStage.Group1_Group2[,'No.']!='')),'Symp Group'], 
										  Gender = diseaseStage.Group1_Group2[which((diseaseStage.Group1_Group2[,'Symp Group'] %in% c(Group1.colorname, Group2.colorname)) & (diseaseStage.Group1_Group2[,'No.']!='')),'Gender'],
										  annotation_legend_param = list(
										  	"Symp Group" = list(
										  		title = "Symp Group", 
										  		at = c(Group1.colorname, Group2.colorname), 
										  		labels = c(Group1.colorname, Group2.colorname)
										  	), 
										  	Gender = list(
										  		title = "Gender", 
										  		at = c("M", "F"), 
										  		labels = c("Male", "Female")
										  	)
										  ), 
										  col = list(
										  	"Symp Group" = diseaseStage.color[c(Group1.coloridx, Group2.coloridx)], 
										  	Gender = c("M" = "royalblue1", "F" = "palevioletred1")
										  ))
			if(ncol(mat3) > 0 & ncol(mat3) <= 10) {				width.inch = 9
			} else if(ncol(mat3) > 10 & ncol(mat3) <= 20) {		width.inch = 10
			} else if(ncol(mat3) > 20 & ncol(mat3) <= 30) {		width.inch = 11
			} else if(ncol(mat3) > 30 & ncol(mat3) <= 40) {		width.inch = 12
			} else if(ncol(mat3) > 40 & ncol(mat3) <= 50) {		width.inch = 13
			} else if(ncol(mat3) > 50 & ncol(mat3) <= 60) {		width.inch = 14
			} else if(ncol(mat3) > 60) {						width.inch = 15
			}
			if(nrow(mat3) > 0 & nrow(mat3) <= 3) {				height.inch = 3
			} else if(nrow(mat3) > 3 & nrow(mat3) <= 6) {		height.inch = 4
			} else if(nrow(mat3) > 6 & nrow(mat3) <= 12) {		height.inch = 5
			} else if(nrow(mat3) > 12 & nrow(mat3) <= 15) {		height.inch = 6
			} else if(nrow(mat3) > 15 & nrow(mat3) <= 18) {		height.inch = 7
			} else if(nrow(mat3) > 18 & nrow(mat3) <= 24) {		height.inch = 8
			} else if(nrow(mat3) > 24 & nrow(mat3) <= 28) {		height.inch = 9
			} else if(nrow(mat3) > 28 & nrow(mat3) <= 32) {		height.inch = 10
			} else if(nrow(mat3) > 32 & nrow(mat3) <= 36) {		height.inch = 11
			} else if(nrow(mat3) > 36 & nrow(mat3) <= 40) {		height.inch = 12
			} else if(nrow(mat3) > 40 & nrow(mat3) <= 50) {		height.inch = 13
			} else if(nrow(mat3) > 50 & nrow(mat3) <= 60) {		height.inch = 14
			} else if(nrow(mat3) > 60 & nrow(mat3) <= 70) {		height.inch = 15
			} else if(nrow(mat3) > 70 & nrow(mat3) <= 80) {		height.inch = 16
			} else if(nrow(mat3) > 80 & nrow(mat3) <= 90) {		height.inch = 18
			} else if(nrow(mat3) > 90 & nrow(mat3) <= 100) {	height.inch = 20
			} else if(nrow(mat3) > 100 & nrow(mat3) <= 110) {	height.inch = 22
			} else if(nrow(mat3) > 110 & nrow(mat3) <= 120) {	height.inch = 24
			} else if(nrow(mat3) > 120) {						height.inch = 25
			}
			png(paste('./Data/02_significant/piRNA/',Group[1],' vs ',Group[2],'/','Heatmap','/',heatmap.mix.plot,'.png',sep=''), width=width.inch, height=height.inch, units="in", res=600)
			ht3 = Heatmap(mat3, 
						  name = "normalized", 
						  column_title = paste(Group[1],' vs ',Group[2],sep=''), 
						  column_title_gp = gpar(fontsize = 20, fontface = "bold"), 
						  
						  cluster_rows = TRUE, 
						  clustering_distance_rows = "euclidean", 
						  clustering_method_rows = "ward.D2", 
						  show_row_dend = TRUE, 
						  row_dend_reorder = TRUE, 
						  
						  cluster_columns = TRUE, 
						  clustering_distance_columns = "euclidean", 
						  clustering_method_columns = "ward.D2", 
						  show_column_dend = TRUE, 
						  column_dend_reorder = TRUE, 
						  
						  top_annotation = column_ha)
			draw(ht3, heatmap_legend_side = "left", annotation_legend_side = "right")
			dev.off()
		}
		
		# data format
		mat = mat3
		
		if(nrow(mat) > 0) {
			samples_serial_num = sapply(strsplit(colnames(mat), '_'), FUN=function(X) X[1])
			samples_id = sapply(strsplit(colnames(mat), '_'), FUN=function(X) X[2])
			
			Group1.colidx.new = unlist(sapply(diseaseStage.ALL[eval(parse(text=paste('diseaseStage.', Group[1], '.idx', sep=''))),'Sample ID'], FUN=function(X) which(samples_id %in% X)))
			Group2.colidx.new = unlist(sapply(diseaseStage.ALL[eval(parse(text=paste('diseaseStage.', Group[2], '.idx', sep=''))),'Sample ID'], FUN=function(X) which(samples_id %in% X)))
			
			sampleValue = mat
			sampleGroup = sapply(samples_id, FUN=function(X) diseaseStage.ALL[which(diseaseStage.ALL[,'Sample ID'] %in% X),'Symp Group'])
			sampleGroup[which(sampleGroup == 'PD-MCI')] = 'PDMCI'
			
			sampleValueGroup = as.data.frame(t(rbind(sampleValue, sampleGroup)), stringsAsFactors=FALSE)
			colnames(sampleValueGroup)[ncol(sampleValueGroup)] = c('Group')
			sampleValueGroup[, -which(colnames(sampleValueGroup) == 'Group')] = apply(sampleValueGroup[, -which(colnames(sampleValueGroup) == 'Group'), drop=FALSE], 2, as.numeric)
			sampleValueGroup[,'Group'] = as.factor(sampleValueGroup[,'Group'])
			
			if(Group[1] == 'PD' | Group[2] == 'PD') {
				PD.subtype.idx = c(which(levels(sampleValueGroup[,'Group'])=='PDND'), which(levels(sampleValueGroup[,'Group'])=='PDMCI'), which(levels(sampleValueGroup[,'Group'])=='PDD'))
				if(Group[1] == 'PD') {
					if(Group[2] == 'HC') {
						sampleValueGroup[,'Group'] = factor(sampleValueGroup[,'Group'], levels(sampleValueGroup[,'Group'])[c(PD.subtype.idx, which(levels(sampleValueGroup[,'Group'])=='HC'))])
					} else if(Group[2] == 'MSA') {
						sampleValueGroup[,'Group'] = factor(sampleValueGroup[,'Group'], levels(sampleValueGroup[,'Group'])[c(PD.subtype.idx, which(levels(sampleValueGroup[,'Group'])=='MSA'))])
					} else if(Group[2] == 'PDND' | Group[2] == 'PDMCI' | Group[2] == 'PDD') {
						sampleValueGroup[,'Group'] = factor(sampleValueGroup[,'Group'], levels(sampleValueGroup[,'Group'])[PD.subtype.idx])
					}
				}
				if(Group[2] == 'PD') {
					if(Group[1] == 'HC') {
						sampleValueGroup[,'Group'] = factor(sampleValueGroup[,'Group'], levels(sampleValueGroup[,'Group'])[c(which(levels(sampleValueGroup[,'Group'])=='HC'), PD.subtype.idx)])
					} else if(Group[1] == 'MSA') {
						sampleValueGroup[,'Group'] = factor(sampleValueGroup[,'Group'], levels(sampleValueGroup[,'Group'])[c(which(levels(sampleValueGroup[,'Group'])=='MSA'), PD.subtype.idx)])
					} else if(Group[1] == 'PDND' | Group[1] == 'PDMCI' | Group[1] == 'PDD') {
						sampleValueGroup[,'Group'] = factor(sampleValueGroup[,'Group'], levels(sampleValueGroup[,'Group'])[PD.subtype.idx])
					}
				}
			} else {
				sampleValueGroup[,'Group'] = factor(sampleValueGroup[,'Group'], levels(sampleValueGroup[,'Group'])[c(which(levels(sampleValueGroup[,'Group'])==Group[1]), 
																													 which(levels(sampleValueGroup[,'Group'])==Group[2]))])
			}
			
			# --- Hierarchical Clustering ---------------------------------------------------------------------
			gc()	# Garbage Collection (for release memory)
			
			mainDir <- paste('./Data/02_significant/piRNA/',Group[1],' vs ',Group[2],sep='')
			subDir <- 'Dendrogram'
			if(file.exists(file.path(mainDir, subDir))) {
			} else {
				dir.create(file.path(mainDir, subDir))
			}
			
			# https://www.jamleecute.com/hierarchical-clustering-%E9%9A%8E%E5%B1%A4%E5%BC%8F%E5%88%86%E7%BE%A4/
			# http://rpubs.com/skydome20/R-Note9-Clustering
			# https://cran.r-project.org/web/packages/dendextend/vignettes/introduction.html
			# https://cran.r-project.org/web/packages/dendextend/vignettes/FAQ.html
			
			#install.packages("factoextra")		# factoextra: Extract and Visualize the Results of Multivariate Data Analyses
			#install.packages("dendextend")		# dendextend: Extending 'dendrogram' Functionality in R
			#suppressPackageStartupMessages(library(factoextra))
			suppressPackageStartupMessages(library(dendextend))
			
			df <- sampleValueGroup[, -which(colnames(sampleValueGroup)=='Group'), drop=FALSE]
			
			E.dist <- dist(df, method = "euclidean")
			h.E.Ward.cluster <- hclust(E.dist, method = "ward.D2")
			
			# Determining optimal number of clusters
			#fviz_nbclust(df, FUN = hcut, method = "wss") + geom_vline(xintercept = 2, linetype = 2)	# Elbow Method
			#fviz_nbclust(df, FUN = hcut, method = "silhouette")										# Average Silhouette Method
			
			#cut.h.cluster <- cutree(h.E.Ward.cluster, k = 2)
			#table(cut.h.cluster, sampleValueGroup[,'Group'])
			
			#fviz_cluster(list(data = df, cluster = cut.h.cluster))
			
			dend <- as.dendrogram(h.E.Ward.cluster)
			
			diseaseStage.darkcolor <- c('darkgreen','indianred','red','darkred','black')
			names(diseaseStage.darkcolor) = diseaseStage
			names(diseaseStage.darkcolor)[which(names(diseaseStage.darkcolor)=='PD-MCI')] = 'PDMCI'
			
			colors_to_use <- sapply(sampleValueGroup[,'Group'], FUN=function(X) diseaseStage.darkcolor[which(names(diseaseStage.darkcolor) %in% X)])
			
			X = cbind(Gender = sapply(samples_id, FUN=function(X) diseaseStage.ALL[which(diseaseStage.ALL[,'Sample ID'] %in% X),'Gender']), colors_to_use)
			rownames(X) = rownames(sampleValueGroup)
			X[which(X[,'Gender']=="M"),'Gender'] = "royalblue1"
			X[which(X[,'Gender']=="F"),'Gender'] = "palevioletred1"
			
			colors_to_use <- colors_to_use[order.dendrogram(dend)]
			
			labels_colors(dend) <- colors_to_use
			#dend <- assign_values_to_leaves_edgePar(dend=dend, value = colors_to_use, edgePar = "col")
			#dend <- dend %>% set("leaves_pch", 19) %>% set("leaves_col", colors_to_use)
			
			if(ncol(mat) > 0 & ncol(mat) <= 10) {				width.inch = 9
			} else if(ncol(mat) > 10 & ncol(mat) <= 20) {		width.inch = 10
			} else if(ncol(mat) > 20 & ncol(mat) <= 30) {		width.inch = 11
			} else if(ncol(mat) > 30 & ncol(mat) <= 40) {		width.inch = 12
			} else if(ncol(mat) > 40 & ncol(mat) <= 50) {		width.inch = 13
			} else if(ncol(mat) > 50 & ncol(mat) <= 60) {		width.inch = 14
			} else if(ncol(mat) > 60) {							width.inch = 15
			}
			
			png(paste('./Data/02_significant/piRNA/',Group[1],' vs ',Group[2],'/','Dendrogram','/',dendrogram.plot,'.png',sep=''), width=width.inch, height=6, units="in", res=600)
			
			par(mar = c(8, 4, 2, 5))
			par(xpd = TRUE)
			
			plot(dend, main = paste(Group[1],' vs ',Group[2],sep=''))
			colored_bars(colors = X, dend = dend, rowLabels = c("Gender","Symp Group"))
			
			legend("topright", bty = "n", inset = c(-0.075,0), legend = levels(sampleValueGroup[,'Group']), fill = sapply(levels(sampleValueGroup[,'Group']), FUN=function(X) diseaseStage.darkcolor[which(names(diseaseStage.darkcolor) %in% X)]))
			names(diseaseStage.darkcolor)[which(names(diseaseStage.darkcolor)=='PDMCI')] = 'PD-MCI'
			
			dev.off()
			
			# --- PCA (Principal Component Analysis) ----------------------------------------------------------
			gc()	# Garbage Collection (for release memory)
			
			mainDir <- paste('./Data/02_significant/piRNA/',Group[1],' vs ',Group[2],sep='')
			subDir <- 'PCA'
			if(file.exists(file.path(mainDir, subDir))) {
			} else {
				dir.create(file.path(mainDir, subDir))
			}
			
			# https://cran.r-project.org/web/packages/ggfortify/vignettes/plot_pca.html
			
			#install.packages("ggfortify")		# ggfortify: Data Visualization Tools for Statistical Analysis Results
			#install.packages("factoextra")		# factoextra: Extract and Visualize the Results of Multivariate Data Analyses
			suppressPackageStartupMessages(library(ggfortify))
			suppressPackageStartupMessages(library(factoextra))
			
			# https://www.bioconductor.org/packages/release/bioc/vignettes/mixOmics/inst/doc/vignette.html
			# http://blog.sciencenet.cn/blog-3406804-1163309.html
			# http://mixomics.org/graphics/sample-plots/plotindiv/
			
			if (!requireNamespace("BiocManager", quietly = TRUE))
				install.packages("BiocManager")	# BiocManager: Installing and Managing Bioconductor Packages
			#BiocManager::install("mixOmics")	# mixOmics: Omics Data Integration Project
			suppressPackageStartupMessages(library(mixOmics))
			
			#install.packages("ggplot2")		# ggplot2: Create Elegant Data Visualisations Using the Grammar of Graphics
			#install.packages("rgl")			# rgl: 3D Visualization Using OpenGL
			suppressPackageStartupMessages(library(ggplot2))
			library(rgl)
			
			df <- sampleValueGroup[, -which(colnames(sampleValueGroup)=='Group'), drop=FALSE]
			group <- sampleValueGroup[,'Group']
			
			if(ncol(df) >= 2) {
				# /// prcomp {stats} ///
				res.pca <- prcomp(df)					# prcomp {stats}
				#summary(res.pca)
				#str(res.pca)
				
				#png(paste('./Data/02_significant/piRNA/',Group[1],' vs ',Group[2],'/','PCA','/',pca2D.plot,'.png',sep=''), width=10, height=10, units="in", res=600)
				
				pic <- autoplot(res.pca, data = sampleValueGroup, colour = 'Group', shape = TRUE, label.size = 3)
				pic <- pic + 
					labs(title = paste(Group[1],' vs ',Group[2],sep=''), fill = "Symp Group", color = "Symp Group") + 
					theme(axis.text.x = element_text(color = "black", size = 24, margin = margin(10,0,0,0)), 
						  axis.text.y = element_text(color = "black", size = 24, margin = margin(0,10,0,0)), 
						  axis.title.x = element_text(color = "black", size = 26, margin = margin(30,0,0,0)), 
						  axis.title.y = element_text(color = "black", size = 26, margin = margin(0,30,0,0)), 
						  plot.title = element_text(color = "black", size = 28, margin = margin(0,0,30,0), hjust = 0.5), 
						  legend.text = element_text(color = "black", size = 18), 
						  legend.title = element_text(color = "black", size = 18))
				#print(pic)
				
				#dev.off()
				
				#ggsave(paste('./Data/02_significant/piRNA/',Group[1],' vs ',Group[2],'/','PCA','/',pca2D.plot,'.png',sep=''), pic, width=10, height=10, units="in", dpi=600, bg="transparent")
				
				# https://rpkgs.datanovia.com/factoextra/reference/fviz_pca.html
				#fviz_pca_ind(res.pca, geom.ind = "point", repel = TRUE, col.ind = group, legend.title = "Symp Group", ggtheme = theme_minimal())
				#fviz_pca_var(res.pca, col.var = "contrib", gradient.cols = c("white","blue","red"), ggtheme = theme_minimal())
				#fviz_pca_biplot(res.pca, label = "var", habillage = group, addEllipses = TRUE, ellipse.level = 0.95, ggtheme = theme_minimal())
				
				# http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/119-pca-in-r-using-ade4-quick-scripts/
				#fviz_eig(res.pca)
				#fviz_contrib(res.pca, choice = "var", axes = 1, top = 10)	# Contributions of variables to PC1
				#fviz_contrib(res.pca, choice = "var", axes = 2, top = 10)	# Contributions of variables to PC2
				#fviz_contrib(res.pca, choice = "var", axes = 1:2, top = 10)
				
				# /// pca {mixOmics} ///
				if(ncol(df) == 2) {
					pca.res <- mixOmics::pca(df, ncomp = 2)	# pca {mixOmics}
				} else if(ncol(df) >= 3) {
					pca.res <- mixOmics::pca(df, ncomp = 3)	# pca {mixOmics}
				}
				#summary(pca.res)
				#str(pca.res)
				
				png(paste('./Data/02_significant/piRNA/',Group[1],' vs ',Group[2],'/','PCA','/',pca2D.plot,'.png',sep=''), width=10, height=10, units="in", res=600)
				plotIndiv(pca.res, group = group, ind.names = FALSE, style = 'ggplot2', legend = TRUE, legend.title = "Symp Group", cex = 3, size.title = 28, size.xlabel = 26, size.ylabel = 26, size.axis = 24, size.legend.title = 18, size.legend = 18, title = paste(Group[1],' vs ',Group[2],sep=''))
				#plotIndiv(pca.res, group = group, ind.names = FALSE, style = 'graphics', legend = TRUE, legend.title = "Symp Group", cex = 2, size.title = 3, size.xlabel = 2, size.ylabel = 2, size.axis = 1.5, size.legend.title = 1.5, size.legend = 1.5, title = paste(Group[1],' vs ',Group[2],sep=''))
				dev.off()
				
				if(ncol(df) >= 3) {
					plotIndiv(pca.res, group = group, ind.names = TRUE, style = '3d', legend = TRUE, legend.title = "Symp Group", size.legend.title = 1.5, size.legend = 1.5, title = paste(Group[1],' vs ',Group[2],sep=''))
					snapshot3d(paste('./Data/02_significant/piRNA/',Group[1],' vs ',Group[2],'/','PCA','/',pca3D.plot,'.png',sep=''))
					rgl.close()
				}
			}
			
			# --- PLS-DA (Partial Least Squares Discriminant Analysis) ----------------------------------------
			gc()	# Garbage Collection (for release memory)
			
			mainDir <- paste('./Data/02_significant/piRNA/',Group[1],' vs ',Group[2],sep='')
			subDir <- 'PLSDA'
			if(file.exists(file.path(mainDir, subDir))) {
			} else {
				dir.create(file.path(mainDir, subDir))
			}
			
			# https://www.bioconductor.org/packages/release/bioc/vignettes/mixOmics/inst/doc/vignette.html
			# http://blog.sciencenet.cn/blog-3406804-1163309.html
			# http://mixomics.org/graphics/sample-plots/plotindiv/
			
			if (!requireNamespace("BiocManager", quietly = TRUE))
				install.packages("BiocManager")	# BiocManager: Installing and Managing Bioconductor Packages
			#BiocManager::install("mixOmics")	# mixOmics: Omics Data Integration Project
			suppressPackageStartupMessages(library(mixOmics))
			
			#install.packages("ggplot2")		# ggplot2: Create Elegant Data Visualisations Using the Grammar of Graphics
			#install.packages("rgl")			# rgl: 3D Visualization Using OpenGL
			suppressPackageStartupMessages(library(ggplot2))
			library(rgl)
			
			df <- sampleValueGroup[, -which(colnames(sampleValueGroup)=='Group'), drop=FALSE]
			group <- sampleValueGroup[,'Group']
			
			if(ncol(df) >= 3) {
				# /// plsda {mixOmics} ///
				plsda_result <- mixOmics::plsda(df, group, ncomp = 3)	# plsda {mixOmics}
				#names(plsda_result)
				
				#plsda_result$explained_variance$X
				#plsda_result$variates$X
				#plsda_result$loadings$X
				plsda_result_eig <- {plsda_result$explained_variance$X}[1:2]
				sample_site <- data.frame(plsda_result$variates)[1:2]
				sample_site$names <- rownames(sample_site)
				names(sample_site)[1:2] <- c('plsda1', 'plsda2')
				#write.table(sample_site, 'plsda_sample.txt', row.names = FALSE, sep = '\t', quote = FALSE)
				
				png(paste('./Data/02_significant/piRNA/',Group[1],' vs ',Group[2],'/','PLSDA','/',plsda2D.plot,'.png',sep=''), width=10, height=10, units="in", res=600)
				plotIndiv(plsda_result, ind.names = FALSE, style = 'ggplot2', legend = TRUE, legend.title = "Symp Group", cex = 3, size.title = 28, size.xlabel = 26, size.ylabel = 26, size.axis = 24, size.legend.title = 18, size.legend = 18, title = paste(Group[1],' vs ',Group[2],sep=''))
				#plotIndiv(plsda_result, ind.names = FALSE, style = 'graphics', legend = TRUE, legend.title = "Symp Group", cex = 2, size.title = 3, size.xlabel = 2, size.ylabel = 2, size.axis = 1.5, size.legend.title = 1.5, size.legend = 1.5, title = paste(Group[1],' vs ',Group[2],sep=''))
				dev.off()
				
				plotIndiv(plsda_result, ind.names = TRUE, style = '3d', legend = TRUE, legend.title = "Symp Group", size.legend.title = 1.5, size.legend = 1.5, title = paste(Group[1],' vs ',Group[2],sep=''))
				snapshot3d(paste('./Data/02_significant/piRNA/',Group[1],' vs ',Group[2],'/','PLSDA','/',plsda3D.plot,'.png',sep=''))
				rgl.close()
				
				plsda_plot <- ggplot(sample_site, aes(plsda1, plsda2, color = group, label = names)) + 
					geom_point(size = 1.5, alpha = 0.6) + 
					stat_ellipse(show.legend = FALSE) + 
					scale_color_manual(values = c('#1D7ACC', '#F67433', '#00815F')) + 
					theme(panel.grid = element_line(color = 'grey50'), panel.background = element_rect(color = 'black', fill = 'transparent')) + 
					theme(legend.title = element_blank(), legend.key = element_rect(fill = 'transparent')) + 
					labs(x = paste('PLS-DA axis1 ( explained variance ', round(100 * plsda_result_eig[1], 2), '% )', sep = ''), y = paste('PLS-DA axis2 ( explained variance ', round(100 * plsda_result_eig[2], 2), '% )', sep = ''))
				#ggsave('plsda_plot.pdf', plsda_plot, width = 6, height = 5)
				#ggsave('plsda_plot.png', plsda_plot, width = 6, height = 5)
			}
		}
	}
	
	# --- write xlsx file by openxlsx package ---------------------------------------------------------
	wb <- openxlsx::createWorkbook(paste('./Data/02_significant/piRNA/',Group[1],' vs ',Group[2],'/','piRNA_datatable ','(',Group[1],' vs ',Group[2],')','.xlsx',sep=''))
	modifyBaseFont(wb, fontSize=12, fontColour="black", fontName="Times New Roman")
	
	for(s in 1:10)
	{
		if(s == 1) {
			sheet = paste('sample information','(',Group[1],',',Group[2],')',sep='')
			sheetData = diseaseStage.Group1_Group2
		} else if(s == 2) {
			sheet = paste('piRNA(normalized)','(',Group[1],',',Group[2],')',sep='')
			sheetData = piRNA.normalizedData.Comp
		} else if(s == 3) {
			sheet = paste('Specific piRNA in ',Group[1],sep='')
			sheetData = piRNA.normalizedData.Comp[Group1.piRNA.rowidx, , drop=FALSE]
		} else if(s == 4) {
			sheet = paste('Specific piRNA in ',Group[2],sep='')
			sheetData = piRNA.normalizedData.Comp[Group2.piRNA.rowidx, , drop=FALSE]
		} else if(s == 5) {
			sheet = paste('Common piRNA in ',Group[1],' & ',Group[2],sep='')
			sheetData = piRNA.normalizedData.Comp.stat
			sheetData = sheetData[order(sheetData[,'pvalue']), , drop=FALSE]
		} else if(s == 6) {
			sheet = paste('pvalue',' < ',pvalue.cutoff,sep='')
			sheetData = piRNA.normalizedData.Comp.stat[which(piRNA.normalizedData.Comp.stat[,'pvalue'] < pvalue.cutoff), , drop=FALSE]
			sheetData = sheetData[order(sheetData[,'pvalue']), , drop=FALSE]
		} else if(s == 7) {
			sheet = paste('pvalue',' >= ',pvalue.cutoff,sep='')
			sheetData = piRNA.normalizedData.Comp.stat[which(piRNA.normalizedData.Comp.stat[,'pvalue'] >= pvalue.cutoff), , drop=FALSE]
			sheetData = sheetData[order(sheetData[,'pvalue']), , drop=FALSE]
		} else if(s == 8) {
			sheet = paste(Group[1],' > ',Group[2],sep='')
			sheetData = piRNA.normalizedData.Comp.UPinGroup1
		} else if(s == 9) {
			sheet = paste(Group[1],' < ',Group[2],sep='')
			sheetData = piRNA.normalizedData.Comp.UPinGroup2
		} else if(s == 10) {
			sheet = 'Not available piRNA'
			sheetData = piRNA.normalizedData.Comp[NA.piRNA.rowidx, , drop=FALSE]
		}
		
		if(!is.null(sheetData))
		{
			if(nrow(sheetData) > 0)
			{
				#na.idx = which(is.na(sheetData), arr.ind=TRUE)
				#if(nrow(na.idx) > 0)
				#{
				#	for(ri in 1:nrow(na.idx))
				#	{
				#		sheetData[na.idx[ri,'row'],na.idx[ri,'col']] = ''
				#	}
				#}
				
				addWorksheet(wb, sheet)
				writeData(wb, sheet, sheetData, colNames=TRUE, rowNames=FALSE, keepNA=FALSE)
				
				setRowHeights(wb, sheet, rows=c(1:(1+nrow(sheetData))), heights=16.5)
				if(length(grep("sample information", sheet)) == 0) {
					if(length(grep("common", sheet, ignore.case=TRUE)) > 0 | length(grep("pvalue", sheet, ignore.case=TRUE)) > 0 | length(grep(">", sheet, ignore.case=TRUE)) > 0 | length(grep("<", sheet, ignore.case=TRUE)) > 0) {
						setColWidths(wb, sheet, cols=c(1:ncol(sheetData)), widths=c(10,30,17.5,rep(12.5,ncol(sheetData)-3)))
						freezePane(wb, sheet, firstActiveRow=2, firstActiveCol=4+4)
					} else {
						setColWidths(wb, sheet, cols=c(1:ncol(sheetData)), widths=c(10,30,17.5,rep(12.5,ncol(sheetData)-3)))
						freezePane(wb, sheet, firstActiveRow=2, firstActiveCol=4)
					}
				} else if(length(grep("sample information", sheet)) > 0) {
					setColWidths(wb, sheet, cols=c(1:ncol(sheetData)), widths=c(7.5,rep(12.5,ncol(sheetData)-4),45,15,12.5))
					freezePane(wb, sheet, firstActiveRow=2)
				}
				
				style <- createStyle(halign="left", valign="center")
				addStyle(wb, sheet, style, rows=c(1:(1+nrow(sheetData))), cols=c(1:ncol(sheetData)), gridExpand=TRUE)
				style <- createStyle(textDecoration="bold")
				addStyle(wb, sheet, style, rows=1, cols=c(1:ncol(sheetData)), stack=TRUE)
				
				if(length(grep("sample information", sheet)) == 0) {
					right_align.idx = which(! colnames(sheetData) %in% c('type','Name','Accession_ID','Drives_from','locus_GRCh38','locus_hg19'))
					if(length(right_align.idx) > 0) {
						style <- createStyle(halign="right", valign="center")
						addStyle(wb, sheet, style, rows=c(1:(1+nrow(sheetData))), cols=right_align.idx, gridExpand=TRUE, stack=TRUE)
					}
				} else if(length(grep("sample information", sheet)) > 0) {
					center_align.idx = which(colnames(sheetData) %in% c('Group','Symp Group','Sample ID','Gender'))
					if(length(center_align.idx) > 0) {
						style <- createStyle(halign="center", valign="center")
						addStyle(wb, sheet, style, rows=c(1:(1+nrow(sheetData))), cols=center_align.idx, gridExpand=TRUE, stack=TRUE)
					}
					
					right_align.idx = which(colnames(sheetData) %in% c('Birth Year','Onset Age','Age','Onset Period'))
					if(length(right_align.idx) > 0) {
						style <- createStyle(halign="right", valign="center")
						addStyle(wb, sheet, style, rows=c(1:(1+nrow(sheetData))), cols=right_align.idx, gridExpand=TRUE, stack=TRUE)
					}
				}
			}
		}
	}
	
	if(length(wb$sheet_names) > 0) {
		openxlsx::saveWorkbook(wb, file=paste('./Data/02_significant/piRNA/',Group[1],' vs ',Group[2],'/','piRNA_datatable ','(',Group[1],' vs ',Group[2],')','.xlsx',sep=''), overwrite=TRUE)
	}
	
	# --- write xlsx file by openxlsx package ---------------------------------------------------------
	if(m == 1) {
		wb2 <- openxlsx::createWorkbook("./Data/02_significant/piRNA/piRNA_datatable (pairwise comparison).xlsx")
		modifyBaseFont(wb2, fontSize=12, fontColour="black", fontName="Times New Roman")
		start.idx = 0
	} else {
		start.idx = 1
	}
	
	for(s in start.idx:2)
	{
		if(s == 0) {
			sheet = "sample information"
			sheetData = diseaseStage.ALL
		} else {
			if(s == 1) {
				sheet = paste(Group[1],' > ',Group[2],sep='')
				sheetData = piRNA.normalizedData.Comp.UPinGroup1
			} else if(s == 2) {
				sheet = paste(Group[1],' < ',Group[2],sep='')
				sheetData = piRNA.normalizedData.Comp.UPinGroup2
			}
		}
		
		if(!is.null(sheetData))
		{
			if(nrow(sheetData) > 0)
			{
				#na.idx = which(is.na(sheetData), arr.ind=TRUE)
				#if(nrow(na.idx) > 0)
				#{
				#	for(ri in 1:nrow(na.idx))
				#	{
				#		sheetData[na.idx[ri,'row'],na.idx[ri,'col']] = ''
				#	}
				#}
				
				addWorksheet(wb2, sheet)
				writeData(wb2, sheet, sheetData, colNames=TRUE, rowNames=FALSE, keepNA=FALSE)
				
				setRowHeights(wb2, sheet, rows=c(1:(1+nrow(sheetData))), heights=16.5)
				if(length(grep("sample information", sheet)) == 0) {
					if(length(grep("common", sheet, ignore.case=TRUE)) > 0 | length(grep("pvalue", sheet, ignore.case=TRUE)) > 0 | length(grep(">", sheet, ignore.case=TRUE)) > 0 | length(grep("<", sheet, ignore.case=TRUE)) > 0) {
						setColWidths(wb2, sheet, cols=c(1:ncol(sheetData)), widths=c(10,30,17.5,rep(12.5,ncol(sheetData)-3)))
						freezePane(wb2, sheet, firstActiveRow=2, firstActiveCol=4+4)
					} else {
						setColWidths(wb2, sheet, cols=c(1:ncol(sheetData)), widths=c(10,30,17.5,rep(12.5,ncol(sheetData)-3)))
						freezePane(wb2, sheet, firstActiveRow=2, firstActiveCol=4)
					}
				} else if(length(grep("sample information", sheet)) > 0) {
					setColWidths(wb2, sheet, cols=c(1:ncol(sheetData)), widths=c(7.5,rep(12.5,ncol(sheetData)-4),45,15,12.5))
					freezePane(wb2, sheet, firstActiveRow=2)
				}
				
				style <- createStyle(halign="left", valign="center")
				addStyle(wb2, sheet, style, rows=c(1:(1+nrow(sheetData))), cols=c(1:ncol(sheetData)), gridExpand=TRUE)
				style <- createStyle(textDecoration="bold")
				addStyle(wb2, sheet, style, rows=1, cols=c(1:ncol(sheetData)), stack=TRUE)
				
				if(length(grep("sample information", sheet)) == 0) {
					right_align.idx = which(! colnames(sheetData) %in% c('type','Name','Accession_ID','Drives_from','locus_GRCh38','locus_hg19'))
					if(length(right_align.idx) > 0) {
						style <- createStyle(halign="right", valign="center")
						addStyle(wb2, sheet, style, rows=c(1:(1+nrow(sheetData))), cols=right_align.idx, gridExpand=TRUE, stack=TRUE)
					}
				} else if(length(grep("sample information", sheet)) > 0) {
					center_align.idx = which(colnames(sheetData) %in% c('Group','Symp Group','Sample ID','Gender'))
					if(length(center_align.idx) > 0) {
						style <- createStyle(halign="center", valign="center")
						addStyle(wb2, sheet, style, rows=c(1:(1+nrow(sheetData))), cols=center_align.idx, gridExpand=TRUE, stack=TRUE)
					}
					
					right_align.idx = which(colnames(sheetData) %in% c('Birth Year','Onset Age','Age','Onset Period'))
					if(length(right_align.idx) > 0) {
						style <- createStyle(halign="right", valign="center")
						addStyle(wb2, sheet, style, rows=c(1:(1+nrow(sheetData))), cols=right_align.idx, gridExpand=TRUE, stack=TRUE)
					}
				}
			}
		}
	}
	
	if(m == ncol(combn.idx)) {
		if(length(wb2$sheet_names) > 0) {
			openxlsx::saveWorkbook(wb2, file="./Data/02_significant/piRNA/piRNA_datatable (pairwise comparison).xlsx", overwrite=TRUE)
		}
	}
}

