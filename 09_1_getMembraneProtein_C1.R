rm(list=ls())	# remove all variables in workspace

setwd("D:/Exo")

Sample <- "C1"

mainDir <- paste('./Data/',Sample,sep="")
subDir <- "09_membraneProtein"
if (file.exists(file.path(mainDir, subDir))){
} else {
	dir.create(file.path(mainDir, subDir))
}
mainDir <- paste('./Data/',Sample,'/09_membraneProtein',sep="")
subDir <- "subset"
if (file.exists(file.path(mainDir, subDir))){
} else {
	dir.create(file.path(mainDir, subDir))
}


#install.packages("xlsx")				# xlsx: Read, write, format Excel 2007 and Excel 97/2000/XP/2003 files
#Sys.setenv(JAVA_HOME='C:/Program Files (x86)/Java/jdk1.8.0_171/jre')	# for 32-bit (change directory when the version is changed)
#Sys.setenv(JAVA_HOME='C:/Program Files/Java/jdk1.8.0_171/jre')			# for 64-bit (change directory when the version is changed)
#options(java.parameters = "-Xmx2000m")	# https://stackoverflow.com/questions/21937640/handling-java-lang-outofmemoryerror-when-writing-to-excel-from-r
#options(java.parameters = "-Xmx2g")	# https://stackoverflow.com/questions/7963393/out-of-memory-error-java-when-using-r-and-xlconnect-package
library(xlsx)

# Importing a big xlsx file into R? https://stackoverflow.com/questions/19147884/importing-a-big-xlsx-file-into-r/43118530#43118530
# Building R for Windows: https://cran.r-project.org/bin/windows/Rtools/
#install.packages("openxlsx")			# openxlsx: Read, Write and Edit XLSX Files
library(openxlsx)


if(Sample == "C1"){
	batch <- c('20150828_C1_1','20150828_C1_2','20150828_C1_3')
}else if(Sample == "C1_8R"){
	batch <- c('20141121_sev_C2','20150311_Sev_exosome','20150324_exosome1','20150324_exosome2','20150428_MinShun_100X','20150428_MinShun_100X_2','C1_150','C1_5th')
}else if(Sample == "C1_ALL"){
	batch <- c('20141121_sev_C2','20150311_Sev_exosome','20150324_exosome1','20150324_exosome2','20150428_MinShun_100X','20150428_MinShun_100X_2','C1_150','C1_5th','20150828_C1_1','20150828_C1_2','20150828_C1_3')
}else if(Sample == "E8"){
	batch <- c('20150828_E8_1','20150828_E8_2','20150828_E8_3')
}else{
	print(paste('Please provide filename of Sample ',Sample,sep=""))
}


# --- keyword of membrane-related proteins ---
membrane <- c('membrane protein','membrane','transmembrane protein','transmembrane')
surface <- c('surface')
tetraspanin <- c('tetraspanin','tspan')
CD <- c('CD','CD antigen')
transport <- c('transporter','transport')
channel <- c('channel')
integrin <- c('integrin')
lipidraft <- c('lipid raft','flotillin','caveolin')
membrane.category <- c(membrane, surface, tetraspanin, CD, transport, channel, integrin, lipidraft)

membrane.FP <- c('Basement membrane','Basement-membrane')	# false-positive keyword
CD.FP <- c('Cdc','cds','CDV','CD2-associated')				# false-positive keyword
integrin.FP <- c('integrin interactor')						# false-positive keyword


# === subset ====================================================================================================================
folder <- c('high','low')

for(r in 1:length(batch))
{
	if (file.exists(paste('./Data/',Sample,'/04_drug/subset/n_',r,sep=""))) {
		
		mainDir <- paste('./Data/',Sample,'/09_membraneProtein/subset',sep="")
		subDir <- paste('n_',r,sep="")
		if (file.exists(file.path(mainDir, subDir))){
		} else {
			dir.create(file.path(mainDir, subDir))
		}
		
		for(f in 1:length(folder))
		{
			if (file.exists(paste('./Data/',Sample,'/04_drug/subset/n_',r,'/',folder[f],sep=""))) {
				
				if (file.exists(paste('./Data/',Sample,'/04_drug/subset/n_',r,'/',folder[f],'/',Sample,'_exo_UniProt_drug.xlsx',sep=""))) {
					
					# === input =======================================================================================
					#input <- read.csv(paste('./Data/',Sample,'/04_drug/subset/n_',r,'/',folder[f],'/',Sample,'_exo_UniProt_drug.csv',sep=""), header=TRUE, stringsAsFactors=FALSE)
					input <- read.xlsx(paste('./Data/',Sample,'/04_drug/subset/n_',r,'/',folder[f],'/',Sample,'_exo_UniProt_drug.xlsx',sep=""), sheet='UniProt_drug', colNames=TRUE, rowNames=FALSE, check.names=FALSE)
					mouse.Protein.names <- as.character(input$mouse.Protein.names)
					
					
					# === search for membrane-related proteins by name ================================================
					category.matched.idx <- c()
					category.matched.value <- c()
					#category.matched.value.levels <- c()
					matched.idx <- c()
					
					# --- write xlsx file by openxlsx package ---------------------------------------------------------
					# Importing a big xlsx file into R? https://stackoverflow.com/questions/19147884/importing-a-big-xlsx-file-into-r/43118530#43118530
					# Building R for Windows: https://cran.r-project.org/bin/windows/Rtools/
					#Sys.setenv("R_ZIPCMD" = "C:/Rtools/bin/zip.exe")	# path to zip.exe
					
					wb <- openxlsx::createWorkbook(paste('./Data/',Sample,'/09_membraneProtein/subset/n_',r,'/',folder[f],'/',Sample,'_exo_UniProt_membraneRelated','.xlsx',sep=""))
					modifyBaseFont(wb, fontSize=12, fontColour="black", fontName="Times New Roman")
					
					wb2 <- openxlsx::createWorkbook(paste('./Data/',Sample,'/09_membraneProtein/',Sample,'_exo_UniProt_membraneRelated','.xlsx',sep=""))
					modifyBaseFont(wb2, fontSize=12, fontColour="black", fontName="Times New Roman")
					
					for(j in 1:length(membrane.category))
					{
						category.matched.idx <- grep(membrane.category[j], mouse.Protein.names, ignore.case=TRUE)
						category.matched.value <- grep(membrane.category[j], mouse.Protein.names, ignore.case=TRUE, value=TRUE)
						#category.matched.value.levels <- levels(as.factor(category.matched.value))
						#category.matched.value.levels = sort(category.matched.value.levels)
						#print(category.matched.idx)
						#print(length(category.matched.idx))
						#print(category.matched.value)
						#print(length(category.matched.value))
						#print(category.matched.value.levels)
						#print(length(category.matched.value.levels))
						
						if(length(category.matched.idx) > 0)
						{
							# /// remove false-positive matching (start) ///
							FP <- c()
							if (membrane.category[j] %in% membrane) {
								FP <- membrane.FP
							} else if (membrane.category[j] %in% CD) {
								FP <- CD.FP
							} else if (membrane.category[j] %in% integrin) {
								FP <- integrin.FP
							}
							
							if(length(FP) > 0)
							{
								FP.comb.matched.idx <- c()
								FP.comb.matched.value <- c()
								for(k in 1:length(FP))
								{
									FP.matched.idx <- grep(FP[k], mouse.Protein.names[category.matched.idx], ignore.case=TRUE)
									FP.matched.value <- grep(FP[k], mouse.Protein.names[category.matched.idx], ignore.case=TRUE, value=TRUE)
									
									if(length(FP.matched.idx) > 0)
									{
										#print(FP.matched.idx)
										#print(FP.matched.value)
										
										FP.comb.matched.idx <- c(FP.comb.matched.idx, FP.matched.idx)
										FP.comb.matched.value <- c(FP.comb.matched.value, FP.matched.value)
									}
								}
								
								if(length(FP.comb.matched.idx) > 0)
								{
									#print(FP.comb.matched.idx)
									#print(FP.comb.matched.value)
									
									category.matched.idx <- category.matched.idx[-FP.comb.matched.idx]
									category.matched.value <- category.matched.value[-FP.comb.matched.idx]
									#print(category.matched.idx)
									#print(category.matched.value)
								}
							}
							# /// remove false-positive matching (end) ///
							
							if(length(category.matched.idx) > 0)
							{
								mainDir <- paste('./Data/',Sample,'/09_membraneProtein/subset/n_',r,sep="")
								subDir <- folder[f]
								if (file.exists(file.path(mainDir, subDir))){
								} else {
									dir.create(file.path(mainDir, subDir))
								}
								
								category.matched.info <- input[category.matched.idx,,drop=FALSE]
								#write.csv(category.matched.info, paste('./Data/',Sample,'/09_membraneProtein/subset/n_',r,'/',folder[f],'/',Sample,'_exo_UniProt_membraneRelated','_',gsub(' ','',membrane.category[j]),'.csv',sep=""), quote=TRUE, row.names=FALSE, na="")
								#if(r == 1 & folder[f] == 'high')	write.csv(category.matched.info, paste('./Data/',Sample,'/09_membraneProtein/',Sample,'_exo_UniProt_membraneRelated','_',gsub(' ','',membrane.category[j]),'.csv',sep=""), quote=TRUE, row.names=FALSE, na="")
								
								sheet = membrane.category[j]
								sheetData = category.matched.info
								
								if(!is.null(sheetData))
								{
									if(nrow(sheetData) > 0)
									{
										addWorksheet(wb, sheet)
										writeData(wb, sheet, sheetData, colNames=TRUE, rowNames=FALSE, keepNA=FALSE)
										setColWidths(wb, sheet, cols=c(1:ncol(sheetData)), widths=12.5)
										setRowHeights(wb, sheet, rows=c(1:(1+nrow(sheetData))), heights=16.5)
										freezePane(wb, sheet, firstActiveRow=2)
										style <- createStyle(halign="left", valign="center")
										addStyle(wb, sheet, style, rows=c(1:(1+nrow(sheetData))), cols=c(1:ncol(sheetData)), gridExpand=TRUE)
										style <- createStyle(textDecoration="bold")
										addStyle(wb, sheet, style, rows=1, cols=c(1:ncol(sheetData)), stack=TRUE)
										
										right_align.idx = which(colnames(sheetData) %in% c('mouse.Frequency'))
										if(length(right_align.idx) > 0) {
											style <- createStyle(halign="right", valign="center")
											addStyle(wb, sheet, style, rows=c(2:(1+nrow(sheetData))), cols=right_align.idx, gridExpand=TRUE, stack=TRUE)
										}
										
										if(r == 1 & folder[f] == 'high')
										{
											addWorksheet(wb2, sheet)
											writeData(wb2, sheet, sheetData, colNames=TRUE, rowNames=FALSE, keepNA=FALSE)
											setColWidths(wb2, sheet, cols=c(1:ncol(sheetData)), widths=12.5)
											setRowHeights(wb2, sheet, rows=c(1:(1+nrow(sheetData))), heights=16.5)
											freezePane(wb2, sheet, firstActiveRow=2)
											style <- createStyle(halign="left", valign="center")
											addStyle(wb2, sheet, style, rows=c(1:(1+nrow(sheetData))), cols=c(1:ncol(sheetData)), gridExpand=TRUE)
											style <- createStyle(textDecoration="bold")
											addStyle(wb2, sheet, style, rows=1, cols=c(1:ncol(sheetData)), stack=TRUE)
											
											right_align.idx = which(colnames(sheetData) %in% c('mouse.Frequency'))
											if(length(right_align.idx) > 0) {
												style <- createStyle(halign="right", valign="center")
												addStyle(wb2, sheet, style, rows=c(2:(1+nrow(sheetData))), cols=right_align.idx, gridExpand=TRUE, stack=TRUE)
											}
										}
									}
								}
								
								matched.idx <- c(matched.idx, category.matched.idx)
							}
						}
					}
					
					if(length(matched.idx) > 0)
					{
						unique.matched.idx <- sort(unique(matched.idx))
						if(length(unique.matched.idx) > 0)
						{
							unique.matched.info <- input[unique.matched.idx,,drop=FALSE]
							#write.csv(unique.matched.info, paste('./Data/',Sample,'/09_membraneProtein/subset/n_',r,'/',folder[f],'/',Sample,'_exo_UniProt_membraneRelated','.csv',sep=""), quote=TRUE, row.names=FALSE, na="")
							#if(r == 1 & folder[f] == 'high')	write.csv(unique.matched.info, paste('./Data/',Sample,'/09_membraneProtein/',Sample,'_exo_UniProt_membraneRelated','.csv',sep=""), quote=TRUE, row.names=FALSE, na="")
							
							unique.matched.info.name = unique.matched.info[,c('mouse.UniProt','mouse.Gene.ID','mouse.MGI.ID','mouse.MGI.Symbol','mouse.Gene.names','mouse.Protein.names')]
							#write.csv(unique.matched.info.name, paste('./Data/',Sample,'/09_membraneProtein/subset/n_',r,'/',folder[f],'/',Sample,'_exo_UniProt_membraneRelated','.name','.csv',sep=""), quote=TRUE, row.names=FALSE, na="")
							#if(r == 1 & folder[f] == 'high')	write.csv(unique.matched.info.name, paste('./Data/',Sample,'/09_membraneProtein/',Sample,'_exo_UniProt_membraneRelated','.name','.csv',sep=""), quote=TRUE, row.names=FALSE, na="")
							
							for(s in 1:2)
							{
								if(s == 1) {
									sheet = paste(Sample,'_exo_membraneRelated',sep="")
									sheetData = unique.matched.info
								} else if(s == 2) {
									sheet = paste(Sample,'_exo_membraneRelated','.','name',sep="")
									sheetData = unique.matched.info.name
								}
								
								if(!is.null(sheetData))
								{
									if(nrow(sheetData) > 0)
									{
										addWorksheet(wb, sheet)
										writeData(wb, sheet, sheetData, colNames=TRUE, rowNames=FALSE, keepNA=FALSE)
										setColWidths(wb, sheet, cols=c(1:ncol(sheetData)), widths=12.5)
										setRowHeights(wb, sheet, rows=c(1:(1+nrow(sheetData))), heights=16.5)
										freezePane(wb, sheet, firstActiveRow=2)
										style <- createStyle(halign="left", valign="center")
										addStyle(wb, sheet, style, rows=c(1:(1+nrow(sheetData))), cols=c(1:ncol(sheetData)), gridExpand=TRUE)
										style <- createStyle(textDecoration="bold")
										addStyle(wb, sheet, style, rows=1, cols=c(1:ncol(sheetData)), stack=TRUE)
										
										right_align.idx = which(colnames(sheetData) %in% c('mouse.Frequency'))
										if(length(right_align.idx) > 0) {
											style <- createStyle(halign="right", valign="center")
											addStyle(wb, sheet, style, rows=c(2:(1+nrow(sheetData))), cols=right_align.idx, gridExpand=TRUE, stack=TRUE)
										}
										
										if(r == 1 & folder[f] == 'high')
										{
											addWorksheet(wb2, sheet)
											writeData(wb2, sheet, sheetData, colNames=TRUE, rowNames=FALSE, keepNA=FALSE)
											setColWidths(wb2, sheet, cols=c(1:ncol(sheetData)), widths=12.5)
											setRowHeights(wb2, sheet, rows=c(1:(1+nrow(sheetData))), heights=16.5)
											freezePane(wb2, sheet, firstActiveRow=2)
											style <- createStyle(halign="left", valign="center")
											addStyle(wb2, sheet, style, rows=c(1:(1+nrow(sheetData))), cols=c(1:ncol(sheetData)), gridExpand=TRUE)
											style <- createStyle(textDecoration="bold")
											addStyle(wb2, sheet, style, rows=1, cols=c(1:ncol(sheetData)), stack=TRUE)
											
											right_align.idx = which(colnames(sheetData) %in% c('mouse.Frequency'))
											if(length(right_align.idx) > 0) {
												style <- createStyle(halign="right", valign="center")
												addStyle(wb2, sheet, style, rows=c(2:(1+nrow(sheetData))), cols=right_align.idx, gridExpand=TRUE, stack=TRUE)
											}
										}
									}
								}
							}
						}
					}
					
					sheetName = c(paste(Sample,'_exo_membraneRelated',sep=""), paste(Sample,'_exo_membraneRelated','.','name',sep=""), membrane.category)
					worksheetOrder(wb) = unlist(sapply(sheetName, FUN=function(N) which(names(wb) %in% unlist(N))))
					if(length(wb$sheet_names) > 0) {
						openxlsx::saveWorkbook(wb, file=paste('./Data/',Sample,'/09_membraneProtein/subset/n_',r,'/',folder[f],'/',Sample,'_exo_UniProt_membraneRelated','.xlsx',sep=""), overwrite=TRUE)
					}
					
					if(r == 1 & folder[f] == 'high')
					{
						sheetName = c(paste(Sample,'_exo_membraneRelated',sep=""), paste(Sample,'_exo_membraneRelated','.','name',sep=""), membrane.category)
						worksheetOrder(wb2) = unlist(sapply(sheetName, FUN=function(N) which(names(wb2) %in% unlist(N))))
						if(length(wb2$sheet_names) > 0) {
							openxlsx::saveWorkbook(wb2, file=paste('./Data/',Sample,'/09_membraneProtein/',Sample,'_exo_UniProt_membraneRelated','.xlsx',sep=""), overwrite=TRUE)
						}
					}
					
					# === read and write to reformat sheetOrder =======================================================
					for(s in 1:2)
					{
						if(s == 1) {
							pwd = paste('./Data/',Sample,'/09_membraneProtein/subset/n_',r,'/',folder[f],'/',Sample,'_exo_UniProt_membraneRelated','.xlsx',sep="")
						} else if(s == 2) {
							pwd = paste('./Data/',Sample,'/09_membraneProtein/',Sample,'_exo_UniProt_membraneRelated','.xlsx',sep="")
						}
						
						if(file.exists(pwd))
						{
							# --- write xlsx file by openxlsx package -----------------------------------------------------------------------------
							# Importing a big xlsx file into R? https://stackoverflow.com/questions/19147884/importing-a-big-xlsx-file-into-r/43118530#43118530
							# Building R for Windows: https://cran.r-project.org/bin/windows/Rtools/
							#Sys.setenv("R_ZIPCMD" = "C:/Rtools/bin/zip.exe")	# path to zip.exe
							
							wb <- openxlsx::createWorkbook(pwd)
							modifyBaseFont(wb, fontSize=12, fontColour="black", fontName="Times New Roman")
							
							for(sh in 1:length(names(getSheets(xlsx::loadWorkbook(pwd)))))
							{
								sheet = names(getSheets(xlsx::loadWorkbook(pwd)))[sh]
								
								sheetData <- xlsx::read.xlsx(pwd, sheetName=sheet, startRow=1, header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
								
								if(!is.null(sheetData))
								{
									if(nrow(sheetData) > 0)
									{
										addWorksheet(wb, sheet)
										writeData(wb, sheet, sheetData, colNames=TRUE, rowNames=FALSE, keepNA=FALSE)
										
										setColWidths(wb, sheet, cols=c(1:ncol(sheetData)), widths=12.5)
										setRowHeights(wb, sheet, rows=c(1:(1+nrow(sheetData))), heights=16.5)
										freezePane(wb, sheet, firstActiveRow=2)
										
										style <- createStyle(halign="left", valign="center")
										addStyle(wb, sheet, style, rows=c(1:(1+nrow(sheetData))), cols=c(1:ncol(sheetData)), gridExpand=TRUE)
										style <- createStyle(textDecoration="bold")
										addStyle(wb, sheet, style, rows=1, cols=c(1:ncol(sheetData)), stack=TRUE)
										
										right_align.idx = which(colnames(sheetData) %in% c('mouse.Frequency'))
										if(length(right_align.idx) > 0) {
											style <- createStyle(halign="right", valign="center")
											addStyle(wb, sheet, style, rows=c(2:(1+nrow(sheetData))), cols=right_align.idx, gridExpand=TRUE, stack=TRUE)
										}
									}
								}
							}
							
							if(length(wb$sheet_names) > 0) {
								openxlsx::saveWorkbook(wb, file=pwd, overwrite=TRUE)
							}
						}
					}
				}
			}
		}
	}
}

