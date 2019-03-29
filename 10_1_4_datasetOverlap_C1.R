rm(list=ls())	# remove all variables in workspace

setwd("D:/Exo")

Sample <- "C1"

mainDir <- paste('./Data/',Sample,sep="")
subDir <- "10_datasetOverlap"
if (file.exists(file.path(mainDir, subDir))){
} else {
	dir.create(file.path(mainDir, subDir))
}
mainDir <- paste('./Data/',Sample,'/10_datasetOverlap',sep="")
subDir <- "Data"
if (file.exists(file.path(mainDir, subDir))){
} else {
	dir.create(file.path(mainDir, subDir))
}
mainDir <- paste('./Data/',Sample,'/10_datasetOverlap/Data',sep="")
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
#library(xlsx)

# Importing a big xlsx file into R? https://stackoverflow.com/questions/19147884/importing-a-big-xlsx-file-into-r/43118530#43118530
# Building R for Windows: https://cran.r-project.org/bin/windows/Rtools/
#install.packages("openxlsx")			# openxlsx: Read, Write and Edit XLSX Files
library(openxlsx)

#install.packages("VennDiagram")		# VennDiagram: Generate High-Resolution Venn and Euler Plots
library(VennDiagram)


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

C1.batch.num = length(c('20150828_C1_1','20150828_C1_2','20150828_C1_3'))
C1_8R.batch.num = length(c('20141121_sev_C2','20150311_Sev_exosome','20150324_exosome1','20150324_exosome2','20150428_MinShun_100X','20150428_MinShun_100X_2','C1_150','C1_5th'))
C1_ALL.batch.num = length(c('20141121_sev_C2','20150311_Sev_exosome','20150324_exosome1','20150324_exosome2','20150428_MinShun_100X','20150428_MinShun_100X_2','C1_150','C1_5th','20150828_C1_1','20150828_C1_2','20150828_C1_3'))
E8.batch.num = length(c('20150828_E8_1','20150828_E8_2','20150828_E8_3'))


# === input Top100 recorded in Microvesicle database ============================================================================
# --- ExoCarta --------------------------------------------------------------------------------------------------------
ExoCarta.ProteinNumber_Type <- c('Top100','All')

# /// mouse format ///
#ExoCarta_Top100Protein.GeneSymbol2MGI_table.IdentificationNumber <- read.csv(paste('./Data/',Sample,'/10_datasetOverlap/Database/mmu/ExoCarta_Top100.Protein_GeneSymbol2MGI.table_IdentificationNumber.csv',sep=""), header=TRUE, stringsAsFactors=FALSE)
ExoCarta_Top100Protein.GeneSymbol2MGI_table.IdentificationNumber <- read.xlsx(paste('./Data/',Sample,'/10_datasetOverlap/Database/mmu/ExoCarta_Top100.Protein_GeneSymbol2MGI.table.xlsx',sep=""), sheet='ExoCarta_Top100 (with number)', colNames=TRUE, rowNames=FALSE, check.names=FALSE)
ExoCarta_Top100.X <- as.character(ExoCarta_Top100Protein.GeneSymbol2MGI_table.IdentificationNumber[,1])
ExoCarta_Top100.mgi_symbol <- as.character(ExoCarta_Top100Protein.GeneSymbol2MGI_table.IdentificationNumber[,'mgi_symbol'])
ExoCarta_Top100.gene_symbol = ExoCarta_Top100.mgi_symbol
ExoCarta_Top100.Number_of_times_identified <- as.character(ExoCarta_Top100Protein.GeneSymbol2MGI_table.IdentificationNumber[,'Number_of_times_identified'])
ExoCarta_Top100.mgi_id <- as.character(ExoCarta_Top100Protein.GeneSymbol2MGI_table.IdentificationNumber[,'mgi_id'])
ExoCarta_Top100.entrezgene <- as.character(ExoCarta_Top100Protein.GeneSymbol2MGI_table.IdentificationNumber[,'entrezgene'])
ExoCarta_Top100.uniprotswissprot <- as.character(ExoCarta_Top100Protein.GeneSymbol2MGI_table.IdentificationNumber[,'uniprotswissprot'])
ExoCarta_Top100.description <- as.character(ExoCarta_Top100Protein.GeneSymbol2MGI_table.IdentificationNumber[,'description'])
ExoCarta_Top100.Protein_names <- as.character(ExoCarta_Top100Protein.GeneSymbol2MGI_table.IdentificationNumber[,'Protein_names'])
#length(ExoCarta_Top100.entrezgene)

#ExoCarta_AllProtein.Gene2MGI_table <- read.csv(paste('./Data/',Sample,'/10_datasetOverlap/Database/mmu/ExoCarta_All.Protein_Gene2MGI.table.csv',sep=""), header=TRUE, stringsAsFactors=FALSE)
ExoCarta_AllProtein.Gene2MGI_table <- read.xlsx(paste('./Data/',Sample,'/10_datasetOverlap/Database/mmu/ExoCarta_All.Protein_Gene2MGI.table.xlsx',sep=""), sheet='ExoCarta_All', colNames=TRUE, rowNames=FALSE, check.names=FALSE)
ExoCarta_All.X <- as.character(ExoCarta_AllProtein.Gene2MGI_table[,1])
ExoCarta_All.mgi_symbol <- as.character(ExoCarta_AllProtein.Gene2MGI_table[,'mgi_symbol'])
ExoCarta_All.gene_symbol = ExoCarta_All.mgi_symbol
ExoCarta_All.mgi_id <- as.character(ExoCarta_AllProtein.Gene2MGI_table[,'mgi_id'])
ExoCarta_All.entrezgene <- as.character(ExoCarta_AllProtein.Gene2MGI_table[,'entrezgene'])
ExoCarta_All.uniprotswissprot <- as.character(ExoCarta_AllProtein.Gene2MGI_table[,'uniprotswissprot'])
ExoCarta_All.description <- as.character(ExoCarta_AllProtein.Gene2MGI_table[,'description'])
ExoCarta_All.Protein_names <- as.character(ExoCarta_AllProtein.Gene2MGI_table[,'Protein_names'])
#length(ExoCarta_All.entrezgene)


# === subset ====================================================================================================================
folder <- c('high','low')

for(r in 1:length(batch))
{
	if (file.exists(paste('./Data/',Sample,'/03_homologene/subset/n_',r,sep=""))) {
		
		mainDir <- paste('./Data/',Sample,'/10_datasetOverlap/Data/subset',sep="")
		subDir <- paste('n_',r,sep="")
		if (file.exists(file.path(mainDir, subDir))){
		} else {
			dir.create(file.path(mainDir, subDir))
		}
		
		for(f in 1:length(folder))
		{
			dataset.Sample1 <- c('C1','C1_8R','C1_ALL')
			dataset.Sample2 <- c('E8')
			dataset <- c(dataset.Sample1, dataset.Sample2)
			
			for(i in 1:length(dataset))
			{
				homologene <- c()
				homologene.mouse.UniProt <- c()
				homologene.mouse.GI <- c()
				homologene.mouse.Gene.ID <- c()
				homologene.mouse.MGI.ID <- c()
				homologene.mouse.MGI.Symbol <- c()
				homologene.mouse.Gene.names <- c()
				homologene.mouse.Protein.names <- c()
				homologene.mouse.Frequency <- c()
				homologene.HID <- c()
				homologene.human.UniProt <- c()
				homologene.human.GI <- c()
				homologene.human.Gene.ID <- c()
				homologene.human.HGNC.ID <- c()
				homologene.human.HGNC.Symbol <- c()
				homologene.human.Gene.names <- c()
				homologene.human.Protein.names <- c()
				
				assign(paste(dataset[i], '.', 'homologene', sep=""), homologene)
				
				assign(paste(dataset[i], '.', 'homologene.mouse.UniProt', sep=""), homologene.mouse.UniProt)
				assign(paste(dataset[i], '.', 'homologene.mouse.GI', sep=""), homologene.mouse.GI)
				assign(paste(dataset[i], '.', 'homologene.mouse.Gene.ID', sep=""), homologene.mouse.Gene.ID)
				assign(paste(dataset[i], '.', 'homologene.mouse.MGI.ID', sep=""), homologene.mouse.MGI.ID)
				assign(paste(dataset[i], '.', 'homologene.mouse.MGI.Symbol', sep=""), homologene.mouse.MGI.Symbol)
				assign(paste(dataset[i], '.', 'homologene.mouse.Gene.names', sep=""), homologene.mouse.Gene.names)
				assign(paste(dataset[i], '.', 'homologene.mouse.Protein.names', sep=""), homologene.mouse.Protein.names)
				assign(paste(dataset[i], '.', 'homologene.mouse.Frequency', sep=""), homologene.mouse.Frequency)
				
				assign(paste(dataset[i], '.', 'homologene.HID', sep=""), homologene.HID)
				
				assign(paste(dataset[i], '.', 'homologene.human.UniProt', sep=""), homologene.human.UniProt)
				assign(paste(dataset[i], '.', 'homologene.human.GI', sep=""), homologene.human.GI)
				assign(paste(dataset[i], '.', 'homologene.human.Gene.ID', sep=""), homologene.human.Gene.ID)
				assign(paste(dataset[i], '.', 'homologene.human.HGNC.ID', sep=""), homologene.human.HGNC.ID)
				assign(paste(dataset[i], '.', 'homologene.human.HGNC.Symbol', sep=""), homologene.human.HGNC.Symbol)
				assign(paste(dataset[i], '.', 'homologene.human.Gene.names', sep=""), homologene.human.Gene.names)
				assign(paste(dataset[i], '.', 'homologene.human.Protein.names', sep=""), homologene.human.Protein.names)
				
				
				dataset.batch.num = eval(parse(text=paste(dataset[i], '.', 'batch.num', sep="")))
				
				if(r > dataset.batch.num & folder[f] == 'low') {
					new_r = 1
					new_folder = 'high'
				} else {
					new_r = r
					new_folder = folder[f]
				}
				
				if (file.exists(paste('./Data/',dataset[i],'/03_homologene/subset/n_',new_r,'/',new_folder,sep=""))) {
					
					if (file.exists(paste('./Data/',dataset[i],'/03_homologene/subset/n_',new_r,'/',new_folder,'/',dataset[i],'_exo_UniProt_homologene.xlsx',sep=""))) {
						
						# === input GeneID ==============================================================================================================
						#homologene <- read.csv(paste('./Data/',dataset[i],'/03_homologene/subset/n_',new_r,'/',new_folder,'/',dataset[i],'_exo_UniProt_homologene.csv',sep=""), header=TRUE, stringsAsFactors=FALSE)
						homologene <- read.xlsx(paste('./Data/',dataset[i],'/03_homologene/subset/n_',new_r,'/',new_folder,'/',dataset[i],'_exo_UniProt_homologene.xlsx',sep=""), sheet='homologene (mouse2human)', colNames=TRUE, rowNames=FALSE, check.names=FALSE)
						
						homologene.mouse.UniProt <- as.character(homologene$mouse.UniProt)
						homologene.mouse.GI <- as.character(homologene$mouse.GI)
						homologene.mouse.Gene.ID <- as.character(homologene$mouse.Gene.ID)
						homologene.mouse.MGI.ID <- as.character(homologene$mouse.MGI.ID)
						homologene.mouse.MGI.Symbol <- as.character(homologene$mouse.MGI.Symbol)
						homologene.mouse.Gene.names <- as.character(homologene$mouse.Gene.names)
						homologene.mouse.Protein.names <- as.character(homologene$mouse.Protein.names)
						homologene.mouse.Frequency <- as.character(homologene$mouse.Frequency)
						homologene.mouse.Gene.ID[which(is.na(homologene.mouse.Gene.ID))] = ''
						
						homologene.HID <- as.character(homologene$HID)
						
						homologene.human.UniProt <- as.character(homologene$human.UniProt)
						homologene.human.GI <- as.character(homologene$human.GI)
						homologene.human.Gene.ID <- as.character(homologene$human.Gene.ID)
						homologene.human.HGNC.ID <- as.character(homologene$human.HGNC.ID)
						homologene.human.HGNC.Symbol <- as.character(homologene$human.HGNC.Symbol)
						homologene.human.Gene.names <- as.character(homologene$human.Gene.names)
						homologene.human.Protein.names <- as.character(homologene$human.Protein.names)
						homologene.human.Gene.ID[which(is.na(homologene.human.Gene.ID))] = ''
						
						#length(homologene.mouse.Gene.ID)
						
						
						assign(paste(dataset[i], '.', 'homologene', sep=""), homologene)
						
						assign(paste(dataset[i], '.', 'homologene.mouse.UniProt', sep=""), homologene.mouse.UniProt)
						assign(paste(dataset[i], '.', 'homologene.mouse.GI', sep=""), homologene.mouse.GI)
						assign(paste(dataset[i], '.', 'homologene.mouse.Gene.ID', sep=""), homologene.mouse.Gene.ID)
						assign(paste(dataset[i], '.', 'homologene.mouse.MGI.ID', sep=""), homologene.mouse.MGI.ID)
						assign(paste(dataset[i], '.', 'homologene.mouse.MGI.Symbol', sep=""), homologene.mouse.MGI.Symbol)
						assign(paste(dataset[i], '.', 'homologene.mouse.Gene.names', sep=""), homologene.mouse.Gene.names)
						assign(paste(dataset[i], '.', 'homologene.mouse.Protein.names', sep=""), homologene.mouse.Protein.names)
						assign(paste(dataset[i], '.', 'homologene.mouse.Frequency', sep=""), homologene.mouse.Frequency)
						
						assign(paste(dataset[i], '.', 'homologene.HID', sep=""), homologene.HID)
						
						assign(paste(dataset[i], '.', 'homologene.human.UniProt', sep=""), homologene.human.UniProt)
						assign(paste(dataset[i], '.', 'homologene.human.GI', sep=""), homologene.human.GI)
						assign(paste(dataset[i], '.', 'homologene.human.Gene.ID', sep=""), homologene.human.Gene.ID)
						assign(paste(dataset[i], '.', 'homologene.human.HGNC.ID', sep=""), homologene.human.HGNC.ID)
						assign(paste(dataset[i], '.', 'homologene.human.HGNC.Symbol', sep=""), homologene.human.HGNC.Symbol)
						assign(paste(dataset[i], '.', 'homologene.human.Gene.names', sep=""), homologene.human.Gene.names)
						assign(paste(dataset[i], '.', 'homologene.human.Protein.names', sep=""), homologene.human.Protein.names)
					}
				}
			}
			
			
			if (file.exists(paste('./Data/',Sample,'/03_homologene/subset/n_',r,'/',folder[f],sep=""))) {
				
				if (file.exists(paste('./Data/',Sample,'/03_homologene/subset/n_',r,'/',folder[f],'/',Sample,'_exo_UniProt_homologene.xlsx',sep=""))) {
					
					mainDir <- paste('./Data/',Sample,'/10_datasetOverlap/Data/subset/n_',r,sep="")
					subDir <- folder[f]
					if (file.exists(file.path(mainDir, subDir))){
					} else {
						dir.create(file.path(mainDir, subDir))
					}
					
					
					# === overlap ===================================================================================================================
					if(Sample == 'E8') {
						Sample.CompareMaxNum = 2
					} else {
						Sample.CompareMaxNum = 3
					}
					
					for(n in 1:Sample.CompareMaxNum)		# dataset number
					{
						mainDir <- paste('./Data/',Sample,'/10_datasetOverlap/Data/subset/n_',r,'/',folder[f],sep="")
						subDir <- paste(n,'_dataset',sep="")
						if (file.exists(file.path(mainDir, subDir))){
						} else {
							dir.create(file.path(mainDir, subDir))
						}
						
						if(n == 1)
						{
							mainDir <- paste('./Data/',Sample,'/10_datasetOverlap/Data/subset/n_',r,'/',folder[f],'/',n,'_dataset',sep="")
							subDir <- "Microvesicle"
							if (file.exists(file.path(mainDir, subDir))){
							} else {
								dir.create(file.path(mainDir, subDir))
							}
							mainDir <- paste('./Data/',Sample,'/10_datasetOverlap/Data/subset/n_',r,'/',folder[f],'/',n,'_dataset/Microvesicle',sep="")
							subDir <- "ExoCarta"
							if (file.exists(file.path(mainDir, subDir))){
							} else {
								dir.create(file.path(mainDir, subDir))
							}
						}
						
						
						combn.idx <- combn(length(dataset), n)
						#print(combn.idx)
						#print(dim(combn.idx))
						
						for(i in 1:dim(combn.idx)[2])		# dataset combination
						{
							combnName <- paste(dataset[combn.idx[,i]], collapse='.')
							#print(combnName)
							
							if(n == 1) {		# n = 1
								
								Name1 <- dataset[combn.idx[1,i]]
								Name1.homologene.mouse.Gene.ID = eval(parse(text=paste(Name1, '.', 'homologene.mouse.Gene.ID', sep="")))
								Name1.homologene.mouse.UniProt = eval(parse(text=paste(Name1, '.', 'homologene.mouse.UniProt', sep="")))
								Name1.homologene.mouse.MGI.ID = eval(parse(text=paste(Name1, '.', 'homologene.mouse.MGI.ID', sep="")))
								Name1.homologene.mouse.MGI.Symbol = eval(parse(text=paste(Name1, '.', 'homologene.mouse.MGI.Symbol', sep="")))
								Name1.homologene.mouse.Gene.names = eval(parse(text=paste(Name1, '.', 'homologene.mouse.Gene.names', sep="")))
								Name1.homologene.mouse.Protein.names = eval(parse(text=paste(Name1, '.', 'homologene.mouse.Protein.names', sep="")))
								input1 = Name1.homologene.mouse.Gene.ID[which(Name1.homologene.mouse.Gene.ID != '')]
								input1.idx = which(Name1.homologene.mouse.Gene.ID != '')
								
								if(length(input1) == 0) {
									# Do nothing
								} else {
									if(Name1 == Sample)
									{
										# === compare to database (ExoCarta) ==================================================================
										for(p in 1:length(ExoCarta.ProteinNumber_Type))
										{
											ExoCarta_Protein.entrezgene = eval(parse(text=paste('ExoCarta_', ExoCarta.ProteinNumber_Type[p], '.entrezgene', sep="")))
											ExoCarta_Protein.uniprotswissprot = eval(parse(text=paste('ExoCarta_', ExoCarta.ProteinNumber_Type[p], '.uniprotswissprot', sep="")))
											ExoCarta_Protein.mgi_id = eval(parse(text=paste('ExoCarta_', ExoCarta.ProteinNumber_Type[p], '.mgi_id', sep="")))
											ExoCarta_Protein.gene_symbol = eval(parse(text=paste('ExoCarta_', ExoCarta.ProteinNumber_Type[p], '.gene_symbol', sep="")))
											ExoCarta_Protein.description = eval(parse(text=paste('ExoCarta_', ExoCarta.ProteinNumber_Type[p], '.description', sep="")))
											ExoCarta_Protein.Protein_names = eval(parse(text=paste('ExoCarta_', ExoCarta.ProteinNumber_Type[p], '.Protein_names', sep="")))
											
											input0 = ExoCarta_Protein.entrezgene
											
											input.union <- union(input0,input1)
											input.intersect <- intersect(input0,input1)
											input.setdiff.0 <- setdiff(input0,input1)
											input.setdiff.1 <- setdiff(input1,input0)
											#print(length(input.union))
											#print(length(input.intersect))
											#print(length(input.setdiff.0))
											#print(length(input.setdiff.1))
											
											
											prename = paste('./Data/',Sample,'/10_datasetOverlap/Data/subset/n_',r,'/',folder[f],'/',n,'_dataset/Microvesicle/ExoCarta/','ExoCarta_',ExoCarta.ProteinNumber_Type[p],'.Protein','_vs_',combnName,'.exosome', sep="")
											
											wb <- openxlsx::createWorkbook(paste(prename,'.xlsx',sep=""))
											modifyBaseFont(wb, fontSize=12, fontColour="black", fontName="Times New Roman")
											
											overlapType <- c('union','intersect','setdiff.0','setdiff.1')
											for(o in 1:length(overlapType))
											{
												input.overlapType = eval(parse(text=paste('input', '.', overlapType[o], sep="")))
												
												if(length(input.overlapType) > 0)
												{
													# --- txt -----------------------------------------------------------------------------------------
													write.table(input.overlapType, paste(prename,'_',overlapType[o],'.txt',sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE)
													
													# --- csv -----------------------------------------------------------------------------------------
													input.overlapType.Sample1.idx = sapply(input.overlapType, FUN=function(X) which(Name1.homologene.mouse.Gene.ID %in% X))
													input.overlapType.ExoCarta.idx = sapply(input.overlapType, FUN=function(X) which(ExoCarta_Protein.entrezgene %in% X))
													uniprot_id.all <- c()
													mgi_id.all <- c()
													gene_symbol.all <- c()
													gene_name.all <- c()
													protein_name.all <- c()
													for(j in 1:length(input.overlapType))
													{
														Sample1.HitNumber = length(unlist(input.overlapType.Sample1.idx[j]))
														ExoCarta.HitNumber = length(unlist(input.overlapType.ExoCarta.idx[j]))
														
														if(Sample1.HitNumber == 1) {
															uniprot_id = Name1.homologene.mouse.UniProt[unlist(input.overlapType.Sample1.idx[j])]
															mgi_id = Name1.homologene.mouse.MGI.ID[unlist(input.overlapType.Sample1.idx[j])]
															gene_symbol = Name1.homologene.mouse.MGI.Symbol[unlist(input.overlapType.Sample1.idx[j])]
															gene_name = Name1.homologene.mouse.Gene.names[unlist(input.overlapType.Sample1.idx[j])]
															protein_name = Name1.homologene.mouse.Protein.names[unlist(input.overlapType.Sample1.idx[j])]
														} else if(Sample1.HitNumber == 0) {
															if(ExoCarta.HitNumber == 1) {
																uniprot_id = ExoCarta_Protein.uniprotswissprot[unlist(input.overlapType.ExoCarta.idx[j])]
																mgi_id = ExoCarta_Protein.mgi_id[unlist(input.overlapType.ExoCarta.idx[j])]
																gene_symbol = ExoCarta_Protein.gene_symbol[unlist(input.overlapType.ExoCarta.idx[j])]
																gene_name = ExoCarta_Protein.description[unlist(input.overlapType.ExoCarta.idx[j])]
																protein_name = ExoCarta_Protein.Protein_names[unlist(input.overlapType.ExoCarta.idx[j])]
															} else if(ExoCarta.HitNumber == 0) {
																uniprot_id = ''
																mgi_id = ''
																gene_symbol = ''
																gene_name = ''
																protein_name = ''
																cat('[Warning] mouse.Gene.ID:',input.overlapType[j],' does not have any id mapping to \"Name1.homologene.mouse.Gene.ID (',dataset[combn.idx[1,i]],'.homologene.mouse.Gene.ID)\" and ExoCarta_Protein.entrezgene.\n', sep="")
															} else if(ExoCarta.HitNumber > 1) {
																uniprot_id = paste(unique(ExoCarta_Protein.uniprotswissprot[unlist(input.overlapType.ExoCarta.idx[j])]), collapse="|")
																mgi_id = paste(unique(ExoCarta_Protein.mgi_id[unlist(input.overlapType.ExoCarta.idx[j])]), collapse="|")
																gene_symbol = paste(unique(ExoCarta_Protein.gene_symbol[unlist(input.overlapType.ExoCarta.idx[j])]), collapse="|")
																gene_name = paste(unique(ExoCarta_Protein.description[unlist(input.overlapType.ExoCarta.idx[j])]), collapse="|")
																protein_name = paste(unique(ExoCarta_Protein.Protein_names[unlist(input.overlapType.ExoCarta.idx[j])]), collapse="|")
																
																if(length(unique(ExoCarta_Protein.uniprotswissprot[unlist(input.overlapType.ExoCarta.idx[j])])) > 1)
																cat('[Note] mouse.Gene.ID:',input.overlapType[j],' has more than two id mapping to \"ExoCarta_Protein.entrezgene\".\n', sep="")
															}
														} else if(Sample1.HitNumber > 1) {
															uniprot_id = paste(unique(Name1.homologene.mouse.UniProt[unlist(input.overlapType.Sample1.idx[j])]), collapse="|")
															mgi_id = paste(unique(Name1.homologene.mouse.MGI.ID[unlist(input.overlapType.Sample1.idx[j])]), collapse="|")
															gene_symbol = paste(unique(Name1.homologene.mouse.MGI.Symbol[unlist(input.overlapType.Sample1.idx[j])]), collapse="|")
															gene_name = paste(unique(Name1.homologene.mouse.Gene.names[unlist(input.overlapType.Sample1.idx[j])]), collapse="|")
															protein_name = paste(unique(Name1.homologene.mouse.Protein.names[unlist(input.overlapType.Sample1.idx[j])]), collapse="|")
															
															if(length(unique(Name1.homologene.mouse.UniProt[unlist(input.overlapType.Sample1.idx[j])])) > 1)
															cat('[Note] mouse.Gene.ID:',input.overlapType[j],' has more than two id mapping to \"Name1.homologene.mouse.Gene.ID (',dataset[combn.idx[1,i]],'.homologene.mouse.Gene.ID)\".\n', sep="")
														}
														
														uniprot_id.all = c(uniprot_id.all, uniprot_id)
														mgi_id.all = c(mgi_id.all, mgi_id)
														gene_symbol.all = c(gene_symbol.all, gene_symbol)
														gene_name.all = c(gene_name.all, gene_name)
														protein_name.all = c(protein_name.all, protein_name)
													}
													input.overlapType.info = cbind(input.overlapType, uniprot_id.all, mgi_id.all, gene_symbol.all, gene_name.all, protein_name.all)
													colnames(input.overlapType.info) = c('mouse.Gene.ID','mouse.UniProt','mouse.MGI.ID','mouse.MGI.Symbol','mouse.Gene.names','mouse.Protein.names')
													#write.table(input.overlapType.info, paste(prename,'_',overlapType[o],'.csv',sep=""), quote=TRUE, na="", sep=",", row.names=FALSE)
													
													# --- xlsx ----------------------------------------------------------------------------------------
													sheet = overlapType[o]
													sheetData = input.overlapType.info
													
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
														}
													}
												}
											}
											
											if(length(wb$sheet_names) > 0) {
												openxlsx::saveWorkbook(wb, file=paste(prename,'.xlsx',sep=""), overwrite=TRUE)
											}
											
											# --- png -----------------------------------------------------------------------------------------
											### Plot venn diagram of 2 datasets ###
											venn.plot <- draw.pairwise.venn(
											area1 = length(input0),								# The size of the first set
											area2 = length(input1),								# The size of the second set
											
											cross.area = length(input.intersect),				# The size of the intersection between the sets
											
											category = c(paste('ExoCarta_',ExoCarta.ProteinNumber_Type[p],' ',paste('(',length(input0),')',sep=""),sep=""), paste(dataset[combn.idx[1,i]],'_exosome',' ',paste('(',length(input1),')',sep=""),sep="")),	# category name (number) of the sets
											cat.pos = c(-15,15),								# the positions (in degrees) of the category names along the circles
											cat.dist = c(0.02,0.02),							# the distances (in npc units) of the category names from the edges of the circles
											cat.cex = 3,										# the size of the category names
											cat.col = c("black","black"),						# the colours of the category names
											#fill = c("gray15","gray30"),						# the colours of the circles' areas
											fill = c("white","white"),							# the colours of the circles' areas
											#alpha = rep(1,2),									# the alpha transparency of the circles' areas
											#alpha = c(0.4,0.6),									# the alpha transparency of the circles' areas
											#lty = "blank",										# the line dash pattern of the circles' circumferences
											cex = 3												# the size of the areas' labels
											)
											
											png(paste(prename,'_vennDiagram.png',sep=""), width=14, height=14, units="in", res=600)
											grid.draw(venn.plot)
											dev.off()
											
											grid.draw(venn.plot);
											grid.newpage();							# clear plotting area for new plots
										}
									}
								}
							} else if(n == 2) {		# n = 2
								
								Name1 <- dataset[combn.idx[1,i]]
								Name1.homologene.mouse.Gene.ID = eval(parse(text=paste(Name1, '.', 'homologene.mouse.Gene.ID', sep="")))
								Name1.homologene.mouse.UniProt = eval(parse(text=paste(Name1, '.', 'homologene.mouse.UniProt', sep="")))
								Name1.homologene.mouse.MGI.ID = eval(parse(text=paste(Name1, '.', 'homologene.mouse.MGI.ID', sep="")))
								Name1.homologene.mouse.MGI.Symbol = eval(parse(text=paste(Name1, '.', 'homologene.mouse.MGI.Symbol', sep="")))
								Name1.homologene.mouse.Gene.names = eval(parse(text=paste(Name1, '.', 'homologene.mouse.Gene.names', sep="")))
								Name1.homologene.mouse.Protein.names = eval(parse(text=paste(Name1, '.', 'homologene.mouse.Protein.names', sep="")))
								input1 = Name1.homologene.mouse.Gene.ID[which(Name1.homologene.mouse.Gene.ID != '')]
								input1.idx = which(Name1.homologene.mouse.Gene.ID != '')
								
								Name2 <- dataset[combn.idx[2,i]]
								Name2.homologene.mouse.Gene.ID = eval(parse(text=paste(Name2, '.', 'homologene.mouse.Gene.ID', sep="")))
								Name2.homologene.mouse.UniProt = eval(parse(text=paste(Name2, '.', 'homologene.mouse.UniProt', sep="")))
								Name2.homologene.mouse.MGI.ID = eval(parse(text=paste(Name2, '.', 'homologene.mouse.MGI.ID', sep="")))
								Name2.homologene.mouse.MGI.Symbol = eval(parse(text=paste(Name2, '.', 'homologene.mouse.MGI.Symbol', sep="")))
								Name2.homologene.mouse.Gene.names = eval(parse(text=paste(Name2, '.', 'homologene.mouse.Gene.names', sep="")))
								Name2.homologene.mouse.Protein.names = eval(parse(text=paste(Name2, '.', 'homologene.mouse.Protein.names', sep="")))
								input2 = Name2.homologene.mouse.Gene.ID[which(Name2.homologene.mouse.Gene.ID != '')]
								input2.idx = which(Name2.homologene.mouse.Gene.ID != '')
								
								if(length(input1) == 0 & length(input2) == 0) {
									# Do nothing
								} else {
									RUN = 0
									if(Sample == 'E8') {
										if((Name1=='E8') | (Name2=='E8')) {
											RUN = 1
										} else {
											RUN = 0
										}
									} else {
										RUN = 1
									}
									
									if(RUN)
									{
										mainDir <- paste('./Data/',Sample,'/10_datasetOverlap/Data/subset/n_',r,'/',folder[f],'/',n,'_dataset',sep="")
										subDir <- combnName
										if (file.exists(file.path(mainDir, subDir))){
										} else {
											dir.create(file.path(mainDir, subDir))
										}
										
										
										# === pairwise comparison =================================================================================
										input.union <- union(input1,input2)
										input.intersect <- intersect(input1,input2)
										input.setdiff.1 <- setdiff(input1,input2)
										input.setdiff.2 <- setdiff(input2,input1)
										#print(length(input.union))
										#print(length(input.intersect))
										#print(length(input.setdiff.1))
										#print(length(input.setdiff.2))
										
										
										prename = paste('./Data/',Sample,'/10_datasetOverlap/Data/subset/n_',r,'/',folder[f],'/',n,'_dataset/',combnName,'/', combnName, '.exosome', sep="")
										
										wb <- openxlsx::createWorkbook(paste(prename,'.xlsx',sep=""))
										modifyBaseFont(wb, fontSize=12, fontColour="black", fontName="Times New Roman")
										
										overlapType <- c('union','intersect','setdiff.1','setdiff.2')
										for(o in 1:length(overlapType))
										{
											input.overlapType = eval(parse(text=paste('input', '.', overlapType[o], sep="")))
											
											if(length(input.overlapType) > 0)
											{
												# --- txt -----------------------------------------------------------------------------------------
												write.table(input.overlapType, paste(prename,'_',overlapType[o],'.txt',sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE)
												
												# --- csv -----------------------------------------------------------------------------------------
												input.overlapType.Sample1.idx = sapply(input.overlapType, FUN=function(X) which(Name1.homologene.mouse.Gene.ID %in% X))
												input.overlapType.Sample2.idx = sapply(input.overlapType, FUN=function(X) which(Name2.homologene.mouse.Gene.ID %in% X))
												uniprot_id.all <- c()
												mgi_id.all <- c()
												gene_symbol.all <- c()
												gene_name.all <- c()
												protein_name.all <- c()
												for(j in 1:length(input.overlapType))
												{
													Sample1.HitNumber = length(unlist(input.overlapType.Sample1.idx[j]))
													Sample2.HitNumber = length(unlist(input.overlapType.Sample2.idx[j]))
													
													if(Sample1.HitNumber == 1) {
														uniprot_id = Name1.homologene.mouse.UniProt[unlist(input.overlapType.Sample1.idx[j])]
														mgi_id = Name1.homologene.mouse.MGI.ID[unlist(input.overlapType.Sample1.idx[j])]
														gene_symbol = Name1.homologene.mouse.MGI.Symbol[unlist(input.overlapType.Sample1.idx[j])]
														gene_name = Name1.homologene.mouse.Gene.names[unlist(input.overlapType.Sample1.idx[j])]
														protein_name = Name1.homologene.mouse.Protein.names[unlist(input.overlapType.Sample1.idx[j])]
													} else if(Sample1.HitNumber == 0) {
														if(Sample2.HitNumber == 1) {
															uniprot_id = Name2.homologene.mouse.UniProt[unlist(input.overlapType.Sample2.idx[j])]
															mgi_id = Name2.homologene.mouse.MGI.ID[unlist(input.overlapType.Sample2.idx[j])]
															gene_symbol = Name2.homologene.mouse.MGI.Symbol[unlist(input.overlapType.Sample2.idx[j])]
															gene_name = Name2.homologene.mouse.Gene.names[unlist(input.overlapType.Sample2.idx[j])]
															protein_name = Name2.homologene.mouse.Protein.names[unlist(input.overlapType.Sample2.idx[j])]
														} else if(Sample2.HitNumber == 0) {
															uniprot_id = ''
															mgi_id = ''
															gene_symbol = ''
															gene_name = ''
															protein_name = ''
															cat('[Warning] mouse.Gene.ID:',input.overlapType[j],' does not have any id mapping to all homologene.mouse.Gene.ID.\n', sep="")
														} else if(Sample2.HitNumber > 1) {
															uniprot_id = paste(unique(Name2.homologene.mouse.UniProt[unlist(input.overlapType.Sample2.idx[j])]), collapse="|")
															mgi_id = paste(unique(Name2.homologene.mouse.MGI.ID[unlist(input.overlapType.Sample2.idx[j])]), collapse="|")
															gene_symbol = paste(unique(Name2.homologene.mouse.MGI.Symbol[unlist(input.overlapType.Sample2.idx[j])]), collapse="|")
															gene_name = paste(unique(Name2.homologene.mouse.Gene.names[unlist(input.overlapType.Sample2.idx[j])]), collapse="|")
															protein_name = paste(unique(Name2.homologene.mouse.Protein.names[unlist(input.overlapType.Sample2.idx[j])]), collapse="|")
															
															if(length(unique(Name2.homologene.mouse.UniProt[unlist(input.overlapType.Sample2.idx[j])])) > 1)
															cat('[Note] mouse.Gene.ID:',input.overlapType[j],' has more than two id mapping to \"Name2.homologene.mouse.Gene.ID (',dataset[combn.idx[2,i]],'.homologene.mouse.Gene.ID)\".\n', sep="")
														}
													} else if(Sample1.HitNumber > 1) {
														uniprot_id = paste(unique(Name1.homologene.mouse.UniProt[unlist(input.overlapType.Sample1.idx[j])]), collapse="|")
														mgi_id = paste(unique(Name1.homologene.mouse.MGI.ID[unlist(input.overlapType.Sample1.idx[j])]), collapse="|")
														gene_symbol = paste(unique(Name1.homologene.mouse.MGI.Symbol[unlist(input.overlapType.Sample1.idx[j])]), collapse="|")
														gene_name = paste(unique(Name1.homologene.mouse.Gene.names[unlist(input.overlapType.Sample1.idx[j])]), collapse="|")
														protein_name = paste(unique(Name1.homologene.mouse.Protein.names[unlist(input.overlapType.Sample1.idx[j])]), collapse="|")
														
														if(length(unique(Name1.homologene.mouse.UniProt[unlist(input.overlapType.Sample1.idx[j])])) > 1)
														cat('[Note] mouse.Gene.ID:',input.overlapType[j],' has more than two id mapping to \"Name1.homologene.mouse.Gene.ID (',dataset[combn.idx[1,i]],'.homologene.mouse.Gene.ID)\".\n', sep="")
													}
													
													uniprot_id.all = c(uniprot_id.all, uniprot_id)
													mgi_id.all = c(mgi_id.all, mgi_id)
													gene_symbol.all = c(gene_symbol.all, gene_symbol)
													gene_name.all = c(gene_name.all, gene_name)
													protein_name.all = c(protein_name.all, protein_name)
												}
												input.overlapType.info = cbind(input.overlapType, uniprot_id.all, mgi_id.all, gene_symbol.all, gene_name.all, protein_name.all)
												colnames(input.overlapType.info) = c('mouse.Gene.ID','mouse.UniProt','mouse.MGI.ID','mouse.MGI.Symbol','mouse.Gene.names','mouse.Protein.names')
												#write.table(input.overlapType.info, paste(prename,'_',overlapType[o],'.csv',sep=""), quote=TRUE, na="", sep=",", row.names=FALSE)
												
												# --- xlsx ----------------------------------------------------------------------------------------
												sheet = overlapType[o]
												sheetData = input.overlapType.info
												
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
													}
												}
											}
										}
										
										if(length(wb$sheet_names) > 0) {
											openxlsx::saveWorkbook(wb, file=paste(prename,'.xlsx',sep=""), overwrite=TRUE)
										}
										
										# --- png -----------------------------------------------------------------------------------------
										### Plot venn diagram of 2 datasets ###
										venn.plot <- draw.pairwise.venn(
										area1 = length(input1),								# The size of the first set
										area2 = length(input2),								# The size of the second set
										
										cross.area = length(input.intersect),				# The size of the intersection between the sets
										
										category = c(paste(dataset[combn.idx[1,i]],'_exosome',' ',paste('(',length(input1),')',sep=""),sep=""), paste(dataset[combn.idx[2,i]],'_exosome',' ',paste('(',length(input2),')',sep=""),sep="")),	# category name (number) of the sets
										cat.pos = c(-15,15),								# the positions (in degrees) of the category names along the circles
										cat.dist = c(0.02,0.02),							# the distances (in npc units) of the category names from the edges of the circles
										cat.cex = 3,										# the size of the category names
										cat.col = c("black","black"),						# the colours of the category names
										#fill = c("gray30","gray45"),						# the colours of the circles' areas
										fill = c("white","white"),							# the colours of the circles' areas
										#alpha = rep(1,2),									# the alpha transparency of the circles' areas
										#alpha = c(0.6,0.8),									# the alpha transparency of the circles' areas
										#lty = "blank",										# the line dash pattern of the circles' circumferences
										cex = 3												# the size of the areas' labels
										)
										
										png(paste(prename,'_vennDiagram.png',sep=""), width=14, height=14, units="in", res=600)
										grid.draw(venn.plot)
										dev.off()
										
										grid.draw(venn.plot);
										grid.newpage();							# clear plotting area for new plots
										
										
										# === compare to database (ExoCarta) ==================================================================
										for(p in 1:length(ExoCarta.ProteinNumber_Type))
										{
											ExoCarta_Protein.entrezgene = eval(parse(text=paste('ExoCarta_', ExoCarta.ProteinNumber_Type[p], '.entrezgene', sep="")))
											ExoCarta_Protein.uniprotswissprot = eval(parse(text=paste('ExoCarta_', ExoCarta.ProteinNumber_Type[p], '.uniprotswissprot', sep="")))
											ExoCarta_Protein.mgi_id = eval(parse(text=paste('ExoCarta_', ExoCarta.ProteinNumber_Type[p], '.mgi_id', sep="")))
											ExoCarta_Protein.gene_symbol = eval(parse(text=paste('ExoCarta_', ExoCarta.ProteinNumber_Type[p], '.gene_symbol', sep="")))
											ExoCarta_Protein.description = eval(parse(text=paste('ExoCarta_', ExoCarta.ProteinNumber_Type[p], '.description', sep="")))
											ExoCarta_Protein.Protein_names = eval(parse(text=paste('ExoCarta_', ExoCarta.ProteinNumber_Type[p], '.Protein_names', sep="")))
											
											input0 = ExoCarta_Protein.entrezgene
											
											input.union <- Reduce(union, list(input0,input1,input2))
											input.intersect <- Reduce(intersect, list(input0,input1,input2))		#intersect(intersect(input0,input1),input2)
											#print(length(input.union))
											#print(length(input.intersect))
											
											input.intersect.01 <- intersect(input0,input1)
											input.intersect.02 <- intersect(input0,input2)
											input.intersect.12 <- intersect(input1,input2)
											#print(length(input.intersect.01))
											#print(length(input.intersect.02))
											#print(length(input.intersect.12))
											
											input.intersect.01_2 <- setdiff(input.intersect.01,input.intersect)
											input.intersect.02_1 <- setdiff(input.intersect.02,input.intersect)
											input.intersect.12_0 <- setdiff(input.intersect.12,input.intersect)
											#print(length(input.intersect.01_2))
											#print(length(input.intersect.02_1))
											#print(length(input.intersect.12_0))
											
											input.setdiff.0 <- Reduce(setdiff, list(input0,input1,input2))
											input.setdiff.1 <- Reduce(setdiff, list(input1,input0,input2))
											input.setdiff.2 <- Reduce(setdiff, list(input2,input0,input1))
											#print(length(input.setdiff.0))
											#print(length(input.setdiff.1))
											#print(length(input.setdiff.2))
											
											
											prename = paste('./Data/',Sample,'/10_datasetOverlap/Data/subset/n_',r,'/',folder[f],'/',n,'_dataset/',combnName,'/','ExoCarta_',ExoCarta.ProteinNumber_Type[p],'.Protein','_vs_',combnName,'.exosome', sep="")
											
											wb <- openxlsx::createWorkbook(paste(prename,'.xlsx',sep=""))
											modifyBaseFont(wb, fontSize=12, fontColour="black", fontName="Times New Roman")
											
											overlapType <- c('union','intersect','intersect.01','intersect.02','intersect.12','intersect.01_2','intersect.02_1','intersect.12_0','setdiff.0','setdiff.1','setdiff.2')
											for(o in 1:length(overlapType))
											{
												input.overlapType = eval(parse(text=paste('input', '.', overlapType[o], sep="")))
												
												if(length(input.overlapType) > 0)
												{
													# --- txt -----------------------------------------------------------------------------------------
													write.table(input.overlapType, paste(prename,'_',overlapType[o],'.txt',sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE)
													
													# --- csv -----------------------------------------------------------------------------------------
													input.overlapType.Sample1.idx = sapply(input.overlapType, FUN=function(X) which(Name1.homologene.mouse.Gene.ID %in% X))
													input.overlapType.Sample2.idx = sapply(input.overlapType, FUN=function(X) which(Name2.homologene.mouse.Gene.ID %in% X))
													input.overlapType.ExoCarta.idx = sapply(input.overlapType, FUN=function(X) which(ExoCarta_Protein.entrezgene %in% X))
													uniprot_id.all <- c()
													mgi_id.all <- c()
													gene_symbol.all <- c()
													gene_name.all <- c()
													protein_name.all <- c()
													for(j in 1:length(input.overlapType))
													{
														Sample1.HitNumber = length(unlist(input.overlapType.Sample1.idx[j]))
														Sample2.HitNumber = length(unlist(input.overlapType.Sample2.idx[j]))
														ExoCarta.HitNumber = length(unlist(input.overlapType.ExoCarta.idx[j]))
														
														if(Sample1.HitNumber == 1) {
															uniprot_id = Name1.homologene.mouse.UniProt[unlist(input.overlapType.Sample1.idx[j])]
															mgi_id = Name1.homologene.mouse.MGI.ID[unlist(input.overlapType.Sample1.idx[j])]
															gene_symbol = Name1.homologene.mouse.MGI.Symbol[unlist(input.overlapType.Sample1.idx[j])]
															gene_name = Name1.homologene.mouse.Gene.names[unlist(input.overlapType.Sample1.idx[j])]
															protein_name = Name1.homologene.mouse.Protein.names[unlist(input.overlapType.Sample1.idx[j])]
														} else if(Sample1.HitNumber == 0) {
															if(Sample2.HitNumber == 1) {
																uniprot_id = Name2.homologene.mouse.UniProt[unlist(input.overlapType.Sample2.idx[j])]
																mgi_id = Name2.homologene.mouse.MGI.ID[unlist(input.overlapType.Sample2.idx[j])]
																gene_symbol = Name2.homologene.mouse.MGI.Symbol[unlist(input.overlapType.Sample2.idx[j])]
																gene_name = Name2.homologene.mouse.Gene.names[unlist(input.overlapType.Sample2.idx[j])]
																protein_name = Name2.homologene.mouse.Protein.names[unlist(input.overlapType.Sample2.idx[j])]
															} else if(Sample2.HitNumber == 0) {
																if(ExoCarta.HitNumber == 1) {
																	uniprot_id = ExoCarta_Protein.uniprotswissprot[unlist(input.overlapType.ExoCarta.idx[j])]
																	mgi_id = ExoCarta_Protein.mgi_id[unlist(input.overlapType.ExoCarta.idx[j])]
																	gene_symbol = ExoCarta_Protein.gene_symbol[unlist(input.overlapType.ExoCarta.idx[j])]
																	gene_name = ExoCarta_Protein.description[unlist(input.overlapType.ExoCarta.idx[j])]
																	protein_name = ExoCarta_Protein.Protein_names[unlist(input.overlapType.ExoCarta.idx[j])]
																} else if(ExoCarta.HitNumber == 0) {
																	uniprot_id = ''
																	mgi_id = ''
																	gene_symbol = ''
																	gene_name = ''
																	protein_name = ''
																	cat('[Warning] mouse.Gene.ID:',input.overlapType[j],' does not have any id mapping to all homologene.mouse.Gene.ID and ExoCarta_Protein.entrezgene.\n', sep="")
																} else if(ExoCarta.HitNumber > 1) {
																	uniprot_id = paste(unique(ExoCarta_Protein.uniprotswissprot[unlist(input.overlapType.ExoCarta.idx[j])]), collapse="|")
																	mgi_id = paste(unique(ExoCarta_Protein.mgi_id[unlist(input.overlapType.ExoCarta.idx[j])]), collapse="|")
																	gene_symbol = paste(unique(ExoCarta_Protein.gene_symbol[unlist(input.overlapType.ExoCarta.idx[j])]), collapse="|")
																	gene_name = paste(unique(ExoCarta_Protein.description[unlist(input.overlapType.ExoCarta.idx[j])]), collapse="|")
																	protein_name = paste(unique(ExoCarta_Protein.Protein_names[unlist(input.overlapType.ExoCarta.idx[j])]), collapse="|")
																	
																	if(length(unique(ExoCarta_Protein.uniprotswissprot[unlist(input.overlapType.ExoCarta.idx[j])])) > 1)
																	cat('[Note] mouse.Gene.ID:',input.overlapType[j],' has more than two id mapping to \"ExoCarta_Protein.entrezgene\".\n', sep="")
																}
															} else if(Sample2.HitNumber > 1) {
																uniprot_id = paste(unique(Name2.homologene.mouse.UniProt[unlist(input.overlapType.Sample2.idx[j])]), collapse="|")
																mgi_id = paste(unique(Name2.homologene.mouse.MGI.ID[unlist(input.overlapType.Sample2.idx[j])]), collapse="|")
																gene_symbol = paste(unique(Name2.homologene.mouse.MGI.Symbol[unlist(input.overlapType.Sample2.idx[j])]), collapse="|")
																gene_name = paste(unique(Name2.homologene.mouse.Gene.names[unlist(input.overlapType.Sample2.idx[j])]), collapse="|")
																protein_name = paste(unique(Name2.homologene.mouse.Protein.names[unlist(input.overlapType.Sample2.idx[j])]), collapse="|")
																
																if(length(unique(Name2.homologene.mouse.UniProt[unlist(input.overlapType.Sample2.idx[j])])) > 1)
																cat('[Note] mouse.Gene.ID:',input.overlapType[j],' has more than two id mapping to \"Name2.homologene.mouse.Gene.ID (',dataset[combn.idx[2,i]],'.homologene.mouse.Gene.ID)\".\n', sep="")
															}
														} else if(Sample1.HitNumber > 1) {
															uniprot_id = paste(unique(Name1.homologene.mouse.UniProt[unlist(input.overlapType.Sample1.idx[j])]), collapse="|")
															mgi_id = paste(unique(Name1.homologene.mouse.MGI.ID[unlist(input.overlapType.Sample1.idx[j])]), collapse="|")
															gene_symbol = paste(unique(Name1.homologene.mouse.MGI.Symbol[unlist(input.overlapType.Sample1.idx[j])]), collapse="|")
															gene_name = paste(unique(Name1.homologene.mouse.Gene.names[unlist(input.overlapType.Sample1.idx[j])]), collapse="|")
															protein_name = paste(unique(Name1.homologene.mouse.Protein.names[unlist(input.overlapType.Sample1.idx[j])]), collapse="|")
															
															if(length(unique(Name1.homologene.mouse.UniProt[unlist(input.overlapType.Sample1.idx[j])])) > 1)
															cat('[Note] mouse.Gene.ID:',input.overlapType[j],' has more than two id mapping to \"Name1.homologene.mouse.Gene.ID (',dataset[combn.idx[1,i]],'.homologene.mouse.Gene.ID)\".\n', sep="")
														}
														
														uniprot_id.all = c(uniprot_id.all, uniprot_id)
														mgi_id.all = c(mgi_id.all, mgi_id)
														gene_symbol.all = c(gene_symbol.all, gene_symbol)
														gene_name.all = c(gene_name.all, gene_name)
														protein_name.all = c(protein_name.all, protein_name)
													}
													input.overlapType.info = cbind(input.overlapType, uniprot_id.all, mgi_id.all, gene_symbol.all, gene_name.all, protein_name.all)
													colnames(input.overlapType.info) = c('mouse.Gene.ID','mouse.UniProt','mouse.MGI.ID','mouse.MGI.Symbol','mouse.Gene.names','mouse.Protein.names')
													#write.table(input.overlapType.info, paste(prename,'_',overlapType[o],'.csv',sep=""), quote=TRUE, na="", sep=",", row.names=FALSE)
													
													# --- xlsx ----------------------------------------------------------------------------------------
													sheet = overlapType[o]
													sheetData = input.overlapType.info
													
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
														}
													}
												}
											}
											
											if(length(wb$sheet_names) > 0) {
												openxlsx::saveWorkbook(wb, file=paste(prename,'.xlsx',sep=""), overwrite=TRUE)
											}
											
											# --- png -----------------------------------------------------------------------------------------
											### Plot venn diagram of 3 datasets ###
											venn.plot <- draw.triple.venn(
											area1 = length(input0),			# The size of the first set
											area2 = length(input1),			# The size of the second set
											area3 = length(input2),			# The size of the third set
											
											n12 = length(input.intersect.01),		# The size of the intersection between the first and the second set
											n23 = length(input.intersect.12),		# The size of the intersection between the second and the third set
											n13 = length(input.intersect.02),		# The size of the intersection between the first and the third set
											n123 = length(input.intersect),			# The size of the intersection between all three sets
											
											category = c(paste('ExoCarta_',ExoCarta.ProteinNumber_Type[p],' ',paste('(',length(input0),')',sep=""),sep=""), paste(dataset[combn.idx[1,i]],'_exosome',' ',paste('(',length(input1),')',sep=""),sep=""), paste(dataset[combn.idx[2,i]],'_exosome',' ',paste('(',length(input2),')',sep=""),sep="")),	# category name (number) of the sets
											cat.pos = c(-15,15,180),				# the positions (in degrees) of the category names along the circles
											cat.dist = c(0.03,0.03,0.02),			# the distances (in npc units) of the category names from the edges of the circles
											cat.cex = 3,							# the size of the category names
											cat.col = c("black","black","black"),	# the colours of the category names
											#fill = c("gray15","gray30","gray45"),	# the colours of the circles' areas
											fill = c("white","white","white"),		# the colours of the circles' areas
											#alpha = rep(1,3),						# the alpha transparency of the circles' areas
											#alpha = c(0.4,0.6,0.8),					# the alpha transparency of the circles' areas
											#lty = "blank",							# the line dash pattern of the circles' circumferences
											cex = 3									# the size of the areas' labels
											)
											
											png(paste(prename,'_vennDiagram.png',sep=""), width=14, height=14, units="in", res=600)
											grid.draw(venn.plot)
											dev.off()
											
											grid.draw(venn.plot);
											grid.newpage();							# clear plotting area for new plots
										}
									}
								}
							} else if(n == 3) {		# n = 3
								
								Name1 <- dataset[combn.idx[1,i]]
								Name1.homologene.mouse.Gene.ID = eval(parse(text=paste(Name1, '.', 'homologene.mouse.Gene.ID', sep="")))
								Name1.homologene.mouse.UniProt = eval(parse(text=paste(Name1, '.', 'homologene.mouse.UniProt', sep="")))
								Name1.homologene.mouse.MGI.ID = eval(parse(text=paste(Name1, '.', 'homologene.mouse.MGI.ID', sep="")))
								Name1.homologene.mouse.MGI.Symbol = eval(parse(text=paste(Name1, '.', 'homologene.mouse.MGI.Symbol', sep="")))
								Name1.homologene.mouse.Gene.names = eval(parse(text=paste(Name1, '.', 'homologene.mouse.Gene.names', sep="")))
								Name1.homologene.mouse.Protein.names = eval(parse(text=paste(Name1, '.', 'homologene.mouse.Protein.names', sep="")))
								input1 = Name1.homologene.mouse.Gene.ID[which(Name1.homologene.mouse.Gene.ID != '')]
								input1.idx = which(Name1.homologene.mouse.Gene.ID != '')
								
								Name2 <- dataset[combn.idx[2,i]]
								Name2.homologene.mouse.Gene.ID = eval(parse(text=paste(Name2, '.', 'homologene.mouse.Gene.ID', sep="")))
								Name2.homologene.mouse.UniProt = eval(parse(text=paste(Name2, '.', 'homologene.mouse.UniProt', sep="")))
								Name2.homologene.mouse.MGI.ID = eval(parse(text=paste(Name2, '.', 'homologene.mouse.MGI.ID', sep="")))
								Name2.homologene.mouse.MGI.Symbol = eval(parse(text=paste(Name2, '.', 'homologene.mouse.MGI.Symbol', sep="")))
								Name2.homologene.mouse.Gene.names = eval(parse(text=paste(Name2, '.', 'homologene.mouse.Gene.names', sep="")))
								Name2.homologene.mouse.Protein.names = eval(parse(text=paste(Name2, '.', 'homologene.mouse.Protein.names', sep="")))
								input2 = Name2.homologene.mouse.Gene.ID[which(Name2.homologene.mouse.Gene.ID != '')]
								input2.idx = which(Name2.homologene.mouse.Gene.ID != '')
								
								Name3 <- dataset[combn.idx[3,i]]
								Name3.homologene.mouse.Gene.ID = eval(parse(text=paste(Name3, '.', 'homologene.mouse.Gene.ID', sep="")))
								Name3.homologene.mouse.UniProt = eval(parse(text=paste(Name3, '.', 'homologene.mouse.UniProt', sep="")))
								Name3.homologene.mouse.MGI.ID = eval(parse(text=paste(Name3, '.', 'homologene.mouse.MGI.ID', sep="")))
								Name3.homologene.mouse.MGI.Symbol = eval(parse(text=paste(Name3, '.', 'homologene.mouse.MGI.Symbol', sep="")))
								Name3.homologene.mouse.Gene.names = eval(parse(text=paste(Name3, '.', 'homologene.mouse.Gene.names', sep="")))
								Name3.homologene.mouse.Protein.names = eval(parse(text=paste(Name3, '.', 'homologene.mouse.Protein.names', sep="")))
								input3 = Name3.homologene.mouse.Gene.ID[which(Name3.homologene.mouse.Gene.ID != '')]
								input3.idx = which(Name3.homologene.mouse.Gene.ID != '')
								
								if(length(input1) == 0 & length(input2) == 0 & length(input3) == 0) {
									# Do nothing
								} else {
									if((Name1!='E8') & (Name2!='E8') & (Name3!='E8'))
									{
										mainDir <- paste('./Data/',Sample,'/10_datasetOverlap/Data/subset/n_',r,'/',folder[f],'/',n,'_dataset',sep="")
										subDir <- combnName
										if (file.exists(file.path(mainDir, subDir))){
										} else {
											dir.create(file.path(mainDir, subDir))
										}
										
										
										# === triple comparison ===============================================================================
										input.union <- Reduce(union, list(input1,input2,input3))
										input.intersect <- Reduce(intersect, list(input1,input2,input3))		#intersect(intersect(input1,input2),input3)
										#print(length(input.union))
										#print(length(input.intersect))
										
										input.intersect.12 <- intersect(input1,input2)
										input.intersect.13 <- intersect(input1,input3)
										input.intersect.23 <- intersect(input2,input3)
										#print(length(input.intersect.12))
										#print(length(input.intersect.13))
										#print(length(input.intersect.23))
										
										input.intersect.12_3 <- setdiff(input.intersect.12,input.intersect)
										input.intersect.13_2 <- setdiff(input.intersect.13,input.intersect)
										input.intersect.23_1 <- setdiff(input.intersect.23,input.intersect)
										#print(length(input.intersect.12_3))
										#print(length(input.intersect.13_2))
										#print(length(input.intersect.23_1))
										
										input.setdiff.1 <- Reduce(setdiff, list(input1,input2,input3))
										input.setdiff.2 <- Reduce(setdiff, list(input2,input1,input3))
										input.setdiff.3 <- Reduce(setdiff, list(input3,input1,input2))
										#print(length(input.setdiff.1))
										#print(length(input.setdiff.2))
										#print(length(input.setdiff.3))
										
										
										prename = paste('./Data/',Sample,'/10_datasetOverlap/Data/subset/n_',r,'/',folder[f],'/',n,'_dataset/',combnName,'/', combnName, '.exosome', sep="")
										
										wb <- openxlsx::createWorkbook(paste(prename,'.xlsx',sep=""))
										modifyBaseFont(wb, fontSize=12, fontColour="black", fontName="Times New Roman")
										
										overlapType <- c('union','intersect','intersect.12','intersect.13','intersect.23','intersect.12_3','intersect.13_2','intersect.23_1','setdiff.1','setdiff.2','setdiff.3')
										for(o in 1:length(overlapType))
										{
											input.overlapType = eval(parse(text=paste('input', '.', overlapType[o], sep="")))
											
											if(length(input.overlapType) > 0)
											{
												# --- txt -----------------------------------------------------------------------------------------
												write.table(input.overlapType, paste(prename,'_',overlapType[o],'.txt',sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE)
												
												# --- csv -----------------------------------------------------------------------------------------
												input.overlapType.Sample1.idx = sapply(input.overlapType, FUN=function(X) which(Name1.homologene.mouse.Gene.ID %in% X))
												input.overlapType.Sample2.idx = sapply(input.overlapType, FUN=function(X) which(Name2.homologene.mouse.Gene.ID %in% X))
												input.overlapType.Sample3.idx = sapply(input.overlapType, FUN=function(X) which(Name3.homologene.mouse.Gene.ID %in% X))
												uniprot_id.all <- c()
												mgi_id.all <- c()
												gene_symbol.all <- c()
												gene_name.all <- c()
												protein_name.all <- c()
												for(j in 1:length(input.overlapType))
												{
													Sample1.HitNumber = length(unlist(input.overlapType.Sample1.idx[j]))
													Sample2.HitNumber = length(unlist(input.overlapType.Sample2.idx[j]))
													Sample3.HitNumber = length(unlist(input.overlapType.Sample3.idx[j]))
													
													if(Sample1.HitNumber == 1) {
														uniprot_id = Name1.homologene.mouse.UniProt[unlist(input.overlapType.Sample1.idx[j])]
														mgi_id = Name1.homologene.mouse.MGI.ID[unlist(input.overlapType.Sample1.idx[j])]
														gene_symbol = Name1.homologene.mouse.MGI.Symbol[unlist(input.overlapType.Sample1.idx[j])]
														gene_name = Name1.homologene.mouse.Gene.names[unlist(input.overlapType.Sample1.idx[j])]
														protein_name = Name1.homologene.mouse.Protein.names[unlist(input.overlapType.Sample1.idx[j])]
													} else if(Sample1.HitNumber == 0) {
														if(Sample2.HitNumber == 1) {
															uniprot_id = Name2.homologene.mouse.UniProt[unlist(input.overlapType.Sample2.idx[j])]
															mgi_id = Name2.homologene.mouse.MGI.ID[unlist(input.overlapType.Sample2.idx[j])]
															gene_symbol = Name2.homologene.mouse.MGI.Symbol[unlist(input.overlapType.Sample2.idx[j])]
															gene_name = Name2.homologene.mouse.Gene.names[unlist(input.overlapType.Sample2.idx[j])]
															protein_name = Name2.homologene.mouse.Protein.names[unlist(input.overlapType.Sample2.idx[j])]
														} else if(Sample2.HitNumber == 0) {
															if(Sample3.HitNumber == 1) {
																uniprot_id = Name3.homologene.mouse.UniProt[unlist(input.overlapType.Sample3.idx[j])]
																mgi_id = Name3.homologene.mouse.MGI.ID[unlist(input.overlapType.Sample3.idx[j])]
																gene_symbol = Name3.homologene.mouse.MGI.Symbol[unlist(input.overlapType.Sample3.idx[j])]
																gene_name = Name3.homologene.mouse.Gene.names[unlist(input.overlapType.Sample3.idx[j])]
																protein_name = Name3.homologene.mouse.Protein.names[unlist(input.overlapType.Sample3.idx[j])]
															} else if(Sample3.HitNumber == 0) {
																uniprot_id = ''
																mgi_id = ''
																gene_symbol = ''
																gene_name = ''
																protein_name = ''
																cat('[Warning] mouse.Gene.ID:',input.overlapType[j],' does not have any id mapping to all homologene.mouse.Gene.ID.\n', sep="")
															} else if(Sample3.HitNumber > 1) {
																uniprot_id = paste(unique(Name3.homologene.mouse.UniProt[unlist(input.overlapType.Sample3.idx[j])]), collapse="|")
																mgi_id = paste(unique(Name3.homologene.mouse.MGI.ID[unlist(input.overlapType.Sample3.idx[j])]), collapse="|")
																gene_symbol = paste(unique(Name3.homologene.mouse.MGI.Symbol[unlist(input.overlapType.Sample3.idx[j])]), collapse="|")
																gene_name = paste(unique(Name3.homologene.mouse.Gene.names[unlist(input.overlapType.Sample3.idx[j])]), collapse="|")
																protein_name = paste(unique(Name3.homologene.mouse.Protein.names[unlist(input.overlapType.Sample3.idx[j])]), collapse="|")
																
																if(length(unique(Name3.homologene.mouse.UniProt[unlist(input.overlapType.Sample3.idx[j])])) > 1)
																cat('[Note] mouse.Gene.ID:',input.overlapType[j],' has more than two id mapping to \"Name3.homologene.mouse.Gene.ID (',dataset[combn.idx[3,i]],'.homologene.mouse.Gene.ID)\".\n', sep="")
															}
														} else if(Sample2.HitNumber > 1) {
															uniprot_id = paste(unique(Name2.homologene.mouse.UniProt[unlist(input.overlapType.Sample2.idx[j])]), collapse="|")
															mgi_id = paste(unique(Name2.homologene.mouse.MGI.ID[unlist(input.overlapType.Sample2.idx[j])]), collapse="|")
															gene_symbol = paste(unique(Name2.homologene.mouse.MGI.Symbol[unlist(input.overlapType.Sample2.idx[j])]), collapse="|")
															gene_name = paste(unique(Name2.homologene.mouse.Gene.names[unlist(input.overlapType.Sample2.idx[j])]), collapse="|")
															protein_name = paste(unique(Name2.homologene.mouse.Protein.names[unlist(input.overlapType.Sample2.idx[j])]), collapse="|")
															
															if(length(unique(Name2.homologene.mouse.UniProt[unlist(input.overlapType.Sample2.idx[j])])) > 1)
															cat('[Note] mouse.Gene.ID:',input.overlapType[j],' has more than two id mapping to \"Name2.homologene.mouse.Gene.ID (',dataset[combn.idx[2,i]],'.homologene.mouse.Gene.ID)\".\n', sep="")
														}
													} else if(Sample1.HitNumber > 1) {
														uniprot_id = paste(unique(Name1.homologene.mouse.UniProt[unlist(input.overlapType.Sample1.idx[j])]), collapse="|")
														mgi_id = paste(unique(Name1.homologene.mouse.MGI.ID[unlist(input.overlapType.Sample1.idx[j])]), collapse="|")
														gene_symbol = paste(unique(Name1.homologene.mouse.MGI.Symbol[unlist(input.overlapType.Sample1.idx[j])]), collapse="|")
														gene_name = paste(unique(Name1.homologene.mouse.Gene.names[unlist(input.overlapType.Sample1.idx[j])]), collapse="|")
														protein_name = paste(unique(Name1.homologene.mouse.Protein.names[unlist(input.overlapType.Sample1.idx[j])]), collapse="|")
														
														if(length(unique(Name1.homologene.mouse.UniProt[unlist(input.overlapType.Sample1.idx[j])])) > 1)
														cat('[Note] mouse.Gene.ID:',input.overlapType[j],' has more than two id mapping to \"Name1.homologene.mouse.Gene.ID (',dataset[combn.idx[1,i]],'.homologene.mouse.Gene.ID)\".\n', sep="")
													}
													
													uniprot_id.all = c(uniprot_id.all, uniprot_id)
													mgi_id.all = c(mgi_id.all, mgi_id)
													gene_symbol.all = c(gene_symbol.all, gene_symbol)
													gene_name.all = c(gene_name.all, gene_name)
													protein_name.all = c(protein_name.all, protein_name)
												}
												input.overlapType.info = cbind(input.overlapType, uniprot_id.all, mgi_id.all, gene_symbol.all, gene_name.all, protein_name.all)
												colnames(input.overlapType.info) = c('mouse.Gene.ID','mouse.UniProt','mouse.MGI.ID','mouse.MGI.Symbol','mouse.Gene.names','mouse.Protein.names')
												#write.table(input.overlapType.info, paste(prename,'_',overlapType[o],'.csv',sep=""), quote=TRUE, na="", sep=",", row.names=FALSE)
												
												# --- xlsx ----------------------------------------------------------------------------------------
												sheet = overlapType[o]
												sheetData = input.overlapType.info
												
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
													}
												}
											}
										}
										
										if(length(wb$sheet_names) > 0) {
											openxlsx::saveWorkbook(wb, file=paste(prename,'.xlsx',sep=""), overwrite=TRUE)
										}
										
										# --- png -----------------------------------------------------------------------------------------
										### Plot venn diagram of 3 datasets ###
										venn.plot <- draw.triple.venn(
										area1 = length(input1),			# The size of the first set
										area2 = length(input2),			# The size of the second set
										area3 = length(input3),			# The size of the third set
										
										n12 = length(input.intersect.12),		# The size of the intersection between the first and the second set
										n23 = length(input.intersect.23),		# The size of the intersection between the second and the third set
										n13 = length(input.intersect.13),		# The size of the intersection between the first and the third set
										n123 = length(input.intersect),			# The size of the intersection between all three sets
										
										category = c(paste(dataset[combn.idx[1,i]],'_exosome',' ',paste('(',length(input1),')',sep=""),sep=""), paste(dataset[combn.idx[2,i]],'_exosome',' ',paste('(',length(input2),')',sep=""),sep=""), paste(dataset[combn.idx[3,i]],'_exosome',' ',paste('(',length(input3),')',sep=""),sep="")),	# category name (number) of the sets
										cat.pos = c(-15,15,180),				# the positions (in degrees) of the category names along the circles
										cat.dist = c(0.03,0.03,0.02),			# the distances (in npc units) of the category names from the edges of the circles
										cat.cex = 3,							# the size of the category names
										cat.col = c("black","black","black"),	# the colours of the category names
										#fill = c("gray30","gray45","gray60"),	# the colours of the circles' areas
										fill = c("white","white","white"),		# the colours of the circles' areas
										#alpha = rep(1,3),						# the alpha transparency of the circles' areas
										#alpha = c(0.6,0.8,1),					# the alpha transparency of the circles' areas
										#lty = "blank",							# the line dash pattern of the circles' circumferences
										cex = 3									# the size of the areas' labels
										)
										
										png(paste(prename,'_vennDiagram.png',sep=""), width=14, height=14, units="in", res=600)
										grid.draw(venn.plot)
										dev.off()
										
										grid.draw(venn.plot);
										grid.newpage();							# clear plotting area for new plots
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

