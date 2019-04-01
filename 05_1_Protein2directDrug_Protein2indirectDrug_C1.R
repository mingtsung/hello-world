rm(list=ls())	# remove all variables in workspace

setwd("D:/Exo")

Sample <- "pcMSC_N230_0"

# /// pcMSC ///
if(Sample == "pcMSC_N230_0")				{	batch <- c('P3_3rd_181022','P4_4th_181022','P4_5th_181022')
}else if(Sample == "pcMSC_N230")			{	batch <- c('pcMSC-1_181102','pcMSC-2_181102')
}else if(Sample == "pcMSC_N226")			{	batch <- c('P14_B2-2_190103','P14_B3-1_190103','P14_B4_190103')
}else if(Sample == "pcMSC")					{	batch <- c('pcMSC-1_181102','pcMSC-2_181102','P14_B2-2_190103','P14_B3-1_190103','P14_B4_190103')

# /// pcMSC_Glucose ///
}else if(Sample == "pcMSC_Glucose_5mM_0")	{	batch <- c('5_190117')
}else if(Sample == "pcMSC_Glucose_25mM_0")	{	batch <- c('25_190117')
}else if(Sample == "pcMSC_Glucose_25mM")	{	batch <- c('Glucos_190312')

# /// pcMSC_TNFalpha ///
}else if(Sample == "pcMSC_TNFalpha")		{	batch <- c('TNFa_190312')

# /// SB ///
}else if(Sample == "SB")					{	batch <- c('SB_190103')

# /// Mascot ///
#}else if(Sample == "pcMSC_N230_Mascot")		{	batch <- c('20181115_F059869_pcMSC-1_181102_260proteins','20181115_F059870_pcMSC-2_181102_112proteins')
#}else if(Sample == "pcMSC_N226_Mascot")		{	batch <- c('20190108_F060331_P14_B2-2_190103_69proteins','20190108_F060332_P14_B3-1_190103_155proteins','20190108_F060333_P14_B4_190103_33proteins')
#}else if(Sample == "SB_Mascot")				{	batch <- c('20190108_F060334_SB_190103_316proteins')

}else										{	print(paste('Please provide filename of Sample ',Sample,sep=""))
}


mainDir <- paste('./Data/',Sample,sep="")
subDir <- "05_protein2drug"
if (file.exists(file.path(mainDir, subDir))){
} else {
	dir.create(file.path(mainDir, subDir))
}
mainDir <- paste('./Data/',Sample,'/05_protein2drug',sep="")
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


# === Cytoscape Attribute (.noa) ====================================================================================================================
# --- search for human id ---
human.idmapping <- read.table("./Database/UniProt/idmapping/by_organism/HUMAN_9606_idmapping_selected.tab", sep="\t", quote="", header=FALSE, stringsAsFactors=FALSE)
human.idmapping.UniProt <- human.idmapping[[1]]
human.idmapping.EntryName <- human.idmapping[[2]]
human.idmapping.GeneID <- human.idmapping[[3]]
human.idmapping.RefSeq <- human.idmapping[[4]]
human.idmapping.GI <- human.idmapping[[5]]
human.UniProt.attribute <- cbind(human.idmapping.UniProt,'human_UniProt')

# --- search for mouse id ---
mouse.idmapping <- read.table("./Database/UniProt/idmapping/by_organism/MOUSE_10090_idmapping_selected.tab", sep="\t", quote="", header=FALSE, stringsAsFactors=FALSE)
mouse.idmapping.UniProt <- mouse.idmapping[[1]]
mouse.idmapping.EntryName <- mouse.idmapping[[2]]
mouse.idmapping.GeneID <- mouse.idmapping[[3]]
mouse.idmapping.RefSeq <- mouse.idmapping[[4]]
mouse.idmapping.GI <- mouse.idmapping[[5]]
mouse.UniProt.attribute <- cbind(mouse.idmapping.UniProt,'mouse_UniProt')

UniProt.attribute <- rbind(human.UniProt.attribute,mouse.UniProt.attribute)

# --- search for Drug id ---
drug_links <- read.csv("./Database/DrugBank/External_Links/External_Drug_Links/drugbank_all_drug_links.csv", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
drug_links.DrugBankID <- trimws(drug_links[,'DrugBank ID'])
DrugBank.attribute <- cbind(drug_links.DrugBankID,'DrugBank')

attribute <- rbind(UniProt.attribute,DrugBank.attribute)
attribute = unique(attribute)
write.table(attribute,paste('./Data/',Sample,'/05_protein2drug/','UniProt','_','DrugBank','_','attribute.noa',sep=""),quote=FALSE,sep="\t",row.names=FALSE,col.names=c('ID','Type'))

#all_drug_links <- read.csv("./Database/DrugBank/External_Links/External_Drug_Links/drugbank_all_drug_links.csv", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
approved_drug_links <- read.csv("./Database/DrugBank/External_Links/External_Drug_Links/drugbank_approved_drug_links.csv", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
#small_molecule_drug_links <- read.csv("./Database/DrugBank/External_Links/External_Drug_Links/drugbank_small_molecule_drug_links.csv", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
#biotech_drug_links <- read.csv("./Database/DrugBank/External_Links/External_Drug_Links/drugbank_biotech_drug_links.csv", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
experimental_drug_links <- read.csv("./Database/DrugBank/External_Links/External_Drug_Links/drugbank_experimental_drug_links.csv", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
nutraceutical_drug_links <- read.csv("./Database/DrugBank/External_Links/External_Drug_Links/drugbank_nutraceutical_drug_links.csv", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
illicit_drug_links <- read.csv("./Database/DrugBank/External_Links/External_Drug_Links/drugbank_illicit_drug_links.csv", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
withdrawn_drug_links <- read.csv("./Database/DrugBank/External_Links/External_Drug_Links/drugbank_withdrawn_drug_links.csv", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
investigational_drug_links <- read.csv("./Database/DrugBank/External_Links/External_Drug_Links/drugbank_investigational_drug_links.csv", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)


# === subset ====================================================================================================================
folder <- c('high','low')

for(r in 1:length(batch))
{
	if (file.exists(paste('./Data/',Sample,'/04_drug/subset/n_',r,sep=""))) {
		
		mainDir <- paste('./Data/',Sample,'/05_protein2drug/subset',sep="")
		subDir <- paste('n_',r,sep="")
		if (file.exists(file.path(mainDir, subDir))){
		} else {
			dir.create(file.path(mainDir, subDir))
		}
		
		for(f in 1:length(folder))
		{
			if (file.exists(paste('./Data/',Sample,'/04_drug/subset/n_',r,'/',folder[f],sep=""))) {
				
				if (file.exists(paste('./Data/',Sample,'/04_drug/subset/n_',r,'/',folder[f],'/',Sample,'_exo_UniProt_drug.xlsx',sep=""))) {
					
					mainDir <- paste('./Data/',Sample,'/05_protein2drug/subset/n_',r,sep="")
					subDir <- folder[f]
					if (file.exists(file.path(mainDir, subDir))){
					} else {
						dir.create(file.path(mainDir, subDir))
					}
					
					
					# === input =======================================================================================
					#input <- read.csv(paste('./Data/',Sample,'/04_drug/subset/n_',r,'/',folder[f],'/',Sample,'_exo_UniProt_drug.csv',sep=""), header=TRUE, stringsAsFactors=FALSE)
					input <- read.xlsx(paste('./Data/',Sample,'/04_drug/subset/n_',r,'/',folder[f],'/',Sample,'_exo_UniProt_drug.xlsx',sep=""), sheet='UniProt_drug', colNames=TRUE, rowNames=FALSE, check.names=FALSE)
					human.UniProt <- as.character(input$human.UniProt)
					human.GI <- as.character(input$human.GI)
					human.Gene.ID <- as.character(input$human.Gene.ID)
					human.HGNC.ID <- as.character(input$human.HGNC.ID)
					human.HGNC.Symbol <- as.character(input$human.HGNC.Symbol)
					HID <- as.character(input$HID)
					mouse.MGI.Symbol <- as.character(input$mouse.MGI.Symbol)
					mouse.MGI.ID <- as.character(input$mouse.MGI.ID)
					mouse.Gene.ID <- as.character(input$mouse.Gene.ID)
					mouse.GI <- as.character(input$mouse.GI)
					mouse.UniProt <- as.character(input$mouse.UniProt)
					human.Gene.names <- as.character(input$human.Gene.names)
					human.Protein.names <- as.character(input$human.Protein.names)
					mouse.Gene.names <- as.character(input$mouse.Gene.names)
					mouse.Protein.names <- as.character(input$mouse.Protein.names)
					human.Frequency <- as.character(input$human.Frequency)
					human.DrugBank <- as.character(input$human.DrugBank)
					human.UniProt.interaction <- as.character(input$human.UniProt.interaction)
					human.UniProt.interaction.drug <- as.character(input$human.UniProt.interaction.drug)
					mouse.DrugBank <- as.character(input$mouse.DrugBank)
					mouse.UniProt.interaction <- as.character(input$mouse.UniProt.interaction)
					mouse.UniProt.interaction.drug <- as.character(input$mouse.UniProt.interaction.drug)
					
					
					# === protein information =========================================================================
					# --- order by frequency ---
					human.protein.list.frequency.order <- cbind(human.UniProt,mouse.UniProt,human.Protein.names,mouse.Protein.names,human.DrugBank,mouse.DrugBank,human.Frequency)
					#dim(human.protein.list.frequency.order)
					#write.table(human.protein.list.frequency.order[,c('human.UniProt','mouse.UniProt','human.Protein.names','mouse.Protein.names','human.Frequency')],paste('./Data/',Sample,'/05_protein2drug/subset/n_',r,'/',folder[f],'/',Sample,'_exo_protein.name.csv',sep=""),quote=TRUE,na="",sep=",",row.names=FALSE)
					#if(r == 1 & folder[f] == 'high')	write.table(human.protein.list.frequency.order[,c('human.UniProt','mouse.UniProt','human.Protein.names','mouse.Protein.names','human.Frequency')],paste('./Data/',Sample,'/05_protein2drug/',Sample,'_exo_protein.name.csv',sep=""),quote=TRUE,na="",sep=",",row.names=FALSE)
					
					# --- order by frequency + mouse homolog filter (with mouse homolog) ---
					human.protein.list.with.mouse.homolog <- human.protein.list.frequency.order[which(human.protein.list.frequency.order[,'mouse.UniProt']!=""),,drop=FALSE]
					#dim(human.protein.list.with.mouse.homolog)
					#write.table(human.protein.list.with.mouse.homolog[,c('human.UniProt','mouse.UniProt','human.Protein.names','mouse.Protein.names','human.Frequency')],paste('./Data/',Sample,'/05_protein2drug/subset/n_',r,'/',folder[f],'/',Sample,'_exo_protein.name_homolog.csv',sep=""),quote=TRUE,na="",sep=",",row.names=FALSE)
					#if(r == 1 & folder[f] == 'high')	write.table(human.protein.list.with.mouse.homolog[,c('human.UniProt','mouse.UniProt','human.Protein.names','mouse.Protein.names','human.Frequency')],paste('./Data/',Sample,'/05_protein2drug/',Sample,'_exo_protein.name_homolog.csv',sep=""),quote=TRUE,na="",sep=",",row.names=FALSE)
					
					# --- order by frequency + mouse homolog filter (without mouse homolog) ---
					human.protein.list.without.mouse.homolog <- human.protein.list.frequency.order[which(human.protein.list.frequency.order[,'mouse.UniProt']==""),,drop=FALSE]
					#dim(human.protein.list.without.mouse.homolog)
					#write.table(human.protein.list.without.mouse.homolog[,c('human.UniProt','mouse.UniProt','human.Protein.names','mouse.Protein.names','human.Frequency')],paste('./Data/',Sample,'/05_protein2drug/subset/n_',r,'/',folder[f],'/',Sample,'_exo_protein.name_no.homolog.csv',sep=""),quote=TRUE,na="",sep=",",row.names=FALSE)
					#if(r == 1 & folder[f] == 'high')	write.table(human.protein.list.without.mouse.homolog[,c('human.UniProt','mouse.UniProt','human.Protein.names','mouse.Protein.names','human.Frequency')],paste('./Data/',Sample,'/05_protein2drug/',Sample,'_exo_protein.name_no.homolog.csv',sep=""),quote=TRUE,na="",sep=",",row.names=FALSE)
					
					# --- write xlsx file by openxlsx package ---------------------------------------------------------
					# Importing a big xlsx file into R? https://stackoverflow.com/questions/19147884/importing-a-big-xlsx-file-into-r/43118530#43118530
					# Building R for Windows: https://cran.r-project.org/bin/windows/Rtools/
					#Sys.setenv("R_ZIPCMD" = "C:/Rtools/bin/zip.exe")	# path to zip.exe
					
					wb <- openxlsx::createWorkbook(paste('./Data/',Sample,'/05_protein2drug/subset/n_',r,'/',folder[f],'/',Sample,'_exo_protein.name.xlsx',sep=""))
					modifyBaseFont(wb, fontSize=12, fontColour="black", fontName="Times New Roman")
					
					wb2 <- openxlsx::createWorkbook(paste('./Data/',Sample,'/05_protein2drug/',Sample,'_exo_protein.name.xlsx',sep=""))
					modifyBaseFont(wb2, fontSize=12, fontColour="black", fontName="Times New Roman")
					
					for(s in 1:3)
					{
						if(s == 1) {
							sheet = 'protein name'
							sheetData = human.protein.list.frequency.order[,c('human.UniProt','mouse.UniProt','human.Protein.names','mouse.Protein.names','human.Frequency'),drop=FALSE]
						} else if(s == 2) {
							sheet = 'human with mouse homolog'
							sheetData = human.protein.list.with.mouse.homolog[,c('human.UniProt','mouse.UniProt','human.Protein.names','mouse.Protein.names','human.Frequency'),drop=FALSE]
						} else if(s == 3) {
							sheet = 'human without mouse homolog'
							sheetData = human.protein.list.without.mouse.homolog[,c('human.UniProt','mouse.UniProt','human.Protein.names','mouse.Protein.names','human.Frequency'),drop=FALSE]
						}
						
						if(!is.null(sheetData))
						{
							if(nrow(sheetData) > 0)
							{
								addWorksheet(wb, sheet)
								writeData(wb, sheet, sheetData, colNames=TRUE, rowNames=FALSE, keepNA=FALSE)
								setColWidths(wb, sheet, cols=c(1:ncol(sheetData)), widths=c(rep(15,2),rep(20,2),rep(15,ncol(sheetData)-(2+2))))
								setRowHeights(wb, sheet, rows=c(1:(1+nrow(sheetData))), heights=16.5)
								freezePane(wb, sheet, firstActiveRow=2)
								style <- createStyle(halign="left", valign="center")
								addStyle(wb, sheet, style, rows=c(1:(1+nrow(sheetData))), cols=c(1:ncol(sheetData)), gridExpand=TRUE)
								style <- createStyle(textDecoration="bold")
								addStyle(wb, sheet, style, rows=1, cols=c(1:ncol(sheetData)), stack=TRUE)
								
								right_align.idx = which(colnames(sheetData) %in% c('human.Frequency'))
								if(length(right_align.idx) > 0) {
									style <- createStyle(halign="right", valign="center")
									addStyle(wb, sheet, style, rows=c(2:(1+nrow(sheetData))), cols=right_align.idx, gridExpand=TRUE, stack=TRUE)
								}
								
								if(r == 1 & folder[f] == 'high')
								{
									addWorksheet(wb2, sheet)
									writeData(wb2, sheet, sheetData, colNames=TRUE, rowNames=FALSE, keepNA=FALSE)
									setColWidths(wb2, sheet, cols=c(1:ncol(sheetData)), widths=c(rep(15,2),rep(20,2),rep(15,ncol(sheetData)-(2+2))))
									setRowHeights(wb2, sheet, rows=c(1:(1+nrow(sheetData))), heights=16.5)
									freezePane(wb2, sheet, firstActiveRow=2)
									style <- createStyle(halign="left", valign="center")
									addStyle(wb2, sheet, style, rows=c(1:(1+nrow(sheetData))), cols=c(1:ncol(sheetData)), gridExpand=TRUE)
									style <- createStyle(textDecoration="bold")
									addStyle(wb2, sheet, style, rows=1, cols=c(1:ncol(sheetData)), stack=TRUE)
									
									right_align.idx = which(colnames(sheetData) %in% c('human.Frequency'))
									if(length(right_align.idx) > 0) {
										style <- createStyle(halign="right", valign="center")
										addStyle(wb2, sheet, style, rows=c(2:(1+nrow(sheetData))), cols=right_align.idx, gridExpand=TRUE, stack=TRUE)
									}
								}
							}
						}
					}
					
					if(length(wb$sheet_names) > 0) {
						openxlsx::saveWorkbook(wb, file=paste('./Data/',Sample,'/05_protein2drug/subset/n_',r,'/',folder[f],'/',Sample,'_exo_protein.name.xlsx',sep=""), overwrite=TRUE)
					}
					
					if(r == 1 & folder[f] == 'high')
					{
						if(length(wb2$sheet_names) > 0) {
							openxlsx::saveWorkbook(wb2, file=paste('./Data/',Sample,'/05_protein2drug/',Sample,'_exo_protein.name.xlsx',sep=""), overwrite=TRUE)
						}
					}
					
					
					# --- order by frequency + mouse homolog filter (with mouse homolog) + human/mouse drug filter (with drug) ---
					#mouse.homolog.idx <- which(human.protein.list.frequency.order[,'mouse.UniProt']!="")
					#mouse.homolog.withdrug.idx <- which((human.protein.list.frequency.order[,'mouse.UniProt']!="") & (human.protein.list.frequency.order[,'mouse.DrugBank']!=""))
					#mouse.homolog.withoutdrug.idx <- which((human.protein.list.frequency.order[,'mouse.UniProt']!="") & (human.protein.list.frequency.order[,'mouse.DrugBank']==""))
					#mouse.homolog.idx
					#length(mouse.homolog.idx)
					#mouse.homolog.withdrug.idx
					#length(mouse.homolog.withdrug.idx)
					#mouse.homolog.withoutdrug.idx
					#length(mouse.homolog.withoutdrug.idx)
					#if((length(mouse.homolog.withdrug.idx) + length(mouse.homolog.withoutdrug.idx)) == length(mouse.homolog.idx))
					#{
					#	print('Please check mouse homolog protein list.')
					#	print(paste('The number of mouse homolog protein is ',length(mouse.homolog.idx),sep=""))
					#	print(paste('The number of mouse homolog protein with drug is ',length(mouse.homolog.withdrug.idx),sep=""))
					#	print(paste('The number of mouse homolog protein without drug is ',length(mouse.homolog.withoutdrug.idx),sep=""))
					#}
					
					#sum(which(human.protein.list.frequency.order[,'mouse.DrugBank']!="") %in% which(human.protein.list.frequency.order[,'mouse.UniProt']!="")) == length(which(human.protein.list.frequency.order[,'mouse.DrugBank']!=""))
					human.drug.idx <- which(((!is.na(human.protein.list.frequency.order[,'human.DrugBank'])) & (human.protein.list.frequency.order[,'human.DrugBank']!="")))
					mouse.drug.idx <- which(((!is.na(human.protein.list.frequency.order[,'mouse.UniProt'])) & (human.protein.list.frequency.order[,'mouse.UniProt']!="")) & ((!is.na(human.protein.list.frequency.order[,'mouse.DrugBank'])) & (human.protein.list.frequency.order[,'mouse.DrugBank']!="")))
					drug.idx <- unique(c(human.drug.idx,mouse.drug.idx))
					#length(human.drug.idx)
					#length(mouse.drug.idx)
					#length(drug.idx)
					
					human.drug.protein.list.frequency.order <- human.protein.list.frequency.order[human.drug.idx,,drop=FALSE]
					mouse.drug.protein.list.frequency.order <- human.protein.list.frequency.order[mouse.drug.idx,,drop=FALSE]
					drug.protein.list.frequency.order <- human.protein.list.frequency.order[drug.idx,,drop=FALSE]
					#dim(human.drug.protein.list.frequency.order)
					#dim(mouse.drug.protein.list.frequency.order)
					#dim(drug.protein.list.frequency.order)
					
					if(nrow(drug.protein.list.frequency.order) > 0)
					{
						# --- order by drug number ---
						#length(drug.protein.list.frequency.order[,'human.DrugBank'])
						human.withDrug.idx <- c()
						human.drug.num <- c()
						for(j in 1:length(drug.protein.list.frequency.order[,'human.DrugBank']))
						{
							if((!is.na(drug.protein.list.frequency.order[,'human.DrugBank'][j])) & (drug.protein.list.frequency.order[,'human.DrugBank'][j]!="")){
								human.withDrug.idx <- c(human.withDrug.idx,j)
								human.drug.num <- c(human.drug.num,length(strsplit(drug.protein.list.frequency.order[,'human.DrugBank'],";")[[j]]))
							}else{
								human.drug.num <- c(human.drug.num,0)
							}	
						}
						#human.withDrug.idx
						#length(human.withDrug.idx)
						#human.drug.num
						#length(human.drug.num)
						sorted.human.drug.num <- sort(human.drug.num,decreasing=TRUE)
						sorted.human.drug.num.order <- order(human.drug.num,decreasing=TRUE)
						#sorted.human.drug.num
						#sorted.human.drug.num.order
						
						#length(drug.protein.list.frequency.order[,'mouse.DrugBank'])
						mouse.withDrug.idx <- c()
						mouse.drug.num <- c()
						for(j in 1:length(drug.protein.list.frequency.order[,'mouse.DrugBank']))
						{
							if((!is.na(drug.protein.list.frequency.order[,'mouse.DrugBank'][j])) & (drug.protein.list.frequency.order[,'mouse.DrugBank'][j]!="")){
								mouse.withDrug.idx <- c(mouse.withDrug.idx,j)
								mouse.drug.num <- c(mouse.drug.num,length(strsplit(drug.protein.list.frequency.order[,'mouse.DrugBank'],";")[[j]]))
							}else{
								mouse.drug.num <- c(mouse.drug.num,0)
							}	
						}
						#mouse.withDrug.idx
						#length(mouse.withDrug.idx)
						#mouse.drug.num
						#length(mouse.drug.num)
						sorted.mouse.drug.num <- sort(mouse.drug.num,decreasing=TRUE)
						sorted.mouse.drug.num.order <- order(mouse.drug.num,decreasing=TRUE)
						#sorted.mouse.drug.num
						#sorted.mouse.drug.num.order
						
						drug.protein.list.drugnum.order <- c()
						if(exists("human.withDrug.idx"))
						{
							if(length(human.withDrug.idx) > 0)
							{
								human.drug.protein.list.drugnum.order <- c()
								for(j in 1:length(which(sorted.human.drug.num.order %in% human.withDrug.idx)))
								{
									human.drug.protein.list.drugnum.order <- rbind(human.drug.protein.list.drugnum.order,drug.protein.list.frequency.order[sorted.human.drug.num.order[j],])
									drug.protein.list.drugnum.order <- rbind(drug.protein.list.drugnum.order,drug.protein.list.frequency.order[sorted.human.drug.num.order[j],])
								}
								#print(human.drug.protein.list.drugnum.order)
							}
						}
						if(exists("mouse.withDrug.idx"))
						{
							if(length(mouse.withDrug.idx) > 0)
							{
								mouse.drug.protein.list.drugnum.order <- c()
								for(j in 1:length(which(sorted.mouse.drug.num.order %in% mouse.withDrug.idx)))
								{
									mouse.drug.protein.list.drugnum.order <- rbind(mouse.drug.protein.list.drugnum.order,drug.protein.list.frequency.order[sorted.mouse.drug.num.order[j],])
								}
								#print(mouse.drug.protein.list.drugnum.order)
								
								for(j in 1:length(sorted.mouse.drug.num.order))
								{
									if(sorted.mouse.drug.num.order[j] %in% mouse.withDrug.idx[which(!mouse.withDrug.idx %in% human.withDrug.idx)])
									{
										drug.protein.list.drugnum.order <- rbind(drug.protein.list.drugnum.order,drug.protein.list.frequency.order[sorted.mouse.drug.num.order[j],])
									}
								}
							}
						}
						#drug.protein.list.drugnum.order
						#dim(drug.protein.list.drugnum.order)
						
						order_name <- c('Frequency','DrugNum')
						for(o in 1:length(order_name))
						{
							if(order_name[o] == 'Frequency') {
								drug.protein.list.order = drug.protein.list.frequency.order
							} else if(order_name[o] == 'DrugNum') {
								drug.protein.list.order = drug.protein.list.drugnum.order
							}
							
							# --- remove duplicated mouse protein id ---
							duplicated.drug.mouseprotein.idx <- which(duplicated(drug.protein.list.order[,'mouse.UniProt']))
							nonduplicated.drug.mouseprotein.idx <- which(!duplicated(drug.protein.list.order[,'mouse.UniProt']))
							#duplicated.drug.mouseprotein.idx
							#nonduplicated.drug.mouseprotein.idx
							
							duplicated.drug.mouseprotein.id <- unique(drug.protein.list.order[,'mouse.UniProt'][which(duplicated(drug.protein.list.order[,'mouse.UniProt']))])
							nonduplicated.drug.mouseprotein.id <- drug.protein.list.order[,'mouse.UniProt'][which(!duplicated(drug.protein.list.order[,'mouse.UniProt']))]
							#duplicated.drug.mouseprotein.id
							#nonduplicated.drug.mouseprotein.id
							#length(duplicated.drug.mouseprotein.id)
							#length(nonduplicated.drug.mouseprotein.id)
							
							drug.uniqueprotein.list.order <- c()
							drug.uniqueprotein.list.order.withDrug <- c()
							for(j in 1:dim(drug.protein.list.order)[1])
							{
								if(drug.protein.list.order[j,'mouse.UniProt'] %in% duplicated.drug.mouseprotein.id) {
									nonduplicated.drug.humanprotein.id <- drug.protein.list.order[which(drug.protein.list.order[,'mouse.UniProt'] %in% drug.protein.list.order[j,'mouse.UniProt'])[1],'human.UniProt']
									nonduplicated.drug.humanprotein.name <- drug.protein.list.order[which(drug.protein.list.order[,'mouse.UniProt'] %in% drug.protein.list.order[j,'mouse.UniProt'])[1],'human.Protein.names']
									if(length(which(drug.protein.list.order[,'mouse.UniProt'] %in% drug.protein.list.order[j,'mouse.UniProt'])) >= 2)
									{
										for(k in 2:length(which(drug.protein.list.order[,'mouse.UniProt'] %in% drug.protein.list.order[j,'mouse.UniProt'])))
										{
											nonduplicated.drug.humanprotein.id <- paste(nonduplicated.drug.humanprotein.id, drug.protein.list.order[which(drug.protein.list.order[,'mouse.UniProt'] %in% drug.protein.list.order[j,'mouse.UniProt'])[k],'human.UniProt'], sep=' / ')
											nonduplicated.drug.humanprotein.name <- paste(nonduplicated.drug.humanprotein.name, drug.protein.list.order[which(drug.protein.list.order[,'mouse.UniProt'] %in% drug.protein.list.order[j,'mouse.UniProt'])[k],'human.Protein.names'], sep=' / ')
										}
									}
									drug.uniqueprotein.list.order <- rbind(drug.uniqueprotein.list.order, cbind(nonduplicated.drug.humanprotein.id, drug.protein.list.order[j,'mouse.UniProt'], nonduplicated.drug.humanprotein.name, drug.protein.list.order[j,'mouse.Protein.names'], drug.protein.list.order[j,'human.Frequency']))
									drug.uniqueprotein.list.order.withDrug <- rbind(drug.uniqueprotein.list.order.withDrug, cbind(nonduplicated.drug.humanprotein.id, drug.protein.list.order[j,'mouse.UniProt'], nonduplicated.drug.humanprotein.name, drug.protein.list.order[j,'mouse.Protein.names'], drug.protein.list.order[j,'human.DrugBank'], drug.protein.list.order[j,'mouse.DrugBank'], drug.protein.list.order[j,'human.Frequency']))
								} else {
									drug.uniqueprotein.list.order <- rbind(drug.uniqueprotein.list.order, drug.protein.list.order[j,c('human.UniProt','mouse.UniProt','human.Protein.names','mouse.Protein.names','human.Frequency')])
									drug.uniqueprotein.list.order.withDrug <- rbind(drug.uniqueprotein.list.order.withDrug, drug.protein.list.order[j,])
								}
							}
							drug.uniqueprotein.list.order <- drug.uniqueprotein.list.order[nonduplicated.drug.mouseprotein.idx,,drop=FALSE]
							drug.uniqueprotein.list.order.withDrug <- drug.uniqueprotein.list.order.withDrug[nonduplicated.drug.mouseprotein.idx,,drop=FALSE]
							colnames(drug.uniqueprotein.list.order) = c('human.UniProt','mouse.UniProt','human.Protein.names','mouse.Protein.names','human.Frequency')
							colnames(drug.uniqueprotein.list.order.withDrug) = c('human.UniProt','mouse.UniProt','human.Protein.names','mouse.Protein.names','human.DrugBank','mouse.DrugBank','human.Frequency')
							
							human.DrugBank.number = sapply(strsplit(drug.uniqueprotein.list.order.withDrug[,'human.DrugBank'],';'), FUN=function(X) unique(ifelse(length(unlist(X))==0, 0, length(unlist(X)))))
							mouse.DrugBank.number = sapply(strsplit(drug.uniqueprotein.list.order.withDrug[,'mouse.DrugBank'],';'), FUN=function(X) unique(ifelse(length(unlist(X))==0, 0, length(unlist(X)))))
							
							drug.uniqueprotein.list.order.withDrug = cbind(drug.uniqueprotein.list.order.withDrug, human.DrugBank.number, mouse.DrugBank.number)
							drug.uniqueprotein.list.order.withDrug = drug.uniqueprotein.list.order.withDrug[,c('human.Frequency','human.UniProt','human.Protein.names','human.DrugBank','human.DrugBank.number','mouse.UniProt','mouse.Protein.names','mouse.DrugBank','mouse.DrugBank.number'),drop=FALSE]
							
							if(order_name[o] == 'Frequency') {
								drug.uniqueprotein.list.order.withDrug.reorder_idx = order(as.integer(drug.uniqueprotein.list.order.withDrug[,'human.Frequency']),as.integer(drug.uniqueprotein.list.order.withDrug[,'human.DrugBank.number']),as.integer(drug.uniqueprotein.list.order.withDrug[,'mouse.DrugBank.number']),decreasing=TRUE)
							} else if(order_name[o] == 'DrugNum') {
								drug.uniqueprotein.list.order.withDrug.reorder_idx = order(as.integer(drug.uniqueprotein.list.order.withDrug[,'human.DrugBank.number']),as.integer(drug.uniqueprotein.list.order.withDrug[,'mouse.DrugBank.number']),as.integer(drug.uniqueprotein.list.order.withDrug[,'human.Frequency']),decreasing=TRUE)
							}
							
							drug.uniqueprotein.list.order.withDrug.reorder = drug.uniqueprotein.list.order.withDrug[drug.uniqueprotein.list.order.withDrug.reorder_idx,,drop=FALSE]
							
							#write.table(drug.uniqueprotein.list.order.withDrug.reorder,paste('./Data/',Sample,'/05_protein2drug/subset/n_',r,'/',folder[f],'/',Sample,'_exo_protein2drug_by',order_name[o],'.csv',sep=""),quote=TRUE,na="",sep=",",row.names=FALSE)
							#if(r == 1 & folder[f] == 'high')	write.table(drug.uniqueprotein.list.order.withDrug.reorder,paste('./Data/',Sample,'/05_protein2drug/',Sample,'_exo_protein2drug_by',order_name[o],'.csv',sep=""),quote=TRUE,na="",sep=",",row.names=FALSE)
							
							
							species <- c('human','mouse')
							
							for(sp in 1:length(species))
							{
								protein2drug.submatrix.new.all <- c()
								
								for(j in 1:dim(drug.uniqueprotein.list.order.withDrug.reorder)[1])
								{
									if(as.numeric(drug.uniqueprotein.list.order.withDrug.reorder[j,paste(species[sp],'.','DrugBank.number',sep="")]) > 0)
									{
										DrugBank_id = unlist(strsplit(drug.uniqueprotein.list.order.withDrug.reorder[j,paste(species[sp],'.','DrugBank',sep="")],';'))
										
										protein2drug.submatrix <- c()
										for(d in 1:length(DrugBank_id))
										{
											drug_id = DrugBank_id[d]
											drug_links.idx = which(drug_links[,'DrugBank ID'] %in% drug_id)
											if(length(drug_links.idx) > 0) {
												drug_links.info = drug_links[drug_links.idx,,drop=FALSE]
												drug_links.info = data.frame(lapply(drug_links.info, trimws), check.names=FALSE, stringsAsFactors=FALSE)
												if(length(drug_links.idx) != 1) {
													cat('[Note] "',drug_id,'" has multiple records',' (',length(drug_links.idx),') in "./Database/DrugBank/External_Links/External_Drug_Links/drugbank_all_drug_links.csv"','.\n', sep="")
												}
											} else {
												drug_links.info = matrix(data='', nrow=1, ncol=ncol(drug_links), byrow=TRUE)
												colnames(drug_links.info) = colnames(drug_links)
												drug_links.info[,'DrugBank ID'] = DrugBank_id[d]
											}
											
											Drug_Group = c('Approved','Experimental','Nutraceutical','Illicit','Withdrawn','Investigational')
											Drug_Group.drug_links <- c()
											for(dg in 1:length(Drug_Group))
											{
												assign(paste(tolower(Drug_Group[dg]),'_','drug_links.idx',sep=""), which(eval(parse(text=paste(tolower(Drug_Group[dg]),'_','drug_links',sep="")))[,'DrugBank ID'] %in% drug_id))
												drug_links.idx = eval(parse(text=paste(tolower(Drug_Group[dg]),'_','drug_links.idx',sep="")))
												
												if(length(drug_links.idx) > 0) {
													Drug_Group.drug_links <- c(Drug_Group.drug_links, Drug_Group[dg])
												}
											}
											if(length(Drug_Group.drug_links) > 0) {
												Drug_Group.drug_links = paste(Drug_Group.drug_links, collapse=", ")
											} else {
												Drug_Group.drug_links = ''
											}
											
											protein2drug.submatrix <- rbind(protein2drug.submatrix, cbind(drug_id, drug_links.info, Drug_Group.drug_links))
										}
										if(nrow(protein2drug.submatrix) != as.numeric(drug.uniqueprotein.list.order.withDrug.reorder[j,paste(species[sp],'.','DrugBank.number',sep="")])) {
											cat('[Warning] Drug number of "protein2drug.submatrix" of "',drug.uniqueprotein.list.order.withDrug.reorder[j,paste(species[sp],'.','UniProt',sep="")],'" (',nrow(protein2drug.submatrix),') is not equal to its "',paste(species[sp],'.','DrugBank.number',sep=""),'" (',as.numeric(drug.uniqueprotein.list.order.withDrug.reorder[j,paste(species[sp],'.','DrugBank.number',sep="")]),')','.\n', sep="")
										}
										
										colnames(protein2drug.submatrix)[ncol(protein2drug.submatrix)] = 'Drug Group'
										protein2drug.submatrix = protein2drug.submatrix[,-which(colnames(protein2drug.submatrix)=='DrugBank ID'),drop=FALSE]
										colnames(protein2drug.submatrix)[1] = 'DrugBank ID'
										rownames(protein2drug.submatrix) = c(1:nrow(protein2drug.submatrix))
										
										protein2drug.submatrix = protein2drug.submatrix[,c(c('DrugBank ID','Name','Drug Type','Drug Group','CAS Number'),setdiff(colnames(protein2drug.submatrix),c('DrugBank ID','Name','Drug Type','Drug Group','CAS Number'))),drop=FALSE]
										
										if(!is.null(protein2drug.submatrix))
										{
											if(nrow(protein2drug.submatrix) > 0)
											{
												na.idx = which(is.na(protein2drug.submatrix), arr.ind=TRUE)
												if(nrow(na.idx) > 0)
												{
													for(ri in 1:nrow(na.idx))
													{
														protein2drug.submatrix[na.idx[ri,'row'],na.idx[ri,'col']] = ''
													}
												}
											}
										}
										
										if(!is.null(protein2drug.submatrix))
										{
											if(nrow(protein2drug.submatrix) > 0)
											{
												Title.INFO = as.character(drug.uniqueprotein.list.order.withDrug.reorder[j,c(paste(species[sp],'.','UniProt',sep=""),paste(species[sp],'.','Protein.names',sep=""))])
												Title.INFO = Title.INFO[which(!is.na(Title.INFO) & Title.INFO != '')]
												Title.INFO = paste(Title.INFO, collapse="||")
												
												Title = c(paste('> ', Title.INFO, sep=""), rep(NA,ncol(protein2drug.submatrix)))
												Header = c(NA, colnames(protein2drug.submatrix))
												SubMatrix = as.matrix(cbind(rownames(protein2drug.submatrix), protein2drug.submatrix))
												BlankLine = c(NA, rep(NA,ncol(protein2drug.submatrix)))
												protein2drug.submatrix.new = rbind(Title, Header, SubMatrix, BlankLine)
												
												protein2drug.submatrix.new.all <- rbind(protein2drug.submatrix.new.all, protein2drug.submatrix.new)
											}
										}
									}
								}
								
								assign(paste(species[sp],'protein2drug.submatrix.new.all',sep=""), protein2drug.submatrix.new.all)
							}
							
							
							# --- write xlsx file by openxlsx package ---------------------------------------------------------
							# Importing a big xlsx file into R? https://stackoverflow.com/questions/19147884/importing-a-big-xlsx-file-into-r/43118530#43118530
							# Building R for Windows: https://cran.r-project.org/bin/windows/Rtools/
							#Sys.setenv("R_ZIPCMD" = "C:/Rtools/bin/zip.exe")	# path to zip.exe
							
							if(order_name[o] == 'Frequency') {
								run = 1
							} else if(order_name[o] == 'DrugNum') {
								run = 2
							}
							
							for(s0 in 1:run)
							{
								if(s0 == 1) {
									wb <- openxlsx::createWorkbook(paste('./Data/',Sample,'/05_protein2drug/subset/n_',r,'/',folder[f],'/',Sample,'_exo_protein2drug_by',order_name[o],'.xlsx',sep=""))
									modifyBaseFont(wb, fontSize=12, fontColour="black", fontName="Times New Roman")
									
									wb2 <- openxlsx::createWorkbook(paste('./Data/',Sample,'/05_protein2drug/',Sample,'_exo_protein2drug_by',order_name[o],'.xlsx',sep=""))
									modifyBaseFont(wb2, fontSize=12, fontColour="black", fontName="Times New Roman")
								} else if(s0 == 2) {
									wb <- openxlsx::createWorkbook(paste('./Data/',Sample,'/05_protein2drug/subset/n_',r,'/',folder[f],'/',Sample,'_exo_protein2drug','.xlsx',sep=""))
									modifyBaseFont(wb, fontSize=12, fontColour="black", fontName="Times New Roman")
									
									wb2 <- openxlsx::createWorkbook(paste('./Data/',Sample,'/05_protein2drug/',Sample,'_exo_protein2drug','.xlsx',sep=""))
									modifyBaseFont(wb2, fontSize=12, fontColour="black", fontName="Times New Roman")
								}
								
								for(s in 1:3)
								{
									if(s == 1) {
										sheet = 'protein2drug'
										if(s0 == 1) {
											sheetData = drug.uniqueprotein.list.order.withDrug.reorder
										} else if(s0 == 2) {
											sheetData = drug.uniqueprotein.list.order.withDrug.reorder[,-which(colnames(drug.uniqueprotein.list.order.withDrug.reorder)=='human.Frequency'),drop=FALSE]
										}
									} else if(s == 2) {
										sheet = 'human protein2drug'
										sheetData = humanprotein2drug.submatrix.new.all
									} else if(s == 3) {
										sheet = 'mouse protein2drug'
										sheetData = mouseprotein2drug.submatrix.new.all
									}
									
									if(!is.null(sheetData))
									{
										if(nrow(sheetData) > 0)
										{
											if(s == 1) {
												addWorksheet(wb, sheet)
												writeData(wb, sheet, sheetData, colNames=TRUE, rowNames=FALSE, keepNA=FALSE)
												setColWidths(wb, sheet, cols=c(1:ncol(sheetData)), widths=c(rep(15,ncol(sheetData)-(4*2)),c(15,20,15,15),c(15,20,15,15)))
												setRowHeights(wb, sheet, rows=c(1:(1+nrow(sheetData))), heights=16.5)
												freezePane(wb, sheet, firstActiveRow=2)
												style <- createStyle(halign="left", valign="center")
												addStyle(wb, sheet, style, rows=c(1:(1+nrow(sheetData))), cols=c(1:ncol(sheetData)), gridExpand=TRUE)
												style <- createStyle(textDecoration="bold")
												addStyle(wb, sheet, style, rows=1, cols=c(1:ncol(sheetData)), stack=TRUE)
												
												right_align.idx = which(colnames(sheetData) %in% c('human.Frequency','human.DrugBank.number','mouse.DrugBank.number'))
												if(length(right_align.idx) > 0) {
													style <- createStyle(halign="right", valign="center")
													addStyle(wb, sheet, style, rows=c(2:(1+nrow(sheetData))), cols=right_align.idx, gridExpand=TRUE, stack=TRUE)
												}
											} else if(s == 2 | s == 3) {
												addWorksheet(wb, sheet)
												writeData(wb, sheet, sheetData, colNames=FALSE, rowNames=FALSE, keepNA=FALSE)
												setColWidths(wb, sheet, cols=c(1:ncol(sheetData)), widths=c(8.5,12.5,c(30,17.5,22.5),15,rep(17.5,2),rep(20,2),rep(12.5,ncol(sheetData)-(1+1+3+1+2+2))))
												setRowHeights(wb, sheet, rows=c(1:nrow(sheetData)), heights=16.5)
												freezePane(wb, sheet, firstActiveCol=6)
												style <- createStyle(halign="left", valign="center")
												addStyle(wb, sheet, style, rows=c(1:nrow(sheetData)), cols=c(1:ncol(sheetData)), gridExpand=TRUE)
											}
											
											if(r == 1 & folder[f] == 'high')
											{
												if(s == 1) {
													addWorksheet(wb2, sheet)
													writeData(wb2, sheet, sheetData, colNames=TRUE, rowNames=FALSE, keepNA=FALSE)
													setColWidths(wb2, sheet, cols=c(1:ncol(sheetData)), widths=c(rep(15,ncol(sheetData)-(4*2)),c(15,20,15,15),c(15,20,15,15)))
													setRowHeights(wb2, sheet, rows=c(1:(1+nrow(sheetData))), heights=16.5)
													freezePane(wb2, sheet, firstActiveRow=2)
													style <- createStyle(halign="left", valign="center")
													addStyle(wb2, sheet, style, rows=c(1:(1+nrow(sheetData))), cols=c(1:ncol(sheetData)), gridExpand=TRUE)
													style <- createStyle(textDecoration="bold")
													addStyle(wb2, sheet, style, rows=1, cols=c(1:ncol(sheetData)), stack=TRUE)
													
													right_align.idx = which(colnames(sheetData) %in% c('human.Frequency','human.DrugBank.number','mouse.DrugBank.number'))
													if(length(right_align.idx) > 0) {
														style <- createStyle(halign="right", valign="center")
														addStyle(wb2, sheet, style, rows=c(2:(1+nrow(sheetData))), cols=right_align.idx, gridExpand=TRUE, stack=TRUE)
													}
												} else if(s == 2 | s == 3) {
													addWorksheet(wb2, sheet)
													writeData(wb2, sheet, sheetData, colNames=FALSE, rowNames=FALSE, keepNA=FALSE)
													setColWidths(wb2, sheet, cols=c(1:ncol(sheetData)), widths=c(8.5,12.5,c(30,17.5,22.5),15,rep(17.5,2),rep(20,2),rep(12.5,ncol(sheetData)-(1+1+3+1+2+2))))
													setRowHeights(wb2, sheet, rows=c(1:nrow(sheetData)), heights=16.5)
													freezePane(wb2, sheet, firstActiveCol=6)
													style <- createStyle(halign="left", valign="center")
													addStyle(wb2, sheet, style, rows=c(1:nrow(sheetData)), cols=c(1:ncol(sheetData)), gridExpand=TRUE)
												}
											}
										}
									}
								}
								
								if(s0 == 1) {
									if(length(wb$sheet_names) > 0) {
										openxlsx::saveWorkbook(wb, file=paste('./Data/',Sample,'/05_protein2drug/subset/n_',r,'/',folder[f],'/',Sample,'_exo_protein2drug_by',order_name[o],'.xlsx',sep=""), overwrite=TRUE)
									}
									
									if(r == 1 & folder[f] == 'high')
									{
										if(length(wb2$sheet_names) > 0) {
											openxlsx::saveWorkbook(wb2, file=paste('./Data/',Sample,'/05_protein2drug/',Sample,'_exo_protein2drug_by',order_name[o],'.xlsx',sep=""), overwrite=TRUE)
										}
									}
								} else if(s0 == 2) {
									if(length(wb$sheet_names) > 0) {
										openxlsx::saveWorkbook(wb, file=paste('./Data/',Sample,'/05_protein2drug/subset/n_',r,'/',folder[f],'/',Sample,'_exo_protein2drug','.xlsx',sep=""), overwrite=TRUE)
									}
									
									if(r == 1 & folder[f] == 'high')
									{
										if(length(wb2$sheet_names) > 0) {
											openxlsx::saveWorkbook(wb2, file=paste('./Data/',Sample,'/05_protein2drug/',Sample,'_exo_protein2drug','.xlsx',sep=""), overwrite=TRUE)
										}
									}
								}
							}
						}
					}
					
					
					# === interaction =================================================================================
					write.table(attribute,paste('./Data/',Sample,'/05_protein2drug/subset/n_',r,'/',folder[f],'/','UniProt','_','DrugBank','_','attribute.noa',sep=""),quote=FALSE,sep="\t",row.names=FALSE,col.names=c('ID','Type'))
					
					human.UniProt = sapply(strsplit(human.UniProt,'\\|'), FUN=function(X) X[1])
					human.UniProt[which(is.na(human.UniProt))] = ''
					mouse.UniProt = sapply(strsplit(mouse.UniProt,'\\|'), FUN=function(X) X[1])
					mouse.UniProt[which(is.na(mouse.UniProt))] = ''
					
					species <- c('human','mouse')
					
					for(sp in 1:length(species))
					{
						if(species[sp] == 'human') {
							species.UniProt = human.UniProt
							species.DrugBank = human.DrugBank
							species.UniProt.interaction = human.UniProt.interaction
							species.UniProt.interaction.drug = human.UniProt.interaction.drug
						} else if(species[sp] == 'mouse') {
							species.UniProt = mouse.UniProt
							species.DrugBank = mouse.DrugBank
							species.UniProt.interaction = mouse.UniProt.interaction
							species.UniProt.interaction.drug = mouse.UniProt.interaction.drug
						}
						
						
						# --- human UniProt -> mouse UniProt (human homolog) / mouse UniProt -> human UniProt (mouse homolog) ---
						species.homologene.list <- c()
						for(j in 1:dim(input)[1])
						{
							#print(j)
							if(human.UniProt[j]!="")
							{
								if(mouse.UniProt[j]!="")
								{
									if(species[sp] == 'human') {
										species.homologene <- cbind(human.UniProt[j],'homologene',mouse.UniProt[j])
									} else if(species[sp] == 'mouse') {
										species.homologene <- cbind(mouse.UniProt[j],'homologene',human.UniProt[j])
									}
									#print(species.homologene)
									
									species.homologene.list <- rbind(species.homologene.list,species.homologene)
								}
							}
						}
						#print(species.homologene.list)
						species.homologene.list = unique(species.homologene.list)
						if(species[sp] == 'human') {
							colnames(species.homologene.list) = c('human','homologene','mouse')
							human.homologene.list = species.homologene.list
							write.table(species.homologene.list,paste('./Data/',Sample,'/05_protein2drug/subset/n_',r,'/',folder[f],'/',Sample,'_exo_homologene','_','human2mouse','.sif',sep=""),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
						} else if(species[sp] == 'mouse') {
							colnames(species.homologene.list) = c('mouse','homologene','human')
							mouse.homologene.list = species.homologene.list
							write.table(species.homologene.list,paste('./Data/',Sample,'/05_protein2drug/subset/n_',r,'/',folder[f],'/',Sample,'_exo_homologene','_','mouse2human','.sif',sep=""),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
						}
						
						
						# --- species UniProt to species DrugBank (protein2drug) ---
						species.drug.list <- c()
						for(j in 1:dim(input)[1])
						{
							#print(j)
							if(species.UniProt[j]!="")
							{
								if((!is.na(species.DrugBank[j])) & (species.DrugBank[j]!=""))
								{
									for(k in 1:length(strsplit(species.DrugBank[j],";")[[1]]))
									{
										species.drug <- cbind(species.UniProt[j],paste(species[sp],'_','protein2drug',sep=""),strsplit(species.DrugBank[j],";")[[1]][k])
										#print(species.drug)
										
										species.drug.list <- rbind(species.drug.list,species.drug)
									}
								}
							}
						}
						#print(species.drug.list)
						species.drug.list = unique(species.drug.list)
						if(species[sp] == 'human') {
							human.drug.list = species.drug.list
						} else if(species[sp] == 'mouse') {
							mouse.drug.list = species.drug.list
						}
						write.table(species.drug.list,paste('./Data/',Sample,'/05_protein2drug/subset/n_',r,'/',folder[f],'/',Sample,'_exo_',species[sp],'_','protein2drug','.sif',sep=""),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
						
						
						# --- species UniProt to species Interaction UniProt (protein2protein) ---
						species.interaction.list <- c()
						for(j in 1:dim(input)[1])
						{
							#print(j)
							if(species.UniProt[j]!="")
							{
								if((!is.na(species.UniProt.interaction[j])) & (species.UniProt.interaction[j]!=""))
								{
									for(k in 1:length(strsplit(species.UniProt.interaction[j],"; ")[[1]]))
									{
										if(strsplit(species.UniProt.interaction[j],"; ")[[1]][k] == "Itself")
										{
											species.interaction <- cbind(species.UniProt[j],paste(species[sp],'_','protein2protein',sep=""),species.UniProt[j])
											#print(species.interaction)
										}
										else
										{
											species.interaction <- cbind(species.UniProt[j],paste(species[sp],'_','protein2protein',sep=""),strsplit(species.UniProt.interaction[j],"; ")[[1]][k])
											#print(species.interaction)
										}
										
										species.interaction.list <- rbind(species.interaction.list,species.interaction)
									}
								}
							}
						}
						#print(species.interaction.list)
						species.interaction.list = unique(species.interaction.list)
						if(species[sp] == 'human') {
							human.interaction.list = species.interaction.list
						} else if(species[sp] == 'mouse') {
							mouse.interaction.list = species.interaction.list
						}
						write.table(species.interaction.list,paste('./Data/',Sample,'/05_protein2drug/subset/n_',r,'/',folder[f],'/',Sample,'_exo_',species[sp],'_','protein2protein','.sif',sep=""),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
						
						
						# --- species UniProt to species Interaction UniProt drug (protein2protein2drug) ---
						species.interaction.drug.list <- c()
						for(j in 1:dim(input)[1])
						{
							#print(j)
							if(species.UniProt[j]!="")
							{
								if((!is.na(species.UniProt.interaction.drug[j])) & (species.UniProt.interaction.drug[j]!=""))
								{
									for(m in 1:length(strsplit(strsplit(species.UniProt.interaction.drug[j],"; ")[[1]],"\\(")))
									{
										species.interaction.target.id <- strsplit(strsplit(species.UniProt.interaction.drug[j],"; ")[[1]],"\\(")[[m]][1]
										
										if(species.interaction.target.id == "Itself") {
											species.interaction.target.id <- species.UniProt[j]
											#print(species.interaction.target.id)
										} else {
											species.interaction.target.id <- strsplit(strsplit(species.UniProt.interaction.drug[j],"; ")[[1]],"\\(")[[m]][1]
											#print(species.interaction.target.id)
										}
										
										for(n in 1:length(strsplit(gsub("[[:punct:]]",replacement=" ",strsplit(strsplit(species.UniProt.interaction.drug[j],"; ")[[1]],"\\(")[[m]][2])," ")[[1]]))
										{
											species.interaction.drug.id <- strsplit(gsub("[[:punct:]]",replacement=" ",strsplit(strsplit(species.UniProt.interaction.drug[j],"; ")[[1]],"\\(")[[m]][2])," ")[[1]][n]
											
											species.interaction.drug <- cbind(species.interaction.target.id,paste(species[sp],'_','protein2protein2drug',sep=""),species.interaction.drug.id)
											#print(species.interaction.drug)
											
											species.interaction.drug.list <- rbind(species.interaction.drug.list,species.interaction.drug)
										}
									}
								}
							}
						}
						#print(species.interaction.drug.list)
						species.interaction.drug.list = unique(species.interaction.drug.list)
						if(species[sp] == 'human') {
							human.interaction.drug.list = species.interaction.drug.list
						} else if(species[sp] == 'mouse') {
							mouse.interaction.drug.list = species.interaction.drug.list
						}
						write.table(species.interaction.drug.list,paste('./Data/',Sample,'/05_protein2drug/subset/n_',r,'/',folder[f],'/',Sample,'_exo_',species[sp],'_','protein2protein2drug','.sif',sep=""),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
						
						
						# --- species interaction with interaction drug (protein2protein_drugsubset) ---
						species.interaction.subset.idx <- which(species.interaction.list[,3] %in% species.interaction.drug.list[!duplicated(species.interaction.drug.list[,1]),1])
						species.interaction.subset <- species.interaction.list[species.interaction.subset.idx,,drop=FALSE]
						species.interaction.subset = unique(species.interaction.subset)
						if(species[sp] == 'human') {
							human.interaction.subset = species.interaction.subset
						} else if(species[sp] == 'mouse') {
							mouse.interaction.subset = species.interaction.subset
						}
						write.table(species.interaction.subset,paste('./Data/',Sample,'/05_protein2drug/subset/n_',r,'/',folder[f],'/',Sample,'_exo_',species[sp],'_','protein2protein_drugsubset','.sif',sep=""),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
						
						
						# --- matched homologene (homologene subset) ---
						species.homologene.drug.subset <- species.homologene.list[sort(unique(which(species.homologene.list[,species[sp]] %in% species.drug.list[!duplicated(species.drug.list[,1]),1]))),,drop=FALSE]
						species.homologene.interaction.subset <- species.homologene.list[sort(unique(which(species.homologene.list[,species[sp]] %in% species.interaction.subset[!duplicated(species.interaction.subset[,1]),1]))),,drop=FALSE]
						species.homologene.subset <- species.homologene.list[sort(unique(c(which(species.homologene.list[,species[sp]] %in% species.drug.list[!duplicated(species.drug.list[,1]),1]),which(species.homologene.list[,species[sp]] %in% species.interaction.subset[!duplicated(species.interaction.subset[,1]),1])))),,drop=FALSE]
						if(species[sp] == 'human') {
							human.homologene.drug.subset = species.homologene.drug.subset
							human.homologene.interaction.subset = species.homologene.interaction.subset
							human.homologene.subset = species.homologene.subset
						} else if(species[sp] == 'mouse') {
							mouse.homologene.drug.subset = species.homologene.drug.subset
							mouse.homologene.interaction.subset = species.homologene.interaction.subset
							mouse.homologene.subset = species.homologene.subset
						}
						
						
						# --- merge (species with drug, interaction drug, and both) ---
						if(species[sp] == 'human') {
							write.table(unique(species.drug.list),paste('./Data/',Sample,'/05_protein2drug/subset/n_',r,'/',folder[f],'/',Sample,'_exo_',species[sp],'_directDrug.sif',sep=""),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
							write.table(unique(rbind(species.interaction.subset,species.interaction.drug.list)),paste('./Data/',Sample,'/05_protein2drug/subset/n_',r,'/',folder[f],'/',Sample,'_exo_',species[sp],'_indirectDrug.sif',sep=""),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
							write.table(unique(rbind(species.drug.list,species.interaction.subset,species.interaction.drug.list)),paste('./Data/',Sample,'/05_protein2drug/subset/n_',r,'/',folder[f],'/',Sample,'_exo_',species[sp],'_Drug.sif',sep=""),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
						} else if(species[sp] == 'mouse') {
							write.table(unique(rbind(species.homologene.drug.subset,species.drug.list)),paste('./Data/',Sample,'/05_protein2drug/subset/n_',r,'/',folder[f],'/',Sample,'_exo_',species[sp],'_directDrug.sif',sep=""),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
							write.table(unique(rbind(species.homologene.interaction.subset,species.interaction.subset,species.interaction.drug.list)),paste('./Data/',Sample,'/05_protein2drug/subset/n_',r,'/',folder[f],'/',Sample,'_exo_',species[sp],'_indirectDrug.sif',sep=""),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
							write.table(unique(rbind(species.homologene.subset,species.drug.list,species.interaction.subset,species.interaction.drug.list)),paste('./Data/',Sample,'/05_protein2drug/subset/n_',r,'/',folder[f],'/',Sample,'_exo_',species[sp],'_Drug.sif',sep=""),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
						}
					}
					
					
					# === human and mouse =========================================================
					# --- Merge human and mouse (with drug) ---
					human2mouse.homologene.drug.subset = unique(rbind(human.homologene.drug.subset[,c('human','homologene','mouse')],mouse.homologene.drug.subset[,c('human','homologene','mouse')]))
					write.table(unique(rbind(human.drug.list,human2mouse.homologene.drug.subset,mouse.drug.list)),paste('./Data/',Sample,'/05_protein2drug/subset/n_',r,'/',folder[f],'/',Sample,'_exo_human2mouse_directDrug.sif',sep=""),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
					
					# --- Merge human and mouse (with interaction drug) ---
					human2mouse.homologene.interaction.subset = unique(rbind(human.homologene.interaction.subset[,c('human','homologene','mouse')],mouse.homologene.interaction.subset[,c('human','homologene','mouse')]))
					write.table(unique(rbind(human.interaction.subset,human.interaction.drug.list,human2mouse.homologene.interaction.subset,mouse.interaction.subset,mouse.interaction.drug.list)),paste('./Data/',Sample,'/05_protein2drug/subset/n_',r,'/',folder[f],'/',Sample,'_exo_human2mouse_indirectDrug.sif',sep=""),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
					
					# --- Merge human and mouse (with drug and interaction drug) ---
					human2mouse.homologene.subset = unique(rbind(human.homologene.subset[,c('human','homologene','mouse')],mouse.homologene.subset[,c('human','homologene','mouse')]))
					write.table(unique(rbind(human.drug.list,human.interaction.subset,human.interaction.drug.list,human2mouse.homologene.subset,mouse.drug.list,mouse.interaction.subset,mouse.interaction.drug.list)),paste('./Data/',Sample,'/05_protein2drug/subset/n_',r,'/',folder[f],'/',Sample,'_exo_human2mouse_Drug.sif',sep=""),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
				}
			}
		}
	}
}

