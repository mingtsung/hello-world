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
subDir <- "Database"
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


# === id mapping ================================================================================================================
#biomaRt: Interface to BioMart databases (e.g. Ensembl, COSMIC, Wormbase and Gramene) (https://www.bioconductor.org/packages/release/bioc/html/biomaRt.html)
#source("https://bioconductor.org/biocLite.R")		#source("http://bioconductor.org/biocLite.R")
#biocLite("biomaRt")								# for id conversion
library(biomaRt)
#ensembl = useMart(biomart="ensembl");	listDatasets(ensembl);
ensembl.hsa = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
ensembl.mmu = useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")
ensembl.rno = useMart(biomart="ensembl", dataset="rnorvegicus_gene_ensembl")

#BridgeDb | standardized access to gene, protein and metabolite identifier mapping services. (http://www.bridgedb.org/)
#BridgeDbR: Code for using BridgeDb identifier mapping framework from within R (https://www.bioconductor.org/packages/release/bioc/html/BridgeDbR.html)
#source("https://bioconductor.org/biocLite.R")		#source("http://bioconductor.org/biocLite.R")
#biocLite("BridgeDbR")								# for id conversion
library(BridgeDbR)
#code.hsa = getOrganismCode("Homo sapiens")
#code.mmu = getOrganismCode("Mus musculus")
#code.rno = getOrganismCode("Rattus norvegicus")
#getBridgeNames()	# Hs_Derby_Ensembl_89.bridge ; Mm_Derby_Ensembl_89.bridge ; Rn_Derby_Ensembl_89.bridge
#dbLocation.hsa <- getDatabase("Homo sapiens",location=getwd())			# http://bridgedb.org/data/gene_database/Hs_Derby_Ensembl_89.bridge
#dbLocation.mmu <- getDatabase("Mus musculus",location=getwd())			# http://bridgedb.org/data/gene_database/Mm_Derby_Ensembl_89.bridge
#dbLocation.rno <- getDatabase("Rattus norvegicus",location=getwd())		# http://bridgedb.org/data/gene_database/Rn_Derby_Ensembl_89.bridge
location.hsa <- './Database/BridgeDb/gene_database/Hs_Derby_Ensembl_89.bridge'
location.mmu <- './Database/BridgeDb/gene_database/Mm_Derby_Ensembl_89.bridge'
location.rno <- './Database/BridgeDb/gene_database/Rn_Derby_Ensembl_89.bridge'
mapper.hsa <- loadDatabase(location.hsa)
mapper.mmu <- loadDatabase(location.mmu)
mapper.rno <- loadDatabase(location.rno)
#Usage: map(mapper, source, identifier, target)
#	mapper = loaded BridgeDb identifier mapper
#	source = system code of the data source		# http://www.bridgedb.org/documentation/system-codes/
#	identifier = identifier to be converted
#	target = system code of the target data source (optional)
#map(mapper.hsa, "L", "196410", "X")		#Entrez Gene identifier (system code: L); Affy identifiers (system code: X)


# === search for mouse id =======================================================================================================
mouse.idmapping <- read.table("./Database/UniProt/idmapping/by_organism/MOUSE_10090_idmapping_selected.tab", sep="\t", quote="", header=FALSE, stringsAsFactors=FALSE)
mouse.idmapping.UniProt <- mouse.idmapping[[1]]
mouse.idmapping.EntryName <- mouse.idmapping[[2]]
mouse.idmapping.GeneID <- mouse.idmapping[[3]]
mouse.idmapping.RefSeq <- mouse.idmapping[[4]]
mouse.idmapping.GI <- mouse.idmapping[[5]]

mouse.UniProt.reviewed <- read.csv("./Database/UniProt/Mouse/uniprot_Mus_musculus_10090_reviewed.csv", header=TRUE, stringsAsFactors=FALSE)
mouse.UniProt.reviewed.Entry <- mouse.UniProt.reviewed$Entry
mouse.UniProt.reviewed.EntryName <- mouse.UniProt.reviewed$Entry.name
mouse.UniProt.reviewed.ProteinNames <- mouse.UniProt.reviewed$Protein.names
mouse.UniProt.reviewed.GeneID <- sub(";$",replacement="",mouse.UniProt.reviewed$Cross.reference..GeneID.)
mouse.UniProt.reviewed.GeneID[is.na(mouse.UniProt.reviewed.GeneID)] <- ''
mouse.UniProt.reviewed.DrugBank <- sub(";$",replacement="",mouse.UniProt.reviewed$Cross.reference..DrugBank.)
mouse.UniProt.reviewed.DrugBank[is.na(mouse.UniProt.reviewed.DrugBank)] <- ''
mouse.UniProt.reviewed.Interaction <- mouse.UniProt.reviewed$Interacts.with

mouse.UniProt.unreviewed <- read.csv("./Database/UniProt/Mouse/uniprot_Mus_musculus_10090_unreviewed.csv", header=TRUE, stringsAsFactors=FALSE)
mouse.UniProt.unreviewed.Entry <- mouse.UniProt.unreviewed$Entry
mouse.UniProt.unreviewed.EntryName <- mouse.UniProt.unreviewed$Entry.name
mouse.UniProt.unreviewed.ProteinNames <- mouse.UniProt.unreviewed$Protein.names
mouse.UniProt.unreviewed.GeneID <- sub(";$",replacement="",mouse.UniProt.unreviewed$Cross.reference..GeneID.)
mouse.UniProt.unreviewed.GeneID[is.na(mouse.UniProt.unreviewed.GeneID)] <- ''
mouse.UniProt.unreviewed.DrugBank <- sub(";$",replacement="",mouse.UniProt.unreviewed$Cross.reference..DrugBank.)
mouse.UniProt.unreviewed.DrugBank[is.na(mouse.UniProt.unreviewed.DrugBank)] <- ''
mouse.UniProt.unreviewed.Interaction <- mouse.UniProt.unreviewed$Interacts.with


# === search for human id =======================================================================================================
human.idmapping <- read.table("./Database/UniProt/idmapping/by_organism/HUMAN_9606_idmapping_selected.tab", sep="\t", quote="", header=FALSE, stringsAsFactors=FALSE)
human.idmapping.UniProt <- human.idmapping[[1]]
human.idmapping.EntryName <- human.idmapping[[2]]
human.idmapping.GeneID <- human.idmapping[[3]]
human.idmapping.RefSeq <- human.idmapping[[4]]
human.idmapping.GI <- human.idmapping[[5]]

human.UniProt.reviewed <- read.csv("./Database/UniProt/Human/uniprot_Homo_sapiens_9606_reviewed.csv", header=TRUE, stringsAsFactors=FALSE)
human.UniProt.reviewed.Entry <- human.UniProt.reviewed$Entry
human.UniProt.reviewed.EntryName <- human.UniProt.reviewed$Entry.name
human.UniProt.reviewed.ProteinNames <- human.UniProt.reviewed$Protein.names
human.UniProt.reviewed.GeneID <- sub(";$",replacement="",human.UniProt.reviewed$Cross.reference..GeneID.)
human.UniProt.reviewed.GeneID[is.na(human.UniProt.reviewed.GeneID)] <- ''
human.UniProt.reviewed.DrugBank <- sub(";$",replacement="",human.UniProt.reviewed$Cross.reference..DrugBank.)
human.UniProt.reviewed.DrugBank[is.na(human.UniProt.reviewed.DrugBank)] <- ''
human.UniProt.reviewed.Interaction <- human.UniProt.reviewed$Interacts.with

human.UniProt.unreviewed <- read.csv("./Database/UniProt/Human/uniprot_Homo_sapiens_9606_unreviewed.csv", header=TRUE, stringsAsFactors=FALSE)
human.UniProt.unreviewed.Entry <- human.UniProt.unreviewed$Entry
human.UniProt.unreviewed.EntryName <- human.UniProt.unreviewed$Entry.name
human.UniProt.unreviewed.ProteinNames <- human.UniProt.unreviewed$Protein.names
human.UniProt.unreviewed.GeneID <- sub(";$",replacement="",human.UniProt.unreviewed$Cross.reference..GeneID.)
human.UniProt.unreviewed.GeneID[is.na(human.UniProt.unreviewed.GeneID)] <- ''
human.UniProt.unreviewed.DrugBank <- sub(";$",replacement="",human.UniProt.unreviewed$Cross.reference..DrugBank.)
human.UniProt.unreviewed.DrugBank[is.na(human.UniProt.unreviewed.DrugBank)] <- ''
human.UniProt.unreviewed.Interaction <- human.UniProt.unreviewed$Interacts.with

human.idmapping.GeneID[which(human.idmapping.UniProt == 'P62861')] = c(2197)					# http://www.uniprot.org/uniprot/P62861; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:3597; https://www.ncbi.nlm.nih.gov/gene/2197
human.UniProt.reviewed.GeneID[which(human.UniProt.reviewed.Entry == 'P62861')] = c(2197)		# http://www.uniprot.org/uniprot/P62861; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:3597; https://www.ncbi.nlm.nih.gov/gene/2197
human.UniProt.unreviewed.GeneID[which(human.UniProt.unreviewed.Entry == 'P62861')] = c(2197)	# http://www.uniprot.org/uniprot/P62861; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:3597; https://www.ncbi.nlm.nih.gov/gene/2197
human.idmapping.GeneID[which(human.idmapping.UniProt == 'O14950')] = c(103910)					# http://www.uniprot.org/uniprot/O14950; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:29827; https://www.ncbi.nlm.nih.gov/gene/103910 (https://www.ncbi.nlm.nih.gov/gene/10627)
human.UniProt.reviewed.GeneID[which(human.UniProt.reviewed.Entry == 'O14950')] = c(103910)		# http://www.uniprot.org/uniprot/O14950; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:29827; https://www.ncbi.nlm.nih.gov/gene/103910 (https://www.ncbi.nlm.nih.gov/gene/10627)
human.UniProt.unreviewed.GeneID[which(human.UniProt.unreviewed.Entry == 'O14950')] = c(103910)	# http://www.uniprot.org/uniprot/O14950; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:29827; https://www.ncbi.nlm.nih.gov/gene/103910 (https://www.ncbi.nlm.nih.gov/gene/10627)


# === input Top100 recorded in Microvesicle database ============================================================================
#Database <- c('ExoCarta','Vesiclepedia','EVpedia')

# --- ExoCarta --------------------------------------------------------------------------------------------------------
ExoCarta_Top100 <- read.table(paste('./Database/Microvesicle/ExoCarta/ExoCarta_top100_protein_details_5.txt',sep=""), sep="\t", quote="", header=TRUE, stringsAsFactors=FALSE)

ExoCarta_Top100.Gene_Symbol <- trimws(ExoCarta_Top100[,'Gene.Symbol'])
ExoCarta_Top100.Number_of_times_identified <- trimws(ExoCarta_Top100[,'Number.of.times.identified'])

# /// human format ///
mainDir <- paste('./Data/',Sample,'/10_datasetOverlap/Database',sep="")
subDir <- "hsa"
if (file.exists(file.path(mainDir, subDir))){
} else {
	dir.create(file.path(mainDir, subDir))
}

ExoCarta_Top100.Gene_Symbol.hsa_format = ExoCarta_Top100.Gene_Symbol

ensembl = ensembl.hsa
#listFilters(ensembl)[,c('name','description')]
#listAttributes(ensembl)[,c('name','description')]

# --- search for IDs by hgnc_symbol ---
a = ExoCarta_Top100.Gene_Symbol
b = ensembl.hsa
ExoCarta_Top100.Gene_Symbol.hgnc_symbol.attributes <- getBM(attributes=c('entrezgene','hgnc_id','uniprotswissprot'), filters='hgnc_symbol', values=a, mart=b)
ExoCarta_Top100.Gene_Symbol.hgnc_symbol.wikigene_description <- getBM(attributes=c('hgnc_id','wikigene_description'), filters='hgnc_symbol', values=a, mart=b)
ExoCarta_Top100.Gene_Symbol.hgnc_symbol.wikigene_description.table = merge(ExoCarta_Top100.Gene_Symbol.hgnc_symbol.wikigene_description, ExoCarta_Top100.Gene_Symbol.hgnc_symbol.attributes, by.x="hgnc_id", by.y="hgnc_id")
ExoCarta_Top100.Gene_Symbol.hgnc_symbol.attributes = ExoCarta_Top100.Gene_Symbol.hgnc_symbol.wikigene_description.table[,c('entrezgene','hgnc_id','uniprotswissprot','wikigene_description'),drop=FALSE]
colnames(ExoCarta_Top100.Gene_Symbol.hgnc_symbol.attributes)[which(colnames(ExoCarta_Top100.Gene_Symbol.hgnc_symbol.attributes)=='wikigene_description')] = 'description'
ExoCarta_Top100.Gene_Symbol.hgnc_symbol.attributes = unique(ExoCarta_Top100.Gene_Symbol.hgnc_symbol.attributes)
ExoCarta_Top100.Gene_Symbol.hgnc_id <- getBM(attributes=c('hgnc_symbol','hgnc_id'), filters='hgnc_symbol', values=a, mart=b)
ExoCarta_Top100.Gene_Symbol.hgnc_id.table = merge(ExoCarta_Top100.Gene_Symbol.hgnc_id, ExoCarta_Top100.Gene_Symbol.hgnc_symbol.attributes, by.x="hgnc_id", by.y="hgnc_id")
ExoCarta_Top100.Gene_Symbol.hgnc_id.table = ExoCarta_Top100.Gene_Symbol.hgnc_id.table[,c('hgnc_symbol','hgnc_id','entrezgene','uniprotswissprot','description')]
ExoCarta_Top100.Gene_Symbol.hgnc_id.table = ExoCarta_Top100.Gene_Symbol.hgnc_id.table[order(match(ExoCarta_Top100.Gene_Symbol.hgnc_id.table[,'hgnc_symbol'], a)),]

# merge multiple uniprotswissprot of single gene
unique.hgnc_id.idx <- which(!duplicated(ExoCarta_Top100.Gene_Symbol.hgnc_id.table[,'hgnc_id']))
unique.hgnc_id <- ExoCarta_Top100.Gene_Symbol.hgnc_id.table[unique.hgnc_id.idx,'hgnc_id']
duplicated.hgnc_id.idx <- which(duplicated(ExoCarta_Top100.Gene_Symbol.hgnc_id.table[,'hgnc_id']))
duplicated.hgnc_id <- ExoCarta_Top100.Gene_Symbol.hgnc_id.table[duplicated.hgnc_id.idx,'hgnc_id']
duplicated.hgnc_id.idx.pair <- unique(sapply(duplicated.hgnc_id, FUN=function(X) which(ExoCarta_Top100.Gene_Symbol.hgnc_id.table[,'hgnc_id'] %in% X)))
duplicated.hgnc_id.idx.uniprotswissprot <- sapply(duplicated.hgnc_id.idx.pair, FUN=function(X) paste(sort(unique(ExoCarta_Top100.Gene_Symbol.hgnc_id.table[X,'uniprotswissprot'][ExoCarta_Top100.Gene_Symbol.hgnc_id.table[X,'uniprotswissprot']!=''])),collapse=';'))
if(length(duplicated.hgnc_id.idx.pair) > 0)
{
	for(i in 1:length(duplicated.hgnc_id.idx.pair))
	{
		ExoCarta_Top100.Gene_Symbol.hgnc_id.table[unlist(duplicated.hgnc_id.idx.pair[i]),'uniprotswissprot'] = duplicated.hgnc_id.idx.uniprotswissprot[i]
	}
}

# uniprotswissprot not found by R package "biomaRt"
empty.uniprotswissprot.idx <- which(ExoCarta_Top100.Gene_Symbol.hgnc_id.table[,'uniprotswissprot']=='')
empty.uniprotswissprot <- ExoCarta_Top100.Gene_Symbol.hgnc_id.table[empty.uniprotswissprot.idx,]

ExoCarta_Top100.Gene_Symbol.hgnc_id.table = unique(ExoCarta_Top100.Gene_Symbol.hgnc_id.table)
#dim(ExoCarta_Top100.Gene_Symbol.hgnc_id.table)

# remove MIR3620 (microRNA 3620) (https://www.ncbi.nlm.nih.gov/gene/100500810)
ExoCarta_Top100.Gene_Symbol.hgnc_id.table = ExoCarta_Top100.Gene_Symbol.hgnc_id.table[which(ExoCarta_Top100.Gene_Symbol.hgnc_id.table[,'entrezgene']!='100500810'),]

rownames(ExoCarta_Top100.Gene_Symbol.hgnc_id.table) = c(1:nrow(ExoCarta_Top100.Gene_Symbol.hgnc_id.table))
#ExoCarta_Top100.Gene_Symbol.hgnc_id.table

ExoCarta_Top100.Gene_Symbol.hgnc_id.table.New = ExoCarta_Top100.Gene_Symbol.hgnc_id.table
#ExoCarta_Top100.Gene_Symbol.hgnc_id.table.New

# Gene_Symbol: attributes not found by R package "biomaRt"
ExoCarta_Top100.Gene_Symbol.hgnc_symbol.NA_attributes <- a[which(! a %in% ExoCarta_Top100.Gene_Symbol.hgnc_id.table[,'hgnc_symbol'])]
if(length(ExoCarta_Top100.Gene_Symbol.hgnc_symbol.NA_attributes) > 0)
{
	manually_curated.hgnc_id.table <- matrix('',length(ExoCarta_Top100.Gene_Symbol.hgnc_symbol.NA_attributes),ncol(ExoCarta_Top100.Gene_Symbol.hgnc_id.table))
	colnames(manually_curated.hgnc_id.table) = colnames(ExoCarta_Top100.Gene_Symbol.hgnc_id.table)
	
	ExoCarta_Top100.Gene_Symbol.hgnc_symbol.NA_attributes.uncurated <- c()
	
	for(i in 1:length(ExoCarta_Top100.Gene_Symbol.hgnc_symbol.NA_attributes))
	{
		# --- Gene type: protein coding -----------------------------------------------------------
		if(ExoCarta_Top100.Gene_Symbol.hgnc_symbol.NA_attributes[i] == 'RAB1A') {	# https://www.ncbi.nlm.nih.gov/gene/5861; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:9758; http://www.uniprot.org/uniprot/P62820
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RAB1A'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:9758'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '5861'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'P62820'
			manually_curated.hgnc_id.table[i,'description'] = 'RAB1A, member RAS oncogene family [Source:HGNC Symbol;Acc:HGNC:9758]'
		} else {
			cat(paste('"',ExoCarta_Top100.Gene_Symbol.hgnc_symbol.NA_attributes[i],'" did not match queried attributes (\'entrezgene\',\'hgnc_id\',\'uniprotswissprot\',\'description\') by using R package "biomaRt".',sep=""),"\n")
			
			ExoCarta_Top100.Gene_Symbol.hgnc_symbol.NA_attributes.uncurated <- c(ExoCarta_Top100.Gene_Symbol.hgnc_symbol.NA_attributes.uncurated, ExoCarta_Top100.Gene_Symbol.hgnc_symbol.NA_attributes[i])
		}
		
		manually_curated.hgnc_id.table[i,'description'] = trimws(manually_curated.hgnc_id.table[i,'description'])
		manually_curated.hgnc_id.table[i,'description'] = sapply(strsplit(manually_curated.hgnc_id.table[i,'description'],'\\|'), FUN=function(U) paste(unique(sapply(strsplit(unlist(U),'\\[Source:'), FUN=function(V) trimws(unlist(V))[1])),collapse="|"))
	}
	
	if(length(ExoCarta_Top100.Gene_Symbol.hgnc_symbol.NA_attributes.uncurated) > 0)
	{
		write.table(ExoCarta_Top100.Gene_Symbol.hgnc_symbol.NA_attributes.uncurated, paste('./Data/',Sample,'/10_datasetOverlap/Database/mmu/ExoCarta_All.Protein_Gene2HGNC_NA.txt',sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE)
	}
	
	ExoCarta_Top100.Gene_Symbol.hgnc_id.table.New = rbind(ExoCarta_Top100.Gene_Symbol.hgnc_id.table, manually_curated.hgnc_id.table)
}

rownames(ExoCarta_Top100.Gene_Symbol.hgnc_id.table.New) = c(1:nrow(ExoCarta_Top100.Gene_Symbol.hgnc_id.table.New))
#ExoCarta_Top100.Gene_Symbol.hgnc_id.table.New

ExoCarta_Top100.Gene_Symbol.hgnc_id.table = ExoCarta_Top100.Gene_Symbol.hgnc_id.table.New
#ExoCarta_Top100.Gene_Symbol.hgnc_id.table

ExoCarta_Top100.Gene_Symbol.hgnc_id.table = ExoCarta_Top100.Gene_Symbol.hgnc_id.table[order(match(ExoCarta_Top100.Gene_Symbol.hgnc_id.table[,'hgnc_symbol'], a)),]
rownames(ExoCarta_Top100.Gene_Symbol.hgnc_id.table) = c(1:nrow(ExoCarta_Top100.Gene_Symbol.hgnc_id.table))
#ExoCarta_Top100.Gene_Symbol.hgnc_id.table

human.UniProt.Entry2ProteinNames = sapply(strsplit(ExoCarta_Top100.Gene_Symbol.hgnc_id.table[,'uniprotswissprot'],";"), FUN=function(X) paste(unique(sapply(unlist(X), FUN=function(Y) ifelse(length(which(human.UniProt.reviewed.Entry %in% unlist(Y)))==0, ifelse(length(which(human.UniProt.unreviewed.Entry %in% unlist(Y)))==0, '', human.UniProt.unreviewed.ProteinNames[which(human.UniProt.unreviewed.Entry %in% unlist(Y))]), human.UniProt.reviewed.ProteinNames[which(human.UniProt.reviewed.Entry %in% unlist(Y))]))),collapse="|"))
ExoCarta_Top100.Gene_Symbol.hgnc_id.table = cbind(ExoCarta_Top100.Gene_Symbol.hgnc_id.table, human.UniProt.Entry2ProteinNames)
colnames(ExoCarta_Top100.Gene_Symbol.hgnc_id.table)[ncol(ExoCarta_Top100.Gene_Symbol.hgnc_id.table)] = 'Protein_names'

#write.csv(ExoCarta_Top100.Gene_Symbol.hgnc_id.table, paste('./Data/',Sample,'/10_datasetOverlap/Database/hsa/ExoCarta_Top100.Protein_GeneSymbol2HGNC.table.csv',sep=""), quote=TRUE, row.names=TRUE)

# add in "Number of times identified" of "Gene Symbol"
ExoCarta_Top100.hsa_format = ExoCarta_Top100
ExoCarta_Top100.hsa_format[,'Gene.Symbol'] = ExoCarta_Top100.Gene_Symbol.hsa_format
ExoCarta_Top100.Gene_Symbol.hgnc_id.table.Number_of_times_identified = merge(ExoCarta_Top100.hsa_format, ExoCarta_Top100.Gene_Symbol.hgnc_id.table, by.x="Gene.Symbol", by.y="hgnc_symbol")
colnames(ExoCarta_Top100.Gene_Symbol.hgnc_id.table.Number_of_times_identified) = c('hgnc_symbol','Number_of_times_identified',colnames(ExoCarta_Top100.Gene_Symbol.hgnc_id.table)[-1])
ExoCarta_Top100.Gene_Symbol.hgnc_id.table.Number_of_times_identified = ExoCarta_Top100.Gene_Symbol.hgnc_id.table.Number_of_times_identified[order(match(ExoCarta_Top100.Gene_Symbol.hgnc_id.table.Number_of_times_identified[,'hgnc_symbol'], a)),]
rownames(ExoCarta_Top100.Gene_Symbol.hgnc_id.table.Number_of_times_identified) = order(match(a, ExoCarta_Top100.Gene_Symbol.hgnc_id.table.Number_of_times_identified[,'hgnc_symbol']))
ExoCarta_Top100.Gene_Symbol.hgnc_id.table.Number_of_times_identified = ExoCarta_Top100.Gene_Symbol.hgnc_id.table.Number_of_times_identified[order(as.numeric(rownames(ExoCarta_Top100.Gene_Symbol.hgnc_id.table.Number_of_times_identified))),]
#ExoCarta_Top100.Gene_Symbol.hgnc_id.table.Number_of_times_identified
#write.csv(ExoCarta_Top100.Gene_Symbol.hgnc_id.table.Number_of_times_identified, paste('./Data/',Sample,'/10_datasetOverlap/Database/hsa/ExoCarta_Top100.Protein_GeneSymbol2HGNC.table_IdentificationNumber.csv',sep=""), quote=TRUE, row.names=TRUE)

# --- write xlsx file by openxlsx package ---------------------------------------------------------
# Importing a big xlsx file into R? https://stackoverflow.com/questions/19147884/importing-a-big-xlsx-file-into-r/43118530#43118530
# Building R for Windows: https://cran.r-project.org/bin/windows/Rtools/
#Sys.setenv("R_ZIPCMD" = "C:/Rtools/bin/zip.exe")	# path to zip.exe

wb <- openxlsx::createWorkbook(paste('./Data/',Sample,'/10_datasetOverlap/Database/hsa/ExoCarta_Top100.Protein_GeneSymbol2HGNC.table.xlsx',sep=""))
modifyBaseFont(wb, fontSize=12, fontColour="black", fontName="Times New Roman")

for(s in 1:2)
{
	if(s == 1) {
		sheet = 'ExoCarta_Top100'
		sheetData = ExoCarta_Top100.Gene_Symbol.hgnc_id.table
	} else if(s == 2) {
		sheet = 'ExoCarta_Top100 (with number)'
		sheetData = ExoCarta_Top100.Gene_Symbol.hgnc_id.table.Number_of_times_identified
	}
	
	if(!is.null(sheetData))
	{
		if(nrow(sheetData) > 0)
		{
			addWorksheet(wb, sheet)
			writeData(wb, sheet, sheetData, colNames=TRUE, rowNames=TRUE, keepNA=FALSE)
			setColWidths(wb, sheet, cols=c(1:(1+ncol(sheetData))), widths=c(8.5,rep(15,ncol(sheetData))))
			setRowHeights(wb, sheet, rows=c(1:(1+nrow(sheetData))), heights=16.5)
			freezePane(wb, sheet, firstActiveRow=2)
			style <- createStyle(halign="left", valign="center")
			addStyle(wb, sheet, style, rows=c(1:(1+nrow(sheetData))), cols=c(1:(1+ncol(sheetData))), gridExpand=TRUE)
			style <- createStyle(textDecoration="bold")
			addStyle(wb, sheet, style, rows=1, cols=c(1:(1+ncol(sheetData))), stack=TRUE)
			
			right_align.idx = which(colnames(sheetData) %in% c('Number_of_times_identified'))
			if(length(right_align.idx) > 0) {
				style <- createStyle(halign="right", valign="center")
				addStyle(wb, sheet, style, rows=c(2:(1+nrow(sheetData))), cols=1+right_align.idx, gridExpand=TRUE, stack=TRUE)
			}
		}
	}
}

if(length(wb$sheet_names) > 0) {
	openxlsx::saveWorkbook(wb, file=paste('./Data/',Sample,'/10_datasetOverlap/Database/hsa/ExoCarta_Top100.Protein_GeneSymbol2HGNC.table.xlsx',sep=""), overwrite=TRUE)
}


# /// mouse format ///
mainDir <- paste('./Data/',Sample,'/10_datasetOverlap/Database',sep="")
subDir <- "mmu"
if (file.exists(file.path(mainDir, subDir))){
} else {
	dir.create(file.path(mainDir, subDir))
}

ExoCarta_Top100.Gene_Symbol.mmu_format <- paste(toupper(substr(ExoCarta_Top100.Gene_Symbol,1,1)), tolower(substr(ExoCarta_Top100.Gene_Symbol,2,nchar(ExoCarta_Top100.Gene_Symbol))), sep="")
ExoCarta_Top100.Gene_Symbol.mmu_format[which(ExoCarta_Top100.Gene_Symbol.mmu_format=='Rab7a')] = 'Rab7'			# Rab7: Also known as Rab7a (https://www.ncbi.nlm.nih.gov/gene/19349)
ExoCarta_Top100.Gene_Symbol.mmu_format.Hist2h4a_idx = which(ExoCarta_Top100.Gene_Symbol.mmu_format=='Hist2h4a')
ExoCarta_Top100.Gene_Symbol.mmu_format[which(ExoCarta_Top100.Gene_Symbol.mmu_format=='Hist2h4a')] = 'Hist1h4b'	# Hist1h4b: Gene Overview MyGene.info: HIST2H4A (http://www.informatics.jax.org/marker/MGI:2448420)

ensembl = ensembl.mmu
#listFilters(ensembl)[,c('name','description')]
#listAttributes(ensembl)[,c('name','description')]

# --- search for IDs by mgi_symbol ---
a = ExoCarta_Top100.Gene_Symbol.mmu_format
b = ensembl.mmu
ExoCarta_Top100.Gene_Symbol.mgi_symbol.attributes <- getBM(attributes=c('entrezgene','mgi_id','uniprotswissprot'), filters='mgi_symbol', values=a, mart=b)
ExoCarta_Top100.Gene_Symbol.mgi_symbol.wikigene_description <- getBM(attributes=c('mgi_id','wikigene_description'), filters='mgi_symbol', values=a, mart=b)
ExoCarta_Top100.Gene_Symbol.mgi_symbol.wikigene_description.table = merge(ExoCarta_Top100.Gene_Symbol.mgi_symbol.wikigene_description, ExoCarta_Top100.Gene_Symbol.mgi_symbol.attributes, by.x="mgi_id", by.y="mgi_id")
ExoCarta_Top100.Gene_Symbol.mgi_symbol.attributes = ExoCarta_Top100.Gene_Symbol.mgi_symbol.wikigene_description.table[,c('entrezgene','mgi_id','uniprotswissprot','wikigene_description'),drop=FALSE]
colnames(ExoCarta_Top100.Gene_Symbol.mgi_symbol.attributes)[which(colnames(ExoCarta_Top100.Gene_Symbol.mgi_symbol.attributes)=='wikigene_description')] = 'description'
ExoCarta_Top100.Gene_Symbol.mgi_symbol.attributes = unique(ExoCarta_Top100.Gene_Symbol.mgi_symbol.attributes)
ExoCarta_Top100.Gene_Symbol.mgi_id <- getBM(attributes=c('mgi_symbol','mgi_id'), filters='mgi_symbol', values=a, mart=b)
ExoCarta_Top100.Gene_Symbol.mgi_id.table = merge(ExoCarta_Top100.Gene_Symbol.mgi_id, ExoCarta_Top100.Gene_Symbol.mgi_symbol.attributes, by.x="mgi_id", by.y="mgi_id")

# https://www.ncbi.nlm.nih.gov/gene/11674; http://www.informatics.jax.org/marker/MGI:87994; http://www.uniprot.org/uniprot/P05064
ExoCarta_Top100.Gene_Symbol.mgi_id.table[which(ExoCarta_Top100.Gene_Symbol.mgi_id.table[,'mgi_symbol']=='Aldoa'),'mgi_id'] = 'MGI:87994'
ExoCarta_Top100.Gene_Symbol.mgi_id.table[which(ExoCarta_Top100.Gene_Symbol.mgi_id.table[,'mgi_symbol']=='Aldoa'),'entrezgene'] = '11674'
ExoCarta_Top100.Gene_Symbol.mgi_id.table[which(ExoCarta_Top100.Gene_Symbol.mgi_id.table[,'mgi_symbol']=='Aldoa'),'uniprotswissprot'] = 'P05064'
ExoCarta_Top100.Gene_Symbol.mgi_id.table[which(ExoCarta_Top100.Gene_Symbol.mgi_id.table[,'mgi_symbol']=='Aldoa'),'description'] = 'aldolase A, fructose-bisphosphate [Source:MGI Symbol;Acc:MGI:87994]'
# https://www.ncbi.nlm.nih.gov/gene/14433; http://www.informatics.jax.org/marker/MGI:95640; http://www.uniprot.org/uniprot/P16858
ExoCarta_Top100.Gene_Symbol.mgi_id.table[which(ExoCarta_Top100.Gene_Symbol.mgi_id.table[,'mgi_symbol']=='Gapdh'),'mgi_id'] = 'MGI:95640'
ExoCarta_Top100.Gene_Symbol.mgi_id.table[which(ExoCarta_Top100.Gene_Symbol.mgi_id.table[,'mgi_symbol']=='Gapdh'),'entrezgene'] = '14433'
ExoCarta_Top100.Gene_Symbol.mgi_id.table[which(ExoCarta_Top100.Gene_Symbol.mgi_id.table[,'mgi_symbol']=='Gapdh'),'uniprotswissprot'] = 'P16858'
ExoCarta_Top100.Gene_Symbol.mgi_id.table[which(ExoCarta_Top100.Gene_Symbol.mgi_id.table[,'mgi_symbol']=='Gapdh'),'description'] = 'glyceraldehyde-3-phosphate dehydrogenase [Source:MGI Symbol;Acc:MGI:95640]'

ExoCarta_Top100.Gene_Symbol.mgi_id.table[,'description'] = trimws(ExoCarta_Top100.Gene_Symbol.mgi_id.table[,'description'])
ExoCarta_Top100.Gene_Symbol.mgi_id.table[,'description'] = sapply(strsplit(ExoCarta_Top100.Gene_Symbol.mgi_id.table[,'description'],'\\|'), FUN=function(U) paste(unique(sapply(strsplit(unlist(U),'\\[Source:'), FUN=function(V) trimws(unlist(V))[1])),collapse="|"))
ExoCarta_Top100.Gene_Symbol.mgi_id.table = unique(ExoCarta_Top100.Gene_Symbol.mgi_id.table)

ExoCarta_Top100.Gene_Symbol.mgi_id.table = ExoCarta_Top100.Gene_Symbol.mgi_id.table[,c('mgi_symbol','mgi_id','entrezgene','uniprotswissprot','description')]
ExoCarta_Top100.Gene_Symbol.mgi_id.table = ExoCarta_Top100.Gene_Symbol.mgi_id.table[order(match(ExoCarta_Top100.Gene_Symbol.mgi_id.table[,'mgi_symbol'], a)),]

# merge multiple uniprotswissprot of single gene
unique.mgi_id.idx <- which(!duplicated(ExoCarta_Top100.Gene_Symbol.mgi_id.table[,'mgi_id']))
unique.mgi_id <- ExoCarta_Top100.Gene_Symbol.mgi_id.table[unique.mgi_id.idx,'mgi_id']
duplicated.mgi_id.idx <- which(duplicated(ExoCarta_Top100.Gene_Symbol.mgi_id.table[,'mgi_id']))
duplicated.mgi_id <- ExoCarta_Top100.Gene_Symbol.mgi_id.table[duplicated.mgi_id.idx,'mgi_id']
duplicated.mgi_id.idx.pair <- unique(sapply(duplicated.mgi_id, FUN=function(X) which(ExoCarta_Top100.Gene_Symbol.mgi_id.table[,'mgi_id'] %in% X)))
duplicated.mgi_id.idx.uniprotswissprot <- sapply(duplicated.mgi_id.idx.pair, FUN=function(X) paste(sort(unique(ExoCarta_Top100.Gene_Symbol.mgi_id.table[X,'uniprotswissprot'][ExoCarta_Top100.Gene_Symbol.mgi_id.table[X,'uniprotswissprot']!=''])),collapse=';'))
if(length(duplicated.mgi_id.idx.pair) > 0)
{
	for(i in 1:length(duplicated.mgi_id.idx.pair))
	{
		ExoCarta_Top100.Gene_Symbol.mgi_id.table[unlist(duplicated.mgi_id.idx.pair[i]),'uniprotswissprot'] = duplicated.mgi_id.idx.uniprotswissprot[i]
	}
}

# uniprotswissprot not found by R package "biomaRt"
empty.uniprotswissprot.idx <- which(ExoCarta_Top100.Gene_Symbol.mgi_id.table[,'uniprotswissprot']=='')
empty.uniprotswissprot <- ExoCarta_Top100.Gene_Symbol.mgi_id.table[empty.uniprotswissprot.idx,]
ExoCarta_Top100.Gene_Symbol.mgi_id.table[which(ExoCarta_Top100.Gene_Symbol.mgi_id.table[,'entrezgene']=='21825'),'uniprotswissprot'] = 'P35441'		# https://www.ncbi.nlm.nih.gov/gene/21825	http://www.informatics.jax.org/marker/MGI:98737		http://www.uniprot.org/uniprot/P35441
ExoCarta_Top100.Gene_Symbol.mgi_id.table[which(ExoCarta_Top100.Gene_Symbol.mgi_id.table[,'entrezgene']=='78388'),'uniprotswissprot'] = 'Q9EQK5'		# https://www.ncbi.nlm.nih.gov/gene/78388	http://www.informatics.jax.org/marker/MGI:1925638	http://www.uniprot.org/uniprot/Q9EQK5

ExoCarta_Top100.Gene_Symbol.mgi_id.table = unique(ExoCarta_Top100.Gene_Symbol.mgi_id.table)
#dim(ExoCarta_Top100.Gene_Symbol.mgi_id.table)

# remove Gapdh-ps15 (https://www.ncbi.nlm.nih.gov/gene/100042025)
ExoCarta_Top100.Gene_Symbol.mgi_id.table = ExoCarta_Top100.Gene_Symbol.mgi_id.table[which(ExoCarta_Top100.Gene_Symbol.mgi_id.table[,'entrezgene']!='100042025'),]

rownames(ExoCarta_Top100.Gene_Symbol.mgi_id.table) = c(1:nrow(ExoCarta_Top100.Gene_Symbol.mgi_id.table))
#ExoCarta_Top100.Gene_Symbol.mgi_id.table
#ExoCarta_Top100.Gene_Symbol.mgi_id.table = rbind(ExoCarta_Top100.Gene_Symbol.mgi_id.table[c(c(1:(which(a == 'Hist1h4b')[1]-1)),c((which(a == 'Hist1h4b')[1]+1):(which(a == 'Hist1h4b')[2]-1))),], ExoCarta_Top100.Gene_Symbol.mgi_id.table[which(a == 'Hist1h4b')[1],], ExoCarta_Top100.Gene_Symbol.mgi_id.table[c((which(a == 'Hist1h4b')[2]):(nrow(ExoCarta_Top100.Gene_Symbol.mgi_id.table))),])
#rownames(ExoCarta_Top100.Gene_Symbol.mgi_id.table) = c(c(1:(which(a == 'Hist1h4b')[1]-1)),c((which(a == 'Hist1h4b')[1]+1):(which(a == 'Hist1h4b')[2]-1)),which(a == 'Hist1h4b')[2],c((which(a == 'Hist1h4b')[2]+1):(nrow(ExoCarta_Top100.Gene_Symbol.mgi_id.table)+1)))
#ExoCarta_Top100.Gene_Symbol.mgi_id.table

ExoCarta_Top100.Gene_Symbol.mgi_id.table.New = ExoCarta_Top100.Gene_Symbol.mgi_id.table
#ExoCarta_Top100.Gene_Symbol.mgi_id.table.New

# Gene_Symbol: attributes not found by R package "biomaRt"
ExoCarta_Top100.Gene_Symbol.mgi_symbol.NA_attributes <- a[which(! a %in% ExoCarta_Top100.Gene_Symbol.mgi_id.table[,'mgi_symbol'])]
if(length(ExoCarta_Top100.Gene_Symbol.mgi_symbol.NA_attributes) > 0)
{
	manually_curated.mgi_id.table <- matrix('',length(ExoCarta_Top100.Gene_Symbol.mgi_symbol.NA_attributes),ncol(ExoCarta_Top100.Gene_Symbol.mgi_id.table))
	colnames(manually_curated.mgi_id.table) = colnames(ExoCarta_Top100.Gene_Symbol.mgi_id.table)
	
	ExoCarta_Top100.Gene_Symbol.mgi_symbol.NA_attributes.uncurated <- c()
	
	for(i in 1:length(ExoCarta_Top100.Gene_Symbol.mgi_symbol.NA_attributes))
	{
		# --- Gene type: protein coding -----------------------------------------------------------
		if(ExoCarta_Top100.Gene_Symbol.mgi_symbol.NA_attributes[i] == 'Rab1a') {	# https://www.ncbi.nlm.nih.gov/gene/19324; http://www.informatics.jax.org/marker/MGI:97842; http://www.uniprot.org/uniprot/P62821
			manually_curated.mgi_id.table[i,'mgi_symbol'] = 'Rab1a'
			manually_curated.mgi_id.table[i,'mgi_id'] = 'MGI:97842'
			manually_curated.mgi_id.table[i,'entrezgene'] = '19324'
			manually_curated.mgi_id.table[i,'uniprotswissprot'] = 'P62821'
			manually_curated.mgi_id.table[i,'description'] = 'RAB1A, member RAS oncogene family [Source:MGI Symbol;Acc:MGI:97842]'
		} else {
			cat(paste('"',ExoCarta_Top100.Gene_Symbol.mgi_symbol.NA_attributes[i],'" did not match queried attributes (\'entrezgene\',\'mgi_id\',\'uniprotswissprot\',\'description\') by using R package "biomaRt".',sep=""),"\n")
			
			ExoCarta_Top100.Gene_Symbol.mgi_symbol.NA_attributes.uncurated <- c(ExoCarta_Top100.Gene_Symbol.mgi_symbol.NA_attributes.uncurated, ExoCarta_Top100.Gene_Symbol.mgi_symbol.NA_attributes[i])
		}
		
		manually_curated.mgi_id.table[i,'description'] = trimws(manually_curated.mgi_id.table[i,'description'])
		manually_curated.mgi_id.table[i,'description'] = sapply(strsplit(manually_curated.mgi_id.table[i,'description'],'\\|'), FUN=function(U) paste(unique(sapply(strsplit(unlist(U),'\\[Source:'), FUN=function(V) trimws(unlist(V))[1])),collapse="|"))
	}
	
	if(length(ExoCarta_Top100.Gene_Symbol.mgi_symbol.NA_attributes.uncurated) > 0)
	{
		write.table(ExoCarta_Top100.Gene_Symbol.mgi_symbol.NA_attributes.uncurated, paste('./Data/',Sample,'/10_datasetOverlap/Database/mmu/ExoCarta_All.Protein_Gene2MGI_NA.txt',sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE)
	}
	
	ExoCarta_Top100.Gene_Symbol.mgi_id.table.New = rbind(ExoCarta_Top100.Gene_Symbol.mgi_id.table, manually_curated.mgi_id.table)
}

rownames(ExoCarta_Top100.Gene_Symbol.mgi_id.table.New) = c(1:nrow(ExoCarta_Top100.Gene_Symbol.mgi_id.table.New))
#ExoCarta_Top100.Gene_Symbol.mgi_id.table.New

ExoCarta_Top100.Gene_Symbol.mgi_id.table = ExoCarta_Top100.Gene_Symbol.mgi_id.table.New
#ExoCarta_Top100.Gene_Symbol.mgi_id.table

ExoCarta_Top100.Gene_Symbol.mgi_id.table = ExoCarta_Top100.Gene_Symbol.mgi_id.table[order(match(ExoCarta_Top100.Gene_Symbol.mgi_id.table[,'mgi_symbol'], a)),]
rownames(ExoCarta_Top100.Gene_Symbol.mgi_id.table) = c(1:nrow(ExoCarta_Top100.Gene_Symbol.mgi_id.table))
#ExoCarta_Top100.Gene_Symbol.mgi_id.table

mouse.UniProt.Entry2ProteinNames = sapply(strsplit(ExoCarta_Top100.Gene_Symbol.mgi_id.table[,'uniprotswissprot'],";"), FUN=function(X) paste(unique(sapply(unlist(X), FUN=function(Y) ifelse(length(which(mouse.UniProt.reviewed.Entry %in% unlist(Y)))==0, ifelse(length(which(mouse.UniProt.unreviewed.Entry %in% unlist(Y)))==0, '', mouse.UniProt.unreviewed.ProteinNames[which(mouse.UniProt.unreviewed.Entry %in% unlist(Y))]), mouse.UniProt.reviewed.ProteinNames[which(mouse.UniProt.reviewed.Entry %in% unlist(Y))]))),collapse="|"))
ExoCarta_Top100.Gene_Symbol.mgi_id.table = cbind(ExoCarta_Top100.Gene_Symbol.mgi_id.table, mouse.UniProt.Entry2ProteinNames)
colnames(ExoCarta_Top100.Gene_Symbol.mgi_id.table)[ncol(ExoCarta_Top100.Gene_Symbol.mgi_id.table)] = 'Protein_names'

#write.csv(ExoCarta_Top100.Gene_Symbol.mgi_id.table, paste('./Data/',Sample,'/10_datasetOverlap/Database/mmu/ExoCarta_Top100.Protein_GeneSymbol2MGI.table.csv',sep=""), quote=TRUE, row.names=TRUE)

# add in "Number of times identified" of "Gene Symbol"
ExoCarta_Top100.mmu_format = ExoCarta_Top100
ExoCarta_Top100.mmu_format[,'Gene.Symbol'] = ExoCarta_Top100.Gene_Symbol.mmu_format
ExoCarta_Top100.Gene_Symbol.mgi_id.table.Number_of_times_identified = merge(ExoCarta_Top100.mmu_format, ExoCarta_Top100.Gene_Symbol.mgi_id.table, by.x="Gene.Symbol", by.y="mgi_symbol")
colnames(ExoCarta_Top100.Gene_Symbol.mgi_id.table.Number_of_times_identified) = c('mgi_symbol','Number_of_times_identified',colnames(ExoCarta_Top100.Gene_Symbol.mgi_id.table)[-1])
ExoCarta_Top100.Gene_Symbol.mgi_id.table.Number_of_times_identified = ExoCarta_Top100.Gene_Symbol.mgi_id.table.Number_of_times_identified[order(match(ExoCarta_Top100.Gene_Symbol.mgi_id.table.Number_of_times_identified[,'mgi_symbol'], a)),]
rownames(ExoCarta_Top100.Gene_Symbol.mgi_id.table.Number_of_times_identified) = order(match(a, ExoCarta_Top100.Gene_Symbol.mgi_id.table.Number_of_times_identified[,'mgi_symbol']))
ExoCarta_Top100.Gene_Symbol.mgi_id.table.Number_of_times_identified = ExoCarta_Top100.Gene_Symbol.mgi_id.table.Number_of_times_identified[order(as.numeric(rownames(ExoCarta_Top100.Gene_Symbol.mgi_id.table.Number_of_times_identified))),]
#ExoCarta_Top100.Gene_Symbol.mgi_id.table.Number_of_times_identified

# merge Hist2h4a and Hist1h4b (Hist1h4b: Gene Overview MyGene.info: HIST2H4A (http://www.informatics.jax.org/marker/MGI:2448420))
#Hist2h4a_Hist1h4b.Number_of_times_identified = sum(ExoCarta_Top100.Gene_Symbol.mgi_id.table.Number_of_times_identified[which(ExoCarta_Top100.Gene_Symbol.mgi_id.table.Number_of_times_identified[,'mgi_symbol']=='Hist1h4b'),'Number_of_times_identified'])
#ExoCarta_Top100.Gene_Symbol.mgi_id.table.Number_of_times_identified[which(ExoCarta_Top100.Gene_Symbol.mgi_id.table.Number_of_times_identified[,'mgi_symbol']=='Hist1h4b'),'Number_of_times_identified'] = Hist2h4a_Hist1h4b.Number_of_times_identified

#ExoCarta_Top100.Gene_Symbol.mgi_id.table.Number_of_times_identified = ExoCarta_Top100.Gene_Symbol.mgi_id.table.Number_of_times_identified[-which(ExoCarta_Top100.Gene_Symbol.mgi_id.table.Number_of_times_identified[,'mgi_symbol']=='Hist1h4b')[which(which(ExoCarta_Top100.Gene_Symbol.mgi_id.table.Number_of_times_identified[,'mgi_symbol']=='Hist1h4b') %in% ExoCarta_Top100.Gene_Symbol.mmu_format.Hist2h4a_idx)],]
ExoCarta_Top100.Gene_Symbol.mgi_id.table.Number_of_times_identified = unique(ExoCarta_Top100.Gene_Symbol.mgi_id.table.Number_of_times_identified)
#ExoCarta_Top100.Gene_Symbol.mgi_id.table.Number_of_times_identified = ExoCarta_Top100.Gene_Symbol.mgi_id.table.Number_of_times_identified[c(which(ExoCarta_Top100.Gene_Symbol.mgi_id.table.Number_of_times_identified[,'Number_of_times_identified'] >= Hist2h4a_Hist1h4b.Number_of_times_identified), which(ExoCarta_Top100.Gene_Symbol.mgi_id.table.Number_of_times_identified[,'Number_of_times_identified'] < Hist2h4a_Hist1h4b.Number_of_times_identified)),]
#ExoCarta_Top100.Gene_Symbol.mgi_id.table.Number_of_times_identified

#write.csv(ExoCarta_Top100.Gene_Symbol.mgi_id.table.Number_of_times_identified, paste('./Data/',Sample,'/10_datasetOverlap/Database/mmu/ExoCarta_Top100.Protein_GeneSymbol2MGI.table_IdentificationNumber.csv',sep=""), quote=TRUE, row.names=TRUE)

# --- write xlsx file by openxlsx package ---------------------------------------------------------
# Importing a big xlsx file into R? https://stackoverflow.com/questions/19147884/importing-a-big-xlsx-file-into-r/43118530#43118530
# Building R for Windows: https://cran.r-project.org/bin/windows/Rtools/
#Sys.setenv("R_ZIPCMD" = "C:/Rtools/bin/zip.exe")	# path to zip.exe

wb <- openxlsx::createWorkbook(paste('./Data/',Sample,'/10_datasetOverlap/Database/mmu/ExoCarta_Top100.Protein_GeneSymbol2MGI.table.xlsx',sep=""))
modifyBaseFont(wb, fontSize=12, fontColour="black", fontName="Times New Roman")

for(s in 1:2)
{
	if(s == 1) {
		sheet = 'ExoCarta_Top100'
		sheetData = ExoCarta_Top100.Gene_Symbol.mgi_id.table
	} else if(s == 2) {
		sheet = 'ExoCarta_Top100 (with number)'
		sheetData = ExoCarta_Top100.Gene_Symbol.mgi_id.table.Number_of_times_identified
	}
	
	if(!is.null(sheetData))
	{
		if(nrow(sheetData) > 0)
		{
			addWorksheet(wb, sheet)
			writeData(wb, sheet, sheetData, colNames=TRUE, rowNames=TRUE, keepNA=FALSE)
			setColWidths(wb, sheet, cols=c(1:(1+ncol(sheetData))), widths=c(8.5,rep(15,ncol(sheetData))))
			setRowHeights(wb, sheet, rows=c(1:(1+nrow(sheetData))), heights=16.5)
			freezePane(wb, sheet, firstActiveRow=2)
			style <- createStyle(halign="left", valign="center")
			addStyle(wb, sheet, style, rows=c(1:(1+nrow(sheetData))), cols=c(1:(1+ncol(sheetData))), gridExpand=TRUE)
			style <- createStyle(textDecoration="bold")
			addStyle(wb, sheet, style, rows=1, cols=c(1:(1+ncol(sheetData))), stack=TRUE)
			
			right_align.idx = which(colnames(sheetData) %in% c('Number_of_times_identified'))
			if(length(right_align.idx) > 0) {
				style <- createStyle(halign="right", valign="center")
				addStyle(wb, sheet, style, rows=c(2:(1+nrow(sheetData))), cols=1+right_align.idx, gridExpand=TRUE, stack=TRUE)
			}
		}
	}
}

if(length(wb$sheet_names) > 0) {
	openxlsx::saveWorkbook(wb, file=paste('./Data/',Sample,'/10_datasetOverlap/Database/mmu/ExoCarta_Top100.Protein_GeneSymbol2MGI.table.xlsx',sep=""), overwrite=TRUE)
}

