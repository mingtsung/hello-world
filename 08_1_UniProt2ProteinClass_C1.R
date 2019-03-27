rm(list=ls())	# remove all variables in workspace

setwd("D:/Exo")

Sample <- "pcMSC_N230"

mainDir <- paste('./Data/',Sample,sep="")
subDir <- "08_UniProt2ProteinClass"
if (file.exists(file.path(mainDir, subDir))){
} else {
	dir.create(file.path(mainDir, subDir))
}
mainDir <- paste('./Data/',Sample,'/08_UniProt2ProteinClass',sep="")
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

#install.packages("gtools")				# gtools: Various R Programming Tools
library(gtools)							# mixedsort(), mixedorder()

#install.packages("ggplot2")			# ggplot2: Create Elegant Data Visualisations Using the Grammar of Graphics
library(ggplot2)

#install.packages("hash")				# hash: Full feature implementation of hash/associated arrays/dictionaries
library(hash)


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


### Code reference ##########################################################################################
#Round values while preserve their rounded sum in R
#https://biostatmatt.com/archives/2902
round_preserve_sum <- function(x, digits = 0) {
	up <- 10 ^ digits
	x <- x * up
	y <- floor(x)
	indices <- tail(order(x-y), round(sum(x)) - sum(y))
	y[indices] <- y[indices] + 1
	y / up
}
#############################################################################################################


if(Sample == "pcMSC_N230"){
	batch <- c('P3_3rd_181022','P4_4th_181022','P4_5th_181022')
}else{
	print(paste('Please provide filename of Sample ',Sample,sep=""))
}


# === PANTHER (Protein ANalysis THrough Evolutionary Relationships) Classification System =======================================
SpeciesName <- c("Homo sapiens","Mus musculus","Rattus norvegicus")
SpeciesName.general <- c("human","mouse","rat")
SpeciesName.code <- c("hsa","mmu","rno")

# --- input Protein Class -------------------------------------------------------------------------------------------------------
Protein_Class <- read.table("./Database/PANTHER/Ontologies/Protein_Class_13.0", sep="\t", quote="", header=FALSE, stringsAsFactors=FALSE)
Protein_class_relationship <- read.table("./Database/PANTHER/Ontologies/Protein_class_relationship", sep="\t", quote="", header=FALSE, stringsAsFactors=FALSE)

Protein_Class.id <- trimws(Protein_Class[,1])
Protein_Class.level <- trimws(Protein_Class[,2])
Protein_Class.name <- trimws(Protein_Class[,3])
Protein_Class.description <- trimws(Protein_Class[,4])

Protein_Class.id2level.hashtable <- hash()
Protein_Class.id2name.hashtable <- hash()
Protein_Class.id2description.hashtable <- hash()
for(hh in 1:length(Protein_Class.id))
{
	Protein_Class.id2level.hashtable[[Protein_Class.id[hh]]] = Protein_Class.level[hh]
	Protein_Class.id2name.hashtable[[Protein_Class.id[hh]]] = Protein_Class.name[hh]
	Protein_Class.id2description.hashtable[[Protein_Class.id[hh]]] = Protein_Class.description[hh]
}

Protein_Class.name2description.hashtable <- hash()
for(hh in 1:length(Protein_Class.name))
{
	Protein_Class.name2description.hashtable[[Protein_Class.name[hh]]] = Protein_Class.description[hh]
}

#fiveNum.idx = which(sapply(strsplit(Protein_Class.level,'\\.'), FUN=function(X) length(X)==5))
#sixNum.idx = which(sapply(strsplit(Protein_Class.level,'\\.'), FUN=function(X) length(X)==6))

Protein_Class.level_0.idx = which(sapply(strsplit(Protein_Class.level,'\\.'), FUN=function(X) sum(X!='00')==1))
Protein_Class.level_1.idx = which(sapply(strsplit(Protein_Class.level,'\\.'), FUN=function(X) sum(X!='00')==2))
Protein_Class.level_2.idx = which(sapply(strsplit(Protein_Class.level,'\\.'), FUN=function(X) sum(X!='00')==3))
Protein_Class.level_3.idx = which(sapply(strsplit(Protein_Class.level,'\\.'), FUN=function(X) sum(X!='00')==4))
Protein_Class.level_4.idx = which(sapply(strsplit(Protein_Class.level,'\\.'), FUN=function(X) sum(X!='00')==5))
Protein_Class.level_5.idx = which(sapply(strsplit(Protein_Class.level,'\\.'), FUN=function(X) sum(X!='00')==6))
#length(Protein_Class.level_0.idx)
#length(Protein_Class.level_1.idx)
#length(Protein_Class.level_2.idx)
#length(Protein_Class.level_3.idx)
#length(Protein_Class.level_4.idx)
#length(Protein_Class.level_5.idx)
#length(Protein_Class.level_0.idx) + length(Protein_Class.level_1.idx) + length(Protein_Class.level_2.idx) + length(Protein_Class.level_3.idx) + length(Protein_Class.level_4.idx) + length(Protein_Class.level_5.idx)
if(length(c(Protein_Class.level_0.idx, Protein_Class.level_1.idx, Protein_Class.level_2.idx, Protein_Class.level_3.idx, Protein_Class.level_4.idx, Protein_Class.level_5.idx)) != nrow(Protein_Class))
{
	cat('[Note] Sum of each level (',length(Protein_Class.level_0.idx),'+',length(Protein_Class.level_1.idx),'+',length(Protein_Class.level_2.idx),'+',length(Protein_Class.level_3.idx),'+',length(Protein_Class.level_4.idx),'+',length(Protein_Class.level_5.idx),'=',length(c(Protein_Class.level_0.idx, Protein_Class.level_1.idx, Protein_Class.level_2.idx, Protein_Class.level_3.idx, Protein_Class.level_4.idx, Protein_Class.level_5.idx)),') is not the same as total protein class number (',nrow(Protein_Class),').\n', sep="")
}

# /// level_0 ///
Protein_Class.level_0 <- Protein_Class[Protein_Class.level_0.idx,,drop=FALSE]
Protein_Class.level_0.id <- Protein_Class.level_0[,1]
Protein_Class.level_0.level <- Protein_Class.level_0[,2]
Protein_Class.level_0.name <- Protein_Class.level_0[,3]
Protein_Class.level_0.description <- Protein_Class.level_0[,4]

Protein_Class.level_0.id2level.hashtable <- hash()
Protein_Class.level_0.id2name.hashtable <- hash()
Protein_Class.level_0.id2description.hashtable <- hash()
for(hh in 1:length(Protein_Class.level_0.id))
{
	Protein_Class.level_0.id2level.hashtable[[Protein_Class.level_0.id[hh]]] = Protein_Class.level_0.level[hh]
	Protein_Class.level_0.id2name.hashtable[[Protein_Class.level_0.id[hh]]] = Protein_Class.level_0.name[hh]
	Protein_Class.level_0.id2description.hashtable[[Protein_Class.level_0.id[hh]]] = Protein_Class.level_0.description[hh]
}

Protein_Class.level_0.name2description.hashtable <- hash()
for(hh in 1:length(Protein_Class.level_0.name))
{
	Protein_Class.level_0.name2description.hashtable[[Protein_Class.level_0.name[hh]]] = Protein_Class.level_0.description[hh]
}

# /// level_1 ///
#Protein_Class.level_1.idx <- which(Protein_Class.level %in% grep('.00.00.00.00',grep('.00.00.00',Protein_Class.level[fiveNum.idx],value=TRUE,invert=FALSE),value=TRUE,invert=TRUE))
Protein_Class.level_1 <- Protein_Class[Protein_Class.level_1.idx,,drop=FALSE]
Protein_Class.level_1.id <- Protein_Class.level_1[,1]
Protein_Class.level_1.level <- Protein_Class.level_1[,2]
Protein_Class.level_1.name <- Protein_Class.level_1[,3]
Protein_Class.level_1.description <- Protein_Class.level_1[,4]

Protein_Class.level_1.id2level.hashtable <- hash()
Protein_Class.level_1.id2name.hashtable <- hash()
Protein_Class.level_1.id2description.hashtable <- hash()
for(hh in 1:length(Protein_Class.level_1.id))
{
	Protein_Class.level_1.id2level.hashtable[[Protein_Class.level_1.id[hh]]] = Protein_Class.level_1.level[hh]
	Protein_Class.level_1.id2name.hashtable[[Protein_Class.level_1.id[hh]]] = Protein_Class.level_1.name[hh]
	Protein_Class.level_1.id2description.hashtable[[Protein_Class.level_1.id[hh]]] = Protein_Class.level_1.description[hh]
}

Protein_Class.level_1.name2description.hashtable <- hash()
for(hh in 1:length(Protein_Class.level_1.name))
{
	Protein_Class.level_1.name2description.hashtable[[Protein_Class.level_1.name[hh]]] = Protein_Class.level_1.description[hh]
}

# /// level_2 ///
#Protein_Class.level_2.idx <- which(Protein_Class.level %in% grep('.00.00.00',grep('.00.00',Protein_Class.level[fiveNum.idx],value=TRUE,invert=FALSE),value=TRUE,invert=TRUE))
Protein_Class.level_2 <- Protein_Class[Protein_Class.level_2.idx,,drop=FALSE]
Protein_Class.level_2.id <- Protein_Class.level_2[,1]
Protein_Class.level_2.level <- Protein_Class.level_2[,2]
Protein_Class.level_2.name <- Protein_Class.level_2[,3]
Protein_Class.level_2.description <- Protein_Class.level_2[,4]

Protein_Class.level_2.id2level.hashtable <- hash()
Protein_Class.level_2.id2name.hashtable <- hash()
Protein_Class.level_2.id2description.hashtable <- hash()
for(hh in 1:length(Protein_Class.level_2.id))
{
	Protein_Class.level_2.id2level.hashtable[[Protein_Class.level_2.id[hh]]] = Protein_Class.level_2.level[hh]
	Protein_Class.level_2.id2name.hashtable[[Protein_Class.level_2.id[hh]]] = Protein_Class.level_2.name[hh]
	Protein_Class.level_2.id2description.hashtable[[Protein_Class.level_2.id[hh]]] = Protein_Class.level_2.description[hh]
}

Protein_Class.level_2.name2description.hashtable <- hash()
for(hh in 1:length(Protein_Class.level_2.name))
{
	Protein_Class.level_2.name2description.hashtable[[Protein_Class.level_2.name[hh]]] = Protein_Class.level_2.description[hh]
}

# /// level_3 ///
#Protein_Class.level_3.idx <- which(Protein_Class.level %in% grep('.00.00',grep('.00',Protein_Class.level[fiveNum.idx],value=TRUE,invert=FALSE),value=TRUE,invert=TRUE))
Protein_Class.level_3 <- Protein_Class[Protein_Class.level_3.idx,,drop=FALSE]
Protein_Class.level_3.id <- Protein_Class.level_3[,1]
Protein_Class.level_3.level <- Protein_Class.level_3[,2]
Protein_Class.level_3.name <- Protein_Class.level_3[,3]
Protein_Class.level_3.description <- Protein_Class.level_3[,4]

Protein_Class.level_3.id2level.hashtable <- hash()
Protein_Class.level_3.id2name.hashtable <- hash()
Protein_Class.level_3.id2description.hashtable <- hash()
for(hh in 1:length(Protein_Class.level_3.id))
{
	Protein_Class.level_3.id2level.hashtable[[Protein_Class.level_3.id[hh]]] = Protein_Class.level_3.level[hh]
	Protein_Class.level_3.id2name.hashtable[[Protein_Class.level_3.id[hh]]] = Protein_Class.level_3.name[hh]
	Protein_Class.level_3.id2description.hashtable[[Protein_Class.level_3.id[hh]]] = Protein_Class.level_3.description[hh]
}

Protein_Class.level_3.name2description.hashtable <- hash()
for(hh in 1:length(Protein_Class.level_3.name))
{
	Protein_Class.level_3.name2description.hashtable[[Protein_Class.level_3.name[hh]]] = Protein_Class.level_3.description[hh]
}

# /// level_4 ///
Protein_Class.level_4 <- Protein_Class[Protein_Class.level_4.idx,,drop=FALSE]
Protein_Class.level_4.id <- Protein_Class.level_4[,1]
Protein_Class.level_4.level <- Protein_Class.level_4[,2]
Protein_Class.level_4.name <- Protein_Class.level_4[,3]
Protein_Class.level_4.description <- Protein_Class.level_4[,4]

Protein_Class.level_4.id2level.hashtable <- hash()
Protein_Class.level_4.id2name.hashtable <- hash()
Protein_Class.level_4.id2description.hashtable <- hash()
for(hh in 1:length(Protein_Class.level_4.id))
{
	Protein_Class.level_4.id2level.hashtable[[Protein_Class.level_4.id[hh]]] = Protein_Class.level_4.level[hh]
	Protein_Class.level_4.id2name.hashtable[[Protein_Class.level_4.id[hh]]] = Protein_Class.level_4.name[hh]
	Protein_Class.level_4.id2description.hashtable[[Protein_Class.level_4.id[hh]]] = Protein_Class.level_4.description[hh]
}

Protein_Class.level_4.name2description.hashtable <- hash()
for(hh in 1:length(Protein_Class.level_4.name))
{
	Protein_Class.level_4.name2description.hashtable[[Protein_Class.level_4.name[hh]]] = Protein_Class.level_4.description[hh]
}

# /// level_5 ///
Protein_Class.level_5 <- Protein_Class[Protein_Class.level_5.idx,,drop=FALSE]
Protein_Class.level_5.id <- Protein_Class.level_5[,1]
Protein_Class.level_5.level <- Protein_Class.level_5[,2]
Protein_Class.level_5.name <- Protein_Class.level_5[,3]
Protein_Class.level_5.description <- Protein_Class.level_5[,4]

Protein_Class.level_5.id2level.hashtable <- hash()
Protein_Class.level_5.id2name.hashtable <- hash()
Protein_Class.level_5.id2description.hashtable <- hash()
for(hh in 1:length(Protein_Class.level_5.id))
{
	Protein_Class.level_5.id2level.hashtable[[Protein_Class.level_5.id[hh]]] = Protein_Class.level_5.level[hh]
	Protein_Class.level_5.id2name.hashtable[[Protein_Class.level_5.id[hh]]] = Protein_Class.level_5.name[hh]
	Protein_Class.level_5.id2description.hashtable[[Protein_Class.level_5.id[hh]]] = Protein_Class.level_5.description[hh]
}

Protein_Class.level_5.name2description.hashtable <- hash()
for(hh in 1:length(Protein_Class.level_5.name))
{
	Protein_Class.level_5.name2description.hashtable[[Protein_Class.level_5.name[hh]]] = Protein_Class.level_5.description[hh]
}

Protein_class_relationship.childid <- trimws(Protein_class_relationship[,1])
Protein_class_relationship.childname <- trimws(Protein_class_relationship[,2])
Protein_class_relationship.parentid <- trimws(Protein_class_relationship[,3])
Protein_class_relationship.parentname <- trimws(Protein_class_relationship[,4])
Protein_class_relationship.childorder <- trimws(Protein_class_relationship[,5])

# --- input sequence_classifications --------------------------------------------------------------------------------------------
for(sp in 1:length(SpeciesName.general))
{
	PTHR13 <- read.table(paste('./Database/PANTHER/sequence_classifications/13.1/PTHR13.1_',SpeciesName.general[sp],sep=""), sep="\t", quote="", header=FALSE, stringsAsFactors=FALSE, comment.char="")
	
	# './Database/PANTHER/sequence_classifications/13.1/README'
	PTHR13.GeneID <- trimws(PTHR13[,1])			# 1) Gene Identifier (format:  organism|gene id source:gene id|protein id source:protein id)
	#PTHR13.ProteinID <- trimws(PTHR13[,2])		# 2) Protein ID - currently empty. The protein ids can be retrieved from above
	PTHR13.SubfamilyID <- trimws(PTHR13[,3])	# 3) PANTHER SF ID ?for example, PTHR12213:SF6.  ":SF" indicates the subfamily ID
	PTHR13.FamilyName <- trimws(PTHR13[,4])		# 4) PANTHER Family Name
	PTHR13.SubfamilyName <- trimws(PTHR13[,5])	# 5) PANTHER Subfamily Name
	PTHR13.MF <- trimws(PTHR13[,6])				# 6) PANTHER Molecular function*
	PTHR13.BP <- trimws(PTHR13[,7])				# 7) PANTHER Biological process*
	PTHR13.CC <- trimws(PTHR13[,8])				# 8) Cellular components*:  PANTHER GO slim cellular component terms assigned to families and subfamilies
	PTHR13.ProteinClass <- trimws(PTHR13[,9])	# 9) Protein class*:  PANTHER protein class terms assigned to families and subfamilies
	PTHR13.Pathway <- trimws(PTHR13[,10])		# 10) Pathway**
	
	if(SpeciesName.general[sp] == "human") {
		PTHR13_human.GeneID <- PTHR13.GeneID
		#PTHR13_human.ProteinID <- PTHR13.ProteinID
		PTHR13_human.SubfamilyID <- PTHR13.SubfamilyID
		PTHR13_human.FamilyName <- PTHR13.FamilyName
		PTHR13_human.SubfamilyName <- PTHR13.SubfamilyName
		PTHR13_human.MF <- PTHR13.MF
		PTHR13_human.BP <- PTHR13.BP
		PTHR13_human.CC <- PTHR13.CC
		PTHR13_human.ProteinClass <- PTHR13.ProteinClass
		PTHR13_human.Pathway <- PTHR13.Pathway
	} else if(SpeciesName.general[sp] == "mouse") {
		PTHR13_mouse.GeneID <- PTHR13.GeneID
		#PTHR13_mouse.ProteinID <- PTHR13.ProteinID
		PTHR13_mouse.SubfamilyID <- PTHR13.SubfamilyID
		PTHR13_mouse.FamilyName <- PTHR13.FamilyName
		PTHR13_mouse.SubfamilyName <- PTHR13.SubfamilyName
		PTHR13_mouse.MF <- PTHR13.MF
		PTHR13_mouse.BP <- PTHR13.BP
		PTHR13_mouse.CC <- PTHR13.CC
		PTHR13_mouse.ProteinClass <- PTHR13.ProteinClass
		PTHR13_mouse.Pathway <- PTHR13.Pathway
	} else if(SpeciesName.general[sp] == "rat") {
		PTHR13_rat.GeneID <- PTHR13.GeneID
		#PTHR13_rat.ProteinID <- PTHR13.ProteinID
		PTHR13_rat.SubfamilyID <- PTHR13.SubfamilyID
		PTHR13_rat.FamilyName <- PTHR13.FamilyName
		PTHR13_rat.SubfamilyName <- PTHR13.SubfamilyName
		PTHR13_rat.MF <- PTHR13.MF
		PTHR13_rat.BP <- PTHR13.BP
		PTHR13_rat.CC <- PTHR13.CC
		PTHR13_rat.ProteinClass <- PTHR13.ProteinClass
		PTHR13_rat.Pathway <- PTHR13.Pathway
	}
}

# +++ human +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
PTHR13_human.GeneID.species <- sapply(strsplit(PTHR13_human.GeneID,'\\|'), FUN=function(X) X[1])
PTHR13_human.GeneID.UniProtKB <- gsub('UniProtKB=','',sapply(strsplit(PTHR13_human.GeneID,'\\|'), FUN=function(X) X[3]))

PTHR13_human.GeneID.HGNCid <- gsub('HGNC=','HGNC:',sapply(strsplit(PTHR13_human.GeneID,'\\|'), FUN=function(X) X[2]))
PTHR13_human.GeneID.HGNCid.NotHGNC.idx <- grep('HGNC:',PTHR13_human.GeneID.HGNCid,invert=TRUE,value=FALSE)
PTHR13_human.GeneID.HGNCid.NotHGNC.id <- grep('HGNC:',PTHR13_human.GeneID.HGNCid,invert=TRUE,value=TRUE)

PTHR13_human.GeneID.EnsemblGene.idx <- grep('Ensembl=',PTHR13_human.GeneID.HGNCid,invert=FALSE,value=FALSE)
PTHR13_human.GeneID.EnsemblGene.id <- grep('Ensembl=',PTHR13_human.GeneID.HGNCid,invert=FALSE,value=TRUE)
PTHR13_human.GeneID.EnsemblGene <- gsub('Ensembl=','',PTHR13_human.GeneID.EnsemblGene.id)

PTHR13_human.GeneID.NCBIGeneID.idx <- grep('GeneID=',PTHR13_human.GeneID.HGNCid,invert=FALSE,value=FALSE)
PTHR13_human.GeneID.NCBIGeneID.id <- grep('GeneID=',PTHR13_human.GeneID.HGNCid,invert=FALSE,value=TRUE)
PTHR13_human.GeneID.NCBIGeneID <- gsub('GeneID=','',PTHR13_human.GeneID.NCBIGeneID.id)

PTHR13_human.GeneID.NCBIGeneSymbol.idx <- grep('Gene=',PTHR13_human.GeneID.HGNCid,invert=FALSE,value=FALSE)
PTHR13_human.GeneID.NCBIGeneSymbol.id <- grep('Gene=',PTHR13_human.GeneID.HGNCid,invert=FALSE,value=TRUE)
PTHR13_human.GeneID.NCBIGeneSymbol <- gsub('Gene=','',PTHR13_human.GeneID.NCBIGeneSymbol.id)

# Note: Still need to check
if(length(PTHR13_human.GeneID.EnsemblGene) + length(PTHR13_human.GeneID.NCBIGeneID) + length(PTHR13_human.GeneID.NCBIGeneSymbol) != length(PTHR13_human.GeneID.HGNCid.NotHGNC.id))
{
#	cat('Please check Gene Identifier (the first column) in \"',paste('./Database/PANTHER/sequence_classifications/13.1/PTHR13.1_human',sep=""),'\" file.','\n','Non-HGNC id is not only \'Ensembl=\', \'GeneID=\', \'Gene=\'.','\n',sep="")
}

ensembl = ensembl.hsa
#listFilters(ensembl)[,c('name','description')]
#listAttributes(ensembl)[,c('name','description')]

a = PTHR13_human.GeneID.EnsemblGene
b = ensembl.hsa
PTHR13_human.GeneID.EnsemblGene.attributes <- getBM(attributes=c('ensembl_gene_id','entrezgene','uniprotswissprot'), filters='ensembl_gene_id', values=a, mart=b)
PTHR13_human.GeneID.EnsemblGene.attributes = unique(PTHR13_human.GeneID.EnsemblGene.attributes)
PTHR13_human.GeneID.EnsemblGene2HGNC <- getBM(attributes=c('ensembl_gene_id','hgnc_id','hgnc_symbol'), filters='ensembl_gene_id', values=a, mart=b)
PTHR13_human.GeneID.EnsemblGene2HGNC.table = merge(PTHR13_human.GeneID.EnsemblGene2HGNC, PTHR13_human.GeneID.EnsemblGene.attributes, by.x="ensembl_gene_id", by.y="ensembl_gene_id")
if(length(which(is.na(PTHR13_human.GeneID.EnsemblGene2HGNC.table[,'entrezgene']))) > 0)		PTHR13_human.GeneID.EnsemblGene2HGNC.table[which(is.na(PTHR13_human.GeneID.EnsemblGene2HGNC.table[,'entrezgene'])),'entrezgene'] = ''
PTHR13_human.GeneID.EnsemblGene2HGNC.table.rm_idx = which(PTHR13_human.GeneID.EnsemblGene2HGNC.table[,'hgnc_id']=='' & PTHR13_human.GeneID.EnsemblGene2HGNC.table[,'hgnc_symbol']=='' & PTHR13_human.GeneID.EnsemblGene2HGNC.table[,'entrezgene']=='' & PTHR13_human.GeneID.EnsemblGene2HGNC.table[,'uniprotswissprot']=='')
if(length(PTHR13_human.GeneID.EnsemblGene2HGNC.table.rm_idx) > 0)		PTHR13_human.GeneID.EnsemblGene2HGNC.table = PTHR13_human.GeneID.EnsemblGene2HGNC.table[-PTHR13_human.GeneID.EnsemblGene2HGNC.table.rm_idx,,drop=FALSE]
PTHR13_human.GeneID.EnsemblGene2HGNC.table = unique(PTHR13_human.GeneID.EnsemblGene2HGNC.table)
PTHR13_human.GeneID.EnsemblGene2HGNC.subtable = unique(PTHR13_human.GeneID.EnsemblGene2HGNC.table[,c('ensembl_gene_id','hgnc_id')])
PTHR13_human.GeneID.EnsemblGene2HGNC.subtable.rm_idx = which(PTHR13_human.GeneID.EnsemblGene2HGNC.subtable[,'hgnc_id']=='')
if(length(PTHR13_human.GeneID.EnsemblGene2HGNC.subtable.rm_idx) > 0)	PTHR13_human.GeneID.EnsemblGene2HGNC.subtable = PTHR13_human.GeneID.EnsemblGene2HGNC.subtable[-PTHR13_human.GeneID.EnsemblGene2HGNC.subtable.rm_idx,,drop=FALSE]
PTHR13_human.GeneID.EnsemblGene2HGNC.subtable = unique(PTHR13_human.GeneID.EnsemblGene2HGNC.subtable)
if(length(PTHR13_human.GeneID.EnsemblGene2HGNC.subtable[,'ensembl_gene_id']) > 0)
{
	PTHR13_human.GeneID.EnsemblGene2HGNC.subtable[,'ensembl_gene_id'] = paste('Ensembl=',PTHR13_human.GeneID.EnsemblGene2HGNC.subtable[,'ensembl_gene_id'],sep="")
	
	for(i in 1:length(PTHR13_human.GeneID.EnsemblGene2HGNC.subtable[,'ensembl_gene_id']))
	{
		PTHR13_human.GeneID.HGNCid[which(PTHR13_human.GeneID.HGNCid %in% PTHR13_human.GeneID.EnsemblGene2HGNC.subtable[i,'ensembl_gene_id'])] = PTHR13_human.GeneID.EnsemblGene2HGNC.subtable[i,'hgnc_id']
	}
}

a = PTHR13_human.GeneID.NCBIGeneID
b = ensembl.hsa
PTHR13_human.GeneID.NCBIGeneID.attributes <- getBM(attributes=c('entrezgene','ensembl_gene_id','uniprotswissprot'), filters='entrezgene', values=a, mart=b)
PTHR13_human.GeneID.NCBIGeneID.attributes = unique(PTHR13_human.GeneID.NCBIGeneID.attributes)
PTHR13_human.GeneID.NCBIGeneID2HGNC <- getBM(attributes=c('entrezgene','hgnc_id','hgnc_symbol'), filters='entrezgene', values=a, mart=b)
PTHR13_human.GeneID.NCBIGeneID2HGNC.table = merge(PTHR13_human.GeneID.NCBIGeneID2HGNC, PTHR13_human.GeneID.NCBIGeneID.attributes, by.x="entrezgene", by.y="entrezgene")
if(length(which(is.na(PTHR13_human.GeneID.NCBIGeneID2HGNC.table[,'entrezgene']))) > 0)		PTHR13_human.GeneID.NCBIGeneID2HGNC.table[which(is.na(PTHR13_human.GeneID.NCBIGeneID2HGNC.table[,'entrezgene'])),'entrezgene'] = ''
PTHR13_human.GeneID.NCBIGeneID2HGNC.table.rm_idx = which(PTHR13_human.GeneID.NCBIGeneID2HGNC.table[,'hgnc_id']=='' & PTHR13_human.GeneID.NCBIGeneID2HGNC.table[,'hgnc_symbol']=='' & PTHR13_human.GeneID.NCBIGeneID2HGNC.table[,'ensembl_gene_id']=='' & PTHR13_human.GeneID.NCBIGeneID2HGNC.table[,'uniprotswissprot']=='')
if(length(PTHR13_human.GeneID.NCBIGeneID2HGNC.table.rm_idx) > 0)		PTHR13_human.GeneID.NCBIGeneID2HGNC.table = PTHR13_human.GeneID.NCBIGeneID2HGNC.table[-PTHR13_human.GeneID.NCBIGeneID2HGNC.table.rm_idx,,drop=FALSE]
PTHR13_human.GeneID.NCBIGeneID2HGNC.table = unique(PTHR13_human.GeneID.NCBIGeneID2HGNC.table)
PTHR13_human.GeneID.NCBIGeneID2HGNC.subtable = unique(PTHR13_human.GeneID.NCBIGeneID2HGNC.table[,c('entrezgene','hgnc_id')])
PTHR13_human.GeneID.NCBIGeneID2HGNC.subtable.rm_idx = which(PTHR13_human.GeneID.NCBIGeneID2HGNC.subtable[,'hgnc_id']=='')
if(length(PTHR13_human.GeneID.NCBIGeneID2HGNC.subtable.rm_idx) > 0)	PTHR13_human.GeneID.NCBIGeneID2HGNC.subtable = PTHR13_human.GeneID.NCBIGeneID2HGNC.subtable[-PTHR13_human.GeneID.NCBIGeneID2HGNC.subtable.rm_idx,,drop=FALSE]
PTHR13_human.GeneID.NCBIGeneID2HGNC.subtable = unique(PTHR13_human.GeneID.NCBIGeneID2HGNC.subtable)
if(length(PTHR13_human.GeneID.NCBIGeneID2HGNC.subtable[,'entrezgene']) > 0)
{
	PTHR13_human.GeneID.NCBIGeneID2HGNC.subtable[,'entrezgene'] = paste('GeneID=',PTHR13_human.GeneID.NCBIGeneID2HGNC.subtable[,'entrezgene'],sep="")
	
	for(i in 1:length(PTHR13_human.GeneID.NCBIGeneID2HGNC.subtable[,'entrezgene']))
	{
		PTHR13_human.GeneID.HGNCid[which(PTHR13_human.GeneID.HGNCid %in% PTHR13_human.GeneID.NCBIGeneID2HGNC.subtable[i,'entrezgene'])] = PTHR13_human.GeneID.NCBIGeneID2HGNC.subtable[i,'hgnc_id']
	}
}

a = PTHR13_human.GeneID.NCBIGeneSymbol
b = ensembl.hsa
PTHR13_human.GeneID.NCBIGeneSymbol.attributes <- getBM(attributes=c('hgnc_symbol','ensembl_gene_id','uniprotswissprot'), filters='hgnc_symbol', values=a, mart=b)
PTHR13_human.GeneID.NCBIGeneSymbol.attributes = unique(PTHR13_human.GeneID.NCBIGeneSymbol.attributes)
PTHR13_human.GeneID.NCBIGeneSymbol2HGNC <- getBM(attributes=c('hgnc_symbol','hgnc_id','entrezgene'), filters='hgnc_symbol', values=a, mart=b)
PTHR13_human.GeneID.NCBIGeneSymbol2HGNC.table = merge(PTHR13_human.GeneID.NCBIGeneSymbol2HGNC, PTHR13_human.GeneID.NCBIGeneSymbol.attributes, by.x="hgnc_symbol", by.y="hgnc_symbol")
if(length(which(is.na(PTHR13_human.GeneID.NCBIGeneSymbol2HGNC.table[,'entrezgene']))) > 0)	PTHR13_human.GeneID.NCBIGeneSymbol2HGNC.table[which(is.na(PTHR13_human.GeneID.NCBIGeneSymbol2HGNC.table[,'entrezgene'])),'entrezgene'] = ''
PTHR13_human.GeneID.NCBIGeneSymbol2HGNC.table.rm_idx = which(PTHR13_human.GeneID.NCBIGeneSymbol2HGNC.table[,'hgnc_id']=='' & PTHR13_human.GeneID.NCBIGeneSymbol2HGNC.table[,'entrezgene']=='' & PTHR13_human.GeneID.NCBIGeneSymbol2HGNC.table[,'ensembl_gene_id']=='' & PTHR13_human.GeneID.NCBIGeneSymbol2HGNC.table[,'uniprotswissprot']=='')
if(length(PTHR13_human.GeneID.NCBIGeneSymbol2HGNC.table.rm_idx) > 0)		PTHR13_human.GeneID.NCBIGeneSymbol2HGNC.table = PTHR13_human.GeneID.NCBIGeneSymbol2HGNC.table[-PTHR13_human.GeneID.NCBIGeneSymbol2HGNC.table.rm_idx,,drop=FALSE]
PTHR13_human.GeneID.NCBIGeneSymbol2HGNC.table = unique(PTHR13_human.GeneID.NCBIGeneSymbol2HGNC.table)
PTHR13_human.GeneID.NCBIGeneSymbol2HGNC.subtable = unique(PTHR13_human.GeneID.NCBIGeneSymbol2HGNC.table[,c('hgnc_symbol','hgnc_id')])
PTHR13_human.GeneID.NCBIGeneSymbol2HGNC.subtable.rm_idx = which(PTHR13_human.GeneID.NCBIGeneSymbol2HGNC.subtable[,'hgnc_id']=='')
if(length(PTHR13_human.GeneID.NCBIGeneSymbol2HGNC.subtable.rm_idx) > 0)	PTHR13_human.GeneID.NCBIGeneSymbol2HGNC.subtable = PTHR13_human.GeneID.NCBIGeneSymbol2HGNC.subtable[-PTHR13_human.GeneID.NCBIGeneSymbol2HGNC.subtable.rm_idx,,drop=FALSE]
PTHR13_human.GeneID.NCBIGeneSymbol2HGNC.subtable = unique(PTHR13_human.GeneID.NCBIGeneSymbol2HGNC.subtable)
if(length(PTHR13_human.GeneID.NCBIGeneSymbol2HGNC.subtable[,'hgnc_symbol']) > 0)
{
	PTHR13_human.GeneID.NCBIGeneSymbol2HGNC.subtable[,'hgnc_symbol'] = paste('Gene=',PTHR13_human.GeneID.NCBIGeneSymbol2HGNC.subtable[,'hgnc_symbol'],sep="")
	
	for(i in 1:length(PTHR13_human.GeneID.NCBIGeneSymbol2HGNC.subtable[,'hgnc_symbol']))
	{
		PTHR13_human.GeneID.HGNCid[which(PTHR13_human.GeneID.HGNCid %in% PTHR13_human.GeneID.NCBIGeneSymbol2HGNC.subtable[i,'hgnc_symbol'])] = PTHR13_human.GeneID.NCBIGeneSymbol2HGNC.subtable[i,'hgnc_id']
	}
}

PTHR13_human.GeneID.HGNCid.NotHGNC.idx <- grep('HGNC:',PTHR13_human.GeneID.HGNCid,invert=TRUE,value=FALSE)
PTHR13_human.GeneID.HGNCid.NotHGNC.id <- grep('HGNC:',PTHR13_human.GeneID.HGNCid,invert=TRUE,value=TRUE)
#PTHR13_human.GeneID.HGNCid.NotHGNC.idx
#PTHR13_human.GeneID.HGNCid.NotHGNC.id
#length(PTHR13_human.GeneID.HGNCid.NotHGNC.idx)
#length(PTHR13_human.GeneID.HGNCid.NotHGNC.id)

# Note: not manually check yet
if(length(PTHR13_human.GeneID.HGNCid.NotHGNC.idx) > 0)
{
	for(i in 1:length(PTHR13_human.GeneID.HGNCid.NotHGNC.idx))
	{
		if(PTHR13_human.GeneID.HGNCid[PTHR13_human.GeneID.HGNCid.NotHGNC.idx[i]] == '') {
			PTHR13_human.GeneID.HGNCid[PTHR13_human.GeneID.HGNCid.NotHGNC.idx[i]] = ''
		} else {
			if(PTHR13_human.GeneID.HGNCid[PTHR13_human.GeneID.HGNCid.NotHGNC.idx[i]] == '') {
				#PTHR13_human.GeneID.HGNCid[PTHR13_human.GeneID.HGNCid.NotHGNC.idx[i]] = ''
			} else {
#				print(cat(paste('"',PTHR13_human.GeneID[PTHR13_human.GeneID.HGNCid.NotHGNC.idx[i]],'" does not have HGNC id in "./Database/PANTHER/sequence_classifications/13.1/PTHR13.1_human" file.\n',sep="")))
			}
		}
	}
}

PTHR13_human.GeneID.HGNCid.NotHGNC.idx <- grep('HGNC:',PTHR13_human.GeneID.HGNCid,invert=TRUE,value=FALSE)
PTHR13_human.GeneID.HGNCid.NotHGNC.id <- grep('HGNC:',PTHR13_human.GeneID.HGNCid,invert=TRUE,value=TRUE)
#PTHR13_human.GeneID.HGNCid.NotHGNC.idx
#PTHR13_human.GeneID.HGNCid.NotHGNC.id
#length(PTHR13_human.GeneID.HGNCid.NotHGNC.idx)
#length(PTHR13_human.GeneID.HGNCid.NotHGNC.id)

# +++ mouse +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
PTHR13_mouse.GeneID.species <- sapply(strsplit(PTHR13_mouse.GeneID,'\\|'), FUN=function(X) X[1])
PTHR13_mouse.GeneID.UniProtKB <- gsub('UniProtKB=','',sapply(strsplit(PTHR13_mouse.GeneID,'\\|'), FUN=function(X) X[3]))

PTHR13_mouse.GeneID.MGIid <- gsub('MGI=MGI=','MGI:',sapply(strsplit(PTHR13_mouse.GeneID,'\\|'), FUN=function(X) X[2]))
PTHR13_mouse.GeneID.MGIid.NotMGI.idx <- grep('MGI:',PTHR13_mouse.GeneID.MGIid,invert=TRUE,value=FALSE)
PTHR13_mouse.GeneID.MGIid.NotMGI.id <- grep('MGI:',PTHR13_mouse.GeneID.MGIid,invert=TRUE,value=TRUE)

PTHR13_mouse.GeneID.EnsemblGene.idx <- grep('Ensembl=',PTHR13_mouse.GeneID.MGIid,invert=FALSE,value=FALSE)
PTHR13_mouse.GeneID.EnsemblGene.id <- grep('Ensembl=',PTHR13_mouse.GeneID.MGIid,invert=FALSE,value=TRUE)
PTHR13_mouse.GeneID.EnsemblGene <- gsub('Ensembl=','',PTHR13_mouse.GeneID.EnsemblGene.id)

PTHR13_mouse.GeneID.NCBIGeneID.idx <- grep('GeneID=',PTHR13_mouse.GeneID.MGIid,invert=FALSE,value=FALSE)
PTHR13_mouse.GeneID.NCBIGeneID.id <- grep('GeneID=',PTHR13_mouse.GeneID.MGIid,invert=FALSE,value=TRUE)
PTHR13_mouse.GeneID.NCBIGeneID <- gsub('GeneID=','',PTHR13_mouse.GeneID.NCBIGeneID.id)

PTHR13_mouse.GeneID.NCBIGeneSymbol.idx <- grep('Gene=',PTHR13_mouse.GeneID.MGIid,invert=FALSE,value=FALSE)
PTHR13_mouse.GeneID.NCBIGeneSymbol.id <- grep('Gene=',PTHR13_mouse.GeneID.MGIid,invert=FALSE,value=TRUE)
PTHR13_mouse.GeneID.NCBIGeneSymbol <- gsub('Gene=','',PTHR13_mouse.GeneID.NCBIGeneSymbol.id)


if(length(PTHR13_mouse.GeneID.EnsemblGene) + length(PTHR13_mouse.GeneID.NCBIGeneID) + length(PTHR13_mouse.GeneID.NCBIGeneSymbol) != length(PTHR13_mouse.GeneID.MGIid.NotMGI.id))
{
	cat('Please check Gene Identifier (the first column) in \"',paste('./Database/PANTHER/sequence_classifications/13.1/PTHR13.1_mouse',sep=""),'\" file.','\n','Non-MGI id is not only \'Ensembl=\', \'GeneID=\', \'Gene=\'.','\n',sep="")
}

ensembl = ensembl.mmu
#listFilters(ensembl)[,c('name','description')]
#listAttributes(ensembl)[,c('name','description')]

a = PTHR13_mouse.GeneID.EnsemblGene
b = ensembl.mmu
PTHR13_mouse.GeneID.EnsemblGene.attributes <- getBM(attributes=c('ensembl_gene_id','entrezgene','uniprotswissprot'), filters='ensembl_gene_id', values=a, mart=b)
PTHR13_mouse.GeneID.EnsemblGene.attributes = unique(PTHR13_mouse.GeneID.EnsemblGene.attributes)
PTHR13_mouse.GeneID.EnsemblGene2MGI <- getBM(attributes=c('ensembl_gene_id','mgi_id','mgi_symbol'), filters='ensembl_gene_id', values=a, mart=b)
PTHR13_mouse.GeneID.EnsemblGene2MGI.table = merge(PTHR13_mouse.GeneID.EnsemblGene2MGI, PTHR13_mouse.GeneID.EnsemblGene.attributes, by.x="ensembl_gene_id", by.y="ensembl_gene_id")
if(length(which(is.na(PTHR13_mouse.GeneID.EnsemblGene2MGI.table[,'entrezgene']))) > 0)		PTHR13_mouse.GeneID.EnsemblGene2MGI.table[which(is.na(PTHR13_mouse.GeneID.EnsemblGene2MGI.table[,'entrezgene'])),'entrezgene'] = ''
PTHR13_mouse.GeneID.EnsemblGene2MGI.table.rm_idx = which(PTHR13_mouse.GeneID.EnsemblGene2MGI.table[,'mgi_id']=='' & PTHR13_mouse.GeneID.EnsemblGene2MGI.table[,'mgi_symbol']=='' & PTHR13_mouse.GeneID.EnsemblGene2MGI.table[,'entrezgene']=='' & PTHR13_mouse.GeneID.EnsemblGene2MGI.table[,'uniprotswissprot']=='')
if(length(PTHR13_mouse.GeneID.EnsemblGene2MGI.table.rm_idx) > 0)		PTHR13_mouse.GeneID.EnsemblGene2MGI.table = PTHR13_mouse.GeneID.EnsemblGene2MGI.table[-PTHR13_mouse.GeneID.EnsemblGene2MGI.table.rm_idx,,drop=FALSE]
PTHR13_mouse.GeneID.EnsemblGene2MGI.table = unique(PTHR13_mouse.GeneID.EnsemblGene2MGI.table)
PTHR13_mouse.GeneID.EnsemblGene2MGI.subtable = unique(PTHR13_mouse.GeneID.EnsemblGene2MGI.table[,c('ensembl_gene_id','mgi_id')])
PTHR13_mouse.GeneID.EnsemblGene2MGI.subtable.rm_idx = which(PTHR13_mouse.GeneID.EnsemblGene2MGI.subtable[,'mgi_id']=='')
if(length(PTHR13_mouse.GeneID.EnsemblGene2MGI.subtable.rm_idx) > 0)	PTHR13_mouse.GeneID.EnsemblGene2MGI.subtable = PTHR13_mouse.GeneID.EnsemblGene2MGI.subtable[-PTHR13_mouse.GeneID.EnsemblGene2MGI.subtable.rm_idx,,drop=FALSE]
PTHR13_mouse.GeneID.EnsemblGene2MGI.subtable = unique(PTHR13_mouse.GeneID.EnsemblGene2MGI.subtable)
if(length(PTHR13_mouse.GeneID.EnsemblGene2MGI.subtable[,'ensembl_gene_id']) > 0)
{
	PTHR13_mouse.GeneID.EnsemblGene2MGI.subtable[,'ensembl_gene_id'] = paste('Ensembl=',PTHR13_mouse.GeneID.EnsemblGene2MGI.subtable[,'ensembl_gene_id'],sep="")
	
	for(i in 1:length(PTHR13_mouse.GeneID.EnsemblGene2MGI.subtable[,'ensembl_gene_id']))
	{
		PTHR13_mouse.GeneID.MGIid[which(PTHR13_mouse.GeneID.MGIid %in% PTHR13_mouse.GeneID.EnsemblGene2MGI.subtable[i,'ensembl_gene_id'])] = PTHR13_mouse.GeneID.EnsemblGene2MGI.subtable[i,'mgi_id']
	}
}

a = PTHR13_mouse.GeneID.NCBIGeneID
b = ensembl.mmu
PTHR13_mouse.GeneID.NCBIGeneID.attributes <- getBM(attributes=c('entrezgene','ensembl_gene_id','uniprotswissprot'), filters='entrezgene', values=a, mart=b)
PTHR13_mouse.GeneID.NCBIGeneID.attributes = unique(PTHR13_mouse.GeneID.NCBIGeneID.attributes)
PTHR13_mouse.GeneID.NCBIGeneID2MGI <- getBM(attributes=c('entrezgene','mgi_id','mgi_symbol'), filters='entrezgene', values=a, mart=b)
PTHR13_mouse.GeneID.NCBIGeneID2MGI.table = merge(PTHR13_mouse.GeneID.NCBIGeneID2MGI, PTHR13_mouse.GeneID.NCBIGeneID.attributes, by.x="entrezgene", by.y="entrezgene")
if(length(which(is.na(PTHR13_mouse.GeneID.NCBIGeneID2MGI.table[,'entrezgene']))) > 0)		PTHR13_mouse.GeneID.NCBIGeneID2MGI.table[which(is.na(PTHR13_mouse.GeneID.NCBIGeneID2MGI.table[,'entrezgene'])),'entrezgene'] = ''
PTHR13_mouse.GeneID.NCBIGeneID2MGI.table.rm_idx = which(PTHR13_mouse.GeneID.NCBIGeneID2MGI.table[,'mgi_id']=='' & PTHR13_mouse.GeneID.NCBIGeneID2MGI.table[,'mgi_symbol']=='' & PTHR13_mouse.GeneID.NCBIGeneID2MGI.table[,'ensembl_gene_id']=='' & PTHR13_mouse.GeneID.NCBIGeneID2MGI.table[,'uniprotswissprot']=='')
if(length(PTHR13_mouse.GeneID.NCBIGeneID2MGI.table.rm_idx) > 0)		PTHR13_mouse.GeneID.NCBIGeneID2MGI.table = PTHR13_mouse.GeneID.NCBIGeneID2MGI.table[-PTHR13_mouse.GeneID.NCBIGeneID2MGI.table.rm_idx,,drop=FALSE]
PTHR13_mouse.GeneID.NCBIGeneID2MGI.table = unique(PTHR13_mouse.GeneID.NCBIGeneID2MGI.table)
PTHR13_mouse.GeneID.NCBIGeneID2MGI.subtable = unique(PTHR13_mouse.GeneID.NCBIGeneID2MGI.table[,c('entrezgene','mgi_id')])
PTHR13_mouse.GeneID.NCBIGeneID2MGI.subtable.rm_idx = which(PTHR13_mouse.GeneID.NCBIGeneID2MGI.subtable[,'mgi_id']=='')
if(length(PTHR13_mouse.GeneID.NCBIGeneID2MGI.subtable.rm_idx) > 0)	PTHR13_mouse.GeneID.NCBIGeneID2MGI.subtable = PTHR13_mouse.GeneID.NCBIGeneID2MGI.subtable[-PTHR13_mouse.GeneID.NCBIGeneID2MGI.subtable.rm_idx,,drop=FALSE]
PTHR13_mouse.GeneID.NCBIGeneID2MGI.subtable = unique(PTHR13_mouse.GeneID.NCBIGeneID2MGI.subtable)
if(length(PTHR13_mouse.GeneID.NCBIGeneID2MGI.subtable[,'entrezgene']) > 0)
{
	PTHR13_mouse.GeneID.NCBIGeneID2MGI.subtable[,'entrezgene'] = paste('GeneID=',PTHR13_mouse.GeneID.NCBIGeneID2MGI.subtable[,'entrezgene'],sep="")
	
	for(i in 1:length(PTHR13_mouse.GeneID.NCBIGeneID2MGI.subtable[,'entrezgene']))
	{
		PTHR13_mouse.GeneID.MGIid[which(PTHR13_mouse.GeneID.MGIid %in% PTHR13_mouse.GeneID.NCBIGeneID2MGI.subtable[i,'entrezgene'])] = PTHR13_mouse.GeneID.NCBIGeneID2MGI.subtable[i,'mgi_id']
	}
}

a = PTHR13_mouse.GeneID.NCBIGeneSymbol
b = ensembl.mmu
PTHR13_mouse.GeneID.NCBIGeneSymbol.attributes <- getBM(attributes=c('mgi_symbol','ensembl_gene_id','uniprotswissprot'), filters='mgi_symbol', values=a, mart=b)
PTHR13_mouse.GeneID.NCBIGeneSymbol.attributes = unique(PTHR13_mouse.GeneID.NCBIGeneSymbol.attributes)
PTHR13_mouse.GeneID.NCBIGeneSymbol2MGI <- getBM(attributes=c('mgi_symbol','mgi_id','entrezgene'), filters='mgi_symbol', values=a, mart=b)
PTHR13_mouse.GeneID.NCBIGeneSymbol2MGI.table = merge(PTHR13_mouse.GeneID.NCBIGeneSymbol2MGI, PTHR13_mouse.GeneID.NCBIGeneSymbol.attributes, by.x="mgi_symbol", by.y="mgi_symbol")
if(length(which(is.na(PTHR13_mouse.GeneID.NCBIGeneSymbol2MGI.table[,'entrezgene']))) > 0)	PTHR13_mouse.GeneID.NCBIGeneSymbol2MGI.table[which(is.na(PTHR13_mouse.GeneID.NCBIGeneSymbol2MGI.table[,'entrezgene'])),'entrezgene'] = ''
PTHR13_mouse.GeneID.NCBIGeneSymbol2MGI.table.rm_idx = which(PTHR13_mouse.GeneID.NCBIGeneSymbol2MGI.table[,'mgi_id']=='' & PTHR13_mouse.GeneID.NCBIGeneSymbol2MGI.table[,'entrezgene']=='' & PTHR13_mouse.GeneID.NCBIGeneSymbol2MGI.table[,'ensembl_gene_id']=='' & PTHR13_mouse.GeneID.NCBIGeneSymbol2MGI.table[,'uniprotswissprot']=='')
if(length(PTHR13_mouse.GeneID.NCBIGeneSymbol2MGI.table.rm_idx) > 0)		PTHR13_mouse.GeneID.NCBIGeneSymbol2MGI.table = PTHR13_mouse.GeneID.NCBIGeneSymbol2MGI.table[-PTHR13_mouse.GeneID.NCBIGeneSymbol2MGI.table.rm_idx,,drop=FALSE]
PTHR13_mouse.GeneID.NCBIGeneSymbol2MGI.table = unique(PTHR13_mouse.GeneID.NCBIGeneSymbol2MGI.table)
PTHR13_mouse.GeneID.NCBIGeneSymbol2MGI.subtable = unique(PTHR13_mouse.GeneID.NCBIGeneSymbol2MGI.table[,c('mgi_symbol','mgi_id')])
PTHR13_mouse.GeneID.NCBIGeneSymbol2MGI.subtable.rm_idx = which(PTHR13_mouse.GeneID.NCBIGeneSymbol2MGI.subtable[,'mgi_id']=='')
if(length(PTHR13_mouse.GeneID.NCBIGeneSymbol2MGI.subtable.rm_idx) > 0)	PTHR13_mouse.GeneID.NCBIGeneSymbol2MGI.subtable = PTHR13_mouse.GeneID.NCBIGeneSymbol2MGI.subtable[-PTHR13_mouse.GeneID.NCBIGeneSymbol2MGI.subtable.rm_idx,,drop=FALSE]
PTHR13_mouse.GeneID.NCBIGeneSymbol2MGI.subtable = unique(PTHR13_mouse.GeneID.NCBIGeneSymbol2MGI.subtable)
if(length(PTHR13_mouse.GeneID.NCBIGeneSymbol2MGI.subtable[,'mgi_symbol']) > 0)
{
	PTHR13_mouse.GeneID.NCBIGeneSymbol2MGI.subtable[,'mgi_symbol'] = paste('Gene=',PTHR13_mouse.GeneID.NCBIGeneSymbol2MGI.subtable[,'mgi_symbol'],sep="")
	
	for(i in 1:length(PTHR13_mouse.GeneID.NCBIGeneSymbol2MGI.subtable[,'mgi_symbol']))
	{
		PTHR13_mouse.GeneID.MGIid[which(PTHR13_mouse.GeneID.MGIid %in% PTHR13_mouse.GeneID.NCBIGeneSymbol2MGI.subtable[i,'mgi_symbol'])] = PTHR13_mouse.GeneID.NCBIGeneSymbol2MGI.subtable[i,'mgi_id']
	}
}

PTHR13_mouse.GeneID.MGIid.NotMGI.idx <- grep('MGI:',PTHR13_mouse.GeneID.MGIid,invert=TRUE,value=FALSE)
PTHR13_mouse.GeneID.MGIid.NotMGI.id <- grep('MGI:',PTHR13_mouse.GeneID.MGIid,invert=TRUE,value=TRUE)
#PTHR13_mouse.GeneID.MGIid.NotMGI.idx
#PTHR13_mouse.GeneID.MGIid.NotMGI.id
#length(PTHR13_mouse.GeneID.MGIid.NotMGI.idx)
#length(PTHR13_mouse.GeneID.MGIid.NotMGI.id)

# Note: not manually check for PTHR13.1_mouse yet (but already check for PTHR11.0_mouse)
if(length(PTHR13_mouse.GeneID.MGIid.NotMGI.idx) > 0)
{
	for(i in 1:length(PTHR13_mouse.GeneID.MGIid.NotMGI.idx))
	{
		if(PTHR13_mouse.GeneID.MGIid[PTHR13_mouse.GeneID.MGIid.NotMGI.idx[i]] == 'GeneID=102637099') {						# https://www.ncbi.nlm.nih.gov/gene/102637099
			PTHR13_mouse.GeneID.MGIid[PTHR13_mouse.GeneID.MGIid.NotMGI.idx[i]] = 'MGI:5621380'								# http://www.informatics.jax.org/marker/MGI:5621380
		} else if(PTHR13_mouse.GeneID.MGIid[PTHR13_mouse.GeneID.MGIid.NotMGI.idx[i]] == 'Gene=Mup8') {						# https://www.ncbi.nlm.nih.gov/gene/100041687
			PTHR13_mouse.GeneID.MGIid[PTHR13_mouse.GeneID.MGIid.NotMGI.idx[i]] = 'MGI:3709619'								# http://www.informatics.jax.org/marker/MGI:3709619
		} else if(PTHR13_mouse.GeneID.MGIid[PTHR13_mouse.GeneID.MGIid.NotMGI.idx[i]] == 'Gene=Mup11') {						# https://www.ncbi.nlm.nih.gov/gene/100039028
			PTHR13_mouse.GeneID.MGIid[PTHR13_mouse.GeneID.MGIid.NotMGI.idx[i]] = 'MGI:3709617'								# http://www.informatics.jax.org/marker/MGI:3709617
		} else if(PTHR13_mouse.GeneID.MGIid[PTHR13_mouse.GeneID.MGIid.NotMGI.idx[i]] == 'Gene=V1ra11') {					# https://www.ncbi.nlm.nih.gov/gene/?term=V1ra11 -> Vmn1r45 (Also known as: V1r2; V1RA9; V1ra2; mV1R2; V1ra11) (http://www.uniprot.org/uniprot/Q8R2E6)
			PTHR13_mouse.GeneID.MGIid[PTHR13_mouse.GeneID.MGIid.NotMGI.idx[i]] = 'MGI:1333762'								# http://www.informatics.jax.org/marker/MGI:1333762
		} else if(PTHR13_mouse.GeneID.MGIid[PTHR13_mouse.GeneID.MGIid.NotMGI.idx[i]] == 'Ensembl=ENSMUSG00000076940') {		# http://asia.ensembl.org/Mus_musculus/Gene/Summary?db=core;g=ENSMUSG00000076940;r=16:19260401-19260859
			PTHR13_mouse.GeneID.MGIid[PTHR13_mouse.GeneID.MGIid.NotMGI.idx[i]] = 'MGI:99548'								# http://www.informatics.jax.org/marker/MGI:99548
		} else if(PTHR13_mouse.GeneID.MGIid[PTHR13_mouse.GeneID.MGIid.NotMGI.idx[i]] == 'Ensembl=ENSMUSG00000078673') {		# http://asia.ensembl.org/Mus_musculus/Gene/Summary?db=core;g=ENSMUSG00000078673;r=4:61778324-61782269
			PTHR13_mouse.GeneID.MGIid[PTHR13_mouse.GeneID.MGIid.NotMGI.idx[i]] = 'MGI:3705235'								# http://www.informatics.jax.org/marker/MGI:3705235
		} else if(PTHR13_mouse.GeneID.MGIid[PTHR13_mouse.GeneID.MGIid.NotMGI.idx[i]] == 'Ensembl=ENSMUSG00000078675') {		# http://asia.ensembl.org/Mus_musculus/Gene/Summary?db=core;g=ENSMUSG00000078675;r=4:61515592-61519531
			PTHR13_mouse.GeneID.MGIid[PTHR13_mouse.GeneID.MGIid.NotMGI.idx[i]] = 'MGI:3780250'								# http://www.informatics.jax.org/marker/MGI:3780250
		} else if(PTHR13_mouse.GeneID.MGIid[PTHR13_mouse.GeneID.MGIid.NotMGI.idx[i]] == 'Ensembl=ENSMUSG00000078686') {		# http://asia.ensembl.org/Mus_musculus/Gene/Summary?db=core;g=ENSMUSG00000078686;r=4:60418046-60421952
			PTHR13_mouse.GeneID.MGIid[PTHR13_mouse.GeneID.MGIid.NotMGI.idx[i]] = 'MGI:3782918'								# http://www.informatics.jax.org/marker/MGI:3782918
		} else if(PTHR13_mouse.GeneID.MGIid[PTHR13_mouse.GeneID.MGIid.NotMGI.idx[i]] == 'Ensembl=ENSMUSG00000102805') {		# http://asia.ensembl.org/Mus_musculus/Gene/Summary?db=core;g=ENSMUSG00000102805;r=3:84496426-85887518
			PTHR13_mouse.GeneID.MGIid[PTHR13_mouse.GeneID.MGIid.NotMGI.idx[i]] = 'MGI:5610468'								# http://www.informatics.jax.org/marker/MGI:5610468
		} else if(PTHR13_mouse.GeneID.MGIid[PTHR13_mouse.GeneID.MGIid.NotMGI.idx[i]] == 'Ensembl=ENSMUSG00000104824') {		# http://asia.ensembl.org/Mus_musculus/Gene/Summary?db=core;g=ENSMUSG00000104824;r=5:25935791-25942381;t=ENSMUST00000196214
			PTHR13_mouse.GeneID.MGIid[PTHR13_mouse.GeneID.MGIid.NotMGI.idx[i]] = 'MGI:5435018'								# http://www.informatics.jax.org/marker/MGI:5435018
		} else if(PTHR13_mouse.GeneID.MGIid[PTHR13_mouse.GeneID.MGIid.NotMGI.idx[i]] == 'Ensembl=ENSMUSG00000105204') {		# http://asia.ensembl.org/Mus_musculus/Gene/Summary?db=core;g=ENSMUSG00000105204;r=3:89081939-89101099;t=ENSMUST00000200659
			PTHR13_mouse.GeneID.MGIid[PTHR13_mouse.GeneID.MGIid.NotMGI.idx[i]] = 'MGI:5663875'								# http://www.informatics.jax.org/marker/MGI:5663875
		} else if(PTHR13_mouse.GeneID.MGIid[PTHR13_mouse.GeneID.MGIid.NotMGI.idx[i]] == 'Ensembl=ENSMUSG00000105630') {		# http://asia.ensembl.org/Mus_musculus/Gene/Summary?db=core;g=ENSMUSG00000105630;r=6:68880178-68880674;t=ENSMUST00000198756
			PTHR13_mouse.GeneID.MGIid[PTHR13_mouse.GeneID.MGIid.NotMGI.idx[i]] = 'MGI:5662680'								# http://www.informatics.jax.org/marker/MGI:5662680
		} else if(PTHR13_mouse.GeneID.MGIid[PTHR13_mouse.GeneID.MGIid.NotMGI.idx[i]] == 'Ensembl=ENSMUSG00000105835') {		# http://asia.ensembl.org/Mus_musculus/Gene/Summary?db=core;g=ENSMUSG00000105835;r=5:65427427-65492833;t=ENSMUST00000196121
			PTHR13_mouse.GeneID.MGIid[PTHR13_mouse.GeneID.MGIid.NotMGI.idx[i]] = 'MGI:5663689'								# http://www.informatics.jax.org/marker/MGI:5663689
		} else if(PTHR13_mouse.GeneID.MGIid[PTHR13_mouse.GeneID.MGIid.NotMGI.idx[i]] == 'Ensembl=ENSMUSG00000105993') {		# http://asia.ensembl.org/Mus_musculus/Gene/Summary?db=core;g=ENSMUSG00000105993;r=X:42489691-42491524;t=ENSMUST00000197237
			PTHR13_mouse.GeneID.MGIid[PTHR13_mouse.GeneID.MGIid.NotMGI.idx[i]] = 'MGI:5663474'								# http://www.informatics.jax.org/marker/MGI:5663474
		} else if(PTHR13_mouse.GeneID.MGIid[PTHR13_mouse.GeneID.MGIid.NotMGI.idx[i]] == 'Ensembl=ENSMUSG00000106445') {		# http://asia.ensembl.org/Mus_musculus/Gene/Summary?db=core;g=ENSMUSG00000106445;r=5:15524456-15529222;t=ENSMUST00000196384
			PTHR13_mouse.GeneID.MGIid[PTHR13_mouse.GeneID.MGIid.NotMGI.idx[i]] = 'MGI:5434545'								# http://www.informatics.jax.org/marker/MGI:5434545
		} else if(PTHR13_mouse.GeneID.MGIid[PTHR13_mouse.GeneID.MGIid.NotMGI.idx[i]] == 'Ensembl=ENSMUSG00000106627') {		# http://asia.ensembl.org/Mus_musculus/Gene/Summary?db=core;g=ENSMUSG00000106627;r=5:25965961-25972663;t=ENSMUST00000200447
			PTHR13_mouse.GeneID.MGIid[PTHR13_mouse.GeneID.MGIid.NotMGI.idx[i]] = 'MGI:5435035'								# http://www.informatics.jax.org/marker/MGI:5435035
		} else {
			if(PTHR13_mouse.GeneID.MGIid[PTHR13_mouse.GeneID.MGIid.NotMGI.idx[i]] == 'Gene=Clul1') {						# http://www.uniprot.org/uniprot/Q3ZRW6
				#PTHR13_mouse.GeneID.MGIid[PTHR13_mouse.GeneID.MGIid.NotMGI.idx[i]] = ''									# http://www.informatics.jax.org/sequence/key/34934360
			} else if(PTHR13_mouse.GeneID.MGIid[PTHR13_mouse.GeneID.MGIid.NotMGI.idx[i]] == 'Gene=Pol') {					# http://www.uniprot.org/uniprot/P10400
				#PTHR13_mouse.GeneID.MGIid[PTHR13_mouse.GeneID.MGIid.NotMGI.idx[i]] = ''									# http://www.informatics.jax.org/sequence/key/34942720
			} else if(PTHR13_mouse.GeneID.MGIid[PTHR13_mouse.GeneID.MGIid.NotMGI.idx[i]] == 'Gene=Smoktcr') {				# http://www.uniprot.org/uniprot/A2KF29
				#PTHR13_mouse.GeneID.MGIid[PTHR13_mouse.GeneID.MGIid.NotMGI.idx[i]] = ''									# http://www.informatics.jax.org/sequence/key/34945114
			} else if(PTHR13_mouse.GeneID.MGIid[PTHR13_mouse.GeneID.MGIid.NotMGI.idx[i]] == 'Ensembl=ENSMUSG00000094329') {	# http://asia.ensembl.org/Mus_musculus/Gene/Idhistory?g=ENSMUSG00000094329 (http://www.uniprot.org/uniprot/A0A075B673; http://www.uniprot.org/uniprot/A0A075B673?version=*)
				#PTHR13_mouse.GeneID.MGIid[PTHR13_mouse.GeneID.MGIid.NotMGI.idx[i]] = ''
			} else if(PTHR13_mouse.GeneID.MGIid[PTHR13_mouse.GeneID.MGIid.NotMGI.idx[i]] == 'Ensembl=ENSMUSG00000096877') {	# http://asia.ensembl.org/Mus_musculus/Gene/Idhistory?g=ENSMUSG00000096877 (http://www.uniprot.org/uniprot/A0A075B673; http://www.uniprot.org/uniprot/A0A075B673?version=*)
				#PTHR13_mouse.GeneID.MGIid[PTHR13_mouse.GeneID.MGIid.NotMGI.idx[i]] = ''
			} else {
#				print(cat(paste('"',PTHR13_mouse.GeneID[PTHR13_mouse.GeneID.MGIid.NotMGI.idx[i]],'" does not have MGI id in "./Database/PANTHER/sequence_classifications/13.1/PTHR13.1_mouse" file.\n',sep="")))
			}
		}
	}
}

PTHR13_mouse.GeneID.MGIid.NotMGI.idx <- grep('MGI:',PTHR13_mouse.GeneID.MGIid,invert=TRUE,value=FALSE)
PTHR13_mouse.GeneID.MGIid.NotMGI.id <- grep('MGI:',PTHR13_mouse.GeneID.MGIid,invert=TRUE,value=TRUE)
#PTHR13_mouse.GeneID.MGIid.NotMGI.idx
#PTHR13_mouse.GeneID.MGIid.NotMGI.id
#length(PTHR13_mouse.GeneID.MGIid.NotMGI.idx)
#length(PTHR13_mouse.GeneID.MGIid.NotMGI.id)


# === subset ====================================================================================================================
folder <- c('high','low')

for(r in 1:length(batch))
{
	if (file.exists(paste('./Data/',Sample,'/03_homologene/subset/n_',r,sep=""))) {
		
		mainDir <- paste('./Data/',Sample,'/08_UniProt2ProteinClass/subset',sep="")
		subDir <- paste('n_',r,sep="")
		if (file.exists(file.path(mainDir, subDir))){
		} else {
			dir.create(file.path(mainDir, subDir))
		}
		
		for(f in 1:length(folder))
		{
			if (file.exists(paste('./Data/',Sample,'/03_homologene/subset/n_',r,'/',folder[f],sep=""))) {
				
				if (file.exists(paste('./Data/',Sample,'/03_homologene/subset/n_',r,'/',folder[f],'/',Sample,'_exo_UniProt_homologene.xlsx',sep=""))) {
					
					mainDir <- paste('./Data/',Sample,'/08_UniProt2ProteinClass/subset/n_',r,sep="")
					subDir <- folder[f]
					if (file.exists(file.path(mainDir, subDir))){
					} else {
						dir.create(file.path(mainDir, subDir))
					}
					
					
					# === input GeneID ==============================================================================================================
					#homologene <- read.csv(paste('./Data/',Sample,'/03_homologene/subset/n_',r,'/',folder[f],'/',Sample,'_exo_UniProt_homologene.csv',sep=""), header=TRUE, stringsAsFactors=FALSE)
					homologene <- read.xlsx(paste('./Data/',Sample,'/03_homologene/subset/n_',r,'/',folder[f],'/',Sample,'_exo_UniProt_homologene.xlsx',sep=""), sheet='homologene (human2mouse)', colNames=TRUE, rowNames=FALSE, check.names=FALSE)
					
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
					
					
					# === search for ProteinClass of MS-identified mouse proteins and their human homolog proteins ==================================
					for(sp in 1:length(SpeciesName.general))
					{
						if(SpeciesName.general[sp] != "rat")
						{
							mainDir <- paste('./Data/',Sample,'/08_UniProt2ProteinClass/subset/n_',r,'/',folder[f],sep="")
							subDir <- SpeciesName.general[sp]
							if (file.exists(file.path(mainDir, subDir))){
							} else {
								dir.create(file.path(mainDir, subDir))
							}
							
							
							if(SpeciesName.general[sp] == "human") {
								SpeciesDB = 'HGNC'
								PTHR13_species.GeneID.UniProtKB = PTHR13_human.GeneID.UniProtKB
								PTHR13_species.GeneID.SpeciesDB = PTHR13_human.GeneID.HGNCid
								PTHR13_species.ProteinClass = PTHR13_human.ProteinClass
								homologene.species.UniProt = homologene.human.UniProt
								homologene.species.SpeciesDB_ID = homologene.human.HGNC.ID
								homologene.species.SpeciesDB_Symbol = homologene.human.HGNC.Symbol
								homologene.species.Gene_ID = homologene.human.Gene.ID
								homologene.species.Protein_names = homologene.human.Protein.names
							} else if(SpeciesName.general[sp] == "mouse") {
								SpeciesDB = 'MGI'
								PTHR13_species.GeneID.UniProtKB = PTHR13_mouse.GeneID.UniProtKB
								PTHR13_species.GeneID.SpeciesDB = PTHR13_mouse.GeneID.MGIid
								PTHR13_species.ProteinClass = PTHR13_mouse.ProteinClass
								homologene.species.UniProt = homologene.mouse.UniProt
								homologene.species.SpeciesDB_ID = homologene.mouse.MGI.ID
								homologene.species.SpeciesDB_Symbol = homologene.mouse.MGI.Symbol
								homologene.species.Gene_ID = homologene.mouse.Gene.ID
								homologene.species.Protein_names = homologene.mouse.Protein.names
							}
							
							
							mat_colnames = c(paste(SpeciesName.general[sp],'.UniProt',sep=""), paste(SpeciesName.general[sp],'.',SpeciesDB,'.ID',sep=""), paste(SpeciesName.general[sp],'.Gene.ID',sep=""), paste(SpeciesName.general[sp],'.',SpeciesDB,'.Symbol',sep=""), paste(SpeciesName.general[sp],'.Protein.names',sep=""))
							
							
							# *** obsolete entry -> new entry (start) *****************************************************
							homologene.species.UniProt.New_UniProt = homologene.species.UniProt
							
							# obsolete entry -> new entry: UniProtKB - P62158 (CALM_HUMAN) -> On May 10, 2017 this entry became obsolete. It can now be found as secondary accession in P0DP23, P0DP24 and P0DP25.
							homologene.species.UniProt.New_UniProt.idx = which(homologene.species.UniProt.New_UniProt == 'P0DP23' | homologene.species.UniProt.New_UniProt == 'P0DP24' | homologene.species.UniProt.New_UniProt == 'P0DP25')
							if(length(homologene.species.UniProt.New_UniProt.idx) > 0)		homologene.species.UniProt.New_UniProt[homologene.species.UniProt.New_UniProt.idx] = 'P62158'
							
							# obsolete entry -> new entry: UniProtKB - P62204 (CALM_MOUSE) -> On May 10, 2017 this entry became obsolete. It can now be found as secondary accession in P0DP26, P0DP27 and P0DP28.
							homologene.species.UniProt.New_UniProt.idx = which(homologene.species.UniProt.New_UniProt == 'P0DP26' | homologene.species.UniProt.New_UniProt == 'P0DP27' | homologene.species.UniProt.New_UniProt == 'P0DP28')
							if(length(homologene.species.UniProt.New_UniProt.idx) > 0)		homologene.species.UniProt.New_UniProt[homologene.species.UniProt.New_UniProt.idx] = 'P62204'
							
							# *** obsolete entry -> new entry (end) *******************************************************
							
							
							PTHR13_species.GeneID.UniProtKB.homologene_species_UniProt.idx <- sapply(homologene.species.UniProt.New_UniProt, FUN=function(X) which(PTHR13_species.GeneID.UniProtKB %in% unique(unlist(strsplit(X,'\\|')))))
							PTHR13_species.GeneID.SpeciesDB.homologene_species_SpeciesDB_ID.idx <- sapply(homologene.species.SpeciesDB_ID, FUN=function(X) which(PTHR13_species.GeneID.SpeciesDB %in% unique(unlist(strsplit(X,'\\|')))))
							
							species.UniProt.ProteinClass_id.mat <- matrix('',length(homologene.species.UniProt),5+6,byrow=TRUE)
							species.UniProt.ProteinClass_level.mat <- matrix('',length(homologene.species.UniProt),5+6,byrow=TRUE)
							species.UniProt.ProteinClass_name.mat <- matrix('',length(homologene.species.UniProt),5+6,byrow=TRUE)
							species.UniProt.ProteinClass_description.mat <- matrix('',length(homologene.species.UniProt),5+6,byrow=TRUE)
							colnames(species.UniProt.ProteinClass_id.mat) = c(mat_colnames,'Level.0','Level.1','Level.2','Level.3','Level.4','Level.5')
							colnames(species.UniProt.ProteinClass_level.mat) = c(mat_colnames,'Level.0','Level.1','Level.2','Level.3','Level.4','Level.5')
							colnames(species.UniProt.ProteinClass_name.mat) = c(mat_colnames,'Level.0','Level.1','Level.2','Level.3','Level.4','Level.5')
							colnames(species.UniProt.ProteinClass_description.mat) = c(mat_colnames,'Level.0','Level.1','Level.2','Level.3','Level.4','Level.5')
							
							homologene.species.UniProt.NoProteinClass <- c()
							homologene.species.UniProt.NA <- c()
							
							Level.submatrix.all.new.all <- c()
							
							for(i in 1:length(homologene.species.UniProt))
							{
								# *** manually check for the missed ProteinClass (start) **************************************
								addinProteinClass <- c()
								
								if(homologene.species.SpeciesDB_Symbol[i] %in% c(c('Anxa1','Anxa2','Anxa3','Anxa4','Anxa5','Anxa6','Anxa7','Anxa8'),toupper(c('Anxa1','Anxa2','Anxa3','Anxa4','Anxa5','Anxa6','Anxa7','Anxa8')))) {		# Protein family of annexins
									addinProteinClass = c('annexin')
								} else if(homologene.species.SpeciesDB_Symbol[i] %in% c(c('Itga3','Itga4','Itga5','Itga6','Itgav','Itgb1'),toupper(c('Itga3','Itga4','Itga5','Itga6','Itgav','Itgb1')))) {								# Protein family of integrins
									addinProteinClass = c('integrin')
								} else if(homologene.species.SpeciesDB_Symbol[i] %in% c(c('Rps6','Rps8','Rps16','Rps18','Rps25','Rpsa','Rpl4','Rpl14','Rpl17','Rpl18','Rpl18a'),toupper(c('Rps6','Rps8','Rps16','Rps18','Rps25','Rpsa','Rpl4','Rpl14','Rpl17','Rpl18','Rpl18a')))) {	# Protein family of ribosomal proteins
									addinProteinClass = c('ribosomal protein')
								} else if(homologene.species.SpeciesDB_Symbol[i] %in% c('Hnrnpa3',toupper('Hnrnpa3'))) {	# 229279|MGI:1917171|Q8BG05|Heterogeneous nuclear ribonucleoprotein A3 (hnRNP A3)
									addinProteinClass = c('ribonucleoprotein')
								} else if(homologene.species.SpeciesDB_Symbol[i] %in% c('Cdh13',toupper('Cdh13'))) {		# 12554|MGI:99551|Q9WTR5|Cadherin-13 (Heart cadherin) (H-cadherin) (Truncated cadherin) (T-cad) (T-cadherin)
									addinProteinClass = c('cadherin')
								} else if(homologene.species.SpeciesDB_Symbol[i] %in% c('H2afz',toupper('H2afz'))) {		# 51788|MGI:1888388|P0C0S6|Histone H2A.Z (H2A/z)
									addinProteinClass = c('histone')
								} else if(homologene.species.SpeciesDB_Symbol[i] %in% c('Rbm3',toupper('Rbm3'))) {			# 19652|MGI:1099460|O89086|RNA-binding protein 3 (RNA-binding motif protein 3)
									addinProteinClass = c('RNA binding protein')
								} else if(homologene.species.SpeciesDB_Symbol[i] %in% c('Aldoa',toupper('Aldoa'))) {		# 11674|MGI:87994|P05064|Fructose-bisphosphate aldolase A (EC 4.1.2.13) (Aldolase 1) (Muscle-type aldolase)
									addinProteinClass = c('aldolase')
								} else if(homologene.species.SpeciesDB_Symbol[i] %in% c('Usp5',toupper('Usp5'))) {			# 22225|MGI:1347343|P56399|Ubiquitin carboxyl-terminal hydrolase 5 (EC 3.4.19.12) (Deubiquitinating enzyme 5) (Isopeptidase T) (Ubiquitin thioesterase 5) (Ubiquitin-specific-processing protease 5)
									addinProteinClass = c('hydrolase')
								} else if(homologene.species.SpeciesDB_Symbol[i] %in% c('Tpi1',toupper('Tpi1'))) {			# 21991|MGI:98797|P17751|Triosephosphate isomerase (TIM) (EC 5.3.1.1) (Triose-phosphate isomerase)
									addinProteinClass = c('isomerase')
								} else if(homologene.species.SpeciesDB_Symbol[i] %in% c('Nme1',toupper('Nme1'))) {			# 18102|MGI:97355|P15532|Nucleoside diphosphate kinase A (NDK A) (NDP kinase A) (EC 2.7.4.6) (Metastasis inhibition factor NM23) (NDPK-A) (Tumor metastatic process-associated protein) (nm23-M1)
									addinProteinClass = c('kinase')
								} else if(homologene.species.SpeciesDB_Symbol[i] %in% c('Pkm',toupper('Pkm'))) {			# 18746|MGI:97591|P52480|Pyruvate kinase PKM (EC 2.7.1.40) (Pyruvate kinase muscle isozyme)
									addinProteinClass = c('kinase')
								} else if(homologene.species.SpeciesDB_Symbol[i] %in% c('Pgam1',toupper('Pgam1'))) {		# 18648|MGI:97552|Q9DBJ1|Phosphoglycerate mutase 1 (EC 5.4.2.11) (EC 5.4.2.4) (BPG-dependent PGAM 1) (Phosphoglycerate mutase isozyme B) (PGAM-B)
									addinProteinClass = c('mutase')
								} else if(homologene.species.SpeciesDB_Symbol[i] %in% c('Gpx4',toupper('Gpx4'))) {			# 625249|MGI:104767|O70325|Phospholipid hydroperoxide glutathione peroxidase, mitochondrial (PHGPx) (EC 1.11.1.12) (Glutathione peroxidase 4) (GPx-4) (GSHPx-4)
									addinProteinClass = c('peroxidase')
								} else if(homologene.species.SpeciesDB_Symbol[i] %in% c('Smpdl3b',toupper('Smpdl3b'))) {	# 100340|MGI:1916022|P58242|Acid sphingomyelinase-like phosphodiesterase 3b (ASM-like phosphodiesterase 3b) (EC 3.1.4.-)
									addinProteinClass = c('phosphodiesterase')
								} else if(homologene.species.SpeciesDB_Symbol[i] %in% c('Pld2',toupper('Pld2'))) {			# 18806|MGI:892877|P97813|Phospholipase D2 (PLD 2) (mPLD2) (EC 3.1.4.4) (Choline phosphatase 2) (PLD1C) (Phosphatidylcholine-hydrolyzing phospholipase D2)
									addinProteinClass = c('phospholipase')
								} else if(homologene.species.SpeciesDB_Symbol[i] %in% c('Dusp19',toupper('Dusp19'))) {		# 68082|MGI:1915332|Q8K4T5|Dual specificity protein phosphatase 19 (EC 3.1.3.16) (EC 3.1.3.48) (Protein phosphatase SKRP1) (Stress-activated protein kinase pathway-regulating phosphatase 1)
									addinProteinClass = c('protein phosphatase')
								} else if(homologene.species.SpeciesDB_Symbol[i] %in% c('Snrnp200',toupper('Snrnp200'))) {	# 320632|MGI:2444401|Q6P4T2|U5 small nuclear ribonucleoprotein 200 kDa helicase (EC 3.6.4.13) (BRR2 homolog) (U5 snRNP-specific 200 kDa protein) (U5-200KD)
									addinProteinClass = c('RNA helicase')
								} else if(homologene.species.SpeciesDB_Symbol[i] %in% c('Gstp1',toupper('Gstp1'))) {		# 14870|MGI:95865|P19157|Glutathione S-transferase P 1 (Gst P1) (EC 2.5.1.18) (GST YF-YF) (GST class-pi) (GST-piB) (Preadipocyte growth factor)
									addinProteinClass = c('transferase')
								} else if(homologene.species.SpeciesDB_Symbol[i] %in% c('Epha2',toupper('Epha2'))) {		# 13836|MGI:95278|Q03145|Ephrin type-A receptor 2 (EC 2.7.10.1) (Epithelial cell kinase) (Tyrosine-protein kinase receptor ECK) (Tyrosine-protein kinase receptor MPK-5) (Tyrosine-protein kinase receptor SEK-2)
									addinProteinClass = c('tyrosine protein kinase receptor')
								} else if(homologene.species.SpeciesDB_Symbol[i] %in% c('Rnf123',toupper('Rnf123'))) {		# 84585|MGI:2148796|Q5XPI3|E3 ubiquitin-protein ligase RNF123 (EC 2.3.2.27) (Kip1 ubiquitination-promoting complex protein 1) (RING finger protein 123) (RING-type E3 ubiquitin transferase RNF123)
									addinProteinClass = c('ubiquitin-protein ligase')
								}
								# *** manually check for the missed ProteinClass (end) ****************************************
								
								
								# /// search by UniProt ///
								#homologene.species.UniProt.GeneID <- mixedsort(unique(PTHR13_species.GeneID[unlist(PTHR13_species.GeneID.UniProtKB.homologene_species_UniProt.idx[i])]))
								#homologene.species.UniProt.ProteinID <- mixedsort(unique(PTHR13_species.ProteinID[unlist(PTHR13_species.GeneID.UniProtKB.homologene_species_UniProt.idx[i])]))
								#homologene.species.UniProt.SubfamilyID <- mixedsort(unique(PTHR13_species.SubfamilyID[unlist(PTHR13_species.GeneID.UniProtKB.homologene_species_UniProt.idx[i])]))
								#homologene.species.UniProt.FamilyName <- mixedsort(unique(PTHR13_species.FamilyName[unlist(PTHR13_species.GeneID.UniProtKB.homologene_species_UniProt.idx[i])]))
								#homologene.species.UniProt.SubfamilyName <- mixedsort(unique(PTHR13_species.SubfamilyName[unlist(PTHR13_species.GeneID.UniProtKB.homologene_species_UniProt.idx[i])]))
								#homologene.species.UniProt.MF <- mixedsort(unique(PTHR13_species.MF[unlist(PTHR13_species.GeneID.UniProtKB.homologene_species_UniProt.idx[i])]))
								#homologene.species.UniProt.BP <- mixedsort(unique(PTHR13_species.BP[unlist(PTHR13_species.GeneID.UniProtKB.homologene_species_UniProt.idx[i])]))
								#homologene.species.UniProt.CC <- mixedsort(unique(PTHR13_species.CC[unlist(PTHR13_species.GeneID.UniProtKB.homologene_species_UniProt.idx[i])]))
								homologene.species.UniProt.ProteinClass <- mixedsort(unique(PTHR13_species.ProteinClass[unlist(PTHR13_species.GeneID.UniProtKB.homologene_species_UniProt.idx[i])]))
								if(length(addinProteinClass) != 0)		homologene.species.UniProt.ProteinClass = paste(unique(c(unique(unlist(strsplit(homologene.species.UniProt.ProteinClass,';'))), unlist(sapply(addinProteinClass, FUN=function(X) paste(X,'#',unique(Protein_Class.id[which(Protein_Class.name == X)]),sep=""))))), collapse=';')
								#homologene.species.UniProt.Pathway <- mixedsort(unique(PTHR13_species.Pathway[unlist(PTHR13_species.GeneID.UniProtKB.homologene_species_UniProt.idx[i])]))
								#print(homologene.species.UniProt.GeneID)
								#print(homologene.species.UniProt.ProteinID)
								#print(homologene.species.UniProt.SubfamilyID)
								#print(homologene.species.UniProt.FamilyName)
								#print(homologene.species.UniProt.SubfamilyName)
								#print(homologene.species.UniProt.MF)
								#print(homologene.species.UniProt.BP)
								#print(homologene.species.UniProt.CC)
								#print(homologene.species.UniProt.ProteinClass)
								#print(homologene.species.UniProt.Pathway)
								
								# /// search by MGI id ///
								#homologene.species.SpeciesDB_ID.GeneID <- mixedsort(unique(PTHR13_species.GeneID[unlist(PTHR13_species.GeneID.SpeciesDB.homologene_species_SpeciesDB_ID.idx[i])]))
								#homologene.species.SpeciesDB_ID.ProteinID <- mixedsort(unique(PTHR13_species.ProteinID[unlist(PTHR13_species.GeneID.SpeciesDB.homologene_species_SpeciesDB_ID.idx[i])]))
								#homologene.species.SpeciesDB_ID.SubfamilyID <- mixedsort(unique(PTHR13_species.SubfamilyID[unlist(PTHR13_species.GeneID.SpeciesDB.homologene_species_SpeciesDB_ID.idx[i])]))
								#homologene.species.SpeciesDB_ID.FamilyName <- mixedsort(unique(PTHR13_species.FamilyName[unlist(PTHR13_species.GeneID.SpeciesDB.homologene_species_SpeciesDB_ID.idx[i])]))
								#homologene.species.SpeciesDB_ID.SubfamilyName <- mixedsort(unique(PTHR13_species.SubfamilyName[unlist(PTHR13_species.GeneID.SpeciesDB.homologene_species_SpeciesDB_ID.idx[i])]))
								#homologene.species.SpeciesDB_ID.MF <- mixedsort(unique(PTHR13_species.MF[unlist(PTHR13_species.GeneID.SpeciesDB.homologene_species_SpeciesDB_ID.idx[i])]))
								#homologene.species.SpeciesDB_ID.BP <- mixedsort(unique(PTHR13_species.BP[unlist(PTHR13_species.GeneID.SpeciesDB.homologene_species_SpeciesDB_ID.idx[i])]))
								#homologene.species.SpeciesDB_ID.CC <- mixedsort(unique(PTHR13_species.CC[unlist(PTHR13_species.GeneID.SpeciesDB.homologene_species_SpeciesDB_ID.idx[i])]))
								homologene.species.SpeciesDB_ID.ProteinClass <- mixedsort(unique(PTHR13_species.ProteinClass[unlist(PTHR13_species.GeneID.SpeciesDB.homologene_species_SpeciesDB_ID.idx[i])]))
								if(length(addinProteinClass) != 0)		homologene.species.SpeciesDB_ID.ProteinClass = paste(unique(c(unique(unlist(strsplit(homologene.species.SpeciesDB_ID.ProteinClass,';'))), unlist(sapply(addinProteinClass, FUN=function(X) paste(X,'#',unique(Protein_Class.id[which(Protein_Class.name == X)]),sep=""))))), collapse=';')
								#homologene.species.SpeciesDB_ID.Pathway <- mixedsort(unique(PTHR13_species.Pathway[unlist(PTHR13_species.GeneID.SpeciesDB.homologene_species_SpeciesDB_ID.idx[i])]))
								#print(homologene.species.SpeciesDB_ID.GeneID)
								#print(homologene.species.SpeciesDB_ID.ProteinID)
								#print(homologene.species.SpeciesDB_ID.SubfamilyID)
								#print(homologene.species.SpeciesDB_ID.FamilyName)
								#print(homologene.species.SpeciesDB_ID.SubfamilyName)
								#print(homologene.species.SpeciesDB_ID.MF)
								#print(homologene.species.SpeciesDB_ID.BP)
								#print(homologene.species.SpeciesDB_ID.CC)
								#print(homologene.species.SpeciesDB_ID.ProteinClass)
								#print(homologene.species.SpeciesDB_ID.Pathway)
								
								# === ProteinClass ===
								species.UniProt.ProteinClass_id.mat[i,paste(SpeciesName.general[sp],'.UniProt',sep="")] = homologene.species.UniProt[i]
								species.UniProt.ProteinClass_id.mat[i,paste(SpeciesName.general[sp],'.',SpeciesDB,'.ID',sep="")] = homologene.species.SpeciesDB_ID[i]
								species.UniProt.ProteinClass_id.mat[i,paste(SpeciesName.general[sp],'.Gene.ID',sep="")] = homologene.species.Gene_ID[i]
								species.UniProt.ProteinClass_id.mat[i,paste(SpeciesName.general[sp],'.',SpeciesDB,'.Symbol',sep="")] = homologene.species.SpeciesDB_Symbol[i]
								species.UniProt.ProteinClass_id.mat[i,paste(SpeciesName.general[sp],'.Protein.names',sep="")] = homologene.species.Protein_names[i]
								species.UniProt.ProteinClass_level.mat[i,paste(SpeciesName.general[sp],'.UniProt',sep="")] = homologene.species.UniProt[i]
								species.UniProt.ProteinClass_level.mat[i,paste(SpeciesName.general[sp],'.',SpeciesDB,'.ID',sep="")] = homologene.species.SpeciesDB_ID[i]
								species.UniProt.ProteinClass_level.mat[i,paste(SpeciesName.general[sp],'.Gene.ID',sep="")] = homologene.species.Gene_ID[i]
								species.UniProt.ProteinClass_level.mat[i,paste(SpeciesName.general[sp],'.',SpeciesDB,'.Symbol',sep="")] = homologene.species.SpeciesDB_Symbol[i]
								species.UniProt.ProteinClass_level.mat[i,paste(SpeciesName.general[sp],'.Protein.names',sep="")] = homologene.species.Protein_names[i]
								species.UniProt.ProteinClass_name.mat[i,paste(SpeciesName.general[sp],'.UniProt',sep="")] = homologene.species.UniProt[i]
								species.UniProt.ProteinClass_name.mat[i,paste(SpeciesName.general[sp],'.',SpeciesDB,'.ID',sep="")] = homologene.species.SpeciesDB_ID[i]
								species.UniProt.ProteinClass_name.mat[i,paste(SpeciesName.general[sp],'.Gene.ID',sep="")] = homologene.species.Gene_ID[i]
								species.UniProt.ProteinClass_name.mat[i,paste(SpeciesName.general[sp],'.',SpeciesDB,'.Symbol',sep="")] = homologene.species.SpeciesDB_Symbol[i]
								species.UniProt.ProteinClass_name.mat[i,paste(SpeciesName.general[sp],'.Protein.names',sep="")] = homologene.species.Protein_names[i]
								species.UniProt.ProteinClass_description.mat[i,paste(SpeciesName.general[sp],'.UniProt',sep="")] = homologene.species.UniProt[i]
								species.UniProt.ProteinClass_description.mat[i,paste(SpeciesName.general[sp],'.',SpeciesDB,'.ID',sep="")] = homologene.species.SpeciesDB_ID[i]
								species.UniProt.ProteinClass_description.mat[i,paste(SpeciesName.general[sp],'.Gene.ID',sep="")] = homologene.species.Gene_ID[i]
								species.UniProt.ProteinClass_description.mat[i,paste(SpeciesName.general[sp],'.',SpeciesDB,'.Symbol',sep="")] = homologene.species.SpeciesDB_Symbol[i]
								species.UniProt.ProteinClass_description.mat[i,paste(SpeciesName.general[sp],'.Protein.names',sep="")] = homologene.species.Protein_names[i]
								
								hit.case = 0
								if(length(homologene.species.UniProt.ProteinClass) != 0) {	# search by UniProt
									if(length(which(homologene.species.UniProt.ProteinClass != '')) != 0) {
										hit.case = 1
										ProteinClass.name_id <- unique(unlist(strsplit(homologene.species.UniProt.ProteinClass,';')))
									} else {	# search by MGI id
										if(length(homologene.species.SpeciesDB_ID.ProteinClass) != 0) {
											if(length(which(homologene.species.SpeciesDB_ID.ProteinClass != '')) != 0) {
												hit.case = 2
												ProteinClass.name_id <- unique(unlist(strsplit(homologene.species.SpeciesDB_ID.ProteinClass,';')))
											} else {
												#print(paste(homologene.species.UniProt[i],' is in PANTHER (Protein ANalysis THrough Evolutionary Relationships) Classification System, but is not classified into any protein class.',sep=""))
												homologene.species.UniProt.NoProteinClass <- rbind(homologene.species.UniProt.NoProteinClass, cbind(i,homologene.species.UniProt[i],homologene.species.SpeciesDB_ID[i],homologene.species.Gene_ID[i],homologene.species.SpeciesDB_Symbol[i],homologene.species.Protein_names[i]))
											}
										} else {
											#print(paste(homologene.species.UniProt[i],' is in PANTHER (Protein ANalysis THrough Evolutionary Relationships) Classification System, but is not classified into any protein class.',sep=""))
											homologene.species.UniProt.NoProteinClass <- rbind(homologene.species.UniProt.NoProteinClass, cbind(i,homologene.species.UniProt[i],homologene.species.SpeciesDB_ID[i],homologene.species.Gene_ID[i],homologene.species.SpeciesDB_Symbol[i],homologene.species.Protein_names[i]))
										}
									}
								} else {	# search by MGI id
									if(length(homologene.species.SpeciesDB_ID.ProteinClass) != 0) {
										if(length(which(homologene.species.SpeciesDB_ID.ProteinClass != '')) != 0) {
											hit.case = 3
											ProteinClass.name_id <- unique(unlist(strsplit(homologene.species.SpeciesDB_ID.ProteinClass,';')))
										} else {
											#print(paste(homologene.species.UniProt[i],' is not in PANTHER (Protein ANalysis THrough Evolutionary Relationships) Classification System.',sep=""))
											homologene.species.UniProt.NoProteinClass <- rbind(homologene.species.UniProt.NoProteinClass, cbind(i,homologene.species.UniProt[i],homologene.species.SpeciesDB_ID[i],homologene.species.Gene_ID[i],homologene.species.SpeciesDB_Symbol[i],homologene.species.Protein_names[i]))
										}
									} else {
										#print(paste(homologene.species.UniProt[i],' is not in PANTHER (Protein ANalysis THrough Evolutionary Relationships) Classification System.',sep=""))
										homologene.species.UniProt.NA <- rbind(homologene.species.UniProt.NA, cbind(i,homologene.species.UniProt[i],homologene.species.SpeciesDB_ID[i],homologene.species.Gene_ID[i],homologene.species.SpeciesDB_Symbol[i],homologene.species.Protein_names[i]))
									}
								}
								
								#############################################################################################################################
								if(hit.case != 0)
								{
									ProteinClass.name <- sapply(sapply(ProteinClass.name_id, FUN=function(X) strsplit(X,'#')), FUN=function(Y) Y[1])
									ProteinClass.id <- sapply(sapply(ProteinClass.name_id, FUN=function(X) strsplit(X,'#')), FUN=function(Y) Y[2])
									
									ProteinClass.id.mappingidx <- sort(which(Protein_Class.id %in% ProteinClass.id))
									ProteinClass.id.mappingid <- Protein_Class.id[ProteinClass.id.mappingidx]
									ProteinClass.id.mappinglevel <- Protein_Class.level[ProteinClass.id.mappingidx]
									ProteinClass.id.mappingname <- Protein_Class.name[ProteinClass.id.mappingidx]
									ProteinClass.id.mappingdescription <- Protein_Class.description[ProteinClass.id.mappingidx]
									#print(ProteinClass.id.mappingidx)
									
									# +++ add in parent level id ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
									ProteinClass.id.mappinglevel.parent <- c()
									for(j in 1:length(ProteinClass.id.mappinglevel))
									{
										ProteinClass.id.mappinglevel.split = unlist(strsplit(ProteinClass.id.mappinglevel[j],'\\.'))
										
										if(max(which(ProteinClass.id.mappinglevel.split != '00')) >= 3)
										{
											for(k in max(which(ProteinClass.id.mappinglevel.split != '00')):3)
											{
												ProteinClass.id.mappinglevel.parent <- c(ProteinClass.id.mappinglevel.parent, paste(c(ProteinClass.id.mappinglevel.split[1:(k-1)], '00', rep('00',length(ProteinClass.id.mappinglevel.split)-k)), collapse='.'))
											}
										}
									}
									ProteinClass.id.mappinglevel.parent = unique(ProteinClass.id.mappinglevel.parent)
									#print(ProteinClass.id.mappinglevel.parent)
									
									ProteinClass.id.mappinglevel.all <- c(ProteinClass.id.mappinglevel, ProteinClass.id.mappinglevel.parent)
									ProteinClass.id.mappinglevel.all = unique(ProteinClass.id.mappinglevel.all)
									#print(ProteinClass.id.mappinglevel.all)
									
									ProteinClass.id.mappingidx <- sort(which(Protein_Class.level %in% ProteinClass.id.mappinglevel.all))
									ProteinClass.id.mappingid <- Protein_Class.id[ProteinClass.id.mappingidx]
									ProteinClass.id.mappinglevel <- Protein_Class.level[ProteinClass.id.mappingidx]
									ProteinClass.id.mappingname <- Protein_Class.name[ProteinClass.id.mappingidx]
									ProteinClass.id.mappingdescription <- Protein_Class.description[ProteinClass.id.mappingidx]
									#print(ProteinClass.id.mappingidx)
									# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
									
									Level_0.idx = which(sapply(strsplit(Protein_Class[ProteinClass.id.mappingidx,2],'\\.'), FUN=function(X) sum(X!='00')==1))
									Level_1.idx = which(sapply(strsplit(Protein_Class[ProteinClass.id.mappingidx,2],'\\.'), FUN=function(X) sum(X!='00')==2))
									Level_2.idx = which(sapply(strsplit(Protein_Class[ProteinClass.id.mappingidx,2],'\\.'), FUN=function(X) sum(X!='00')==3))
									Level_3.idx = which(sapply(strsplit(Protein_Class[ProteinClass.id.mappingidx,2],'\\.'), FUN=function(X) sum(X!='00')==4))
									Level_4.idx = which(sapply(strsplit(Protein_Class[ProteinClass.id.mappingidx,2],'\\.'), FUN=function(X) sum(X!='00')==5))
									Level_5.idx = which(sapply(strsplit(Protein_Class[ProteinClass.id.mappingidx,2],'\\.'), FUN=function(X) sum(X!='00')==6))
									
									species.UniProt.ProteinClass_id.mat[i,5+1] = paste(unique(Protein_Class[Protein_Class.level_0.idx,1]),collapse=';')
									species.UniProt.ProteinClass_id.mat[i,5+2] = paste(unique(Protein_Class[ProteinClass.id.mappingidx,1][Level_1.idx]),collapse=';')
									species.UniProt.ProteinClass_id.mat[i,5+3] = paste(unique(Protein_Class[ProteinClass.id.mappingidx,1][Level_2.idx]),collapse=';')
									species.UniProt.ProteinClass_id.mat[i,5+4] = paste(unique(Protein_Class[ProteinClass.id.mappingidx,1][Level_3.idx]),collapse=';')
									species.UniProt.ProteinClass_id.mat[i,5+5] = paste(unique(Protein_Class[ProteinClass.id.mappingidx,1][Level_4.idx]),collapse=';')
									species.UniProt.ProteinClass_id.mat[i,5+6] = paste(unique(Protein_Class[ProteinClass.id.mappingidx,1][Level_5.idx]),collapse=';')
									
									species.UniProt.ProteinClass_level.mat[i,5+1] = paste(unique(Protein_Class[Protein_Class.level_0.idx,2]),collapse=';')
									species.UniProt.ProteinClass_level.mat[i,5+2] = paste(unique(Protein_Class[ProteinClass.id.mappingidx,2][Level_1.idx]),collapse=';')
									species.UniProt.ProteinClass_level.mat[i,5+3] = paste(unique(Protein_Class[ProteinClass.id.mappingidx,2][Level_2.idx]),collapse=';')
									species.UniProt.ProteinClass_level.mat[i,5+4] = paste(unique(Protein_Class[ProteinClass.id.mappingidx,2][Level_3.idx]),collapse=';')
									species.UniProt.ProteinClass_level.mat[i,5+5] = paste(unique(Protein_Class[ProteinClass.id.mappingidx,2][Level_4.idx]),collapse=';')
									species.UniProt.ProteinClass_level.mat[i,5+6] = paste(unique(Protein_Class[ProteinClass.id.mappingidx,2][Level_5.idx]),collapse=';')
									
									species.UniProt.ProteinClass_name.mat[i,5+1] = paste(unique(Protein_Class[Protein_Class.level_0.idx,3]),collapse=';')
									species.UniProt.ProteinClass_name.mat[i,5+2] = paste(unique(Protein_Class[ProteinClass.id.mappingidx,3][Level_1.idx]),collapse=';')
									species.UniProt.ProteinClass_name.mat[i,5+3] = paste(unique(Protein_Class[ProteinClass.id.mappingidx,3][Level_2.idx]),collapse=';')
									species.UniProt.ProteinClass_name.mat[i,5+4] = paste(unique(Protein_Class[ProteinClass.id.mappingidx,3][Level_3.idx]),collapse=';')
									species.UniProt.ProteinClass_name.mat[i,5+5] = paste(unique(Protein_Class[ProteinClass.id.mappingidx,3][Level_4.idx]),collapse=';')
									species.UniProt.ProteinClass_name.mat[i,5+6] = paste(unique(Protein_Class[ProteinClass.id.mappingidx,3][Level_5.idx]),collapse=';')
									
									species.UniProt.ProteinClass_description.mat[i,5+1] = paste(unique(Protein_Class[Protein_Class.level_0.idx,4]),collapse=';')
									species.UniProt.ProteinClass_description.mat[i,5+2] = paste(unique(Protein_Class[ProteinClass.id.mappingidx,4][Level_1.idx]),collapse='|')
									species.UniProt.ProteinClass_description.mat[i,5+3] = paste(unique(Protein_Class[ProteinClass.id.mappingidx,4][Level_2.idx]),collapse='|')
									species.UniProt.ProteinClass_description.mat[i,5+4] = paste(unique(Protein_Class[ProteinClass.id.mappingidx,4][Level_3.idx]),collapse='|')
									species.UniProt.ProteinClass_description.mat[i,5+5] = paste(unique(Protein_Class[ProteinClass.id.mappingidx,4][Level_4.idx]),collapse='|')
									species.UniProt.ProteinClass_description.mat[i,5+6] = paste(unique(Protein_Class[ProteinClass.id.mappingidx,4][Level_5.idx]),collapse='|')
									
									
									# --- Level.submatrix -------------------------------------------------------------------------------
									Level_0.name = unique(Protein_Class[Protein_Class.level_0.idx,3])
									Level_1.name = unique(Protein_Class[ProteinClass.id.mappingidx,3][Level_1.idx])
									Level_2.name = unique(Protein_Class[ProteinClass.id.mappingidx,3][Level_2.idx])
									Level_3.name = unique(Protein_Class[ProteinClass.id.mappingidx,3][Level_3.idx])
									Level_4.name = unique(Protein_Class[ProteinClass.id.mappingidx,3][Level_4.idx])
									Level_5.name = unique(Protein_Class[ProteinClass.id.mappingidx,3][Level_5.idx])
									
									Level.submatrix.all <- c()
									
									for(Lv in 0:5)
									{
										Level_n.name = eval(parse(text=paste('Level_', Lv, '.name', sep="")))
										Protein_Class.level_n.id2level.hashtable = eval(parse(text=paste('Protein_Class.level_', Lv, '.id2level.hashtable', sep="")))
										Protein_Class.level_n.id2name.hashtable = eval(parse(text=paste('Protein_Class.level_', Lv, '.id2name.hashtable', sep="")))
										Protein_Class.level_n.id2description.hashtable = eval(parse(text=paste('Protein_Class.level_', Lv, '.id2description.hashtable', sep="")))
										
										if(length(Level_n.name) > 0)
										{
											Level_n.submatrix.all <- c()
											
											for(hh in 1:length(Level_n.name))
											{
												if(has.key(Level_n.name[hh],invert(Protein_Class.level_n.id2name.hashtable)))
												{
													Level_n.name2id = as.character(invert(Protein_Class.level_n.id2name.hashtable)[[Level_n.name[hh]]])
													
													if(has.key(Level_n.name2id,Protein_Class.level_n.id2level.hashtable)) {
														Level_n.name2id2level = as.character(Protein_Class.level_n.id2level.hashtable[[Level_n.name2id]])
													} else {
														Level_n.name2id2level = ''
													}
													
													if(has.key(Level_n.name2id,Protein_Class.level_n.id2name.hashtable)) {
														Level_n.name2id2name = as.character(Protein_Class.level_n.id2name.hashtable[[Level_n.name2id]])
														if(Level_n.name2id2name != Level_n.name[hh]) {
															cat('[Warning] value (name) of "Protein_Class.level_',Lv,'.id2name.hashtable','" (',Level_n.name2id2name,') is not the same as name of "Level.',Lv,'" (',Level_n.name[hh],').\n', sep="")
														}
													} else {
														Level_n.name2id2name = ''
													}
													
													if(has.key(Level_n.name2id,Protein_Class.level_n.id2description.hashtable)) {
														Level_n.name2id2description = as.character(Protein_Class.level_n.id2description.hashtable[[Level_n.name2id]])
													} else {
														Level_n.name2id2description = ''
													}
													
													Level_n.submatrix = data.frame(Lv, Level_n.name2id2name, Level_n.name2id, Level_n.name2id2level, Level_n.name2id2description, check.names=FALSE, stringsAsFactors=FALSE)
													
													Level_n.submatrix.all <- rbind(Level_n.submatrix.all, Level_n.submatrix)
												}
											}
											
											Level.submatrix.all <- rbind(Level.submatrix.all, Level_n.submatrix.all)
										}
									}
									
									if(!is.null(Level.submatrix.all))
									{
										if(nrow(Level.submatrix.all) > 0)
										{
											colnames(Level.submatrix.all) = c('Level','Protein class name','Protein class id','Protein class level','Protein class description')
											rownames(Level.submatrix.all) = c(1:nrow(Level.submatrix.all))
										}
									}
									
									if(!is.null(Level.submatrix.all))
									{
										if(nrow(Level.submatrix.all) > 0)
										{
											na.idx = which(is.na(Level.submatrix.all), arr.ind=TRUE)
											if(nrow(na.idx) > 0)
											{
												for(ri in 1:nrow(na.idx))
												{
													Level.submatrix.all[na.idx[ri,'row'],na.idx[ri,'col']] = ''
												}
											}
										}
									}
									
									if(!is.null(Level.submatrix.all))
									{
										if(nrow(Level.submatrix.all) > 0)
										{
											Title.INFO = c(homologene.species.UniProt[i], homologene.species.SpeciesDB_ID[i], homologene.species.Gene_ID[i], homologene.species.SpeciesDB_Symbol[i], homologene.species.Protein_names[i])
											Title.INFO = Title.INFO[which(!is.na(Title.INFO) & Title.INFO != '')]
											Title.INFO = paste(Title.INFO, collapse="||")
											
											Title = c(paste('> ', Title.INFO, sep=""), rep(NA,ncol(Level.submatrix.all)))
											Header = c(NA, colnames(Level.submatrix.all))
											SubMatrix = as.matrix(cbind(rownames(Level.submatrix.all), Level.submatrix.all))
											BlankLine = c(NA, rep(NA,ncol(Level.submatrix.all)))
											Level.submatrix.all.new = rbind(Title, Header, SubMatrix, BlankLine)
											
											Level.submatrix.all.new.all <- rbind(Level.submatrix.all.new.all, Level.submatrix.all.new)
										}
									}
								}
							}
							
							species.UniProt.ProteinClass_id.mat.empty.idx = which(species.UniProt.ProteinClass_id.mat[,paste(SpeciesName.general[sp],'.UniProt',sep="")] == '')
							if(length(species.UniProt.ProteinClass_id.mat.empty.idx) > 0) {
								species.UniProt.ProteinClass_id.mat.output = species.UniProt.ProteinClass_id.mat[-species.UniProt.ProteinClass_id.mat.empty.idx,,drop=FALSE]
							} else {
								species.UniProt.ProteinClass_id.mat.output = species.UniProt.ProteinClass_id.mat
							}
							
							species.UniProt.ProteinClass_level.mat.empty.idx = which(species.UniProt.ProteinClass_level.mat[,paste(SpeciesName.general[sp],'.UniProt',sep="")] == '')
							if(length(species.UniProt.ProteinClass_level.mat.empty.idx) > 0) {
								species.UniProt.ProteinClass_level.mat.output = species.UniProt.ProteinClass_level.mat[-species.UniProt.ProteinClass_level.mat.empty.idx,,drop=FALSE]
							} else {
								species.UniProt.ProteinClass_level.mat.output = species.UniProt.ProteinClass_level.mat
							}
							
							species.UniProt.ProteinClass_name.mat.empty.idx = which(species.UniProt.ProteinClass_name.mat[,paste(SpeciesName.general[sp],'.UniProt',sep="")] == '')
							if(length(species.UniProt.ProteinClass_name.mat.empty.idx) > 0) {
								species.UniProt.ProteinClass_name.mat.output = species.UniProt.ProteinClass_name.mat[-species.UniProt.ProteinClass_name.mat.empty.idx,,drop=FALSE]
							} else {
								species.UniProt.ProteinClass_name.mat.output = species.UniProt.ProteinClass_name.mat
							}
							
							species.UniProt.ProteinClass_description.mat.empty.idx = which(species.UniProt.ProteinClass_description.mat[,paste(SpeciesName.general[sp],'.UniProt',sep="")] == '')
							if(length(species.UniProt.ProteinClass_description.mat.empty.idx) > 0) {
								species.UniProt.ProteinClass_description.mat.output = species.UniProt.ProteinClass_description.mat[-species.UniProt.ProteinClass_description.mat.empty.idx,,drop=FALSE]
							} else {
								species.UniProt.ProteinClass_description.mat.output = species.UniProt.ProteinClass_description.mat
							}
							
							#write.csv(species.UniProt.ProteinClass_id.mat.output, paste('./Data/',Sample,'/08_UniProt2ProteinClass/subset/n_',r,'/',folder[f],'/',SpeciesName.general[sp],'/',Sample,'_exo_',SpeciesName.general[sp],'_UniProt_ProteinClass_id.csv',sep=""), quote=TRUE, row.names=FALSE)
							#write.csv(species.UniProt.ProteinClass_level.mat.output, paste('./Data/',Sample,'/08_UniProt2ProteinClass/subset/n_',r,'/',folder[f],'/',SpeciesName.general[sp],'/',Sample,'_exo_',SpeciesName.general[sp],'_UniProt_ProteinClass_level.csv',sep=""), quote=TRUE, row.names=FALSE)
							#write.csv(species.UniProt.ProteinClass_name.mat.output, paste('./Data/',Sample,'/08_UniProt2ProteinClass/subset/n_',r,'/',folder[f],'/',SpeciesName.general[sp],'/',Sample,'_exo_',SpeciesName.general[sp],'_UniProt_ProteinClass_name.csv',sep=""), quote=TRUE, row.names=FALSE)
							#write.csv(species.UniProt.ProteinClass_description.mat.output, paste('./Data/',Sample,'/08_UniProt2ProteinClass/subset/n_',r,'/',folder[f],'/',SpeciesName.general[sp],'/',Sample,'_exo_',SpeciesName.general[sp],'_UniProt_ProteinClass_description.csv',sep=""), quote=TRUE, row.names=FALSE)
							
							#if(r == 1 & folder[f] == 'high')
							#{
							#	write.csv(species.UniProt.ProteinClass_id.mat.output, paste('./Data/',Sample,'/08_UniProt2ProteinClass/',Sample,'_exo_',SpeciesName.general[sp],'_UniProt_ProteinClass_id.csv',sep=""), quote=TRUE, row.names=FALSE)
							#	write.csv(species.UniProt.ProteinClass_level.mat.output, paste('./Data/',Sample,'/08_UniProt2ProteinClass/',Sample,'_exo_',SpeciesName.general[sp],'_UniProt_ProteinClass_level.csv',sep=""), quote=TRUE, row.names=FALSE)
							#	write.csv(species.UniProt.ProteinClass_name.mat.output, paste('./Data/',Sample,'/08_UniProt2ProteinClass/',Sample,'_exo_',SpeciesName.general[sp],'_UniProt_ProteinClass_name.csv',sep=""), quote=TRUE, row.names=FALSE)
							#	write.csv(species.UniProt.ProteinClass_description.mat.output, paste('./Data/',Sample,'/08_UniProt2ProteinClass/',Sample,'_exo_',SpeciesName.general[sp],'_UniProt_ProteinClass_description.csv',sep=""), quote=TRUE, row.names=FALSE)
							#}
							
							# --- write xlsx file by openxlsx package ---------------------------------------------------------
							# Importing a big xlsx file into R? https://stackoverflow.com/questions/19147884/importing-a-big-xlsx-file-into-r/43118530#43118530
							# Building R for Windows: https://cran.r-project.org/bin/windows/Rtools/
							#Sys.setenv("R_ZIPCMD" = "C:/Rtools/bin/zip.exe")	# path to zip.exe
							
							wb <- openxlsx::createWorkbook(paste('./Data/',Sample,'/08_UniProt2ProteinClass/subset/n_',r,'/',folder[f],'/',SpeciesName.general[sp],'/',Sample,'_exo_',SpeciesName.general[sp],'_UniProt_ProteinClass.xlsx',sep=""))
							modifyBaseFont(wb, fontSize=12, fontColour="black", fontName="Times New Roman")
							
							wb2 <- openxlsx::createWorkbook(paste('./Data/',Sample,'/08_UniProt2ProteinClass/',Sample,'_exo_',SpeciesName.general[sp],'_UniProt_ProteinClass.xlsx',sep=""))
							modifyBaseFont(wb2, fontSize=12, fontColour="black", fontName="Times New Roman")
							
							for(s in 1:4)
							{
								if(s == 1) {
									sheet = 'ProteinClass_id'
									sheetData = species.UniProt.ProteinClass_id.mat.output
								} else if(s == 2) {
									sheet = 'ProteinClass_level'
									sheetData = species.UniProt.ProteinClass_level.mat.output
								} else if(s == 3) {
									sheet = 'ProteinClass_name'
									sheetData = species.UniProt.ProteinClass_name.mat.output
								} else if(s == 4) {
									sheet = 'ProteinClass_description'
									sheetData = species.UniProt.ProteinClass_description.mat.output
								}
								
								if(!is.null(sheetData))
								{
									if(nrow(sheetData) > 0)
									{
										addWorksheet(wb, sheet)
										writeData(wb, sheet, sheetData, colNames=TRUE, rowNames=FALSE, keepNA=FALSE)
										setColWidths(wb, sheet, cols=c(1:ncol(sheetData)), widths=c(rep(15,3),rep(20,2),rep(12.5,ncol(sheetData)-(3+2))))
										setRowHeights(wb, sheet, rows=c(1:(1+nrow(sheetData))), heights=16.5)
										freezePane(wb, sheet, firstActiveRow=2, firstActiveCol=6)
										style <- createStyle(halign="left", valign="center")
										addStyle(wb, sheet, style, rows=c(1:(1+nrow(sheetData))), cols=c(1:ncol(sheetData)), gridExpand=TRUE)
										style <- createStyle(textDecoration="bold")
										addStyle(wb, sheet, style, rows=1, cols=c(1:ncol(sheetData)), stack=TRUE)
										
										if(r == 1 & folder[f] == 'high')
										{
											addWorksheet(wb2, sheet)
											writeData(wb2, sheet, sheetData, colNames=TRUE, rowNames=FALSE, keepNA=FALSE)
											setColWidths(wb2, sheet, cols=c(1:ncol(sheetData)), widths=c(rep(15,3),rep(20,2),rep(12.5,ncol(sheetData)-(3+2))))
											setRowHeights(wb2, sheet, rows=c(1:(1+nrow(sheetData))), heights=16.5)
											freezePane(wb2, sheet, firstActiveRow=2, firstActiveCol=6)
											style <- createStyle(halign="left", valign="center")
											addStyle(wb2, sheet, style, rows=c(1:(1+nrow(sheetData))), cols=c(1:ncol(sheetData)), gridExpand=TRUE)
											style <- createStyle(textDecoration="bold")
											addStyle(wb2, sheet, style, rows=1, cols=c(1:ncol(sheetData)), stack=TRUE)
										}
									}
								}
							}
							
							if(length(wb$sheet_names) > 0) {
								openxlsx::saveWorkbook(wb, file=paste('./Data/',Sample,'/08_UniProt2ProteinClass/subset/n_',r,'/',folder[f],'/',SpeciesName.general[sp],'/',Sample,'_exo_',SpeciesName.general[sp],'_UniProt_ProteinClass.xlsx',sep=""), overwrite=TRUE)
							}
							
							if(r == 1 & folder[f] == 'high')
							{
								if(length(wb2$sheet_names) > 0) {
									openxlsx::saveWorkbook(wb2, file=paste('./Data/',Sample,'/08_UniProt2ProteinClass/',Sample,'_exo_',SpeciesName.general[sp],'_UniProt_ProteinClass.xlsx',sep=""), overwrite=TRUE)
								}
							}
							
							
							if(!is.null(homologene.species.UniProt.NoProteinClass))
							{
								colnames(homologene.species.UniProt.NoProteinClass) = c('index',mat_colnames)
								#write.csv(homologene.species.UniProt.NoProteinClass, paste('./Data/',Sample,'/08_UniProt2ProteinClass/subset/n_',r,'/',folder[f],'/',SpeciesName.general[sp],'/',Sample,'_exo_',SpeciesName.general[sp],'_UniProt_NoProteinClass.csv',sep=""), quote=TRUE, row.names=FALSE)
								#if(r == 1 & folder[f] == 'high')	write.csv(homologene.species.UniProt.NoProteinClass, paste('./Data/',Sample,'/08_UniProt2ProteinClass/',Sample,'_exo_',SpeciesName.general[sp],'_UniProt_NoProteinClass.csv',sep=""), quote=TRUE, row.names=FALSE)
							}
							if(!is.null(homologene.species.UniProt.NA))
							{
								colnames(homologene.species.UniProt.NA) = c('index',mat_colnames)
								#write.csv(homologene.species.UniProt.NA, paste('./Data/',Sample,'/08_UniProt2ProteinClass/subset/n_',r,'/',folder[f],'/',SpeciesName.general[sp],'/',Sample,'_exo_',SpeciesName.general[sp],'_UniProt_NAinPANTHER.csv',sep=""), quote=TRUE, row.names=FALSE)
								#if(r == 1 & folder[f] == 'high')	write.csv(homologene.species.UniProt.NA, paste('./Data/',Sample,'/08_UniProt2ProteinClass/',Sample,'_exo_',SpeciesName.general[sp],'_UniProt_NAinPANTHER.csv',sep=""), quote=TRUE, row.names=FALSE)
							}
							
							# --- write xlsx file by openxlsx package ---------------------------------------------------------
							# Importing a big xlsx file into R? https://stackoverflow.com/questions/19147884/importing-a-big-xlsx-file-into-r/43118530#43118530
							# Building R for Windows: https://cran.r-project.org/bin/windows/Rtools/
							#Sys.setenv("R_ZIPCMD" = "C:/Rtools/bin/zip.exe")	# path to zip.exe
							
							wb <- openxlsx::createWorkbook(paste('./Data/',Sample,'/08_UniProt2ProteinClass/subset/n_',r,'/',folder[f],'/',SpeciesName.general[sp],'/',Sample,'_exo_',SpeciesName.general[sp],'_UniProt_NoProteinClass.xlsx',sep=""))
							modifyBaseFont(wb, fontSize=12, fontColour="black", fontName="Times New Roman")
							
							wb2 <- openxlsx::createWorkbook(paste('./Data/',Sample,'/08_UniProt2ProteinClass/',Sample,'_exo_',SpeciesName.general[sp],'_UniProt_NoProteinClass.xlsx',sep=""))
							modifyBaseFont(wb2, fontSize=12, fontColour="black", fontName="Times New Roman")
							
							for(s in 1:2)
							{
								if(s == 1) {
									sheet = 'No Protein Class'
									sheetData = homologene.species.UniProt.NoProteinClass
								} else if(s == 2) {
									sheet = 'NA in PANTHER'
									sheetData = homologene.species.UniProt.NA
								}
								
								if(!is.null(sheetData))
								{
									if(nrow(sheetData) > 0)
									{
										addWorksheet(wb, sheet)
										writeData(wb, sheet, sheetData, colNames=TRUE, rowNames=FALSE, keepNA=FALSE)
										setColWidths(wb, sheet, cols=c(1:ncol(sheetData)), widths=c(8.5,rep(15,3),rep(20,2),rep(12.5,ncol(sheetData)-(1+3+2))))
										setRowHeights(wb, sheet, rows=c(1:(1+nrow(sheetData))), heights=16.5)
										freezePane(wb, sheet, firstActiveRow=2)
										style <- createStyle(halign="left", valign="center")
										addStyle(wb, sheet, style, rows=c(1:(1+nrow(sheetData))), cols=c(1:ncol(sheetData)), gridExpand=TRUE)
										style <- createStyle(textDecoration="bold")
										addStyle(wb, sheet, style, rows=1, cols=c(1:ncol(sheetData)), stack=TRUE)
										
										if(r == 1 & folder[f] == 'high')
										{
											addWorksheet(wb2, sheet)
											writeData(wb2, sheet, sheetData, colNames=TRUE, rowNames=FALSE, keepNA=FALSE)
											setColWidths(wb2, sheet, cols=c(1:ncol(sheetData)), widths=c(8.5,rep(15,3),rep(20,2),rep(12.5,ncol(sheetData)-(1+3+2))))
											setRowHeights(wb2, sheet, rows=c(1:(1+nrow(sheetData))), heights=16.5)
											freezePane(wb2, sheet, firstActiveRow=2)
											style <- createStyle(halign="left", valign="center")
											addStyle(wb2, sheet, style, rows=c(1:(1+nrow(sheetData))), cols=c(1:ncol(sheetData)), gridExpand=TRUE)
											style <- createStyle(textDecoration="bold")
											addStyle(wb2, sheet, style, rows=1, cols=c(1:ncol(sheetData)), stack=TRUE)
										}
									}
								}
							}
							
							if(length(wb$sheet_names) > 0) {
								openxlsx::saveWorkbook(wb, file=paste('./Data/',Sample,'/08_UniProt2ProteinClass/subset/n_',r,'/',folder[f],'/',SpeciesName.general[sp],'/',Sample,'_exo_',SpeciesName.general[sp],'_UniProt_NoProteinClass.xlsx',sep=""), overwrite=TRUE)
							}
							
							if(r == 1 & folder[f] == 'high')
							{
								if(length(wb2$sheet_names) > 0) {
									openxlsx::saveWorkbook(wb2, file=paste('./Data/',Sample,'/08_UniProt2ProteinClass/',Sample,'_exo_',SpeciesName.general[sp],'_UniProt_NoProteinClass.xlsx',sep=""), overwrite=TRUE)
								}
							}
							
							
							# === Protein Class (Level 0-5) to MS-identified mouse proteins and their human homolog proteins ================================
							# --- write xlsx file by openxlsx package ---------------------------------------------------------
							# Importing a big xlsx file into R? https://stackoverflow.com/questions/19147884/importing-a-big-xlsx-file-into-r/43118530#43118530
							# Building R for Windows: https://cran.r-project.org/bin/windows/Rtools/
							#Sys.setenv("R_ZIPCMD" = "C:/Rtools/bin/zip.exe")	# path to zip.exe
							
							for(s in 1:3)
							{
								if(s == 1) {
									assign(paste('prename',s,'1',sep=""), paste('./Data/',Sample,'/08_UniProt2ProteinClass/subset/n_',r,'/',folder[f],'/',SpeciesName.general[sp],'/',Sample,'_exo_',SpeciesName.general[sp],'_','ProteinClass',sep=""))
									assign(paste('prename',s,'2',sep=""), paste('./Data/',Sample,'/08_UniProt2ProteinClass/',Sample,'_exo_',SpeciesName.general[sp],'_','ProteinClass',sep=""))
								} else if(s == 2) {
									assign(paste('prename',s,'1',sep=""), paste('./Data/',Sample,'/08_UniProt2ProteinClass/subset/n_',r,'/',folder[f],'/',SpeciesName.general[sp],'/',Sample,'_exo_',SpeciesName.general[sp],'_','ProteinClass','_LevelNameSort',sep=""))
									assign(paste('prename',s,'2',sep=""), paste('./Data/',Sample,'/08_UniProt2ProteinClass/',Sample,'_exo_',SpeciesName.general[sp],'_','ProteinClass','_LevelNameSort',sep=""))
								} else if(s == 3) {
									assign(paste('prename',s,'1',sep=""), paste('./Data/',Sample,'/08_UniProt2ProteinClass/subset/n_',r,'/',folder[f],'/',SpeciesName.general[sp],'/',Sample,'_exo_',SpeciesName.general[sp],'_','ProteinClass','_ProteinNumSort',sep=""))
									assign(paste('prename',s,'2',sep=""), paste('./Data/',Sample,'/08_UniProt2ProteinClass/',Sample,'_exo_',SpeciesName.general[sp],'_','ProteinClass','_ProteinNumSort',sep=""))
								}
								
								assign(paste('wb',s,'1',sep=""), openxlsx::createWorkbook(paste(eval(parse(text=paste('prename',s,'1',sep=""))),'.xlsx',sep="")))
								modifyBaseFont(eval(parse(text=paste('wb',s,'1',sep=""))), fontSize=12, fontColour="black", fontName="Times New Roman")
								
								assign(paste('wb',s,'2',sep=""), openxlsx::createWorkbook(paste(eval(parse(text=paste('prename',s,'2',sep=""))),'.xlsx',sep="")))
								modifyBaseFont(eval(parse(text=paste('wb',s,'2',sep=""))), fontSize=12, fontColour="black", fontName="Times New Roman")
								
								for(s0 in 1:2)
								{
									if(s0 == 1) {
										sheet = 'Protein2ProteinClass'
										sheetData = species.UniProt.ProteinClass_name.mat.output
									} else if(s0 == 2) {
										sheet = 'Protein2ProteinClass_name'
										sheetData = Level.submatrix.all.new.all
									}
									
									if(!is.null(sheetData))
									{
										if(nrow(sheetData) > 0)
										{
											if(s0 == 1) {
												addWorksheet(eval(parse(text=paste('wb',s,'1',sep=""))), sheet)
												writeData(eval(parse(text=paste('wb',s,'1',sep=""))), sheet, sheetData, colNames=TRUE, rowNames=FALSE, keepNA=FALSE)
												setColWidths(eval(parse(text=paste('wb',s,'1',sep=""))), sheet, cols=c(1:ncol(sheetData)), widths=c(rep(15,3),rep(20,2),rep(12.5,ncol(sheetData)-(3+2))))
												setRowHeights(eval(parse(text=paste('wb',s,'1',sep=""))), sheet, rows=c(1:(1+nrow(sheetData))), heights=16.5)
												freezePane(eval(parse(text=paste('wb',s,'1',sep=""))), sheet, firstActiveRow=2, firstActiveCol=6)
												style <- createStyle(halign="left", valign="center")
												addStyle(eval(parse(text=paste('wb',s,'1',sep=""))), sheet, style, rows=c(1:(1+nrow(sheetData))), cols=c(1:ncol(sheetData)), gridExpand=TRUE)
												style <- createStyle(textDecoration="bold")
												addStyle(eval(parse(text=paste('wb',s,'1',sep=""))), sheet, style, rows=1, cols=c(1:ncol(sheetData)), stack=TRUE)
											} else if(s0 == 2) {
												addWorksheet(eval(parse(text=paste('wb',s,'1',sep=""))), sheet)
												writeData(eval(parse(text=paste('wb',s,'1',sep=""))), sheet, sheetData, colNames=FALSE, rowNames=FALSE, keepNA=FALSE)
												setColWidths(eval(parse(text=paste('wb',s,'1',sep=""))), sheet, cols=c(1:ncol(sheetData)), widths=c(8.5,12.5,40,rep(17.5,ncol(sheetData)-(1+1+1+1)),40))
												setRowHeights(eval(parse(text=paste('wb',s,'1',sep=""))), sheet, rows=c(1:nrow(sheetData)), heights=16.5)
												style <- createStyle(halign="left", valign="center")
												addStyle(eval(parse(text=paste('wb',s,'1',sep=""))), sheet, style, rows=c(1:nrow(sheetData)), cols=c(1:ncol(sheetData)), gridExpand=TRUE)
											}
											
											if(r == 1 & folder[f] == 'high')
											{
												if(s0 == 1) {
													addWorksheet(eval(parse(text=paste('wb',s,'2',sep=""))), sheet)
													writeData(eval(parse(text=paste('wb',s,'2',sep=""))), sheet, sheetData, colNames=TRUE, rowNames=FALSE, keepNA=FALSE)
													setColWidths(eval(parse(text=paste('wb',s,'2',sep=""))), sheet, cols=c(1:ncol(sheetData)), widths=c(rep(15,3),rep(20,2),rep(12.5,ncol(sheetData)-(3+2))))
													setRowHeights(eval(parse(text=paste('wb',s,'2',sep=""))), sheet, rows=c(1:(1+nrow(sheetData))), heights=16.5)
													freezePane(eval(parse(text=paste('wb',s,'2',sep=""))), sheet, firstActiveRow=2, firstActiveCol=6)
													style <- createStyle(halign="left", valign="center")
													addStyle(eval(parse(text=paste('wb',s,'2',sep=""))), sheet, style, rows=c(1:(1+nrow(sheetData))), cols=c(1:ncol(sheetData)), gridExpand=TRUE)
													style <- createStyle(textDecoration="bold")
													addStyle(eval(parse(text=paste('wb',s,'2',sep=""))), sheet, style, rows=1, cols=c(1:ncol(sheetData)), stack=TRUE)
												} else if(s0 == 2) {
													addWorksheet(eval(parse(text=paste('wb',s,'2',sep=""))), sheet)
													writeData(eval(parse(text=paste('wb',s,'2',sep=""))), sheet, sheetData, colNames=FALSE, rowNames=FALSE, keepNA=FALSE)
													setColWidths(eval(parse(text=paste('wb',s,'2',sep=""))), sheet, cols=c(1:ncol(sheetData)), widths=c(8.5,12.5,40,rep(17.5,ncol(sheetData)-(1+1+1+1)),40))
													setRowHeights(eval(parse(text=paste('wb',s,'2',sep=""))), sheet, rows=c(1:nrow(sheetData)), heights=16.5)
													style <- createStyle(halign="left", valign="center")
													addStyle(eval(parse(text=paste('wb',s,'2',sep=""))), sheet, style, rows=c(1:nrow(sheetData)), cols=c(1:ncol(sheetData)), gridExpand=TRUE)
												}
											}
										}
									}
								}
							}
							
							
							for(Lv in 0:5)
							{
								mainDir <- paste('./Data/',Sample,'/08_UniProt2ProteinClass/subset/n_',r,'/',folder[f],'/',SpeciesName.general[sp],sep="")
								subDir <- paste('Level',Lv,sep="")
								if (file.exists(file.path(mainDir, subDir))){
								} else {
									dir.create(file.path(mainDir, subDir))
								}
								
								###############################################################################################
								Protein_Class.level.data = eval(parse(text=paste('Protein_Class.level_', Lv, sep="")))
								#Protein_Class.level.data
								
								###############################################################################################
								# /// LevelIDSort ///
								Protein_Class.level.data.species <- matrix('',nrow(Protein_Class.level.data),5,byrow=TRUE)
								colnames(Protein_Class.level.data.species) = c(paste('Level.',Lv,sep=""),paste(SpeciesName.general[sp],'.UniProt',sep=""),paste(SpeciesName.general[sp],'.Protein.names',sep=""),paste(SpeciesName.general[sp],'.Protein.number',sep=""),paste(SpeciesName.general[sp],'.Protein.percent',sep=""))
								
								Protein_Class.level.data.species.Protein.number.vector <- c()
								for(i in 1:nrow(Protein_Class.level.data))
								{
									Protein_Class.level.data.species.idx <- c()
									
									for(j in 1:nrow(species.UniProt.ProteinClass_name.mat))
									{
										if(Protein_Class.level.data[i,3] %in% unique(unlist(strsplit(species.UniProt.ProteinClass_name.mat[j,paste('Level.',Lv,sep="")],';'))))
										{
											Protein_Class.level.data.species.idx <- c(Protein_Class.level.data.species.idx, j)
										}
									}
									
									Protein_Class.level.data.species.UniProt = paste(as.matrix(homologene)[Protein_Class.level.data.species.idx,c(paste(SpeciesName.general[sp],'.UniProt',sep=""))],collapse=';')
									Protein_Class.level.data.species.Protein.names = paste(as.matrix(homologene)[Protein_Class.level.data.species.idx,c(paste(SpeciesName.general[sp],'.Protein.names',sep=""))],collapse='|')
									Protein_Class.level.data.species.Protein.number = length(Protein_Class.level.data.species.idx)
									Protein_Class.level.data.species.Protein.number.vector <- c(Protein_Class.level.data.species.Protein.number.vector, Protein_Class.level.data.species.Protein.number)
									
									Protein_Class.level.data.species[i,paste('Level.',Lv,sep="")] = Protein_Class.level.data[i,3]
									Protein_Class.level.data.species[i,paste(SpeciesName.general[sp],'.UniProt',sep="")] = Protein_Class.level.data.species.UniProt
									Protein_Class.level.data.species[i,paste(SpeciesName.general[sp],'.Protein.names',sep="")] = Protein_Class.level.data.species.Protein.names
									Protein_Class.level.data.species[i,paste(SpeciesName.general[sp],'.Protein.number',sep="")] = Protein_Class.level.data.species.Protein.number
								}
								
								for(i in 1:nrow(Protein_Class.level.data))
								{
									if(sum(Protein_Class.level.data.species.Protein.number.vector) != 0) {
										Protein_Class.level.data.species[i,paste(SpeciesName.general[sp],'.Protein.percent',sep="")] = as.numeric(as.character(round_preserve_sum((Protein_Class.level.data.species.Protein.number.vector[i]/sum(Protein_Class.level.data.species.Protein.number.vector))*100,1)))
									} else {
										Protein_Class.level.data.species[i,paste(SpeciesName.general[sp],'.Protein.percent',sep="")] = 0
									}
								}
								
								#write.csv(Protein_Class.level.data.species, paste('./Data/',Sample,'/08_UniProt2ProteinClass/subset/n_',r,'/',folder[f],'/',SpeciesName.general[sp],'/',paste('Level',Lv,sep=""),'/',Sample,'_exo_',SpeciesName.general[sp],'_','ProteinClass_Level',Lv,'.csv',sep=""), quote=TRUE, row.names=FALSE)
								
								if(nrow(Protein_Class.level.data.species) > 0)
								{
									Protein_Class.level.data.species.submatrix.new.all <- c()
									
									sif2UniProt <- c()
									sif2SpeciesDB <- c()
									sif2Gene <- c()
									sif2GeneName <- c()
									sif2ProteinName <- c()
									
									for(m in 1:nrow(Protein_Class.level.data.species))
									{
										if(Protein_Class.level.data.species[m,paste(SpeciesName.general[sp],'.Protein.number',sep="")] != 0)
										{
											class_name = Protein_Class.level.data.species[m,paste('Level.',Lv,sep="")]
											class_name.idx <- which(Protein_Class.level.data[,3] %in% class_name)
											class_id = Protein_Class.level.data[class_name.idx,1]
											class_level = Protein_Class.level.data[class_name.idx,2]
											class_description = Protein_Class.level.data[class_name.idx,4]
											
											protein.idx = sapply(unique(unlist(strsplit(Protein_Class.level.data.species[m,paste(SpeciesName.general[sp],'.UniProt',sep="")],';'))), FUN=function(X) which(species.UniProt.ProteinClass_name.mat[,paste(SpeciesName.general[sp],'.UniProt',sep="")] %in% X)[1])
											protein_info = species.UniProt.ProteinClass_name.mat[protein.idx,,drop=FALSE]
											
											Protein_Class.level.data.species.submatrix <- protein_info
											colnames(Protein_Class.level.data.species.submatrix) = colnames(protein_info)
											rownames(Protein_Class.level.data.species.submatrix) = c(1:nrow(Protein_Class.level.data.species.submatrix))
											
											if(!is.null(Protein_Class.level.data.species.submatrix))
											{
												if(nrow(Protein_Class.level.data.species.submatrix) > 0)
												{
													na.idx = which(is.na(Protein_Class.level.data.species.submatrix), arr.ind=TRUE)
													if(nrow(na.idx) > 0)
													{
														for(ri in 1:nrow(na.idx))
														{
															Protein_Class.level.data.species.submatrix[na.idx[ri,'row'],na.idx[ri,'col']] = ''
														}
													}
												}
											}
											
											if(!is.null(Protein_Class.level.data.species.submatrix))
											{
												if(nrow(Protein_Class.level.data.species.submatrix) > 0)
												{
													Title.INFO = c(class_name, class_id, class_level, class_description)
													Title.INFO = Title.INFO[which(!is.na(Title.INFO) & Title.INFO != '')]
													Title.INFO = paste(Title.INFO, collapse="|")
													
													Title = c(paste('> ', Title.INFO, sep=""), rep(NA,ncol(Protein_Class.level.data.species.submatrix)))
													Header = c(NA, colnames(Protein_Class.level.data.species.submatrix))
													SubMatrix = as.matrix(cbind(rownames(Protein_Class.level.data.species.submatrix), Protein_Class.level.data.species.submatrix))
													BlankLine = c(NA, rep(NA,ncol(Protein_Class.level.data.species.submatrix)))
													Protein_Class.level.data.species.submatrix.new = rbind(Title, Header, SubMatrix, BlankLine)
													
													Protein_Class.level.data.species.submatrix.new.all <- rbind(Protein_Class.level.data.species.submatrix.new.all, Protein_Class.level.data.species.submatrix.new)
													
													
													# --- Cytoscape -----------------------------------------------------------------------
													Source <- class_name
													
													Relationship2UniProt <- paste('ProteinClass','2','UniProt',sep="")
													Target2UniProt <- Protein_Class.level.data.species.submatrix[,paste(SpeciesName.general[sp],'.UniProt',sep="")]
													Target2UniProt.NotEmpty.idx = which(Target2UniProt!='')
													if(length(Target2UniProt.NotEmpty.idx) > 0) {
														Target2UniProt.NotEmpty = unique(Target2UniProt[Target2UniProt.NotEmpty.idx])
														sif2UniProt <- rbind(sif2UniProt, cbind(Source, Relationship2UniProt, Target2UniProt.NotEmpty))
													}
													
													Relationship2SpeciesDB <- paste('ProteinClass','2',SpeciesDB,sep="")
													Target2SpeciesDB <- Protein_Class.level.data.species.submatrix[,paste(SpeciesName.general[sp],'.',SpeciesDB,'.ID',sep="")]
													Target2SpeciesDB.NotEmpty.idx = which(Target2SpeciesDB!='')
													if(length(Target2SpeciesDB.NotEmpty.idx) > 0) {
														Target2SpeciesDB.NotEmpty = unique(Target2SpeciesDB[Target2SpeciesDB.NotEmpty.idx])
														sif2SpeciesDB <- rbind(sif2SpeciesDB, cbind(Source, Relationship2SpeciesDB, Target2SpeciesDB.NotEmpty))
													}
													
													Relationship2Gene <- paste('ProteinClass','2','Gene',sep="")
													Target2Gene <- Protein_Class.level.data.species.submatrix[,paste(SpeciesName.general[sp],'.Gene.ID',sep="")]
													Target2Gene.NotEmpty.idx = which(Target2Gene!='')
													if(length(Target2Gene.NotEmpty.idx) > 0) {
														Target2Gene.NotEmpty = unique(Target2Gene[Target2Gene.NotEmpty.idx])
														sif2Gene <- rbind(sif2Gene, cbind(Source, Relationship2Gene, Target2Gene.NotEmpty))
													}
													
													Relationship2GeneName <- paste('ProteinClass','2','GeneName',sep="")
													Target2GeneName <- Protein_Class.level.data.species.submatrix[,paste(SpeciesName.general[sp],'.',SpeciesDB,'.Symbol',sep="")]
													Target2GeneName.NotEmpty.idx = which(Target2GeneName!='')
													if(length(Target2GeneName.NotEmpty.idx) > 0) {
														Target2GeneName.NotEmpty = unique(Target2GeneName[Target2GeneName.NotEmpty.idx])
														sif2GeneName <- rbind(sif2GeneName, cbind(Source, Relationship2GeneName, Target2GeneName.NotEmpty))
													}
													
													Relationship2ProteinName <- paste('ProteinClass','2','ProteinName',sep="")
													Target2ProteinName <- Protein_Class.level.data.species.submatrix[,paste(SpeciesName.general[sp],'.Protein.names',sep="")]
													Target2ProteinName.NotEmpty.idx = which(Target2ProteinName!='')
													if(length(Target2ProteinName.NotEmpty.idx) > 0) {
														Target2ProteinName.NotEmpty = unique(Target2ProteinName[Target2ProteinName.NotEmpty.idx])
														sif2ProteinName <- rbind(sif2ProteinName, cbind(Source, Relationship2ProteinName, Target2ProteinName.NotEmpty))
													}
												}
											}
										}
									}
									
									DB <- c('UniProt',SpeciesDB,'Gene','GeneName','ProteinName')
									for(db in 1:length(DB))
									{
										if(DB[db] == 'UniProt') {
											sif = sif2UniProt
										} else if(DB[db] == SpeciesDB) {
											sif = sif2SpeciesDB
										} else if(DB[db] == 'Gene') {
											sif = sif2Gene
										} else if(DB[db] == 'GeneName') {
											sif = sif2GeneName
										} else if(DB[db] == 'ProteinName') {
											sif = sif2ProteinName
										}
										
										if(!is.null(sif))
										{
											sif[,1] = sapply(sif[,1], FUN=function(X) unlist(strsplit(X," "))[1])
											sif = unique(sif)
											write.table(sif, paste('./Data/',Sample,'/08_UniProt2ProteinClass/subset/n_',r,'/',folder[f],'/',SpeciesName.general[sp],'/',paste('Level',Lv,sep=""),'/',Sample,'_exo_',SpeciesName.general[sp],'_','ProteinClass_Level',Lv,'_',DB[db],'.sif',sep=""), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
											
											Source <- mixedsort(unique(sif[,1]))
											Target <- mixedsort(unique(sif[,3]))
											noa <- rbind(cbind(Source,'ProteinClass'), cbind(Target,DB[db]))
											colnames(noa) <- c('node','attribute')
											noa = unique(noa)
											#print(noa)
											write.table(noa, paste('./Data/',Sample,'/08_UniProt2ProteinClass/subset/n_',r,'/',folder[f],'/',SpeciesName.general[sp],'/',paste('Level',Lv,sep=""),'/',Sample,'_exo_',SpeciesName.general[sp],'_','ProteinClass_Level',Lv,'_',DB[db],'.noa',sep=""), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
										}
									}
								}
								
								
								# /// LevelNameSort ///
								Protein_Class.level.data.species.LevelNameSort <- Protein_Class.level.data.species[order(Protein_Class.level.data.species[,paste('Level.',Lv,sep="")], decreasing=FALSE),,drop=FALSE]
								#write.csv(Protein_Class.level.data.species.LevelNameSort, paste('./Data/',Sample,'/08_UniProt2ProteinClass/subset/n_',r,'/',folder[f],'/',SpeciesName.general[sp],'/',paste('Level',Lv,sep=""),'/',Sample,'_exo_',SpeciesName.general[sp],'_','ProteinClass_Level',Lv,'_LevelNameSort.csv',sep=""), quote=TRUE, row.names=FALSE)
								
								if(nrow(Protein_Class.level.data.species.LevelNameSort) > 0)
								{
									Protein_Class.level.data.species.LevelNameSort.submatrix.new.all <- c()
									
									for(m in 1:nrow(Protein_Class.level.data.species.LevelNameSort))
									{
										if(Protein_Class.level.data.species.LevelNameSort[m,paste(SpeciesName.general[sp],'.Protein.number',sep="")] != 0)
										{
											class_name = Protein_Class.level.data.species.LevelNameSort[m,paste('Level.',Lv,sep="")]
											class_name.idx <- which(Protein_Class.level.data[,3] %in% class_name)
											class_id = Protein_Class.level.data[class_name.idx,1]
											class_level = Protein_Class.level.data[class_name.idx,2]
											class_description = Protein_Class.level.data[class_name.idx,4]
											
											protein.idx = sapply(unique(unlist(strsplit(Protein_Class.level.data.species.LevelNameSort[m,paste(SpeciesName.general[sp],'.UniProt',sep="")],';'))), FUN=function(X) which(species.UniProt.ProteinClass_name.mat[,paste(SpeciesName.general[sp],'.UniProt',sep="")] %in% X)[1])
											protein_info = species.UniProt.ProteinClass_name.mat[protein.idx,,drop=FALSE]
											
											Protein_Class.level.data.species.LevelNameSort.submatrix <- protein_info
											colnames(Protein_Class.level.data.species.LevelNameSort.submatrix) = colnames(protein_info)
											rownames(Protein_Class.level.data.species.LevelNameSort.submatrix) = c(1:nrow(Protein_Class.level.data.species.LevelNameSort.submatrix))
											
											if(!is.null(Protein_Class.level.data.species.LevelNameSort.submatrix))
											{
												if(nrow(Protein_Class.level.data.species.LevelNameSort.submatrix) > 0)
												{
													na.idx = which(is.na(Protein_Class.level.data.species.LevelNameSort.submatrix), arr.ind=TRUE)
													if(nrow(na.idx) > 0)
													{
														for(ri in 1:nrow(na.idx))
														{
															Protein_Class.level.data.species.LevelNameSort.submatrix[na.idx[ri,'row'],na.idx[ri,'col']] = ''
														}
													}
												}
											}
											
											if(!is.null(Protein_Class.level.data.species.LevelNameSort.submatrix))
											{
												if(nrow(Protein_Class.level.data.species.LevelNameSort.submatrix) > 0)
												{
													Title.INFO = c(class_name, class_id, class_level, class_description)
													Title.INFO = Title.INFO[which(!is.na(Title.INFO) & Title.INFO != '')]
													Title.INFO = paste(Title.INFO, collapse="|")
													
													Title = c(paste('> ', Title.INFO, sep=""), rep(NA,ncol(Protein_Class.level.data.species.LevelNameSort.submatrix)))
													Header = c(NA, colnames(Protein_Class.level.data.species.LevelNameSort.submatrix))
													SubMatrix = as.matrix(cbind(rownames(Protein_Class.level.data.species.LevelNameSort.submatrix), Protein_Class.level.data.species.LevelNameSort.submatrix))
													BlankLine = c(NA, rep(NA,ncol(Protein_Class.level.data.species.LevelNameSort.submatrix)))
													Protein_Class.level.data.species.LevelNameSort.submatrix.new = rbind(Title, Header, SubMatrix, BlankLine)
													
													Protein_Class.level.data.species.LevelNameSort.submatrix.new.all <- rbind(Protein_Class.level.data.species.LevelNameSort.submatrix.new.all, Protein_Class.level.data.species.LevelNameSort.submatrix.new)
												}
											}
										}
									}
								}
								
								
								# /// ProteinNumSort ///
								Protein_Class.level.data.species.ProteinNumSort <- Protein_Class.level.data.species[order(as.numeric(as.character(Protein_Class.level.data.species[,paste(SpeciesName.general[sp],'.Protein.number',sep="")])), decreasing=TRUE),,drop=FALSE]
								#write.csv(Protein_Class.level.data.species.ProteinNumSort, paste('./Data/',Sample,'/08_UniProt2ProteinClass/subset/n_',r,'/',folder[f],'/',SpeciesName.general[sp],'/',paste('Level',Lv,sep=""),'/',Sample,'_exo_',SpeciesName.general[sp],'_','ProteinClass_Level',Lv,'_ProteinNumSort.csv',sep=""), quote=TRUE, row.names=FALSE)
								
								if(nrow(Protein_Class.level.data.species.ProteinNumSort) > 0)
								{
									if(sum(as.numeric(as.character(Protein_Class.level.data.species.ProteinNumSort[,paste(SpeciesName.general[sp],'.Protein.number',sep="")]))) != 0)
									{
										mainDir <- paste('./Data/',Sample,'/08_UniProt2ProteinClass/subset/n_',r,'/',folder[f],'/',SpeciesName.general[sp],'/',paste('Level',Lv,sep=""),sep="")
										subDir <- 'id'
										if (file.exists(file.path(mainDir, subDir))){
										} else {
											dir.create(file.path(mainDir, subDir))
										}
										mainDir <- paste('./Data/',Sample,'/08_UniProt2ProteinClass/subset/n_',r,'/',folder[f],'/',SpeciesName.general[sp],'/',paste('Level',Lv,sep=""),'/id',sep="")
										subDir <- 'UniProt'
										if (file.exists(file.path(mainDir, subDir))){
										} else {
											dir.create(file.path(mainDir, subDir))
										}
										mainDir <- paste('./Data/',Sample,'/08_UniProt2ProteinClass/subset/n_',r,'/',folder[f],'/',SpeciesName.general[sp],'/',paste('Level',Lv,sep=""),'/id',sep="")
										subDir <- SpeciesDB
										if (file.exists(file.path(mainDir, subDir))){
										} else {
											dir.create(file.path(mainDir, subDir))
										}
										mainDir <- paste('./Data/',Sample,'/08_UniProt2ProteinClass/subset/n_',r,'/',folder[f],'/',SpeciesName.general[sp],'/',paste('Level',Lv,sep=""),'/id',sep="")
										subDir <- 'Gene'
										if (file.exists(file.path(mainDir, subDir))){
										} else {
											dir.create(file.path(mainDir, subDir))
										}
										
										mainDir <- paste('./Data/',Sample,'/08_UniProt2ProteinClass/subset/n_',r,'/',folder[f],'/',SpeciesName.general[sp],'/',paste('Level',Lv,sep=""),sep="")
										subDir <- 'name'
										if (file.exists(file.path(mainDir, subDir))){
										} else {
											dir.create(file.path(mainDir, subDir))
										}
										mainDir <- paste('./Data/',Sample,'/08_UniProt2ProteinClass/subset/n_',r,'/',folder[f],'/',SpeciesName.general[sp],'/',paste('Level',Lv,sep=""),'/name',sep="")
										subDir <- 'Protein'
										if (file.exists(file.path(mainDir, subDir))){
										} else {
											dir.create(file.path(mainDir, subDir))
										}
									}
									
									
									Protein_Class.level.data.species.ProteinNumSort.submatrix.new.all <- c()
									
									for(m in 1:nrow(Protein_Class.level.data.species.ProteinNumSort))
									{
										if(Protein_Class.level.data.species.ProteinNumSort[m,paste(SpeciesName.general[sp],'.Protein.number',sep="")] != 0)
										{
											class_name = Protein_Class.level.data.species.ProteinNumSort[m,paste('Level.',Lv,sep="")]
											class_name.idx <- which(Protein_Class.level.data[,3] %in% class_name)
											class_id = Protein_Class.level.data[class_name.idx,1]
											class_level = Protein_Class.level.data[class_name.idx,2]
											class_description = Protein_Class.level.data[class_name.idx,4]
											
											protein.idx = sapply(unique(unlist(strsplit(Protein_Class.level.data.species.ProteinNumSort[m,paste(SpeciesName.general[sp],'.UniProt',sep="")],';'))), FUN=function(X) which(species.UniProt.ProteinClass_name.mat[,paste(SpeciesName.general[sp],'.UniProt',sep="")] %in% X)[1])
											protein_info = species.UniProt.ProteinClass_name.mat[protein.idx,,drop=FALSE]
											
											Protein_Class.level.data.species.ProteinNumSort.submatrix <- protein_info
											colnames(Protein_Class.level.data.species.ProteinNumSort.submatrix) = colnames(protein_info)
											rownames(Protein_Class.level.data.species.ProteinNumSort.submatrix) = c(1:nrow(Protein_Class.level.data.species.ProteinNumSort.submatrix))
											
											if(!is.null(Protein_Class.level.data.species.ProteinNumSort.submatrix))
											{
												if(nrow(Protein_Class.level.data.species.ProteinNumSort.submatrix) > 0)
												{
													na.idx = which(is.na(Protein_Class.level.data.species.ProteinNumSort.submatrix), arr.ind=TRUE)
													if(nrow(na.idx) > 0)
													{
														for(ri in 1:nrow(na.idx))
														{
															Protein_Class.level.data.species.ProteinNumSort.submatrix[na.idx[ri,'row'],na.idx[ri,'col']] = ''
														}
													}
												}
											}
											
											if(!is.null(Protein_Class.level.data.species.ProteinNumSort.submatrix))
											{
												if(nrow(Protein_Class.level.data.species.ProteinNumSort.submatrix) > 0)
												{
													Title.INFO = c(class_name, class_id, class_level, class_description)
													Title.INFO = Title.INFO[which(!is.na(Title.INFO) & Title.INFO != '')]
													Title.INFO = paste(Title.INFO, collapse="|")
													
													Title = c(paste('> ', Title.INFO, sep=""), rep(NA,ncol(Protein_Class.level.data.species.ProteinNumSort.submatrix)))
													Header = c(NA, colnames(Protein_Class.level.data.species.ProteinNumSort.submatrix))
													SubMatrix = as.matrix(cbind(rownames(Protein_Class.level.data.species.ProteinNumSort.submatrix), Protein_Class.level.data.species.ProteinNumSort.submatrix))
													BlankLine = c(NA, rep(NA,ncol(Protein_Class.level.data.species.ProteinNumSort.submatrix)))
													Protein_Class.level.data.species.ProteinNumSort.submatrix.new = rbind(Title, Header, SubMatrix, BlankLine)
													
													Protein_Class.level.data.species.ProteinNumSort.submatrix.new.all <- rbind(Protein_Class.level.data.species.ProteinNumSort.submatrix.new.all, Protein_Class.level.data.species.ProteinNumSort.submatrix.new)
												}
											}
											
											
											# /// print id out for id uploading onto DB search ///
											write.table(Protein_Class.level.data.species.ProteinNumSort.submatrix[which(Protein_Class.level.data.species.ProteinNumSort.submatrix[,paste(SpeciesName.general[sp],'.UniProt',sep="")]!=''),paste(SpeciesName.general[sp],'.UniProt',sep="")], paste('./Data/',Sample,'/08_UniProt2ProteinClass/subset/n_',r,'/',folder[f],'/',SpeciesName.general[sp],'/',paste('Level',Lv,sep=""),'/id/UniProt/',nrow(Protein_Class.level.data.species.ProteinNumSort.submatrix),'_',Sample,'_exo_',SpeciesName.general[sp],'_','PC_L',Lv,'_',gsub('/','.',Protein_Class.level.data.species.ProteinNumSort[m,paste('Level.',Lv,sep="")]),'_UniProt.txt',sep=""), sep=",", quote=FALSE, row.names=FALSE, col.names=FALSE, append=FALSE)
											write.table(Protein_Class.level.data.species.ProteinNumSort.submatrix[which(Protein_Class.level.data.species.ProteinNumSort.submatrix[,paste(SpeciesName.general[sp],'.',SpeciesDB,'.ID',sep="")]!=''),paste(SpeciesName.general[sp],'.',SpeciesDB,'.ID',sep="")], paste('./Data/',Sample,'/08_UniProt2ProteinClass/subset/n_',r,'/',folder[f],'/',SpeciesName.general[sp],'/',paste('Level',Lv,sep=""),'/id/',SpeciesDB,'/',nrow(Protein_Class.level.data.species.ProteinNumSort.submatrix),'_',Sample,'_exo_',SpeciesName.general[sp],'_','PC_L',Lv,'_',gsub('/','.',Protein_Class.level.data.species.ProteinNumSort[m,paste('Level.',Lv,sep="")]),'_',SpeciesDB,'.txt',sep=""), sep=",", quote=FALSE, row.names=FALSE, col.names=FALSE, append=FALSE)
											write.table(Protein_Class.level.data.species.ProteinNumSort.submatrix[which(Protein_Class.level.data.species.ProteinNumSort.submatrix[,paste(SpeciesName.general[sp],'.Gene.ID',sep="")]!=''),paste(SpeciesName.general[sp],'.Gene.ID',sep="")], paste('./Data/',Sample,'/08_UniProt2ProteinClass/subset/n_',r,'/',folder[f],'/',SpeciesName.general[sp],'/',paste('Level',Lv,sep=""),'/id/Gene/',nrow(Protein_Class.level.data.species.ProteinNumSort.submatrix),'_',Sample,'_exo_',SpeciesName.general[sp],'_','PC_L',Lv,'_',gsub('/','.',Protein_Class.level.data.species.ProteinNumSort[m,paste('Level.',Lv,sep="")]),'_Gene.txt',sep=""), sep=",", quote=FALSE, row.names=FALSE, col.names=FALSE, append=FALSE)
											
											# /// print name out for view ///
											write.table(mixedsort(Protein_Class.level.data.species.ProteinNumSort.submatrix[which(Protein_Class.level.data.species.ProteinNumSort.submatrix[,paste(SpeciesName.general[sp],'.Protein.names',sep="")]!=''),paste(SpeciesName.general[sp],'.Protein.names',sep="")]), paste('./Data/',Sample,'/08_UniProt2ProteinClass/subset/n_',r,'/',folder[f],'/',SpeciesName.general[sp],'/',paste('Level',Lv,sep=""),'/name/Protein/',nrow(Protein_Class.level.data.species.ProteinNumSort.submatrix),'_',Sample,'_exo_',SpeciesName.general[sp],'_','PC_L',Lv,'_',gsub('/','.',Protein_Class.level.data.species.ProteinNumSort[m,paste('Level.',Lv,sep="")]),'_ProteinName.txt',sep=""), sep=",", quote=FALSE, row.names=FALSE, col.names=FALSE, append=FALSE)
										}
									}
								}
								
								
								# --- write xlsx file by openxlsx package -----------------------------------------------------------------------------
								# Importing a big xlsx file into R? https://stackoverflow.com/questions/19147884/importing-a-big-xlsx-file-into-r/43118530#43118530
								# Building R for Windows: https://cran.r-project.org/bin/windows/Rtools/
								#Sys.setenv("R_ZIPCMD" = "C:/Rtools/bin/zip.exe")	# path to zip.exe
								
								for(s in 1:3)
								{
									if(s == 1) {
										prename = paste('./Data/',Sample,'/08_UniProt2ProteinClass/subset/n_',r,'/',folder[f],'/',SpeciesName.general[sp],'/',paste('Level',Lv,sep=""),'/',Sample,'_exo_',SpeciesName.general[sp],'_','ProteinClass_Level',Lv,sep="")
									} else if(s == 2) {
										prename = paste('./Data/',Sample,'/08_UniProt2ProteinClass/subset/n_',r,'/',folder[f],'/',SpeciesName.general[sp],'/',paste('Level',Lv,sep=""),'/',Sample,'_exo_',SpeciesName.general[sp],'_','ProteinClass_Level',Lv,'_LevelNameSort',sep="")
									} else if(s == 3) {
										prename = paste('./Data/',Sample,'/08_UniProt2ProteinClass/subset/n_',r,'/',folder[f],'/',SpeciesName.general[sp],'/',paste('Level',Lv,sep=""),'/',Sample,'_exo_',SpeciesName.general[sp],'_','ProteinClass_Level',Lv,'_ProteinNumSort',sep="")
									}
									
									wb <- openxlsx::createWorkbook(paste(prename,'.xlsx',sep=""))
									modifyBaseFont(wb, fontSize=12, fontColour="black", fontName="Times New Roman")
									
									# --- sheet=1 -------------------------------------------------------------------------------------
									if(s == 1) {
										sheet = paste('Level','.',Lv,sep="")
										sheetData = Protein_Class.level.data.species
									} else if(s == 2) {
										sheet = paste('Level','.',Lv,sep="")
										sheetData = Protein_Class.level.data.species.LevelNameSort
									} else if(s == 3) {
										sheet = paste('Level','.',Lv,sep="")
										sheetData = Protein_Class.level.data.species.ProteinNumSort
									}
									
									if(!is.null(sheetData))
									{
										if(nrow(sheetData) > 0)
										{
											addWorksheet(wb, sheet)
											writeData(wb, sheet, sheetData, colNames=TRUE, rowNames=FALSE, keepNA=FALSE)
											
											setColWidths(wb, sheet, cols=c(1:ncol(sheetData)), widths=c(25,rep(20,ncol(sheetData)-1)))
											setRowHeights(wb, sheet, rows=c(1:(1+nrow(sheetData))), heights=16.5)
											freezePane(wb, sheet, firstActiveRow=2, firstActiveCol=2)
											
											style <- createStyle(halign="left", valign="center")
											addStyle(wb, sheet, style, rows=c(1:(1+nrow(sheetData))), cols=c(1:ncol(sheetData)), gridExpand=TRUE)
											style <- createStyle(textDecoration="bold")
											addStyle(wb, sheet, style, rows=1, cols=c(1:ncol(sheetData)), stack=TRUE)
											
											right_align.idx = which(colnames(sheetData) %in% c(paste(SpeciesName.general[sp],'.Protein.number',sep=""),paste(SpeciesName.general[sp],'.Protein.percent',sep="")))
											if(length(right_align.idx) > 0) {
												style <- createStyle(halign="right", valign="center")
												addStyle(wb, sheet, style, rows=c(2:(1+nrow(sheetData))), cols=right_align.idx, gridExpand=TRUE, stack=TRUE)
											}
											
											# -----------------------------------------------------------------------------------
											
											addWorksheet(eval(parse(text=paste('wb',s,'1',sep=""))), sheet)
											writeData(eval(parse(text=paste('wb',s,'1',sep=""))), sheet, sheetData, colNames=TRUE, rowNames=FALSE, keepNA=FALSE)
											
											setColWidths(eval(parse(text=paste('wb',s,'1',sep=""))), sheet, cols=c(1:ncol(sheetData)), widths=c(25,rep(20,ncol(sheetData)-1)))
											setRowHeights(eval(parse(text=paste('wb',s,'1',sep=""))), sheet, rows=c(1:(1+nrow(sheetData))), heights=16.5)
											freezePane(eval(parse(text=paste('wb',s,'1',sep=""))), sheet, firstActiveRow=2, firstActiveCol=2)
											
											style <- createStyle(halign="left", valign="center")
											addStyle(eval(parse(text=paste('wb',s,'1',sep=""))), sheet, style, rows=c(1:(1+nrow(sheetData))), cols=c(1:ncol(sheetData)), gridExpand=TRUE)
											style <- createStyle(textDecoration="bold")
											addStyle(eval(parse(text=paste('wb',s,'1',sep=""))), sheet, style, rows=1, cols=c(1:ncol(sheetData)), stack=TRUE)
											
											right_align.idx = which(colnames(sheetData) %in% c(paste(SpeciesName.general[sp],'.Protein.number',sep=""),paste(SpeciesName.general[sp],'.Protein.percent',sep="")))
											if(length(right_align.idx) > 0) {
												style <- createStyle(halign="right", valign="center")
												addStyle(eval(parse(text=paste('wb',s,'1',sep=""))), sheet, style, rows=c(2:(1+nrow(sheetData))), cols=right_align.idx, gridExpand=TRUE, stack=TRUE)
											}
											
											if(r == 1 & folder[f] == 'high')
											{
												addWorksheet(eval(parse(text=paste('wb',s,'2',sep=""))), sheet)
												writeData(eval(parse(text=paste('wb',s,'2',sep=""))), sheet, sheetData, colNames=TRUE, rowNames=FALSE, keepNA=FALSE)
												
												setColWidths(eval(parse(text=paste('wb',s,'2',sep=""))), sheet, cols=c(1:ncol(sheetData)), widths=c(25,rep(20,ncol(sheetData)-1)))
												setRowHeights(eval(parse(text=paste('wb',s,'2',sep=""))), sheet, rows=c(1:(1+nrow(sheetData))), heights=16.5)
												freezePane(eval(parse(text=paste('wb',s,'2',sep=""))), sheet, firstActiveRow=2, firstActiveCol=2)
												
												style <- createStyle(halign="left", valign="center")
												addStyle(eval(parse(text=paste('wb',s,'2',sep=""))), sheet, style, rows=c(1:(1+nrow(sheetData))), cols=c(1:ncol(sheetData)), gridExpand=TRUE)
												style <- createStyle(textDecoration="bold")
												addStyle(eval(parse(text=paste('wb',s,'2',sep=""))), sheet, style, rows=1, cols=c(1:ncol(sheetData)), stack=TRUE)
												
												right_align.idx = which(colnames(sheetData) %in% c(paste(SpeciesName.general[sp],'.Protein.number',sep=""),paste(SpeciesName.general[sp],'.Protein.percent',sep="")))
												if(length(right_align.idx) > 0) {
													style <- createStyle(halign="right", valign="center")
													addStyle(eval(parse(text=paste('wb',s,'2',sep=""))), sheet, style, rows=c(2:(1+nrow(sheetData))), cols=right_align.idx, gridExpand=TRUE, stack=TRUE)
												}
											}
										}
									}
									
									# --- sheet=2 -------------------------------------------------------------------------------------
									if(s == 1) {
										sheet = paste('Level','.',Lv,'_','ProteinClass2Protein',sep="")
										sheetData = Protein_Class.level.data.species.submatrix.new.all
									} else if(s == 2) {
										sheet = paste('Level','.',Lv,'_','ProteinClass2Protein',sep="")
										sheetData = Protein_Class.level.data.species.LevelNameSort.submatrix.new.all
									} else if(s == 3) {
										sheet = paste('Level','.',Lv,'_','ProteinClass2Protein',sep="")
										sheetData = Protein_Class.level.data.species.ProteinNumSort.submatrix.new.all
									}
									
									if(!is.null(sheetData))
									{
										if(nrow(sheetData) > 0)
										{
											addWorksheet(wb, sheet)
											writeData(wb, sheet, sheetData, colNames=FALSE, rowNames=FALSE, keepNA=FALSE)
											
											setColWidths(wb, sheet, cols=c(1:ncol(sheetData)), widths=c(8.5,rep(15,3),rep(20,2),rep(12.5,ncol(sheetData)-(1+3+2))))
											setRowHeights(wb, sheet, rows=c(1:nrow(sheetData)), heights=16.5)
											freezePane(wb, sheet, firstActiveCol=7)
											
											style <- createStyle(halign="left", valign="center")
											addStyle(wb, sheet, style, rows=c(1:nrow(sheetData)), cols=c(1:ncol(sheetData)), gridExpand=TRUE)
											
											# -----------------------------------------------------------------------------------
											
											addWorksheet(eval(parse(text=paste('wb',s,'1',sep=""))), sheet)
											writeData(eval(parse(text=paste('wb',s,'1',sep=""))), sheet, sheetData, colNames=FALSE, rowNames=FALSE, keepNA=FALSE)
											
											setColWidths(eval(parse(text=paste('wb',s,'1',sep=""))), sheet, cols=c(1:ncol(sheetData)), widths=c(8.5,rep(15,3),rep(20,2),rep(12.5,ncol(sheetData)-(1+3+2))))
											setRowHeights(eval(parse(text=paste('wb',s,'1',sep=""))), sheet, rows=c(1:nrow(sheetData)), heights=16.5)
											freezePane(eval(parse(text=paste('wb',s,'1',sep=""))), sheet, firstActiveCol=7)
											
											style <- createStyle(halign="left", valign="center")
											addStyle(eval(parse(text=paste('wb',s,'1',sep=""))), sheet, style, rows=c(1:nrow(sheetData)), cols=c(1:ncol(sheetData)), gridExpand=TRUE)
											
											if(r == 1 & folder[f] == 'high')
											{
												addWorksheet(eval(parse(text=paste('wb',s,'2',sep=""))), sheet)
												writeData(eval(parse(text=paste('wb',s,'2',sep=""))), sheet, sheetData, colNames=FALSE, rowNames=FALSE, keepNA=FALSE)
												
												setColWidths(eval(parse(text=paste('wb',s,'2',sep=""))), sheet, cols=c(1:ncol(sheetData)), widths=c(8.5,rep(15,3),rep(20,2),rep(12.5,ncol(sheetData)-(1+3+2))))
												setRowHeights(eval(parse(text=paste('wb',s,'2',sep=""))), sheet, rows=c(1:nrow(sheetData)), heights=16.5)
												freezePane(eval(parse(text=paste('wb',s,'2',sep=""))), sheet, firstActiveCol=7)
												
												style <- createStyle(halign="left", valign="center")
												addStyle(eval(parse(text=paste('wb',s,'2',sep=""))), sheet, style, rows=c(1:nrow(sheetData)), cols=c(1:ncol(sheetData)), gridExpand=TRUE)
											}
										}
									}
									
									if(length(wb$sheet_names) > 0) {
										openxlsx::saveWorkbook(wb, file=paste(prename,'.xlsx',sep=""), overwrite=TRUE)
									}
								}
								
								
								# --- plotting ------------------------------------------------------------------------------------------------------
								number.cutoff = 5
								percent.cutoff = 5
								
								df = data.frame(group=Protein_Class.level.data.species.ProteinNumSort[,paste('Level.',Lv,sep="")], value=as.numeric(as.character(Protein_Class.level.data.species.ProteinNumSort[,paste(SpeciesName.general[sp],'.Protein.number',sep="")])), percent=as.numeric(as.character(Protein_Class.level.data.species.ProteinNumSort[,paste(SpeciesName.general[sp],'.Protein.percent',sep="")])))
								#df$group = factor(df$group, levels=df$group)
								df$value = as.numeric(as.character(df$value))
								df$percent = as.numeric(as.character(df$percent))
								
								if(Lv == 0 | Lv == 1) {
									df.new = df
								} else {
									if(Lv == 2) {
										number.cutoff = 5
										percent.cutoff = 1
									} else if(Lv == 3) {
										number.cutoff = 5
										percent.cutoff = 5
									} else if(Lv == 4 | Lv == 5) {
										number.cutoff = 1
										percent.cutoff = 5
									}
									
									if(length(which(df$value < number.cutoff)) > 0) {
										df.new = rbind(df[which(df$value >= number.cutoff),,drop=FALSE], data.frame(cbind(group='Others', value=sum(df[which(df$value < number.cutoff),'value']), percent=sum(df[which(df$value < number.cutoff),'percent']))))
										if(df.new[which(df.new[,'group']=='Others'),'percent']==0 | df.new[which(df.new[,'group']=='Others'),'percent']=='NaN')	df.new = df.new[-which(df.new[,'group']=='Others'),,drop=FALSE]
									} else {
										df.new = df
									}
								}
								
								if(nrow(df.new) != 0)
								{
									df.new$group = factor(df.new$group, levels=df.new$group)
									df.new$value = as.numeric(as.character(df.new$value))
									df.new$percent = as.numeric(as.character(df.new$percent))
									
									if(sum(df.new$percent) != 0)
									{
										#pic <- ggplot(df.new, aes(x="", y=percent, fill=group)) + geom_bar(width=1, stat="identity")
										#pic <- pic + geom_col(position=position_stack(reverse = TRUE), width=1)
										
										#if(Lv == 0) {
										#	pic <- ggplot(df.new, aes(x=group, y=percent, fill=group)) + geom_bar(width=df.new$percent/100, stat="identity")
										#} else if(Lv == 1) {
										#	pic <- ggplot(df.new, aes(x=group, y=percent, fill=group)) + geom_bar(width=df.new$percent/15, stat="identity")
										#} else if(Lv == 2) {
										#	pic <- ggplot(df.new, aes(x=group, y=percent, fill=group)) + geom_bar(width=df.new$percent/20, stat="identity")
										#} else if(Lv == 3) {
										#	pic <- ggplot(df.new, aes(x=group, y=percent, fill=group)) + geom_bar(width=df.new$percent/25, stat="identity")
										#} else if(Lv == 4) {
										#	pic <- ggplot(df.new, aes(x=group, y=percent, fill=group)) + geom_bar(width=df.new$percent/30, stat="identity")
										#}
										pic <- ggplot(df.new, aes(x=group, y=percent, fill=group)) + geom_bar(width=df.new$percent/max(df.new$percent), stat="identity")
										
										pic <- pic + coord_flip()
										pic <- pic + scale_x_discrete(limits=rev(df.new$group))
										
										#pic <- pic + scale_fill_grey()
										pic <- pic + scale_fill_grey(start=0.2, end=0.6)
										
										pic <- pic + theme_bw()
										pic <- pic + theme(panel.grid.major = element_blank())
										pic <- pic + theme(panel.grid.minor = element_blank())
										pic <- pic + theme(panel.border = element_blank())
										
										pic <- pic + guides(fill=guide_legend(title='protein class'))
										pic <- pic + theme(legend.title = element_text(color = "black", size = 16))
										pic <- pic + theme(legend.text = element_text(color = "black", size = 14))
										pic <- pic + theme(legend.position="none")
										
										x.pos <- rep(1.0,length(df.new$percent))
										#x.pos <- rep(1.0,length(which(df.new$percent >= percent.cutoff)))
										#if(length(which(df.new$percent < percent.cutoff)) > 0)
										#{
										#	x.pos.shift_dist = 0.025
										#	
										#	for(k in 1:length(which(df.new$percent < percent.cutoff)))
										#	{
										#		x.pos.shift_pos = 1.0 + (x.pos.shift_dist * k)
										#		
										#		#if(k%%2 == 0) {
										#		#	x.pos.shift_pos = 1.0 + (x.pos.shift_dist * 2)
										#		#} else if(k%%2 == 1) {
										#		#	x.pos.shift_pos = 1.0 + (x.pos.shift_dist * 1)
										#		#}
										#		
										#		x.pos <- c(x.pos, x.pos.shift_pos)
										#	}
										#}
										#x.pos
										
										df.new.breaks <- df.new$percent[1]/2
										if(nrow(df.new) >= 2) {
											for(k in 2:nrow(df.new))
											{
												df.new.breaks <- c(df.new.breaks, sum(df.new$percent[1:(k-1)]) + df.new$percent[k]/2)
											}
										}
										#df.new.breaks
										
										#pic <- pic + geom_text(aes(x=x.pos, y=df.new.breaks, label=paste(formatC(as.numeric(as.character(round_preserve_sum(df.new$percent,2))),format='f',digits=2),'%',sep='')), size=2.5, color="white")
										
										if(Lv == 0 | Lv == 5) {
											y.pos.shift_dist = 5
										} else if(Lv == 1) {
											y.pos.shift_dist = (0.8+0.95)/2
											if(Sample=='C1')			y.pos.shift_dist = 0.8
											if(Sample=='E8')			y.pos.shift_dist = 0.95
											if(Sample=='C1_8R')			y.pos.shift_dist = 0.85
											if(Sample=='C1_ALL')		y.pos.shift_dist = 0.85
											if(Sample=='pcMSC_N230')	y.pos.shift_dist = 0.95
										} else if(Lv == 2) {
											y.pos.shift_dist = (1.5+3.0)/2
											if(Sample=='C1')			y.pos.shift_dist = 3.0
											if(Sample=='E8')			y.pos.shift_dist = 3.0
											if(Sample=='C1_8R')			y.pos.shift_dist = 1.5
											if(Sample=='C1_ALL')		y.pos.shift_dist = 1.5
											if(Sample=='pcMSC_N230')	y.pos.shift_dist = 3.0
										} else if(Lv == 3) {
											y.pos.shift_dist = (1.7+3.8)/2
											if(Sample=='C1')			y.pos.shift_dist = 3.4
											if(Sample=='E8')			y.pos.shift_dist = 3.8
											if(Sample=='C1_8R')			y.pos.shift_dist = 1.7
											if(Sample=='C1_ALL')		y.pos.shift_dist = 1.7
											if(Sample=='pcMSC_N230')	y.pos.shift_dist = 3.8
										} else if(Lv == 4) {
											y.pos.shift_dist = (1.8+2.6)/2
											if(Sample=='C1')			y.pos.shift_dist = 2.6
											if(Sample=='E8')			y.pos.shift_dist = 2.2
											if(Sample=='C1_8R')			y.pos.shift_dist = 1.8
											if(Sample=='C1_ALL')		y.pos.shift_dist = 1.8
											if(Sample=='pcMSC_N230')	y.pos.shift_dist = 2.6
										}
										pic <- pic + geom_text(aes(y=(df.new$percent)+y.pos.shift_dist, label=paste(formatC(as.numeric(as.character(round_preserve_sum(df.new$percent,1))),format='f',digits=1),'%',sep='')), size=5, color="black")
										
										#pic <- pic + scale_y_continuous(breaks=df.new.breaks, labels=df.new$group, position="right")
										#pic <- pic + scale_y_reverse(breaks=df.new.breaks, labels=df.new$group, position="right")
										
										pic <- pic + theme(axis.text.x = element_text(color = "black", size = 14, hjust = 1))
										pic <- pic + theme(axis.text.y = element_text(color = "black", size = 14, vjust = 0.3))
										
										#pic <- pic + coord_polar(theta="y", start=0)
										
										pic <- pic + theme(axis.title.x = element_blank())
										pic <- pic + theme(axis.title.y = element_blank())
										pic <- pic + theme(axis.text.x = element_blank())
										#pic <- pic + theme(axis.text.y = element_blank())
										pic <- pic + theme(axis.ticks.x = element_blank())
										pic <- pic + theme(axis.ticks.y = element_blank())
										
										pic <- pic + theme(panel.background = element_rect(fill="transparent", colour=NA))
										pic <- pic + theme(plot.background = element_rect(fill="transparent", colour=NA))
										pic <- pic + theme(legend.background = element_rect(fill="transparent", colour=NA))
										pic <- pic + theme(legend.box.background = element_rect(fill="transparent", colour=NA))
										pic <- pic + theme(legend.key = element_rect(fill="transparent", colour=NA))
										
										#pic
										
										
										if(nrow(df.new) == 1) {
											ggsave(paste('./Data/',Sample,'/08_UniProt2ProteinClass/subset/n_',r,'/',folder[f],'/',SpeciesName.general[sp],'/',paste('Level',Lv,sep=""),'/',Sample,'_exo_',SpeciesName.general[sp],'_','ProteinClass_Level',Lv,'_ProteinNumSort.png',sep=""), pic, width=12, height=3/5, units="in", dpi=600, bg="transparent")
											ggsave(paste('./Data/',Sample,'/08_UniProt2ProteinClass/subset/n_',r,'/',folder[f],'/',SpeciesName.general[sp],'/',Sample,'_exo_',SpeciesName.general[sp],'_','ProteinClass_Level',Lv,'_ProteinNumSort.png',sep=""), pic, width=12, height=3/5, units="in", dpi=600, bg="transparent")
											if(r == 1 & folder[f] == 'high')	ggsave(paste('./Data/',Sample,'/08_UniProt2ProteinClass/',Sample,'_exo_',SpeciesName.general[sp],'_','ProteinClass_Level',Lv,'_ProteinNumSort.png',sep=""), pic, width=12, height=3/5, units="in", dpi=600, bg="transparent")
										} else if(nrow(df.new) > 1 & nrow(df.new) <= 5) {
											ggsave(paste('./Data/',Sample,'/08_UniProt2ProteinClass/subset/n_',r,'/',folder[f],'/',SpeciesName.general[sp],'/',paste('Level',Lv,sep=""),'/',Sample,'_exo_',SpeciesName.general[sp],'_','ProteinClass_Level',Lv,'_ProteinNumSort.png',sep=""), pic, width=12, height=5/5, units="in", dpi=600, bg="transparent")
											ggsave(paste('./Data/',Sample,'/08_UniProt2ProteinClass/subset/n_',r,'/',folder[f],'/',SpeciesName.general[sp],'/',Sample,'_exo_',SpeciesName.general[sp],'_','ProteinClass_Level',Lv,'_ProteinNumSort.png',sep=""), pic, width=12, height=5/5, units="in", dpi=600, bg="transparent")
											if(r == 1 & folder[f] == 'high')	ggsave(paste('./Data/',Sample,'/08_UniProt2ProteinClass/',Sample,'_exo_',SpeciesName.general[sp],'_','ProteinClass_Level',Lv,'_ProteinNumSort.png',sep=""), pic, width=12, height=5/5, units="in", dpi=600, bg="transparent")
										} else if(nrow(df.new) > 5 & nrow(df.new) <= 10) {
											ggsave(paste('./Data/',Sample,'/08_UniProt2ProteinClass/subset/n_',r,'/',folder[f],'/',SpeciesName.general[sp],'/',paste('Level',Lv,sep=""),'/',Sample,'_exo_',SpeciesName.general[sp],'_','ProteinClass_Level',Lv,'_ProteinNumSort.png',sep=""), pic, width=12, height=10/5, units="in", dpi=600, bg="transparent")
											ggsave(paste('./Data/',Sample,'/08_UniProt2ProteinClass/subset/n_',r,'/',folder[f],'/',SpeciesName.general[sp],'/',Sample,'_exo_',SpeciesName.general[sp],'_','ProteinClass_Level',Lv,'_ProteinNumSort.png',sep=""), pic, width=12, height=10/5, units="in", dpi=600, bg="transparent")
											if(r == 1 & folder[f] == 'high')	ggsave(paste('./Data/',Sample,'/08_UniProt2ProteinClass/',Sample,'_exo_',SpeciesName.general[sp],'_','ProteinClass_Level',Lv,'_ProteinNumSort.png',sep=""), pic, width=12, height=10/5, units="in", dpi=600, bg="transparent")
										} else if(nrow(df.new) > 10 & nrow(df.new) <= 15) {
											ggsave(paste('./Data/',Sample,'/08_UniProt2ProteinClass/subset/n_',r,'/',folder[f],'/',SpeciesName.general[sp],'/',paste('Level',Lv,sep=""),'/',Sample,'_exo_',SpeciesName.general[sp],'_','ProteinClass_Level',Lv,'_ProteinNumSort.png',sep=""), pic, width=12, height=15/5, units="in", dpi=600, bg="transparent")
											ggsave(paste('./Data/',Sample,'/08_UniProt2ProteinClass/subset/n_',r,'/',folder[f],'/',SpeciesName.general[sp],'/',Sample,'_exo_',SpeciesName.general[sp],'_','ProteinClass_Level',Lv,'_ProteinNumSort.png',sep=""), pic, width=12, height=15/5, units="in", dpi=600, bg="transparent")
											if(r == 1 & folder[f] == 'high')	ggsave(paste('./Data/',Sample,'/08_UniProt2ProteinClass/',Sample,'_exo_',SpeciesName.general[sp],'_','ProteinClass_Level',Lv,'_ProteinNumSort.png',sep=""), pic, width=12, height=15/5, units="in", dpi=600, bg="transparent")
										} else if(nrow(df.new) > 15 & nrow(df.new) <= 20) {
											ggsave(paste('./Data/',Sample,'/08_UniProt2ProteinClass/subset/n_',r,'/',folder[f],'/',SpeciesName.general[sp],'/',paste('Level',Lv,sep=""),'/',Sample,'_exo_',SpeciesName.general[sp],'_','ProteinClass_Level',Lv,'_ProteinNumSort.png',sep=""), pic, width=12, height=20/5, units="in", dpi=600, bg="transparent")
											ggsave(paste('./Data/',Sample,'/08_UniProt2ProteinClass/subset/n_',r,'/',folder[f],'/',SpeciesName.general[sp],'/',Sample,'_exo_',SpeciesName.general[sp],'_','ProteinClass_Level',Lv,'_ProteinNumSort.png',sep=""), pic, width=12, height=20/5, units="in", dpi=600, bg="transparent")
											if(r == 1 & folder[f] == 'high')	ggsave(paste('./Data/',Sample,'/08_UniProt2ProteinClass/',Sample,'_exo_',SpeciesName.general[sp],'_','ProteinClass_Level',Lv,'_ProteinNumSort.png',sep=""), pic, width=12, height=20/5, units="in", dpi=600, bg="transparent")
										} else if(nrow(df.new) > 20 & nrow(df.new) <= 25) {
											ggsave(paste('./Data/',Sample,'/08_UniProt2ProteinClass/subset/n_',r,'/',folder[f],'/',SpeciesName.general[sp],'/',paste('Level',Lv,sep=""),'/',Sample,'_exo_',SpeciesName.general[sp],'_','ProteinClass_Level',Lv,'_ProteinNumSort.png',sep=""), pic, width=12, height=25/5, units="in", dpi=600, bg="transparent")
											ggsave(paste('./Data/',Sample,'/08_UniProt2ProteinClass/subset/n_',r,'/',folder[f],'/',SpeciesName.general[sp],'/',Sample,'_exo_',SpeciesName.general[sp],'_','ProteinClass_Level',Lv,'_ProteinNumSort.png',sep=""), pic, width=12, height=25/5, units="in", dpi=600, bg="transparent")
											if(r == 1 & folder[f] == 'high')	ggsave(paste('./Data/',Sample,'/08_UniProt2ProteinClass/',Sample,'_exo_',SpeciesName.general[sp],'_','ProteinClass_Level',Lv,'_ProteinNumSort.png',sep=""), pic, width=12, height=25/5, units="in", dpi=600, bg="transparent")
										} else if(nrow(df.new) > 25 & nrow(df.new) <= 30) {
											ggsave(paste('./Data/',Sample,'/08_UniProt2ProteinClass/subset/n_',r,'/',folder[f],'/',SpeciesName.general[sp],'/',paste('Level',Lv,sep=""),'/',Sample,'_exo_',SpeciesName.general[sp],'_','ProteinClass_Level',Lv,'_ProteinNumSort.png',sep=""), pic, width=12, height=30/5, units="in", dpi=600, bg="transparent")
											ggsave(paste('./Data/',Sample,'/08_UniProt2ProteinClass/subset/n_',r,'/',folder[f],'/',SpeciesName.general[sp],'/',Sample,'_exo_',SpeciesName.general[sp],'_','ProteinClass_Level',Lv,'_ProteinNumSort.png',sep=""), pic, width=12, height=30/5, units="in", dpi=600, bg="transparent")
											if(r == 1 & folder[f] == 'high')	ggsave(paste('./Data/',Sample,'/08_UniProt2ProteinClass/',Sample,'_exo_',SpeciesName.general[sp],'_','ProteinClass_Level',Lv,'_ProteinNumSort.png',sep=""), pic, width=12, height=30/5, units="in", dpi=600, bg="transparent")
										} else if(nrow(df.new) > 30 & nrow(df.new) <= 35) {
											ggsave(paste('./Data/',Sample,'/08_UniProt2ProteinClass/subset/n_',r,'/',folder[f],'/',SpeciesName.general[sp],'/',paste('Level',Lv,sep=""),'/',Sample,'_exo_',SpeciesName.general[sp],'_','ProteinClass_Level',Lv,'_ProteinNumSort.png',sep=""), pic, width=12, height=35/5, units="in", dpi=600, bg="transparent")
											ggsave(paste('./Data/',Sample,'/08_UniProt2ProteinClass/subset/n_',r,'/',folder[f],'/',SpeciesName.general[sp],'/',Sample,'_exo_',SpeciesName.general[sp],'_','ProteinClass_Level',Lv,'_ProteinNumSort.png',sep=""), pic, width=12, height=35/5, units="in", dpi=600, bg="transparent")
											if(r == 1 & folder[f] == 'high')	ggsave(paste('./Data/',Sample,'/08_UniProt2ProteinClass/',Sample,'_exo_',SpeciesName.general[sp],'_','ProteinClass_Level',Lv,'_ProteinNumSort.png',sep=""), pic, width=12, height=35/5, units="in", dpi=600, bg="transparent")
										} else {
											ggsave(paste('./Data/',Sample,'/08_UniProt2ProteinClass/subset/n_',r,'/',folder[f],'/',SpeciesName.general[sp],'/',paste('Level',Lv,sep=""),'/',Sample,'_exo_',SpeciesName.general[sp],'_','ProteinClass_Level',Lv,'_ProteinNumSort.png',sep=""), pic, width=12, height=nrow(df.new)/5, units="in", dpi=600, bg="transparent")
											ggsave(paste('./Data/',Sample,'/08_UniProt2ProteinClass/subset/n_',r,'/',folder[f],'/',SpeciesName.general[sp],'/',Sample,'_exo_',SpeciesName.general[sp],'_','ProteinClass_Level',Lv,'_ProteinNumSort.png',sep=""), pic, width=12, height=nrow(df.new)/5, units="in", dpi=600, bg="transparent")
											if(r == 1 & folder[f] == 'high')	ggsave(paste('./Data/',Sample,'/08_UniProt2ProteinClass/',Sample,'_exo_',SpeciesName.general[sp],'_','ProteinClass_Level',Lv,'_ProteinNumSort.png',sep=""), pic, width=12, height=nrow(df.new)/5, units="in", dpi=600, bg="transparent")
										}
									}
								}
							}
							
							
							for(s in 1:3)
							{
								if(length(eval(parse(text=paste('wb',s,'1',sep="")))$sheet_names) > 0) {
									openxlsx::saveWorkbook(eval(parse(text=paste('wb',s,'1',sep=""))), file=paste(eval(parse(text=paste('prename',s,'1',sep=""))),'.xlsx',sep=""), overwrite=TRUE)
								}
								
								if(r == 1 & folder[f] == 'high')
								{
									if(length(eval(parse(text=paste('wb',s,'2',sep="")))$sheet_names) > 0) {
										openxlsx::saveWorkbook(eval(parse(text=paste('wb',s,'2',sep=""))), file=paste(eval(parse(text=paste('prename',s,'2',sep=""))),'.xlsx',sep=""), overwrite=TRUE)
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

