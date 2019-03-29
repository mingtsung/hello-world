rm(list=ls())	# remove all variables in workspace

setwd("D:/Exo")

Sample <- "pcMSC_N230"

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


# === input All proteins recorded in Microvesicle database ======================================================================
#Database <- c('ExoCarta','Vesiclepedia','EVpedia')

# --- ExoCarta --------------------------------------------------------------------------------------------------------
ExoCarta_Protein_mRNA <- read.table(paste('./Database/Microvesicle/ExoCarta/EXOCARTA_PROTEIN_MRNA_DETAILS_5.txt',sep=""), sep="\t", quote="", header=TRUE, stringsAsFactors=FALSE)

# https://www.ncbi.nlm.nih.gov/gene/81822 [Rattus norvegicus (Norway rat)] -> EXPERIMENT ID: 34 -> PUBMED ID: 19190083 -> https://www.ncbi.nlm.nih.gov/pubmed/19190083 -> Characterization of exosome-like vesicles released from human tracheobronchial ciliated epithelium: a possible role in innate defense -> human -> Homo sapiens
ExoCarta_Protein_mRNA[which(trimws(ExoCarta_Protein_mRNA[,'SPECIES']) == 'Homo sapiens' & (trimws(ExoCarta_Protein_mRNA[,'ENTREZ.GENE.ID']) == '11461' | trimws(ExoCarta_Protein_mRNA[,'ENTREZ.GENE.ID']) == '81822' | trimws(ExoCarta_Protein_mRNA[,'ENTREZ.GENE.ID']) == '280979' | trimws(ExoCarta_Protein_mRNA[,'ENTREZ.GENE.ID']) == '443052')),'GENE.SYMBOL'] = 'ACTB'
ExoCarta_Protein_mRNA[which(trimws(ExoCarta_Protein_mRNA[,'SPECIES']) == 'Homo sapiens' & (trimws(ExoCarta_Protein_mRNA[,'ENTREZ.GENE.ID']) == '11461' | trimws(ExoCarta_Protein_mRNA[,'ENTREZ.GENE.ID']) == '81822' | trimws(ExoCarta_Protein_mRNA[,'ENTREZ.GENE.ID']) == '280979' | trimws(ExoCarta_Protein_mRNA[,'ENTREZ.GENE.ID']) == '443052')),'ENTREZ.GENE.ID'] = '60'

# https://www.ncbi.nlm.nih.gov/gene/287022 [Bos taurus (cattle)] -> EXPERIMENT ID: 34 -> PUBMED ID: 19190083 -> https://www.ncbi.nlm.nih.gov/pubmed/19190083 -> Characterization of exosome-like vesicles released from human tracheobronchial ciliated epithelium: a possible role in innate defense -> human -> Homo sapiens
ExoCarta_Protein_mRNA[which(trimws(ExoCarta_Protein_mRNA[,'SPECIES']) == 'Homo sapiens' & trimws(ExoCarta_Protein_mRNA[,'ENTREZ.GENE.ID']) == '287022'),'GENE.SYMBOL'] = 'YWHAZ'
ExoCarta_Protein_mRNA[which(trimws(ExoCarta_Protein_mRNA[,'SPECIES']) == 'Homo sapiens' & trimws(ExoCarta_Protein_mRNA[,'ENTREZ.GENE.ID']) == '287022'),'ENTREZ.GENE.ID'] = '7534'

# EXPERIMENT ID: 79,80,81 -> PUBMED ID: 20458337 -> https://www.ncbi.nlm.nih.gov/pubmed/20458337 -> MHC class II-associated proteins in B-cell exosomes and potential functional implications for exosome biogenesis -> Supplementary Table 1a (icb201064x3.xls) -> entrez number:29823 -> protein description:MYL6B;MYL6 Isoform Non-muscle of Myosin light polypeptide 6 ; gene symbol:MYL6B -> https://www.ncbi.nlm.nih.gov/gene/140465 -> GeneID:140465 ; Official Symbol:MYL6B
ExoCarta_Protein_mRNA[which(trimws(ExoCarta_Protein_mRNA[,'SPECIES']) == 'Homo sapiens' & trimws(ExoCarta_Protein_mRNA[,'ENTREZ.GENE.ID']) == '29823'),'GENE.SYMBOL'] = 'MYL6B'
ExoCarta_Protein_mRNA[which(trimws(ExoCarta_Protein_mRNA[,'SPECIES']) == 'Homo sapiens' & trimws(ExoCarta_Protein_mRNA[,'ENTREZ.GENE.ID']) == '29823'),'ENTREZ.GENE.ID'] = '140465'

# https://www.ncbi.nlm.nih.gov/gene/14252 [Mus musculus (house mouse)] -> EXPERIMENT ID: 65 -> PUBMED ID: 19415654 -> https://www.ncbi.nlm.nih.gov/pubmed/19415654 -> Proteomics of MUC1-containing lipid rafts from plasma membranes and exosomes of human breast carcinoma cells MCF-7 -> human -> Homo sapiens
ExoCarta_Protein_mRNA[which(trimws(ExoCarta_Protein_mRNA[,'SPECIES']) == 'Homo sapiens' & trimws(ExoCarta_Protein_mRNA[,'ENTREZ.GENE.ID']) == '14252' & trimws(ExoCarta_Protein_mRNA[,'EXPERIMENT.ID']) == '65'),'GENE.SYMBOL'] = 'FLOT2'
ExoCarta_Protein_mRNA[which(trimws(ExoCarta_Protein_mRNA[,'SPECIES']) == 'Homo sapiens' & trimws(ExoCarta_Protein_mRNA[,'ENTREZ.GENE.ID']) == '14252' & trimws(ExoCarta_Protein_mRNA[,'EXPERIMENT.ID']) == '65'),'ENTREZ.GENE.ID'] = '2319'

# EXPERIMENT ID: 130 -> PUBMED ID: 18309083 -> https://www.ncbi.nlm.nih.gov/pubmed/18309083 -> Ceramide triggers budding of exosome vesicles into multivesicular endosomes -> mouse -> Mus musculus
ExoCarta_Protein_mRNA[which(trimws(ExoCarta_Protein_mRNA[,'SPECIES']) == 'Homo sapiens' & trimws(ExoCarta_Protein_mRNA[,'EXPERIMENT.ID']) == '130'),'SPECIES'] = 'Mus musculus'

# https://www.ncbi.nlm.nih.gov/gene/301252 [Rattus norvegicus (Norway rat)] -> EXPERIMENT ID: 39 -> PUBMED ID: 15210972 -> https://www.ncbi.nlm.nih.gov/pubmed/15210972 -> Cells release prions in association with exosomes -> mouse -> Mus musculus
ExoCarta_Protein_mRNA[which(trimws(ExoCarta_Protein_mRNA[,'SPECIES']) == 'Mus musculus' & trimws(ExoCarta_Protein_mRNA[,'ENTREZ.GENE.ID']) == '301252'),'GENE.SYMBOL'] = 'Hsp90ab1'
ExoCarta_Protein_mRNA[which(trimws(ExoCarta_Protein_mRNA[,'SPECIES']) == 'Mus musculus' & trimws(ExoCarta_Protein_mRNA[,'ENTREZ.GENE.ID']) == '301252'),'ENTREZ.GENE.ID'] = '15516'

# https://www.ncbi.nlm.nih.gov/gene/51552 [Homo sapiens (human)] -> EXPERIMENT ID: 73 -> PUBMED ID: 18494037 -> https://www.ncbi.nlm.nih.gov/pubmed/18494037 -> Difference gel electrophoresis analysis of Ras-transformed fibroblast cell-derived exosomes -> murine fibroblast NIH3T3 cells -> Mus musculus, mouse
ExoCarta_Protein_mRNA[which(trimws(ExoCarta_Protein_mRNA[,'SPECIES']) == 'Mus musculus' & trimws(ExoCarta_Protein_mRNA[,'ENTREZ.GENE.ID']) == '51552'),'GENE.SYMBOL'] = 'Rab14'
ExoCarta_Protein_mRNA[which(trimws(ExoCarta_Protein_mRNA[,'SPECIES']) == 'Mus musculus' & trimws(ExoCarta_Protein_mRNA[,'ENTREZ.GENE.ID']) == '51552'),'ENTREZ.GENE.ID'] = '68365'

ExoCarta_Protein_mRNA.hsa <- ExoCarta_Protein_mRNA[which(trimws(ExoCarta_Protein_mRNA[,'SPECIES']) == 'Homo sapiens'),]
ExoCarta_Protein_mRNA.mmu <- ExoCarta_Protein_mRNA[which(trimws(ExoCarta_Protein_mRNA[,'SPECIES']) == 'Mus musculus'),]
ExoCarta_Protein_mRNA.rno <- ExoCarta_Protein_mRNA[which(trimws(ExoCarta_Protein_mRNA[,'SPECIES']) == 'Rattus norvegicus'),]
#dim(ExoCarta_Protein_mRNA.hsa)
#dim(ExoCarta_Protein_mRNA.mmu)
#dim(ExoCarta_Protein_mRNA.rno)


# /// human format ///
# --- protein ---
ExoCarta_Protein.hsa <- ExoCarta_Protein_mRNA.hsa[which(trimws(ExoCarta_Protein_mRNA.hsa[,'CONTENT.TYPE']) == 'protein'),]
#dim(ExoCarta_Protein.hsa)
ExoCarta_Protein.hsa.CONTENT_ID <- trimws(ExoCarta_Protein.hsa[,'CONTENT.ID'])
ExoCarta_Protein.hsa.CONTENT_TYPE <- trimws(ExoCarta_Protein.hsa[,'CONTENT.TYPE'])
ExoCarta_Protein.hsa.ENTREZ_GENE_ID <- trimws(ExoCarta_Protein.hsa[,'ENTREZ.GENE.ID'])
ExoCarta_Protein.hsa.GENE_SYMBOL <- trimws(ExoCarta_Protein.hsa[,'GENE.SYMBOL'])
ExoCarta_Protein.hsa.SPECIES <- trimws(ExoCarta_Protein.hsa[,'SPECIES'])
ExoCarta_Protein.hsa.EXPERIMENT_ID <- trimws(ExoCarta_Protein.hsa[,'EXPERIMENT.ID'])
ExoCarta_Protein.hsa.METHODS <- trimws(ExoCarta_Protein.hsa[,'METHODS'])
#length(ExoCarta_Protein.hsa.ENTREZ_GENE_ID)
ExoCarta_Protein.hsa.ENTREZ_GENE_ID.unique = unique(ExoCarta_Protein.hsa.ENTREZ_GENE_ID)
#length(ExoCarta_Protein.hsa.ENTREZ_GENE_ID.unique)

# remove non-protein id (Gene type: ncRNA)
notprotein.idx <- c()
notprotein.idx = c(notprotein.idx, which(ExoCarta_Protein.hsa.ENTREZ_GENE_ID.unique == '57653'))	# https://www.ncbi.nlm.nih.gov/gene/57653	http://www.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000255036;r=9:97238497-97377287
notprotein.idx = c(notprotein.idx, which(ExoCarta_Protein.hsa.ENTREZ_GENE_ID.unique == '442497'))	# https://www.ncbi.nlm.nih.gov/gene/442497	http://www.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000249574;r=7:379359-382712;t=ENST00000515213
ExoCarta_Protein.hsa.ENTREZ_GENE_ID.unique = ExoCarta_Protein.hsa.ENTREZ_GENE_ID.unique[-notprotein.idx]
ExoCarta_Protein.hsa.ENTREZ_GENE_ID.unique = unique(ExoCarta_Protein.hsa.ENTREZ_GENE_ID.unique)
#length(ExoCarta_Protein.hsa.ENTREZ_GENE_ID.unique)

# --- mrna ---
ExoCarta_mRNA.hsa <- ExoCarta_Protein_mRNA.hsa[which(trimws(ExoCarta_Protein_mRNA.hsa[,'CONTENT.TYPE']) == 'mRNA' | trimws(ExoCarta_Protein_mRNA.hsa[,'CONTENT.TYPE']) == 'mrna'),]
#dim(ExoCarta_mRNA.hsa)
ExoCarta_mRNA.hsa.CONTENT_ID <- trimws(ExoCarta_mRNA.hsa[,'CONTENT.ID'])
ExoCarta_mRNA.hsa.CONTENT_TYPE <- trimws(ExoCarta_mRNA.hsa[,'CONTENT.TYPE'])
ExoCarta_mRNA.hsa.ENTREZ_GENE_ID <- trimws(ExoCarta_mRNA.hsa[,'ENTREZ.GENE.ID'])
ExoCarta_mRNA.hsa.GENE_SYMBOL <- trimws(ExoCarta_mRNA.hsa[,'GENE.SYMBOL'])
ExoCarta_mRNA.hsa.SPECIES <- trimws(ExoCarta_mRNA.hsa[,'SPECIES'])
ExoCarta_mRNA.hsa.EXPERIMENT_ID <- trimws(ExoCarta_mRNA.hsa[,'EXPERIMENT.ID'])
ExoCarta_mRNA.hsa.METHODS <- trimws(ExoCarta_mRNA.hsa[,'METHODS'])
#length(ExoCarta_mRNA.hsa.ENTREZ_GENE_ID)
ExoCarta_mRNA.hsa.ENTREZ_GENE_ID.unique = unique(ExoCarta_mRNA.hsa.ENTREZ_GENE_ID)
#length(ExoCarta_mRNA.hsa.ENTREZ_GENE_ID.unique)

# /// mouse format ///
# --- protein ---
ExoCarta_Protein.mmu <- ExoCarta_Protein_mRNA.mmu[which(trimws(ExoCarta_Protein_mRNA.mmu[,'CONTENT.TYPE']) == 'protein'),]
#dim(ExoCarta_Protein.mmu)
ExoCarta_Protein.mmu.CONTENT_ID <- trimws(ExoCarta_Protein.mmu[,'CONTENT.ID'])
ExoCarta_Protein.mmu.CONTENT_TYPE <- trimws(ExoCarta_Protein.mmu[,'CONTENT.TYPE'])
ExoCarta_Protein.mmu.ENTREZ_GENE_ID <- trimws(ExoCarta_Protein.mmu[,'ENTREZ.GENE.ID'])
ExoCarta_Protein.mmu.GENE_SYMBOL <- trimws(ExoCarta_Protein.mmu[,'GENE.SYMBOL'])
ExoCarta_Protein.mmu.SPECIES <- trimws(ExoCarta_Protein.mmu[,'SPECIES'])
ExoCarta_Protein.mmu.EXPERIMENT_ID <- trimws(ExoCarta_Protein.mmu[,'EXPERIMENT.ID'])
ExoCarta_Protein.mmu.METHODS <- trimws(ExoCarta_Protein.mmu[,'METHODS'])
#length(ExoCarta_Protein.mmu.ENTREZ_GENE_ID)
ExoCarta_Protein.mmu.ENTREZ_GENE_ID.unique = unique(ExoCarta_Protein.mmu.ENTREZ_GENE_ID)
#length(ExoCarta_Protein.mmu.ENTREZ_GENE_ID.unique)

# --- mrna ---
ExoCarta_mRNA.mmu <- ExoCarta_Protein_mRNA.mmu[which(trimws(ExoCarta_Protein_mRNA.mmu[,'CONTENT.TYPE']) == 'mRNA' | trimws(ExoCarta_Protein_mRNA.mmu[,'CONTENT.TYPE']) == 'mrna'),]
#dim(ExoCarta_mRNA.mmu)
ExoCarta_mRNA.mmu.CONTENT_ID <- trimws(ExoCarta_mRNA.mmu[,'CONTENT.ID'])
ExoCarta_mRNA.mmu.CONTENT_TYPE <- trimws(ExoCarta_mRNA.mmu[,'CONTENT.TYPE'])
ExoCarta_mRNA.mmu.ENTREZ_GENE_ID <- trimws(ExoCarta_mRNA.mmu[,'ENTREZ.GENE.ID'])
ExoCarta_mRNA.mmu.GENE_SYMBOL <- trimws(ExoCarta_mRNA.mmu[,'GENE.SYMBOL'])
ExoCarta_mRNA.mmu.SPECIES <- trimws(ExoCarta_mRNA.mmu[,'SPECIES'])
ExoCarta_mRNA.mmu.EXPERIMENT_ID <- trimws(ExoCarta_mRNA.mmu[,'EXPERIMENT.ID'])
ExoCarta_mRNA.mmu.METHODS <- trimws(ExoCarta_mRNA.mmu[,'METHODS'])
#length(ExoCarta_mRNA.mmu.ENTREZ_GENE_ID)
ExoCarta_mRNA.mmu.ENTREZ_GENE_ID.unique = unique(ExoCarta_mRNA.mmu.ENTREZ_GENE_ID)
#length(ExoCarta_mRNA.mmu.ENTREZ_GENE_ID.unique)


# === id mapping ================================================================================================================
# /// human format ///
mainDir <- paste('./Data/',Sample,'/10_datasetOverlap/Database',sep="")
subDir <- "hsa"
if (file.exists(file.path(mainDir, subDir))){
} else {
	dir.create(file.path(mainDir, subDir))
}

ensembl = ensembl.hsa
#listFilters(ensembl)[,c('name','description')]
#listAttributes(ensembl)[,c('name','description')]

# --- search for IDs by entrezgene ---
a = ExoCarta_Protein.hsa.ENTREZ_GENE_ID.unique
b = ensembl.hsa
ExoCarta_Protein.entrezgene.attributes <- getBM(attributes=c('hgnc_symbol','hgnc_id','uniprotswissprot'), filters='entrezgene', values=a, mart=b)
ExoCarta_Protein.entrezgene.wikigene_description <- getBM(attributes=c('hgnc_symbol','wikigene_description'), filters='entrezgene', values=a, mart=b)
ExoCarta_Protein.entrezgene.wikigene_description.table = merge(ExoCarta_Protein.entrezgene.wikigene_description, ExoCarta_Protein.entrezgene.attributes, by.x="hgnc_symbol", by.y="hgnc_symbol")
ExoCarta_Protein.entrezgene.attributes = ExoCarta_Protein.entrezgene.wikigene_description.table[,c('hgnc_symbol','hgnc_id','uniprotswissprot','wikigene_description'),drop=FALSE]
colnames(ExoCarta_Protein.entrezgene.attributes)[which(colnames(ExoCarta_Protein.entrezgene.attributes)=='wikigene_description')] = 'description'
ExoCarta_Protein.entrezgene.attributes = unique(ExoCarta_Protein.entrezgene.attributes)
ExoCarta_Protein.hgnc_id <- getBM(attributes=c('entrezgene','hgnc_id'), filters='entrezgene', values=a, mart=b)
ExoCarta_Protein.hgnc_id.table = merge(ExoCarta_Protein.hgnc_id, ExoCarta_Protein.entrezgene.attributes, by.x="hgnc_id", by.y="hgnc_id")

# https://www.ncbi.nlm.nih.gov/gene/226; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:414; http://www.uniprot.org/uniprot/P04075
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='226'),'hgnc_symbol'] = 'ALDOA'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='226'),'hgnc_id'] = 'HGNC:414'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='226'),'uniprotswissprot'] = 'P04075'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='226'),'description'] = 'aldolase, fructose-bisphosphate A [Source:HGNC Symbol;Acc:HGNC:414]'
# https://www.ncbi.nlm.nih.gov/gene/640; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:1057; http://www.uniprot.org/uniprot/P51451
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='640'),'hgnc_symbol'] = 'BLK'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='640'),'hgnc_id'] = 'HGNC:1057'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='640'),'uniprotswissprot'] = 'P51451'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='640'),'description'] = 'BLK proto-oncogene, Src family tyrosine kinase [Source:HGNC Symbol;Acc:HGNC:1057]'
# https://www.ncbi.nlm.nih.gov/gene/720; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:1323; http://www.uniprot.org/uniprot/P0C0L4
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='720'),'hgnc_symbol'] = 'C4A'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='720'),'hgnc_id'] = 'HGNC:1323'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='720'),'uniprotswissprot'] = 'P0C0L4'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='720'),'description'] = 'complement C4A (Rodgers blood group) [Source:HGNC Symbol;Acc:HGNC:1323]'
# https://www.ncbi.nlm.nih.gov/gene/1409; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:2388; http://www.uniprot.org/uniprot/P02489
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='1409'),'hgnc_symbol'] = 'CRYAA'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='1409'),'hgnc_id'] = 'HGNC:2388'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='1409'),'uniprotswissprot'] = 'P02489'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='1409'),'description'] = 'crystallin alpha A [Source:HGNC Symbol;Acc:HGNC:2388]'
# https://www.ncbi.nlm.nih.gov/gene/1565; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:2625; http://www.uniprot.org/uniprot/P10635
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='1565'),'hgnc_symbol'] = 'CYP2D6'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='1565'),'hgnc_id'] = 'HGNC:2625'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='1565'),'uniprotswissprot'] = 'P10635'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='1565'),'description'] = 'cytochrome P450 family 2 subfamily D member 6 [Source:HGNC Symbol;Acc:HGNC:2625]'
# https://www.ncbi.nlm.nih.gov/gene/1652; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:2732; http://www.uniprot.org/uniprot/P30046
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='1652'),'hgnc_symbol'] = 'DDT'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='1652'),'hgnc_id'] = 'HGNC:2732'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='1652'),'uniprotswissprot'] = 'P30046'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='1652'),'description'] = 'D-dopachrome tautomerase [Source:HGNC Symbol;Acc:HGNC:2732]'
# https://www.ncbi.nlm.nih.gov/gene/2222; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:3629; http://www.uniprot.org/uniprot/P37268
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='2222'),'hgnc_symbol'] = 'FDFT1'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='2222'),'hgnc_id'] = 'HGNC:3629'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='2222'),'uniprotswissprot'] = 'P37268'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='2222'),'description'] = 'farnesyl-diphosphate farnesyltransferase 1 [Source:HGNC Symbol;Acc:HGNC:3629]'
# https://www.ncbi.nlm.nih.gov/gene/2632; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:4180; http://www.uniprot.org/uniprot/Q04446
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='2632'),'hgnc_symbol'] = 'GBE1'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='2632'),'hgnc_id'] = 'HGNC:4180'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='2632'),'uniprotswissprot'] = 'Q04446'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='2632'),'description'] = '1,4-alpha-glucan branching enzyme 1 [Source:HGNC Symbol;Acc:HGNC:4180]'
# https://www.ncbi.nlm.nih.gov/gene/2657; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:4214; http://www.uniprot.org/uniprot/P27539
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='2657'),'hgnc_symbol'] = 'GDF1'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='2657'),'hgnc_id'] = 'HGNC:4214'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='2657'),'uniprotswissprot'] = 'P27539'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='2657'),'description'] = 'growth differentiation factor 1 [Source:HGNC Symbol;Acc:HGNC:4214]'
# https://www.ncbi.nlm.nih.gov/gene/3123; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:4948
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='3123'),'hgnc_symbol'] = 'HLA-DRB1'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='3123'),'hgnc_id'] = 'HGNC:4948'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='3123'),'uniprotswissprot'] = 'P01911;P01912;P04229;P20039;Q29974;Q30134;Q30167;Q5Y7A7;Q95IE3;Q9GIY3'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='3123'),'description'] = 'major histocompatibility complex, class II, DR beta 1 [Source:HGNC Symbol;Acc:HGNC:4948]'
# https://www.ncbi.nlm.nih.gov/gene/3125; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:4951; http://www.uniprot.org/uniprot/P79483
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='3125'),'hgnc_symbol'] = 'HLA-DRB3'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='3125'),'hgnc_id'] = 'HGNC:4951'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='3125'),'uniprotswissprot'] = 'P79483'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='3125'),'description'] = 'major histocompatibility complex, class II, DR beta 3 [Source:HGNC Symbol;Acc:HGNC:4951]'
# https://www.ncbi.nlm.nih.gov/gene/3481; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:5466; http://www.uniprot.org/uniprot/P01344
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='3481'),'hgnc_symbol'] = 'IGF2'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='3481'),'hgnc_id'] = 'HGNC:5466'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='3481'),'uniprotswissprot'] = 'P01344'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='3481'),'description'] = 'insulin like growth factor 2 [Source:HGNC Symbol;Acc:HGNC:5466]'
# https://www.ncbi.nlm.nih.gov/gene/4681; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:7650; http://www.uniprot.org/uniprot/P41271
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='4681'),'hgnc_symbol'] = 'NBL1'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='4681'),'hgnc_id'] = 'HGNC:7650'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='4681'),'uniprotswissprot'] = 'P41271'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='4681'),'description'] = 'NBL1, DAN family BMP antagonist [Source:HGNC Symbol;Acc:HGNC:7650]'
# https://www.ncbi.nlm.nih.gov/gene/5430; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:9187; http://www.uniprot.org/uniprot/P24928
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='5430'),'hgnc_symbol'] = 'POLR2A'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='5430'),'hgnc_id'] = 'HGNC:9187'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='5430'),'uniprotswissprot'] = 'P24928'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='5430'),'description'] = 'RNA polymerase II subunit A [Source:HGNC Symbol;Acc:HGNC:9187]'
# https://www.ncbi.nlm.nih.gov/gene/6203; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:10442; http://www.uniprot.org/uniprot/P46781
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='6203'),'hgnc_symbol'] = 'RPS9'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='6203'),'hgnc_id'] = 'HGNC:10442'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='6203'),'uniprotswissprot'] = 'P46781'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='6203'),'description'] = 'ribosomal protein S9 [Source:HGNC Symbol;Acc:HGNC:10442]'
# https://www.ncbi.nlm.nih.gov/gene/6693; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:11249; http://www.uniprot.org/uniprot/P16150
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='6693'),'hgnc_symbol'] = 'SPN'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='6693'),'hgnc_id'] = 'HGNC:11249'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='6693'),'uniprotswissprot'] = 'P16150'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='6693'),'description'] = 'sialophorin [Source:HGNC Symbol;Acc:HGNC:11249]'
# https://www.ncbi.nlm.nih.gov/gene/6905; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:11582; http://www.uniprot.org/uniprot/Q15813
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='6905'),'hgnc_symbol'] = 'TBCE'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='6905'),'hgnc_id'] = 'HGNC:11582'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='6905'),'uniprotswissprot'] = 'Q15813'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='6905'),'description'] = 'tubulin folding cofactor E [Source:HGNC Symbol;Acc:HGNC:11582]'
# https://www.ncbi.nlm.nih.gov/gene/7375; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:12627; http://www.uniprot.org/uniprot/Q13107
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='7375'),'hgnc_symbol'] = 'USP4'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='7375'),'hgnc_id'] = 'HGNC:12627'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='7375'),'uniprotswissprot'] = 'Q13107'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='7375'),'description'] = 'ubiquitin specific peptidase 4 [Source:HGNC Symbol;Acc:HGNC:12627]'
# https://www.ncbi.nlm.nih.gov/gene/8411; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:3185; http://www.uniprot.org/uniprot/Q15075
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='8411'),'hgnc_symbol'] = 'EEA1'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='8411'),'hgnc_id'] = 'HGNC:3185'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='8411'),'uniprotswissprot'] = 'Q15075'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='8411'),'description'] = 'early endosome antigen 1 [Source:HGNC Symbol;Acc:HGNC:3185]'
# https://www.ncbi.nlm.nih.gov/gene/8471; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:6128; http://www.uniprot.org/uniprot/O14654
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='8471'),'hgnc_symbol'] = 'IRS4'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='8471'),'hgnc_id'] = 'HGNC:6128'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='8471'),'uniprotswissprot'] = 'O14654'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='8471'),'description'] = 'insulin receptor substrate 4 [Source:HGNC Symbol;Acc:HGNC:6128]'
# https://www.ncbi.nlm.nih.gov/gene/8519; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:5412; http://www.uniprot.org/uniprot/P13164
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='8519'),'hgnc_symbol'] = 'IFITM1'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='8519'),'hgnc_id'] = 'HGNC:5412'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='8519'),'uniprotswissprot'] = 'P13164'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='8519'),'description'] = 'interferon induced transmembrane protein 1 [Source:HGNC Symbol;Acc:HGNC:5412]'
# https://www.ncbi.nlm.nih.gov/gene/8566; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:8819; http://www.uniprot.org/uniprot/O00764
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='8566'),'hgnc_symbol'] = 'PDXK'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='8566'),'hgnc_id'] = 'HGNC:8819'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='8566'),'uniprotswissprot'] = 'O00764'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='8566'),'description'] = 'pyridoxal kinase [Source:HGNC Symbol;Acc:HGNC:8819]'
# https://www.ncbi.nlm.nih.gov/gene/8857; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:13572; http://www.uniprot.org/uniprot/Q9Y6R7
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='8857'),'hgnc_symbol'] = 'FCGBP'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='8857'),'hgnc_id'] = 'HGNC:13572'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='8857'),'uniprotswissprot'] = 'Q9Y6R7'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='8857'),'description'] = 'Fc fragment of IgG binding protein [Source:HGNC Symbol;Acc:HGNC:13572]'
# https://www.ncbi.nlm.nih.gov/gene/8878; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:11280; http://www.uniprot.org/uniprot/Q13501
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='8878'),'hgnc_symbol'] = 'SQSTM1'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='8878'),'hgnc_id'] = 'HGNC:11280'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='8878'),'uniprotswissprot'] = 'Q13501'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='8878'),'description'] = 'sequestosome 1 [Source:HGNC Symbol;Acc:HGNC:11280]'
# https://www.ncbi.nlm.nih.gov/gene/9138; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:681; http://www.uniprot.org/uniprot/Q92888
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='9138'),'hgnc_symbol'] = 'ARHGEF1'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='9138'),'hgnc_id'] = 'HGNC:681'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='9138'),'uniprotswissprot'] = 'Q92888'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='9138'),'description'] = 'Rho guanine nucleotide exchange factor 1 [Source:HGNC Symbol;Acc:HGNC:681]'
# https://www.ncbi.nlm.nih.gov/gene/9414; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:11828; http://www.uniprot.org/uniprot/Q9UDY2
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='9414'),'hgnc_symbol'] = 'TJP2'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='9414'),'hgnc_id'] = 'HGNC:11828'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='9414'),'uniprotswissprot'] = 'Q9UDY2'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='9414'),'description'] = 'tight junction protein 2 [Source:HGNC Symbol;Acc:HGNC:11828]'
# https://www.ncbi.nlm.nih.gov/gene/9554; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:10700; http://www.uniprot.org/uniprot/O75396
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='9554'),'hgnc_symbol'] = 'SEC22B'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='9554'),'hgnc_id'] = 'HGNC:10700'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='9554'),'uniprotswissprot'] = 'O75396'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='9554'),'description'] = 'SEC22 homolog B, vesicle trafficking protein (gene/pseudogene) [Source:HGNC Symbol;Acc:HGNC:10700]'
# https://www.ncbi.nlm.nih.gov/gene/9659; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:15580; http://www.uniprot.org/uniprot/Q5VU43
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='9659'),'hgnc_symbol'] = 'PDE4DIP'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='9659'),'hgnc_id'] = 'HGNC:15580'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='9659'),'uniprotswissprot'] = 'Q5VU43'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='9659'),'description'] = 'phosphodiesterase 4D interacting protein [Source:HGNC Symbol;Acc:HGNC:15580]'
# https://www.ncbi.nlm.nih.gov/gene/9782; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:6912; http://www.uniprot.org/uniprot/P43243
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='9782'),'hgnc_symbol'] = 'MATR3'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='9782'),'hgnc_id'] = 'HGNC:6912'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='9782'),'uniprotswissprot'] = 'P43243'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='9782'),'description'] = 'matrin 3 [Source:HGNC Symbol;Acc:HGNC:6912]'
# https://www.ncbi.nlm.nih.gov/gene/10061; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:71; http://www.uniprot.org/uniprot/Q9UG63
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='10061'),'hgnc_symbol'] = 'ABCF2'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='10061'),'hgnc_id'] = 'HGNC:71'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='10061'),'uniprotswissprot'] = 'Q9UG63'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='10061'),'description'] = 'ATP binding cassette subfamily F member 2 [Source:HGNC Symbol;Acc:HGNC:71]'
# https://www.ncbi.nlm.nih.gov/gene/10436; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:16912; http://www.uniprot.org/uniprot/Q92979
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='10436'),'hgnc_symbol'] = 'EMG1'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='10436'),'hgnc_id'] = 'HGNC:16912'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='10436'),'uniprotswissprot'] = 'Q92979'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='10436'),'description'] = 'EMG1, N1-specific pseudouridine methyltransferase [Source:HGNC Symbol;Acc:HGNC:16912]'
# https://www.ncbi.nlm.nih.gov/gene/10715; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:14253; http://www.uniprot.org/uniprot/P27544
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='10715'),'hgnc_symbol'] = 'CERS1'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='10715'),'hgnc_id'] = 'HGNC:14253'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='10715'),'uniprotswissprot'] = 'P27544'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='10715'),'description'] = 'ceramide synthase 1 [Source:HGNC Symbol;Acc:HGNC:14253]'
# https://www.ncbi.nlm.nih.gov/gene/11046; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:20799; http://www.uniprot.org/uniprot/Q76EJ3
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='11046'),'hgnc_symbol'] = 'SLC35D2'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='11046'),'hgnc_id'] = 'HGNC:20799'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='11046'),'uniprotswissprot'] = 'Q76EJ3'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='11046'),'description'] = 'solute carrier family 35 member D2 [Source:HGNC Symbol;Acc:HGNC:20799]'
# https://www.ncbi.nlm.nih.gov/gene/11163; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:8051; http://www.uniprot.org/uniprot/Q9NZJ9
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='11163'),'hgnc_symbol'] = 'NUDT4'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='11163'),'hgnc_id'] = 'HGNC:8051'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='11163'),'uniprotswissprot'] = 'Q9NZJ9'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='11163'),'description'] = 'nudix hydrolase 4 [Source:HGNC Symbol;Acc:HGNC:8051]'
# https://www.ncbi.nlm.nih.gov/gene/11272; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:18020; http://www.uniprot.org/uniprot/Q16378
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='11272'),'hgnc_symbol'] = 'PRR4'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='11272'),'hgnc_id'] = 'HGNC:18020'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='11272'),'uniprotswissprot'] = 'Q16378'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='11272'),'description'] = 'proline rich 4 [Source:HGNC Symbol;Acc:HGNC:18020]'
# https://www.ncbi.nlm.nih.gov/gene/23370; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:17090; http://www.uniprot.org/uniprot/Q6ZSZ5
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='23370'),'hgnc_symbol'] = 'ARHGEF18'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='23370'),'hgnc_id'] = 'HGNC:17090'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='23370'),'uniprotswissprot'] = 'Q6ZSZ5'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='23370'),'description'] = 'Rho/Rac guanine nucleotide exchange factor 18 [Source:HGNC Symbol;Acc:HGNC:17090]'
# https://www.ncbi.nlm.nih.gov/gene/23475; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:9755; http://www.uniprot.org/uniprot/Q15274
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='23475'),'hgnc_symbol'] = 'QPRT'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='23475'),'hgnc_id'] = 'HGNC:9755'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='23475'),'uniprotswissprot'] = 'Q15274'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='23475'),'description'] = 'quinolinate phosphoribosyltransferase [Source:HGNC Symbol;Acc:HGNC:9755]'
# https://www.ncbi.nlm.nih.gov/gene/23499; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:13664; http://www.uniprot.org/uniprot/Q9UPN3
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='23499'),'hgnc_symbol'] = 'MACF1'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='23499'),'hgnc_id'] = 'HGNC:13664'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='23499'),'uniprotswissprot'] = 'Q9UPN3'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='23499'),'description'] = 'microtubule-actin crosslinking factor 1 [Source:HGNC Symbol;Acc:HGNC:13664]'
# https://www.ncbi.nlm.nih.gov/gene/23637; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:17155; http://www.uniprot.org/uniprot/Q9Y3P9
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='23637'),'hgnc_symbol'] = 'RABGAP1'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='23637'),'hgnc_id'] = 'HGNC:17155'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='23637'),'uniprotswissprot'] = 'Q9Y3P9'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='23637'),'description'] = 'RAB GTPase activating protein 1 [Source:HGNC Symbol;Acc:HGNC:17155]'
# https://www.ncbi.nlm.nih.gov/gene/26121; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:15446; http://www.uniprot.org/uniprot/Q8WWY3
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='26121'),'hgnc_symbol'] = 'PRPF31'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='26121'),'hgnc_id'] = 'HGNC:15446'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='26121'),'uniprotswissprot'] = 'Q8WWY3'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='26121'),'description'] = 'pre-mRNA processing factor 31 [Source:HGNC Symbol;Acc:HGNC:15446]'
# https://www.ncbi.nlm.nih.gov/gene/26585; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:2001; http://www.uniprot.org/uniprot/O60565
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='26585'),'hgnc_symbol'] = 'GREM1'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='26585'),'hgnc_id'] = 'HGNC:2001'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='26585'),'uniprotswissprot'] = 'O60565'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='26585'),'description'] = 'gremlin 1, DAN family BMP antagonist [Source:HGNC Symbol;Acc:HGNC:2001]'
# https://www.ncbi.nlm.nih.gov/gene/51206; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:14388; http://www.uniprot.org/uniprot/Q9HCN6
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='51206'),'hgnc_symbol'] = 'GP6'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='51206'),'hgnc_id'] = 'HGNC:14388'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='51206'),'uniprotswissprot'] = 'Q9HCN6'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='51206'),'description'] = 'glycoprotein VI platelet [Source:HGNC Symbol;Acc:HGNC:14388]'
# https://www.ncbi.nlm.nih.gov/gene/56924; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:16061; http://www.uniprot.org/uniprot/Q9NQU5
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='56924'),'hgnc_symbol'] = 'PAK6'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='56924'),'hgnc_id'] = 'HGNC:16061'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='56924'),'uniprotswissprot'] = 'Q9NQU5'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='56924'),'description'] = 'p21 (RAC1) activated kinase 6 [Source:HGNC Symbol;Acc:HGNC:16061]'
# https://www.ncbi.nlm.nih.gov/gene/57535; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:29618; http://www.uniprot.org/uniprot/Q6UXG2
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='57535'),'hgnc_symbol'] = 'KIAA1324'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='57535'),'hgnc_id'] = 'HGNC:29618'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='57535'),'uniprotswissprot'] = 'Q6UXG2'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='57535'),'description'] = 'KIAA1324 [Source:HGNC Symbol;Acc:HGNC:29618]'
# https://www.ncbi.nlm.nih.gov/gene/80006; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:25828; http://www.uniprot.org/uniprot/A5PLN9
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='80006'),'hgnc_symbol'] = 'TRAPPC13'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='80006'),'hgnc_id'] = 'HGNC:25828'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='80006'),'uniprotswissprot'] = 'A5PLN9'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='80006'),'description'] = 'trafficking protein particle complex 13 [Source:HGNC Symbol;Acc:HGNC:25828]'
# https://www.ncbi.nlm.nih.gov/gene/83986; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:14163; http://www.uniprot.org/uniprot/Q9H0X4
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='83986'),'hgnc_symbol'] = 'FAM234A'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='83986'),'hgnc_id'] = 'HGNC:14163'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='83986'),'uniprotswissprot'] = 'Q9H0X4'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='83986'),'description'] = 'family with sequence similarity 234 member A [Source:HGNC Symbol;Acc:HGNC:14163]'
# https://www.ncbi.nlm.nih.gov/gene/84858; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:23589; http://www.uniprot.org/uniprot/Q96F45
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='84858'),'hgnc_symbol'] = 'ZNF503'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='84858'),'hgnc_id'] = 'HGNC:23589'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='84858'),'uniprotswissprot'] = 'Q96F45'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='84858'),'description'] = 'zinc finger protein 503 [Source:HGNC Symbol;Acc:HGNC:23589]'
# https://www.ncbi.nlm.nih.gov/gene/84876; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:25896; http://www.uniprot.org/uniprot/Q96D31
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='84876'),'hgnc_symbol'] = 'ORAI1'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='84876'),'hgnc_id'] = 'HGNC:25896'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='84876'),'uniprotswissprot'] = 'Q96D31'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='84876'),'description'] = 'ORAI calcium release-activated calcium modulator 1 [Source:HGNC Symbol;Acc:HGNC:25896]'
# https://www.ncbi.nlm.nih.gov/gene/150786; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:30272; http://www.uniprot.org/uniprot/Q53S08
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='150786'),'hgnc_symbol'] = 'RAB6D'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='150786'),'hgnc_id'] = 'HGNC:30272'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='150786'),'uniprotswissprot'] = 'Q53S08'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='150786'),'description'] = 'RAB6D, member RAS oncogene family [Source:HGNC Symbol;Acc:HGNC:30272]'
# https://www.ncbi.nlm.nih.gov/gene/203569; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:31804; http://www.uniprot.org/uniprot/Q7Z2X7
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='203569'),'hgnc_symbol'] = 'PAGE2'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='203569'),'hgnc_id'] = 'HGNC:31804'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='203569'),'uniprotswissprot'] = 'Q7Z2X7'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='203569'),'description'] = 'PAGE family member 2 [Source:HGNC Symbol;Acc:HGNC:31804]'
# https://www.ncbi.nlm.nih.gov/gene/246181; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:24056; http://www.uniprot.org/uniprot/Q8NHP1
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='246181'),'hgnc_symbol'] = 'AKR7L'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='246181'),'hgnc_id'] = 'HGNC:24056'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='246181'),'uniprotswissprot'] = 'Q8NHP1'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='246181'),'description'] = 'aldo-keto reductase family 7 like (gene/pseudogene) [Source:HGNC Symbol;Acc:HGNC:24056]'
# https://www.ncbi.nlm.nih.gov/gene/283820; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:22652; http://www.uniprot.org/uniprot/Q5JPE7
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='283820'),'hgnc_symbol'] = 'NOMO2'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='283820'),'hgnc_id'] = 'HGNC:22652'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='283820'),'uniprotswissprot'] = 'Q5JPE7'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='283820'),'description'] = 'NODAL modulator 2 [Source:HGNC Symbol;Acc:HGNC:22652]'
# https://www.ncbi.nlm.nih.gov/gene/337873; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:20516; http://www.uniprot.org/uniprot/Q6DN03
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='337873'),'hgnc_symbol'] = 'HIST2H2BC'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='337873'),'hgnc_id'] = 'HGNC:20516'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='337873'),'uniprotswissprot'] = 'Q6DN03'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='337873'),'description'] = 'histone cluster 2 H2B family member c (pseudogene) [Source:HGNC Symbol;Acc:HGNC:20516]'
# https://www.ncbi.nlm.nih.gov/gene/440563; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:48813; http://www.uniprot.org/uniprot/B2RXH8
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='440563'),'hgnc_symbol'] = 'HNRNPCL2'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='440563'),'hgnc_id'] = 'HGNC:48813'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='440563'),'uniprotswissprot'] = 'B2RXH8'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='440563'),'description'] = 'heterogeneous nuclear ribonucleoprotein C-like 2 [Source:HGNC Symbol;Acc:HGNC:48813]'
# https://www.ncbi.nlm.nih.gov/gene/441459; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:23644; http://www.uniprot.org/uniprot/A2A2Z9
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='441459'),'hgnc_symbol'] = 'ANKRD18B'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='441459'),'hgnc_id'] = 'HGNC:23644'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='441459'),'uniprotswissprot'] = 'A2A2Z9'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='441459'),'description'] = 'ankyrin repeat domain 18B [Source:HGNC Symbol;Acc:HGNC:23644]'
# https://www.ncbi.nlm.nih.gov/gene/653145; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:546; http://www.uniprot.org/uniprot/P13928
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='653145'),'hgnc_symbol'] = 'ANXA8'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='653145'),'hgnc_id'] = 'HGNC:546'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='653145'),'uniprotswissprot'] = 'P13928'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='653145'),'description'] = 'annexin A8 [Source:HGNC Symbol;Acc:HGNC:546]'
# https://www.ncbi.nlm.nih.gov/gene/100996720
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='100996720'),'hgnc_symbol'] = 'LOC100996720'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='100996720'),'hgnc_id'] = ''
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='100996720'),'uniprotswissprot'] = ''
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='100996720'),'description'] = 'uncharacterized LOC100996720'
# https://www.ncbi.nlm.nih.gov/gene/101060301; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:51333; http://www.uniprot.org/uniprot/P0DMR1
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='101060301'),'hgnc_symbol'] = 'HNRNPCL4'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='101060301'),'hgnc_id'] = 'HGNC:51333'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='101060301'),'uniprotswissprot'] = 'P0DMR1'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='101060301'),'description'] = 'heterogeneous nuclear ribonucleoprotein C-like 4 [Source:HGNC Symbol;Acc:HGNC:51333]'
# https://www.ncbi.nlm.nih.gov/gene/102724334; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:4762; http://www.uniprot.org/uniprot/P57053
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='102724334'),'hgnc_symbol'] = 'H2BFS'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='102724334'),'hgnc_id'] = 'HGNC:4762'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='102724334'),'uniprotswissprot'] = 'P57053'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='102724334'),'description'] = 'H2B histone family member S [Source:HGNC Symbol;Acc:HGNC:4762]'
# https://www.ncbi.nlm.nih.gov/gene/102724631; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:51240; http://www.uniprot.org/uniprot/A0JP26
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='102724631'),'hgnc_symbol'] = 'POTEB3'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='102724631'),'hgnc_id'] = 'HGNC:51240'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='102724631'),'uniprotswissprot'] = 'A0JP26'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='102724631'),'description'] = 'POTE ankyrin domain family member B3 [Source:HGNC Symbol;Acc:HGNC:51240]'
# https://www.ncbi.nlm.nih.gov/gene/102724652; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:2388; http://www.uniprot.org/uniprot/P02489
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='102724652'),'hgnc_symbol'] = 'CRYAA'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='102724652'),'hgnc_id'] = 'HGNC:2388'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='102724652'),'uniprotswissprot'] = 'P02489'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='102724652'),'description'] = 'crystallin alpha A [Source:HGNC Symbol;Acc:HGNC:2388]'

# https://www.ncbi.nlm.nih.gov/gene/276; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:474; http://www.uniprot.org/uniprot/P04745
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='276'),'hgnc_symbol'] = 'AMY1A'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='276'),'hgnc_id'] = 'HGNC:474'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='276'),'uniprotswissprot'] = 'P04745'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='276'),'description'] = 'amylase, alpha 1A (salivary) [Source:HGNC Symbol;Acc:HGNC:474]'
# https://www.ncbi.nlm.nih.gov/gene/277; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:475; http://www.uniprot.org/uniprot/P04745
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='277'),'hgnc_symbol'] = 'AMY1B'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='277'),'hgnc_id'] = 'HGNC:475'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='277'),'uniprotswissprot'] = 'P04745'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='277'),'description'] = 'amylase, alpha 1B (salivary) [Source:HGNC Symbol;Acc:HGNC:475]'
# https://www.ncbi.nlm.nih.gov/gene/278; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:476; http://www.uniprot.org/uniprot/P04745
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='278'),'hgnc_symbol'] = 'AMY1C'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='278'),'hgnc_id'] = 'HGNC:476'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='278'),'uniprotswissprot'] = 'P04745'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='278'),'description'] = 'amylase, alpha 1C (salivary) [Source:HGNC Symbol;Acc:HGNC:476]'
# https://www.ncbi.nlm.nih.gov/gene/279; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:477; http://www.uniprot.org/uniprot/P04746
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='279'),'hgnc_symbol'] = 'AMY2A'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='279'),'hgnc_id'] = 'HGNC:477'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='279'),'uniprotswissprot'] = 'P04746'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='279'),'description'] = 'amylase, alpha 2A (pancreatic) [Source:HGNC Symbol;Acc:HGNC:477]'
# https://www.ncbi.nlm.nih.gov/gene/280; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:478; http://www.uniprot.org/uniprot/P19961
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='280'),'hgnc_symbol'] = 'AMY2B'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='280'),'hgnc_id'] = 'HGNC:478'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='280'),'uniprotswissprot'] = 'P19961'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='280'),'description'] = 'amylase, alpha 2B (pancreatic) [Source:HGNC Symbol;Acc:HGNC:478]'

# https://www.ncbi.nlm.nih.gov/gene/801; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:1442; http://www.uniprot.org/uniprot/P0DP23
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='801'),'hgnc_symbol'] = 'CALM1'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='801'),'hgnc_id'] = 'HGNC:1442'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='801'),'uniprotswissprot'] = 'P0DP23'	# P62158 (CALM_HUMAN): Demerged into P0DP23 (CALM1_HUMAN), P0DP24 (CALM2_HUMAN) and P0DP25 (CALM3_HUMAN).
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='801'),'description'] = 'calmodulin 1 [Source:HGNC Symbol;Acc:HGNC:1442]'
# https://www.ncbi.nlm.nih.gov/gene/805; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:1445; http://www.uniprot.org/uniprot/P0DP24
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='805'),'hgnc_symbol'] = 'CALM2'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='805'),'hgnc_id'] = 'HGNC:1445'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='805'),'uniprotswissprot'] = 'P0DP24'	# P62158 (CALM_HUMAN): Demerged into P0DP23 (CALM1_HUMAN), P0DP24 (CALM2_HUMAN) and P0DP25 (CALM3_HUMAN).
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='805'),'description'] = 'calmodulin 2 [Source:HGNC Symbol;Acc:HGNC:1445]'
# https://www.ncbi.nlm.nih.gov/gene/808; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:1449; http://www.uniprot.org/uniprot/P0DP25
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='808'),'hgnc_symbol'] = 'CALM3'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='808'),'hgnc_id'] = 'HGNC:1449'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='808'),'uniprotswissprot'] = 'P0DP25'	# P62158 (CALM_HUMAN): Demerged into P0DP23 (CALM1_HUMAN), P0DP24 (CALM2_HUMAN) and P0DP25 (CALM3_HUMAN).
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='808'),'description'] = 'calmodulin 3 [Source:HGNC Symbol;Acc:HGNC:1449]'

# https://www.ncbi.nlm.nih.gov/gene/2678; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:4250; http://www.uniprot.org/uniprot/P19440
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='2678'),'hgnc_symbol'] = 'GGT1'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='2678'),'hgnc_id'] = 'HGNC:4250'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='2678'),'uniprotswissprot'] = 'P19440'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='2678'),'description'] = 'gamma-glutamyltransferase 1 [Source:HGNC Symbol;Acc:HGNC:4250]'
# https://www.ncbi.nlm.nih.gov/gene/728441; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:4251; http://www.uniprot.org/uniprot/P36268
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='728441'),'hgnc_symbol'] = 'GGT2'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='728441'),'hgnc_id'] = 'HGNC:4251'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='728441'),'uniprotswissprot'] = 'P36268'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='728441'),'description'] = 'gamma-glutamyltransferase 2 [Source:HGNC Symbol;Acc:HGNC:4251]'

# https://www.ncbi.nlm.nih.gov/gene/6606; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:11117; http://www.uniprot.org/uniprot/Q16637
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='6606'),'hgnc_symbol'] = 'SMN1'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='6606'),'hgnc_id'] = 'HGNC:11117'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='6606'),'uniprotswissprot'] = 'Q16637'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='6606'),'description'] = 'survival of motor neuron 1, telomeric [Source:HGNC Symbol;Acc:HGNC:11117]'
# https://www.ncbi.nlm.nih.gov/gene/6607; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:11118; http://www.uniprot.org/uniprot/Q16637
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='6607'),'hgnc_symbol'] = 'SMN2'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='6607'),'hgnc_id'] = 'HGNC:11118'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='6607'),'uniprotswissprot'] = 'Q16637'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='6607'),'description'] = 'survival of motor neuron 2, centromeric [Source:HGNC Symbol;Acc:HGNC:11118]'

# https://www.ncbi.nlm.nih.gov/gene/2212; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:3616; http://www.uniprot.org/uniprot/P12318
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='2212'),'hgnc_symbol'] = 'FCGR2A'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='2212'),'hgnc_id'] = 'HGNC:3616'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='2212'),'uniprotswissprot'] = 'P12318'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='2212'),'description'] = 'Fc fragment of IgG receptor IIa [Source:HGNC Symbol;Acc:HGNC:3616]'
# https://www.ncbi.nlm.nih.gov/gene/9103; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:15626; http://www.uniprot.org/uniprot/P31995
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='9103'),'hgnc_symbol'] = 'FCGR2C'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='9103'),'hgnc_id'] = 'HGNC:15626'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='9103'),'uniprotswissprot'] = 'P31995'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='9103'),'description'] = 'Fc fragment of IgG receptor IIc (gene/pseudogene) [Source:HGNC Symbol;Acc:HGNC:15626]'

# https://www.ncbi.nlm.nih.gov/gene/3115; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:4940; http://www.uniprot.org/uniprot/P04440
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='3115'),'hgnc_symbol'] = 'HLA-DPB1'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='3115'),'hgnc_id'] = 'HGNC:4940'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='3115'),'uniprotswissprot'] = 'P04440'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='3115'),'description'] = 'major histocompatibility complex, class II, DP beta 1 [Source:HGNC Symbol;Acc:HGNC:4940]'
# https://www.ncbi.nlm.nih.gov/gene/3117; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:4942; http://www.uniprot.org/uniprot/P01909
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='3117'),'hgnc_symbol'] = 'HLA-DQA1'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='3117'),'hgnc_id'] = 'HGNC:4942'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='3117'),'uniprotswissprot'] = 'P01909'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='3117'),'description'] = 'major histocompatibility complex, class II, DQ alpha 1 [Source:HGNC Symbol;Acc:HGNC:4942]'

# https://www.ncbi.nlm.nih.gov/gene/6232; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:10416; http://www.uniprot.org/uniprot/P42677
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='6232'),'hgnc_symbol'] = 'RPS27'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='6232'),'hgnc_id'] = 'HGNC:10416'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='6232'),'uniprotswissprot'] = 'P42677'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='6232'),'description'] = 'ribosomal protein S27 [Source:HGNC Symbol;Acc:HGNC:10416]'
# https://www.ncbi.nlm.nih.gov/gene/51065; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:18476; http://www.uniprot.org/uniprot/Q71UM5
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='51065'),'hgnc_symbol'] = 'RPS27L'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='51065'),'hgnc_id'] = 'HGNC:18476'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='51065'),'uniprotswissprot'] = 'Q71UM5'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='51065'),'description'] = 'ribosomal protein S27 like [Source:HGNC Symbol;Acc:HGNC:18476]'

# https://www.ncbi.nlm.nih.gov/gene/7307; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:12453; http://www.uniprot.org/uniprot/Q01081
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='7307'),'hgnc_symbol'] = 'U2AF1'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='7307'),'hgnc_id'] = 'HGNC:12453'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='7307'),'uniprotswissprot'] = 'Q01081'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='7307'),'description'] = 'U2 small nuclear RNA auxiliary factor 1 [Source:HGNC Symbol;Acc:HGNC:12453]'
# https://www.ncbi.nlm.nih.gov/gene/102724594; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:51830; http://www.uniprot.org/uniprot/P0DN76
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='102724594'),'hgnc_symbol'] = 'U2AF1L5'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='102724594'),'hgnc_id'] = 'HGNC:51830'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='102724594'),'uniprotswissprot'] = 'P0DN76'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='102724594'),'description'] = 'U2 small nuclear RNA auxiliary factor 1 like 5 [Source:HGNC Symbol;Acc:HGNC:51830]'

# https://www.ncbi.nlm.nih.gov/gene/552900; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:29488; http://www.uniprot.org/uniprot/Q9H3K6
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='552900'),'hgnc_symbol'] = 'BOLA2'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='552900'),'hgnc_id'] = 'HGNC:29488'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='552900'),'uniprotswissprot'] = 'Q9H3K6'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='552900'),'description'] = 'bolA family member 2 [Source:HGNC Symbol;Acc:HGNC:29488]'
# https://www.ncbi.nlm.nih.gov/gene/654483; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:32479; http://www.uniprot.org/uniprot/Q9H3K6
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='654483'),'hgnc_symbol'] = 'BOLA2B'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='654483'),'hgnc_id'] = 'HGNC:32479'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='654483'),'uniprotswissprot'] = 'Q9H3K6'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='654483'),'description'] = 'bolA family member 2B [Source:HGNC Symbol;Acc:HGNC:32479]'

# https://www.ncbi.nlm.nih.gov/gene/3963; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:6568; http://www.uniprot.org/uniprot/P47929
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='3963'),'hgnc_symbol'] = 'LGALS7'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='3963'),'hgnc_id'] = 'HGNC:6568'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='3963'),'uniprotswissprot'] = 'P47929'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='3963'),'description'] = 'galectin 7 [Source:HGNC Symbol;Acc:HGNC:6568]'
# https://www.ncbi.nlm.nih.gov/gene/653499; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:34447; http://www.uniprot.org/uniprot/P47929
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='653499'),'hgnc_symbol'] = 'LGALS7B'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='653499'),'hgnc_id'] = 'HGNC:34447'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='653499'),'uniprotswissprot'] = 'P47929'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='653499'),'description'] = 'galectin 7B [Source:HGNC Symbol;Acc:HGNC:34447]'

# https://www.ncbi.nlm.nih.gov/gene/200315; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:17343; http://www.uniprot.org/uniprot/P31941
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='200315'),'hgnc_symbol'] = 'APOBEC3A'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='200315'),'hgnc_id'] = 'HGNC:17343'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='200315'),'uniprotswissprot'] = 'P31941'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='200315'),'description'] = 'apolipoprotein B mRNA editing enzyme catalytic subunit 3A [Source:HGNC Symbol;Acc:HGNC:17343]'
# https://www.ncbi.nlm.nih.gov/gene/100913187; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:44196; http://www.uniprot.org/uniprot/P31941
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='100913187'),'hgnc_symbol'] = 'APOBEC3A_B'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='100913187'),'hgnc_id'] = 'HGNC:44196'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='100913187'),'uniprotswissprot'] = 'P31941'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='100913187'),'description'] = 'APOBEC3A and APOBEC3B deletion hybrid [Source:HGNC Symbol;Acc:HGNC:44196]'

# https://www.ncbi.nlm.nih.gov/gene/875; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:1550; http://www.uniprot.org/uniprot/P35520
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='875'),'hgnc_symbol'] = 'CBS'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='875'),'hgnc_id'] = 'HGNC:1550'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='875'),'uniprotswissprot'] = 'P35520'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='875'),'description'] = 'cystathionine-beta-synthase [Source:HGNC Symbol;Acc:HGNC:1550]'
# https://www.ncbi.nlm.nih.gov/gene/102724560; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:51829; http://www.uniprot.org/uniprot/P0DN79
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='102724560'),'hgnc_symbol'] = 'CBSL'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='102724560'),'hgnc_id'] = 'HGNC:51829'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='102724560'),'uniprotswissprot'] = 'P0DN79'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='102724560'),'description'] = 'cystathionine-beta-synthase like [Source:HGNC Symbol;Acc:HGNC:51829]'

# https://www.ncbi.nlm.nih.gov/gene/6231; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:10414; http://www.uniprot.org/uniprot/P62854
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='6231'),'hgnc_symbol'] = 'RPS26'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='6231'),'hgnc_id'] = 'HGNC:10414'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='6231'),'uniprotswissprot'] = 'P62854'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='6231'),'description'] = 'ribosomal protein S26 [Source:HGNC Symbol;Acc:HGNC:10414]'
# https://www.ncbi.nlm.nih.gov/gene/101929876; http://www.uniprot.org/uniprot/P62854
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='101929876'),'hgnc_symbol'] = 'LOC101929876'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='101929876'),'hgnc_id'] = ''
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='101929876'),'uniprotswissprot'] = 'P62854'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='101929876'),'description'] = '40S ribosomal protein S26'

# https://www.ncbi.nlm.nih.gov/gene/23042; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:28995; http://www.uniprot.org/uniprot/Q6P996
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='23042'),'hgnc_symbol'] = 'PDXDC1'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='23042'),'hgnc_id'] = 'HGNC:28995'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='23042'),'uniprotswissprot'] = 'Q6P996'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='23042'),'description'] = 'pyridoxal dependent decarboxylase domain containing 1 [Source:HGNC Symbol;Acc:HGNC:28995]'
# https://www.ncbi.nlm.nih.gov/gene/102724985; http://www.uniprot.org/uniprot/Q6P996
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='102724985'),'hgnc_symbol'] = 'LOC102724985'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='102724985'),'hgnc_id'] = ''
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='102724985'),'uniprotswissprot'] = 'Q6P996'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='102724985'),'description'] = 'pyridoxal-dependent decarboxylase domain-containing protein 1'

# https://www.ncbi.nlm.nih.gov/gene/55718; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:30347; http://www.uniprot.org/uniprot/Q9NVU0
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='55718'),'hgnc_symbol'] = 'POLR3E'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='55718'),'hgnc_id'] = 'HGNC:30347'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='55718'),'uniprotswissprot'] = 'Q9NVU0'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='55718'),'description'] = 'RNA polymerase III subunit E [Source:HGNC Symbol;Acc:HGNC:30347]'
# https://www.ncbi.nlm.nih.gov/gene/101060521; http://www.uniprot.org/uniprot/Q9NVU0
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='101060521'),'hgnc_symbol'] = 'LOC101060521'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='101060521'),'hgnc_id'] = ''
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='101060521'),'uniprotswissprot'] = 'Q9NVU0'
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='101060521'),'description'] = 'DNA-directed RNA polymerase III subunit RPC5'

ExoCarta_Protein.hgnc_id.table[,'description'] = trimws(ExoCarta_Protein.hgnc_id.table[,'description'])
ExoCarta_Protein.hgnc_id.table[,'description'] = sapply(strsplit(ExoCarta_Protein.hgnc_id.table[,'description'],'\\|'), FUN=function(U) paste(unique(sapply(strsplit(unlist(U),'\\[Source:'), FUN=function(V) trimws(unlist(V))[1])),collapse="|"))
ExoCarta_Protein.hgnc_id.table = unique(ExoCarta_Protein.hgnc_id.table)

ExoCarta_Protein.hgnc_id.table = ExoCarta_Protein.hgnc_id.table[,c('hgnc_symbol','hgnc_id','entrezgene','uniprotswissprot','description')]
ExoCarta_Protein.hgnc_id.table = ExoCarta_Protein.hgnc_id.table[order(match(ExoCarta_Protein.hgnc_id.table[,'entrezgene'], a)),]

# merge multiple uniprotswissprot of single gene
unique.hgnc_id.idx <- which(!duplicated(ExoCarta_Protein.hgnc_id.table[,'hgnc_id']))
unique.hgnc_id <- ExoCarta_Protein.hgnc_id.table[unique.hgnc_id.idx,'hgnc_id']
duplicated.hgnc_id.idx <- which(duplicated(ExoCarta_Protein.hgnc_id.table[,'hgnc_id']))
duplicated.hgnc_id <- ExoCarta_Protein.hgnc_id.table[duplicated.hgnc_id.idx,'hgnc_id']
duplicated.hgnc_id.idx.pair <- unique(sapply(duplicated.hgnc_id, FUN=function(X) which(ExoCarta_Protein.hgnc_id.table[,'hgnc_id'] %in% X)))
duplicated.hgnc_id.idx.uniprotswissprot <- sapply(duplicated.hgnc_id.idx.pair, FUN=function(X) paste(sort(unique(ExoCarta_Protein.hgnc_id.table[X,'uniprotswissprot'][ExoCarta_Protein.hgnc_id.table[X,'uniprotswissprot']!=''])),collapse=';'))
if(length(duplicated.hgnc_id.idx.pair) > 0)
{
	for(i in 1:length(duplicated.hgnc_id.idx.pair))
	{
		ExoCarta_Protein.hgnc_id.table[unlist(duplicated.hgnc_id.idx.pair[i]),'uniprotswissprot'] = duplicated.hgnc_id.idx.uniprotswissprot[i]
	}
}

# uniprotswissprot not found by R package "biomaRt"
empty.uniprotswissprot.idx <- which(ExoCarta_Protein.hgnc_id.table[,'uniprotswissprot']=='')
empty.uniprotswissprot <- ExoCarta_Protein.hgnc_id.table[empty.uniprotswissprot.idx,]

ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='820'),'uniprotswissprot'] = 'P49913'				# https://www.ncbi.nlm.nih.gov/gene/820			https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:1472	http://www.uniprot.org/uniprot/P49913
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='828'),'uniprotswissprot'] = 'Q13938'				# https://www.ncbi.nlm.nih.gov/gene/828			https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:1487	http://www.uniprot.org/uniprot/Q13938
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='984'),'uniprotswissprot'] = 'P21127'				# https://www.ncbi.nlm.nih.gov/gene/984			https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:1729	http://www.uniprot.org/uniprot/P21127
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='1056'),'uniprotswissprot'] = 'P19835'				# https://www.ncbi.nlm.nih.gov/gene/1056		https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:1848	http://www.uniprot.org/uniprot/P19835
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='1731'),'uniprotswissprot'] = 'Q8WYJ6'				# https://www.ncbi.nlm.nih.gov/gene/1731		https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:2879	http://www.uniprot.org/uniprot/Q8WYJ6
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='4925'),'uniprotswissprot'] = 'P80303'				# https://www.ncbi.nlm.nih.gov/gene/4925		https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:8044	http://www.uniprot.org/uniprot/P80303
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='5175'),'uniprotswissprot'] = 'P16284'				# https://www.ncbi.nlm.nih.gov/gene/5175		https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:8823	http://www.uniprot.org/uniprot/P16284
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='6854'),'uniprotswissprot'] = 'Q92777'				# https://www.ncbi.nlm.nih.gov/gene/6854		https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:11495	http://www.uniprot.org/uniprot/Q92777
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='7873'),'uniprotswissprot'] = 'P55145'				# https://www.ncbi.nlm.nih.gov/gene/7873		https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:15461	http://www.uniprot.org/uniprot/P55145
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='10189'),'uniprotswissprot'] = 'Q86V81'				# https://www.ncbi.nlm.nih.gov/gene/10189		https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:19071	http://www.uniprot.org/uniprot/Q86V81
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='11117'),'uniprotswissprot'] = 'Q9Y6C2'				# https://www.ncbi.nlm.nih.gov/gene/11117		https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:19880	http://www.uniprot.org/uniprot/Q9Y6C2
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='23145'),'uniprotswissprot'] = 'A2VEC9'				# https://www.ncbi.nlm.nih.gov/gene/23145		https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:21998	http://www.uniprot.org/uniprot/A2VEC9
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='23380'),'uniprotswissprot'] = 'O75044'				# https://www.ncbi.nlm.nih.gov/gene/23380		https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:19751	http://www.uniprot.org/uniprot/O75044
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='51411'),'uniprotswissprot'] = 'Q9UBW5'				# https://www.ncbi.nlm.nih.gov/gene/51411		https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:1053	http://www.uniprot.org/uniprot/Q9UBW5
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='51497'),'uniprotswissprot'] = 'Q8IXH7'				# https://www.ncbi.nlm.nih.gov/gene/51497		https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:15934	http://www.uniprot.org/uniprot/Q8IXH7
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='51667'),'uniprotswissprot'] = 'Q9Y5A7'				# https://www.ncbi.nlm.nih.gov/gene/51667		https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:17623	http://www.uniprot.org/uniprot/Q9Y5A7
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='54798'),'uniprotswissprot'] = 'Q6V1P9'				# https://www.ncbi.nlm.nih.gov/gene/54798		https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:23111	http://www.uniprot.org/uniprot/Q6V1P9
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='55742'),'uniprotswissprot'] = 'Q9NVD7'				# https://www.ncbi.nlm.nih.gov/gene/55742		https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:14652	http://www.uniprot.org/uniprot/Q9NVD7
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='80852'),'uniprotswissprot'] = 'Q9C0E4'				# https://www.ncbi.nlm.nih.gov/gene/80852		https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:23841	http://www.uniprot.org/uniprot/Q9C0E4
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='83481'),'uniprotswissprot'] = 'P58107'				# https://www.ncbi.nlm.nih.gov/gene/83481		https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:15577	http://www.uniprot.org/uniprot/P58107
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='83882'),'uniprotswissprot'] = 'Q9H1Z9'				# https://www.ncbi.nlm.nih.gov/gene/83882		https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:29942	http://www.uniprot.org/uniprot/Q9H1Z9
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='283450'),'uniprotswissprot'] = 'Q9Y4D8'			# https://www.ncbi.nlm.nih.gov/gene/283450		https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:26611	http://www.uniprot.org/uniprot/Q9Y4D8
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='375791'),'uniprotswissprot'] = 'A8MQ03'			# https://www.ncbi.nlm.nih.gov/gene/375791		https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:30529	http://www.uniprot.org/uniprot/A8MQ03
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='389898'),'uniprotswissprot'] = 'Q5JXB2'			# https://www.ncbi.nlm.nih.gov/gene/389898		https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:31710	http://www.uniprot.org/uniprot/Q5JXB2
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='654364'),'uniprotswissprot'] = 'P22392'			# https://www.ncbi.nlm.nih.gov/gene/654364		https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:33531	http://www.uniprot.org/uniprot/P22392
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='100133941'),'uniprotswissprot'] = 'P25063'			# https://www.ncbi.nlm.nih.gov/gene/100133941	https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:1645	http://www.uniprot.org/uniprot/P25063
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='100526664'),'uniprotswissprot'] = 'O60449'			# https://www.ncbi.nlm.nih.gov/gene/100526664	https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:38828	http://www.uniprot.org/uniprot/O60449
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='100526693'),'uniprotswissprot'] = 'A0A0A6YYG9'		# https://www.ncbi.nlm.nih.gov/gene/100526693	https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:38830	http://www.uniprot.org/uniprot/A0A0A6YYG9
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='100526737'),'uniprotswissprot'] = 'Q96PK6'			# https://www.ncbi.nlm.nih.gov/gene/100526737	https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:38840	http://www.uniprot.org/uniprot/Q96PK6
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='100526767'),'uniprotswissprot'] = 'Q9Y3E7'			# https://www.ncbi.nlm.nih.gov/gene/100526767	https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:38847	http://www.uniprot.org/uniprot/Q9Y3E7
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='100526842'),'uniprotswissprot'] = 'A0A0A0MRF8'		# https://www.ncbi.nlm.nih.gov/gene/100526842	https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:44661	http://www.uniprot.org/uniprot/A0A0A0MRF8
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='100529144'),'uniprotswissprot'] = 'A0A0A6YYL4'		# https://www.ncbi.nlm.nih.gov/gene/100529144	https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:44424	http://www.uniprot.org/uniprot/A0A0A6YYL4
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='100534599'),'uniprotswissprot'] = 'Q9ULR0'			# https://www.ncbi.nlm.nih.gov/gene/100534599	https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:42969	http://www.uniprot.org/uniprot/Q9ULR0
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='387522'),'uniprotswissprot'] = 'A5PLL7;Q13404'		# https://www.ncbi.nlm.nih.gov/gene/387522		https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:33521	http://www.uniprot.org/uniprot/A5PLL7 + http://www.uniprot.org/uniprot/Q13404
ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']=='100302736'),'uniprotswissprot'] = 'Q9Y3B3;Q86XR7'	# https://www.ncbi.nlm.nih.gov/gene/100302736	https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:33945	http://www.uniprot.org/uniprot/Q9Y3B3 + http://www.uniprot.org/uniprot/Q86XR7

ExoCarta_Protein.hgnc_id.table = unique(ExoCarta_Protein.hgnc_id.table)
#dim(ExoCarta_Protein.hgnc_id.table)

# remove MIR3620 (microRNA 3620) (https://www.ncbi.nlm.nih.gov/gene/100500810)
ExoCarta_Protein.hgnc_id.table = ExoCarta_Protein.hgnc_id.table[which(ExoCarta_Protein.hgnc_id.table[,'entrezgene']!='100500810'),]

rownames(ExoCarta_Protein.hgnc_id.table) = c(1:nrow(ExoCarta_Protein.hgnc_id.table))
ExoCarta_Protein.hgnc_id.table.New = ExoCarta_Protein.hgnc_id.table
#ExoCarta_Protein.hgnc_id.table
#ExoCarta_Protein.hgnc_id.table.New

#human.UniProt.Entry2ProteinNames = sapply(strsplit(ExoCarta_Protein.hgnc_id.table[,'uniprotswissprot'],";"), FUN=function(X) paste(unique(sapply(unlist(X), FUN=function(Y) ifelse(length(which(human.UniProt.reviewed.Entry %in% unlist(Y)))==0, ifelse(length(which(human.UniProt.unreviewed.Entry %in% unlist(Y)))==0, '', human.UniProt.unreviewed.ProteinNames[which(human.UniProt.unreviewed.Entry %in% unlist(Y))]), human.UniProt.reviewed.ProteinNames[which(human.UniProt.reviewed.Entry %in% unlist(Y))]))),collapse="|"))
#ExoCarta_Protein.hgnc_id.table = cbind(ExoCarta_Protein.hgnc_id.table, human.UniProt.Entry2ProteinNames)
#colnames(ExoCarta_Protein.hgnc_id.table)[ncol(ExoCarta_Protein.hgnc_id.table)] = 'Protein_names'

#write.csv(ExoCarta_Protein.hgnc_id.table, paste('./Data/',Sample,'/10_datasetOverlap/Database/hsa/ExoCarta_All.Protein_Gene2HGNC.table.csv',sep=""), quote=TRUE, row.names=TRUE)

# entrezgene: attributes not found by R package "biomaRt"
ExoCarta_Protein.entrezgene.NA_attributes <- a[which(! a %in% ExoCarta_Protein.hgnc_id.table[,'entrezgene'])]
if(length(ExoCarta_Protein.entrezgene.NA_attributes) > 0)
{
	manually_curated.hgnc_id.table <- matrix('',length(ExoCarta_Protein.entrezgene.NA_attributes),ncol(ExoCarta_Protein.hgnc_id.table))
	colnames(manually_curated.hgnc_id.table) = colnames(ExoCarta_Protein.hgnc_id.table)
	
	ExoCarta_Protein.entrezgene.NA_attributes.uncurated <- c()
	
	for(i in 1:length(ExoCarta_Protein.entrezgene.NA_attributes))
	{
		# --- Gene type: protein coding -----------------------------------------------------------
		if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '3492') {				# https://www.ncbi.nlm.nih.gov/gene/3492; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:5477
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'IGH'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:5477'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '3492'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'immunoglobulin heavy locus [Source:HGNC Symbol;Acc:HGNC:5477]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '3535') {			# https://www.ncbi.nlm.nih.gov/gene/3535; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:5853
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'IGL'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:5853'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '3535'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'immunoglobulin lambda locus [Source:HGNC Symbol;Acc:HGNC:5853]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '50802') {		# https://www.ncbi.nlm.nih.gov/gene/50802; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:5715
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'IGK'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:5715'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '50802'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'immunoglobulin kappa locus [Source:HGNC Symbol;Acc:HGNC:5715]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '388524') {		# https://www.ncbi.nlm.nih.gov/gene/388524; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:36809
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RPSAP58'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:36809'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '388524'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'ribosomal protein SA pseudogene 58 [Source:HGNC Symbol;Acc:HGNC:36809]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '567') {			# https://www.ncbi.nlm.nih.gov/gene/567; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:914; http://www.uniprot.org/uniprot/P61769
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'B2M'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:914'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '567'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'P61769'
			manually_curated.hgnc_id.table[i,'description'] = 'beta-2-microglobulin [Source:HGNC Symbol;Acc:HGNC:914]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '2219') {			# https://www.ncbi.nlm.nih.gov/gene/2219; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:3623; http://www.uniprot.org/uniprot/O00602
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'FCN1'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:3623'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '2219'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'O00602'
			manually_curated.hgnc_id.table[i,'description'] = 'ficolin 1 [Source:HGNC Symbol;Acc:HGNC:3623]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '2632') {			# https://www.ncbi.nlm.nih.gov/gene/2632; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:4180; http://www.uniprot.org/uniprot/Q04446
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'GBE1'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:4180'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '2632'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'Q04446'
			manually_curated.hgnc_id.table[i,'description'] = '1,4-alpha-glucan branching enzyme 1 [Source:HGNC Symbol;Acc:HGNC:4180]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '2781') {			# https://www.ncbi.nlm.nih.gov/gene/2781; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:4395; http://www.uniprot.org/uniprot/P19086
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'GNAZ'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:4395'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '2781'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'P19086'
			manually_curated.hgnc_id.table[i,'description'] = 'G protein subunit alpha z [Source:HGNC Symbol;Acc:HGNC:4395]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '3866') {			# https://www.ncbi.nlm.nih.gov/gene/3866; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:6421; http://www.uniprot.org/uniprot/P19012
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'KRT15'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:6421'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '3866'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'P19012'
			manually_curated.hgnc_id.table[i,'description'] = 'keratin 15 [Source:HGNC Symbol;Acc:HGNC:6421]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '3934') {			# https://www.ncbi.nlm.nih.gov/gene/3934; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:6526; http://www.uniprot.org/uniprot/P80188
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'LCN2'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:6526'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '3934'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'P80188'
			manually_curated.hgnc_id.table[i,'description'] = 'lipocalin 2 [Source:HGNC Symbol;Acc:HGNC:6526]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '5430') {			# https://www.ncbi.nlm.nih.gov/gene/5430; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:9187; http://www.uniprot.org/uniprot/P24928
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'POLR2A'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:9187'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '5430'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'P24928'
			manually_curated.hgnc_id.table[i,'description'] = 'RNA polymerase II subunit A [Source:HGNC Symbol;Acc:HGNC:9187]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '6125') {			# https://www.ncbi.nlm.nih.gov/gene/6125; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:10360; http://www.uniprot.org/uniprot/P46777
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RPL5'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:10360'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '6125'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'P46777'
			manually_curated.hgnc_id.table[i,'description'] = 'ribosomal protein L5 [Source:HGNC Symbol;Acc:HGNC:10360]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '6202') {			# https://www.ncbi.nlm.nih.gov/gene/6202; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:10441; http://www.uniprot.org/uniprot/P62241
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RPS8'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:10441'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '6202'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'P62241'
			manually_curated.hgnc_id.table[i,'description'] = 'ribosomal protein S8 [Source:HGNC Symbol;Acc:HGNC:10441]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '6207') {			# https://www.ncbi.nlm.nih.gov/gene/6207; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:10386; http://www.uniprot.org/uniprot/P62277
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RPS13'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:10386'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '6207'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'P62277'
			manually_curated.hgnc_id.table[i,'description'] = 'ribosomal protein S13 [Source:HGNC Symbol;Acc:HGNC:10386]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '6627') {			# https://www.ncbi.nlm.nih.gov/gene/6627; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:11152; http://www.uniprot.org/uniprot/P09661
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'SNRPA1'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:11152'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '6627'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'P09661'
			manually_curated.hgnc_id.table[i,'description'] = 'small nuclear ribonucleoprotein polypeptide A\' [Source:HGNC Symbol;Acc:HGNC:11152]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '8351') {			# https://www.ncbi.nlm.nih.gov/gene/8351; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:4767; http://www.uniprot.org/uniprot/P68431
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'HIST1H3D'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:4767'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '8351'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'P68431'
			manually_curated.hgnc_id.table[i,'description'] = 'histone cluster 1 H3 family member d [Source:HGNC Symbol;Acc:HGNC:4767]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '8411') {			# https://www.ncbi.nlm.nih.gov/gene/8411; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:3185; http://www.uniprot.org/uniprot/Q15075
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'EEA1'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:3185'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '8411'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'Q15075'
			manually_curated.hgnc_id.table[i,'description'] = 'early endosome antigen 1 [Source:HGNC Symbol;Acc:HGNC:3185]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '8519') {			# https://www.ncbi.nlm.nih.gov/gene/8519; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:5412; http://www.uniprot.org/uniprot/P13164
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'IFITM1'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:5412'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '8519'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'P13164'
			manually_curated.hgnc_id.table[i,'description'] = 'interferon induced transmembrane protein 1 [Source:HGNC Symbol;Acc:HGNC:5412]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '9573') {			# https://www.ncbi.nlm.nih.gov/gene/9573; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:4218; http://www.uniprot.org/uniprot/Q9NR23
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'GDF3'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:4218'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '9573'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'Q9NR23'
			manually_curated.hgnc_id.table[i,'description'] = 'growth differentiation factor 3 [Source:HGNC Symbol;Acc:HGNC:4218]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '11054') {		# https://www.ncbi.nlm.nih.gov/gene/11054; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:15768; http://www.uniprot.org/uniprot/Q9NZT2
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'OGFR'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:15768'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '11054'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'Q9NZT2'
			manually_curated.hgnc_id.table[i,'description'] = 'opioid growth factor receptor [Source:HGNC Symbol;Acc:HGNC:15768]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '55766') {		# https://www.ncbi.nlm.nih.gov/gene/55766; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:14456; http://www.uniprot.org/uniprot/Q9BTM1
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'H2AFJ'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:14456'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '55766'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'Q9BTM1'
			manually_curated.hgnc_id.table[i,'description'] = 'H2A histone family member J [Source:HGNC Symbol;Acc:HGNC:14456]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '57535') {		# https://www.ncbi.nlm.nih.gov/gene/57535; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:29618; http://www.uniprot.org/uniprot/Q6UXG2
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'KIAA1324'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:29618'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '57535'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'Q6UXG2'
			manually_curated.hgnc_id.table[i,'description'] = 'KIAA1324 [Source:HGNC Symbol;Acc:HGNC:29618]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '80086') {		# https://www.ncbi.nlm.nih.gov/gene/80086; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:18637; http://www.uniprot.org/uniprot/Q9H853
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'TUBA4B'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:18637'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '80086'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'Q9H853'
			manually_curated.hgnc_id.table[i,'description'] = 'tubulin alpha 4b [Source:HGNC Symbol;Acc:HGNC:18637]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '84876') {		# https://www.ncbi.nlm.nih.gov/gene/84876; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:25896; http://www.uniprot.org/uniprot/Q96D31
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'ORAI1'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:25896'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '84876'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'Q96D31'
			manually_curated.hgnc_id.table[i,'description'] = 'ORAI calcium release-activated calcium modulator 1 [Source:HGNC Symbol;Acc:HGNC:25896]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '113189') {		# https://www.ncbi.nlm.nih.gov/gene/113189; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:24464; http://www.uniprot.org/uniprot/Q8NCH0
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'CHST14'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:24464'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '113189'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'Q8NCH0'
			manually_curated.hgnc_id.table[i,'description'] = 'carbohydrate sulfotransferase 14 [Source:HGNC Symbol;Acc:HGNC:24464]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '114827') {		# https://www.ncbi.nlm.nih.gov/gene/114827; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:29408; http://www.uniprot.org/uniprot/B1AJZ9
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'FHAD1'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:29408'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '114827'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'B1AJZ9'
			manually_curated.hgnc_id.table[i,'description'] = 'forkhead associated phosphopeptide binding domain 1 [Source:HGNC Symbol;Acc:HGNC:29408]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '119467') {		# https://www.ncbi.nlm.nih.gov/gene/119467; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:20795; http://www.uniprot.org/uniprot/Q8NCR9
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'CLRN3'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:20795'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '119467'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'Q8NCR9'
			manually_curated.hgnc_id.table[i,'description'] = 'clarin 3 [Source:HGNC Symbol;Acc:HGNC:20795]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '151790') {		# https://www.ncbi.nlm.nih.gov/gene/151790; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:26587; http://www.uniprot.org/uniprot/Q8IV35
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'WDR49'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:26587'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '151790'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'Q8IV35'
			manually_curated.hgnc_id.table[i,'description'] = 'WD repeat domain 49 [Source:HGNC Symbol;Acc:HGNC:26587]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '246181') {		# https://www.ncbi.nlm.nih.gov/gene/246181; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:24056; http://www.uniprot.org/uniprot/Q8NHP1
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'AKR7L'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:24056'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '246181'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'Q8NHP1'
			manually_curated.hgnc_id.table[i,'description'] = 'aldo-keto reductase family 7 like (gene/pseudogene) [Source:HGNC Symbol;Acc:HGNC:24056]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '375248') {		# https://www.ncbi.nlm.nih.gov/gene/375248; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:24079; http://www.uniprot.org/uniprot/A6QL64
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'ANKRD36'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:24079'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '375248'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'A6QL64'
			manually_curated.hgnc_id.table[i,'description'] = 'ankyrin repeat domain 36 [Source:HGNC Symbol;Acc:HGNC:24079]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '414061') {		# https://www.ncbi.nlm.nih.gov/gene/414061; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:32397; http://www.uniprot.org/uniprot/Q8WWF6
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'DNAJB3'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:32397'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '414061'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'Q8WWF6'
			manually_curated.hgnc_id.table[i,'description'] = 'DnaJ heat shock protein family (Hsp40) member B3 [Source:HGNC Symbol;Acc:HGNC:32397]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '441459') {		# https://www.ncbi.nlm.nih.gov/gene/441459; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:23644; http://www.uniprot.org/uniprot/A2A2Z9
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'ANKRD18B'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:23644'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '441459'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'A2A2Z9'
			manually_curated.hgnc_id.table[i,'description'] = 'ankyrin repeat domain 18B [Source:HGNC Symbol;Acc:HGNC:23644]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '100996720') {	# https://www.ncbi.nlm.nih.gov/gene/100996720
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'LOC100996720'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = ''
			manually_curated.hgnc_id.table[i,'entrezgene'] = '100996720'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'uncharacterized LOC100996720'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '102724187') {	# https://www.ncbi.nlm.nih.gov/gene/102724187 -> https://www.ncbi.nlm.nih.gov/gene/100134938; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:37278; http://www.uniprot.org/uniprot/B0FP48
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'UPK3BL1'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:37278'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '100134938'				# discontinued on 1-Jul-2015: This record was replaced with Gene ID: 100134938
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'B0FP48'
			manually_curated.hgnc_id.table[i,'description'] = 'uroplakin 3B like 1 [Source:HGNC Symbol;Acc:HGNC:37278]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '642441') {		# https://www.ncbi.nlm.nih.gov/gene/642441 -> discontinued on 9-Jun-2016: This record has been withdrawn by NCBI because the model on which it was based was not predicted in a later annotation.
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'LOC642441'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = ''
			manually_curated.hgnc_id.table[i,'entrezgene'] = '642441'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'uncharacterized LOC642441'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '101060521') {	# https://www.ncbi.nlm.nih.gov/gene/101060521
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'LOC101060521'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = ''
			manually_curated.hgnc_id.table[i,'entrezgene'] = '101060521'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'Q9NVU0'
			manually_curated.hgnc_id.table[i,'description'] = 'DNA-directed RNA polymerase III subunit RPC5'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '101929876') {	# https://www.ncbi.nlm.nih.gov/gene/101929876
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'LOC101929876'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = ''
			manually_curated.hgnc_id.table[i,'entrezgene'] = '101929876'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'P62854'
			manually_curated.hgnc_id.table[i,'description'] = '40S ribosomal protein S26'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '102724844') {	# https://www.ncbi.nlm.nih.gov/gene/102724844 -> discontinued on 9-Jun-2016: This record has been withdrawn by NCBI because the model on which it was based was not predicted in a later annotation.
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'LOC102724844'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = ''
			manually_curated.hgnc_id.table[i,'entrezgene'] = '102724844'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'immunoglobulin superfamily member 3-like'
		# --- Gene type: other --------------------------------------------------------------------
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '440786') {		# https://www.ncbi.nlm.nih.gov/gene/440786
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'LOC440786'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = ''
			manually_curated.hgnc_id.table[i,'entrezgene'] = '440786'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'uncharacterized LOC440786'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '652070') {		# https://www.ncbi.nlm.nih.gov/gene/652070
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'SCFV'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = ''
			manually_curated.hgnc_id.table[i,'entrezgene'] = '652070'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'single-chain Fv fragment'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '388560') {		# https://www.ncbi.nlm.nih.gov/gene/388560; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:35201
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'ERVFRD-2'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:35201'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '388560'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'endogenous retrovirus group FRD member 2 [Source:HGNC Symbol;Acc:HGNC:35201]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '3493') {			# https://www.ncbi.nlm.nih.gov/gene/3493; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:5478; http://www.uniprot.org/uniprot/P01876
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'IGHA1'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:5478'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '3493'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'P01876'
			manually_curated.hgnc_id.table[i,'description'] = 'immunoglobulin heavy constant alpha 1 [Source:HGNC Symbol;Acc:HGNC:5478]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '3494') {			# https://www.ncbi.nlm.nih.gov/gene/3494; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:5479; http://www.uniprot.org/uniprot/P01877
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'IGHA2'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:5479'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '3494'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'P01877'
			manually_curated.hgnc_id.table[i,'description'] = 'immunoglobulin heavy constant alpha 2 (A2m marker) [Source:HGNC Symbol;Acc:HGNC:5479]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '3500') {			# https://www.ncbi.nlm.nih.gov/gene/3500; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:5525; http://www.uniprot.org/uniprot/P01857
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'IGHG1'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:5525'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '3500'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'P01857'
			manually_curated.hgnc_id.table[i,'description'] = 'immunoglobulin heavy constant gamma 1 (G1m marker) [Source:HGNC Symbol;Acc:HGNC:5525]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '3501') {			# https://www.ncbi.nlm.nih.gov/gene/3501; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:5526; http://www.uniprot.org/uniprot/P01859
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'IGHG2'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:5526'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '3501'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'P01859'
			manually_curated.hgnc_id.table[i,'description'] = 'immunoglobulin heavy constant gamma 2 (G2m marker) [Source:HGNC Symbol;Acc:HGNC:5526]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '3502') {			# https://www.ncbi.nlm.nih.gov/gene/3502; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:5527; http://www.uniprot.org/uniprot/P01860
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'IGHG3'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:5527'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '3502'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'P01860'
			manually_curated.hgnc_id.table[i,'description'] = 'immunoglobulin heavy constant gamma 3 (G3m marker) [Source:HGNC Symbol;Acc:HGNC:5527]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '3503') {			# https://www.ncbi.nlm.nih.gov/gene/3503; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:5528; http://www.uniprot.org/uniprot/P01861
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'IGHG4'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:5528'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '3503'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'P01861'
			manually_curated.hgnc_id.table[i,'description'] = 'immunoglobulin heavy constant gamma 4 (G4m marker) [Source:HGNC Symbol;Acc:HGNC:5528]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '3507') {			# https://www.ncbi.nlm.nih.gov/gene/3507; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:5541; http://www.uniprot.org/uniprot/P01871
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'IGHM'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:5541'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '3507'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'P01871'
			manually_curated.hgnc_id.table[i,'description'] = 'immunoglobulin heavy constant mu [Source:HGNC Symbol;Acc:HGNC:5541]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '3537') {			# https://www.ncbi.nlm.nih.gov/gene/3537; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:5855; http://www.uniprot.org/uniprot/P0CG04
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'IGLC1'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:5855'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '3537'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'P0CG04'
			manually_curated.hgnc_id.table[i,'description'] = 'immunoglobulin lambda constant 1 [Source:HGNC Symbol;Acc:HGNC:5855]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '3538') {			# https://www.ncbi.nlm.nih.gov/gene/3538; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:5856; http://www.uniprot.org/uniprot/P0DOY2
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'IGLC2'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:5856'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '3538'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'P0DOY2'
			manually_curated.hgnc_id.table[i,'description'] = 'immunoglobulin lambda constant 2 [Source:HGNC Symbol;Acc:HGNC:5856]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '3539') {			# https://www.ncbi.nlm.nih.gov/gene/3539; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:5857; http://www.uniprot.org/uniprot/P0DOY3
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'IGLC3'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:5857'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '3539'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'P0DOY3'
			manually_curated.hgnc_id.table[i,'description'] = 'immunoglobulin lambda constant 3 (Kern-Oz+ marker) [Source:HGNC Symbol;Acc:HGNC:5857]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '28299') {		# https://www.ncbi.nlm.nih.gov/gene/28299; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:5741; http://www.uniprot.org/uniprot/P01602
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'IGKV1-5'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:5741'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '28299'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'P01602'
			manually_curated.hgnc_id.table[i,'description'] = 'immunoglobulin kappa variable 1-5 [Source:HGNC Symbol;Acc:HGNC:5741]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '28450') {		# https://www.ncbi.nlm.nih.gov/gene/28450; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:5580; http://www.uniprot.org/uniprot/P01762
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'IGHV3-11'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:5580'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '28450'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'P01762'
			manually_curated.hgnc_id.table[i,'description'] = 'immunoglobulin heavy variable 3-11 (gene/pseudogene) [Source:HGNC Symbol;Acc:HGNC:5580]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '28452') {		# https://www.ncbi.nlm.nih.gov/gene/28452; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:5620; http://www.uniprot.org/uniprot/P01780
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'IGHV3-7'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:5620'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '28452'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'P01780'
			manually_curated.hgnc_id.table[i,'description'] = 'immunoglobulin heavy variable 3-7 [Source:HGNC Symbol;Acc:HGNC:5620]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '28786') {		# https://www.ncbi.nlm.nih.gov/gene/28786; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:5919; http://www.uniprot.org/uniprot/A0A075B6K6
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'IGLV4-3'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:5919'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '28786'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'A0A075B6K6'
			manually_curated.hgnc_id.table[i,'description'] = 'immunoglobulin lambda variable 4-3 [Source:HGNC Symbol;Acc:HGNC:5919]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '28796') {		# https://www.ncbi.nlm.nih.gov/gene/28796; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:5905; http://www.uniprot.org/uniprot/P80748
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'IGLV3-21'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:5905'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '28796'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'P80748'
			manually_curated.hgnc_id.table[i,'description'] = 'immunoglobulin lambda variable 3-21 [Source:HGNC Symbol;Acc:HGNC:5905]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '28816') {		# https://www.ncbi.nlm.nih.gov/gene/28816; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:5887; http://www.uniprot.org/uniprot/P01706
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'IGLV2-11'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:5887'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '28816'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'P01706'
			manually_curated.hgnc_id.table[i,'description'] = 'immunoglobulin lambda variable 2-11 [Source:HGNC Symbol;Acc:HGNC:5887]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '28823') {		# https://www.ncbi.nlm.nih.gov/gene/28823; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:5879; http://www.uniprot.org/uniprot/P01699
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'IGLV1-44'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:5879'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '28823'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'P01699'
			manually_curated.hgnc_id.table[i,'description'] = 'immunoglobulin lambda variable 1-44 [Source:HGNC Symbol;Acc:HGNC:5879]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '28825') {		# https://www.ncbi.nlm.nih.gov/gene/28825; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:5877; http://www.uniprot.org/uniprot/P01703
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'IGLV1-40'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:5877'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '28825'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'P01703'
			manually_curated.hgnc_id.table[i,'description'] = 'immunoglobulin lambda variable 1-40 [Source:HGNC Symbol;Acc:HGNC:5877]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '28834') {		# https://www.ncbi.nlm.nih.gov/gene/28834; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:5861; http://www.uniprot.org/uniprot/A0M8Q6
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'IGLC7'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:5861'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '28834'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'A0M8Q6'
			manually_curated.hgnc_id.table[i,'description'] = 'immunoglobulin lambda constant 7 [Source:HGNC Symbol;Acc:HGNC:5861]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '28875') {		# https://www.ncbi.nlm.nih.gov/gene/28875; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:5824; http://www.uniprot.org/uniprot/A0A087WSY6
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'IGKV3D-15'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:5824'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '28875'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'A0A087WSY6'
			manually_curated.hgnc_id.table[i,'description'] = 'immunoglobulin kappa variable 3D-15 (gene/pseudogene) [Source:HGNC Symbol;Acc:HGNC:5824]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '28908') {		# https://www.ncbi.nlm.nih.gov/gene/28908; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:5834; http://www.uniprot.org/uniprot/P06312
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'IGKV4-1'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:5834'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '28908'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'P06312'
			manually_curated.hgnc_id.table[i,'description'] = 'immunoglobulin kappa variable 4-1 [Source:HGNC Symbol;Acc:HGNC:5834]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '28912') {		# https://www.ncbi.nlm.nih.gov/gene/28912; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:5817; http://www.uniprot.org/uniprot/P01619
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'IGKV3-20'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:5817'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '28912'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'P01619'
			manually_curated.hgnc_id.table[i,'description'] = 'immunoglobulin kappa variable 3-20 [Source:HGNC Symbol;Acc:HGNC:5817]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '28923') {		# https://www.ncbi.nlm.nih.gov/gene/28923; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:5781; http://www.uniprot.org/uniprot/A0A0C4DH68
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'IGKV2-24'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:5781'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '28923'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'A0A0C4DH68'
			manually_curated.hgnc_id.table[i,'description'] = 'immunoglobulin kappa variable 2-24 [Source:HGNC Symbol;Acc:HGNC:5781]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '3514') {			# https://www.ncbi.nlm.nih.gov/gene/3514; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:5716; http://www.uniprot.org/uniprot/P01834 (http://www.uniprot.org/uniprot/P01601)
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'IGKC'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:5716'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '3514'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'P01834'
			manually_curated.hgnc_id.table[i,'description'] = 'immunoglobulin kappa constant [Source:HGNC Symbol;Acc:HGNC:5716]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '28396') {		# https://www.ncbi.nlm.nih.gov/gene/28396; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:5649; http://www.uniprot.org/uniprot/P0DP07 (http://www.uniprot.org/uniprot/A0A087WSY4)
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'IGHV4-31'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:5649'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '28396'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'P0DP07'
			manually_curated.hgnc_id.table[i,'description'] = 'immunoglobulin heavy variable 4-31 [Source:HGNC Symbol;Acc:HGNC:5649]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '64006') {		# https://www.ncbi.nlm.nih.gov/gene/64006; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:13915; http://www.uniprot.org/uniprot/Q69384 (http://www.uniprot.org/uniprot/Q69383)
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'ERVK-6'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:13915'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '64006'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'Q69384'
			manually_curated.hgnc_id.table[i,'description'] = 'endogenous retrovirus group K member 6, envelope [Source:HGNC Symbol;Acc:HGNC:13915]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '100130100') {	# https://www.ncbi.nlm.nih.gov/gene/100130100 -> discontinued on 9-Jun-2016: This record has been withdrawn by NCBI because the model on which it was based was not predicted in a later annotation.
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'LOC100130100'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = ''
			manually_curated.hgnc_id.table[i,'entrezgene'] = '100130100'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'Ig kappa chain V-I region Walker-like'
		# --- Gene type: pseudo -------------------------------------------------------------------
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '400750') {		# https://www.ncbi.nlm.nih.gov/gene/400750
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'LOC400750'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = ''
			manually_curated.hgnc_id.table[i,'entrezgene'] = '400750'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'heat shock protein family A (Hsp70) member 5 pseudogene'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '441241') {		# https://www.ncbi.nlm.nih.gov/gene/441241
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'LOC441241'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = ''
			manually_curated.hgnc_id.table[i,'entrezgene'] = '441241'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'vitamin K epoxide reductase complex subunit 1 like 1 pseudogene'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '646127') {		# https://www.ncbi.nlm.nih.gov/gene/646127
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'LOC646127'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = ''
			manually_curated.hgnc_id.table[i,'entrezgene'] = '646127'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'telomeric repeat-binding factor 1 pseudogene'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '728533') {		# https://www.ncbi.nlm.nih.gov/gene/728533
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'LOC728533'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = ''
			manually_curated.hgnc_id.table[i,'entrezgene'] = '728533'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'IST1 homolog'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '100129492') {	# https://www.ncbi.nlm.nih.gov/gene/100129492
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'LOC100129492'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = ''
			manually_curated.hgnc_id.table[i,'entrezgene'] = '100129492'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'small nuclear ribonucleoprotein Sm D1 pseudogene'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '100133211') {	# https://www.ncbi.nlm.nih.gov/gene/100133211
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'LOC100133211'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = ''
			manually_curated.hgnc_id.table[i,'entrezgene'] = '100133211'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'related RAS viral (r-ras) oncogene homolog 2 pseudogene'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '102724737') {	# https://www.ncbi.nlm.nih.gov/gene/102724737
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'LOC102724737'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = ''
			manually_curated.hgnc_id.table[i,'entrezgene'] = '102724737'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = '40S ribosomal protein S8 pseudogene'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '643752') {		# https://www.ncbi.nlm.nih.gov/gene/643752; http://www.uniprot.org/uniprot/A6NIZ1
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RAP1BL'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = ''
			manually_curated.hgnc_id.table[i,'entrezgene'] = '643752'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'A6NIZ1'
			manually_curated.hgnc_id.table[i,'description'] = 'RAP1B, member of RAS oncogene family pseudogene'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '245') {			# https://www.ncbi.nlm.nih.gov/gene/245; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:432
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'ALOX12P2'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:432'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '245'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'arachidonate 12-lipoxygenase pseudogene 2 [Source:HGNC Symbol;Acc:HGNC:432]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '303') {			# https://www.ncbi.nlm.nih.gov/gene/303; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:538
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'ANXA2P1'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:538'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '303'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'annexin A2 pseudogene 1 [Source:HGNC Symbol;Acc:HGNC:538]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '56604') {		# https://www.ncbi.nlm.nih.gov/gene/56604; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:12413
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'TUBB7P'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:12413'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '56604'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'tubulin beta 7 pseudogene [Source:HGNC Symbol;Acc:HGNC:12413]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '56969') {		# https://www.ncbi.nlm.nih.gov/gene/56969; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:25345
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RPL23AP32'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:25345'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '56969'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'ribosomal protein L23a pseudogene 32 [Source:HGNC Symbol;Acc:HGNC:25345]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '92755') {		# https://www.ncbi.nlm.nih.gov/gene/92755; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:12414
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'TUBBP1'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:12414'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '92755'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'tubulin beta pseudogene 1 [Source:HGNC Symbol;Acc:HGNC:12414]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '113157') {		# https://www.ncbi.nlm.nih.gov/gene/113157; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:17960
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RPLP0P2'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:17960'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '113157'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'ribosomal protein lateral stalk subunit P0 pseudogene 2 [Source:HGNC Symbol;Acc:HGNC:17960]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '122589') {		# https://www.ncbi.nlm.nih.gov/gene/122589; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:23558
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RPLP0P3'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:23558'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '122589'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'ribosomal protein lateral stalk subunit P0 pseudogene 3 [Source:HGNC Symbol;Acc:HGNC:23558]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '128192') {		# https://www.ncbi.nlm.nih.gov/gene/128192; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:53659
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'PPIAP35'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:53659'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '128192'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'peptidylprolyl isomerase A pseudogene 35 [Source:HGNC Symbol;Acc:HGNC:53659]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '149501') {		# https://www.ncbi.nlm.nih.gov/gene/149501; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:39879
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'KRT8P45'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:39879'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '149501'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'keratin 8 pseudogene 45 [Source:HGNC Symbol;Acc:HGNC:39879]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '220885') {		# https://www.ncbi.nlm.nih.gov/gene/220885; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:31464
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RPSAP15'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:31464'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '220885'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'ribosomal protein SA pseudogene 15 [Source:HGNC Symbol;Acc:HGNC:31464]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '221438') {		# https://www.ncbi.nlm.nih.gov/gene/221438; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:30808
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'TREML5P'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:30808'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '221438'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'triggering receptor expressed on myeloid cells like 5, pseudogene [Source:HGNC Symbol;Acc:HGNC:30808]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '221547') {		# https://www.ncbi.nlm.nih.gov/gene/221547; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:21631
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RANP1'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:21631'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '221547'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'RAN, member RAS oncogene family pseudogene 1 [Source:HGNC Symbol;Acc:HGNC:21631]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '283412') {		# https://www.ncbi.nlm.nih.gov/gene/283412; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:36608
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RPL29P12'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:36608'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '283412'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'ribosomal protein L29 pseudogene 12 [Source:HGNC Symbol;Acc:HGNC:36608]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '283523') {		# https://www.ncbi.nlm.nih.gov/gene/283523; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:39686
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'TERF1P5'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:39686'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '283523'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'telomeric repeat binding factor 1 pseudogene 5 [Source:HGNC Symbol;Acc:HGNC:39686]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '284764') {		# https://www.ncbi.nlm.nih.gov/gene/284764; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:4000
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'FTLP3'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:4000'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '284764'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'ferritin light chain pseudogene 3 [Source:HGNC Symbol;Acc:HGNC:4000]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '286444') {		# https://www.ncbi.nlm.nih.gov/gene/286444; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:36140
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RPS2P55'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:36140'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '286444'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'ribosomal protein S2 pseudogene 55 [Source:HGNC Symbol;Acc:HGNC:36140]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '339781') {		# https://www.ncbi.nlm.nih.gov/gene/339781; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:33387
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'KRT18P19'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:33387'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '339781'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'keratin 18 pseudogene 19 [Source:HGNC Symbol;Acc:HGNC:33387]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '343184') {		# https://www.ncbi.nlm.nih.gov/gene/343184; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:36350
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RPS2P11'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:36350'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '343184'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'ribosomal protein S2 pseudogene 11 [Source:HGNC Symbol;Acc:HGNC:36350]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '345041') {		# https://www.ncbi.nlm.nih.gov/gene/345041; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:5266
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'HSPD1P5'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:5266'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '345041'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'heat shock protein family D (Hsp60) member 1 pseudogene 5 [Source:HGNC Symbol;Acc:HGNC:5266]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '376693') {		# https://www.ncbi.nlm.nih.gov/gene/376693; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:36423
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RPS10P7'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:36423'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '376693'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'ribosomal protein S10 pseudogene 7 [Source:HGNC Symbol;Acc:HGNC:36423]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '387867') {		# https://www.ncbi.nlm.nih.gov/gene/387867; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:6504
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RPSAP12'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:6504'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '387867'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'ribosomal protein SA pseudogene 12 [Source:HGNC Symbol;Acc:HGNC:6504]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '388339') {		# https://www.ncbi.nlm.nih.gov/gene/388339; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:35967
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RPS18P12'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:35967'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '388339'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'ribosomal protein S18 pseudogene 12 [Source:HGNC Symbol;Acc:HGNC:35967]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '388707') {		# https://www.ncbi.nlm.nih.gov/gene/388707; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:36106
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RPSAP18'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:36106'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '388707'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'ribosomal protein SA pseudogene 18 [Source:HGNC Symbol;Acc:HGNC:36106]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '389141') {		# https://www.ncbi.nlm.nih.gov/gene/389141; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:36680
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RPSAP29'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:36680'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '389141'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'ribosomal protein SA pseudogene 29 [Source:HGNC Symbol;Acc:HGNC:36680]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '390183') {		# https://www.ncbi.nlm.nih.gov/gene/390183; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:36522
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RPS4XP13'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:36522'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '390183'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'ribosomal protein S4X pseudogene 13 [Source:HGNC Symbol;Acc:HGNC:36522]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '390601') {		# https://www.ncbi.nlm.nih.gov/gene/390601; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:33363
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'KRT8P9'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:33363'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '390601'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'keratin 8 pseudogene 9 [Source:HGNC Symbol;Acc:HGNC:33363]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '391777') {		# https://www.ncbi.nlm.nih.gov/gene/391777; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:36057
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RPS4XP6'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:36057'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '391777'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'ribosomal protein S4X pseudogene 6 [Source:HGNC Symbol;Acc:HGNC:36057]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '391833') {		# https://www.ncbi.nlm.nih.gov/gene/391833; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:36775
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RPS10P11'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:36775'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '391833'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'ribosomal protein S10 pseudogene 11 [Source:HGNC Symbol;Acc:HGNC:36775]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '399942') {		# https://www.ncbi.nlm.nih.gov/gene/399942; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:14531
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'TUBAP2'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:14531'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '399942'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'tubulin alpha pseudogene 2 [Source:HGNC Symbol;Acc:HGNC:14531]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '400389') {		# https://www.ncbi.nlm.nih.gov/gene/400389; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:35969
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RPL12P35'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:35969'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '400389'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'ribosomal protein L12 pseudogene 35 [Source:HGNC Symbol;Acc:HGNC:35969]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '400578') {		# https://www.ncbi.nlm.nih.gov/gene/400578; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:37807
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'KRT16P2'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:37807'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '400578'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'keratin 16 pseudogene 2 [Source:HGNC Symbol;Acc:HGNC:37807]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '400963') {		# https://www.ncbi.nlm.nih.gov/gene/400963; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:35764
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RPS2P17'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:35764'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '400963'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'ribosomal protein S2 pseudogene 17 [Source:HGNC Symbol;Acc:HGNC:35764]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '401331') {		# https://www.ncbi.nlm.nih.gov/gene/401331; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:44185
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RASA4CP'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:44185'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '401331'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'RAS p21 protein activator 4C, pseudogene [Source:HGNC Symbol;Acc:HGNC:44185]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '401817') {		# https://www.ncbi.nlm.nih.gov/gene/401817; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:35814
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RPS10P22'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:35814'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '401817'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'ribosomal protein S10 pseudogene 22 [Source:HGNC Symbol;Acc:HGNC:35814]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '439953') {		# https://www.ncbi.nlm.nih.gov/gene/439953; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:44962
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'PPIAP31'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:44962'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '439953'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'peptidylprolyl isomerase A pseudogene 31 [Source:HGNC Symbol;Acc:HGNC:44962]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '440176') {		# https://www.ncbi.nlm.nih.gov/gene/440176; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:23538
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RPL12P6'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:23538'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '440176'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'ribosomal protein L12 pseudogene 6 [Source:HGNC Symbol;Acc:HGNC:23538]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '440589') {		# https://www.ncbi.nlm.nih.gov/gene/440589; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:36207
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RPS2P8'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:36207'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '440589'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'ribosomal protein S2 pseudogene 8 [Source:HGNC Symbol;Acc:HGNC:36207]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '440917') {		# https://www.ncbi.nlm.nih.gov/gene/440917; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:49439
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'YWHAEP5'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:49439'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '440917'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'tyrosine 3-monooxygenase/tryptophan 5-monooxygenase activation protein epsilon pseudogene 5 [Source:HGNC Symbol;Acc:HGNC:49439]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '440991') {		# https://www.ncbi.nlm.nih.gov/gene/440991; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:36700
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RPS3P3'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:36700'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '440991'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'ribosomal protein S3 pseudogene 3 [Source:HGNC Symbol;Acc:HGNC:36700]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '441876') {		# https://www.ncbi.nlm.nih.gov/gene/441876; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:36462
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RPS16P1'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:36462'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '441876'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'ribosomal protein S16 pseudogene 1 [Source:HGNC Symbol;Acc:HGNC:36462]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '442308') {		# https://www.ncbi.nlm.nih.gov/gene/442308; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:42175
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'TUBBP6'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:42175'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '442308'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'tubulin beta class I pseudogene 6 [Source:HGNC Symbol;Acc:HGNC:42175]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '541611') {		# https://www.ncbi.nlm.nih.gov/gene/541611; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:32540
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'HSP90AB6P'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:32540'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '541611'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'heat shock protein 90 alpha family class B member 6, pseudogene [Source:HGNC Symbol;Acc:HGNC:32540]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '643300') {		# https://www.ncbi.nlm.nih.gov/gene/643300; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:35133
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'HSPD1P1'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:35133'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '643300'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'heat shock protein family D (Hsp60) member 1 pseudogene 1 [Source:HGNC Symbol;Acc:HGNC:35133]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '643358') {		# https://www.ncbi.nlm.nih.gov/gene/643358; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:36855
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RPS27AP16'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:36855'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '643358'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'ribosomal protein S27a pseudogene 16 [Source:HGNC Symbol;Acc:HGNC:36855]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '643531') {		# https://www.ncbi.nlm.nih.gov/gene/643531; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:35717
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RPL29P26'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:35717'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '643531'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'ribosomal protein L29 pseudogene 26 [Source:HGNC Symbol;Acc:HGNC:35717]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '643617') {		# https://www.ncbi.nlm.nih.gov/gene/643617; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:31462
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RPSAP8'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:31462'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '643617'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'ribosomal protein SA pseudogene 8 [Source:HGNC Symbol;Acc:HGNC:31462]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '643751') {		# https://www.ncbi.nlm.nih.gov/gene/643751; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:44430
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'CDC42P6'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:44430'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '643751'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'cell division cycle 42 pseudogene 6 [Source:HGNC Symbol;Acc:HGNC:44430]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '644464') {		# https://www.ncbi.nlm.nih.gov/gene/644464; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:36760
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RPSAP61'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:36760'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '644464'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'ribosomal protein SA pseudogene 61 [Source:HGNC Symbol;Acc:HGNC:36760]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '644745') {		# https://www.ncbi.nlm.nih.gov/gene/644745; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:35146
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'HSPD1P4'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:35146'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '644745'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'heat shock protein family D (Hsp60) member 1 pseudogene 4 [Source:HGNC Symbol;Acc:HGNC:35146]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '645018') {		# https://www.ncbi.nlm.nih.gov/gene/645018; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:35992
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RPS2P20'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:35992'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '645018'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'ribosomal protein S2 pseudogene 20 [Source:HGNC Symbol;Acc:HGNC:35992]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '645548') {		# https://www.ncbi.nlm.nih.gov/gene/645548; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:5267
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'HSPD1P6'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:5267'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '645548'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'heat shock protein family D (Hsp60) member 1 pseudogene 6 [Source:HGNC Symbol;Acc:HGNC:5267]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '646282') {		# https://www.ncbi.nlm.nih.gov/gene/646282; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:911
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'AZGP1P1'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:911'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '646282'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'alpha-2-glycoprotein 1, zinc-binding pseudogene 1 [Source:HGNC Symbol;Acc:HGNC:911]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '646316') {		# https://www.ncbi.nlm.nih.gov/gene/646316; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:38499
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'TERF1P3'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:38499'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '646316'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'telomeric repeat binding factor 1 pseudogene 3 [Source:HGNC Symbol;Acc:HGNC:38499]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '646359') {		# https://www.ncbi.nlm.nih.gov/gene/646359; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:38110
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'TERF1P2'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:38110'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '646359'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'telomeric repeat binding factor 1 pseudogene 2 [Source:HGNC Symbol;Acc:HGNC:38110]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '646785') {		# https://www.ncbi.nlm.nih.gov/gene/646785; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:36374
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RPS10P13'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:36374'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '646785'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'ribosomal protein S10 pseudogene 13 [Source:HGNC Symbol;Acc:HGNC:36374]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '646875') {		# https://www.ncbi.nlm.nih.gov/gene/646875; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:16070
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RPL12P2'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:16070'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '646875'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'ribosomal protein L12 pseudogene 2 [Source:HGNC Symbol;Acc:HGNC:16070]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '646949') {		# https://www.ncbi.nlm.nih.gov/gene/646949; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:35554
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RPL23P6'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:35554'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '646949'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'ribosomal protein L23 pseudogene 6 [Source:HGNC Symbol;Acc:HGNC:35554]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '647000') {		# https://www.ncbi.nlm.nih.gov/gene/647000; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:12415
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'TUBBP2'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:12415'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '647000'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'tubulin beta pseudogene 2 [Source:HGNC Symbol;Acc:HGNC:12415]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '647099') {		# https://www.ncbi.nlm.nih.gov/gene/647099; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:35942
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RPL23AP42'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:35942'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '647099'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'ribosomal protein L23a pseudogene 42 [Source:HGNC Symbol;Acc:HGNC:35942]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '647285') {		# https://www.ncbi.nlm.nih.gov/gene/647285; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:35729
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RPL29P9'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:35729'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '647285'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'ribosomal protein L29 pseudogene 9 [Source:HGNC Symbol;Acc:HGNC:35729]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '650901') {		# https://www.ncbi.nlm.nih.gov/gene/650901; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:36546
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RPS2P12'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:36546'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '650901'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'ribosomal protein S2 pseudogene 12 [Source:HGNC Symbol;Acc:HGNC:36546]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '653162') {		# https://www.ncbi.nlm.nih.gov/gene/653162; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:31463
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RPSAP9'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:31463'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '653162'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'ribosomal protein SA pseudogene 9 [Source:HGNC Symbol;Acc:HGNC:31463]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '653214') {		# https://www.ncbi.nlm.nih.gov/gene/653214; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:17236
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'PPIAP22'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:17236'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '653214'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'peptidylprolyl isomerase A pseudogene 22 [Source:HGNC Symbol;Acc:HGNC:17236]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '653232') {		# https://www.ncbi.nlm.nih.gov/gene/653232; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:21538
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RPL15P3'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:21538'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '653232'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'ribosomal protein L15 pseudogene 3 [Source:HGNC Symbol;Acc:HGNC:21538]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '653553') {		# https://www.ncbi.nlm.nih.gov/gene/653553; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:5251
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'HSPB1P1'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:5251'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '653553'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'heat shock protein family B (small) member 1 pseudogene 1 [Source:HGNC Symbol;Acc:HGNC:5251]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '728002') {		# https://www.ncbi.nlm.nih.gov/gene/728002; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:36852
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RPL15P17'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:36852'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '728002'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'ribosomal protein L15 pseudogene 17 [Source:HGNC Symbol;Acc:HGNC:36852]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '728088') {		# https://www.ncbi.nlm.nih.gov/gene/728088; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:36515
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RPL15P18'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:36515'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '728088'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'ribosomal protein L15 pseudogene 18 [Source:HGNC Symbol;Acc:HGNC:36515]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '728576') {		# https://www.ncbi.nlm.nih.gov/gene/728576; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:36808
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RPL15P7'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:36808'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '728576'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'ribosomal protein L15 pseudogene 7 [Source:HGNC Symbol;Acc:HGNC:36808]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '728590') {		# https://www.ncbi.nlm.nih.gov/gene/728590; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:36126
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RPS27AP11'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:36126'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '728590'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'ribosomal protein S27a pseudogene 11 [Source:HGNC Symbol;Acc:HGNC:36126]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '728641') {		# https://www.ncbi.nlm.nih.gov/gene/728641; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:31070
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'FABP5P7'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:31070'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '728641'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'fatty acid binding protein 5 pseudogene 7 [Source:HGNC Symbol;Acc:HGNC:31070]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '728791') {		# https://www.ncbi.nlm.nih.gov/gene/728791; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:31364
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RPS10P4'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:31364'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '728791'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'ribosomal protein S10 pseudogene 4 [Source:HGNC Symbol;Acc:HGNC:31364]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '728979') {		# https://www.ncbi.nlm.nih.gov/gene/728979; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:35947
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RPL10AP9'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:35947'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '728979'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'ribosomal protein L10a pseudogene 9 [Source:HGNC Symbol;Acc:HGNC:35947]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '729500') {		# https://www.ncbi.nlm.nih.gov/gene/729500; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:36495
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RPL12P14'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:36495'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '729500'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'ribosomal protein L12 pseudogene 14 [Source:HGNC Symbol;Acc:HGNC:36495]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '729659') {		# https://www.ncbi.nlm.nih.gov/gene/729659; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:10491
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'S100A11P1'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:10491'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '729659'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'S100 calcium binding protein A11 pseudogene 1 [Source:HGNC Symbol;Acc:HGNC:10491]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '729679') {		# https://www.ncbi.nlm.nih.gov/gene/729679; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:36422
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RPS2P51'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:36422'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '729679'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'ribosomal protein S2 pseudogene 51 [Source:HGNC Symbol;Acc:HGNC:36422]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '729682') {		# https://www.ncbi.nlm.nih.gov/gene/729682; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:33697
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'KRT17P3'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:33697'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '729682'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'keratin 17 pseudogene 3 [Source:HGNC Symbol;Acc:HGNC:33697]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '729708') {		# https://www.ncbi.nlm.nih.gov/gene/729708; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:35449
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'TPI1P1'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:35449'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '729708'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'triosephosphate isomerase 1 pseudogene 1 [Source:HGNC Symbol;Acc:HGNC:35449]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '729903') {		# https://www.ncbi.nlm.nih.gov/gene/729903; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:35739
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RPS16P10'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:35739'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '729903'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'ribosomal protein S16 pseudogene 10 [Source:HGNC Symbol;Acc:HGNC:35739]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '730029') {		# https://www.ncbi.nlm.nih.gov/gene/730029; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:36508
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RPSAP19'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:36508'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '730029'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'ribosomal protein SA pseudogene 19 [Source:HGNC Symbol;Acc:HGNC:36508]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '100128936') {	# https://www.ncbi.nlm.nih.gov/gene/100128936; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:36504
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RPL10AP6'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:36504'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '100128936'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'ribosomal protein L10a pseudogene 6 [Source:HGNC Symbol;Acc:HGNC:36504]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '100129982') {	# https://www.ncbi.nlm.nih.gov/gene/100129982; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:36248
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RPL12P19'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:36248'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '100129982'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'ribosomal protein L12 pseudogene 19 [Source:HGNC Symbol;Acc:HGNC:36248]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '100130446') {	# https://www.ncbi.nlm.nih.gov/gene/100130446; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:35732
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RPS27AP12'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:35732'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '100130446'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'ribosomal protein S27a pseudogene 12 [Source:HGNC Symbol;Acc:HGNC:35732]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '100130562') {	# https://www.ncbi.nlm.nih.gov/gene/100130562; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:31386
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RPS2P5'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:31386'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '100130562'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'ribosomal protein S2 pseudogene 5 [Source:HGNC Symbol;Acc:HGNC:31386]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '100130624') {	# https://www.ncbi.nlm.nih.gov/gene/100130624; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:36933
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RPL15P22'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:36933'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '100130624'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'ribosomal protein L15 pseudogene 22 [Source:HGNC Symbol;Acc:HGNC:36933]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '100131713') {	# https://www.ncbi.nlm.nih.gov/gene/100131713; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:36905
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RPL29P11'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:36905'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '100131713'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'ribosomal protein L29 pseudogene 11 [Source:HGNC Symbol;Acc:HGNC:36905]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '100131863') {	# https://www.ncbi.nlm.nih.gov/gene/100131863; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:36287
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RPS18P5'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:36287'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '100131863'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'ribosomal protein S18 pseudogene 5 [Source:HGNC Symbol;Acc:HGNC:36287]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '100132795') {	# https://www.ncbi.nlm.nih.gov/gene/100132795; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:36597
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RPL12P32'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:36597'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '100132795'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'ribosomal protein L12 pseudogene 32 [Source:HGNC Symbol;Acc:HGNC:36597]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '304') {			# https://www.ncbi.nlm.nih.gov/gene/304; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:539; http://www.uniprot.org/uniprot/A6NMY6
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'ANXA2P2'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:539'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '304'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'A6NMY6'
			manually_curated.hgnc_id.table[i,'description'] = 'annexin A2 pseudogene 2 [Source:HGNC Symbol;Acc:HGNC:539]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '2679') {			# https://www.ncbi.nlm.nih.gov/gene/2679; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:4252; http://www.uniprot.org/uniprot/A6NGU5
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'GGT3P'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:4252'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '2679'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'A6NGU5'
			manually_curated.hgnc_id.table[i,'description'] = 'gamma-glutamyltransferase 3 pseudogene [Source:HGNC Symbol;Acc:HGNC:4252]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '3136') {			# https://www.ncbi.nlm.nih.gov/gene/3136; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:4965; http://www.uniprot.org/uniprot/P01893
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'HLA-H'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:4965'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '3136'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'P01893'
			manually_curated.hgnc_id.table[i,'description'] = 'major histocompatibility complex, class I, H (pseudogene) [Source:HGNC Symbol;Acc:HGNC:4965]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '3311') {			# https://www.ncbi.nlm.nih.gov/gene/3311; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:5240; http://www.uniprot.org/uniprot/P48741
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'HSPA7'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:5240'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '3311'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'P48741'
			manually_curated.hgnc_id.table[i,'description'] = 'heat shock protein family A (Hsp70) member 7 [Source:HGNC Symbol;Acc:HGNC:5240]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '3323') {			# https://www.ncbi.nlm.nih.gov/gene/3323; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:5255; http://www.uniprot.org/uniprot/Q58FG1
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'HSP90AA4P'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:5255'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '3323'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'Q58FG1'
			manually_curated.hgnc_id.table[i,'description'] = 'heat shock protein 90 alpha family class A member 4, pseudogene [Source:HGNC Symbol;Acc:HGNC:5255]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '3324') {			# https://www.ncbi.nlm.nih.gov/gene/3324; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:5256; http://www.uniprot.org/uniprot/Q14568
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'HSP90AA2P'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:5256'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '3324'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'Q14568'
			manually_curated.hgnc_id.table[i,'description'] = 'heat shock protein 90 alpha family class A member 2, pseudogene [Source:HGNC Symbol;Acc:HGNC:5256]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '3327') {			# https://www.ncbi.nlm.nih.gov/gene/3327; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:5259; http://www.uniprot.org/uniprot/Q58FF7
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'HSP90AB3P'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:5259'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '3327'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'Q58FF7'
			manually_curated.hgnc_id.table[i,'description'] = 'heat shock protein 90 alpha family class B member 3, pseudogene [Source:HGNC Symbol;Acc:HGNC:5259]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '3542') {			# https://www.ncbi.nlm.nih.gov/gene/3542; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:5860; http://www.uniprot.org/uniprot/P0CF74
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'IGLC6'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:5860'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '3542'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'P0CF74'
			manually_curated.hgnc_id.table[i,'description'] = 'immunoglobulin lambda constant 6 (gene/pseudogene) [Source:HGNC Symbol;Acc:HGNC:5860]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '144106') {		# https://www.ncbi.nlm.nih.gov/gene/144106; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:18556; http://www.uniprot.org/uniprot/Q8NFI4
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'ST13P5'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:18556'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '144106'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'Q8NFI4'
			manually_curated.hgnc_id.table[i,'description'] = 'ST13, Hsp70 interacting protein pseudogene 5 [Source:HGNC Symbol;Acc:HGNC:18556]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '145165') {		# https://www.ncbi.nlm.nih.gov/gene/145165; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:18487; http://www.uniprot.org/uniprot/Q8IZP2
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'ST13P4'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:18487'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '145165'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'Q8IZP2'
			manually_curated.hgnc_id.table[i,'description'] = 'ST13, Hsp70 interacting protein pseudogene 4 [Source:HGNC Symbol;Acc:HGNC:18487]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '158078') {		# https://www.ncbi.nlm.nih.gov/gene/158078; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:3200; http://www.uniprot.org/uniprot/Q5VTE0
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'EEF1A1P5'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:3200'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '158078'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'Q5VTE0'
			manually_curated.hgnc_id.table[i,'description'] = 'eukaryotic translation elongation factor 1 alpha 1 pseudogene 5 [Source:HGNC Symbol;Acc:HGNC:3200]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '220717') {		# https://www.ncbi.nlm.nih.gov/gene/220717; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:36404; http://www.uniprot.org/uniprot/Q8NHW5
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RPLP0P6'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:36404'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '220717'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'Q8NHW5'
			manually_curated.hgnc_id.table[i,'description'] = 'ribosomal protein lateral stalk subunit P0 pseudogene 6 [Source:HGNC Symbol;Acc:HGNC:36404]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '283458') {		# https://www.ncbi.nlm.nih.gov/gene/283458; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:31358; http://www.uniprot.org/uniprot/O60361
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'NME2P1'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:31358'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '283458'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'O60361'
			manually_curated.hgnc_id.table[i,'description'] = 'NME/NM23 nucleoside diphosphate kinase 2 pseudogene 1 [Source:HGNC Symbol;Acc:HGNC:31358]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '286310') {		# https://www.ncbi.nlm.nih.gov/gene/286310; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:23412; http://www.uniprot.org/uniprot/Q5VSP4
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'LCN1P1'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:23412'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '286310'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'Q5VSP4'
			manually_curated.hgnc_id.table[i,'description'] = 'lipocalin 1 pseudogene 1 [Source:HGNC Symbol;Acc:HGNC:23412]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '337873') {		# https://www.ncbi.nlm.nih.gov/gene/337873; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:20516; http://www.uniprot.org/uniprot/Q6DN03
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'HIST2H2BC'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:20516'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '337873'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'Q6DN03'
			manually_curated.hgnc_id.table[i,'description'] = 'histone cluster 2 H2B family member c (pseudogene) [Source:HGNC Symbol;Acc:HGNC:20516]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '391634') {		# https://www.ncbi.nlm.nih.gov/gene/391634; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:32537; http://www.uniprot.org/uniprot/Q58FF8
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'HSP90AB2P'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:32537'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '391634'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'Q58FF8'
			manually_curated.hgnc_id.table[i,'description'] = 'heat shock protein 90 alpha family class B member 2, pseudogene [Source:HGNC Symbol;Acc:HGNC:32537]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '440915') {		# https://www.ncbi.nlm.nih.gov/gene/440915; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:30182; http://www.uniprot.org/uniprot/Q9BYX7
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'POTEKP'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:30182'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '440915'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'Q9BYX7'
			manually_curated.hgnc_id.table[i,'description'] = 'POTE ankyrin domain family member K, pseudogene [Source:HGNC Symbol;Acc:HGNC:30182]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '441400') {		# https://www.ncbi.nlm.nih.gov/gene/441400; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:23683; http://www.uniprot.org/uniprot/Q92928
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'RAB1C'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:23683'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '441400'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'Q92928'
			manually_curated.hgnc_id.table[i,'description'] = 'RAB1C, member RAS oncogene family pseudogene [Source:HGNC Symbol;Acc:HGNC:23683]'
		# --- Gene type: ncRNA --------------------------------------------------------------------
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '100996740') {	# https://www.ncbi.nlm.nih.gov/gene/100996740
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'LOC100996740'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = ''
			manually_curated.hgnc_id.table[i,'entrezgene'] = '100996740'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'uncharacterized LOC100996740'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '101060400') {	# https://www.ncbi.nlm.nih.gov/gene/101060400
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'LOC101060400'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = ''
			manually_curated.hgnc_id.table[i,'entrezgene'] = '101060400'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'uncharacterized LOC101060400'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '284889') {		# https://www.ncbi.nlm.nih.gov/gene/284889; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:27669
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'MIF-AS1'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:27669'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '284889'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'MIF antisense RNA 1 [Source:HGNC Symbol;Acc:HGNC:27669]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '677779') {		# https://www.ncbi.nlm.nih.gov/gene/677779; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:32675
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'LINC00488'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:32675'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '677779'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'long intergenic non-protein coding RNA 488 [Source:HGNC Symbol;Acc:HGNC:32675]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '79667') {		# https://www.ncbi.nlm.nih.gov/gene/79667; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:25796; http://www.uniprot.org/uniprot/Q9H8V8
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'KLF3-AS1'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:25796'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '79667'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'Q9H8V8'
			manually_curated.hgnc_id.table[i,'description'] = 'KLF3 antisense RNA 1 [Source:HGNC Symbol;Acc:HGNC:25796]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '493913') {		# https://www.ncbi.nlm.nih.gov/gene/493913; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:35152; http://www.uniprot.org/uniprot/Q5QFB9
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'PAPPA-AS1'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:35152'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '493913'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = 'Q5QFB9'
			manually_curated.hgnc_id.table[i,'description'] = 'PAPPA antisense RNA 1 [Source:HGNC Symbol;Acc:HGNC:35152]'
		# --- Gene type: snoRNA -------------------------------------------------------------------
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '619499') {		# https://www.ncbi.nlm.nih.gov/gene/619499; https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:32617
			manually_curated.hgnc_id.table[i,'hgnc_symbol'] = 'SNORA27'
			manually_curated.hgnc_id.table[i,'hgnc_id'] = 'HGNC:32617'
			manually_curated.hgnc_id.table[i,'entrezgene'] = '619499'
			manually_curated.hgnc_id.table[i,'uniprotswissprot'] = ''
			manually_curated.hgnc_id.table[i,'description'] = 'small nucleolar RNA, H/ACA box 27 [Source:HGNC Symbol;Acc:HGNC:32617]'
		} else {
			cat(paste('entrezgene "',ExoCarta_Protein.entrezgene.NA_attributes[i],'" did not match queried attributes (\'hgnc_symbol\',\'hgnc_id\',\'uniprotswissprot\',\'description\') by using R package "biomaRt".',sep=""),"\n")
			
			ExoCarta_Protein.entrezgene.NA_attributes.uncurated <- c(ExoCarta_Protein.entrezgene.NA_attributes.uncurated, ExoCarta_Protein.entrezgene.NA_attributes[i])
		}
		
		manually_curated.hgnc_id.table[i,'description'] = trimws(manually_curated.hgnc_id.table[i,'description'])
		manually_curated.hgnc_id.table[i,'description'] = sapply(strsplit(manually_curated.hgnc_id.table[i,'description'],'\\|'), FUN=function(U) paste(unique(sapply(strsplit(unlist(U),'\\[Source:'), FUN=function(V) trimws(unlist(V))[1])),collapse="|"))
	}
	
	if(length(ExoCarta_Protein.entrezgene.NA_attributes.uncurated) > 0)
	{
		write.table(ExoCarta_Protein.entrezgene.NA_attributes.uncurated, paste('./Data/',Sample,'/10_datasetOverlap/Database/hsa/ExoCarta_All.Protein_Gene2HGNC_NA.txt',sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE)
	}
	
	ExoCarta_Protein.hgnc_id.table.New = rbind(ExoCarta_Protein.hgnc_id.table, manually_curated.hgnc_id.table)
}

rownames(ExoCarta_Protein.hgnc_id.table.New) = c(1:nrow(ExoCarta_Protein.hgnc_id.table.New))
#ExoCarta_Protein.hgnc_id.table.New

human.UniProt.Entry2ProteinNames = sapply(strsplit(ExoCarta_Protein.hgnc_id.table.New[,'uniprotswissprot'],";"), FUN=function(X) paste(unique(sapply(unlist(X), FUN=function(Y) ifelse(length(which(human.UniProt.reviewed.Entry %in% unlist(Y)))==0, ifelse(length(which(human.UniProt.unreviewed.Entry %in% unlist(Y)))==0, '', human.UniProt.unreviewed.ProteinNames[which(human.UniProt.unreviewed.Entry %in% unlist(Y))]), human.UniProt.reviewed.ProteinNames[which(human.UniProt.reviewed.Entry %in% unlist(Y))]))),collapse="|"))
ExoCarta_Protein.hgnc_id.table.New = cbind(ExoCarta_Protein.hgnc_id.table.New, human.UniProt.Entry2ProteinNames)
colnames(ExoCarta_Protein.hgnc_id.table.New)[ncol(ExoCarta_Protein.hgnc_id.table.New)] = 'Protein_names'

ExoCarta_Protein.hgnc_id.table.New = ExoCarta_Protein.hgnc_id.table.New[which(rowSums(ExoCarta_Protein.hgnc_id.table.New=='')!=ncol(ExoCarta_Protein.hgnc_id.table.New)),,drop=FALSE]
rownames(ExoCarta_Protein.hgnc_id.table.New) = c(1:nrow(ExoCarta_Protein.hgnc_id.table.New))

#write.csv(ExoCarta_Protein.hgnc_id.table.New, paste('./Data/',Sample,'/10_datasetOverlap/Database/hsa/ExoCarta_All.Protein_Gene2HGNC.table.csv',sep=""), quote=TRUE, row.names=TRUE)

# --- write xlsx file by openxlsx package ---------------------------------------------------------
# Importing a big xlsx file into R? https://stackoverflow.com/questions/19147884/importing-a-big-xlsx-file-into-r/43118530#43118530
# Building R for Windows: https://cran.r-project.org/bin/windows/Rtools/
#Sys.setenv("R_ZIPCMD" = "C:/Rtools/bin/zip.exe")	# path to zip.exe

wb <- openxlsx::createWorkbook(paste('./Data/',Sample,'/10_datasetOverlap/Database/hsa/ExoCarta_All.Protein_Gene2HGNC.table.xlsx',sep=""))
modifyBaseFont(wb, fontSize=12, fontColour="black", fontName="Times New Roman")

sheet = 'ExoCarta_All'
sheetData = ExoCarta_Protein.hgnc_id.table.New

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
	}
}

if(length(wb$sheet_names) > 0) {
	openxlsx::saveWorkbook(wb, file=paste('./Data/',Sample,'/10_datasetOverlap/Database/hsa/ExoCarta_All.Protein_Gene2HGNC.table.xlsx',sep=""), overwrite=TRUE)
}


# /// mouse format ///
mainDir <- paste('./Data/',Sample,'/10_datasetOverlap/Database',sep="")
subDir <- "mmu"
if (file.exists(file.path(mainDir, subDir))){
} else {
	dir.create(file.path(mainDir, subDir))
}

ensembl = ensembl.mmu
#listFilters(ensembl)[,c('name','description')]
#listAttributes(ensembl)[,c('name','description')]

# --- search for IDs by entrezgene ---
a = ExoCarta_Protein.mmu.ENTREZ_GENE_ID.unique
b = ensembl.mmu
ExoCarta_Protein.entrezgene.attributes <- getBM(attributes=c('mgi_symbol','mgi_id','uniprotswissprot'), filters='entrezgene', values=a, mart=b)
ExoCarta_Protein.entrezgene.wikigene_description <- getBM(attributes=c('mgi_symbol','wikigene_description'), filters='entrezgene', values=a, mart=b)
ExoCarta_Protein.entrezgene.wikigene_description.table = merge(ExoCarta_Protein.entrezgene.wikigene_description, ExoCarta_Protein.entrezgene.attributes, by.x="mgi_symbol", by.y="mgi_symbol")
ExoCarta_Protein.entrezgene.attributes = ExoCarta_Protein.entrezgene.wikigene_description.table[,c('mgi_symbol','mgi_id','uniprotswissprot','wikigene_description'),drop=FALSE]
colnames(ExoCarta_Protein.entrezgene.attributes)[which(colnames(ExoCarta_Protein.entrezgene.attributes)=='wikigene_description')] = 'description'
ExoCarta_Protein.entrezgene.attributes = unique(ExoCarta_Protein.entrezgene.attributes)
ExoCarta_Protein.mgi_id <- getBM(attributes=c('entrezgene','mgi_id'), filters='entrezgene', values=a, mart=b)
ExoCarta_Protein.mgi_id.table = merge(ExoCarta_Protein.mgi_id, ExoCarta_Protein.entrezgene.attributes, by.x="mgi_id", by.y="mgi_id")
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='76219'),'uniprotswissprot'] = 'C0HK79'		# https://www.ncbi.nlm.nih.gov/gene/76219	http://www.uniprot.org/uniprot/C0HK79
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='76976'),'uniprotswissprot'] = 'C0HK80'		# https://www.ncbi.nlm.nih.gov/gene/76976	http://www.uniprot.org/uniprot/C0HK80

# https://www.ncbi.nlm.nih.gov/gene/11674; http://www.informatics.jax.org/marker/MGI:87994; http://www.uniprot.org/uniprot/P05064
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='11674'),'mgi_symbol'] = 'Aldoa'
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='11674'),'mgi_id'] = 'MGI:87994'
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='11674'),'uniprotswissprot'] = 'P05064'
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='11674'),'description'] = 'aldolase A, fructose-bisphosphate [Source:MGI Symbol;Acc:MGI:87994]'
# https://www.ncbi.nlm.nih.gov/gene/11946; http://www.informatics.jax.org/marker/MGI:88115; http://www.uniprot.org/uniprot/Q03265
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='11946'),'mgi_symbol'] = 'Atp5a1'
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='11946'),'mgi_id'] = 'MGI:88115'
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='11946'),'uniprotswissprot'] = 'Q03265'
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='11946'),'description'] = 'ATP synthase, H+ transporting, mitochondrial F1 complex, alpha subunit 1 [Source:MGI Symbol;Acc:MGI:88115]'
# https://www.ncbi.nlm.nih.gov/gene/14433; http://www.informatics.jax.org/marker/MGI:95640; http://www.uniprot.org/uniprot/P16858
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='14433'),'mgi_symbol'] = 'Gapdh'
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='14433'),'mgi_id'] = 'MGI:95640'
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='14433'),'uniprotswissprot'] = 'P16858'
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='14433'),'description'] = 'glyceraldehyde-3-phosphate dehydrogenase [Source:MGI Symbol;Acc:MGI:95640]'
# https://www.ncbi.nlm.nih.gov/gene/14980; http://www.informatics.jax.org/marker/MGI:95912; http://www.uniprot.org/uniprot/P01897
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='14980'),'mgi_symbol'] = 'H2-L'
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='14980'),'mgi_id'] = 'MGI:95912'
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='14980'),'uniprotswissprot'] = 'P01897'
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='14980'),'description'] = 'histocompatibility 2, D region locus L [Source:MGI Symbol;Acc:MGI:95912]'
# https://www.ncbi.nlm.nih.gov/gene/15129; http://www.informatics.jax.org/marker/MGI:96021; http://www.uniprot.org/uniprot/P02088
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='15129'),'mgi_symbol'] = 'Hbb-b1'
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='15129'),'mgi_id'] = 'MGI:96021'
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='15129'),'uniprotswissprot'] = 'P02088'
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='15129'),'description'] = 'hemoglobin, beta adult major chain [Source:MGI Symbol;Acc:MGI:96021]'
# https://www.ncbi.nlm.nih.gov/gene/15130; http://www.informatics.jax.org/marker/MGI:96022; http://www.uniprot.org/uniprot/P02089
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='15130'),'mgi_symbol'] = 'Hbb-b2'
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='15130'),'mgi_id'] = 'MGI:96022'
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='15130'),'uniprotswissprot'] = 'P02089'
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='15130'),'description'] = 'hemoglobin, beta adult minor chain [Source:MGI Symbol;Acc:MGI:96022]'
# https://www.ncbi.nlm.nih.gov/gene/15289; http://www.informatics.jax.org/marker/MGI:96113; http://www.uniprot.org/uniprot/P63158
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='15289'),'mgi_symbol'] = 'Hmgb1'
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='15289'),'mgi_id'] = 'MGI:96113'
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='15289'),'uniprotswissprot'] = 'P63158'
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='15289'),'description'] = 'high mobility group box 1 [Source:MGI Symbol;Acc:MGI:96113]'
# https://www.ncbi.nlm.nih.gov/gene/19324; http://www.informatics.jax.org/marker/MGI:97842; http://www.uniprot.org/uniprot/P62821
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='19324'),'mgi_symbol'] = 'Rab1a'
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='19324'),'mgi_id'] = 'MGI:97842'
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='19324'),'uniprotswissprot'] = 'P62821'
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='19324'),'description'] = 'RAB1A, member RAS oncogene family [Source:MGI Symbol;Acc:MGI:97842]'
# https://www.ncbi.nlm.nih.gov/gene/19942; http://www.informatics.jax.org/marker/MGI:98036; http://www.uniprot.org/uniprot/P61358
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='19942'),'mgi_symbol'] = 'Rpl27'
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='19942'),'mgi_id'] = 'MGI:98036'
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='19942'),'uniprotswissprot'] = 'P61358'
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='19942'),'description'] = 'ribosomal protein L27 [Source:MGI Symbol;Acc:MGI:98036]'
# https://www.ncbi.nlm.nih.gov/gene/19944; http://www.informatics.jax.org/marker/MGI:99687; http://www.uniprot.org/uniprot/P47915
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='19944'),'mgi_symbol'] = 'Rpl29'
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='19944'),'mgi_id'] = 'MGI:99687'
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='19944'),'uniprotswissprot'] = 'P47915'
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='19944'),'description'] = 'ribosomal protein L29 [Source:MGI Symbol;Acc:MGI:99687]'
# https://www.ncbi.nlm.nih.gov/gene/20104; http://www.informatics.jax.org/marker/MGI:98159; http://www.uniprot.org/uniprot/P62754
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='20104'),'mgi_symbol'] = 'Rps6'
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='20104'),'mgi_id'] = 'MGI:98159'
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='20104'),'uniprotswissprot'] = 'P62754'
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='20104'),'description'] = 'ribosomal protein S6 [Source:MGI Symbol;Acc:MGI:98159]'
# https://www.ncbi.nlm.nih.gov/gene/20116; http://www.informatics.jax.org/marker/MGI:98166; http://www.uniprot.org/uniprot/P62242
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='20116'),'mgi_symbol'] = 'Rps8'
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='20116'),'mgi_id'] = 'MGI:98166'
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='20116'),'uniprotswissprot'] = 'P62242'
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='20116'),'description'] = 'ribosomal protein S8 [Source:MGI Symbol;Acc:MGI:98166]'
# https://www.ncbi.nlm.nih.gov/gene/22186; http://www.informatics.jax.org/marker/MGI:98887; http://www.uniprot.org/uniprot/P62984
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='22186'),'mgi_symbol'] = 'Uba52'
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='22186'),'mgi_id'] = 'MGI:98887'
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='22186'),'uniprotswissprot'] = 'P62984'
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='22186'),'description'] = 'ubiquitin A-52 residue ribosomal protein fusion product 1 [Source:MGI Symbol;Acc:MGI:98887]'
# https://www.ncbi.nlm.nih.gov/gene/72157; http://www.informatics.jax.org/marker/MGI:97565; http://www.uniprot.org/uniprot/Q9D0F9
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='72157'),'mgi_symbol'] = 'Pgm2'
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='72157'),'mgi_id'] = 'MGI:97565'
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='72157'),'uniprotswissprot'] = 'Q9D0F9'
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='72157'),'description'] = 'phosphoglucomutase 2 [Source:MGI Symbol;Acc:MGI:97565]'
# https://www.ncbi.nlm.nih.gov/gene/100862072; http://www.informatics.jax.org/marker/MGI:5434806; http://www.uniprot.org/uniprot/P58022
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='100862072'),'mgi_symbol'] = 'Gm21451'
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='100862072'),'mgi_id'] = 'MGI:5434806'
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='100862072'),'uniprotswissprot'] = 'P58022'
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='100862072'),'description'] = 'predicted gene, 21451 [Source:MGI Symbol;Acc:MGI:5434806]'
# https://www.ncbi.nlm.nih.gov/gene/100862258; http://www.informatics.jax.org/marker/MGI:5434951; http://www.uniprot.org/uniprot/P63158
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='100862258'),'mgi_symbol'] = 'Gm21596'
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='100862258'),'mgi_id'] = 'MGI:5434951'
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='100862258'),'uniprotswissprot'] = 'P63158'
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='100862258'),'description'] = 'predicted gene, 21596 [Source:MGI Symbol;Acc:MGI:5434951]'
# https://www.ncbi.nlm.nih.gov/gene/101488212; http://www.informatics.jax.org/marker/MGI:5439444; http://www.uniprot.org/uniprot/Q8VD58
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='101488212'),'mgi_symbol'] = 'Evi2'
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='101488212'),'mgi_id'] = 'MGI:5439444'
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='101488212'),'uniprotswissprot'] = 'Q8VD58'
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='101488212'),'description'] = 'ecotropic viral integration site 2 [Source:MGI Symbol;Acc:MGI:5439444]'

# https://www.ncbi.nlm.nih.gov/gene/15018; http://www.informatics.jax.org/marker/MGI:95936; http://www.uniprot.org/uniprot/P14429
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='15018'),'mgi_symbol'] = 'H2-Q7'
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='15018'),'mgi_id'] = 'MGI:95936'
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='15018'),'uniprotswissprot'] = 'P14429'
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='15018'),'description'] = 'histocompatibility 2, Q region locus 7 [Source:MGI Symbol;Acc:MGI:95936]'
# https://www.ncbi.nlm.nih.gov/gene/110558; http://www.informatics.jax.org/marker/MGI:95938; http://www.uniprot.org/uniprot/P14431
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='110558'),'mgi_symbol'] = 'H2-Q9'
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='110558'),'mgi_id'] = 'MGI:95938'
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='110558'),'uniprotswissprot'] = 'P14431'
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='110558'),'description'] = 'histocompatibility 2, Q region locus 9 [Source:MGI Symbol;Acc:MGI:95938]'

ExoCarta_Protein.mgi_id.table[,'description'] = trimws(ExoCarta_Protein.mgi_id.table[,'description'])
ExoCarta_Protein.mgi_id.table[,'description'] = sapply(strsplit(ExoCarta_Protein.mgi_id.table[,'description'],'\\|'), FUN=function(U) paste(unique(sapply(strsplit(unlist(U),'\\[Source:'), FUN=function(V) trimws(unlist(V))[1])),collapse="|"))
ExoCarta_Protein.mgi_id.table = unique(ExoCarta_Protein.mgi_id.table)

ExoCarta_Protein.mgi_id.table = ExoCarta_Protein.mgi_id.table[,c('mgi_symbol','mgi_id','entrezgene','uniprotswissprot','description')]
ExoCarta_Protein.mgi_id.table = ExoCarta_Protein.mgi_id.table[order(match(ExoCarta_Protein.mgi_id.table[,'entrezgene'], a)),]

# merge multiple uniprotswissprot of single gene
unique.mgi_id.idx <- which(!duplicated(ExoCarta_Protein.mgi_id.table[,'mgi_id']))
unique.mgi_id <- ExoCarta_Protein.mgi_id.table[unique.mgi_id.idx,'mgi_id']
duplicated.mgi_id.idx <- which(duplicated(ExoCarta_Protein.mgi_id.table[,'mgi_id']))
duplicated.mgi_id <- ExoCarta_Protein.mgi_id.table[duplicated.mgi_id.idx,'mgi_id']
duplicated.mgi_id.idx.pair <- unique(sapply(duplicated.mgi_id, FUN=function(X) which(ExoCarta_Protein.mgi_id.table[,'mgi_id'] %in% X)))
duplicated.mgi_id.idx.uniprotswissprot <- sapply(duplicated.mgi_id.idx.pair, FUN=function(X) paste(sort(unique(ExoCarta_Protein.mgi_id.table[X,'uniprotswissprot'][ExoCarta_Protein.mgi_id.table[X,'uniprotswissprot']!=''])),collapse=';'))
if(length(duplicated.mgi_id.idx.pair) > 0)
{
	for(i in 1:length(duplicated.mgi_id.idx.pair))
	{
		ExoCarta_Protein.mgi_id.table[unlist(duplicated.mgi_id.idx.pair[i]),'uniprotswissprot'] = duplicated.mgi_id.idx.uniprotswissprot[i]
	}
}

# uniprotswissprot not found by R package "biomaRt"
empty.uniprotswissprot.idx <- which(ExoCarta_Protein.mgi_id.table[,'uniprotswissprot']=='')
empty.uniprotswissprot <- ExoCarta_Protein.mgi_id.table[empty.uniprotswissprot.idx,]

ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='11803'),'uniprotswissprot'] = 'Q03157'		# https://www.ncbi.nlm.nih.gov/gene/11803	http://www.informatics.jax.org/marker/MGI:88046		http://www.uniprot.org/uniprot/Q03157
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='12304'),'uniprotswissprot'] = 'P08003'		# https://www.ncbi.nlm.nih.gov/gene/12304	http://www.informatics.jax.org/marker/MGI:104864	http://www.uniprot.org/uniprot/P08003
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='12332'),'uniprotswissprot'] = 'P24452'		# https://www.ncbi.nlm.nih.gov/gene/12332	http://www.informatics.jax.org/marker/MGI:1098259	http://www.uniprot.org/uniprot/P24452
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='12340'),'uniprotswissprot'] = 'P47753'		# https://www.ncbi.nlm.nih.gov/gene/12340	http://www.informatics.jax.org/marker/MGI:106227	http://www.uniprot.org/uniprot/P47753
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='12505'),'uniprotswissprot'] = 'P15379'		# https://www.ncbi.nlm.nih.gov/gene/12505	http://www.informatics.jax.org/marker/MGI:88338		http://www.uniprot.org/uniprot/P15379
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='13211'),'uniprotswissprot'] = 'O70133'		# https://www.ncbi.nlm.nih.gov/gene/13211	http://www.informatics.jax.org/marker/MGI:108177	http://www.uniprot.org/uniprot/O70133
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='14325'),'uniprotswissprot'] = 'P29391'		# https://www.ncbi.nlm.nih.gov/gene/14325	http://www.informatics.jax.org/marker/MGI:95589		http://www.uniprot.org/uniprot/P29391
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='14337'),'uniprotswissprot'] = 'P49945'		# https://www.ncbi.nlm.nih.gov/gene/14337	http://www.informatics.jax.org/marker/MGI:95590		http://www.uniprot.org/uniprot/P49945
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='14805'),'uniprotswissprot'] = 'Q60934'		# https://www.ncbi.nlm.nih.gov/gene/14805	http://www.informatics.jax.org/marker/MGI:95814		http://www.uniprot.org/uniprot/Q60934
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='14991'),'uniprotswissprot'] = 'Q31093'		# https://www.ncbi.nlm.nih.gov/gene/14991	http://www.informatics.jax.org/marker/MGI:95915		http://www.uniprot.org/uniprot/Q31093
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='15122'),'uniprotswissprot'] = 'P01942'		# https://www.ncbi.nlm.nih.gov/gene/15122	http://www.informatics.jax.org/marker/MGI:96015		http://www.uniprot.org/uniprot/P01942
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='15135'),'uniprotswissprot'] = 'P02104'		# https://www.ncbi.nlm.nih.gov/gene/15135	http://www.informatics.jax.org/marker/MGI:96027		http://www.uniprot.org/uniprot/P02104
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='15525'),'uniprotswissprot'] = 'Q61316'		# https://www.ncbi.nlm.nih.gov/gene/15525	http://www.informatics.jax.org/marker/MGI:1342292	http://www.uniprot.org/uniprot/Q61316
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='15530'),'uniprotswissprot'] = 'Q05793'		# https://www.ncbi.nlm.nih.gov/gene/15530	http://www.informatics.jax.org/marker/MGI:96257		http://www.uniprot.org/uniprot/Q05793
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='16401'),'uniprotswissprot'] = 'Q00651'		# https://www.ncbi.nlm.nih.gov/gene/16401	http://www.informatics.jax.org/marker/MGI:96603		http://www.uniprot.org/uniprot/Q00651
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='16404'),'uniprotswissprot'] = 'Q61738'		# https://www.ncbi.nlm.nih.gov/gene/16404	http://www.informatics.jax.org/marker/MGI:102700	http://www.uniprot.org/uniprot/Q61738
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='16409'),'uniprotswissprot'] = 'P05555'		# https://www.ncbi.nlm.nih.gov/gene/16409	http://www.informatics.jax.org/marker/MGI:96607		http://www.uniprot.org/uniprot/P05555
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='16414'),'uniprotswissprot'] = 'P11835'		# https://www.ncbi.nlm.nih.gov/gene/16414	http://www.informatics.jax.org/marker/MGI:96611		http://www.uniprot.org/uniprot/P11835
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='16425'),'uniprotswissprot'] = 'Q61703'		# https://www.ncbi.nlm.nih.gov/gene/16425	http://www.informatics.jax.org/marker/MGI:96619		http://www.uniprot.org/uniprot/Q61703
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='16451'),'uniprotswissprot'] = 'P52332'		# https://www.ncbi.nlm.nih.gov/gene/16451	http://www.informatics.jax.org/marker/MGI:96628		http://www.uniprot.org/uniprot/P52332
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='16452'),'uniprotswissprot'] = 'Q62120'		# https://www.ncbi.nlm.nih.gov/gene/16452	http://www.informatics.jax.org/marker/MGI:96629		http://www.uniprot.org/uniprot/Q62120
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='16661'),'uniprotswissprot'] = 'P02535'		# https://www.ncbi.nlm.nih.gov/gene/16661	http://www.informatics.jax.org/marker/MGI:96685		http://www.uniprot.org/uniprot/P02535
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='16665'),'uniprotswissprot'] = 'Q61414'		# https://www.ncbi.nlm.nih.gov/gene/16665	http://www.informatics.jax.org/marker/MGI:96689		http://www.uniprot.org/uniprot/Q61414
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='16688'),'uniprotswissprot'] = 'Q9Z331'		# https://www.ncbi.nlm.nih.gov/gene/16688	http://www.informatics.jax.org/marker/MGI:1333768	http://www.uniprot.org/uniprot/Q9Z331
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='16728'),'uniprotswissprot'] = 'P11627'		# https://www.ncbi.nlm.nih.gov/gene/16728	http://www.informatics.jax.org/marker/MGI:96721		http://www.uniprot.org/uniprot/P11627
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='16854'),'uniprotswissprot'] = 'P16110'		# https://www.ncbi.nlm.nih.gov/gene/16854	http://www.informatics.jax.org/marker/MGI:96778		http://www.uniprot.org/uniprot/P16110
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='16912'),'uniprotswissprot'] = 'P28076'		# https://www.ncbi.nlm.nih.gov/gene/16912	http://www.informatics.jax.org/marker/MGI:1346526	http://www.uniprot.org/uniprot/P28076
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='16971'),'uniprotswissprot'] = 'Q91ZX7'		# https://www.ncbi.nlm.nih.gov/gene/16971	http://www.informatics.jax.org/marker/MGI:96828		http://www.uniprot.org/uniprot/Q91ZX7
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='17228'),'uniprotswissprot'] = 'P21844'		# https://www.ncbi.nlm.nih.gov/gene/17228	http://www.informatics.jax.org/marker/MGI:96941		http://www.uniprot.org/uniprot/P21844
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='17441'),'uniprotswissprot'] = 'Q61885'		# https://www.ncbi.nlm.nih.gov/gene/17441	http://www.informatics.jax.org/marker/MGI:97435		http://www.uniprot.org/uniprot/Q61885
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='18007'),'uniprotswissprot'] = 'P97798'		# https://www.ncbi.nlm.nih.gov/gene/18007	http://www.informatics.jax.org/marker/MGI:1097159	http://www.uniprot.org/uniprot/P97798
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='18040'),'uniprotswissprot'] = 'P08553'		# https://www.ncbi.nlm.nih.gov/gene/18040	http://www.informatics.jax.org/marker/MGI:97314		http://www.uniprot.org/uniprot/P08553
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='18176'),'uniprotswissprot'] = 'P08556'		# https://www.ncbi.nlm.nih.gov/gene/18176	http://www.informatics.jax.org/marker/MGI:97376		http://www.uniprot.org/uniprot/P08556
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='18596'),'uniprotswissprot'] = 'P05622'		# https://www.ncbi.nlm.nih.gov/gene/18596	http://www.informatics.jax.org/marker/MGI:97531		http://www.uniprot.org/uniprot/P05622
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='18770'),'uniprotswissprot'] = 'P53657'		# https://www.ncbi.nlm.nih.gov/gene/18770	http://www.informatics.jax.org/marker/MGI:97604		http://www.uniprot.org/uniprot/P53657
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='18787'),'uniprotswissprot'] = 'P22777'		# https://www.ncbi.nlm.nih.gov/gene/18787	http://www.informatics.jax.org/marker/MGI:97608		http://www.uniprot.org/uniprot/P22777
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='19933'),'uniprotswissprot'] = 'O09167'		# https://www.ncbi.nlm.nih.gov/gene/19933	http://www.informatics.jax.org/marker/MGI:1278340	http://www.uniprot.org/uniprot/O09167
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='20042'),'uniprotswissprot'] = 'P63323'		# https://www.ncbi.nlm.nih.gov/gene/20042	http://www.informatics.jax.org/marker/MGI:98105		http://www.uniprot.org/uniprot/P63323
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='20103'),'uniprotswissprot'] = 'P97461'		# https://www.ncbi.nlm.nih.gov/gene/20103	http://www.informatics.jax.org/marker/MGI:1097682	http://www.uniprot.org/uniprot/P97461
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='20514'),'uniprotswissprot'] = 'P51912'		# https://www.ncbi.nlm.nih.gov/gene/20514	http://www.informatics.jax.org/marker/MGI:105305	http://www.uniprot.org/uniprot/P51912
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='20846'),'uniprotswissprot'] = 'P42225'		# https://www.ncbi.nlm.nih.gov/gene/20846	http://www.informatics.jax.org/marker/MGI:103063	http://www.uniprot.org/uniprot/P42225
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='21814'),'uniprotswissprot'] = 'O88393'		# https://www.ncbi.nlm.nih.gov/gene/21814	http://www.informatics.jax.org/marker/MGI:104637	http://www.uniprot.org/uniprot/O88393
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='21825'),'uniprotswissprot'] = 'P35441'		# https://www.ncbi.nlm.nih.gov/gene/21825	http://www.informatics.jax.org/marker/MGI:98737		http://www.uniprot.org/uniprot/P35441
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='22284'),'uniprotswissprot'] = 'P70398'		# https://www.ncbi.nlm.nih.gov/gene/22284	http://www.informatics.jax.org/marker/MGI:894681	http://www.uniprot.org/uniprot/P70398
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='24088'),'uniprotswissprot'] = 'Q9QUN7'		# https://www.ncbi.nlm.nih.gov/gene/24088	http://www.informatics.jax.org/marker/MGI:1346060	http://www.uniprot.org/uniprot/Q9QUN7
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='26364'),'uniprotswissprot'] = 'Q9Z0M6'		# https://www.ncbi.nlm.nih.gov/gene/26364	http://www.informatics.jax.org/marker/MGI:1347095	http://www.uniprot.org/uniprot/Q9Z0M6
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='26921'),'uniprotswissprot'] = 'P97820'		# https://www.ncbi.nlm.nih.gov/gene/26921	http://www.informatics.jax.org/marker/MGI:1349394	http://www.uniprot.org/uniprot/P97820
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='27374'),'uniprotswissprot'] = 'Q8CIG8'		# https://www.ncbi.nlm.nih.gov/gene/27374	http://www.informatics.jax.org/marker/MGI:1351645	http://www.uniprot.org/uniprot/Q8CIG8
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='30806'),'uniprotswissprot'] = 'P57110'		# https://www.ncbi.nlm.nih.gov/gene/30806	http://www.informatics.jax.org/marker/MGI:1353468	http://www.uniprot.org/uniprot/P57110
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='50908'),'uniprotswissprot'] = 'Q8CG14'		# https://www.ncbi.nlm.nih.gov/gene/50908	http://www.informatics.jax.org/marker/MGI:1355312	http://www.uniprot.org/uniprot/Q8CG14
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='52118'),'uniprotswissprot'] = 'Q8K094'		# https://www.ncbi.nlm.nih.gov/gene/52118	http://www.informatics.jax.org/marker/MGI:107741	http://www.uniprot.org/uniprot/Q8K094
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='53612'),'uniprotswissprot'] = 'O88384'		# https://www.ncbi.nlm.nih.gov/gene/53612	http://www.informatics.jax.org/marker/MGI:1855688	http://www.uniprot.org/uniprot/O88384
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='54217'),'uniprotswissprot'] = 'P47964'		# https://www.ncbi.nlm.nih.gov/gene/54217	http://www.informatics.jax.org/marker/MGI:1860603	http://www.uniprot.org/uniprot/P47964
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='56087'),'uniprotswissprot'] = 'D3YYQ8'		# https://www.ncbi.nlm.nih.gov/gene/56087	http://www.informatics.jax.org/marker/MGI:1860299	http://www.uniprot.org/uniprot/D3YYQ8
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='56752'),'uniprotswissprot'] = 'Q9JLJ2'		# https://www.ncbi.nlm.nih.gov/gene/56752	http://www.informatics.jax.org/marker/MGI:1861622	http://www.uniprot.org/uniprot/Q9JLJ2
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='66588'),'uniprotswissprot'] = 'Q9DBP5'		# https://www.ncbi.nlm.nih.gov/gene/66588	http://www.informatics.jax.org/marker/MGI:1913838	http://www.uniprot.org/uniprot/Q9DBP5
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='66656'),'uniprotswissprot'] = 'P57776'		# https://www.ncbi.nlm.nih.gov/gene/66656	http://www.informatics.jax.org/marker/MGI:1913906	http://www.uniprot.org/uniprot/P57776
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='67373'),'uniprotswissprot'] = 'Q9CPN9'		# https://www.ncbi.nlm.nih.gov/gene/67373	http://www.informatics.jax.org/marker/MGI:1914623	http://www.uniprot.org/uniprot/Q9CPN9
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='69748'),'uniprotswissprot'] = 'Q571I9'		# https://www.ncbi.nlm.nih.gov/gene/69748	http://www.informatics.jax.org/marker/MGI:1916998	http://www.uniprot.org/uniprot/Q571I9
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='70699'),'uniprotswissprot'] = 'B9EJ54'		# https://www.ncbi.nlm.nih.gov/gene/70699	http://www.informatics.jax.org/marker/MGI:2141625	http://www.uniprot.org/uniprot/B9EJ54
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='71824'),'uniprotswissprot'] = 'B9EHI3'		# https://www.ncbi.nlm.nih.gov/gene/71824	http://www.informatics.jax.org/marker/MGI:1919074	http://www.uniprot.org/uniprot/B9EHI3
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='71853'),'uniprotswissprot'] = 'Q922R8'		# https://www.ncbi.nlm.nih.gov/gene/71853	http://www.informatics.jax.org/marker/MGI:1919103	http://www.uniprot.org/uniprot/Q922R8
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='78388'),'uniprotswissprot'] = 'Q9EQK5'		# https://www.ncbi.nlm.nih.gov/gene/78388	http://www.informatics.jax.org/marker/MGI:1925638	http://www.uniprot.org/uniprot/Q9EQK5
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='80910'),'uniprotswissprot'] = 'Q8CIM5'		# https://www.ncbi.nlm.nih.gov/gene/80910	http://www.informatics.jax.org/marker/MGI:1934129	http://www.uniprot.org/uniprot/Q8CIM5
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='99543'),'uniprotswissprot'] = 'Q8BK62'		# https://www.ncbi.nlm.nih.gov/gene/99543	http://www.informatics.jax.org/marker/MGI:1914877	http://www.uniprot.org/uniprot/Q8BK62
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='101706'),'uniprotswissprot'] = 'E9Q7G0'		# https://www.ncbi.nlm.nih.gov/gene/101706	http://www.informatics.jax.org/marker/MGI:2443665	http://www.uniprot.org/uniprot/E9Q7G0
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='106572'),'uniprotswissprot'] = 'Q921E2'		# https://www.ncbi.nlm.nih.gov/gene/106572	http://www.informatics.jax.org/marker/MGI:1914603	http://www.uniprot.org/uniprot/Q921E2
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='114228'),'uniprotswissprot'] = 'Q9Z1R9'		# https://www.ncbi.nlm.nih.gov/gene/114228	http://www.informatics.jax.org/marker/MGI:98839		http://www.uniprot.org/uniprot/Q9Z1R9
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='140559'),'uniprotswissprot'] = 'Q8R366'		# https://www.ncbi.nlm.nih.gov/gene/140559	http://www.informatics.jax.org/marker/MGI:2154090	http://www.uniprot.org/uniprot/Q8R366
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='210293'),'uniprotswissprot'] = 'Q8BZN6'		# https://www.ncbi.nlm.nih.gov/gene/210293	http://www.informatics.jax.org/marker/MGI:2146320	http://www.uniprot.org/uniprot/Q8BZN6
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='214290'),'uniprotswissprot'] = 'Q5BLK4'		# https://www.ncbi.nlm.nih.gov/gene/214290	http://www.informatics.jax.org/marker/MGI:2387179	http://www.uniprot.org/uniprot/Q5BLK4
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='215384'),'uniprotswissprot'] = 'E9Q9C6'		# https://www.ncbi.nlm.nih.gov/gene/215384	http://www.informatics.jax.org/marker/MGI:2444336	http://www.uniprot.org/uniprot/E9Q9C6
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='225348'),'uniprotswissprot'] = 'Q3TAQ9'		# https://www.ncbi.nlm.nih.gov/gene/225348	http://www.informatics.jax.org/marker/MGI:1917819	http://www.uniprot.org/uniprot/Q3TAQ9
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='226519'),'uniprotswissprot'] = 'P02468'		# https://www.ncbi.nlm.nih.gov/gene/226519	http://www.informatics.jax.org/marker/MGI:99914		http://www.uniprot.org/uniprot/P02468
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='239673'),'uniprotswissprot'] = 'E9Q1Z0'		# https://www.ncbi.nlm.nih.gov/gene/239673	http://www.informatics.jax.org/marker/MGI:3045312	http://www.uniprot.org/uniprot/E9Q1Z0
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='319163'),'uniprotswissprot'] = 'Q8CGP4'		# https://www.ncbi.nlm.nih.gov/gene/319163	http://www.informatics.jax.org/marker/MGI:2448285	http://www.uniprot.org/uniprot/Q8CGP4
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='319195'),'uniprotswissprot'] = 'Q9CPR4'		# https://www.ncbi.nlm.nih.gov/gene/319195	http://www.informatics.jax.org/marker/MGI:2448270	http://www.uniprot.org/uniprot/Q9CPR4
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='319482'),'uniprotswissprot'] = 'E9PVG8'		# https://www.ncbi.nlm.nih.gov/gene/319482	http://www.informatics.jax.org/marker/MGI:2442118	http://www.uniprot.org/uniprot/E9PVG8
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='353204'),'uniprotswissprot'] = 'A6ZI46'		# https://www.ncbi.nlm.nih.gov/gene/353204	http://www.informatics.jax.org/marker/MGI:2447811	http://www.uniprot.org/uniprot/A6ZI46
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='404710'),'uniprotswissprot'] = 'F8VQ29'		# https://www.ncbi.nlm.nih.gov/gene/404710	http://www.informatics.jax.org/marker/MGI:3028642	http://www.uniprot.org/uniprot/F8VQ29
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='623131'),'uniprotswissprot'] = 'B2RW88'		# https://www.ncbi.nlm.nih.gov/gene/623131	http://www.informatics.jax.org/marker/MGI:3648539	http://www.uniprot.org/uniprot/B2RW88
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='670464'),'uniprotswissprot'] = 'B1AQB0'		# https://www.ncbi.nlm.nih.gov/gene/670464	http://www.informatics.jax.org/marker/MGI:3652177	http://www.uniprot.org/uniprot/B1AQB0
ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']=='670496'),'uniprotswissprot'] = 'A2A4L9'		# https://www.ncbi.nlm.nih.gov/gene/670496	http://www.informatics.jax.org/marker/MGI:3650329	http://www.uniprot.org/uniprot/A2A4L9

ExoCarta_Protein.mgi_id.table = unique(ExoCarta_Protein.mgi_id.table)
#dim(ExoCarta_Protein.mgi_id.table)

# remove Gapdh-ps15 (https://www.ncbi.nlm.nih.gov/gene/100042025)
ExoCarta_Protein.mgi_id.table = ExoCarta_Protein.mgi_id.table[which(ExoCarta_Protein.mgi_id.table[,'entrezgene']!='100042025'),]

rownames(ExoCarta_Protein.mgi_id.table) = c(1:nrow(ExoCarta_Protein.mgi_id.table))
ExoCarta_Protein.mgi_id.table.New = ExoCarta_Protein.mgi_id.table
#ExoCarta_Protein.mgi_id.table
#ExoCarta_Protein.mgi_id.table.New

#mouse.UniProt.Entry2ProteinNames = sapply(strsplit(ExoCarta_Protein.mgi_id.table[,'uniprotswissprot'],";"), FUN=function(X) paste(unique(sapply(unlist(X), FUN=function(Y) ifelse(length(which(mouse.UniProt.reviewed.Entry %in% unlist(Y)))==0, ifelse(length(which(mouse.UniProt.unreviewed.Entry %in% unlist(Y)))==0, '', mouse.UniProt.unreviewed.ProteinNames[which(mouse.UniProt.unreviewed.Entry %in% unlist(Y))]), mouse.UniProt.reviewed.ProteinNames[which(mouse.UniProt.reviewed.Entry %in% unlist(Y))]))),collapse="|"))
#ExoCarta_Protein.mgi_id.table = cbind(ExoCarta_Protein.mgi_id.table, mouse.UniProt.Entry2ProteinNames)
#colnames(ExoCarta_Protein.mgi_id.table)[ncol(ExoCarta_Protein.mgi_id.table)] = 'Protein_names'

#write.csv(ExoCarta_Protein.mgi_id.table, paste('./Data/',Sample,'/10_datasetOverlap/Database/mmu/ExoCarta_All.Protein_Gene2MGI.table.csv',sep=""), quote=TRUE, row.names=TRUE)

# entrezgene: attributes not found by R package "biomaRt"
ExoCarta_Protein.entrezgene.NA_attributes <- a[which(! a %in% ExoCarta_Protein.mgi_id.table[,'entrezgene'])]
if(length(ExoCarta_Protein.entrezgene.NA_attributes) > 0)
{
	manually_curated.mgi_id.table <- matrix('',length(ExoCarta_Protein.entrezgene.NA_attributes),ncol(ExoCarta_Protein.mgi_id.table))
	colnames(manually_curated.mgi_id.table) = colnames(ExoCarta_Protein.mgi_id.table)
	
	ExoCarta_Protein.entrezgene.NA_attributes.uncurated <- c()
	
	for(i in 1:length(ExoCarta_Protein.entrezgene.NA_attributes))
	{
		# --- Gene type: protein coding -----------------------------------------------------------
		if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '12466') {				# https://www.ncbi.nlm.nih.gov/gene/12466; http://www.informatics.jax.org/marker/MGI:107943; http://www.uniprot.org/uniprot/P80317
			manually_curated.mgi_id.table[i,'mgi_symbol'] = 'Cct6a'
			manually_curated.mgi_id.table[i,'mgi_id'] = 'MGI:107943'
			manually_curated.mgi_id.table[i,'entrezgene'] = '12466'
			manually_curated.mgi_id.table[i,'uniprotswissprot'] = 'P80317'
			manually_curated.mgi_id.table[i,'description'] = 'chaperonin containing Tcp1, subunit 6a (zeta) [Source:MGI Symbol;Acc:MGI:107943]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '14980') {		# https://www.ncbi.nlm.nih.gov/gene/14980; http://www.informatics.jax.org/marker/MGI:95912; http://www.uniprot.org/uniprot/P01897
			manually_curated.mgi_id.table[i,'mgi_symbol'] = 'H2-L'
			manually_curated.mgi_id.table[i,'mgi_id'] = 'MGI:95912'
			manually_curated.mgi_id.table[i,'entrezgene'] = '14980'
			manually_curated.mgi_id.table[i,'uniprotswissprot'] = 'P01897'
			manually_curated.mgi_id.table[i,'description'] = 'histocompatibility 2, D region locus L [Source:MGI Symbol;Acc:MGI:95912]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '15129') {		# https://www.ncbi.nlm.nih.gov/gene/15129; http://www.informatics.jax.org/marker/MGI:96021; http://www.uniprot.org/uniprot/P02088
			manually_curated.mgi_id.table[i,'mgi_symbol'] = 'Hbb-b1'
			manually_curated.mgi_id.table[i,'mgi_id'] = 'MGI:96021'
			manually_curated.mgi_id.table[i,'entrezgene'] = '15129'
			manually_curated.mgi_id.table[i,'uniprotswissprot'] = 'P02088'
			manually_curated.mgi_id.table[i,'description'] = 'hemoglobin, beta adult major chain [Source:MGI Symbol;Acc:MGI:96021]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '15130') {		# https://www.ncbi.nlm.nih.gov/gene/15130; http://www.informatics.jax.org/marker/MGI:96022; http://www.uniprot.org/uniprot/P02089
			manually_curated.mgi_id.table[i,'mgi_symbol'] = 'Hbb-b2'
			manually_curated.mgi_id.table[i,'mgi_id'] = 'MGI:96022'
			manually_curated.mgi_id.table[i,'entrezgene'] = '15130'
			manually_curated.mgi_id.table[i,'uniprotswissprot'] = 'P02089'
			manually_curated.mgi_id.table[i,'description'] = 'hemoglobin, beta adult minor chain [Source:MGI Symbol;Acc:MGI:96022]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '16402') {		# https://www.ncbi.nlm.nih.gov/gene/16402; http://www.informatics.jax.org/marker/MGI:96604; http://www.uniprot.org/uniprot/P11688
			manually_curated.mgi_id.table[i,'mgi_symbol'] = 'Itga5'
			manually_curated.mgi_id.table[i,'mgi_id'] = 'MGI:96604'
			manually_curated.mgi_id.table[i,'entrezgene'] = '16402'
			manually_curated.mgi_id.table[i,'uniprotswissprot'] = 'P11688'
			manually_curated.mgi_id.table[i,'description'] = 'integrin alpha 5 (fibronectin receptor alpha) [Source:MGI Symbol;Acc:MGI:96604]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '66914') {		# https://www.ncbi.nlm.nih.gov/gene/66914; http://www.informatics.jax.org/marker/MGI:1914164; http://www.uniprot.org/uniprot/Q9D1C8
			manually_curated.mgi_id.table[i,'mgi_symbol'] = 'Vps28'
			manually_curated.mgi_id.table[i,'mgi_id'] = 'MGI:1914164'
			manually_curated.mgi_id.table[i,'entrezgene'] = '66914'
			manually_curated.mgi_id.table[i,'uniprotswissprot'] = 'Q9D1C8'
			manually_curated.mgi_id.table[i,'description'] = 'vacuolar protein sorting 28 [Source:MGI Symbol;Acc:MGI:1914164]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '67332') {		# https://www.ncbi.nlm.nih.gov/gene/67332; http://www.informatics.jax.org/marker/MGI:1914582; http://www.uniprot.org/uniprot/P62320
			manually_curated.mgi_id.table[i,'mgi_symbol'] = 'Snrpd3'
			manually_curated.mgi_id.table[i,'mgi_id'] = 'MGI:1914582'
			manually_curated.mgi_id.table[i,'entrezgene'] = '67332'
			manually_curated.mgi_id.table[i,'uniprotswissprot'] = 'P62320'
			manually_curated.mgi_id.table[i,'description'] = 'small nuclear ribonucleoprotein D3 [Source:MGI Symbol;Acc:MGI:1914582]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '100042025') {	# https://www.ncbi.nlm.nih.gov/gene/100042025; http://www.informatics.jax.org/marker/MGI:5434255; http://www.uniprot.org/uniprot/P16858
			manually_curated.mgi_id.table[i,'mgi_symbol'] = 'Gapdh-ps15'
			manually_curated.mgi_id.table[i,'mgi_id'] = 'MGI:5434255'
			manually_curated.mgi_id.table[i,'entrezgene'] = '100042025'
			manually_curated.mgi_id.table[i,'uniprotswissprot'] = 'P16858'
			manually_curated.mgi_id.table[i,'description'] = 'glyceraldehyde-3-phosphate dehydrogenase, pseudogene 15 [Source:MGI Symbol;Acc:MGI:5434255]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '100862072') {	# https://www.ncbi.nlm.nih.gov/gene/100862072; http://www.informatics.jax.org/marker/MGI:5434806; http://www.uniprot.org/uniprot/P58022
			manually_curated.mgi_id.table[i,'mgi_symbol'] = 'Gm21451'
			manually_curated.mgi_id.table[i,'mgi_id'] = 'MGI:5434806'
			manually_curated.mgi_id.table[i,'entrezgene'] = '100862072'
			manually_curated.mgi_id.table[i,'uniprotswissprot'] = 'P58022'
			manually_curated.mgi_id.table[i,'description'] = 'predicted gene, 21451 [Source:MGI Symbol;Acc:MGI:5434806]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '100862258') {	# https://www.ncbi.nlm.nih.gov/gene/100862258; http://www.informatics.jax.org/marker/MGI:5434951; http://www.uniprot.org/uniprot/P63158
			manually_curated.mgi_id.table[i,'mgi_symbol'] = 'Gm21596'
			manually_curated.mgi_id.table[i,'mgi_id'] = 'MGI:5434951'
			manually_curated.mgi_id.table[i,'entrezgene'] = '100862258'
			manually_curated.mgi_id.table[i,'uniprotswissprot'] = 'P63158'
			manually_curated.mgi_id.table[i,'description'] = 'predicted gene, 21596 [Source:MGI Symbol;Acc:MGI:5434951]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '102641229') {	# https://www.ncbi.nlm.nih.gov/gene/102641229; http://www.uniprot.org/uniprot/P62806
			manually_curated.mgi_id.table[i,'mgi_symbol'] = 'LOC102641229'
			manually_curated.mgi_id.table[i,'mgi_id'] = ''
			manually_curated.mgi_id.table[i,'entrezgene'] = '102641229'
			manually_curated.mgi_id.table[i,'uniprotswissprot'] = 'P62806'				# http://exocarta.org/gene_summary?gene_id=102641229
			manually_curated.mgi_id.table[i,'description'] = 'histone H4 [Source:NCBI gene;Acc:102641229]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '102641332') {	# https://www.ncbi.nlm.nih.gov/gene/102641332; http://www.uniprot.org/uniprot/P62827
			manually_curated.mgi_id.table[i,'mgi_symbol'] = 'LOC102641332'
			manually_curated.mgi_id.table[i,'mgi_id'] = ''
			manually_curated.mgi_id.table[i,'entrezgene'] = '102641332'
			manually_curated.mgi_id.table[i,'uniprotswissprot'] = 'P62827'				# http://exocarta.org/gene_summary?gene_id=102641332
			manually_curated.mgi_id.table[i,'description'] = 'GTP-binding nuclear protein Ran [Source:NCBI gene;Acc:102641332]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '102641678') {	# https://www.ncbi.nlm.nih.gov/gene/102641678; http://www.uniprot.org/uniprot/Q9DB20
			manually_curated.mgi_id.table[i,'mgi_symbol'] = 'LOC102641678'
			manually_curated.mgi_id.table[i,'mgi_id'] = ''
			manually_curated.mgi_id.table[i,'entrezgene'] = '102641678'
			manually_curated.mgi_id.table[i,'uniprotswissprot'] = 'Q9DB20'				# http://exocarta.org/gene_summary?gene_id=102641678
			manually_curated.mgi_id.table[i,'description'] = 'ATP synthase subunit O, mitochondrial [Source:NCBI gene;Acc:102641678]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '102642689') {	# https://www.ncbi.nlm.nih.gov/gene/102642689; http://www.uniprot.org/uniprot/P14206
			manually_curated.mgi_id.table[i,'mgi_symbol'] = 'LOC102642689'
			manually_curated.mgi_id.table[i,'mgi_id'] = ''
			manually_curated.mgi_id.table[i,'entrezgene'] = '102642689'
			manually_curated.mgi_id.table[i,'uniprotswissprot'] = 'P14206'				# http://exocarta.org/gene_summary?gene_id=102642689
			manually_curated.mgi_id.table[i,'description'] = '40S ribosomal protein SA [Source:NCBI gene;Acc:102642689]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '102642938') {	# https://www.ncbi.nlm.nih.gov/gene/102642938; http://www.uniprot.org/uniprot/O88569
			manually_curated.mgi_id.table[i,'mgi_symbol'] = 'LOC102642938'
			manually_curated.mgi_id.table[i,'mgi_id'] = ''
			manually_curated.mgi_id.table[i,'entrezgene'] = '102642938'
			manually_curated.mgi_id.table[i,'uniprotswissprot'] = 'O88569'				# http://exocarta.org/gene_summary?gene_id=102642938
			manually_curated.mgi_id.table[i,'description'] = 'heterogeneous nuclear ribonucleoproteins A2/B1 [Source:NCBI gene;Acc:102642938]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '102643269') {	# https://www.ncbi.nlm.nih.gov/gene/102643269; http://www.uniprot.org/uniprot/P35550
			manually_curated.mgi_id.table[i,'mgi_symbol'] = 'LOC102643269'
			manually_curated.mgi_id.table[i,'mgi_id'] = ''
			manually_curated.mgi_id.table[i,'entrezgene'] = '102643269'
			manually_curated.mgi_id.table[i,'uniprotswissprot'] = 'P35550'				# http://exocarta.org/gene_summary?gene_id=102643269
			manually_curated.mgi_id.table[i,'description'] = 'rRNA 2\'-O-methyltransferase fibrillarin [Source:NCBI gene;Acc:102643269]'
		# --- Gene type: other --------------------------------------------------------------------
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '14352') {		# https://www.ncbi.nlm.nih.gov/gene/14352; http://www.informatics.jax.org/marker/MGI:95598; http://www.uniprot.org/uniprot/P11370
			manually_curated.mgi_id.table[i,'mgi_symbol'] = 'Fv4'
			manually_curated.mgi_id.table[i,'mgi_id'] = 'MGI:95598'
			manually_curated.mgi_id.table[i,'entrezgene'] = '14352'
			manually_curated.mgi_id.table[i,'uniprotswissprot'] = 'P11370'
			manually_curated.mgi_id.table[i,'description'] = 'Friend virus susceptibility 4 [Source:MGI Symbol;Acc:MGI:95598]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '15127') {		# https://www.ncbi.nlm.nih.gov/gene/15127; http://www.informatics.jax.org/marker/MGI:96020
			manually_curated.mgi_id.table[i,'mgi_symbol'] = 'Hbb'
			manually_curated.mgi_id.table[i,'mgi_id'] = 'MGI:96020'
			manually_curated.mgi_id.table[i,'entrezgene'] = '15127'
			manually_curated.mgi_id.table[i,'uniprotswissprot'] = ''
			manually_curated.mgi_id.table[i,'description'] = 'hemoglobin beta chain complex [Source:MGI Symbol;Acc:MGI:96020]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '17276') {		# https://www.ncbi.nlm.nih.gov/gene/17276; http://www.informatics.jax.org/marker/MGI:107565; http://www.uniprot.org/uniprot/P70355
			manually_curated.mgi_id.table[i,'mgi_symbol'] = 'Mela'
			manually_curated.mgi_id.table[i,'mgi_id'] = 'MGI:107565'
			manually_curated.mgi_id.table[i,'entrezgene'] = '17276'
			manually_curated.mgi_id.table[i,'uniprotswissprot'] = 'P70355'
			manually_curated.mgi_id.table[i,'description'] = 'melanoma antigen [Source:MGI Symbol;Acc:MGI:107565]'
		# --- Gene type: pseudo -------------------------------------------------------------------
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '14337') {		# https://www.ncbi.nlm.nih.gov/gene/14337; http://www.informatics.jax.org/marker/MGI:95590
			manually_curated.mgi_id.table[i,'mgi_symbol'] = 'Ftl2-ps'
			manually_curated.mgi_id.table[i,'mgi_id'] = 'MGI:95590'
			manually_curated.mgi_id.table[i,'entrezgene'] = '14337'
			manually_curated.mgi_id.table[i,'uniprotswissprot'] = ''
			manually_curated.mgi_id.table[i,'description'] = 'ferritin light polypeptide 2, pseudogene [Source:MGI Symbol;Acc:MGI:95590]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '103324') {		# https://www.ncbi.nlm.nih.gov/gene/103324; http://www.informatics.jax.org/marker/MGI:3645521
			manually_curated.mgi_id.table[i,'mgi_symbol'] = 'Gm4735'
			manually_curated.mgi_id.table[i,'mgi_id'] = 'MGI:3645521'
			manually_curated.mgi_id.table[i,'entrezgene'] = '103324'
			manually_curated.mgi_id.table[i,'uniprotswissprot'] = ''
			manually_curated.mgi_id.table[i,'description'] = 'predicted gene 4735 [Source:MGI Symbol;Acc:MGI:3645521]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '380687') {		# https://www.ncbi.nlm.nih.gov/gene/380687; http://www.informatics.jax.org/marker/MGI:3779464
			manually_curated.mgi_id.table[i,'mgi_symbol'] = 'Gm5138'
			manually_curated.mgi_id.table[i,'mgi_id'] = 'MGI:3779464'
			manually_curated.mgi_id.table[i,'entrezgene'] = '380687'
			manually_curated.mgi_id.table[i,'uniprotswissprot'] = ''
			manually_curated.mgi_id.table[i,'description'] = 'predicted gene 5138 [Source:MGI Symbol;Acc:MGI:3779464]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '433464') {		# https://www.ncbi.nlm.nih.gov/gene/433464; http://www.informatics.jax.org/marker/MGI:3651689
			manually_curated.mgi_id.table[i,'mgi_symbol'] = 'Gm13991'
			manually_curated.mgi_id.table[i,'mgi_id'] = 'MGI:3651689'
			manually_curated.mgi_id.table[i,'entrezgene'] = '433464'
			manually_curated.mgi_id.table[i,'uniprotswissprot'] = ''
			manually_curated.mgi_id.table[i,'description'] = 'predicted gene 13991 [Source:MGI Symbol;Acc:MGI:3651689]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '433594') {		# https://www.ncbi.nlm.nih.gov/gene/433594; http://www.informatics.jax.org/marker/MGI:3643179
			manually_curated.mgi_id.table[i,'mgi_symbol'] = 'Gm5537'
			manually_curated.mgi_id.table[i,'mgi_id'] = 'MGI:3643179'
			manually_curated.mgi_id.table[i,'entrezgene'] = '433594'
			manually_curated.mgi_id.table[i,'uniprotswissprot'] = ''
			manually_curated.mgi_id.table[i,'description'] = 'predicted gene 5537 [Source:MGI Symbol;Acc:MGI:3643179]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '433941') {		# https://www.ncbi.nlm.nih.gov/gene/433941; http://www.informatics.jax.org/marker/MGI:3804971
			manually_curated.mgi_id.table[i,'mgi_symbol'] = 'Gm5561'
			manually_curated.mgi_id.table[i,'mgi_id'] = 'MGI:3804971'
			manually_curated.mgi_id.table[i,'entrezgene'] = '433941'
			manually_curated.mgi_id.table[i,'uniprotswissprot'] = ''
			manually_curated.mgi_id.table[i,'description'] = 'predicted gene 5561 [Source:MGI Symbol;Acc:MGI:3804971]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '546298') {		# https://www.ncbi.nlm.nih.gov/gene/546298; http://www.informatics.jax.org/marker/MGI:3705640
			manually_curated.mgi_id.table[i,'mgi_symbol'] = 'Rps2-ps13'
			manually_curated.mgi_id.table[i,'mgi_id'] = 'MGI:3705640'
			manually_curated.mgi_id.table[i,'entrezgene'] = '546298'
			manually_curated.mgi_id.table[i,'uniprotswissprot'] = ''
			manually_curated.mgi_id.table[i,'description'] = 'ribosomal protein S2, pseudogene 13 [Source:MGI Symbol;Acc:MGI:3705640]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '619900') {		# https://www.ncbi.nlm.nih.gov/gene/619900; http://www.informatics.jax.org/marker/MGI:3643614
			manually_curated.mgi_id.table[i,'mgi_symbol'] = 'Rps27a-ps2'
			manually_curated.mgi_id.table[i,'mgi_id'] = 'MGI:3643614'
			manually_curated.mgi_id.table[i,'entrezgene'] = '619900'
			manually_curated.mgi_id.table[i,'uniprotswissprot'] = ''
			manually_curated.mgi_id.table[i,'description'] = 'ribosomal protein S27A, pseudogene 2 [Source:MGI Symbol;Acc:MGI:3643614]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '621100') {		# https://www.ncbi.nlm.nih.gov/gene/621100; http://www.informatics.jax.org/marker/MGI:3646174
			manually_curated.mgi_id.table[i,'mgi_symbol'] = 'Rpl27-ps3'
			manually_curated.mgi_id.table[i,'mgi_id'] = 'MGI:3646174'
			manually_curated.mgi_id.table[i,'entrezgene'] = '621100'
			manually_curated.mgi_id.table[i,'uniprotswissprot'] = ''
			manually_curated.mgi_id.table[i,'description'] = 'ribosomal protein L27, pseudogene 3 [Source:MGI Symbol;Acc:MGI:3646174]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '625298') {		# https://www.ncbi.nlm.nih.gov/gene/625298; http://www.informatics.jax.org/marker/MGI:3647042
			manually_curated.mgi_id.table[i,'mgi_symbol'] = 'Rps13-ps1'
			manually_curated.mgi_id.table[i,'mgi_id'] = 'MGI:3647042'
			manually_curated.mgi_id.table[i,'entrezgene'] = '625298'
			manually_curated.mgi_id.table[i,'uniprotswissprot'] = ''
			manually_curated.mgi_id.table[i,'description'] = 'ribosomal protein S13, pseudogene 1 [Source:MGI Symbol;Acc:MGI:3647042]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '631069') {		# https://www.ncbi.nlm.nih.gov/gene/631069; http://www.informatics.jax.org/marker/MGI:3647164 (Error Page: No Marker Found)
			manually_curated.mgi_id.table[i,'mgi_symbol'] = 'Gm7053'
			manually_curated.mgi_id.table[i,'mgi_id'] = 'MGI:3647164'					# Error Page: No Marker Found
			manually_curated.mgi_id.table[i,'entrezgene'] = '631069'					# This record has been withdrawn by NCBI staff. This record represented a gene that is not currently annotated by NCBI.
			manually_curated.mgi_id.table[i,'uniprotswissprot'] = ''
			manually_curated.mgi_id.table[i,'description'] = 'predicted gene 7053 [Source:MGI Symbol;Acc:MGI:3647164]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '633082') {		# https://www.ncbi.nlm.nih.gov/gene/633082; http://www.informatics.jax.org/marker/MGI:3650637
			manually_curated.mgi_id.table[i,'mgi_symbol'] = 'Gm12164'
			manually_curated.mgi_id.table[i,'mgi_id'] = 'MGI:3650637'
			manually_curated.mgi_id.table[i,'entrezgene'] = '633082'
			manually_curated.mgi_id.table[i,'uniprotswissprot'] = ''
			manually_curated.mgi_id.table[i,'description'] = 'predicted gene 12164 [Source:MGI Symbol;Acc:MGI:3650637]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '668838') {		# https://www.ncbi.nlm.nih.gov/gene/668838; http://www.informatics.jax.org/marker/MGI:3643262; http://www.uniprot.org/uniprot/C5H0E8
			manually_curated.mgi_id.table[i,'mgi_symbol'] = 'Gm9392'
			manually_curated.mgi_id.table[i,'mgi_id'] = 'MGI:3643262'
			manually_curated.mgi_id.table[i,'entrezgene'] = '668838'
			manually_curated.mgi_id.table[i,'uniprotswissprot'] = 'C5H0E8'
			manually_curated.mgi_id.table[i,'description'] = 'predicted pseudogene 9392 [Source:MGI Symbol;Acc:MGI:3643262]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '100040053') {	# https://www.ncbi.nlm.nih.gov/gene/100040053; http://www.informatics.jax.org/marker/MGI:3780741
			manually_curated.mgi_id.table[i,'mgi_symbol'] = 'Gm2574'
			manually_curated.mgi_id.table[i,'mgi_id'] = 'MGI:3780741'
			manually_curated.mgi_id.table[i,'entrezgene'] = '100040053'
			manually_curated.mgi_id.table[i,'uniprotswissprot'] = ''
			manually_curated.mgi_id.table[i,'description'] = 'predicted pseudogene 2574 [Source:MGI Symbol;Acc:MGI:3780741]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '100041325') {	# https://www.ncbi.nlm.nih.gov/gene/100041325; http://www.informatics.jax.org/marker/MGI:3781450
			manually_curated.mgi_id.table[i,'mgi_symbol'] = 'Gm3272'
			manually_curated.mgi_id.table[i,'mgi_id'] = 'MGI:3781450'
			manually_curated.mgi_id.table[i,'entrezgene'] = '100041325'
			manually_curated.mgi_id.table[i,'uniprotswissprot'] = ''
			manually_curated.mgi_id.table[i,'description'] = 'predicted pseudogene 3272 [Source:MGI Symbol;Acc:MGI:3781450]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '100041342') {	# https://www.ncbi.nlm.nih.gov/gene/100041342; http://www.informatics.jax.org/marker/MGI:3650720
			manually_curated.mgi_id.table[i,'mgi_symbol'] = 'Gm12537'
			manually_curated.mgi_id.table[i,'mgi_id'] = 'MGI:3650720'
			manually_curated.mgi_id.table[i,'entrezgene'] = '100041342'
			manually_curated.mgi_id.table[i,'uniprotswissprot'] = ''
			manually_curated.mgi_id.table[i,'description'] = 'predicted gene 12537 [Source:MGI Symbol;Acc:MGI:3650720]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '100042349') {	# https://www.ncbi.nlm.nih.gov/gene/100042349; http://www.informatics.jax.org/marker/MGI:3708724
			manually_curated.mgi_id.table[i,'mgi_symbol'] = 'Gm10359'
			manually_curated.mgi_id.table[i,'mgi_id'] = 'MGI:3708724'
			manually_curated.mgi_id.table[i,'entrezgene'] = '100042349'
			manually_curated.mgi_id.table[i,'uniprotswissprot'] = ''
			manually_curated.mgi_id.table[i,'description'] = 'predicted gene 10359 [Source:MGI Symbol;Acc:MGI:3708724]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '100042746') {	# https://www.ncbi.nlm.nih.gov/gene/100042746; http://www.informatics.jax.org/marker/MGI:3650772
			manually_curated.mgi_id.table[i,'mgi_symbol'] = 'Gm12033'
			manually_curated.mgi_id.table[i,'mgi_id'] = 'MGI:3650772'
			manually_curated.mgi_id.table[i,'entrezgene'] = '100042746'
			manually_curated.mgi_id.table[i,'uniprotswissprot'] = ''
			manually_curated.mgi_id.table[i,'description'] = 'predicted gene 12033 [Source:MGI Symbol;Acc:MGI:3650772]'
		} else if(ExoCarta_Protein.entrezgene.NA_attributes[i] == '102634437') {	# https://www.ncbi.nlm.nih.gov/gene/102634437; http://www.informatics.jax.org/marker/MGI:3650241
			manually_curated.mgi_id.table[i,'mgi_symbol'] = 'Ywhaq-ps3'
			manually_curated.mgi_id.table[i,'mgi_id'] = 'MGI:3650241'
			manually_curated.mgi_id.table[i,'entrezgene'] = '102634437'
			manually_curated.mgi_id.table[i,'uniprotswissprot'] = ''
			manually_curated.mgi_id.table[i,'description'] = 'tyrosine 3-monooxygenase/tryptophan 5-monooxygenase activation protein theta, pseudogene 3 [Source:MGI Symbol;Acc:MGI:3650241]'
		} else {
			cat(paste('entrezgene "',ExoCarta_Protein.entrezgene.NA_attributes[i],'" did not match queried attributes (\'mgi_symbol\',\'mgi_id\',\'uniprotswissprot\',\'description\') by using R package "biomaRt".',sep=""),"\n")
			
			ExoCarta_Protein.entrezgene.NA_attributes.uncurated <- c(ExoCarta_Protein.entrezgene.NA_attributes.uncurated, ExoCarta_Protein.entrezgene.NA_attributes[i])
		}
		
		manually_curated.mgi_id.table[i,'description'] = trimws(manually_curated.mgi_id.table[i,'description'])
		manually_curated.mgi_id.table[i,'description'] = sapply(strsplit(manually_curated.mgi_id.table[i,'description'],'\\|'), FUN=function(U) paste(unique(sapply(strsplit(unlist(U),'\\[Source:'), FUN=function(V) trimws(unlist(V))[1])),collapse="|"))
	}
	
	if(length(ExoCarta_Protein.entrezgene.NA_attributes.uncurated) > 0)
	{
		write.table(ExoCarta_Protein.entrezgene.NA_attributes.uncurated, paste('./Data/',Sample,'/10_datasetOverlap/Database/mmu/ExoCarta_All.Protein_Gene2MGI_NA.txt',sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE)
	}
	
	ExoCarta_Protein.mgi_id.table.New = rbind(ExoCarta_Protein.mgi_id.table, manually_curated.mgi_id.table)
}

rownames(ExoCarta_Protein.mgi_id.table.New) = c(1:nrow(ExoCarta_Protein.mgi_id.table.New))
#ExoCarta_Protein.mgi_id.table.New

mouse.UniProt.Entry2ProteinNames = sapply(strsplit(ExoCarta_Protein.mgi_id.table.New[,'uniprotswissprot'],";"), FUN=function(X) paste(unique(sapply(unlist(X), FUN=function(Y) ifelse(length(which(mouse.UniProt.reviewed.Entry %in% unlist(Y)))==0, ifelse(length(which(mouse.UniProt.unreviewed.Entry %in% unlist(Y)))==0, '', mouse.UniProt.unreviewed.ProteinNames[which(mouse.UniProt.unreviewed.Entry %in% unlist(Y))]), mouse.UniProt.reviewed.ProteinNames[which(mouse.UniProt.reviewed.Entry %in% unlist(Y))]))),collapse="|"))
ExoCarta_Protein.mgi_id.table.New = cbind(ExoCarta_Protein.mgi_id.table.New, mouse.UniProt.Entry2ProteinNames)
colnames(ExoCarta_Protein.mgi_id.table.New)[ncol(ExoCarta_Protein.mgi_id.table.New)] = 'Protein_names'

ExoCarta_Protein.mgi_id.table.New = ExoCarta_Protein.mgi_id.table.New[which(rowSums(ExoCarta_Protein.mgi_id.table.New=='')!=ncol(ExoCarta_Protein.mgi_id.table.New)),,drop=FALSE]
rownames(ExoCarta_Protein.mgi_id.table.New) = c(1:nrow(ExoCarta_Protein.mgi_id.table.New))

#write.csv(ExoCarta_Protein.mgi_id.table.New, paste('./Data/',Sample,'/10_datasetOverlap/Database/mmu/ExoCarta_All.Protein_Gene2MGI.table.csv',sep=""), quote=TRUE, row.names=TRUE)

# --- write xlsx file by openxlsx package ---------------------------------------------------------
# Importing a big xlsx file into R? https://stackoverflow.com/questions/19147884/importing-a-big-xlsx-file-into-r/43118530#43118530
# Building R for Windows: https://cran.r-project.org/bin/windows/Rtools/
#Sys.setenv("R_ZIPCMD" = "C:/Rtools/bin/zip.exe")	# path to zip.exe

wb <- openxlsx::createWorkbook(paste('./Data/',Sample,'/10_datasetOverlap/Database/mmu/ExoCarta_All.Protein_Gene2MGI.table.xlsx',sep=""))
modifyBaseFont(wb, fontSize=12, fontColour="black", fontName="Times New Roman")

sheet = 'ExoCarta_All'
sheetData = ExoCarta_Protein.mgi_id.table.New

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
	}
}

if(length(wb$sheet_names) > 0) {
	openxlsx::saveWorkbook(wb, file=paste('./Data/',Sample,'/10_datasetOverlap/Database/mmu/ExoCarta_All.Protein_Gene2MGI.table.xlsx',sep=""), overwrite=TRUE)
}

