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

#install.packages("gtools")				# gtools: Various R Programming Tools
library(gtools)							# mixedsort(), mixedorder()


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


# === input Top100 recorded in Microvesicle database ============================================================================
#Database <- c('ExoCarta','Vesiclepedia','EVpedia')

# --- EVpedia ---------------------------------------------------------------------------------------------------------
# /// mouse format ///
mainDir <- paste('./Data/',Sample,'/10_datasetOverlap/Database',sep="")
subDir <- "mmu"
if (file.exists(file.path(mainDir, subDir))){
} else {
	dir.create(file.path(mainDir, subDir))
}

EVpedia_Top100.Protein_Mouse <- read.csv(paste('./Database/Microvesicle/EVpedia/Top100_EVmarkers/Protein_Mouse.csv',sep=""), header=TRUE, stringsAsFactors=FALSE)

EVpedia_Top100.Protein_Mouse.Index <- trimws(EVpedia_Top100.Protein_Mouse[,'Index'])
EVpedia_Top100.Protein_Mouse.UniProt_accession <- trimws(EVpedia_Top100.Protein_Mouse[,'UniProt.accession'])
EVpedia_Top100.Protein_Mouse.UniProt_name <- trimws(EVpedia_Top100.Protein_Mouse[,'UniProt.name'])
EVpedia_Top100.Protein_Mouse.Status <- trimws(EVpedia_Top100.Protein_Mouse[,'Status'])
EVpedia_Top100.Protein_Mouse.Protein_name <- trimws(EVpedia_Top100.Protein_Mouse[,'Protein.name'])
EVpedia_Top100.Protein_Mouse.Gene_symbol <- trimws(EVpedia_Top100.Protein_Mouse[,'Gene.symbol'])
EVpedia_Top100.Protein_Mouse.Species <- trimws(EVpedia_Top100.Protein_Mouse[,'Species'])
EVpedia_Top100.Protein_Mouse.Identification_number <- trimws(EVpedia_Top100.Protein_Mouse[,'Identification.number'])

# --- add in missed Gene symbol or change '0' to '' ---
for(i in 1:nrow(EVpedia_Top100.Protein_Mouse))
{
	if(EVpedia_Top100.Protein_Mouse.Gene_symbol[i] == '0') {
		if(EVpedia_Top100.Protein_Mouse.Protein_name[i] != 'Deleted.') {
			if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'P01806') {				# http://www.uniprot.org/uniprot/P01806; Igh-VX24 (http://www.informatics.jax.org/marker/MGI:96492); https://www.ncbi.nlm.nih.gov/gene/195176
				EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = 'Igh-VX24'
			} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'P01837') {		# http://www.uniprot.org/uniprot/P01837; Igkc (http://www.informatics.jax.org/marker/MGI:96495); https://www.ncbi.nlm.nih.gov/gene/16071
				EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = 'Igkc'
			} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q3TYS2') {		# http://www.uniprot.org/uniprot/Q3TYS2; BC017643 (https://www.ncbi.nlm.nih.gov/gene/217370); http://www.informatics.jax.org/marker/MGI:2384959; Mm.27182. (https://www.ncbi.nlm.nih.gov/UniGene/clust.cgi?ORG=Mm&CID=27182); Mm.490051.
				EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = 'BC017643'
			} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q3U6N9') {		# http://www.uniprot.org/uniprot/Q3U6N9; 1110038F14Rik (https://www.ncbi.nlm.nih.gov/gene/117171); http://www.informatics.jax.org/marker/MGI:2152337; https://www.ncbi.nlm.nih.gov/UniGene/clust.cgi?ORG=Mm&CID=322968; https://www.ncbi.nlm.nih.gov/UniGene/clust.cgi?ORG=Mm&CID=491343; https://www.ncbi.nlm.nih.gov/UniGene/clust.cgi?ORG=Mm&CID=480608
				EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = '1110038F14Rik'
			} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q3UE31') {		# http://www.uniprot.org/uniprot/Q3UE31; 5031439G07Rik (https://www.ncbi.nlm.nih.gov/gene/223739); http://www.informatics.jax.org/marker/MGI:2444899; https://www.ncbi.nlm.nih.gov/UniGene/clust.cgi?ORG=Mm&CID=323925
				EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = '5031439G07Rik'
			} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q3USZ8') {		# http://www.uniprot.org/uniprot/Q3USZ8; 1190002N15Rik (https://www.ncbi.nlm.nih.gov/gene/68861); http://www.informatics.jax.org/marker/MGI:1916111; https://www.ncbi.nlm.nih.gov/UniGene/clust.cgi?ORG=Mm&CID=258746
				EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = '1190002N15Rik'
			} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q3UZV7') {		# http://www.uniprot.org/uniprot/Q3UZV7; 9330182L06Rik (https://www.ncbi.nlm.nih.gov/gene/231014); http://www.informatics.jax.org/marker/MGI:2443264; https://www.ncbi.nlm.nih.gov/UniGene/clust.cgi?ORG=Mm&CID=130063
				EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = '9330182L06Rik'
			} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q5XFZ0') {		# http://www.uniprot.org/uniprot/Q5XFZ0; 2700062C07Rik (https://www.ncbi.nlm.nih.gov/gene/68046); http://www.informatics.jax.org/marker/MGI:1915296; https://www.ncbi.nlm.nih.gov/UniGene/clust.cgi?ORG=Mm&CID=276574
				EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = '2700062C07Rik'
			} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q66JV7') {		# http://www.uniprot.org/uniprot/Q66JV7; 9930012K11Rik (https://www.ncbi.nlm.nih.gov/gene/268759); http://www.informatics.jax.org/marker/MGI:2145726; https://www.ncbi.nlm.nih.gov/UniGene/clust.cgi?ORG=Mm&CID=12148
				EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = '9930012K11Rik'
			} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q8BGA7') {		# http://www.uniprot.org/uniprot/Q8BGA7; E130309D02Rik (https://www.ncbi.nlm.nih.gov/gene/231868); http://www.informatics.jax.org/marker/MGI:2442621; https://www.ncbi.nlm.nih.gov/UniGene/clust.cgi?ORG=Mm&CID=440105
				EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = 'E130309D02Rik'
			} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q8BWQ6') {		# http://www.uniprot.org/uniprot/Q8BWQ6; 9030624J02Rik (https://www.ncbi.nlm.nih.gov/gene/71517); http://www.informatics.jax.org/marker/MGI:1918767; https://www.ncbi.nlm.nih.gov/UniGene/clust.cgi?ORG=Mm&CID=213298
				EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = '9030624J02Rik'
			} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q8C3W1') {		# http://www.uniprot.org/uniprot/Q8C3W1; 2310022B05Rik (https://www.ncbi.nlm.nih.gov/gene/69551); http://www.informatics.jax.org/marker/MGI:1916801; https://www.ncbi.nlm.nih.gov/UniGene/clust.cgi?ORG=Mm&CID=261920
				EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = '2310022B05Rik'
			} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q8C708') {		# http://www.uniprot.org/uniprot/Q8C708; AI467606 (https://www.ncbi.nlm.nih.gov/gene/101602); http://www.informatics.jax.org/marker/MGI:2141979; https://www.ncbi.nlm.nih.gov/UniGene/clust.cgi?ORG=Mm&CID=490357
				EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = 'AI467606'
			} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q91V76') {		# http://www.uniprot.org/uniprot/Q91V76; 4931406C07Rik (https://www.ncbi.nlm.nih.gov/gene/70984); http://www.informatics.jax.org/marker/MGI:1918234; https://www.ncbi.nlm.nih.gov/UniGene/clust.cgi?ORG=Mm&CID=261842
				EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = '4931406C07Rik'
			} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q9CPS8') {		# http://www.uniprot.org/uniprot/Q9CPS8; 1700019D03Rik (https://www.ncbi.nlm.nih.gov/gene/67080); http://www.informatics.jax.org/marker/MGI:1914330; https://www.ncbi.nlm.nih.gov/UniGene/clust.cgi?ORG=Mm&CID=236287
				EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = '1700019D03Rik'
			} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q9CQT9') {		# http://www.uniprot.org/uniprot/Q9CQT9; 1110008F13Rik (https://www.ncbi.nlm.nih.gov/gene/67388); http://www.informatics.jax.org/marker/MGI:1914638; https://www.ncbi.nlm.nih.gov/UniGene/clust.cgi?ORG=Mm&CID=29264
				EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = '1110008F13Rik'
			} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q9CR34') {		# http://www.uniprot.org/uniprot/Q9CR34; 4921536K21Rik (https://www.ncbi.nlm.nih.gov/gene/67430); http://www.informatics.jax.org/marker/MGI:1914680; https://www.ncbi.nlm.nih.gov/UniGene/clust.cgi?ORG=Mm&CID=347950
				EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = '4921536K21Rik'
			} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q9CXL3') {		# http://www.uniprot.org/uniprot/Q9CXL3; 3110082I17Rik (https://www.ncbi.nlm.nih.gov/gene/73212); http://www.informatics.jax.org/marker/MGI:1920462; https://www.ncbi.nlm.nih.gov/UniGene/clust.cgi?ORG=Mm&CID=84294
				EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = '3110082I17Rik'
			} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q9CYI0') {		# http://www.uniprot.org/uniprot/Q9CYI0; 5730455P16Rik (https://www.ncbi.nlm.nih.gov/gene/70591); http://www.informatics.jax.org/marker/MGI:1917841; https://www.ncbi.nlm.nih.gov/UniGene/clust.cgi?ORG=Mm&CID=86860
				EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = '5730455P16Rik'
			} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q9CYS6') {		# http://www.uniprot.org/uniprot/Q9CYS6; 2810459M11Rik (https://www.ncbi.nlm.nih.gov/gene/72792); http://www.informatics.jax.org/marker/MGI:1920042; https://www.ncbi.nlm.nih.gov/UniGene/clust.cgi?ORG=Mm&CID=54981
				EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = '2810459M11Rik'
			} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q9D2Q3') {		# http://www.uniprot.org/uniprot/Q9D2Q3; 2310057M21Rik (https://www.ncbi.nlm.nih.gov/gene/68277); http://www.informatics.jax.org/marker/MGI:1915527; https://www.ncbi.nlm.nih.gov/UniGene/clust.cgi?ORG=Mm&CID=121122
				EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = '2310057M21Rik'
			} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q9D7E4') {		# http://www.uniprot.org/uniprot/Q9D7E4; 2310011J03Rik (https://www.ncbi.nlm.nih.gov/gene/66374); http://www.informatics.jax.org/marker/MGI:1913624; https://www.ncbi.nlm.nih.gov/UniGene/clust.cgi?ORG=Mm&CID=383824
				EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = '2310011J03Rik'
			} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q9D937') {		# http://www.uniprot.org/uniprot/Q9D937; 1810009A15Rik (https://www.ncbi.nlm.nih.gov/gene/66276); http://www.informatics.jax.org/marker/MGI:1913526; https://www.ncbi.nlm.nih.gov/UniGene/clust.cgi?ORG=Mm&CID=27503
				EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = '1810009A15Rik'
			} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q8BH86') {		# http://www.uniprot.org/uniprot/Q8BH86; Dglucy (https://www.ncbi.nlm.nih.gov/gene/217830); http://www.informatics.jax.org/marker/MGI:2444813; 9030617O03Rik (https://www.ncbi.nlm.nih.gov/UniGene/clust.cgi?ORG=Mm&CID=218590)
				EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = 'Dglucy 9030617O03Rik'
			} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q9CQE8') {		# http://www.uniprot.org/uniprot/Q9CQE8; Rtraf (https://www.ncbi.nlm.nih.gov/gene/68045); http://www.informatics.jax.org/marker/MGI:1915295; 2700060E02Rik (https://www.ncbi.nlm.nih.gov/UniGene/clust.cgi?ORG=Mm&CID=21932)
				EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = 'Rtraf 2700060E02Rik'
			} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q9CYZ6') {		# http://www.uniprot.org/uniprot/Q9CYZ6; Rex1bd (https://www.ncbi.nlm.nih.gov/gene/66462); http://www.informatics.jax.org/marker/MGI:1913712; 2810428I15Rik (https://www.ncbi.nlm.nih.gov/UniGene/clust.cgi?ORG=Mm&CID=28242)
				EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = 'Rex1bd 2810428I15Rik'
			} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q9D735') {		# http://www.uniprot.org/uniprot/Q9D735; Trir (https://www.ncbi.nlm.nih.gov/gene/68544); http://www.informatics.jax.org/marker/MGI:1922833; 2310036O22Rik (https://www.ncbi.nlm.nih.gov/UniGene/clust.cgi?ORG=Mm&CID=196005)
				EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = 'Trir 2310036O22Rik'
			} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q6PIU9') {		# http://www.uniprot.org/uniprot/Q6PIU9; Gm38495 (https://www.ncbi.nlm.nih.gov/gene/102637099); http://www.informatics.jax.org/marker/MGI:5621380; Aak1 (https://www.ncbi.nlm.nih.gov/UniGene/clust.cgi?ORG=Mm&CID=221038)
				EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = 'Gm38495 Aak1'
			} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q3THJ6') {		# http://www.uniprot.org/uniprot/Q3THJ6; Ribosomal protein L10 (Rpl10) (https://www.ncbi.nlm.nih.gov/UniGene/clust.cgi?ORG=Mm&CID=100113; https://www.ncbi.nlm.nih.gov/UniGene/clust.cgi?ORG=Mm&CID=315052; https://www.ncbi.nlm.nih.gov/UniGene/clust.cgi?ORG=Mm&CID=352201); https://www.ncbi.nlm.nih.gov/UniGene/clust.cgi?ORG=Mm&CID=472707
				EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = 'Rpl10'
			} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'P01878') {		# http://www.uniprot.org/uniprot/P01878; UniGene (Mm.246497.; Mm.303767.; Mm.304472.; Mm.305376.; Mm.313447.; Mm.313465.; Mm.313480.; Mm.313488.; Mm.335723.; Mm.342177.; Mm.342187.; Mm.351752.; Mm.390491.; Mm.390523.; Mm.390541.; Mm.423002.; Mm.423936.; Mm.425757.; Mm.425779.; Mm.427704.; Mm.435577.; Mm.436242.; Mm.436336.; Mm.436337.; Mm.458003.; Mm.475073.)
				EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = ''
			} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'P01630') {		# http://www.uniprot.org/uniprot/P01630; Immunoglobulin kappa chain complex (Igk) (https://www.ncbi.nlm.nih.gov/UniGene/clust.cgi?ORG=Mm&CID=333124)
				EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = 'Igk'
			} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'P01631') {		# http://www.uniprot.org/uniprot/P01631; Immunoglobulin kappa chain complex (Igk) (https://www.ncbi.nlm.nih.gov/UniGene/clust.cgi?ORG=Mm&CID=333124); Immunoglobulin kappa variable 1-117 (Igkv1-117) (https://www.ncbi.nlm.nih.gov/UniGene/clust.cgi?ORG=Mm&CID=489779)
				EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = 'Igk Igkv1-117'
			} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'P01851') {		# http://www.uniprot.org/uniprot/P01851; T cell receptor beta, joining region (Tcrb-J) (https://www.ncbi.nlm.nih.gov/UniGene/clust.cgi?ORG=Mm&CID=333026)
				EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = 'Tcrb-J'
			} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'P01915') {		# http://www.uniprot.org/uniprot/P01915; Histocompatibility 2, class II antigen E beta (H2-Eb1) (https://www.ncbi.nlm.nih.gov/UniGene/clust.cgi?ORG=Mm&CID=22564)
				EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = 'H2-Eb1'
			} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'P03987') {		# http://www.uniprot.org/uniprot/P03987; Immunoglobulin heavy constant mu (Ighm) (https://www.ncbi.nlm.nih.gov/UniGene/clust.cgi?ORG=Mm&CID=342177)
				EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = 'Ighm'
			} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q61649') {		# http://www.uniprot.org/uniprot/Q61649; Hemoglobin alpha, adult chain 2 (Hba-a2) (https://www.ncbi.nlm.nih.gov/UniGene/clust.cgi?ORG=Mm&CID=196110)
				EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = 'Hba-a2'
			} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q8BZ90') {		# http://www.uniprot.org/uniprot/Q8BZ90; Family with sequence similarity 13, member C (Fam13c) (https://www.ncbi.nlm.nih.gov/UniGene/clust.cgi?ORG=Mm&CID=477813)
				EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = 'Fam13c'
			} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q8C847') {		# http://www.uniprot.org/uniprot/Q8C847; Galactosidase, beta 1 (Glb1) (https://www.ncbi.nlm.nih.gov/UniGene/clust.cgi?ORG=Mm&CID=290516)
				EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = 'Glb1'
			} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q9Z1Z7') {		# http://www.uniprot.org/uniprot/Q9Z1Z7; Fibronectin 1 (Fn1) (https://www.ncbi.nlm.nih.gov/UniGene/clust.cgi?ORG=Mm&CID=193099)
				EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = 'Fn1'
			} else {
				if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'P01635') {			# http://www.uniprot.org/uniprot/P01635; Igkv12-41 (https://www.ncbi.nlm.nih.gov/gene/619960); http://www.informatics.jax.org/marker/MGI:4439772
					# Gene: N/A; # Protein: Ig kappa chain V-V region K2
					EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = 'Igkv12-41'
				} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'P01645') {	# http://www.uniprot.org/uniprot/P01645
					# Gene: N/A; # Protein: Ig kappa chain V-V region HP 93G7
				} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'P01646') {	# http://www.uniprot.org/uniprot/P01646
					# Gene: N/A; # Protein: Ig kappa chain V-V region HP 123E6
				} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'P01660') {	# http://www.uniprot.org/uniprot/P01660
					# Gene: N/A; # Protein: Ig kappa chain V-III region PC 3741/TEPC 111
				} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'P01663') {	# http://www.uniprot.org/uniprot/P01663
					# Gene: N/A; # Protein: Ig kappa chain V-III region PC 4050
				} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'P01668') {	# http://www.uniprot.org/uniprot/P01668
					# Gene: N/A; # Protein: Ig kappa chain V-III region PC 7210
				} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'P01670') {	# http://www.uniprot.org/uniprot/P01670
					# Gene: N/A; # Protein: Ig kappa chain V-III region PC 6684
				} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'P01750') {	# http://www.uniprot.org/uniprot/P01750
					# Gene: N/A; # Protein: Ig heavy chain V region 102
				} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'P01786') {	# http://www.uniprot.org/uniprot/P01786
					# Gene: N/A; # Protein: Ig heavy chain V region MOPC 47A
				} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'P01794') {	# http://www.uniprot.org/uniprot/P01794
					# Gene: N/A; # Protein: Ig heavy chain V region HPCG14
				} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'P06330') {	# http://www.uniprot.org/uniprot/P06330
					# Gene: N/A; # Protein: Ig heavy chain V region AC38 205.12
				} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'P10404') {	# http://www.uniprot.org/uniprot/P10404
					# Gene: N/A; # Protein: MLV-related proviral Env polyprotein
				} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'A0N8M6') {	# http://www.uniprot.org/uniprot/A0N8M6
					# Gene: N/A; # Protein: Submitted name: H4a-3 coding region
				} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'A2N1N0') {	# http://www.uniprot.org/uniprot/A2N1N0
					# Gene: N/A; # Protein: Submitted name: B cell antigen receptor
				} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'B8YPG2') {	# http://www.uniprot.org/uniprot/B8YPG2
					# Gene: N/A; # Protein: Submitted name: Envelope glycoprotein 52
				} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q1KYM2') {	# http://www.uniprot.org/uniprot/Q1KYM2
					# Gene: N/A; # Protein: Submitted name: Gag-pro-pol polyprotein
				} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q61915') {	# http://www.uniprot.org/uniprot/Q61915
					# Gene: N/A; # Protein: Submitted name: Mtv-9 protein
				} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q61975') {	# http://www.uniprot.org/uniprot/Q61975
					# Gene: N/A; # Protein: Submitted name: Uncharacterized protein
				} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q6YIY0') {	# http://www.uniprot.org/uniprot/Q6YIY0
					# Gene: N/A; # Protein: Submitted name: Gag protein
				} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q8C9G6') {	# http://www.uniprot.org/uniprot/Q8C9G6
					# Gene: N/A; # Protein: Submitted name: Uncharacterized protein
				} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q9D2M4') {	# http://www.uniprot.org/uniprot/Q9D2M4
					# Gene: N/A; # Protein: Submitted name: Uncharacterized protein
				} else {
					cat(paste('[Warning] Please check Gene symbol of UniProt accession "',EVpedia_Top100.Protein_Mouse.UniProt_accession[i],'". Gene symbol might be wrong.\n',sep=""))
				}
				EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = ''
			}
			EVpedia_Top100.Protein_Mouse.Gene_symbol[i] = EVpedia_Top100.Protein_Mouse[i,'Gene.symbol']
		} else {
			EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = ''
			EVpedia_Top100.Protein_Mouse.Gene_symbol[i] = EVpedia_Top100.Protein_Mouse[i,'Gene.symbol']
		}
	}
}

# --- add in missed Gene name(s) ---
for(i in 1:nrow(EVpedia_Top100.Protein_Mouse))
{
	if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'E9PZY8') {				# http://www.uniprot.org/uniprot/E9PZY8; Virma (https://www.ncbi.nlm.nih.gov/gene/66185)
		EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = paste(trimws(EVpedia_Top100.Protein_Mouse[i,'Gene.symbol']),' Virma',sep="")
	} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'P17182') {		# http://www.uniprot.org/uniprot/P17182; Eno1b (https://www.ncbi.nlm.nih.gov/gene/433182)
		EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = paste(trimws(EVpedia_Top100.Protein_Mouse[i,'Gene.symbol']),' Eno1b',sep="")
	} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'P52480') {		# http://www.uniprot.org/uniprot/P52480; Pkm (https://www.ncbi.nlm.nih.gov/gene/18746)
		EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = paste(trimws(EVpedia_Top100.Protein_Mouse[i,'Gene.symbol']),' Pkm',sep="")
	} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'P62821') {		# http://www.uniprot.org/uniprot/P62821; Rab1a (https://www.ncbi.nlm.nih.gov/gene/19324)
		EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = paste(trimws(EVpedia_Top100.Protein_Mouse[i,'Gene.symbol']),' Rab1a',sep="")
	} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'P68040') {		# http://www.uniprot.org/uniprot/P68040; Rack1 (https://www.ncbi.nlm.nih.gov/gene/14694)
		EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = paste(trimws(EVpedia_Top100.Protein_Mouse[i,'Gene.symbol']),' Rack1',sep="")
	} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'P97351') {		# http://www.uniprot.org/uniprot/P97351; Rps3a1 (https://www.ncbi.nlm.nih.gov/gene/20091)
		EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = paste(trimws(EVpedia_Top100.Protein_Mouse[i,'Gene.symbol']),' Rps3a1',sep="")
	} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'P97825') {		# http://www.uniprot.org/uniprot/P97825; Jpt1 (https://www.ncbi.nlm.nih.gov/gene/15374)
		EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = paste(trimws(EVpedia_Top100.Protein_Mouse[i,'Gene.symbol']),' Jpt1',sep="")
	} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q3UMU9') {		# http://www.uniprot.org/uniprot/Q3UMU9; Hdgfl2 (https://www.ncbi.nlm.nih.gov/gene/15193)
		EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = paste(trimws(EVpedia_Top100.Protein_Mouse[i,'Gene.symbol']),' Hdgfl2',sep="")
	} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q61549') {		# http://www.uniprot.org/uniprot/Q61549; Adgre1 (https://www.ncbi.nlm.nih.gov/gene/13733)
		EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = paste(trimws(EVpedia_Top100.Protein_Mouse[i,'Gene.symbol']),' Adgre1',sep="")
	} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q64281') {		# http://www.uniprot.org/uniprot/Q64281; Lilrb4a (https://www.ncbi.nlm.nih.gov/gene/14728)
		EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = paste(trimws(EVpedia_Top100.Protein_Mouse[i,'Gene.symbol']),' Lilrb4a',sep="")
	} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q6A0A9') {		# http://www.uniprot.org/uniprot/Q6A0A9; Fam120a (https://www.ncbi.nlm.nih.gov/gene/218236)
		EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = paste(trimws(EVpedia_Top100.Protein_Mouse[i,'Gene.symbol']),' Fam120a',sep="")
	} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q6XQH0') {		# http://www.uniprot.org/uniprot/Q6XQH0; Gal3st2b (https://www.ncbi.nlm.nih.gov/gene/100041596)
		EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = paste(trimws(EVpedia_Top100.Protein_Mouse[i,'Gene.symbol']),' Gal3st2b',sep="")
	} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q7TSY8') {		# http://www.uniprot.org/uniprot/Q7TSY8; Sgo2a (https://www.ncbi.nlm.nih.gov/gene/68549)
		EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = paste(trimws(EVpedia_Top100.Protein_Mouse[i,'Gene.symbol']),' Sgo2a',sep="")
	} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q8BFQ8') {		# http://www.uniprot.org/uniprot/Q8BFQ8; Gatd1 (https://www.ncbi.nlm.nih.gov/gene/213350)
		EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = paste(trimws(EVpedia_Top100.Protein_Mouse[i,'Gene.symbol']),' Gatd1',sep="")
	} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q8BUI3') {		# http://www.uniprot.org/uniprot/Q8BUI3; Lrwd1 (https://www.ncbi.nlm.nih.gov/gene/71735)
		EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = paste(trimws(EVpedia_Top100.Protein_Mouse[i,'Gene.symbol']),' Lrwd1',sep="")
	} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q8CFG8') {		# http://www.uniprot.org/uniprot/Q8CFG8; C1s2 (https://www.ncbi.nlm.nih.gov/gene/317677)
		EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = paste(trimws(EVpedia_Top100.Protein_Mouse[i,'Gene.symbol']),' C1s2',sep="")
	} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q8CHH5') {		# http://www.uniprot.org/uniprot/Q8CHH5; Bicral (https://www.ncbi.nlm.nih.gov/gene/210982)
		EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = paste(trimws(EVpedia_Top100.Protein_Mouse[i,'Gene.symbol']),' Bicral',sep="")
	} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q8K3W0') {		# http://www.uniprot.org/uniprot/Q8K3W0; Babam2 (https://www.ncbi.nlm.nih.gov/gene/107976)
		EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = paste(trimws(EVpedia_Top100.Protein_Mouse[i,'Gene.symbol']),' Babam2',sep="")
	} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q9CV28') {		# http://www.uniprot.org/uniprot/Q9CV28; Mindy3 (https://www.ncbi.nlm.nih.gov/gene/66960)
		EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = paste(trimws(EVpedia_Top100.Protein_Mouse[i,'Gene.symbol']),' Mindy3',sep="")
	} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q9D1N4') {		# http://www.uniprot.org/uniprot/Q9D1N4; Mymk (https://www.ncbi.nlm.nih.gov/gene/66139)
		EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = paste(trimws(EVpedia_Top100.Protein_Mouse[i,'Gene.symbol']),' Mymk',sep="")
	} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q9ERE8') {		# http://www.uniprot.org/uniprot/Q9ERE8; Tlnrd1 (https://www.ncbi.nlm.nih.gov/gene/80889)
		EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = paste(trimws(EVpedia_Top100.Protein_Mouse[i,'Gene.symbol']),' Tlnrd1',sep="")
	} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q9ERY9') {		# http://www.uniprot.org/uniprot/Q9ERY9; Erg28 (https://www.ncbi.nlm.nih.gov/gene/58520)
		EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = paste(trimws(EVpedia_Top100.Protein_Mouse[i,'Gene.symbol']),' Erg28',sep="")
	} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q9R112') {		# http://www.uniprot.org/uniprot/Q9R112; Sqor (https://www.ncbi.nlm.nih.gov/gene/59010)
		EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = paste(trimws(EVpedia_Top100.Protein_Mouse[i,'Gene.symbol']),' Sqor',sep="")
	} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q9Z1R4') {		# http://www.uniprot.org/uniprot/Q9Z1R4; D17H6S53E (https://www.ncbi.nlm.nih.gov/gene/114585)
		EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = paste(trimws(EVpedia_Top100.Protein_Mouse[i,'Gene.symbol']),' D17H6S53E',sep="")
	} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q9Z222') {		# http://www.uniprot.org/uniprot/Q9Z222; B3gnt2 (https://www.ncbi.nlm.nih.gov/gene/53625)
		EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = paste(trimws(EVpedia_Top100.Protein_Mouse[i,'Gene.symbol']),' B3gnt2',sep="")
	} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q8QZV7') {		# http://www.uniprot.org/uniprot/Q8QZV7; Ints13 (https://www.ncbi.nlm.nih.gov/gene/71177)
		EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = paste(trimws(EVpedia_Top100.Protein_Mouse[i,'Gene.symbol']),' Ints13',sep="")
	} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q8R3P6') {		# http://www.uniprot.org/uniprot/Q8R3P6; Ints14 (https://www.ncbi.nlm.nih.gov/gene/69882)
		EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = paste(trimws(EVpedia_Top100.Protein_Mouse[i,'Gene.symbol']),' Ints14',sep="")
	} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'O54724') {		# http://www.uniprot.org/uniprot/O54724; Cavin1 (https://www.ncbi.nlm.nih.gov/gene/19285)
		EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = paste(trimws(EVpedia_Top100.Protein_Mouse[i,'Gene.symbol']),' Cavin1',sep="")
	} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q63918') {		# http://www.uniprot.org/uniprot/Q63918; Cavin2 (https://www.ncbi.nlm.nih.gov/gene/20324)
		EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = paste(trimws(EVpedia_Top100.Protein_Mouse[i,'Gene.symbol']),' Cavin2',sep="")
	} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q8VE91') {		# http://www.uniprot.org/uniprot/Q8VE91; Retreg1 (https://www.ncbi.nlm.nih.gov/gene/66270)
		EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = paste(trimws(EVpedia_Top100.Protein_Mouse[i,'Gene.symbol']),' Retreg1',sep="")
	} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q6NS82') {		# http://www.uniprot.org/uniprot/Q6NS82; Retreg2 (https://www.ncbi.nlm.nih.gov/gene/227298)
		EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = paste(trimws(EVpedia_Top100.Protein_Mouse[i,'Gene.symbol']),' Retreg2',sep="")
	} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q9CQV4') {		# http://www.uniprot.org/uniprot/Q9CQV4; Retreg3 (https://www.ncbi.nlm.nih.gov/gene/67998)
		EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = paste(trimws(EVpedia_Top100.Protein_Mouse[i,'Gene.symbol']),' Retreg3',sep="")
	} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'E9Q1Z0') {		# http://www.uniprot.org/uniprot/E9Q1Z0; Krt90 (https://www.ncbi.nlm.nih.gov/gene/239673)		# http://www.uniprot.org/uniprot/E9Q1Z0 -> Not found by "biomaRt", but could be found by "BridgeDbR" R package.
		EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = paste(trimws(EVpedia_Top100.Protein_Mouse[i,'Gene.symbol']),' Krt90',sep="")
	} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'P50247') {		# http://www.uniprot.org/uniprot/P50247; Gm4737 (https://www.ncbi.nlm.nih.gov/gene/11615)
		EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = paste(trimws(EVpedia_Top100.Protein_Mouse[i,'Gene.symbol']),' Gm4737',sep="")
	} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q80SU7') {		# http://www.uniprot.org/uniprot/Q80SU7; Gm4070 (https://www.ncbi.nlm.nih.gov/gene/100042856)
		EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = paste(trimws(EVpedia_Top100.Protein_Mouse[i,'Gene.symbol']),' Gm4070',sep="")
	} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q69ZQ1') {		# http://www.uniprot.org/uniprot/Q69ZQ1; AI464131 (https://www.ncbi.nlm.nih.gov/gene/329828)
		EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = paste(trimws(EVpedia_Top100.Protein_Mouse[i,'Gene.symbol']),' AI464131',sep="")
	} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q6PDI5') {		# http://www.uniprot.org/uniprot/Q6PDI5; AI314180 (https://www.ncbi.nlm.nih.gov/gene/230249); http://www.informatics.jax.org/marker/MGI:2140220
		EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = paste(trimws(EVpedia_Top100.Protein_Mouse[i,'Gene.symbol']),' AI314180',sep="")
	} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q4PJX1') {		# http://www.uniprot.org/uniprot/Q4PJX1; BC003331 (https://www.ncbi.nlm.nih.gov/gene/226499)
		EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = paste(trimws(EVpedia_Top100.Protein_Mouse[i,'Gene.symbol']),' BC003331',sep="")
	} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'A2AAE1') {		# http://www.uniprot.org/uniprot/A2AAE1; 4932438A13Rik (https://www.ncbi.nlm.nih.gov/gene/229227)
		EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = paste(trimws(EVpedia_Top100.Protein_Mouse[i,'Gene.symbol']),' 4932438A13Rik',sep="")
	} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'P56379') {		# http://www.uniprot.org/uniprot/P56379; 2010107E04Rik (https://www.ncbi.nlm.nih.gov/gene/70257)
		EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = paste(trimws(EVpedia_Top100.Protein_Mouse[i,'Gene.symbol']),' 2010107E04Rik',sep="")
	} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q8C6C9') {		# http://www.uniprot.org/uniprot/Q8C6C9; 2310057J18Rik (https://www.ncbi.nlm.nih.gov/gene/67719)
		EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = paste(trimws(EVpedia_Top100.Protein_Mouse[i,'Gene.symbol']),' 2310057J18Rik',sep="")
	} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q8CB14') {		# http://www.uniprot.org/uniprot/Q8CB14; 4930503L19Rik (https://www.ncbi.nlm.nih.gov/gene/269033)
		EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = paste(trimws(EVpedia_Top100.Protein_Mouse[i,'Gene.symbol']),' 4930503L19Rik',sep="")
	} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q8R404') {		# http://www.uniprot.org/uniprot/Q8R404; 2410015M20Rik (https://www.ncbi.nlm.nih.gov/gene/224904)
		EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = paste(trimws(EVpedia_Top100.Protein_Mouse[i,'Gene.symbol']),' 2410015M20Rik',sep="")
	} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q8VC42') {		# http://www.uniprot.org/uniprot/Q8VC42; 3110002H16Rik (https://www.ncbi.nlm.nih.gov/gene/76482)
		EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = paste(trimws(EVpedia_Top100.Protein_Mouse[i,'Gene.symbol']),' 3110002H16Rik',sep="")
	} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q91X21') {		# http://www.uniprot.org/uniprot/Q91X21; 2510039O18Rik (https://www.ncbi.nlm.nih.gov/gene/77034)
		EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = paste(trimws(EVpedia_Top100.Protein_Mouse[i,'Gene.symbol']),' 2510039O18Rik',sep="")
	} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q99LX5') {		# http://www.uniprot.org/uniprot/Q99LX5; 2310033P09Rik (https://www.ncbi.nlm.nih.gov/gene/67862)
		EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = paste(trimws(EVpedia_Top100.Protein_Mouse[i,'Gene.symbol']),' 2310033P09Rik',sep="")
	} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q9CQD4') {		# http://www.uniprot.org/uniprot/Q9CQD4; 2610002M06Rik (https://www.ncbi.nlm.nih.gov/gene/67028)
		EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = paste(trimws(EVpedia_Top100.Protein_Mouse[i,'Gene.symbol']),' 2610002M06Rik',sep="")
	} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'P51807') {		# http://www.uniprot.org/uniprot/P51807; Dynlt1a (https://www.ncbi.nlm.nih.gov/gene/100310872); Dynlt1b (https://www.ncbi.nlm.nih.gov/gene/21648); Dynlt1c (https://www.ncbi.nlm.nih.gov/gene/100040563); Dynlt1f (https://www.ncbi.nlm.nih.gov/gene/100040531)
		EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = paste(trimws(EVpedia_Top100.Protein_Mouse[i,'Gene.symbol']),' Dynlt1a Dynlt1b Dynlt1c Dynlt1f',sep="")
	} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'Q8BR63') {		# http://www.uniprot.org/uniprot/Q8BR63; Fam177a (https://www.ncbi.nlm.nih.gov/gene/73385); 1700047I17Rik2 (https://www.ncbi.nlm.nih.gov/gene/100101807)
		EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = paste(trimws(EVpedia_Top100.Protein_Mouse[i,'Gene.symbol']),' Fam177a 1700047I17Rik2',sep="")
	} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'P10853') {		# http://www.uniprot.org/uniprot/P10853; Hist1h2bq (https://www.ncbi.nlm.nih.gov/gene/665596), Hist1h2br (https://www.ncbi.nlm.nih.gov/gene/665622)
		EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = paste(trimws(EVpedia_Top100.Protein_Mouse[i,'Gene.symbol']),'; Hist1h2bq','; Hist1h2br',sep="")
	} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'P22752') {		# http://www.uniprot.org/uniprot/P22752; Hist1h2ap (https://www.ncbi.nlm.nih.gov/gene/319171)
		EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = paste(trimws(EVpedia_Top100.Protein_Mouse[i,'Gene.symbol']),'; Hist1h2ap',sep="")
	} else if(EVpedia_Top100.Protein_Mouse.UniProt_accession[i] == 'P62806') {		# http://www.uniprot.org/uniprot/P62806; Hist1h4n (https://www.ncbi.nlm.nih.gov/gene/319161)
		EVpedia_Top100.Protein_Mouse[i,'Gene.symbol'] = paste(trimws(EVpedia_Top100.Protein_Mouse[i,'Gene.symbol']),'; Hist1h4n',sep="")
	}
	
	EVpedia_Top100.Protein_Mouse.Gene_symbol[i] = EVpedia_Top100.Protein_Mouse[i,'Gene.symbol']
}

EVpedia_Top100.Protein_Mouse.Index <- trimws(EVpedia_Top100.Protein_Mouse[,'Index'])
EVpedia_Top100.Protein_Mouse.UniProt_accession <- trimws(EVpedia_Top100.Protein_Mouse[,'UniProt.accession'])
EVpedia_Top100.Protein_Mouse.UniProt_name <- trimws(EVpedia_Top100.Protein_Mouse[,'UniProt.name'])
EVpedia_Top100.Protein_Mouse.Status <- trimws(EVpedia_Top100.Protein_Mouse[,'Status'])
EVpedia_Top100.Protein_Mouse.Protein_name <- trimws(EVpedia_Top100.Protein_Mouse[,'Protein.name'])
EVpedia_Top100.Protein_Mouse.Gene_symbol <- trimws(EVpedia_Top100.Protein_Mouse[,'Gene.symbol'])
EVpedia_Top100.Protein_Mouse.Species <- trimws(EVpedia_Top100.Protein_Mouse[,'Species'])
EVpedia_Top100.Protein_Mouse.Identification_number <- trimws(EVpedia_Top100.Protein_Mouse[,'Identification.number'])


# === id mapping ================================================================================================================
ensembl = ensembl.mmu
#listFilters(ensembl)[,c('name','description')]
#listAttributes(ensembl)[,c('name','description')]

# --- search for IDs by uniprotswissprot ---
a = EVpedia_Top100.Protein_Mouse.UniProt_accession
b = ensembl.mmu

EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes <- getBM(attributes=c('entrezgene','mgi_id','mgi_symbol'), filters='uniprotswissprot', values=a, mart=b)
EVpedia_Top100.Protein_Mouse.uniprotswissprot.wikigene_description <- getBM(attributes=c('mgi_symbol','wikigene_description'), filters='uniprotswissprot', values=a, mart=b)
EVpedia_Top100.Protein_Mouse.uniprotswissprot.wikigene_description.table = merge(EVpedia_Top100.Protein_Mouse.uniprotswissprot.wikigene_description, EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes, by.x="mgi_symbol", by.y="mgi_symbol")
EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes = EVpedia_Top100.Protein_Mouse.uniprotswissprot.wikigene_description.table[,c('entrezgene','mgi_id','mgi_symbol','wikigene_description'),drop=FALSE]
colnames(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes)[which(colnames(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes)=='wikigene_description')] = 'description'
EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes = unique(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes)

EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:95640'),'entrezgene'] = c(14433)		# http://www.informatics.jax.org/marker/MGI:95640; https://www.ncbi.nlm.nih.gov/gene/14433
EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:95655'),'entrezgene'] = c(14451)		# http://www.informatics.jax.org/marker/MGI:95655; https://www.ncbi.nlm.nih.gov/gene/14451
EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:96604'),'entrezgene'] = c(16402)		# http://www.informatics.jax.org/marker/MGI:96604; https://www.ncbi.nlm.nih.gov/gene/16402
EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:96904'),'entrezgene'] = c(17113)		# http://www.informatics.jax.org/marker/MGI:96904; https://www.ncbi.nlm.nih.gov/gene/17113
EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:97298'),'entrezgene'] = c(18000)		# http://www.informatics.jax.org/marker/MGI:97298; https://www.ncbi.nlm.nih.gov/gene/18000
EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:97562'),'entrezgene'] = c(668435)		# http://www.informatics.jax.org/marker/MGI:97562; https://www.ncbi.nlm.nih.gov/gene/668435
EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:97843'),'entrezgene'] = c(19339)		# http://www.informatics.jax.org/marker/MGI:97843; https://www.ncbi.nlm.nih.gov/gene/19339
EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:98036'),'entrezgene'] = c(19942)		# http://www.informatics.jax.org/marker/MGI:98036; https://www.ncbi.nlm.nih.gov/gene/19942
EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:98159'),'entrezgene'] = c(20104)		# http://www.informatics.jax.org/marker/MGI:98159; https://www.ncbi.nlm.nih.gov/gene/20104
EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:102701'),'entrezgene'] = c(14728)		# http://www.informatics.jax.org/marker/MGI:102701; https://www.ncbi.nlm.nih.gov/gene/14728
EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:104681'),'entrezgene'] = c(15239)		# http://www.informatics.jax.org/marker/MGI:104681; https://www.ncbi.nlm.nih.gov/gene/15239
EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:107686'),'entrezgene'] = c(17992)		# http://www.informatics.jax.org/marker/MGI:107686; https://www.ncbi.nlm.nih.gov/gene/17992
EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:109279'),'entrezgene'] = c(18115)		# http://www.informatics.jax.org/marker/MGI:109279; https://www.ncbi.nlm.nih.gov/gene/18115
EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:109620'),'entrezgene'] = c(11877)		# http://www.informatics.jax.org/marker/MGI:109620; https://www.ncbi.nlm.nih.gov/gene/11877
EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:894645'),'entrezgene'] = c(18693)		# http://www.informatics.jax.org/marker/MGI:894645; https://www.ncbi.nlm.nih.gov/gene/18693
EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:1343176'),'entrezgene'] = c(21855)		# http://www.informatics.jax.org/marker/MGI:1343176; https://www.ncbi.nlm.nih.gov/gene/21855
EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:1346871'),'entrezgene'] = c(26400)		# http://www.informatics.jax.org/marker/MGI:1346871; https://www.ncbi.nlm.nih.gov/gene/26400
EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:1347049'),'entrezgene'] = c(26372)		# http://www.informatics.jax.org/marker/MGI:1347049; https://www.ncbi.nlm.nih.gov/gene/26372
EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:1349183'),'entrezgene'] = c(18829)		# http://www.informatics.jax.org/marker/MGI:1349183; https://www.ncbi.nlm.nih.gov/gene/18829
EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:1913412'),'entrezgene'] = c(66162)		# http://www.informatics.jax.org/marker/MGI:1913412; https://www.ncbi.nlm.nih.gov/gene/66162
EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:1913712'),'entrezgene'] = c(66462)		# http://www.informatics.jax.org/marker/MGI:1913712; https://www.ncbi.nlm.nih.gov/gene/66462
EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:1914582'),'entrezgene'] = c(67332)		# http://www.informatics.jax.org/marker/MGI:1914582; https://www.ncbi.nlm.nih.gov/gene/67332
EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:1915755'),'entrezgene'] = c(68505)		# http://www.informatics.jax.org/marker/MGI:1915755; https://www.ncbi.nlm.nih.gov/gene/68505
EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:1915896'),'entrezgene'] = c(68646)		# http://www.informatics.jax.org/marker/MGI:1915896; https://www.ncbi.nlm.nih.gov/gene/68646
EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:1915964'),'entrezgene'] = c(102115)	# http://www.informatics.jax.org/marker/MGI:1915964; https://www.ncbi.nlm.nih.gov/gene/102115
EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:1919016'),'entrezgene'] = c(71766)		# http://www.informatics.jax.org/marker/MGI:1919016; https://www.ncbi.nlm.nih.gov/gene/71766
EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:1919094'),'entrezgene'] = c(71844)		# http://www.informatics.jax.org/marker/MGI:1919094; https://www.ncbi.nlm.nih.gov/gene/71844
EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:1919999'),'entrezgene'] = c(72749)		# http://www.informatics.jax.org/marker/MGI:1919999; https://www.ncbi.nlm.nih.gov/gene/72749
EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:1920199'),'entrezgene'] = c(72949)		# http://www.informatics.jax.org/marker/MGI:1920199; https://www.ncbi.nlm.nih.gov/gene/72949
EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:1920994'),'entrezgene'] = c(73744)		# http://www.informatics.jax.org/marker/MGI:1920994; https://www.ncbi.nlm.nih.gov/gene/73744
EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:2137226'),'entrezgene'] = c(75398)		# http://www.informatics.jax.org/marker/MGI:2137226; https://www.ncbi.nlm.nih.gov/gene/75398
EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:2443816'),'entrezgene'] = c(231834)	# http://www.informatics.jax.org/marker/MGI:2443816; https://www.ncbi.nlm.nih.gov/gene/231834
EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:2678374'),'entrezgene'] = c(233073)	# http://www.informatics.jax.org/marker/MGI:2678374; https://www.ncbi.nlm.nih.gov/gene/233073
EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:2684957'),'entrezgene'] = c(227656)	# http://www.informatics.jax.org/marker/MGI:2684957; https://www.ncbi.nlm.nih.gov/gene/227656
EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:2685142'),'entrezgene'] = c(380924)	# http://www.informatics.jax.org/marker/MGI:2685142; https://www.ncbi.nlm.nih.gov/gene/380924
EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:3039607'),'entrezgene'] = c(228812)	# http://www.informatics.jax.org/marker/MGI:3039607; https://www.ncbi.nlm.nih.gov/gene/228812
EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:3641855'),'entrezgene'] = c(100041480)	# http://www.informatics.jax.org/marker/MGI:3641855; https://www.ncbi.nlm.nih.gov/gene/100041480
EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:3643110'),'entrezgene'] = c(666609)	# http://www.informatics.jax.org/marker/MGI:3643110; https://www.ncbi.nlm.nih.gov/gene/666609
EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:3644269'),'entrezgene'] = c(317677)	# http://www.informatics.jax.org/marker/MGI:3644269; https://www.ncbi.nlm.nih.gov/gene/317677
EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:3645078'),'entrezgene'] = c(665989)	# http://www.informatics.jax.org/marker/MGI:3645078; https://www.ncbi.nlm.nih.gov/gene/665989
EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:3646174'),'entrezgene'] = c(621100)	# http://www.informatics.jax.org/marker/MGI:3646174; https://www.ncbi.nlm.nih.gov/gene/621100
EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:3646298'),'entrezgene'] = c(620648)	# http://www.informatics.jax.org/marker/MGI:3646298; https://www.ncbi.nlm.nih.gov/gene/620648
EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:3646907'),'entrezgene'] = c(625193)	# http://www.informatics.jax.org/marker/MGI:3646907; https://www.ncbi.nlm.nih.gov/gene/625193
EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:3704198'),'entrezgene'] = c(100039026)	# http://www.informatics.jax.org/marker/MGI:3704198; https://www.ncbi.nlm.nih.gov/gene/100039026
EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:3704288'),'entrezgene'] = c(100043712)	# http://www.informatics.jax.org/marker/MGI:3704288; https://www.ncbi.nlm.nih.gov/gene/100043712
EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:3704308'),'entrezgene'] = c(100042777)	# http://www.informatics.jax.org/marker/MGI:3704308; https://www.ncbi.nlm.nih.gov/gene/100042777
EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:3704312'),'entrezgene'] = c(100042405)	# http://www.informatics.jax.org/marker/MGI:3704312; https://www.ncbi.nlm.nih.gov/gene/100042405
EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:3704336'),'entrezgene'] = c(100043346)	# http://www.informatics.jax.org/marker/MGI:3704336; https://www.ncbi.nlm.nih.gov/gene/100043346
EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:3704359'),'entrezgene'] = c(100042179)	# http://www.informatics.jax.org/marker/MGI:3704359; https://www.ncbi.nlm.nih.gov/gene/100042179
EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:3704480'),'entrezgene'] = c(100043906)	# http://www.informatics.jax.org/marker/MGI:3704480; https://www.ncbi.nlm.nih.gov/gene/100043906
EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:3704493'),'entrezgene'] = c(672195)	# http://www.informatics.jax.org/marker/MGI:3704493; https://www.ncbi.nlm.nih.gov/gene/672195
EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:3780170'),'entrezgene'] = c(100038991)	# http://www.informatics.jax.org/marker/MGI:3780170; https://www.ncbi.nlm.nih.gov/gene/100038991
EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:3782538'),'entrezgene'] = c(100043313)	# http://www.informatics.jax.org/marker/MGI:3782538; https://www.ncbi.nlm.nih.gov/gene/100043313
EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:4439772'),'entrezgene'] = c(619960)	# http://www.informatics.jax.org/marker/MGI:4439772; https://www.ncbi.nlm.nih.gov/gene/619960
EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:4439824'),'entrezgene'] = c(674018)	# http://www.informatics.jax.org/marker/MGI:4439824; https://www.ncbi.nlm.nih.gov/gene/674018
EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:5012260'),'entrezgene'] = c(100504125)	# http://www.informatics.jax.org/marker/MGI:5012260; https://www.ncbi.nlm.nih.gov/gene/100504125
EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:5306917'),'entrezgene'] = c(100316903)	# http://www.informatics.jax.org/marker/MGI:5306917; https://www.ncbi.nlm.nih.gov/gene/100316903
EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:5313163'),'entrezgene'] = c(105980075)	# http://www.informatics.jax.org/marker/MGI:5313163; https://www.ncbi.nlm.nih.gov/gene/105980075
#EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:5141853'),'entrezgene'] = c()			# http://www.informatics.jax.org/marker/MGI:5141853
#EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:5547771'),'entrezgene'] = c()			# http://www.informatics.jax.org/marker/MGI:5547771
#EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[which(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes[,'mgi_id']=='MGI:5579341'),'entrezgene'] = c()			# http://www.informatics.jax.org/marker/MGI:5579341

EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes = unique(EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes)
EVpedia_Top100.Protein_Mouse.mgi_id <- getBM(attributes=c('uniprotswissprot','mgi_id'), filters='uniprotswissprot', values=a, mart=b)
EVpedia_Top100.Protein_Mouse.mgi_id.table = merge(EVpedia_Top100.Protein_Mouse.mgi_id, EVpedia_Top100.Protein_Mouse.uniprotswissprot.attributes, by.x="mgi_id", by.y="mgi_id")

# http://www.uniprot.org/uniprot/D6RG99; http://www.informatics.jax.org/marker/MGI:1315196; https://www.ncbi.nlm.nih.gov/gene/14356
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='D6RG99'),'entrezgene'] = c(14356)
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='D6RG99'),'mgi_id'] = 'MGI:1315196'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='D6RG99'),'mgi_symbol'] = 'Timm10b'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='D6RG99'),'description'] = 'translocase of inner mitochondrial membrane 10B [Source:MGI Symbol;Acc:MGI:1315196]'
# http://www.uniprot.org/uniprot/O70152; http://www.informatics.jax.org/marker/MGI:1330239; https://www.ncbi.nlm.nih.gov/gene/13480
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='O70152'),'entrezgene'] = c(13480)
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='O70152'),'mgi_id'] = 'MGI:1330239'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='O70152'),'mgi_symbol'] = 'Dpm1'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='O70152'),'description'] = 'dolichol-phosphate (beta-D) mannosyltransferase 1 [Source:MGI Symbol;Acc:MGI:1330239]'
# http://www.uniprot.org/uniprot/P09411; http://www.informatics.jax.org/marker/MGI:97555; https://www.ncbi.nlm.nih.gov/gene/18655
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='P09411'),'entrezgene'] = c(18655)
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='P09411'),'mgi_id'] = 'MGI:97555'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='P09411'),'mgi_symbol'] = 'Pgk1'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='P09411'),'description'] = 'phosphoglycerate kinase 1 [Source:MGI Symbol;Acc:MGI:97555]'
# http://www.uniprot.org/uniprot/P09602; http://www.informatics.jax.org/marker/MGI:96136; https://www.ncbi.nlm.nih.gov/gene/15331
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='P09602'),'entrezgene'] = c(15331)
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='P09602'),'mgi_id'] = 'MGI:96136'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='P09602'),'mgi_symbol'] = 'Hmgn2'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='P09602'),'description'] = 'high mobility group nucleosomal binding domain 2 [Source:MGI Symbol;Acc:MGI:96136]'
# http://www.uniprot.org/uniprot/P17095; http://www.informatics.jax.org/marker/MGI:96160; https://www.ncbi.nlm.nih.gov/gene/15361
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='P17095'),'entrezgene'] = c(15361)
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='P17095'),'mgi_id'] = 'MGI:96160'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='P17095'),'mgi_symbol'] = 'Hmga1'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='P17095'),'description'] = 'high mobility group AT-hook 1 [Source:MGI Symbol;Acc:MGI:96160]'
# http://www.uniprot.org/uniprot/P17182; http://www.informatics.jax.org/marker/MGI:95393; https://www.ncbi.nlm.nih.gov/gene/13806
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='P17182'),'entrezgene'] = c(13806)
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='P17182'),'mgi_id'] = 'MGI:95393'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='P17182'),'mgi_symbol'] = 'Eno1'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='P17182'),'description'] = 'enolase 1, alpha non-neuron [Source:MGI Symbol;Acc:MGI:95393]'
# http://www.uniprot.org/uniprot/P30999; http://www.informatics.jax.org/marker/MGI:105100; https://www.ncbi.nlm.nih.gov/gene/12388
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='P30999'),'entrezgene'] = c(12388)
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='P30999'),'mgi_id'] = 'MGI:105100'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='P30999'),'mgi_symbol'] = 'Ctnnd1'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='P30999'),'description'] = 'catenin (cadherin associated protein), delta 1 [Source:MGI Symbol;Acc:MGI:105100]'
# http://www.uniprot.org/uniprot/P50247; http://www.informatics.jax.org/marker/MGI:87968; https://www.ncbi.nlm.nih.gov/gene/269378
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='P50247'),'entrezgene'] = c(269378)
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='P50247'),'mgi_id'] = 'MGI:87968'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='P50247'),'mgi_symbol'] = 'Ahcy'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='P50247'),'description'] = 'S-adenosylhomocysteine hydrolase [Source:MGI Symbol;Acc:MGI:87968]'
# http://www.uniprot.org/uniprot/P52293; http://www.informatics.jax.org/marker/MGI:103561; https://www.ncbi.nlm.nih.gov/gene/16647
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='P52293'),'entrezgene'] = c(16647)
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='P52293'),'mgi_id'] = 'MGI:103561'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='P52293'),'mgi_symbol'] = 'Kpna2'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='P52293'),'description'] = 'karyopherin (importin) alpha 2 [Source:MGI Symbol;Acc:MGI:103561]'
# http://www.uniprot.org/uniprot/P61358; http://www.informatics.jax.org/marker/MGI:98036; https://www.ncbi.nlm.nih.gov/gene/19942
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='P61358'),'entrezgene'] = c(19942)
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='P61358'),'mgi_id'] = 'MGI:98036'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='P61358'),'mgi_symbol'] = 'Rpl27'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='P61358'),'description'] = 'ribosomal protein L27 [Source:MGI Symbol;Acc:MGI:98036]'
# http://www.uniprot.org/uniprot/P62309; http://www.informatics.jax.org/marker/MGI:1915261; https://www.ncbi.nlm.nih.gov/gene/68011
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='P62309'),'entrezgene'] = c(68011)
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='P62309'),'mgi_id'] = 'MGI:1915261'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='P62309'),'mgi_symbol'] = 'Snrpg'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='P62309'),'description'] = 'small nuclear ribonucleoprotein polypeptide G [Source:MGI Symbol;Acc:MGI:1915261]'
# http://www.uniprot.org/uniprot/P62897; http://www.informatics.jax.org/marker/MGI:88578; https://www.ncbi.nlm.nih.gov/gene/13063
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='P62897'),'entrezgene'] = c(13063)
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='P62897'),'mgi_id'] = 'MGI:88578'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='P62897'),'mgi_symbol'] = 'Cycs'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='P62897'),'description'] = 'cytochrome c, somatic [Source:MGI Symbol;Acc:MGI:88578]'
# http://www.uniprot.org/uniprot/P62984; http://www.informatics.jax.org/marker/MGI:98887; https://www.ncbi.nlm.nih.gov/gene/22186
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='P62984'),'entrezgene'] = c(22186)
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='P62984'),'mgi_id'] = 'MGI:98887'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='P62984'),'mgi_symbol'] = 'Uba52'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='P62984'),'description'] = 'ubiquitin A-52 residue ribosomal protein fusion product 1 [Source:MGI Symbol;Acc:MGI:98887]'
# http://www.uniprot.org/uniprot/P84089; http://www.informatics.jax.org/marker/MGI:108089; https://www.ncbi.nlm.nih.gov/gene/13877
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='P84089'),'entrezgene'] = c(13877)
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='P84089'),'mgi_id'] = 'MGI:108089'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='P84089'),'mgi_symbol'] = 'Erh'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='P84089'),'description'] = 'enhancer of rudimentary homolog (Drosophila) [Source:MGI Symbol;Acc:MGI:108089]'
# http://www.uniprot.org/uniprot/Q60737; http://www.informatics.jax.org/marker/MGI:88543; https://www.ncbi.nlm.nih.gov/gene/12995
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q60737'),'entrezgene'] = c(12995)
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q60737'),'mgi_id'] = 'MGI:88543'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q60737'),'mgi_symbol'] = 'Csnk2a1'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q60737'),'description'] = 'casein kinase 2, alpha 1 polypeptide [Source:MGI Symbol;Acc:MGI:88543]'
# http://www.uniprot.org/uniprot/Q64281; http://www.informatics.jax.org/marker/MGI:102701; https://www.ncbi.nlm.nih.gov/gene/14728
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q64281'),'entrezgene'] = c(14728)
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q64281'),'mgi_id'] = 'MGI:102701'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q64281'),'mgi_symbol'] = 'Lilrb4a'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q64281'),'description'] = 'leukocyte immunoglobulin-like receptor, subfamily B, member 4A [Source:MGI Symbol;Acc:MGI:102701]'
# http://www.uniprot.org/uniprot/Q6PB93; http://www.informatics.jax.org/marker/MGI:894694; https://www.ncbi.nlm.nih.gov/gene/108148
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q6PB93'),'entrezgene'] = c(108148)
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q6PB93'),'mgi_id'] = 'MGI:894694'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q6PB93'),'mgi_symbol'] = 'Galnt2'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q6PB93'),'description'] = 'polypeptide N-acetylgalactosaminyltransferase 2 [Source:MGI Symbol;Acc:MGI:894694]'
# http://www.uniprot.org/uniprot/Q6ZWU9; http://www.informatics.jax.org/marker/MGI:1888676; https://www.ncbi.nlm.nih.gov/gene/57294
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q6ZWU9'),'entrezgene'] = c(57294)
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q6ZWU9'),'mgi_id'] = 'MGI:1888676'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q6ZWU9'),'mgi_symbol'] = 'Rps27'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q6ZWU9'),'description'] = 'ribosomal protein S27 [Source:MGI Symbol;Acc:MGI:1888676]'
# http://www.uniprot.org/uniprot/Q6ZWV3; http://www.informatics.jax.org/marker/MGI:105943; https://www.ncbi.nlm.nih.gov/gene/110954
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q6ZWV3'),'entrezgene'] = c(110954)
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q6ZWV3'),'mgi_id'] = 'MGI:105943'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q6ZWV3'),'mgi_symbol'] = 'Rpl10'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q6ZWV3'),'description'] = 'ribosomal protein L10 [Source:MGI Symbol;Acc:MGI:105943]'
# http://www.uniprot.org/uniprot/Q6ZWV7; http://www.informatics.jax.org/marker/MGI:1913739; https://www.ncbi.nlm.nih.gov/gene/66489
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q6ZWV7'),'entrezgene'] = c(66489)
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q6ZWV7'),'mgi_id'] = 'MGI:1913739'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q6ZWV7'),'mgi_symbol'] = 'Rpl35'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q6ZWV7'),'description'] = 'ribosomal protein L35 [Source:MGI Symbol;Acc:MGI:1913739]'
# http://www.uniprot.org/uniprot/Q6ZWY8; http://www.informatics.jax.org/marker/MGI:109146; https://www.ncbi.nlm.nih.gov/gene/19240
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q6ZWY8'),'entrezgene'] = c(19240)
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q6ZWY8'),'mgi_id'] = 'MGI:109146'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q6ZWY8'),'mgi_symbol'] = 'Tmsb10'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q6ZWY8'),'description'] = 'thymosin, beta 10 [Source:MGI Symbol;Acc:MGI:109146]'
# http://www.uniprot.org/uniprot/Q80WG5; http://www.informatics.jax.org/marker/MGI:2652847; https://www.ncbi.nlm.nih.gov/gene/241296
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q80WG5'),'entrezgene'] = c(241296)
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q80WG5'),'mgi_id'] = 'MGI:2652847'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q80WG5'),'mgi_symbol'] = 'Lrrc8a'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q80WG5'),'description'] = 'leucine rich repeat containing 8A [Source:MGI Symbol;Acc:MGI:2652847]'
# http://www.uniprot.org/uniprot/Q8BGJ9; http://www.informatics.jax.org/marker/MGI:2678374; https://www.ncbi.nlm.nih.gov/gene/233073
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q8BGJ9'),'entrezgene'] = c(233073)
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q8BGJ9'),'mgi_id'] = 'MGI:2678374'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q8BGJ9'),'mgi_symbol'] = 'U2af1l4'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q8BGJ9'),'description'] = 'U2 small nuclear RNA auxiliary factor 1-like 4 [Source:MGI Symbol;Acc:MGI:2678374]'
# http://www.uniprot.org/uniprot/Q8C854; http://www.informatics.jax.org/marker/MGI:104592; https://www.ncbi.nlm.nih.gov/gene/17876
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q8C854'),'entrezgene'] = c(17876)
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q8C854'),'mgi_id'] = 'MGI:104592'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q8C854'),'mgi_symbol'] = 'Myef2'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q8C854'),'description'] = 'myelin basic protein expression factor 2, repressor [Source:MGI Symbol;Acc:MGI:104592]'
# http://www.uniprot.org/uniprot/Q8CE90; http://www.informatics.jax.org/marker/MGI:1346871; https://www.ncbi.nlm.nih.gov/gene/26400
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q8CE90'),'entrezgene'] = c(26400)
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q8CE90'),'mgi_id'] = 'MGI:1346871'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q8CE90'),'mgi_symbol'] = 'Map2k7'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q8CE90'),'description'] = 'mitogen-activated protein kinase kinase 7 [Source:MGI Symbol;Acc:MGI:1346871]'
# http://www.uniprot.org/uniprot/Q8R4R6; http://www.informatics.jax.org/marker/MGI:1916732; https://www.ncbi.nlm.nih.gov/gene/69482
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q8R4R6'),'entrezgene'] = c(69482)
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q8R4R6'),'mgi_id'] = 'MGI:1916732'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q8R4R6'),'mgi_symbol'] = 'Nup35'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q8R4R6'),'description'] = 'nucleoporin 35 [Source:MGI Symbol;Acc:MGI:1916732]'
# http://www.uniprot.org/uniprot/Q921G6; http://www.informatics.jax.org/marker/MGI:1917193; https://www.ncbi.nlm.nih.gov/gene/231798
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q921G6'),'entrezgene'] = c(231798)
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q921G6'),'mgi_id'] = 'MGI:1917193'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q921G6'),'mgi_symbol'] = 'Lrch4'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q921G6'),'description'] = 'leucine-rich repeats and calponin homology (CH) domain containing 4 [Source:MGI Symbol;Acc:MGI:1917193]'
# http://www.uniprot.org/uniprot/Q99LN9; http://www.informatics.jax.org/marker/MGI:1915964; https://www.ncbi.nlm.nih.gov/gene/102115
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q99LN9'),'entrezgene'] = c(102115)
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q99LN9'),'mgi_id'] = 'MGI:1915964'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q99LN9'),'mgi_symbol'] = 'Dohh'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q99LN9'),'description'] = 'deoxyhypusine hydroxylase/monooxygenase [Source:MGI Symbol;Acc:MGI:1915964]'
# http://www.uniprot.org/uniprot/Q9CQH7; http://www.informatics.jax.org/marker/MGI:1915312; https://www.ncbi.nlm.nih.gov/gene/70533
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q9CQH7'),'entrezgene'] = c(70533)
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q9CQH7'),'mgi_id'] = 'MGI:1915312'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q9CQH7'),'mgi_symbol'] = 'Btf3l4'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q9CQH7'),'description'] = 'basic transcription factor 3-like 4 [Source:MGI Symbol;Acc:MGI:1915312]'
# http://www.uniprot.org/uniprot/Q9CQV1; http://www.informatics.jax.org/marker/MGI:1913699; https://www.ncbi.nlm.nih.gov/gene/66449
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q9CQV1'),'entrezgene'] = c(66449)
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q9CQV1'),'mgi_id'] = 'MGI:1913699'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q9CQV1'),'mgi_symbol'] = 'Pam16'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q9CQV1'),'description'] = 'presequence translocase-asssociated motor 16 homolog (S. cerevisiae) [Source:MGI Symbol;Acc:MGI:1913699]'
# http://www.uniprot.org/uniprot/Q9CW46; http://www.informatics.jax.org/marker/MGI:1919016; https://www.ncbi.nlm.nih.gov/gene/71766
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q9CW46'),'entrezgene'] = c(71766)
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q9CW46'),'mgi_id'] = 'MGI:1919016'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q9CW46'),'mgi_symbol'] = 'Raver1'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q9CW46'),'description'] = 'ribonucleoprotein, PTB-binding 1 [Source:MGI Symbol;Acc:MGI:1919016]'
# http://www.uniprot.org/uniprot/Q9D1C8; http://www.informatics.jax.org/marker/MGI:1914164; https://www.ncbi.nlm.nih.gov/gene/66914
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q9D1C8'),'entrezgene'] = c(66914)
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q9D1C8'),'mgi_id'] = 'MGI:1914164'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q9D1C8'),'mgi_symbol'] = 'Vps28'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q9D1C8'),'description'] = 'vacuolar protein sorting 28 [Source:MGI Symbol;Acc:MGI:1914164]'
# http://www.uniprot.org/uniprot/Q9D1J3; http://www.informatics.jax.org/marker/MGI:1913368; https://www.ncbi.nlm.nih.gov/gene/66118
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q9D1J3'),'entrezgene'] = c(66118)
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q9D1J3'),'mgi_id'] = 'MGI:1913368'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q9D1J3'),'mgi_symbol'] = 'Sarnp'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q9D1J3'),'description'] = 'SAP domain containing ribonucleoprotein [Source:MGI Symbol;Acc:MGI:1913368]'
# http://www.uniprot.org/uniprot/Q9D868; http://www.informatics.jax.org/marker/MGI:106499; https://www.ncbi.nlm.nih.gov/gene/66101
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q9D868'),'entrezgene'] = c(66101)
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q9D868'),'mgi_id'] = 'MGI:106499'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q9D868'),'mgi_symbol'] = 'Ppih'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q9D868'),'description'] = 'peptidyl prolyl isomerase H [Source:MGI Symbol;Acc:MGI:106499]'
# http://www.uniprot.org/uniprot/Q9R1Q7; http://www.informatics.jax.org/marker/MGI:1298382; https://www.ncbi.nlm.nih.gov/gene/18824
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q9R1Q7'),'entrezgene'] = c(18824)
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q9R1Q7'),'mgi_id'] = 'MGI:1298382'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q9R1Q7'),'mgi_symbol'] = 'Plp2'
EVpedia_Top100.Protein_Mouse.mgi_id.table[which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']=='Q9R1Q7'),'description'] = 'proteolipid protein 2 [Source:MGI Symbol;Acc:MGI:1298382]'

EVpedia_Top100.Protein_Mouse.mgi_id.table[,'description'] = trimws(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'description'])
EVpedia_Top100.Protein_Mouse.mgi_id.table[,'description'] = sapply(strsplit(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'description'],'\\|'), FUN=function(U) paste(unique(sapply(strsplit(unlist(U),'\\[Source:'), FUN=function(V) trimws(unlist(V))[1])),collapse="|"))
EVpedia_Top100.Protein_Mouse.mgi_id.table = unique(EVpedia_Top100.Protein_Mouse.mgi_id.table)

EVpedia_Top100.Protein_Mouse.mgi_id.table = EVpedia_Top100.Protein_Mouse.mgi_id.table[order(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot'],EVpedia_Top100.Protein_Mouse.mgi_id.table[,'mgi_symbol'],decreasing=FALSE),]
EVpedia_Top100.Protein_Mouse.mgi_id.table = EVpedia_Top100.Protein_Mouse.mgi_id.table[,c('mgi_symbol','mgi_id','entrezgene','uniprotswissprot','description')]
EVpedia_Top100.Protein_Mouse.mgi_id.table = EVpedia_Top100.Protein_Mouse.mgi_id.table[order(match(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot'], a)),]

# --- uniprotswissprot not found by R package "biomaRt" ---
query.uniprotswissprot = EVpedia_Top100.Protein_Mouse.UniProt_accession
subject.uniprotswissprot = EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']
query.uniprotswissprot.not_in_biomaRt.idx = which(! query.uniprotswissprot %in% subject.uniprotswissprot)
query.uniprotswissprot.not_in_biomaRt = query.uniprotswissprot[query.uniprotswissprot.not_in_biomaRt.idx]

if(length(query.uniprotswissprot.not_in_biomaRt) > 0)
{
	query.uniprotswissprot.not_in_biomaRt.mgi_id.table <- c()
	query.uniprotswissprot.not_in_biomaRt_BridgeDbR <- c()
	for(i in 1:length(query.uniprotswissprot.not_in_biomaRt))
	{
		UniprotTrEMBL2MGI = map(mapper.mmu, "S", query.uniprotswissprot.not_in_biomaRt[i], "M")	#Uniprot-TrEMBL (system code: S); MGI (system code: M)
		#print(UniprotTrEMBL2MGI)
		
		if(query.uniprotswissprot.not_in_biomaRt[i] == 'Q54AC6') {			# http://www.uniprot.org/uniprot/Q54AC6; https://www.ncbi.nlm.nih.gov/gene/12856; http://www.informatics.jax.org/marker/MGI:1333806
			UniprotTrEMBL2MGI = 'MGI:1333806'
		} else if(query.uniprotswissprot.not_in_biomaRt[i] == 'A2BGI9') {	# http://www.uniprot.org/uniprot/A2BGI9; http://www.informatics.jax.org/marker/MGI:106499 (http://www.informatics.jax.org/marker/MGI:3645078)
			UniprotTrEMBL2MGI = 'MGI:106499'
		} else if(query.uniprotswissprot.not_in_biomaRt[i] == 'B0V3P3') {	# http://www.uniprot.org/uniprot/B0V3P3; http://www.informatics.jax.org/marker/MGI:1337130 (http://www.informatics.jax.org/marker/MGI:5663201)
			UniprotTrEMBL2MGI = 'MGI:1337130'
		} else if(query.uniprotswissprot.not_in_biomaRt[i] == 'Q9CXU4') {	# http://www.uniprot.org/uniprot/Q9CXU4; https://www.ncbi.nlm.nih.gov/gene/53600; http://www.informatics.jax.org/marker/MGI:1858317 (http://www.informatics.jax.org/marker/MGI:3704362)
			UniprotTrEMBL2MGI = 'MGI:1858317'
		} else if(query.uniprotswissprot.not_in_biomaRt[i] == 'E9PY90') {	# http://www.uniprot.org/uniprot/E9PY90; https://www.ncbi.nlm.nih.gov/gene/338320; http://www.informatics.jax.org/marker/MGI:2159614 (http://www.informatics.jax.org/marker/MGI:1346056)
			UniprotTrEMBL2MGI = 'MGI:2159614'
		}
		#print(UniprotTrEMBL2MGI)
		
		if(length(UniprotTrEMBL2MGI) > 0) {
			a = UniprotTrEMBL2MGI
			b = ensembl.mmu
			query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes <- getBM(attributes=c('mgi_symbol','mgi_id','entrezgene'), filters='mgi_id', values=a, mart=b)
			#print(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes)
			
			query.uniprotswissprot.not_in_biomaRt.mgi_id.wikigene_description <- getBM(attributes=c('mgi_symbol','wikigene_description'), filters='mgi_id', values=a, mart=b)
			query.uniprotswissprot.not_in_biomaRt.mgi_id.wikigene_description.table = merge(query.uniprotswissprot.not_in_biomaRt.mgi_id.wikigene_description, query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes, by.x="mgi_symbol", by.y="mgi_symbol")
			query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes = query.uniprotswissprot.not_in_biomaRt.mgi_id.wikigene_description.table[,c('mgi_symbol','mgi_id','entrezgene','wikigene_description'),drop=FALSE]
			colnames(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes)[which(colnames(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes)=='wikigene_description')] = 'description'
			query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes = unique(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes)
			#print(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes)
			
			if(nrow(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes) > 0)
			{
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:87994'),'entrezgene'] = c(11674)			# http://www.informatics.jax.org/marker/MGI:87994; https://www.ncbi.nlm.nih.gov/gene/11674
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:95655'),'entrezgene'] = c(14451)			# http://www.informatics.jax.org/marker/MGI:95655; https://www.ncbi.nlm.nih.gov/gene/14451
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:96607'),'entrezgene'] = c(16409)			# http://www.informatics.jax.org/marker/MGI:96607; https://www.ncbi.nlm.nih.gov/gene/16409
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:97298'),'entrezgene'] = c(18000)			# http://www.informatics.jax.org/marker/MGI:97298; https://www.ncbi.nlm.nih.gov/gene/18000
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:97843'),'entrezgene'] = c(19339)			# http://www.informatics.jax.org/marker/MGI:97843; https://www.ncbi.nlm.nih.gov/gene/19339
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:99705'),'entrezgene'] = c(621054)		# http://www.informatics.jax.org/marker/MGI:99705; https://www.ncbi.nlm.nih.gov/gene/621054
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:107686'),'entrezgene'] = c(17992)		# http://www.informatics.jax.org/marker/MGI:107686; https://www.ncbi.nlm.nih.gov/gene/17992
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:1335093'),'entrezgene'] = c(19054)		# http://www.informatics.jax.org/marker/MGI:1335093; https://www.ncbi.nlm.nih.gov/gene/19054
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:1336201'),'entrezgene'] = c(20610)		# http://www.informatics.jax.org/marker/MGI:1336201; https://www.ncbi.nlm.nih.gov/gene/20610
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:1347049'),'entrezgene'] = c(26372)		# http://www.informatics.jax.org/marker/MGI:1347049; https://www.ncbi.nlm.nih.gov/gene/26372
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:1913712'),'entrezgene'] = c(66462)		# http://www.informatics.jax.org/marker/MGI:1913712; https://www.ncbi.nlm.nih.gov/gene/66462
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:1915896'),'entrezgene'] = c(68646)		# http://www.informatics.jax.org/marker/MGI:1915896; https://www.ncbi.nlm.nih.gov/gene/68646
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:1915964'),'entrezgene'] = c(102115)		# http://www.informatics.jax.org/marker/MGI:1915964; https://www.ncbi.nlm.nih.gov/gene/102115
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:1920994'),'entrezgene'] = c(73744)		# http://www.informatics.jax.org/marker/MGI:1920994; https://www.ncbi.nlm.nih.gov/gene/73744
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:2144831'),'entrezgene'] = c(100041194)	# http://www.informatics.jax.org/marker/MGI:2144831; https://www.ncbi.nlm.nih.gov/gene/100041194
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:2152971'),'entrezgene'] = c(116837)		# http://www.informatics.jax.org/marker/MGI:2152971; https://www.ncbi.nlm.nih.gov/gene/116837
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:2159728'),'entrezgene'] = c(171429)		# http://www.informatics.jax.org/marker/MGI:2159728; https://www.ncbi.nlm.nih.gov/gene/171429
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:2385199'),'entrezgene'] = c(230500)		# http://www.informatics.jax.org/marker/MGI:2385199; https://www.ncbi.nlm.nih.gov/gene/230500
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:2443183'),'entrezgene'] = c(338366)		# http://www.informatics.jax.org/marker/MGI:2443183; https://www.ncbi.nlm.nih.gov/gene/338366
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:2652858'),'entrezgene'] = c(269211)		# http://www.informatics.jax.org/marker/MGI:2652858; https://www.ncbi.nlm.nih.gov/gene/269211
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:3641908'),'entrezgene'] = c(100042734)	# http://www.informatics.jax.org/marker/MGI:3641908; https://www.ncbi.nlm.nih.gov/gene/100042734
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:3642192'),'entrezgene'] = c(100039748)	# http://www.informatics.jax.org/marker/MGI:3642192; https://www.ncbi.nlm.nih.gov/gene/100039748
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:3642279'),'entrezgene'] = c(100041265)	# http://www.informatics.jax.org/marker/MGI:3642279; https://www.ncbi.nlm.nih.gov/gene/100041265
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:3642298'),'entrezgene'] = c(100039740)	# http://www.informatics.jax.org/marker/MGI:3642298; https://www.ncbi.nlm.nih.gov/gene/100039740
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:3642334'),'entrezgene'] = c(100040551)	# http://www.informatics.jax.org/marker/MGI:3642334; https://www.ncbi.nlm.nih.gov/gene/100040551
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:3642408'),'entrezgene'] = c(100039281)	# http://www.informatics.jax.org/marker/MGI:3642408; https://www.ncbi.nlm.nih.gov/gene/100039281
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:3642754'),'entrezgene'] = c(100039782)	# http://www.informatics.jax.org/marker/MGI:3642754; https://www.ncbi.nlm.nih.gov/gene/100039782
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:3642825'),'entrezgene'] = c(100042864)	# http://www.informatics.jax.org/marker/MGI:3642825; https://www.ncbi.nlm.nih.gov/gene/100042864
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:3643173'),'entrezgene'] = c(432466)		# http://www.informatics.jax.org/marker/MGI:3643173; https://www.ncbi.nlm.nih.gov/gene/432466
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:3643769'),'entrezgene'] = c(667759)		# http://www.informatics.jax.org/marker/MGI:3643769; https://www.ncbi.nlm.nih.gov/gene/667759
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:3645078'),'entrezgene'] = c(665989)		# http://www.informatics.jax.org/marker/MGI:3645078; https://www.ncbi.nlm.nih.gov/gene/665989
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:3645172'),'entrezgene'] = c(668559)		# http://www.informatics.jax.org/marker/MGI:3645172; https://www.ncbi.nlm.nih.gov/gene/668559
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:3645809'),'entrezgene'] = c(667073)		# http://www.informatics.jax.org/marker/MGI:3645809; https://www.ncbi.nlm.nih.gov/gene/667073
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:3646032'),'entrezgene'] = c(667728)		# http://www.informatics.jax.org/marker/MGI:3646032; https://www.ncbi.nlm.nih.gov/gene/667728
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:3646561'),'entrezgene'] = c(631906)		# http://www.informatics.jax.org/marker/MGI:3646561; https://www.ncbi.nlm.nih.gov/gene/631906
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:3647118'),'entrezgene'] = c(434375)		# http://www.informatics.jax.org/marker/MGI:3647118; https://www.ncbi.nlm.nih.gov/gene/434375
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:3647789'),'entrezgene'] = c(432502)		# http://www.informatics.jax.org/marker/MGI:3647789; https://www.ncbi.nlm.nih.gov/gene/432502
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:3647835'),'entrezgene'] = c(100503759)	# http://www.informatics.jax.org/marker/MGI:3647835; https://www.ncbi.nlm.nih.gov/gene/100503759
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:3647964'),'entrezgene'] = c(666974)		# http://www.informatics.jax.org/marker/MGI:3647964; https://www.ncbi.nlm.nih.gov/gene/666974
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:3648085'),'entrezgene'] = c(266459)		# http://www.informatics.jax.org/marker/MGI:3648085; https://www.ncbi.nlm.nih.gov/gene/266459
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:3648259'),'entrezgene'] = c(627557)		# http://www.informatics.jax.org/marker/MGI:3648259; https://www.ncbi.nlm.nih.gov/gene/627557
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:3648528'),'entrezgene'] = c(382096)		# http://www.informatics.jax.org/marker/MGI:3648528; https://www.ncbi.nlm.nih.gov/gene/382096
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:3648910'),'entrezgene'] = c(545208)		# http://www.informatics.jax.org/marker/MGI:3648910; https://www.ncbi.nlm.nih.gov/gene/545208
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:3649357'),'entrezgene'] = c(435692)		# http://www.informatics.jax.org/marker/MGI:3649357; https://www.ncbi.nlm.nih.gov/gene/435692
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:3649411'),'entrezgene'] = c(100040340)	# http://www.informatics.jax.org/marker/MGI:3649411; https://www.ncbi.nlm.nih.gov/gene/100040340
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:3649769'),'entrezgene'] = c(100502754)	# http://www.informatics.jax.org/marker/MGI:3649769; https://www.ncbi.nlm.nih.gov/gene/100502754
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:3687012'),'entrezgene'] = c(436522)		# http://www.informatics.jax.org/marker/MGI:3687012; https://www.ncbi.nlm.nih.gov/gene/436522
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:3702053'),'entrezgene'] = c(668119)		# http://www.informatics.jax.org/marker/MGI:3702053; https://www.ncbi.nlm.nih.gov/gene/668119
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:3702055'),'entrezgene'] = c(668100)		# http://www.informatics.jax.org/marker/MGI:3702055; https://www.ncbi.nlm.nih.gov/gene/668100
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:3702064'),'entrezgene'] = c(668113)		# http://www.informatics.jax.org/marker/MGI:3702064; https://www.ncbi.nlm.nih.gov/gene/668113
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:3702066'),'entrezgene'] = c(545432)		# http://www.informatics.jax.org/marker/MGI:3702066; https://www.ncbi.nlm.nih.gov/gene/545432
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:3702068'),'entrezgene'] = c(668107)		# http://www.informatics.jax.org/marker/MGI:3702068; https://www.ncbi.nlm.nih.gov/gene/668107
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:3702070'),'entrezgene'] = c(668115)		# http://www.informatics.jax.org/marker/MGI:3702070; https://www.ncbi.nlm.nih.gov/gene/668115
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:3702072'),'entrezgene'] = c(668096)		# http://www.informatics.jax.org/marker/MGI:3702072; https://www.ncbi.nlm.nih.gov/gene/668096
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:3704322'),'entrezgene'] = c(100042812)	# http://www.informatics.jax.org/marker/MGI:3704322; https://www.ncbi.nlm.nih.gov/gene/100042812
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:3704362'),'entrezgene'] = c(100039174)	# http://www.informatics.jax.org/marker/MGI:3704362; https://www.ncbi.nlm.nih.gov/gene/100039174
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:3704479'),'entrezgene'] = c(15181)		# http://www.informatics.jax.org/marker/MGI:3704479; https://www.ncbi.nlm.nih.gov/gene/15181
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:3704493'),'entrezgene'] = c(672195)		# http://www.informatics.jax.org/marker/MGI:3704493; https://www.ncbi.nlm.nih.gov/gene/672195
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:3704503'),'entrezgene'] = c(100034727)	# http://www.informatics.jax.org/marker/MGI:3704503; https://www.ncbi.nlm.nih.gov/gene/100034727
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:3708114'),'entrezgene'] = c(666488)		# http://www.informatics.jax.org/marker/MGI:3708114; https://www.ncbi.nlm.nih.gov/gene/666488
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:3708621'),'entrezgene'] = c(100039316)	# http://www.informatics.jax.org/marker/MGI:3708621; https://www.ncbi.nlm.nih.gov/gene/100039316
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:3779903'),'entrezgene'] = c(102638448)	# http://www.informatics.jax.org/marker/MGI:3779903; https://www.ncbi.nlm.nih.gov/gene/102638448
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:3781727'),'entrezgene'] = c(100041859)	# http://www.informatics.jax.org/marker/MGI:3781727; https://www.ncbi.nlm.nih.gov/gene/100041859
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:3782011'),'entrezgene'] = c(100042427)	# http://www.informatics.jax.org/marker/MGI:3782011; https://www.ncbi.nlm.nih.gov/gene/100042427
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:3822538'),'entrezgene'] = c(414125)		# http://www.informatics.jax.org/marker/MGI:3822538; https://www.ncbi.nlm.nih.gov/gene/414125
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:4937303'),'entrezgene'] = c(108168363)	# http://www.informatics.jax.org/marker/MGI:4937303; https://www.ncbi.nlm.nih.gov/gene/108168363
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:5439454'),'entrezgene'] = c(108168782)	# http://www.informatics.jax.org/marker/MGI:5439454; https://www.ncbi.nlm.nih.gov/gene/108168782
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:3641718'),'entrezgene'] = NA				# http://www.informatics.jax.org/marker/MGI:3641718
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:4937062'),'entrezgene'] = NA				# http://www.informatics.jax.org/marker/MGI:4937062
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:5141855'),'entrezgene'] = NA				# http://www.informatics.jax.org/marker/MGI:5141855
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:5141890'),'entrezgene'] = NA				# http://www.informatics.jax.org/marker/MGI:5141890
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:5141896'),'entrezgene'] = NA				# http://www.informatics.jax.org/marker/MGI:5141896
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:5141906'),'entrezgene'] = NA				# http://www.informatics.jax.org/marker/MGI:5141906
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:5141968'),'entrezgene'] = NA				# http://www.informatics.jax.org/marker/MGI:5141968
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:5141983'),'entrezgene'] = NA				# http://www.informatics.jax.org/marker/MGI:5141983
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:5141986'),'entrezgene'] = NA				# http://www.informatics.jax.org/marker/MGI:5141986
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:5142012'),'entrezgene'] = NA				# http://www.informatics.jax.org/marker/MGI:5142012
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:5313118'),'entrezgene'] = NA				# http://www.informatics.jax.org/marker/MGI:5313118
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:5439456'),'entrezgene'] = NA				# http://www.informatics.jax.org/marker/MGI:5439456
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:5504144'),'entrezgene'] = NA				# http://www.informatics.jax.org/marker/MGI:5504144
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:5547772'),'entrezgene'] = NA				# http://www.informatics.jax.org/marker/MGI:5547772
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:5547784'),'entrezgene'] = NA				# http://www.informatics.jax.org/marker/MGI:5547784
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:5662811'),'entrezgene'] = NA				# http://www.informatics.jax.org/marker/MGI:5662811
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:5663201'),'entrezgene'] = NA				# http://www.informatics.jax.org/marker/MGI:5663201
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[which(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes[,'mgi_id']=='MGI:5804914'),'entrezgene'] = NA				# http://www.informatics.jax.org/marker/MGI:5804914
				#print(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes)
				
				query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes = unique(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes)
				#print(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes)
				
				query.uniprotswissprot.not_in_biomaRt.mgi_id = cbind(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes, query.uniprotswissprot.not_in_biomaRt[i])
				colnames(query.uniprotswissprot.not_in_biomaRt.mgi_id) = c(colnames(query.uniprotswissprot.not_in_biomaRt.mgi_id.attributes), 'uniprotswissprot')
				query.uniprotswissprot.not_in_biomaRt.mgi_id = query.uniprotswissprot.not_in_biomaRt.mgi_id[,c('mgi_symbol','mgi_id','entrezgene','uniprotswissprot','description')]
				query.uniprotswissprot.not_in_biomaRt.mgi_id.table = rbind(query.uniprotswissprot.not_in_biomaRt.mgi_id.table, query.uniprotswissprot.not_in_biomaRt.mgi_id)
			}
		} else {
			query.uniprotswissprot.not_in_biomaRt_BridgeDbR <- c(query.uniprotswissprot.not_in_biomaRt_BridgeDbR, query.uniprotswissprot.not_in_biomaRt[i])
		}
	}
	#query.uniprotswissprot.not_in_biomaRt.mgi_id.table
	#query.uniprotswissprot.not_in_biomaRt_BridgeDbR
	
	EVpedia_Top100.Protein_Mouse.mgi_id.table = rbind(EVpedia_Top100.Protein_Mouse.mgi_id.table, query.uniprotswissprot.not_in_biomaRt.mgi_id.table)
}

# --- re-order ---
a = EVpedia_Top100.Protein_Mouse.UniProt_accession
EVpedia_Top100.Protein_Mouse.mgi_id.table = EVpedia_Top100.Protein_Mouse.mgi_id.table[order(match(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot'], a)),]

# --- add in entrezgene, missed description of entrezgene, and correct mgi_id, mgi_symbol, description, and uniprotswissprot of entrezgene ---
for(i in 1:nrow(EVpedia_Top100.Protein_Mouse.mgi_id.table))
{
	if(is.na(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'entrezgene'])) {
		if(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_id'] == 'MGI:87994') {
			EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'entrezgene'] = c(11674)				# Aldoa (http://www.informatics.jax.org/marker/MGI:87994); https://www.ncbi.nlm.nih.gov/gene/11674
		} else if(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_id'] == 'MGI:96607') {		# http://www.uniprot.org/uniprot/G5E8F1 -> Not found by "biomaRt", but could be found by "BridgeDbR" R package.
			EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'entrezgene'] = c(16409)				# Itgam (http://www.informatics.jax.org/marker/MGI:96607); https://www.ncbi.nlm.nih.gov/gene/16409
		} else if(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_id'] == 'MGI:104681') {
			EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'entrezgene'] = c(15239)				# Hgs (http://www.informatics.jax.org/marker/MGI:104681); https://www.ncbi.nlm.nih.gov/gene/15239
		} else if(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_id'] == 'MGI:109279') {
			EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'entrezgene'] = c(18115)				# Nnt (http://www.informatics.jax.org/marker/MGI:109279); https://www.ncbi.nlm.nih.gov/gene/18115
		} else if(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_id'] == 'MGI:894645') {
			EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'entrezgene'] = c(18693)				# Pick1 (http://www.informatics.jax.org/marker/MGI:894645); https://www.ncbi.nlm.nih.gov/gene/18693
		} else if(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_id'] == 'MGI:1096365') {
			EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'entrezgene'] = c(19188)				# Psme2 (http://www.informatics.jax.org/marker/MGI:1096365); https://www.ncbi.nlm.nih.gov/gene/19188
		} else if(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_id'] == 'MGI:1919308') {
			EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'entrezgene'] = c(72058)				# Igsf5 (http://www.informatics.jax.org/marker/MGI:1919308); https://www.ncbi.nlm.nih.gov/gene/72058
		} else if(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_id'] == 'MGI:1928344') {
			EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'entrezgene'] = c(56433)				# Vps29 (http://www.informatics.jax.org/marker/MGI:1928344); https://www.ncbi.nlm.nih.gov/gene/56433
		} else if(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_id'] == 'MGI:3649356') {
			EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'entrezgene'] = c(666586)				# Gm11808 (http://www.informatics.jax.org/marker/MGI:3649356); https://www.ncbi.nlm.nih.gov/gene/666586
		} else {
			if(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_id'] == 'MGI:3641718') {		# http://www.informatics.jax.org/marker/MGI:3641718
				EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'entrezgene'] = NA
			} else if(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_id'] == 'MGI:4937062') {	# http://www.informatics.jax.org/marker/MGI:4937062
				EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'entrezgene'] = NA
			} else if(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_id'] == 'MGI:5141855') {	# http://www.informatics.jax.org/marker/MGI:5141855
				EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'entrezgene'] = NA
			} else if(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_id'] == 'MGI:5141890') {	# http://www.informatics.jax.org/marker/MGI:5141890
				EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'entrezgene'] = NA
			} else if(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_id'] == 'MGI:5141896') {	# http://www.informatics.jax.org/marker/MGI:5141896
				EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'entrezgene'] = NA
			} else if(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_id'] == 'MGI:5141906') {	# http://www.informatics.jax.org/marker/MGI:5141906
				EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'entrezgene'] = NA
			} else if(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_id'] == 'MGI:5141968') {	# http://www.informatics.jax.org/marker/MGI:5141968
				EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'entrezgene'] = NA
			} else if(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_id'] == 'MGI:5141983') {	# http://www.informatics.jax.org/marker/MGI:5141983
				EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'entrezgene'] = NA
			} else if(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_id'] == 'MGI:5141986') {	# http://www.informatics.jax.org/marker/MGI:5141986
				EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'entrezgene'] = NA
			} else if(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_id'] == 'MGI:5142012') {	# http://www.informatics.jax.org/marker/MGI:5142012
				EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'entrezgene'] = NA
			} else if(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_id'] == 'MGI:5313118') {	# http://www.informatics.jax.org/marker/MGI:5313118
				EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'entrezgene'] = NA
			} else if(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_id'] == 'MGI:5439456') {	# http://www.informatics.jax.org/marker/MGI:5439456
				EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'entrezgene'] = NA
			} else if(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_id'] == 'MGI:5504144') {	# http://www.informatics.jax.org/marker/MGI:5504144
				EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'entrezgene'] = NA
			} else if(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_id'] == 'MGI:5547772') {	# http://www.informatics.jax.org/marker/MGI:5547772
				EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'entrezgene'] = NA
			} else if(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_id'] == 'MGI:5547784') {	# http://www.informatics.jax.org/marker/MGI:5547784
				EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'entrezgene'] = NA
			} else if(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_id'] == 'MGI:5662811') {	# http://www.informatics.jax.org/marker/MGI:5662811
				EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'entrezgene'] = NA
			} else if(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_id'] == 'MGI:5663201') {	# http://www.informatics.jax.org/marker/MGI:5663201
				EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'entrezgene'] = NA
			} else if(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_id'] == 'MGI:5804914') {	# http://www.informatics.jax.org/marker/MGI:5804914
				EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'entrezgene'] = NA
			} else if(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_id'] == 'MGI:6121499') {	# http://www.informatics.jax.org/marker/MGI:6121499
				EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'entrezgene'] = NA
			} else if(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_id'] == 'MGI:6121530') {	# http://www.informatics.jax.org/marker/MGI:6121530
				EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'entrezgene'] = NA
			} else if(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_id'] == 'MGI:6121629') {	# http://www.informatics.jax.org/marker/MGI:6121629
				EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'entrezgene'] = NA
			} else {
				cat(paste('[Warning] Please check entrezgene of mgi_id "',EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_id'],'" of uniprotswissprot "',EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'uniprotswissprot'],'". Entrezgene is NA (Not Available) by R package "biomaRt".\n',sep=""))
			}
		}
	} else {
		if(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'entrezgene'] == c(22668)) {				# https://www.ncbi.nlm.nih.gov/gene/22668; http://www.informatics.jax.org/marker/MGI:1095403; http://www.uniprot.org/uniprot/Q64213
			EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_id'] = 'MGI:1095403'
			EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_symbol'] = 'Sf1'
			EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'description'] = 'splicing factor 1 [Source:MGI Symbol;Acc:MGI:1095403]'
			EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'uniprotswissprot'] = 'Q64213'
		} else if(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'entrezgene'] == c(70472)) {		# https://www.ncbi.nlm.nih.gov/gene/70472; http://www.informatics.jax.org/marker/MGI:1917722; http://www.uniprot.org/uniprot/Q8CDM1, http://www.uniprot.org/uniprot/G3X963
			EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_id'] = 'MGI:1917722'
			EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_symbol'] = 'Atad2'
			EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'description'] = 'ATPase family, AAA domain containing 2 [Source:MGI Symbol;Acc:MGI:1917722]'
			EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'uniprotswissprot'] = paste(c('Q8CDM1','G3X963'),collapse="|")
		} else if(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'entrezgene'] == c(320207)) {			# https://www.ncbi.nlm.nih.gov/gene/320207; http://www.informatics.jax.org/marker/MGI:2443588; http://www.uniprot.org/uniprot/Q5SW28
			EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_id'] = 'MGI:2443588'
			EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_symbol'] = 'Pik3r5'
			EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'description'] = 'phosphoinositide-3-kinase regulatory subunit 5 [Source:MGI Symbol;Acc:MGI:2443588]'
			EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'uniprotswissprot'] = 'Q5SW28'
		} else if(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'entrezgene'] == c(104709)) {			# https://www.ncbi.nlm.nih.gov/gene/104709; http://www.informatics.jax.org/marker/MGI:2144613; http://www.uniprot.org/uniprot/Q3U6Q4
			EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_id'] = 'MGI:2144613'
			EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_symbol'] = 'Pik3r6'
			EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'description'] = 'phosphoinositide-3-kinase regulatory subunit 6 [Source:MGI Symbol;Acc:MGI:2144613]'	# wrong name in MGI: phosphoinositide-3-kinase regulatory subunit 5 -> correct name in UniProt: Phosphoinositide 3-kinase regulatory subunit 6
			EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'uniprotswissprot'] = 'Q3U6Q4'
		} else if(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'entrezgene'] == c(99889)) {		# https://www.ncbi.nlm.nih.gov/gene/99889; http://www.informatics.jax.org/marker/MGI:1277120
			EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_id'] = 'MGI:1277120'
			EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_symbol'] = 'Arfip1'
			EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'description'] = 'ADP-ribosylation factor interacting protein 1 [Source:MGI Symbol;Acc:MGI:1277120]'
		} else if(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'entrezgene'] == c(619960)) {		# https://www.ncbi.nlm.nih.gov/gene/619960; http://www.informatics.jax.org/marker/MGI:4439772
			EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_id'] = 'MGI:4439772'
			EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_symbol'] = 'Igkv12-41'
			EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'description'] = 'immunoglobulin kappa chain variable 12-41 [Source:MGI Symbol;Acc:MGI:4439772]'
		} else if(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'entrezgene'] == c(545432)) {		# https://www.ncbi.nlm.nih.gov/gene/545432; http://www.informatics.jax.org/marker/MGI:3702066
			EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_id'] = 'MGI:3702066'
			EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_symbol'] = 'Gm13695'
			EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'description'] = 'predicted gene 13695 [Source:MGI Symbol;Acc:MGI:3702066]'
		} else if(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'entrezgene'] == c(668096)) {		# https://www.ncbi.nlm.nih.gov/gene/668096; http://www.informatics.jax.org/marker/MGI:3702072
			EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_id'] = 'MGI:3702072'
			EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_symbol'] = 'Gm13698'
			EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'description'] = 'predicted gene 13698 [Source:MGI Symbol;Acc:MGI:3702072]'
		} else if(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'entrezgene'] == c(668100)) {		# https://www.ncbi.nlm.nih.gov/gene/668100; http://www.informatics.jax.org/marker/MGI:3702055
			EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_id'] = 'MGI:3702055'
			EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_symbol'] = 'Gm13693'
			EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'description'] = 'predicted gene 13693 [Source:MGI Symbol;Acc:MGI:3702055]'
		} else if(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'entrezgene'] == c(668107)) {		# https://www.ncbi.nlm.nih.gov/gene/668107; http://www.informatics.jax.org/marker/MGI:3702068
			EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_id'] = 'MGI:3702068'
			EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_symbol'] = 'Gm13696'
			EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'description'] = 'predicted gene 13696 [Source:MGI Symbol;Acc:MGI:3702068]'
		} else if(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'entrezgene'] == c(668113)) {		# https://www.ncbi.nlm.nih.gov/gene/668113; http://www.informatics.jax.org/marker/MGI:3702064
			EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_id'] = 'MGI:3702064'
			EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_symbol'] = 'Gm13694'
			EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'description'] = 'predicted gene 13694 [Source:MGI Symbol;Acc:MGI:3702064]'
		} else if(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'entrezgene'] == c(668115)) {		# https://www.ncbi.nlm.nih.gov/gene/668115; http://www.informatics.jax.org/marker/MGI:3702070
			EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_id'] = 'MGI:3702070'
			EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_symbol'] = 'Gm13697'
			EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'description'] = 'predicted gene 13697 [Source:MGI Symbol;Acc:MGI:3702070]'
		} else if(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'entrezgene'] == c(668119)) {		# https://www.ncbi.nlm.nih.gov/gene/668119; http://www.informatics.jax.org/marker/MGI:3702053
			EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_id'] = 'MGI:3702053'
			EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_symbol'] = 'Gm13691'
			EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'description'] = 'predicted gene 13691 [Source:MGI Symbol;Acc:MGI:3702053]'
		}
		
		EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'description'] = trimws(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'description'])
		EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'description'] = sapply(strsplit(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'description'],'\\|'), FUN=function(U) paste(unique(sapply(strsplit(unlist(U),'\\[Source:'), FUN=function(V) trimws(unlist(V))[1])),collapse="|"))
	}
}

# --- add in missed description of mgi_id ---
for(i in 1:nrow(EVpedia_Top100.Protein_Mouse.mgi_id.table))
{
	if(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_id'] == 'MGI:5141983') {				# http://www.informatics.jax.org/marker/MGI:5141983
		EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'entrezgene'] = NA
		EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_symbol'] = 'Gm20518'
		EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'description'] = 'predicted gene 20518 [Source:MGI Symbol;Acc:MGI:5141983]'
	} else if(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_id'] == 'MGI:5141986') {			# http://www.informatics.jax.org/marker/MGI:5141986
		EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'entrezgene'] = NA
		EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_symbol'] = 'Gm20521'
		EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'description'] = 'predicted gene 20521 [Source:MGI Symbol;Acc:MGI:5141986]'
	} else if(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_id'] == 'MGI:5504144') {			# http://www.informatics.jax.org/marker/MGI:5504144
		EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'entrezgene'] = NA
		EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_symbol'] = 'Gm27029'
		EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'description'] = 'predicted gene, 27029 [Source:MGI Symbol;Acc:MGI:5504144]'
	} else if(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_id'] == 'MGI:5547772') {			# http://www.informatics.jax.org/marker/MGI:5547772
		EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'entrezgene'] = NA
		EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_symbol'] = 'Gm28036'
		EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'description'] = 'predicted gene, 28036 [Source:MGI Symbol;Acc:MGI:5547772]'
	} else if(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_id'] == 'MGI:5547784') {			# http://www.informatics.jax.org/marker/MGI:5547784
		EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'entrezgene'] = NA
		EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_symbol'] = 'Gm28048'
		EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'description'] = 'predicted gene, 28048 [Source:MGI Symbol;Acc:MGI:5547784]'
	} else if(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_id'] == 'MGI:5804914') {			# http://www.informatics.jax.org/marker/MGI:5804914
		EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'entrezgene'] = NA
		EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_symbol'] = 'Gm45799'
		EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'description'] = 'predicted gene 45799 [Source:MGI Symbol;Acc:MGI:5804914]'
	}
	
	EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'description'] = trimws(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'description'])
	EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'description'] = sapply(strsplit(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'description'],'\\|'), FUN=function(U) paste(unique(sapply(strsplit(unlist(U),'\\[Source:'), FUN=function(V) trimws(unlist(V))[1])),collapse="|"))
}

# --- revise wrong entrezgene ---
for(i in 1:nrow(EVpedia_Top100.Protein_Mouse.mgi_id.table))
{
	if(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'uniprotswissprot'] == 'P62889') {			# http://www.uniprot.org/uniprot/P62889; Rpl30 (https://www.ncbi.nlm.nih.gov/gene/19946)
		EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'entrezgene'] = c(19946)
	} else if(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'uniprotswissprot'] == 'Q64281') {	# http://www.uniprot.org/uniprot/Q64281; Lilrb4a (https://www.ncbi.nlm.nih.gov/gene/14728)
		EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'entrezgene'] = c(14728)
	} else if(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'uniprotswissprot'] == 'P84244') {	# http://www.uniprot.org/uniprot/P84244
		if(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_id'] == 'MGI:1097686') {			# http://www.informatics.jax.org/marker/MGI:1097686; https://www.ncbi.nlm.nih.gov/gene/15078
			EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'entrezgene'] = c(15078)
		} else if(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_id'] == 'MGI:1101768') {		# http://www.informatics.jax.org/marker/MGI:1101768; https://www.ncbi.nlm.nih.gov/gene/15081
			EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'entrezgene'] = c(15081)
		}
	} else if(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_id'] == 'MGI:95896') {			# http://www.informatics.jax.org/marker/MGI:95896; H2-D1 (https://www.ncbi.nlm.nih.gov/gene/14964)
		EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'entrezgene'] = c(14964)
	}
}

# --- replace Synonyms to Symbol (Official Symbol) ---
for(i in 1:nrow(EVpedia_Top100.Protein_Mouse.mgi_id.table))
{
	if(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_id'] == 'MGI:1922471') {				# http://www.informatics.jax.org/marker/MGI:1922471; https://www.ncbi.nlm.nih.gov/gene/75221; http://www.uniprot.org/uniprot/Q99KK7
		EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_symbol'] = 'Dpp3'
	} else if(EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_id'] == 'MGI:2443973') {			# http://www.informatics.jax.org/marker/MGI:2443973; https://www.ncbi.nlm.nih.gov/gene/231991; http://www.uniprot.org/uniprot/Q8K1L0
		EVpedia_Top100.Protein_Mouse.mgi_id.table[i,'mgi_symbol'] = 'Creb5'
	}
}

# --- Gene symbol re-mapping ---
EVpedia_Top100.Protein_Mouse.Gene_symbol.idx <- sapply(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'mgi_symbol'], FUN=function(Y) which(sapply(strsplit(EVpedia_Top100.Protein_Mouse.Gene_symbol,'\\; | '), FUN=function(X) Y %in% X)))
EVpedia_Top100.Protein_Mouse.Gene_symbol.idx.length <- sapply(EVpedia_Top100.Protein_Mouse.Gene_symbol.idx, FUN=function(Z) as.numeric(length(Z)))

EVpedia_Top100.Protein_Mouse.Gene_symbol.match_0.idx = which(EVpedia_Top100.Protein_Mouse.Gene_symbol.idx.length == 0)
EVpedia_Top100.Protein_Mouse.Gene_symbol.match_1.idx = which(EVpedia_Top100.Protein_Mouse.Gene_symbol.idx.length == 1)
EVpedia_Top100.Protein_Mouse.Gene_symbol.match_2.idx = which(EVpedia_Top100.Protein_Mouse.Gene_symbol.idx.length == 2)
EVpedia_Top100.Protein_Mouse.Gene_symbol.match_3.idx = which(EVpedia_Top100.Protein_Mouse.Gene_symbol.idx.length == 3)
EVpedia_Top100.Protein_Mouse.Gene_symbol.match_4.idx = which(EVpedia_Top100.Protein_Mouse.Gene_symbol.idx.length == 4)
EVpedia_Top100.Protein_Mouse.Gene_symbol.match_5.idx = which(EVpedia_Top100.Protein_Mouse.Gene_symbol.idx.length == 5)
EVpedia_Top100.Protein_Mouse.Gene_symbol.match_6.idx = which(EVpedia_Top100.Protein_Mouse.Gene_symbol.idx.length > 5)
#length(EVpedia_Top100.Protein_Mouse.Gene_symbol.match_0.idx)
#length(EVpedia_Top100.Protein_Mouse.Gene_symbol.match_1.idx)
#length(EVpedia_Top100.Protein_Mouse.Gene_symbol.match_2.idx)
#length(EVpedia_Top100.Protein_Mouse.Gene_symbol.match_3.idx)
#length(EVpedia_Top100.Protein_Mouse.Gene_symbol.match_4.idx)
#length(EVpedia_Top100.Protein_Mouse.Gene_symbol.match_5.idx)
#length(EVpedia_Top100.Protein_Mouse.Gene_symbol.match_6.idx)
#length(EVpedia_Top100.Protein_Mouse.Gene_symbol.match_0.idx) + length(EVpedia_Top100.Protein_Mouse.Gene_symbol.match_1.idx) + length(EVpedia_Top100.Protein_Mouse.Gene_symbol.match_2.idx) + length(EVpedia_Top100.Protein_Mouse.Gene_symbol.match_3.idx) + length(EVpedia_Top100.Protein_Mouse.Gene_symbol.match_4.idx) + length(EVpedia_Top100.Protein_Mouse.Gene_symbol.match_5.idx) + length(EVpedia_Top100.Protein_Mouse.Gene_symbol.match_6.idx)
if((length(EVpedia_Top100.Protein_Mouse.Gene_symbol.match_0.idx) + length(EVpedia_Top100.Protein_Mouse.Gene_symbol.match_1.idx) + length(EVpedia_Top100.Protein_Mouse.Gene_symbol.match_2.idx) + length(EVpedia_Top100.Protein_Mouse.Gene_symbol.match_3.idx) + length(EVpedia_Top100.Protein_Mouse.Gene_symbol.match_4.idx) + length(EVpedia_Top100.Protein_Mouse.Gene_symbol.match_5.idx) + length(EVpedia_Top100.Protein_Mouse.Gene_symbol.match_6.idx)) != nrow(EVpedia_Top100.Protein_Mouse.mgi_id.table))
{
	cat('[Warning] There is/are "Gene symbol" with more than 3 time identities by mapping between variable "strsplit(EVpedia_Top100.Protein_Mouse.Gene_symbol,\' \')" and variable "EVpedia_Top100.Protein_Mouse.mgi_id.table[,\'mgi_symbol\']".\n')
}

# --- remove incorrect gene mapping ---
if(length(EVpedia_Top100.Protein_Mouse.Gene_symbol.match_0.idx) > 0)
{
	#EVpedia_Top100.Protein_Mouse.mgi_id.table[EVpedia_Top100.Protein_Mouse.Gene_symbol.match_0.idx,]
	
	remove.idx <- c()
	for(i in 1:length(EVpedia_Top100.Protein_Mouse.Gene_symbol.match_0.idx))
	{
		EVpedia_Top100.Protein_Mouse.Gene_symbol.match_0.uniprotswissprot = EVpedia_Top100.Protein_Mouse.mgi_id.table[EVpedia_Top100.Protein_Mouse.Gene_symbol.match_0.idx[i],'uniprotswissprot']
		EVpedia_Top100.Protein_Mouse.Gene_symbol.match_0.entrezgene = EVpedia_Top100.Protein_Mouse.mgi_id.table[EVpedia_Top100.Protein_Mouse.Gene_symbol.match_0.idx[i],'entrezgene']
		EVpedia_Top100.Protein_Mouse.Gene_symbol.match_0.mgi_symbol = EVpedia_Top100.Protein_Mouse.mgi_id.table[EVpedia_Top100.Protein_Mouse.Gene_symbol.match_0.idx[i],'mgi_symbol']
		#print(EVpedia_Top100.Protein_Mouse.Gene_symbol.match_0.uniprotswissprot)
		#print(EVpedia_Top100.Protein_Mouse.Gene_symbol.match_0.entrezgene)
		#print(EVpedia_Top100.Protein_Mouse.Gene_symbol.match_0.mgi_symbol)
		
		if(EVpedia_Top100.Protein_Mouse.Gene_symbol.match_0.uniprotswissprot=='E9PZY8' & EVpedia_Top100.Protein_Mouse.Gene_symbol.match_0.entrezgene!='66185') {		# http://www.uniprot.org/uniprot/E9PZY8; Virma (https://www.ncbi.nlm.nih.gov/gene/66185)
			remove.idx <- c(remove.idx, EVpedia_Top100.Protein_Mouse.Gene_symbol.match_0.idx[i])
		} else if(EVpedia_Top100.Protein_Mouse.Gene_symbol.match_0.uniprotswissprot=='P09411' & EVpedia_Top100.Protein_Mouse.Gene_symbol.match_0.entrezgene!='18655') {	# http://www.uniprot.org/uniprot/P09411; Pgk1 (https://www.ncbi.nlm.nih.gov/gene/18655)
			remove.idx <- c(remove.idx, EVpedia_Top100.Protein_Mouse.Gene_symbol.match_0.idx[i])
		} else if(EVpedia_Top100.Protein_Mouse.Gene_symbol.match_0.uniprotswissprot=='P61358' & EVpedia_Top100.Protein_Mouse.Gene_symbol.match_0.entrezgene!='19942') {	# http://www.uniprot.org/uniprot/P61358; Rpl27 (https://www.ncbi.nlm.nih.gov/gene/19942); http://www.informatics.jax.org/marker/MGI:98036
			remove.idx <- c(remove.idx, EVpedia_Top100.Protein_Mouse.Gene_symbol.match_0.idx[i])
		} else if(EVpedia_Top100.Protein_Mouse.Gene_symbol.match_0.uniprotswissprot=='P62984' & EVpedia_Top100.Protein_Mouse.Gene_symbol.match_0.entrezgene!='22186') {	# http://www.uniprot.org/uniprot/P62984; Uba52 (https://www.ncbi.nlm.nih.gov/gene/22186)
			remove.idx <- c(remove.idx, EVpedia_Top100.Protein_Mouse.Gene_symbol.match_0.idx[i])
		} else if(EVpedia_Top100.Protein_Mouse.Gene_symbol.match_0.uniprotswissprot=='Q60737' & EVpedia_Top100.Protein_Mouse.Gene_symbol.match_0.entrezgene!='12995') {	# http://www.uniprot.org/uniprot/Q60737; Csnk2a1 (https://www.ncbi.nlm.nih.gov/gene/12995)
			remove.idx <- c(remove.idx, EVpedia_Top100.Protein_Mouse.Gene_symbol.match_0.idx[i])
		}
		
		if(EVpedia_Top100.Protein_Mouse.Gene_symbol.match_0.uniprotswissprot=='Q8VD58' & EVpedia_Top100.Protein_Mouse.Gene_symbol.match_0.mgi_symbol!='Evi2b') {		# http://www.uniprot.org/uniprot/Q8VD58; Evi2b (https://www.ncbi.nlm.nih.gov/gene/216984) / Evi2 (https://www.ncbi.nlm.nih.gov/gene/101488212)
			remove.idx <- c(remove.idx, EVpedia_Top100.Protein_Mouse.Gene_symbol.match_0.idx[i])
		}
	}
	
	if(length(remove.idx) > 0)
	{
		EVpedia_Top100.Protein_Mouse.mgi_id.table = EVpedia_Top100.Protein_Mouse.mgi_id.table[-remove.idx,]
	}
}

# --- merge multiple mgi_symbol, mgi_id, entrezgene, description of single uniprotswissprot ---
unique.uniprotswissprot.idx <- which(!duplicated(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']))
unique.uniprotswissprot <- EVpedia_Top100.Protein_Mouse.mgi_id.table[unique.uniprotswissprot.idx,'uniprotswissprot']
duplicated.uniprotswissprot.idx <- which(duplicated(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot']))
duplicated.uniprotswissprot <- EVpedia_Top100.Protein_Mouse.mgi_id.table[duplicated.uniprotswissprot.idx,'uniprotswissprot']
duplicated.uniprotswissprot.idx.pair <- unique(sapply(duplicated.uniprotswissprot[which(duplicated.uniprotswissprot!='')], FUN=function(X) which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'uniprotswissprot'] %in% X), simplify=FALSE))

duplicated.uniprotswissprot.idx.mgi_symbol <- sapply(duplicated.uniprotswissprot.idx.pair, FUN=function(X) paste(unique(unlist(strsplit(as.character(EVpedia_Top100.Protein_Mouse.mgi_id.table[X,'mgi_symbol'][EVpedia_Top100.Protein_Mouse.mgi_id.table[X,'mgi_symbol']!='']),'\\|'))),collapse="|"))
duplicated.uniprotswissprot.idx.mgi_id <- sapply(duplicated.uniprotswissprot.idx.pair, FUN=function(X) paste(unique(unlist(strsplit(as.character(EVpedia_Top100.Protein_Mouse.mgi_id.table[X,'mgi_id'][EVpedia_Top100.Protein_Mouse.mgi_id.table[X,'mgi_id']!='']),'\\|'))),collapse="|"))
duplicated.uniprotswissprot.idx.entrezgene <- sapply(duplicated.uniprotswissprot.idx.pair, FUN=function(X) paste(unique(unlist(strsplit(as.character(EVpedia_Top100.Protein_Mouse.mgi_id.table[X,'entrezgene'][EVpedia_Top100.Protein_Mouse.mgi_id.table[X,'entrezgene']!='']),'\\|'))),collapse="|"))
duplicated.uniprotswissprot.idx.uniprotswissprot <- sapply(duplicated.uniprotswissprot.idx.pair, FUN=function(X) paste(unique(unlist(strsplit(as.character(EVpedia_Top100.Protein_Mouse.mgi_id.table[X,'uniprotswissprot'][EVpedia_Top100.Protein_Mouse.mgi_id.table[X,'uniprotswissprot']!='']),'\\|'))),collapse="|"))
duplicated.uniprotswissprot.idx.description <- sapply(duplicated.uniprotswissprot.idx.pair, FUN=function(X) paste(unique(unlist(strsplit(as.character(EVpedia_Top100.Protein_Mouse.mgi_id.table[X,'description'][EVpedia_Top100.Protein_Mouse.mgi_id.table[X,'description']!='']),'\\|'))),collapse="|"))

if(length(duplicated.uniprotswissprot.idx.pair) > 0)
{
	for(i in 1:length(duplicated.uniprotswissprot.idx.pair))
	{
		EVpedia_Top100.Protein_Mouse.mgi_id.table[unlist(duplicated.uniprotswissprot.idx.pair[i]),'mgi_symbol'] = duplicated.uniprotswissprot.idx.mgi_symbol[i]
		EVpedia_Top100.Protein_Mouse.mgi_id.table[unlist(duplicated.uniprotswissprot.idx.pair[i]),'mgi_id'] = duplicated.uniprotswissprot.idx.mgi_id[i]
		EVpedia_Top100.Protein_Mouse.mgi_id.table[unlist(duplicated.uniprotswissprot.idx.pair[i]),'entrezgene'] = duplicated.uniprotswissprot.idx.entrezgene[i]
		EVpedia_Top100.Protein_Mouse.mgi_id.table[unlist(duplicated.uniprotswissprot.idx.pair[i]),'uniprotswissprot'] = duplicated.uniprotswissprot.idx.uniprotswissprot[i]
		EVpedia_Top100.Protein_Mouse.mgi_id.table[unlist(duplicated.uniprotswissprot.idx.pair[i]),'description'] = duplicated.uniprotswissprot.idx.description[i]
	}
}

# --- merge multiple uniprotswissprot of single mgi_id ---
unique.mgi_id.idx <- which(!duplicated(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'mgi_id']))
unique.mgi_id <- EVpedia_Top100.Protein_Mouse.mgi_id.table[unique.mgi_id.idx,'mgi_id']
duplicated.mgi_id.idx <- which(duplicated(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'mgi_id']))
duplicated.mgi_id <- EVpedia_Top100.Protein_Mouse.mgi_id.table[duplicated.mgi_id.idx,'mgi_id']
duplicated.mgi_id.idx.pair <- unique(sapply(duplicated.mgi_id[which(duplicated.mgi_id!='')], FUN=function(X) which(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'mgi_id'] %in% X), simplify=FALSE))

duplicated.mgi_id.idx.mgi_symbol <- sapply(duplicated.mgi_id.idx.pair, FUN=function(X) paste(unique(unlist(strsplit(as.character(EVpedia_Top100.Protein_Mouse.mgi_id.table[X,'mgi_symbol'][EVpedia_Top100.Protein_Mouse.mgi_id.table[X,'mgi_symbol']!='']),'\\|'))),collapse="|"))
duplicated.mgi_id.idx.mgi_id <- sapply(duplicated.mgi_id.idx.pair, FUN=function(X) paste(unique(unlist(strsplit(as.character(EVpedia_Top100.Protein_Mouse.mgi_id.table[X,'mgi_id'][EVpedia_Top100.Protein_Mouse.mgi_id.table[X,'mgi_id']!='']),'\\|'))),collapse="|"))
duplicated.mgi_id.idx.entrezgene <- sapply(duplicated.mgi_id.idx.pair, FUN=function(X) paste(unique(unlist(strsplit(as.character(EVpedia_Top100.Protein_Mouse.mgi_id.table[X,'entrezgene'][EVpedia_Top100.Protein_Mouse.mgi_id.table[X,'entrezgene']!='']),'\\|'))),collapse="|"))
duplicated.mgi_id.idx.uniprotswissprot <- sapply(duplicated.mgi_id.idx.pair, FUN=function(X) paste(unique(unlist(strsplit(as.character(EVpedia_Top100.Protein_Mouse.mgi_id.table[X,'uniprotswissprot'][EVpedia_Top100.Protein_Mouse.mgi_id.table[X,'uniprotswissprot']!='']),'\\|'))),collapse="|"))
duplicated.mgi_id.idx.description <- sapply(duplicated.mgi_id.idx.pair, FUN=function(X) paste(unique(unlist(strsplit(as.character(EVpedia_Top100.Protein_Mouse.mgi_id.table[X,'description'][EVpedia_Top100.Protein_Mouse.mgi_id.table[X,'description']!='']),'\\|'))),collapse="|"))

if(length(duplicated.mgi_id.idx.pair) > 0)
{
	for(i in 1:length(duplicated.mgi_id.idx.pair))
	{
		EVpedia_Top100.Protein_Mouse.mgi_id.table[unlist(duplicated.mgi_id.idx.pair[i]),'mgi_symbol'] = duplicated.mgi_id.idx.mgi_symbol[i]
		EVpedia_Top100.Protein_Mouse.mgi_id.table[unlist(duplicated.mgi_id.idx.pair[i]),'mgi_id'] = duplicated.mgi_id.idx.mgi_id[i]
		EVpedia_Top100.Protein_Mouse.mgi_id.table[unlist(duplicated.mgi_id.idx.pair[i]),'entrezgene'] = duplicated.mgi_id.idx.entrezgene[i]
		EVpedia_Top100.Protein_Mouse.mgi_id.table[unlist(duplicated.mgi_id.idx.pair[i]),'uniprotswissprot'] = duplicated.mgi_id.idx.uniprotswissprot[i]
		EVpedia_Top100.Protein_Mouse.mgi_id.table[unlist(duplicated.mgi_id.idx.pair[i]),'description'] = duplicated.mgi_id.idx.description[i]
	}
}

EVpedia_Top100.Protein_Mouse.mgi_id.table = unique(EVpedia_Top100.Protein_Mouse.mgi_id.table)
#dim(EVpedia_Top100.Protein_Mouse.mgi_id.table)

rownames(EVpedia_Top100.Protein_Mouse.mgi_id.table) = c(1:nrow(EVpedia_Top100.Protein_Mouse.mgi_id.table))
#EVpedia_Top100.Protein_Mouse.mgi_id.table


# === combine "id table" and "raw input data" ===================================================================================
EVpedia_Top100.Protein_Mouse.Gene_symbol.idx <- c()
for(i in 1:length(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'mgi_symbol']))
{
	name = EVpedia_Top100.Protein_Mouse.mgi_id.table[,'mgi_symbol'][[i]]
	value = paste(sort(unique(as.numeric(unlist(sapply(unlist(strsplit(EVpedia_Top100.Protein_Mouse.mgi_id.table[,'mgi_symbol'][[i]],'\\|')), FUN=function(Y) which(sapply(strsplit(EVpedia_Top100.Protein_Mouse.Gene_symbol,'\\; | '), FUN=function(X) Y %in% X))))))), collapse="|")
	names(value) = name
	
	EVpedia_Top100.Protein_Mouse.Gene_symbol.idx <- c(EVpedia_Top100.Protein_Mouse.Gene_symbol.idx, value)
}
#EVpedia_Top100.Protein_Mouse.Gene_symbol.idx

EVpedia_Top100.Protein_Mouse.Index.New = sapply(EVpedia_Top100.Protein_Mouse.Gene_symbol.idx, FUN=function(X) paste(EVpedia_Top100.Protein_Mouse.Index[as.numeric(unlist(strsplit(X,'\\|')))],collapse="|"))
EVpedia_Top100.Protein_Mouse.UniProt_accession.New = sapply(EVpedia_Top100.Protein_Mouse.Gene_symbol.idx, FUN=function(X) paste(EVpedia_Top100.Protein_Mouse.UniProt_accession[as.numeric(unlist(strsplit(X,'\\|')))],collapse="|"))
EVpedia_Top100.Protein_Mouse.UniProt_name.New = sapply(EVpedia_Top100.Protein_Mouse.Gene_symbol.idx, FUN=function(X) paste(EVpedia_Top100.Protein_Mouse.UniProt_name[as.numeric(unlist(strsplit(X,'\\|')))],collapse="|"))
EVpedia_Top100.Protein_Mouse.Status.New = sapply(EVpedia_Top100.Protein_Mouse.Gene_symbol.idx, FUN=function(X) paste(EVpedia_Top100.Protein_Mouse.Status[as.numeric(unlist(strsplit(X,'\\|')))],collapse="|"))
EVpedia_Top100.Protein_Mouse.Protein_name.New = sapply(EVpedia_Top100.Protein_Mouse.Gene_symbol.idx, FUN=function(X) paste(EVpedia_Top100.Protein_Mouse.Protein_name[as.numeric(unlist(strsplit(X,'\\|')))],collapse="|"))
EVpedia_Top100.Protein_Mouse.Gene_symbol.New = sapply(EVpedia_Top100.Protein_Mouse.Gene_symbol.idx, FUN=function(X) paste(EVpedia_Top100.Protein_Mouse.Gene_symbol[as.numeric(unlist(strsplit(X,'\\|')))],collapse="|"))
EVpedia_Top100.Protein_Mouse.Species.New = unique(unlist(sapply(EVpedia_Top100.Protein_Mouse.Gene_symbol.idx, FUN=function(X) unique(unlist(EVpedia_Top100.Protein_Mouse.Species[as.numeric(unlist(strsplit(X,'\\|')))])))))
EVpedia_Top100.Protein_Mouse.Identification_number.New = sapply(EVpedia_Top100.Protein_Mouse.Gene_symbol.idx, FUN=function(X) paste(EVpedia_Top100.Protein_Mouse.Identification_number[as.numeric(unlist(strsplit(X,'\\|')))],collapse="|"))
EVpedia_Top100.Protein_Mouse.Identification_number.New.Sum = sapply(EVpedia_Top100.Protein_Mouse.Identification_number.New, FUN=function(X) sum(as.numeric(unlist(strsplit(X,'\\|')))))

EVpedia_Top100.Protein_Mouse.New = cbind(EVpedia_Top100.Protein_Mouse.Index.New, EVpedia_Top100.Protein_Mouse.UniProt_accession.New, EVpedia_Top100.Protein_Mouse.UniProt_name.New, EVpedia_Top100.Protein_Mouse.Status.New, EVpedia_Top100.Protein_Mouse.Protein_name.New, EVpedia_Top100.Protein_Mouse.Gene_symbol.New, EVpedia_Top100.Protein_Mouse.Species.New, EVpedia_Top100.Protein_Mouse.Identification_number.New, EVpedia_Top100.Protein_Mouse.Identification_number.New.Sum)
colnames(EVpedia_Top100.Protein_Mouse.New) = c(gsub('\\.','_',colnames(EVpedia_Top100.Protein_Mouse)), paste(gsub('\\.','_',colnames(EVpedia_Top100.Protein_Mouse))[length(gsub('\\.','_',colnames(EVpedia_Top100.Protein_Mouse)))],'_New',sep=""))
#EVpedia_Top100.Protein_Mouse.New

EVpedia_Top100.Protein_Mouse.mgi_id.table.New = cbind(EVpedia_Top100.Protein_Mouse.mgi_id.table, EVpedia_Top100.Protein_Mouse.New)
for(i in 1:ncol(EVpedia_Top100.Protein_Mouse.mgi_id.table.New))
{
	EVpedia_Top100.Protein_Mouse.mgi_id.table.New[,i] = as.character(EVpedia_Top100.Protein_Mouse.mgi_id.table.New[,i])
}
#EVpedia_Top100.Protein_Mouse.mgi_id.table.New
#dim(EVpedia_Top100.Protein_Mouse.mgi_id.table.New)
#write.csv(EVpedia_Top100.Protein_Mouse.mgi_id.table.New, paste('./Data/',Sample,'/10_datasetOverlap/Database/mmu/EVpedia_Top.Protein_UniProt2MGI.table_IdentificationNumber.csv',sep=""), quote=TRUE, row.names=TRUE)

# --- check for duplicated Index ---
query.Index = as.numeric(EVpedia_Top100.Protein_Mouse.Index)
subject.Index = as.numeric(unlist(strsplit(as.character(EVpedia_Top100.Protein_Mouse.mgi_id.table.New[,'Index']),'\\|')))
subject.Index.duplicated = subject.Index[which(duplicated(subject.Index))]

if(length(subject.Index.duplicated) > 0)
{
	for(i in 1:length(subject.Index.duplicated))
	{
		subject.Index.duplicated.idx = which(sapply(strsplit(as.character(EVpedia_Top100.Protein_Mouse.mgi_id.table.New[,'Index']),'\\|'), FUN=function(X) which(X %in% subject.Index.duplicated[i])) > 0)
		#EVpedia_Top100.Protein_Mouse.mgi_id.table.New[subject.Index.duplicated.idx,]
		
		subject.Index.duplicated.idx.queryIndex = as.numeric(unique(unlist(strsplit(as.character(EVpedia_Top100.Protein_Mouse.mgi_id.table.New[subject.Index.duplicated.idx,'Index']),'\\|'))))
		
		EVpedia_Top100.Protein_Mouse.mgi_id.table.New[subject.Index.duplicated.idx,'mgi_symbol'] = paste(unique(unlist(strsplit(as.character(EVpedia_Top100.Protein_Mouse.mgi_id.table.New[subject.Index.duplicated.idx,'mgi_symbol']),'\\|'))),collapse="|")
		EVpedia_Top100.Protein_Mouse.mgi_id.table.New[subject.Index.duplicated.idx,'mgi_id'] = paste(unique(unlist(strsplit(as.character(EVpedia_Top100.Protein_Mouse.mgi_id.table.New[subject.Index.duplicated.idx,'mgi_id']),'\\|'))),collapse="|")
		EVpedia_Top100.Protein_Mouse.mgi_id.table.New[subject.Index.duplicated.idx,'entrezgene'] = paste(unique(unlist(strsplit(as.character(EVpedia_Top100.Protein_Mouse.mgi_id.table.New[subject.Index.duplicated.idx,'entrezgene']),'\\|'))),collapse="|")
		EVpedia_Top100.Protein_Mouse.mgi_id.table.New[subject.Index.duplicated.idx,'uniprotswissprot'] = paste(unique(unlist(strsplit(as.character(EVpedia_Top100.Protein_Mouse.mgi_id.table.New[subject.Index.duplicated.idx,'uniprotswissprot']),'\\|'))),collapse="|")
		EVpedia_Top100.Protein_Mouse.mgi_id.table.New[subject.Index.duplicated.idx,'description'] = paste(unique(unlist(strsplit(as.character(EVpedia_Top100.Protein_Mouse.mgi_id.table.New[subject.Index.duplicated.idx,'description']),'\\|'))),collapse="|")
		
		EVpedia_Top100.Protein_Mouse.mgi_id.table.New[subject.Index.duplicated.idx,'Index'] = paste(subject.Index.duplicated.idx.queryIndex,collapse="|")
		
		EVpedia_Top100.Protein_Mouse.Index.idx = which(EVpedia_Top100.Protein_Mouse.Index %in% subject.Index.duplicated.idx.queryIndex)
		EVpedia_Top100.Protein_Mouse.mgi_id.table.New[subject.Index.duplicated.idx,'UniProt_accession'] = paste(EVpedia_Top100.Protein_Mouse.UniProt_accession[EVpedia_Top100.Protein_Mouse.Index.idx],collapse="|")
		EVpedia_Top100.Protein_Mouse.mgi_id.table.New[subject.Index.duplicated.idx,'UniProt_name'] = paste(EVpedia_Top100.Protein_Mouse.UniProt_name[EVpedia_Top100.Protein_Mouse.Index.idx],collapse="|")
		EVpedia_Top100.Protein_Mouse.mgi_id.table.New[subject.Index.duplicated.idx,'Status'] = paste(EVpedia_Top100.Protein_Mouse.Status[EVpedia_Top100.Protein_Mouse.Index.idx],collapse="|")
		EVpedia_Top100.Protein_Mouse.mgi_id.table.New[subject.Index.duplicated.idx,'Protein_name'] = paste(EVpedia_Top100.Protein_Mouse.Protein_name[EVpedia_Top100.Protein_Mouse.Index.idx],collapse="|")
		EVpedia_Top100.Protein_Mouse.mgi_id.table.New[subject.Index.duplicated.idx,'Gene_symbol'] = paste(EVpedia_Top100.Protein_Mouse.Gene_symbol[EVpedia_Top100.Protein_Mouse.Index.idx],collapse="|")
		EVpedia_Top100.Protein_Mouse.mgi_id.table.New[subject.Index.duplicated.idx,'Species'] = unique(EVpedia_Top100.Protein_Mouse.Species[EVpedia_Top100.Protein_Mouse.Index.idx])
		EVpedia_Top100.Protein_Mouse.mgi_id.table.New[subject.Index.duplicated.idx,'Identification_number'] = paste(as.numeric(EVpedia_Top100.Protein_Mouse.Identification_number[EVpedia_Top100.Protein_Mouse.Index.idx]),collapse="|")
		EVpedia_Top100.Protein_Mouse.mgi_id.table.New[subject.Index.duplicated.idx,'Identification_number_New'] = sum(as.numeric(EVpedia_Top100.Protein_Mouse.Identification_number[EVpedia_Top100.Protein_Mouse.Index.idx]))
		
		#EVpedia_Top100.Protein_Mouse.mgi_id.table.New[subject.Index.duplicated.idx,]
	}
}

EVpedia_Top100.Protein_Mouse.mgi_id.table.New = unique(EVpedia_Top100.Protein_Mouse.mgi_id.table.New)
#EVpedia_Top100.Protein_Mouse.mgi_id.table.New
#dim(EVpedia_Top100.Protein_Mouse.mgi_id.table.New)
#write.csv(EVpedia_Top100.Protein_Mouse.mgi_id.table.New, paste('./Data/',Sample,'/10_datasetOverlap/Database/mmu/EVpedia_Top.Protein_UniProt2MGI.table_IdentificationNumber.csv',sep=""), quote=TRUE, row.names=TRUE)

# --- re-check for duplicated Index again ---
query.Index = as.numeric(EVpedia_Top100.Protein_Mouse.Index)
subject.Index = as.numeric(unlist(strsplit(as.character(EVpedia_Top100.Protein_Mouse.mgi_id.table.New[,'Index']),'\\|')))
subject.Index.duplicated = subject.Index[which(duplicated(subject.Index))]

if(length(subject.Index.duplicated) > 0)
{
	cat(paste('Please check for duplicated Index of variable "EVpedia_Top100.Protein_Mouse.mgi_id.table.New[,\'Index\']". It still has duplicated Index.\n',sep=""))
	print(subject.Index.duplicated)
}

# --- check for uniqueness of each column ---
checked.duplicated.idx <- c()
for(i in 1:ncol(EVpedia_Top100.Protein_Mouse.mgi_id.table.New))
{
	if(! colnames(EVpedia_Top100.Protein_Mouse.mgi_id.table.New)[i] %in% c('Status','Species','Identification_number','Identification_number_New'))
	{
		column.row_duplicated_idx = which(duplicated(EVpedia_Top100.Protein_Mouse.mgi_id.table.New[,colnames(EVpedia_Top100.Protein_Mouse.mgi_id.table.New)[i]]))
		#print(column.row_duplicated_idx)
		
		if(length(column.row_duplicated_idx) > 0)
		{
			column.row_duplicated_value = EVpedia_Top100.Protein_Mouse.mgi_id.table.New[column.row_duplicated_idx,colnames(EVpedia_Top100.Protein_Mouse.mgi_id.table.New)[i]]
			#print(column.row_duplicated_value)
			
			removed.checked_duplicated_name = c(NA,'NA','',
												'Uncharacterized protein','Uncharacterized protein|Uncharacterized protein','Uncharacterized protein (Fragment)','Putative uncharacterized protein',
												'Olfactory receptor',
												'Histone H2A',
												'Spectrin beta chain','Spectrin beta chain|Spectrin beta chain',
												'60S ribosomal protein L29',
												'Phospholipid-transporting ATPase (EC 3.6.3.1)',
												'Diacylglycerol kinase (DAG kinase) (EC 2.7.1.107)',
												'ribosomal protein S6 kinase, polypeptide 2')
			
			column.row_duplicated_value.idx = which(!column.row_duplicated_value %in% removed.checked_duplicated_name)
			#print(column.row_duplicated_value.idx)
			
			column.row_duplicated_idx = column.row_duplicated_idx[column.row_duplicated_value.idx]
			column.row_duplicated_value = column.row_duplicated_value[column.row_duplicated_value.idx]
			#print(column.row_duplicated_idx)
			#print(column.row_duplicated_value)
			
			if(length(column.row_duplicated_value) > 0)
			{
				cat(paste('Please check variable "EVpedia_Top100.Protein_Mouse.mgi_id.table.New[,\'',colnames(EVpedia_Top100.Protein_Mouse.mgi_id.table.New)[i],'\']".',' Column "',colnames(EVpedia_Top100.Protein_Mouse.mgi_id.table.New)[i],'" contains duplicated items.\n',sep=""))
				checked.duplicated.idx <- c(checked.duplicated.idx, column.row_duplicated_idx)
			}
		}
	}
}

if(length(checked.duplicated.idx) > 0)
{
	checked.duplicated.idx = sort(unique(checked.duplicated.idx))
	print(EVpedia_Top100.Protein_Mouse.mgi_id.table.New[checked.duplicated.idx,])
}

rownames(EVpedia_Top100.Protein_Mouse.mgi_id.table.New) = c(1:nrow(EVpedia_Top100.Protein_Mouse.mgi_id.table.New))
#EVpedia_Top100.Protein_Mouse.mgi_id.table.New

EVpedia_Top100.Protein_Mouse.mgi_id.table.New = EVpedia_Top100.Protein_Mouse.mgi_id.table.New[mixedorder(EVpedia_Top100.Protein_Mouse.mgi_id.table.New[,'Identification_number_New'],decreasing=TRUE),]
rownames(EVpedia_Top100.Protein_Mouse.mgi_id.table.New) = c(1:nrow(EVpedia_Top100.Protein_Mouse.mgi_id.table.New))
#EVpedia_Top100.Protein_Mouse.mgi_id.table.New

#write.csv(EVpedia_Top100.Protein_Mouse.mgi_id.table.New, paste('./Data/',Sample,'/10_datasetOverlap/Database/mmu/EVpedia_Top.Protein_UniProt2MGI.table_IdentificationNumber.csv',sep=""), quote=TRUE, row.names=TRUE)

# === uniprotswissprot not found by R package "biomaRt" and "BridgeDbR" ========================================================
query.Index = as.numeric(EVpedia_Top100.Protein_Mouse.Index)
#length(query.Index)
query.Index.sorted = sort(unique(query.Index))

subject.Index = as.numeric(unlist(strsplit(as.character(EVpedia_Top100.Protein_Mouse.mgi_id.table.New[,'Index']),'\\|')))
#length(subject.Index)
subject.Index.sorted = sort(unique(subject.Index))

unmapping.Index = which(! query.Index.sorted %in% subject.Index.sorted)
#length(unmapping.Index)
unmapping.Index.idx = which(EVpedia_Top100.Protein_Mouse.Index %in% unmapping.Index)
unmapping.Index.table = cbind('', '', '', '', '', EVpedia_Top100.Protein_Mouse.Index[unmapping.Index.idx], EVpedia_Top100.Protein_Mouse.UniProt_accession[unmapping.Index.idx], EVpedia_Top100.Protein_Mouse.UniProt_name[unmapping.Index.idx], EVpedia_Top100.Protein_Mouse.Status[unmapping.Index.idx], EVpedia_Top100.Protein_Mouse.Protein_name[unmapping.Index.idx], EVpedia_Top100.Protein_Mouse.Gene_symbol[unmapping.Index.idx], EVpedia_Top100.Protein_Mouse.Species[unmapping.Index.idx], EVpedia_Top100.Protein_Mouse.Identification_number[unmapping.Index.idx], EVpedia_Top100.Protein_Mouse.Identification_number[unmapping.Index.idx])
colnames(unmapping.Index.table) = colnames(EVpedia_Top100.Protein_Mouse.mgi_id.table.New)
rownames(unmapping.Index.table) = c(1:nrow(unmapping.Index.table))
#unmapping.Index.table

EVpedia_Top100.Protein_Mouse.mgi_id.table.New.withUnmappingProtein = rbind(EVpedia_Top100.Protein_Mouse.mgi_id.table.New, unmapping.Index.table)
rownames(EVpedia_Top100.Protein_Mouse.mgi_id.table.New.withUnmappingProtein) = c(1:nrow(EVpedia_Top100.Protein_Mouse.mgi_id.table.New.withUnmappingProtein))
#EVpedia_Top100.Protein_Mouse.mgi_id.table.New.withUnmappingProtein

#write.csv(EVpedia_Top100.Protein_Mouse.mgi_id.table.New.withUnmappingProtein, paste('./Data/',Sample,'/10_datasetOverlap/Database/mmu/EVpedia_Top.Protein_UniProt2MGI.table_IdentificationNumber_withUnmappingProtein.csv',sep=""), quote=TRUE, row.names=TRUE)

# --- write xlsx file by openxlsx package ---------------------------------------------------------
# Importing a big xlsx file into R? https://stackoverflow.com/questions/19147884/importing-a-big-xlsx-file-into-r/43118530#43118530
# Building R for Windows: https://cran.r-project.org/bin/windows/Rtools/
#Sys.setenv("R_ZIPCMD" = "C:/Rtools/bin/zip.exe")	# path to zip.exe

wb <- openxlsx::createWorkbook(paste('./Data/',Sample,'/10_datasetOverlap/Database/mmu/EVpedia_Top.Protein_UniProt2MGI.table_IdentificationNumber.xlsx',sep=""))
modifyBaseFont(wb, fontSize=12, fontColour="black", fontName="Times New Roman")

for(s in 1:2)
{
	if(s == 1) {
		sheet = 'EVpedia_Top (mapped)'
		sheetData = EVpedia_Top100.Protein_Mouse.mgi_id.table.New
	} else if(s == 2) {
		sheet = 'EVpedia_Top (unmapped)'
		sheetData = unmapping.Index.table
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
			
			right_align.idx = which(colnames(sheetData) %in% c('Identification_number','Identification_number_New'))
			if(length(right_align.idx) > 0) {
				style <- createStyle(halign="right", valign="center")
				addStyle(wb, sheet, style, rows=c(2:(1+nrow(sheetData))), cols=1+right_align.idx, gridExpand=TRUE, stack=TRUE)
			}
		}
	}
}

if(length(wb$sheet_names) > 0) {
	openxlsx::saveWorkbook(wb, file=paste('./Data/',Sample,'/10_datasetOverlap/Database/mmu/EVpedia_Top.Protein_UniProt2MGI.table_IdentificationNumber.xlsx',sep=""), overwrite=TRUE)
}

