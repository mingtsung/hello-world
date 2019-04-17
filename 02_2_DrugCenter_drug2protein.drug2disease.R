# --- drug center ---------------------------------------------------------------------------------
drug.center <- c()
for(m in 1:length(drug))
{
	# /// drug2protein ///
	if(DRUG == 1) {
		protein.subset <- sort(unique(as.character(ProteinDrugLink[as.character(ProteinDrugLink[,3]) %in% drug[m],1])))
		if(exists('RemovedProteinDrugLink')) {	# add removed protein(s)
			protein.subset <- c(protein.subset, sort(unique(as.character(RemovedProteinDrugLink[as.character(RemovedProteinDrugLink[,3]) %in% drug[m],1]))))
		}
		protein.subset <- sort(unique(protein.subset))
		#protein.subset <- protein.subset[! protein.subset %in% remove.protein]	# remove specific protein(s)
	} else {
		protein.subset <- c()
	}
	
	protein.list <- ''
	if(length(protein.subset) > 0) {
		for(k in 1:length(protein.subset))
		{
			protein_id = protein.subset[k]
			
			uniprotId.idx = which(sapply(trimws(human_geneID_DisGeNET$c2.uniprotId), FUN=function(X) protein_id %in% unlist(strsplit(unlist(X),';'))))
			
			gene.geneId <- unique(as.character(trimws(human_geneID_DisGeNET[uniprotId.idx,,drop=FALSE]$c2.geneId)))
			gene.symbol <- unique(as.character(trimws(human_geneID_DisGeNET[uniprotId.idx,,drop=FALSE]$c2.symbol)))
			gene.description <- unique(as.character(trimws(human_geneID_DisGeNET[uniprotId.idx,,drop=FALSE]$c2.description)))
			gene.pantherName <- unique(as.character(trimws(human_geneID_DisGeNET[uniprotId.idx,,drop=FALSE]$c2.pantherName)))
			
			if(length(gene.geneId) > 1)				{	cat('Please check gene geneId. "',protein_id,'" in human_geneID_DisGeNET_',Func[i],'.txt matches to multiple gene geneId','.\n', sep="")			}
			if(length(gene.symbol) > 1)				{	cat('Please check gene symbol. "',protein_id,'" in human_geneID_DisGeNET_',Func[i],'.txt matches to multiple gene symbol','.\n', sep="")			}
			if(length(gene.description) > 1)		{	cat('Please check gene description. "',protein_id,'" in human_geneID_DisGeNET_',Func[i],'.txt matches to multiple gene description','.\n', sep="")	}
			if(length(gene.pantherName) > 1)		{	cat('Please check gene pantherName. "',protein_id,'" in human_geneID_DisGeNET_',Func[i],'.txt matches to multiple gene pantherName','.\n', sep="")	}
			
			protein.genedescription = gene.description
			if(length(protein.genedescription) > 1) {	print(paste('Please check protein name. ',protein_id,' matches to multiple protein names.',sep=""))	}
			
			if(length(protein.genedescription) == 1){
				if(protein.genedescription != '') {
					if(k == 1) {
						protein.list <- paste(protein_id,protein.genedescription,sep='~')
					} else {
						protein.list <- paste(protein.list,'|',paste(protein_id,protein.genedescription,sep='~'),sep="")
					}
				} else {
					if(k == 1) {
						protein.list <- paste(protein_id,'N/A',sep='~')
					} else {
						protein.list <- paste(protein.list,'|',paste(protein_id,'N/A',sep='~'),sep="")
					}
				}
			} else if(length(protein.genedescription) == 0){
				if(k == 1) {
					protein.list <- paste(protein_id,'N/A',sep='~')
				} else {
					protein.list <- paste(protein.list,'|',paste(protein_id,'N/A',sep='~'),sep="")
				}
			}
			
			if(k == 1) {
				disease.subset <- sort(unique(as.character(ProteinDiseaseLink[as.character(ProteinDiseaseLink[,1]) %in% protein_id,3])))
			} else {
				disease.subset <- c(disease.subset, sort(unique(as.character(ProteinDiseaseLink[as.character(ProteinDiseaseLink[,1]) %in% protein_id,3]))))
			}
		}
		#print(protein.list)
		
		# /// drug2disease ///
		if(exists('RemovedProteinDiseaseLink')) {	# add removed protein(s)
			disease.subset <- c(disease.subset, sort(unique(as.character(RemovedProteinDiseaseLink[as.character(RemovedProteinDiseaseLink[,1]) %in% protein_id,3]))))
		}
		disease.subset <- sort(unique(disease.subset))
		
		disease.list <- ''
		if(length(disease.subset) > 0) {
			if(disease.level[j] == 'Name') {
				for(k in 1:length(disease.subset))
				{
					disease.diseaseId <- unique(as.character(trimws(human_geneID_DisGeNET[trimws(human_geneID_DisGeNET$c1.name) %in% disease.subset[k],,drop=FALSE]$c1.diseaseId)))
					if(length(disease.diseaseId) > 1)	{	cat('Please check disease name. "',disease.subset[k],'" in human_geneID_DisGeNET_',Func[i],'.txt matches to multiple disease names','.\n', sep="")	}
					
					if(length(disease.diseaseId) == 1){
						if(disease.diseaseId != '') {
							if(k == 1) {
								disease.list <- paste(disease.diseaseId,disease.subset[k],sep='~')
							} else {
								disease.list <- paste(disease.list,'|',paste(disease.diseaseId,disease.subset[k],sep='~'),sep="")
							}
						} else {
							if(k == 1) {
								#disease.list <- paste('N/A',disease.subset[k],sep='~')
								disease.list <- disease.subset[k]
							} else {
								#disease.list <- paste(disease.list,'|',paste('N/A',disease.subset[k],sep='~'),sep="")
								disease.list <- paste(disease.list,'|',disease.subset[k],sep="")
							}
						}
					} else if(length(disease.diseaseId) == 0){
						if(k == 1) {
							#disease.list <- paste('N/A',disease.subset[k],sep='~')
							disease.list <- disease.subset[k]
						} else {
							#disease.list <- paste(disease.list,'|',paste('N/A',disease.subset[k],sep='~'),sep="")
							disease.list <- paste(disease.list,'|',disease.subset[k],sep="")
						}
					}
				}
				#print(disease.list)
			} else {
				for(k in 1:length(disease.subset))
				{
					if(k == 1) {
						disease.list <- disease.subset[k]
					} else {
						disease.list <- paste(disease.list,'|',disease.subset[k],sep="")
					}
				}
				#print(disease.list)
			}
		} else {
			disease.list <- ''
		}
		
		drug.center <- rbind(drug.center, cbind(drug[m], protein.list, length(protein.subset), disease.list, length(disease.subset)))
	} else {
		protein.list <- ''
	}
	
	#drug.center <- rbind(drug.center, cbind(drug[m], protein.list, length(protein.subset), disease.list, length(disease.subset)))
}

# === insert drug name and other information (i.e. drug type, drug group) ===============
drug.name <- c()
drug.type <- c()
drug.group <- c()
for(k in 1:length(drug.center[,1]))
{
	drug.drugid = drug.center[,1][k]
	drug_links.idx = which(drug_links.DrugBankID %in% drug.drugid)
	
	if(length(drug_links.idx) == 0) {
		drugbank_vocabulary.IDs.idx = which(sapply(drugbank_vocabulary.IDs, FUN=function(X) drug.drugid %in% trimws(unlist(strsplit(unlist(X),"\\|")))))
		if(length(drugbank_vocabulary.IDs.idx) > 0) {
			drug.drugid.official = trimws(unlist(strsplit(drugbank_vocabulary.IDs[drugbank_vocabulary.IDs.idx],"\\|")))[1]
			drug_links.idx = which(drug_links.DrugBankID %in% drug.drugid.official)
		}
	}
	
	drug.drugname = as.character(drug_links.Name[drug_links.idx])
	drug.drugtype = as.character(drug_links.DrugType[drug_links.idx])
	
	if(length(drug_links.idx) > 0) {
		if(drug_links[drug_links.idx,'Name'] != drug.drugname) {
			cat('[Warning] "',drug.drugid,'" has drug name',' (',drug.drugname,')',' in "drug.proteinorder"',', but, drug name',' (',drug_links[drug_links.idx,'Name'],')',' in "./Database/DrugBank/External_Links/External_Drug_Links/drugbank_all_drug_links.csv"','.\n', sep="")
		}
		if(drug_links[drug_links.idx,'Drug Type'] != drug.drugtype) {
			cat('[Warning] "',drug.drugid,'" has drug type',' (',drug.drugtype,')',' in "drug.proteinorder"',', but, drug type',' (',drug_links[drug_links.idx,'Drug Type'],')',' in "./Database/DrugBank/External_Links/External_Drug_Links/drugbank_all_drug_links.csv"','.\n', sep="")
		}
		drug_links.info = drug_links[drug_links.idx,,drop=FALSE]
		drug_links.info = data.frame(lapply(drug_links.info, trimws), check.names=FALSE, stringsAsFactors=FALSE)
		if(length(drug_links.idx) != 1) {
			cat('[Note] "',drug.drugid,'" has multiple records',' (',length(drug_links.idx),') in "./Database/DrugBank/External_Links/External_Drug_Links/drugbank_all_drug_links.csv"','.\n', sep="")
		}
	} else {
		drug_links.info = matrix(data='', nrow=1, ncol=ncol(drug_links), byrow=TRUE)
		colnames(drug_links.info) = colnames(drug_links)
		drug_links.info[,'DrugBank ID'] = drug.center[,1][k]
		drug_links.info[,'Name'] = ''
	}
	
	Drug_Group = c('Approved','Experimental','Illicit','Investigational','Nutraceutical','Withdrawn')
	Drug_Group.drug_links <- c()
	for(dg in 1:length(Drug_Group))
	{
		assign(paste(tolower(Drug_Group[dg]),'_','drug_links.idx',sep=""), which(eval(parse(text=paste(tolower(Drug_Group[dg]),'_','drug_links',sep="")))[,'DrugBank ID'] %in% drug.drugid))
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
	
	drug.name <- c(drug.name, drug_links.info[,'Name'])
	drug.type <- c(drug.type, drug_links.info[,'Drug Type'])
	drug.group <- c(drug.group, Drug_Group.drug_links)
}
#print(drug.name)
#print(drug.type)
#print(length(drug.name))
#print(length(drug.type))
drug.center <- cbind(drug.center[,1], drug.name, drug.type, drug.group, drug.center[,-1,drop=FALSE])

if(disease.level[j] == 'DisGeNET') {
	disease.HeaderTerm = 'disease'
} else if(disease.level[j] == 'Name') {
	disease.HeaderTerm = 'disease'
} else if(disease.level[j] == 'DisClassName') {
	disease.HeaderTerm = 'MeSH class'
} else if(disease.level[j] == 'PantherName') {
	disease.HeaderTerm = 'protein class'
}
colnames(drug.center) <- c('drug','drug name','drug type','drug group','protein','protein number',disease.HeaderTerm,paste(disease.HeaderTerm,' ','number',sep=""))
#print(drug.center)

if(DRUG == 1) {
	if(sum(as.numeric(drug.center[,4+2])) != nrow(unique(rbind(ProteinDrugLink, RemovedProteinDrugLink))))
	{
		cat('[Warning] Please check the number of protein. (DrugCenter_drug2protein.drug2disease)','\n')
	}
}

if(dim(drug.center)[1] == 1) {
	drug.proteinorder <- drug.center
	drug.diseaseorder <- drug.center
} else {
	drug.proteinorder <- drug.center[order(as.numeric(drug.center[,4+2]), as.numeric(drug.center[,6+2]), decreasing=TRUE),,drop=FALSE]
	drug.diseaseorder <- drug.center[order(as.numeric(drug.center[,6+2]), as.numeric(drug.center[,4+2]), decreasing=TRUE),,drop=FALSE]
}
#print(drug.proteinorder)
#print(drug.diseaseorder)
#write.table(drug.proteinorder, paste('./Data/',Sample,'/07_disease2protein2drug/subset/n_',r,'/',folder[f],'/','human','/',Func[i],'/',disease.level[j],'/',Sample,'_exo_','human','_geneID','_DisGeNET','_',Func[i],disease.level[j],'_','drug','_','protein','order','.csv',sep=""), quote=TRUE, sep=",", row.names=FALSE, col.names=TRUE)
#write.table(drug.diseaseorder, paste('./Data/',Sample,'/07_disease2protein2drug/subset/n_',r,'/',folder[f],'/','human','/',Func[i],'/',disease.level[j],'/',Sample,'_exo_','human','_geneID','_DisGeNET','_',Func[i],disease.level[j],'_','drug','_',gsub(' ','',disease.HeaderTerm),'order','.csv',sep=""), quote=TRUE, sep=",", row.names=FALSE, col.names=TRUE)

if(DRUG == 1) {
	assign(paste(Func[i],'.',disease.level[j],'.','drug.proteinorder',sep=""), drug.proteinorder)
	assign(paste(Func[i],'.',disease.level[j],'.','drug.diseaseorder',sep=""), drug.diseaseorder)
} else if(DRUG == 0) {
	
} else if(DRUG == -1) {
	cat('[Note] No relation files in "','./Data/',Sample,'/07_disease2protein2drug/subset/n_',r,'/',folder[f],'/','human','/',Func[i],'/',disease.level[j],'".\n', sep="")
}


# --- submatrix -----------------------------------------------------------------------------------
for(s in 1:2)
{
	if(s == 1) {
		drug.order = drug.proteinorder
		assign(paste(Func[i],'.',disease.level[j],'.','drug.proteinorder','.','A2B.submatrix.new.all',sep=""), c())
		assign(paste(Func[i],'.',disease.level[j],'.','drug.proteinorder','.','A2C.submatrix.new.all',sep=""), c())
	} else if(s == 2) {
		drug.order = drug.diseaseorder
		assign(paste(Func[i],'.',disease.level[j],'.','drug.diseaseorder','.','A2B.submatrix.new.all',sep=""), c())
		assign(paste(Func[i],'.',disease.level[j],'.','drug.diseaseorder','.','A2C.submatrix.new.all',sep=""), c())
	}
	
	A2B.submatrix.new.all <- c()
	A2C.submatrix.new.all <- c()
	
	for(k in 1:dim(drug.order)[1])
	{
		if(as.numeric(drug.order[k,'protein number']) > 0)
		{
			# --- protein -------------------------------------------------------------------------------------
			UniProt_id2name = unlist(strsplit(drug.order[k,'protein'],'\\|'))
			UniProt_id = sapply(UniProt_id2name, FUN=function(X) substr(unlist(X), 0, which(unlist(strsplit(unlist(X),''))=='~')[1]-1))
			UniProt_name = sapply(UniProt_id2name, FUN=function(X) substr(unlist(X), which(unlist(strsplit(unlist(X),''))=='~')[1]+1, nchar(unlist(X))))
			
			A2B.submatrix <- c()
			for(d in 1:length(UniProt_id))
			{
				protein_id = UniProt_id[d]
				protein_name = UniProt_name[d]
				
				uniprotId.idx = which(sapply(trimws(human_geneID_DisGeNET$c2.uniprotId), FUN=function(X) protein_id %in% unlist(strsplit(unlist(X),';'))))
				
				gene.geneId <- unique(as.character(trimws(human_geneID_DisGeNET[uniprotId.idx,,drop=FALSE]$c2.geneId)))
				gene.symbol <- unique(as.character(trimws(human_geneID_DisGeNET[uniprotId.idx,,drop=FALSE]$c2.symbol)))
				gene.description <- unique(as.character(trimws(human_geneID_DisGeNET[uniprotId.idx,,drop=FALSE]$c2.description)))
				gene.pantherName <- unique(as.character(trimws(human_geneID_DisGeNET[uniprotId.idx,,drop=FALSE]$c2.pantherName)))
				
				if(length(gene.geneId) > 1)				{	cat('Please check gene geneId. "',protein_id,'" in human_geneID_DisGeNET_',Func[i],'.txt matches to multiple gene geneId','.\n', sep="")			}
				if(length(gene.symbol) > 1)				{	cat('Please check gene symbol. "',protein_id,'" in human_geneID_DisGeNET_',Func[i],'.txt matches to multiple gene symbol','.\n', sep="")			}
				if(length(gene.description) > 1)		{	cat('Please check gene description. "',protein_id,'" in human_geneID_DisGeNET_',Func[i],'.txt matches to multiple gene description','.\n', sep="")	}
				if(length(gene.pantherName) > 1)		{	cat('Please check gene pantherName. "',protein_id,'" in human_geneID_DisGeNET_',Func[i],'.txt matches to multiple gene pantherName','.\n', sep="")	}
				
				if(gene.description != protein_name)	{	cat('[Warning] Please check gene description. "description" of "',protein_id,'" in human_geneID_DisGeNET_',Func[i],'.txt (',gene.description,') is not equal to its "protein name" in "drug.order" (',protein_name,')','.\n', sep="")	}
				
				A2B.submatrix <- rbind(A2B.submatrix, cbind(protein_id, gene.geneId, gene.symbol, protein_name, gene.pantherName))
			}
			colnames(A2B.submatrix) = c('UniProt ID','Gene ID','Gene Symbol','Gene Name','Protein Class')
			rownames(A2B.submatrix) = c(1:nrow(A2B.submatrix))
			
			if(nrow(A2B.submatrix) != as.numeric(drug.order[k,'protein number'])) {
				cat('[Warning] Protein number of "A2B.submatrix" of "',drug.order[k,'drug'],'" (',nrow(A2B.submatrix),') is not equal to its "','protein number','" (',as.numeric(drug.order[k,'protein number']),')','.\n', sep="")
			}
			
			if(!is.null(A2B.submatrix))
			{
				if(nrow(A2B.submatrix) > 0)
				{
					na.idx = which(is.na(A2B.submatrix), arr.ind=TRUE)
					if(nrow(na.idx) > 0)
					{
						for(ri in 1:nrow(na.idx))
						{
							A2B.submatrix[na.idx[ri,'row'],na.idx[ri,'col']] = ''
						}
					}
				}
			}
			
			if(!is.null(A2B.submatrix))
			{
				if(nrow(A2B.submatrix) > 0)
				{
					Title.INFO = as.character(drug.order[k,c('drug','drug name','drug type','drug group')])
					Title.INFO = Title.INFO[which(!is.na(Title.INFO) & Title.INFO != '')]
					Title.INFO = paste(Title.INFO, collapse="||")
					
					Title = c(paste('> ', Title.INFO, sep=""), rep(NA,ncol(A2B.submatrix)))
					Header = c(NA, colnames(A2B.submatrix))
					SubMatrix = as.matrix(cbind(rownames(A2B.submatrix), A2B.submatrix))
					BlankLine = c(NA, rep(NA,ncol(A2B.submatrix)))
					A2B.submatrix.new = rbind(Title, Header, SubMatrix, BlankLine)
					
					A2B.submatrix.new.all <- rbind(A2B.submatrix.new.all, A2B.submatrix.new)
				}
			}
		}
		
		if(as.numeric(drug.order[k,paste(disease.HeaderTerm,' ','number',sep="")]) > 0)
		{
			# --- disease -------------------------------------------------------------------------------------
			if(disease.level[j] == 'Name') {
				MedGen_id2name = unlist(strsplit(drug.order[k,disease.HeaderTerm],'\\|'))
				MedGen_id = sapply(MedGen_id2name, FUN=function(X) substr(unlist(X), 0, which(unlist(strsplit(unlist(X),''))=='~')[1]-1))
				MedGen_name = sapply(MedGen_id2name, FUN=function(X) substr(unlist(X), which(unlist(strsplit(unlist(X),''))=='~')[1]+1, nchar(unlist(X))))
				
				A2C.submatrix <- c()
				for(d in 1:length(MedGen_id))
				{
					disease_id = MedGen_id[d]
					disease_name = MedGen_name[d]
					
					disease.diseaseId <- unique(as.character(trimws(human_geneID_DisGeNET[trimws(human_geneID_DisGeNET$c1.name) %in% disease_name,,drop=FALSE]$c1.diseaseId)))
					disease.MESH <- unique(as.character(trimws(human_geneID_DisGeNET[trimws(human_geneID_DisGeNET$c1.name) %in% disease_name,,drop=FALSE]$c1.MESH)))
					disease.OMIM <- unique(as.character(trimws(human_geneID_DisGeNET[trimws(human_geneID_DisGeNET$c1.name) %in% disease_name,,drop=FALSE]$c1.OMIM)))
					disease.STY <- unique(as.character(trimws(human_geneID_DisGeNET[trimws(human_geneID_DisGeNET$c1.name) %in% disease_name,,drop=FALSE]$c1.STY)))
					disease.diseaseClassName <- unique(as.character(trimws(human_geneID_DisGeNET[trimws(human_geneID_DisGeNET$c1.name) %in% disease_name,,drop=FALSE]$c1.diseaseClassName)))
					disease.type <- unique(as.character(trimws(human_geneID_DisGeNET[trimws(human_geneID_DisGeNET$c1.name) %in% disease_name,,drop=FALSE]$c1.type)))
					
					if(length(disease.diseaseId) > 1)			{	cat('Please check disease diseaseId. "',disease_name,'" in human_geneID_DisGeNET_',Func[i],'.txt matches to multiple disease diseaseId','.\n', sep="")					}
					if(length(disease.MESH) > 1)				{	cat('Please check disease MESH_id. "',disease_name,'" in human_geneID_DisGeNET_',Func[i],'.txt matches to multiple disease MESH_id','.\n', sep="")						}
					if(length(disease.OMIM) > 1)				{	cat('Please check disease OMIM_id. "',disease_name,'" in human_geneID_DisGeNET_',Func[i],'.txt matches to multiple disease OMIM_id','.\n', sep="")						}
					if(length(disease.STY) > 1)					{	cat('Please check disease STY. "',disease_name,'" in human_geneID_DisGeNET_',Func[i],'.txt matches to multiple disease STY','.\n', sep="")								}
					if(length(disease.diseaseClassName) > 1)	{	cat('Please check disease diseaseClassName. "',disease_name,'" in human_geneID_DisGeNET_',Func[i],'.txt matches to multiple disease diseaseClassName','.\n', sep="")	}
					if(length(disease.type) > 1)				{	cat('Please check disease type. "',disease_name,'" in human_geneID_DisGeNET_',Func[i],'.txt matches to multiple disease type','.\n', sep="")							}
					
					if(disease.diseaseId != disease_id)			{	cat('[Warning] Please check disease diseaseId. "diseaseId" of "',disease_name,'" in human_geneID_DisGeNET_',Func[i],'.txt (',disease.diseaseId,') is not equal to its "disease id" in "drug.order" (',disease_id,')','.\n', sep="")	}
					
					A2C.submatrix <- rbind(A2C.submatrix, cbind(disease_id, disease_name, disease.type, disease.STY, disease.diseaseClassName, disease.MESH, disease.OMIM))
				}
				colnames(A2C.submatrix) = c('UMLS CUI','Name','Type','Semantic Type','MeSH Class','MeSH','OMIM')
				rownames(A2C.submatrix) = c(1:nrow(A2C.submatrix))
			} else {
				MedGen_name = unlist(strsplit(drug.order[k,disease.HeaderTerm],'\\|'))
				
				A2C.submatrix <- c()
				for(d in 1:length(MedGen_name))
				{
					disease_name = MedGen_name[d]
					
					A2C.submatrix <- rbind(A2C.submatrix, cbind(disease_name))
				}
				if(disease.level[j] == 'DisGeNET') {
					colnames(A2C.submatrix) = c('Name')
				} else if(disease.level[j] == 'DisClassName') {
					colnames(A2C.submatrix) = c('MeSH Class')
				} else if(disease.level[j] == 'PantherName') {
					colnames(A2C.submatrix) = c('Protein Class')
				}
				rownames(A2C.submatrix) = c(1:nrow(A2C.submatrix))
			}
			
			if(nrow(A2C.submatrix) != as.numeric(drug.order[k,paste(disease.HeaderTerm,' ','number',sep="")])) {
				cat('[Warning] Disease number of "A2C.submatrix" of "',drug.order[k,'drug'],'" (',nrow(A2C.submatrix),') is not equal to its "',paste(disease.HeaderTerm,' ','number',sep=""),'" (',as.numeric(drug.order[k,paste(disease.HeaderTerm,' ','number',sep="")]),')','.\n', sep="")
			}
			
			if(!is.null(A2C.submatrix))
			{
				if(nrow(A2C.submatrix) > 0)
				{
					na.idx = which(is.na(A2C.submatrix), arr.ind=TRUE)
					if(nrow(na.idx) > 0)
					{
						for(ri in 1:nrow(na.idx))
						{
							A2C.submatrix[na.idx[ri,'row'],na.idx[ri,'col']] = ''
						}
					}
				}
			}
			
			if(!is.null(A2C.submatrix))
			{
				if(nrow(A2C.submatrix) > 0)
				{
					Title.INFO = as.character(drug.order[k,c('drug','drug name','drug type','drug group')])
					Title.INFO = Title.INFO[which(!is.na(Title.INFO) & Title.INFO != '')]
					Title.INFO = paste(Title.INFO, collapse="||")
					
					Title = c(paste('> ', Title.INFO, sep=""), rep(NA,ncol(A2C.submatrix)))
					Header = c(NA, colnames(A2C.submatrix))
					SubMatrix = as.matrix(cbind(rownames(A2C.submatrix), A2C.submatrix))
					BlankLine = c(NA, rep(NA,ncol(A2C.submatrix)))
					A2C.submatrix.new = rbind(Title, Header, SubMatrix, BlankLine)
					
					A2C.submatrix.new.all <- rbind(A2C.submatrix.new.all, A2C.submatrix.new)
				}
			}
		}
	}
	
	if(s == 1) {
		assign(paste(Func[i],'.',disease.level[j],'.','drug.proteinorder','.','A2B.submatrix.new.all',sep=""), A2B.submatrix.new.all)
		assign(paste(Func[i],'.',disease.level[j],'.','drug.proteinorder','.','A2C.submatrix.new.all',sep=""), A2C.submatrix.new.all)
	} else if(s == 2) {
		assign(paste(Func[i],'.',disease.level[j],'.','drug.diseaseorder','.','A2B.submatrix.new.all',sep=""), A2B.submatrix.new.all)
		assign(paste(Func[i],'.',disease.level[j],'.','drug.diseaseorder','.','A2C.submatrix.new.all',sep=""), A2C.submatrix.new.all)
	}
}


# --- write xlsx file by openxlsx package ---------------------------------------------------------
# Importing a big xlsx file into R? https://stackoverflow.com/questions/19147884/importing-a-big-xlsx-file-into-r/43118530#43118530
# Building R for Windows: https://cran.r-project.org/bin/windows/Rtools/
#Sys.setenv("R_ZIPCMD" = "C:/Rtools/bin/zip.exe")	# path to zip.exe

for(s in 1:2)
{
	if(s == 1) {
		prename = paste('./Data/',Sample,'/07_disease2protein2drug/subset/n_',r,'/',folder[f],'/','human','/',Func[i],'/',disease.level[j],'/',Sample,'_exo_','human','_geneID','_DisGeNET','_',Func[i],disease.level[j],'_','drug','_','protein','order',sep="")
	} else if(s == 2) {
		prename = paste('./Data/',Sample,'/07_disease2protein2drug/subset/n_',r,'/',folder[f],'/','human','/',Func[i],'/',disease.level[j],'/',Sample,'_exo_','human','_geneID','_DisGeNET','_',Func[i],disease.level[j],'_','drug','_',gsub(' ','',disease.HeaderTerm),'order',sep="")
	}
	
	wb <- openxlsx::createWorkbook(paste(prename,'.xlsx',sep=""))
	modifyBaseFont(wb, fontSize=12, fontColour="black", fontName="Times New Roman")
	
	for(sn in 0:2)
	{
		if(s == 1) {
			if(sn == 0) {
				sheet = paste(Func[i],'_','drug',sep="")
				if(exists(paste(Func[i],'.',disease.level[j],'.','drug.proteinorder',sep=""))) {
					sheetData = eval(parse(text=paste(Func[i],'.',disease.level[j],'.','drug.proteinorder',sep="")))
				} else {
					sheetData <- c()
				}
			} else if(sn == 1) {
				sheet = paste(Func[i],'_','drug','2','protein',sep="")
				if(exists(paste(Func[i],'.',disease.level[j],'.','drug.proteinorder','.','A2B.submatrix.new.all',sep=""))) {
					sheetData = eval(parse(text=paste(Func[i],'.',disease.level[j],'.','drug.proteinorder','.','A2B.submatrix.new.all',sep="")))
				} else {
					sheetData <- c()
				}
			} else if(sn == 2) {
				sheet = paste(Func[i],'_','drug','2',gsub(' ','',disease.HeaderTerm),sep="")
				if(exists(paste(Func[i],'.',disease.level[j],'.','drug.proteinorder','.','A2C.submatrix.new.all',sep=""))) {
					sheetData = eval(parse(text=paste(Func[i],'.',disease.level[j],'.','drug.proteinorder','.','A2C.submatrix.new.all',sep="")))
				} else {
					sheetData <- c()
				}
			}
		} else if(s == 2) {
			if(sn == 0) {
				sheet = paste(Func[i],'_','drug',sep="")
				if(exists(paste(Func[i],'.',disease.level[j],'.','drug.diseaseorder',sep=""))) {
					sheetData = eval(parse(text=paste(Func[i],'.',disease.level[j],'.','drug.diseaseorder',sep="")))
				} else {
					sheetData <- c()
				}
			} else if(sn == 1) {
				sheet = paste(Func[i],'_','drug','2',gsub(' ','',disease.HeaderTerm),sep="")
				if(exists(paste(Func[i],'.',disease.level[j],'.','drug.diseaseorder','.','A2C.submatrix.new.all',sep=""))) {
					sheetData = eval(parse(text=paste(Func[i],'.',disease.level[j],'.','drug.diseaseorder','.','A2C.submatrix.new.all',sep="")))
				} else {
					sheetData <- c()
				}
			} else if(sn == 2) {
				sheet = paste(Func[i],'_','drug','2','protein',sep="")
				if(exists(paste(Func[i],'.',disease.level[j],'.','drug.diseaseorder','.','A2B.submatrix.new.all',sep=""))) {
					sheetData = eval(parse(text=paste(Func[i],'.',disease.level[j],'.','drug.diseaseorder','.','A2B.submatrix.new.all',sep="")))
				} else {
					sheetData <- c()
				}
			}
		}
		
		if(!is.null(sheetData))
		{
			if(nrow(sheetData) > 0)
			{
				if(sn == 0) {
					addWorksheet(wb, sheet)
					writeData(wb, sheet, sheetData, colNames=TRUE, rowNames=FALSE, keepNA=FALSE)
					setColWidths(wb, sheet, cols=c(1:ncol(sheetData)), widths=c(12.5,c(30,17.5,22.5),rep(15,ncol(sheetData)-(1+3))))
					setRowHeights(wb, sheet, rows=c(1:(1+nrow(sheetData))), heights=16.5)
					freezePane(wb, sheet, firstActiveRow=2)
					style <- createStyle(halign="left", valign="center")
					addStyle(wb, sheet, style, rows=c(1:(1+nrow(sheetData))), cols=c(1:ncol(sheetData)), gridExpand=TRUE)
					style <- createStyle(textDecoration="bold")
					addStyle(wb, sheet, style, rows=1, cols=c(1:ncol(sheetData)), stack=TRUE)
					
					right_align.idx = which(colnames(sheetData) %in% c('protein number',paste(disease.HeaderTerm,' ','number',sep=""),'drug number'))
					if(length(right_align.idx) > 0) {
						style <- createStyle(halign="right", valign="center")
						addStyle(wb, sheet, style, rows=c(2:(1+nrow(sheetData))), cols=right_align.idx, gridExpand=TRUE, stack=TRUE)
					}
				} else if(sn == 1) {
					if(s == 1) {
						addWorksheet(wb, sheet)
						writeData(wb, sheet, sheetData, colNames=FALSE, rowNames=FALSE, keepNA=FALSE)
						setColWidths(wb, sheet, cols=c(1:ncol(sheetData)), widths=c(8.5,12.5,rep(12.5,2),30,22.5,rep(15,ncol(sheetData)-(1+1+2+1+1))))
						setRowHeights(wb, sheet, rows=c(1:nrow(sheetData)), heights=16.5)
						style <- createStyle(halign="left", valign="center")
						addStyle(wb, sheet, style, rows=c(1:nrow(sheetData)), cols=c(1:ncol(sheetData)), gridExpand=TRUE)
					} else if(s == 2) {
						addWorksheet(wb, sheet)
						writeData(wb, sheet, sheetData, colNames=FALSE, rowNames=FALSE, keepNA=FALSE)
						if(disease.level[j] == 'Name') {
							setColWidths(wb, sheet, cols=c(1:ncol(sheetData)), widths=c(8.5,12.5,30,12.5,20,30,rep(12.5,ncol(sheetData)-(1+1+1+1+1+1))))
						} else if(disease.level[j] == 'PantherName') {
							setColWidths(wb, sheet, cols=c(1:ncol(sheetData)), widths=c(8.5,rep(22.5,ncol(sheetData)-1)))
						} else {
							setColWidths(wb, sheet, cols=c(1:ncol(sheetData)), widths=c(8.5,rep(30,ncol(sheetData)-1)))
						}
						setRowHeights(wb, sheet, rows=c(1:nrow(sheetData)), heights=16.5)
						style <- createStyle(halign="left", valign="center")
						addStyle(wb, sheet, style, rows=c(1:nrow(sheetData)), cols=c(1:ncol(sheetData)), gridExpand=TRUE)
					}
				} else if(sn == 2) {
					if(s == 1) {
						addWorksheet(wb, sheet)
						writeData(wb, sheet, sheetData, colNames=FALSE, rowNames=FALSE, keepNA=FALSE)
						if(disease.level[j] == 'Name') {
							setColWidths(wb, sheet, cols=c(1:ncol(sheetData)), widths=c(8.5,12.5,30,12.5,20,30,rep(12.5,ncol(sheetData)-(1+1+1+1+1+1))))
						} else if(disease.level[j] == 'PantherName') {
							setColWidths(wb, sheet, cols=c(1:ncol(sheetData)), widths=c(8.5,rep(22.5,ncol(sheetData)-1)))
						} else {
							setColWidths(wb, sheet, cols=c(1:ncol(sheetData)), widths=c(8.5,rep(30,ncol(sheetData)-1)))
						}
						setRowHeights(wb, sheet, rows=c(1:nrow(sheetData)), heights=16.5)
						style <- createStyle(halign="left", valign="center")
						addStyle(wb, sheet, style, rows=c(1:nrow(sheetData)), cols=c(1:ncol(sheetData)), gridExpand=TRUE)
					} else if(s == 2) {
						addWorksheet(wb, sheet)
						writeData(wb, sheet, sheetData, colNames=FALSE, rowNames=FALSE, keepNA=FALSE)
						setColWidths(wb, sheet, cols=c(1:ncol(sheetData)), widths=c(8.5,12.5,rep(12.5,2),30,22.5,rep(15,ncol(sheetData)-(1+1+2+1+1))))
						setRowHeights(wb, sheet, rows=c(1:nrow(sheetData)), heights=16.5)
						style <- createStyle(halign="left", valign="center")
						addStyle(wb, sheet, style, rows=c(1:nrow(sheetData)), cols=c(1:ncol(sheetData)), gridExpand=TRUE)
					}
				}
			}
		}
	}
	
	if(length(wb$sheet_names) > 0) {
		openxlsx::saveWorkbook(wb, file=paste(prename,'.xlsx',sep=""), overwrite=TRUE)
	}
}

