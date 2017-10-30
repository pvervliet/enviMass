# Run homologue series detection
	
	####################################################################################
    measurements <- read.csv(file=file.path(logfile[[1]],"dataframes","measurements"), colClasses = "character");
	#measurements[,names(measurements)=="homologues"]<-"FALSE"
	if(mute(logfile$parameters$prof_select=="TRUE")){
		for_IDs<-measurements[(measurements$include=="TRUE") & (measurements$homologues=="FALSE") & (measurements$profiled!="FALSE"),]$ID
	}else{
		for_IDs<-measurements[(measurements$include=="TRUE") & (measurements$homologues=="FALSE") ,]$ID		
	}
	intstand <- read.table(file=file.path(logfile[[1]],"dataframes","IS.txt"),header=TRUE,sep="\t",colClasses = "character");
	targets <- read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character");
	####################################################################################	
	####################################################################################
	if(length(for_IDs)){
		if(logfile$parameters$homol_units[1]!="FALSE"){
			these <- enviPat::check_chemform(isotopes,strsplit(logfile$parameters$homol_units,",")[[1]])[,3]
			mzfilter <- c(enviPat::check_chemform(isotopes,strsplit(logfile$parameters$homol_units,",")[[1]])[,3] %*% t(1/as.numeric(strsplit(logfile$parameters$homol_charges,",")[[1]])))
			mzfilter <- unique(mzfilter);
			elements <- unique(unlist(sapply(enviMass::check_chemform(isotopes, strsplit(logfile$parameters$homol_units,",")[[1]], get_list = TRUE),names)))
			use_minmz <- (min(mzfilter) - .1)
			use_maxmz <- (max(mzfilter) + .1)
		}else{
			mzfilter <- FALSE
			elements <- unique(as.character(isotopes[,1])[1:295]) #then use all available elements
			use_minmz <- as.numeric(logfile$parameters$homol_minmz)
			use_maxmz <- as.numeric(logfile$parameters$homol_maxmz)		
		}
		for(i in for_IDs){
			if(file.exists(file.path(logfile[[1]], "results", "componentization", "homologues", paste(i, sep = "_")))) file.remove(file.path(logfile[[1]], "results", "componentization", "homologues", paste(i, sep = "_")))
			if(file.exists(file.path(logfile[[1]], "results", "componentization", "homologues", paste("full" , i, sep = "_")))) file.remove(file.path(logfile[[1]], "results", "componentization", "homologues", paste("full", i, sep = "_")))
		}
		if(FALSE){ # for debugging - outside clusters
			for(i in for_IDs) homol_search2_wrap(x = i, logfile) 
		}
		clusterEvalQ(cl = clus,{rm(list=ls()); NULL})
		clusterExport(cl = clus, varlist = c("mzfilter", "elements", "use_minmz", "use_maxmz", "isotopes", "measurements", "intstand", "targets"), envir = environment())
		cluster_results <- clusterApplyLB(cl = clus, 
			x = for_IDs, 
			fun = enviMass:::homol_search2_wrap, 
			logfile = logfile
		)
		clusterEvalQ(cl = clus,{rm(list=ls()); NULL})
	}
	measurements[!is.na(match(measurements$ID, for_IDs)),"homologues"] <- "TRUE"
	write.csv(measurements, file = file.path(logfile[[1]], "dataframes", "measurements"), row.names = FALSE);	
	####################################################################################	
	rm(mzfilter, elements, use_minmz, use_maxmz, measurements)
	####################################################################################	

	
	
	
	