# Run adduct grouping, filewise
    measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");	
	if(mute(logfile$parameters$prof_select=="TRUE")){
		for_IDs<-measurements[(measurements$include=="TRUE") & (measurements$adducts=="FALSE") & (measurements$profiled!="FALSE"),]$ID
	}else{
		for_IDs<-measurements[(measurements$include=="TRUE") & (measurements$adducts=="FALSE") ,]$ID		
	}
	if(length(for_IDs)){
		for(i in for_IDs){
			if(file.exists(file.path(logfile[[1]],"results","componentization","adducts",paste("full",for_file,sep="_")))) 
				file.remove(file.path(logfile[[1]],"results","componentization","adducts",paste("full",for_file,sep="_")))
			if(file.exists(file.path(logfile[[1]],"results","componentization","adducts",paste(for_file,sep="")))) 
				file.remove(file.path(logfile[[1]],"results","componentization","adducts",paste(for_file,sep="")))
		}
		clusterEvalQ(cl = clus,{rm(list=ls()); NULL})
		clusterExport(cl = clus, varlist = c("adducts"), envir = environment())
		cluster_results <- clusterApplyLB(cl = clus, 
			x = for_IDs, 
			fun = enviMass:::adduct_search2_wrap, 
			logfile = logfile,
			measurements = measurements
		)
		clusterEvalQ(cl = clus,{rm(list=ls()); NULL})
	}
	measurements[!is.na(match(measurements$ID,for_IDs)),"adducts"]<-"TRUE"
	write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
	rm(measurements)	