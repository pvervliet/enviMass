######################################################################################################################
# Run isotopologue grouping, filewise & parallel #####################################################################
###################################################################################################################### 	

    measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
	cat("Isotopologue grouping: ")
	load(file.path(logfile[[1]],"dataframes","quantiz"))
	if((quantiz$R_set!=logfile$parameters$resolution) & (quantiz$R_set!="Sciex_all")){
		cat("\n WARNING: seems the quantized data for isotopologue grouping does NOT MATCH your selected resolution! Please resolve this issue.")
	}
	if(mute(logfile$parameters$prof_select=="TRUE")){
		for_IDs<-measurements[(measurements$include=="TRUE") & (measurements$isotopologues=="FALSE") & (measurements$profiled!="FALSE"),]$ID
	}else{
		for_IDs<-measurements[(measurements$include=="TRUE") & (measurements$isotopologues=="FALSE") ,]$ID		
	}
	if(length(for_IDs)){
		clusterEvalQ(cl = clus,{rm(list=ls()); NULL})
		clusterExport(cl = clus, varlist = c("quantiz"), envir = environment())
		cluster_results <- clusterApplyLB(cl = clus, 
			x = for_IDs, 
			fun = pattern_search3_wrap, 
			logfile = logfile
		)
		clusterEvalQ(cl = clus,{rm(list=ls()); NULL})
	}
	measurements[!is.na(match(measurements$ID,for_IDs)),"isotopologues"] <- "TRUE"
	write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
	rm(measurements)
