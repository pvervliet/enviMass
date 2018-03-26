######################################################################################
# Run EIC correlation, filewise ######################################################
######################################################################################

	####################################################################################
    measurements <- read.csv(file = file.path(logfile[[1]],"dataframes","measurements"), colClasses = "character");
	#measurements[,names(measurements)=="EIC_correlation"]<-"FALSE"
	if(mute(logfile$parameters$prof_select == "TRUE")){
		for_IDs <- measurements[(measurements$include == "TRUE") & (measurements$EIC_correlation == "FALSE") & (measurements$profiled != "FALSE"),]$ID
	}else{
		for_IDs <- measurements[(measurements$include == "TRUE") & (measurements$EIC_correlation == "FALSE") ,]$ID		
	}
	####################################################################################	
	do_cor <- TRUE			# omit correlation
	####################################################################################
	
	####################################################################################	
	if(length(for_IDs)){
		clusterEvalQ(cl = clus,{rm(list=ls()); NULL})
		#clusterExport(cl = clus, varlist = c(""), envir = environment())
		#for(x in for_IDs) enviMass:::EIC_correlation_wrap(x, logfile, do_cor) # for debugging
		cluster_results <- clusterApplyLB(cl = clus, 
			x = for_IDs, 
			fun = enviMass:::EIC_correlation_wrap, 
			logfile = logfile,
			do_cor = do_cor
		)
		clusterEvalQ(cl = clus,{rm(list = ls()); NULL})
	}
	####################################################################################	
	
	measurements[!is.na(match(measurements$ID,for_IDs)),"EIC_correlation"] <- "TRUE"
	write.csv(measurements, file = file.path(logfile[[1]],"dataframes","measurements"), row.names = FALSE);		
	rm(measurements)

