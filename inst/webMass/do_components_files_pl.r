###############################################################################################
# Run a file-wise componentization ############################################################
###############################################################################################

	###########################################################################################
	do_isot <- (logfile$workflow[names(logfile$workflow) == "isotopologues"] == "yes")
	do_addu <- (logfile$workflow[names(logfile$workflow) == "adducts"] == "yes")
	do_homol <- (logfile$workflow[names(logfile$workflow) == "homologues"] == "yes")
	
	if( do_isot | do_addu ){ # homol alone not sufficient to run nontarget::combine
	
		#######################################################################################
		measurements <- read.csv(file = file.path(logfile[[1]], "dataframes", "measurements"), colClasses = "character");
		#for_IDs <- measurements$ID
		if(mute(logfile$parameters$prof_select == "TRUE")){
			for_IDs <- measurements[(measurements$include=="TRUE") & (measurements$components_files=="FALSE") & (measurements$profiled!="FALSE"),]$ID
		}else{
			for_IDs <- measurements[(measurements$include=="TRUE") & (measurements$components_files=="FALSE") ,]$ID		
		}	
		#######################################################################################
		if(length(for_IDs)){	
			for(i in for_IDs){
				if(file.exists(file.path(logfile[[1]],"results","componentization","components",paste(i)))){ 
					file.remove(file.path(logfile[[1]],"results","componentization","components",paste(i)))
				}
			}
			if(FALSE){ # for debugging - outside clusters
				for(i in for_IDs) combine2_wrap(x = i, logfile, measurements) 
			}
			clusterEvalQ(cl = clus,{rm(list=ls()); NULL})
			clusterExport(cl = clus, varlist = c("do_isot", "do_addu", "do_homol"), envir = environment())
			cluster_results <- clusterApplyLB(cl = clus, 
				x = for_IDs, 
				fun = enviMass:::combine2_wrap, 
				logfile = logfile,
				measurements = measurements
			)
			clusterEvalQ(cl = clus,{rm(list=ls()); NULL})	
		}
		#######################################################################################	
		measurements[!is.na(match(measurements$ID, for_IDs)), "components_files"] <- "TRUE"
		write.csv(measurements, file = file.path(logfile[[1]], "dataframes", "measurements"), row.names = FALSE);		
		rm(measurements)
		
	}
	
	

	

