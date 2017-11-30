
    ############################################################################
	measurements <- read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
	replic <- (measurements$tag3 [(measurements$tag3 != "FALSE") & (measurements$include == "TRUE")])
    ############################################################################

    ############################################################################   
    # clean all peaklists until replicate intersection is made sample-wise #####
# > BAUSTELLE
    clean_peaklists <- function(x, logfile){
    	
    	########################################################################
    	if(file.exists(file.path(logfile[[1]], "peaklist", x))){
			load( file = file.path(logfile[[1]], "peaklist", x), envir=environment()); 
			peaklist[,colnames(peaklist)=="keep"] <- 1  
			save(peaklist, file = file.path(logfile[[1]], "peaklist", x))
			rm(peaklist)
		}
    	########################################################################

    }
	clusterEvalQ(cl = clus,{rm(list = ls()); NULL})
	cluster_results <- clusterApplyLB(cl = clus, 
		x = measurements$ID, 
		fun = clean_peaklists,
		logfile = logfile
	)
	clusterEvalQ(cl = clus,{rm(list = ls()); NULL})	  
# < BAUSTELLE
    ############################################################################

    ############################################################################
	if(length(replic)){

		########################################################################
		replic<-replic[duplicated(replic)]
		replic<-unique(replic)
		######################################################################  
		clusterEvalQ(cl = clus,{rm(list = ls()); NULL})
		cluster_results <- clusterApplyLB(cl = clus, 
			x = replic, 
			fun = enviMass:::replicates_wrap, 
			logfile = logfile,
			measurements = measurements,
			ppm = logfile$parameters$replicate_ppm,
			mz_tol = as.numeric(logfile$parameters$replicate_dmz),
			rt_tol = as.numeric(logfile$parameters$replicate_delRT),
			int_tol = (10^(as.numeric(logfile$parameters$replicate_IS_dInt)))
		)
		clusterEvalQ(cl = clus,{rm(list = ls()); NULL})	  
		######################################################################
		cluster_results <-unlist(cluster_results)
		#if(any(cluster_results=="done")){
			#measurements$replicates[match(for_IDs_pos[cluster_results=="done"], measurements$ID)]<-"TRUE" # this measurements column does not exist yet!
		#}
		rm(cluster_results, replic)
		######################################################################

	}
    ############################################################################


	
