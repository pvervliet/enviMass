
    ############################################################################
    # retrieve monoisotopic masses for IS ######################################
    ############################################################################
    measurements <- read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
	#measurements$recal<-"FALSE"
	reset_mass <- function( x, logfile){
		for_file <- x
		if(any(objects(envir = as.environment(".GlobalEnv")) == "peaklist")){rm(peaklist, envir = as.environment(".GlobalEnv"))}
		if(any(objects() == "peaklist")){rm(peaklist)}
		load(file = file.path(logfile[[1]], "peaklist", for_file), envir = as.environment(".GlobalEnv"));
		peaklist[,"m/z_corr"] <- peaklist[,"m/z"];
		save(peaklist,file = file.path(logfile[[1]], "peaklist", for_file))
		rm(peaklist)
		return("done")
	}
	if(logfile$parameters$recal_include_pos == "TRUE"){
		for_IDs_pos <- measurements[ (measurements$include == "TRUE") & (measurements$recal == "FALSE") & (measurements$Mode == "positive"),]$ID
	}else{
		for_IDs_pos <- c()
		for_IDs_pos_non <- measurements[ (measurements$include == "TRUE") & (measurements$recal == "FALSE") & (measurements$Mode == "positive"),]$ID
		if(length(for_IDs_pos_non)){ # e.g. if logfile$parameters$recal_include_pos set from TRUE to FALSE
			##################################################################
			clusterEvalQ(cl = clus,{rm(list = ls()); NULL})
			cluster_results <- clusterApplyLB(cl = clus, 
					x = for_IDs_pos_non, 
					fun = reset_mass, 
					logfile = logfile	
			)
			clusterEvalQ(cl = clus,{rm(list = ls()); NULL})	  
			measurements$recal[match(for_IDs_pos_non, measurements$ID)] <- "TRUE"
			##################################################################
		}
	}
	if(logfile$parameters$recal_include_neg == "TRUE"){	
		for_IDs_neg <- measurements[ (measurements$include == "TRUE") & (measurements$recal == "FALSE") & (measurements$Mode == "negative"),]$ID
	}else{
		for_IDs_neg <- c()
		for_IDs_neg_non <- measurements[ (measurements$include == "TRUE") & (measurements$recal == "FALSE") & (measurements$Mode == "negative"),]$ID
		if(length(for_IDs_neg_non)){ # e.g. if logfile$parameters$recal_include_neg set from TRUE to FALSE
			##################################################################
			clusterEvalQ(cl = clus,{rm(list = ls()); NULL})
			cluster_results <- clusterApplyLB(cl = clus, 
					x = for_IDs_neg_non, 
					fun = reset_mass, 
					logfile = logfile	
			)
			clusterEvalQ(cl = clus,{rm(list = ls()); NULL})	  
			measurements$recal[match(for_IDs_neg_non, measurements$ID)] <- "TRUE"
			##################################################################
		}
	}	
    ############################################################################	
	
	
    ############################################################################	
    if(length(for_IDs_pos)){ # positive ########################################
	
		########################################################################
		mz_pos <- c();
		RT_pos <- c();	
		if(logfile$parameters$recal_use_pos == "Internal standards"){
			if(file.exists(file.path(logfile[[1]],"results","intmass_pos_IS"))){
				if(any(objects(envir=as.environment(".GlobalEnv")) == "intmass_pos_IS")){rm(intmass_pos_IS, envir = as.environment(".GlobalEnv"))}
				if(any(objects() == "intmass_pos_IS")){rm(intmass_pos_IS)}
				load(file = file.path(logfile[[1]], "results", "intmass_pos_IS"), envir = as.environment(".GlobalEnv"))                    
				mz_pos <- c(mz_pos,intmass_pos_IS[,1]);
				RT_pos <- c(RT_pos,intmass_pos_IS[,2]);
			}else{stop("\n IS recalibration masses not found, positive mode ... check project / settings?")}
		}
		if(logfile$parameters$recal_use_pos == "Target compounds"){
			if(file.exists(file.path(logfile[[1]], "results", "intmass_pos_target"))){	  
				if(any(objects(envir=as.environment(".GlobalEnv")) == "intmass_pos_target")){rm(intmass_pos_target, envir = as.environment(".GlobalEnv"))}
				if(any(objects() == "intmass_pos_target")){rm(intmass_pos_target)}
				load(file = file.path(logfile[[1]], "results", "intmass_pos_target"), envir = as.environment(".GlobalEnv"))                    
				mz_pos <- c(mz_pos,intmass_pos_target[,1]);
				RT_pos <- c(RT_pos,intmass_pos_target[,2]);
			}else{stop("\n Target recalibration masses not found, positive mode ... check project / settings?")}      
		}
		if(logfile$parameters$recal_use_pos == "both"){
			if(file.exists(file.path(logfile[[1]], "results", "intmass_pos_IS"))){	  
				if(any(objects(envir=as.environment(".GlobalEnv")) == "intmass_pos_IS")){rm(intmass_pos_IS, envir = as.environment(".GlobalEnv"))}
				if(any(objects() == "intmass_pos_IS")){rm(intmass_pos_IS)}
				load(file = file.path(logfile[[1]], "results", "intmass_pos_IS"), envir = as.environment(".GlobalEnv"))                    
				mz_pos <- c(mz_pos,intmass_pos_IS[,1]);
				RT_pos <- c(RT_pos,intmass_pos_IS[,2]);
			}else{stop("\n IS recalibration masses not found, positive mode ... check project / settings?")}		
			if(file.exists(file.path(logfile[[1]], "results", "intmass_pos_target"))){	
				if(any(objects(envir=as.environment(".GlobalEnv")) == "intmass_pos_target")){rm(intmass_pos_target, envir = as.environment(".GlobalEnv"))}
				if(any(objects() == "intmass_pos_target")){rm(intmass_pos_target)}		
				load(file = file.path(logfile[[1]], "results", "intmass_pos_target"), envir = as.environment(".GlobalEnv"))                    
				mz_pos <- c(mz_pos,intmass_pos_target[,1]);
				RT_pos <- c(RT_pos,intmass_pos_target[,2]);
			}else{stop("\n Target recalibration masses not found, positive mode ... check project / settings?")}			
		}
		mz_pos <- c(as.numeric(as.character(mz_pos)));
		RT_pos <- c(as.numeric(as.character(RT_pos)));
		######################################################################  
		clusterEvalQ(cl = clus,{rm(list = ls()); NULL})
		clusterExport(cl = clus, varlist = c("mz_pos", "RT_pos"), envir = environment())
		cluster_results <- clusterApplyLB(cl = clus, 
			x = for_IDs_pos, 
			fun = enviMass:::recalib_wrap, 
			logfile = logfile,
			mz = mz_pos,
			tolmz = as.numeric(logfile$parameters$recal_dmz_pos),
			ppm = as.character(logfile$parameters$recal_ppm_pos),
			ret = RT_pos,
			tolret = as.numeric(logfile$parameters$recal_drt_pos),
			max_recal = as.numeric(logfile$parameters$recal_maxdmz_pos)			
		)
		clusterEvalQ(cl = clus,{rm(list = ls()); NULL})	  
		######################################################################
		cluster_results <- unlist(cluster_results)
		if(any(cluster_results == "done")){
			measurements$recal[match(for_IDs_pos[cluster_results == "done"], measurements$ID)] <- "TRUE"
		}
		rm(cluster_results, mz_pos, RT_pos)
		######################################################################

    }
    ##########################################################################		
	
	##########################################################################	
    if(length(for_IDs_neg)){ # negative ######################################
	
		#######################################################################
		mz_neg <- c();
		RT_neg <- c();	
		if(logfile$parameters$recal_use_neg == "Internal standards"){
			if(file.exists(file.path(logfile[[1]], "results", "intmass_neg_IS"))){
				if(any(objects(envir = as.environment(".GlobalEnv")) == "intmass_neg_IS")){rm(intmass_neg_IS, envir = as.environment(".GlobalEnv"))}
				if(any(objects() == "intmass_neg_IS")){rm(intmass_neg_IS)}
				load(file = file.path(logfile[[1]],"results","intmass_neg_IS"), envir = as.environment(".GlobalEnv"))                    
				mz_neg <- c(mz_neg,intmass_neg_IS[,1]);
				RT_neg <- c(RT_neg,intmass_neg_IS[,2]);
			}else{stop("\n IS recalibration masses not found, negative mode ... check project / settings?")}
		}
		if(logfile$parameters$recal_use_neg == "Target compounds"){
			if(file.exists(file.path(logfile[[1]], "results", "intmass_neg_target"))){	  
				if(any(objects(envir = as.environment(".GlobalEnv")) == "intmass_neg_target")){rm(intmass_neg_target, envir = as.environment(".GlobalEnv"))}
				if(any(objects() == "intmass_neg_target")){rm(intmass_neg_target)}
				load(file = file.path(logfile[[1]],"results","intmass_neg_target"), envir = as.environment(".GlobalEnv"))                    
				mz_neg <- c(mz_neg,intmass_neg_target[,1]);
				RT_neg <- c(RT_neg,intmass_neg_target[,2]);
			}else{stop("\n Target recalibration masses not found, negative mode ... check project / settings?")}      
		}
		if(logfile$parameters$recal_use_neg == "both"){
			if(file.exists(file.path(logfile[[1]], "results", "intmass_neg_IS"))){	  
				if(any(objects(envir = as.environment(".GlobalEnv")) == "intmass_neg_IS")){rm(intmass_neg_IS, envir = as.environment(".GlobalEnv"))}
				if(any(objects() == "intmass_neg_IS")){rm(intmass_neg_IS)}
				load(file = file.path(logfile[[1]], "results", "intmass_neg_IS"), envir = as.environment(".GlobalEnv"))                    
				mz_neg <- c(mz_neg,intmass_neg_IS[,1]);
				RT_neg <- c(RT_neg,intmass_neg_IS[,2]);
			}else{stop("\n IS recalibration masses not found, negative mode ... check project / settings?")}		
			if(file.exists(file.path(logfile[[1]], "results", "intmass_neg_target"))){	
				if(any(objects(envir = as.environment(".GlobalEnv")) == "intmass_neg_target")){rm(intmass_neg_target, envir = as.environment(".GlobalEnv"))}
				if(any(objects() == "intmass_neg_target")){rm(intmass_neg_target)}		
				load(file = file.path(logfile[[1]], "results", "intmass_neg_target"), envir = as.environment(".GlobalEnv"))                    
				mz_neg <- c(mz_neg,intmass_neg_target[,1]);
				RT_neg <- c(RT_neg,intmass_neg_target[,2]);
			}else{stop("\n Target recalibration masses not found, negative mode ... check project / settings?")}			
		}
		mz_neg <- c(as.numeric(as.character(mz_neg)));
		RT_neg <- c(as.numeric(as.character(RT_neg)));
		######################################################################  
		clusterEvalQ(cl = clus,{rm(list = ls()); NULL})
		clusterExport(cl = clus, varlist = c("mz_neg", "RT_neg"), envir = environment())
		cluster_results <- clusterApplyLB(cl = clus, 
			x = for_IDs_neg,
			fun = enviMass:::recalib_wrap, 
			logfile = logfile,
			mz = mz_neg,
			tolmz = as.numeric(logfile$parameters$recal_dmz_neg),
			ppm = as.character(logfile$parameters$recal_ppm_neg),
			ret = RT_neg,
			tolret = as.numeric(logfile$parameters$recal_drt_neg),
			max_recal = as.numeric(logfile$parameters$recal_maxdmz_neg)	
		)
		clusterEvalQ(cl = clus,{rm(list = ls()); NULL})	  
		######################################################################
		cluster_results <- unlist(cluster_results)
		if(any(cluster_results == "done")){
			measurements$recal[match(for_IDs_neg[cluster_results=="done"], measurements$ID)] <- "TRUE"
		}
		rm(cluster_results, mz_neg, RT_neg)
		######################################################################
		
    }
    ##########################################################################		
	
	##########################################################################		
	if(TRUE){ # make dummy plots
		not_for_IDs <- measurements$ID[
			(is.na(match(measurements$ID, c(for_IDs_neg, for_IDs_pos)))) &
			(measurements$recal == "FALSE")
		]
		if(length(not_for_IDs)){
			for(i in 1:length(not_for_IDs)){
				png(filename = file.path(logfile[[1]],"pics",paste("recal_",not_for_IDs[i],sep="")), bg = "white", width = 1100, height= 300)
				plot.new()
				plot.window(xlim=c(0,1),ylim=c(0,1))
				text(0.5,0.5,"File excluded from recalibration.",cex=1)
				dev.off();
			}
		}
	}
	write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);	
	rm(measurements)
	##########################################################################		
	