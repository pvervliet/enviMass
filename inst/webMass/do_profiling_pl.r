	
	############################################################################################
	# REMOVE OLD RESULTS #######################################################################
	if(file.exists(file.path(as.character(logfile[[1]]),"results","profileList_pos"))){file.remove(file.path(as.character(logfile[[1]]),"results","profileList_pos"))}
	if(file.exists(file.path(as.character(logfile[[1]]),"results","profileList_pos_copy"))){file.remove(file.path(as.character(logfile[[1]]),"results","profileList_pos_copy"))}
	if(file.exists(file.path(as.character(logfile[[1]]),"results","profileList_neg"))){file.remove(file.path(as.character(logfile[[1]]),"results","profileList_neg"))}
	if(file.exists(file.path(as.character(logfile[[1]]),"results","profileList_neg_copy"))){file.remove(file.path(as.character(logfile[[1]]),"results","profileList_neg_copy"))}
	if(any(objects(envir = as.environment(".GlobalEnv")) == "peaklist")){rm(peaklist, envir = as.environment(".GlobalEnv"))}
	if(any(objects() == "peaklist")){rm(peaklist)}
	if(any(objects(envir = as.environment(".GlobalEnv")) == "profileList_pos")){rm(profileList_pos, envir = as.environment(".GlobalEnv"))}
	if(any(objects() == "profileList_pos")){rm(profileList_pos)}		
	if(any(objects(envir = as.environment(".GlobalEnv")) == "profileList_neg")){rm(profileList_neg, envir = as.environment(".GlobalEnv"))}
	if(any(objects() == "profileList_neg")){rm(profileList_neg)}
	############################################################################################
	
	############################################################################################
    measurements <- read.csv(file=file.path(logfile[[1]], "dataframes", "measurements"), colClasses = "character");
	measurements <- measurements[measurements[,"include"] == "TRUE",, drop = FALSE]
	if(logfile$parameters$prof_select=="TRUE"){
		measurements<-measurements[measurements[,names(measurements) == "profiled"] == "TRUE",]
	}
	############################################################################################

	############################################################################################		
	if(any(measurements[,"Mode"]=="positive")){
	 
		#########################################################################################
		# set up profileList_pos ################################################################
	    profileList_pos <- list(0)
	    profileList_pos[[1]] <- data.frame(TRUE, FALSE, FALSE, FALSE)    # state
	    colnames(profileList_pos[[1]]) <- c("peaks?", "agglom?", "profiling", "trends?")
	    profileList_pos[[2]] <- 0  # peaks
	    profileList_pos[[3]] <- 0  # datetime
	    profileList_pos[[4]] <- 0  # time
	    profileList_pos[[5]] <- 0  # place
	    profileList_pos[[6]] <- 0  # index_agglom
	    profileList_pos[[7]] <- 0  # index_prof
	    profileList_pos[[8]] <- 0  # parameters
	    profileList_pos[[9]] <- 0  # sample type
	    names(profileList_pos) <- c("state","peaks","datetime","sampleID","place",
			"index_agglom","index_prof","parameters","type")
		for_files<-enviMass:::get_measurement_IDs(
				logfile,
				sets = as.numeric(logfile$parameters$prof_maxfiles),
				ion_mode = "positive",
				until = logfile$parameters$upto_file,
				selective = logfile$parameters$prof_select,
				types = c("sample", "blank", "spiked"),
				places = FALSE,
				check_exist = TRUE
			)
    	profileList_pos[["datetime"]] <- for_files[["datetime"]];
    	profileList_pos[["sampleID"]] <- for_files[["sampleID"]];
		profileList_pos[["place"]] <- for_files[["locus"]];
    	profileList_pos[["type"]] <- for_files[["typus"]];
		#########################################################################################
		# load peaks into profileList_pos #######################################################
		clusterEvalQ(cl = clus,{rm(list=ls()); NULL})
		cluster_results <- clusterApplyLB(cl = clus, 
			x = for_files[["sampleID"]], 
			fun = startprofiles_pl, 
			logfile = logfile,
			frac = FALSE,
			blind_omit=as.logical(logfile$parameters$blind_omit)
		)
		clusterEvalQ(cl = clus,{rm(list=ls()); NULL})
		profileList_pos[["peaks"]] <- do.call(rbind, cluster_results)
		rm(cluster_results)
		profileList_pos[["peaks"]] <- profileList_pos[["peaks"]][order(profileList_pos[["peaks"]][,"m/z"], decreasing = FALSE),]
		#########################################################################################
		if(any(profileList_pos[[2]][,2] == 0)){stop("\n issue in do_profiling: zero intensities detected. Try to rerun the workflow including the peakpicking, using -> Settings -> General -> Reset project including peak picking.")}
		profileList_pos <- agglomer(
			profileList_pos,
			dmass = (as.numeric(logfile$parameters$prof_dmz) + 1),
			ppm = as.logical(as.character(logfile$parameters$prof_ppm)),
			dret = (as.numeric(logfile$parameters$prof_drt) + 10)
		)
		#########################################################################################
		if( (logfile$workflow[names(logfile$workflow)=="replicates"] == "no") & (logfile$parameters$replicates_prof == "no") ){ 				
			profileList_pos <- partcluster_pl(
				profileList = profileList_pos,
				dmass = as.numeric(logfile$parameters$prof_dmz),
				ppm = as.logical(as.character(logfile$parameters$prof_ppm)),
				dret = as.numeric(logfile$parameters$prof_drt),
				replicates = FALSE,
				IDs = FALSE,
				clus = clus
			)
		}else{ # run a profiling in the replicate groups first
			replicates <- measurements[measurements[,"Mode"] == "positive", "tag3"]
			IDs <- measurements[measurements[,"Mode"]=="positive", "ID"]
			profileList_pos <- partcluster_pl(
				profileList = profileList_pos,
				dmass = as.numeric(logfile$parameters$prof_dmz),
				ppm = as.logical(as.character(logfile$parameters$prof_ppm)),
				dret = as.numeric(logfile$parameters$prof_drt),
				replicates = replicates,
				IDs = IDs,
				clus = clus
			)
		}
		########################################################################################
		if(mute(as.logical(logfile$parameters$test))){
			####################################################################################
			# profile IDs correct? #############################################################
			for(i in 1:dim(profileList_pos[["index_prof"]])[1]){
				if(
					!all(profileList_pos[["peaks"]][	
						profileList_pos[["index_prof"]][i, "start_ID"]:profileList_pos[["index_prof"]][i, "end_ID"]
					,"profileIDs"] == i)
				){
					stop("\n Debug partcluster_pl.r at #1")
				}
			}
			####################################################################################	
		}
		#########################################################################################
		profileList_pos <<- profileList_pos
		save(profileList_pos, file = file.path(as.character(logfile[[1]]), "results", "profileList_pos"), compress = FALSE)
		profileList_pos_copy <- profileList_pos
		save(profileList_pos_copy, file = file.path(as.character(logfile[[1]]), "results", "profileList_pos_copy"), compress = FALSE) # used for screening - does not include modifications of downstream compound subtraction		
		rm(profileList_pos_copy)
		links_peaks_pos<-list(); # each entry with 6 lists itself: targets, IS, EIC_correl, isotop, adducts, homol; links_profiles_pos defined in do_components_profiles.r
		save(links_peaks_pos, file = file.path(as.character(logfile[[1]]), "results", "links_peaks_pos"));	
		#########################################################################################
		
	}
	############################################################################################		

	############################################################################################		
	if(any(measurements[,"Mode"]=="negative")){
	
		#########################################################################################
		# set up profileList_neg ################################################################
	    profileList_neg <- list(0)
	    profileList_neg[[1]] <- data.frame(TRUE, FALSE, FALSE, FALSE)    # state
	    colnames(profileList_neg[[1]]) <- c("peaks?", "agglom?", "profiling", "trends?")
	    profileList_neg[[2]] <- 0  # peaks
	    profileList_neg[[3]] <- 0  # datetime
	    profileList_neg[[4]] <- 0  # time
	    profileList_neg[[5]] <- 0  # place
	    profileList_neg[[6]] <- 0  # index_agglom
	    profileList_neg[[7]] <- 0  # index_prof
	    profileList_neg[[8]] <- 0  # parameters
	    profileList_neg[[9]] <- 0  # sample type
	    names(profileList_neg) <- c("state", "peaks", "datetime", "sampleID", "place",
			"index_agglom", "index_prof", "parameters", "type")
		for_files <- enviMass:::get_measurement_IDs(
				logfile,
				sets = as.numeric(logfile$parameters$prof_maxfiles),
				ion_mode = "negative",
				until = logfile$parameters$upto_file,
				selective = logfile$parameters$prof_select,
				types = c("sample","blank","spiked"),
				places = FALSE,
				check_exist = TRUE
			)
    	profileList_neg[["datetime"]] <- for_files[["datetime"]];
    	profileList_neg[["sampleID"]] <- for_files[["sampleID"]];
		profileList_neg[["place"]] <- for_files[["locus"]];
    	profileList_neg[["type"]] <- for_files[["typus"]];
		#########################################################################################
		# load peaks into profileList_neg #######################################################
		clusterEvalQ(cl = clus,{rm(list=ls()); NULL})
		cluster_results <- clusterApplyLB(cl = clus, 
			x = for_files[["sampleID"]], 
			fun = startprofiles_pl, 
			logfile = logfile,
			frac = FALSE,
			blind_omit=as.logical(logfile$parameters$blind_omit)
		)
		clusterEvalQ(cl = clus,{rm(list=ls()); NULL})
		profileList_neg[["peaks"]] <- do.call(rbind,cluster_results)
		profileList_neg[["peaks"]] <- profileList_neg[["peaks"]][order(profileList_neg[["peaks"]][,"m/z"],decreasing=FALSE),]
		#########################################################################################
		if(any(profileList_neg[[2]][,2]==0)){stop("\n issue in do_profiling: zero intensities detected. Try to rerun the workflow including the peakpicking, using -> Settings -> General -> Reset project including peak picking.")}
		profileList_neg <- agglomer(
			profileList_neg,
			dmass = (as.numeric(logfile$parameters$prof_dmz)+1),
			ppm = as.logical(as.character(logfile$parameters$prof_ppm)),
			dret = (as.numeric(logfile$parameters$prof_drt)+10)
		)
		#########################################################################################
		if( (logfile$workflow[names(logfile$workflow)=="replicates"] == "no") & (logfile$parameters$replicates_prof == "no") ){ 				
			profileList_neg <- partcluster_pl(
				profileList = profileList_neg,
				dmass = as.numeric(logfile$parameters$prof_dmz),
				ppm = as.logical(as.character(logfile$parameters$prof_ppm)),
				dret = as.numeric(logfile$parameters$prof_drt),
				replicates = FALSE,
				IDs = FALSE,
				clus = clus
			)
		}else{ # run a profiling in the replicate groups first
			replicates <- measurements[measurements[,"Mode"] == "negative", "tag3"]
			IDs <- measurements[measurements[,"Mode"] == "negative", "ID"]
			profileList_neg <- partcluster_pl(
				profileList = profileList_neg,
				dmass = as.numeric(logfile$parameters$prof_dmz),
				ppm = as.logical(as.character(logfile$parameters$prof_ppm)),
				dret = as.numeric(logfile$parameters$prof_drt),
				replicates = replicates,
				IDs = IDs,
				clus = clus
			)
		}
		########################################################################################
		if(mute(as.logical(logfile$parameters$test))){
			####################################################################################
			# profile IDs correct? #############################################################
			for(i in 1:dim(profileList_pos[["index_prof"]])[1]){
				if(
					!all(profileList_pos[["peaks"]][	
						profileList_pos[["index_prof"]][i, "start_ID"]:profileList_pos[["index_prof"]][i, "end_ID"]
					,"profileIDs"] == i)
				){
					stop("\n Debug partcluster_pl.r at #1")
				}
			}
			####################################################################################	
		}
		#########################################################################################
		profileList_neg<<-profileList_neg
		save(profileList_neg,file=file.path(as.character(logfile[[1]]),"results","profileList_neg"),compress=FALSE);
		profileList_neg_copy<-profileList_neg
		save(profileList_neg_copy,file=file.path(as.character(logfile[[1]]),"results","profileList_neg_copy"),compress=FALSE); # used for screening - does not include modifications of downstream compound subtraction		
		links_peaks_neg<-list(); # each entry with 6 lists itself: targets, IS, EIC_correl, isotop, adducts, homol; links_profiles_neg defined in do_components_profiles.r
		save(links_peaks_neg,file=file.path(as.character(logfile[[1]]),"results","links_peaks_neg"));	
		#########################################################################################
		
	}
	############################################################################################		
