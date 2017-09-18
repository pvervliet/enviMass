########################################################################################
# NEW: IS / target screening results interesection with profiles is now run here - #####
# - not anymore in the do_components_profiles_pl script ################################
########################################################################################

########################################################################################
# REMOVE OLD RESULTS ###################################################################
if(file.exists(file.path(as.character(logfile[[1]]),"results","links_profiles_pos"))){file.remove(file.path(as.character(logfile[[1]]),"results","links_profiles_pos"))}
if(file.exists(file.path(as.character(logfile[[1]]),"results","links_profiles_neg"))){file.remove(file.path(as.character(logfile[[1]]),"results","links_profiles_neg"))}
if(any(objects(envir=as.environment(".GlobalEnv"))=="profileList_pos")){rm(profileList_pos,envir=as.environment(".GlobalEnv"))}
if(any(objects()=="profileList_pos")){rm(profileList_pos)}	
if(any(objects(envir=as.environment(".GlobalEnv"))=="links_peaks_pos")){rm(links_peaks_pos,envir=as.environment(".GlobalEnv"))}
if(any(objects()=="links_peaks_pos")){rm(links_peaks_pos)}				
if(any(objects(envir=as.environment(".GlobalEnv"))=="links_profiles_pos")){rm(links_profiles_pos,envir=as.environment(".GlobalEnv"))}
if(any(objects()=="links_profiles_pos")){rm(links_profiles_pos)}					
if(any(objects(envir=as.environment(".GlobalEnv"))=="profileList_neg")){rm(profileList_neg,envir=as.environment(".GlobalEnv"))}
if(any(objects()=="profileList_neg")){rm(profileList_neg)}	
if(any(objects(envir=as.environment(".GlobalEnv"))=="links_peaks_neg")){rm(links_peaks_neg,envir=as.environment(".GlobalEnv"))}
if(any(objects()=="links_peaks_neg")){rm(links_peaks_neg)}				
if(any(objects(envir=as.environment(".GlobalEnv"))=="links_profiles_neg")){rm(links_profiles_neg,envir=as.environment(".GlobalEnv"))}
if(any(objects()=="links_profiles_neg")){rm(links_profiles_neg)}	
########################################################################################





########################################################################################
# on POSITIVE profiles #################################################################
if(
	file.exists(file.path(logfile[[1]],"results","profileList_pos"))  
){

	####################################################################################	
	if(any(objects(envir=as.environment(".GlobalEnv"))=="profileList_pos")){rm(profileList_pos,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="profileList_pos")){rm(profileList_pos)}	
	if(any(objects(envir=as.environment(".GlobalEnv"))=="links_peaks_pos")){rm(links_peaks_pos,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="links_peaks_pos")){rm(links_peaks_pos)}				
	if(any(objects(envir=as.environment(".GlobalEnv"))=="links_profiles_pos")){rm(links_profiles_pos,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="links_profiles_pos")){rm(links_profiles_pos)}					
	load(file.path(as.character(logfile[[1]]),"results","profileList_pos"),envir=as.environment(".GlobalEnv"));
	# links_peaks_pos<-list(); # each entry with 6 lists itself: targets, IS, EIC_correl, isotop, adducts, homol - defined in do_profiling.r
	load(file.path(as.character(logfile[[1]]),"results","links_peaks_pos"),envir=as.environment(".GlobalEnv"));	
	links_profiles_pos<-list(); # each entry with 6 lists itself: targets, IS, EIC_correl, isotop, adducts, homol
	measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
	use_entries_profiles<-enviMass::find_empty(links_profiles_pos) # also finds gaps
	profileList_pos[["index_prof"]][,"links"]<-0
	with_bar<-FALSE
	####################################################################################

	####################################################################################
	# (1) ANNOTATE TARGET & ISTD SCREENING MACTHES stored in links_peaks_pos ###########
	if(
		(
			(logfile$workflow[names(logfile$workflow)=="subtr"]=="no") |	
			(
				((logfile$workflow[names(logfile$workflow)=="IS_screen"]=="yes") &  (logfile$parameters$subtr_IS!="yes")) |
				((logfile$workflow[names(logfile$workflow)=="target_screen"]=="yes") & (logfile$parameters$subtr_target!="yes"))	
			) 
		) & 
		(length(links_peaks_pos)>0) # anything screened?
	){

		###########################################################################
		if(with_bar){pBar <- txtProgressBar(min = 0, max = dim(profileList_pos[["index_prof"]])[1], style = 3)}
		for(i in 1:dim(profileList_pos[["index_prof"]])[1]){
			if(with_bar){setTxtProgressBar(pBar, i, title = NULL, label = NULL)}
			if(
				any(profileList_pos[["peaks"]][
					(profileList_pos[["index_prof"]][i,"start_ID"]:profileList_pos[["index_prof"]][i,"end_ID"]),"links"
				]!=0)			
			){
			
				###################################################################
				# add a new link to the profile ###################################
				if( profileList_pos[["index_prof"]][i,"links"]==0 ){ 	# establish a new link ...
					if(length(use_entries_profiles)>0){
						at_entry<-use_entries_profiles[1]
						use_entries<-use_entries_profiles[-1]
					}else{
						at_entry<-(length(links_profiles_pos)+1)
					}
					links_profiles_pos[[at_entry]]<-enviMass::new_entry_links_profiles(profileList_pos[["index_prof"]][i,"number_peaks_total"][[1]])
					names(links_profiles_pos)[at_entry]<-as.character(i)
					profileList_pos[["index_prof"]][i,"links"]<-at_entry						
				}else{
					at_entry<-profileList_pos[["index_prof"]][i,"links"]
				}				
				###################################################################
				# search peaks & their links to compounds #########################
				for(j in (profileList_pos[["index_prof"]][i,"start_ID"]:profileList_pos[["index_prof"]][i,"end_ID"])){
					if(profileList_pos[["peaks"]][j,"links"]!=0){
						# add IS link #############################################
						if(	length(links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[2]]) > 0 ){
							for(k in 1:length(links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[2]])){ # if several compound matches exist for this peak
								if( length(links_profiles_pos[[at_entry]][[2]])==0 ){ # make a new entry for profile link to IS
									links_profiles_pos[[at_entry]][[2]]<-
										data.frame(
											links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[2]][[k]][[1]],
											1,
											links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[2]][[k]][[2]],
											stringsAsFactors = FALSE
										)
									names(links_profiles_pos[[at_entry]][[2]])<-c("Compound","Counts","max_score")
								}else{
									at<-which(links_profiles_pos[[at_entry]][[2]][,1] == links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[2]][[k]][[1]])
									if(length(at)>0){ 	# increment existing link ...
										links_profiles_pos[[at_entry]][[2]][at,2]<-(links_profiles_pos[[at_entry]][[2]][at,2]+1)
										if(links_profiles_pos[[at_entry]][[2]][at,2] < links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[2]][[k]][[2]] ){
											links_profiles_pos[[at_entry]][[2]][at,2] <- links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[2]][[k]][[2]]
										}
									}else{	# ... or add a new one?
										links_profiles_pos[[at_entry]][[2]]<-data.frame(
											c(
												links_profiles_pos[[at_entry]][[2]][,1], 
												links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[2]][[k]][[1]]
											),
											c(links_profiles_pos[[at_entry]][[2]][,2],1),
											c(
												links_profiles_pos[[at_entry]][[2]][,3], 
												links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[2]][[k]][[2]]
											),
											stringsAsFactors = FALSE
										)									
										names(links_profiles_pos[[at_entry]][[2]])<-c("Compound","Counts","max_score")
									}
								}
							}	
						}					
						# add target link #########################################
						if(	length(links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[1]])>0 ){
							for(k in 1:length(links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[1]])){ # if several compound matches exist for this peak
								if( length(links_profiles_pos[[at_entry]][[1]])==0 ){ # make a new entry for profile link to IS
									links_profiles_pos[[at_entry]][[1]]<-
										data.frame(
											links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[1]][[k]][[1]],
											1,
											links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[1]][[k]][[2]],
											stringsAsFactors = FALSE
										)
									names(links_profiles_pos[[at_entry]][[1]])<-c("Compound","Counts","max_score")
								}else{
									at<-which(links_profiles_pos[[at_entry]][[1]][,1] == links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[1]][[k]][[1]])
									if(length(at)>0){ 	# increment existing link ...
										links_profiles_pos[[at_entry]][[1]][at,2]<-(links_profiles_pos[[at_entry]][[1]][at,2]+1)
										if(links_profiles_pos[[at_entry]][[1]][at,2] < links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[1]][[k]][[2]] ){
											links_profiles_pos[[at_entry]][[1]][at,2] <- links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[1]][[k]][[2]]
										}
									}else{	# ... or add a new one?
										links_profiles_pos[[at_entry]][[1]]<-data.frame(
											c(
												links_profiles_pos[[at_entry]][[1]][,1], 
												links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[1]][[k]][[1]]
											),
											c(links_profiles_pos[[at_entry]][[1]][,2],1),
											c(
												links_profiles_pos[[at_entry]][[1]][,3], 
												links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[1]][[k]][[2]]
											),
											stringsAsFactors = FALSE
										)									
										names(links_profiles_pos[[at_entry]][[1]])<-c("Compound","Counts","max_score")
									}
								}
							}					
						}
						###########################################################
					}
				}
				###################################################################
			}
		}
		if(with_bar){close(pBar)}
	}
	if(!length(links_profiles_pos)){
		stop("Not enough internal standard (ISTD) screening matches detected to run the ISTD normalization, positive ionization - remove this workflow step?")
	}
	####################################################################################	

	####################################################################################
	# run normalization ################################################################
	# -> screen IS intensity profiles ##################################################
	min_count <- floor(length(profileList_pos[[4]]) * as.numeric(logfile$parameters$ISnorm_percfiles) / 100);
	lis_delint_IS <- list()
	lis_median_IS <- list()
	#lis_RT_IS <- list()
	for(p in 1:length(profileList_pos[["sampleID"]])){
		lis_delint_IS[[p]] <- numeric(0)
		lis_median_IS[[p]] <- numeric(0)
		#lis_RT_IS[[p]] <- numeric(0)
	}
	unused_profile <- rep(TRUE, dim(profileList_pos[["index_prof"]])[1])
	for(i in 1:length(links_profiles_pos)){
		if(!length(links_profiles_pos[[i]]$IS)) next
		these <- which((links_profiles_pos[[i]]$IS$Counts >= min_count) & (links_profiles_pos[[i]]$IS$max_score >= as.numeric(logfile$parameters$ISnorm_score) ) )
		if(!length(these)) next
		at_profile <- as.numeric(names(links_profiles_pos)[i])
		unused_profile[at_profile] <- FALSE
		if(profileList_pos[["index_prof"]][at_profile,"number_peaks_total"] < min_count) next
		for_peaks <- profileList_pos[["index_prof"]][at_profile,"start_ID"]:profileList_pos[["index_prof"]][at_profile,"end_ID"]
		median_intensity <- median( log10(profileList_pos[["peaks"]][for_peaks, "intensity"]) )
		for(j in for_peaks){
			at_file <- match(profileList_pos[["peaks"]][j,"sampleIDs"], profileList_pos[["sampleID"]])
			at_int <- log10(profileList_pos[["peaks"]][j,"intensity"][[1]])
			#at_RT <- profileList_pos[["peaks"]][j,"RT"]
			lis_delint_IS[[at_file]] <- c(lis_delint_IS[[at_file]], (at_int - median_intensity))
			lis_median_IS[[at_file]] <- c(lis_median_IS[[at_file]], median_intensity)
			#lis_RT_IS[[at_file]] <- c(lis_RT_IS[[at_file]], at_RT)
		}	
	}
	# -> screen other profiles #########################################################
	if( logfile$parameters$ISnorm_medblank == "TRUE" || logfile$parameters$ISnorm_medsam == "TRUE" ){
		lis_delint_nb <- list()
		lis_median_nb <- list()
		lis_delint_b <- list()
		lis_median_b <- list()			
		################################################################################
		if(logfile$parameters$ISnorm_medsam == "TRUE"){
			sam_IDs <- profileList_pos[["sampleID"]][profileList_pos[["type"]] != "blank"]
			for(p in 1:length(sam_IDs)){
				lis_delint_nb[[p]] <- numeric(0)
				lis_median_nb[[p]] <- numeric(0)
			}
			names(lis_delint_nb) <- sam_IDs
			names(lis_median_nb) <- sam_IDs
			ids_nb <- order(profileList_pos[["index_prof"]][,"number_peaks_sample"], decreasing = TRUE) 
			ids_nb <- ids_nb[
				( profileList_pos[["index_prof"]][ids_nb, "number_peaks_sample"] > 0 ) & 
				( profileList_pos[["index_prof"]][ids_nb, "number_peaks_blind"] == 0 )
			]
			if( length(ids_nb) ){
				ids_nb <- ids_nb[unused_profile[ids_nb]]
			}
			if(logfile$parameters$ISnorm_usesubsam == "TRUE" & (length(ids_nb) > as.numeric(logfile$parameters$ISnorm_numsam))){
				ids_nb <- ids_nb[1:as.numeric(logfile$parameters$ISnorm_numsam)]
			}
			if( length(ids_nb) ){
				for(i in ids_nb){
					for_peaks <- profileList_pos[["index_prof"]][i,"start_ID"]:profileList_pos[["index_prof"]][i,"end_ID"]
					median_intensity <- median( log10(profileList_pos[["peaks"]][for_peaks, "intensity"]) )
					for(j in for_peaks){
						at_file <- match(profileList_pos[["peaks"]][j, "sampleIDs"], sam_IDs)
						at_int <- log10(profileList_pos[["peaks"]][j, "intensity"][[1]])
						lis_delint_nb[[at_file]] <- c(lis_delint_nb[[at_file]], (at_int - median_intensity))
						lis_median_nb[[at_file]] <- c(lis_median_nb[[at_file]], median_intensity)
					}	
				}
			}
		}
		################################################################################
		if(logfile$parameters$ISnorm_medblank == "TRUE"){
			blank_IDs <- profileList_pos[["sampleID"]][profileList_pos[["type"]] == "blank"]
			for(p in 1:length(profileList_pos[["sampleID"]])){
				lis_delint_b[[p]] <- numeric(0)
				lis_median_b[[p]] <- numeric(0)
			}		
			names(lis_delint_b) <- blank_IDs 
			names(lis_median_b) <- blank_IDs 
			ids_b <- order(profileList_pos[["index_prof"]][, "number_peaks_blind"], decreasing =TRUE) 
			ids_b <- ids_b[( profileList_pos[["index_prof"]][ids_b, "number_peaks_blind"] > 0 )]
			if( length(ids_b) ){
				ids_b <- ids_b[unused_profile[ids_b]]
			}
			if(logfile$parameters$ISnorm_usesubblank == "TRUE" & (length(ids_b) > as.numeric(logfile$parameters$ISnorm_numblank))){
				ids_b <- ids_b[1:as.numeric(logfile$parameters$ISnorm_numblank)]
			}
			if( length(ids_b) ){
				for(i in ids_b){
					for_peaks <- profileList_pos[["index_prof"]][i,"start_ID"]:profileList_pos[["index_prof"]][i,"end_ID"]
					median_intensity <- median( log10(profileList_pos[["peaks"]][for_peaks, "intensity"]) )
					for(j in for_peaks){
						at_file <- match(profileList_pos[["peaks"]][j, "sampleIDs"], sam_IDs)
						at_int <- log10(profileList_pos[["peaks"]][j, "intensity"][[1]])
						lis_delint_b[[at_file]] <- c(lis_delint_b[[at_file]], (at_int - median_intensity))
						lis_median_b[[at_file]] <- c(lis_median_b[[at_file]], median_intensity)
					}	
				}
			}
		}
		################################################################################
	}
	# -> correct intensities - replace in profileList, recalculate mean_int ############
	corfac <- rep(1, length(lis_delint_IS)) # default: no correction, factor =1
	use_corfac <- rep(FALSE, length(lis_delint_IS))
	for(k in 1:length(lis_delint_IS)){
		if( length(lis_delint_IS[[k]]) >= as.numeric(logfile$parameters$ISnorm_numbIS) ){
			corfac[k] <- (10 ^ median(lis_delint_IS[[k]]))
			use_corfac[k] <- TRUE 
		}
	}
	sampleID <- profileList_pos[["sampleID"]];
	corr_intens <- .Call("_enviMass_correct_intens",
						as.numeric(corfac),	  # correction factor
						as.integer(sampleID),       
						as.numeric(profileList_pos[["peaks"]][,"intensity"]), # intensities
						as.integer(profileList_pos[["peaks"]][,"sampleIDs"]),  
						PACKAGE = "enviMass"
					)
	profileList_pos[["peaks"]][,"intensity"] <- corr_intens
	for(k in 1:dim(profileList_pos[["index_prof"]])[1]){
		profileList_pos[["index_prof"]][k,"mean_int"] <- mean(profileList_pos[["peaks"]][(profileList_pos[["index_prof"]][k,"start_ID"]:profileList_pos[["index_prof"]][k,"end_ID"]),"intensity"])
	}
	# -> data structure to store results for later / intermediate plotting #############
	int_norm_ISTD_pos <- list()
	int_norm_ISTD_pos[[1]] <- lis_delint_IS
	int_norm_ISTD_pos[[2]] <- lis_median_IS
	#int_norm_ISTD_pos[[3]] <- lis_RT_IS
	int_norm_ISTD_pos[[4]] <- use_corfac
	int_norm_ISTD_pos[[5]] <- lis_delint_nb
	int_norm_ISTD_pos[[6]] <- lis_median_nb
	int_norm_ISTD_pos[[7]] <- lis_delint_b
	int_norm_ISTD_pos[[8]] <- lis_median_b
	int_norm_ISTD_pos[[9]] <- profileList_pos[["datetime"]]
	int_norm_ISTD_pos[[10]] <- profileList_pos[["type"]]
	int_norm_ISTD_pos[[11]] <- profileList_pos[["sampleID"]]	
	names(int_norm_ISTD_pos) <- c("lis_delint_IS", " lis_median_IS", "lis_RT_IS", "use_corfac", "lis_delint_nb", 
		"lis_median_nb", "lis_delint_b", "lis_median_b", "atPOSIX", "sampletype", "sampleID")
	# -> save data & derive plots ######################################################
	





	####################################################################################

	# -> save data #####################################################################
	save(profileList_pos, file = file.path(as.character(logfile[[1]]), "results", "profileList_pos"));
	save(links_profiles_pos, file = file.path(as.character(logfile[[1]]), "results", "links_profiles_pos"));	
	save(int_norm_ISTD_pos, file = file.path(as.character(logfile[[1]]), "results", "int_norm_ISTD_pos"));
	rm(links_profiles_pos, profileList_pos, int_norm_ISTD_pos)
	####################################################################################

}

########################################################################################





########################################################################################
# on NEGATIVE profiles #################################################################
if(
	file.exists(file.path(logfile[[1]],"results","profileList_neg"))  
){

	####################################################################################	
	if(any(objects(envir=as.environment(".GlobalEnv"))=="profileList_neg")){rm(profileList_neg,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="profileList_neg")){rm(profileList_neg)}	
	if(any(objects(envir=as.environment(".GlobalEnv"))=="links_peaks_neg")){rm(links_peaks_neg,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="links_peaks_neg")){rm(links_peaks_neg)}				
	if(any(objects(envir=as.environment(".GlobalEnv"))=="links_profiles_neg")){rm(links_profiles_neg,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="links_profiles_neg")){rm(links_profiles_neg)}					
	load(file.path(as.character(logfile[[1]]),"results","profileList_neg"),envir=as.environment(".GlobalEnv"));
	# links_peaks_neg<-list(); # each entry with 6 lists itself: targets, IS, EIC_correl, isotop, adducts, homol - defined in do_profiling.r
	load(file.path(as.character(logfile[[1]]),"results","links_peaks_neg"),envir=as.environment(".GlobalEnv"));	
	links_profiles_neg<-list(); # each entry with 6 lists itself: targets, IS, EIC_correl, isotop, adducts, homol
	measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
	use_entries_profiles<-enviMass::find_empty(links_profiles_neg) # also finds gaps
	profileList_neg[["index_prof"]][,"links"]<-0
	with_bar<-FALSE
	####################################################################################

	if(
		(
			(logfile$workflow[names(logfile$workflow)=="subtr"]=="no") |	
			(
				((logfile$workflow[names(logfile$workflow)=="IS_screen"]=="yes") &  (logfile$parameters$subtr_IS!="yes")) |
				((logfile$workflow[names(logfile$workflow)=="target_screen"]=="yes") & (logfile$parameters$subtr_target!="yes"))	
			) 
		) & 
		(length(links_peaks_neg)>0) # anything screened?
	){

		###########################################################################
		if(with_bar){pBar <- txtProgressBar(min = 0, max = dim(profileList_neg[["index_prof"]])[1], style = 3)}
		for(i in 1:dim(profileList_neg[["index_prof"]])[1]){
			if(with_bar){setTxtProgressBar(pBar, i, title = NULL, label = NULL)}
			if(
				any(profileList_neg[["peaks"]][
					(profileList_neg[["index_prof"]][i,"start_ID"]:profileList_neg[["index_prof"]][i,"end_ID"]),"links"
				]!=0)			
			){
			
				###################################################################
				# add a new link to the profile ###################################
				if( profileList_neg[["index_prof"]][i,"links"]==0 ){ 	# establish a new link ...
					if(length(use_entries_profiles)>0){
						at_entry<-use_entries_profiles[1]
						use_entries<-use_entries_profiles[-1]
					}else{
						at_entry<-(length(links_profiles_neg)+1)
					}
					links_profiles_neg[[at_entry]]<-enviMass::new_entry_links_profiles(profileList_neg[["index_prof"]][i,"number_peaks_total"][[1]])
					names(links_profiles_neg)[at_entry]<-as.character(i)
					profileList_neg[["index_prof"]][i,"links"]<-at_entry						
				}else{
					at_entry<-profileList_neg[["index_prof"]][i,"links"]
				}				
				###################################################################
				# search peaks & their links to compounds #########################
				for(j in (profileList_neg[["index_prof"]][i,"start_ID"]:profileList_neg[["index_prof"]][i,"end_ID"])){
					if(profileList_neg[["peaks"]][j,"links"]!=0){
						# add IS link #############################################
						if(	length(links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[2]]) > 0 ){
							for(k in 1:length(links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[2]])){ # if several compound matches exist for this peak
								if( length(links_profiles_neg[[at_entry]][[2]])==0 ){ # make a new entry for profile link to IS
									links_profiles_neg[[at_entry]][[2]]<-
										data.frame(
											links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[2]][[k]][[1]],
											1,
											links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[2]][[k]][[2]],
											stringsAsFactors = FALSE
										)
									names(links_profiles_neg[[at_entry]][[2]])<-c("Compound","Counts","max_score")
								}else{
									at<-which(links_profiles_neg[[at_entry]][[2]][,1] == links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[2]][[k]][[1]])
									if(length(at)>0){ 	# increment existing link ...
										links_profiles_neg[[at_entry]][[2]][at,2]<-(links_profiles_neg[[at_entry]][[2]][at,2]+1)
										if(links_profiles_neg[[at_entry]][[2]][at,2] < links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[2]][[k]][[2]] ){
											links_profiles_neg[[at_entry]][[2]][at,2] <- links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[2]][[k]][[2]]
										}
									}else{	# ... or add a new one?
										links_profiles_neg[[at_entry]][[2]]<-data.frame(
											c(
												links_profiles_neg[[at_entry]][[2]][,1], 
												links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[2]][[k]][[1]]
											),
											c(links_profiles_neg[[at_entry]][[2]][,2],1),
											c(
												links_profiles_neg[[at_entry]][[2]][,3], 
												links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[2]][[k]][[2]]
											),
											stringsAsFactors = FALSE
										)									
										names(links_profiles_neg[[at_entry]][[2]])<-c("Compound","Counts","max_score")
									}
								}
							}	
						}					
						# add target link #########################################
						if(	length(links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[1]])>0 ){
							for(k in 1:length(links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[1]])){ # if several compound matches exist for this peak
								if( length(links_profiles_neg[[at_entry]][[1]])==0 ){ # make a new entry for profile link to IS
									links_profiles_neg[[at_entry]][[1]]<-
										data.frame(
											links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[1]][[k]][[1]],
											1,
											links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[1]][[k]][[2]],
											stringsAsFactors = FALSE
										)
									names(links_profiles_neg[[at_entry]][[1]])<-c("Compound","Counts","max_score")
								}else{
									at<-which(links_profiles_neg[[at_entry]][[1]][,1] == links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[1]][[k]][[1]])
									if(length(at)>0){ 	# increment existing link ...
										links_profiles_neg[[at_entry]][[1]][at,2]<-(links_profiles_neg[[at_entry]][[1]][at,2]+1)
										if(links_profiles_neg[[at_entry]][[1]][at,2] < links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[1]][[k]][[2]] ){
											links_profiles_neg[[at_entry]][[1]][at,2] <- links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[1]][[k]][[2]]
										}
									}else{	# ... or add a new one?
										links_profiles_neg[[at_entry]][[1]]<-data.frame(
											c(
												links_profiles_neg[[at_entry]][[1]][,1], 
												links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[1]][[k]][[1]]
											),
											c(links_profiles_neg[[at_entry]][[1]][,2],1),
											c(
												links_profiles_neg[[at_entry]][[1]][,3], 
												links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[1]][[k]][[2]]
											),
											stringsAsFactors = FALSE
										)									
										names(links_profiles_neg[[at_entry]][[1]])<-c("Compound","Counts","max_score")
									}
								}
							}					
						}
						###########################################################
					}
				}
				###################################################################
			}
		}
		if(with_bar){close(pBar)}
	}
	if(!length(links_profiles_neg)){
		stop("Not enough internal standard (ISTD) screening matches detected to run the ISTD normalization, positive ionization - remove this workflow step?")
	}
	####################################################################################	

	####################################################################################
	# run normalization ################################################################






	####################################################################################

	####################################################################################
	save(profileList_neg,file=file.path(as.character(logfile[[1]]),"results","profileList_neg"));
	save(links_profiles_neg,file=file.path(as.character(logfile[[1]]),"results","links_profiles_neg"));	
	rm(links_profiles_neg,profileList_neg)
	####################################################################################

}

########################################################################################







