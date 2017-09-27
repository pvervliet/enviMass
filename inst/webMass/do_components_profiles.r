##################################################################################
# Profile componentization #######################################################
##################################################################################

##################################################################################
# REMOVE OLD RESULTS #############################################################
if(file.exists(file.path(as.character(logfile[[1]]),"results","links_profiles_pos"))){file.remove(file.path(as.character(logfile[[1]]),"results","links_profiles_pos"))}
if(file.exists(file.path(as.character(logfile[[1]]),"results","links_profiles_neg"))){file.remove(file.path(as.character(logfile[[1]]),"results","links_profiles_neg"))}
#if(file.exists(file.path(as.character(logfile[[1]]),"results","componentization","components_pos"))){file.remove(file.path(as.character(logfile[[1]]),"results","componentization","components_pos"))}
#if(file.exists(file.path(as.character(logfile[[1]]),"results","componentization","components_neg"))){file.remove(file.path(as.character(logfile[[1]]),"results","componentization","components_neg"))}
##################################################################################



##################################################################################
# POSITIVE #######################################################################
if( 
	file.exists(file.path(logfile[[1]],"results","profileList_pos"))  
){

	##############################################################################	
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
	peaks<-profileList_pos[["peaks"]][,c("sampleIDs","peakIDs","profileIDs","RT")] # to retrieve relations with, sampleID, peakID, profileID, RT
	ord<-order(peaks[,"sampleIDs"],peaks[,"peakIDs"],peaks[,"profileIDs"],decreasing=FALSE)
	peaks<-peaks[ord,]
	use_entries_profiles<-enviMass::find_empty(links_profiles_pos) # also finds gaps
	profileList_pos[["index_prof"]][,"links"]<-0
	with_bar<-TRUE
	##############################################################################	

	##############################################################################	
	# (1) ANNOTATE TARGET & ISTD SCREENING MACTHES stored in links_peaks_pos #####
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
		cat("\n Annotation of screening results to profiles")
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
						if(	length(links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[2]])>0 ){
							if(any(duplicated(links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[1]]))){ # From different peak vs. isotopologue pattern combinations
								links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[2]]<-
									links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[2]][!duplicated(links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[2]])]
							}		
							for(k in 1:length(links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[2]])){ # if several compound matches exist for this peak
								if( length(links_profiles_pos[[at_entry]][[2]])==0 ){ # make a new entry for profile link to IS
									links_profiles_pos[[at_entry]][[2]]<-
										data.frame(
											links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[2]][[k]],1,stringsAsFactors = FALSE
										)
									names(links_profiles_pos[[at_entry]][[2]])<-c("Compound","Counts")
								}else{
									at<-which(links_profiles_pos[[at_entry]][[2]][,1]==links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[2]][[k]])
									if(length(at)>0){ 	# increment existing link ...
										links_profiles_pos[[at_entry]][[2]][at,2]<-(links_profiles_pos[[at_entry]][[2]][at,2]+1)
									}else{	# or just add a new one?
										links_profiles_pos[[at_entry]][[2]]<-data.frame(
											c(links_profiles_pos[[at_entry]][[2]][,1], links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[2]][[k]]),
											c(links_profiles_pos[[at_entry]][[2]][,2],1),
											stringsAsFactors = FALSE
										)									
										names(links_profiles_pos[[at_entry]][[2]])<-c("Compound","Counts")
									}
								}
							}	
						}					
						# add target link #########################################
						if(	length(links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[1]])>0 ){
							if(any(duplicated(links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[1]]))){ # From different peak vs. isotopologue pattern combinations
								links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[1]]<-
									links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[1]][!duplicated(links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[1]])]
							}						
							for(k in 1:length(links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[1]])){ # if several compound matches exist for this peak
								if( length(links_profiles_pos[[at_entry]][[1]])==0 ){ # make a new entry for profile link to IS
									links_profiles_pos[[at_entry]][[1]]<-
										data.frame(
											links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[1]][[k]],1,stringsAsFactors = FALSE
										)
									names(links_profiles_pos[[at_entry]][[1]])<-c("Compound","Counts")
								}else{
									at<-which(links_profiles_pos[[at_entry]][[1]][,1]==links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[1]][[k]])
									if(length(at)>0){ 	# increment existing link ...
										links_profiles_pos[[at_entry]][[1]][at,2]<-(links_profiles_pos[[at_entry]][[1]][at,2]+1)
									}else{	# or just add a new one?
										links_profiles_pos[[at_entry]][[1]]<-data.frame(
											c(links_profiles_pos[[at_entry]][[1]][,1], links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[1]][[k]]),
											c(links_profiles_pos[[at_entry]][[1]][,2],1),
											stringsAsFactors = FALSE
										)									
										names(links_profiles_pos[[at_entry]][[1]])<-c("Compound","Counts")
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
		cat(" done.")
	}
	##############################################################################	
	
	##############################################################################	
	# for each profile, filter all non-target relations to other profiles ########
	# (2) INSERT EIC and CO-OCCURRENCE INFORMATION ###############################
	if(logfile$workflow[names(logfile$workflow)=="EIC_correlation"]=="yes"){
		cat("\n Retrieving EIC correlation links ")
		# (2.1) INSERT EIC LINKS #################################################
		##########################################################################
		forIDs<-profileList_pos[["sampleID"]]
		for_files<-list.files(file.path(logfile[[1]],"results","componentization","EIC_corr"))
		keep<-match(forIDs,for_files) # which files are available?
		if(any(is.na(keep))){cat("\n Just note: not all files found in profiles have EIC correlation results (1). \n")}
		TRUE_IDs<-(measurements$ID[measurements$adducts=="TRUE"]) # for files which have run through that step
		keep2<-match(forIDs,for_files) # ensure file availability
		forIDs<-forIDs[!is.na(keep) & !is.na(keep2)]
		if(logfile$parameters$dofile_latest_profcomp=="TRUE"){ # restrict to latest files
			atPOSIX<-profileList_pos[["datetime"]];
			matchID<-profileList_pos[["sampleID"]];
			atdate<-c();attime<-c();
			for(i in 1:length(atPOSIX)){
				atdate<-c(atdate, strsplit(atPOSIX[i]," ")[[1]][1]);
				attime<-c(attime, strsplit(atPOSIX[i]," ")[[1]][2]);
			}
			attime<-as.difftime(attime);
			atdate<-as.Date(atdate, tz="GMT");
			ord<-order(as.numeric(atdate),as.numeric(attime),matchID,decreasing=TRUE);
			matchID<-matchID[ord];
			forIDs<-forIDs[match(forIDs,matchID)]
			if(length(forIDs)>as.numeric(logfile$parameters$numfile_latest_profcomp)){
				forIDs<-forIDs[1:as.numeric(logfile$parameters$numfile_latest_profcomp)]
			}
		}
		not_found1<-0;inserted1<-0
		if(length(forIDs)>0){
			if(with_bar){pBar <- txtProgressBar(min = 0, max = length(forIDs), style = 3)}
			for(i in 1:length(forIDs)){
				if(with_bar){setTxtProgressBar(pBar, i, title = NULL, label = NULL)}
				load(file=file.path(logfile[[1]],"results","componentization","EIC_corr",forIDs[i]))
				#
				#EIC_pairs<-EIC_pairs[EIC_pairs[,4]>=logfile$parameters$EICor_mincor,,drop=FALSE]		
				#
				if(length(EIC_pairs[,1])==0){next}
				# find profiles for first peak
				get1<-cbind(
					rep(as.numeric(forIDs[i]),length(EIC_pairs[,1])),EIC_pairs[,1] # first column sorted correctly
				)
				found1<-enviMass::rows_compare(get1,peaks[,c("sampleIDs","peakIDs")],row_order=FALSE,column_order_a=FALSE,column_order_b=FALSE,get_index=TRUE)
				# find profiles for second peak
				get2<-cbind(
					rep(as.numeric(forIDs[i]),length(EIC_pairs[,2])),EIC_pairs[,2] # second column requires sorting
				)
				found2<-enviMass::rows_compare(get2,peaks[,c("sampleIDs","peakIDs")],row_order=FALSE,column_order_a=TRUE,column_order_b=FALSE,get_index=TRUE)
				for(j in 1:length(found1)){ # insert links
					# insert PROFILE LINKS ########################################
					if(found1[j]==0){not_found1<-(not_found1+1);next} # e.g., peak blind-removed 
					if(found2[j]==0){not_found1<-(not_found1+1);next}		
					inserted1<-(inserted1+1);				
					# (1) insert link to second profile for the first profile
					prof1<-peaks[found1[j],"profileIDs"][[1]]
					prof2<-peaks[found2[j],"profileIDs"][[1]]
					# enable check ... pairs musst have very similar retention times
					if(FALSE){
						cat("\n");
						cat(peaks[found1[j],"RT"]);cat(" - ")
						cat(peaks[found2[j],"RT"])
					}
					if(profileList_pos[["index_prof"]][prof1,"profile_ID"]!=prof1){stop("\nComponentization: debug me, #1!")}				
					if(profileList_pos[["index_prof"]][prof1,"links"]==0){ 	# establish a new link ...
						if(length(use_entries_profiles)>0){
							at_entry_1<-use_entries_profiles[1]
							use_entries<-use_entries_profiles[-1]
						}else{
							at_entry_1<-(length(links_profiles_pos)+1)
						}
						links_profiles_pos[[at_entry_1]]<-enviMass::new_entry_links_profiles(profileList_pos[["index_prof"]][prof1,"number_peaks_total"][[1]])
						names(links_profiles_pos)[at_entry_1]<-as.character(prof1)
						profileList_pos[["index_prof"]][prof1,"links"]<-at_entry_1						
					}else{
						at_entry_1<-profileList_pos[["index_prof"]][prof1,"links"]
					}
					here1<-which(links_profiles_pos[[at_entry_1]][["EIC"]][,"linked profile"]==prof2)
					if(length(here1)==0){			
						links_profiles_pos[[at_entry_1]][["EIC"]]<-rbind(
							links_profiles_pos[[at_entry_1]][["EIC"]], c(prof2,1,0,0,0,0,1,NA)
						)
						here1<-dim(links_profiles_pos[[at_entry_1]][["EIC"]])[1]
					}else{
						links_profiles_pos[[at_entry_1]][["EIC"]][here1,"link counts"]<-(links_profiles_pos[[at_entry_1]][["EIC"]][here1,"link counts"]+1)
					}
					# (2) insert link to first profile for the second profile
					if(profileList_pos[["index_prof"]][prof2,"profile_ID"]!=prof2){stop("\nCross-profile componentization: debug me, EIC_1!")}				
					if(profileList_pos[["index_prof"]][prof2,"links"]==0){ 	# establish a new link ...
						if(length(use_entries_profiles)>0){
							at_entry_2<-use_entries_profiles[1]
							use_entries<-use_entries_profiles[-1]
						}else{
							at_entry_2<-(length(links_profiles_pos)+1)
						}
						links_profiles_pos[[at_entry_2]]<-enviMass::new_entry_links_profiles(profileList_pos[["index_prof"]][prof2,"number_peaks_total"][[1]])		
						names(links_profiles_pos)[at_entry_2]<-as.character(prof2)
						profileList_pos[["index_prof"]][prof2,"links"]<-at_entry_2						
					}else{
						at_entry_2<-profileList_pos[["index_prof"]][prof2,"links"]
					}
					here2<-which(links_profiles_pos[[at_entry_2]][["EIC"]][,"linked profile"]==prof1)
					if(length(here2)==0){
						links_profiles_pos[[at_entry_2]][["EIC"]]<-rbind(
							links_profiles_pos[[at_entry_2]][["EIC"]], c(prof1,1,0,0,0,0,0,NA)
						)
						here2<-dim(links_profiles_pos[[at_entry_2]][["EIC"]])[1]
					}else{
						links_profiles_pos[[at_entry_2]][["EIC"]][here2,"link counts"]<-(links_profiles_pos[[at_entry_2]][["EIC"]][here2,"link counts"]+1)
					}				
					# insert EIC correlation
					if(links_profiles_pos[[at_entry_1]][["EIC"]][here1,"use"]==1){
						if(length(links_profiles_pos[[at_entry_1]]$EIC_cor)<here1){
							links_profiles_pos[[at_entry_1]]$EIC_cor[[here1]]<-EIC_pairs[j,4]
						}else{
							if(is.null(links_profiles_pos[[at_entry_1]]$EIC_cor[[here1]])){
								links_profiles_pos[[at_entry_1]]$EIC_cor[[here1]]<-EIC_pairs[j,4]
							}else{
								links_profiles_pos[[at_entry_1]]$EIC_cor[[here1]]<-c(links_profiles_pos[[at_entry_1]]$EIC_cor[[here1]],EIC_pairs[j,4])
							}
						}
					}else{
						if(length(links_profiles_pos[[at_entry_2]]$EIC_cor)<here2){
							links_profiles_pos[[at_entry_2]]$EIC_cor[[here2]]<-EIC_pairs[j,4]
						}else{
							if(is.null(links_profiles_pos[[at_entry_2]]$EIC_cor[[here2]])){
								links_profiles_pos[[at_entry_2]]$EIC_cor[[here2]]<-EIC_pairs[j,4]
							}else{
								links_profiles_pos[[at_entry_2]]$EIC_cor[[here2]]<-c(links_profiles_pos[[at_entry_2]]$EIC_cor[[here2]],EIC_pairs[j,4])
							}
						}
					}
					# insert PEAK LINKS ##########################################
				}
			}
			if(with_bar){close(pBar)}		
			if(any(objects()=="EIC_pairs")){rm(EIC_pairs)}	
			cat(" done.")
		}
	}
	# (3) INSERT ISOTOPOLOGUE LINKS ##############################################
	if(logfile$workflow[names(logfile$workflow)=="isotopologues"]=="yes"){	
		forIDs<-profileList_pos[["sampleID"]]
		for_files<-list.files(file.path(logfile[[1]],"results","componentization","isotopologues"))
		keep<-match(forIDs,for_files) # which files are available?
		if(any(is.na(keep))){cat("\n Just note: not all files found in profiles have isotopologue links available. \n")}
		TRUE_IDs<-(measurements$ID[measurements$adducts=="TRUE"]) # for files which have run through that step
		keep2<-match(forIDs,for_files)
		forIDs<-forIDs[!is.na(keep) & !is.na(keep2)]
		if(logfile$parameters$dofile_latest_profcomp=="TRUE"){ # restrict to latest files
			atPOSIX<-profileList_pos[["datetime"]];
			matchID<-profileList_pos[["sampleID"]];
			atdate<-c();attime<-c();
			for(i in 1:length(atPOSIX)){
				atdate<-c(atdate, strsplit(atPOSIX[i]," ")[[1]][1]);
				attime<-c(attime, strsplit(atPOSIX[i]," ")[[1]][2]);
			}
			attime<-as.difftime(attime);
			atdate<-as.Date(atdate, tz="GMT");
			ord<-order(as.numeric(atdate),as.numeric(attime),matchID,decreasing=TRUE);
			matchID<-matchID[ord];
			forIDs<-forIDs[match(forIDs,matchID)]
			if(length(forIDs)>as.numeric(logfile$parameters$numfile_latest_profcomp)){
				forIDs<-forIDs[1:as.numeric(logfile$parameters$numfile_latest_profcomp)]
			}
		}
		not_found3<-0;inserted3<-0
		if(length(forIDs)>0){
			cat("\n Retrieving isotopologue links ")
			if(with_bar){pBar <- txtProgressBar(min = 0, max = length(forIDs), style = 3)}
			for(i in 1:length(forIDs)){
				if(with_bar){setTxtProgressBar(pBar, i, title = NULL, label = NULL)}
				load(file=file.path(logfile[[1]],"results","componentization","isotopologues",forIDs[i]))
				if(length(Isot_pairs[,1])==0){next}
				# find profiles for first peak
				get1<-cbind(
					rep(as.numeric(forIDs[i]),length(Isot_pairs[,1])),Isot_pairs[,1]
				)
				found1<-enviMass::rows_compare(get1,peaks[,c("sampleIDs","peakIDs")],row_order=FALSE,column_order_a=FALSE,column_order_b=FALSE,get_index=TRUE)
				# find profiles for second peak
				get2<-cbind(
					rep(as.numeric(forIDs[i]),length(Isot_pairs[,2])),Isot_pairs[,2]
				)
				found2<-enviMass::rows_compare(get2,peaks[,c("sampleIDs","peakIDs")],row_order=FALSE,column_order_a=TRUE,column_order_b=FALSE,get_index=TRUE)
				for(j in 1:length(found1)){ # insert links
					# insert PROFILE LINKS ########################################
					if(found1[j]==0){not_found3<-(not_found3+1);next} # e.g., peak blind-removed 
					if(found2[j]==0){not_found3<-(not_found3+1);next}			
					inserted3<-(inserted3+1);				
					# enable check ... pairs musst have very similar retention times
					if(FALSE){
						cat("\n");
						cat(peaks[found1[j],"RT"]);cat(" - ")
						cat(peaks[found2[j],"RT"])				
					}
					# (1) insert link to second profile for the first profile
					prof1<-peaks[found1[j],"profileIDs"][[1]]
					prof2<-peaks[found2[j],"profileIDs"][[1]]				
					if(profileList_pos[["index_prof"]][prof1,"profile_ID"]!=prof1){stop("\nComponentization: debug me, #1!")}				
					if(profileList_pos[["index_prof"]][prof1,"links"]==0){ 	# establish a new link ...
						if(length(use_entries_profiles)>0){
							at_entry_1<-use_entries_profiles[1]
							use_entries<-use_entries_profiles[-1]
						}else{
							at_entry_1<-(length(links_profiles_pos)+1)
						}
						links_profiles_pos[[at_entry_1]]<-enviMass::new_entry_links_profiles(profileList_pos[["index_prof"]][prof1,"number_peaks_total"][[1]])
						names(links_profiles_pos)[at_entry_1]<-as.character(prof1)
						profileList_pos[["index_prof"]][prof1,"links"]<-at_entry_1						
					}else{
						at_entry_1<-profileList_pos[["index_prof"]][prof1,"links"]
					}
					here1<-which(links_profiles_pos[[at_entry_1]][["isot"]][,"linked profile"]==prof2)
					if(length(here1)==0){
						links_profiles_pos[[at_entry_1]][["isot"]]<-rbind(
							links_profiles_pos[[at_entry_1]][["isot"]], c(prof2,1,0,1,NA)
						)
						here1<-dim(links_profiles_pos[[at_entry_1]][["isot"]])[1]
					}else{
						links_profiles_pos[[at_entry_1]][["isot"]][here1,"link counts"]<-(links_profiles_pos[[at_entry_1]][["isot"]][here1,"link counts"]+1)
					}
					# (2) insert link to first profile for the second profile
					if(profileList_pos[["index_prof"]][prof2,"profile_ID"]!=prof2){stop("\nComponentization: debug me, #1!")}				
					if(profileList_pos[["index_prof"]][prof2,"links"]==0){ 	# establish a new link ...
						if(length(use_entries_profiles)>0){
							at_entry_2<-use_entries_profiles[1]
							use_entries<-use_entries_profiles[-1]
						}else{																		
							at_entry_2<-(length(links_profiles_pos)+1)
						}
						links_profiles_pos[[at_entry_2]]<-enviMass::new_entry_links_profiles(profileList_pos[["index_prof"]][prof2,"number_peaks_total"][[1]])
						names(links_profiles_pos)[at_entry_2]<-as.character(prof2)
						profileList_pos[["index_prof"]][prof2,"links"]<-at_entry_2						
					}else{
						at_entry_2<-profileList_pos[["index_prof"]][prof2,"links"]
					}
					here2<-which(links_profiles_pos[[at_entry_2]][["isot"]][,"linked profile"]==prof1)
					if(length(here2)==0){
						links_profiles_pos[[at_entry_2]][["isot"]]<-rbind(
							links_profiles_pos[[at_entry_2]][["isot"]], c(prof1,1,0,0,NA)
						)
						here2<-dim(links_profiles_pos[[at_entry_2]][["isot"]])[1]
					}else{
						links_profiles_pos[[at_entry_2]][["isot"]][here2,"link counts"]<-(links_profiles_pos[[at_entry_2]][["isot"]][here2,"link counts"]+1)
					}				
					# INSERT ref_1: total number of co-occurences ################
					got_entry<-FALSE
					#if(dim(links_profiles_pos[[at_entry_1]]$EIC)[1]!=0){ # already dealed with during EIC correlation?
					#	here3<-which(links_profiles_pos[[at_entry_1]]$EIC[,"linked profile"]==prof2)
					#	if(length(here3)>0){
					#		got_entry<-TRUE
					#		if(length(here3)>1){
					#			stop("\n DEBUG ME !")
					#		}
					#	}
					#}
					if(got_entry){ # use value of existing entry ...
						links_profiles_pos[[at_entry_1]]$isot[here1,"ref_1"]<-links_profiles_pos[[at_entry_1]]$EIC[here3,"ref_1"]
						links_profiles_pos[[at_entry_2]]$isot[here2,"ref_1"]<-links_profiles_pos[[at_entry_1]]$EIC[here3,"ref_1"]
					}else{ # ... or retrieve them
						these<-profileList_pos[["peaks"]][
							profileList_pos[["index_prof"]][prof1,"start_ID"]:profileList_pos[["index_prof"]][prof1,"end_ID"]
						,"sampleIDs"]
						those<-profileList_pos[["peaks"]][
							profileList_pos[["index_prof"]][prof2,"start_ID"]:profileList_pos[["index_prof"]][prof2,"end_ID"]
						,"sampleIDs"]
						matched<-match(these,those)
						not_NA<-sum(!is.na(matched))
						links_profiles_pos[[at_entry_1]]$isot[here1,"ref_1"]<-not_NA
						links_profiles_pos[[at_entry_2]]$isot[here2,"ref_1"]<-not_NA
					}
					if(any(links_profiles_pos[[at_entry_1]]$isot[,"ref_1"]==0)){stop("\n DEBUG ME ! FOUND_1")}
					if(any(links_profiles_pos[[at_entry_2]]$isot[,"ref_1"]==0)){stop("\n DEBUG ME ! FOUND_2")}
					# insert PEAK LINKS ##########################################
				}
			}
			if(with_bar){close(pBar)}
			if(any(objects()=="Isot_pairs")){rm(Isot_pairs)}
			cat(" done.")
		}
	}
	# (4) INSERT ADDUCT LINKS ####################################################
	if(logfile$workflow[names(logfile$workflow)=="adducts"]=="yes"){
		forIDs<-profileList_pos[["sampleID"]]
		for_files<-list.files(file.path(logfile[[1]],"results","componentization","adducts"))
		keep<-match(forIDs,for_files) # which files are available?
		if(any(is.na(keep))){cat("\n Just note: not all files found in profiles have adduct links available. \n")}
		TRUE_IDs<-(measurements$ID[measurements$adducts=="TRUE"]) # for files which have run through that step
		keep2<-match(forIDs,for_files)
		forIDs<-forIDs[!is.na(keep) & !is.na(keep2)]
		if(logfile$parameters$dofile_latest_profcomp=="TRUE"){ # restrict to latest files
			atPOSIX<-profileList_pos[["datetime"]];
			matchID<-profileList_pos[["sampleID"]];
			atdate<-c();attime<-c();
			for(i in 1:length(atPOSIX)){
				atdate<-c(atdate, strsplit(atPOSIX[i]," ")[[1]][1]);
				attime<-c(attime, strsplit(atPOSIX[i]," ")[[1]][2]);
			}
			attime<-as.difftime(attime);
			atdate<-as.Date(atdate, tz="GMT");
			ord<-order(as.numeric(atdate),as.numeric(attime),matchID,decreasing=TRUE);
			matchID<-matchID[ord];
			forIDs<-forIDs[match(forIDs,matchID)]
			if(length(forIDs)>as.numeric(logfile$parameters$numfile_latest_profcomp)){
				forIDs<-forIDs[1:as.numeric(logfile$parameters$numfile_latest_profcomp)]
			}
		}
		not_found4<-0;inserted4<-0
		if(length(forIDs)>0){
			cat("\n Retrieving adduct links ")
			if(with_bar){pBar <- txtProgressBar(min = 0, max = length(forIDs), style = 3)}
			for(i in 1:length(forIDs)){
				if(with_bar){setTxtProgressBar(pBar, i, title = NULL, label = NULL)}
				load(file=file.path(logfile[[1]],"results","componentization","adducts",forIDs[i]))
				if(length(Adduct_pairs[,1])==0){next}
				# find profiles for first peak
				get1<-cbind(
					rep(as.numeric(forIDs[i]),length(Adduct_pairs[,1])),Adduct_pairs[,1]
				)
				found1<-enviMass::rows_compare(get1,peaks[,c("sampleIDs","peakIDs")],row_order=FALSE,column_order_a=FALSE,column_order_b=FALSE,get_index=TRUE)
				# find profiles for second peak
				get2<-cbind(
					rep(as.numeric(forIDs[i]),length(Adduct_pairs[,2])),Adduct_pairs[,2]
				)
				found2<-enviMass::rows_compare(get2,peaks[,c("sampleIDs","peakIDs")],row_order=FALSE,column_order_a=TRUE,column_order_b=FALSE,get_index=TRUE)
				for(j in 1:length(found1)){ # insert links	
					if(found1[j]==0){not_found4<-(not_found4+1);next} # e.g., peak blind-removed 
					if(found2[j]==0){not_found4<-(not_found4+1);next}		
					inserted4<-(inserted4+1);
					# enable check ... pairs musst have very similar retention times
					if(FALSE){
						cat("\n");
						cat(peaks[found1[j],"RT"]);cat(" - ")
						cat(peaks[found2[j],"RT"])				
					}
					# insert PROFILE LINKS ######################################	
					# (1) insert link to second profile for the first profile
					prof1<-peaks[found1[j],"profileIDs"][[1]]
					prof2<-peaks[found2[j],"profileIDs"][[1]]
					if(profileList_pos[["index_prof"]][prof1,"profile_ID"]!=prof1){stop("\nComponentization: debug me, #1!")}				
					if(profileList_pos[["index_prof"]][prof1,"links"]==0){ 	# establish a new link ...
						if(length(use_entries_profiles)>0){
							at_entry_1<-use_entries_profiles[1]
							use_entries<-use_entries_profiles[-1]
						}else{
							at_entry_1<-(length(links_profiles_pos)+1)
						}
						links_profiles_pos[[at_entry_1]]<-enviMass::new_entry_links_profiles(profileList_pos[["index_prof"]][prof1,"number_peaks_total"][[1]])
						names(links_profiles_pos)[at_entry_1]<-as.character(prof1)
						profileList_pos[["index_prof"]][prof1,"links"]<-at_entry_1						
					}else{
						at_entry_1<-profileList_pos[["index_prof"]][prof1,"links"]
					}
					here1<-which(links_profiles_pos[[at_entry_1]][["adduc"]][,"linked profile"]==prof2)
					if(length(here1)==0){
						links_profiles_pos[[at_entry_1]][["adduc"]]<-rbind(
							links_profiles_pos[[at_entry_1]][["adduc"]], c(prof2,1,0,1,NA)
						)
						here1<-dim(links_profiles_pos[[at_entry_1]][["adduc"]])[1]
					}else{
						links_profiles_pos[[at_entry_1]][["adduc"]][here1,"link counts"]<-(links_profiles_pos[[at_entry_1]][["adduc"]][here1,"link counts"]+1)
					}
					# (2) insert link to first profile for the second profile
					if(profileList_pos[["index_prof"]][prof2,"profile_ID"]!=prof2){stop("\nComponentization: debug me, #1!")}				
					if(profileList_pos[["index_prof"]][prof2,"links"]==0){ 	# establish a new link ...
						if(length(use_entries_profiles)>0){
							at_entry_2<-use_entries_profiles[1]
							use_entries<-use_entries_profiles[-1]
						}else{
							at_entry_2<-(length(links_profiles_pos)+1)
						}
						links_profiles_pos[[at_entry_2]]<-enviMass::new_entry_links_profiles(profileList_pos[["index_prof"]][prof2,"number_peaks_total"][[1]])
						names(links_profiles_pos)[at_entry_2]<-as.character(prof2)
						profileList_pos[["index_prof"]][prof2,"links"]<-at_entry_2						
					}else{
						at_entry_2<-profileList_pos[["index_prof"]][prof2,"links"]
					}
					here2<-which(links_profiles_pos[[at_entry_2]][["adduc"]][,"linked profile"]==prof1)
					if(length(here2)==0){
						links_profiles_pos[[at_entry_2]][["adduc"]]<-rbind(
							links_profiles_pos[[at_entry_2]][["adduc"]], c(prof1,1,0,0,NA)
						)
						here2<-dim(links_profiles_pos[[at_entry_2]][["adduc"]])[1]
					}else{
						links_profiles_pos[[at_entry_2]][["adduc"]][here2,"link counts"]<-(links_profiles_pos[[at_entry_2]][["adduc"]][here2,"link counts"]+1)
					}		
					# INSERT ref_1: total number of co-occurences ################
					got_entry<-FALSE
					#if(dim(links_profiles_pos[[at_entry_1]]$EIC)[1]!=0){ # already dealed with during EIC correlation?
					#	here3<-which(links_profiles_pos[[at_entry_1]]$EIC[,"linked profile"]==prof2)
					#	if(length(here3)>0){
					#		got_entry<-TRUE
					#	}
					#}
					if(got_entry){ # use value of existing entry ...
						links_profiles_pos[[at_entry_1]]$adduc[here1,"ref_1"]<-links_profiles_pos[[at_entry_1]]$EIC[here3,"ref_1"]
						links_profiles_pos[[at_entry_2]]$adduc[here2,"ref_1"]<-links_profiles_pos[[at_entry_1]]$EIC[here3,"ref_1"]
					}else{ # ... or retrieve them
						these<-profileList_pos[["peaks"]][
							profileList_pos[["index_prof"]][prof1,"start_ID"]:profileList_pos[["index_prof"]][prof1,"end_ID"]
						,"sampleIDs"]
						those<-profileList_pos[["peaks"]][
							profileList_pos[["index_prof"]][prof2,"start_ID"]:profileList_pos[["index_prof"]][prof2,"end_ID"]
						,"sampleIDs"]
						matched<-match(these,those)
						not_NA<-sum(!is.na(matched))
						links_profiles_pos[[at_entry_1]]$adduc[here1,"ref_1"]<-not_NA
						links_profiles_pos[[at_entry_2]]$adduc[here2,"ref_1"]<-not_NA
					}
					if(any(links_profiles_pos[[at_entry_1]]$adduc[,"ref_1"]==0)){stop("\n DEBUG ME ! FOUND_3")}
					if(any(links_profiles_pos[[at_entry_2]]$adduc[,"ref_1"]==0)){stop("\n DEBUG ME ! FOUND_4")}
					# insert PEAK LINKS? ########################################
				}
			}
			if(with_bar){close(pBar)}
			if(any(objects()=="Adduct_pairs")){rm(Adduct_pairs)}
			cat(" done.")
		}	
	}
	# (5) INSERT HOMOLOGUE SERIES LINKS ##########################################
	if(logfile$workflow[names(logfile$workflow)=="homologues"]=="yes"){
		forIDs<-profileList_pos[["sampleID"]]
		for_files<-list.files(file.path(logfile[[1]],"results","componentization","homologues"))
		keep<-match(forIDs,for_files) # which files are available?
		if(any(is.na(keep))){cat("\n Just note: not all files found in profiles have homologue links available. \n")}
		TRUE_IDs<-(measurements$ID[measurements$homologues=="TRUE"]) # for files which have run through that step
		keep2<-match(forIDs,for_files)
		forIDs<-forIDs[!is.na(keep) & !is.na(keep2)]
		if(logfile$parameters$dofile_latest_profcomp=="TRUE"){ # restrict to latest files
			atPOSIX<-profileList_pos[["datetime"]];
			matchID<-profileList_pos[["sampleID"]];
			atdate<-c();attime<-c();
			for(i in 1:length(atPOSIX)){
				atdate<-c(atdate, strsplit(atPOSIX[i]," ")[[1]][1]);
				attime<-c(attime, strsplit(atPOSIX[i]," ")[[1]][2]);
			}
			attime<-as.difftime(attime);
			atdate<-as.Date(atdate, tz="GMT");
			ord<-order(as.numeric(atdate),as.numeric(attime),matchID,decreasing=TRUE);
			matchID<-matchID[ord];
			forIDs<-forIDs[match(forIDs,matchID)]
			if(length(forIDs)>as.numeric(logfile$parameters$numfile_latest_profcomp)){
				forIDs<-forIDs[1:as.numeric(logfile$parameters$numfile_latest_profcomp)]
			}
		}
		not_found5<-0;inserted5<-0
		if(length(forIDs)>0){
			cat("\n Retrieving homologue links ")
			if(with_bar){pBar <- txtProgressBar(min = 0, max = length(forIDs), style = 3)}
			for(i in 1:length(forIDs)){
				if(with_bar){setTxtProgressBar(pBar, i, title = NULL, label = NULL)}
				load(file=file.path(logfile[[1]],"results","componentization","homologues",forIDs[i]))	
				if(length(Homol_groups[,1])==0){next}	
				keep_group<-rep(TRUE,max(Homol_groups[,3])) # make sure all peaks in a group are still present
				# find profiles for first peak
				get1<-cbind(
					rep(as.numeric(forIDs[i]),length(Homol_groups[,1])),Homol_groups[,1]
				)
				found1<-enviMass::rows_compare(get1,peaks[,c("sampleIDs","peakIDs")],row_order=FALSE,column_order_a=FALSE,column_order_b=FALSE,get_index=TRUE)
				keep_group[Homol_groups[found1==0,3]]<-FALSE
				if(!any(keep_group)){next} # none left with all peaks still embedded into profiles
				# find profiles for second peak
				get2<-cbind(
					rep(as.numeric(forIDs[i]),length(Homol_groups[,2])),Homol_groups[,2]
				)
				found2<-enviMass::rows_compare(get2,peaks[,c("sampleIDs","peakIDs")],row_order=FALSE,column_order_a=TRUE,column_order_b=FALSE,get_index=TRUE)
				if(!any(keep_group)){next} # none left with all peaks still embedded into profiles
				for(j in 1:length(found1)){ # insert links	
					if(found1[j]==0){not_found5<-(not_found5+1);next} # e.g., peak blind-removed 
					if(found2[j]==0){not_found5<-(not_found5+1);next}
					if(Homol_groups[j,3]!=Homol_groups[j,3]){stop("\n Debug componentization for its homologue relations!")}
					if(!keep_group[Homol_groups[j,3]]){next} # not all peaks for this group were embedded into profiles, e.g., peak blind-removed
					inserted5<-(inserted5+1);
					# enable check ... pairs musst have very similar retention times
					if(FALSE){
						cat("\n");
						cat(peaks[found1[j],"RT"]);cat(" - ")
						cat(peaks[found2[j],"RT"])				
					}
					# insert PROFILE LINKS #############################################	
					# (1) insert link to second profile for the first profile
					prof1<-peaks[found1[j],"profileIDs"][[1]]
					prof2<-peaks[found2[j],"profileIDs"][[1]]
					if(profileList_pos[["index_prof"]][prof1,"profile_ID"]!=prof1){stop("\nComponentization: debug me, #1!")}				
					if(profileList_pos[["index_prof"]][prof1,"links"]==0){ 	# establish a new link ...
						if(length(use_entries_profiles)>0){
							at_entry<-use_entries_profiles[1]
							use_entries<-use_entries_profiles[-1]
						}else{
							at_entry<-(length(links_profiles_pos)+1)
						}
						links_profiles_pos[[at_entry]]<-enviMass::new_entry_links_profiles(profileList_pos[["index_prof"]][prof1,"number_peaks_total"][[1]])
						names(links_profiles_pos)[at_entry]<-as.character(prof1)
						profileList_pos[["index_prof"]][prof1,"links"]<-at_entry						
					}else{
						at_entry<-profileList_pos[["index_prof"]][prof1,"links"]
					}
					here<-which(links_profiles_pos[[at_entry]][["homol"]][,"linked profile"]==prof2)
					if(length(here)==0){
						links_profiles_pos[[at_entry]][["homol"]]<-rbind(
							links_profiles_pos[[at_entry]][["homol"]], c(prof2,1,0)
						)
					}else{
						links_profiles_pos[[at_entry]][["homol"]][here,"link counts"]<-(links_profiles_pos[[at_entry]][["homol"]][here,"link counts"]+1)
					}
					# (2) insert link to first profile for the second profile
					if(profileList_pos[["index_prof"]][prof2,"profile_ID"]!=prof2){stop("\nComponentization: debug me, #1!")}				
					if(profileList_pos[["index_prof"]][prof2,"links"]==0){ 	# establish a new link ...
						if(length(use_entries_profiles)>0){
							at_entry<-use_entries_profiles[1]
							use_entries<-use_entries_profiles[-1]
						}else{
							at_entry<-(length(links_profiles_pos)+1)
						}
						links_profiles_pos[[at_entry]]<-enviMass::new_entry_links_profiles(profileList_pos[["index_prof"]][prof2,"number_peaks_total"][[1]])
						names(links_profiles_pos)[at_entry]<-as.character(prof2)
						profileList_pos[["index_prof"]][prof2,"links"]<-at_entry						
					}else{
						at_entry<-profileList_pos[["index_prof"]][prof2,"links"]
					}
					here<-which(links_profiles_pos[[at_entry]][["homol"]][,"linked profile"]==prof1)
					if(length(here)==0){
						links_profiles_pos[[at_entry]][["homol"]]<-rbind(
							links_profiles_pos[[at_entry]][["homol"]], c(prof1,1,0)
						)
					}else{
						links_profiles_pos[[at_entry]][["homol"]][here,"link counts"]<-(links_profiles_pos[[at_entry]][["homol"]][here,"link counts"]+1)
					}				
				}	
			}
			if(with_bar){close(pBar)}
			if(any(objects()=="Homol_groups")){rm(Homol_groups)}
			cat(" done.")
		}
	}
	##############################################################################	

	##############################################################################	
	# (6) Filter #################################################################
	cut_delRT_EIC<<-NA
	if(
		(logfile$parameters$filter_profcomp_pos=="TRUE") &
		( (logfile$workflow[names(logfile$workflow)=="isotopologues"] == "yes") || (logfile$workflow[names(logfile$workflow)=="adducts"]=="yes") )
	){
		# (6.1) by ISTD - first get their characteristics on delRT and correl. ###
		fil1<-enviMass::analyseA_links_profiles(
				links_profiles = links_profiles_pos, 
				profileList = profileList_pos, 
				min_rat=.7, 	# isot, adduc: "link counts"/"ref_1"
				min_count=.4, 	# isot, adduc:  "ref_1">=(min_count*number_samples)
				perc =.9,
				for_which=logfile$parameters$for_which_profcomp_pos
			)
		# clean isotopologues ####################################################
		#cut_delRT_isot<-median(fil1$delRT_isot)
		cut_delRT_isot<<-boxplot.stats(c(fil1$delRT_isot))$stats[5]
		cut_cor_isot<<-(boxplot.stats(c(fil1$int_cor_isot))$stats[1])
		if(!is.na(cut_delRT_isot)&!is.na(cut_cor_isot)){
			links_profiles_pos<-enviMass::cleanA_links_profiles(
				links_profiles = links_profiles_pos, 
				profileList = profileList_pos,
				cut_delRT_isot = cut_delRT_isot, 
				cut_cor_isot = cut_cor_isot, 
				cut_frac_iso = .85
			)
		}else{cat("\n No isotopologue linkage filtering feasible")}
		# clean adducts ##########################################################
		#cut_delRT_adduc<-median(fil1$delRT_adduc)
		cut_delRT_adduc<<-boxplot.stats(c(fil1$delRT_adduc))$stats[5]
		if(!is.na(cut_delRT_adduc)){
			links_profiles_pos<-enviMass::cleanB_links_profiles( 
				links_profiles = links_profiles_pos, 
				profileList = profileList_pos,
				cut_delRT_adduc = cut_delRT_adduc, 
				cut_frac_adduc = .85
			)
		}else{cat("\n No adduct linkage filtering feasible")}
		# (6.2) by ISTD - check their EIC correlation ###########################
		fil2<-enviMass::analyseB_links_profiles(
				links_profiles = links_profiles_pos,  
				min_count=.4, 	# isot, adduc:  "ref_1">=(min_count*number_samples)
				for_which=logfile$parameters$for_which_profcomp_pos
			)
		use_EIC<-c(fil2$EIC_cor_isot,fil2$EIC_cor_adduc)
		cut_EIC<<-(boxplot.stats(use_EIC)$stats[1])
		cut_delRT_EIC<<-max(cut_delRT_isot,cut_delRT_adduc)
		if(!is.na(cut_EIC)&!is.na(cut_delRT_EIC)){	
			links_profiles_pos<-enviMass::cleanC_links_profiles(
				links_profiles = links_profiles_pos, 
				profileList = profileList_pos,
				cut_EIC = cut_EIC, 
				cut_frac_EIC = .9, 
				cut_delRT_EIC = cut_delRT_EIC
			)
		}else{cat("\n No EIC linkage filtering feasible")}
		# (6.3) Clean lists ######################################################
		for(n in 1:length(links_profiles_pos)){
			is_empty<-enviMass::analyseC_links_profiles(links_profiles_pos, at_entry = n)
			if(is_empty){
				links_profiles_pos[[n]]<-NA
				profileList_pos[["index_prof"]][as.numeric(names(links_profiles_pos)[n]),"links"]<-0
			}
		}
	}
	##############################################################################		

	##############################################################################		
	# (7) build profile components ###############################################
	plot_it<-FALSE
	plot_what<-"profiles"
	#plot_what<-"relations"	
	with_test<-FALSE
	along<-order(profileList_pos[["index_prof"]][,"number_peaks_total"],decreasing=TRUE) # doesn`t matter actually ...
	if( (logfile$parameters$filter_profcomp_pos=="TRUE") & (!is.na(cut_delRT_EIC)) ){
		use_del_RT<<-cut_delRT_EIC
	}else{
		use_del_RT<<-as.numeric(logfile$parameters$corr_del_RT)
	}
	######################
	if(with_bar){pBar <- txtProgressBar(min = 0, max = length(along), style = 3)}
	for(i in 1:length(along)){
		if(with_bar){setTxtProgressBar(pBar, i, title = NULL, label = NULL)}
		if(logfile$parameters$prof_comp_link_only=="TRUE"){
			# (1) get in-/directly related isotopologue profiles #####################
			prof_isot_IDs<-enviMass::get_isotopol(
				profileList=profileList_pos,
				prof_ID=along[i],
				links_profiles=links_profiles_pos,
				min_peaks=as.numeric(logfile$parameters$corr_min_peaks),
				skip_peaks=as.logical(logfile$parameters$corr_skip_peaks),
				min_cor=as.numeric(logfile$parameters$comp_corr),
				with_test=with_test,
				only_direct=FALSE,
				del_RT=use_del_RT,
				#omit_profiles=in_group
				omit_profiles=FALSE
			)		
			# collect adducts for these isotopologue profiles ########################
			prof_adduct_IDs<-c()
			for(j in 1:length(prof_isot_IDs)){
				got_adducts<-enviMass::get_adducts(
					profileList=profileList_pos,
					prof_ID=prof_isot_IDs[j],
					links_profiles=links_profiles_pos,
					min_peaks=as.numeric(logfile$parameters$corr_min_peaks),
					skip_peaks=as.logical(logfile$parameters$corr_skip_peaks),
					min_cor=as.numeric(logfile$parameters$comp_corr),
					with_test=with_test,
					#omit_profiles=in_group
					omit_profiles=FALSE
				)
				prof_adduct_IDs<-c(prof_adduct_IDs,got_adducts)
			}
			prof_all_IDs<-c(prof_isot_IDs,prof_adduct_IDs)
		}else{
			# (3) get all in-/directly related profiles ##############################	
			prof_all_IDs<-enviMass::get_all(
				profileList=profileList_pos,
				prof_ID=along[i],
				links_profiles=links_profiles_pos,
				min_peaks=as.numeric(logfile$parameters$corr_min_peaks),
				skip_peaks=as.logical(logfile$parameters$corr_skip_peaks),
				min_cor=as.numeric(logfile$parameters$comp_corr),
				with_test=with_test,
				only_direct=FALSE,
				del_RT=use_del_RT,
				#omit_profiles=in_group
				omit_profiles=FALSE
			)		
		}
		##############################################################################
		if(plot_it){
			if(length(prof_all_IDs)==1) next
			enviMass::plot_components(
				profileList=profileList_pos,
				prof_IDs=prof_all_IDs,
				links_profiles=links_profiles_pos,
				what=plot_what,
				xlim=FALSE,ylim=FALSE,await_input=TRUE,
				skipit=TRUE,
				min_peaks=as.numeric(logfile$parameters$corr_min_peaks),
				norma=TRUE
			)	
		}	
		##############################################################################
		if(with_test){
			if(!enviMass::check_interlinked(prof_all_IDs,profileList_pos,links_profiles_pos)) stop("Not all interlinked!")
		}
		##############################################################################
		if(length(prof_all_IDs)>1){
			if(with_test){if(prof_all_IDs[1]!=along[i]){stop("\n\nDebug_not_first!")}}	
			at_entry<-profileList_pos[["index_prof"]][prof_all_IDs[1],"links"]
			links_profiles_pos[[at_entry]][["group"]]<-prof_all_IDs[-1]
		}
		##############################################################################		
	}
	if(with_bar){close(pBar)}
	##############################################################################	

	##############################################################################
	# save! ######################################################################
	save(profileList_pos,file=file.path(as.character(logfile[[1]]),"results","profileList_pos"));
	save(links_profiles_pos,file=file.path(as.character(logfile[[1]]),"results","links_profiles_pos"));	
	##############################################################################	
	rm(links_profiles_pos,profileList_pos)

}


##################################################################################
# NEGATIVE #######################################################################
if( 
	file.exists(file.path(logfile[[1]],"results","profileList_neg")) 
){

	##############################################################################	
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
	peaks<-profileList_neg[["peaks"]][,c("sampleIDs","peakIDs","profileIDs","RT")] # to retrieve relations with, sampleID, peakID, profileID, RT
	ord<-order(peaks[,"sampleIDs"],peaks[,"peakIDs"],peaks[,"profileIDs"],decreasing=FALSE)
	peaks<-peaks[ord,]
	use_entries_profiles<-enviMass::find_empty(links_profiles_neg) # also finds gaps
	profileList_neg[["index_prof"]][,"links"]<-0
	with_bar<-TRUE
	##############################################################################	

	##############################################################################	
	# (1) ANNOTATE TARGET & ISTD SCREENING MACTHES stored in links_peaks_neg #####
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
		cat("\n Annotation of screening results to profiles")
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
						if(	length(links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[2]])>0 ){
							for(k in 1:length(links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[2]])){ # if several compound matches exist for this peak
								if( length(links_profiles_neg[[at_entry]][[2]])==0 ){ # make a new entry for profile link to IS
									links_profiles_neg[[at_entry]][[2]]<-
										data.frame(
											links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[2]][[k]],1,stringsAsFactors = FALSE
										)
									names(links_profiles_neg[[at_entry]][[2]])<-c("Compound","Counts")
								}else{
									at<-which(links_profiles_neg[[at_entry]][[2]][,1]==links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[2]][[k]])
									if(length(at)>0){ 	# increment existing link ...
										links_profiles_neg[[at_entry]][[2]][at,2]<-(links_profiles_neg[[at_entry]][[2]][at,2]+1)
									}else{	# or just add a new one?
										links_profiles_neg[[at_entry]][[2]]<-data.frame(
											c(links_profiles_neg[[at_entry]][[2]][,1], links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[2]][[k]]),
											c(links_profiles_neg[[at_entry]][[2]][,2],1),
											stringsAsFactors = FALSE
										)									
										names(links_profiles_neg[[at_entry]][[2]])<-c("Compound","Counts")
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
											links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[1]][[k]],1,stringsAsFactors = FALSE
										)
									names(links_profiles_neg[[at_entry]][[1]])<-c("Compound","Counts")
								}else{
									at<-which(links_profiles_neg[[at_entry]][[1]][,1]==links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[1]][[k]])
									if(length(at)>0){ 	# increment existing link ...
										links_profiles_neg[[at_entry]][[1]][at,2]<-(links_profiles_neg[[at_entry]][[1]][at,2]+1)
									}else{	# or just add a new one?
										links_profiles_neg[[at_entry]][[1]]<-data.frame(
											c(links_profiles_neg[[at_entry]][[1]][,1], links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[1]][[k]]),
											c(links_profiles_neg[[at_entry]][[1]][,2],1),
											stringsAsFactors = FALSE
										)									
										names(links_profiles_neg[[at_entry]][[1]])<-c("Compound","Counts")
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
		cat(" done.")
	}
	##############################################################################	
	
	##############################################################################	
	# for each profile, filter all non-target relations to other profiles ########
	# (2) INSERT EIC and CO-OCCURRENCE INFORMATION ###############################
	if(logfile$workflow[names(logfile$workflow)=="EIC_correlation"]=="yes"){
		cat("\n Retrieving EIC correlation links ")
		# (2.1) INSERT EIC LINKS #################################################
		##########################################################################
		forIDs<-profileList_neg[["sampleID"]]
		for_files<-list.files(file.path(logfile[[1]],"results","componentization","EIC_corr"))
		keep<-match(forIDs,for_files) # which files are available?
		if(any(is.na(keep))){cat("\n Just note: not all files found in profiles have EIC correlation results (1). \n")}
		TRUE_IDs<-(measurements$ID[measurements$adducts=="TRUE"]) # for files which have run through that step
		keep2<-match(forIDs,for_files) # ensure file availability
		forIDs<-forIDs[!is.na(keep) & !is.na(keep2)]
		if(logfile$parameters$dofile_latest_profcomp=="TRUE"){ # restrict to latest files
			atPOSIX<-profileList_neg[["datetime"]];
			matchID<-profileList_neg[["sampleID"]];
			atdate<-c();attime<-c();
			for(i in 1:length(atPOSIX)){
				atdate<-c(atdate, strsplit(atPOSIX[i]," ")[[1]][1]);
				attime<-c(attime, strsplit(atPOSIX[i]," ")[[1]][2]);
			}
			attime<-as.difftime(attime);
			atdate<-as.Date(atdate, tz="GMT");
			ord<-order(as.numeric(atdate),as.numeric(attime),matchID,decreasing=TRUE);
			matchID<-matchID[ord];
			forIDs<-forIDs[match(forIDs,matchID)]
			if(length(forIDs)>as.numeric(logfile$parameters$numfile_latest_profcomp)){
				forIDs<-forIDs[1:as.numeric(logfile$parameters$numfile_latest_profcomp)]
			}
		}
		not_found1<-0;inserted1<-0
		if(length(forIDs)>0){
			if(with_bar){pBar <- txtProgressBar(min = 0, max = length(forIDs), style = 3)}
			for(i in 1:length(forIDs)){
				if(with_bar){setTxtProgressBar(pBar, i, title = NULL, label = NULL)}
				load(file=file.path(logfile[[1]],"results","componentization","EIC_corr",forIDs[i]))
				#
				#EIC_pairs<-EIC_pairs[EIC_pairs[,4]>=logfile$parameters$EICor_mincor,,drop=FALSE]		
				#
				if(length(EIC_pairs[,1])==0){next}
				# find profiles for first peak
				get1<-cbind(
					rep(as.numeric(forIDs[i]),length(EIC_pairs[,1])),EIC_pairs[,1] # first column sorted correctly
				)
				found1<-enviMass::rows_compare(get1,peaks[,c("sampleIDs","peakIDs")],row_order=FALSE,column_order_a=FALSE,column_order_b=FALSE,get_index=TRUE)
				# find profiles for second peak
				get2<-cbind(
					rep(as.numeric(forIDs[i]),length(EIC_pairs[,2])),EIC_pairs[,2] # second column requires sorting
				)
				found2<-enviMass::rows_compare(get2,peaks[,c("sampleIDs","peakIDs")],row_order=FALSE,column_order_a=TRUE,column_order_b=FALSE,get_index=TRUE)
				for(j in 1:length(found1)){ # insert links
					# insert PROFILE LINKS ########################################
					if(found1[j]==0){not_found1<-(not_found1+1);next} # e.g., peak blind-removed 
					if(found2[j]==0){not_found1<-(not_found1+1);next}		
					inserted1<-(inserted1+1);				
					# (1) insert link to second profile for the first profile
					prof1<-peaks[found1[j],"profileIDs"][[1]]
					prof2<-peaks[found2[j],"profileIDs"][[1]]
					# enable check ... pairs musst have very similar retention times
					if(FALSE){
						cat("\n");
						cat(peaks[found1[j],"RT"]);cat(" - ")
						cat(peaks[found2[j],"RT"])
					}
					if(profileList_neg[["index_prof"]][prof1,"profile_ID"]!=prof1){stop("\nComponentization: debug me, #1!")}				
					if(profileList_neg[["index_prof"]][prof1,"links"]==0){ 	# establish a new link ...
						if(length(use_entries_profiles)>0){
							at_entry_1<-use_entries_profiles[1]
							use_entries<-use_entries_profiles[-1]
						}else{
							at_entry_1<-(length(links_profiles_neg)+1)
						}
						links_profiles_neg[[at_entry_1]]<-enviMass::new_entry_links_profiles(profileList_neg[["index_prof"]][prof1,"number_peaks_total"][[1]])
						names(links_profiles_neg)[at_entry_1]<-as.character(prof1)
						profileList_neg[["index_prof"]][prof1,"links"]<-at_entry_1						
					}else{
						at_entry_1<-profileList_neg[["index_prof"]][prof1,"links"]
					}
					here1<-which(links_profiles_neg[[at_entry_1]][["EIC"]][,"linked profile"]==prof2)
					if(length(here1)==0){			
						links_profiles_neg[[at_entry_1]][["EIC"]]<-rbind(
							links_profiles_neg[[at_entry_1]][["EIC"]], c(prof2,1,0,0,0,0,1,NA)
						)
						here1<-dim(links_profiles_neg[[at_entry_1]][["EIC"]])[1]
					}else{
						links_profiles_neg[[at_entry_1]][["EIC"]][here1,"link counts"]<-(links_profiles_neg[[at_entry_1]][["EIC"]][here1,"link counts"]+1)
					}
					# (2) insert link to first profile for the second profile
					if(profileList_neg[["index_prof"]][prof2,"profile_ID"]!=prof2){stop("\nCross-profile componentization: debug me, EIC_1!")}				
					if(profileList_neg[["index_prof"]][prof2,"links"]==0){ 	# establish a new link ...
						if(length(use_entries_profiles)>0){
							at_entry_2<-use_entries_profiles[1]
							use_entries<-use_entries_profiles[-1]
						}else{
							at_entry_2<-(length(links_profiles_neg)+1)
						}
						links_profiles_neg[[at_entry_2]]<-enviMass::new_entry_links_profiles(profileList_neg[["index_prof"]][prof2,"number_peaks_total"][[1]])		
						names(links_profiles_neg)[at_entry_2]<-as.character(prof2)
						profileList_neg[["index_prof"]][prof2,"links"]<-at_entry_2						
					}else{
						at_entry_2<-profileList_neg[["index_prof"]][prof2,"links"]
					}
					here2<-which(links_profiles_neg[[at_entry_2]][["EIC"]][,"linked profile"]==prof1)
					if(length(here2)==0){
						links_profiles_neg[[at_entry_2]][["EIC"]]<-rbind(
							links_profiles_neg[[at_entry_2]][["EIC"]], c(prof1,1,0,0,0,0,0,NA)
						)
						here2<-dim(links_profiles_neg[[at_entry_2]][["EIC"]])[1]
					}else{
						links_profiles_neg[[at_entry_2]][["EIC"]][here2,"link counts"]<-(links_profiles_neg[[at_entry_2]][["EIC"]][here2,"link counts"]+1)
					}				
					# insert EIC correlation
					if(links_profiles_neg[[at_entry_1]][["EIC"]][here1,"use"]==1){
						if(length(links_profiles_neg[[at_entry_1]]$EIC_cor)<here1){
							links_profiles_neg[[at_entry_1]]$EIC_cor[[here1]]<-EIC_pairs[j,4]
						}else{
							if(is.null(links_profiles_neg[[at_entry_1]]$EIC_cor[[here1]])){
								links_profiles_neg[[at_entry_1]]$EIC_cor[[here1]]<-EIC_pairs[j,4]
							}else{
								links_profiles_neg[[at_entry_1]]$EIC_cor[[here1]]<-c(links_profiles_neg[[at_entry_1]]$EIC_cor[[here1]],EIC_pairs[j,4])
							}
						}
					}else{
						if(length(links_profiles_neg[[at_entry_2]]$EIC_cor)<here2){
							links_profiles_neg[[at_entry_2]]$EIC_cor[[here2]]<-EIC_pairs[j,4]
						}else{
							if(is.null(links_profiles_neg[[at_entry_2]]$EIC_cor[[here2]])){
								links_profiles_neg[[at_entry_2]]$EIC_cor[[here2]]<-EIC_pairs[j,4]
							}else{
								links_profiles_neg[[at_entry_2]]$EIC_cor[[here2]]<-c(links_profiles_neg[[at_entry_2]]$EIC_cor[[here2]],EIC_pairs[j,4])
							}
						}
					}
					# insert PEAK LINKS ##########################################
				}
			}
			if(with_bar){close(pBar)}		
			if(any(objects()=="EIC_pairs")){rm(EIC_pairs)}	
			cat(" done.")
		}
	}
	# (3) INSERT ISOTOPOLOGUE LINKS ##############################################
	if(logfile$workflow[names(logfile$workflow)=="isotopologues"]=="yes"){	
		forIDs<-profileList_neg[["sampleID"]]
		for_files<-list.files(file.path(logfile[[1]],"results","componentization","isotopologues"))
		keep<-match(forIDs,for_files) # which files are available?
		if(any(is.na(keep))){cat("\n Just note: not all files found in profiles have isotopologue links available. \n")}
		TRUE_IDs<-(measurements$ID[measurements$adducts=="TRUE"]) # for files which have run through that step
		keep2<-match(forIDs,for_files)
		forIDs<-forIDs[!is.na(keep) & !is.na(keep2)]
		if(logfile$parameters$dofile_latest_profcomp=="TRUE"){ # restrict to latest files
			atPOSIX<-profileList_neg[["datetime"]];
			matchID<-profileList_neg[["sampleID"]];
			atdate<-c();attime<-c();
			for(i in 1:length(atPOSIX)){
				atdate<-c(atdate, strsplit(atPOSIX[i]," ")[[1]][1]);
				attime<-c(attime, strsplit(atPOSIX[i]," ")[[1]][2]);
			}
			attime<-as.difftime(attime);
			atdate<-as.Date(atdate, tz="GMT");
			ord<-order(as.numeric(atdate),as.numeric(attime),matchID,decreasing=TRUE);
			matchID<-matchID[ord];
			forIDs<-forIDs[match(forIDs,matchID)]
			if(length(forIDs)>as.numeric(logfile$parameters$numfile_latest_profcomp)){
				forIDs<-forIDs[1:as.numeric(logfile$parameters$numfile_latest_profcomp)]
			}
		}
		not_found3<-0;inserted3<-0
		if(length(forIDs)>0){
			cat("\n Retrieving isotopologue links ")
			if(with_bar){pBar <- txtProgressBar(min = 0, max = length(forIDs), style = 3)}
			for(i in 1:length(forIDs)){
				if(with_bar){setTxtProgressBar(pBar, i, title = NULL, label = NULL)}
				load(file=file.path(logfile[[1]],"results","componentization","isotopologues",forIDs[i]))
				if(length(Isot_pairs[,1])==0){next}
				# find profiles for first peak
				get1<-cbind(
					rep(as.numeric(forIDs[i]),length(Isot_pairs[,1])),Isot_pairs[,1]
				)
				found1<-enviMass::rows_compare(get1,peaks[,c("sampleIDs","peakIDs")],row_order=FALSE,column_order_a=FALSE,column_order_b=FALSE,get_index=TRUE)
				# find profiles for second peak
				get2<-cbind(
					rep(as.numeric(forIDs[i]),length(Isot_pairs[,2])),Isot_pairs[,2]
				)
				found2<-enviMass::rows_compare(get2,peaks[,c("sampleIDs","peakIDs")],row_order=FALSE,column_order_a=TRUE,column_order_b=FALSE,get_index=TRUE)
				for(j in 1:length(found1)){ # insert links
					# insert PROFILE LINKS ########################################
					if(found1[j]==0){not_found3<-(not_found3+1);next} # e.g., peak blind-removed 
					if(found2[j]==0){not_found3<-(not_found3+1);next}			
					inserted3<-(inserted3+1);				
					# enable check ... pairs musst have very similar retention times
					if(FALSE){
						cat("\n");
						cat(peaks[found1[j],"RT"]);cat(" - ")
						cat(peaks[found2[j],"RT"])				
					}
					# (1) insert link to second profile for the first profile
					prof1<-peaks[found1[j],"profileIDs"][[1]]
					prof2<-peaks[found2[j],"profileIDs"][[1]]				
					if(profileList_neg[["index_prof"]][prof1,"profile_ID"]!=prof1){stop("\nComponentization: debug me, #1!")}				
					if(profileList_neg[["index_prof"]][prof1,"links"]==0){ 	# establish a new link ...
						if(length(use_entries_profiles)>0){
							at_entry_1<-use_entries_profiles[1]
							use_entries<-use_entries_profiles[-1]
						}else{
							at_entry_1<-(length(links_profiles_neg)+1)
						}
						links_profiles_neg[[at_entry_1]]<-enviMass::new_entry_links_profiles(profileList_neg[["index_prof"]][prof1,"number_peaks_total"][[1]])
						names(links_profiles_neg)[at_entry_1]<-as.character(prof1)
						profileList_neg[["index_prof"]][prof1,"links"]<-at_entry_1						
					}else{
						at_entry_1<-profileList_neg[["index_prof"]][prof1,"links"]
					}
					here1<-which(links_profiles_neg[[at_entry_1]][["isot"]][,"linked profile"]==prof2)
					if(length(here1)==0){
						links_profiles_neg[[at_entry_1]][["isot"]]<-rbind(
							links_profiles_neg[[at_entry_1]][["isot"]], c(prof2,1,0,1,NA)
						)
						here1<-dim(links_profiles_neg[[at_entry_1]][["isot"]])[1]
					}else{
						links_profiles_neg[[at_entry_1]][["isot"]][here1,"link counts"]<-(links_profiles_neg[[at_entry_1]][["isot"]][here1,"link counts"]+1)
					}
					# (2) insert link to first profile for the second profile
					if(profileList_neg[["index_prof"]][prof2,"profile_ID"]!=prof2){stop("\nComponentization: debug me, #1!")}				
					if(profileList_neg[["index_prof"]][prof2,"links"]==0){ 	# establish a new link ...
						if(length(use_entries_profiles)>0){
							at_entry_2<-use_entries_profiles[1]
							use_entries<-use_entries_profiles[-1]
						}else{																		
							at_entry_2<-(length(links_profiles_neg)+1)
						}
						links_profiles_neg[[at_entry_2]]<-enviMass::new_entry_links_profiles(profileList_neg[["index_prof"]][prof2,"number_peaks_total"][[1]])
						names(links_profiles_neg)[at_entry_2]<-as.character(prof2)
						profileList_neg[["index_prof"]][prof2,"links"]<-at_entry_2						
					}else{
						at_entry_2<-profileList_neg[["index_prof"]][prof2,"links"]
					}
					here2<-which(links_profiles_neg[[at_entry_2]][["isot"]][,"linked profile"]==prof1)
					if(length(here2)==0){
						links_profiles_neg[[at_entry_2]][["isot"]]<-rbind(
							links_profiles_neg[[at_entry_2]][["isot"]], c(prof1,1,0,0,NA)
						)
						here2<-dim(links_profiles_neg[[at_entry_2]][["isot"]])[1]
					}else{
						links_profiles_neg[[at_entry_2]][["isot"]][here2,"link counts"]<-(links_profiles_neg[[at_entry_2]][["isot"]][here2,"link counts"]+1)
					}				
					# INSERT ref_1: total number of co-occurences ################
					got_entry<-FALSE
					#if(dim(links_profiles_neg[[at_entry_1]]$EIC)[1]!=0){ # already dealed with during EIC correlation?
					#	here3<-which(links_profiles_neg[[at_entry_1]]$EIC[,"linked profile"]==prof2)
					#	if(length(here3)>0){
					#		got_entry<-TRUE
					#		if(length(here3)>1){
					#			stop("\n DEBUG ME !")
					#		}
					#	}
					#}
					if(got_entry){ # use value of existing entry ...
						links_profiles_neg[[at_entry_1]]$isot[here1,"ref_1"]<-links_profiles_neg[[at_entry_1]]$EIC[here3,"ref_1"]
						links_profiles_neg[[at_entry_2]]$isot[here2,"ref_1"]<-links_profiles_neg[[at_entry_1]]$EIC[here3,"ref_1"]
					}else{ # ... or retrieve them
						these<-profileList_neg[["peaks"]][
							profileList_neg[["index_prof"]][prof1,"start_ID"]:profileList_neg[["index_prof"]][prof1,"end_ID"]
						,"sampleIDs"]
						those<-profileList_neg[["peaks"]][
							profileList_neg[["index_prof"]][prof2,"start_ID"]:profileList_neg[["index_prof"]][prof2,"end_ID"]
						,"sampleIDs"]
						matched<-match(these,those)
						not_NA<-sum(!is.na(matched))
						links_profiles_neg[[at_entry_1]]$isot[here1,"ref_1"]<-not_NA
						links_profiles_neg[[at_entry_2]]$isot[here2,"ref_1"]<-not_NA
					}
					if(any(links_profiles_neg[[at_entry_1]]$isot[,"ref_1"]==0)){stop("\n DEBUG ME ! FOUND_1")}
					if(any(links_profiles_neg[[at_entry_2]]$isot[,"ref_1"]==0)){stop("\n DEBUG ME ! FOUND_2")}
					# insert PEAK LINKS ##########################################
				}
			}
			if(with_bar){close(pBar)}
			if(any(objects()=="Isot_pairs")){rm(Isot_pairs)}
			cat(" done.")
		}
	}
	# (4) INSERT ADDUCT LINKS ####################################################
	if(logfile$workflow[names(logfile$workflow)=="adducts"]=="yes"){
		forIDs<-profileList_neg[["sampleID"]]
		for_files<-list.files(file.path(logfile[[1]],"results","componentization","adducts"))
		keep<-match(forIDs,for_files) # which files are available?
		if(any(is.na(keep))){cat("\n Just note: not all files found in profiles have adduct links available. \n")}
		TRUE_IDs<-(measurements$ID[measurements$adducts=="TRUE"]) # for files which have run through that step
		keep2<-match(forIDs,for_files)
		forIDs<-forIDs[!is.na(keep) & !is.na(keep2)]
		if(logfile$parameters$dofile_latest_profcomp=="TRUE"){ # restrict to latest files
			atPOSIX<-profileList_neg[["datetime"]];
			matchID<-profileList_neg[["sampleID"]];
			atdate<-c();attime<-c();
			for(i in 1:length(atPOSIX)){
				atdate<-c(atdate, strsplit(atPOSIX[i]," ")[[1]][1]);
				attime<-c(attime, strsplit(atPOSIX[i]," ")[[1]][2]);
			}
			attime<-as.difftime(attime);
			atdate<-as.Date(atdate, tz="GMT");
			ord<-order(as.numeric(atdate),as.numeric(attime),matchID,decreasing=TRUE);
			matchID<-matchID[ord];
			forIDs<-forIDs[match(forIDs,matchID)]
			if(length(forIDs)>as.numeric(logfile$parameters$numfile_latest_profcomp)){
				forIDs<-forIDs[1:as.numeric(logfile$parameters$numfile_latest_profcomp)]
			}
		}
		not_found4<-0;inserted4<-0
		if(length(forIDs)>0){
			cat("\n Retrieving adduct links ")
			if(with_bar){pBar <- txtProgressBar(min = 0, max = length(forIDs), style = 3)}
			for(i in 1:length(forIDs)){
				if(with_bar){setTxtProgressBar(pBar, i, title = NULL, label = NULL)}
				load(file=file.path(logfile[[1]],"results","componentization","adducts",forIDs[i]))
				if(length(Adduct_pairs[,1])==0){next}
				# find profiles for first peak
				get1<-cbind(
					rep(as.numeric(forIDs[i]),length(Adduct_pairs[,1])),Adduct_pairs[,1]
				)
				found1<-enviMass::rows_compare(get1,peaks[,c("sampleIDs","peakIDs")],row_order=FALSE,column_order_a=FALSE,column_order_b=FALSE,get_index=TRUE)
				# find profiles for second peak
				get2<-cbind(
					rep(as.numeric(forIDs[i]),length(Adduct_pairs[,2])),Adduct_pairs[,2]
				)
				found2<-enviMass::rows_compare(get2,peaks[,c("sampleIDs","peakIDs")],row_order=FALSE,column_order_a=TRUE,column_order_b=FALSE,get_index=TRUE)
				for(j in 1:length(found1)){ # insert links	
					if(found1[j]==0){not_found4<-(not_found4+1);next} # e.g., peak blind-removed 
					if(found2[j]==0){not_found4<-(not_found4+1);next}		
					inserted4<-(inserted4+1);
					# enable check ... pairs musst have very similar retention times
					if(FALSE){
						cat("\n");
						cat(peaks[found1[j],"RT"]);cat(" - ")
						cat(peaks[found2[j],"RT"])				
					}
					# insert PROFILE LINKS ######################################	
					# (1) insert link to second profile for the first profile
					prof1<-peaks[found1[j],"profileIDs"][[1]]
					prof2<-peaks[found2[j],"profileIDs"][[1]]
					if(profileList_neg[["index_prof"]][prof1,"profile_ID"]!=prof1){stop("\nComponentization: debug me, #1!")}				
					if(profileList_neg[["index_prof"]][prof1,"links"]==0){ 	# establish a new link ...
						if(length(use_entries_profiles)>0){
							at_entry_1<-use_entries_profiles[1]
							use_entries<-use_entries_profiles[-1]
						}else{
							at_entry_1<-(length(links_profiles_neg)+1)
						}
						links_profiles_neg[[at_entry_1]]<-enviMass::new_entry_links_profiles(profileList_neg[["index_prof"]][prof1,"number_peaks_total"][[1]])
						names(links_profiles_neg)[at_entry_1]<-as.character(prof1)
						profileList_neg[["index_prof"]][prof1,"links"]<-at_entry_1						
					}else{
						at_entry_1<-profileList_neg[["index_prof"]][prof1,"links"]
					}
					here1<-which(links_profiles_neg[[at_entry_1]][["adduc"]][,"linked profile"]==prof2)
					if(length(here1)==0){
						links_profiles_neg[[at_entry_1]][["adduc"]]<-rbind(
							links_profiles_neg[[at_entry_1]][["adduc"]], c(prof2,1,0,1,NA)
						)
						here1<-dim(links_profiles_neg[[at_entry_1]][["adduc"]])[1]
					}else{
						links_profiles_neg[[at_entry_1]][["adduc"]][here1,"link counts"]<-(links_profiles_neg[[at_entry_1]][["adduc"]][here1,"link counts"]+1)
					}
					# (2) insert link to first profile for the second profile
					if(profileList_neg[["index_prof"]][prof2,"profile_ID"]!=prof2){stop("\nComponentization: debug me, #1!")}				
					if(profileList_neg[["index_prof"]][prof2,"links"]==0){ 	# establish a new link ...
						if(length(use_entries_profiles)>0){
							at_entry_2<-use_entries_profiles[1]
							use_entries<-use_entries_profiles[-1]
						}else{
							at_entry_2<-(length(links_profiles_neg)+1)
						}
						links_profiles_neg[[at_entry_2]]<-enviMass::new_entry_links_profiles(profileList_neg[["index_prof"]][prof2,"number_peaks_total"][[1]])
						names(links_profiles_neg)[at_entry_2]<-as.character(prof2)
						profileList_neg[["index_prof"]][prof2,"links"]<-at_entry_2						
					}else{
						at_entry_2<-profileList_neg[["index_prof"]][prof2,"links"]
					}
					here2<-which(links_profiles_neg[[at_entry_2]][["adduc"]][,"linked profile"]==prof1)
					if(length(here2)==0){
						links_profiles_neg[[at_entry_2]][["adduc"]]<-rbind(
							links_profiles_neg[[at_entry_2]][["adduc"]], c(prof1,1,0,0,NA)
						)
						here2<-dim(links_profiles_neg[[at_entry_2]][["adduc"]])[1]
					}else{
						links_profiles_neg[[at_entry_2]][["adduc"]][here2,"link counts"]<-(links_profiles_neg[[at_entry_2]][["adduc"]][here2,"link counts"]+1)
					}		
					# INSERT ref_1: total number of co-occurences ################
					got_entry<-FALSE
					#if(dim(links_profiles_neg[[at_entry_1]]$EIC)[1]!=0){ # already dealed with during EIC correlation?
					#	here3<-which(links_profiles_neg[[at_entry_1]]$EIC[,"linked profile"]==prof2)
					#	if(length(here3)>0){
					#		got_entry<-TRUE
					#	}
					#}
					if(got_entry){ # use value of existing entry ...
						links_profiles_neg[[at_entry_1]]$adduc[here1,"ref_1"]<-links_profiles_neg[[at_entry_1]]$EIC[here3,"ref_1"]
						links_profiles_neg[[at_entry_2]]$adduc[here2,"ref_1"]<-links_profiles_neg[[at_entry_1]]$EIC[here3,"ref_1"]
					}else{ # ... or retrieve them
						these<-profileList_neg[["peaks"]][
							profileList_neg[["index_prof"]][prof1,"start_ID"]:profileList_neg[["index_prof"]][prof1,"end_ID"]
						,"sampleIDs"]
						those<-profileList_neg[["peaks"]][
							profileList_neg[["index_prof"]][prof2,"start_ID"]:profileList_neg[["index_prof"]][prof2,"end_ID"]
						,"sampleIDs"]
						matched<-match(these,those)
						not_NA<-sum(!is.na(matched))
						links_profiles_neg[[at_entry_1]]$adduc[here1,"ref_1"]<-not_NA
						links_profiles_neg[[at_entry_2]]$adduc[here2,"ref_1"]<-not_NA
					}
					if(any(links_profiles_neg[[at_entry_1]]$adduc[,"ref_1"]==0)){stop("\n DEBUG ME ! FOUND_3")}
					if(any(links_profiles_neg[[at_entry_2]]$adduc[,"ref_1"]==0)){stop("\n DEBUG ME ! FOUND_4")}
					# insert PEAK LINKS? ########################################
				}
			}
			if(with_bar){close(pBar)}
			if(any(objects()=="Adduct_pairs")){rm(Adduct_pairs)}
			cat(" done.")
		}	
	}
	# (5) INSERT HOMOLOGUE SERIES LINKS ##########################################
	if(logfile$workflow[names(logfile$workflow)=="homologues"]=="yes"){
		forIDs<-profileList_neg[["sampleID"]]
		for_files<-list.files(file.path(logfile[[1]],"results","componentization","homologues"))
		keep<-match(forIDs,for_files) # which files are available?
		if(any(is.na(keep))){cat("\n Just note: not all files found in profiles have homologue links available. \n")}
		TRUE_IDs<-(measurements$ID[measurements$homologues=="TRUE"]) # for files which have run through that step
		keep2<-match(forIDs,for_files)
		forIDs<-forIDs[!is.na(keep) & !is.na(keep2)]
		if(logfile$parameters$dofile_latest_profcomp=="TRUE"){ # restrict to latest files
			atPOSIX<-profileList_neg[["datetime"]];
			matchID<-profileList_neg[["sampleID"]];
			atdate<-c();attime<-c();
			for(i in 1:length(atPOSIX)){
				atdate<-c(atdate, strsplit(atPOSIX[i]," ")[[1]][1]);
				attime<-c(attime, strsplit(atPOSIX[i]," ")[[1]][2]);
			}
			attime<-as.difftime(attime);
			atdate<-as.Date(atdate, tz="GMT");
			ord<-order(as.numeric(atdate),as.numeric(attime),matchID,decreasing=TRUE);
			matchID<-matchID[ord];
			forIDs<-forIDs[match(forIDs,matchID)]
			if(length(forIDs)>as.numeric(logfile$parameters$numfile_latest_profcomp)){
				forIDs<-forIDs[1:as.numeric(logfile$parameters$numfile_latest_profcomp)]
			}
		}
		not_found5<-0;inserted5<-0
		if(length(forIDs)>0){
			cat("\n Retrieving homologue links ")
			if(with_bar){pBar <- txtProgressBar(min = 0, max = length(forIDs), style = 3)}
			for(i in 1:length(forIDs)){
				if(with_bar){setTxtProgressBar(pBar, i, title = NULL, label = NULL)}
				load(file=file.path(logfile[[1]],"results","componentization","homologues",forIDs[i]))	
				if(length(Homol_groups[,1])==0){next}	
				keep_group<-rep(TRUE,max(Homol_groups[,3])) # make sure all peaks in a group are still present
				# find profiles for first peak
				get1<-cbind(
					rep(as.numeric(forIDs[i]),length(Homol_groups[,1])),Homol_groups[,1]
				)
				found1<-enviMass::rows_compare(get1,peaks[,c("sampleIDs","peakIDs")],row_order=FALSE,column_order_a=FALSE,column_order_b=FALSE,get_index=TRUE)
				keep_group[Homol_groups[found1==0,3]]<-FALSE
				if(!any(keep_group)){next} # none left with all peaks still embedded into profiles
				# find profiles for second peak
				get2<-cbind(
					rep(as.numeric(forIDs[i]),length(Homol_groups[,2])),Homol_groups[,2]
				)
				found2<-enviMass::rows_compare(get2,peaks[,c("sampleIDs","peakIDs")],row_order=FALSE,column_order_a=TRUE,column_order_b=FALSE,get_index=TRUE)
				if(!any(keep_group)){next} # none left with all peaks still embedded into profiles
				for(j in 1:length(found1)){ # insert links	
					if(found1[j]==0){not_found5<-(not_found5+1);next} # e.g., peak blind-removed 
					if(found2[j]==0){not_found5<-(not_found5+1);next}
					if(Homol_groups[j,3]!=Homol_groups[j,3]){stop("\n Debug componentization for its homologue relations!")}
					if(!keep_group[Homol_groups[j,3]]){next} # not all peaks for this group were embedded into profiles, e.g., peak blind-removed
					inserted5<-(inserted5+1);
					# enable check ... pairs musst have very similar retention times
					if(FALSE){
						cat("\n");
						cat(peaks[found1[j],"RT"]);cat(" - ")
						cat(peaks[found2[j],"RT"])				
					}
					# insert PROFILE LINKS #############################################	
					# (1) insert link to second profile for the first profile
					prof1<-peaks[found1[j],"profileIDs"][[1]]
					prof2<-peaks[found2[j],"profileIDs"][[1]]
					if(profileList_neg[["index_prof"]][prof1,"profile_ID"]!=prof1){stop("\nComponentization: debug me, #1!")}				
					if(profileList_neg[["index_prof"]][prof1,"links"]==0){ 	# establish a new link ...
						if(length(use_entries_profiles)>0){
							at_entry<-use_entries_profiles[1]
							use_entries<-use_entries_profiles[-1]
						}else{
							at_entry<-(length(links_profiles_neg)+1)
						}
						links_profiles_neg[[at_entry]]<-enviMass::new_entry_links_profiles(profileList_neg[["index_prof"]][prof1,"number_peaks_total"][[1]])
						names(links_profiles_neg)[at_entry]<-as.character(prof1)
						profileList_neg[["index_prof"]][prof1,"links"]<-at_entry						
					}else{
						at_entry<-profileList_neg[["index_prof"]][prof1,"links"]
					}
					here<-which(links_profiles_neg[[at_entry]][["homol"]][,"linked profile"]==prof2)
					if(length(here)==0){
						links_profiles_neg[[at_entry]][["homol"]]<-rbind(
							links_profiles_neg[[at_entry]][["homol"]], c(prof2,1,0)
						)
					}else{
						links_profiles_neg[[at_entry]][["homol"]][here,"link counts"]<-(links_profiles_neg[[at_entry]][["homol"]][here,"link counts"]+1)
					}
					# (2) insert link to first profile for the second profile
					if(profileList_neg[["index_prof"]][prof2,"profile_ID"]!=prof2){stop("\nComponentization: debug me, #1!")}				
					if(profileList_neg[["index_prof"]][prof2,"links"]==0){ 	# establish a new link ...
						if(length(use_entries_profiles)>0){
							at_entry<-use_entries_profiles[1]
							use_entries<-use_entries_profiles[-1]
						}else{
							at_entry<-(length(links_profiles_neg)+1)
						}
						links_profiles_neg[[at_entry]]<-enviMass::new_entry_links_profiles(profileList_neg[["index_prof"]][prof2,"number_peaks_total"][[1]])
						names(links_profiles_neg)[at_entry]<-as.character(prof2)
						profileList_neg[["index_prof"]][prof2,"links"]<-at_entry						
					}else{
						at_entry<-profileList_neg[["index_prof"]][prof2,"links"]
					}
					here<-which(links_profiles_neg[[at_entry]][["homol"]][,"linked profile"]==prof1)
					if(length(here)==0){
						links_profiles_neg[[at_entry]][["homol"]]<-rbind(
							links_profiles_neg[[at_entry]][["homol"]], c(prof1,1,0)
						)
					}else{
						links_profiles_neg[[at_entry]][["homol"]][here,"link counts"]<-(links_profiles_neg[[at_entry]][["homol"]][here,"link counts"]+1)
					}				
				}	
			}
			if(with_bar){close(pBar)}
			if(any(objects()=="Homol_groups")){rm(Homol_groups)}
			cat(" done.")
		}
	}
	##############################################################################	

	##############################################################################	
	# (6) Filter #################################################################
	cut_delRT_EIC<<-NA
	if(
		(logfile$parameters$filter_profcomp_neg=="TRUE") &
		( (logfile$workflow[names(logfile$workflow)=="isotopologues"] == "yes") || (logfile$workflow[names(logfile$workflow)=="adducts"]=="yes") )
	){
		# (6.1) by ISTD - first get their characteristics on delRT and correl. ###
		fil1<-enviMass::analyseA_links_profiles(
				links_profiles = links_profiles_neg, 
				profileList = profileList_neg, 
				min_rat=.7, 	# isot, adduc: "link counts"/"ref_1"
				min_count=.4, 	# isot, adduc:  "ref_1">=(min_count*number_samples)
				perc =.9,
				for_which=logfile$parameters$for_which_profcomp_neg
			)
		# clean isotopologues ####################################################
		#cut_delRT_isot<-median(fil1$delRT_isot)
		cut_delRT_isot<<-boxplot.stats(c(fil1$delRT_isot))$stats[5]
		cut_cor_isot<<-(boxplot.stats(c(fil1$int_cor_isot))$stats[1])
		if(!is.na(cut_delRT_isot)&!is.na(cut_cor_isot)){
			links_profiles_neg<-enviMass::cleanA_links_profiles(
				links_profiles = links_profiles_neg, 
				profileList = profileList_neg,
				cut_delRT_isot = cut_delRT_isot, 
				cut_cor_isot = cut_cor_isot, 
				cut_frac_iso = .85
			)
		}else{cat("\n No isotopologue linkage filtering feasible")}
		# clean adducts ##########################################################
		#cut_delRT_adduc<-median(fil1$delRT_adduc)
		cut_delRT_adduc<<-boxplot.stats(c(fil1$delRT_adduc))$stats[5]
		if(!is.na(cut_delRT_adduc)){
			links_profiles_neg<-enviMass::cleanB_links_profiles( 
				links_profiles = links_profiles_neg, 
				profileList = profileList_neg,
				cut_delRT_adduc = cut_delRT_adduc, 
				cut_frac_adduc = .85
			)
		}else{cat("\n No adduct linkage filtering feasible")}
		# (6.2) by ISTD - check their EIC correlation ###########################
		fil2<-enviMass::analyseB_links_profiles(
				links_profiles = links_profiles_neg,  
				min_count=.4, 	# isot, adduc:  "ref_1">=(min_count*number_samples)
				for_which=logfile$parameters$for_which_profcomp_neg
			)
		use_EIC<-c(fil2$EIC_cor_isot,fil2$EIC_cor_adduc)
		cut_EIC<<-(boxplot.stats(use_EIC)$stats[1])
		cut_delRT_EIC<<-max(cut_delRT_isot,cut_delRT_adduc)
		if(!is.na(cut_EIC)&!is.na(cut_delRT_EIC)){	
			links_profiles_neg<-enviMass::cleanC_links_profiles(
				links_profiles = links_profiles_neg, 
				profileList = profileList_neg,
				cut_EIC = cut_EIC, 
				cut_frac_EIC = .9, 
				cut_delRT_EIC = cut_delRT_EIC
			)
		}else{cat("\n No EIC linkage filtering feasible")}
		# (6.3) Clean lists ######################################################
		for(n in 1:length(links_profiles_neg)){
			is_empty<-enviMass::analyseC_links_profiles(links_profiles_neg, at_entry = n)
			if(is_empty){
				links_profiles_neg[[n]]<-NA
				profileList_neg[["index_prof"]][as.numeric(names(links_profiles_neg)[n]),"links"]<-0
			}
		}
	}
	##############################################################################		

	##############################################################################		
	# (7) build profile components ###############################################
	plot_it<-FALSE
	plot_what<-"profiles"
	#plot_what<-"relations"	
	with_test<-FALSE
	along<-order(profileList_neg[["index_prof"]][,"number_peaks_total"],decreasing=TRUE) # doesn`t matter actually ...
	if( (logfile$parameters$filter_profcomp_pos=="TRUE") & (!is.na(cut_delRT_EIC)) ){
		use_del_RT<<-cut_delRT_EIC
	}else{
		use_del_RT<<-as.numeric(logfile$parameters$corr_del_RT)
	}
	######################
	if(with_bar){pBar <- txtProgressBar(min = 0, max = length(along), style = 3)}
	for(i in 1:length(along)){
		if(with_bar){setTxtProgressBar(pBar, i, title = NULL, label = NULL)}
		if(logfile$parameters$prof_comp_link_only=="TRUE"){
			# (1) get in-/directly related isotopologue profiles #################
			prof_isot_IDs<-enviMass::get_isotopol(
				profileList=profileList_neg,
				prof_ID=along[i],
				links_profiles=links_profiles_neg,
				min_peaks=as.numeric(logfile$parameters$corr_min_peaks),
				skip_peaks=as.logical(logfile$parameters$corr_skip_peaks),
				min_cor=as.numeric(logfile$parameters$comp_corr),
				with_test=with_test,
				only_direct=FALSE,
				del_RT=use_del_RT,
				#omit_profiles=in_group
				omit_profiles=FALSE
			)		
			# collect adducts for these isotopologue profiles ####################
			prof_adduct_IDs<-c()
			for(j in 1:length(prof_isot_IDs)){
				got_adducts<-enviMass::get_adducts(
					profileList=profileList_neg,
					prof_ID=prof_isot_IDs[j],
					links_profiles=links_profiles_neg,
					min_peaks=as.numeric(logfile$parameters$corr_min_peaks),
					skip_peaks=as.logical(logfile$parameters$corr_skip_peaks),
					min_cor=as.numeric(logfile$parameters$comp_corr),
					with_test=with_test,
					#omit_profiles=in_group
					omit_profiles=FALSE
				)
				prof_adduct_IDs<-c(prof_adduct_IDs,got_adducts)
			}
			prof_all_IDs<-c(prof_isot_IDs,prof_adduct_IDs)
		}else{
			# (3) get all in-/directly related profiles ##########################
			prof_all_IDs<-enviMass::get_all(
				profileList=profileList_neg,
				prof_ID=along[i],
				links_profiles=links_profiles_neg,
				min_peaks=as.numeric(logfile$parameters$corr_min_peaks),
				skip_peaks=as.logical(logfile$parameters$corr_skip_peaks),
				min_cor=as.numeric(logfile$parameters$comp_corr),
				with_test=with_test,
				only_direct=FALSE,
				del_RT=use_del_RT,
				#omit_profiles=in_group
				omit_profiles=FALSE
			)		
		}
		##########################################################################
		if(plot_it){
			if(length(prof_all_IDs)==1) next
			enviMass::plot_components(
				profileList=profileList_neg,
				prof_IDs=prof_all_IDs,
				links_profiles=links_profiles_neg,
				what=plot_what,
				xlim=FALSE,ylim=FALSE,await_input=TRUE,
				skipit=TRUE,
				min_peaks=as.numeric(logfile$parameters$corr_min_peaks),
				norma=TRUE
			)	
		}	
		##########################################################################
		if(with_test){
			if(!enviMass::check_interlinked(prof_all_IDs,profileList_neg,links_profiles_neg)) stop("Not all interlinked!")
		}
		##########################################################################
		if(length(prof_all_IDs)>1){
			if(with_test){if(prof_all_IDs[1]!=along[i]){stop("\n\nDebug_not_first!")}}	
			at_entry<-profileList_neg[["index_prof"]][prof_all_IDs[1],"links"]
			links_profiles_neg[[at_entry]][["group"]]<-prof_all_IDs[-1]
		}
		##########################################################################	
	}
	if(with_bar){close(pBar)}
	##############################################################################	

	##############################################################################
	# save! ######################################################################
	save(profileList_neg,file=file.path(as.character(logfile[[1]]),"results","profileList_neg"));
	save(links_profiles_neg,file=file.path(as.character(logfile[[1]]),"results","links_profiles_neg"));	
	##############################################################################	
	rm(links_profiles_neg,profileList_neg)

}
