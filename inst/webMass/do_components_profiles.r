# Componentization ##################################################################
###########################################################################################
# REMOVE OLD RESULTS #######################################################################
if(file.exists(file.path(as.character(logfile[[1]]),"results","links_profiles_pos"))){file.remove(file.path(as.character(logfile[[1]]),"results","links_profiles_pos"))}
if(file.exists(file.path(as.character(logfile[[1]]),"results","links_profiles_neg"))){file.remove(file.path(as.character(logfile[[1]]),"results","links_profiles_neg"))}
if(file.exists(file.path(as.character(logfile[[1]]),"results","componentization","components_pos"))){file.remove(file.path(as.character(logfile[[1]]),"results","componentization","components_pos"))}
if(file.exists(file.path(as.character(logfile[[1]]),"results","componentization","components_neg"))){file.remove(file.path(as.character(logfile[[1]]),"results","componentization","components_neg"))}




##################################################################################
# POSITIVE #######################################################################
if( 
	file.exists(file.path(logfile[[1]],"results","profileList_pos")) 
){
#system.time({
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
	use_entries_profiles<-enviMass:::find_empty(links_profiles_pos) # also finds gaps
	profileList_pos[["index_prof"]][,"links"]<-0
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
		pBar <- txtProgressBar(min = 0, max = dim(profileList_pos[["index_prof"]])[1], style = 3)
		for(i in 1:dim(profileList_pos[["index_prof"]])[1]){
			setTxtProgressBar(pBar, i, title = NULL, label = NULL)
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
					links_profiles_pos[[at_entry]]<-enviMass:::new_entry_links_profiles(profileList_pos[["index_prof"]][i,"number_peaks_total"][[1]])
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
		close(pBar)
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
		keep2<-match(forIDs,for_files)
		forIDs<-forIDs[!is.na(keep) & !is.na(keep2)]
		not_found1<-0
		inserted1<-0
		if(length(forIDs)>0){
			pBar <- txtProgressBar(min = 0, max = length(forIDs), style = 3)
			for(i in 1:length(forIDs)){
				setTxtProgressBar(pBar, i, title = NULL, label = NULL)
				load(file=file.path(logfile[[1]],"results","componentization","EIC_corr",forIDs[i]))
				#
				#EIC_pairs<-EIC_pairs[EIC_pairs[,4]>=logfile$parameters$EICor_mincor,,drop=FALSE]		
				#
				if(length(EIC_pairs[,1])==0){next}
				# find profiles for first peak
				get1<-cbind(
					rep(as.numeric(forIDs[i]),length(EIC_pairs[,1])),EIC_pairs[,1] # first column sorted correctly
				)
				found1<-enviMass:::rows_compare(get1,peaks[,c("sampleIDs","peakIDs")],row_order=FALSE,column_order_a=FALSE,column_order_b=FALSE,get_index=TRUE)
				# find profiles for second peak
				get2<-cbind(
					rep(as.numeric(forIDs[i]),length(EIC_pairs[,2])),EIC_pairs[,2] # second column requires sorting
				)
				found2<-enviMass:::rows_compare(get2,peaks[,c("sampleIDs","peakIDs")],row_order=FALSE,column_order_a=TRUE,column_order_b=FALSE,get_index=TRUE)
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
						links_profiles_pos[[at_entry_1]]<-enviMass:::new_entry_links_profiles(profileList_pos[["index_prof"]][prof1,"number_peaks_total"][[1]])
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
						links_profiles_pos[[at_entry_2]]<-enviMass:::new_entry_links_profiles(profileList_pos[["index_prof"]][prof2,"number_peaks_total"][[1]])		
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
			close(pBar)			
			if(any(objects()=="EIC_pairs")){rm(EIC_pairs)}	
		}
		##########################################################################				
		# (2.2) INSERT EIC NON-LINKS #############################################
		if(FALSE){
			forIDs<-profileList_pos[["sampleID"]]
			for_files<-list.files(file.path(logfile[[1]],"results","componentization","EIC_corr"))
			keep<-match(forIDs,for_files) # which files are available?
			if(any(is.na(keep))){cat("\n Just note: not all files found in profiles have EIC correlation results (2). \n")}
			TRUE_IDs<-(measurements$ID[measurements$adducts=="TRUE"]) # for files which have run through that step
			keep2<-match(forIDs,for_files)
			forIDs<-forIDs[!is.na(keep) & !is.na(keep2)]
			not_found2<-0
			inserted2<-0
			if(length(forIDs)>0){
				pBar <- txtProgressBar(min = 0, max = length(forIDs), style = 3)
				for(i in 1:length(forIDs)){
					setTxtProgressBar(pBar, i, title = NULL, label = NULL)
					load(file=file.path(logfile[[1]],"results","componentization","EIC_corr",forIDs[i]))
					EIC_pairs<-EIC_pairs[EIC_pairs[,4]<logfile$parameters$EICor_mincor,,drop=FALSE]		
					if(length(EIC_pairs[,1])==0){next}
					# find profiles for first peak
					get1<-cbind(
						rep(as.numeric(forIDs[i]),length(EIC_pairs[,1])),EIC_pairs[,1] # first column sorted correctly
					)
					found1<-enviMass:::rows_compare(get1,peaks[,c("sampleIDs","peakIDs")],row_order=FALSE,column_order_a=FALSE,column_order_b=FALSE,get_index=TRUE)
					# find profiles for second peak
					get2<-cbind(
						rep(as.numeric(forIDs[i]),length(EIC_pairs[,2])),EIC_pairs[,2] # second column requires sorting
					)
					found2<-enviMass:::rows_compare(get2,peaks[,c("sampleIDs","peakIDs")],row_order=FALSE,column_order_a=TRUE,column_order_b=FALSE,get_index=TRUE)
					for(j in 1:length(found1)){ # insert links
						# insert PROFILE LINKS #######################################
						if(found1[j]==0){not_found2<-(not_found2+1);next} # e.g., peak blind-removed 
						if(found2[j]==0){not_found2<-(not_found2+1);next}			
						# (1) insert link to second profile for the first profile
						prof1<-peaks[found1[j],"profileIDs"][[1]]
						prof2<-peaks[found2[j],"profileIDs"][[1]]
						# profiles have to have positive entries from above
						if(profileList_pos[["index_prof"]][prof1,"links"]==0) next
						if(profileList_pos[["index_prof"]][prof2,"links"]==0) next
						inserted2<-(inserted2+1);
						# enable check ... pairs musst have very similar retention times
						if(FALSE){
							cat("\n");
							cat(peaks[found1[j],"RT"]);cat(" - ")
							cat(peaks[found2[j],"RT"])
						}
						if(profileList_pos[["index_prof"]][prof1,"profile_ID"]!=prof1){stop("\nComponentization: debug me, #1!")}				
						at_entry<-profileList_pos[["index_prof"]][prof1,"links"]
						#if( # no non-link required if no link inserted before
						#	!any(links_profiles_pos[[at_entry]]$EIC[,1]==prof2) 
						#){ next }	
						here<-which(links_profiles_pos[[at_entry]][["EIC"]][,"linked profile"]==prof2)
						links_profiles_pos[[at_entry]][["EIC"]][here,"no-link counts"]<-(links_profiles_pos[[at_entry]][["EIC"]][here,"no-link counts"]+1)
						# (2) insert link to first profile for the second profile
						if(profileList_pos[["index_prof"]][prof2,4]!=prof2){stop("\nComponentization: debug me, #1!")}				
						at_entry<-profileList_pos[["index_prof"]][prof2,"links"] # must exist already - see above no non-link "next"
						# GAPPED?
						here<-which(links_profiles_pos[[at_entry]][["EIC"]][,"linked profile"]==prof1)
						links_profiles_pos[[at_entry]][["EIC"]][here,"no-link counts"]<-(links_profiles_pos[[at_entry]][["EIC"]][here,"no-link counts"]+1)		
						# insert PEAK LINKS ? ########################################
					}
				}
				close(pBar)
				if(any(objects()=="EIC_pairs")){rm(EIC_pairs)}
		}
		}
		##########################################################################				
		# (2.3) INSERT CO-OCCURENCES #############################################
		if(FALSE){
			if(length(forIDs)>0){ # from above
				pBar <- txtProgressBar(min = 0, max = length(forIDs), style = 3)
				for(i in 1:length(links_profiles_pos)){
					setTxtProgressBar(pBar, i, title = NULL, label = NULL)
					if(dim(links_profiles_pos[[i]]$EIC)[1]>0){
						for(m in 1:dim(links_profiles_pos[[i]]$EIC)[1]){
							if(links_profiles_pos[[i]]$EIC[m,c("ref_1")]!=0) next; # done via linked profile
							prof1<-as.numeric(names(links_profiles_pos)[i])
							prof2<-links_profiles_pos[[i]]$EIC[m,"linked profile"]
							del1<-profileList_pos[["index_prof"]][prof1,"number_peaks_total"][[1]]
							del2<-profileList_pos[["index_prof"]][prof2,"number_peaks_total"][[1]]
							summed<-sum(links_profiles_pos[[i]]$EIC[m,c("link counts","no-link counts")])
							at_entry<-profileList_pos[["index_prof"]][prof2,"links"]
							here<-which(links_profiles_pos[[at_entry]][["EIC"]][,"linked profile"]==prof1)
							# which peaks co-occur in the same sample?
							these<-profileList_pos[["peaks"]][
								profileList_pos[["index_prof"]][prof1,"start_ID"]:profileList_pos[["index_prof"]][prof1,"end_ID"]
							,"sampleIDs"]
							those<-profileList_pos[["peaks"]][
								profileList_pos[["index_prof"]][prof2,"start_ID"]:profileList_pos[["index_prof"]][prof2,"end_ID"]
							,"sampleIDs"]
							matched<-match(these,those)
							links_profiles_pos[[i]]$EIC[m,"ref_1"]<-sum(!is.na(matched))
							RT_1<-(profileList_pos[["peaks"]][
									(profileList_pos[["index_prof"]][prof1,"start_ID"]:profileList_pos[["index_prof"]][prof1,"end_ID"])
									[!is.na(matched)]
								,"RT"])
							RT_2<-(profileList_pos[["peaks"]][
									(profileList_pos[["index_prof"]][prof2,"start_ID"]:profileList_pos[["index_prof"]][prof2,"end_ID"])
									[matched[!is.na(matched)]]
								,"RT"])
							del_RT<-abs(RT_1-RT_2)
							links_profiles_pos[[i]]$EIC[m,"ref_2"]<-sum(del_RT<=as.numeric(logfile$parameters$isotop_rttol))
							links_profiles_pos[[i]]$EIC[m,"ref_3"]<-sum(del_RT<=as.numeric(logfile$parameters$adducts_rttol))
							# insert values in other linked profile
							at_entry<-profileList_pos[["index_prof"]][prof2,"links"]
							here<-which(links_profiles_pos[[at_entry]][["EIC"]][,"linked profile"]==prof1)
							links_profiles_pos[[at_entry]]$EIC[here,c("ref_1","ref_2","ref_3")]<-links_profiles_pos[[i]]$EIC[m,c("ref_1","ref_2","ref_3")]
						}
					}
				}
				close(pBar)
			}
		}
		##########################################################################
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
		not_found3<-0
		inserted3<-0
		if(length(forIDs)>0){
			cat("\n Retrieving isotopologue links ")
			pBar <- txtProgressBar(min = 0, max = length(forIDs), style = 3)
			for(i in 1:length(forIDs)){
				setTxtProgressBar(pBar, i, title = NULL, label = NULL)
				load(file=file.path(logfile[[1]],"results","componentization","isotopologues",forIDs[i]))
				if(length(Isot_pairs[,1])==0){next}
				# find profiles for first peak
				get1<-cbind(
					rep(as.numeric(forIDs[i]),length(Isot_pairs[,1])),Isot_pairs[,1]
				)
				found1<-enviMass:::rows_compare(get1,peaks[,c("sampleIDs","peakIDs")],row_order=FALSE,column_order_a=FALSE,column_order_b=FALSE,get_index=TRUE)
				# find profiles for second peak
				get2<-cbind(
					rep(as.numeric(forIDs[i]),length(Isot_pairs[,2])),Isot_pairs[,2]
				)
				found2<-enviMass:::rows_compare(get2,peaks[,c("sampleIDs","peakIDs")],row_order=FALSE,column_order_a=TRUE,column_order_b=FALSE,get_index=TRUE)
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
						links_profiles_pos[[at_entry_1]]<-enviMass:::new_entry_links_profiles(profileList_pos[["index_prof"]][prof1,"number_peaks_total"][[1]])
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
						links_profiles_pos[[at_entry_2]]<-enviMass:::new_entry_links_profiles(profileList_pos[["index_prof"]][prof2,"number_peaks_total"][[1]])
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
			close(pBar)
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
		not_found4<-0
		inserted4<-0
		if(length(forIDs)>0){
			cat("\n Retrieving adduct links ")
			pBar <- txtProgressBar(min = 0, max = length(forIDs), style = 3)
			for(i in 1:length(forIDs)){
				setTxtProgressBar(pBar, i, title = NULL, label = NULL)
				load(file=file.path(logfile[[1]],"results","componentization","adducts",forIDs[i]))
				if(length(Adduct_pairs[,1])==0){next}
				# find profiles for first peak
				get1<-cbind(
					rep(as.numeric(forIDs[i]),length(Adduct_pairs[,1])),Adduct_pairs[,1]
				)
				found1<-enviMass:::rows_compare(get1,peaks[,c("sampleIDs","peakIDs")],row_order=FALSE,column_order_a=FALSE,column_order_b=FALSE,get_index=TRUE)
				# find profiles for second peak
				get2<-cbind(
					rep(as.numeric(forIDs[i]),length(Adduct_pairs[,2])),Adduct_pairs[,2]
				)
				found2<-enviMass:::rows_compare(get2,peaks[,c("sampleIDs","peakIDs")],row_order=FALSE,column_order_a=TRUE,column_order_b=FALSE,get_index=TRUE)
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
						links_profiles_pos[[at_entry_1]]<-enviMass:::new_entry_links_profiles(profileList_pos[["index_prof"]][prof1,"number_peaks_total"][[1]])
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
						links_profiles_pos[[at_entry_2]]<-enviMass:::new_entry_links_profiles(profileList_pos[["index_prof"]][prof2,"number_peaks_total"][[1]])
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
			close(pBar)
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
		not_found5<-0
		inserted5<-0
		if(length(forIDs)>0){
			cat("\n Retrieving homologue links ")
			pBar <- txtProgressBar(min = 0, max = length(forIDs), style = 3)
			for(i in 1:length(forIDs)){
				setTxtProgressBar(pBar, i, title = NULL, label = NULL)
				load(file=file.path(logfile[[1]],"results","componentization","homologues",forIDs[i]))	
				if(length(Homol_groups[,1])==0){next}	
				keep_group<-rep(TRUE,max(Homol_groups[,3])) # make sure all peaks in a group are still present
				# find profiles for first peak
				get1<-cbind(
					rep(as.numeric(forIDs[i]),length(Homol_groups[,1])),Homol_groups[,1]
				)
				found1<-enviMass:::rows_compare(get1,peaks[,c("sampleIDs","peakIDs")],row_order=FALSE,column_order_a=FALSE,column_order_b=FALSE,get_index=TRUE)
				keep_group[Homol_groups[found1==0,3]]<-FALSE
				if(!any(keep_group)){next} # none left with all peaks still embedded into profiles
				# find profiles for second peak
				get2<-cbind(
					rep(as.numeric(forIDs[i]),length(Homol_groups[,2])),Homol_groups[,2]
				)
				found2<-enviMass:::rows_compare(get2,peaks[,c("sampleIDs","peakIDs")],row_order=FALSE,column_order_a=TRUE,column_order_b=FALSE,get_index=TRUE)
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
						links_profiles_pos[[at_entry]]<-enviMass:::new_entry_links_profiles(profileList_pos[["index_prof"]][prof1,"number_peaks_total"][[1]])
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
						links_profiles_pos[[at_entry]]<-enviMass:::new_entry_links_profiles(profileList_pos[["index_prof"]][prof2,"number_peaks_total"][[1]])
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
			close(pBar)
			if(any(objects()=="Homol_groups")){rm(Homol_groups)}
			cat(" done.")
		}
	}
	##############################################################################	

	##############################################################################	
	# (6) filter #################################################################
	# (6.1) by ISTD - first get their characteristics on delRT and correl. #######
	fil1<-enviMass:::analyseA_links_profiles(
			links_profiles=links_profiles_pos, profileList=profileList_pos, 
			min_rat=.7, 	# isot, adduc: "link counts"/"ref_1"
			min_count=.4, 	# isot, adduc:  "ref_1">=(min_count*number_samples)
			for_which="ISTD"
		)
	# clean isotopologues ########################################################
	cut_delRT_isot<-median(fil1$delRT_isot)
	cut_cor_isot<-(boxplot.stats(c(fil1$int_cor_isot))$conf[1])
	links_profiles_pos<-enviMass:::cleanA_links_profiles(
		links_profiles=links_profiles_pos, profileList=profileList_pos,
		cut_delRT_isot = cut_delRT_isot, 
		cut_cor_isot = cut_cor_isot, 
		cut_frac_iso = .9
	)
	# clean adducts ##############################################################
	cut_delRT_adduc<-median(fil1$delRT_adduc)
	links_profiles_pos<-enviMass:::cleanB_links_profiles( 
		links_profiles=links_profiles_pos, profileList=profileList_pos,
		cut_delRT_adduc = cut_delRT_adduc, 
		cut_frac_adduc = .9
	)
	# (6.2) by ISTD - check their EIC correlation ################################
	fil2<-enviMass:::analyseB_links_profiles(
			links_profiles=links_profiles_pos,  
			min_count=.4, 	# isot, adduc:  "ref_1">=(min_count*number_samples)
			for_which="ISTD"
		)
	use_EIC<-c(fil2$EIC_cor_isot,fil2$EIC_cor_adduc)
	cut_EIC<-(boxplot.stats(use_EIC)$stats[1])



#})	

	##############################################################################
	# save! ######################################################################
	save(profileList_pos,file=file.path(as.character(logfile[[1]]),"results","profileList_pos"));
	save(links_profiles_pos,file=file.path(as.character(logfile[[1]]),"results","links_profiles_pos"));	
	##############################################################################	
	rm(links_profiles_pos,profileList_pos,links_peaks_pos)

}






##################################################################################
# Negative #######################################################################
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
	# links_peaks_neg<-list(); # each entry with 6 lists itself: targets, IS, EIC_correl, isotop, adducts, homol
	load(file.path(as.character(logfile[[1]]),"results","links_peaks_neg"));	
	links_profiles_neg<-list(); # each entry with 6 lists itself: targets, IS, EIC_correl, isotop, adducts, homol
	measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
	peaks<-profileList_neg[["peaks"]][,c("sampleIDs","peakIDs","profileIDs","RT")] # to retrieve relations with, sampleID, peakID, profileID, RT
	ord<-order(peaks[,"sampleIDs"],peaks[,"peakIDs"],peaks[,"profileIDs"],decreasing=FALSE)
	peaks<-peaks[ord,]
	use_entries_profiles<-enviMass:::find_empty(links_profiles_neg)
	profileList_neg[["index_prof"]][,"links"]<-0
	##############################################################################	

	##############################################################################	
	# annotate target & IS screening matches stored in links_peaks_neg ###########
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
		cat("\n  Annotation of screening results to profiles")
		for(i in 1:length(profileList_neg[[7]][,1])){
			if(
				any(profileList_neg[[2]][
					(profileList_neg[[7]][i,"start_ID"]:profileList_neg[[7]][i,"end_ID"]),"links"
				]!=0)			
			){
			
				###################################################################
				# add a new link to the profile ###################################
				if( profileList_neg[[7]][i,"links"]==0 ){ 	# establish a new link ...
					if(length(use_entries_profiles)>0){
						at_entry<-use_entries_profiles[1]
						use_entries<-use_entries_profiles[-1]
					}else{
						at_entry<-(length(links_profiles_neg)+1)
					}
					links_profiles_neg[[at_entry]]<-list()
					links_profiles_neg[[at_entry]][[1]]<-list() # target
					links_profiles_neg[[at_entry]][[2]]<-list()	# IS
					links_profiles_neg[[at_entry]][[3]]<-list()	# EIC_correl
						links_profiles_neg[[at_entry]][[3]]<-matrix(ncol=4,nrow=0)
						colnames(links_profiles_neg[[at_entry]][[3]])<-c("linked profile","link counts","ref","no-link counts")
					links_profiles_neg[[at_entry]][[4]]<-list()	# isotop
						links_profiles_neg[[at_entry]][[4]]<-matrix(ncol=3,nrow=0)
						colnames(links_profiles_neg[[at_entry]][[4]])<-c("linked profile","link counts","ref")					
					links_profiles_neg[[at_entry]][[5]]<-list()	# adducts
						links_profiles_neg[[at_entry]][[5]]<-matrix(ncol=3,nrow=0)
						colnames(links_profiles_neg[[at_entry]][[5]])<-c("linked profile","link counts","ref")					
					links_profiles_neg[[at_entry]][[6]]<-list()	# homol		
						links_profiles_neg[[at_entry]][[6]]<-matrix(ncol=3,nrow=0)
						colnames(links_profiles_neg[[at_entry]][[6]])<-c("linked profile","link counts","ref")
					links_profiles_neg[[at_entry]][[7]]<-list()	# group		
					names(links_profiles_neg[[at_entry]])<-c("targ","IS","EIC","isot","adduc","homol","group")
					names(links_profiles_neg)[at_entry]<-as.character(i)
					profileList_neg[[7]][i,"links"]<-at_entry						
				}else{
					at_entry<-profileList_neg[[7]][i,"links"]
				}				
				###################################################################
				# search peaks & their links to compounds #########################
				for(j in (profileList_neg[[7]][i,"start_ID"]:profileList_neg[[7]][i,"end_ID"])){
					if(profileList_neg[[2]][j,"links"]!=0){
						# add IS link #############################################
						if(	length(links_peaks_neg[[profileList_neg[[2]][j,"links"]]][[2]])>0 ){
							for(k in 1:length(links_peaks_neg[[profileList_neg[[2]][j,"links"]]][[2]])){ # if several compound matches exist for this peak
								if( length(links_profiles_neg[[at_entry]][[2]])==0 ){ # make a new entry for profile link to IS
									links_profiles_neg[[at_entry]][[2]]<-
										data.frame(
											links_peaks_neg[[profileList_neg[[2]][j,"links"]]][[2]][[k]],1,stringsAsFactors = FALSE
										)
									names(links_profiles_neg[[at_entry]][[2]])<-c("Compound","Counts")
								}else{
									at<-which(links_profiles_neg[[at_entry]][[2]][,1]==links_peaks_neg[[profileList_neg[[2]][j,"links"]]][[2]][[k]])
									if(length(at)>0){ 	# increment existing link ...
										links_profiles_neg[[at_entry]][[2]][at,2]<-(links_profiles_neg[[at_entry]][[2]][at,2]+1)
									}else{	# or just add a new one?
										links_profiles_neg[[at_entry]][[2]]<-data.frame(
											c(links_profiles_neg[[at_entry]][[2]][,1], links_peaks_neg[[profileList_neg[[2]][j,"links"]]][[2]][[k]]),
											c(links_profiles_neg[[at_entry]][[2]][,2],1),
											stringsAsFactors = FALSE
										)									
										names(links_profiles_neg[[at_entry]][[2]])<-c("Compound","Counts")
									}
								}
							}	
						}					
						# add target link #########################################
						if(	length(links_peaks_neg[[profileList_neg[[2]][j,"links"]]][[1]])>0 ){
							for(k in 1:length(links_peaks_neg[[profileList_neg[[2]][j,"links"]]][[1]])){ # if several compound matches exist for this peak
								if( length(links_profiles_neg[[at_entry]][[1]])==0 ){ # make a new entry for profile link to IS
									links_profiles_neg[[at_entry]][[1]]<-
										data.frame(
											links_peaks_neg[[profileList_neg[[2]][j,"links"]]][[1]][[k]],1,stringsAsFactors = FALSE
										)
									names(links_profiles_neg[[at_entry]][[1]])<-c("Compound","Counts")
								}else{
									at<-which(links_profiles_neg[[at_entry]][[1]][,1]==links_peaks_neg[[profileList_neg[[2]][j,"links"]]][[1]][[k]])
									if(length(at)>0){ 	# increment existing link ...
										links_profiles_neg[[at_entry]][[1]][at,2]<-(links_profiles_neg[[at_entry]][[1]][at,2]+1)
									}else{	# or just add a new one?
										links_profiles_neg[[at_entry]][[1]]<-data.frame(
											c(links_profiles_neg[[at_entry]][[1]][,1], links_peaks_neg[[profileList_neg[[2]][j,"links"]]][[1]][[k]]),
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
	}
	##############################################################################	
	
	##############################################################################	
	# for each profile, filter all relations to other profiles ###################
			
	##############################################################################	
	
	##############################################################################
	# save! ######################################################################
	save(profileList_neg,file=file.path(as.character(logfile[[1]]),"results","profileList_neg"));
	save(links_profiles_neg,file=file.path(as.character(logfile[[1]]),"results","links_profiles_neg"));	
	##############################################################################	
	rm(links_profiles_neg,profileList_neg,links_peaks_neg)

}

