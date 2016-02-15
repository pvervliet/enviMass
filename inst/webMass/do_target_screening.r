	########################################################################################################	
	# clean ! ##############################################################################################
	those<-list.files(file.path(logfile$project_folder,"results","screening"))
	if(length(those)>0){
		for(i in 1:length(those)){
			if(grepl("target",those[i])){ # distinguish with targets
				file.remove(file.path(logfile$project_folder,"results","screening",those[i]))
			}
		}
	}
	########################################################################################################	

	########################################################################################################
	# load available LOD smoothing spline models ###########################################################
	if(file.exists(file=file.path(logfile$project_folder,"results","LOD","LOD_splined"))){
		load(file=file.path(logfile$project_folder,"results","LOD","LOD_splined"));
		do_LOD<-TRUE
	}else{
		do_LOD<-FALSE	
	}		
	########################################################################################################

	########################################################################################################	
	########################################################################################################
	# target screening on positive ionization ##################################################################
	if(any(objects(envir=as.environment(".GlobalEnv"))=="peaklist")){rm(peaklist,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="peaklist")){rm(peaklist)}
	if(any(objects(envir=as.environment(".GlobalEnv"))=="profileList_pos")){rm(profileList_pos,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="profileList_pos")){rm(profileList_pos)}
	if(any(objects(envir=as.environment(".GlobalEnv"))=="pattern_pos_target")){rm(pattern_pos_target,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="pattern_pos_target")){rm(pattern_pos_target)}
	if(any(objects(envir=as.environment(".GlobalEnv"))=="patternRT_pos_target")){rm(patternRT_pos_target,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="patternRT_pos_target")){rm(patternRT_pos_target)}
	if(any(objects(envir=as.environment(".GlobalEnv"))=="patternDelRT_pos_target")){rm(patternDelRT_pos_target,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="patternDelRT_pos_target")){rm(patternDelRT_pos_target)}
	if(
		file.exists(file.path(as.character(logfile[[1]]),"results","profileList_pos")) &
		file.exists(file.path(logfile[[1]],"results","pattern_pos_target"))
	){

		load(file=file.path(as.character(logfile[[1]]),"results","profileList_pos"),envir=as.environment(".GlobalEnv"));	
		load(file=file.path(logfile[[1]],"results","pattern_pos_target"),envir=as.environment(".GlobalEnv"));
		pattern<<-pattern_pos_target;rm(pattern_pos_target,envir=as.environment(".GlobalEnv"));
		load(file=file.path(logfile[[1]],"results","patternRT_pos_target"),envir=as.environment(".GlobalEnv"));
		pattern_RT<<-patternRT_pos_target;rm(patternRT_pos_target,envir=as.environment(".GlobalEnv"));
		load(file=file.path(logfile[[1]],"results","patternDelRT_pos_target"),envir=as.environment(".GlobalEnv"));
		pattern_delRT<<-patternDelRT_pos_target;rm(patternDelRT_pos_target,envir=as.environment(".GlobalEnv"));
		
		mztol<-as.numeric(logfile$parameters$tar_dmz)				# m/z tolerance ...
		ppm<-as.logical(as.character(logfile$parameters$tar_ppm))	# ... given in pppm?
		cutint<-as.numeric(logfile$parameters$tar_intcut)			# Lower intensity threhold
		int_tol<-as.numeric(logfile$parameters$tar_inttol)			# Intensity tolerance %
		RT_tol_outside<-as.numeric(logfile$parameters$tar_drt1)		# RT tolerance of peaks in sample relative to their expected RT [s]
		RT_tol_inside<-as.numeric(logfile$parameters$tar_drt2)		# RT tolerance of peaks within an isotope pattern [s]
		cut_score<-as.numeric(logfile$parameters$tar_w1)	
		
		peaks<-profileList_pos[[7]];
		peaklist<-peaks[,c(14,16,15)];
		# screen centroids
		count_nonmax<-0
		for(i in 1:length(pattern)){
			count_nonmax<-(count_nonmax+
				length(pattern[[i]][,1])
			)
		}
		centro_mass<-rep(0,count_nonmax)
		centro_ID<-rep(0,count_nonmax)
		centro_number<-rep(0,count_nonmax)
		centro_RT<-rep(0,count_nonmax)
		centro_dRT<-rep(0,count_nonmax)
		at_ID<-1
		screen_list<-as.list(rep("FALSE",length(pattern)))
		for(i in 1:length(pattern)){
			n<-length(pattern[[i]][,1])
			centro_mass[at_ID:(at_ID+n-1)]<-pattern[[i]][,1]
			centro_ID[at_ID:(at_ID+n-1)]<-i
			centro_number[at_ID:(at_ID+n-1)]<-(1:n)
			centro_RT[at_ID:(at_ID+n-1)]<-pattern_RT[i]
			centro_dRT[at_ID:(at_ID+n-1)]<-pattern_delRT[i]
			screen_list[[i]]<-as.list(rep("FALSE",n))
			at_ID<-(at_ID+n)
		}
		getit <- search_peak( 
			peaklist, 
			centro_mass, 
			dmz=mztol*4, # precheck
			ppm=ppm, 
			RT=centro_RT, 
			dRT=centro_dRT
		)	
		for(i in 1:length(getit)){ # transfer to a fist list of compoundadduct x centroids
			screen_list[[centro_ID[i]]][[centro_number[i]]]<-getit[i]
		}
		# resort to a full result list: pattern x sample x (centroids,matches) ( = peak index in profileList_pos)
		target_pos_screen_listed<-list()  # default: no match at all
		for(i in 1:length(screen_list)){ # over compound x adduct = i
			if(any(is.na(screen_list[[i]]==FALSE))){
				target_pos_screen_listed[[i]]<-list() # m-level		
				for(j in 1:length(screen_list[[i]])){ # over its centroids = j
					if(screen_list[[i]][[j]]!="FALSE"){ 
						profs<-as.numeric(strsplit(screen_list[[i]][[j]]," / ")[[1]])
						for(k in 1:length(profs)){ # over their matched profile peaks = k
							if(profileList_pos[[7]][profs[k],4]!=profs[k]){cat("\n debug me: profile ID mismatch");stop();} # just a check
							for(m in profileList_pos[[7]][profs[k],1]:profileList_pos[[7]][profs[k],2]){ # over their sample peaks
								delmass<-abs(profileList_pos[[2]][m,1]-pattern[[i]][j,1])		
								if(!ppm){
									if(delmass>mztol){next}
								}else{
									if(delmass*1E6/pattern[[i]][j,1]>mztol){next}
								}
								if(length(target_pos_screen_listed[[i]])<profileList_pos[[2]][m,6][[1]] ){
									target_pos_screen_listed[[i]][[profileList_pos[[2]][m,6][[1]]]]<-matrix(ncol=2,nrow=0)	# sample level
								}
								target_pos_screen_listed[[i]][[profileList_pos[[2]][m,6][[1]]]]<-rbind(
									target_pos_screen_listed[[i]][[profileList_pos[[2]][m,6][[1]]]],c(j,m)
								)
							}							
						}
					}
				}
			}else{
				target_pos_screen_listed[[i]]<-numeric(0)
			}
		}
		# decompose ###########################################################################		
		many<-0
		many_unamb<-0
		res_target_pos_screen<-list()  # default: no match at all
		if(length(target_pos_screen_listed)>0){
			j<-1
			for(i in j:length(target_pos_screen_listed)){ # i - on compound_adduct
				if(length(target_pos_screen_listed[[i]])>0){	
					res_target_pos_screen[[i]]<-list()
					for(m in 1:length(target_pos_screen_listed[[i]])){ # m - sample
						if(length(target_pos_screen_listed[[i]][[m]])>0){
							if(do_LOD){
								with_model<-which(names(LOD_splined)==paste("LOD_",m,sep=""))
								if(length(with_model)>0){						
									use_cutint<-10^(predict(LOD_splined[[with_model]],pattern_RT[i])$y)
								}
							}else{
								use_cutint<-cutint
							}
							combination_matches<-recomb_score(
								cent_peak_mat=target_pos_screen_listed[[i]][[m]],
								pattern_compound=pattern[[i]],
								profileList=profileList_pos,
								LOD=use_cutint,
								RT_tol_inside=RT_tol_inside,
								int_tol=int_tol,
								score_cut=FALSE
							)			
							res_target_pos_screen[[i]][[m]]<-combination_matches
							if(length(combination_matches)>1){many_unamb<-(many_unamb+1)}
							many<-(many+1)
						}
					}
				}else{
					res_target_pos_screen[[i]]<-numeric(0)
				}
			}
			names(res_target_pos_screen)<-names(target_pos_screen_listed)
		}
		# save list ########################################################################################
		save(res_target_pos_screen,file=file.path(logfile$project_folder,"results","screening","res_target_pos_screen"))
		# assemble output table of length(list) ############################################################
		if(length(target_pos_screen_listed)>0){
			measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
			intstand<-read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character");
			results_screen_target_pos<-get_screening_results(
				screened_listed=res_target_pos_screen,
				pattern=pattern,
				profileList=profileList_pos,
				measurements_table=measurements,
				compound_table=intstand,
				cut_score=cut_score
			)
			save(results_screen_target_pos,file=file.path(logfile$project_folder,"results","screening","results_screen_target_pos"))
			rm(measurements,intstand,results_screen_target_pos);		
		}
		####################################################################################################
		rm(getit,target_pos_screen_listed,res_target_pos_screen)
		rm(pattern,pattern_RT,pattern_delRT,envir=as.environment(".GlobalEnv"))
}		
	########################################################################################################
	########################################################################################################
	
			
	########################################################################################################	
	########################################################################################################
	# target screening on negative ionization ##############################################################
	if(any(objects(envir=as.environment(".GlobalEnv"))=="peaklist")){rm(peaklist,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="peaklist")){rm(peaklist)}
	if(any(objects(envir=as.environment(".GlobalEnv"))=="profileList_neg")){rm(profileList_neg,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="profileList_neg")){rm(profileList_neg)}
	if(any(objects(envir=as.environment(".GlobalEnv"))=="pattern_neg_target")){rm(pattern_neg_target,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="pattern_neg_target")){rm(pattern_neg_target)}
	if(any(objects(envir=as.environment(".GlobalEnv"))=="patternRT_neg_target")){rm(patternRT_neg_target,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="patternRT_neg_target")){rm(patternRT_neg_target)}
	if(any(objects(envir=as.environment(".GlobalEnv"))=="patternDelRT_neg_target")){rm(patternDelRT_neg_target,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="patternDelRT_neg_target")){rm(patternDelRT_neg_target)}
	
	if(
		file.exists(file.path(as.character(logfile[[1]]),"results","profileList_neg")) &
		file.exists(file.path(logfile[[1]],"results","pattern_neg_target"))
	){

		load(file=file.path(as.character(logfile[[1]]),"results","profileList_neg"),envir=as.environment(".GlobalEnv"));	
		load(file=file.path(logfile[[1]],"results","pattern_neg_target"),envir=as.environment(".GlobalEnv"));
		pattern<<-pattern_neg_target;rm(pattern_neg_target,envir=as.environment(".GlobalEnv"));
		load(file=file.path(logfile[[1]],"results","patternRT_neg_target"),envir=as.environment(".GlobalEnv"));
		pattern_RT<<-patternRT_neg_target;rm(patternRT_neg_target,envir=as.environment(".GlobalEnv"));
		load(file=file.path(logfile[[1]],"results","patternDelRT_neg_target"),envir=as.environment(".GlobalEnv"));
		pattern_delRT<<-patternDelRT_neg_target;rm(patternDelRT_neg_target,envir=as.environment(".GlobalEnv"));
		
		mztol<-as.numeric(logfile$parameters$tar_dmz)				# m/z tolerance ...
		ppm<-as.logical(as.character(logfile$parameters$tar_ppm))	# ... given in pppm?
		cutint<-as.numeric(logfile$parameters$tar_intcut)			# Lower intensity threhold
		int_tol<-as.numeric(logfile$parameters$tar_inttol)			# Intensity tolerance %
		RT_tol_outside<-as.numeric(logfile$parameters$tar_drt1)		# RT tolerance of peaks in sample relative to their expected RT [s]
		RT_tol_inside<-as.numeric(logfile$parameters$tar_drt2)		# RT tolerance of peaks within an isotope pattern [s]
		cut_score<-as.numeric(logfile$parameters$tar_w1)	
		
		peaks<-profileList_neg[[7]];
		peaklist<-peaks[,c(14,16,15)];
		# screen centroids
		count_nonmax<-0
		for(i in 1:length(pattern)){
			count_nonmax<-(count_nonmax+
				length(pattern[[i]][,1])
			)
		}
		centro_mass<-rep(0,count_nonmax)
		centro_ID<-rep(0,count_nonmax)
		centro_number<-rep(0,count_nonmax)
		centro_RT<-rep(0,count_nonmax)
		centro_dRT<-rep(0,count_nonmax)
		at_ID<-1
		screen_list<-as.list(rep("FALSE",length(pattern)))
		for(i in 1:length(pattern)){
			n<-length(pattern[[i]][,1])
			centro_mass[at_ID:(at_ID+n-1)]<-pattern[[i]][,1]
			centro_ID[at_ID:(at_ID+n-1)]<-i
			centro_number[at_ID:(at_ID+n-1)]<-(1:n)
			centro_RT[at_ID:(at_ID+n-1)]<-pattern_RT[i]
			centro_dRT[at_ID:(at_ID+n-1)]<-pattern_delRT[i]
			screen_list[[i]]<-as.list(rep("FALSE",n))
			at_ID<-(at_ID+n)
		}
		getit <- search_peak(
			peaklist, 
			centro_mass, 
			dmz=mztol*4, # precheck
			ppm=ppm, 
			RT=centro_RT, 
			dRT=centro_dRT
		)	
		for(i in 1:length(getit)){ # transfer to a fist list of compoundadduct x centroids
			screen_list[[centro_ID[i]]][[centro_number[i]]]<-getit[i]
		}
		# resort to a full result list: pattern x sample x (centroids,matches) ( = peak index in profileList_neg)
		target_neg_screen_listed<-list()  # default: no match at all
		for(i in 1:length(screen_list)){ # over compound x adduct = i
			if(any(is.na(screen_list[[i]]==FALSE))){
				target_neg_screen_listed[[i]]<-list() # m-level		
				for(j in 1:length(screen_list[[i]])){ # over its centroids = j
					if(screen_list[[i]][[j]]!="FALSE"){ 
						profs<-as.numeric(strsplit(screen_list[[i]][[j]]," / ")[[1]])
						for(k in 1:length(profs)){ # over their matched profile peaks = k
							if(profileList_neg[[7]][profs[k],4]!=profs[k]){cat("\n debug me: profile ID mismatch");stop();} # just a check
							for(m in profileList_neg[[7]][profs[k],1]:profileList_neg[[7]][profs[k],2]){ # over their sample peaks
								delmass<-abs(profileList_neg[[2]][m,1]-pattern[[i]][j,1])		
								if(!ppm){
									if(delmass>mztol){next}
								}else{
									if(delmass*1E6/pattern[[i]][j,1]>mztol){next}
								}
								if(length(target_neg_screen_listed[[i]])<profileList_neg[[2]][m,6][[1]] ){
									target_neg_screen_listed[[i]][[profileList_neg[[2]][m,6][[1]]]]<-matrix(ncol=2,nrow=0)	# sample level
								}
								target_neg_screen_listed[[i]][[profileList_neg[[2]][m,6][[1]]]]<-rbind(
									target_neg_screen_listed[[i]][[profileList_neg[[2]][m,6][[1]]]],c(j,m)
								)
							}							
						}
					}
				}
			}else{
				target_neg_screen_listed[[i]]<-numeric(0)
			}
		}
		# decompose ###########################################################################		
		many<-0
		many_unamb<-0
		res_target_neg_screen<-list()  # default: no match at all
		if(length(target_neg_screen_listed)>0){
			j<-1
			for(i in j:length(target_neg_screen_listed)){ # i - on compound_adduct
				if(length(target_neg_screen_listed[[i]])>0){	
					res_target_neg_screen[[i]]<-list()
					for(m in 1:length(target_neg_screen_listed[[i]])){ # m - sample
						if(length(target_neg_screen_listed[[i]][[m]])>0){
							if(do_LOD){
								with_model<-which(names(LOD_splined)==paste("LOD_",m,sep=""))
								if(length(with_model)>0){						
									use_cutint<-10^(predict(LOD_splined[[with_model]],pattern_RT[i])$y)
								}
							}else{
								use_cutint<-cutint
							}
							combination_matches<-recomb_score(
								cent_peak_mat=target_neg_screen_listed[[i]][[m]],
								pattern_compound=pattern[[i]],
								profileList=profileList_neg,
								LOD=use_cutint,
								RT_tol_inside=RT_tol_inside,
								int_tol=int_tol,
								score_cut=FALSE
							)
							res_target_neg_screen[[i]][[m]]<-combination_matches
							if(length(combination_matches)>1){many_unamb<-(many_unamb+1)}
							many<-(many+1)
						}
					}
				}else{
					res_target_neg_screen[[i]]<-numeric(0)
				}
			}
			names(res_target_neg_screen)<-names(target_neg_screen_listed)
		}
		# save list ########################################################################################
		save(res_target_neg_screen,file=file.path(logfile$project_folder,"results","screening","res_target_neg_screen"))
		# assemble output table of length(list) ############################################################
		# iterator m is directly equal to the sample ID ####################################################
		if(length(target_neg_screen_listed)>0){
			measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
			intstand<-read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character");
			results_screen_target_neg<-get_screening_results(
				screened_listed=res_target_neg_screen,
				pattern=pattern,
				profileList=profileList_neg,
				measurements_table=measurements,
				compound_table=intstand,
				cut_score=cut_score
			)
			save(results_screen_target_neg,file=file.path(logfile$project_folder,"results","screening","results_screen_target_neg"))
			rm(measurements,intstand,results_screen_target_neg);
		}
		####################################################################################################
		rm(getit,target_neg_screen_listed,res_target_neg_screen)
		rm(pattern,pattern_RT,pattern_delRT,envir=as.environment(".GlobalEnv"))
}		
	########################################################################################################
	########################################################################################################
	
	



