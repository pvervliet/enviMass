
    measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
	measurements<-measurements[measurements[,"include"]=="TRUE",]
	if( 
		any(measurements[,"Mode"]=="positive") & 
		file.exists(file.path(as.character(logfile[[1]]),"results","pattern_pos_IS"))  &
		file.exists(file.path(as.character(logfile[[1]]),"results","profileList_pos"))	&
		as.logical(logfile$parameters$ISnorm_include_pos)		
	){
	
		if(any(objects(envir=as.environment(".GlobalEnv"))=="peaklist")){rm(peaklist,envir=as.environment(".GlobalEnv"))}
		if(any(objects()=="peaklist")){rm(peaklist)}
		if(any(objects(envir=as.environment(".GlobalEnv"))=="profileList_pos")){rm(profileList_pos,envir=as.environment(".GlobalEnv"))}
		if(any(objects()=="profileList_pos")){rm(profileList_pos)}
		if(any(objects(envir=as.environment(".GlobalEnv"))=="pattern_pos_IS")){rm(pattern_pos_IS,envir=as.environment(".GlobalEnv"))}
		if(any(objects()=="pattern_pos_IS")){rm(pattern_pos_IS)}
		if(any(objects(envir=as.environment(".GlobalEnv"))=="patternRT_pos_IS")){rm(patternRT_pos_IS,envir=as.environment(".GlobalEnv"))}
		if(any(objects()=="patternRT_pos_IS")){rm(patternRT_pos_IS)}
		if(any(objects(envir=as.environment(".GlobalEnv"))=="patternDelRT_pos_IS")){rm(patternDelRT_pos_IS,envir=as.environment(".GlobalEnv"))}
		if(any(objects()=="patternDelRT_pos_IS")){rm(patternDelRT_pos_IS)}
		load(file=file.path(as.character(logfile[[1]]),"results","profileList_pos"),envir=as.environment(".GlobalEnv"));	
		load(file=file.path(logfile[[1]],"results","pattern_pos_IS"),envir=as.environment(".GlobalEnv"));
		load(file=file.path(logfile[[1]],"results","patternRT_pos_IS"),envir=as.environment(".GlobalEnv"));
		load(file=file.path(logfile[[1]],"results","patternDelRT_pos_IS"),envir=as.environment(".GlobalEnv"));
		peaks<-profileList_pos[[7]];
		peaklist<-peaks[,c(14,16,15)];
		treshscore<-as.numeric(logfile$parameters$ISnorm_score_pos)
		we1=(1-as.numeric(logfile$parameters$IS_w1))
		we2=(1-we1)
		cat("- screening:");
		screened<-screening(	peaklist, 
								blanklist=FALSE, 
								pattern_pos_IS, 
								RT = patternRT_pos_IS,
								dmz = as.numeric(logfile$parameters$IS_dmz), 
								ppm = logfile$parameters$IS_ppm, 
								dRT = patternDelRT_pos_IS, 
								dRTwithin = as.numeric(logfile$parameters$IS_drt2), 
								dRTblank = FALSE, 
								dInt = as.numeric(logfile$parameters$IS_inttol), 
								Intcut = as.numeric(logfile$parameters$IS_intcut), 
								w1=we1, 
								w2=we2, 	
								w3=0	
		)
		# set matrix to sort & store data from a profile ###################################
		atPOSIX<-profileList_pos[[3]];
		sampletype<-profileList_pos[[9]];
		sampleID<-profileList_pos[[4]];
		# filter out other file types such as spiked ones
		keep<-((sampletype=="sample")|(sampletype=="blank"))
		atPOSIX<-atPOSIX[keep]
		sampletype<-sampletype[keep]
		sampleID<-sampleID[keep]
		#
		atdate<-c();
		attime<-c();
		for(i in 1:length(atPOSIX)){
				  atdate<-c(atdate, strsplit(atPOSIX[i]," ")[[1]][1]);
				  attime<-c(attime, strsplit(atPOSIX[i]," ")[[1]][2]);
		}
		attime<-as.difftime(attime);
		atdate<-as.Date(atdate, tz="GMT");
		ord<-order(as.numeric(atdate),as.numeric(attime),sampleID);
		atPOSIXsort<-atPOSIX[ord];
		atdate<-atdate[ord];
		attime<-attime[ord];
		sampleID<-sampleID[ord];
		sampletype<-sampletype[ord];
		# to retrieve data
		timeset<-matrix(nrow=length(atPOSIX),ncol=5,0);
		for(i in 1:length(sampleID)){
			  if(sampletype[i]=="sample"){
				timeset[i,2]<-as.numeric(sampleID[i]);
			  }
			  if(sampletype[i]=="blank"){
				timeset[i,3]<-as.numeric(sampleID[i]);
			  }
		}
		# screen IS intensity profiles #####################################################
		cat("- on IS profiles");
		min_count<-floor(length(profileList_pos[[4]])*as.numeric(logfile$parameters$ISnorm_percfiles_pos)/100);
		lis_delint_IS<-list(0)
		lis_median_IS<-list(0)
		stillin<-rep(TRUE,length(peaks[,1]))
		for(p in 1:length(timeset[,1])){
			lis_delint_IS[[p]]<-numeric(0)
			lis_median_IS[[p]]<-numeric(0)
		}
		for(i in 1:length(screened)){
			if(length(screened[[i]])!=1){
				for(j in 1:length(screened[[i]][,1])){
					if(  !(as.character(screened[[i]][j,9]))=="NaN"  ){
						if( as.numeric(as.character(screened[[i]][j,9]))>=treshscore ){	
							hits<-strsplit(as.character(screened[[i]][j,2]),"/")[[1]];
							hits<-hits[hits!="-"];
							hits<-as.numeric(hits);
							for(b in 1:length(hits)){
								if( peaks[hits[b],3]>min_count ){
									profID<-as.numeric(peaks[hits[b],4])
									timeset[,4:5]<-0;
									timeset[,c(4,5)] <-.Call("_enviMass_fill_timeset",
														as.numeric(timeset),
														as.numeric(profileList_pos[[2]][(profileList_pos[[7]][profileList_pos[[7]][,4]==profID,1]:profileList_pos[[7]][profileList_pos[[7]][,4]==profID,2]),6]), # sampleIDs
														as.numeric(profileList_pos[[2]][(profileList_pos[[7]][profileList_pos[[7]][,4]==profID,1]:profileList_pos[[7]][profileList_pos[[7]][,4]==profID,2]),2]), # intensities
														as.integer(length(timeset[,1])),
														PACKAGE="enviMass"
									)
									timeset[timeset[,4]!=0,4]<-log10(timeset[timeset[,4]!=0,4])
									timeset[timeset[,5]!=0,5]<-log10(timeset[timeset[,5]!=0,5])
									median_int<-c(timeset[timeset[,4]!=0,4],timeset[timeset[,5]!=0,5])
									median_int<-median(median_int);
									for(k in 1:length(timeset[,1])){
										if(timeset[k,5]!=0){
											lis_delint_IS[[k]]<-c(lis_delint_IS[[k]],(timeset[k,5]-median_int))
											lis_median_IS[[k]]<-c(lis_median_IS[[k]],median_int)									
											next;
										}
										if(timeset[k,4]!=0){
											lis_delint_IS[[k]]<-c(lis_delint_IS[[k]],(timeset[k,4]-median_int))
											lis_median_IS[[k]]<-c(lis_median_IS[[k]],median_int)	
										}
									}
									stillin[hits[b]]<-FALSE									
								}	
							}
						}
					}	
				}
			}
		}
		# screen other profiles ###############################################################
		if( (logfile$parameters$ISnorm_medblank_pos=="TRUE" || logfile$parameters$ISnorm_medsam_pos=="TRUE") ){
			cat("- on blank / non-blank profiles");
			peaks<-peaks[stillin,];
			peaks<-peaks[sample(1:length(peaks[,1]),length(peaks[,1]), replace = FALSE),];
			if(logfile$parameters$ISnorm_medblank_pos=="TRUE"){ # use blank
				count_b<-0
				if( logfile$parameters$ISnorm_usesubblank_pos=="TRUE" ){ # use subsampling
					max_count_b<-as.numeric(logfile$parameters$ISnorm_numblank_pos)
				}else{
					max_count_b<-length(peaks[,1])
				}
			}else{
				count_b<-0
				max_count_b<-0
			}
			if(logfile$parameters$ISnorm_medsam_pos=="TRUE"){ # use nonblank
				count_nb<-0
				if( logfile$parameters$ISnorm_usesubsam_pos=="TRUE" ){ # use subsampling
					max_count_nb<-as.numeric(logfile$parameters$ISnorm_numsam_pos)
				}else{
					max_count_nb<-length(peaks[,1])
				}
			}else{
				count_nb<-0
				max_count_nb<-0
			}
			lis_delint_nb<-list(0)
			lis_median_nb<-list(0)
			lis_delint_b<-list(0)
			lis_median_b<-list(0)			
			stillin<-rep(TRUE,length(peaks[,1]))
			for(p in 1:length(timeset[,1])){
				lis_delint_nb[[p]]<-numeric(0)
				lis_median_nb[[p]]<-numeric(0)
				lis_delint_b[[p]]<-numeric(0)
				lis_median_b[[p]]<-numeric(0)
			}
			for(i in 1:length(peaks[,1])){		
				profID<-as.numeric(peaks[i,4])
				timeset[,4:5]<-0;
				timeset[,c(4,5)] <-.Call("_enviMass_fill_timeset",
					as.numeric(timeset),
					as.numeric(profileList_pos[[2]][(profileList_pos[[7]][profileList_pos[[7]][,4]==profID,1]:profileList_pos[[7]][profileList_pos[[7]][,4]==profID,2]),6]), # sampleIDs
					as.numeric(profileList_pos[[2]][(profileList_pos[[7]][profileList_pos[[7]][,4]==profID,1]:profileList_pos[[7]][profileList_pos[[7]][,4]==profID,2]),2]), # intensities
					as.integer(length(timeset[,1])),
					PACKAGE="enviMass"
				)
				if(	 # for blind #############################################3
					(logfile$parameters$ISnorm_medblank_pos=="TRUE") &&
					(any(timeset[,5]>0)) &&
					(count_b<max_count_b)
				){ 
					timeset[timeset[,4]!=0,4]<-log10(timeset[timeset[,4]!=0,4])
					timeset[timeset[,5]!=0,5]<-log10(timeset[timeset[,5]!=0,5])
					median_int<-c(timeset[timeset[,4]!=0,4],timeset[timeset[,5]!=0,5])
					median_int<-median(median_int);
					for(k in 1:length(timeset[,1])){
						if(timeset[k,5]!=0){
							lis_delint_b[[k]]<-c(lis_delint_b[[k]],(timeset[k,5]-median_int))
							lis_median_b[[k]]<-c(lis_median_b[[k]],median_int)									
							next;
						}
						if(timeset[k,4]!=0){
							lis_delint_b[[k]]<-c(lis_delint_b[[k]],(timeset[k,4]-median_int))
							lis_median_b[[k]]<-c(lis_median_b[[k]],median_int)	
						}
					}
					count_b<-(count_b+1);
					next;
				}
				if(	 # for non-blind ##########################################
					(logfile$parameters$ISnorm_medsam_pos=="TRUE") &&
					(any(timeset[,4]>0)) &&
					(count_nb<max_count_nb)
				){ 
					timeset[timeset[,4]!=0,4]<-log10(timeset[timeset[,4]!=0,4])
					median_int<-median(timeset[timeset[,4]!=0,4]);
					for(k in 1:length(timeset[,1])){
						if(timeset[k,4]!=0){
							lis_delint_nb[[k]]<-c(lis_delint_nb[[k]],(timeset[k,4]-median_int))
							lis_median_nb[[k]]<-c(lis_median_nb[[k]],median_int)									
						}
					}
					count_nb<-(count_nb+1);
				}
				if( # enough subsampling?
					(count_b>=max_count_b) &&
					(count_nb>=max_count_nb) 
				){
					break
				}
			}
		}
		# correct intensities - replace in profileList, recalculate mean_int ###############
		cat("- correcting intensities");
		corfac<-c()
		for(k in 1:length(lis_delint_IS)){
			if( length(lis_delint_IS[[k]])>as.numeric(logfile$parameters$ISnorm_numbIS_pos) ){
				corfac<-c( corfac,10^median(lis_delint_IS[[k]]) )
			}else{			
				corfac<-c(corfac,1)
			}
		}
		corr_intens <- .Call(	"_enviMass_correct_intens",
								as.numeric(corfac),	  # correction factor
								as.integer(sampleID),       
								as.numeric(profileList_pos[[2]][,2]), # intensities
								as.integer(profileList_pos[[2]][,6]),  
								PACKAGE="enviMass"
							)
		profileList_pos[[2]][,2]<-corr_intens
		for(k in 1:length(profileList_pos[[7]][,8])){
			profileList_pos[[7]][k,16]<-mean(profileList_pos[[2]][(profileList_pos[[7]][k,1]:profileList_pos[[7]][k,2]),2])
		}
		save(profileList_pos,file=file.path(as.character(logfile[[1]]),"results","profileList_pos"),compress=FALSE);	
		
# > BAUSTELLE

		# derive plots #####################################################################
		int_norm_ISTD_pos <<- list()
		int_norm_ISTD_pos[[1]] <<- lis_delint_IS
		int_norm_ISTD_pos[[2]] <<- lis_median_IS
		#int_norm_ISTD_pos[[3]] <- lis_RT_IS
		int_norm_ISTD_pos[[4]] <<- use_corfac
		int_norm_ISTD_pos[[5]] <<- list()#lis_delint_nb
		int_norm_ISTD_pos[[6]] <<- list()#lis_median_nb
		int_norm_ISTD_pos[[7]] <<- list()#lis_delint_b
		int_norm_ISTD_pos[[8]] <<- list()#lis_median_b
		int_norm_ISTD_pos[[9]] <<- profileList_pos[["datetime"]]
		int_norm_ISTD_pos[[10]] <<- profileList_pos[["type"]]
		int_norm_ISTD_pos[[11]] <<- profileList_pos[["sampleID"]]	
		names(int_norm_ISTD_pos) <<- c("lis_delint_IS", " lis_median_IS", "lis_RT_IS", "use_corfac", "lis_delint_nb", 
			"lis_median_nb", "lis_delint_b", "lis_median_b", "atPOSIX", "sampletype", "sampleID")
		# -> save data & derive plots ######################################################
		output$int_norm_ISTD_pos_median <- renderPlot({   
			par(mar = c(.2, 4.5, .9, 8))
			enviMass:::plot_ISTD_norm(
				int_norm_ISTD = int_norm_ISTD_pos,
				logfile = logfile,
				what = "normalization"
			)
		},res = 100) 
		output$int_norm_ISTD_pos_counts <- renderPlot({   
			par(mar = c(4.5, 4.5, .9, 8))
			enviMass:::plot_ISTD_norm(
				int_norm_ISTD = int_norm_ISTD_pos,
				logfile = logfile,
				what = "counts"
			)
		},res = 100) 
		save(int_norm_ISTD_pos, file = file.path(as.character(logfile[[1]]), "results", "int_norm_ISTD_pos"));
		####################################################################################

# < BAUSTELLE

	}

	########################################################################################
	if( 
		any(measurements[,"Mode"]=="negative") & 
		file.exists(file.path(as.character(logfile[[1]]),"results","pattern_neg_IS")) &
		file.exists(file.path(as.character(logfile[[1]]),"results","profileList_neg")) &
		as.logical(logfile$parameters$ISnorm_include_neg)		
	){
		if(any(objects(envir=as.environment(".GlobalEnv"))=="peaklist")){rm(peaklist,envir=as.environment(".GlobalEnv"))}
		if(any(objects()=="peaklist")){rm(peaklist)}
		if(any(objects(envir=as.environment(".GlobalEnv"))=="profileList_neg")){rm(profileList_neg,envir=as.environment(".GlobalEnv"))}
		if(any(objects()=="profileList_neg")){rm(profileList_neg)}
		if(any(objects(envir=as.environment(".GlobalEnv"))=="pattern_neg_IS")){rm(pattern_neg_IS,envir=as.environment(".GlobalEnv"))}
		if(any(objects()=="pattern_neg_IS")){rm(pattern_neg_IS)}
		if(any(objects(envir=as.environment(".GlobalEnv"))=="patternRT_neg_IS")){rm(patternRT_neg_IS,envir=as.environment(".GlobalEnv"))}
		if(any(objects()=="patternRT_neg_IS")){rm(patternRT_neg_IS)}
		if(any(objects(envir=as.environment(".GlobalEnv"))=="patternDelRT_neg_IS")){rm(patternDelRT_neg_IS,envir=as.environment(".GlobalEnv"))}
		if(any(objects()=="patternDelRT_neg_IS")){rm(patternDelRT_neg_IS)}
		load(file=file.path(as.character(logfile[[1]]),"results","profileList_neg"),envir=as.environment(".GlobalEnv"));	
		load(file=file.path(logfile[[1]],"results","pattern_neg_IS"),envir=as.environment(".GlobalEnv"));
		load(file=file.path(logfile[[1]],"results","patternRT_neg_IS"),envir=as.environment(".GlobalEnv"));
		load(file=file.path(logfile[[1]],"results","patternDelRT_neg_IS"),envir=as.environment(".GlobalEnv"));
		peaks<-profileList_neg[[7]];
		peaklist<-peaks[,c(14,16,15)];
		treshscore<-as.numeric(logfile$parameters$ISnorm_score_neg)
		we1=(1-as.numeric(logfile$parameters$IS_w1))
		we2=(1-we1)
		cat("- screening:");
		screened<-screening(	peaklist, 
								blanklist=FALSE, 
								pattern_neg_IS, 
								RT = patternRT_neg_IS,
								dmz = as.numeric(logfile$parameters$IS_dmz), 
								ppm = logfile$parameters$IS_ppm, 
								dRT = patternDelRT_neg_IS, 
								dRTwithin = as.numeric(logfile$parameters$IS_drt2), 
								dRTblank = FALSE, 
								dInt = as.numeric(logfile$parameters$IS_inttol), 
								Intcut = as.numeric(logfile$parameters$IS_intcut), 
								w1=we1, 
								w2=we2, 	
								w3=0	
		)
		# set matrix to sort & store data from a profile ###################################
		cat("- on IS profiles")
		atPOSIX<-profileList_neg[[3]];
		sampletype<-profileList_neg[[9]];
		sampleID<-profileList_neg[[4]];
		# filter out other file types such as spiked ones
		keep<-((sampletype=="sample")|(sampletype=="blank"))
		atPOSIX<-atPOSIX[keep]
		sampletype<-sampletype[keep]
		sampleID<-sampleID[keep]
		#
		atdate<-c();
		attime<-c();
		for(i in 1:length(atPOSIX)){
				  atdate<-c(atdate, strsplit(atPOSIX[i]," ")[[1]][1]);
				  attime<-c(attime, strsplit(atPOSIX[i]," ")[[1]][2]);
		}
		attime<-as.difftime(attime);
		atdate<-as.Date(atdate, tz="GMT");
		ord<-order(as.numeric(atdate),as.numeric(attime),sampleID);
		atPOSIXsort<-atPOSIX[ord];
		atdate<-atdate[ord];
		attime<-attime[ord];
		sampleID<-sampleID[ord];
		sampletype<-sampletype[ord];
		# to retrieve data
		timeset<-matrix(nrow=length(atPOSIX),ncol=5,0);
		for(i in 1:length(sampleID)){
			  if(sampletype[i]=="sample"){
				timeset[i,2]<-as.numeric(sampleID[i]);
			  }
			  if(sampletype[i]=="blank"){
				timeset[i,3]<-as.numeric(sampleID[i]);
			  }
		}
		# screen IS intensity profiles #####################################################
		min_count<-floor(length(profileList_neg[[4]])*as.numeric(logfile$parameters$ISnorm_percfiles_neg)/100);
		lis_delint_IS<-list(0)
		lis_median_IS<-list(0)
		stillin<-rep(TRUE,length(peaks[,1]))
		for(p in 1:length(timeset[,1])){
			lis_delint_IS[[p]]<-numeric(0)
			lis_median_IS[[p]]<-numeric(0)
		}
		for(i in 1:length(screened)){
			if(length(screened[[i]])!=1){
				for(j in 1:length(screened[[i]][,1])){
					if(  !(as.character(screened[[i]][j,9]))=="NaN"  ){
						if( as.numeric(as.character(screened[[i]][j,9]))>=treshscore ){
							hits<-strsplit(as.character(screened[[i]][j,2]),"/")[[1]]
							hits<-hits[hits!="-"];
							hits<-as.numeric(hits);
							for(b in 1:length(hits)){
								if( peaks[hits[b],3]>min_count ){
									profID<-as.numeric(peaks[hits[b],4])
									timeset[,4:5]<-0;
									timeset[,c(4,5)] <-.Call("_enviMass_fill_timeset",
														as.numeric(timeset),
														as.numeric(profileList_neg[[2]][(profileList_neg[[7]][profileList_neg[[7]][,4]==profID,1]:profileList_neg[[7]][profileList_neg[[7]][,4]==profID,2]),6]), # sampleIDs
														as.numeric(profileList_neg[[2]][(profileList_neg[[7]][profileList_neg[[7]][,4]==profID,1]:profileList_neg[[7]][profileList_neg[[7]][,4]==profID,2]),2]), # intensities
														as.integer(length(timeset[,1])),
														PACKAGE="enviMass"
									)
									timeset[timeset[,4]!=0,4]<-log10(timeset[timeset[,4]!=0,4])
									timeset[timeset[,5]!=0,5]<-log10(timeset[timeset[,5]!=0,5])
									median_int<-c(timeset[timeset[,4]!=0,4],timeset[timeset[,5]!=0,5])
									median_int<-median(median_int);
									for(k in 1:length(timeset[,1])){
										if(timeset[k,5]!=0){
											lis_delint_IS[[k]]<-c(lis_delint_IS[[k]],(timeset[k,5]-median_int))
											lis_median_IS[[k]]<-c(lis_median_IS[[k]],median_int)									
											next;
										}
										if(timeset[k,4]!=0){
											lis_delint_IS[[k]]<-c(lis_delint_IS[[k]],(timeset[k,4]-median_int))
											lis_median_IS[[k]]<-c(lis_median_IS[[k]],median_int)	
										}
									}
									stillin[hits[b]]<-FALSE									
								}	
							}
						}
					}	
				}
			}
		}
		# screen other profiles ###############################################################
		if( logfile$parameters$ISnorm_medblank_neg=="TRUE" || logfile$parameters$ISnorm_medsam_neg=="TRUE" ){
			cat("- on blank / non-blank profiles");
			peaks<-peaks[stillin,];
			peaks<-peaks[sample(1:length(peaks[,1]),length(peaks[,1]), replace = FALSE),];
			if(logfile$parameters$ISnorm_medblank_neg=="TRUE"){ # use blank
				count_b<-0
				if( logfile$parameters$ISnorm_usesubblank_neg=="TRUE" ){ # use subsampling
					max_count_b<-as.numeric(logfile$parameters$ISnorm_numblank_neg)
				}else{
					max_count_b<-length(peaks[,1])
				}
			}else{
				count_b<-0
				max_count_b<-0
			}
			if(logfile$parameters$ISnorm_medsam_neg=="TRUE"){ # use nonblank
				count_nb<-0
				if( logfile$parameters$ISnorm_usesubsam_neg=="TRUE" ){ # use subsampling
					max_count_nb<-as.numeric(logfile$parameters$ISnorm_numsam_neg)
				}else{
					max_count_nb<-length(peaks[,1])
				}
			}else{
				count_nb<-0
				max_count_nb<-0
			}
			lis_delint_nb<-list(0)
			lis_median_nb<-list(0)
			lis_delint_b<-list(0)
			lis_median_b<-list(0)			
			stillin<-rep(TRUE,length(peaks[,1]))
			for(p in 1:length(timeset[,1])){
				lis_delint_nb[[p]]<-numeric(0)
				lis_median_nb[[p]]<-numeric(0)
				lis_delint_b[[p]]<-numeric(0)
				lis_median_b[[p]]<-numeric(0)
			}
			for(i in 1:length(peaks[,1])){		
				profID<-as.numeric(peaks[i,4])
				timeset[,4:5]<-0;
				timeset[,c(4,5)] <-.Call("_enviMass_fill_timeset",
					as.numeric(timeset),
					as.numeric(profileList_neg[[2]][(profileList_neg[[7]][profileList_neg[[7]][,4]==profID,1]:profileList_neg[[7]][profileList_neg[[7]][,4]==profID,2]),6]), # sampleIDs
					as.numeric(profileList_neg[[2]][(profileList_neg[[7]][profileList_neg[[7]][,4]==profID,1]:profileList_neg[[7]][profileList_neg[[7]][,4]==profID,2]),2]), # intensities
					as.integer(length(timeset[,1])),
					PACKAGE="enviMass"
				)
				if(	 # for blind #############################################3
					(logfile$parameters$ISnorm_medblank_neg=="TRUE") &&
					(any(timeset[,5]>0)) &&
					(count_b<max_count_b)
				){ 
					timeset[timeset[,4]!=0,4]<-log10(timeset[timeset[,4]!=0,4])
					timeset[timeset[,5]!=0,5]<-log10(timeset[timeset[,5]!=0,5])
					median_int<-c(timeset[timeset[,4]!=0,4],timeset[timeset[,5]!=0,5])
					median_int<-median(median_int);
					for(k in 1:length(timeset[,1])){
						if(timeset[k,5]!=0){
							lis_delint_b[[k]]<-c(lis_delint_b[[k]],(timeset[k,5]-median_int))
							lis_median_b[[k]]<-c(lis_median_b[[k]],median_int)									
							next;
						}
						if(timeset[k,4]!=0){
							lis_delint_b[[k]]<-c(lis_delint_b[[k]],(timeset[k,4]-median_int))
							lis_median_b[[k]]<-c(lis_median_b[[k]],median_int)	
						}
					}
					count_b<-(count_b+1);
					next;
				}
				if(	 # for non-blind ##########################################
					(logfile$parameters$ISnorm_medsam_neg=="TRUE") &&
					(any(timeset[,4]>0)) &&
					(count_nb<max_count_nb)
				){ 
					timeset[timeset[,4]!=0,4]<-log10(timeset[timeset[,4]!=0,4])
					median_int<-median(timeset[timeset[,4]!=0,4]);
					for(k in 1:length(timeset[,1])){
						if(timeset[k,4]!=0){
							lis_delint_nb[[k]]<-c(lis_delint_nb[[k]],(timeset[k,4]-median_int))
							lis_median_nb[[k]]<-c(lis_median_nb[[k]],median_int)									
						}
					}
					count_nb<-(count_nb+1);
				}
				if( # enough subsampling?
					(count_b>=max_count_b) &&
					(count_nb>=max_count_nb) 
				){
					break
				}
			}
		}
		# correct intensities - replace in profileList, recalculate mean_int ###############
		cat("- correcting intensities");
		corfac<-c()
		for(k in 1:length(lis_delint_IS)){
			if( length(lis_delint_IS[[k]])>as.numeric(logfile$parameters$ISnorm_numbIS_neg) ){
				corfac<-c( corfac,10^median(lis_delint_IS[[k]]) )
			}else{			
				corfac<-c(corfac,1)
			}
		}
		corr_intens <- .Call("_enviMass_correct_intens",
							as.numeric(corfac),	  # correction factor
							as.integer(sampleID),       
							as.numeric(profileList_neg[[2]][,2]), # intensities
							as.integer(profileList_neg[[2]][,6]),  
							PACKAGE="enviMass"
							)
		profileList_neg[[2]][,2]<-corr_intens
		for(k in 1:length(profileList_neg[[7]][,8])){
			profileList_neg[[7]][k,16]<-mean(profileList_neg[[2]][(profileList_neg[[7]][k,1]:profileList_neg[[7]][k,2]),2])
		}
		save(profileList_neg,file=file.path(as.character(logfile[[1]]),"results","profileList_neg"),compress=FALSE);	
		# derive plots #####################################################################

# > BAUSTELLE

		# derive plots #####################################################################
		int_norm_ISTD_neg <<- list()
		int_norm_ISTD_neg[[1]] <<- lis_delint_IS
		int_norm_ISTD_neg[[2]] <<- lis_median_IS
		#int_norm_ISTD_neg[[3]] <- lis_RT_IS
		int_norm_ISTD_neg[[4]] <<- use_corfac
		int_norm_ISTD_neg[[5]] <<- list()#lis_delint_nb
		int_norm_ISTD_neg[[6]] <<- list()#lis_median_nb
		int_norm_ISTD_neg[[7]] <<- list()#lis_delint_b
		int_norm_ISTD_neg[[8]] <<- list()#lis_median_b
		int_norm_ISTD_neg[[9]] <<- profileList_neg[["datetime"]]
		int_norm_ISTD_neg[[10]] <<- profileList_neg[["type"]]
		int_norm_ISTD_neg[[11]] <<- profileList_neg[["sampleID"]]	
		names(int_norm_ISTD_neg) <<- c("lis_delint_IS", " lis_median_IS", "lis_RT_IS", "use_corfac", "lis_delint_nb", 
			"lis_median_nb", "lis_delint_b", "lis_median_b", "atPOSIX", "sampletype", "sampleID")
		# -> save data & derive plots ######################################################
		output$int_norm_ISTD_neg_median <- renderPlot({   
			par(mar = c(.2, 4.5, .9, 8))
			enviMass:::plot_ISTD_norm(
				int_norm_ISTD = int_norm_ISTD_neg,
				logfile = logfile,
				what = "normalization"
			)
		},res = 100) 
		output$int_norm_ISTD_neg_counts <- renderPlot({   
			par(mar = c(4.5, 4.5, .9, 8))
			enviMass:::plot_ISTD_norm(
				int_norm_ISTD = int_norm_ISTD_neg,
				logfile = logfile,
				what = "counts"
			)
		},res = 100) 
		save(int_norm_ISTD_neg, file = file.path(as.character(logfile[[1]]), "results", "int_norm_ISTD_neg"));
		####################################################################################

# < BAUSTELLE


		####################################################################################

	}



