#######################################################################################################
# POSITIVE IONIZATION #################################################################################
if(any(objects(envir=as.environment(".GlobalEnv"))=="profileList_pos")){rm(profileList_pos,envir=as.environment(".GlobalEnv"))}
if(any(objects()=="profileList_pos")){rm(profileList_pos)}
if(any(objects(envir=as.environment(".GlobalEnv"))=="links_peaks_pos")){rm(links_peaks_pos,envir=as.environment(".GlobalEnv"))}
if(any(objects()=="links_peaks_pos")){rm(links_peaks_pos)}

if(
	file.exists(file.path(as.character(logfile[["project_folder"]]),"results","profileList_pos")) 
){

	load(file=file.path(as.character(logfile[["project_folder"]]),"results","profileList_pos"),envir=as.environment(".GlobalEnv"));	
	###################################################################################################
	# blind peak subtraction ##########################################################################
	if(
		(logfile$parameters$subtr_blind=="yes") &	
		(any(profileList_pos[["peaks"]][,"in_blind"]!=Inf)) # if any blind peaks exist
	){
		keep_peaks<-rep(TRUE,dim(profileList_pos[["peaks"]])[1])
		# sample peaks found in blind
		keep_peaks[profileList_pos[["peaks"]][,"in_blind"]<as.numeric(logfile$parameters$blind_threshold)]<-FALSE
		# blind peaks themselves
		these<-as.numeric(profileList_pos[["sampleID"]][profileList_pos[["type"]]=="blank"])
		#system.time({# replace by match()?		
		if(length(these)>0){
			for(i in 1:length(these)){
				keep_peaks[profileList_pos[["peaks"]][,"sampleIDs"]==these[i]]<-FALSE
			}
		}		
		#})		
		profileList_pos[["peaks"]]<<-profileList_pos[["peaks"]][keep_peaks,]
		cat(paste("\nBlind in profile subtraction, positive: ",round((sum(!keep_peaks)/length(keep_peaks)*100),digits=3),"% of peaks removed\n",sep=""))
	}else{
		cat("\nblind peak subtraction nothing to subtract, positive ioniz.")
	}
	###################################################################################################
	# spike file peak subtraction #####################################################################
	if( logfile$parameters$subtr_spiked=="yes" ){
		these<-as.numeric(profileList_pos[["sampleID"]][profileList_pos[["type"]]=="spiked"])
		if(length(these)>0){
			keep_peaks<-rep(TRUE,dim(profileList_pos[["peaks"]])[1])
			for(i in 1:length(these)){
				keep_peaks[profileList_pos[["peaks"]][,"sampleIDs"]==these[i]]<-FALSE
			}
			profileList_pos[["peaks"]]<<-profileList_pos[["peaks"]][keep_peaks,]
		}else{
			cat("\nNo peaks from files to remove peaks for, positive.")
		}	
	}	
	###################################################################################################
	# IS peak subtraction #############################################################################
	if(
		(logfile$parameters$subtr_IS=="yes") &
		(file.exists(file.path(as.character(logfile[[1]]),"results","links_peaks_pos")))
	){
		load(file=file.path(as.character(logfile[[1]]),"results","links_peaks_pos"));
		if(length(links_peaks_pos)>0){
			keep_peaks<-rep(TRUE,dim(profileList_pos[["peaks"]])[1])
			for(i in 1:length(profileList_pos[["peaks"]][,1])){
				if(profileList_pos[["peaks"]][i,"links"]!=0){
					if(length(links_peaks_pos[[profileList_pos[["peaks"]][i,"links"]]][[2]])>0){
						keep_peaks[i]<-FALSE
						links_peaks_pos[[profileList_pos[["peaks"]][i,"links"]]][[2]]<-list() # = IS entry[[2]] in links_peaks_pos
					}
				}
			}
			profileList_pos[["peaks"]]<<-profileList_pos[["peaks"]][keep_peaks,]
			save(links_peaks_pos,file=file.path(as.character(logfile[[1]]),"results","links_peaks_pos"));
			rm(links_peaks_pos)		
			cat(paste("\nIS subtraction: ",round((sum(!keep_peaks)/length(keep_peaks)*100),digits=3),"% of peaks removed\n",sep=""))
		}else{
			cat("\n IS peak subtraction: nothing to subtract, positive ioniz.")
		}
	}
	###################################################################################################
	# target peak subtraction #########################################################################
	if(
		(logfile$parameters$subtr_target=="yes")  &
		(file.exists(file.path(as.character(logfile[[1]]),"results","links_peaks_pos")))
	){
		load(file=file.path(as.character(logfile[[1]]),"results","links_peaks_pos"));
		if(length(links_peaks_pos)>0){
			keep_peaks<-rep(TRUE,dim(profileList_pos[["peaks"]])[1])
			for(i in 1:length(profileList_pos[["peaks"]][,1])){
				if(profileList_pos[["peaks"]][i,"links"]!=0){
					if(length(links_peaks_pos[[profileList_pos[["peaks"]][i,"links"]]][[1]])>0){
						keep_peaks[i]<-FALSE
						links_peaks_pos[[profileList_pos[["peaks"]][i,"links"]]][[1]]<-list() # = target entry[[1]] in links_peaks_pos
					}
				}
			}
			profileList_pos[["peaks"]]<<-profileList_pos[["peaks"]][keep_peaks,]
			save(links_peaks_pos,file=file.path(as.character(logfile[[1]]),"results","links_peaks_pos"));	
			rm(links_peaks_pos)	
			cat(paste("\nTarget subtraction: ",round((sum(!keep_peaks)/length(keep_peaks)*100),digits=3),"% of peaks removed\n",sep=""))
		}else{
			cat("\n Target peak subtraction: nothing to subtract, positive ioniz.")
		}
	}
	###################################################################################################
	profileList_pos[[2]][,"profileIDs"]<<-0
	profileList_pos[[6]]<<-0  # index_agglom
	profileList_pos[[7]]<<-0  # index_prof
	profileList_pos[[8]]<<-0  # parameters
	profileList_pos[[1]]<<-data.frame(TRUE,FALSE,FALSE,FALSE)    # state
	colnames(profileList_pos[[1]])<<-c("peaks?","agglom?","profiling","trends?")
	###################################################################################################
	profileList_pos<<-agglomer(
						profileList_pos,
						dmass=(as.numeric(logfile$parameters$prof_dmz)+1),
						ppm=as.logical(as.character(logfile$parameters$prof_ppm)),
						dret=(as.numeric(logfile$parameters$prof_drt)+10)
					)
	if(logfile$parameters$replicates_prof=="no"){ 				
		profileList_pos<<-partcluster(
							profileList=profileList_pos,
							dmass=as.numeric(logfile$parameters$prof_dmz),
							ppm=as.logical(as.character(logfile$parameters$prof_ppm)),
							dret=as.numeric(logfile$parameters$prof_drt),
							from=FALSE,
							to=FALSE,
							progbar=logfile$parameters$progressBar,
							plot_it=FALSE,
							replicates=FALSE
						)
	}else{ # run a profiling in the replicate groups first
		measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
		measurements<-measurements[measurements[,"include"]=="TRUE",]
		measurements<-measurements[measurements[,"Mode"]=="positive",]
		replicates<-measurements$tag3
		IDs<-measurements$ID
		profileList_pos<<-partcluster(
							profileList=profileList_pos,
							dmass=as.numeric(logfile$parameters$prof_dmz),
							ppm=as.logical(as.character(logfile$parameters$prof_ppm)),
							dret=as.numeric(logfile$parameters$prof_drt),
							from=FALSE,
							to=FALSE,
							progbar=logfile$parameters$progressBar,
							plot_it=FALSE,
							replicates=replicates,
							IDs
						)
	}
	###################################################################################################
# > BAUSTELLE: clean links_peaks_pos for removed links!

# <
	profileList_pos<<-enviMass:::in_blind(profileList_pos)
	save(profileList_pos,file=file.path(as.character(logfile[["project_folder"]]),"results","profileList_pos"),compress=FALSE);
	profpeaks_pos<<-enviMass:::profiletopeak(profileList_pos,progbar=logfile$parameters$progressBar)		
	profpeaks_pos<<-profpeaks_pos[order(profpeaks_pos[,13],decreasing=TRUE),];
	save(profpeaks_pos,file=file.path(as.character(logfile[["project_folder"]]),"results","profpeaks_pos"));	
	###################################################################################################

}



#######################################################################################################
# NEGATIVE IONIZATION #################################################################################
if(any(objects(envir=as.environment(".GlobalEnv"))=="profileList_neg")){rm(profileList_neg,envir=as.environment(".GlobalEnv"))}
if(any(objects()=="profileList_neg")){rm(profileList_neg)}
if(any(objects(envir=as.environment(".GlobalEnv"))=="links_peaks_neg")){rm(links_peaks_neg,envir=as.environment(".GlobalEnv"))}
if(any(objects()=="links_peaks_neg")){rm(links_peaks_neg)}

if(
	file.exists(file.path(as.character(logfile[["project_folder"]]),"results","profileList_neg")) 
){

	load(file=file.path(as.character(logfile[["project_folder"]]),"results","profileList_neg"),envir=as.environment(".GlobalEnv"));	
	###################################################################################################
	# blind peak subtraction ##########################################################################
	if(
		(logfile$parameters$subtr_blind=="yes") &	
		(any(profileList_neg[["peaks"]][,9]!=Inf)) # if any blind peaks exist
	){
		keep_peaks<-rep(TRUE,dim(profileList_neg[["peaks"]])[1])
		keep_peaks[profileList_pos[["peaks"]][,"in_blind"]<as.numeric(logfile$parameters$blind_threshold)]<-FALSE
		# blind peaks themselves
		these<-as.numeric(profileList_neg[["sampleID"]][profileList_neg[["type"]]=="blank"])
		if(length(these)>0){
			for(i in 1:length(these)){
				keep_peaks[profileList_neg[["peaks"]][,6]==these[i]]<-FALSE
			}
		}		
		profileList_neg[["peaks"]]<<-profileList_neg[["peaks"]][keep_peaks,]
		cat(paste("\nBlind in profile subtraction, negative: ",round((sum(!keep_peaks)/length(keep_peaks)*100),digits=3),"% of peaks removed\n",sep=""))
	}else{
		cat("\nblind peak subtraction: nothing to subtract, negative ioniz.")
	}
	###################################################################################################
	# spike file peak subtraction #####################################################################
	if( logfile$parameters$subtr_spiked=="yes" ){
		these<-as.numeric(profileList_neg[["sampleID"]][profileList_neg[["type"]]=="spiked"])
		if(length(these)>0){
			keep_peaks<-rep(TRUE,dim(profileList_neg[["peaks"]])[1])
			for(i in 1:length(these)){
				keep_peaks[profileList_neg[["peaks"]][,6]==these[i]]<-FALSE
			}
			profileList_neg[["peaks"]]<<-profileList_neg[["peaks"]][keep_peaks,]
		}else{
			cat("\nNo peaks from files to remove peaks for, negative.")
		}	
	}	
	###################################################################################################
	# IS peak subtraction #############################################################################
	if(
		(logfile$parameters$subtr_IS=="yes") &
		(file.exists(file.path(as.character(logfile[["project_folder"]]),"results","links_peaks_neg")))
	){
		load(file=file.path(as.character(logfile[["project_folder"]]),"results","links_peaks_neg"));
		if(length(links_peaks_neg)>0){
			keep_peaks<-rep(TRUE,dim(profileList_neg[["peaks"]])[1])
			for(i in 1:length(profileList_neg[["peaks"]][,1])){
				if(profileList_neg[["peaks"]][i,5]!=0){
					if(length(links_peaks_neg[[profileList_neg[["peaks"]][i,5]]][[2]])>0){
						keep_peaks[i]<-FALSE
						links_peaks_neg[[profileList_neg[["peaks"]][i,5]]][[2]]<-list() # = IS entry[[2]] in links_peaks_neg
					}
				}
			}
			profileList_neg[["peaks"]]<<-profileList_neg[["peaks"]][keep_peaks,]
			save(links_peaks_neg,file=file.path(as.character(logfile[["project_folder"]]),"results","links_peaks_neg"));
			rm(links_peaks_neg)		
			cat(paste("\nIS subtraction: ",round((sum(!keep_peaks)/length(keep_peaks)*100),digits=3),"% of peaks removed\n",sep=""))
		}else{
			cat("\n IS peak subtraction: nothing to subtract, negative ioniz.")
		}
	}
	###################################################################################################
	# target peak subtraction #########################################################################
	if(
		(logfile$parameters$subtr_target=="yes")  &
		(file.exists(file.path(as.character(logfile[["project_folder"]]),"results","links_peaks_neg")))
	){
		load(file=file.path(as.character(logfile[["project_folder"]]),"results","links_peaks_neg"));
		if(length(links_peaks_neg)>0){
			keep_peaks<-rep(TRUE,dim(profileList_neg[["peaks"]])[1])
			for(i in 1:length(profileList_neg[["peaks"]][,1])){
				if(profileList_neg[["peaks"]][i,5]!=0){
					if(length(links_peaks_neg[[profileList_neg[["peaks"]][i,5]]][[1]])>0){
						keep_peaks[i]<-FALSE
						links_peaks_neg[[profileList_neg[["peaks"]][i,5]]][[1]]<-list() # = target entry[[1]] in links_peaks_neg
					}
				}
			}
			profileList_neg[["peaks"]]<<-profileList_neg[["peaks"]][keep_peaks,]
			save(links_peaks_neg,file=file.path(as.character(logfile[[1]]),"results","links_peaks_neg"));	
			rm(links_peaks_neg)	
			cat(paste("\nTarget subtraction: ",round((sum(!keep_peaks)/length(keep_peaks)*100),digits=3),"% of peaks removed\n",sep=""))
		}else{
			cat("\n Target peak subtraction: nothing to subtract, negative ioniz.")
		}
	}
	###################################################################################################
	profileList_neg[[2]][,"profileIDs"]<<-0
	profileList_neg[[6]]<<-0  # index_agglom
	profileList_neg[[7]]<<-0  # index_prof
	profileList_neg[[8]]<<-0  # parameters
	profileList_neg[[1]]<<-data.frame(TRUE,FALSE,FALSE,FALSE)    # state
	colnames(profileList_neg[[1]])<<-c("peaks?","agglom?","profiling","trends?")
	###################################################################################################
	profileList_neg<<-agglomer(
						profileList_neg,
						dmass=(as.numeric(logfile$parameters$prof_dmz)+1),
						ppm=as.logical(as.character(logfile$parameters$prof_ppm)),
						dret=(as.numeric(logfile$parameters$prof_drt)+10)
					)
	if(logfile$parameters$replicates_prof=="no"){ 				
		profileList_neg<<-partcluster(
							profileList=profileList_neg,
							dmass=as.numeric(logfile$parameters$prof_dmz),
							ppm=as.logical(as.character(logfile$parameters$prof_ppm)),
							dret=as.numeric(logfile$parameters$prof_drt),
							from=FALSE,
							to=FALSE,
							progbar=logfile$parameters$progressBar,
							plot_it=FALSE,
							replicates=FALSE
						)
	}else{ # run a profiling in the replicate groups first
		measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
		measurements<-measurements[measurements[,"include"]=="TRUE",]
		measurements<-measurements[measurements[,"Mode"]=="negative",]
		replicates<-measurements$tag3
		IDs<-measurements$ID
		profileList_neg<<-partcluster(
							profileList=profileList_neg,
							dmass=as.numeric(logfile$parameters$prof_dmz),
							ppm=as.logical(as.character(logfile$parameters$prof_ppm)),
							dret=as.numeric(logfile$parameters$prof_drt),
							from=FALSE,
							to=FALSE,
							progbar=logfile$parameters$progressBar,
							plot_it=FALSE,
							replicates=replicates,
							IDs
						)
	}
	###################################################################################################
# > BAUSTELLE: clean links_peaks_neg for removed peak links!

# <
	profileList_neg<<-enviMass:::in_blind(profileList_neg)
	save(profileList_neg,file=file.path(as.character(logfile[["project_folder"]]),"results","profileList_neg"),compress=FALSE);
	profpeaks_neg<<-enviMass:::profiletopeak(profileList_neg,progbar=logfile$parameters$progressBar)		
	profpeaks_neg<<-profpeaks_neg[order(profpeaks_neg[,13],decreasing=TRUE),];
	save(profpeaks_neg,file=file.path(as.character(logfile[["project_folder"]]),"results","profpeaks_neg"));	
	###################################################################################################

}

if(any(ls()=="logfile")){stop("\n illegal logfile detected in do_subtr.r!")}
########################################################################################################