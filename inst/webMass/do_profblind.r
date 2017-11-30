# insert value for each profile: median above blind intensity ratio


##################################################################################
# POSITIVE #######################################################################
if( 
	file.exists(file.path(logfile[[1]],"results","profileList_pos")) 
){

	##############################################################################
	if(any(objects(envir=as.environment(".GlobalEnv")) == "profileList_pos")){rm(profileList_pos, envir = as.environment(".GlobalEnv"))}
	if(any(objects()=="profileList_pos")){rm(profileList_pos)}	
	load(file.path(as.character(logfile[[1]]),"results","profileList_pos"), envir = as.environment(".GlobalEnv"));
	##############################################################################
	len <- dim(profileList_pos[["index_prof"]])[1]
	profileList_pos[["index_prof"]][,"above_blind?"] <- Inf
	for(n in 1:len){
		these <- profileList_pos[["peaks"]][
				profileList_pos[["index_prof"]][n,"start_ID"]:profileList_pos[["index_prof"]][n,"end_ID"]
			,"in_blind"]	
		if(any(these != Inf)){
			#a<-median(these[these!=Inf])
			a <- mean(these[these != Inf])
			profileList_pos[["index_prof"]][n,"above_blind?"] <- a
		}else{
			profileList_pos[["index_prof"]][n,"above_blind?"] <- Inf
		}
	}
	##############################################################################
	if(as.logical(logfile$parameters$test)){
		##########################################################################
		# profile IDs correct? ###################################################
		for(i in 1:dim(profileList_pos[["index_prof"]])[1]){
			if(
				!all(profileList_pos[["peaks"]][	
					profileList_pos[["index_prof"]][i, "start_ID"]:profileList_pos[["index_prof"]][i, "end_ID"]
				,"profileIDs"] == i)
			){
				stop("\n Debug do_IS_normaliz.r at #2")
			}
		}
		##########################################################################
	}
	##############################################################################
	save(profileList_pos,file=file.path(as.character(logfile[[1]]),"results","profileList_pos"));
	rm(profileList_pos,envir=as.environment(".GlobalEnv"))
	##############################################################################
	
}
##################################################################################

##################################################################################
# NEGATIVE #######################################################################
if( 
	file.exists(file.path(logfile[[1]],"results","profileList_neg")) 
){

	##############################################################################
	if(any(objects(envir=as.environment(".GlobalEnv")) == "profileList_neg")){rm(profileList_pos, envir = as.environment(".GlobalEnv"))}
	if(any(objects() == "profileList_neg")){rm(profileList_pos)}	
	load(file.path(as.character(logfile[[1]]),"results","profileList_neg"), envir = as.environment(".GlobalEnv"));
	##############################################################################
	len <- dim(profileList_neg[["index_prof"]])[1]
	profileList_neg[["index_prof"]][,"above_blind?"] <- Inf
	for(n in 1:len){
		these<-profileList_neg[["peaks"]][
				profileList_neg[["index_prof"]][n,"start_ID"]:profileList_neg[["index_prof"]][n,"end_ID"]
			,"in_blind"]	
		if(any(these != Inf)){
			#a<-median(these[these!=Inf])
			a <- mean(these[these != Inf])			
			profileList_neg[["index_prof"]][n,"above_blind?"] <- a
		}else{
			profileList_neg[["index_prof"]][n,"above_blind?"] <- Inf
		}
	}
	##############################################################################
	if(as.logical(logfile$parameters$test)){
		##########################################################################
		# profile IDs correct? ###################################################
		for(i in 1:dim(profileList_neg[["index_prof"]])[1]){
			if(
				!all(profileList_neg[["peaks"]][	
					profileList_neg[["index_prof"]][i, "start_ID"]:profileList_neg[["index_prof"]][i, "end_ID"]
				,"profileIDs"] == i)
			){
				stop("\n Debug do_IS_normaliz.r at #2")
			}
		}
		##########################################################################
	}
	##############################################################################
	save(profileList_neg, file = file.path(as.character(logfile[[1]]),"results","profileList_neg"));
	rm(profileList_neg, envir = as.environment(".GlobalEnv"))
	##############################################################################
	
}
##################################################################################


