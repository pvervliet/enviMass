# insert value for each profile: median above blind intensity ratio


##################################################################################
# POSITIVE #######################################################################
if( 
	file.exists(file.path(logfile[[1]],"results","profileList_pos")) 
){

	##############################################################################
	if(any(objects(envir=as.environment(".GlobalEnv"))=="profileList_pos")){rm(profileList_pos,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="profileList_pos")){rm(profileList_pos)}	
	load(file.path(as.character(logfile[[1]]),"results","profileList_pos"),envir=as.environment(".GlobalEnv"));
	##############################################################################
	profileList_pos[["index_prof"]][,"above_blind?"]<-Inf
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
	if(any(objects(envir=as.environment(".GlobalEnv"))=="profileList_neg")){rm(profileList_neg,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="profileList_neg")){rm(profileList_neg)}	
	load(file.path(as.character(logfile[[1]]),"results","profileList_neg"),envir=as.environment(".GlobalEnv"));
	##############################################################################
	profileList_neg[["index_prof"]][,"above_blind?"]<-Inf
	##############################################################################
	save(profileList_neg,file=file.path(as.character(logfile[[1]]),"results","profileList_neg"));
	rm(profileList_neg,envir=as.environment(".GlobalEnv"))
	##############################################################################
	
}
##################################################################################


