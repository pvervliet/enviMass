##########
if(file.exists(file=file.path(as.character(logfile[[1]]),"results","profpeaks_pos"))){
	file.remove(file=file.path(as.character(logfile[[1]]),"results","profpeaks_pos"))
}
if(file.exists(file=file.path(as.character(logfile[[1]]),"results","profpeaks_neg"))){
	file.remove(file=file.path(as.character(logfile[[1]]),"results","profpeaks_neg"))
}
if(file.exists(file=file.path(as.character(logfile[[1]]),"results","profileList_pos"))){
	file.remove(file=file.path(as.character(logfile[[1]]),"results","profileList_pos"))
}  
if(file.exists(file=file.path(as.character(logfile[[1]]),"results","profileList_neg"))){
	file.remove(file=file.path(as.character(logfile[[1]]),"results","profileList_neg"))
}  
if(file.exists(file=file.path(as.character(logfile[[1]]),"results","links_peaks_pos"))){
	file.remove(file=file.path(as.character(logfile[[1]]),"results","links_peaks_pos"))
}
if(file.exists(file=file.path(as.character(logfile[[1]]),"results","links_peaks_neg"))){
	file.remove(file=file.path(as.character(logfile[[1]]),"results","links_peaks_neg"))
}
if(file.exists(file=file.path(as.character(logfile[[1]]),"results","links_profiles _pos"))){
	file.remove(file=file.path(as.character(logfile[[1]]),"results","links_profiles_pos"))
}
if(file.exists(file=file.path(as.character(logfile[[1]]),"results","links_profiles_neg"))){
	file.remove(file=file.path(as.character(logfile[[1]]),"results","links_profiles_neg"))
}
##########
if(any(objects(envir=as.environment(".GlobalEnv"))=="profpeaks2")){
	rm(profpeaks2,envir=as.environment(".GlobalEnv"))
}
if(any(objects(envir=as.environment(".GlobalEnv"))=="profpeaks3")){
	rm(profpeaks3,envir=as.environment(".GlobalEnv"))
}
if(any(objects(envir=as.environment(".GlobalEnv"))=="profpeaks_pos")){
	rm(profpeaks_pos,envir=as.environment(".GlobalEnv"))
}
if(any(objects(envir=as.environment(".GlobalEnv"))=="profpeaks_neg")){
	rm(profpeaks_neg,envir=as.environment(".GlobalEnv"))
}
if(any(objects(envir=as.environment(".GlobalEnv"))=="profileList_pos")){
	rm(profileList_pos,envir=as.environment(".GlobalEnv"))
} 
if(any(objects(envir=as.environment(".GlobalEnv"))=="profileList_neg")){
	rm(profileList_neg,envir=as.environment(".GlobalEnv"))
} 
if(any(objects(envir=as.environment(".GlobalEnv"))=="profileList")){
	rm(profileList,envir=as.environment(".GlobalEnv"))
} 
##########
if(any(objects()=="profpeaks2")){
	rm(profpeaks2)
}
if(any(objects()=="profpeaks3")){
	rm(profpeaks3)
}
if(any(objects()=="profpeaks_pos")){
	rm(profpeaks_pos)
}
if(any(objects()=="profpeaks_neg")){
	rm(profpeaks_neg)
}
if(any(objects()=="profileList_pos")){
	rm(profileList_pos)
} 
if(any(objects()=="profileList_neg")){
	rm(profileList_neg)
} 
if(any(objects()=="profileList")){
	rm(profileList)
} 
##########