VGMYK9509S <- function(
logfile,
Ion_mode_profiles,
...
){
if(any(objects(envir = as.environment(".GlobalEnv")) == "profileList_pos")){rm(profileList_pos, envir = as.environment(".GlobalEnv"))}
if(any(objects() == "profileList_pos")){rm(profileList_pos)}
if(any(objects(envir = as.environment(".GlobalEnv")) == "profileList_neg")){rm(profileList_neg, envir = as.environment(".GlobalEnv"))}
if(any(objects() == "profileList_neg")){rm(profileList_neg)}
if(any(objects(envir = as.environment(".GlobalEnv")) == "profileList_pos_cal")){rm(profileList_pos_cal, envir = as.environment(".GlobalEnv"))}
if(any(objects() == "profileList_pos_cal")){rm(profileList_pos_cal)}
if(any(objects(envir = as.environment(".GlobalEnv")) == "profileList_neg_cal")){rm(profileList_neg_cal, envir = as.environment(".GlobalEnv"))}
if(any(objects() == "profileList_neg_cal")){rm(profileList_neg_cal)}
if(any(objects(envir = as.environment(".GlobalEnv")) == "links_profiles_pos")){rm(links_profiles_pos, envir = as.environment(".GlobalEnv"))}
if(any(objects() == "links_profiles_pos")){rm(links_profiles_pos)}
if(any(objects(envir = as.environment(".GlobalEnv")) == "links_profiles_neg")){rm(links_profiles_neg, envir = as.environment(".GlobalEnv"))}
if(any(objects() == "links_profiles_neg")){rm(links_profiles_neg)}
if(
file.exists(file.path(logfile$project_folder,"results","profileList_pos")) &
!any(objects(envir = as.environment(".GlobalEnv")) == "no_load") &
(Ion_mode_profiles =="positive")
){
load(file=file.path(as.character(logfile$project_folder),"results","profileList_pos"),envir=as.environment(".GlobalEnv"), verbose=TRUE);
assign("profileList",profileList_pos,envir=as.environment(".GlobalEnv"));
if(file.exists(file.path(logfile[[1]], "results", "links_profiles_pos"))){
load(file = file.path(as.character(logfile[[1]]),"results","links_profiles_pos"), envir = as.environment(".GlobalEnv"), verbose=TRUE);
assign("links_profiles", links_profiles_pos, envir = as.environment(".GlobalEnv"));
}
}
if(
file.exists(file.path(logfile$project_folder,"quantification","profileList_pos_cal")) &
!any(objects(envir = as.environment(".GlobalEnv")) == "no_load")
){
load(file=file.path(as.character(logfile$project_folder),"quantification","profileList_pos_cal"),envir=as.environment(".GlobalEnv"), verbose=TRUE);
}
if(
file.exists(file.path(logfile$project_folder,"results","profileList_neg")) &
!any(objects(envir = as.environment(".GlobalEnv")) == "no_load") &
(Ion_mode_profiles == "negative")
){
load(file=file.path(as.character(logfile$project_folder),"results","profileList_neg"),envir=as.environment(".GlobalEnv"), verbose=TRUE);
assign("profileList",profileList_neg,envir=as.environment(".GlobalEnv"));
if(file.exists(file.path(logfile[[1]], "results", "links_profiles_neg"))){
load(file = file.path(as.character(logfile[[1]]),"results","links_profiles_neg"), envir = as.environment(".GlobalEnv"), verbose=TRUE);
assign("links_profiles", links_profiles_neg, envir = as.environment(".GlobalEnv"));
}
}
if(
file.exists(file.path(logfile$project_folder,"quantification","profileList_neg_cal")) &
!any(objects(envir = as.environment(".GlobalEnv")) == "no_load")
){
load(file=file.path(as.character(logfile$project_folder),"quantification","profileList_neg_cal"),envir=as.environment(".GlobalEnv"), verbose=TRUE);
}
return("Workflow R objects updated")
}
