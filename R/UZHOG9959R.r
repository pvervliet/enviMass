UZHOG9959R <- function(
profileList,
prof_ID,
links_profiles,
min_peaks=3,
skip_peaks=FALSE,
min_cor=.9,
with_test=FALSE,
omit_profiles=FALSE
){
if(omit_profiles[1]!="FALSE"){
if(length(omit_profiles)!=length(profileList[["index_prof"]][,"links"])){stop("omit_profiles not equal to number of profiles")}
}
if(with_test){if(prof_ID>length(profileList[["index_prof"]][,"profile_ID"])){stop("\n Debug get_isotopol _1!")}}
if(with_test){if(profileList[["index_prof"]][prof_ID,"profile_ID"]!=prof_ID){stop("\n Debug get_isotopol _2!")}}
if(profileList[["index_prof"]][prof_ID,"links"]==0){
return(c())
}
in_link<-profileList[["index_prof"]][prof_ID,"links"][[1]]
if(length(links_profiles[[in_link]]$adduc[,1])>0){
if(skip_peaks){
those<-which(
((links_profiles[[in_link]]$adduc[,"correl"]/1000)>=min_cor) &
(links_profiles[[in_link]]$adduc[,"ref_1"]>=min_peaks)
)
if(length(those)>0){
those<-links_profiles[[in_link]]$adduc[those,"linked profile"]
if(omit_profiles[1]!="FALSE"){
those<-those[omit_profiles[those]==0]
}
}
}else{
those<-which(
((links_profiles[[in_link]]$adduc[,"correl"]/1000)>=min_cor) |
(links_profiles[[in_link]]$adduc[,"ref_1"]<min_peaks)
)
if(length(those)>0){
those<-links_profiles[[in_link]]$adduc[those,"linked profile"]
if(omit_profiles[1]!="FALSE"){
those<-those[omit_profiles[those]==0]
}
}
}
if(length(those)==0){
prof_adduc_IDs<-c()
}else{
prof_adduc_IDs<-those
}
}else{
return(c())
}
if(with_test){
if(any(profileList[["index_prof"]][prof_adduc_IDs,"profile_ID"]!=prof_adduc_IDs)){
stop("\n Debug get_isotopol _3!")
}
}
return(prof_adduc_IDs)
}
