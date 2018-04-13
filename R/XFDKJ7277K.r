XFDKJ7277K <- function(
links_profiles,
profileList,
logfile,
min_rat=.7,
min_count=.4,
perc =.9,
for_which="ISTD"
){
if(logfile$parameters$dofile_latest_profcomp == "TRUE"){
num_sam <- as.numeric(logfile$parameters$numfile_latest_profcomp) * min_count
}else{
num_sam <- floor(length(profileList[["sampleID"]]) * min_count)
}
if(for_which == "ISTD"){
for_profs <- which(unlist(lapply(lapply(links_profiles, `[[`, 2),length))!=0)
}
if(for_which == "target"){
for_profs <- which(unlist(lapply(lapply(links_profiles, `[[`, 1),length))!=0)
}
if(for_which == "both"){
for_profs1 <- which(unlist(lapply(lapply(links_profiles, `[[`, 1),length))!=0)
for_profs2 <- which(unlist(lapply(lapply(links_profiles, `[[`, 2),length))!=0)
for_profs <- unique(c(for_profs1,for_profs2))
}
if(for_which == "all"){
for_profs <- (1:length(links_profiles))
}
delRT_isot <- c()
delRT_adduc <- c()
int_cor_isot <- c()
int_cor_adduc <- c()
for(n in 1:length(for_profs)){
if(links_profiles[[for_profs[n]]]$total < num_sam) next
prof1 <- as.numeric(names(links_profiles)[for_profs[n]])
if(dim(links_profiles[[for_profs[n]]]$isot)[1] > 0){
for(m in 1:dim(links_profiles[[for_profs[n]]]$isot)[1]){
if((links_profiles[[for_profs[n]]]$isot[m,"link counts"]/links_profiles[[for_profs[n]]]$isot[m,"ref_1"]) < min_rat) next;
if(links_profiles[[for_profs[n]]]$isot[m,"ref_1"] < num_sam) next;
if(links_profiles[[for_profs[n]]]$isot[m,"use"] == 0) next;
prof2<-links_profiles[[for_profs[n]]]$isot[m,"linked profile"]
these<-profileList[["peaks"]][
profileList[["index_prof"]][prof1,"start_ID"]:profileList[["index_prof"]][prof1,"end_ID"]
,"sampleIDs"]
those<-profileList[["peaks"]][
profileList[["index_prof"]][prof2,"start_ID"]:profileList[["index_prof"]][prof2,"end_ID"]
,"sampleIDs"]
matched<-match(these,those)
RT_1<-(profileList[["peaks"]][
(profileList[["index_prof"]][prof1,"start_ID"]:profileList[["index_prof"]][prof1,"end_ID"])
[!is.na(matched)]
,"RT"])
RT_2<-(profileList[["peaks"]][
(profileList[["index_prof"]][prof2,"start_ID"]:profileList[["index_prof"]][prof2,"end_ID"])
[matched[!is.na(matched)]]
,"RT"])
dRT<-abs(RT_1-RT_2)
delRT_isot<-c(delRT_isot,quantile(dRT, perc)[[1]])
int_1<-(profileList[["peaks"]][
(profileList[["index_prof"]][prof1,"start_ID"]:profileList[["index_prof"]][prof1,"end_ID"])
[!is.na(matched)]
,"intensity"])
int_2<-(profileList[["peaks"]][
(profileList[["index_prof"]][prof2,"start_ID"]:profileList[["index_prof"]][prof2,"end_ID"])
[matched[!is.na(matched)]]
,"intensity"])
int_cor_isot<-c(int_cor_isot,cor(int_1,int_2))
}
}
if(dim(links_profiles[[for_profs[n]]]$adduc)[1]>0){
for(m in 1:dim(links_profiles[[for_profs[n]]]$adduc)[1]){
if((links_profiles[[for_profs[n]]]$adduc[m,"link counts"]/links_profiles[[for_profs[n]]]$adduc[m,"ref_1"])<min_rat) next;
if(links_profiles[[for_profs[n]]]$adduc[m,"ref_1"]<num_sam) next;
if(links_profiles[[for_profs[n]]]$adduc[m,"use"]==0) next;
prof2<-links_profiles[[for_profs[n]]]$adduc[m,"linked profile"]
these<-profileList[["peaks"]][
profileList[["index_prof"]][prof1,"start_ID"]:profileList[["index_prof"]][prof1,"end_ID"]
,"sampleIDs"]
those<-profileList[["peaks"]][
profileList[["index_prof"]][prof2,"start_ID"]:profileList[["index_prof"]][prof2,"end_ID"]
,"sampleIDs"]
matched<-match(these,those)
RT_1<-(profileList[["peaks"]][
(profileList[["index_prof"]][prof1,"start_ID"]:profileList[["index_prof"]][prof1,"end_ID"])
[!is.na(matched)]
,"RT"])
RT_2<-(profileList[["peaks"]][
(profileList[["index_prof"]][prof2,"start_ID"]:profileList[["index_prof"]][prof2,"end_ID"])
[matched[!is.na(matched)]]
,"RT"])
dRT<-abs(RT_1-RT_2)
delRT_adduc<-c(delRT_adduc,quantile(dRT, perc)[[1]])
int_1<-(profileList[["peaks"]][
(profileList[["index_prof"]][prof1,"start_ID"]:profileList[["index_prof"]][prof1,"end_ID"])
[!is.na(matched)]
,"intensity"])
int_2<-(profileList[["peaks"]][
(profileList[["index_prof"]][prof2,"start_ID"]:profileList[["index_prof"]][prof2,"end_ID"])
[matched[!is.na(matched)]]
,"intensity"])
int_cor_adduc<-c(int_cor_adduc,cor(int_1,int_2))
}
}
}
listed<-list()
if(length(delRT_isot)) listed[[1]]<-delRT_isot else listed[[1]]<-NA
if(length(delRT_adduc)) listed[[2]]<-delRT_adduc else listed[[2]]<-NA
if(length(int_cor_isot)) listed[[3]]<-int_cor_isot else listed[[3]]<-NA
if(length(int_cor_adduc)) listed[[4]]<-int_cor_adduc else listed[[4]]<-NA
names(listed)<-c("delRT_isot","delRT_adduc","int_cor_isot","int_cor_adduc")
return(listed)
}
