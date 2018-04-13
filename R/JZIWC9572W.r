JZIWC9572W <- function(
links_profiles,
profileList,
cut_delRT_adduc = 5,
cut_frac_adduc = .9
){
for_profs<-(1:length(links_profiles))
found<-0
removed<-0
for(n in 1:length(for_profs)){
prof1<-as.numeric(names(links_profiles)[for_profs[n]])
if(dim(links_profiles[[for_profs[n]]]$adduc)[1]>0){
keep<-rep(TRUE,dim(links_profiles[[for_profs[n]]]$adduc)[1])
found<-(found+length(keep))
for(m in 1:dim(links_profiles[[for_profs[n]]]$adduc)[1]){
if(links_profiles[[for_profs[n]]]$adduc[m,"use"]==0) next;
prof2<-links_profiles[[for_profs[n]]]$adduc[m,"linked profile"]
num_peaks<-links_profiles[[for_profs[n]]]$adduc[m,"ref_1"]
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
got_frac<-(sum(dRT<cut_delRT_adduc)/length(dRT))
if(got_frac<cut_frac_adduc){
keep[m]<-FALSE
next;
}
if(num_peaks>2){
int_1<-(profileList[["peaks"]][
(profileList[["index_prof"]][prof1,"start_ID"]:profileList[["index_prof"]][prof1,"end_ID"])
[!is.na(matched)]
,"intensity"])
int_2<-(profileList[["peaks"]][
(profileList[["index_prof"]][prof2,"start_ID"]:profileList[["index_prof"]][prof2,"end_ID"])
[matched[!is.na(matched)]]
,"intensity"])
int_cor<-cor(int_1,int_2)
if(!is.na(int_cor)){
int_cor<-as.integer(round(int_cor*1000))
links_profiles[[for_profs[n]]]$adduc[m,"correl"]<-int_cor
at_entry_2<-profileList[["index_prof"]][prof2,"links"]
at<-which(links_profiles[[at_entry_2]]$adduc[,"linked profile"]==prof1)
links_profiles[[at_entry_2]]$adduc[at,"correl"]<-int_cor
}
}
}
if(any(!keep)){
removed<-(removed+sum(keep==FALSE))
those<-links_profiles[[for_profs[n]]]$adduc[!keep,"linked profile"]
for(prof2 in those){
at_entry_2<-profileList[["index_prof"]][prof2,"links"]
at<-which(links_profiles[[at_entry_2]]$adduc[,"linked profile"]==prof1)
links_profiles[[at_entry_2]]$adduc<-links_profiles[[at_entry_2]]$adduc[-at,,drop=FALSE]
}
links_profiles[[for_profs[n]]]$adduc<-
links_profiles[[for_profs[n]]]$adduc[keep,,drop=FALSE]
}
}
}
perce<-as.character(round(removed/found,digits=3)*100)
cat("\n");cat(paste0(perce, "% of ", found," adduct links filtered."));cat("\n");
return(links_profiles)
}
