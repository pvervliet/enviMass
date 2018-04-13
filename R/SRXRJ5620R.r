SRXRJ5620R <- function(
profileList,
prof_ID,
links_profiles,
min_peaks=3,
skip_peaks=FALSE,
min_cor=.9,
with_test=FALSE,
only_direct=FALSE,
del_RT=30,
omit_profiles=FALSE
){
prof_isot_IDs<-c()
if(omit_profiles[1]!="FALSE"){
if(length(omit_profiles)!=length(profileList[["index_prof"]][,"links"])){stop("omit_profiles not equal to number of profiles")}
}
if(with_test){if(prof_ID>length(profileList[["index_prof"]][,"profile_ID"])){stop("\n Debug get_isotopol _1!")}}
if(with_test){if(profileList[["index_prof"]][prof_ID,"profile_ID"]!=prof_ID){stop("\n Debug get_isotopol _2!")}}
if(profileList[["index_prof"]][prof_ID,"links"]==0){
return(prof_ID)
}
in_link<-profileList[["index_prof"]][prof_ID,"links"][[1]]
inter_links<-c()
not_inter_links<-c()
if(length(links_profiles[[in_link]]$isot[,1])>0){
if(skip_peaks){
those<-which(
((links_profiles[[in_link]]$isot[,"correl"]/1000)>=min_cor) &
(links_profiles[[in_link]]$isot[,"ref_1"]>=min_peaks)
)
if(length(those)>0){
those<-links_profiles[[in_link]]$isot[those,"linked profile"]
if(omit_profiles[1]!="FALSE"){
those<-those[omit_profiles[those]==0]
}
}
not_those<-which(
((links_profiles[[in_link]]$isot[,"correl"]/1000)<min_cor) |
(links_profiles[[in_link]]$isot[,"ref_1"]<min_peaks)
)
if(length(not_those)>0){
not_those<-links_profiles[[in_link]]$isot[not_those,"linked profile"]
}
}else{
those<-which(
((links_profiles[[in_link]]$isot[,"correl"]/1000)>=min_cor) |
(links_profiles[[in_link]]$isot[,"ref_1"]<min_peaks)
)
if(length(those)>0){
those<-links_profiles[[in_link]]$isot[those,"linked profile"]
if(omit_profiles[1]!="FALSE"){
those<-those[omit_profiles[those]==0]
}
}
not_those<-which(
((links_profiles[[in_link]]$isot[,"correl"]/1000)<min_cor) &
(links_profiles[[in_link]]$isot[,"ref_1"]>min_peaks)
)
if(length(not_those)>0){
not_those<-links_profiles[[in_link]]$isot[not_those,"linked profile"]
}
}
if(length(those)==0){
return(prof_ID)
}else{
inter_links<-c(inter_links,those)
not_inter_links<-c(not_inter_links,not_those)
if(only_direct){return(c(prof_ID,inter_links))}
}
}else{
return(prof_ID)
}
prof_isot_IDs<-c(prof_ID)
if(skip_peaks & (profileList[["index_prof"]][prof_ID,"number_peaks_total"][[1]]<min_peaks)){
return(prof_isot_IDs)
}
while(length(inter_links)>0){
in_link<-profileList[["index_prof"]][inter_links[1],"links"][[1]]
if(length(links_profiles[[in_link]]$isot[,1])>0){
those<-links_profiles[[in_link]]$isot[,"linked profile"]
those<-those[is.na(match(those,c(inter_links,prof_isot_IDs)))]
if(length(those)>0){
those<-those[is.na(match(those,c(not_inter_links)))]
}
if(length(those)>0 & omit_profiles[1]!="FALSE"){
those<-those[omit_profiles[those]==0]
}
if(length(those)>0 & del_RT!=FALSE){
those<-those[abs(profileList[["index_prof"]][prof_ID,"mean_RT"]-profileList[["index_prof"]][those,"mean_RT"])<=del_RT]
}
if(length(those)>0){
this<-profileList[["peaks"]][profileList[["index_prof"]][prof_ID,"start_ID"]:profileList[["index_prof"]][prof_ID,"end_ID"],"sampleIDs"]
del1<-profileList[["index_prof"]][prof_ID,"number_peaks_total"][[1]]
keep<-rep(FALSE,length(those))
for(k in 1:length(those)){
if(skip_peaks){
del2<-profileList[["index_prof"]][those[k],"number_peaks_total"][[1]]
if(del2<min_peaks){next}
that<-profileList[["peaks"]][profileList[["index_prof"]][those[k],"start_ID"]:profileList[["index_prof"]][those[k],"end_ID"],"sampleIDs"]
matched<-match(this,that)
if(sum(!is.na(matched)) < min_peaks) next
if(sum(!is.na(matched)) < 2) next
int1 <- ((profileList[["peaks"]][profileList[["index_prof"]][prof_ID,"start_ID"]:profileList[["index_prof"]][prof_ID,"end_ID"],"intensity"])[!is.na(matched)])
int2 <- ((profileList[["peaks"]][profileList[["index_prof"]][those[k],"start_ID"]:profileList[["index_prof"]][those[k],"end_ID"],"intensity"])[matched[!is.na(matched)]])
correl <- cor(int1, int2)
if(!is.na(correl)){
if(correl >= min_cor){
keep[k] <- TRUE
}
}
}else{
del2 <- profileList[["index_prof"]][those[k], "number_peaks_total"][[1]]
if(del2 < min_peaks){
keep[k] <- TRUE
next
}
that <- profileList[["peaks"]][profileList[["index_prof"]][those[k],"start_ID"]:profileList[["index_prof"]][those[k],"end_ID"],"sampleIDs"]
matched <- match(this,that)
if(sum(!is.na(matched)) < min_peaks){
keep[k] <- TRUE
next
}
if(sum(!is.na(matched)) < 2) next
int1 <- ((profileList[["peaks"]][profileList[["index_prof"]][prof_ID,"start_ID"]:profileList[["index_prof"]][prof_ID,"end_ID"],"intensity"])[!is.na(matched)])
int2 <- ((profileList[["peaks"]][profileList[["index_prof"]][those[k],"start_ID"]:profileList[["index_prof"]][those[k],"end_ID"],"intensity"])[matched[!is.na(matched)]])
correl <- cor(int1,int2)
if(!is.na(correl)){
if(correl >= min_cor){
keep[k]<-TRUE
}
}
}
}
if(length(those[keep])>0){
inter_links<-c(inter_links,those[keep])
}
if(length(those[!keep])>0){
not_inter_links<-c(not_inter_links,those[!keep])
}
}
}
prof_isot_IDs<-c(prof_isot_IDs,inter_links[1])
inter_links<-inter_links[-1]
}
if(with_test){
if(any(profileList[["index_prof"]][prof_isot_IDs,"profile_ID"]!=prof_isot_IDs)){
stop("\n Debug get_isotopol _3!")
}
}
return(prof_isot_IDs)
}
