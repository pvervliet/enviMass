RDQZL3486I <- function(
prof_IDs,
profileList,
links_profiles
){
if(length(prof_IDs)==1){return(TRUE)}
found<-rep(FALSE,length(prof_IDs));
started<-prof_IDs[1]
while(length(started)>0){
in_link<-profileList[["index_prof"]][started[1],"links"][[1]]
if(length(links_profiles[[in_link]]$isot[,1])>0){
those1<-links_profiles[[in_link]]$isot[,"linked profile"]
}else{those1<-c()}
if(length(links_profiles[[in_link]]$adduc[,1])>0){
those2<-links_profiles[[in_link]]$adduc[,"linked profile"]
}else{those2<-c()}
if(length(links_profiles[[in_link]]$EIC[,1])>0){
those3<-links_profiles[[in_link]]$EIC[,"linked profile"]
}else{those3<-c()}
those<-c(those1,those2,those3)
those<-unique(those)
those<-those[!is.na(match(those,prof_IDs))]
if(length(those)>0){
those<-those[found[match(those,prof_IDs)]==FALSE]
those<-those[!is.na(those)]
started<-c(started,those)
started<-unique(started)
}
found[prof_IDs==started[1]]<-TRUE;
started<-started[-1]
}
if(all(found)){
return(TRUE)
}else{
return(FALSE)
}
}
