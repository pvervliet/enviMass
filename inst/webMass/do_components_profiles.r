if(file.exists(file.path(as.character(logfile[[1]]),"results","links_profiles_pos"))){file.remove(file.path(as.character(logfile[[1]]),"results","links_profiles_pos"))}
if(file.exists(file.path(as.character(logfile[[1]]),"results","links_profiles_neg"))){file.remove(file.path(as.character(logfile[[1]]),"results","links_profiles_neg"))}
if(
file.exists(file.path(logfile[[1]],"results","profileList_pos"))
){
if(any(objects(envir=as.environment(".GlobalEnv"))=="profileList_pos")){rm(profileList_pos,envir=as.environment(".GlobalEnv"))}
if(any(objects()=="profileList_pos")){rm(profileList_pos)}
if(any(objects(envir=as.environment(".GlobalEnv"))=="links_peaks_pos")){rm(links_peaks_pos,envir=as.environment(".GlobalEnv"))}
if(any(objects()=="links_peaks_pos")){rm(links_peaks_pos)}
if(any(objects(envir=as.environment(".GlobalEnv"))=="links_profiles_pos")){rm(links_profiles_pos,envir=as.environment(".GlobalEnv"))}
if(any(objects()=="links_profiles_pos")){rm(links_profiles_pos)}
load(file.path(as.character(logfile[[1]]),"results","profileList_pos"),envir=as.environment(".GlobalEnv"));
load(file.path(as.character(logfile[[1]]),"results","links_peaks_pos"), envir = as.environment(".GlobalEnv"));
assign("links_profiles_pos", list(), envir = as.environment(".GlobalEnv"))
measurements <- read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
peaks <- profileList_pos[["peaks"]][,c("sampleIDs","peakIDs","profileIDs","RT")]
ord <- order(peaks[,"sampleIDs"], peaks[,"peakIDs"], peaks[,"profileIDs"], decreasing = FALSE)
peaks <- peaks[ord,]
use_entries_profiles <- enviMass::VOLZN1972A(links_profiles_pos)
profileList_pos[["index_prof"]][,"links"] <<- 0
with_bar<-TRUE
if(
(
(logfile$workflow[names(logfile$workflow)=="subtr"]=="no") |
(
((logfile$workflow[names(logfile$workflow)=="IS_screen"]=="yes") &  (logfile$parameters$subtr_IS!="yes")) |
((logfile$workflow[names(logfile$workflow)=="target_screen"]=="yes") & (logfile$parameters$subtr_target!="yes"))
)
) &
(length(links_peaks_pos)>0)
){
cat("\n Annotation of screening results to profiles")
if(with_bar){pBar <- txtProgressBar(min = 0, max = dim(profileList_pos[["index_prof"]])[1], style = 3)}
for(i in 1:dim(profileList_pos[["index_prof"]])[1]){
if(with_bar){setTxtProgressBar(pBar, i, title = NULL, label = NULL)}
if(
any(profileList_pos[["peaks"]][
(profileList_pos[["index_prof"]][i,"start_ID"]:profileList_pos[["index_prof"]][i,"end_ID"]),"links"
]!=0)
){
if( profileList_pos[["index_prof"]][i,"links"]==0 ){
if(length(use_entries_profiles)>0){
at_entry<-use_entries_profiles[1]
use_entries<-use_entries_profiles[-1]
}else{
at_entry<-(length(links_profiles_pos)+1)
}
links_profiles_pos[[at_entry]] <<- enviMass::RHHPS2808T(profileList_pos[["index_prof"]][i,"number_peaks_total"][[1]])
names(links_profiles_pos)[at_entry] <<- as.character(i)
profileList_pos[["index_prof"]][i,"links"] <<- at_entry
}else{
at_entry<-profileList_pos[["index_prof"]][i,"links"]
}
for(j in (profileList_pos[["index_prof"]][i,"start_ID"]:profileList_pos[["index_prof"]][i, "end_ID"])){
if(profileList_pos[["peaks"]][j,"links"] != 0){
if(	length(links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[2]]) > 0 ){
for(k in 1:length(links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[2]])){
if( length(links_profiles_pos[[at_entry]][[2]]) == 0 ){
links_profiles_pos[[at_entry]][[2]] <<-
data.frame(
links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[2]][[k]][[1]],
1,
links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[2]][[k]][[2]],
stringsAsFactors = FALSE
)
names(links_profiles_pos[[at_entry]][[2]]) <<- c("Compound","Counts","max_score")
}else{
at<-which(links_profiles_pos[[at_entry]][[2]][,1] == links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[2]][[k]][[1]])
if(length(at) > 0){
links_profiles_pos[[at_entry]][[2]][at,2] <<- (links_profiles_pos[[at_entry]][[2]][at,2]+1)
if(links_profiles_pos[[at_entry]][[2]][at,2] < links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[2]][[k]][[2]] ){
links_profiles_pos[[at_entry]][[2]][at,2] <<- links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[2]][[k]][[2]]
}
}else{
links_profiles_pos[[at_entry]][[2]] <<- data.frame(
c(
links_profiles_pos[[at_entry]][[2]][,1],
links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[2]][[k]][[1]]
),
c(links_profiles_pos[[at_entry]][[2]][,2],1),
c(
links_profiles_pos[[at_entry]][[2]][,3],
links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[2]][[k]][[2]]
),
stringsAsFactors = FALSE
)
names(links_profiles_pos[[at_entry]][[2]]) <<- c("Compound","Counts","max_score")
}
}
}
}
if(	length(links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[1]]) > 0 ){
for(k in 1:length(links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[1]])){
if( length(links_profiles_pos[[at_entry]][[1]])==0 ){
links_profiles_pos[[at_entry]][[1]] <<-
data.frame(
links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[1]][[k]][[1]],
1,
links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[1]][[k]][[2]],
stringsAsFactors = FALSE
)
names(links_profiles_pos[[at_entry]][[1]]) <<- c("Compound","Counts","max_score")
}else{
at <- which(links_profiles_pos[[at_entry]][[1]][,1] == links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[1]][[k]][[1]])
if(length(at) > 0){
links_profiles_pos[[at_entry]][[1]][at,2] <<- (links_profiles_pos[[at_entry]][[1]][at,2]+1)
if(links_profiles_pos[[at_entry]][[1]][at,2] < links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[1]][[k]][[2]] ){
links_profiles_pos[[at_entry]][[1]][at,2] <<- links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[1]][[k]][[2]]
}
}else{
links_profiles_pos[[at_entry]][[1]] <<- data.frame(
c(
links_profiles_pos[[at_entry]][[1]][,1],
links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[1]][[k]][[1]]
),
c(links_profiles_pos[[at_entry]][[1]][,2],1),
c(
links_profiles_pos[[at_entry]][[1]][,3],
links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[1]][[k]][[2]]
),
stringsAsFactors = FALSE
)
names(links_profiles_pos[[at_entry]][[1]]) <<- c("Compound","Counts","max_score")
}
}
}
}
}
}
}
}
if(with_bar){close(pBar)}
cat(" done.")
}
if(logfile$workflow[names(logfile$workflow) == "EIC_correlation"] == "yes"){
cat("\n Retrieving EIC correlation links ")
forIDs <- profileList_pos[["sampleID"]]
for_files <- list.files(file.path(logfile[[1]],"results","componentization","EIC_corr"))
keep <- match(forIDs,for_files)
if(any(is.na(keep))){cat("\n Just note: not all files found in profiles have EIC correlation results (1). \n")}
TRUE_IDs <- (measurements$ID[measurements$EIC_correlation == "TRUE"])
keep2 <- match(forIDs, TRUE_IDs)
forIDs <- forIDs[!is.na(keep) & !is.na(keep2)]
if(logfile$parameters$dofile_latest_profcomp == "TRUE"){
atPOSIX <- profileList_pos[["datetime"]];
matchID <- profileList_pos[["sampleID"]];
atdate <- c(); attime <- c();
for(i in 1:length(atPOSIX)){
atdate <- c(atdate, strsplit(atPOSIX[i]," ")[[1]][1]);
attime <- c(attime, strsplit(atPOSIX[i]," ")[[1]][2]);
}
attime <- as.difftime(attime);
atdate <- as.Date(atdate, tz="GMT");
ord <- order(as.numeric(atdate), as.numeric(attime), matchID, decreasing=TRUE);
matchID <- matchID[ord];
forIDs <- forIDs[match(forIDs, matchID)]
if(length(forIDs) > as.numeric(logfile$parameters$numfile_latest_profcomp)){
forIDs <- forIDs[1:as.numeric(logfile$parameters$numfile_latest_profcomp)]
}
forIDs <- forIDs[!is.na(forIDs)]
}
forIDs <- as.numeric(forIDs)
not_found1 <- 0; inserted1 <- 0
if(length(forIDs) > 0){
if(with_bar){pBar <- txtProgressBar(min = 0, max = length(forIDs), style = 3)}
for(i in 1:length(forIDs)){
if(with_bar){setTxtProgressBar(pBar, i, title = NULL, label = NULL)}
load(file = file.path(logfile[[1]], "results", "componentization", "EIC_corr", as.character(forIDs[i])))
if(length(EIC_pairs[,1]) == 0){next}
get1 <- cbind(
rep(forIDs[i], length(EIC_pairs[,1])), EIC_pairs[,1]
)
found1 <- enviMass::GSZEL3799Y(get1, peaks[,c("sampleIDs","peakIDs")], row_order = FALSE, column_order_a = FALSE, column_order_b = FALSE, get_index = TRUE)
get2 <- cbind(
rep(forIDs[i],length(EIC_pairs[,2])), EIC_pairs[,2]
)
found2 <- enviMass::GSZEL3799Y(get2, peaks[,c("sampleIDs", "peakIDs")], row_order = FALSE, column_order_a = TRUE, column_order_b = FALSE, get_index = TRUE)
for(j in 1:length(found1)){
if(found1[j] == 0){not_found1 <- (not_found1 + 1); next}
if(found2[j] == 0){not_found1 <- (not_found1 + 1); next}
inserted1 <- (inserted1+1);
prof1 <- peaks[found1[j], "profileIDs"][[1]]
prof2 <- peaks[found2[j], "profileIDs"][[1]]
if(FALSE){
cat("\n");
cat(peaks[found1[j],"RT"]); cat(" - ")
cat(peaks[found2[j],"RT"])
}
if(profileList_pos[["index_prof"]][prof1,"profile_ID"] != prof1){stop("\nComponentization: debug me, _1!")}
if(profileList_pos[["index_prof"]][prof1,"links"] == 0){
if(length(use_entries_profiles) > 0){
at_entry_1 <- use_entries_profiles[1]
use_entries <- use_entries_profiles[-1]
}else{
at_entry_1 <- (length(links_profiles_pos)+1)
}
links_profiles_pos[[at_entry_1]] <<- enviMass::RHHPS2808T(profileList_pos[["index_prof"]][prof1,"number_peaks_total"][[1]])
names(links_profiles_pos)[at_entry_1] <<- as.character(prof1)
profileList_pos[["index_prof"]][prof1,"links"] <<- at_entry_1
}else{
at_entry_1 <- profileList_pos[["index_prof"]][prof1,"links"]
}
here1 <- which(links_profiles_pos[[at_entry_1]][["EIC"]][,"linked profile"] == prof2)
if(length(here1) == 0){
links_profiles_pos[[at_entry_1]][["EIC"]] <<- rbind(
links_profiles_pos[[at_entry_1]][["EIC"]], c(prof2,1,0,0,0,0,1,NA)
)
here1 <- dim(links_profiles_pos[[at_entry_1]][["EIC"]])[1]
}else{
links_profiles_pos[[at_entry_1]][["EIC"]][here1,"link counts"] <<- (links_profiles_pos[[at_entry_1]][["EIC"]][here1,"link counts"] + 1)
}
if(profileList_pos[["index_prof"]][prof2,"profile_ID"]!=prof2){stop("\nCross-profile componentization: debug me, EIC_1!")}
if(profileList_pos[["index_prof"]][prof2,"links"]==0){
if(length(use_entries_profiles)>0){
at_entry_2<-use_entries_profiles[1]
use_entries<-use_entries_profiles[-1]
}else{
at_entry_2<-(length(links_profiles_pos)+1)
}
links_profiles_pos[[at_entry_2]] <<- enviMass::RHHPS2808T(profileList_pos[["index_prof"]][prof2,"number_peaks_total"][[1]])
names(links_profiles_pos)[at_entry_2] <<- as.character(prof2)
profileList_pos[["index_prof"]][prof2,"links"] <<- at_entry_2
}else{
at_entry_2<-profileList_pos[["index_prof"]][prof2,"links"]
}
here2<-which(links_profiles_pos[[at_entry_2]][["EIC"]][,"linked profile"]==prof1)
if(length(here2)==0){
links_profiles_pos[[at_entry_2]][["EIC"]] <<- rbind(
links_profiles_pos[[at_entry_2]][["EIC"]], c(prof1,1,0,0,0,0,0,NA)
)
here2<-dim(links_profiles_pos[[at_entry_2]][["EIC"]])[1]
}else{
links_profiles_pos[[at_entry_2]][["EIC"]][here2,"link counts"] <<- (links_profiles_pos[[at_entry_2]][["EIC"]][here2,"link counts"]+1)
}
if(links_profiles_pos[[at_entry_1]][["EIC"]][here1,"use"]==1){
if(length(links_profiles_pos[[at_entry_1]]$EIC_cor) < here1){
links_profiles_pos[[at_entry_1]]$EIC_cor[[here1]] <<- EIC_pairs[j,4]
}else{
if(is.null(links_profiles_pos[[at_entry_1]]$EIC_cor[[here1]])){
links_profiles_pos[[at_entry_1]]$EIC_cor[[here1]] <<- EIC_pairs[j,4]
}else{
links_profiles_pos[[at_entry_1]]$EIC_cor[[here1]] <<- c(links_profiles_pos[[at_entry_1]]$EIC_cor[[here1]],EIC_pairs[j,4])
}
}
}else{
if(length(links_profiles_pos[[at_entry_2]]$EIC_cor) < here2){
links_profiles_pos[[at_entry_2]]$EIC_cor[[here2]] <<- EIC_pairs[j,4]
}else{
if(is.null(links_profiles_pos[[at_entry_2]]$EIC_cor[[here2]])){
links_profiles_pos[[at_entry_2]]$EIC_cor[[here2]] <<- EIC_pairs[j,4]
}else{
links_profiles_pos[[at_entry_2]]$EIC_cor[[here2]] <<- c(links_profiles_pos[[at_entry_2]]$EIC_cor[[here2]],EIC_pairs[j,4])
}
}
}
}
}
if(with_bar){close(pBar)}
if(any(objects()=="EIC_pairs")){rm(EIC_pairs)}
cat(" done.")
}
}
if(logfile$workflow[names(logfile$workflow) == "isotopologues"] == "yes"){
forIDs <- profileList_pos[["sampleID"]]
for_files <- list.files(file.path(logfile[[1]], "results", "componentization", "isotopologues"))
keep <- match(forIDs, for_files)
if(any(is.na(keep))){cat("\n Just note: not all files found in profiles have isotopologue links available. \n")}
TRUE_IDs <- (measurements$ID[measurements$isotopologues == "TRUE"])
keep2 <- match(forIDs, TRUE_IDs)
forIDs <- forIDs[!is.na(keep) & !is.na(keep2)]
if(logfile$parameters$dofile_latest_profcomp == "TRUE"){
atPOSIX <- profileList_pos[["datetime"]];
matchID <- profileList_pos[["sampleID"]];
atdate <- c(); attime <- c();
for(i in 1:length(atPOSIX)){
atdate <- c(atdate, strsplit(atPOSIX[i]," ")[[1]][1]);
attime <- c(attime, strsplit(atPOSIX[i]," ")[[1]][2]);
}
attime <- as.difftime(attime);
atdate <- as.Date(atdate, tz = "GMT");
ord <- order(as.numeric(atdate), as.numeric(attime), matchID, decreasing = TRUE);
matchID <- matchID[ord];
forIDs <- forIDs[match(forIDs, matchID)]
if(length(forIDs) > as.numeric(logfile$parameters$numfile_latest_profcomp)){
forIDs <- forIDs[1:as.numeric(logfile$parameters$numfile_latest_profcomp)]
}
forIDs <- forIDs[!is.na(forIDs)]
}
forIDs <- as.numeric(forIDs)
not_found3 <- 0; inserted3 <- 0
if(length(forIDs) > 0){
cat("\n Retrieving isotopologue links ")
if(with_bar){pBar <- txtProgressBar(min = 0, max = length(forIDs), style = 3)}
for(i in 1:length(forIDs)){
if(with_bar){setTxtProgressBar(pBar, i, title = NULL, label = NULL)}
load(file = file.path(logfile[[1]], "results", "componentization", "isotopologues", as.character(forIDs[i])))
if(length(Isot_pairs[,1]) == 0){next}
get1 <- cbind(
rep(forIDs[i], length(Isot_pairs[,1])), Isot_pairs[,1]
)
found1 <- enviMass::GSZEL3799Y(get1, peaks[,c("sampleIDs", "peakIDs")], row_order = FALSE, column_order_a = FALSE, column_order_b = FALSE, get_index = TRUE)
get2 <- cbind(
rep(forIDs[i], length(Isot_pairs[,2])), Isot_pairs[,2]
)
found2 <- enviMass::GSZEL3799Y(get2, peaks[,c("sampleIDs", "peakIDs")], row_order = FALSE, column_order_a = TRUE, column_order_b = FALSE, get_index = TRUE)
for(j in 1:length(found1)){
if(found1[j] == 0){not_found3 <- (not_found3 + 1); next}
if(found2[j] == 0){not_found3 <- (not_found3 + 1); next}
inserted3 <- (inserted3 + 1);
if(FALSE){
cat("\n");
cat(peaks[found1[j], "RT"]); cat(" - ");
cat(peaks[found2[j], "RT"])
}
prof1 <- peaks[found1[j], "profileIDs"][[1]]
prof2 <- peaks[found2[j], "profileIDs"][[1]]
if(profileList_pos[["index_prof"]][prof1, "profile_ID"] != prof1){stop("\nComponentization: debug me, _1!")}
if(profileList_pos[["index_prof"]][prof1, "links"] == 0){
if(length(use_entries_profiles) > 0){
at_entry_1 <- use_entries_profiles[1]
use_entries <- use_entries_profiles[-1]
}else{
at_entry_1 <- (length(links_profiles_pos) + 1)
}
links_profiles_pos[[at_entry_1]] <<- enviMass::RHHPS2808T(profileList_pos[["index_prof"]][prof1, "number_peaks_total"][[1]])
names(links_profiles_pos)[at_entry_1] <<- as.character(prof1)
profileList_pos[["index_prof"]][prof1, "links"] <<- at_entry_1
}else{
at_entry_1<-profileList_pos[["index_prof"]][prof1, "links"]
}
here1 <- which(links_profiles_pos[[at_entry_1]][["isot"]][, "linked profile"] == prof2)
if(length(here1)==0){
links_profiles_pos[[at_entry_1]][["isot"]] <<- rbind(
links_profiles_pos[[at_entry_1]][["isot"]], c(prof2,1,0,1,NA)
)
here1 <- dim(links_profiles_pos[[at_entry_1]][["isot"]])[1]
is_new1 <- TRUE
}else{
links_profiles_pos[[at_entry_1]][["isot"]][here1, "link counts"] <<- (links_profiles_pos[[at_entry_1]][["isot"]][here1,"link counts"] + 1)
is_new1 <- FALSE
}
if(profileList_pos[["index_prof"]][prof2, "profile_ID"] != prof2){stop("\nComponentization: debug me, _1!")}
if(profileList_pos[["index_prof"]][prof2, "links"] == 0){
if(length(use_entries_profiles) > 0){
at_entry_2 <- use_entries_profiles[1]
use_entries <- use_entries_profiles[-1]
}else{
at_entry_2 <- (length(links_profiles_pos)+1)
}
links_profiles_pos[[at_entry_2]] <<- enviMass::RHHPS2808T(profileList_pos[["index_prof"]][prof2,"number_peaks_total"][[1]])
names(links_profiles_pos)[at_entry_2] <<- as.character(prof2)
profileList_pos[["index_prof"]][prof2,"links"] <<- at_entry_2
}else{
at_entry_2 <- profileList_pos[["index_prof"]][prof2,"links"]
}
here2<-which(links_profiles_pos[[at_entry_2]][["isot"]][,"linked profile"] == prof1)
if(length(here2) == 0){
links_profiles_pos[[at_entry_2]][["isot"]] <<- rbind(
links_profiles_pos[[at_entry_2]][["isot"]], c(prof1,1,0,0,NA)
)
here2 <- dim(links_profiles_pos[[at_entry_2]][["isot"]])[1]
is_new2 <- TRUE
}else{
links_profiles_pos[[at_entry_2]][["isot"]][here2,"link counts"] <<- (links_profiles_pos[[at_entry_2]][["isot"]][here2,"link counts"]+1)
is_new2 <- FALSE
}
if(is_new1 | is_new2){
these <- profileList_pos[["peaks"]][
profileList_pos[["index_prof"]][prof1,"start_ID"]:profileList_pos[["index_prof"]][prof1,"end_ID"]
,"sampleIDs"]
these <- these[!is.na(match(these,forIDs))]
those <- profileList_pos[["peaks"]][
profileList_pos[["index_prof"]][prof2,"start_ID"]:profileList_pos[["index_prof"]][prof2,"end_ID"]
,"sampleIDs"]
those <- those[!is.na(match(these,forIDs))]
matched <- match(these,those)
not_NA <- sum(!is.na(matched))
links_profiles_pos[[at_entry_1]]$isot[here1,"ref_1"] <<- not_NA
links_profiles_pos[[at_entry_2]]$isot[here2,"ref_1"] <<- not_NA
}
if(any(links_profiles_pos[[at_entry_1]]$isot[,"ref_1"] == 0)){stop("\n DEBUG ME ! FOUND_1")}
if(any(links_profiles_pos[[at_entry_2]]$isot[,"ref_1"] == 0)){stop("\n DEBUG ME ! FOUND_2")}
}
}
if(with_bar){close(pBar)}
if(any(objects()=="Isot_pairs")){rm(Isot_pairs)}
cat(" done.")
}
}
if(logfile$workflow[names(logfile$workflow)=="adducts"] == "yes"){
forIDs<-profileList_pos[["sampleID"]]
for_files<-list.files(file.path(logfile[[1]],"results","componentization","adducts"))
keep<-match(forIDs,for_files)
if(any(is.na(keep))){cat("\n Just note: not all files found in profiles have adduct links available. \n")}
TRUE_IDs <- (measurements$ID[measurements$adducts=="TRUE"])
keep2 <- match(forIDs, TRUE_IDs)
forIDs<-forIDs[!is.na(keep) & !is.na(keep2)]
if(logfile$parameters$dofile_latest_profcomp=="TRUE"){
atPOSIX<-profileList_pos[["datetime"]];
matchID<-profileList_pos[["sampleID"]];
atdate<-c();attime<-c();
for(i in 1:length(atPOSIX)){
atdate<-c(atdate, strsplit(atPOSIX[i]," ")[[1]][1]);
attime<-c(attime, strsplit(atPOSIX[i]," ")[[1]][2]);
}
attime<-as.difftime(attime);
atdate<-as.Date(atdate, tz = "GMT");
ord<-order(as.numeric(atdate),as.numeric(attime),matchID,decreasing=TRUE);
matchID<-matchID[ord];
forIDs<-forIDs[match(forIDs,matchID)]
if(length(forIDs)>as.numeric(logfile$parameters$numfile_latest_profcomp)){
forIDs<-forIDs[1:as.numeric(logfile$parameters$numfile_latest_profcomp)]
}
forIDs <- forIDs[!is.na(forIDs)]
}
forIDs <- as.numeric(forIDs)
not_found4 <- 0; inserted4 <- 0
if(length(forIDs) > 0){
cat("\n Retrieving adduct links ")
if(with_bar){pBar <- txtProgressBar(min = 0, max = length(forIDs), style = 3)}
for(i in 1:length(forIDs)){
if(with_bar){setTxtProgressBar(pBar, i, title = NULL, label = NULL)}
load(file=file.path(logfile[[1]],"results","componentization","adducts",as.character(forIDs[i])))
if(length(Adduct_pairs[,1])==0){next}
get1<-cbind(
rep(forIDs[i],length(Adduct_pairs[,1])),Adduct_pairs[,1]
)
found1<-enviMass::GSZEL3799Y(get1,peaks[,c("sampleIDs","peakIDs")],row_order=FALSE,column_order_a=FALSE,column_order_b=FALSE,get_index=TRUE)
get2<-cbind(
rep(forIDs[i],length(Adduct_pairs[,2])),Adduct_pairs[,2]
)
found2<-enviMass::GSZEL3799Y(get2,peaks[,c("sampleIDs","peakIDs")],row_order=FALSE,column_order_a=TRUE,column_order_b=FALSE,get_index=TRUE)
for(j in 1:length(found1)){
if(found1[j]==0){not_found4<-(not_found4+1);next}
if(found2[j]==0){not_found4<-(not_found4+1);next}
inserted4<-(inserted4+1);
if(FALSE){
cat("\n");
cat(peaks[found1[j],"RT"]);cat(" - ")
cat(peaks[found2[j],"RT"])
}
prof1<-peaks[found1[j],"profileIDs"][[1]]
prof2<-peaks[found2[j],"profileIDs"][[1]]
if(profileList_pos[["index_prof"]][prof1,"profile_ID"]!=prof1){stop("\nComponentization: debug me, _1!")}
if(profileList_pos[["index_prof"]][prof1,"links"]==0){
if(length(use_entries_profiles)>0){
at_entry_1<-use_entries_profiles[1]
use_entries<-use_entries_profiles[-1]
}else{
at_entry_1<-(length(links_profiles_pos)+1)
}
links_profiles_pos[[at_entry_1]] <<- enviMass::RHHPS2808T(profileList_pos[["index_prof"]][prof1,"number_peaks_total"][[1]])
names(links_profiles_pos)[at_entry_1] <<- as.character(prof1)
profileList_pos[["index_prof"]][prof1,"links"] <<- at_entry_1
}else{
at_entry_1<-profileList_pos[["index_prof"]][prof1,"links"]
}
here1<-which(links_profiles_pos[[at_entry_1]][["adduc"]][,"linked profile"]==prof2)
if(length(here1)==0){
links_profiles_pos[[at_entry_1]][["adduc"]] <<- rbind(
links_profiles_pos[[at_entry_1]][["adduc"]], c(prof2,1,0,1,NA)
)
here1<-dim(links_profiles_pos[[at_entry_1]][["adduc"]])[1]
is_new1 <- TRUE
}else{
links_profiles_pos[[at_entry_1]][["adduc"]][here1,"link counts"] <<- (links_profiles_pos[[at_entry_1]][["adduc"]][here1,"link counts"]+1)
is_new1 <- FALSE
}
if(profileList_pos[["index_prof"]][prof2,"profile_ID"]!=prof2){stop("\nComponentization: debug me, _1!")}
if(profileList_pos[["index_prof"]][prof2,"links"]==0){
if(length(use_entries_profiles)>0){
at_entry_2<-use_entries_profiles[1]
use_entries<-use_entries_profiles[-1]
}else{
at_entry_2<-(length(links_profiles_pos)+1)
}
links_profiles_pos[[at_entry_2]] <<- enviMass::RHHPS2808T(profileList_pos[["index_prof"]][prof2,"number_peaks_total"][[1]])
names(links_profiles_pos)[at_entry_2] <<- as.character(prof2)
profileList_pos[["index_prof"]][prof2,"links"] <<- at_entry_2
}else{
at_entry_2 <- profileList_pos[["index_prof"]][prof2,"links"]
}
here2 <- which(links_profiles_pos[[at_entry_2]][["adduc"]][,"linked profile"]==prof1)
if(length(here2)==0){
links_profiles_pos[[at_entry_2]][["adduc"]] <<- rbind(
links_profiles_pos[[at_entry_2]][["adduc"]], c(prof1,1,0,0,NA)
)
here2 <- dim(links_profiles_pos[[at_entry_2]][["adduc"]])[1]
is_new2 <- TRUE
}else{
links_profiles_pos[[at_entry_2]][["adduc"]][here2,"link counts"] <<- (links_profiles_pos[[at_entry_2]][["adduc"]][here2,"link counts"]+1)
is_new2 <- FALSE
}
if(is_new1 | is_new2){
these<-profileList_pos[["peaks"]][
profileList_pos[["index_prof"]][prof1,"start_ID"]:profileList_pos[["index_prof"]][prof1,"end_ID"]
,"sampleIDs"]
these <- these[!is.na(match(these,forIDs))]
those<-profileList_pos[["peaks"]][
profileList_pos[["index_prof"]][prof2,"start_ID"]:profileList_pos[["index_prof"]][prof2,"end_ID"]
,"sampleIDs"]
those <- those[!is.na(match(those,forIDs))]
matched <- match(these,those)
not_NA <- sum(!is.na(matched))
links_profiles_pos[[at_entry_1]]$adduc[here1,"ref_1"] <<- not_NA
links_profiles_pos[[at_entry_2]]$adduc[here2,"ref_1"] <<- not_NA
}
if(any(links_profiles_pos[[at_entry_1]]$adduc[,"ref_1"] == 0)){stop("\n DEBUG ME ! FOUND_3")}
if(any(links_profiles_pos[[at_entry_2]]$adduc[,"ref_1"] == 0)){stop("\n DEBUG ME ! FOUND_4")}
}
}
if(with_bar){close(pBar)}
if(any(objects()=="Adduct_pairs")){rm(Adduct_pairs)}
cat(" done.")
}
}
if(logfile$workflow[names(logfile$workflow)=="homologues"] == "yes"){
forIDs<-profileList_pos[["sampleID"]]
for_files<-list.files(file.path(logfile[[1]],"results","componentization","homologues"))
keep<-match(forIDs,for_files)
if(any(is.na(keep))){cat("\n Just note: not all files found in profiles have homologue links available. \n")}
TRUE_IDs <- (measurements$ID[measurements$homologues=="TRUE"])
keep2 <- match(forIDs, TRUE_IDs)
forIDs <- forIDs[!is.na(keep) & !is.na(keep2)]
if(logfile$parameters$dofile_latest_profcomp=="TRUE"){
atPOSIX<-profileList_pos[["datetime"]];
matchID<-profileList_pos[["sampleID"]];
atdate<-c();attime<-c();
for(i in 1:length(atPOSIX)){
atdate<-c(atdate, strsplit(atPOSIX[i]," ")[[1]][1]);
attime<-c(attime, strsplit(atPOSIX[i]," ")[[1]][2]);
}
attime<-as.difftime(attime);
atdate<-as.Date(atdate, tz="GMT");
ord<-order(as.numeric(atdate),as.numeric(attime),matchID,decreasing=TRUE);
matchID<-matchID[ord];
forIDs<-forIDs[match(forIDs,matchID)]
if(length(forIDs)>as.numeric(logfile$parameters$numfile_latest_profcomp)){
forIDs<-forIDs[1:as.numeric(logfile$parameters$numfile_latest_profcomp)]
}
forIDs <- forIDs[!is.na(forIDs)]
}
forIDs <- as.numeric(forIDs)
not_found5<-0;inserted5<-0
if(length(forIDs)>0){
cat("\n Retrieving homologue links ")
if(with_bar){pBar <- txtProgressBar(min = 0, max = length(forIDs), style = 3)}
for(i in 1:length(forIDs)){
if(with_bar){setTxtProgressBar(pBar, i, title = NULL, label = NULL)}
load(file=file.path(logfile[[1]],"results","componentization","homologues",as.character(forIDs[i])))
if(length(Homol_groups[,1])==0){next}
keep_group<-rep(TRUE,max(Homol_groups[,3]))
get1<-cbind(
rep(forIDs[i],length(Homol_groups[,1])),Homol_groups[,1]
)
found1<-enviMass::GSZEL3799Y(get1,peaks[,c("sampleIDs","peakIDs")],row_order=FALSE,column_order_a=FALSE,column_order_b=FALSE,get_index=TRUE)
keep_group[Homol_groups[found1==0,3]]<-FALSE
if(!any(keep_group)){next}
get2<-cbind(
rep(forIDs[i],length(Homol_groups[,2])),Homol_groups[,2]
)
found2<-enviMass::GSZEL3799Y(get2,peaks[,c("sampleIDs","peakIDs")],row_order=FALSE,column_order_a=TRUE,column_order_b=FALSE,get_index=TRUE)
if(!any(keep_group)){next}
for(j in 1:length(found1)){
if(found1[j]==0){not_found5<-(not_found5+1);next}
if(found2[j]==0){not_found5<-(not_found5+1);next}
if(Homol_groups[j,3]!=Homol_groups[j,3]){stop("\n Debug componentization for its homologue relations!")}
if(!keep_group[Homol_groups[j,3]]){next}
inserted5<-(inserted5+1);
if(FALSE){
cat("\n");
cat(peaks[found1[j],"RT"]);cat(" - ")
cat(peaks[found2[j],"RT"])
}
prof1<-peaks[found1[j],"profileIDs"][[1]]
prof2<-peaks[found2[j],"profileIDs"][[1]]
if(profileList_pos[["index_prof"]][prof1,"profile_ID"]!=prof1){stop("\nComponentization: debug me, _1!")}
if(profileList_pos[["index_prof"]][prof1,"links"]==0){
if(length(use_entries_profiles)>0){
at_entry<-use_entries_profiles[1]
use_entries<-use_entries_profiles[-1]
}else{
at_entry<-(length(links_profiles_pos)+1)
}
links_profiles_pos[[at_entry]] <<- enviMass::RHHPS2808T(profileList_pos[["index_prof"]][prof1,"number_peaks_total"][[1]])
names(links_profiles_pos)[at_entry] <<- as.character(prof1)
profileList_pos[["index_prof"]][prof1,"links"] <<- at_entry
}else{
at_entry <- profileList_pos[["index_prof"]][prof1,"links"]
}
here <- which(links_profiles_pos[[at_entry]][["homol"]][,"linked profile"]==prof2)
if(length(here)==0){
links_profiles_pos[[at_entry]][["homol"]] <<- rbind(
links_profiles_pos[[at_entry]][["homol"]], c(prof2,1,0)
)
}else{
links_profiles_pos[[at_entry]][["homol"]][here,"link counts"] <<- (links_profiles_pos[[at_entry]][["homol"]][here,"link counts"]+1)
}
if(profileList_pos[["index_prof"]][prof2,"profile_ID"]!=prof2){stop("\nComponentization: debug me, _1!")}
if(profileList_pos[["index_prof"]][prof2,"links"]==0){
if(length(use_entries_profiles)>0){
at_entry<-use_entries_profiles[1]
use_entries<-use_entries_profiles[-1]
}else{
at_entry<-(length(links_profiles_pos)+1)
}
links_profiles_pos[[at_entry]] <<- enviMass::RHHPS2808T(profileList_pos[["index_prof"]][prof2,"number_peaks_total"][[1]])
names(links_profiles_pos)[at_entry] <<- as.character(prof2)
profileList_pos[["index_prof"]][prof2,"links"] <<- at_entry
}else{
at_entry <- profileList_pos[["index_prof"]][prof2,"links"]
}
here <- which(links_profiles_pos[[at_entry]][["homol"]][,"linked profile"]==prof1)
if(length(here)==0){
links_profiles_pos[[at_entry]][["homol"]] <<- rbind(
links_profiles_pos[[at_entry]][["homol"]], c(prof1,1,0)
)
}else{
links_profiles_pos[[at_entry]][["homol"]][here,"link counts"] <<- (links_profiles_pos[[at_entry]][["homol"]][here,"link counts"]+1)
}
}
}
if(with_bar){close(pBar)}
if(any(objects()=="Homol_groups")){rm(Homol_groups)}
cat(" done.")
}
}
cut_delRT_EIC <<- NA
if(
(logfile$parameters$filter_profcomp_pos == "TRUE") &
( 	(logfile$workflow[names(logfile$workflow) == "isotopologues"] == "yes") ||
(logfile$workflow[names(logfile$workflow)=="adducts"]=="yes")	)
){
anaA <- enviMass:::XFDKJ7277K(
links_profiles = links_profiles_pos,
profileList = profileList_pos,
logfile = logfile,
min_rat = .7,
min_count = .4,
perc = .9,
for_which = logfile$parameters$for_which_profcomp_pos
)
cut_delRT_isot <<- boxplot.stats(c(anaA$delRT_isot))$stats[5]
cut_cor_isot <<- (boxplot.stats(c(anaA$int_cor_isot))$stats[1])
if(!is.na(cut_delRT_isot) & !is.na(cut_cor_isot)){
links_profiles_pos <<- enviMass::ZQEZK7587Y(
links_profiles = links_profiles_pos,
profileList = profileList_pos,
cut_delRT_isot = cut_delRT_isot,
cut_cor_isot = cut_cor_isot,
cut_frac_iso = .85
)
}else{
cat("\n No isotopologue linkage filtering feasible")
links_profiles_pos <<- enviMass::ZQEZK7587Y(
links_profiles = links_profiles_pos,
profileList = profileList_pos,
cut_delRT_isot = Inf,
cut_cor_isot = -Inf,
cut_frac_iso = 0
)
}
cut_delRT_adduc <<- boxplot.stats(c(anaA$delRT_adduc))$stats[5]
if(!is.na(cut_delRT_adduc)){
links_profiles_pos <<- enviMass::JZIWC9572W(
links_profiles = links_profiles_pos,
profileList = profileList_pos,
cut_delRT_adduc = cut_delRT_adduc,
cut_frac_adduc = .85
)
}else{
cat("\n No adduct linkage filtering feasible")
links_profiles_pos <<- enviMass::JZIWC9572W(
links_profiles = links_profiles_pos,
profileList = profileList_pos,
cut_delRT_adduc = Inf,
cut_frac_adduc = 0
)
}
anaB <- enviMass::FDGDT7166D(
links_profiles = links_profiles_pos,
profileList = profileList_pos,
min_count = .4,
for_which = logfile$parameters$for_which_profcomp_pos
)
use_EIC <- c(anaB$EIC_cor_isot, anaB$EIC_cor_adduc)
cut_EIC <<- (boxplot.stats(use_EIC)$stats[1])
cut_delRT_EIC <<- max(cut_delRT_isot,cut_delRT_adduc)
if(!is.na(cut_EIC) & !is.na(cut_delRT_EIC)){
links_profiles_pos <<- enviMass::UMNYY9939J(
links_profiles = links_profiles_pos,
profileList = profileList_pos,
cut_EIC = cut_EIC,
cut_frac_EIC = .9,
cut_delRT_EIC = cut_delRT_EIC
)
}else{
cat("\n No EIC linkage filtering feasible")
links_profiles_pos <<- enviMass::UMNYY9939J(
links_profiles = links_profiles_pos,
profileList = profileList_pos,
cut_EIC = 0,
cut_frac_EIC = 0,
cut_delRT_EIC = Inf
)
}
for(n in 1:length(links_profiles_pos)){
is_empty <- enviMass::XMSFU1875I(links_profiles_pos, at_entry = n)
if(is_empty){
links_profiles_pos[[n]] <<- NA
profileList_pos[["index_prof"]][as.numeric(names(links_profiles_pos)[n]),"links"] <<- 0
}
}
}else{
links_profiles_pos <<- enviMass::ZQEZK7587Y(
links_profiles = links_profiles_pos,
profileList = profileList_pos,
cut_delRT_isot = Inf,
cut_cor_isot = -Inf,
cut_frac_iso = 0
)
links_profiles_pos <<- enviMass::JZIWC9572W(
links_profiles = links_profiles_pos,
profileList = profileList_pos,
cut_delRT_adduc = Inf,
cut_frac_adduc = 0
)
links_profiles_pos <<- enviMass::UMNYY9939J(
links_profiles = links_profiles_pos,
profileList = profileList_pos,
cut_EIC = 0,
cut_frac_EIC = 0,
cut_delRT_EIC = Inf
)
}
plot_it <- FALSE
plot_what <- "profiles"
with_test <- mute(as.logical(logfile$parameters$test))
along <- order(profileList_pos[["index_prof"]][,"number_peaks_total"], decreasing = TRUE)
if( (logfile$parameters$filter_profcomp_pos == "TRUE") & (!is.na(cut_delRT_EIC)) ){
use_del_RT <<- cut_delRT_EIC
}else{
use_del_RT <<- as.numeric(logfile$parameters$corr_del_RT)
}
if(with_bar){pBar <- txtProgressBar(min = 0, max = length(along), style = 3)}
found <- 0
for(i in 1:length(along)){
if(with_bar){setTxtProgressBar(pBar, i, title = NULL, label = NULL)}
if( logfile$parameters$prof_comp_link_only == "TRUE" ){
prof_isot_IDs <- enviMass::SRXRJ5620R(
profileList = profileList_pos,
prof_ID = along[i],
links_profiles = links_profiles_pos,
min_peaks = as.numeric(logfile$parameters$corr_min_peaks),
skip_peaks = as.logical(logfile$parameters$corr_skip_peaks),
min_cor = as.numeric(logfile$parameters$comp_corr),
with_test = with_test,
only_direct = FALSE,
del_RT = use_del_RT,
omit_profiles = FALSE
)
prof_adduct_IDs <- c()
for(j in 1:length(prof_isot_IDs)){
got_adducts <- enviMass::UZHOG9959R(
profileList = profileList_pos,
prof_ID = prof_isot_IDs[j],
links_profiles = links_profiles_pos,
min_peaks = as.numeric(logfile$parameters$corr_min_peaks),
skip_peaks = as.logical(logfile$parameters$corr_skip_peaks),
min_cor = as.numeric(logfile$parameters$comp_corr),
with_test = with_test,
omit_profiles = FALSE
)
prof_adduct_IDs <- c(prof_adduct_IDs, got_adducts)
}
prof_all_IDs <- c(prof_isot_IDs, prof_adduct_IDs)
}else{
prof_all_IDs <- enviMass::CNPJL3793Y(
profileList = profileList_pos,
prof_ID = along[i],
links_profiles = links_profiles_pos,
min_peaks = as.numeric(logfile$parameters$corr_min_peaks),
skip_peaks = as.logical(logfile$parameters$corr_skip_peaks),
min_cor = as.numeric(logfile$parameters$comp_corr),
with_test = with_test,
only_direct = FALSE,
del_RT = use_del_RT,
omit_profiles = FALSE
)
}
if(plot_it){
if(length(prof_all_IDs) == 1) next
enviMass::GIEOK5253D(
profileList = profileList_pos,
prof_IDs = prof_all_IDs,
links_profiles = links_profiles_pos,
what = plot_what,
xlim = FALSE, ylim = FALSE, await_input = TRUE,
skipit = TRUE,
min_peaks = as.numeric(logfile$parameters$corr_min_peaks),
norma = TRUE
)
}
if(with_test){
if(!enviMass::RDQZL3486I(prof_all_IDs, profileList_pos, links_profiles_pos)) stop("Not all interlinked!")
}
if(length(prof_all_IDs) > 1){
found <- (found + 1)
if(with_test){if(prof_all_IDs[1] != along[i]){stop("\n\nDebug_not_first!")}}
at_entry <- profileList_pos[["index_prof"]][prof_all_IDs[1], "links"]
links_profiles_pos[[at_entry]][["group"]] <<- prof_all_IDs[-1]
}
}
if(with_bar){close(pBar)}
save(profileList_pos,file=file.path(as.character(logfile[[1]]),"results","profileList_pos"));
save(links_profiles_pos,file=file.path(as.character(logfile[[1]]),"results","links_profiles_pos"));
if(any(ls()=="profileList_pos")){stop("\n illegal profileList_pos detected _1 in do_components_profiles_pl!")}
if(any(ls()=="links_profiles_pos")){stop("\n illegal links_profiles_pos detected _1 in do_components_profiles_pl!")}
rm(links_profiles_pos, profileList_pos, envir=as.environment(".GlobalEnv"))
}
if(
file.exists(file.path(logfile[[1]],"results","profileList_neg"))
){
if(any(objects(envir=as.environment(".GlobalEnv"))=="profileList_neg")){rm(profileList_neg,envir=as.environment(".GlobalEnv"))}
if(any(objects()=="profileList_neg")){rm(profileList_neg)}
if(any(objects(envir=as.environment(".GlobalEnv"))=="links_peaks_neg")){rm(links_peaks_neg,envir=as.environment(".GlobalEnv"))}
if(any(objects()=="links_peaks_neg")){rm(links_peaks_neg)}
if(any(objects(envir=as.environment(".GlobalEnv"))=="links_profiles_neg")){rm(links_profiles_neg,envir=as.environment(".GlobalEnv"))}
if(any(objects()=="links_profiles_neg")){rm(links_profiles_neg)}
load(file.path(as.character(logfile[[1]]),"results","profileList_neg"),envir=as.environment(".GlobalEnv"));
load(file.path(as.character(logfile[[1]]),"results","links_peaks_neg"),envir=as.environment(".GlobalEnv"));
assign("links_profiles_neg", list(), envir = as.environment(".GlobalEnv"))
measurements <- read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
peaks <- profileList_neg[["peaks"]][,c("sampleIDs","peakIDs","profileIDs","RT")]
ord <- order(peaks[,"sampleIDs"],peaks[,"peakIDs"],peaks[,"profileIDs"],decreasing=FALSE)
peaks <- peaks[ord,]
use_entries_profiles <- enviMass::VOLZN1972A(links_profiles_neg)
profileList_neg[["index_prof"]][,"links"] <<- 0
with_bar <- TRUE
if(
(
(logfile$workflow[names(logfile$workflow)=="subtr"]=="no") |
(
((logfile$workflow[names(logfile$workflow)=="IS_screen"]=="yes") &  (logfile$parameters$subtr_IS!="yes")) |
((logfile$workflow[names(logfile$workflow)=="target_screen"]=="yes") & (logfile$parameters$subtr_target!="yes"))
)
) &
(length(links_peaks_neg)>0)
){
cat("\n Annotation of screening results to profiles")
if(with_bar){pBar <- txtProgressBar(min = 0, max = dim(profileList_neg[["index_prof"]])[1], style = 3)}
for(i in 1:dim(profileList_neg[["index_prof"]])[1]){
if(with_bar){setTxtProgressBar(pBar, i, title = NULL, label = NULL)}
if(
any(profileList_neg[["peaks"]][
(profileList_neg[["index_prof"]][i,"start_ID"]:profileList_neg[["index_prof"]][i,"end_ID"]),"links"
]!=0)
){
if( profileList_neg[["index_prof"]][i,"links"]==0 ){
if(length(use_entries_profiles)>0){
at_entry<-use_entries_profiles[1]
use_entries<-use_entries_profiles[-1]
}else{
at_entry<-(length(links_profiles_neg)+1)
}
links_profiles_neg[[at_entry]]<<-enviMass::RHHPS2808T(profileList_neg[["index_prof"]][i,"number_peaks_total"][[1]])
names(links_profiles_neg)[at_entry]<<-as.character(i)
profileList_neg[["index_prof"]][i,"links"]<<-at_entry
}else{
at_entry<-profileList_neg[["index_prof"]][i,"links"]
}
for(j in (profileList_neg[["index_prof"]][i,"start_ID"]:profileList_neg[["index_prof"]][i,"end_ID"])){
if(profileList_neg[["peaks"]][j,"links"]!=0){
if(	length(links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[2]]) > 0 ){
for(k in 1:length(links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[2]])){
if( length(links_profiles_neg[[at_entry]][[2]])==0 ){
links_profiles_neg[[at_entry]][[2]]<<-
data.frame(
links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[2]][[k]][[1]],
1,
links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[2]][[k]][[2]],
stringsAsFactors = FALSE
)
names(links_profiles_neg[[at_entry]][[2]])<<-c("Compound","Counts","max_score")
}else{
at<-which(links_profiles_neg[[at_entry]][[2]][,1] == links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[2]][[k]][[1]])
if(length(at)>0){
links_profiles_neg[[at_entry]][[2]][at,2]<<-(links_profiles_neg[[at_entry]][[2]][at,2]+1)
if(links_profiles_neg[[at_entry]][[2]][at,2] < links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[2]][[k]][[2]] ){
links_profiles_neg[[at_entry]][[2]][at,2] <<- links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[2]][[k]][[2]]
}
}else{
links_profiles_neg[[at_entry]][[2]]<<-data.frame(
c(
links_profiles_neg[[at_entry]][[2]][,1],
links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[2]][[k]][[1]]
),
c(links_profiles_neg[[at_entry]][[2]][,2],1),
c(
links_profiles_neg[[at_entry]][[2]][,3],
links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[2]][[k]][[2]]
),
stringsAsFactors = FALSE
)
names(links_profiles_neg[[at_entry]][[2]])<<-c("Compound","Counts","max_score")
}
}
}
}
if(	length(links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[1]])>0 ){
for(k in 1:length(links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[1]])){
if( length(links_profiles_neg[[at_entry]][[1]])==0 ){
links_profiles_neg[[at_entry]][[1]]<<-
data.frame(
links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[1]][[k]][[1]],
1,
links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[1]][[k]][[2]],
stringsAsFactors = FALSE
)
names(links_profiles_neg[[at_entry]][[1]])<<-c("Compound","Counts","max_score")
}else{
at<-which(links_profiles_neg[[at_entry]][[1]][,1] == links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[1]][[k]][[1]])
if(length(at)>0){
links_profiles_neg[[at_entry]][[1]][at,2]<<-(links_profiles_neg[[at_entry]][[1]][at,2]+1)
if(links_profiles_neg[[at_entry]][[1]][at,2] < links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[1]][[k]][[2]] ){
links_profiles_neg[[at_entry]][[1]][at,2] <<- links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[1]][[k]][[2]]
}
}else{
links_profiles_neg[[at_entry]][[1]]<<-data.frame(
c(
links_profiles_neg[[at_entry]][[1]][,1],
links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[1]][[k]][[1]]
),
c(links_profiles_neg[[at_entry]][[1]][,2],1),
c(
links_profiles_neg[[at_entry]][[1]][,3],
links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[1]][[k]][[2]]
),
stringsAsFactors = FALSE
)
names(links_profiles_neg[[at_entry]][[1]])<<-c("Compound","Counts","max_score")
}
}
}
}
}
}
}
}
if(with_bar){close(pBar)}
cat(" done.")
}
if(logfile$workflow[names(logfile$workflow) == "EIC_correlation"] == "yes"){
cat("\n Retrieving EIC correlation links ")
forIDs <- profileList_neg[["sampleID"]]
for_files <- list.files(file.path(logfile[[1]],"results","componentization","EIC_corr"))
keep <- match(forIDs,for_files)
if(any(is.na(keep))){cat("\n Just note: not all files found in profiles have EIC correlation results (1). \n")}
TRUE_IDs <- (measurements$ID[measurements$EIC_correlation=="TRUE"])
keep2 <- match(forIDs, TRUE_IDs)
forIDs <- forIDs[!is.na(keep) & !is.na(keep2)]
if(logfile$parameters$dofile_latest_profcomp == "TRUE"){
atPOSIX<-profileList_neg[["datetime"]];
matchID<-profileList_neg[["sampleID"]];
atdate<-c();attime<-c();
for(i in 1:length(atPOSIX)){
atdate <- c(atdate, strsplit(atPOSIX[i]," ")[[1]][1]);
attime <- c(attime, strsplit(atPOSIX[i]," ")[[1]][2]);
}
attime <- as.difftime(attime);
atdate <- as.Date(atdate, tz="GMT");
ord <- order(as.numeric(atdate), as.numeric(attime), matchID, decreasing=TRUE);
matchID <- matchID[ord];
forIDs <- forIDs[match(forIDs, matchID)]
if(length(forIDs) > as.numeric(logfile$parameters$numfile_latest_profcomp)){
forIDs <- forIDs[1:as.numeric(logfile$parameters$numfile_latest_profcomp)]
}
forIDs <- forIDs[!is.na(forIDs)]
}
forIDs <- as.numeric(forIDs)
not_found1 <- 0; inserted1 <- 0
if(length(forIDs) > 0){
if(with_bar){pBar <- txtProgressBar(min = 0, max = length(forIDs), style = 3)}
for(i in 1:length(forIDs)){
if(with_bar){setTxtProgressBar(pBar, i, title = NULL, label = NULL)}
load(file = file.path(logfile[[1]], "results", "componentization", "EIC_corr", as.character(forIDs[i])))
if(length(EIC_pairs[,1]) == 0){next}
get1 <- cbind(
rep(forIDs[i], length(EIC_pairs[,1])), EIC_pairs[,1]
)
found1<-enviMass::GSZEL3799Y(get1,peaks[,c("sampleIDs","peakIDs")],row_order=FALSE,column_order_a=FALSE,column_order_b=FALSE,get_index=TRUE)
get2<-cbind(
rep(forIDs[i],length(EIC_pairs[,2])),EIC_pairs[,2]
)
found2<-enviMass::GSZEL3799Y(get2,peaks[,c("sampleIDs","peakIDs")],row_order=FALSE,column_order_a=TRUE,column_order_b=FALSE,get_index=TRUE)
for(j in 1:length(found1)){
if(found1[j]==0){not_found1<-(not_found1+1); next}
if(found2[j]==0){not_found1<-(not_found1+1); next}
inserted1<-(inserted1+1);
prof1<-peaks[found1[j],"profileIDs"][[1]]
prof2<-peaks[found2[j],"profileIDs"][[1]]
if(FALSE){
cat("\n");
cat(peaks[found1[j],"RT"]);cat(" - ")
cat(peaks[found2[j],"RT"])
}
if(profileList_neg[["index_prof"]][prof1,"profile_ID"]!=prof1){stop("\nComponentization: debug me, _1!")}
if(profileList_neg[["index_prof"]][prof1,"links"]==0){
if(length(use_entries_profiles)>0){
at_entry_1<-use_entries_profiles[1]
use_entries<-use_entries_profiles[-1]
}else{
at_entry_1<-(length(links_profiles_neg)+1)
}
links_profiles_neg[[at_entry_1]]<<-enviMass::RHHPS2808T(profileList_neg[["index_prof"]][prof1,"number_peaks_total"][[1]])
names(links_profiles_neg)[at_entry_1]<<-as.character(prof1)
profileList_neg[["index_prof"]][prof1,"links"]<<-at_entry_1
}else{
at_entry_1<-profileList_neg[["index_prof"]][prof1,"links"]
}
here1<-which(links_profiles_neg[[at_entry_1]][["EIC"]][,"linked profile"]==prof2)
if(length(here1)==0){
links_profiles_neg[[at_entry_1]][["EIC"]]<<-rbind(
links_profiles_neg[[at_entry_1]][["EIC"]], c(prof2,1,0,0,0,0,1,NA)
)
here1<-dim(links_profiles_neg[[at_entry_1]][["EIC"]])[1]
}else{
links_profiles_neg[[at_entry_1]][["EIC"]][here1,"link counts"]<<-(links_profiles_neg[[at_entry_1]][["EIC"]][here1,"link counts"]+1)
}
if(profileList_neg[["index_prof"]][prof2,"profile_ID"]!=prof2){stop("\nCross-profile componentization: debug me, EIC_1!")}
if(profileList_neg[["index_prof"]][prof2,"links"]==0){
if(length(use_entries_profiles)>0){
at_entry_2<-use_entries_profiles[1]
use_entries<-use_entries_profiles[-1]
}else{
at_entry_2<-(length(links_profiles_neg)+1)
}
links_profiles_neg[[at_entry_2]]<<-enviMass::RHHPS2808T(profileList_neg[["index_prof"]][prof2,"number_peaks_total"][[1]])
names(links_profiles_neg)[at_entry_2]<<-as.character(prof2)
profileList_neg[["index_prof"]][prof2,"links"]<<-at_entry_2
}else{
at_entry_2<-profileList_neg[["index_prof"]][prof2,"links"]
}
here2<-which(links_profiles_neg[[at_entry_2]][["EIC"]][,"linked profile"]==prof1)
if(length(here2)==0){
links_profiles_neg[[at_entry_2]][["EIC"]]<<-rbind(
links_profiles_neg[[at_entry_2]][["EIC"]], c(prof1,1,0,0,0,0,0,NA)
)
here2<-dim(links_profiles_neg[[at_entry_2]][["EIC"]])[1]
}else{
links_profiles_neg[[at_entry_2]][["EIC"]][here2,"link counts"]<<-(links_profiles_neg[[at_entry_2]][["EIC"]][here2,"link counts"]+1)
}
if(links_profiles_neg[[at_entry_1]][["EIC"]][here1,"use"] == 1){
if(length(links_profiles_neg[[at_entry_1]]$EIC_cor) < here1){
links_profiles_neg[[at_entry_1]]$EIC_cor[[here1]] <<- EIC_pairs[j,4]
}else{
if(is.null(links_profiles_neg[[at_entry_1]]$EIC_cor[[here1]])){
links_profiles_neg[[at_entry_1]]$EIC_cor[[here1]] <<- EIC_pairs[j,4]
}else{
links_profiles_neg[[at_entry_1]]$EIC_cor[[here1]] <<- c(links_profiles_neg[[at_entry_1]]$EIC_cor[[here1]],EIC_pairs[j,4])
}
}
}else{
if(length(links_profiles_neg[[at_entry_2]]$EIC_cor) < here2){
links_profiles_neg[[at_entry_2]]$EIC_cor[[here2] ]<<- EIC_pairs[j,4]
}else{
if(is.null(links_profiles_neg[[at_entry_2]]$EIC_cor[[here2]])){
links_profiles_neg[[at_entry_2]]$EIC_cor[[here2]] <<- EIC_pairs[j,4]
}else{
links_profiles_neg[[at_entry_2]]$EIC_cor[[here2]] <<- c(links_profiles_neg[[at_entry_2]]$EIC_cor[[here2]],EIC_pairs[j,4])
}
}
}
}
}
if(with_bar){close(pBar)}
if(any(objects()=="EIC_pairs")){rm(EIC_pairs)}
cat(" done.")
}
}
if(logfile$workflow[names(logfile$workflow) == "isotopologues"] == "yes"){
forIDs <- profileList_neg[["sampleID"]]
for_files <- list.files(file.path(logfile[[1]], "results", "componentization", "isotopologues"))
keep <- match(forIDs, for_files)
if(any(is.na(keep))){cat("\n Just note: not all files found in profiles have isotopologue links available. \n")}
TRUE_IDs <- (measurements$ID[measurements$isotopologues == "TRUE"])
keep2 <- match(forIDs, TRUE_IDs)
forIDs <- forIDs[!is.na(keep) & !is.na(keep2)]
if(logfile$parameters$dofile_latest_profcomp == "TRUE"){
atPOSIX <- profileList_neg[["datetime"]];
matchID <- profileList_neg[["sampleID"]];
atdate <- c(); attime <- c();
for(i in 1:length(atPOSIX)){
atdate <- c(atdate, strsplit(atPOSIX[i]," ")[[1]][1]);
attime <- c(attime, strsplit(atPOSIX[i]," ")[[1]][2]);
}
attime <- as.difftime(attime);
atdate <- as.Date(atdate, tz = "GMT");
ord <- order(as.numeric(atdate), as.numeric(attime), matchID, decreasing = TRUE);
matchID <- matchID[ord];
forIDs <- forIDs[match(forIDs, matchID)]
if(length(forIDs) > as.numeric(logfile$parameters$numfile_latest_profcomp)){
forIDs <- forIDs[1:as.numeric(logfile$parameters$numfile_latest_profcomp)]
}
forIDs <- forIDs[!is.na(forIDs)]
}
forIDs <- as.numeric(forIDs)
not_found3 <- 0; inserted3 <- 0
if(length(forIDs) > 0){
cat("\n Retrieving isotopologue links ")
if(with_bar){pBar <- txtProgressBar(min = 0, max = length(forIDs), style = 3)}
for(i in 1:length(forIDs)){
if(with_bar){setTxtProgressBar(pBar, i, title = NULL, label = NULL)}
load(file = file.path(logfile[[1]], "results", "componentization", "isotopologues", as.character(forIDs[i])))
if(length(Isot_pairs[,1]) == 0){next}
get1 <- cbind(
rep(forIDs[i], length(Isot_pairs[,1])), Isot_pairs[,1]
)
found1 <- enviMass::GSZEL3799Y(get1, peaks[,c("sampleIDs", "peakIDs")], row_order = FALSE, column_order_a = FALSE, column_order_b = FALSE, get_index = TRUE)
get2 <- cbind(
rep(forIDs[i], length(Isot_pairs[,2])), Isot_pairs[,2]
)
found2 <- enviMass::GSZEL3799Y(get2, peaks[,c("sampleIDs", "peakIDs")], row_order = FALSE, column_order_a = TRUE, column_order_b = FALSE, get_index = TRUE)
for(j in 1:length(found1)){
if(found1[j] == 0){not_found3 <- (not_found3 + 1); next}
if(found2[j] == 0){not_found3 <- (not_found3 + 1); next}
inserted3 <- (inserted3 + 1);
if(FALSE){
cat("\n");
cat(peaks[found1[j], "RT"]); cat(" - ");
cat(peaks[found2[j], "RT"])
}
prof1 <- peaks[found1[j], "profileIDs"][[1]]
prof2 <- peaks[found2[j], "profileIDs"][[1]]
if(profileList_neg[["index_prof"]][prof1, "profile_ID"] != prof1){stop("\nComponentization: debug me, _1!")}
if(profileList_neg[["index_prof"]][prof1, "links"] == 0){
if(length(use_entries_profiles) > 0){
at_entry_1 <- use_entries_profiles[1]
use_entries <- use_entries_profiles[-1]
}else{
at_entry_1 <- (length(links_profiles_neg) + 1)
}
links_profiles_neg[[at_entry_1]] <<- enviMass::RHHPS2808T(profileList_neg[["index_prof"]][prof1, "number_peaks_total"][[1]])
names(links_profiles_neg)[at_entry_1] <<- as.character(prof1)
profileList_neg[["index_prof"]][prof1, "links"] <<- at_entry_1
}else{
at_entry_1<-profileList_neg[["index_prof"]][prof1, "links"]
}
here1 <- which(links_profiles_neg[[at_entry_1]][["isot"]][, "linked profile"] == prof2)
if(length(here1)==0){
links_profiles_neg[[at_entry_1]][["isot"]] <<- rbind(
links_profiles_neg[[at_entry_1]][["isot"]], c(prof2,1,0,1,NA)
)
here1 <- dim(links_profiles_neg[[at_entry_1]][["isot"]])[1]
is_new1 <- TRUE
}else{
links_profiles_neg[[at_entry_1]][["isot"]][here1, "link counts"] <<- (links_profiles_neg[[at_entry_1]][["isot"]][here1,"link counts"] + 1)
is_new1 <- FALSE
}
if(profileList_neg[["index_prof"]][prof2, "profile_ID"] != prof2){stop("\nComponentization: debug me, _1!")}
if(profileList_neg[["index_prof"]][prof2, "links"] == 0){
if(length(use_entries_profiles) > 0){
at_entry_2 <- use_entries_profiles[1]
use_entries <- use_entries_profiles[-1]
}else{
at_entry_2 <- (length(links_profiles_neg)+1)
}
links_profiles_neg[[at_entry_2]] <<- enviMass::RHHPS2808T(profileList_neg[["index_prof"]][prof2,"number_peaks_total"][[1]])
names(links_profiles_neg)[at_entry_2] <<- as.character(prof2)
profileList_neg[["index_prof"]][prof2,"links"] <<- at_entry_2
}else{
at_entry_2 <- profileList_neg[["index_prof"]][prof2,"links"]
}
here2<-which(links_profiles_neg[[at_entry_2]][["isot"]][,"linked profile"] == prof1)
if(length(here2) == 0){
links_profiles_neg[[at_entry_2]][["isot"]] <<- rbind(
links_profiles_neg[[at_entry_2]][["isot"]], c(prof1,1,0,0,NA)
)
here2 <- dim(links_profiles_neg[[at_entry_2]][["isot"]])[1]
is_new2 <- TRUE
}else{
links_profiles_neg[[at_entry_2]][["isot"]][here2,"link counts"] <<- (links_profiles_neg[[at_entry_2]][["isot"]][here2,"link counts"]+1)
is_new2 <- FALSE
}
if(is_new1 | is_new2){
these <- profileList_neg[["peaks"]][
profileList_neg[["index_prof"]][prof1,"start_ID"]:profileList_neg[["index_prof"]][prof1,"end_ID"]
,"sampleIDs"]
these <- these[!is.na(match(these,forIDs))]
those <- profileList_neg[["peaks"]][
profileList_neg[["index_prof"]][prof2,"start_ID"]:profileList_neg[["index_prof"]][prof2,"end_ID"]
,"sampleIDs"]
those <- those[!is.na(match(these,forIDs))]
matched <- match(these,those)
not_NA <- sum(!is.na(matched))
links_profiles_neg[[at_entry_1]]$isot[here1,"ref_1"] <<- not_NA
links_profiles_neg[[at_entry_2]]$isot[here2,"ref_1"] <<- not_NA
}
if(any(links_profiles_neg[[at_entry_1]]$isot[,"ref_1"] == 0)){stop("\n DEBUG ME ! FOUND_1")}
if(any(links_profiles_neg[[at_entry_2]]$isot[,"ref_1"] == 0)){stop("\n DEBUG ME ! FOUND_2")}
}
}
if(with_bar){close(pBar)}
if(any(objects()=="Isot_pairs")){rm(Isot_pairs)}
cat(" done.")
}
}
if(logfile$workflow[names(logfile$workflow)=="adducts"]=="yes"){
forIDs<-profileList_neg[["sampleID"]]
for_files<-list.files(file.path(logfile[[1]],"results","componentization","adducts"))
keep<-match(forIDs,for_files)
if(any(is.na(keep))){cat("\n Just note: not all files found in profiles have adduct links available. \n")}
TRUE_IDs <- (measurements$ID[measurements$adducts=="TRUE"])
keep2 <- match(forIDs, TRUE_IDs)
forIDs<-forIDs[!is.na(keep) & !is.na(keep2)]
if(logfile$parameters$dofile_latest_profcomp=="TRUE"){
atPOSIX<-profileList_neg[["datetime"]];
matchID<-profileList_neg[["sampleID"]];
atdate<-c();attime<-c();
for(i in 1:length(atPOSIX)){
atdate<-c(atdate, strsplit(atPOSIX[i]," ")[[1]][1]);
attime<-c(attime, strsplit(atPOSIX[i]," ")[[1]][2]);
}
attime<-as.difftime(attime);
atdate<-as.Date(atdate, tz = "GMT");
ord<-order(as.numeric(atdate),as.numeric(attime),matchID,decreasing=TRUE);
matchID<-matchID[ord];
forIDs<-forIDs[match(forIDs,matchID)]
if(length(forIDs)>as.numeric(logfile$parameters$numfile_latest_profcomp)){
forIDs<-forIDs[1:as.numeric(logfile$parameters$numfile_latest_profcomp)]
}
forIDs <- forIDs[!is.na(forIDs)]
}
forIDs <- as.numeric(forIDs)
not_found4 <- 0; inserted4 <- 0
if(length(forIDs) > 0){
cat("\n Retrieving adduct links ")
if(with_bar){pBar <- txtProgressBar(min = 0, max = length(forIDs), style = 3)}
for(i in 1:length(forIDs)){
if(with_bar){setTxtProgressBar(pBar, i, title = NULL, label = NULL)}
load(file=file.path(logfile[[1]],"results","componentization","adducts",as.character(forIDs[i])))
if(length(Adduct_pairs[,1])==0){next}
get1<-cbind(
rep(forIDs[i],length(Adduct_pairs[,1])),Adduct_pairs[,1]
)
found1<-enviMass::GSZEL3799Y(get1,peaks[,c("sampleIDs","peakIDs")],row_order=FALSE,column_order_a=FALSE,column_order_b=FALSE,get_index=TRUE)
get2<-cbind(
rep(forIDs[i],length(Adduct_pairs[,2])),Adduct_pairs[,2]
)
found2<-enviMass::GSZEL3799Y(get2,peaks[,c("sampleIDs","peakIDs")],row_order=FALSE,column_order_a=TRUE,column_order_b=FALSE,get_index=TRUE)
for(j in 1:length(found1)){
if(found1[j]==0){not_found4<-(not_found4+1);next}
if(found2[j]==0){not_found4<-(not_found4+1);next}
inserted4<-(inserted4+1);
if(FALSE){
cat("\n");
cat(peaks[found1[j],"RT"]);cat(" - ")
cat(peaks[found2[j],"RT"])
}
prof1<-peaks[found1[j],"profileIDs"][[1]]
prof2<-peaks[found2[j],"profileIDs"][[1]]
if(profileList_neg[["index_prof"]][prof1,"profile_ID"]!=prof1){stop("\nComponentization: debug me, _1!")}
if(profileList_neg[["index_prof"]][prof1,"links"]==0){
if(length(use_entries_profiles)>0){
at_entry_1<-use_entries_profiles[1]
use_entries<-use_entries_profiles[-1]
}else{
at_entry_1<-(length(links_profiles_neg)+1)
}
links_profiles_neg[[at_entry_1]] <<- enviMass::RHHPS2808T(profileList_neg[["index_prof"]][prof1,"number_peaks_total"][[1]])
names(links_profiles_neg)[at_entry_1] <<- as.character(prof1)
profileList_neg[["index_prof"]][prof1,"links"] <<- at_entry_1
}else{
at_entry_1<-profileList_neg[["index_prof"]][prof1,"links"]
}
here1<-which(links_profiles_neg[[at_entry_1]][["adduc"]][,"linked profile"]==prof2)
if(length(here1)==0){
links_profiles_neg[[at_entry_1]][["adduc"]] <<- rbind(
links_profiles_neg[[at_entry_1]][["adduc"]], c(prof2,1,0,1,NA)
)
here1<-dim(links_profiles_neg[[at_entry_1]][["adduc"]])[1]
is_new1 <- TRUE
}else{
links_profiles_neg[[at_entry_1]][["adduc"]][here1,"link counts"] <<- (links_profiles_neg[[at_entry_1]][["adduc"]][here1,"link counts"]+1)
is_new1 <- FALSE
}
if(profileList_neg[["index_prof"]][prof2,"profile_ID"]!=prof2){stop("\nComponentization: debug me, _1!")}
if(profileList_neg[["index_prof"]][prof2,"links"]==0){
if(length(use_entries_profiles)>0){
at_entry_2<-use_entries_profiles[1]
use_entries<-use_entries_profiles[-1]
}else{
at_entry_2<-(length(links_profiles_neg)+1)
}
links_profiles_neg[[at_entry_2]] <<- enviMass::RHHPS2808T(profileList_neg[["index_prof"]][prof2,"number_peaks_total"][[1]])
names(links_profiles_neg)[at_entry_2] <<- as.character(prof2)
profileList_neg[["index_prof"]][prof2,"links"] <<- at_entry_2
}else{
at_entry_2<-profileList_neg[["index_prof"]][prof2,"links"]
}
here2<-which(links_profiles_neg[[at_entry_2]][["adduc"]][,"linked profile"]==prof1)
if(length(here2)==0){
links_profiles_neg[[at_entry_2]][["adduc"]]<<-rbind(
links_profiles_neg[[at_entry_2]][["adduc"]], c(prof1,1,0,0,NA)
)
here2<-dim(links_profiles_neg[[at_entry_2]][["adduc"]])[1]
is_new2 <- TRUE
}else{
links_profiles_neg[[at_entry_2]][["adduc"]][here2,"link counts"]<<-(links_profiles_neg[[at_entry_2]][["adduc"]][here2,"link counts"]+1)
is_new2 <- FALSE
}
if(is_new1 | is_new2){
these<-profileList_neg[["peaks"]][
profileList_neg[["index_prof"]][prof1,"start_ID"]:profileList_neg[["index_prof"]][prof1,"end_ID"]
,"sampleIDs"]
these <- these[!is.na(match(these,forIDs))]
those<-profileList_neg[["peaks"]][
profileList_neg[["index_prof"]][prof2,"start_ID"]:profileList_neg[["index_prof"]][prof2,"end_ID"]
,"sampleIDs"]
those <- those[!is.na(match(those,forIDs))]
matched <- match(these,those)
not_NA <- sum(!is.na(matched))
links_profiles_neg[[at_entry_1]]$adduc[here1,"ref_1"]<<-not_NA
links_profiles_neg[[at_entry_2]]$adduc[here2,"ref_1"]<<-not_NA
}
if(any(links_profiles_neg[[at_entry_1]]$adduc[,"ref_1"]==0)){stop("\n DEBUG ME ! FOUND_3")}
if(any(links_profiles_neg[[at_entry_2]]$adduc[,"ref_1"]==0)){stop("\n DEBUG ME ! FOUND_4")}
}
}
if(with_bar){close(pBar)}
if(any(objects()=="Adduct_pairs")){rm(Adduct_pairs)}
cat(" done.")
}
}
if(logfile$workflow[names(logfile$workflow)=="homologues"]=="yes"){
forIDs<-profileList_neg[["sampleID"]]
for_files<-list.files(file.path(logfile[[1]],"results","componentization","homologues"))
keep <- match(forIDs, for_files)
if(any(is.na(keep))){cat("\n Just note: not all files found in profiles have homologue links available. \n")}
TRUE_IDs <- (measurements$ID[measurements$homologues == "TRUE"])
keep2<-match(forIDs, TRUE_IDs)
forIDs <- forIDs[!is.na(keep) & !is.na(keep2)]
if(logfile$parameters$dofile_latest_profcomp == "TRUE"){
atPOSIX<-profileList_neg[["datetime"]];
matchID<-profileList_neg[["sampleID"]];
atdate<-c();attime<-c();
for(i in 1:length(atPOSIX)){
atdate<-c(atdate, strsplit(atPOSIX[i]," ")[[1]][1]);
attime<-c(attime, strsplit(atPOSIX[i]," ")[[1]][2]);
}
attime<-as.difftime(attime);
atdate<-as.Date(atdate, tz="GMT");
ord<-order(as.numeric(atdate),as.numeric(attime),matchID,decreasing=TRUE);
matchID<-matchID[ord];
forIDs<-forIDs[match(forIDs,matchID)]
if(length(forIDs)>as.numeric(logfile$parameters$numfile_latest_profcomp)){
forIDs<-forIDs[1:as.numeric(logfile$parameters$numfile_latest_profcomp)]
}
forIDs <- forIDs[!is.na(forIDs)]
}
not_found5<-0;inserted5<-0
if(length(forIDs)>0){
cat("\n Retrieving homologue links ")
if(with_bar){pBar <- txtProgressBar(min = 0, max = length(forIDs), style = 3)}
for(i in 1:length(forIDs)){
if(with_bar){setTxtProgressBar(pBar, i, title = NULL, label = NULL)}
load(file=file.path(logfile[[1]],"results","componentization","homologues",forIDs[i]))
if(length(Homol_groups[,1])==0){next}
keep_group<-rep(TRUE,max(Homol_groups[,3]))
get1<-cbind(
rep(as.numeric(forIDs[i]),length(Homol_groups[,1])),Homol_groups[,1]
)
found1<-enviMass::GSZEL3799Y(get1,peaks[,c("sampleIDs","peakIDs")],row_order=FALSE,column_order_a=FALSE,column_order_b=FALSE,get_index=TRUE)
keep_group[Homol_groups[found1==0,3]]<-FALSE
if(!any(keep_group)){next}
get2<-cbind(
rep(as.numeric(forIDs[i]),length(Homol_groups[,2])),Homol_groups[,2]
)
found2<-enviMass::GSZEL3799Y(get2,peaks[,c("sampleIDs","peakIDs")],row_order=FALSE,column_order_a=TRUE,column_order_b=FALSE,get_index=TRUE)
if(!any(keep_group)){next}
for(j in 1:length(found1)){
if(found1[j]==0){not_found5<-(not_found5+1);next}
if(found2[j]==0){not_found5<-(not_found5+1);next}
if(Homol_groups[j,3]!=Homol_groups[j,3]){stop("\n Debug componentization for its homologue relations!")}
if(!keep_group[Homol_groups[j,3]]){next}
inserted5<-(inserted5+1);
if(FALSE){
cat("\n");
cat(peaks[found1[j],"RT"]);cat(" - ")
cat(peaks[found2[j],"RT"])
}
prof1<-peaks[found1[j],"profileIDs"][[1]]
prof2<-peaks[found2[j],"profileIDs"][[1]]
if(profileList_neg[["index_prof"]][prof1,"profile_ID"]!=prof1){stop("\nComponentization: debug me, _1!")}
if(profileList_neg[["index_prof"]][prof1,"links"]==0){
if(length(use_entries_profiles)>0){
at_entry<-use_entries_profiles[1]
use_entries<-use_entries_profiles[-1]
}else{
at_entry<-(length(links_profiles_neg)+1)
}
links_profiles_neg[[at_entry]]<<-enviMass::RHHPS2808T(profileList_neg[["index_prof"]][prof1,"number_peaks_total"][[1]])
names(links_profiles_neg)[at_entry]<<-as.character(prof1)
profileList_neg[["index_prof"]][prof1,"links"]<<-at_entry
}else{
at_entry<-profileList_neg[["index_prof"]][prof1,"links"]
}
here<-which(links_profiles_neg[[at_entry]][["homol"]][,"linked profile"]==prof2)
if(length(here)==0){
links_profiles_neg[[at_entry]][["homol"]]<<-rbind(
links_profiles_neg[[at_entry]][["homol"]], c(prof2,1,0)
)
}else{
links_profiles_neg[[at_entry]][["homol"]][here,"link counts"]<<-(links_profiles_neg[[at_entry]][["homol"]][here,"link counts"]+1)
}
if(profileList_neg[["index_prof"]][prof2,"profile_ID"]!=prof2){stop("\nComponentization: debug me, _1!")}
if(profileList_neg[["index_prof"]][prof2,"links"]==0){
if(length(use_entries_profiles)>0){
at_entry<-use_entries_profiles[1]
use_entries<-use_entries_profiles[-1]
}else{
at_entry<-(length(links_profiles_neg)+1)
}
links_profiles_neg[[at_entry]]<<-enviMass::RHHPS2808T(profileList_neg[["index_prof"]][prof2,"number_peaks_total"][[1]])
names(links_profiles_neg)[at_entry]<<-as.character(prof2)
profileList_neg[["index_prof"]][prof2,"links"]<<-at_entry
}else{
at_entry<-profileList_neg[["index_prof"]][prof2,"links"]
}
here<-which(links_profiles_neg[[at_entry]][["homol"]][,"linked profile"]==prof1)
if(length(here)==0){
links_profiles_neg[[at_entry]][["homol"]]<<-rbind(
links_profiles_neg[[at_entry]][["homol"]], c(prof1,1,0)
)
}else{
links_profiles_neg[[at_entry]][["homol"]][here,"link counts"]<<-(links_profiles_neg[[at_entry]][["homol"]][here,"link counts"]+1)
}
}
}
if(with_bar){close(pBar)}
if(any(objects()=="Homol_groups")){rm(Homol_groups)}
cat(" done.")
}
}
cut_delRT_EIC<<-NA
if(
(logfile$parameters$filter_profcomp_neg == "TRUE") &
( (logfile$workflow[names(logfile$workflow) == "isotopologues"] == "yes") || (logfile$workflow[names(logfile$workflow)=="adducts"]=="yes") )
){
fil1<-enviMass::XFDKJ7277K(
links_profiles = links_profiles_neg,
profileList = profileList_neg,
logfile = logfile,
min_rat=.7,
min_count=.4,
perc =.9,
for_which=logfile$parameters$for_which_profcomp_neg
)
cut_delRT_isot<<-boxplot.stats(c(fil1$delRT_isot))$stats[5]
cut_cor_isot<<-(boxplot.stats(c(fil1$int_cor_isot))$stats[1])
if(!is.na(cut_delRT_isot)&!is.na(cut_cor_isot)){
links_profiles_neg <<- enviMass::ZQEZK7587Y(
links_profiles = links_profiles_neg,
profileList = profileList_neg,
cut_delRT_isot = cut_delRT_isot,
cut_cor_isot = cut_cor_isot,
cut_frac_iso = .85
)
}else{
cat("\n No isotopologue linkage filtering feasible")
links_profiles_neg <<- enviMass::ZQEZK7587Y(
links_profiles = links_profiles_neg,
profileList = profileList_neg,
cut_delRT_isot = Inf,
cut_cor_isot = -Inf,
cut_frac_iso = 0
)
}
cut_delRT_adduc<<-boxplot.stats(c(fil1$delRT_adduc))$stats[5]
if(!is.na(cut_delRT_adduc)){
links_profiles_neg <<- enviMass::JZIWC9572W(
links_profiles = links_profiles_neg,
profileList = profileList_neg,
cut_delRT_adduc = cut_delRT_adduc,
cut_frac_adduc = .85
)
}else{
cat("\n No adduct linkage filtering feasible")
links_profiles_neg <<- enviMass::JZIWC9572W(
links_profiles = links_profiles_neg,
profileList = profileList_neg,
cut_delRT_adduc = Inf,
cut_frac_adduc = 0
)
}
fil2<-enviMass::FDGDT7166D(
links_profiles = links_profiles_neg,
profileList = profileList_neg,
min_count=.4,
for_which=logfile$parameters$for_which_profcomp_neg
)
use_EIC<-c(fil2$EIC_cor_isot,fil2$EIC_cor_adduc)
cut_EIC<<-(boxplot.stats(use_EIC)$stats[1])
cut_delRT_EIC<<-max(cut_delRT_isot,cut_delRT_adduc)
if(!is.na(cut_EIC)&!is.na(cut_delRT_EIC)){
links_profiles_neg <<- enviMass::UMNYY9939J(
links_profiles = links_profiles_neg,
profileList = profileList_neg,
cut_EIC = cut_EIC,
cut_frac_EIC = .9,
cut_delRT_EIC = cut_delRT_EIC
)
}else{
cat("\n No EIC linkage filtering feasible")
links_profiles_neg <<- enviMass::UMNYY9939J(
links_profiles = links_profiles_neg,
profileList = profileList_neg,
cut_EIC = 0,
cut_frac_EIC = 0,
cut_delRT_EIC = Inf
)
}
for(n in 1:length(links_profiles_neg)){
is_empty<-enviMass::XMSFU1875I(links_profiles_neg, at_entry = n)
if(is_empty){
links_profiles_neg[[n]] <<- NA
profileList_neg[["index_prof"]][as.numeric(names(links_profiles_neg)[n]),"links"] <<- 0
}
}
}else{
links_profiles_neg <<- enviMass::ZQEZK7587Y(
links_profiles = links_profiles_neg,
profileList = profileList_neg,
cut_delRT_isot = Inf,
cut_cor_isot = -Inf,
cut_frac_iso = 0
)
links_profiles_neg <<- enviMass::JZIWC9572W(
links_profiles = links_profiles_neg,
profileList = profileList_neg,
cut_delRT_adduc = Inf,
cut_frac_adduc = 0
)
links_profiles_neg <<- enviMass::UMNYY9939J(
links_profiles = links_profiles_neg,
profileList = profileList_neg,
cut_EIC = 0,
cut_frac_EIC = 0,
cut_delRT_EIC = Inf
)
}
plot_it <- FALSE
plot_what <- "profiles"
with_test <- mute(as.logical(logfile$parameters$test))
along <- order(profileList_neg[["index_prof"]][,"number_peaks_total"], decreasing = TRUE)
if( (logfile$parameters$filter_profcomp_neg == "TRUE") & (!is.na(cut_delRT_EIC)) ){
use_del_RT <<- cut_delRT_EIC
}else{
use_del_RT <<- as.numeric(logfile$parameters$corr_del_RT)
}
if(with_bar){pBar <- txtProgressBar(min = 0, max = length(along), style = 3)}
for(i in 1:length(along)){
if(with_bar){setTxtProgressBar(pBar, i, title = NULL, label = NULL)}
if(logfile$parameters$prof_comp_link_only=="TRUE"){
prof_isot_IDs<-enviMass::SRXRJ5620R(
profileList=profileList_neg,
prof_ID=along[i],
links_profiles=links_profiles_neg,
min_peaks=as.numeric(logfile$parameters$corr_min_peaks),
skip_peaks=as.logical(logfile$parameters$corr_skip_peaks),
min_cor=as.numeric(logfile$parameters$comp_corr),
with_test=with_test,
only_direct=FALSE,
del_RT=use_del_RT,
omit_profiles=FALSE
)
prof_adduct_IDs<-c()
for(j in 1:length(prof_isot_IDs)){
got_adducts<-enviMass::UZHOG9959R(
profileList=profileList_neg,
prof_ID=prof_isot_IDs[j],
links_profiles=links_profiles_neg,
min_peaks=as.numeric(logfile$parameters$corr_min_peaks),
skip_peaks=as.logical(logfile$parameters$corr_skip_peaks),
min_cor=as.numeric(logfile$parameters$comp_corr),
with_test=with_test,
omit_profiles=FALSE
)
prof_adduct_IDs<-c(prof_adduct_IDs,got_adducts)
}
prof_all_IDs<-c(prof_isot_IDs,prof_adduct_IDs)
}else{
prof_all_IDs<-enviMass::CNPJL3793Y(
profileList=profileList_neg,
prof_ID=along[i],
links_profiles=links_profiles_neg,
min_peaks=as.numeric(logfile$parameters$corr_min_peaks),
skip_peaks=as.logical(logfile$parameters$corr_skip_peaks),
min_cor=as.numeric(logfile$parameters$comp_corr),
with_test=with_test,
only_direct=FALSE,
del_RT=use_del_RT,
omit_profiles=FALSE
)
}
if(plot_it){
if(length(prof_all_IDs)==1) next
enviMass::GIEOK5253D(
profileList=profileList_neg,
prof_IDs=prof_all_IDs,
links_profiles=links_profiles_neg,
what=plot_what,
xlim=FALSE,ylim=FALSE,await_input=TRUE,
skipit=TRUE,
min_peaks=as.numeric(logfile$parameters$corr_min_peaks),
norma=TRUE
)
}
if(with_test){
if(!enviMass::RDQZL3486I(prof_all_IDs,profileList_neg,links_profiles_neg)) stop("Not all interlinked!")
}
if(length(prof_all_IDs)>1){
if(with_test){if(prof_all_IDs[1]!=along[i]){stop("\n\nDebug_not_first!")}}
at_entry<-profileList_neg[["index_prof"]][prof_all_IDs[1],"links"]
links_profiles_neg[[at_entry]][["group"]] <<- prof_all_IDs[-1]
}
}
if(with_bar){close(pBar)}
save(profileList_neg,file=file.path(as.character(logfile[[1]]),"results","profileList_neg"));
save(links_profiles_neg,file=file.path(as.character(logfile[[1]]),"results","links_profiles_neg"));
if(any(ls()=="profileList_neg")){stop("\n illegal profileList_neg detected _1 in do_IS_normaliz.r!")}
if(any(ls()=="links_profiles_neg")){stop("\n illegal links_profiles_neg detected _1 in do_IS_normaliz.r!")}
rm(links_profiles_neg, profileList_neg, envir=as.environment(".GlobalEnv"))
}
