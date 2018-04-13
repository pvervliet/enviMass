if(file.exists(file.path(as.character(logfile[[1]]),"results","links_profiles_pos"))){file.remove(file.path(as.character(logfile[[1]]),"results","links_profiles_pos"))}
if(file.exists(file.path(as.character(logfile[[1]]),"results","links_profiles_neg"))){file.remove(file.path(as.character(logfile[[1]]),"results","links_profiles_neg"))}
if(file.exists(file.path(as.character(logfile[[1]]),"results","int_norm_ISTD_pos"))){file.remove(file.path(as.character(logfile[[1]]),"results","int_norm_ISTD_pos"))}
if(file.exists(file.path(as.character(logfile[[1]]),"results","int_norm_ISTD_neg"))){file.remove(file.path(as.character(logfile[[1]]),"results","int_norm_ISTD_neg"))}
if(any(objects(envir=as.environment(".GlobalEnv"))=="profileList_pos")){rm(profileList_pos,envir=as.environment(".GlobalEnv"))}
if(any(objects()=="profileList_pos")){rm(profileList_pos)}
if(any(objects(envir=as.environment(".GlobalEnv"))=="links_peaks_pos")){rm(links_peaks_pos,envir=as.environment(".GlobalEnv"))}
if(any(objects()=="links_peaks_pos")){rm(links_peaks_pos)}
if(any(objects(envir=as.environment(".GlobalEnv"))=="links_profiles_pos")){rm(links_profiles_pos,envir=as.environment(".GlobalEnv"))}
if(any(objects()=="links_profiles_pos")){rm(links_profiles_pos)}
if(any(objects(envir=as.environment(".GlobalEnv"))=="profileList_neg")){rm(profileList_neg,envir=as.environment(".GlobalEnv"))}
if(any(objects()=="profileList_neg")){rm(profileList_neg)}
if(any(objects(envir=as.environment(".GlobalEnv"))=="links_peaks_neg")){rm(links_peaks_neg,envir=as.environment(".GlobalEnv"))}
if(any(objects()=="links_peaks_neg")){rm(links_peaks_neg)}
if(any(objects(envir=as.environment(".GlobalEnv"))=="links_profiles_neg")){rm(links_profiles_neg,envir=as.environment(".GlobalEnv"))}
if(any(objects()=="links_profiles_neg")){rm(links_profiles_neg)}
if(
file.exists(file.path(logfile[[1]],"results","profileList_pos")) &
(as.logical(logfile$parameters$ISnorm_include_pos) | logfile$workflow[names(logfile$workflow) == "components_profiles"] == "yes")
){
if(any(objects(envir=as.environment(".GlobalEnv"))=="profileList_pos")){rm(profileList_pos,envir=as.environment(".GlobalEnv"))}
if(any(objects()=="profileList_pos")){rm(profileList_pos)}
if(any(objects(envir=as.environment(".GlobalEnv"))=="links_peaks_pos")){rm(links_peaks_pos,envir=as.environment(".GlobalEnv"))}
if(any(objects()=="links_peaks_pos")){rm(links_peaks_pos)}
if(any(objects(envir=as.environment(".GlobalEnv"))=="links_profiles_pos")){rm(links_profiles_pos,envir=as.environment(".GlobalEnv"))}
if(any(objects()=="links_profiles_pos")){rm(links_profiles_pos)}
load(file.path(as.character(logfile[[1]]),"results","profileList_pos"),envir=as.environment(".GlobalEnv"));
load(file.path(as.character(logfile[[1]]),"results","links_peaks_pos"),envir=as.environment(".GlobalEnv"));
assign("links_profiles_pos", list(), envir = as.environment(".GlobalEnv"))
measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
use_entries_profiles<-enviMass::VOLZN1972A(links_profiles_pos)
profileList_pos[["index_prof"]][,"links"] <<- 0
with_bar<-FALSE
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
at_entry <- profileList_pos[["index_prof"]][i,"links"]
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
}
if(!length(links_profiles_pos)){
stop("Certainly not enough internal standard (ISTD) screening matches detected to run the ISTD normalization, positive ionization - remove this workflow step?")
}
if(as.logical(logfile$parameters$ISnorm_include_pos)){
if(as.logical(logfile$parameters$screen_IS_restrict)){
min_count <- min(c(
floor(as.numeric(logfile$parameters$screen_IS_restrict_many) * as.numeric(logfile$parameters$ISnorm_percfiles_pos) / 100),
floor(length(profileList_pos[[4]]) * as.numeric(logfile$parameters$ISnorm_percfiles_pos) / 100)
))
}else{
min_count <- floor(length(profileList_pos[[4]]) * as.numeric(logfile$parameters$ISnorm_percfiles_pos) / 100);
}
lis_delint_IS <- list()
lis_median_IS <- list()
samIDs <- profileList_pos[["sampleID"]]
atPOSIX <- profileList_pos[["datetime"]]
type <- profileList_pos[["type"]]
ord <- order(atPOSIX, type, decreasing = FALSE)
samIDs <- samIDs[ord]
type <- type[ord]
atPOSIX <- atPOSIX[ord]
for(p in 1:length(profileList_pos[["sampleID"]])){
lis_delint_IS[[p]] <- numeric(0)
lis_median_IS[[p]] <- numeric(0)
}
names(lis_delint_IS) <- samIDs
names(lis_median_IS) <- samIDs
unused_profile <- rep(TRUE, dim(profileList_pos[["index_prof"]])[1])
for(i in 1:length(links_profiles_pos)){
if(!length(links_profiles_pos[[i]]$IS)) next
these <- which((links_profiles_pos[[i]]$IS$Counts >= min_count) & (links_profiles_pos[[i]]$IS$max_score >= as.numeric(logfile$parameters$ISnorm_score_pos) ) )
if(!length(these)) next
at_profile <- as.numeric(names(links_profiles_pos)[i])
unused_profile[at_profile] <- FALSE
if(profileList_pos[["index_prof"]][at_profile,"number_peaks_total"] < min_count) next
for_peaks <- profileList_pos[["index_prof"]][at_profile,"start_ID"]:profileList_pos[["index_prof"]][at_profile,"end_ID"]
median_intensity <- median( log10(profileList_pos[["peaks"]][for_peaks, "intensity"]) )
for(j in for_peaks){
at_file <- match(profileList_pos[["peaks"]][j,"sampleIDs"], samIDs)
at_int <- log10(profileList_pos[["peaks"]][j,"intensity"][[1]])
lis_delint_IS[[at_file]] <- c(lis_delint_IS[[at_file]], (at_int - median_intensity))
lis_median_IS[[at_file]] <- c(lis_median_IS[[at_file]], median_intensity)
}
}
lis_delint_nb <- list()
lis_median_nb <- list()
lis_delint_b <- list()
lis_median_b <- list()
if( logfile$parameters$ISnorm_medblank_pos == "TRUE" || logfile$parameters$ISnorm_medsam_pos == "TRUE" ){
if(logfile$parameters$ISnorm_medsam_pos == "TRUE"){
for(p in 1:length(profileList_pos[["sampleID"]])){
lis_delint_nb[[p]] <- numeric(0)
lis_median_nb[[p]] <- numeric(0)
}
names(lis_delint_nb) <- samIDs
names(lis_median_nb) <- samIDs
ids_nb <- order(profileList_pos[["index_prof"]][,"number_peaks_sample"], decreasing = TRUE)
ids_nb <- ids_nb[
( profileList_pos[["index_prof"]][ids_nb, "number_peaks_sample"] > 0 ) &
( profileList_pos[["index_prof"]][ids_nb, "number_peaks_blind"] == 0 )
]
if( length(ids_nb) ){
ids_nb <- ids_nb[unused_profile[ids_nb]]
}
if(logfile$parameters$ISnorm_usesubsam_pos == "TRUE" & (length(ids_nb) > as.numeric(logfile$parameters$ISnorm_numsam_pos))){
ids_nb <- sample(ids_nb, size = as.numeric(logfile$parameters$ISnorm_numsam_pos), replace = FALSE)
}
if( length(ids_nb) ){
for(i in ids_nb){
for_peaks <- profileList_pos[["index_prof"]][i,"start_ID"]:profileList_pos[["index_prof"]][i,"end_ID"]
median_intensity <- median( log10(profileList_pos[["peaks"]][for_peaks, "intensity"]) )
for(j in for_peaks){
at_file <- match(profileList_pos[["peaks"]][j, "sampleIDs"], samIDs)
at_int <- log10(profileList_pos[["peaks"]][j, "intensity"][[1]])
lis_delint_nb[[at_file]] <- c(lis_delint_nb[[at_file]], (at_int - median_intensity))
lis_median_nb[[at_file]] <- c(lis_median_nb[[at_file]], median_intensity)
}
}
}
}
if((logfile$parameters$ISnorm_medblank_pos == "TRUE") & (any(profileList_pos[["type"]] == "blank"))){
for(p in 1:length(profileList_pos[["sampleID"]])){
lis_delint_b[[p]] <- numeric(0)
lis_median_b[[p]] <- numeric(0)
}
names(lis_delint_b) <- samIDs
names(lis_median_b) <- samIDs
ids_b <- order(profileList_pos[["index_prof"]][, "number_peaks_blind"], decreasing =TRUE)
ids_b <- ids_b[( profileList_pos[["index_prof"]][ids_b, "number_peaks_blind"] > 0 )]
if( length(ids_b) ){
ids_b <- ids_b[unused_profile[ids_b]]
}
if(logfile$parameters$ISnorm_usesubblank_pos == "TRUE" & (length(ids_b) > as.numeric(logfile$parameters$ISnorm_numblank_pos))){
ids_b <- sample(ids_b, size = as.numeric(logfile$parameters$ISnorm_numblank_pos), replace = FALSE)
}
if( length(ids_b) ){
for(i in ids_b){
for_peaks <- profileList_pos[["index_prof"]][i,"start_ID"]:profileList_pos[["index_prof"]][i,"end_ID"]
median_intensity <- median( log10(profileList_pos[["peaks"]][for_peaks, "intensity"]) )
for(j in for_peaks){
at_file <- match(profileList_pos[["peaks"]][j, "sampleIDs"], samIDs)
at_int <- log10(profileList_pos[["peaks"]][j, "intensity"][[1]])
lis_delint_b[[at_file]] <- c(lis_delint_b[[at_file]], (at_int - median_intensity))
lis_median_b[[at_file]] <- c(lis_median_b[[at_file]], median_intensity)
}
}
}
}
}
corfac <- rep(1, length(lis_delint_IS))
use_corfac <- rep(FALSE, length(lis_delint_IS))
for(k in 1:length(lis_delint_IS)){
if( length(lis_delint_IS[[k]]) >= as.numeric(logfile$parameters$ISnorm_numbIS_pos) ){
corfac[k] <- (10 ^ median(lis_delint_IS[[k]]))
use_corfac[k] <- TRUE
}
}
sampleID <- profileList_pos[["sampleID"]];
corr_intens <- .Call("_enviMass_correct_intens",
as.numeric(corfac),
as.integer(sampleID),
as.numeric(profileList_pos[["peaks"]][,"intensity"]),
as.integer(profileList_pos[["peaks"]][,"sampleIDs"]),
PACKAGE = "enviMass"
)
profileList_pos[["peaks"]][,"intensity"] <<- corr_intens
for(k in 1:dim(profileList_pos[["index_prof"]])[1]){
profileList_pos[["index_prof"]][k,"mean_int"] <<- mean(profileList_pos[["peaks"]][(profileList_pos[["index_prof"]][k,"start_ID"]:profileList_pos[["index_prof"]][k,"end_ID"]),"intensity"])
}
int_norm_ISTD_pos <- list()
int_norm_ISTD_pos[[1]] <- lis_delint_IS
int_norm_ISTD_pos[[2]] <- lis_median_IS
int_norm_ISTD_pos[[4]] <- use_corfac
int_norm_ISTD_pos[[5]] <- lis_delint_nb
int_norm_ISTD_pos[[6]] <- lis_median_nb
int_norm_ISTD_pos[[7]] <- lis_delint_b
int_norm_ISTD_pos[[8]] <- lis_median_b
int_norm_ISTD_pos[[9]] <- atPOSIX
int_norm_ISTD_pos[[10]] <- type
int_norm_ISTD_pos[[11]] <- samIDs
names(int_norm_ISTD_pos) <- c("lis_delint_IS", " lis_median_IS", "lis_RT_IS", "use_corfac", "lis_delint_nb",
"lis_median_nb", "lis_delint_b", "lis_median_b", "atPOSIX", "sampletype", "sampleID")
output$int_norm_ISTD_pos_median <- renderPlot({
par(mar = c(.2, 4.5, .9, 8))
enviMass:::TICFG4083H(
int_norm_ISTD = int_norm_ISTD_pos,
logfile = logfile,
what = "normalization"
)
},res = 100)
output$int_norm_ISTD_pos_counts <- renderPlot({
par(mar = c(4.5, 4.5, .9, 8))
enviMass:::TICFG4083H(
int_norm_ISTD = int_norm_ISTD_pos,
logfile = logfile,
what = "counts"
)
},res = 100)
if(as.logical(logfile$parameters$test)){
for(i in 1:dim(profileList_pos[["index_prof"]])[1]){
if(
!all(profileList_pos[["peaks"]][
profileList_pos[["index_prof"]][i, "start_ID"]:profileList_pos[["index_prof"]][i, "end_ID"]
,"profileIDs"] == i)
){
stop("\n Debug do_IS_normaliz.r at _2")
}
}
}
save(int_norm_ISTD_pos, file = file.path(as.character(logfile[[1]]), "results", "int_norm_ISTD_pos"));
rm(int_norm_ISTD_pos, envir = as.environment(".GlobalEnv"))
}
save(profileList_pos, file = file.path(as.character(logfile[[1]]), "results", "profileList_pos"));
save(links_profiles_pos, file = file.path(as.character(logfile[[1]]), "results", "links_profiles_pos"));
if(any(ls()=="profileList_pos")){stop("\n illegal profileList_pos detected _1 in do_IS_normaliz.r!")}
if(any(ls()=="links_profiles_pos")){stop("\n illegal links_profiles_pos detected _1 in do_IS_normaliz.r!")}
rm(links_profiles_pos, profileList_pos, envir=as.environment(".GlobalEnv"))
}
if(
file.exists(file.path(logfile[[1]],"results","profileList_neg")) &
(as.logical(logfile$parameters$ISnorm_include_neg) | logfile$workflow[names(logfile$workflow) == "components_profiles"] == "yes")
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
use_entries_profiles <- enviMass::VOLZN1972A(links_profiles_neg)
profileList_neg[["index_prof"]][,"links"] <<- 0
with_bar<-FALSE
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
links_profiles_neg[[at_entry]] <<- enviMass::RHHPS2808T(profileList_neg[["index_prof"]][i,"number_peaks_total"][[1]])
names(links_profiles_neg)[at_entry] <<- as.character(i)
profileList_neg[["index_prof"]][i,"links"] <<- at_entry
}else{
at_entry<-profileList_neg[["index_prof"]][i,"links"]
}
for(j in (profileList_neg[["index_prof"]][i,"start_ID"]:profileList_neg[["index_prof"]][i,"end_ID"])){
if(profileList_neg[["peaks"]][j,"links"]!=0){
if(	length(links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[2]]) > 0 ){
for(k in 1:length(links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[2]])){
if( length(links_profiles_neg[[at_entry]][[2]])==0 ){
links_profiles_neg[[at_entry]][[2]] <<-
data.frame(
links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[2]][[k]][[1]],
1,
links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[2]][[k]][[2]],
stringsAsFactors = FALSE
)
names(links_profiles_neg[[at_entry]][[2]]) <<- c("Compound","Counts","max_score")
}else{
at<-which(links_profiles_neg[[at_entry]][[2]][,1] == links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[2]][[k]][[1]])
if(length(at)>0){
links_profiles_neg[[at_entry]][[2]][at,2] <<- (links_profiles_neg[[at_entry]][[2]][at,2]+1)
if(links_profiles_neg[[at_entry]][[2]][at,2] < links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[2]][[k]][[2]] ){
links_profiles_neg[[at_entry]][[2]][at,2] <<- links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[2]][[k]][[2]]
}
}else{
links_profiles_neg[[at_entry]][[2]] <<- data.frame(
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
names(links_profiles_neg[[at_entry]][[2]]) <<- c("Compound","Counts","max_score")
}
}
}
}
if(	length(links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[1]])>0 ){
for(k in 1:length(links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[1]])){
if( length(links_profiles_neg[[at_entry]][[1]]) == 0){
links_profiles_neg[[at_entry]][[1]] <<-
data.frame(
links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[1]][[k]][[1]],
1,
links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[1]][[k]][[2]],
stringsAsFactors = FALSE
)
names(links_profiles_neg[[at_entry]][[1]]) <<- c("Compound","Counts","max_score")
}else{
at <- which(links_profiles_neg[[at_entry]][[1]][,1] == links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[1]][[k]][[1]])
if(length(at)>0){
links_profiles_neg[[at_entry]][[1]][at,2] <<- (links_profiles_neg[[at_entry]][[1]][at,2]+1)
if(links_profiles_neg[[at_entry]][[1]][at,2] < links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[1]][[k]][[2]] ){
links_profiles_neg[[at_entry]][[1]][at,2] <<- links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[1]][[k]][[2]]
}
}else{
links_profiles_neg[[at_entry]][[1]] <<- data.frame(
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
names(links_profiles_neg[[at_entry]][[1]]) <<- c("Compound","Counts","max_score")
}
}
}
}
}
}
}
}
if(with_bar){close(pBar)}
}
if(!length(links_profiles_neg)){
stop("Certainly not enough internal standard (ISTD) screening matches detected to run the ISTD normalization, negative ionization - remove this workflow step?")
}
if(as.logical(logfile$parameters$ISnorm_include_neg)){
if(as.logical(logfile$parameters$screen_IS_restrict)){
min_count <- min(c(
floor(as.numeric(logfile$parameters$screen_IS_restrict_many) * as.numeric(logfile$parameters$ISnorm_percfiles_neg) / 100),
floor(length(profileList_neg[[4]]) * as.numeric(logfile$parameters$ISnorm_percfiles_neg) / 100)
))
}else{
min_count <- floor(length(profileList_neg[[4]]) * as.numeric(logfile$parameters$ISnorm_percfiles_neg) / 100);
}
lis_delint_IS <- list()
lis_median_IS <- list()
samIDs <- profileList_neg[["sampleID"]]
atPOSIX <- profileList_neg[["datetime"]]
type <- profileList_neg[["type"]]
ord <- order(atPOSIX, type, decreasing = FALSE)
samIDs <- samIDs[ord]
type <- type[ord]
atPOSIX <- atPOSIX[ord]
for(p in 1:length(profileList_neg[["sampleID"]])){
lis_delint_IS[[p]] <- numeric(0)
lis_median_IS[[p]] <- numeric(0)
}
names(lis_delint_IS) <- samIDs
names(lis_median_IS) <- samIDs
unused_profile <- rep(TRUE, dim(profileList_neg[["index_prof"]])[1])
for(i in 1:length(links_profiles_neg)){
if(!length(links_profiles_neg[[i]]$IS)) next
these <- which((links_profiles_neg[[i]]$IS$Counts >= min_count) & (links_profiles_neg[[i]]$IS$max_score >= as.numeric(logfile$parameters$ISnorm_score_neg) ) )
if(!length(these)) next
at_profile <- as.numeric(names(links_profiles_neg)[i])
unused_profile[at_profile] <- FALSE
if(profileList_neg[["index_prof"]][at_profile,"number_peaks_total"] < min_count) next
for_peaks <- profileList_neg[["index_prof"]][at_profile,"start_ID"]:profileList_neg[["index_prof"]][at_profile,"end_ID"]
median_intensity <- median( log10(profileList_neg[["peaks"]][for_peaks, "intensity"]) )
for(j in for_peaks){
at_file <- match(profileList_neg[["peaks"]][j,"sampleIDs"], samIDs)
at_int <- log10(profileList_neg[["peaks"]][j,"intensity"][[1]])
lis_delint_IS[[at_file]] <- c(lis_delint_IS[[at_file]], (at_int - median_intensity))
lis_median_IS[[at_file]] <- c(lis_median_IS[[at_file]], median_intensity)
}
}
lis_delint_nb <- list()
lis_median_nb <- list()
lis_delint_b <- list()
lis_median_b <- list()
if( logfile$parameters$ISnorm_medblank_neg == "TRUE" || logfile$parameters$ISnorm_medsam_neg == "TRUE" ){
if(logfile$parameters$ISnorm_medsam_neg == "TRUE"){
for(p in 1:length(profileList_neg[["sampleID"]])){
lis_delint_nb[[p]] <- numeric(0)
lis_median_nb[[p]] <- numeric(0)
}
names(lis_delint_nb) <- samIDs
names(lis_median_nb) <- samIDs
ids_nb <- order(profileList_neg[["index_prof"]][,"number_peaks_sample"], decreasing = TRUE)
ids_nb <- ids_nb[
( profileList_neg[["index_prof"]][ids_nb, "number_peaks_sample"] > 0 ) &
( profileList_neg[["index_prof"]][ids_nb, "number_peaks_blind"] == 0 )
]
if( length(ids_nb) ){
ids_nb <- ids_nb[unused_profile[ids_nb]]
}
if(logfile$parameters$ISnorm_usesubsam_neg == "TRUE" & (length(ids_nb) > as.numeric(logfile$parameters$ISnorm_numsam_neg))){
ids_nb <- sample(ids_nb, size = as.numeric(logfile$parameters$ISnorm_numsam_neg), replace = FALSE)
}
if( length(ids_nb) ){
for(i in ids_nb){
for_peaks <- profileList_neg[["index_prof"]][i,"start_ID"]:profileList_neg[["index_prof"]][i,"end_ID"]
median_intensity <- median( log10(profileList_neg[["peaks"]][for_peaks, "intensity"]) )
for(j in for_peaks){
at_file <- match(profileList_neg[["peaks"]][j, "sampleIDs"], samIDs)
at_int <- log10(profileList_neg[["peaks"]][j, "intensity"][[1]])
lis_delint_nb[[at_file]] <- c(lis_delint_nb[[at_file]], (at_int - median_intensity))
lis_median_nb[[at_file]] <- c(lis_median_nb[[at_file]], median_intensity)
}
}
}
}
if((logfile$parameters$ISnorm_medblank_neg == "TRUE") & (any(profileList_pos[["type"]] == "blank"))){
for(p in 1:length(profileList_neg[["sampleID"]])){
lis_delint_b[[p]] <- numeric(0)
lis_median_b[[p]] <- numeric(0)
}
names(lis_delint_b) <- samIDs
names(lis_median_b) <- samIDs
ids_b <- order(profileList_neg[["index_prof"]][, "number_peaks_blind"], decreasing =TRUE)
ids_b <- ids_b[( profileList_neg[["index_prof"]][ids_b, "number_peaks_blind"] > 0 )]
if( length(ids_b) ){
ids_b <- ids_b[unused_profile[ids_b]]
}
if(logfile$parameters$ISnorm_usesubblank_neg == "TRUE" & (length(ids_b) > as.numeric(logfile$parameters$ISnorm_numblank_neg))){
ids_b <- sample(ids_b, size = as.numeric(logfile$parameters$ISnorm_numblank_neg), replace = FALSE)
}
if( length(ids_b) ){
for(i in ids_b){
for_peaks <- profileList_neg[["index_prof"]][i,"start_ID"]:profileList_neg[["index_prof"]][i,"end_ID"]
median_intensity <- median( log10(profileList_neg[["peaks"]][for_peaks, "intensity"]) )
for(j in for_peaks){
at_file <- match(profileList_neg[["peaks"]][j, "sampleIDs"], samIDs)
at_int <- log10(profileList_neg[["peaks"]][j, "intensity"][[1]])
lis_delint_b[[at_file]] <- c(lis_delint_b[[at_file]], (at_int - median_intensity))
lis_median_b[[at_file]] <- c(lis_median_b[[at_file]], median_intensity)
}
}
}
}
}
corfac <- rep(1, length(lis_delint_IS))
use_corfac <- rep(FALSE, length(lis_delint_IS))
for(k in 1:length(lis_delint_IS)){
if( length(lis_delint_IS[[k]]) >= as.numeric(logfile$parameters$ISnorm_numbIS_neg) ){
corfac[k] <- (10 ^ median(lis_delint_IS[[k]]))
use_corfac[k] <- TRUE
}
}
sampleID <- profileList_neg[["sampleID"]];
corr_intens <- .Call("_enviMass_correct_intens",
as.numeric(corfac),
as.integer(sampleID),
as.numeric(profileList_neg[["peaks"]][,"intensity"]),
as.integer(profileList_neg[["peaks"]][,"sampleIDs"]),
PACKAGE = "enviMass"
)
profileList_neg[["peaks"]][,"intensity"] <<- corr_intens
for(k in 1:dim(profileList_neg[["index_prof"]])[1]){
profileList_neg[["index_prof"]][k,"mean_int"] <<- mean(profileList_neg[["peaks"]][(profileList_neg[["index_prof"]][k,"start_ID"]:profileList_neg[["index_prof"]][k,"end_ID"]),"intensity"])
}
int_norm_ISTD_neg <<- list()
int_norm_ISTD_neg[[1]] <<- lis_delint_IS
int_norm_ISTD_neg[[2]] <<- lis_median_IS
int_norm_ISTD_neg[[4]] <<- use_corfac
int_norm_ISTD_neg[[5]] <<- lis_delint_nb
int_norm_ISTD_neg[[6]] <<- lis_median_nb
int_norm_ISTD_neg[[7]] <<- lis_delint_b
int_norm_ISTD_neg[[8]] <<- lis_median_b
int_norm_ISTD_neg[[9]] <<- atPOSIX
int_norm_ISTD_neg[[10]] <<- type
int_norm_ISTD_neg[[11]] <<- samIDs
names(int_norm_ISTD_neg) <<- c("lis_delint_IS", " lis_median_IS", "lis_RT_IS", "use_corfac", "lis_delint_nb",
"lis_median_nb", "lis_delint_b", "lis_median_b", "atPOSIX", "sampletype", "sampleID")
output$int_norm_ISTD_neg_median <- renderPlot({
par(mar = c(.2, 4.5, .9, 8))
enviMass:::TICFG4083H(
int_norm_ISTD = int_norm_ISTD_neg,
logfile = logfile,
what = "normalization"
)
},res = 100)
output$int_norm_ISTD_neg_counts <- renderPlot({
par(mar = c(4.5, 4.5, .9, 8))
enviMass:::TICFG4083H(
int_norm_ISTD = int_norm_ISTD_neg,
logfile = logfile,
what = "counts"
)
},res = 100)
if(as.logical(logfile$parameters$test)){
for(i in 1:dim(profileList_neg[["index_prof"]])[1]){
if(
!all(profileList_neg[["peaks"]][
profileList_neg[["index_prof"]][i, "start_ID"]:profileList_neg[["index_prof"]][i, "end_ID"]
,"profileIDs"] == i)
){
stop("\n Debug do_IS_normaliz.r at _4")
}
}
}
save(int_norm_ISTD_neg, file = file.path(as.character(logfile[[1]]), "results", "int_norm_ISTD_neg"));
rm(int_norm_ISTD_neg, envir = as.environment(".GlobalEnv"))
}
save(profileList_neg,file=file.path(as.character(logfile[[1]]),"results","profileList_neg"));
save(links_profiles_neg,file=file.path(as.character(logfile[[1]]),"results","links_profiles_neg"));
if(any(ls()=="profileList_neg")){stop("\n illegal profileList_neg detected _1 in do_IS_normaliz.r!")}
if(any(ls()=="links_profiles_neg")){stop("\n illegal links_profiles_neg detected _1 in do_IS_normaliz.r!")}
rm(links_profiles_neg, profileList_neg, envir=as.environment(".GlobalEnv"))
}
