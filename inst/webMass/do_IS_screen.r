those<-list.files(file.path(logfile$project_folder,"results","screening"))
if(length(those)>0){
for(i in 1:length(those)){
if(grepl("IS",those[i])){
file.remove(file.path(logfile$project_folder,"results","screening",those[i]))
}
}
}
if(any(objects(envir=as.environment(".GlobalEnv"))=="LOD_splined")){rm(LOD_splined,envir=as.environment(".GlobalEnv"))}
if(any(objects()=="LOD_splined")){rm(LOD_splined)}
if(file.exists(file=file.path(logfile$project_folder,"results","LOD","LOD_splined"))){
load(file=file.path(logfile$project_folder,"results","LOD","LOD_splined"));
do_LOD<-TRUE
}else{
do_LOD<-FALSE
}
those_objects<-c(
"peaklist",
"profileList_pos",
"links_peaks_pos",
"pattern_pos_IS",
"patternRT_pos_IS",
"patternDelRT_pos_IS",
"pattern",
"pattern_delRT"
)
for(i in 1:length(those_objects)){
if(exists(those_objects[i],envir=as.environment(".GlobalEnv"))){rm(list=(those_objects[i]),envir=as.environment(".GlobalEnv"))}
if(exists(those_objects[i])){rm(list=(those_objects[i]))}
}
if(
file.exists(file.path(as.character(logfile[[1]]),"results","profileList_pos")) &
file.exists(file.path(logfile[[1]],"results","pattern_pos_IS"))
){
load(file=file.path(as.character(logfile[[1]]),"results","profileList_pos"),envir=as.environment(".GlobalEnv"));
load(file=file.path(as.character(logfile[[1]]),"results","profileList_pos_copy"),envir=as.environment(".GlobalEnv"));
if(dim(profileList_pos_copy[["peaks"]])[1]!=dim(profileList_pos[["peaks"]])[1]){
stop("Differing profileList dimensions_3 - report this issue immediately!")
}
load(file=file.path(as.character(logfile[[1]]),"results","links_peaks_pos"),envir=as.environment(".GlobalEnv"));
load(file=file.path(logfile[[1]],"results","pattern_pos_IS"),envir=as.environment(".GlobalEnv"));
pattern<<-pattern_pos_IS;rm(pattern_pos_IS,envir=as.environment(".GlobalEnv"));
load(file=file.path(logfile[[1]],"results","patternRT_pos_IS"),envir=as.environment(".GlobalEnv"));
pattern_RT<<-patternRT_pos_IS;rm(patternRT_pos_IS,envir=as.environment(".GlobalEnv"));
load(file=file.path(logfile[[1]],"results","patternDelRT_pos_IS"),envir=as.environment(".GlobalEnv"));
pattern_delRT<<-patternDelRT_pos_IS;rm(patternDelRT_pos_IS,envir=as.environment(".GlobalEnv"));
mztol<-as.numeric(logfile$parameters$IS_dmz)
ppm<-as.logical(as.character(logfile$parameters$IS_ppm))
cutint<-as.numeric(logfile$parameters$IS_intcut)
int_tol<-as.numeric(logfile$parameters$IS_inttol)
RT_tol_inside<-as.numeric(logfile$parameters$IS_drt2)
mztol_prof <- as.numeric(logfile$parameters$prof_dmz)
cut_score<-as.numeric(logfile$parameters$IS_w1)
measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
maxID <- max(as.numeric(measurements[,"ID"]))
if(logfile$parameters$screen_IS_restrict == "TRUE"){
measurements <- measurements[measurements[,"Mode"]=="positive",,drop=FALSE]
measurements <- measurements[(measurements[,"Type"]=="sample" | measurements[,"Type"]=="blank" | measurements[,"Type"]=="spiked" ),,drop=FALSE]
starttime <- as.difftime(measurements[,"Time"]);
startdate <- as.Date(measurements[,"Date"], tz="GMT");
numstart <- (as.numeric(startdate)+as.numeric(starttime/(24*60*60)))
if( length(numstart) > as.numeric(logfile$parameters$screen_IS_restrict_many) ){
retain_sample <- rep(FALSE, maxID)
retain_sample[
as.numeric(measurements[
(order(numstart,decreasing=TRUE)[1:as.numeric(logfile$parameters$screen_IS_restrict_many)])
,"ID"])
] <- TRUE
}else{
retain_sample <- rep(TRUE, maxID)
}
}else{
retain_sample <- rep(TRUE, maxID)
}
rm(measurements)
peaks<-profileList_pos[["index_prof"]];
peaklist<-peaks[,c("mean_mz","mean_int","mean_RT")];
count_nonmax<-0
for(i in 1:length(pattern)){
count_nonmax<-(count_nonmax+
length(pattern[[i]][,1])
)
}
centro_mass<-rep(0,count_nonmax)
centro_ID<-rep(0,count_nonmax)
centro_maxpeak<-rep(FALSE,count_nonmax)
centro_number<-rep(0,count_nonmax)
centro_RT<-rep(0,count_nonmax)
centro_dRT<-rep(0,count_nonmax)
at_ID<-1
screen_list<-as.list(rep("FALSE",length(pattern)))
for(i in 1:length(pattern)){
n<-length(pattern[[i]][,1])
centro_mass[at_ID:(at_ID+n-1)]<-pattern[[i]][,1]
centro_ID[at_ID:(at_ID+n-1)]<-i
centro_maxpeak[at_ID:(at_ID+n-1)]<-(pattern[[i]][,2]==max(pattern[[i]][,2]))
centro_number[at_ID:(at_ID+n-1)]<-(1:n)
centro_RT[at_ID:(at_ID+n-1)]<-pattern_RT[i]
centro_dRT[at_ID:(at_ID+n-1)]<-pattern_delRT[i]
screen_list[[i]]<-as.list(rep("FALSE",n))
at_ID<-(at_ID+n)
}
getit <- BYRIB9287R(
peaklist,
centro_mass,
dmz = (mztol + 2 * mztol_prof),
ppm = ppm,
RT = centro_RT,
dRT = (centro_dRT + as.numeric(logfile$parameters$prof_drt))
)
if( as.character(logfile$parameters$screen_IS_maxonly)=="TRUE" ){
getit[!centro_maxpeak]<-"FALSE"
}
for(i in 1:length(getit)){
screen_list[[centro_ID[i]]][[centro_number[i]]]<-getit[i]
}
IS_pos_screen_listed<-list()
set_ID<-seq(1:length(profileList_pos[[4]]))
for(i in 1:length(screen_list)){
if(any(is.na(screen_list[[i]]==FALSE))){
IS_pos_screen_listed[[i]]<-list()
for(j in 1:length(screen_list[[i]])){
if(screen_list[[i]][[j]]!="FALSE"){
profs<-as.numeric(strsplit(screen_list[[i]][[j]]," / ")[[1]])
for(k in 1:length(profs)){
if(profileList_pos[[7]][profs[k],4] != profs[k]){cat("\n debug me: profile ID mismatch");stop();} # a check
for(m in profileList_pos[[7]][profs[k],1]:profileList_pos[[7]][profs[k],2]){
if(retain_sample[profileList_pos[[2]][m,"sampleIDs"]] == FALSE) next
delmass <- abs(profileList_pos[[2]][m,1] - pattern[[i]][j,1])
if(!ppm){
if(delmass > (mztol / 1000)) next
}else{
if((delmass * 1E6 / pattern[[i]][j,1]) > mztol) next
}
delRT <- abs(profileList_pos[[2]][m,3] - pattern_RT[i])
if(delRT > pattern_delRT[i]) next
at_ID<-set_ID[profileList_pos[[4]]==as.character(profileList_pos[[2]][m,6])]
if(length(IS_pos_screen_listed[[i]])<at_ID){
IS_pos_screen_listed[[i]][[at_ID]]<-matrix(ncol=2,nrow=0)
}else{
if(length(IS_pos_screen_listed[[i]][[at_ID]])==0){
IS_pos_screen_listed[[i]][[at_ID]]<-matrix(ncol=2,nrow=0)
}
}
IS_pos_screen_listed[[i]][[at_ID]]<-rbind(
IS_pos_screen_listed[[i]][[at_ID]],c(j,m)
)
colnames(IS_pos_screen_listed[[i]][[at_ID]])<-c(as.character(profileList_pos[[4]][at_ID]),"")
}
}
}
}
}else{
IS_pos_screen_listed[[i]]<-numeric(0)
}
}
names(IS_pos_screen_listed)<-names(pattern)
if( logfile$parameters$screen_IS_cutit=="TRUE" ){
use_score_cut<-TRUE;
score_cut<-cut_score
}else{
use_score_cut<-FALSE;
score_cut<-0
}
many<-0
many_unamb<-0
res_IS_pos_screen<-list()
if(length(IS_pos_screen_listed)>0){
for(i in 1:length(IS_pos_screen_listed)){
if(length(IS_pos_screen_listed[[i]])>0){
res_IS_pos_screen[[i]]<-list()
for(m in 1:length(IS_pos_screen_listed[[i]])){
at_ID<-set_ID[profileList_pos[[4]]==colnames(IS_pos_screen_listed[[i]][[m]])[1]]
if(length(IS_pos_screen_listed[[i]][[m]])>0){
if(do_LOD){
with_model<-which(names(LOD_splined)==paste("LOD_",colnames(IS_pos_screen_listed[[i]][[m]])[1],sep=""))
if(length(with_model)>0){
at_RT <- profileList_pos[["peaks"]][IS_pos_screen_listed[[i]][[m]][1,2], 3]
use_cutint<-10^(predict(LOD_splined[[with_model]], at_RT)$y)
}else{
cat("\n Missing LOD model; using default intensity threshold. Debug?")
use_cutint<-cutint;
}
}else{
use_cutint<-cutint
}
combination_matches <- FSFDE6485E(
cent_peak_mat = IS_pos_screen_listed[[i]][[m]],
pattern_compound = pattern[[i]],
peaks = profileList_pos[["peaks"]],
LOD = use_cutint,
RT_tol_inside = RT_tol_inside,
int_tol = int_tol,
use_score_cut = use_score_cut,
score_cut = score_cut,
plot_it = FALSE,
verbose = FALSE,
RT_seperate = TRUE
)
for(k in 1:length(combination_matches)){
combination_matches[[k]][[10]]<-colnames(IS_pos_screen_listed[[i]][[m]])[1]
names(combination_matches[[k]])[10]<-"file_ID"
}
res_IS_pos_screen[[i]][[at_ID]] <- combination_matches
names(res_IS_pos_screen[[i]])[[at_ID]] <- combination_matches[[k]][[10]]
if(length(combination_matches) > 1){many_unamb<-(many_unamb+1)}
many<-(many+1)
}
}
}else{
res_IS_pos_screen[[i]]<-numeric(0)
}
}
names(res_IS_pos_screen)<-names(IS_pos_screen_listed)
}
save(res_IS_pos_screen,file=file.path(logfile$project_folder,"results","screening","res_IS_pos_screen"))
if( length(IS_pos_screen_listed)>0 ){
measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
measurements<-measurements[measurements[,"include"]=="TRUE",]
if(logfile$parameters$prof_select=="TRUE"){
measurements<-measurements[measurements[,names(measurements)=="profiled"]=="TRUE",]
}
intstand<-read.table(file=file.path(logfile[[1]],"dataframes","IS.txt"),header=TRUE,sep="\t",colClasses = "character");
results_screen_IS_pos<-enviMass::QGFEB5722K(
screened_listed=res_IS_pos_screen,
pattern=pattern,
at_RT=pattern_RT,
profileList=profileList_pos,
measurements_table=measurements,
compound_table=intstand,
cut_score=cut_score
)
results_screen_IS_pos[[1]]<-results_screen_IS_pos[[1]][,-c(11,12)]
results_screen_IS_pos[[2]]<-results_screen_IS_pos[[2]][,-c(8,9)]
use_entries<-VOLZN1972A(links_peaks_pos)
for(i in 1:length(res_IS_pos_screen)){
if(length(res_IS_pos_screen[[i]])>0){
for(m in 1:length(res_IS_pos_screen[[i]])){
if(length(res_IS_pos_screen[[i]][[m]])>0){
for(k in 1:length(res_IS_pos_screen[[i]][[m]])){
local_score<-0
if(!is.na(res_IS_pos_screen[[i]][[m]][[k]]$score_1)){
local_score<-(local_score+res_IS_pos_screen[[i]][[m]][[k]]$score_1)
}
if( (local_score>=1) || (is.na(res_IS_pos_screen[[i]][[m]][[k]]$score_1)) ){
if(!is.na(res_IS_pos_screen[[i]][[m]][[k]]$score_2)){
local_score<-(local_score+res_IS_pos_screen[[i]][[m]][[k]]$score_2)
}
}
if(local_score>=cut_score){
for(a in 1:length(res_IS_pos_screen[[i]][[m]][[k]]$Peaks[,2]) ){
if(profileList_pos[[2]][res_IS_pos_screen[[i]][[m]][[k]]$Peaks[a,2],5]==0){
if(length(use_entries)>0){
at_entry<-use_entries[1]
use_entries<-use_entries[-1]
}else{
at_entry<-(length(links_peaks_pos)+1)
}
links_peaks_pos[[at_entry]]<-list()
links_peaks_pos[[at_entry]][[1]]<-list()
links_peaks_pos[[at_entry]][[2]]<-list()
links_peaks_pos[[at_entry]][[3]]<-list()
links_peaks_pos[[at_entry]][[4]]<-list()
links_peaks_pos[[at_entry]][[5]]<-list()
links_peaks_pos[[at_entry]][[6]]<-list()
profileList_pos[[2]][res_IS_pos_screen[[i]][[m]][[k]]$Peaks[a,2],5]<<-at_entry
profileList_pos_copy[[2]][res_IS_pos_screen[[i]][[m]][[k]]$Peaks[a,2],5]<<-at_entry
links_peaks_pos[[at_entry]][[2]][[1]] <- list()
links_peaks_pos[[at_entry]][[2]][[1]][[1]] <- names(pattern)[i]
links_peaks_pos[[at_entry]][[2]][[1]][[2]] <- local_score
}else{
at_entry<-profileList_pos[[2]][res_IS_pos_screen[[i]][[m]][[k]]$Peaks[a,2],5]
at_list<-(length(links_peaks_pos[[at_entry]][[2]])+1)
links_peaks_pos[[at_entry]][[2]][[at_list]] <- list()
links_peaks_pos[[at_entry]][[2]][[at_list]][[1]] <- names(pattern)[i]
links_peaks_pos[[at_entry]][[2]][[at_list]][[2]] <- local_score
}
}
}
}
}
}
}
}
save(profileList_pos,file=file.path(as.character(logfile[[1]]),"results","profileList_pos"),compress=FALSE)
save(profileList_pos_copy,file=file.path(as.character(logfile[[1]]),"results","profileList_pos_copy"),compress=FALSE)
save(links_peaks_pos,file=file.path(as.character(logfile[[1]]),"results","links_peaks_pos"))
save(results_screen_IS_pos,file=file.path(logfile$project_folder,"results","screening","results_screen_IS_pos"))
rm(measurements,intstand,results_screen_IS_pos);
}
rm(getit,IS_pos_screen_listed,res_IS_pos_screen)
rm(pattern,pattern_RT,pattern_delRT,profileList_pos,profileList_pos_copy,envir=as.environment(".GlobalEnv"))
}
those_objects<-c(
"peaklist",
"profileList_neg",
"links_peaks_neg",
"pattern_neg_IS",
"patternRT_neg_IS",
"patternDelRT_neg_IS",
"pattern",
"pattern_delRT"
)
for(i in 1:length(those_objects)){
if(exists(those_objects[i],envir=as.environment(".GlobalEnv"))){rm(list=(those_objects[i]),envir=as.environment(".GlobalEnv"))}
if(exists(those_objects[i])){rm(list=(those_objects[i]))}
}
if(
file.exists(file.path(as.character(logfile[[1]]),"results","profileList_neg")) &
file.exists(file.path(logfile[[1]],"results","pattern_neg_IS"))
){
load(file=file.path(as.character(logfile[[1]]),"results","profileList_neg"),envir=as.environment(".GlobalEnv"));
load(file=file.path(as.character(logfile[[1]]),"results","profileList_neg_copy"),envir=as.environment(".GlobalEnv"));
if(dim(profileList_neg_copy[["peaks"]])[1]!=dim(profileList_neg[["peaks"]])[1]){
stop("Differing profileList dimensions_4 - report this issue immediately!")
}
load(file=file.path(as.character(logfile[[1]]),"results","links_peaks_neg"),envir=as.environment(".GlobalEnv"));
load(file=file.path(logfile[[1]],"results","pattern_neg_IS"),envir=as.environment(".GlobalEnv"));
pattern<<-pattern_neg_IS;rm(pattern_neg_IS,envir=as.environment(".GlobalEnv"));
load(file=file.path(logfile[[1]],"results","patternRT_neg_IS"),envir=as.environment(".GlobalEnv"));
pattern_RT<<-patternRT_neg_IS;rm(patternRT_neg_IS,envir=as.environment(".GlobalEnv"));
load(file=file.path(logfile[[1]],"results","patternDelRT_neg_IS"),envir=as.environment(".GlobalEnv"));
pattern_delRT<<-patternDelRT_neg_IS;rm(patternDelRT_neg_IS,envir=as.environment(".GlobalEnv"));
if(FALSE){
use_that<-names(pattern)=="30n_M-H_none_none_none"
pattern_RT<<-pattern_RT[use_that]
pattern_delRT<<-pattern_delRT[use_that]
pattern<<-pattern[use_that]
}
mztol<-as.numeric(logfile$parameters$IS_dmz)
ppm<-as.logical(as.character(logfile$parameters$IS_ppm))
cutint<-as.numeric(logfile$parameters$IS_intcut)
int_tol<-as.numeric(logfile$parameters$IS_inttol)
RT_tol_inside<-as.numeric(logfile$parameters$IS_drt2)
mztol_prof <- as.numeric(logfile$parameters$prof_dmz)
cut_score<-as.numeric(logfile$parameters$IS_w1)
measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
maxID <- max(as.numeric(measurements[,"ID"]))
if(logfile$parameters$screen_IS_restrict == "TRUE"){
measurements <- measurements[measurements[,"Mode"]=="negative",,drop=FALSE]
measurements <- measurements[(measurements[,"Type"]=="sample" | measurements[,"Type"]=="blank" | measurements[,"Type"]=="spiked" ),,drop=FALSE]
starttime <- as.difftime(measurements[,"Time"]);
startdate <- as.Date(measurements[,"Date"], tz="GMT");
numstart <- (as.numeric(startdate)+as.numeric(starttime/(24*60*60)))
if( length(numstart) > as.numeric(logfile$parameters$screen_IS_restrict_many) ){
retain_sample <- rep(FALSE, maxID)
retain_sample[
as.numeric(measurements[
(order(numstart,decreasing=TRUE)[1:as.numeric(logfile$parameters$screen_IS_restrict_many)])
,"ID"])
] <- TRUE
}else{
retain_sample <- rep(TRUE, maxID)
}
}else{
retain_sample <- rep(TRUE, maxID)
}
rm(measurements)
peaks<-profileList_neg[["index_prof"]];
peaklist<-peaks[,c("mean_mz","mean_int","mean_RT")];
count_nonmax<-0
for(i in 1:length(pattern)){
count_nonmax<-(count_nonmax+
length(pattern[[i]][,1])
)
}
centro_mass<-rep(0,count_nonmax)
centro_ID<-rep(0,count_nonmax)
centro_maxpeak<-rep(FALSE,count_nonmax)
centro_number<-rep(0,count_nonmax)
centro_RT<-rep(0,count_nonmax)
centro_dRT<-rep(0,count_nonmax)
at_ID<-1
screen_list<-as.list(rep("FALSE",length(pattern)))
for(i in 1:length(pattern)){
n<-length(pattern[[i]][,1])
centro_mass[at_ID:(at_ID+n-1)]<-pattern[[i]][,1]
centro_ID[at_ID:(at_ID+n-1)]<-i
centro_maxpeak[at_ID:(at_ID+n-1)]<-(pattern[[i]][,2]==max(pattern[[i]][,2]))
centro_number[at_ID:(at_ID+n-1)]<-(1:n)
centro_RT[at_ID:(at_ID+n-1)]<-pattern_RT[i]
centro_dRT[at_ID:(at_ID+n-1)]<-pattern_delRT[i]
screen_list[[i]]<-as.list(rep("FALSE",n))
at_ID<-(at_ID+n)
}
getit <- BYRIB9287R(
peaklist,
centro_mass,
dmz = (mztol + 2 * mztol_prof),
ppm = ppm,
RT = centro_RT,
dRT = (centro_dRT + as.numeric(logfile$parameters$prof_drt))
)
if( as.character(logfile$parameters$screen_IS_maxonly)=="TRUE" ){
getit[!centro_maxpeak]<-"FALSE"
}
for(i in 1:length(getit)){
screen_list[[centro_ID[i]]][[centro_number[i]]]<-getit[i]
}
IS_neg_screen_listed<-list()
set_ID<-seq(1:length(profileList_neg[[4]]))
for(i in 1:length(screen_list)){
if(any(is.na(screen_list[[i]]==FALSE))){
IS_neg_screen_listed[[i]]<-list()
for(j in 1:length(screen_list[[i]])){
if(screen_list[[i]][[j]]!="FALSE"){
profs<-as.numeric(strsplit(screen_list[[i]][[j]]," / ")[[1]])
for(k in 1:length(profs)){
if(profileList_neg[[7]][profs[k],4]!=profs[k]){cat("\n debug me: profile ID mismatch");stop();} # a check
for(m in profileList_neg[[7]][profs[k],1]:profileList_neg[[7]][profs[k],2]){
if(retain_sample[profileList_neg[[2]][m,"sampleIDs"]]==FALSE){next}
delmass<-abs(profileList_neg[[2]][m,1]-pattern[[i]][j,1])
if(!ppm){
if(delmass>(mztol/1000)){next}
}else{
if((delmass*1E6/pattern[[i]][j,1])>mztol){next}
}
delRT <- abs(profileList_neg[[2]][m,3] - pattern_RT[i])
if(delRT > pattern_delRT[i]) next
at_ID<-set_ID[profileList_neg[[4]]==as.character(profileList_neg[[2]][m,6])]
if(length(IS_neg_screen_listed[[i]])<at_ID){
IS_neg_screen_listed[[i]][[at_ID]]<-matrix(ncol=2,nrow=0)
}else{
if(length(IS_neg_screen_listed[[i]][[at_ID]])==0){
IS_neg_screen_listed[[i]][[at_ID]]<-matrix(ncol=2,nrow=0)
}
}
IS_neg_screen_listed[[i]][[at_ID]]<-rbind(
IS_neg_screen_listed[[i]][[at_ID]],c(j,m)
)
colnames(IS_neg_screen_listed[[i]][[at_ID]])<-c(as.character(profileList_neg[[4]][at_ID]),"")
}
}
}
}
}else{
IS_neg_screen_listed[[i]]<-numeric(0)
}
}
names(IS_neg_screen_listed)<-names(pattern)
if( logfile$parameters$screen_IS_cutit=="TRUE" ){
use_score_cut<-TRUE;
score_cut<-cut_score
}else{
use_score_cut<-FALSE;
score_cut<-0
}
many<-0
many_unamb<-0
res_IS_neg_screen<-list()
if(length(IS_neg_screen_listed)>0){
for(i in 1:length(IS_neg_screen_listed)){
if(length(IS_neg_screen_listed[[i]])>0){
res_IS_neg_screen[[i]]<-list()
for(m in 1:length(IS_neg_screen_listed[[i]])){
at_ID<-set_ID[profileList_neg[[4]]==colnames(IS_neg_screen_listed[[i]][[m]])[1]]
if(length(IS_neg_screen_listed[[i]][[m]])>0){
if(do_LOD){
with_model<-which(names(LOD_splined)==paste("LOD_",colnames(IS_neg_screen_listed[[i]][[m]])[1],sep=""))
if(length(with_model)>0){
at_RT <- profileList_neg[["peaks"]][IS_neg_screen_listed[[i]][[m]][1,2], 3]
use_cutint<-10^(predict(LOD_splined[[with_model]], at_RT)$y)
}else{
cat("\n Missing LOD model; using default intensity threshold. Debug?")
use_cutint<-cutint;
}
}else{
use_cutint<-cutint
}
combination_matches <- FSFDE6485E(
cent_peak_mat = IS_neg_screen_listed[[i]][[m]],
pattern_compound = pattern[[i]],
peaks = profileList_neg[["peaks"]],
LOD = use_cutint,
RT_tol_inside = RT_tol_inside,
int_tol = int_tol,
use_score_cut = use_score_cut,
score_cut = score_cut,
plot_it = FALSE,
verbose = FALSE,
RT_seperate = TRUE
)
for(k in 1:length(combination_matches)){
combination_matches[[k]][[10]]<-colnames(IS_neg_screen_listed[[i]][[m]])[1]
names(combination_matches[[k]])[10]<-"file_ID"
}
res_IS_neg_screen[[i]][[at_ID]]<-combination_matches
names(res_IS_neg_screen[[i]])[[at_ID]]<-combination_matches[[k]][[10]]
if(length(combination_matches)>1){many_unamb<-(many_unamb+1)}
many<-(many+1)
}
}
}else{
res_IS_neg_screen[[i]]<-numeric(0)
}
}
names(res_IS_neg_screen)<-names(IS_neg_screen_listed)
}
save(res_IS_neg_screen,file=file.path(logfile$project_folder,"results","screening","res_IS_neg_screen"))
if( length(IS_neg_screen_listed)>0 ){
measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
measurements<-measurements[measurements[,"include"]=="TRUE",]
if(logfile$parameters$prof_select=="TRUE"){
measurements<-measurements[measurements[,names(measurements)=="profiled"]=="TRUE",]
}
intstand<-read.table(file=file.path(logfile[[1]],"dataframes","IS.txt"),header=TRUE,sep="\t",colClasses = "character");
results_screen_IS_neg<-enviMass::QGFEB5722K(
screened_listed=res_IS_neg_screen,
pattern=pattern,
at_RT=pattern_RT,
profileList=profileList_neg,
measurements_table=measurements,
compound_table=intstand,
cut_score=cut_score
)
results_screen_IS_neg[[1]]<-results_screen_IS_neg[[1]][,-c(11,12)]
results_screen_IS_neg[[2]]<-results_screen_IS_neg[[2]][,-c(8,9)]
use_entries<-VOLZN1972A(links_peaks_neg)
for(i in 1:length(res_IS_neg_screen)){
if(length(res_IS_neg_screen[[i]])>0){
for(m in 1:length(res_IS_neg_screen[[i]])){
if(length(res_IS_neg_screen[[i]][[m]])>0){
for(k in 1:length(res_IS_neg_screen[[i]][[m]])){
local_score<-0
if(!is.na(res_IS_neg_screen[[i]][[m]][[k]]$score_1)){
local_score<-(local_score+res_IS_neg_screen[[i]][[m]][[k]]$score_1)
}
if( (local_score>=1) || (is.na(res_IS_neg_screen[[i]][[m]][[k]]$score_1)) ){
if(!is.na(res_IS_neg_screen[[i]][[m]][[k]]$score_2)){
local_score<-(local_score+res_IS_neg_screen[[i]][[m]][[k]]$score_2)
}
}
if(local_score>=cut_score){
for(a in 1:length(res_IS_neg_screen[[i]][[m]][[k]]$Peaks[,2]) ){
if(profileList_neg[[2]][res_IS_neg_screen[[i]][[m]][[k]]$Peaks[a,2],5]==0){
if(length(use_entries)>0){
at_entry<-use_entries[1]
use_entries<-use_entries[-1]
}else{
at_entry<-(length(links_peaks_neg)+1)
}
links_peaks_neg[[at_entry]]<-list()
links_peaks_neg[[at_entry]][[1]]<-list()
links_peaks_neg[[at_entry]][[2]]<-list()
links_peaks_neg[[at_entry]][[3]]<-list()
links_peaks_neg[[at_entry]][[4]]<-list()
links_peaks_neg[[at_entry]][[5]]<-list()
links_peaks_neg[[at_entry]][[6]]<-list()
profileList_neg[[2]][res_IS_neg_screen[[i]][[m]][[k]]$Peaks[a,2],5]<<-at_entry
profileList_neg_copy[[2]][res_IS_neg_screen[[i]][[m]][[k]]$Peaks[a,2],5]<<-at_entry
links_peaks_neg[[at_entry]][[2]][[1]] <- list()
links_peaks_neg[[at_entry]][[2]][[1]][[1]] <- names(pattern)[i]
links_peaks_neg[[at_entry]][[2]][[1]][[2]] <- local_score
}else{
at_entry<-profileList_neg[[2]][res_IS_neg_screen[[i]][[m]][[k]]$Peaks[a,2],5]
at_list<-(length(links_peaks_neg[[at_entry]][[2]])+1)
links_peaks_neg[[at_entry]][[2]][[at_list]] <- list()
links_peaks_neg[[at_entry]][[2]][[at_list]][[1]] <- names(pattern)[i]
links_peaks_neg[[at_entry]][[2]][[at_list]][[2]] <- local_score
}
}
}
}
}
}
}
}
save(profileList_neg,file=file.path(as.character(logfile[[1]]),"results","profileList_neg"),compress=FALSE)
save(profileList_neg_copy,file=file.path(as.character(logfile[[1]]),"results","profileList_neg_copy"),compress=FALSE)
save(links_peaks_neg,file=file.path(as.character(logfile[[1]]),"results","links_peaks_neg"))
save(results_screen_IS_neg,file=file.path(logfile$project_folder,"results","screening","results_screen_IS_neg"))
rm(measurements,intstand,results_screen_IS_neg);
}
rm(getit,IS_neg_screen_listed,res_IS_neg_screen)
rm(pattern,pattern_RT,pattern_delRT,profileList_neg,profileList_neg_copy,envir=as.environment(".GlobalEnv"))
}
