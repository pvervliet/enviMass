QTGGY3365L <- function(
x,
logfile,
measurements,
...
){
for_file <- x
for_mode <- measurements[measurements$ID == for_file, "Mode"]
if(
do_isot & file.exists(file.path(logfile[[1]], "results", "componentization", "isotopologues", paste("full", for_file, sep = "_")))
){
load(file.path(logfile[[1]],"results","componentization","isotopologues", paste("full", for_file, sep = "_")))
do_pattern <- TRUE
}else{
pattern <- FALSE
do_pattern <- FALSE
}
if(
do_addu & file.exists(file.path(logfile[[1]], "results", "componentization","adducts",paste("full", for_file,sep = "_")))
){
load(file.path(logfile[[1]], "results", "componentization", "adducts", paste("full", for_file, sep = "_")))
do_adduct <- TRUE
}else{
adduct<-FALSE
do_adduct <- FALSE
}
if(
do_homol & file.exists(file.path(logfile[[1]], "results", "componentization", "homologues", paste("full", for_file, sep = "_")))
){
load(file.path(logfile[[1]],"results", "componentization", "homologues", paste("full", for_file, sep = "_")))
}else{
homol <- FALSE
}
if(!do_adduct & !do_pattern){
return(paste0("No adduct, isotopologue or HS results available for file with ID ", for_file))
}
component <- enviMass::FENUW0294V(
pattern,
adduct,
homol,
rules = c(FALSE, FALSE, FALSE),
dont = FALSE
)
component[[1]][,17] <- as.character(component[[1]][,17])
component[[1]][,18] <- as.character(component[[1]][,18])
named <- colnames(component[[1]])
component[[1]] <- cbind(
component[[1]],
rep("-",dim(component[[1]])[1]),
rep("-",dim(component[[1]])[1]),
rep(0,dim(component[[1]])[1]),
rep(0,dim(component[[1]])[1]),
rep(0,dim(component[[1]])[1]),
rep(0,dim(component[[1]])[1]),
rep(0,dim(component[[1]])[1]),
rep(0,dim(component[[1]])[1]),
rep(0,dim(component[[1]])[1])
)
component[[1]][,19] <- as.character(component[[1]][,19])
component[[1]][,20] <- as.character(component[[1]][,20])
colnames(component[[1]]) <- c(named,
"Target peaks","ISTD peaks","Total peak number","Blind peak number",
"Monois. peak ID |","Monois. m/z |","Monois. RT |","Monois. int. |","Monois. sample/blind int. ratio"
)
if( logfile$workflow[names(logfile$workflow)=="target_screen"] == "yes" ){
targets<-read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character")
if(for_mode == "positive"){
if(exists("profileList_pos_copy")){rm("profileList_pos_copy")}
if(exists("links_peaks_pos")){rm("links_peaks_pos")}
load(file = file.path(as.character(logfile[[1]]),"results","profileList_pos_copy"), envir = environment())
load(file = file.path(as.character(logfile[[1]]),"results","links_peaks_pos"), envir = environment())
for_peaks <- which(
profileList_pos_copy[["peaks"]][,"sampleIDs"] == for_file
)
peaks <- profileList_pos_copy[["peaks"]][for_peaks,]
rm(for_peaks)
if(dim(peaks)[1] > 0){
peaks<-peaks[order(peaks[,"peakIDs"],decreasing=FALSE),]
for(i in 1:dim(component[["Components"]])[1]){
collect_comp<-c()
collect_peak<-c()
if(component[["Components"]][i,"ID pattern peaks |"] != "-"){
those <- as.numeric(strsplit(component[["Components"]][i, "ID pattern peaks |"], ",")[[1]])
these <- as.numeric(those)
these <- match(these, peaks[,"peakIDs"])
these <- these[!is.na(these)]
if(length(these)){
for(n in 1:length(these) ){
if(peaks[these[n], "links"] != 0){
got <- unlist(links_peaks_pos[[peaks[these[n], "links"]]]$target)
if(!length(got)) next
got <- got[seq(1, length(got), 2)]
collect_comp <- c(collect_comp, got)
collect_peak <- c(collect_peak, rep(these[n], length(got)))
}
}
}
}
if(component[["Components"]][i,"ID adduct peaks |"] != "-"){
those <- as.numeric(strsplit(component[["Components"]][i,"ID adduct peaks |"],",")[[1]])
these <- as.numeric(those)
these <- match(these, peaks[,"peakIDs"])
these <- these[!is.na(these)]
if(length(these)){
for(n in 1:length(these) ){
if(peaks[these[n], "links"] != 0){
got <- unlist(links_peaks_pos[[peaks[these[n], "links"]]]$target)
if(!length(got)) next
got <- got[seq(1, length(got), 2)]
collect_comp <- c(collect_comp, got)
collect_peak <- c(collect_peak, rep(these[n], length(got)))
}
}
}
}
if(component[["Components"]][i,"ID interfering peaks |"] != "-"){
those <- as.numeric(strsplit(component[["Components"]][i, "ID interfering peaks |"], ",")[[1]])
these <- as.numeric(those)
these <- match(these,peaks[,"peakIDs"])
these <- these[!is.na(these)]
if(length(these)){
for(n in 1:length(these) ){
if(peaks[these[n],"links"] != 0){
got <- unlist(links_peaks_pos[[peaks[these[n], "links"]]]$target)
if(!length(got)) next
got <- got[seq(1, length(got), 2)]
collect_comp <- c(collect_comp, got)
collect_peak <- c(collect_peak, rep(these[n], length(got)))
}
}
}
}
if(length(collect_comp)>0){
got_comp<-unique(collect_comp)
collect_all<-""
for(n in 1:length(got_comp)){
tar <- strsplit(got_comp[n], "_")[[1]][1:2]
tar_name <- targets[targets[,"ID"] == tar[1], "Name"]
tar <- paste(tar_name, tar[2], sep = "_")
collect_all<-paste0(c(collect_all,
paste(tar,
paste0(unique(collect_peak[collect_comp==got_comp[n]]),collapse=","),
sep=": "
)),
collapse=" / "
)
}
if(nchar(collect_all)>3){
collect_all<-substr(collect_all,4,nchar(collect_all))
}
component[[1]][i,"Target peaks"]<-collect_all
}
}
rm(collect_comp, collect_peak)
}
rm(profileList_pos_copy, links_peaks_pos)
}
if(for_mode=="negative"){
if(exists("profileList_neg_copy",envir=as.environment(".GlobalEnv"))){rm("profileList_neg_copy",envir=as.environment(".GlobalEnv"))}
if(exists("profileList_neg_copy")){rm("profileList_neg_copy")}
if(exists("links_peaks_neg",envir=as.environment(".GlobalEnv"))){rm("links_peaks_neg",envir=as.environment(".GlobalEnv"))}
if(exists("links_peaks_neg")){rm("links_peaks_neg")}
load(file=file.path(as.character(logfile[[1]]),"results","profileList_neg_copy"), envir = environment())
load(file=file.path(as.character(logfile[[1]]),"results","links_peaks_neg"), envir = environment())
for_peaks <- which(
profileList_neg_copy[["peaks"]][,"sampleIDs"] == for_file
)
peaks <- profileList_neg_copy[["peaks"]][for_peaks,]
rm(for_peaks)
if(dim(peaks)[1] > 0){
peaks <- peaks[order(peaks[,"peakIDs"],decreasing=FALSE),]
for(i in 1:dim(component[[1]])[1]){
collect_comp <- c()
collect_peak <- c()
if(component[[1]][i,3] != "-"){
those <- as.numeric(strsplit(component[[1]][i,3],",")[[1]])
these <- as.numeric(those)
these <- match(these, peaks[,"peakIDs"])
these <- these[!is.na(these)]
if(length(these)){
for(n in 1:length(these) ){
if(peaks[these[n],"links"] != 0){
got <- unlist(links_peaks_neg[[peaks[these[n],"links"]]]$target)
if(!length(got)) next
got <- got[seq(1,length(got),2)]
collect_comp <- c(collect_comp, got)
collect_peak <- c(collect_peak, rep(these[n], length(got)))
}
}
}
}
if(component[[1]][i,5] != "-"){
those<-as.numeric(strsplit(component[[1]][i,5],",")[[1]])
these<-as.numeric(those)
these<-match(these,peaks[,"peakIDs"])
these <- these[!is.na(these)]
if(length(these)){
for(n in 1:length(these) ){
if(peaks[these[n],"links"]!=0){
got<-unlist(links_peaks_neg[[peaks[these[n],"links"]]]$target)
if(!length(got)) next
got <- got[seq(1,length(got),2)]
collect_comp<-c(collect_comp,got)
collect_peak<-c(collect_peak,rep(these[n],length(got)))
}
}
}
}
if(component[[1]][i,7]!="-"){
those<-as.numeric(strsplit(component[[1]][i,7],",")[[1]])
these<-as.numeric(those)
these<-match(these,peaks[,"peakIDs"])
these <- these[!is.na(these)]
if(length(these)){
for(n in 1:length(these) ){
if(peaks[these[n],"links"]!=0){
got<-unlist(links_peaks_neg[[peaks[these[n],"links"]]]$target)
if(!length(got)) next
got <- got[seq(1,length(got),2)]
collect_comp<-c(collect_comp,got)
collect_peak<-c(collect_peak,rep(these[n],length(got)))
}
}
}
}
if(length(collect_comp)>0){
got_comp<-unique(collect_comp)
collect_all<-""
for(n in 1:length(got_comp)){
tar <- strsplit(got_comp[n], "_")[[1]][1:2]
tar_name <- targets[targets[,"ID"] == tar[1], "Name"]
tar <- paste(tar_name, tar[2], sep = "_")
collect_all<-paste0(c(collect_all,
paste(tar,
paste0(unique(collect_peak[collect_comp==got_comp[n]]),collapse=","),
sep=": "
)),
collapse=" / "
)
}
if(nchar(collect_all)>3){
collect_all<-substr(collect_all,4,nchar(collect_all))
}
component[[1]][i,"Target peaks"]<-collect_all
}
}
rm(collect_comp, collect_peak)
}
rm(profileList_neg_copy, links_peaks_neg)
}
}
if(logfile$workflow[names(logfile$workflow)=="IS_screen"]=="yes"){
intstand<-read.table(file=file.path(logfile[[1]],"dataframes","IS.txt"),header=TRUE,sep="\t",colClasses = "character");
if(for_mode=="positive"){
if(exists("profileList_pos_copy")){rm("profileList_pos_copy")}
if(exists("links_peaks_pos")){rm("links_peaks_pos")}
load(file=file.path(as.character(logfile[[1]]),"results","profileList_pos_copy"), envir= environment())
load(file=file.path(as.character(logfile[[1]]),"results","links_peaks_pos"), envir= environment())
for_peaks<-which(
profileList_pos_copy[["peaks"]][,"sampleIDs"]==for_file
)
peaks<-profileList_pos_copy[["peaks"]][for_peaks,]
rm(for_peaks)
if(dim(peaks)[1]>0){
peaks<-peaks[order(peaks[,"peakIDs"],decreasing=FALSE),]
for(i in 1:dim(component[[1]])[1]){
collect_comp<-c()
collect_peak<-c()
if(component[[1]][i,3]!="-"){
those<-as.numeric(strsplit(component[[1]][i,3],",")[[1]])
these<-as.numeric(those)
these<-match(these,peaks[,"peakIDs"])
these <- these[!is.na(these)]
if(length(these)){
for(n in 1:length(these) ){
if(peaks[these[n],"links"]!=0){
got<-unlist(links_peaks_pos[[peaks[these[n],"links"]]]$IS)
if(!length(got)) next
got <- got[seq(1,length(got),2)]
collect_comp<-c(collect_comp,got)
collect_peak<-c(collect_peak,rep(these[n],length(got)))
}
}
}
}
if(component[[1]][i,5]!="-"){
those<-as.numeric(strsplit(component[[1]][i,5],",")[[1]])
these<-as.numeric(those)
these<-match(these,peaks[,"peakIDs"])
these <- these[!is.na(these)]
if(length(these)){
for(n in 1:length(these) ){
if(peaks[these[n],"links"]!=0){
got<-unlist(links_peaks_pos[[peaks[these[n],"links"]]]$IS)
if(!length(got)) next
got <- got[seq(1,length(got),2)]
collect_comp<-c(collect_comp,got)
collect_peak<-c(collect_peak,rep(these[n],length(got)))
}
}
}
}
if(component[[1]][i,7]!="-"){
those<-as.numeric(strsplit(component[[1]][i,7],",")[[1]])
these<-as.numeric(those)
these<-match(these,peaks[,"peakIDs"])
these <- these[!is.na(these)]
if(length(these)){
for(n in 1:length(these) ){
if(peaks[these[n],"links"]!=0){
got<-unlist(links_peaks_pos[[peaks[these[n],"links"]]]$IS)
if(!length(got)) next
got <- got[seq(1,length(got),2)]
collect_comp<-c(collect_comp,got)
collect_peak<-c(collect_peak,rep(these[n],length(got)))
}
}
}
}
if(length(collect_comp)>0){
got_comp<-unique(collect_comp)
collect_all<-""
for(n in 1:length(got_comp)){
tar<-strsplit(got_comp[n],"_")[[1]][1:2]
tar_name<-intstand[intstand[,"ID"]==tar[1],"Name"]
tar<-paste(tar_name,tar[2],sep="_")
collect_all<-paste0(c(collect_all,
paste(tar,
paste0(unique(collect_peak[collect_comp==got_comp[n]]),collapse=","),
sep=": "
)),
collapse=" / "
)
}
if(nchar(collect_all)>3){
collect_all<-substr(collect_all,4,nchar(collect_all))
}
component[[1]][i,"ISTD peaks"]<-collect_all
}
}
rm(collect_comp, collect_peak)
}
rm(profileList_pos_copy, links_peaks_pos)
}
if(for_mode=="negative"){
if(exists("profileList_neg_copy")){rm(profileList_neg_copy)}
if(exists("links_peaks_neg")){rm(links_peaks_neg)}
load(file=file.path(as.character(logfile[[1]]),"results","profileList_neg_copy"), envir = environment())
load(file=file.path(as.character(logfile[[1]]),"results","links_peaks_neg"), envir = environment())
for_peaks<-which(
profileList_neg_copy[["peaks"]][,"sampleIDs"]==for_file
)
peaks<-profileList_neg_copy[["peaks"]][for_peaks,]
rm(for_peaks)
if(dim(peaks)[1]>0){
peaks<-peaks[order(peaks[,"peakIDs"],decreasing=FALSE),]
for(i in 1:dim(component[[1]])[1]){
collect_comp<-c()
collect_peak<-c()
if(component[[1]][i,3]!="-"){
those<-as.numeric(strsplit(component[[1]][i,3],",")[[1]])
these<-as.numeric(those)
these<-match(these,peaks[,"peakIDs"])
these <- these[!is.na(these)]
if(length(these)){
for(n in 1:length(these) ){
if(peaks[these[n],"links"]!=0){
got<-unlist(links_peaks_neg[[peaks[these[n],"links"]]]$IS)
if(!length(got)) next
got <- got[seq(1,length(got),2)]
collect_comp<-c(collect_comp,got)
collect_peak<-c(collect_peak,rep(these[n],length(got)))
}
}
}
}
if(component[[1]][i,5]!="-"){
those<-as.numeric(strsplit(component[[1]][i,5],",")[[1]])
these<-as.numeric(those)
these<-match(these,peaks[,"peakIDs"])
these <- these[!is.na(these)]
if(length(these)){
for(n in 1:length(these) ){
if(peaks[these[n],"links"]!=0){
got<-unlist(links_peaks_neg[[peaks[these[n],"links"]]]$IS)
if(!length(got)) next
got <- got[seq(1,length(got),2)]
collect_comp<-c(collect_comp,got)
collect_peak<-c(collect_peak,rep(these[n],length(got)))
}
}
}
}
if(component[[1]][i,7]!="-"){
those<-as.numeric(strsplit(component[[1]][i,7],",")[[1]])
these<-as.numeric(those)
these<-match(these,peaks[,"peakIDs"])
these <- these[!is.na(these)]
if(length(these)){
for(n in 1:length(these) ){
if(peaks[these[n],"links"]!=0){
got<-unlist(links_peaks_neg[[peaks[these[n],"links"]]]$IS)
if(!length(got)) next
got <- got[seq(1,length(got),2)]
collect_comp<-c(collect_comp,got)
collect_peak<-c(collect_peak,rep(these[n],length(got)))
}
}
}
}
if(length(collect_comp)>0){
got_comp<-unique(collect_comp)
collect_all<-""
for(n in 1:length(got_comp)){
tar<-strsplit(got_comp[n],"_")[[1]][1:2]
tar_name<-intstand[intstand[,"ID"]==tar[1],"Name"]
tar<-paste(tar_name,tar[2],sep="_")
collect_all<-paste0(c(collect_all,
paste(tar,
paste0(unique(collect_peak[collect_comp==got_comp[n]]),collapse=","),
sep=": "
)),
collapse=" / "
)
}
if(nchar(collect_all)>3){
collect_all<-substr(collect_all,4,nchar(collect_all))
}
component[[1]][i,"ISTD peaks"]<-collect_all
}
}
rm(collect_comp, collect_peak)
}
rm(profileList_neg_copy, links_peaks_neg)
}
}
if(exists("peaklist", envir = as.environment(".GlobalEnv"))){rm("peaklist", envir = as.environment(".GlobalEnv"))}
if(exists("peaklist")){rm("peaklist")}
load(file = file.path(logfile[[1]], "peaklist", as.character(for_file)), envir = environment());
peaklist <- peaklist[order(peaklist[,10], decreasing = FALSE),]
for(i in 1:dim(component[["Components"]])[1]){
those <- as.numeric(strsplit(component[["Components"]][i, "ID pattern peaks |"], ",")[[1]])
those <- match(those,peaklist[,"peak_ID"])
those <- those[order(peaklist[those,"m/z"])]
those <- those[1]
component[["Components"]][i,"Monois. peak ID |"] <- (peaklist[those,"peak_ID"])
component[["Components"]][i,"Monois. m/z |"] <- (peaklist[those,"m/z_corr"])
component[["Components"]][i,"Monois. RT |"] <- (peaklist[those,"RT_corr"])
component[["Components"]][i,"Monois. int. |"] <- (peaklist[those,"int_corr"])
component[["Components"]][i,"Monois. sample/blind int. ratio"] <- (peaklist[those,"keep_2"])
}
if(logfile$workflow[names(logfile$workflow) == "blind"] == "yes"){
for(i in 1:dim(component[[1]])[1]){
len_tot <- 0
len_blind <- 0
if(component[["Components"]][i, "ID pattern peaks |"] != "-"){
those <- as.numeric(strsplit(component[[1]][i, "ID pattern peaks |"], ",")[[1]])
these <- as.numeric(those)
these <- match(these, peaklist[,"peak_ID"])
len_tot <- (len_tot + length(these))
for(j in 1:length(these)){
if(peaklist[these[j],"keep_2"] < as.numeric(logfile$parameters$blind_threshold)){
those[j] <- paste0(those[j],"*")
len_blind <- (len_blind+1)
}
}
component[["Components"]][i,"ID pattern peaks |"] <- paste0(those,collapse=",")
}
if(component[["Components"]][i,"ID adduct peaks |"] != "-"){
those <- as.numeric(strsplit(component[[1]][i,"ID adduct peaks |"], ",")[[1]])
these <- as.numeric(those)
these <- match(these, peaklist[,"peak_ID"])
len_tot <- (len_tot + length(these))
for(j in 1:length(these)){
if(peaklist[these[j],"keep_2"] < as.numeric(logfile$parameters$blind_threshold)){
those[j] <- paste0(those[j],"*")
len_blind <- (len_blind+1)
}
}
component[["Components"]][i, "ID adduct peaks |"] <- paste0(those,collapse=",")
}
if(component[["Components"]][i,"ID interfering peaks |"]!="-"){
those<-as.numeric(strsplit(component[[1]][i,"ID interfering peaks |"],",")[[1]])
these<-as.numeric(those)
these<-match(these,peaklist[,"peak_ID"])
len_tot<-(len_tot+length(these))
for(j in 1:length(these)){
if(peaklist[these[j],"keep_2"]<as.numeric(logfile$parameters$blind_threshold)){
those[j]<-paste0(those[j],"*")
len_blind<-(len_blind+1)
}
}
component[["Components"]][i,"ID interfering peaks |"]<-paste0(those,collapse=",")
}
component[["Components"]][i,"Total peak number"]<-as.character(len_tot)
component[["Components"]][i,"Blind peak number"]<-as.character(len_blind)
}
rm(peaklist)
}
if(
(logfile$parameters$do_atom_bounds_components == "TRUE") &
(length(logfile$parameters$atom_bounds_components))
){
num_elem <- length(mute(logfile$parameters$atom_bounds_components))
num_atom <- matrix(nrow = dim(component[["Components"]])[1], ncol = num_elem, 0)
for(i in 1:dim(component[["Components"]])[1]){
use_peaks <- as.numeric(strsplit(gsub("*", "", component[["Components"]][i, "ID pattern peaks |"], fixed = TRUE), ",")[[1]])
if(component[["Components"]][i, "z"] == "-"){
use_charges <- 1
}else{
use_charges <- as.numeric(strsplit(component[["Components"]][i, "z"], "/")[[1]])
use_charges <- use_charges[order(use_charges, decreasing = FALSE)]
}
if(length(component[["pattern peak list"]]) > 1){
these <- match(use_peaks, component[["pattern peak list"]][,"peak ID"])
use_masses <- component[["pattern peak list"]][these, "m/z_corr"]
use_intensities <- component[["pattern peak list"]][these, "int_corr"]
use_RT <- mean(component[["pattern peak list"]][these, "RT_corr"])
}else{
these <- match(use_peaks, component[["adduct peak list"]][,"peak ID"])
use_masses <- component[["adduct peak list"]][these, "m/z_corr"]
use_intensities <- component[["adduct peak list"]][these, "int_corr"]
use_RT <- mean(component[["adduct peak list"]][these, "RT_corr"])
}
if(do_LOD){
with_model <- which(names(LOD_splined) == paste("LOD_", for_file, sep = ""))
if(length(with_model) > 0){
use_cutint <- 10^(predict(LOD_splined[[with_model]], use_RT)$y)
}else{
cat("\n Missing LOD model; using default intensity threshold. Debug?")
use_cutint <- as.numeric(logfile$parameters$IS_intcut)
}
}else{
use_cutint <- as.numeric(logfile$parameters$IS_intcut)
}
atom_counts <- enviMass:::QZTGM3353W(
masses = use_masses,
intensities = use_intensities,
elements = logfile$parameters$atom_bounds_components,
dmz = rep(30, num_elem),
ppm = TRUE,
charges = use_charges,
isotopes,
int_cut = use_cutint,
inttol = (as.numeric(logfile$parameter$isotop_inttol) / 100),
use_C = FALSE,
must_peak = TRUE
)
num_atom[i,] <- do.call(pmax, data.frame(t(atom_counts)))
}
named <- colnames(component[["Components"]])
component[["Components"]] <- cbind(component[["Components"]], num_atom)
names(component[["Components"]]) <- c(named, paste0("max_atom_", logfile$parameters$atom_bounds_components))
}
save(component, file = (file.path(logfile[[1]],"results","componentization","components",paste(for_file))))
return(NULL);
}
