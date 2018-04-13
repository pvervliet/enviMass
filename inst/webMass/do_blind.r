measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
measurements<-measurements[measurements[,"include"]=="TRUE",]
IDs<-(measurements[,"ID"])
incl<-(measurements[,"include"])
filetypus<-(measurements[,"Type"])
ionmode<-(measurements[,"Mode"])
atdate<-(measurements[,"Date"])
attime<-(measurements[,"Time"])
attime2<-as.difftime(attime);
atdate<-as.Date(atdate, tz="GMT");
sampleID<-(measurements[,"ID"])
old_samplewise<-(measurements[,"blind"])
new_samplewise<-old_samplewise
ord<-order(as.numeric(atdate),as.numeric(attime2),filetypus,sampleID);
ppm<-logfile$parameters$blind_ppm
dmz<-as.numeric(logfile$parameters$blind_dmz)
dRT<-as.numeric(logfile$parameters$blind_drt)
int_ratio<-as.numeric(logfile$parameters$blind_threshold)
if(FALSE){
ppm<-TRUE
dmz<-3
dRT<-30
int_ratio<-10
}
for(i in 1:length(IDs)){
if(incl[i]=="FALSE"){next}
if(filetypus[i]=="blank"){next}
if(old_samplewise[i]=="TRUE"){next}
load(file=file.path(logfile[[1]],"peaklist",as.character(IDs[i])),envir=as.environment(".GlobalEnv"),verbose=FALSE);
peaklist[,"keep_2"]<-Inf
save(peaklist,file=file.path(logfile[[1]],"peaklist",as.character(IDs[i])))
rm(peaklist)
}
if((logfile$parameters$subtract_pos_bydate=="TRUE") || (logfile$parameters$subtract_neg_bydate=="TRUE")){
blank_ID_last<-"FALSE"
for(i in 1:length(ord)){
#if(sampleID[ord[i]]=="167"){stop()}
if((logfile$parameters$subtract_pos_bydate=="FALSE") & (ionmode[ord[i]]=="positive")){next}
if((logfile$parameters$subtract_neg_bydate=="FALSE") & (ionmode[ord[i]]=="negative")){next}
if(old_samplewise[ord[i]]=="TRUE"){next}
if(incl[ord[i]]=="FALSE"){next}
if(filetypus[ord[i]]=="sample"){
sam_ID<-sampleID[ord[i]]
found_blank<-FALSE
if(i>1){
for(j in (i-1):1){
if((filetypus[ord[j]]=="blank") & (ionmode[ord[i]]==ionmode[ord[j]])){
blank_ID<-sampleID[ord[j]]
found_blank<-TRUE
break;
}
}
}
if(!found_blank){
next;
}
if(blank_ID!=blank_ID_last){
load(file=file.path(logfile[[1]],"peaklist",as.character(blank_ID)),verbose=FALSE);
peaks_blank<-peaklist[,c("m/z_corr","int_corr","RT_corr")]
rm(peaklist)
}
load(file=file.path(logfile[[1]],"peaklist",as.character(sam_ID)),verbose=FALSE);
peaks_sample<-peaklist[,c("m/z_corr","int_corr","RT_corr")]
getit <- BYRIB9287R(
peaklist=peaks_blank,
mz=peaks_sample[,"m/z_corr"],
dmz=(dmz*2),
ppm=ppm,
RT=peaks_sample[,"RT_corr"],
dRT=dRT,
onlymax=TRUE,
int_ratio=int_ratio,
int=peaks_sample[,"int_corr"],
get_matches=FALSE,
get_ratio=TRUE
)
which_affected<-(peaklist[,"keep_2"]>getit)
peaklist[which_affected,"keep_2"]<-getit[which_affected]
save(peaklist,file=file.path(logfile[[1]],"peaklist",as.character(sam_ID)))
blank_ID_last<-blank_ID
cat(paste("\n",
round(sum(peaklist[,"keep_2"]==Inf)/length(peaklist[,1])*100,digits=1),
" % of ",
length(peaklist[,1]),
" peaks not found in blank/blind files (files ",
sam_ID," vs. ",blank_ID,", ",ionmode[ord[i]],", by date & time)."
,sep=""))
rm(peaklist);
new_samplewise[ord[i]]<-"TRUE"
}
}
}
if( (logfile$parameters$subtract_pos_byfile=="TRUE") & any(logfile$Positive_subtraction_files!="FALSE") ){
selec_pos<-logfile$Positive_subtraction_files
selec_pos<-selec_pos[selec_pos!="FALSE"]
for(i in 1:length(IDs)){
if(any(measurements[,"ID"]==IDs[i])){
if( filetypus[measurements[,"ID"]==IDs[i]]=="sample" &  ionmode[measurements[,"ID"]==IDs[i]]=="positive" ){
if(old_samplewise[measurements[,"ID"]==IDs[i]]=="TRUE"){next}
if(incl[measurements[,"ID"]==IDs[i]]=="FALSE"){next}
load(file=file.path(logfile[[1]],"peaklist",as.character(IDs[i])),verbose=FALSE);
sam_peaklist<-peaklist;rm(peaklist);
for(j in 1:length(selec_pos)){
subID<-strsplit(selec_pos[j]," - ")[[1]][1]
load(file=file.path(logfile[[1]],"peaklist",as.character(subID)),verbose=FALSE)
peaks_blank<-peaklist[,c("m/z_corr","int_corr","RT_corr")];rm(peaklist);
getit <- BYRIB9287R(
peaklist=peaks_blank,
mz=sam_peaklist[,"m/z_corr"],
dmz=(dmz*2),
ppm=ppm,
RT=sam_peaklist[,"RT_corr"],
dRT=dRT,
onlymax=TRUE,
int_ratio=int_ratio,
int=sam_peaklist[,"int_corr"],
get_matches=FALSE,
get_ratio=TRUE
)
which_affected<-(sam_peaklist[,"keep_2"]>getit)
sam_peaklist[which_affected,"keep_2"]<-getit[which_affected]
rm(peaks_blank)
}
peaklist<-sam_peaklist
save(peaklist,file=file.path(logfile[[1]],"peaklist",as.character(IDs[i])))
cat(paste("\n",
round(sum(peaklist[,colnames(peaklist)=="keep_2"]==Inf)/length(peaklist[,1])*100,digits=1),
" % of ",
length(peaklist[,1]),
" peaks not found in blank/blind file(s) (selective, file ",
as.character(IDs[i]),"). "
,sep=""))
rm(peaklist,sam_peaklist);
new_samplewise[measurements[,"ID"]==IDs[i]]<-"TRUE"
}
}else{
cat("\n Orphaned peaklist detected - from an older workflow run?")
}
}
}
if( (logfile$parameters$subtract_neg_byfile=="TRUE") & any(logfile$Negative_subtraction_files!="FALSE") ){
selec_neg<-logfile$Negative_subtraction_files
selec_neg<-selec_neg[selec_neg!="FALSE"]
for(i in 1:length(IDs)){
if(any(measurements[,"ID"]==IDs[i])){
if(filetypus[measurements[,"ID"]==IDs[i]]=="sample" &  ionmode[measurements[,"ID"]==IDs[i]]=="negative"){
if(old_samplewise[measurements[,"ID"]==IDs[i]]=="TRUE"){next}
if(incl[measurements[,"ID"]==IDs[i]]=="FALSE"){next}
load(file=file.path(logfile[[1]],"peaklist",as.character(IDs[i])),verbose=FALSE);
sam_peaklist<-peaklist;rm(peaklist);
for(j in 1:length(selec_neg)){
subID<-strsplit(selec_neg[j]," - ")[[1]][1]
load(file=file.path(logfile[[1]],"peaklist",as.character(subID)),verbose=FALSE)
peaks_blank<-peaklist[,c("m/z_corr","int_corr","RT_corr")];rm(peaklist);
getit <- BYRIB9287R(
peaklist=peaks_blank,
mz=sam_peaklist[,"m/z_corr"],
dmz=(dmz*2),
ppm=ppm,
RT=sam_peaklist[,"RT_corr"],
dRT=dRT,
onlymax=TRUE,
int_ratio=int_ratio,
int=sam_peaklist[,"int_corr"],
get_matches=FALSE,
get_ratio=TRUE
)
which_affected<-(sam_peaklist[,"keep_2"]>getit)
sam_peaklist[which_affected,"keep_2"]<-getit[which_affected]
rm(peaks_blank)
}
peaklist<-sam_peaklist
save(peaklist,file=file.path(logfile[[1]],"peaklist",as.character(IDs[i])))
cat(paste("\n",
round(sum(peaklist[,colnames(peaklist)=="keep_2"]==Inf)/length(peaklist[,1])*100,digits=1),
" % of ",
length(peaklist[,1]),
" peaks not found in blank/blind file(s) (selective, file ",
as.character(IDs[i]),"). "
,sep=""))
rm(peaklist,sam_peaklist);
new_samplewise[measurements[,"ID"]==IDs[i]]<-"TRUE"
}
}else{
cat("\n Orphaned peaklist detected - from an older workflow run?")
}
}
}
measurements[,"blind"]<-new_samplewise
write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
rm(measurements)
