if( file.exists(file.path(logfile[[1]],"results","profileList_pos")) ){
if(any(objects(envir=as.environment(".GlobalEnv"))=="profileList_pos")){rm(profileList_pos,envir=as.environment(".GlobalEnv"))}
if(any(objects()=="profileList_pos")){rm(profileList_pos)}
load(file.path(as.character(logfile[[1]]),"results","profileList_pos"),envir=as.environment(".GlobalEnv"))
lags <- as.numeric(strsplit(logfile$parameters$trend_lags,",")[[1]])
if(logfile$parameters$trend_blind=="yes"){ blindsub <- TRUE }else{ blindsub <- FALSE}
blindfold <- as.numeric(logfile$parameters$blind_threshold)
from <- FALSE
to <- FALSE
threshold <- as.numeric(logfile$parameters$trend_threshold)
notrend <- FALSE
omit_trend <- FALSE
if(!profileList_pos[[ "state"]][[3]]){stop("profileList_pos not profiled; aborted.")}
if(!is.numeric(lags)){stop("lags argument must be numeric; aborted.")}
if(!is.logical(notrend)){stop("notrend must be logical.")}
if(!is.logical(omit_trend)){stop("notrend must be logical.")}
atPOSIX<-profileList_pos[["datetime"]];
sampletype<-profileList_pos[["type"]];
sampleID<-profileList_pos[["sampleID"]];
keep<-((sampletype=="sample")|(sampletype=="blank"))
atPOSIX<-atPOSIX[keep]
sampletype<-sampletype[keep]
sampleID<-sampleID[keep]
atdate<-c();attime<-c();
for(i in 1:length(atPOSIX)){
atdate<-c(atdate, strsplit(atPOSIX[i]," ")[[1]][1]);
attime<-c(attime, strsplit(atPOSIX[i]," ")[[1]][2]);
}
attime<-as.difftime(attime);
atdate<-as.Date(atdate, tz="GMT");
ord<-order(as.numeric(atdate),as.numeric(attime),sampleID);
atPOSIXsort<-atPOSIX[ord];
atdate<-atdate[ord];
attime<-attime[ord];
sampleID<-sampleID[ord];
sampletype<-sampletype[ord];
timeset<-matrix(nrow=length(atPOSIX),ncol=(5+(4*length(lags))),0);
for(i in 1:length(sampleID)){
if(sampletype[i]=="sample"){
timeset[i,2]<-as.numeric(sampleID[i]);
}
if(sampletype[i]=="blank"){
timeset[i,3]<-as.numeric(sampleID[i]);
}
}
numtime<-(as.numeric(atdate)+as.numeric(attime/(24*60*60)))
colnames(timeset)<-c("above blank?","sampleID","blankID","sample_int","blank_int",rep("lag_int",length(lags)),rep("del_int",length(lags)),rep("max_time",length(lags)),rep("blind_int",length(lags)))
leng<-max(seq(1,length(timeset[,1]),1)[timeset[,2]!=0])
latestID<-timeset[leng,2][[1]]
if(!omit_trend){
if(any(lags>(max(numtime)-min(numtime)+1))){
lags<-lags[lags<=(max(numtime)-min(numtime)+1)]
cat("WARNING: at least one lag longer than covered time period - omitted!\n")
if(length(lags)==0){
stop("...no lags left; aborted.")
}
}
}
if(max(diff(numtime[timeset[,2]!=0],lag=1))>max(lags)){
warning("\n At least one gap in sample time series larger than largest lag!")
}
profileList_pos[["index_prof"]][,5:7]<-0;
along <- seq(1,dim(profileList_pos[["index_prof"]])[1],1)
size_MB <- (as.numeric(object.size(profileList_pos[["peaks"]])) / 1048576)
for_split <- (size_MB / 2)
for_split <- round( dim(profileList_pos[["index_prof"]])[1] / for_split )
along <- split(along, ceiling(seq_along(along) / for_split))
clus_peaks <- list()
clus_profs <- list()
for(i in 1:length(along)){
clus_peaks[[i]] <- profileList_pos[["peaks"]][
(profileList_pos[["index_prof"]][along[[i]][1], "start_ID"]):(profileList_pos[["index_prof"]][along[[i]][length(along[[i]])], "end_ID"])
,,drop=FALSE]
clus_profs[[i]] <- profileList_pos[["index_prof"]][along[[i]][1]:along[[i]][length(along[[i]])],, drop = FALSE]
start_at <- (clus_profs[[i]][1, "start_ID"][[1]] - 1)
clus_profs[[i]][,1] <- (clus_profs[[i]][,1] - start_at)
clus_profs[[i]][,2] <- (clus_profs[[i]][,2] - start_at)
}
clusterEvalQ(cl = clus,{rm(list=ls()); gc(verbose=FALSE); NULL})
clusterExport(cl = clus,
varlist = c("timeset", "lags", "threshold", "notrend", "omit_trend", "blindsub", "blindfold", "numtime", "latestID", "leng"),
envir = environment())
cluster_results <- clusterMap(
cl = clus,
fun = enviMass:::SOPHZ8293U,
peaks = clus_peaks,
index_prof = clus_profs,
RECYCLE = TRUE,
SIMPLIFY = FALSE,
USE.NAMES = FALSE,
.scheduling = c("dynamic")
)
clusterEvalQ(cl = clus,{rm(list = ls()); gc(verbose = FALSE); NULL})
for(i in 1:length(cluster_results)){
profileList_pos[["index_prof"]][
cluster_results[[i]][1,"profile_ID"]:cluster_results[[i]][dim(cluster_results[[i]])[1],"profile_ID"]
,c("deltaint_newest", "deltaint_global", "absolute_mean_dev", "newest_intensity")] <-
cluster_results[[i]][,c("deltaint_newest", "deltaint_global", "absolute_mean_dev", "newest_intensity")]
}
rm(cluster_results, clus_peaks, clus_profs)
profileList_pos[["state"]][[4]] <- TRUE;
profileList_pos <<- profileList_pos
save(profileList_pos,file=file.path(as.character(logfile[[1]]),"results","profileList_pos"),compress=FALSE);
png(filename = file.path(as.character(logfile[[1]]),"pics","boxprofile_pos"), width = 800, bg = "white")
enviMass::KXNZN5280X(profileList_pos)
dev.off()
expr4p <- list(src=file.path(logfile[[1]],"pics","boxprofile_pos"))
output$boxprofile<-renderImage(expr4p, deleteFile = FALSE)
if(isolate(input$Ion_mode)=="positive"){
profileList <<- profileList_pos;
}
}
if( file.exists(file.path(logfile[[1]],"results","profileList_neg")) ){
if(any(objects(envir=as.environment(".GlobalEnv"))=="profileList_neg")){rm(profileList_neg,envir=as.environment(".GlobalEnv"))}
if(any(objects()=="profileList_neg")){rm(profileList_neg)}
load(file.path(as.character(logfile[[1]]),"results","profileList_neg"),envir=as.environment(".GlobalEnv"))
lags <- as.numeric(strsplit(logfile$parameters$trend_lags,",")[[1]])
if(logfile$parameters$trend_blind=="yes"){ blindsub <- TRUE }else{ blindsub <- FALSE}
blindfold <- as.numeric(logfile$parameters$blind_threshold)
from <- FALSE
to <- FALSE
threshold <- as.numeric(logfile$parameters$trend_threshold)
notrend <- FALSE
omit_trend <- FALSE
if(!profileList_neg[[ "state"]][[3]]){stop("profileList_neg not profiled; aborted.")}
if(!is.numeric(lags)){stop("lags argument must be numeric; aborted.")}
if(!is.logical(notrend)){stop("notrend must be logical.")}
if(!is.logical(omit_trend)){stop("notrend must be logical.")}
atPOSIX<-profileList_neg[["datetime"]];
sampletype<-profileList_neg[["type"]];
sampleID<-profileList_neg[["sampleID"]];
keep<-((sampletype=="sample")|(sampletype=="blank"))
atPOSIX<-atPOSIX[keep]
sampletype<-sampletype[keep]
sampleID<-sampleID[keep]
atdate<-c();attime<-c();
for(i in 1:length(atPOSIX)){
atdate<-c(atdate, strsplit(atPOSIX[i]," ")[[1]][1]);
attime<-c(attime, strsplit(atPOSIX[i]," ")[[1]][2]);
}
attime<-as.difftime(attime);
atdate<-as.Date(atdate, tz="GMT");
ord<-order(as.numeric(atdate),as.numeric(attime),sampleID);
atPOSIXsort<-atPOSIX[ord];
atdate<-atdate[ord];
attime<-attime[ord];
sampleID<-sampleID[ord];
sampletype<-sampletype[ord];
timeset<-matrix(nrow=length(atPOSIX),ncol=(5+(4*length(lags))),0);
for(i in 1:length(sampleID)){
if(sampletype[i]=="sample"){
timeset[i,2]<-as.numeric(sampleID[i]);
}
if(sampletype[i]=="blank"){
timeset[i,3]<-as.numeric(sampleID[i]);
}
}
numtime<-(as.numeric(atdate)+as.numeric(attime/(24*60*60)))
colnames(timeset)<-c("above blank?","sampleID","blankID","sample_int","blank_int",rep("lag_int",length(lags)),rep("del_int",length(lags)),rep("max_time",length(lags)),rep("blind_int",length(lags)))
leng<-max(seq(1,length(timeset[,1]),1)[timeset[,2]!=0])
latestID<-timeset[leng,2][[1]]
if(!omit_trend){
if(any(lags>(max(numtime)-min(numtime)+1))){
lags<-lags[lags<=(max(numtime)-min(numtime)+1)]
cat("WARNING: at least one lag longer than covered time period - omitted!\n")
if(length(lags)==0){
stop("...no lags left; aborted.")
}
}
}
if(max(diff(numtime[timeset[,2]!=0],lag=1))>max(lags)){
warning("\n At least one gap in sample time series larger than largest lag!")
}
profileList_neg[["index_prof"]][,5:7]<-0;
along <- seq(1,dim(profileList_neg[["index_prof"]])[1],1)
size_MB <- (as.numeric(object.size(profileList_neg[["peaks"]])) / 1048576)
for_split <- (size_MB / 20)
for_split <- round( dim(profileList_neg[["index_prof"]])[1] / for_split )
along <- split(along, ceiling(seq_along(along)/for_split))
clus_peaks<-list()
clus_profs<-list()
for(i in 1:length(along)){
clus_peaks[[i]]<-profileList_neg[["peaks"]][
(profileList_neg[["index_prof"]][along[[i]][1],"start_ID"]):(profileList_neg[["index_prof"]][along[[i]][length(along[[i]])],"end_ID"])
,,drop=FALSE]
clus_profs[[i]]<-profileList_neg[["index_prof"]][along[[i]][1]:along[[i]][length(along[[i]])],,drop=FALSE]
start_at<-(clus_profs[[i]][1,"start_ID"][[1]]-1)
clus_profs[[i]][,1]<-(clus_profs[[i]][,1]-start_at)
clus_profs[[i]][,2]<-(clus_profs[[i]][,2]-start_at)
}
clusterEvalQ(cl = clus,{rm(list=ls()); gc(verbose=FALSE); NULL})
clusterExport(cl = clus,
varlist = c("timeset", "lags", "threshold", "notrend", "omit_trend", "blindsub", "blindfold", "numtime", "latestID", "leng"),
envir = environment())
cluster_results <- clusterMap(
cl = clus,
fun = enviMass:::SOPHZ8293U,
peaks = clus_peaks,
index_prof = clus_profs,
RECYCLE = TRUE,
SIMPLIFY = FALSE,
USE.NAMES = FALSE,
.scheduling = c("dynamic")
)
clusterEvalQ(cl = clus,{rm(list=ls()); gc(verbose=FALSE); NULL})
for(i in 1:length(cluster_results)){
profileList_neg[["index_prof"]][
cluster_results[[i]][1,"profile_ID"]:cluster_results[[i]][dim(cluster_results[[i]])[1],"profile_ID"]
,c("deltaint_newest", "deltaint_global", "absolute_mean_dev", "newest_intensity")] <-
cluster_results[[i]][,c("deltaint_newest", "deltaint_global", "absolute_mean_dev", "newest_intensity")]
}
rm(cluster_results, clus_peaks, clus_profs)
profileList_neg[["state"]][[4]] <- TRUE;
profileList_neg<<-profileList_neg
save(profileList_neg,file=file.path(as.character(logfile[[1]]),"results","profileList_neg"),compress=FALSE);
png(filename = file.path(as.character(logfile[[1]]),"pics","boxprofile_neg"), width = 800, bg = "white")
enviMass::KXNZN5280X(profileList_neg)
dev.off()
expr4p<-list(src=file.path(logfile[[1]],"pics","boxprofile_neg"))
output$boxprofile<-renderImage(expr4p, deleteFile = FALSE)
if(isolate(input$Ion_mode)=="negative"){
profileList<<-profileList_neg;
}
}
