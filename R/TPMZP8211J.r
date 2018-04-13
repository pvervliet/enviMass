TPMZP8211J <- function(
profileList,
from=FALSE,
to=FALSE,
progbar=FALSE,
blindsub=TRUE,
blindfold=100,
lags=c(5,14),
threshold=3,
notrend=FALSE,
omit_trend=FALSE
){
if(!profileList[[ "state"]][[3]]){stop("profileList not profiled; aborted.")}
if(!from){m=1}else{m=from}
if(!to){n=length(profileList[["index_prof"]][,1])}else{n=to}
if(blindsub!=FALSE){if(!is.numeric(blindfold) || (blindfold<0)){stop("Invalid blindfold argument; aborted.")}}
if(blindsub!=FALSE){subit=1;subrat=blindfold;}else{subit=2;subrat=0;}
if(!is.numeric(lags)){stop("lags argument must be numeric; aborted.")}
if(!is.logical(notrend)){stop("notrend must be logical.")}
if(!is.logical(omit_trend)){stop("notrend must be logical.")}
atPOSIX<-profileList[["datetime"]];
sampletype<-profileList[["type"]];
sampleID<-profileList[["sampleID"]];
keep<-((sampletype=="sample")|(sampletype=="blank"))
atPOSIX<-atPOSIX[keep]
sampletype<-sampletype[keep]
sampleID<-sampleID[keep]
atdate<-c();
attime<-c();
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
profileList[["index_prof"]][,5:7]<-0;
if(progbar==TRUE){  prog<-winProgressBar("Extract intensity differences...",min=m,max=n);
setWinProgressBar(prog, 0, title = "Extract intensity differences...", label = NULL);}
for(k in m:n){
if(progbar==TRUE){setWinProgressBar(prog, k, title = "Extract intensity differences...", label = NULL)}
if(profileList[["index_prof"]][k,"number_peaks_total"]>1){
timeset[,4:length(timeset[1,])]<-0;
timeset[,c(4,5)] <-.Call("_enviMass_fill_timeset",
as.numeric(timeset),
as.numeric(profileList[["peaks"]][(profileList[["index_prof"]][k,"start_ID"]:profileList[[7]][k,"end_ID"]),"sampleIDs"]),
as.numeric(profileList[["peaks"]][(profileList[["index_prof"]][k,"start_ID"]:profileList[[7]][k,"end_ID"]),"intensity"]),
as.integer(length(timeset[,1])),
PACKAGE="enviMass"
)
if(any(timeset[,4]>0)){
what<-1
if(!omit_trend){
that<-.Call("_enviMass_meandel",
as.numeric(timeset),
as.integer(subit),
as.numeric(subrat),
as.numeric(numtime),
as.integer(what),
as.numeric(lags),
as.numeric(threshold),
as.integer(notrend),
PACKAGE="enviMass"
)
}
if( (what!=1) & (!omit_trend)){
that<-matrix(nrow=length(atPOSIX),ncol=(5+(4*length(lags))),that);
colnames(that)<-c("above blank?","sampleID","blankID","sample_int","blank_int",rep("lag_int",length(lags)),rep("del_int",length(lags)),rep("max_time",length(lags)),rep("blind_int",length(lags)))
plot.new();
plot.window(xlim=c(min(numtime),max(numtime)),ylim=c(min(that[that[,2]!=0,4]),max(that[that[,2]!=0,4])));
box();axis(1);axis(2);
title(xlab="Time",ylab="Intensity")
points(numtime[that[,2]!=0],that[that[,2]!=0,4],col="red",type="l");
for(i in 1:length(lags)){
points(numtime[that[,2]!=0],that[that[,2]!=0,(5+i)],col="darkgrey",type="l");
points(that[that[,(5+i+(length(lags)*2))]!=0,(5+i+(length(lags)*2))],that[that[,(5+i+(length(lags)*2))]!=0,(5+i+length(lags))],col="blue",pch=19);
}
stop(" YOU wanted the smoothed series...\n")
}else{
if(!omit_trend){
profileList[["index_prof"]][k,"deltaint_newest"]<-max(that[4,]);
profileList[["index_prof"]][k,"deltaint_global"]<-max(that[5,]);
profileList[["index_prof"]][k,"absolute_mean_dev"]<-max(that[3,]);
}
if(any(timeset[,5]>0)){
profileList[["index_prof"]][k,"in_blind?"]<-1
profileList[["index_prof"]][k,"number_peaks_blind"]<-length(timeset[timeset[,5]!=0,5])
profileList[["index_prof"]][k,"mean_int_blind"]<-mean(timeset[timeset[,5]!=0,5])
profileList[["index_prof"]][k,"max_int_blind"]<-max(timeset[timeset[,5]!=0,5])
}else{
profileList[["index_prof"]][k,"in_blind?"]<-0
profileList[["index_prof"]][k,"number_peaks_blind"]<-0
profileList[["index_prof"]][k,"mean_int_blind"]<-0
profileList[["index_prof"]][k,"max_int_blind"]<-0
}
profileList[["index_prof"]][k,"number_peaks_sample"]<-length(timeset[timeset[,4]!=0,4])
profileList[["index_prof"]][k,"mean_int_sample"]<-mean(timeset[timeset[,4]!=0,4])
profileList[["index_prof"]][k,"max_int_sample"]<-max(timeset[timeset[,4]!=0,4])
}
}else{
profileList[["index_prof"]][k,"in_blind?"]<-1
profileList[["index_prof"]][k,"number_peaks_sample"]<-0
profileList[["index_prof"]][k,"number_peaks_blind"]<-length(timeset[timeset[,5]!=0,5])
profileList[["index_prof"]][k,"mean_int_sample"]<-0
profileList[["index_prof"]][k,"max_int_sample"]<-0
profileList[["index_prof"]][k,"mean_int_blind"]<-mean(timeset[timeset[,5]!=0,5])
profileList[["index_prof"]][k,"max_int_blind"]<-max(timeset[timeset[,5]!=0,5])
}
profileList[["index_prof"]][k,"newest_intensity"]<-timeset[leng,4]
}else{
profileList[["index_prof"]][k,"absolute_mean_dev"]<-0;
if( any( timeset[,3]== (profileList[[2]][profileList[["index_prof"]][k,1],6]) ) ){
profileList[["index_prof"]][k,"in_blind?"]<-1
profileList[["index_prof"]][k,"number_peaks_sample"]<-0
profileList[["index_prof"]][k,"number_peaks_blind"]<-1
profileList[["index_prof"]][k,"mean_int_sample"]<-0
profileList[["index_prof"]][k,"max_int_sample"]<-0
profileList[["index_prof"]][k,"mean_int_blind"]<-(profileList[[2]][profileList[[7]][k,1],2])
profileList[["index_prof"]][k,"max_int_blind"]<-(profileList[[2]][profileList[[7]][k,1],2])
}else{
profileList[["index_prof"]][k,"in_blind?"]<-0
profileList[["index_prof"]][k,"deltaint_global"]<-(profileList[[2]][profileList[["index_prof"]][k,1],2])
if(any(profileList[[2]][profileList[["index_prof"]][k,1],6]==latestID)){
profileList[["index_prof"]][k,"deltaint_newest"]<-(profileList[["peaks"]][profileList[["index_prof"]][k,1],2])
profileList[["index_prof"]][k,"newest_intensity"]<-(profileList[["peaks"]][profileList[["index_prof"]][k,1],2])
}
profileList[["index_prof"]][k,"number_peaks_sample"]<-1
profileList[["index_prof"]][k,"number_peaks_blind"]<-0
profileList[["index_prof"]][k,"mean_int_sample"]<-(profileList[["peaks"]][profileList[["index_prof"]][k,1],2])
profileList[["index_prof"]][k,"max_int_sample"]<-(profileList[["peaks"]][profileList[["index_prof"]][k,1],2])
profileList[["index_prof"]][k,"mean_int_blind"]<-0
profileList[["index_prof"]][k,"max_int_blind"]<-0
}
}
}
if(progbar==TRUE){close(prog);}
profileList[[1]][[4]]<-TRUE;
return(profileList)
}
