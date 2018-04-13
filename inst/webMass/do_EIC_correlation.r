cat("EIC correlations: ")
measurements <- read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
EICor_plot_it <- FALSE
do_cor <- TRUE
for(b in 1:length(measurements[,"ID"])){
if(
(measurements[b,"include"]=="TRUE") &
(measurements[b,"EIC_correlation"]=="FALSE")
){
if( (mute(logfile$parameters$prof_select=="TRUE")) & (measurements$profiled[b]=="FALSE") ){
cat("\n Skip file - not included in profile building.");next;
}
cat(paste("\n Doing ",as.character(b)," of ",as.character(length(measurements[,"ID"]))," files: ",sep=""))
cat("loading - ")
for_file<-measurements[b,1]
if( file.exists(file.path(logfile[[1]],"results","componentization","EIC_corr",for_file) ) ){
file.remove(file.path(logfile[[1]],"results","componentization","EIC_corr",for_file) )
}
load(file=file.path(logfile[[1]],"peaklist",as.character(for_file)));
load(file=file.path(logfile[[1]],"MSlist",as.character(for_file)));
cat("pairing - ")
ord <- order(peaklist[,"RT_corr"],decreasing=FALSE)
peaklist <- peaklist[ord,,drop=FALSE]
use_peak <- MSlist[[7]][peaklist[,10],3]>=as.numeric(logfile$parameters$EICor_minpeaks)
at <- 1
paired <- matrix(ncol=4,nrow=1E7,0)
len <- dim(paired)[1]
for(i in 1:(length(peaklist[,1])-1)){
if(!use_peak[i]){next}
for(j in (i+1):length(peaklist[,1])){
if(!use_peak[j]){next}
RTdif<-(peaklist[j,14]-peaklist[i,14])
if(RTdif<=as.numeric(logfile$parameters$EICor_delRT)){
paired[at,1]<-i
paired[at,2]<-j
paired[at,3]<-RTdif
at<-(at+1)
if(at>len){
cat(".")
paired2<-matrix(ncol=4,nrow=1E7,0)
paired<-rbind(paired,paired2)
len<-dim(paired)[1]
}
}else{
break;
}
}
}
paired<-paired[1:(at-1),,drop=FALSE]
if(all(paired[,1]==0)){
cat("nothing found; aborted.");next;
}
cat("correlating - ")
skipped<-0
for(i in 1:length(paired[,1])){
if(i>1){
if(paired[i,1]!=paired[i-1,1]){
PEAK_ID_1<-peaklist[paired[i,1],"peak_ID"]
start_1<-MSlist[[7]][PEAK_ID_1,1]
end_1<-MSlist[[7]][PEAK_ID_1,2]
del_1<-MSlist[[7]][PEAK_ID_1,3]
}
}else{
PEAK_ID_1<-peaklist[paired[i,1],10]
start_1<-MSlist[[7]][PEAK_ID_1,1]
end_1<-MSlist[[7]][PEAK_ID_1,2]
del_1<-MSlist[[7]][PEAK_ID_1,3]
}
PEAK_ID_2<-peaklist[paired[i,2],10]
start_2<-MSlist[[7]][PEAK_ID_2,1]
end_2<-MSlist[[7]][PEAK_ID_2,2]
del_2<-MSlist[[7]][PEAK_ID_2,3]
if(del_1<=del_2){
those<-match(MSlist[[4]][[2]][start_1:end_1,3],MSlist[[4]][[2]][start_2:end_2,3])
if(sum(!is.na(those))<as.numeric(logfile$parameters$EICor_minpeaks)){
paired[i,4]<-(-2000)
skipped<-(skipped+1);
next;
}
if(do_cor){
paired[i,4]<-cor(
MSlist[[4]][[2]][start_1:end_1,2][!is.na(those)],
MSlist[[4]][[2]][start_2:end_2,2][those[!is.na(those)]],
method="spearman"
)
}
}else{
those<-match(MSlist[[4]][[2]][start_2:end_2,3],MSlist[[4]][[2]][start_1:end_1,3])
if(sum(!is.na(those))<as.numeric(logfile$parameters$EICor_minpeaks)){
paired[i,4]<-(-2000)
skipped<-(skipped+1);
next;
}
if(do_cor){
paired[i,4]<-cor(
MSlist[[4]][[2]][start_2:end_2,2][!is.na(those)],
MSlist[[4]][[2]][start_1:end_1,2][those[!is.na(those)]],
method="spearman"
)
}
}
if(EICor_plot_it){
if(paired[i,4]>=EICor_mincor){
split.screen(c(2,1))
xlim<-c(
min(c(MSlist[[4]][[2]][start_1:end_1,3],MSlist[[4]][[2]][start_2:end_2,3])),
max(c(MSlist[[4]][[2]][start_1:end_1,3],MSlist[[4]][[2]][start_2:end_2,3]))
)
par(mar=c(3.5,3.5,.1,.1))
screen(1)
plot(
MSlist[[4]][[2]][start_1:end_1,3],
MSlist[[4]][[2]][start_1:end_1,2],
type="h",xlim=xlim,ylab="Intensity",xlab=""
)
screen(2)
par(mar=c(3.5,3.5,.1,.1))
plot(
MSlist[[4]][[2]][start_2:end_2,3],
MSlist[[4]][[2]][start_2:end_2,2],
type="h",xlim=xlim,ylab="Intensity",xlab="RT"
)
close.screen(n, all.screens = TRUE)
Sys.sleep(.3)
}
}
}
cat("filtering: ")
EIC_pairs <- paired[paired[,4] >= 0,,drop=FALSE]
EIC_pairs[,1] <- peaklist[EIC_pairs[,1],"peak_ID"]
EIC_pairs[,2] <- peaklist[EIC_pairs[,2],"peak_ID"]
if(length(EIC_pairs[,1])==0){
cat("nothing found; aborted.");next;
}
those<-(EIC_pairs[,1]>EIC_pairs[,2])
if(any(those)){
EIC_pairs[those,]<-EIC_pairs[those,c(2,1,3,4), drop = FALSE]
}
EIC_pairs<-EIC_pairs[order(EIC_pairs[,1],EIC_pairs[,2], decreasing = FALSE),, drop = FALSE]
save(EIC_pairs,file=file.path(logfile[[1]],"results","componentization","EIC_corr", for_file))
cat(as.character(round((sum(EIC_pairs[,4] >= as.numeric(logfile$parameters$EICor_mincor)) / length(paired[,1]) * 100), digits = 2)));
cat("% = ");
cat(as.character(length(EIC_pairs[,1])));
cat(" pairs corrrelated - ");
rm(paired,EIC_pairs,peaklist,MSlist);
measurements[b, names(measurements)=="EIC_correlation"] <- "TRUE"
write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
cat("done.")
}else{
cat("\n EIC correlation done before.")
}
}
rm(measurements)
