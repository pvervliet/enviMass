measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
measurements_incl<-measurements[measurements[,"include"]=="TRUE",,drop=FALSE]
leng<-length(measurements_incl[,"include"])
if(leng>0){
for(i in 1:leng){
if(measurements_incl[i,"isotopologues"]=="TRUE"){next}
measurements[measurements[,"ID"]==measurements_incl[i,"ID"],"isotopologues"]<-"TRUE"
}
}
write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
rm(measurements,measurements_incl)
remo<-list.files(file.path(logfile[[1]],"results","componentization","isotopologues"))
if(length(remo)>0){
for(i in 1:length(remo)){
file.remove(file.path(logfile[[1]],"results","componentization","isotopologues",remo[i]))
}
}
