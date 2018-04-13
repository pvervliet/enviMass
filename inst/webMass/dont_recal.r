measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
measurements_incl<-measurements[measurements[,"include"]=="TRUE",,drop=FALSE]
leng<-length(measurements_incl[,"include"])
if(leng>0){
for(i in 1:leng){
if(measurements_incl[i,"recal"]=="TRUE"){next}
if(any(objects(envir=as.environment(".GlobalEnv"))=="peaklist")){rm(peaklist,envir=as.environment(".GlobalEnv"))}
if(any(objects()=="peaklist")){rm(peaklist)}
load(file=file.path(logfile[[1]],"peaklist",as.character(measurements_incl[i,"ID"])),envir=as.environment(".GlobalEnv"));
peaklist[,"m/z_corr"] <- peaklist[,"m/z"];
save(peaklist,file=file.path(logfile[[1]],"peaklist",as.character(measurements_incl[i,"ID"])))
if(any(objects(envir=as.environment(".GlobalEnv"))=="peaklist")){rm(peaklist,envir=as.environment(".GlobalEnv"))}
if(any(objects()=="peaklist")){rm(peaklist)}
measurements[measurements[,"ID"]==measurements_incl[i,"ID"],"recal"]<-"TRUE"
}
}
write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
rm(measurements,measurements_incl)
recal_gams<-list.files(file.path(logfile[[1]],"results","recalibration"))
if(length(recal_gams)>0){
for(i in 1:length(recal_gams)){
file.remove(file.path(logfile[[1]],"results","recalibration",recal_gams[i]))
}
}
path=file.path(logfile[[1]],"pics","recal_none")
png(filename = path, bg = "white")
plot.new();plot.window(xlim=c(0,1),ylim=c(0,1));text(0.5,0.5,"nothing selected \n or not available",cex=1)
dev.off()
exprrec<-list(src=path)
output$recal_pic<-renderImage(exprrec, deleteFile = FALSE);
output$peakhist_pic<-renderImage(exprrec, deleteFile = FALSE);
