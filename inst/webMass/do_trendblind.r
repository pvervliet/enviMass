if( file.exists(file.path(logfile[[1]],"results","profileList_pos")) ){
if(any(objects(envir=as.environment(".GlobalEnv"))=="profileList_pos")){rm(profileList_pos,envir=as.environment(".GlobalEnv"))}
if(any(objects()=="profileList_pos")){rm(profileList_pos)}
load(file.path(as.character(logfile[[1]]),"results","profileList_pos"),envir=as.environment(".GlobalEnv"))
lagit<-as.numeric(strsplit(logfile$parameters$trend_lags,",")[[1]])
if(logfile$parameters$trend_blind=="yes"){
blindsubtract<-TRUE
}else{
blindsubtract<-FALSE
}
profileList_pos<-TPMZP8211J(
profileList=profileList_pos,
from=FALSE,
to=FALSE,
progbar=logfile$parameters$progressBar,
blindsub=blindsubtract,
blindfold=as.numeric(logfile$parameters$blind_threshold),
lags=lagit,
threshold=as.numeric(logfile$parameters$trend_threshold),
notrend=as.logical(logfile$parameters$notrend)
)
profileList_pos<<-profileList_pos
save(profileList_pos,file=file.path(as.character(logfile[[1]]),"results","profileList_pos"),compress=FALSE);
png(filename = file.path(as.character(logfile[[1]]),"pics","boxprofile_pos"), width = 800, bg = "white")
enviMass::KXNZN5280X(profileList_pos)
dev.off()
expr4p<-list(src=file.path(logfile[[1]],"pics","boxprofile_pos"))
output$boxprofile<-renderImage(expr4p, deleteFile = FALSE)
if(isolate(input$Ion_mode)=="positive"){
profileList<<-profileList_pos;
}
}
if( file.exists(file.path(logfile[[1]],"results","profileList_neg")) ){
if(any(objects(envir=as.environment(".GlobalEnv"))=="profileList_neg")){rm(profileList_neg,envir=as.environment(".GlobalEnv"))}
if(any(objects()=="profileList_neg")){rm(profileList_neg)}
load(file.path(as.character(logfile[[1]]),"results","profileList_neg"),envir=as.environment(".GlobalEnv"))
lagit<-as.numeric(strsplit(logfile$parameters$trend_lags,",")[[1]])
if(logfile$parameters$trend_blind=="yes"){
blindsubtract<-TRUE
}else{
blindsubtract<-FALSE
}
profileList_neg<-TPMZP8211J(
profileList_neg,
from=FALSE,
to=FALSE,
progbar=logfile$parameters$progressBar,
blindsub=blindsubtract,
blindfold=as.numeric(logfile$parameters$blind_threshold),
lags=lagit,
threshold=as.numeric(logfile$parameters$trend_threshold),
notrend=as.logical(logfile$parameters$notrend)
)
profileList_neg<<-profileList_neg
save(profileList_neg,file=file.path(as.character(logfile[[1]]),"results","profileList_neg"),compress=FALSE);
png(filename = file.path(as.character(logfile[[1]]),"pics","boxprofile_neg"), width = 800, bg = "white")
enviMass::KXNZN5280X(profileList_neg)
dev.off()
expr4n<-list(src=file.path(logfile[[1]],"pics","boxprofile_neg"))
output$boxprofile<-renderImage(expr4n, deleteFile = FALSE)
if(isolate(input$Ion_mode)=="negative"){
profileList<<-profileList_neg;
}
}
