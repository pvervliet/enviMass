measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
cat("Isotopologue grouping: ")
if(!file.exists(file.path(logfile[[1]], "dataframes", "quantiz"))){stop("\nNo isotopologue space available.")}
load(file.path(logfile[[1]],"dataframes","quantiz"))
if((quantiz$R_set!=logfile$parameters$resolution) & (quantiz$R_set!="Sciex_all")){
cat("\n WARNING: seems the quantized data for isotopologue grouping does NOT MATCH your selected resolution! Please resolve this issue.")
}
for(b in 1:length(measurements[,"ID"])){
if(
(measurements[b,"include"]=="TRUE") &
(measurements[b,"isotopologues"]=="FALSE")
){
if( (mute(logfile$parameters$prof_select=="TRUE")) & (measurements$profiled[b]=="FALSE") ){
cat("\n Skip file - not included in profile building.");next;
}
cat(paste("\n Doing ",as.character(b)," of ",as.character(length(measurements[,"ID"]))," files: ",sep=""))
cat("loading - ")
for_file <- measurements[b,"ID"]
if( file.exists(file.path(logfile[[1]],"results","componentization","isotopologues",for_file) ) ){
file.remove(file.path(logfile[[1]],"results","componentization","isotopologues",for_file) )
}
if( file.exists(file.path(logfile[[1]],"results","componentization","isotopologues",paste("full",for_file,sep="_")) ) ){
file.remove(file.path(logfile[[1]],"results","componentization","isotopologues",paste("full",for_file,sep="_")) )
}
load(file=file.path(logfile[[1]],"peaklist",as.character(for_file)));
peaklist <- peaklist[order(peaklist[,10],decreasing=FALSE),]
if((logfile$workflow[names(logfile$workflow)=="EIC_correlation"]=="yes") & TRUE){
if(file.exists(file.path(logfile[[1]],"results","componentization","EIC_corr",for_file))){
load(file.path(logfile[[1]],"results","componentization","EIC_corr",for_file))
exclude<-EIC_pairs[
EIC_pairs[,4]<as.numeric(logfile$parameters$EICor_mincor)
,1:2,drop=FALSE]
rm(EIC_pairs)
if(length(exclude[,1])==0){
exclude<-FALSE
}else{
cat("with exclusion - ")
}
}else{
exclude<-FALSE
}
}else{
exclude<-FALSE
}
cat("grouping - ")
peaklist2 <- as.data.frame(peaklist[peaklist[,"keep"] == 1, c("m/z_corr","int_corr","RT_corr","peak_ID")])
if(logfile$parameters$isotop_ppm=="TRUE"){
use_mztol <- as.numeric(logfile$parameters$isotop_mztol)
}else{
use_mztol <- (as.numeric(logfile$parameters$isotop_mztol)/1000)
}
pattern <- try(
enviMass::LWBSL2520N(
peaklist=peaklist2[,c("m/z_corr","int_corr","RT_corr","peak_ID")],
quantiz,
mztol=use_mztol,
ppm=logfile$parameters$isotop_ppm,
inttol=(as.numeric(logfile$parameters$isotop_inttol)/100),
rttol=as.numeric(logfile$parameters$isotop_rttol),
use_isotopes=FALSE,
use_charges=logfile$parameters$isotop_use_charges,
use_marker=TRUE,
quick=TRUE,
isotopes,
exclude
)
)
if(class(pattern) == "try-error"){
cat("\n Isotopologue detection failed - adapt parameters?");
next;
}
if(length(pattern[["Pairs"]][,1]) == 0){
cat("\n No adduct relations detected");
next;
}
Isot_pairs <- pattern[["Pairs"]]
pattern[["Pairs"]] <- 0
those <- (Isot_pairs[,1] > Isot_pairs[,2])
if(any(those)){
Isot_pairs[those,] <- Isot_pairs[those, c(2,1), drop = FALSE]
}
Isot_pairs <- Isot_pairs[order(Isot_pairs[,1], Isot_pairs[,2], decreasing = FALSE),, drop = FALSE]
save(Isot_pairs, file = (file.path(logfile[[1]],"results","componentization","isotopologues",paste(for_file,sep="_"))))
save(pattern, file = (file.path(logfile[[1]],"results","componentization","isotopologues",paste("full",for_file,sep="_"))))
measurements[b,"isotopologues"] <- "TRUE"
write.csv(measurements, file = file.path(logfile[[1]],"dataframes","measurements"), row.names = FALSE);
rm(peaklist, peaklist2, pattern, Isot_pairs)
cat("done.")
}else{
cat("\n Isotopologues detected before or file not included.")
}
}
rm(quantiz,measurements)
