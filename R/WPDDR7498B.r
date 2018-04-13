WPDDR7498B <- function(
x,
logfile,
frac = FALSE,
blind_omit = FALSE
){
if(file.exists(file.path(logfile[[1]],"peaklist",x))){
load(file=file.path(logfile[[1]],"peaklist",x),verbose=FALSE)
}else{
return(NULL)
}
if(blind_omit == TRUE){
peaklist <- peaklist[(peaklist[,colnames(peaklist) == "keep_2"] >= as.numeric(logfile$parameters$blind_threshold)),, drop = FALSE]
}
peaklist <- peaklist[(peaklist[,colnames(peaklist) == "keep"] == 1),, drop = FALSE]
if(length(peaklist[,1]) == 0){stop(paste0("\n startprofiles_pl: emtpy peaklist detected for file_ID: ", x))}
if(frac != FALSE){
peaklist <- peaklist[order(peaklist[,2], decreasing = TRUE),];
that <- floor(length(peaklist[,1]) * frac)
}else{
that <- length(peaklist[,1])
}
use_columns <- c("m/z_corr", "int_corr", "RT_corr", "peak_ID")
peaks <- as.matrix(cbind( peaklist[1:that, use_columns],
rep(0, that), rep(as.numeric(x), that), rep(0, that), rep(0, that), peaklist[1:that, colnames(peaklist) == "keep_2"])
);
colnames(peaks) <- c("m/z", "intensity", "RT", "peakIDs", "links", "sampleIDs", "partitionIDs", "profileIDs", "in_blind")
return(peaks)
}
