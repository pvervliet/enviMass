SOPHZ8293U <- function(
peaks,
index_prof,
...
){
if(blindsub != FALSE){subit = 1; subrat = blindfold;}else{subit = 2; subrat = 0;}
if(blindsub != FALSE){if(!is.numeric(blindfold) || (blindfold < 0)){stop("Invalid blindfold argument; aborted.")}}
m = 1
n = dim(index_prof)[1]
for(k in m:n){
if(index_prof[k,"number_peaks_total"] > 1){
timeset[,4:length(timeset[1,])] <- 0;
timeset[,c(4,5)] <- .Call("_enviMass_fill_timeset",
as.numeric(timeset),
as.numeric(peaks[(index_prof[k, "start_ID"]:index_prof[k, "end_ID"]), "sampleIDs"]),
as.numeric(peaks[(index_prof[k, "start_ID"]:index_prof[k, "end_ID"]), "intensity"]),
as.integer(length(timeset[,1])),
PACKAGE="enviMass"
)
if(any(timeset[,4] > 0)){
what<-1
if(!omit_trend){
that <- .Call("_enviMass_meandel",
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
if(!omit_trend){
index_prof[k,"deltaint_newest"] <- max(that[4,]);
index_prof[k,"deltaint_global"] <- max(that[5,]);
index_prof[k,"absolute_mean_dev"] <- max(that[3,]);
}
}
index_prof[k,"newest_intensity"] <- timeset[leng, 4]
}else{
index_prof[k, "absolute_mean_dev"] <- 0;
if( any( timeset[,3] == (peaks[index_prof[k,1], 6]) ) ){
}else{
index_prof[k, "deltaint_global"] <- (peaks[index_prof[k,1], 2])
if(any(peaks[index_prof[k,1], 6] == latestID)){
index_prof[k, "deltaint_newest"] <- (peaks[index_prof[k,1], 2])
index_prof[k, "newest_intensity"] <- (peaks[index_prof[k,1], 2])
}
}
}
}
return(index_prof)
}
