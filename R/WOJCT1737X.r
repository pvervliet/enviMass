WOJCT1737X <- function(
x,
logfile,
do_cor
){
for_file <- x
if( file.exists(file.path(logfile[[1]], "results", "componentization", "EIC_corr", for_file) ) ){
file.remove(file.path(logfile[[1]], "results", "componentization", "EIC_corr", for_file) )
}
load(file = file.path(logfile[[1]], "peaklist", as.character(for_file)));
load(file = file.path(logfile[[1]], "MSlist", as.character(for_file)));
ord <- order(peaklist[,"RT_corr"], decreasing = FALSE)
peaklist <- peaklist[ord,,drop = FALSE]
use_peak <- MSlist[[7]][peaklist[,10],3] >= as.numeric(logfile$parameters$EICor_minpeaks)
at <- 1
paired <- matrix(ncol = 4, nrow = 1E7, 0)
len <- dim(paired)[1]
for(i in 1:(length(peaklist[,1]) - 1)){
if(!use_peak[i]){next}
for(j in (i+1):length(peaklist[,1])){
if(!use_peak[j]){next}
RTdif <- (peaklist[j, 14] - peaklist[i, 14])
if(RTdif <= as.numeric(logfile$parameters$EICor_delRT)){
paired[at, 1] <- i
paired[at, 2] <- j
paired[at, 3] <- RTdif
at <- (at + 1)
if(at > len){
paired2 <- matrix(ncol = 4, nrow = 1E7, 0)
paired <- rbind(paired, paired2)
len <- dim(paired)[1]
}
}else{
break;
}
}
}
paired <- paired[1:(at-1),, drop = FALSE]
if(all(paired[,1]==0)){
return("nothing found; aborted.");
}
skipped <- 0
for(i in 1:length(paired[,1])){
if(i > 1){
if(paired[i,1] != paired[i-1, 1]){
PEAK_ID_1 <- peaklist[paired[i, 1], "peak_ID"]
start_1 <- MSlist[[7]][PEAK_ID_1, 1]
end_1 <- MSlist[[7]][PEAK_ID_1, 2]
del_1 <- MSlist[[7]][PEAK_ID_1, 3]
}
}else{
PEAK_ID_1 <- peaklist[paired[i,1], 10]
start_1 <- MSlist[[7]][PEAK_ID_1, 1]
end_1 <- MSlist[[7]][PEAK_ID_1, 2]
del_1 <- MSlist[[7]][PEAK_ID_1, 3]
}
PEAK_ID_2 <- peaklist[paired[i,2], 10]
start_2 <- MSlist[[7]][PEAK_ID_2, 1]
end_2 <- MSlist[[7]][PEAK_ID_2, 2]
del_2 <- MSlist[[7]][PEAK_ID_2, 3]
if(del_1 <= del_2){
those <- match(MSlist[[4]][[2]][start_1:end_1, 3],MSlist[[4]][[2]][start_2:end_2, 3])
if(sum(!is.na(those)) < as.numeric(logfile$parameters$EICor_minpeaks)){
paired[i,4] <- (-2000)
skipped <- (skipped + 1);
next;
}
if(do_cor){
paired[i,4] <- cor(
MSlist[[4]][[2]][start_1:end_1,2][!is.na(those)],
MSlist[[4]][[2]][start_2:end_2,2][those[!is.na(those)]],
method = "spearman"
)
}
}else{
those <- match(MSlist[[4]][[2]][start_2:end_2,3],MSlist[[4]][[2]][start_1:end_1,3])
if(sum(!is.na(those)) < as.numeric(logfile$parameters$EICor_minpeaks)){
paired[i,4] <- (-2000)
skipped <- (skipped + 1);
next;
}
if(do_cor){
paired[i,4] <- cor(
MSlist[[4]][[2]][start_2:end_2,2][!is.na(those)],
MSlist[[4]][[2]][start_1:end_1,2][those[!is.na(those)]],
method = "spearman"
)
}
}
}
EIC_pairs <- paired[paired[,4] >= 0,, drop = FALSE]
EIC_pairs[,1] <- peaklist[EIC_pairs[,1],"peak_ID"]
EIC_pairs[,2] <- peaklist[EIC_pairs[,2],"peak_ID"]
if(length(EIC_pairs[,1]) == 0){
return("nothing found; aborted.");
}
those <- (EIC_pairs[,1] > EIC_pairs[,2])
if(any(those)){
EIC_pairs[those,] <- EIC_pairs[those, c(2,1,3,4), drop = FALSE]
}
EIC_pairs <- EIC_pairs[order(EIC_pairs[,1],EIC_pairs[,2], decreasing = FALSE),, drop = FALSE]
save(EIC_pairs, file = file.path(logfile[[1]],"results","componentization","EIC_corr", for_file))
rm(paired, EIC_pairs, peaklist, MSlist);
return("done")
}
