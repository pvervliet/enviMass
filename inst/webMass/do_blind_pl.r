measurements <- read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
for_IDs <- measurements$ID[
(measurements$include == "TRUE") & (measurements$blind == "FALSE") & (measurements$Type != "blank")
]
if(length(for_IDs)){
if(any(objects()=="peaklist")){rm(peaklist)}
clusterEvalQ(cl = clus,{rm(list = ls()); NULL})
cluster_results <- clusterApplyLB(cl = clus,
x = for_IDs,
fun = enviMass:::BTYMI6379F,
logfile = logfile,
measurements = measurements,
ppm = logfile$parameters$blind_ppm,
dmz = as.numeric(logfile$parameters$blind_dmz),
dRT = as.numeric(logfile$parameters$blind_drt),
int_ratio = as.numeric(logfile$parameters$blind_threshold)
)
clusterEvalQ(cl = clus,{rm(list = ls()); NULL})
cluster_results <- unlist (cluster_results)
if(any(cluster_results == "done")){
measurements$blind[match(for_IDs[cluster_results == "done"], measurements$ID)] <- "TRUE"
}
write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
rm(measurements, cluster_results)
}
