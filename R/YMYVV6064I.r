YMYVV6064I <- function(
peaks,
agglom,
...
){
m = 1
n = dim(agglom)[1]
startat<-c(0);
for(k in m:n){
if(agglom[k,"number_peaks"]>1){
delmz<-(max(peaks[agglom[k,"start_ID"]:agglom[k,"end_ID"],"m/z"]) - min(peaks[agglom[k,"start_ID"]:agglom[k,"end_ID"],"m/z"]))
if(ppm){
delmz <- (delmz / mean(peaks[agglom[k,"start_ID"]:agglom[k,"end_ID"],"m/z"]) * 1E6)
}else{
delmz <- (delmz / 1000)
}
delRT <- (max(peaks[agglom[k,"start_ID"]:agglom[k,"end_ID"],"RT"]) - min(peaks[agglom[k,"start_ID"]:agglom[k,"end_ID"],"RT"]))
if(
(delmz > (dmass*2)) ||
(delRT > dret) ||
(any(duplicated(peaks[agglom[k,"start_ID"]:agglom[k,"end_ID"],"sampleIDs"])))
){
if(!do_replicates){
those <- (agglom[k,1] : agglom[k,2])
clusters <- extractProfiles(
peaks = peaks[those, c("m/z", "intensity", "RT", "sampleIDs")],
in_order = order(peaks[those, "intensity"], decreasing = TRUE),
dmass = dmass,
ppm = ppm,
dret = dret
)
clusters <- (clusters + startat);
peaks[those, 8] <- clusters;
peaks[those,] <- (peaks[those,][order(clusters, decreasing = FALSE),]);
startat <- max(clusters)
}else{
those <- (agglom[k,1] : agglom[k,2])
pregroup <- rep(0, length(those))
with_replic <- replic_ID[those]
if(any(with_replic > 0)){
for_replic <- unique(with_replic)
for_replic <- for_replic[for_replic != 0]
startit <- 0
for(i in 1:length(for_replic)){
these <- those[with_replic == for_replic[i]]
if(length(these) == 1) next
clusters <- extractProfiles(
peaks = peaks[these, c("m/z", "intensity", "RT", "sampleIDs")],
in_order = order(peaks[these, "intensity"], decreasing = TRUE),
dmass = dmass,
ppm = ppm,
dret = dret
)
clusters <- (clusters + startit)
pregroup[which(with_replic == for_replic[i])] <- clusters
startit <- max(clusters)
}
}
clusters <- extractProfiles_replicates(
peaks = peaks[those, c("m/z", "intensity", "RT", "sampleIDs")],
in_order = order(peaks[those, "intensity"], decreasing = TRUE),
dmass = dmass,
ppm = ppm,
dret = dret,
pregroup = pregroup
)
clusters <- (clusters + startat);
peaks[those, "profileIDs"] <- clusters ;
peaks[those,] <- (peaks[those,][order(clusters, decreasing = FALSE),]);
startat <- max(clusters)
}
}else{
startat<-(startat+1);
peaks[agglom[k,"start_ID"]:agglom[k,"end_ID"],"profileIDs"]<-startat
}
}else{
startat<-(startat+1);
peaks[agglom[k,"start_ID"]:agglom[k,"end_ID"],"profileIDs"]<-startat
}
}
return(peaks)
}
