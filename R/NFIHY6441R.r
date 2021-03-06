NFIHY6441R <- function(
clus_centroids,
...
){
clus_centroids[, "measureID"] <- seq(1, length(clus_centroids[, "measureID"]), 1);
if(dim(clus_centroids)[1] >= minpeak){
out1 <- .Call("gapfill",
as.numeric(clus_centroids[, "RT"]),
as.numeric(clus_centroids[, "intensity"]),
as.integer(order(clus_centroids[, "RT"], decreasing = FALSE)),
as.numeric(clus_centroids[, "m/z"]),
as.numeric(clus_centroids[, "measureID"]),
as.numeric(Retens),
as.numeric(drtfill),
PACKAGE = "enviPick"
)
out1 <- matrix(out1, ncol = 10);
colnames(out1) <- c("m/z", "intens", "RT", "index", "intens_filt", "1pick", "pickcrit", "baseline", "intens_corr", "2pick");
out1[,5] <- out1[,2];
out2 <- .Call("pickpeak",
as.numeric(out1),
as.numeric(drtsmall),
as.numeric(drttotal),
as.integer(minpeak),
as.integer(recurs),
as.numeric(weight),
as.numeric(SB),
as.numeric(SN),
as.numeric(minint),
as.numeric(maxint),
as.integer(ended),
as.integer(2),
PACKAGE = "enviPick"
)
out2 <- matrix(out2, ncol = 10);
out2 <- out2[out2[,4] != 0,]
out3 <- rep(0, length(clus_centroids[, "RT"]))
out3[out2[,4]] <- out2[,10]
return(out3)
}else{
return(rep(0, dim(clus_centroids)[1]))
}
}
