HUESR4183N <- function(
profileList,
dmass = 3,
ppm = TRUE,
dret = 60,
replicates = FALSE,
IDs = FALSE,
clus
){
if(!profileList[[1]][[2]]){stop("run agglom first on that profileList; aborted.")}
if(!is.numeric(dmass)){stop("dmass must be numeric; aborted.")}
if(!is.numeric(dret)){stop("dret must be numeric; aborted.")}
if(!is.logical(ppm)){stop("ppm must be logical; aborted.")}
startat <- c(0);
do_replicates<-FALSE
replic_ID <- "FALSE"
if(any(replicates!="FALSE")){
replic <- replicates[duplicated(replicates)]
replic <- unique(replic)
replic <- replic[replic != "FALSE"]
if(length(replic) > 0){
if(length(replicates) != length(IDs)){
stop("\n replicates vector longer than file ID vector")
}
do_replicates <- TRUE
}
if(do_replicates){
IDs_rep <- match(replicates, replic, nomatch = 0)
replic_ID <- IDs_rep[match(profileList[["peaks"]][,"sampleIDs"], IDs)]
}
}
m = 1
n = length(profileList[["index_agglom"]][,1])
along<-seq(1,dim(profileList[["index_agglom"]])[1],1)
size_MB <- (as.numeric(object.size(profileList[["peaks"]])) / 1048576)
for_split <- (size_MB / 2)
for_split <- round( dim(profileList[["index_agglom"]])[1] / for_split )
along <- split(along, ceiling(seq_along(along) / for_split) )
clus_peaks <- list()
clus_agglom <- list()
for(i in 1:length(along)){
clus_peaks[[i]] <- profileList[["peaks"]][
(profileList[["index_agglom"]][along[[i]][1],"start_ID"]) : (profileList[["index_agglom"]][along[[i]][length(along[[i]])],"end_ID"])
,,drop=FALSE]
clus_agglom[[i]] <- profileList[["index_agglom"]][along[[i]][1] : along[[i]][length(along[[i]])],,drop=FALSE]
start_at <- (clus_agglom[[i]][1,"start_ID"][[1]]-1)
clus_agglom[[i]][,1] <- (clus_agglom[[i]][,1]-start_at)
clus_agglom[[i]][,2] <- (clus_agglom[[i]][,2]-start_at)
}
clusterEvalQ(cl = clus,{rm(list = ls()); gc(verbose=FALSE); NULL})
clusterExport(cl = clus,
varlist = c("dmass", "ppm", "dret", "do_replicates", "replic_ID"),
envir = environment())
cluster_results <- clusterMap(
cl = clus,
fun = enviMass:::YMYVV6064I,
peaks = clus_peaks,
agglom = clus_agglom,
RECYCLE = TRUE,
SIMPLIFY = FALSE,
USE.NAMES = FALSE,
.scheduling = c("dynamic")
)
clusterEvalQ(cl = clus,{rm(list = ls()); gc(verbose = FALSE); NULL})
if(length(cluster_results) > 1){
for(i in 2:length(cluster_results)){
cluster_results[[i]][,"profileIDs"] <- (
cluster_results[[i]][,"profileIDs"] + max(cluster_results[[i - 1]][,"profileIDs"])
)
}
}
profileList[["peaks"]] <- do.call(rbind, cluster_results)
rm(cluster_results, clus_agglom, clus_peaks)
index <- .Call("_enviMass_indexed",
as.integer(profileList[["peaks"]][, "profileIDs"]),
as.integer(max(profileList[["peaks"]][, "profileIDs"])),
as.integer(26),
PACKAGE="enviMass"
)
index <- index[index[,1] != 0,]
index[,4] <- seq(length(index[,4]))
colnames(index) <- c(
"start_ID",
"end_ID",
"number_peaks_total",
"profile_ID",
"deltaint_newest",
"deltaint_global",
"absolute_mean_dev",
"in_blind?",
"above_blind?",
"number_peaks_sample",
"number_peaks_blind",
"mean_int_sample",
"mean_int_blind",
"mean_mz",
"mean_RT",
"mean_int",
"newest_intensity",
"links",
"component",
"max_int_sample",
"max_int_blind",
"max_int",
"var_mz",
"min_RT",
"max_RT",
"Mass defect"
)
profileList[[7]] <- index
m = 1
n = dim(profileList[["index_prof"]])[1]
mean_mz <- rep(0, n)
mean_RT <- rep(0, n)
mean_int <- rep(0, n)
mean_int_blind <- rep(0, n)
mean_int_sample <- rep(0, n)
max_int <- rep(0, n)
max_int_blind <- rep(0, n)
max_int_sample <- rep(0, n)
count_blind <- rep(0, n)
count_sam <- rep(0, n)
min_RT <- rep(0, n)
max_RT <- rep(0, n)
var_mz <- rep(0, n)
Mass_defect <- rep(0, n)
sample_IDs <- as.numeric(profileList[["sampleID"]])
sample_types <- profileList[["type"]]
for(k in m:n){
sub_mat <- profileList[["peaks"]][profileList[["index_prof"]][k,1]:profileList[["index_prof"]][k, 2], c(1:3, 6), drop = FALSE]
len <- dim(sub_mat)[1]
mean_mass <- (sum(sub_mat[,1])/len)
mean_mz[k] <- mean_mass
mean_RT[k] <- (sum(sub_mat[,3])/len)
mean_int[k] <- (sum(sub_mat[,2])/len)
max_int[k] <- max(sub_mat[,2])
min_RT[k] <- min(sub_mat[,3])
max_RT[k] <- max(sub_mat[,3])
var_mz <- (sum((sub_mat[,1] - mean_mass)^2) / (len - 1))
Mass_defect[k] <- (round(mean_mass) - mean_mass)
are_blind <- sub_mat[,4] %in% (sample_IDs[sample_types == "blank"])
count_blind[k] <- sum(are_blind)
count_sam[k] <- (len - count_blind[k])
if(count_blind[k] > 0){
mean_int_blind[k] <- (sum(sub_mat[are_blind, "intensity"]) / (len - 1))
max_int_blind[k] <- max(sub_mat[are_blind, "intensity"])
}
if(count_sam[k] > 0){
mean_int_sample[k] <- (sum(sub_mat[!are_blind, "intensity"]) / (len - 1))
max_int_sample[k] <- max(sub_mat[!are_blind, "intensity"])
}
}
profileList[["index_prof"]][,"mean_mz"] <- mean_mz
profileList[["index_prof"]][,"mean_RT"] <- mean_RT
profileList[["index_prof"]][,"mean_int"] <- mean_int
profileList[["index_prof"]][,"mean_int_blind"] <- mean_int_blind
profileList[["index_prof"]][,"mean_int_sample"] <- mean_int_sample
profileList[["index_prof"]][,"max_int"] <- max_int
profileList[["index_prof"]][,"max_int_sample"] <- max_int_sample
profileList[["index_prof"]][,"max_int_blind"] <- max_int_blind
profileList[["index_prof"]][,"in_blind?"] <- (count_blind > 0)
profileList[["index_prof"]][,"number_peaks_blind"] <- count_blind
profileList[["index_prof"]][,"number_peaks_sample"] <- count_sam
profileList[["index_prof"]][,"min_RT"] <- min_RT
profileList[["index_prof"]][,"max_RT"] <- max_RT
profileList[["index_prof"]][,"Mass defect"] <- Mass_defect
profileList[["index_prof"]][,"var_mz"] <- var_mz
profileList[["state"]][[3]] <- TRUE
return(profileList)
}
