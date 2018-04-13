logfile$parameters$parallel
those<-list(0)
those[[1]]<-file.path(logfile[[1]],"results","pattern_pos_IS")
those[[2]]<-file.path(logfile[[1]],"results","patternRT_pos_IS")
those[[3]]<-file.path(logfile[[1]],"results","patternDelRT_pos_IS")
those[[4]]<-file.path(logfile[[1]],"results","pattern_neg_IS")
those[[5]]<-file.path(logfile[[1]],"results","patternRT_neg_IS")
those[[6]]<-file.path(logfile[[1]],"results","patternDelRT_neg_IS")
those[[7]]<-file.path(logfile[[1]],"results","pattern_pos_target")
those[[8]]<-file.path(logfile[[1]],"results","patternRT_pos_target")
those[[9]]<-file.path(logfile[[1]],"results","patternDelRT_pos_target")
those[[10]]<-file.path(logfile[[1]],"results","pattern_neg_target")
those[[11]]<-file.path(logfile[[1]],"results","patternRT_neg_target")
those[[12]]<-file.path(logfile[[1]],"results","patternDelRT_neg_target")
for(n in 1:length(those)){
if(file.exists(those[[n]])){
file.remove(those[[n]])
}
}
rm(those)
intstand <- read.table(file = file.path(logfile[[1]], "dataframes", "IS.txt"), header = TRUE, sep = "\t", colClasses = "character");
if(length(intstand$Formula[intstand$Formula != "-"]) > 0){
index <- (1:dim(intstand)[1])
clusterEvalQ(cl = clus,{rm(list = ls()); NULL})
clusterExport(cl = clus, varlist = c("adducts", "isotopes", "resolution_list"), envir = environment())
cluster_results <- clusterApplyLB(cl = clus,
x = index,
fun = enviMass:::QBQLS8277J,
compound_table = intstand,
logfile = logfile
)
clusterEvalQ(cl = clus,{rm(list = ls()); NULL})
pattern_pos_IS <- list()
names_pos <- c()
patternRT_pos_IS <- c()
patternDelRT_pos_IS <- c()
mz_pos_IS <- c()
RT_pos_IS <- c()
at_pos <- 1
pattern_neg_IS <- list()
names_neg <- c()
patternRT_neg_IS <- c()
patternDelRT_neg_IS <- c()
mz_neg_IS <- c()
RT_neg_IS <- c()
at_neg <- 1
for(i in 1:length(index)){
len <- (length(cluster_results[[i]]) -1)
RT <- (as.numeric(intstand[i,"RT"]) * 60)
if(intstand[i,"RT_tolerance"] != "FALSE"){
delRT <- (as.numeric(intstand[i,"RT_tolerance"]) * 60)
}else{
delRT <- as.numeric(logfile$parameters$IS_drt1)
}
with_mode <- intstand[i,"ion_mode"]
use_for_recal <- as.logical(intstand[i,"use_for_recalibration"])
restrict_adduct <- as.logical(intstand[i,"restrict_adduct"])
take_all <- cluster_results[[i]][[len + 1]]
for(j in 1:len){
if(cluster_results[[i]][[j]][1] == "nothing"){next}
if(with_mode == "positive"){
pattern_pos_IS[[at_pos]] <- cluster_results[[i]][[j]]
names_pos <- c(names_pos, names(cluster_results[[i]])[j])
at_pos <- (at_pos + 1)
patternRT_pos_IS <- c(patternRT_pos_IS, RT)
patternDelRT_pos_IS <- c(patternDelRT_pos_IS, delRT)
if(use_for_recal){
if(restrict_adduct & (j != 1)){next}
mz_pos_IS <- c(mz_pos_IS, cluster_results[[i]][[j]][which.max(cluster_results[[i]][[j]][,2]),1][[1]] )
RT_pos_IS <- c(RT_pos_IS, RT)
}
}
if(with_mode == "negative"){
pattern_neg_IS[[at_neg]] <- cluster_results[[i]][[j]]
names_neg <- c(names_neg, names(cluster_results[[i]])[j])
at_neg <- (at_neg + 1)
patternRT_neg_IS <- c(patternRT_neg_IS, RT)
patternDelRT_neg_IS <- c(patternDelRT_neg_IS, delRT)
if(use_for_recal){
if(restrict_adduct & (j != 1)){next}
mz_neg_IS <- c(mz_neg_IS, cluster_results[[i]][[j]][which.max(cluster_results[[i]][[j]][,2]),1][[1]] )
RT_neg_IS <- c(RT_neg_IS, RT)
}
}
}
}
if(length(pattern_pos_IS)){
names(pattern_pos_IS) <- names_pos
save(pattern_pos_IS, file = file.path(logfile[[1]], "results", "pattern_pos_IS"))
save(patternRT_pos_IS, file = file.path(logfile[[1]], "results", "patternRT_pos_IS"))
save(patternDelRT_pos_IS, file = file.path(logfile[[1]], "results", "patternDelRT_pos_IS"))
}
if(length(mz_pos_IS)){
intmass_pos_IS <- data.frame(mz_pos_IS, RT_pos_IS)
names(intmass_pos_IS) <- c("m/z","RT")
save(intmass_pos_IS, file = file.path(logfile[[1]], "results", "intmass_pos_IS"))
}
if(length(pattern_neg_IS)){
names(pattern_neg_IS) <- names_neg
save(pattern_neg_IS, file = file.path(logfile[[1]], "results", "pattern_neg_IS"))
save(patternRT_neg_IS, file = file.path(logfile[[1]], "results", "patternRT_neg_IS"))
save(patternDelRT_neg_IS, file = file.path(logfile[[1]], "results", "patternDelRT_neg_IS"))
}
if(length(mz_neg_IS)){
intmass_neg_IS <- data.frame(mz_neg_IS, RT_neg_IS)
names(intmass_neg_IS) <- c("m/z","RT")
save(intmass_neg_IS, file = file.path(logfile[[1]],"results","intmass_neg_IS"))
}
rm(cluster_results,	pattern_pos_IS, names_pos,patternRT_pos_IS, patternDelRT_pos_IS,mz_pos_IS,
RT_pos_IS, at_pos, pattern_neg_IS, names_neg, patternRT_neg_IS, patternDelRT_neg_IS, mz_neg_IS, RT_neg_IS, at_neg)
}
targets <- read.table(file=file.path(logfile[[1]], "dataframes", "targets.txt"), header = TRUE, sep = "\t", colClasses = "character")
if(length(targets$Formula[targets$Formula != "-"]) > 0){
index <- (1:dim(targets)[1])
clusterEvalQ(cl = clus,{rm(list = ls()); NULL})
clusterExport(cl = clus, varlist = c("adducts", "isotopes", "resolution_list"), envir = environment())
cluster_results <- clusterApplyLB(cl = clus,
x = index,
fun = enviMass:::QBQLS8277J,
compound_table = targets,
logfile = logfile
)
clusterEvalQ(cl = clus,{rm(list = ls()); NULL})
pattern_pos_target <- list()
names_pos <- c()
patternRT_pos_target <- c()
patternDelRT_pos_target <- c()
mz_pos_target <- c()
RT_pos_target <- c()
at_pos <- 1
pattern_neg_target <- list()
names_neg <- c()
patternRT_neg_target <- c()
patternDelRT_neg_target <- c()
mz_neg_target <- c()
RT_neg_target <- c()
at_neg <- 1
for(i in 1:length(index)){
len <- (length(cluster_results[[i]]) - 1)
RT <- (as.numeric(targets[i,"RT"]) * 60)
if(targets[i,"RT_tolerance"] != "FALSE"){
delRT <- (as.numeric(targets[i,"RT_tolerance"]) * 60)
}else{
delRT <- as.numeric(logfile$parameters$tar_drt1)
}
with_mode <- targets[i,"ion_mode"]
use_for_recal <- as.logical(targets[i, "use_for_recalibration"])
restrict_adduct <- as.logical(targets[i, "restrict_adduct"])
take_all <- cluster_results[[i]][[len+1]]
for(j in 1:len){
if(cluster_results[[i]][[j]][1] == "nothing"){next}
if(with_mode == "positive"){
pattern_pos_target[[at_pos]] <- cluster_results[[i]][[j]]
names_pos <- c(names_pos, names(cluster_results[[i]])[j])
at_pos <- (at_pos + 1)
patternRT_pos_target <- c(patternRT_pos_target, RT)
patternDelRT_pos_target <- c(patternDelRT_pos_target, delRT)
if(use_for_recal){
if(restrict_adduct & (j != 1)){next}
mz_pos_target <- c(mz_pos_target, cluster_results[[i]][[j]][which.max(cluster_results[[i]][[j]][,2]),1][[1]] )
RT_pos_target <- c(RT_pos_target, RT)
}
}
if(with_mode == "negative"){
pattern_neg_target[[at_neg]] <- cluster_results[[i]][[j]]
names_neg <- c(names_neg, names(cluster_results[[i]])[j])
at_neg <- (at_neg + 1)
patternRT_neg_target <- c(patternRT_neg_target, RT)
patternDelRT_neg_target <- c(patternDelRT_neg_target, delRT)
if(use_for_recal){
if(restrict_adduct & (j != 1)){next}
mz_neg_target <- c(mz_neg_target, cluster_results[[i]][[j]][which.max(cluster_results[[i]][[j]][,2]),1][[1]] )
RT_neg_target <- c(RT_neg_target, RT)
}
}
}
}
if(length(pattern_pos_target)){
names(pattern_pos_target) <- names_pos
save(pattern_pos_target, file = file.path(logfile[[1]], "results", "pattern_pos_target"))
save(patternRT_pos_target, file = file.path(logfile[[1]], "results", "patternRT_pos_target"))
save(patternDelRT_pos_target, file = file.path(logfile[[1]], "results", "patternDelRT_pos_target"))
}
if(length(mz_pos_target)){
intmass_pos_target <- data.frame(mz_pos_target, RT_pos_target)
names(intmass_pos_target) <- c("m/z","RT")
save(intmass_pos_target, file = file.path(logfile[[1]], "results", "intmass_pos_target"))
}
if(length(pattern_neg_target)){
names(pattern_neg_target) <- names_neg
save(pattern_neg_target, file = file.path(logfile[[1]], "results", "pattern_neg_target"))
save(patternRT_neg_target, file = file.path(logfile[[1]], "results", "patternRT_neg_target"))
save(patternDelRT_neg_target, file = file.path(logfile[[1]], "results", "patternDelRT_neg_target"))
}
if(length(mz_neg_target)){
intmass_neg_target <- data.frame(mz_neg_target, RT_neg_target)
names(intmass_neg_target) <- c("m/z","RT")
save(intmass_neg_target, file = file.path(logfile[[1]], "results", "intmass_neg_target"))
}
rm(pattern_pos_target, names_pos, patternRT_pos_target, patternDelRT_pos_target, mz_pos_target, RT_pos_target,
at_pos, pattern_neg_target, names_neg, patternRT_neg_target, patternDelRT_neg_target, mz_neg_target, RT_neg_target, at_neg)
}
