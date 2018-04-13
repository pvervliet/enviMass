do_isot <- (logfile$workflow[names(logfile$workflow) == "isotopologues"] == "yes")
do_addu <- (logfile$workflow[names(logfile$workflow) == "adducts"] == "yes")
do_homol <- (logfile$workflow[names(logfile$workflow) == "homologues"] == "yes")
if( do_isot | do_addu ){
measurements <- read.csv(file = file.path(logfile[[1]], "dataframes", "measurements"), colClasses = "character");
if(mute(logfile$parameters$prof_select == "TRUE")){
for_IDs <- measurements[(measurements$include == "TRUE") & (measurements$components_files == "FALSE") & (measurements$profiled != "FALSE"),]$ID
}else{
for_IDs <- measurements[(measurements$include == "TRUE") & (measurements$components_files == "FALSE") ,]$ID
}
if(length(for_IDs)){
for(i in for_IDs){
if(file.exists(file.path(logfile[[1]], "results", "componentization", "components", paste(i)))){
file.remove(file.path(logfile[[1]], "results", "componentization", "components", paste(i)))
}
}
if(FALSE){
for(i in for_IDs) enviMass:::QTGGY3365L(x = i, logfile, measurements)
}
clusterEvalQ(cl = clus,{rm(list=ls()); NULL})
if(
(logfile$parameters$do_atom_bounds_components == "TRUE") &
(length(logfile$parameters$atom_bounds_components))
){
if(any(objects(envir=as.environment(".GlobalEnv"))=="LOD_splined")){rm(LOD_splined, envir = as.environment(".GlobalEnv"))}
if(any(objects() == "LOD_splined")){rm(LOD_splined)}
if(file.exists(file = file.path(logfile$project_folder,"results","LOD","LOD_splined"))){
load(file = file.path(logfile$project_folder,"results","LOD","LOD_splined"))
do_LOD <- TRUE
clusterExport(cl = clus, varlist = c("do_LOD", "LOD_splined", "isotopes"), envir = environment())
}else{
do_LOD <- FALSE
clusterExport(cl = clus, varlist = c("do_LOD", "isotopes"), envir = environment())
}
}
clusterExport(cl = clus, varlist = c("do_isot", "do_addu", "do_homol"), envir = environment())
cluster_results <- clusterApplyLB(cl = clus,
x = for_IDs,
fun = enviMass:::QTGGY3365L,
logfile = logfile,
measurements = measurements
)
clusterEvalQ(cl = clus,{rm(list=ls()); NULL})
}
measurements[!is.na(match(measurements$ID, for_IDs)), "components_files"] <- "TRUE"
write.csv(measurements, file = file.path(logfile[[1]], "dataframes", "measurements"), row.names = FALSE);
rm(measurements)
}
