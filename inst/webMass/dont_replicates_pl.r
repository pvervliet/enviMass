
##########################################################################
measurements <- read.csv(file = file.path(logfile[[1]], "dataframes", "measurements"), colClasses = "character");
for_IDs <- measurements$ID[
	(measurements$include == "TRUE") & (measurements$blind == "FALSE") & (measurements$tag3 != "FALSE")
]
###############################################################################



##########################################################################
if(length(for_IDs)){

	######################################################################
	reset_replic <- function(x, logfile){
		for_file <- x
		if(any(objects(envir = as.environment(".GlobalEnv")) == "peaklist")){rm(peaklist, envir = as.environment(".GlobalEnv"))}
		if(any(objects() == "peaklist")){rm(peaklist)}
		load(file = file.path(logfile[[1]], "peaklist", for_file), envir = as.environment(".GlobalEnv"));
		keep <- rep(Inf, length(peaklist[,1])) 
		peaklist[,colnames(peaklist)=="keep"] <- keep_2
		save(peaklist,file = file.path(logfile[[1]], "peaklist", for_file))
		rm(peaklist)
		return("done")
	}
	######################################################################
	clusterEvalQ(cl = clus,{rm(list = ls()); NULL})
	cluster_results <- clusterApplyLB(cl = clus, 
			x = for_IDs, 
			fun = reset_replic, 
			logfile = logfile	
	)
	clusterEvalQ(cl = clus,{rm(list = ls()); NULL})	  
	rm(cluster_results)
	######################################################################
	measurements$blind[match(for_IDs, measurements$ID)] <- "TRUE"
	write.csv(measurements, file = file.path(logfile[[1]], "dataframes", "measurements"), row.names = FALSE);

}
##########################################################################
rm(measurements)


