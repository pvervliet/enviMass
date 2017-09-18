
##########################################################################
measurements <- read.csv(file = file.path(logfile[[1]], "dataframes", "measurements"), colClasses = "character");
for_IDs <- measurements[ (measurements$include == "TRUE") & (measurements$recal == "FALSE") ,]$ID
##########################################################################
	
##########################################################################
if(length(for_IDs)){

	######################################################################
	reset_mass <- function( x, logfile){
		for_file <- x
		if(any(objects(envir = as.environment(".GlobalEnv")) == "peaklist")){rm(peaklist, envir = as.environment(".GlobalEnv"))}
		if(any(objects() == "peaklist")){rm(peaklist)}
		load(file = file.path(logfile[[1]], "peaklist", for_file), envir = as.environment(".GlobalEnv"));
		peaklist[,"m/z_corr"] <- peaklist[,"m/z"];
		save(peaklist,file = file.path(logfile[[1]], "peaklist", for_file))
		rm(peaklist)
		return("done")
	}
	######################################################################
	clusterEvalQ(cl = clus,{rm(list = ls()); NULL})
	cluster_results <- clusterApplyLB(cl = clus, 
			x = for_IDs, 
			fun = reset_mass, 
			logfile = logfile	
	)
	clusterEvalQ(cl = clus,{rm(list = ls()); NULL})	  
	rm(cluster_results)
	######################################################################
	measurements$recal[match(for_IDs, measurements$ID)] <- "TRUE"
	write.csv(measurements, file = file.path(logfile[[1]], "dataframes", "measurements"), row.names = FALSE);

}
##########################################################################
recal_gams<-list.files(file.path(logfile[[1]],"results","recalibration"))
if(length(recal_gams)>0){
	for(i in 1:length(recal_gams)){
		file.remove(file.path(logfile[[1]],"results","recalibration",recal_gams[i])) 
	}
}
path=file.path(logfile[[1]],"pics","recal_none")
png(filename = path, bg = "white")
plot.new();plot.window(xlim=c(0,1),ylim=c(0,1));text(0.5,0.5,"nothing selected \n or not available",cex=1)
dev.off()
exprrec<-list(src=path)
output$recal_pic<-renderImage(exprrec, deleteFile = FALSE);		
output$peakhist_pic<-renderImage(exprrec, deleteFile = FALSE);	
##########################################################################
rm(measurements)




