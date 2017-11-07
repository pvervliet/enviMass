
measurements <- read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
#measurements$LOD <- "FALSE"
for_IDs <- measurements$ID[(measurements$include == "TRUE") & (measurements$LOD == "FALSE")]


if(length(for_IDs)){

	#################################################################################
	clusterEvalQ(cl = clus,{rm(list = ls()); NULL})
	cluster_results <- clusterApplyLB(cl = clus, 
		x = for_IDs, 
		fun = enviMass:::LOD_wrap, 
		logfile = logfile
	)
	clusterEvalQ(cl = clus,{rm(list = ls()); NULL})	  
	#################################################################################
	if(file.exists(file.path(logfile$project_folder, "results", "LOD", "LOD_splined"))){ # load existing model
		load(file.path(logfile$project_folder, "results", "LOD", "LOD_splined"))
		# clean list
		if(length(LOD_splined)){ # list entry indice == as.numeric(file ID)
			for(i in 1:length(LOD_splined)){
				if(is.na(match(as.character(i), measurements$ID))){
					if(length(LOD_splined) >= i){
						LOD_splined[[i]] <- NULL 
						names(LOD_splined)[i] <- ""
					}
				}
			}
		}
	}else{
		LOD_splined<-list()
	}
	if(length(cluster_results)){
		for(i in 1:length(cluster_results)){
			if(cluster_results[[i]][1] != "nothing"){
				assign(paste("LOD_", for_IDs[i], sep=""), cluster_results[[i]]);
				LOD_splined[[as.numeric(for_IDs[i])]] <- get(paste("LOD_", for_IDs[i], sep=""))
				names(LOD_splined)[as.numeric(for_IDs[i])] <- paste("LOD_",for_IDs[i], sep="")
				weg <- paste("LOD_", for_IDs[i], sep="")
				rm(weg)
			}else{
				if(length(LOD_splined)){
					if(length(LOD_splined) >= as.numeric(for_IDs[i])){
						LOD_splined[[as.numeric(for_IDs[i])]] <- NULL 
						names(LOD_splined)[as.numeric(for_IDs[i])] <- ""
					}
				}
			}
		}
	}
	save(LOD_splined, file = file.path(logfile$project_folder, "results", "LOD", "LOD_splined"))
	rm(LOD_splined)
	#################################################################################
	cluster_results <-unlist(cluster_results)
	if(any(cluster_results != "nothing")){
		measurements$LOD[match(for_IDs[cluster_results=="done"], measurements$ID)]<-"TRUE"
	}
	write.csv(measurements, file = file.path(logfile[[1]], "dataframes", "measurements"), row.names=FALSE);
	rm(measurements, cluster_results)
	#################################################################################

}


