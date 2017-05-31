
	measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
	measurements_incl<-measurements[measurements[,"include"]=="TRUE",,drop=FALSE]
	leng<-length(measurements_incl[,"include"])
	if(leng>0){
		for(i in 1:leng){
			if(measurements_incl[i,"LOD"]=="TRUE"){next}
			measurements[measurements[,"ID"]==measurements_incl[i,"ID"],"LOD"]<-"TRUE"
		}
	}
	write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
	rm(measurements,measurements_incl)


	# delete old LOD gams & results ####################################
	those<-list.files(file.path(logfile$project_folder,"results","LOD"))
	if(length(those)>0){
		for(i in 1:length(those)){
			file.remove(file.path(logfile$project_folder,"results","LOD",those[i]))
		}
	}


