	


	measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
	measurements_incl<-measurements[measurements[,"include"]=="TRUE",,drop=FALSE]
	leng<-length(measurements_incl[,"include"])
	if(leng>0){
		for(i in 1:leng){
			if(measurements_incl[i,"blind"]=="TRUE"){next}
			if(measurements_incl[i,"Type"]=="blank"){next}
			if(any(objects(envir=as.environment(".GlobalEnv"))=="peaklist")){rm(peaklist,envir=as.environment(".GlobalEnv"))}
			if(any(objects()=="peaklist")){rm(peaklist)}
			load(file=file.path(logfile[[1]],"peaklist",as.character(measurements_incl[i,"ID"])),envir=as.environment(".GlobalEnv"));
			peaklist[,colnames(peaklist)=="keep_2"] <- Inf
			save(peaklist,file=file.path(logfile[[1]],"peaklist",as.character(measurements_incl[i,"ID"])))
			if(any(objects(envir=as.environment(".GlobalEnv"))=="peaklist")){rm(peaklist,envir=as.environment(".GlobalEnv"))}
			if(any(objects()=="peaklist")){rm(peaklist)}
			# although not done strictly file-wise, mark normalization:
			measurements[measurements[,"ID"]==measurements_incl[i,"ID"],"blind"]<-"TRUE"
		}
	}
	write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
	rm(measurements,measurements_incl)




