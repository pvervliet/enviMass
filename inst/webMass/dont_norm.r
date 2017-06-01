
	cat("Inserting intensities \n");
	use_int<-logfile$parameters$peak_which_intensity
	measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
	measurements_incl<-measurements[measurements[,"include"]=="TRUE",,drop=FALSE]
	leng<-length(measurements_incl[,"include"])
	if(leng>0){
		for(i in 1:leng){
			if(measurements_incl[i,"norm"]=="TRUE"){next}
			if(any(objects(envir=as.environment(".GlobalEnv"))=="peaklist")){rm(peaklist,envir=as.environment(".GlobalEnv"))}
			if(any(objects()=="peaklist")){rm(peaklist)}
			load(file=file.path(logfile[[1]],"peaklist",as.character(measurements_incl[i,"ID"])),envir=as.environment(".GlobalEnv"));
			if(use_int=="maximum"){
				peaklist[,"int_corr"]<-peaklist[,"max_int"];
			}else{
				peaklist[,"int_corr"]<-peaklist[,"sum_int"];				
			}
			save(peaklist,file=file.path(logfile[[1]],"peaklist",as.character(measurements_incl[i,"ID"])))
			if(any(objects(envir=as.environment(".GlobalEnv"))=="peaklist")){rm(peaklist,envir=as.environment(".GlobalEnv"))}
			if(any(objects()=="peaklist")){rm(peaklist)}
			# although not done strictly file-wise, mark normalization:
			measurements[measurements[,"ID"]==measurements_incl[i,"ID"],"norm"]<-"TRUE"
		}
	}
	write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
	rm(measurements,measurements_incl)

	path=file.path(logfile[[1]],"pics","int_distr_pos")
		png(filename = path, bg = "white")
		plot.new();plot.window(xlim=c(0,1),ylim=c(0,1));text(0.5,0.5,"nothing selected \n or not available",cex=1)
		dev.off()
		expr3p<-list(src=path)
		output$pic_int_distr_pos<-renderImage(expr3p, deleteFile = FALSE)
		path=file.path(logfile[[1]],"pics","int_distr_neg")
	png(filename = path, bg = "white")
		plot.new();plot.window(xlim=c(0,1),ylim=c(0,1));text(0.5,0.5,"nothing selected \n or not available",cex=1)
		dev.off()
		expr3n<-list(src=path)
		output$pic_int_distr_neg<-renderImage(expr3n, deleteFile = FALSE)	


