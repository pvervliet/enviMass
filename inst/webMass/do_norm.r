	
	use_int <- logfile$parameters$peak_which_intensity
	measurements <- read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
	measurements_incl <- measurements[measurements[,"include"]=="TRUE",,drop=FALSE]
	# positive #################################################################
	if(any(measurements_incl[,"Mode"]=="positive")){
		measurements_pos<-measurements_incl[measurements_incl[,"Mode"]=="positive",,drop=FALSE]
		leng<-length(measurements_pos[,8])
		meanint<-c();
		maxint<-c();
		minint<-c();
		for(i in 1:leng){
			if(any(objects(envir=as.environment(".GlobalEnv"))=="peaklist")){rm(peaklist,envir=as.environment(".GlobalEnv"))}
			if(any(objects()=="peaklist")){rm(peaklist)}
			load(file=file.path(logfile[[1]],"peaklist",as.character(measurements_pos[i,1])),envir=as.environment(".GlobalEnv"));
			if(use_int=="maximum"){
				meanint<-c(meanint,median(peaklist[,"max_int"]));
				maxint<-c(maxint,max(peaklist[,"max_int"]));
				minint<-c(minint,min(peaklist[,"max_int"]));
			}else{ # = area
				meanint<-c(meanint,median(peaklist[,"sum_int"]));
				maxint<-c(maxint,max(peaklist[,"sum_int"]));
				minint<-c(minint,min(peaklist[,"sum_int"]));
			}
			rm(peaklist,envir=as.environment(".GlobalEnv"))
		}
		atmean<-median(meanint)
		png(filename = file.path(logfile[[1]],"pics","int_distr_pos"), bg = "white",width = 680, height = 480)
		plot.new();
		plot.window(xlim=c(0,(leng+1)),ylim=c(log10(min(minint)),log10(max(maxint))))
		title(xlab="Measurement sequence (w/o outliers)",ylab="log10 Intensity",main="sample:green / blank:blue / other:grey")
		box();axis(1);axis(2);
		for(i in 1:leng){
			if(any(objects(envir=as.environment(".GlobalEnv"))=="peaklist")){rm(peaklist,envir=as.environment(".GlobalEnv"))}
			if(any(objects()=="peaklist")){rm(peaklist)}
			load(file=file.path(logfile[[1]],"peaklist",as.character(measurements_pos[i,"ID"])),envir=as.environment(".GlobalEnv"));
			doneit<-FALSE
			if(measurements_pos[i,"Type"]=="sample"){
				abline(h=log10(atmean),col="red",lty=2)
				if(use_int=="max_int"){
					boxplot(log10(peaklist[,"max_int"]),add=TRUE,at=i,cex=0.5,col="green")
				}else{
					boxplot(log10(peaklist[,"sum_int"]),add=TRUE,at=i,cex=0.5,col="green")				
				}
				doneit<-TRUE
			}
			if(measurements_pos[i,"Type"]=="blank"){
				abline(h=log10(atmean),col="red",lty=2)
				if(use_int=="max_int"){
					boxplot(log10(peaklist[,"max_int"]),add=TRUE,at=i,cex=0.5,col="blue")
				}else{
					boxplot(log10(peaklist[,"sum_int"]),add=TRUE,at=i,cex=0.5,col="blue")				
				}
				doneit<-TRUE
			}
			if(doneit==FALSE){
				abline(h=log10(atmean),col="red",lty=2)
				if(use_int=="max_int"){
					boxplot(log10(peaklist[,"max_int"]),add=TRUE,at=i,cex=0.5,col="lightgrey")
				}else{
					boxplot(log10(peaklist[,"sum_int"]),add=TRUE,at=i,cex=0.5,col="lightgrey")				
				}				
				doneit<-TRUE
			}
			if(use_int=="max_int"){
				peaklist[,"int_corr"]<-c(peaklist[,"max_int"]/median(peaklist[,"max_int"])*atmean)
			}else{
				peaklist[,"int_corr"]<-c(peaklist[,"sum_int"]/median(peaklist[,"sum_int"])*atmean)			
			}
			save(peaklist,file=file.path(logfile[[1]],"peaklist",as.character(measurements_pos[i,"ID"])))
			if(any(objects(envir=as.environment(".GlobalEnv"))=="peaklist")){rm(peaklist,envir=as.environment(".GlobalEnv"))}
			if(any(objects()=="peaklist")){rm(peaklist)}
			# although not done strictly file-wise, mark normalization:
			measurements[measurements[,"ID"]==measurements_pos[i,"ID"],"norm"]<-"TRUE"
		}
		dev.off()
		expr1p<-list(src=file.path(logfile[[1]],"pics","int_distr_pos"))
		output$pic_int_distr_pos<-renderImage(expr1p, deleteFile = FALSE)
		rm(measurements_pos)
		write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
	}	
	# negative #################################################################	
	if(any(measurements_incl[,"Mode"]=="negative")){
		measurements_neg<-measurements_incl[measurements_incl[,"Mode"]=="negative",,drop=FALSE]
		leng<-length(measurements_neg[,8])
		meanint<-c();
		maxint<-c();
		minint<-c();
		for(i in 1:leng){
			if(any(objects(envir=as.environment(".GlobalEnv"))=="peaklist")){rm(peaklist,envir=as.environment(".GlobalEnv"))}
			if(any(objects()=="peaklist")){rm(peaklist)}
			load(file=file.path(logfile[[1]],"peaklist",as.character(measurements_neg[i,"ID"])),envir=as.environment(".GlobalEnv"));
			if(use_int=="maximum"){
				meanint<-c(meanint,median(peaklist[,"max_int"]));
				maxint<-c(maxint,max(peaklist[,"max_int"]));
				minint<-c(minint,min(peaklist[,"max_int"]));
			}else{	
				meanint<-c(meanint,median(peaklist[,"sum_int"]));
				maxint<-c(maxint,max(peaklist[,"sum_int"]));
				minint<-c(minint,min(peaklist[,"sum_int"]));			
			}	
			rm(peaklist,envir=as.environment(".GlobalEnv"))
		}
		atmean<-median(meanint)
		png(filename = file.path(logfile[[1]],"pics","int_distr_neg"), bg = "white",width = 680, height = 480)
		plot.new();
		plot.window(xlim=c(0,(leng+1)),ylim=c(log10(min(minint)),log10(max(maxint))))
		title(xlab="Measurement sequence (w/o outliers)",ylab="log10 Intensity",main="sample:green / blank:blue / other:grey")
		box();axis(1);axis(2);
		for(i in 1:leng){
			if(any(objects(envir=as.environment(".GlobalEnv"))=="peaklist")){rm(peaklist,envir=as.environment(".GlobalEnv"))}
			if(any(objects()=="peaklist")){rm(peaklist)}
			load(file=file.path(logfile[[1]],"peaklist",as.character(measurements_neg[i,"ID"])),envir=as.environment(".GlobalEnv"));
			doneit<-FALSE
			if(measurements_neg[i,"Type"]=="sample"){
				abline(h=log10(atmean),col="red",lty=2)
				if(use_int=="max_int"){
					boxplot(log10(peaklist[,"max_int"]),add=TRUE,at=i,cex=0.5,col="green")
				}else{
					boxplot(log10(peaklist[,"sum_int"]),add=TRUE,at=i,cex=0.5,col="green")
				}
				doneit<-TRUE
			}
			if(measurements_neg[i,"Type"]=="blank"){
				abline(h=log10(atmean),col="red",lty=2)
				if(use_int=="max_int"){
					boxplot(log10(peaklist[,"max_int"]),add=TRUE,at=i,cex=0.5,col="blue")
				}else{
					boxplot(log10(peaklist[,"sum_int"]),add=TRUE,at=i,cex=0.5,col="blue")
				}
				doneit<-TRUE
			}
			if(doneit==FALSE){
				abline(h=log10(atmean),col="red",lty=2)
				if(use_int=="max_int"){
					boxplot(log10(peaklist[,"max_int"]),add=TRUE,at=i,cex=0.5,col="lightgrey")
				}else{
					boxplot(log10(peaklist[,"sum_int"]),add=TRUE,at=i,cex=0.5,col="lightgrey")
				}
				doneit<-TRUE
			}
			if(use_int=="max_int"){
				peaklist[,"int_corr"]<-c(peaklist[,"max_int"]/median(peaklist[,"max_int"])*atmean)			
			}else{
				peaklist[,"int_corr"]<-c(peaklist[,"sum_int"]/median(peaklist[,"sum_int"])*atmean)
			}
			save(peaklist,file=file.path(logfile[[1]],"peaklist",as.character(measurements_neg[i,"ID"])))
			if(any(objects(envir=as.environment(".GlobalEnv"))=="peaklist")){rm(peaklist,envir=as.environment(".GlobalEnv"))}
			if(any(objects()=="peaklist")){rm(peaklist)}
			# although not done strictly file-wise, mark normalization:
			measurements[measurements[,"ID"]==measurements_neg[i,"ID"],"norm"]<-"TRUE"
		}
		dev.off()
		expr1n<-list(src=file.path(logfile[[1]],"pics","int_distr_neg"))
		output$pic_int_distr_neg<-renderImage(expr1n, deleteFile = FALSE)
		rm(measurements_neg)
		write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
	}	
	rm(measurements,measurements_incl)
	############################################################################
