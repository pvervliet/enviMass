# Run adduct grouping, filewise

    measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
	cat("Adduct grouping: ")
	for(b in 1:length(measurements[,"ID"])){
		if( 
			(measurements[b,"include"]=="TRUE") & 			# included?
			(measurements[b,"adducts"]=="FALSE")  	# not yet done
		){ 

			##########################################################################
			# exclude files that do not end up in profiles ###########################
			if( (mute(logfile$parameters$prof_select=="TRUE")) & (measurements$profiled[b]=="FALSE") ){
				cat("\n Skip file - not included in profile building.");next;
			}
			##########################################################################
			cat(paste("\n Doing ",as.character(b)," of ",as.character(length(measurements[,"ID"]))," files: ",sep=""))		
			##########################################################################
			# LOAD FILES & REMOVE OLD RESULTS ########################################
			cat("loading - ")
			for_file<-measurements[b,"ID"]
			if( file.exists(file.path(logfile[[1]],"results","componentization","adducts",for_file) ) ){
				file.remove(file.path(logfile[[1]],"results","componentization","adducts",for_file) )
			}			
			if( file.exists(file.path(logfile[[1]],"results","componentization","adducts",paste("full",for_file,sep="_")) ) ){
				file.remove(file.path(logfile[[1]],"results","componentization","adducts",paste("full",for_file,sep="_")) )
			}			
			# Peaklist
			load(file=file.path(logfile[[1]],"peaklist",as.character(for_file)));   
			peaklist<-peaklist[order(peaklist[,"peak_ID"],decreasing=FALSE),] # match with IDs - for saving pattern; IDs are retrieved for pairs seperately
			# EIC pairs
			if((logfile$workflow[names(logfile$workflow)=="EIC_correlation"]=="yes") & TRUE){ # load EIC correlation results - removed, also in EIC -> isot. depends matrix!
				if(file.exists(file.path(logfile[[1]],"results","componentization","EIC_corr",for_file))){
					load(file.path(logfile[[1]],"results","componentization","EIC_corr",for_file))
					exclude<-EIC_pairs[
						EIC_pairs[,4]<as.numeric(logfile$parameters$EICor_mincor) 
					,1:2,drop=FALSE]
					rm(EIC_pairs)
					if(length(exclude[,1])==0){
						exclude<-FALSE
					}else{
						cat("with exclusion - ")
					}
				}else{
					exclude<-FALSE
				}
			}else{
				exclude<-FALSE
			}
			##########################################################################	
			cat("grouping - ")			
			if(measurements[b,names(measurements)=="Mode"]=="positive"){
				with_adducts<-logfile$adducts_pos_group
				with_mode<-"positive"
			}else{
				with_adducts<-logfile$adducts_neg_group
				with_mode<-"negative"
			}
			if(length(with_adducts)<2){
				cat("\n Not enough adducts for this ionization mode specified - skipped ...")
				next;
			}			
			peaklist2 <- as.data.frame(peaklist[peaklist[,"keep"]==1, c("m/z_corr","int_corr","RT_corr","peak_ID")])	
			if(logfile$parameters$adducts_ppm=="TRUE"){
				use_mztol<-as.numeric(logfile$parameters$adducts_mztol)
			}else{ # mmu
				use_mztol<-(as.numeric(logfile$parameters$adducts_mztol)/1000)
			}
			adduct<-try(
				enviMass::adduct_search2( # dont name it "adducts" -> conflict
					peaklist=peaklist2[,c("m/z_corr","int_corr","RT_corr","peak_ID")], 
					adducts, 
					rttol = as.numeric(logfile$parameters$adducts_rttol), 
					mztol = use_mztol,
					ppm = logfile$parameters$adducts_ppm, 
					use_adducts = with_adducts, 
					ion_mode = with_mode,
					exclude
				)
			)
			if(class(adduct)=="try-error"){
				cat("\n Adduct detection failed - adpat parameters?");
				next;
			}				
			if(length(adduct[["Pairs"]][,1])==0){
				cat("\n No adduct relations detected");
				next;
			}
			Adduct_pairs <- adduct[["Pairs"]][,c(1,2)]
			those<-(Adduct_pairs[,1]>Adduct_pairs[,2])
			if(any(those)){
				Adduct_pairs[those,]<-Adduct_pairs[those,c(2,1)]
			}
			Adduct_pairs<-Adduct_pairs[order(Adduct_pairs[,1],Adduct_pairs[,2],decreasing=FALSE),]
			save(Adduct_pairs,file=(file.path(logfile[[1]],"results","componentization","adducts",paste(for_file,sep=""))))
			adduct[["Pairs"]]<-0
			save(adduct,file=(file.path(logfile[[1]],"results","componentization","adducts",paste("full",for_file,sep="_"))))
			rm(peaklist,peaklist2,those,Adduct_pairs,with_adducts,with_mode,for_file,adduct)
			##########################################################################	
			measurements[b,"adducts"]<-"TRUE"
			write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
			cat("done.")
			##########################################################################		
		}else{
			cat("\n Adducts extracted before or file not included.")
		}
	}
	rm(measurements)
	