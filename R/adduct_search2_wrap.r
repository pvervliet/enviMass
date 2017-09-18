adduct_search2_wrap <-function(
	x,
	logfile,
	measurements
){

	##########################################################################
	# LOAD FILES & REMOVE OLD RESULTS ########################################
	cat("loading - ")
	for_file<-x
	if( file.exists(file.path(logfile[[1]],"results","componentization","adducts",for_file) ) ){
		file.remove(file.path(logfile[[1]],"results","componentization","adducts",for_file) )
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
	if(measurements[which(measurements$ID==for_file),names(measurements)=="Mode"]=="positive"){
		with_adducts<-logfile$adducts_pos_group
		with_mode<-"positive"
	}else{
		with_adducts<-logfile$adducts_neg_group
		with_mode<-"negative"
	}
	if(length(with_adducts)<2){
		return("\n Not enough adducts for this ionization mode specified - skipped ...")
	}			
	peaklist2<-as.data.frame(peaklist[peaklist[,"keep"]==1,c("m/z_corr","int_corr","RT_corr","peak_ID")])	
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
		return("\n Adduct detection failed - adapt parameters?");
	}				
	if(length(adduct[["Pairs"]][,1])==0){
		return("\n No adduct relations detected");
	}
	Adduct_pairs<-adduct[["Pairs"]][,c(1,2)]
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
	##########################################################################		
	return(NULL)

}








