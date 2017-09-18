#' @title Start new enviMass workflow project
#'
#' @description \code{newproject} initializes an enviMass project
#'
#' @param pro_name Character string of project name.
#' @param pro_dir Character string. Path to a valid folder that will contains the project folder.
#' @param IS Dataframe containing internal standard (IS) compounds
#' @param targets Dataframe containing target and suspect compounds
#'
#' @details enviMass workflow function. Creates a project folder, various subfolders, a logfile, dummy outputs, etc
#' 

pattern_search3_wrap<-function(
	x,
	logfile,
	...
){


	##########################################################################
	# LOAD FILES & REMOVE OLD RESULTS ########################################
	cat("loading - ")
	for_file<-x
	if( file.exists(file.path(logfile[[1]],"results","componentization","isotopologues",for_file) ) ){
		file.remove(file.path(logfile[[1]],"results","componentization","isotopologues",for_file) )
	}			
	load(file=file.path(logfile[[1]],"peaklist",as.character(for_file))); # Peaklist  
	peaklist<-peaklist[order(peaklist[,10],decreasing=FALSE),] # match with IDs - for saving pattern; IDs are retrieved for pairs seperately
	##########################################################################
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
	peaklist2<-as.data.frame(peaklist[peaklist[,"keep"]==1,c("m/z_corr","int_corr","RT_corr","peak_ID")])	
	if(logfile$parameters$isotop_ppm=="TRUE"){
		use_mztol<-as.numeric(logfile$parameters$isotop_mztol)
	}else{ # mmu
		use_mztol<-(as.numeric(logfile$parameters$isotop_mztol)/1000)
	}
	pattern<-try(
		enviMass::pattern_search3(
			peaklist=peaklist2[,c("m/z_corr","int_corr","RT_corr","peak_ID")],
			quantiz,
			mztol=use_mztol,
			ppm=logfile$parameters$isotop_ppm,
			inttol=(as.numeric(logfile$parameters$isotop_inttol)/100),
			rttol=as.numeric(logfile$parameters$isotop_rttol),
			use_isotopes=FALSE,
			use_charges=logfile$parameters$isotop_use_charges,
			use_marker=TRUE,
			quick=TRUE,
			isotopes,
			exclude
		)
	)
	if(class(pattern)=="try-error"){
		return("\n Isotopologue detection failed - adapt parameters?");
	}				
	if(length(pattern[["Pairs"]][,1])==0){
		return("\n No adduct relations detected");
	}
	Isot_pairs<-pattern[["Pairs"]]
	pattern[["Pairs"]]<-0
	those<-(Isot_pairs[,1]>Isot_pairs[,2])
	if(any(those)){
		Isot_pairs[those,]<-Isot_pairs[those,c(2,1)]
	}
	Isot_pairs<-Isot_pairs[order(Isot_pairs[,1],Isot_pairs[,2],decreasing=FALSE),]				
	save(Isot_pairs,file=(file.path(logfile[[1]],"results","componentization","isotopologues",paste(for_file,sep="_"))))
	save(pattern,file=(file.path(logfile[[1]],"results","componentization","isotopologues",paste("full",for_file,sep="_"))))
	##########################################################################	
	rm(peaklist,peaklist2,pattern,Isot_pairs)			
	cat("done.")
	##########################################################################		
 	return("done");
  
}
