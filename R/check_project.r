#' @title Check enviMass project
#'
#' @description \code{} runs several cconsistency checks on a enviMass project
#'
#' @param logfile enviMass project logfile
#' @param isotopes Isotope list
#' @param adducts Adducts list from package enviPat
#' @param skipcheck Logical. Should project check be skipped?
#' @param ignorefiles Logical. Skip checks involving .mzXML or MSlist files
#' @param write_tables Logical. If TRUE, corrected chemical formulas are written into tables (in check_compounds) 
#'
#' @return Character string with either specific error message or a "Project consistent" message.
#' 
#' @details enviMass workflow function; run before further calculations are started in the workflow.
#' 

check_project <- function(
	isotopes,
	adducts,
	skipcheck = FALSE,
	ignorefiles = FALSE,
	write_tables = FALSE,
	...
){
	
	say <- "Project consistent"
	if(skipcheck){
		return(say);
	}
	if(any(ls() == "logfile")){stop("\n illegal logfile detected #1 in check_project.r!")}
	###############################################################################	
	# on adduct table #############################################################
	formula_add <- adducts$Formula_add[adducts$Formula_add != "FALSE"]
	if(length(formula_add)){
		formula_add2 <- enviPat::check_chemform(isotopes, formula_add )
		if(any(formula_add2$warning)){
			say <- paste0("Wrong Formula_add in adduct table found: ",
				paste(formula_add2$new_formula[formula_add2$warning], collapse = ", "),
				". Please revise!"
			)	
		}else{ # insert corrected formulas: e.g., missing numbers
			k <- which(formula_add2$new_formula != formula_add)
			if(length(k)){
				adducts$Formula_add[adducts$Formula_add != "FALSE"] <<- formula_add2$new_formula
			}
		}
	}
	formula_ded <- adducts$Formula_ded[adducts$Formula_ded != "FALSE"]
	if(length(formula_ded)){
		formula_ded2 <- enviPat::check_chemform(isotopes, formula_ded )
		if(any(formula_ded2$warning)){
			say <- paste0("Wrong Formula_ded in adduct table found: ",
				paste(formula_ded2$new_formula[formula_ded2$warning], collapse = ", "),
				". Please revise!"
			)	
		}else{ # insert corrected formulas: e.g., missing numbers
			k <- which(formula_ded2$new_formula != formula_ded)
			if(length(k)){
				adducts$Formula_ded[adducts$Formula_ded != "FALSE"] <<- formula_ded2$new_formula
			}
		}
	}	
	rm(formula_ded, formula_add)
	###############################################################################	
	# on project folders ##########################################################
	if(!all(names(logfile) == c(
		"project_folder",
		"Tasks_to_redo",
		"summary",
		"PW MSconvert path",
		"parameters",
		"workflow",
		"adducts_pos",
		"adducts_neg",
		"isotopes",
		"version",
		"workflow_depend",
		"workflow_must",
		"Positive_subtraction_files",
		"Negative_subtraction_files",
		"adducts_pos_group",
		"adducts_neg_group",
		"UI_options",
		"method_setup"))
	){
		say <- "Your logfile is corrupted_1,  project debug required!"
	}
	###############################################################################
	# on parallelization: #########################################################
	if( (logfile$parameters$parallel == "TRUE") & (logfile$parameters$parallel_restrict == "TRUE") ){
		available_num_cores <- detectCores(all.tests = FALSE, logical = TRUE)
		if(available_num_cores < as.numeric(logfile$parameters$parallel_cores)){
			say <- "Selected number of cores/threads to restrict to for parallel processing exceeds available number for your computer.
				Please revise settings under Settings -> General -> Multi-core processing."
		}
	}
	###############################################################################
	# wrong upstream "must not" executions? #######################################
	must <- logfile[[12]]
	for(i in 1:length(must[1,])){
		for(j in 1:length(must[,i])){	
			if(must[j,i] == -1){
				if(logfile$workflow[names(logfile$workflow) == colnames(must)[i]] == "yes"){
					if(logfile$workflow[names(logfile$workflow) == rownames(must)[j]] == "yes"){
						say <- paste("wokflow step",names(logfile$workflow)[i],"excludes",names(logfile$workflow)[j],"- adapt workflow settings!")
					}
				}
			}		
		}
	}
	##############################################################################
	# directories available? ##################################################### 
	if(!file.exists(file.path(logfile[[1]], "files"))& ignorefiles=="FALSE"){say <- "files directory missing!"}
	if(!file.exists(file.path(logfile[[1]], "MSlist"))& ignorefiles=="FALSE"){say <- "MSlist directory missing!"}  
	if(!file.exists(file.path(logfile[[1]], "MSraw"))& ignorefiles=="FALSE"){say <- "MSraw directory missing!"} 
	if(!file.exists(file.path(logfile[[1]], "features"))){say <- "features directory missing!"}
	if(!file.exists(file.path(logfile[[1]], "results"))){say <- "results directory missing!"}
	if(!file.exists(file.path(logfile[[1]], "results", "screening"))){say <- "results/screening directory missing!"} 
	if(!file.exists(file.path(logfile[[1]], "quantification"))){say <- "results/quantification directory missing!"} 
	if(!file.exists(file.path(logfile[[1]], "results", "LOD"))){say <- "results/LOD directory missing!"} 
	if(!file.exists(file.path(logfile[[1]], "results", "recalibration"))){say <- "results/recalibration directory missing!"} 	
	if(!file.exists(file.path(logfile[[1]], "dataframes"))){say <- "dataframes directory missing!"}
	if(!file.exists(file.path(logfile[[1]], "pics"))){say <- "pics directory missing!"}
	if(!file.exists(file.path(logfile[[1]], "exports"))){say <- "exports directory missing!"}  
	##############################################################################
	# compounds available & ok? ##################################################
	intstand_check <- read.table(file = file.path(logfile[[1]], "dataframes", "IS.txt"), header = TRUE, sep = "\t", colClasses = "character", blank.lines.skip = TRUE);
	targets_check <- read.table(file = file.path(logfile[[1]], "dataframes", "targets.txt"), header = TRUE, sep = "\t", colClasses = "character", blank.lines.skip = TRUE);
	say1 <- enviMass::check_compounds(intstand_check, targets_check, isotopes, adducts, logfile, write_tables = TRUE)
	if(say1 != "Project consistent"){say <- say1}
	if(any(ls() == "logfile")){stop("\n illegal logfile detected #2 in check_project.r!")}		
	# adducts available? #########################################################
	check_adducts_pos <- logfile[["adducts_pos"]]
	check_adducts_pos <- check_adducts_pos[check_adducts_pos != "FALSE"]
	if(length(check_adducts_pos)){	
		found_add <- match(check_adducts_pos, adducts[adducts[,"Ion_mode"] == "positive", "Name"])
		if(any(is.na(found_add))){
			say <- paste("Invalid screening adducts found for positive mode, i.e., this/these saved adducts are not present in the Settings -> Screening adduct table: \n",
						paste(check_adducts_pos[is.na(found_add)], collapse = ", "),
						" \nEither provide a valid adduct table containing this adduct or select new adducts in the Screening Settings and press Apply",
						sep = "")
		}
	}
	check_adducts_neg <- logfile[["adducts_neg"]]
	check_adducts_neg <- check_adducts_neg[check_adducts_neg != "FALSE"]
	if(length(check_adducts_neg)){	
		found_add <- match(check_adducts_neg, adducts[adducts[,"Ion_mode"] == "negative", "Name"])
		if(any(is.na(found_add))){
			say <- paste("Invalid screening adducts found for negative mode, i.e., this/these saved adducts are not present in the Settings -> Screening adduct table: \n",
						paste(check_adducts_pos[is.na(found_add)], collapse = ", "),
						" \nEither provide a valid adduct table containing this adduct or select new adducts in the Screening Settings and press Apply",
						sep = "")
		}
	}
	check_adducts_pos <- logfile[["adducts_pos_group"]]
	check_adducts_pos <- check_adducts_pos[check_adducts_pos != "FALSE"]
	if(length(check_adducts_pos)){	
		found_add <- match(check_adducts_pos, adducts[adducts[,"Ion_mode"] == "positive", "Name"])
		if(any(is.na(found_add))){
			say <- paste("Invalid screening adducts found for positive mode, i.e., this/these saved adducts are not present in the Settings -> Screening adduct table: \n",
						paste(check_adducts_pos[is.na(found_add)], collapse = ", "),
						" \nEither provide a valid adduct table containing this adduct or select new adducts in the Screening Settings and press Apply",
						sep = "")
		}
	}
	check_adducts_neg <- logfile[["adducts_neg_group"]]
	check_adducts_neg <- check_adducts_neg[check_adducts_neg != "FALSE"]
	if(length(check_adducts_neg)){	
		found_add <- match(check_adducts_neg, adducts[adducts[,"Ion_mode"] == "negative", "Name"])
		if(any(is.na(found_add))){
			say <- paste("Invalid screening adducts found for negative mode, i.e., this/these saved adducts are not present in the Settings -> Screening adduct table: \n",
						paste(check_adducts_pos[is.na(found_add)], collapse = ", "),
						" \nEither provide a valid adduct table containing this adduct or select new adducts in the Screening Settings and press Apply",
						sep = "")
		}
	}
	# enough compounds for recalibration available? ##############################
	if(logfile$workflow[names(logfile$workflow) == "recal"] == "yes"){
		# check for positive mode
		if(logfile$parameters$recal_include_pos == "TRUE"){
			if(logfile$parameters$recal_use_pos == "Internal standards"){
				IS <- read.table(file=file.path(logfile[[1]],"dataframes","IS.txt"),header=TRUE,sep="\t",colClasses = "character");
				IS <- IS[IS[,"ion_mode"] == "positive",, drop = FALSE]
				if(length(IS[IS[,"use_for_recalibration"] == "TRUE", 1]) < 10){
					say <- "Not enough internal standards available for mass recalibration in positive mode ... revise, maybe exlude mass recalibration for the positive mode only (Settings -> Recalibration)?"   
				}
			}
			if(logfile$parameters$recal_use_pos == "Target compounds"){
				targets <- read.table(file = file.path(logfile[[1]], "dataframes","targets.txt"), header = TRUE, sep = "\t", colClasses = "character");
				targets <- targets[targets[,"ion_mode"] == "positive",, drop = FALSE]
				if(length(targets[targets[,"use_for_recalibration"] == "TRUE", 1]) < 10){
					say <- "Not enough target compounds available for mass recalibration in positive mode ... revise, maybe exlude mass recalibration for the positive mode only (Settings -> Recalibration)?"    
				}
			}
			if(logfile$parameters$recal_use_pos == "both"){
				IS <- read.table(file = file.path(logfile[[1]], "dataframes", "IS.txt"), header = TRUE, sep = "\t", colClasses = "character");	  
				IS <- IS[IS[,"ion_mode"] == "positive",, drop = FALSE]
				a <- length(IS[IS[,8] == "TRUE",1])
				targets <- read.table(file = file.path(logfile[[1]], "dataframes", "targets.txt"), header = TRUE, sep = "\t", colClasses = "character");  
				targets <- targets[targets[,"ion_mode"] == "positive",, drop = FALSE]
				b <- length(targets[targets[,9] == "TRUE", 1])
				if((a + b) < 10){
					say <- "Not enough target compounds + internal standards available for mass recalibration in positive mode ... revise, maybe exlude mass recalibration for the positive mode only (Settings -> Recalibration)?"    
				}
			}
		}
		# check for negative mode
		if(logfile$parameters$recal_include_neg=="TRUE"){		
			if(logfile$parameters$recal_use_neg=="Internal standards"){
				IS<-read.table(file=file.path(logfile[[1]],"dataframes","IS.txt"),header=TRUE,sep="\t",colClasses = "character");
				IS<-IS[IS[,"ion_mode"]=="negative",,drop=FALSE]
				if(length(IS[IS[,"use_for_recalibration"]=="TRUE",1])<10){
					say <- "Not enough internal standards available for mass recalibration in negative mode ... revise, maybe exlude mass recalibration for the negative mode only (Settings -> Recalibration)?"    
				}
			}
			if(logfile$parameters$recal_use_neg=="Target compounds"){
				targets <- read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character");
				targets <- targets[targets[,"ion_mode"] == "negative",, drop = FALSE]
				if(length(targets[targets[,"use_for_recalibration"] == "TRUE", 1]) < 10){
					say <- "Not enough target compounds available for mass recalibration in negative mode ... revise, maybe exlude mass recalibration for the negative mode only (Settings -> Recalibration)?"    
				}
			}
			if(logfile$parameters$recal_use_neg=="both"){
				IS <- read.table(file = file.path(logfile[[1]], "dataframes", "IS.txt"), header = TRUE, sep = "\t", colClasses = "character");	  
				IS <- IS[IS[,"ion_mode"] == "positive",, drop = FALSE]
				a <- length(IS[IS[,8] == "TRUE", 1])
				targets <- read.table(file = file.path(logfile[[1]], "dataframes", "targets.txt"), header = TRUE, sep = "\t", colClasses = "character");  
				targets <- targets[targets[,"ion_mode"] == "negative",, drop = FALSE]
				b <- length(targets[targets[,9] == "TRUE", 1])
				if((a + b) < 10){
					say <- "Not enough target compounds + internal standards available for mass recalibration in negative mode ... revise, maybe exlude mass recalibration for the negative mode only (Settings -> Recalibration)?"    
				}
			}	
		}
	}
	##############################################################################
	# parameters ok? #############################################################
	# (1) on trend time lags #####################################################
	if(logfile$workflow[names(logfile$workflow) == "trendblind"] == "yes"){
		lags<-as.numeric(strsplit(as.character(logfile$parameters$trend_lags), ",")[[1]])
		if(any(is.na(lags))){say< - "Invalid trend lags - have you used comma separated numerics?"}
		measurements<-read.csv(file = file.path(logfile[[1]], "dataframes", "measurements"), colClasses = "character");
		if(!any(measurements[,"include"] == "TRUE")){
			say <- "No file included into workflow; nothing to be calculated."
			return(say)
		}
		measurements <- measurements[measurements[,"include"] == "TRUE",]
		dated <- measurements[,"Date"]
		timed <- measurements[,"Time"]
		datetime <- c()
		if (length(timed) > 0) {
			for(i in 1:length(timed)){
				datetime<-c(datetime, paste(dated[i], timed[i], "CET", sep = " "))
			}
			atPOSIX <- as.POSIXct(datetime);
			atPOSIX <- as.numeric(atPOSIX)
			if(min(lags) > (((max(atPOSIX) - min(atPOSIX)) / (24 * 60 * 60)) + 1)){say <- "Trend lags longer than time span of the measurements ... abort"}
		}
		rm(lags); ##############################################################################
		if(
			(logfile$workflow[names(logfile$workflow)=="profiling"]=="yes") &
			( !any(measurements[,names(measurements)=="profiled"]=="TRUE") ) &
			(logfile$parameters$prof_select=="TRUE")
		){
			say<-"Workflow option profiling and settings->profile restriction enabled, but no file included as TRUE in measurement table for profiling."
		}
	}
	if(logfile$parameters$screen_IS_restrict=="TRUE"){
		if(is.na(as.numeric(logfile$parameters$screen_IS_restrict_many)) | (as.numeric(logfile$parameters$screen_IS_restrict_many)<1)){
			say<-"ISTD screening: number of latest files to include invalid - must be >=1. Please revise!"
		}
	}
	if(logfile$parameters$screen_target_restrict=="TRUE"){
		if(is.na(as.numeric(logfile$parameters$screen_target_restrict_many)) | (as.numeric(logfile$parameters$screen_target_restrict_many)<1)){
			say<-"Target screening: number of latest files to include invalid - must be >=1. Please revise!"
		}
	}
	if((logfile$workflow[names(logfile$workflow)=="homologues"]=="yes") & (logfile$parameters$homol_units[1]!="FALSE")){
		these<-enviPat::check_chemform(isotopes,strsplit(logfile$parameters$homol_units,",")[[1]])
		if(any(these[,1]=="TRUE")){
			say<-"Invalid homologue units defined - please revise (empty spaces? not comma-seperated?)"
		}
	}
	# data sets ok? ##############################################################
	filed<-list.files(file.path(logfile[[1]],"files"))
	if(!length(filed) & ignorefiles == "FALSE"){say <- "No files available!"}
	measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"), colClasses = "character")
	if(length(names(measurements)) != 28){
		say <- "Measurement table seems corrupted. Have you made any updates recently? Please report this issue!"
	}
	if(any(measurements[,"Mode"] == "positive")){
		if(all(is.na(match(measurements[measurements[,"Mode"] == "positive","Type"],c("sample","blank"))))){
			say <- "No blanks or sample files included - hence, nothing to do. Add at least one to continue, positive mode."
		}	
	}
	if(any(measurements[,"Mode"] == "negative")){
		if(all(is.na(match(measurements[measurements[,"Mode"] == "negative","Type"],c("sample","blank"))))){
			say <- "No blanks or sample files included - hence, nothing to do. Add at least one to continue, negative mode."
		}	
	}
	##############################################################################
	# blind peak subtraction enabled, but no blind files selected? ###############
	if(logfile$workflow[names(logfile$workflow) == "blind"] == "yes"){
		if(!any(measurements[,"Type"] == "blank")){
			say <- "You have included the blind/blank subtraction in the workflow - but there are no blind/blank files. Please revise!"
		}
		if(
			(logfile$parameters$subtract_pos_byfile == "TRUE") & 
			all(logfile$Positive_subtraction_files == "FALSE") 
		){
			say <- "You have included the blind/blank subtraction in the workflow, and have chosen the blind settings option Additional non-sample files to subtract (positive ionization). But you missed to mark any files underneath that selection. Please revise!"
		}
		if(
			(logfile$parameters$subtract_neg_byfile == "TRUE") & 
			all(logfile$Negative_subtraction_files == "FALSE") 
		){
			say <- "You have included the blind/blank subtraction in the workflow, and have chosen the blind settings option Additional non-sample files to subtract (negative ionization). But you missed to mark any files underneath that selection. Please revise!"
		}		
	}
	# check: do calibration?
	if( 
		(logfile$workflow[names(logfile$workflow)=="calibration"]=="yes") &
		!any(measurements[,"Type"]=="calibration")
	){
		say<-"You have included the calibration in the workflow - but there are no calibration files. Please revise!"
	}
	# adduct grouping - but only one adduct selected?
	if( 
		(logfile$workflow[names(logfile$workflow)=="adducts"]=="yes") 
	){
		if( (any(measurements[,"Mode"]=="positive")) & (length(logfile$adducts_pos_group)<2) ){
			say<-"You want to group nontarget adducts with the workflow - but have <2 adducts specified in the Settings-> Componentization -> Adducts (positive mode). Please revise!"
		}
		if( (any(measurements[,"Mode"]=="negative")) & (length(logfile$adducts_neg_group)<2) ){
			say<-"You want to group nontarget adducts with the workflow - but have <2 adducts specified in the Settings-> Componentization -> Adducts (negative mode). Please revise!"
		}
	}
	if(length(measurements[measurements[,"ID"]!="-",1,drop=FALSE])==0){
		say<-"No files available."; 
		return(say);
	}
	if(any(duplicated(measurements[,"ID"]))){
		say<-paste("Duplicated file IDs.",measurements[duplicated(measurements[,"ID"]),1],"Revise!")
	}
	if(any(is.na(as.numeric(measurements[,"ID"])))){
		these<-which(is.na(as.numeric(measurements[,"ID"])))
		say<-paste("Non-numeric file IDs (",measurements[these,1],"). Revise entry in or delete entry from the enviMass file table!",sep="")
	}  
	measurements_ID<-measurements[,"ID"]
	measurements_ID<-paste(measurements_ID,".mzXML",sep="")
	if(any(match(measurements_ID,filed,nomatch=0)==0) & ignorefiles=="FALSE"){
		these<-which(match(measurements_ID,filed,nomatch=0)==0)
		say<-paste("Missing mzXML file for files with ID: ", paste(these,collapse=", "),". Revise - best delete the concerned file from enviMass file table and reload it.",sep="")
	}
	# check date  & time formats
	a<-try({as.Date(measurements[,"Date"], tz="GMT")},silent=TRUE)
	if(any(class(a)=="try-error" | is.na(a))){
		these<-which(class(a)=="try-error" | is.na(a))
		these<-measurements[these,1]
		say<-paste("Invalid date format found for file(s) with ID(s) ",
		paste(these,collapse=", "),". Please revise concerned file(s) in the files tab!",sep="")
	}
	b<-try({as.Date(paste("2014-03-13",measurements[,"Time"]),format= "%Y-%m-%d %H:%M:%S", tz="GMT")},silent=TRUE)  
	if(any(class(b)=="try-error" | is.na(b))){
		these<-which(class(b)=="try-error" | is.na(b))
		these<-measurements[these,1]
		say<-paste("Invalid time format found for file(s) with ID(s) ",
		paste(these,collapse=", "),". Please revise concerned file(s) in the files tab (- if there are any loaded so far)!",sep="")
	}
  	# check calibration files
  measurements_cal<-measurements[measurements[,"Type"]=="calibration",,drop=FALSE] 
  if(length(measurements_cal[,1])>0){
	  a<-try({as.Date(c(measurements_cal[,"date_end"]), tz="GMT")},silent=TRUE)
	  if(any(class(a)=="try-error" | is.na(a))){
			these<-which(class(a)=="try-error" | is.na(a))
			these<-measurements_cal[these,1]
			say<-paste("Invalid end date format found for calibration file(s) with ID(s) ",
			paste(these,collapse=", "),". Please revise concerned calibration file(s) in the files tab!",sep="")
	  }
	  b<-try({as.Date(paste("2014-03-13",measurements_cal[,"time_end"]),format= "%Y-%m-%d %H:%M:%S", tz="GMT")},silent=TRUE)  
	  if(any(class(b)=="try-error" | is.na(b))){
			these<-which(class(b)=="try-error" | is.na(b))
			these<-measurements_cal[these,1]
			say<-paste("Invalid time format found for calibration file(s) with ID(s) ",
			paste(these,collapse=", "),". Please revise concerned calibration file(s) in the files tab!",sep="")
	  }
	  if(any(is.na(as.numeric(measurements_cal$tag1)))){
		say<-paste("Invalid non-numeric concentration (tag1) found for calibration file with ID:",
			measurements_cal[which(is.na(as.numeric(measurements_cal$tag1))),"ID"][1],"- please revise",
			sep=" ")
	  }

  }
	if(logfile$workflow[names(logfile$workflow)=="quantification"]=="yes" & any(measurements[,"Type"]=="calibration")){
		# no period overlaps! ######################################################
		# -> positive 
		cal_files<-measurements[(measurements[,"Mode"]=="positive")&(measurements[,"Type"]=="calibration"),,drop=FALSE]
		if(length(cal_files[,1])>0){
			# do calibration sets overlap in periods?
			cal_files2<-unique(cal_files[,c("tag3","Date","Time","date_end","time_end"),drop=FALSE])
			starttime<-as.difftime(cal_files2[,3]);startdate<-as.Date(cal_files2[,2], tz="GMT");
			numstart<-(as.numeric(startdate)+as.numeric(starttime/(24*60*60)))		
			endtime<-as.difftime(cal_files2[,5]);enddate<-as.Date(cal_files2[,4], tz="GMT");
			numend<-(as.numeric(enddate)+as.numeric(endtime/(24*60*60)))
			if(length(starttime)>1){			
				for(i in 1:(length(starttime)-1)){
					for(j in (i+1):length(starttime)){		
						if(
							(numstart[j]<=numend[i])&
							(numstart[i]<=numend[j])
						){
							say<-paste("Time periods of calibration set ",cal_files2[i,"tag2"]," and ",cal_files2[j,"tag2"]," (positive mode) overlap. Time periods of different calibration files sets must not overlap - please revise!",sep="")
						}
					}
				}
			}
			for(i in 1:length(starttime)){			
				if(numstart[i]>=numend[i]){
						say<-paste("Start Date/Time >= end Date/Time for calibration set ",cal_files2[i,"tag2"]," (positive mode). Please revise!",sep="")				
				}
			}
			# do all files in one calibration set have identical start & end times?		
			tags2<-unique(cal_files[,"tag2"])
			for(i in 1:length(tags2)){
				if(length(unique(cal_files[cal_files[,"tag2"]==tags2[i],"Date"]))>1){
					say<-paste("Positive mode calibration files of set ",tags2[i]," must have identical start dates - but they do not. Please revise!",sep="")
				}
				if(length(unique(cal_files[cal_files[,"tag2"]==tags2[i],"Time"]))>1){
					say<-paste("Positive mode calibration files of set ",tags2[i]," must have identical start times - but they do not. Please revise!",sep="")
				}		
				if(length(unique(cal_files[cal_files[,"tag2"]==tags2[i],"date_end"]))>1){
					say<-paste("Positive mode calibration files of set ",tags2[i]," must have identical end dates - but they do not. Please revise!",sep="")
				}
				if(length(unique(cal_files[cal_files[,"tag2"]==tags2[i],"time_end"]))>1){
					say<-paste("Positive mode calibration files of set ",tags2[i]," must have identical end times - but they do not. Please revise!",sep="")
				}		
			}
		}	
		# -> negative
		cal_files<-measurements[(measurements[,"Mode"]=="negative")&(measurements[,"Type"]=="calibration"),,drop=FALSE]	
		if(length(cal_files[,1])>0){
			# do calibration sets overlap in periods?
			cal_files2<-unique(cal_files[,c("tag3","Date","Time","date_end","time_end"),drop=FALSE])
			starttime<-as.difftime(cal_files2[,3]);startdate<-as.Date(cal_files2[,2], tz="GMT");
			numstart<-(as.numeric(startdate)+as.numeric(starttime/(24*60*60)))		
			endtime<-as.difftime(cal_files2[,5]);enddate<-as.Date(cal_files2[,4], tz="GMT");
			numend<-(as.numeric(enddate)+as.numeric(endtime/(24*60*60)))	
			if(length(starttime)>1){
				for(i in 1:(length(starttime)-1)){
					for(j in (i+1):length(starttime)){		
						if(
							(numstart[j]<=numend[i])&
							(numstart[i]<=numend[j])
						){
							say<-paste("Time periods of calibration set ",cal_files2[i,"tag2"]," and ",cal_files2[j,"tag2"]," (negative mode) overlap. Time periods of different calibration files sets must not overlap - please revise!",sep="")
						}
					}
				}
			}
			for(i in 1:length(starttime)){			
				if(numstart[i]>=numend[i]){
						say<-paste("Start Date/Time >= end Date time for calibration set ",cal_files2[i,"tag2"]," (negative mode) overlap. Please revise!",sep="")				
				}
			}
			# do all files in one calibration set have identical start & end times?		
			tags2<-unique(cal_files[,"tag2"])
			for(i in 1:length(tags2)){
				if(length(unique(cal_files[cal_files[,"tag2"]==tags2[i],"Date"]))>1){
					say<-paste("Negative mode calibration files of set ",tags2[i]," must have identical start dates - but they do not. Please revise!",sep="")
				}
				if(length(unique(cal_files[cal_files[,"tag2"]==tags2[i],"Time"]))>1){
					say<-paste("Negative mode calibration files of set ",tags2[i]," must have identical start times - but they do not. Please revise!",sep="")
				}		
				if(length(unique(cal_files[cal_files[,"tag2"]==tags2[i],"date_end"]))>1){
					say<-paste("Negative mode calibration files of set ",tags2[i]," must have identical end dates - but they do not. Please revise!",sep="")
				}
				if(length(unique(cal_files[cal_files[,"tag2"]==tags2[i],"time_end"]))>1){
					say<-paste("Negative mode calibration files of set ",tags2[i]," must have identical end times - but they do not. Please revise!",sep="")
				}		
			}
		}	
	}
	measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
	if(logfile$workflow[names(logfile$workflow)=="quantification"]=="yes" & !any(measurements[,"Type"]=="calibration")){
		say<-"You want to run a quantification but you have not even calibration files in your project? Either add such files or remove the Quantification step from your project!"
	} 
	if(logfile$workflow[names(logfile$workflow)=="calibration"]=="yes" & !any(measurements[,"Type"]=="calibration")){
		say<-"You want to run a calibration but you have no calibration files in your project? Either add such files or remove the Calibration step from your project!"
	} 
	if(logfile$workflow[names(logfile$workflow)=="recovery"]=="yes" & !any(measurements[,"Type"]=="calibration")){
		say<-"You want to run a recovery but you have not even calibration files in your project for a quantification? Either add such files or remove the Recovery step from your project!"
	}	   
	# check spiked files
	measurements_spiked<-measurements[ (measurements[,"Type"]=="spiked" & measurements[,"include"]=="TRUE"),,drop=FALSE]   
	if(
		(length(measurements_spiked[,"ID"]) > 0) & 
		(logfile$workflow[names(logfile$workflow) == "recovery"] == "yes")
	){  
		these_pos<-which(is.na(match(
			measurements_spiked[measurements_spiked$Mode=="positive",]$tag2,
			measurements[measurements$Mode=="positive","ID"]))
		)
		if(length(these_pos)>0){
			these_pos<-measurements_spiked[measurements_spiked$Mode=="positive","ID"][these_pos]
		}
		these_neg<-which(is.na(match(
			measurements_spiked[measurements_spiked$Mode=="negative",]$tag2,
			measurements[measurements$Mode=="negative","ID"]))
		)
		if(length(these_neg)>0){	
			these_neg<-measurements_spiked[measurements_spiked$Mode=="negative","ID"][these_neg]
		}
		if(length(these_pos)>0 || length(these_neg)>0){
			say<-paste("Invalid file IDs (tag2) to subtract from for spiked file(s) with ID(s) ",
			paste(c(these_pos,these_neg),collapse=", "),". Please revise concerned spiked file(s) in the files tab!",sep="")	
		}
	}
	##############################################################################
	# progress bar? ##############################################################
	if(interactive() && !.Platform$OS.type == "windows" && .Platform$GUI == "Rgui" && logfile[[5]][21]=="TRUE"){
		say<-"Disable the progress bar in the Settings General Tab; works only under Windows OS"
	}
	##############################################################################
	# Isotopologue grouping - quantiz data set available? ######################## 
	# Download from www.envimass.ch ##############################################
	redo_load_quantiz<-FALSE
	if( # quantiz file in project? check if it can be loaded and agrees with selected instrument&resolution
		(logfile$workflow[names(logfile$workflow)=="isotopologues"]=="yes") &
		file.exists(file.path(logfile[[1]],"dataframes","quantiz") )
	){
		load_quantiz<-try(load(file.path(logfile[[1]],"dataframes","quantiz")))  
		if(class(load_quantiz)=="try-error"){
			redo_load_quantiz<-TRUE
		}else{
			if(quantiz$R_set!="Sciex_all"){
				if(quantiz$R_set!=logfile$parameters$resolution){
					redo_load_quantiz<-TRUE
				}
			}else{ # Sciex instruments all set by one quantiz simulation
				if(!grepl("Sciex",logfile$parameters$resolution)){
					redo_load_quantiz<-TRUE
				}
			}
		}
	}
	if( # quantiz not available or not correct?
		(logfile$workflow[names(logfile$workflow)=="isotopologues"]=="yes") &
		( !file.exists(file.path(logfile[[1]],"dataframes","quantiz") ) || redo_load_quantiz)
	){
		avail<-c(
			"OrbitrapXL,Velos,VelosPro_R60000@400",
			"Q-Exactive,ExactivePlus_280K@200",
			"Q-Exactive,ExactivePlus_R140000@200",
			"Q-Exactive,ExactivePlus_R70000@200",
			"Sciex_TripleTOF5600_R25000@200",
			"Sciex_TripleTOF6600_R25000@200",
			"Sciex_QTOFX500R_R25000@200",
			"Agilent_QTOF6550_low_extended2GHz_highRes",
			"Agilent_QTOF6550_low_highRes4GHz_highRes"
		)
		if(any(avail == logfile$parameters$resolution)){ # available on www.envimass.ch
			if(logfile$parameters$resolution == "OrbitrapXL,Velos,VelosPro_R60000@400"){
				get_url <- "http://www.looscomputing.ch/eng/enviMass/inputs/quantiz/OrbitrapXL_Velos_VelosPro_R60K@400/quantiz"
			}
			if(logfile$parameters$resolution=="Q-Exactive,ExactivePlus_280K@200"){
				get_url <- "http://www.looscomputing.ch/eng/enviMass/inputs/quantiz/Q-Exactive_ExactivePlus_280K@200/quantiz"
			}			
			if(logfile$parameters$resolution=="Q-Exactive,ExactivePlus_R140000@200"){
				get_url <- "http://www.looscomputing.ch/eng/enviMass/inputs/quantiz/Q-Exactive_ExactivePlus_R140K@200/quantiz"
			}				
			if(logfile$parameters$resolution=="Q-Exactive,ExactivePlus_R70000@200"){
				get_url <- "http://www.looscomputing.ch/eng/enviMass/inputs/quantiz/Q-Exactive_ExactivePlus_R70K@200/quantiz"
			}				
			if(logfile$parameters$resolution=="Sciex_TripleTOF5600_R25000@200"){
				get_url <- "http://www.looscomputing.ch/eng/enviMass/inputs/quantiz/Sciex_all/quantiz"
			}
			if(logfile$parameters$resolution=="Sciex_TripleTOF6600_R25000@200"){
				get_url <- "http://www.looscomputing.ch/eng/enviMass/inputs/quantiz/Sciex_all/quantiz"
			}
			if(logfile$parameters$resolution=="Sciex_QTOFX500R_R25000@200"){
				get_url <- "http://www.looscomputing.ch/eng/enviMass/inputs/quantiz/Sciex_all/quantiz"
			}
			if(logfile$parameters$resolution=="Agilent_QTOF6550_low_extended2GHz_highRes"){
				get_url <- "http://www.looscomputing.ch/eng/enviMass/inputs/quantiz/Agilent_QTOF6550_low_extended2GHz_highRes/quantiz"
			}			
			if(logfile$parameters$resolution=="Agilent_QTOF6550_low_highRes4GHz_highRes"){
				get_url <- "http://www.looscomputing.ch/eng/enviMass/inputs/quantiz/Agilent_QTOF6550_low_highRes4GHz_highRes/quantiz"
			}	
			dest_file <- file.path(logfile[[1]], "dataframes", "quantiz")
			url_quantiz <- try(download.file(url=get_url, destfile = dest_file, mode = "wb"))
			if(class(url_quantiz) == "try-error"){
				cat("\n Download of missing isotopologue space failed.")
				say <- "Unable to download missing isotopologue space. Please either check your internet connection and retry or
				proceed manually as described on www.enviMass.ch -> Data input -> Download available isotopologue spaces." 
			}else{
				cat("\n Download of missing isotopologue space completed.\n")			  
				load_quantiz <- try(load(file.path(logfile[[1]], "dataframes", "quantiz")) )
				if(class(load_quantiz) == "try-error"){
					cat("\n Loading of missing isotopologue space failed.\n")	
					say <- "Loading failure of downloaded isotopologue space. Please proceed manually as described on www.enviMass.ch -> Data input -> Download available isotopologue spaces."
				}
			}
		}else{	# not available on www.wnvimass.ch
			say <- "The isotopologue grouping step is enabled in the workflow. We regret that no data set for the isotopologue space for the selected instrument and resolution is at present available for enviMass. Please remove the isotopologue grouping step from the workflow." 
		} 
	}
	###############################################################################
	# Homologue series detection ################################################## 
	if(logfile$workflow[names(logfile$workflow) == "homologues"] == "yes"){
		if(logfile$parameters$external$homol_units[1] != "FALSE"){
			these <- enviPat::check_chemform(isotopes, logfile$parameters$external$homol_units)[,1] 
			if(any(these != "FALSE")){
				say <- "Invalid chemical formulas for predefined homologue series units found - please revise" 
			}  
		}
	} 
	###############################################################################
	# Method setup ################################################################
	if( as.logical(logfile$parameters$method_use) & is.character(logfile$method_setup) ){
		say <- "Method setup is enabeled but no valid existing method is available. Either disable the method setup or define a method in the project Settings."
	}
	##############################################################################
	if(any(ls()=="logfile")){stop("\n illegal logfile detected #2 in check_project.r!")}
	return(say);
	##############################################################################  
}
