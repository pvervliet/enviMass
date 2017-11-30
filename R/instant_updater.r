#' @title Update (very) old projects
#'
#' @description Update (very) old projects by rewriting instead of sequential updating
#'
#' @param project_path path to old project
#' @param project_path_new path to build new/updated project in
#' @param skip logical Check version number & skip if current
#' 
#' @details enviMass workflow function to rewrite an old into a new project while updating the version
#' 

#project_path <- "F:/PART_1/MS PROJECTS/RUeS/enviMass"
#project_path_new <- "F:/PART_1/MS PROJECTS/RUeS"

instant_updater <-
function(
	project_path,
	project_path_new,
	skip = TRUE
){

	new_version <- as.character(packageVersion("enviMass"))
	if(!file.exists(file.path(project_path, "logfile.emp"))) stop("\nInvalid project_path - abort.")
	say_path <- enviMass::check_path(pro_name = "updated_project", pro_dir = project_path_new)
	if(say_path!="Project path ok") stop("\n Problem with project_path_new: project already exists, the specified path is invalid or you lack permissions.")
	########################################################################################
	# make copies ##########################################################################
	load(file.path(project_path, "logfile.emp"), envir = as.environment(".GlobalEnv"))
	if(logfile$version == new_version & skip) return("\nProject build under latest version - update skipped")
	logfile$project_folder <- project_path
	save(logfile, file = file.path(as.character(logfile[[1]]),"logfile_updater_copy.emp"))
	if(file.exists(file.path(as.character(logfile[[1]]), "dataframes", "IS.txt"))){
		file.copy(
			from = file.path(as.character(logfile[[1]]), "dataframes", "IS.txt"),
			to = file.path(as.character(logfile[[1]]), "dataframes", "IS_updater_copy.txt")
		)
	}
	if(file.exists(file.path(as.character(logfile[[1]]), "dataframes", "targets.txt"))){
		file.copy(
			from = file.path(as.character(logfile[[1]]), "dataframes", "targets.txt"),
			to = file.path(as.character(logfile[[1]]), "dataframes", "targets_updater_copy.txt")
		)
	}
	if(file.exists(file.path(as.character(logfile[[1]]), "dataframes", "measurements.txt"))){
		file.copy(
			from = file.path(as.character(logfile[[1]]), "dataframes", "measurements.txt"),
			to = file.path(as.character(logfile[[1]]), "dataframes", "measurements_updater_copy.txt")
		)
	}
	########################################################################################
	# make a new project inside the old one ################################################
	if(!exists("IS", envir = as.environment(".GlobalEnv"))){data(IS, package = "enviMass")}
	if(!exists("targets", envir = as.environment(".GlobalEnv"))){data(targets, package = "enviMass")}
	workflow_depend <- try(read.table( # for library(enviMass)
		file = file.path(find.package("enviMass"), "webMass", "workflow_depend")		
	), silent = TRUE)
	if(class(workflow_depend) == "try-error"){ # devtools case
		workflow_depend <- try(read.table( # for library(enviMass)
			file = file.path(find.package("enviMass"), "inst", "webMass", "workflow_depend")		
		))
	}
	if(class(workflow_depend) == "try-error") stop("\n Cannot access workflow_depend")
	workflow_depend <- as.matrix(workflow_depend)
	workflow_must <- try(read.table( # for library(enviMass)
		file = file.path(find.package("enviMass"), "webMass", "workflow_must")		
	), silent = TRUE)
	if(class(workflow_must) == "try-error"){ # devtools case
		workflow_must <- try(read.table( # for library(enviMass)
			file = file.path(find.package("enviMass"), "inst", "webMass", "workflow_must")		
		))
	}
	if(class(workflow_must) == "try-error") stop("\n Cannot access workflow_must")
	workflow_must <- as.matrix(workflow_must)
	enviMass::newproject(
		pro_name = "updated_project", 
		pro_dir = project_path_new, 
		IS, targets,
		workflow_depend = workflow_depend,
		workflow_must = workflow_must,
		skript_check = FALSE
	)
	########################################################################################
	# update logfile & reset to peakpicking ################################################
	logfile_old <- logfile; rm(logfile)
	load(file.path(project_path_new, "updated_project", "logfile.emp"), envir = as.environment(".GlobalEnv"))
	# on parameters
	for(i in 1:length(logfile$parameters)){
		if(names(logfile$parameters)[i] == "") next # for any empty entries
		if(names(logfile$parameters)[i] == "external") next # for external parameters
		that <- which(names(logfile_old$parameters) == names(logfile$parameters)[i])
		if(!length(that)) next
		logfile$parameters[i] <- logfile_old$parameters[that]
	}
	# on workflow options - all nodes to be run by default (->  logfile$summary)
	for(i in 1:length(logfile$workflow)){	
		that <- which(names(logfile_old$workflow) == names(logfile$workflow)[i])
		if(!length(that)) next	
		logfile$workflow[i] <- logfile_old$workflow[that]
	}
	# subtraction files
	if(any(names(logfile_old) == "Positive_subtraction_files")) logfile$Positive_subtraction_files <- logfile_old$Positive_subtraction_files
	if(any(names(logfile_old) == "Negative_subtraction_files")) logfile$Negative_subtraction_files <- logfile_old$Negative_subtraction_files
	# on adducts
	if(any(names(logfile_old) == "adducts_pos_group")) logfile$adducts_pos_group <- logfile_old$adducts_pos_group
	if(any(names(logfile_old) == "adducts_neg_group")) logfile$adducts_neg_group <- logfile_old$adducts_neg_group
	if(any(names(logfile_old) == "adducts_pos")) logfile$adducts_pos <- logfile_old$adducts_pos	
	if(any(names(logfile_old) == "adducts_neg")) logfile$adducts_neg <- logfile_old$adducts_neg
	# PW path
	if(any(names(logfile_old) == "PW MSconvert path")){ 
		logfile[names(logfile) == "PW MSconvert path"] <- logfile_old[names(logfile_old) == "PW MSconvert path"]
	}
	save(logfile, file = file.path(as.character(logfile[["project_folder"]]), "logfile.emp"));	
	########################################################################################
	# update ISTD table ####################################################################
	if(file.exists(file.path(as.character(logfile_old[[1]]), "dataframes", "IS.txt"))){
		IS_old <- read.table(file=file.path(logfile_old[[1]],"dataframes","IS.txt"), header = TRUE, sep = "\t", colClasses = "character");
		len <- dim(IS_old)[1]	
		IS <- IS[1,, drop = FALSE]
		IS1 <- IS
		for(i in 2:len) IS <- rbind(IS, IS1) # know this is brute
		IS$ID <- IS_old$ID
		IS$Name <- IS_old$Name
		IS$Formula <- IS_old$Formula
		IS$RT <- IS_old$RT
		for(j in 5:dim(IS)[2]){
			if(!any(names(IS)[j] == names(IS_old))) next 
			that <- which(names(IS_old) == names(IS)[j])
			IS[,j] <- IS_old[,that]
		}
		write.table(IS, file = file.path(as.character(logfile[[1]]), "dataframes", "IS.txt"), row.names = FALSE, sep = "\t", quote = FALSE)      	  	
		rm(IS, IS_old)
	}
	########################################################################################	
	# update targets table #################################################################
	if(file.exists(file.path(as.character(logfile_old[[1]]), "dataframes", "targets.txt"))){
		targets_old <- read.table(file=file.path(logfile_old[[1]],"dataframes","targets.txt"), header = TRUE, sep = "\t", colClasses = "character");
		len <- dim(targets_old)[1]	
		targets <- targets[1,, drop = FALSE]
		targets1 <- targets
		for(i in 2:len) targets <- rbind(targets, targets1) # know this is brute
		targets$ID <- targets_old$ID
		targets$Name <- targets_old$Name
		targets$Formula <- targets_old$Formula
		targets$RT <- targets_old$RT
		for(j in 5:dim(targets)[2]){
			if(!any(names(targets)[j] == names(targets_old))) next 
			that <- which(names(targets_old) == names(targets)[j])
			targets[,j] <- targets_old[,that]
		}
		write.table(targets, file = file.path(as.character(logfile[[1]]), "dataframes", "targets.txt"), row.names = FALSE, sep = "\t", quote = FALSE)      	  	
		rm(targets, targets_old)
	}	
	########################################################################################
	# update measurements ##################################################################
	if(file.exists(file.path(as.character(logfile_old[[1]]), "dataframes", "measurements"))){	
		measurements <- read.csv(file = file.path(logfile[[1]], "dataframes", "measurements"), colClasses = "character");
		measurements[1,] <- "FALSE"
		measurements_old <- read.csv(file = file.path(logfile_old[[1]], "dataframes", "measurements"), colClasses = "character");
		# check whether old .mzXML files exist
		mzXML_files <- list.files(file.path(logfile_old[[1]], "files"))
		if(any(is.na(match(mzXML_files, paste0(measurements_old$ID, ".mzXML"))))) stop("\n Missing mzXML files - cannot proceed update - abort")
		len <- dim(measurements_old)[1]	
		measurements1 <- measurements[1,,drop = FALSE]
		for(i in 2:len) measurements <- rbind(measurements, measurements1) # know this is brute		
		measurements$ID <- measurements_old$ID
		measurements$Name <- measurements_old$Name
		measurements$Type <- measurements_old$Type
		measurements$Mode <- measurements_old$Mode
		measurements$Place <- measurements_old$Place
		measurements$Date <- measurements_old$Date
		measurements$Time <- measurements_old$Time
		measurements$include <- measurements_old$include
		measurements$copied <- measurements_old$copied
		measurements$tag1 <- measurements_old$tag1		
		measurements$tag2 <- measurements_old$tag2		
		measurements$tag3 <- measurements_old$tag3
		if(any(names(measurements_old) == "ID_2")) measurements$ID_2 <- measurements_old$ID_2
		if(any(names(measurements_old) == "date_end")) measurements$date_end <- measurements_old$date_end		
		if(any(names(measurements_old) == "time_end")) measurements$time_end <- measurements_old$time_end		
		# copy files #######################################################################
		for(i in 1:len){
			file.copy(
				from = file.path(logfile_old[[1]], "files", paste0(measurements[i,"ID"], ".mzXML")),
				to = file.path(logfile[[1]], "files", paste0(measurements[i,"ID"], ".mzXML"))
			)
		}
		write.csv(measurements, file = file.path(logfile[[1]],"dataframes","measurements"), row.names = FALSE);		
	}else{
		stop("\n Where is the measurements table of the project to be updated? - abort.")
	}
	######################################################################################## 
	rm(logfile, logfile_old)
	return("\n Instant update completed - open & rerun project to check if the update was successfull!")
	
}



