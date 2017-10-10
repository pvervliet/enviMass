#' @title Reset workflow functionalities
#'
#' @export
#'
#' @description Reset downstream workflow functionalities
#'
#' @param down Name of node (logfile$Tasks_to_redo) to be altered
#' @param added Name of node (logfile$Tasks_to_redo) to be added - additionally
#' @param except Name of node (logfile$Tasks_to_redo) to be excluded - dangerous
#' @param check_node TRUE = resetting only evaluate when the concerned node is enabled in the workflow ...
#' @param down_TF Vector with node execution types: c("TRUE"), c("FALSE") or c("TRUE","FALSE") ... and it matches do_ or dont_
#' @param single_file TRUE = surpress resetting of all file-wise handlers to FALSE for all affected nodes?
#' 
#' @details enviMass workflow function. 
#' 

workflow_set<-function(
		down,
		added=FALSE,
		except=FALSE,
		check_node=FALSE,
		down_TF=c("TRUE","FALSE"),
		single_file=FALSE,
		...
	){

	########################################################################################	
	if(any(ls()=="logfile")){stop("\n illegal logfile detected #1 in workflow_set.r!")}
	if(any(is.na(match(down_TF,c("TRUE","FALSE"))))){stop("\n Wrong down_TF detectected - debug me!")}
	########################################################################################
	# check inputs #########################################################################
	if(!is.logical(added) & !is.logical(except)){
		if((any(!is.na(match(added,except)))) || (any(!is.na(match(except,added))))){
			stop("workflow_set: added or except nodes but not both for a node.")
		}
	}
	if(is.na(match(down,names(logfile$Tasks_to_redo)))){
		stop(paste("workflow_set: unknown down argument",down))
	}
	if(length(down)>1){
		stop("workflow_set: from which node down? Please specify only ONE node at a time.")
	}
	if(any(is.na(match(logfile$workflow,c("yes","no"))))){
		stop("workflow_set: missing yes/no in logfile$workflow found - please debug!")	
	}
	if(any(is.na(match(logfile$Tasks_to_redo,c("TRUE","FALSE"))))){
		stop("workflow_set: missing yes/no in logfile$workflow found - please debug!")	
	}	
	########################################################################################

	########################################################################################
	# only evaluate if affected script matches node exec (do_ vs. dont_) type - ############
	# otherwise evaluated once workflow exec settings change anyway! #######################
	# file handler - if node will later be included it must be re-exec. ###
	if(check_node){
		if(length(down_TF) == 1){
			if(
				((logfile$workflow[names(logfile$workflow)==down] == "yes") & (down_TF == "FALSE")) |
				((logfile$workflow[names(logfile$workflow)==down] == "no") & (down_TF == "TRUE"))
			){
# - REMOVE:
				#measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
				#if(any(names(measurements)==down)){
				#	measurements[,names(measurements)==down]<-FALSE;
				#	write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
				#	cat(" ** ")
				#}
				#rm(measurements)
# <
				return(NULL);
			}
		}	
	}
	########################################################################################
	depend <- logfile$workflow_depend
	diag(depend) <- 1
	must <- logfile$workflow_must
	########################################################################################
	# retrieve Tasks to redo ###############################################################
	# redefined dependencies:
	if(!is.logical(added)){
		for(i in 1:length(added)){
			depend[rownames(depend) == added[i],colnames(depend) == down] <- 1
		}
	}
	# direct dependencies:
	work_stream <- rownames(depend)[
		(depend[,colnames(depend)==down] == 1) |
		(depend[,colnames(depend)==down] == 2) |
		(depend[,colnames(depend)==down] == 3)
	]
	# collect indirect downstream dependencies
	doit <- TRUE
	while(doit){
		doit <- FALSE
		new_stream <- work_stream
		for(i in 1:length(work_stream)){
			new_nodes <- rownames(depend)[
				(depend[,colnames(depend) == work_stream[i]] == 1) |
				(depend[,colnames(depend) == work_stream[i]] == 2)			
			] # i.e., omit indirect downstream dependencies == 3
			new_nodes <- new_nodes[is.na(match(new_nodes, new_stream))]
			if(length(new_nodes) > 0){
				new_stream <- c(new_stream, new_nodes)
				doit <- TRUE
			}
		}
		work_stream <- unique(new_stream)
	}
	# explicitly excluded dependencies:
	if(!is.logical(except)){
		if(any(work_stream == except)){
			those <- match(work_stream,except)	
			work_stream <- work_stream[is.na(those)]
		}
	}
	########################################################################################
	
	########################################################################################
	# UPDATE $workflow #####################################################################
	# adapt chained requirements if set to "no" ############################################
	########################################################################################
	doit <- TRUE
	while(doit){ # if chained musts exist
		doit <- FALSE	
		for(i in 1:length(must[1,])){		# over columns
			if(logfile$workflow[names(logfile$workflow) == (colnames(must)[i])] == "yes"){
				for(j in 1:length(must[,1])){	# over rows
					if(must[j,i] > 0){ # upstream (1) & downstream (2) musts
						if(logfile$workflow[names(logfile$workflow) == (rownames(must)[j])] != "yes"){
							cat("\n", colnames(must)[i], " requires execution of ", rownames(must)[j])
							logfile$workflow[names(logfile$workflow) == (rownames(must)[j])] <<- "yes"
							doit <- TRUE
						}
					}
				}				
			}
		}
	}
	########################################################################################	
	
	########################################################################################
	# UPDATE $Tasks_to_redo ################################################################
	########################################################################################
	measurements <- read.csv(file = file.path(logfile[[1]], "dataframes", "measurements"), colClasses = "character");
	cat("\n"); cat("\nAdapt dependent nodes: ");
	while(length(work_stream) > 0){
		####################################################################################
# - REMOVE: if( (logfile$workflow[names(logfile$workflow)==work_stream[1]]=="yes") || (!check_node) ){
			cat(work_stream[1]); cat(" - ")
			logfile$Tasks_to_redo[names(logfile$Tasks_to_redo) == work_stream[1]] <<- TRUE;			
			if(	# reset single files
				(any(names(measurements) == work_stream[1])) & (!single_file)
			){
				measurements[,names(measurements) == work_stream[1]] <- FALSE;
				cat(" * ")
			}
# - REMOVE: }	
		####################################################################################
		work_stream<-work_stream[-1]
	}
	write.csv(measurements, file = file.path(logfile[[1]], "dataframes", "measurements"), row.names = FALSE);
	rm(measurements);
	########################################################################################
		
	########################################################################################
	save(logfile, file = file.path(as.character(logfile[[1]]), "logfile.emp"));
    logfile <<- logfile;
	if(any(ls() == "logfile")){stop("\n illegal logfile detected #2 in workflow_set.r!")}
	########################################################################################	
				
}
