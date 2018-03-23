			
			
			
			
if(length(logfile$comparisons)){
	at_comparisons <- which(names(logfile$comparisons) != "")
	if(length(at_comparisons)){
	
		if(any(objects(envir = as.environment(".GlobalEnv")) == "profileList_pos")){rm(profileList_pos, envir = as.environment(".GlobalEnv"))}
		if(any(objects() == "profileList_pos")){rm(profileList_pos)}	
		if(any(objects(envir = as.environment(".GlobalEnv")) == "profileList_neg")){rm(profileList_neg, envir = as.environment(".GlobalEnv"))}
		if(any(objects() == "profileList_neg")){rm(profileList_neg)}	
		measurements <- read.csv(file=file.path(logfile[[1]],"dataframes","measurements"), colClasses = "character");
		comparisons_pos <- c()		
		comparisons_neg <- c()
		for(i in at_comparisons){
			if(!logfile$comparisons[[i]][[1]]) next; # comparison defined, but not used in project
			# parse comparison, find file entries, get ionization mode
			parsed <- logfile$comparisons[[i]][[2]]
			parsed <- strsplit(parsed, "\"")[[1]]
			is_file <- rep(0, length(parsed))
			for(j in 1:length(is_file)){
				that_file <- which(parsed[j] == measurements$ID)
				if(!length(that_file) & is.na(as.numeric(parsed[j]))) next;
				if(!length(that_file) & !is.na(as.numeric(parsed[j]))) stop("\n NOTE: problem in comparison: invalid file ID detected") # invalid file IDs?
				is_file[j] <- that_file	
			}	
			for_mode <- unique(measurements$Mode[is_file])
			# check if specified files are all from the same ionization mode
			if(length(for_mode) > 1) stop(paste("\n NOTE: the comparison ", names(logfile$comparisons)[i], " involves positive and negative file IDs - please revise."))
			# POSITIVE IONIZATION comparisons #####################################################
			if(for_mode == "positive"){
				load(file.path(as.character(logfile[[1]]),"results","profileList_pos"), envir = as.environment(".GlobalEnv"))
				sampleID <- profileList_pos[["sampleID"]];
				# filter out other file types such as spiked ones
				timeset <- matrix(nrow = length(sampleID), ncol = 5,0);
				timeset[,2] <- as.numeric(sampleID)
				# check if specified files are all included	
				if(any(is.na(match(measurements$ID[is_file[is_file != 0]], sampleID)))){
					prob_file <- which(is.na(match(measurements$ID[is_file[is_file != 0]], sampleID)))
					stop(paste("\nNOTE: problem in comparison ", names(logfile$comparisons)[i], " for file IDs: ", measurements$ID[is_file[is_file != 0]][prob_file])," - not found in profiles; revise!")
				}
				# insert parsed timeset row entries to evaluate
				for(j in 1:length(is_file)){
					if(is_file[j] == 0) next;
					at_timeset <- which(timeset[,2] == as.numeric(measurements$ID[is_file[j]]))
					parsed[j] <- paste0( 
					"timeset[", at_timeset, ",4]")
				}	
				parsed_expression <- paste0(parsed, collapse = "")
				parsed_expression <- gsub("AND", "&&", parsed_expression)
				parsed_expression <- gsub("OR", "||", parsed_expression)				
				parsed_expression <- gsub("MUST", "0<", parsed_expression)	
				parsed_expression <- gsub("NOT", "0==", parsed_expression)	
				parsed_expression <- parse(text = parsed_expression)
				m <- 1
				n <- dim(profileList_pos[["index_prof"]])[1]
				store_comparison <- rep(0, n)
				for(k in m:n){
					timeset[,4:length(timeset[1,])] <- 0;
					timeset[,c(4,5)] <- .Call("_enviMass_fill_timeset",
											as.numeric(timeset),
											as.numeric(profileList_pos[["peaks"]][(profileList_pos[["index_prof"]][k,"start_ID"]:profileList_pos[[7]][k,"end_ID"]),"sampleIDs"]), 
											as.numeric(profileList_pos[["peaks"]][(profileList_pos[["index_prof"]][k,"start_ID"]:profileList_pos[[7]][k,"end_ID"]),"intensity"]), 
											as.integer(length(timeset[,1])),
											PACKAGE="enviMass"
										)	
					outcome <- eval(parsed_expression)
					if(is.logical(outcome)) outcome <- as.numeric(outcome)
					if(length(outcome) > 1) stop(paste("\nNOTE: problem in comparison ", names(logfile$comparisons)[i], ": produces several instead of one single value. Please revise!"))
					if(is.numeric(outcome)){
						store_comparison[k] <- outcome
					}else{
						stop("\nNOTE: problem in comparison - non-logical or non-numeric result observed; revise!")
					}
				}
				named <- colnames(profileList_pos[["index_prof"]])
				profileList_pos[["index_prof"]] <- cbind(profileList_pos[["index_prof"]], store_comparison)
				colnames(profileList_pos[["index_prof"]]) <- c(named, names(logfile$comparisons)[i])
				save(profileList_pos, file = file.path(as.character(logfile[[1]]),"results","profileList_pos"), compress = FALSE);
				comparisons_pos <- c(comparisons_pos, names(logfile$comparisons)[i])
				if(any(objects(envir = as.environment(".GlobalEnv")) == "profileList_pos")){rm(profileList_pos, envir = as.environment(".GlobalEnv"))}
				if(any(objects() == "profileList_pos")){rm(profileList_pos)}	
			}
			# NEGATIVE IONIZATION comparisons #####################################################
			if(for_mode == "negative"){
				load(file.path(as.character(logfile[[1]]),"results","profileList_neg"), envir = as.environment(".GlobalEnv"))
				sampleID <- profileList_neg[["sampleID"]];
				# filter out other file types such as spiked ones
				timeset <- matrix(nrow = length(sampleID), ncol = 5,0);
				timeset[,2] <- as.numeric(sampleID)
				# check if specified files are all included	
				if(any(is.na(match(measurements$ID[is_file[is_file != 0]], sampleID)))){
					prob_file <- which(is.na(match(measurements$ID[is_file[is_file != 0]], sampleID)))
					stop(paste("\nNOTE: problem in comparison ", names(logfile$comparisons)[i], " for file IDs: ", measurements$ID[is_file[is_file != 0]][prob_file])," - not found in profiles; revise!")
				}
				# insert parsed timeset row entries to evaluate
				for(j in 1:length(is_file)){
					if(is_file[j] == 0) next;
					at_timeset <- which(timeset[,2] == as.numeric(measurements$ID[is_file[j]]))
					parsed[j] <- paste0( 
					"timeset[", at_timeset, ",4]")
				}	
				parsed_expression <- paste0(parsed, collapse = "")
				parsed_expression <- gsub("AND", "&&", parsed_expression)
				parsed_expression <- gsub("OR", "||", parsed_expression)				
				parsed_expression <- gsub("MUST", "0<", parsed_expression)	
				parsed_expression <- gsub("NOT", "0==", parsed_expression)	
				parsed_expression <- parse(text = parsed_expression)
				m <- 1
				n <- dim(profileList_neg[["index_prof"]])[1]
				store_comparison <- rep(0, n)
				for(k in m:n){
					timeset[,4:length(timeset[1,])] <- 0;
					timeset[,c(4,5)] <- .Call("_enviMass_fill_timeset",
											as.numeric(timeset),
											as.numeric(profileList_neg[["peaks"]][(profileList_neg[["index_prof"]][k,"start_ID"]:profileList_neg[[7]][k,"end_ID"]),"sampleIDs"]), 
											as.numeric(profileList_neg[["peaks"]][(profileList_neg[["index_prof"]][k,"start_ID"]:profileList_neg[[7]][k,"end_ID"]),"intensity"]), 
											as.integer(length(timeset[,1])),
											PACKAGE="enviMass"
										)	
					outcome <- eval(parsed_expression)
					if(is.logical(outcome)) outcome <- as.numeric(outcome)
					if(length(outcome) > 1) stop(paste("\nNOTE: problem in comparison ", names(logfile$comparisons)[i], ": produces several instead of one single value. Please revise!"))
					if(is.numeric(outcome)){
						store_comparison[k] <- outcome
					}else{
						stop("\nNOTE: problem in comparison - non-logical or non-numeric result observed; revise!")
					}
				}
				named <- colnames(profileList_neg[["index_prof"]])
				profileList_neg[["index_prof"]] <- cbind(profileList_neg[["index_prof"]], store_comparison)
				colnames(profileList_neg[["index_prof"]]) <- c(named, names(logfile$comparisons)[i])
				save(profileList_neg, file = file.path(as.character(logfile[[1]]),"results","profileList_neg"), compress = FALSE);
				comparisons_neg <- c(comparisons_neg, names(logfile$comparisons)[i])
				if(any(objects(envir = as.environment(".GlobalEnv")) == "profileList_neg")){rm(profileList_neg, envir = as.environment(".GlobalEnv"))}
				if(any(objects() == "profileList_neg")){rm(pprofileList_neg)}	
			}
		}
		if( (isolate(input$Ion_mode) == "positive") & (length(comparisons_pos)) ){
			updateSelectInput(session, "filterProf_comparison", choices = c("None", comparisons_pos), selected = "None")
		}
		if( (isolate(input$Ion_mode) == "negative") & (length(comparisons_neg)) ){
			updateSelectInput(session, "filterProf_comparison", choices = c("None", comparisons_neg), selected = "None")
		}
	}
}














