##############################################################################
# ON PROFILES ################################################################
##############################################################################

##############################################################################
# PROFILE FITERING (1) - update results for changes in ion mode selection ####
##############################################################################
maincalc3<-reactive({

	input$Ion_mode
	##########################################################################
	if( (isolate(init$a) == "TRUE") & (isolate(input$Ion_mode) == "positive") ){
		exprprofnorm_pos <- list(src=file.path(logfile[[1]],"pics","profnorm_pos"))
		#output$profnorm<-renderImage(exprprofnorm_pos, deleteFile = FALSE)
		exprprofcount_pos <- list(src = file.path(logfile[[1]], "pics", "profcount_pos"))
		output$profcount <- renderImage(exprprofcount_pos, deleteFile = FALSE)
		if(any(objects(envir = as.environment(".GlobalEnv")) == "profileList")){ rm(profileList, inherits = TRUE) }
		if(file.exists(file.path(logfile[[1]], "results", "profileList_pos"))){
			if(any(objects(envir = as.environment(".GlobalEnv")) == "profileList_pos")){rm(profileList_pos, envir = as.environment(".GlobalEnv"))}
			if(any(objects() == "profileList_pos")){rm(profileList_pos)}				
			load(file = file.path(as.character(logfile[[1]]),"results","profileList_pos"), envir = as.environment(".GlobalEnv"), verbose=TRUE);
			assign("profileList", profileList_pos, envir = as.environment(".GlobalEnv"));
		}else{
			return("No profiles available for this ionization mode")
		}
		if(file.exists(file.path(logfile[[1]], "results", "links_profiles_pos"))){
			if(any(objects(envir = as.environment(".GlobalEnv")) == "links_profiles_pos")){rm(links_profiles_pos,envir=as.environment(".GlobalEnv"))}
			if(any(objects() == "links_profiles_pos")){rm(links_profiles_pos)}				
			load(file = file.path(as.character(logfile[[1]]),"results","links_profiles_pos"), envir = as.environment(".GlobalEnv"), verbose=TRUE);
			assign("links_profiles", links_profiles_pos, envir = as.environment(".GlobalEnv"));
		}				
		expr4p <- list(src = file.path(logfile[[1]],"pics","boxprofile_pos"))
		output$boxprofile <- renderImage(expr4p, deleteFile = FALSE)		
		isolate(init$b <<- (init$b+1))
		# load compounds as well #############################################
		if(any(objects(envir = as.environment(".GlobalEnv")) == "compounds_targets")){ rm(compounds_targets, inherits = TRUE) }
		if(any(objects(envir = as.environment(".GlobalEnv")) == "compounds_IS")){ rm(compounds_IS, inherits = TRUE) }
		compounds_IS <<- read.table(file=file.path(logfile[[1]],"dataframes","IS.txt"),header=TRUE,sep="\t",colClasses = "character")
		compounds_IS <<- compounds_IS[compounds_IS$ion_mode == "positive",,drop = FALSE]
		if(!dim(compounds_IS)[1]) rm(compounds_IS, envir = as.environment(".GlobalEnv"))
		compounds_targets <<- read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character")
		compounds_targets <<- compounds_targets[compounds_targets$ion_mode == "positive",,drop = FALSE]
		if(!dim(compounds_targets)[1]) rm(compounds_targets, envir = as.environment(".GlobalEnv"))
		######################################################################
		if(any(objects() == "profileList")){stop("illegal profileList found, #1");}
		if(any(objects() == "profpeaks2")){stop("illegal profpeaks2 found, #1");}
		if(any(objects() == "profpeaks3")){stop("illegal profpeaks3 found, #1");}
		if(any(objects() == "compounds_targets")){stop("illegal compounds_targets found, #1");}
		if(any(objects() == "compounds_IS")){stop("illegal compounds_IS found, #1");}
		return("Select ionization (switch to negative):\n")
	}
	##########################################################################
	if( (isolate(init$a) == "TRUE") &  (isolate(input$Ion_mode) == "negative") ){
		exprprofnorm_neg <- list(src=file.path(logfile[[1]],"pics","profnorm_neg"))
		output$profnorm <- renderImage(exprprofnorm_neg, deleteFile = FALSE)
		exprprofcount_neg <- list(src=file.path(logfile[[1]],"pics","profcount_neg"))
		output$profcount <- renderImage(exprprofcount_neg, deleteFile = FALSE)
		if(any(objects(envir = as.environment(".GlobalEnv")) == "profileList")){ rm(profileList, inherits = TRUE) }
		if(file.exists(file.path(logfile[[1]], "results", "profileList_neg"))){
			if(any(objects(envir = as.environment(".GlobalEnv")) == "profileList_neg")){rm(profileList_neg, envir = as.environment(".GlobalEnv"))}
			if(any(objects() == "profileList_neg")){rm(profileList_neg)}				
			load(file = file.path(as.character(logfile[[1]]),"results","profileList_neg"), envir = as.environment(".GlobalEnv"), verbose=TRUE);
			assign("profileList", profileList_neg, envir = as.environment(".GlobalEnv"));
		}else{
			return("No profiles available for this ionization mode")			
		}
		if(file.exists(file.path(logfile[[1]],"results","links_profiles_neg"))){
			if(any(objects(envir=as.environment(".GlobalEnv"))=="links_profiles_neg")){rm(links_profiles_neg,envir=as.environment(".GlobalEnv"))}
			if(any(objects()=="links_profiles_neg")){rm(links_profiles_neg)}				
			load(file=file.path(as.character(logfile[[1]]),"results","links_profiles_neg"),envir=as.environment(".GlobalEnv"), verbose=TRUE);
			assign("links_profiles",links_profiles_neg,envir=as.environment(".GlobalEnv"));
		}	
		expr4n<-list(src=file.path(logfile[[1]],"pics","boxprofile_neg"))
		output$boxprofile<-renderImage(expr4n, deleteFile = FALSE)	
		isolate(init$b<<-(init$b+1))
		# load compounds as well #############################################
		if(any(objects(envir = as.environment(".GlobalEnv")) == "compounds_targets")){ rm(compounds_targets, inherits = TRUE) }
		if(any(objects(envir = as.environment(".GlobalEnv")) == "compounds_IS")){ rm(compounds_IS, inherits = TRUE) }
		compounds_IS <<- read.table(file=file.path(logfile[[1]],"dataframes","IS.txt"),header=TRUE,sep="\t",colClasses = "character")
		compounds_IS <<- compounds_IS[compounds_IS$ion_mode == "negative",,drop = FALSE]
		if(!dim(compounds_IS)[1]) rm(compounds_IS, envir = as.environment(".GlobalEnv"))
		compounds_targets <<- read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character")
		compounds_targets <<- compounds_targets[compounds_targets$ion_mode == "negative",,drop = FALSE]
		if(!dim(compounds_targets)[1]) rm(compounds_targets, envir = as.environment(".GlobalEnv"))
		######################################################################
		if(any(objects() == "profileList")){stop("illegal profileList found, #1");}
		if(any(objects() == "profpeaks2")){stop("illegal profpeaks2 found, #1");}
		if(any(objects() == "profpeaks3")){stop("illegal profpeaks3 found, #1");}
		if(any(objects() == "compounds_targets")){stop("illegal compounds_targets found, #1");}
		if(any(objects() == "compounds_IS")){stop("illegal compounds_IS found, #1");}
		return("Select ionization (switch to positive):\n")	
	}
	##########################################################################
})
output$had_ion<-renderText(paste(maincalc3())) 
##############################################################################

##############################################################################
# PROFILE FITERING (2) - Sort and filter the profile list ####################
##############################################################################
maincalc6 <- reactive({

	init$a # in project?
	init$b # number of calculations - update profpeaks2 after each such
	#cat(" \n ... ");print(isolate(init$b));cat(" ... ")
	input$filterProf_maxmass
	input$filterProf_minmass
	input$filterProf_minrt
	input$filterProf_maxrt
	input$filterProf_meanblind
	input$filterProf_notblind
	input$filterProf_sort
	# input$filterProf_count <- REMOVE!
	input$filterProf_medianblind
	input$filterProf_medianblind_value
	input$filterProf_minMD
	input$filterProf_maxMD
	input$filterProf_components
	input$search_profile_compound
	cat("\n profileList filtered and sorted_1")
    if( 
		(isolate(init$a)=="TRUE") & 
		(any(objects(envir=as.environment(".GlobalEnv"))=="profileList")) & 
		!is.na(isolate(input$filterProf_minmass)) & 
		!is.na(isolate(input$filterProf_maxmass)) & 
		!is.na(isolate(input$filterProf_minrt)) & 
		!is.na(isolate(input$filterProf_maxrt)) &
		!is.na(isolate(input$filterProf_minMD)) &		
		!is.na(isolate(input$filterProf_maxMD))		
	){

		cat("\n profileList filtered and sorted_2")		
		if(any(objects(envir = as.environment(".GlobalEnv")) == "profpeaks2")){rm(profpeaks2,envir=as.environment(".GlobalEnv"))}
		if(any(objects() == "profpeaks2")){rm(profpeaks2)}		
		assign("profpeaks2", profileList[["index_prof"]], envir = as.environment(".GlobalEnv"));
		if(any(objects() == "profileList")){stop("illegal profileList found, #3");}
		if(any(objects() == "profpeaks2")){stop("illegal profpeaks2 found, #3");}
		if(any(objects() == "profpeaks3")){stop("illegal profpeaks3 found, #3");}		
		###################################################################################################		
		profpeaks2 <<- profpeaks2[profpeaks2[,"mean_mz"] >= isolate(input$filterProf_minmass),,drop = FALSE]
		profpeaks2 <<- profpeaks2[profpeaks2[,"mean_mz"] <= isolate(input$filterProf_maxmass),,drop = FALSE]
		profpeaks2 <<- profpeaks2[profpeaks2[,"mean_RT"] >= isolate(input$filterProf_minrt),,drop = FALSE]
		profpeaks2 <<- profpeaks2[profpeaks2[,"mean_RT"] <= isolate(input$filterProf_maxrt),,drop = FALSE]
		profpeaks2 <<- profpeaks2[profpeaks2[,"Mass defect"] >= isolate(input$filterProf_minMD),,drop = FALSE]
		profpeaks2 <<- profpeaks2[profpeaks2[,"Mass defect"] <= isolate(input$filterProf_maxMD),,drop = FALSE]		
		###################################################################################################		
		if(isolate(input$filterProf_medianblind) == "yes"){
			profpeaks2 <<- profpeaks2[(profpeaks2[,"above_blind?"] >= as.numeric(isolate(input$filterProf_medianblind_value))),, drop = FALSE] 
		}
		if(isolate(input$filterProf_notblind) == "yes"){
			profpeaks2 <<- profpeaks2[profpeaks2[,"in_blind?"]==0,, drop = FALSE] # not in blind, = profileList[[7]][k,8]
		}
		###################################################################################################		
		# compound searches ###############################################################################
		if( 
			isolate(input$search_profile_compound != "Type in compound name/ID or erase any entry to clear filter") &
			isolate(input$search_profile_compound != "") &
			((any(objects(envir=as.environment(".GlobalEnv")) == "compounds_targets")) | (any(objects(envir=as.environment(".GlobalEnv")) == "compounds_IS")))
		){
			cat("\n Searching compound profiles")
			if(
				(logfile$workflow[names(logfile$workflow)=="components_profiles"]=="yes") &
				(logfile$summary[logfile$summary[,1]=="components_profiles",2]=="TRUE") &
				(any(objects(envir=as.environment(".GlobalEnv"))=="links_profiles")) &
				(dim(profpeaks2)[1] > 0)
			){
			
				keep_IDs <<- rep(FALSE, dim(profpeaks2)[1])
				# find compounds by search string
				# on ISTDs ################################################################################				
				if(any(objects(envir = as.environment(".GlobalEnv")) == "compounds_IS")){
					search_compound_ID <<- c()
					this <- which(compounds_IS$ID == isolate(input$search_profile_compound))
					if(length(this)) search_compound_ID <<- c(search_compound_ID, compounds_IS$ID[this])
					this <- which(compounds_IS$Name == isolate(input$search_profile_compound))
					if(length(this)) search_compound_ID <<- c(search_compound_ID, compounds_IS$ID[this])
					len <- nchar(search_compound_ID)
					if(length(search_compound_ID)){
						for(m in 1:dim(profpeaks2)[1]){
							if(profpeaks2[m, "links"] == 0) next
							if(length(links_profiles[[profpeaks2[m, "links"]]]$IS) == 0) next
							for(mm in 1:dim(links_profiles[[profpeaks2[m, "links"]]]$IS)[1]){
								if(substr(links_profiles[[profpeaks2[m, "links"]]]$IS[mm, "Compound"], 1, len) == search_compound_ID) keep_IDs[m] <<- TRUE
							}
						}
					}
				}
				# on targets ##############################################################################
				if(any(objects(envir = as.environment(".GlobalEnv")) == "compounds_targets")){
					search_compound_ID <<- c()
					this <- which(compounds_targets$ID == isolate(input$search_profile_compound))
					if(length(this)) search_compound_ID <<- c(search_compound_ID, compounds_targets$ID[this])
					this <- which(compounds_targets$Name == isolate(input$search_profile_compound))
					if(length(this)) search_compound_ID <<- c(search_compound_ID, compounds_targets$ID[this])
					len <- nchar(search_compound_ID)
					if(length(search_compound_ID)){
						for(m in 1:dim(profpeaks2)[1]){
							if(profpeaks2[m, "links"] == 0) next
							if(length(links_profiles[[profpeaks2[m, "links"]]]$targ) == 0) next
							for(mm in 1:dim(links_profiles[[profpeaks2[m, "links"]]]$targ)[1]){
								if(substr(links_profiles[[profpeaks2[m, "links"]]]$targ[mm, "Compound"], 1, len) == search_compound_ID) keep_IDs[m] <<- TRUE
							}
						}
					}
				}
				###########################################################################################
				if(all(!keep_IDs)){
					shinytoastr::toastr_warning("No profile match for this compound found - compound filter omitted!", 
						title = "Profile filtering:", closeButton = TRUE, position = c("top-center"), timeOut = 0);					
				}else{
					profpeaks2 <<- profpeaks2[keep_IDs,,drop = FALSE]
				}
			}else{
				shinytoastr::toastr_warning("You need to enable the workflow profile componentization to filter for specific compound profiles. Disable the compound search filter or enable the mentioned workflow step and recalculate!", 
					title = "Profile filtering:", closeButton = TRUE, position = c("top-center"), timeOut = 0);
			}
		}
		###################################################################################################	
		if(isolate(input$filterProf_sort) == "ID (increasing)"){
			profpeaks2 <<- profpeaks2[order(profpeaks2[,"profile_ID"], decreasing = FALSE),,drop = FALSE]
			sort_by <<- "profile_ID"; sort_by_decreasing <<- "FALSE"
		}
		if(isolate(input$filterProf_sort) == "mean m/z (increasing)"){
			profpeaks2 <<- profpeaks2[order(profpeaks2[,"mean_mz"], decreasing = FALSE),,drop = FALSE]
			sort_by <<- "mean_mz"; sort_by_decreasing <<- "FALSE"
		}
		if(isolate(input$filterProf_sort) == "mean m/z (decreasing)"){
			profpeaks2 <<- profpeaks2[order(profpeaks2[,"mean_mz"], decreasing = TRUE),,drop = FALSE]
			sort_by <<- "mean_mz"; sort_by_decreasing <<- "TRUE"
		}		
		if(isolate(input$filterProf_sort) == "mean RT (increasing)"){
			profpeaks2<<-profpeaks2[order(profpeaks2[,"mean_RT"], decreasing = FALSE),,drop = FALSE]
			sort_by<<-"mean_RT"; sort_by_decreasing <<- "FALSE"
		}	
		if(isolate(input$filterProf_sort) == "mean RT (decreasing)"){
			profpeaks2 <<- profpeaks2[order(profpeaks2[,"mean_RT"], decreasing = TRUE),,drop = FALSE]
			sort_by <<- "mean_RT"; sort_by_decreasing <<- "TRUE"
		}
		if(isolate(input$filterProf_sort) == "minimum RT (decreasing)"){		
			profpeaks2 <<- profpeaks2[order(profpeaks2[,"min_RT"], decreasing = TRUE),,drop = FALSE]	
			sort_by <<- "min_RT"; sort_by_decreasing <<- "TRUE"			
		}	
		if(isolate(input$filterProf_sort) == "maximum RT (decreasing)"){		
			profpeaks2 <<- profpeaks2[order(profpeaks2[,"max_RT" ], decreasing = TRUE),,drop = FALSE]	
			sort_by <<- "max_RT"; sort_by_decreasing <<- "TRUE"			
		}
		if(isolate(input$filterProf_sort) == "maximum overall intensity (decreasing)"){
			profpeaks2<<-profpeaks2[order(profpeaks2[,"max_int"], decreasing = TRUE),,drop = FALSE]
			sort_by <<- "max_int"; sort_by_decreasing <<- "TRUE"
		}
		if(isolate(input$filterProf_sort) == "maximum intensity in samples (decreasing, zeros removed)"){		
			profpeaks2 <<- profpeaks2[order(profpeaks2[,"max_int_sample"],decreasing=TRUE),,drop = FALSE]	
			profpeaks2 <<- profpeaks2[profpeaks2[,"max_int_sample"] != 0,,drop = FALSE]			
			sort_by <<- "max_int_sample"; sort_by_decreasing <<- "TRUE"
		}
		if(isolate(input$filterProf_sort) == "maximum intensity in blanks/blinds (decreasing, zeros removed)"){		
			profpeaks2 <<- profpeaks2[order(profpeaks2[,"max_int_blind"], decreasing = TRUE),, drop = FALSE]	
			profpeaks2 <<- profpeaks2[profpeaks2[,"max_int_blind"] != 0,, drop = FALSE]			
			sort_by <<- "max_int_blind"; sort_by_decreasing <<- "TRUE"
		}
		if(isolate(input$filterProf_sort) == "mean intensity (decreasing)"){
			profpeaks2 <<- profpeaks2[order(profpeaks2[,"mean_int"], decreasing = TRUE),, drop = FALSE]
			sort_by <<- "mean_int"; sort_by_decreasing <<- "TRUE"
		}
		if(isolate(input$filterProf_sort) == "mean intensity in samples (decreasing, zeros removed)"){		
			profpeaks2 <<- profpeaks2[order(profpeaks2[,"mean_int_sample"], decreasing = TRUE),,drop = FALSE]	
			profpeaks2 <<- profpeaks2[profpeaks2[,"mean_int_sample"] != 0,,drop = FALSE]
			sort_by <<- "mean_int_sample"; sort_by_decreasing <<- "TRUE"
		}
		if(isolate(input$filterProf_sort) == "mean intensity in blanks/blinds (decreasing, zeros removed)"){		
			profpeaks2 <<- profpeaks2[order(profpeaks2[,"mean_int_blind"],decreasing = TRUE),, drop = FALSE]	
			profpeaks2 <<- profpeaks2[profpeaks2[,"mean_int_blind"] != 0,,drop = FALSE]
			sort_by <<- "mean_int_blind"; sort_by_decreasing <<- "TRUE"
		}
		if(isolate(input$filterProf_sort) == "past & current trend intensity (decreasing)"){		
			profpeaks2 <<- profpeaks2[order(profpeaks2[,"deltaint_global"],decreasing = TRUE),, drop = FALSE]
			profpeaks2 <<- profpeaks2[profpeaks2[,"deltaint_global"] != 0,, drop = FALSE]
			sort_by <<- "deltaint_global"; sort_by_decreasing <<- "TRUE"
		}
		if(isolate(input$filterProf_sort) == "current trend intensity (decreasing)"){
			profpeaks2 <<- profpeaks2[order(profpeaks2[,"deltaint_newest"], decreasing = TRUE),, drop = FALSE]
			profpeaks2 <<- profpeaks2[profpeaks2[,"deltaint_newest"] != 0,,drop = FALSE]
			sort_by <<- "deltaint_newest"; sort_by_decreasing <<- "TRUE"
		}
		if(isolate(input$filterProf_sort)=="total peak number (decreasing)"){
			profpeaks2 <<- profpeaks2[order(profpeaks2[,"number_peaks_total"], profpeaks2[,"mean_mz"], decreasing = TRUE),, drop = FALSE]
			sort_by <<- "number_peaks_total";sort_by_decreasing <<- "TRUE"
		}			
		if(isolate(input$filterProf_sort) == "peak number in samples (decreasing, zeros removed)"){		
			profpeaks2 <<- profpeaks2[order(profpeaks2[,"number_peaks_sample"], profpeaks2[,"mean_mz"], decreasing = TRUE),, drop = FALSE]	
			profpeaks2 <<- profpeaks2[profpeaks2[,"number_peaks_sample"]!=0,,drop = FALSE]		
			sort_by <<-"number_peaks_sample";sort_by_decreasing <<- "TRUE"
		}
		if(isolate(input$filterProf_sort) == "peak number in blanks/blinds (decreasing, zeros removed)"){		
			profpeaks2 <<- profpeaks2[order(profpeaks2[,"number_peaks_blind"], profpeaks2[,"mean_mz"], decreasing = TRUE),,drop = FALSE]	
			profpeaks2 <<- profpeaks2[profpeaks2[,"number_peaks_blind"]!=0,,drop = FALSE]	
			sort_by <<- "number_peaks_blind"; sort_by_decreasing <<- "TRUE"
		}
		if(isolate(input$filterProf_sort) == "mass defect (increasing)"){		
			profpeaks2 <<- profpeaks2[order(profpeaks2[,"Mass defect"], decreasing = FALSE),,drop = FALSE]	
			sort_by <<- "Mass defect"; sort_by_decreasing <<- "FALSE"				
		}
		if(isolate(input$filterProf_sort) == "mass defect (decreasing)"){		
			profpeaks2 <<-profpeaks2[order(profpeaks2[,"Mass defect"], decreasing = TRUE),,drop = FALSE]		
			sort_by <<-"Mass defect"; sort_by_decreasing <<- "TRUE"
		}		
		if(isolate(input$filterProf_sort) == "median sample above blind intensity (decreasing)"){		
			profpeaks2 <<- profpeaks2[order(profpeaks2[,"in_blind?"], decreasing = TRUE),,drop = FALSE]		
			sort_by <<- "in_blind?"; sort_by_decreasing <<- "TRUE"
		}
		###################################################################################################		
		if(
			(isolate(input$filterProf_components)=="TRUE") &
			(any(objects(envir=as.environment(".GlobalEnv"))=="profpeaks2")) 			
		){
			if(
				(logfile$workflow[names(logfile$workflow) == "components_profiles"] == "yes") &
				(logfile$summary[logfile$summary[,1] == "components_profiles", 2] == "TRUE") &
				(any(objects(envir=as.environment(".GlobalEnv")) == "links_profiles")) 
			){
				keep_IDs <<- enviMass::analyseE_links_profiles(
					profileList_index = profpeaks2, 
					links_profiles, 
					sort_what = sort_by, 
					sort_decreasing = sort_by_decreasing,
					use_profile = NULL, 
					with_bar = FALSE, 
					return_excl = FALSE		
				)
				if(length(keep_IDs) & all(!is.na(keep_IDs))){
					profpeaks2 <<- profpeaks2[match(keep_IDs, profpeaks2[,"profile_ID"]),, drop = FALSE]
				}else{
					stop("\nDebug server_obs_res_meas.r @ 1")
				}
			}else{
				shinytoastr::toastr_warning("You need to enable the workflow profile componentization to filter 'lower-ranked profiles with redundant intensity patterns'. Disable this filter or enable the mentioned workflow step and recalculate!", 
					title = "Profile filtering:", closeButton = TRUE, position = c("top-center"), timeOut = 0);
			}
		}
		###################################################################################################	


		
		
		
		
		
		###################################################################################################			
		if(dim(profpeaks2)[1] > 0){
			atit1 <- sum(profpeaks2[,"number_peaks_total"]) 
			output$atprof1<-renderText({ atit1 })
			atit2 <- dim(profpeaks2)[1]
			#output$atprof2<-renderText({ atit2 })
			atit3 <- sum(profpeaks2[,"number_peaks_blind"] > 0)
			output$atprof3 <- renderText({ atit3 })
			atit4 <- sum(profpeaks2[,"deltaint_global"] > 0)
			output$atprof4 <- renderText({ atit4 })
			atit5 <- sum(profpeaks2[,"deltaint_newest"] > 0)
			output$atprof5 <- renderText({ atit5 })
			# intensity histogram
			path = file.path(logfile[[1]], "pics", "profilehisto.png");
                png(filename = path, bg = "white", width = 600);
                plot_profiles_intensity_histograms(	mean_intensities=profpeaks2[,"mean_int"],
                                                    max_intensities=profpeaks2[,"max_int"],
                                                    past_incidents=profpeaks2[,"deltaint_global"],
                                                    current_incidents=profpeaks2[,"deltaint_newest"]);
            dev.off();
			expr6 <- list(src=file.path(logfile[[1]],"pics","profilehisto.png"));
			output$profilehisto <- renderImage(expr6, deleteFile = FALSE);
			# generate output table ############################################
			#if( (length(profpeaks2[,1]) > isolate(input$filterProf_count))  &  !is.na(isolate(input$filterProf_count)) ){
			#	profpeaks2 <<- profpeaks2[1:isolate(input$filterProf_count),,drop=FALSE]
			#}
			profpeaks3 <<- as.data.frame(profpeaks2)
			profpeaks3[,"mean_mz" ] <<- round(profpeaks3[,"mean_mz"],digits=6)
			profpeaks3[,"mean_RT" ] <<- round(profpeaks3[,"mean_RT"],digits=1)
			profpeaks3[,"min_RT" ] <<- round(profpeaks3[,"min_RT"],digits=1)
			profpeaks3[,"max_RT" ] <<- round(profpeaks3[,"max_RT"],digits=1)
			profpeaks3[,"mean_int"] <<- round(log10(profpeaks3[,"mean_int"]),digits=3)
			profpeaks3[,"max_int"] <<- round(log10(profpeaks3[,"max_int"]),digits=3)
			profpeaks3[,"in_blind?"] <<- as.integer(profpeaks3[,"in_blind?"])
			profpeaks3[,"above_blind?"] <<- round(profpeaks3[,"above_blind?"],digits=4)
			profpeaks3[,"var_mz"] <<- round(log10(profpeaks3[,"var_mz"]),digits=3)
			profpeaks3[,"profile_ID"] <<- as.character(profpeaks3[,"profile_ID"])
			profpeaks3[,"number_peaks_total"] <<- as.integer(profpeaks3[,"number_peaks_total"])
			profpeaks3[,"deltaint_global"] <<- round(log10(profpeaks3[,"deltaint_global"]),digits=2)
			profpeaks3[,"deltaint_newest"] <<- round(log10(profpeaks3[,"deltaint_newest"]),digits=2)
			profpeaks3 <<- profpeaks3[,c("profile_ID","number_peaks_total","deltaint_global","deltaint_newest","mean_mz",
				"mean_int","mean_RT","mean_RT","max_int","in_blind?","above_blind?","var_mz","min_RT","max_RT")]
			profpeaks3[,8] <<- round((profpeaks3[,8] / 60), digits = 2)
			output$allproftable <- DT::renderDataTable(
				DT::datatable(profpeaks3,
					colnames=c("profile ID","number of peaks","log10 global trend intensity","log10 current trend intensity","mean m/z","log10 mean intensity", 
						"mean RT [s]","mean RT [min]","log10 maximum Intensity","in blind?","above blind?","log10 m/z variance","minimum RT [s]","maximum RT [s]"),
					rownames = FALSE,
					filter = 'top',
                    selection = list(mode = 'single', target = 'row'),
					extensions = c('Buttons','FixedHeader','ColReorder'),
					options = list(
						lengthMenu = list(c(25, 50, 100, 200, -1), list('25', '50', '100', '200', 'All')),
						fixedHeader = FALSE,
						ordering = T,
						dom = 'Blfrtip',
						buttons = c('excel', 'csv','colvis'),#buttons = c('excel', 'pdf', 'print', 'csv'),
						scrollX = TRUE#,
						#colReorder = TRUE
					)
				),
				server = TRUE
			)
			updateNumericInput(session, "profID", value = 0);
			updateNumericInput(session, "profentry", value = 0);
			return(as.character(atit2));
		}else{
#profpeaks3 <<- as.data.frame(profpeaks2)
			shinytoastr::toastr_warning("No profiles remaining with these filter settings", 
				title = "Profile filtering:", closeButton = TRUE, position = c("top-center"), timeOut = 0);
			output$atprof1 <- renderText({ "0" })
			#output$atprof2<-renderText({ "0" })
			output$atprof3 <- renderText({ "0" })
			output$atprof4 <- renderText({ "0" })
			output$atprof5 <- renderText({ "0" })
			path = file.path(logfile[[1]],"pics","profilehisto.png");
			png(filename = path, bg = "white", width = 1100);
				plot.new()
				plot.window(xlim = c(0,1), ylim = c(0,1))
				text(0.5, 0.5, labels = "0 profiles for these filter settings", cex = 1.8, col = "red")
			dev.off()
			expr6 <- list(src = file.path(logfile[[1]],"pics","profilehisto.png"))
			output$profilehisto <- renderImage(expr6, deleteFile = FALSE);
			output$allproftable <- renderText("No profiles left")
			updateNumericInput(session,"profID",value = 0)
			updateNumericInput(session,"profentry",value = 0)
			return("0")
		}
	}else{
		if( isolate(init$a)=="TRUE" ){
			cat("\n No profiles available\n")
			output$atprof1 <- renderText({"0"}) # now used as reactive output
			#output$atprof2<-renderText({"0"})
			output$atprof3 <- renderText({"0"})
			output$atprof4 <- renderText({"0"})
			output$atprof5 <- renderText({"0"})	
			output$allproftable <- renderText("No profiles available")
			path = file.path(logfile[[1]],"pics","profilehisto.png")
			png(filename = path, bg = "white", width = 1100)
				plot.new()
				plot.window(xlim = c(0,1), ylim = c(0,1))
				text(0.5, 0.5, labels = "0 profiles for these filter settings", cex = 1.8, col = "red")
			dev.off()
			expr6 <- list(src=file.path(logfile[[1]],"pics","profilehisto.png"))
			output$profilehisto <- renderImage(expr6, deleteFile = FALSE)
			updateNumericInput(session, "profID", value = 0)
			updateNumericInput(session, "profentry", value = 0)
			return("0")
		}
	}
	if(any(objects() == "profileList")){stop("illegal profileList found, #4");}
	if(any(objects() == "profpeaks2")){stop("illegal profpeaks2 found, #4");}
	if(any(objects() == "profpeaks3")){stop("illegal profpeaks3 found, #4");}
})	
output$atprof2 <- renderText(paste(maincalc6()))
output$prof_number <- renderText(paste(maincalc6())) 
##############################################################################

##############################################################################
# update results per profilepeak list entry index ############################
##############################################################################
observe({
	input$profentry
	init$b
	if(	
		(isolate(init$a)=="TRUE") &  
		!is.na(isolate(input$profentry)) & 
		(isolate(input$profentry) != 0) & 
		(any(objects(envir=as.environment(".GlobalEnv")) == "profileList")) & 
		(any(objects(envir=as.environment(".GlobalEnv")) == "profpeaks2"))
	){

		if( 
			(as.numeric(isolate(input$profentry)) <= length(profpeaks2[,1])) & 
			(as.numeric(isolate(input$profentry)) > 0)
		){
			if(any( profileList[["index_prof"]][,"profile_ID"] == as.numeric(profpeaks2[isolate(input$profentry), "profile_ID"]))){
				updateNumericInput(session, "profID", value = as.numeric(as.character(profpeaks2[isolate(input$profentry), "profile_ID"])))
			}else{
				updateNumericInput(session, "profID", value = 0)
				output$timeprofile <- renderPlot({	
					plot.new()
					plot.window(xlim = c(0, 1), ylim = c(0, 1))
					text(0.5, 0.5, labels="Invalid list entry", cex = 1.8, col = "red")
				})			
				output$oneproftable<-DT::renderDataTable(cbind("No date available",""));
			}
		}else{
			updateNumericInput(session, "profID", value = 0)
			output$timeprofile <- renderPlot({	
				plot.new()
				plot.window(xlim = c(0, 1), ylim = c(0, 1))
				text(0.5, 0.5, labels="Invalid list entry", cex = 1.8, col = "red")
			})			
			output$oneproftable<-DT::renderDataTable(cbind("No date available",""));
		}
	}else{
		if(isolate(init$a) == "TRUE"){
			updateNumericInput(session, "profID", value = 0)
			output$timeprofile <- renderPlot({	
				plot.new()
				plot.window(xlim = c(0, 1), ylim = c(0, 1))
				text(0.5, 0.5, labels = "Nothing to plot - invalid ID", cex = 1.8, col = "red")
			})			
			output$oneproftable<-DT::renderDataTable(cbind("No date available",""));
		}
	}
})	
##############################################################################

##############################################################################
# update results for individual profileIDs ###################################
##############################################################################
ranges_timeprofile <- reactiveValues(x = NULL, y = NULL)
observe({
	input$profID
	init$b
	cat("\n plotting profile _1")
	if( 
		(isolate(init$a) == "TRUE") &  
		(!is.na(isolate(input$profID))) & 
		(isolate(input$profID) != 0) & 
		any(objects(envir = as.environment(".GlobalEnv")) == "profileList") & 
		any(objects(envir = as.environment(".GlobalEnv")) == "profpeaks2") 	
	){
		cat("\n plotting profile _2")
		if(any(objects() == "profileList")){stop("illegal profileList found, #5");}
		if(any(profileList[["index_prof"]][,"profile_ID"] == as.numeric(isolate(input$profID)))){
			
			cat("\n plotting profile with ID "); cat(as.numeric(isolate(input$profID)))
			prof_plot_ID1 <<- as.numeric(isolate(input$profID))
			at_entry <<- profileList[["index_prof"]][
				match(prof_plot_ID1,  profileList[["index_prof"]][,"profile_ID"])
			,"links"]
			##################################################################
			if(logfile$parameters$trend_blind == "yes"){
				blindsubtract <- TRUE
			}else{
				blindsubtract <- FALSE
			}
			lagit <- as.numeric(strsplit(logfile$parameters$trend_lags, ",")[[1]])		
			if(isolate(input$prof_log) == "yes"){
				logscaled <- TRUE
			}else{
				logscaled <- FALSE
			}
			output$timeprofile <- renderPlot({			
				plotaprofile(
					profileList,
					profileID = prof_plot_ID1,
					logint = logscaled,
					blindsub = blindsubtract,
					blindfold = as.numeric(logfile$parameters$blind_threshold),
					lags = lagit,
					threshold = as.numeric(logfile$parameters$trend_threshold),
					ranges_x = ranges_timeprofile$x,
					ranges_y = ranges_timeprofile$y,
					return_data = FALSE
				)
				if((!is.null(ranges_timeprofile$x)) || (!is.null(ranges_timeprofile$y))){
					mtext("Now zoomed in", side = 3, col = "gray")
				}
			})			
			##################################################################
			# target / ISTD matches? #########################################
			if(
				(at_entry != 0) &
				any(objects(envir = as.environment(".GlobalEnv")) == "links_profiles")
			){
				# on targets #################################################
				if(
					(length(links_profiles[[at_entry]]$targ) !=  0) &
					(any(objects(envir = as.environment(".GlobalEnv")) == "compounds_targets")) 
				){
					targ_out <- c()
					for(k in 1:dim(links_profiles[[at_entry]]$targ)[1]){
						at_tar <- strsplit(links_profiles[[at_entry]]$targ[k, 1], "_")[[1]][1]
						at_add <- strsplit(links_profiles[[at_entry]]$targ[k, 1], "_")[[1]][2]
						at_tar <- compounds_targets[match(at_tar, compounds_targets$ID), "Name"]
						targ_out <- c(targ_out, 
							paste0(at_tar, " (", at_add , ", peak counts: ",  links_profiles[[at_entry]]$targ[k, 2], ", max score: ", round(links_profiles[[at_entry]]$targ[k, 3], digits = 2),")")
						)		
					}
					targ_out <- paste(targ_out, collapse = ",")
					targ_out <- paste("Target matches: ", targ_out, sep = "")
					output$prof_targets <- renderText(targ_out)
				}else{
					output$prof_targets <- renderText("Target matches: none") 
				}
				# on ISTDs ###################################################
				if(
					(length(links_profiles[[at_entry]]$IS) !=  0) &
					(any(objects(envir = as.environment(".GlobalEnv")) == "compounds_IS"))
				){
					ISTD_out <- c()
					for(k in 1:dim(links_profiles[[at_entry]]$IS)[1]){
						at_ISTD <- strsplit(links_profiles[[at_entry]]$IS[k, 1], "_")[[1]][1]
						at_add <- strsplit(links_profiles[[at_entry]]$IS[k, 1], "_")[[1]][2]
						at_ISTD <- compounds_IS[match(at_ISTD, compounds_IS$ID), "Name"]
						ISTD_out <- c(ISTD_out, 
							paste0(at_ISTD, " (", at_add , ", peak counts: ",  links_profiles[[at_entry]]$IS[k, 2], ", max score: ", round(links_profiles[[at_entry]]$IS[k, 3], digits = 2),")")
						)		
					}
					ISTD_out <- paste(ISTD_out, collapse = ",")
					ISTD_out <- paste("ISTD matches:", ISTD_out, sep = "")
					output$prof_ISTD  <- renderText(ISTD_out)
				}else{
					output$prof_ISTD  <- renderText("ISTD matches: none") 
				}
				##############################################################
			}else{
				output$prof_targets <- renderText("Target matches: none") 			
				output$prof_ISTD <- renderText("ISTD matches: none") 			
			}
			##################################################################
			# plot similar profiles? #########################################
			plot_similar_profiles <- FALSE
			if(
				(at_entry != 0) &
				any(objects(envir=as.environment(".GlobalEnv")) == "links_profiles")
			){
				if(length(links_profiles[[at_entry]]$group) > 0){
					prof_plot_IDs <<- c(prof_plot_ID1, links_profiles[[at_entry]]$group)
					prof_plot_IDs <<- unique(prof_plot_IDs)
					plot_similar_profiles <- TRUE
				}else{
					prof_plot_IDs <<- prof_plot_ID1				
				}
			}else{
				prof_plot_IDs <<- prof_plot_ID1
			}
			output$similar_profiles_plot <- renderPlot({
				if(plot_similar_profiles){
						enviMass::plot_components(
							profileList = profileList,
							prof_IDs = prof_plot_IDs,
							links_profiles = links_profiles,
							what = "profiles",
							xlim = FALSE,  ylim = FALSE,  await_input = FALSE,
							skipit = TRUE,
							min_peaks = NULL,
							norma = TRUE
						)	
				}else{
					if(logfile$workflow[names(logfile$workflow) == "components_profiles"] == "yes"){
						plot.new(); plot.window(xlim = c(0, 1), ylim = c(0, 1))
						text(0.5, 0.5, labels = "Nothing to plot - no other similar profiles available", cex = 1.8, col = "red")
					}else{	
						plot.new(); plot.window(xlim = c(0, 1), ylim = c(0, 1))
						text(0.5, 0.5, labels = "No profile (cross-file) componentization included in the workflow", cex = 1.8, col = "red")						
					}
				}
			})
			#output$similar_profiles_relations <- renderPlot({	
			#	if(plot_similar_profiles){
			#			enviMass::plot_components(
			#				profileList=profileList,
			#				prof_IDs=prof_plot_IDs,
			#				links_profiles=links_profiles,
			#				what="relations",
			#				xlim=FALSE,ylim=FALSE,await_input=FALSE,
			#				skipit=TRUE,
			#				min_peaks=NULL,
			#				norma=TRUE
			#			)	
			#	}else{
			#		plot.new()
			#		plot.window(xlim=c(0,1),ylim=c(0,1))
			#		text(0.5,0.5,labels="Nothing to plot - no other similar profiles available",cex=1.8,col="red")
			#	}
			#})
			use_prof_IDs <- match(prof_plot_IDs, profileList[["index_prof"]][,"profile_ID"])
			use_prof_IDs <- use_prof_IDs[!is.na(use_prof_IDs)]
			similar_profiles_tab <- profileList[["index_prof"]][use_prof_IDs,c("profile_ID","mean_mz","mean_int","max_int_sample","mean_RT",
				"Mass defect","above_blind?","number_peaks_sample","number_peaks_blind","number_peaks_total"),drop=FALSE]
			similar_profiles_tab <- as.data.frame(similar_profiles_tab)
			similar_profiles_tab[,"mean_mz" ] <- format(similar_profiles_tab[,"mean_mz" ],digits=5)
			similar_profiles_tab[,"mean_int"] <- format(similar_profiles_tab[,"mean_int"],scientific=TRUE,digits=2)
			similar_profiles_tab[,"max_int_sample"] <- format(similar_profiles_tab[,"max_int_sample"],scientific=TRUE,digits=2)
			similar_profiles_tab[,"mean_RT"] <- format(similar_profiles_tab[,"mean_RT"],digits=1)
			similar_profiles_tab[,"Mass defect"] <- format(similar_profiles_tab[,"Mass defect"],digits=3)
			similar_profiles_tab[,"above_blind?"] <- format(similar_profiles_tab[,"above_blind?"],digits=3)
			similar_profiles_tab[,"number_peaks_sample"] <- as.integer(similar_profiles_tab[,"number_peaks_sample"])
			similar_profiles_tab[,"number_peaks_blind"] <- as.integer(similar_profiles_tab[,"number_peaks_blind"])
			similar_profiles_tab[,"number_peaks_total"] <- as.integer(similar_profiles_tab[,"number_peaks_total"])
			output$similar_profiles_table <- DT::renderDataTable({
				DT::datatable(similar_profiles_tab,
					colnames = c("Profile ID","Mean m/z","Mean int. overall","Max int. samples","Mean RT",
						"Mass defect","Median int. ratio blind vs. samples","Peak number samples","Peak number blinds","Peak number total"),
					rownames = FALSE,
					filter = 'top',
					selection = list(mode = 'single', target = 'row'),
					extensions = c('Buttons','FixedHeader','ColReorder'),
					options = list(
						lengthMenu = c(15, 30, 50, 100),
						fixedHeader = FALSE,
						ordering = T,
						dom = 'Blfrtip',
						buttons = c('excel', 'colvis'),#buttons = c('excel', 'pdf', 'print', 'csv'),
						scrollX = TRUE,
						colReorder = TRUE
					)						
				)			
			});
			##################################################################
			peakTable <<- plotaprofile( # peakTable will also be required further below for EIC extraction
					profileList,
					profileID = prof_plot_ID1,
					plotit = FALSE,
					return_data = TRUE
				)
			peakTable <- cbind(seq(1, dim(peakTable)[1], 1), peakTable)
			output$oneproftable <- DT::renderDataTable({
				DT::datatable(
					peakTable,
					colnames = c("", "Date", "Time", "ID sample", "Intensity sample", "ID blind", "Intensity blind", "ID peak"),
					extensions = c('Buttons'),
					rownames = FALSE,					
					options = list(
						lengthMenu = c(50, 150, 300),
						ordering = T,
						dom = 'Blfrtip',
						buttons = c('excel', 'colvis')
					)					
				)
			});
			updateNumericInput(session,"profpeakID",value = 0);
			path=file.path(logfile[[1]], "pics", "massdens.png");
			png(filename = path, bg = "white", width = 550, height = 200);			
				plot.new()
				plot.window(xlim = c(0, 1), ylim = c(0, 1))
				text(0.5, 0.5, labels = "Waiting...", cex = 1.8, col = "red")
			dev.off();
			expr_massdens<-list(src=file.path(logfile[[1]],"pics","massdens.png"));
			output$massdens<-renderImage(expr_massdens, deleteFile = FALSE);				
			path=file.path(logfile[[1]],"pics","mass_int.png");
			png(filename = path, bg = "white", width = 550,height=300);	
				plot.new()
			dev.off();
			expr_mass_int<-list(src=file.path(logfile[[1]],"pics","mass_int.png"));
			output$massint<-renderImage(expr_mass_int, deleteFile = FALSE);	
		}else{
			output$timeprofile <- renderPlot({	
				plot.new()
				plot.window(xlim=c(0,1),ylim=c(0,1))
				text(0.5,0.5,labels="Nothing to plot - invalid ID",cex=1.8,col="red")
			})
			output$oneproftable<-DT::renderDataTable(cbind("No date available",""));
			updateNumericInput(session,"profpeakID",value = 0);		
			path=file.path(logfile[[1]],"pics","massdens.png");
			png(filename = path, bg = "white", width = 550,height=200);			
				plot.new()
				plot.window(xlim=c(0,1),ylim=c(0,1))
				text(0.5,0.5,labels="Not available",cex=1.8,col="red")
			dev.off();
			expr_massdens<-list(src=file.path(logfile[[1]],"pics","massdens.png"));
			output$massdens<-renderImage(expr_massdens, deleteFile = FALSE);				
			path=file.path(logfile[[1]],"pics","mass_int.png");
			png(filename = path, bg = "white", width = 550,height=300);	
				plot.new()
			dev.off();
			expr_mass_int<-list(src=file.path(logfile[[1]],"pics","mass_int.png"));
			output$massint<-renderImage(expr_mass_int, deleteFile = FALSE);	
		}
	}else{
		if(isolate(init$a)=="TRUE"){
			output$timeprofile <- renderPlot({	
				plot.new()
				plot.window(xlim=c(0,1),ylim=c(0,1))
				text(0.5,0.5,labels="Nothing to plot - invalid ID",cex=1.8,col="red")
			})			
			output$oneproftable<-DT::renderDataTable(cbind("No date available",""));
			updateNumericInput(session,"profpeakID",value = 0);
			path=file.path(logfile[[1]],"pics","massdens.png");
			png(filename = path, bg = "white", width = 550,height=200);			
				plot.new()
				plot.window(xlim=c(0,1),ylim=c(0,1))
				text(0.5,0.5,labels="Not available",cex=1.8,col="red")
			dev.off();
			expr_massdens<-list(src=file.path(logfile[[1]],"pics","massdens.png"));
			output$massdens<-renderImage(expr_massdens, deleteFile = FALSE);			
			path=file.path(logfile[[1]],"pics","mass_int.png");
			png(filename = path, bg = "white", width = 550,height=300);	
				plot.new()
			dev.off();
			expr_mass_int<-list(src=file.path(logfile[[1]],"pics","mass_int.png"));
			output$massint<-renderImage(expr_mass_int, deleteFile = FALSE);	
		}
	}
})	
##############################################################################

##############################################################################
# PLOT ZOOM & TABLE SELECTION ################################################
# When a double-click happens, check if there's a brush on the plot.
# If so, zoom to the brush bounds; if not, reset the zoom.
observeEvent(input$timeprofile_dblclick, { # - N
	if(verbose){cat("\n in N")}
    brush <- input$timeprofile_brush
    if (!is.null(brush)) {
		ranges_timeprofile$x <- c(brush$xmin, brush$xmax)
		ranges_timeprofile$y <- c(brush$ymin, brush$ymax)
    } else {
		ranges_timeprofile$x <- NULL
		ranges_timeprofile$y <- NULL
    }
})
##############################################################################

##############################################################################
# get EICs for individual profiles ###########################################
##############################################################################
ranges_profiles <- reactiveValues(RTchrom = FALSE,  intchrom = FALSE)

observe({ # seconds <-> minutes switch when zoomed ###########################
	input$profile_EIC_time
	if(isolate(init$a) == "TRUE" & isolate(ranges_profiles$RTchrom[1] != FALSE)){
		if(isolate(input$profile_EIC_time) == "minutes"){
			isolate(ranges_profiles$RTchrom <- (ranges_profiles$RTchrom / 60))
		}
		if(isolate(input$profile_EIC_time) == "seconds"){
			isolate(ranges_profiles$RTchrom <- (ranges_profiles$RTchrom * 60))
		}
	}
})
observeEvent(input$profile_EIC_dblclick, { 
	brush <- isolate(input$profile_EIC_brush)
    if (!is.null(brush)) {
        cat("\n Zoom in_1")
        isolate(ranges_profiles$RTchrom <- c(brush$xmin, brush$xmax))
        isolate(ranges_profiles$intchrom <- c(brush$ymin, brush$ymax))
    } else {
        cat("\n Zoom out full_1")
        isolate(ranges_profiles$RTchrom_ <- FALSE)
        isolate(ranges_profiles$intchrom <- FALSE)
    }
})
observeEvent(input$profile_EIC_click, { # NOTE: brushing already triggers a click -> use brush with delay=0, which embeds the slower click
    cat("\n Zoom out part_1_a")
    brush <- isolate(input$profile_EIC_brush)
    if (is.null(brush)) {
        cat("\n Zoom out part_1_b")
        if(isolate(ranges_profiles$intchrom[1])!=FALSE){
            old_range_mass<-abs(isolate(ranges_profiles$intchrom[2]-ranges_profiles$intchrom[1]))
            isolate(ranges_profiles$intchrom[1]<-ranges_profiles$intchrom[1]-.3*old_range_mass)
            isolate(ranges_profiles$intchrom[2]<-ranges_profiles$intchrom[2]+.3*old_range_mass)
        }
        if(isolate(ranges_profiles$RTchrom[1])!=FALSE){
		    old_range_RT<-abs(isolate(ranges_profiles$RTchrom[2]-ranges_profiles$RTchrom[1]))
		    isolate(ranges_profiles$RTchrom[1]<-ranges_profiles$RTchrom[1]-.1*old_range_RT)
		    isolate(ranges_profiles$RTchrom[2]<-ranges_profiles$RTchrom[2]+.1*old_range_RT)
	    }
    }else{
        cat("\n Doing hover - and nothing")
    }   
})
maincalc4<-reactive({
	input$profpeakID
	if( 	
		(isolate(init$a)=="TRUE") & 
		(!is.na(isolate(input$profID))) & 
		(isolate(input$profpeakID)>0) & 
		(any(objects(envir=as.environment(".GlobalEnv"))=="peakTable")) &
		(any(objects(envir=as.environment(".GlobalEnv"))=="profileList")) 
	){
		if( (isolate(input$profpeakID)<=length(peakTable[,1])) & (any(profileList[[7]][,4]==as.numeric(isolate(input$profID)))) ){
			# positioning plot ###############################################
			timed<-as.POSIXct(paste(peakTable[,1],peakTable[,2],sep=" "))
			output$profile_position <- renderPlot({
				par_old<-par(mar=c(1,.3,.3,.3))				
				plot.new()
				plot.window(xlim=c(min(timed),max(timed)),ylim=c(0,max(max(as.numeric(as.character(peakTable[,4]))),max(as.numeric(as.character(peakTable[,6]))))))
				abline(v=timed[as.numeric(isolate(input$profpeakID))],col="darkgrey",lwd=5)
				box();
				points(timed[peakTable[,3]!=0],peakTable[peakTable[,3]!=0,4],type="l",col="darkgreen");
				points(timed[peakTable[,5]!=0],peakTable[peakTable[,5]!=0,6],type="l",col="red");
				par(par_old);
			})			
			# EIC plot ########################################################
			if(any(objects(envir=as.environment(".GlobalEnv"))=="MSlist")){rm(MSlist,envir=as.environment(".GlobalEnv"))}
			if(any(objects()=="MSlist")){rm(MSlist)}				
			if(any(objects(envir=as.environment(".GlobalEnv"))=="EIC_ID")){rm(EIC_ID,envir=as.environment(".GlobalEnv"))}
			if(any(objects()=="EIC_ID")){rm(EIC_ID)}				
			fileID<-peakTable[,3]
			fileID[fileID=="0"]<-peakTable[peakTable[,3]=="0",5]
			if(any(peakTable[as.numeric(isolate(input$profpeakID)),7]!=0)){
				use_file_ID<-fileID[as.numeric(isolate(input$profpeakID))]
	        	if(any(objects(envir=as.environment(".GlobalEnv"))=="MSlist")){
	        		if(any(names(MSlist)=="File_ID")){
	        			if(MSlist[["File_ID"]]!=as.character(use_file_ID)){ # File_ID does not match
							load(file.path(logfile[[1]],"MSlist",as.character(use_file_ID)),envir=as.environment(".GlobalEnv")) 
	        			}
	        		}else{ # available MSlist not with File_ID yet
						load(file.path(logfile[[1]],"MSlist",as.character(use_file_ID)),envir=as.environment(".GlobalEnv"))  
						MSlist[["File_ID"]]<-as.character(use_file_ID)
	        		}
	        	}else{ # no MSlist in GlobalEnv
					load(file.path(logfile[[1]],"MSlist",as.character(use_file_ID)),envir=as.environment(".GlobalEnv"))  
	        	}			
				output$profile_EIC  <- renderPlot({
		            par(mar=c(4.5,4,.8,.8))
		            enviMass:::plotchromat(
		              	MSlist,
		              	peakIDs=as.numeric(peakTable[as.numeric(isolate(input$profpeakID)),7]),
		              	RTlim=ranges_profiles$RTchrom,
		              	Intlim=ranges_profiles$intchrom,
		              	set_RT=input$profile_EIC_time,
		              	normalize=FALSE,
		              	chromat_full=input$profile_EIC_type
		            );
		        },res=100)			
				return(
					paste("= sample file ID: ",as.character(fileID[as.numeric(isolate(input$profpeakID))])," (",as.character(timed[as.numeric(isolate(input$profpeakID))]),")" )
				);		
				#}					
			}else{
			output$profile_EIC <- renderPlot({
				plot.new()
				plot.window(xlim=c(0,1),ylim=c(0,1))
				text(0.5,0.5,labels="No sample peak available",cex=1.8,col="red")
			})				
			return("No peak or EIC available");
			}
		}else{
			output$profile_position <- renderPlot({
				plot.new()
				plot.window(xlim=c(0,1),ylim=c(0,1))
				text(0.5,0.5,labels="Selection out of range",cex=1.8,col="red")
			})			
			output$profile_EIC <- renderPlot({
				plot.new()
				plot.window(xlim=c(0,1),ylim=c(0,1))
				text(0.5,0.5,labels="No sample peak available",cex=1.8,col="red")
			})				
			return("No peak or EIC available");
		}
	}else{
		if( (isolate(init$a)=="TRUE") ){
			output$profile_position <- renderPlot({
				plot.new()
				plot.window(xlim=c(0,1),ylim=c(0,1))
				text(0.5,0.5,labels="Selection out of range",cex=1.8,col="red")
			})			
			output$profile_EIC <- renderPlot({
				plot.new()
				plot.window(xlim=c(0,1),ylim=c(0,1))
				text(0.5,0.5,labels="No sample peak available",cex=1.8,col="red")
			})				
			return("No peak or EIC available");
		}
	}
})
output$prof_peak_text<-renderText(paste(maincalc4())) 
##############################################################################
	
##############################################################################
# get mass estimates for individual profiles #################################
##############################################################################
maincalc5 <- reactive({
	input$dens_mass
	if(	(isolate(input$dens_mass)) &
		(isolate(init$a) == "TRUE") & 
		(isolate(input$profID) != 0) & 
		(any(objects(envir=as.environment(".GlobalEnv")) == "profileList")) 
	){
		if(any(profileList[[7]][,4] == as.numeric(isolate(input$profID)))){
			######################################################################
			cat("\n kernel density ...")
			if(isolate(input$use_weight)=="yes"){use_weights<-TRUE}else{use_weights<-FALSE}
			path=file.path(logfile[[1]],"pics","massdens.png");
			png(filename = path, bg = "white", width = 550,height=400);	
				getmass <- mass_dens(
						profileList,
						profileID = as.numeric(isolate(input$profID)),
						bootstrap = TRUE,
						boot_size = as.numeric(isolate(input$boot_size)),
						use_weights
				)
			dev.off();
			expr_massdens <- list(src=file.path(logfile[[1]],"pics","massdens.png"));
			output$massdens <- renderImage(expr_massdens, deleteFile = FALSE);	
			path = file.path(logfile[[1]],"pics","mass_int.png");
			png(filename = path, bg = "white", width = 550, height = 300);	
					mass_int(
						profileList,
						profileID = as.numeric(isolate(input$profID))
					)
			dev.off();
			expr_mass_int <- list(src = file.path(logfile[[1]],"pics","mass_int.png"));
			output$massint <- renderImage(expr_mass_int, deleteFile = FALSE);	
			return(getmass);
			######################################################################
		}else{
			return("Not available");
		}
	}else{
		######################################################################
		cat("\n kernel density ...")
		path = file.path(logfile[[1]],"pics","massdens.png");
		png(filename = path, bg = "white", width = 550, height = 200);			
			plot.new()
			plot.window(xlim = c(0, 1), ylim = c(0,1))
			text(0.5, 0.5, labels = "Not available", cex = 1.8, col = "red")
		dev.off();
		expr_massdens <- list(src = file.path(logfile[[1]], "pics", "massdens.png"));
		output$massdens <- renderImage(expr_massdens, deleteFile = FALSE);		
		return("...");
		######################################################################
	}
})
output$prof_mass <- renderText(paste(maincalc5())) 
##############################################################################

  
  
  
  
  
  
