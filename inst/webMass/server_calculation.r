mainchecked<-reactive({
    input$Check
    if(input$Check){
    	if( (isolate(input$ignore_large_files) == "TRUE") || (logfile$parameters$is_example == "TRUE")){ ignorefiles <- TRUE }else{ ignorefiles <- FALSE }
		say<-enviMass::check_project(isotopes, adducts, skipcheck = isolate(input$do_project_check), ignorefiles = ignorefiles, write_tables = FALSE);
		output$dowhat <<- renderText(say)
		if(say == "Project consistent"){
			shinytoastr::toastr_success(say, title = "Project check message:");
			cat("Project consistent\n");
			return("Project consistent\n");
		}else{
			cat("Project inconsistent\n");
			shinytoastr::toastr_error(say, title = "Project check message:", closeButton = TRUE, position = c("top-center"), timeOut = 0);
			#shinyjs::info(say);
			return("Project inconsistent\n");		
		}
	}
})
output$had_checked <- renderText(paste(mainchecked()))  

observe({ # Set calculation counter "do_flow"
    input$Calc
    if(input$Calc){
      do_flow <<- 0
	  time_start <<- Sys.time()
    }
})

#observe({ # Run calculations
maincalc <- reactive({
    input$Calc
    if(input$Calc){
		if(!any(objects(envir = as.environment(".GlobalEnv")) == "do_flow")){ # just in case the above has not reacted (happened) ...
			do_flow <<- 0
			time_start <<- Sys.time()	
		}
		if(do_flow == 0){	# check only once, initially at do_flow==0! really?
			enviMass::reset_selections(session)
    		if( (isolate(input$ignore_large_files) == "TRUE") || (logfile$parameters$is_example == "TRUE")){ ignorefiles <- TRUE }else{ ignorefiles <- FALSE }
			say <- enviMass::check_project(
				isotopes, adducts, 
				skipcheck = isolate(input$do_project_check), 
				ignorefiles = ignorefiles, 
				write_tables = FALSE); # because of write_tables = TRUE only here, this check must remain here!
			output$dowhat <<- renderText(say)
			enviMass::reset_selections(session)
		}else{
			say <- "Project consistent"
		}
      if(say == "Project consistent"){
	  
		if(any(ls() == "logfile")){stop("\n illegal logfile detected #1 in server_calculation.r!")}
        ########################################################################
        # restart logfile$summary & mark data availability #####################      
        if(do_flow == 0){

			####################################################################
			# set timer table ##################################################
			calTim <<- data.frame(
				logfile$summary[,1],
				rep(NA, length(logfile$summary[,1])),
				row.names=NULL
			)
			names(calTim) <<- c("Step","Time")
			save(calTim,file = file.path(as.character(logfile[[1]]), "exports", "calTim"))
			####################################################################
			if(logfile$parameters$parallel == "TRUE"){
				cat("\nStarting worker sessions ... ")
				if(logfile$parameters$parallel_restrict == "FALSE"){
					num_cores <- detectCores(all.tests = FALSE, logical = TRUE)
				}else{
					num_cores <- as.numeric(logfile$parameters$parallel_cores)
				}
				clus <<- makeCluster(getOption("cl.cores", num_cores))
				#	stopCluster(clus)
				clusterEvalQ(cl = clus,{library(enviMass, verbose = FALSE); NULL})
				cat("done.\n")				
			}
			####################################################################
			closeAlert(session, alertId = "a3")
			logfile$summary[1,2] <<- c(TRUE); # peakpicking
			save(logfile, file = file.path(as.character(logfile[[1]]), "logfile.emp"));
			summa <<- logfile$summary
			summa[,2] <<- "..."
			output$summa_html <- renderText(enviMass::summary_html(summa));
        }
		########################################################################
		# run calculations #####################################################
        if(do_flow > 0 & do_flow <= length(logfile$summary[,1])){	
			at_node <- logfile$summary[do_flow,1]
			if(logfile$parameters$parallel != "TRUE"){
				at_script_do <- paste("do_", at_node, ".R", sep = "")
				at_script_dont <- paste("dont_", at_node, ".R", sep = "")
				if(!file.exists(at_script_do)){
					at_script_do <- FALSE
				}	
				if(!file.exists(at_script_dont)){
					at_script_dont <- FALSE
				}			
			}else{
				at_script_do <- paste("do_", at_node, "_pl", ".R", sep = "") # use parallelized script ... 
				if(!file.exists(at_script_do)){ # ... if it exists
					at_script_do <- paste("do_", at_node, ".R", sep = "")
					if(!file.exists(at_script_do)){
						at_script_do <- FALSE
					}
				}		
				at_script_dont <- paste("dont_", at_node, "_pl", ".R", sep = "") # use parallelized script ... 
				if(!file.exists(at_script_dont)){ # ... if it exists
					at_script_dont <- paste("dont_", at_node, ".R", sep = "")
					if(!file.exists(at_script_dont)){
						at_script_dont <- FALSE
					}
				}						
			}
			time_start_part <- Sys.time()
			try_flow <- try({
				enviMass::workflow_node(
					at_node, at_node, at_node, at_node,
					path_do = at_script_do,
					path_undo = at_script_dont,
					session, output, input
				)  		
			})	
			time_diff_part <- as.numeric(difftime(Sys.time(), time_start_part, units = "mins"))
			calTim[do_flow,2] <<- round(time_diff_part,digits=4)
			write.table(calTim, file = file.path(as.character(logfile[[1]]), "exports", "calTim"))
			if(class(try_flow) == "try-error"){
				do_flow <<- 1000
				try_flow_message <- paste0("Workflow problem encountered at project node ", at_node, ". Revise settings or report the problem. Details: ", try_flow[1]);
				shinytoastr::toastr_warning(try_flow_message, title = "Project check message:", closeButton = TRUE, position = c("top-center"), timeOut = 0);
				if( (logfile$parameters$parallel == "TRUE") & (any(objects() == "clus")) ){
					cat("\n Closing worker sessions ... ")
					stopCluster(clus)
					if(any(objects() == "clus")){rm(clus)}
					if(any(objects(envir = as.environment(".GlobalEnv")) == "clus")){rm(clus, envir = as.environment(".GlobalEnv"))}
					cat("done. ")				
				}						
			}
		}
		########################################################################
        # make function reiterate ##############################################
        if(do_flow == (length(logfile$summary[,1]) + 1)){
			# set all Tasks_to_redo to FALSE ###################################
			logfile$Tasks_to_redo[1:length(logfile$Tasks_to_redo)] <<- "FALSE"	
			save(logfile, file = file.path(as.character(logfile[[1]]), "logfile.emp"));  
			# clean .GlobalEnv and reload results ############################## 
			enviMass:::workflow_objects(
				logfile,
				Ion_mode_profiles = isolate(input$Ion_mode)
			)			
			####################################################################
        }
        do_flow <<- (do_flow+1);
		if(do_flow == (length(logfile$summary[,1]) + 2)){
			output$summa_html <<- renderText(enviMass::summary_html(logfile$summary));
			measurements <- read.csv(file = file.path(logfile[[1]], "dataframes", "measurements"), colClasses = "character")
			measurements_tab <- measurements # because measurements[,"ID"] used further down!
			measurements_tab[,1] <- as.numeric(measurements_tab[,1])
			measurements_tab <- DT::datatable(
				measurements_tab[,c("ID","Name","Type","Mode","Place","Date","Time","include","profiled","tag1",
					"tag2","tag3","date_end","time_end","ID_2")],
				extensions = c('Buttons','ColReorder','FixedHeader'),
				rownames = FALSE,
				options = list(
					lengthMenu = list(c(25, 50, 100, -1), list('25', '50', '100', 'All')),
					fixedHeader = TRUE,
					ordering = T,
					dom = 'Blfrtip',
					buttons = c('excel', 'csv', 'colvis'),#buttons = c('excel', 'pdf', 'print', 'csv'),
					scrollX = TRUE,
					colReorder = TRUE
				)					
			)
			output$measurements <- DT::renderDataTable(measurements_tab); 		
		}
        if(do_flow < (length(logfile$summary[,1]) + 3)){
			invalidateLater(1, session=NULL)
			cat("Calculating...");
			return("Calculating...")
		}else{		
			if(do_flow < 1000){
				output$summa_html <<- renderText(enviMass::summary_html(logfile$summary));
				isolate(init$b <- (init$b + 1))
				if(any(ls() == "logfile")){stop("\n illegal logfile detected #2 in server_calculation.r!")}
				cat("Calculations completed, with a \n")
				time_diff <- (Sys.time() - time_start)
				print(time_diff)
				cat("\n")
				if(logfile$parameters$parallel == "TRUE"){
					cat("\n Closing worker sessions ... ")
					stopCluster(clus)
					if(any(objects() == "clus")){rm(clus)}
					if(any(objects(envir = as.environment(".GlobalEnv")) == "clus")){rm(clus, envir = as.environment(".GlobalEnv"))}
					cat("done. ")				
				}								
				return("\n Calculations completed \n")
			}else{
				cat("\n")
				return(try_flow_message);
				cat(try_flow_message);				
			}	
		}
        ########################################################################
      }else{
        cat("Project inconsistent\n");
		#shinyjs::info(say);
		shinytoastr::toastr_error(say, title = "Project check message:", closeButton = TRUE, position = c("top-center"), timeOut = 0);
		return("Project inconsistent\n");
      }
    }
})
output$had_calculated <- renderText(paste(maincalc()))  
