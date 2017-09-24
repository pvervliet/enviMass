	measurements <- read.csv(file = file.path(logfile[[1]], "dataframes", "measurements"), colClasses = "character");
	#measurements$qc <- "FALSE"
    for_IDs <- measurements[(measurements$include == "TRUE") & (measurements$qc == "FALSE") ,]$ID	
    not_for_IDs <- measurements[(measurements$include == "FALSE"),]$ID		
	
	if(file.exists(file.path(logfile[[1]], "results", "int_distrib"))){
		load(file.path(logfile[[1]], "results", "int_distrib"))
	}else{
		int_distrib <- list()
	}
	######################################################################################################
	
	######################################################################################################	
	if(length(not_for_IDs)){ # clean old entries
		for(i in not_for_IDs){
			at_entry <- as.numeric(i)
			int_distrib[[at_entry]] <- list()		
		}
	}	
	if(length(for_IDs)){
		iles<-seq(0, 1, .01)
		for(i in for_IDs){
		
			##############################################################################################
			load(file = file.path(logfile[[1]], "peaklist", i));   
			at_entry <- as.numeric(i)
			int_distrib[[at_entry]] <- list()
			int_distrib[[at_entry]]$ionization <- measurements[measurements$ID == i,"Mode"]
			int_distrib[[at_entry]]$peak_count <- dim(peaklist)[1]
			int_distrib[[at_entry]]$quantiles <- quantile(log10(peaklist[,"max_int"]), iles)
			int_distrib[[at_entry]]$boxplot <- boxplot(log10(peaklist[,"max_int"]), plot = FALSE)
			rm(peaklist)
			##############################################################################################
			
		}
	}
	save(int_distrib, file = file.path(logfile[[1]], "results", "int_distrib"))
	int_distrib <<- int_distrib
	######################################################################################################

	
	######################################################################################################	
	# plot positive ######################################################################################
	if( any(measurements[(measurements$include == "TRUE"),]$Mode == "positive") ){
		# boxplots #######################################################################################
		output$int_box_pos <- renderPlot({   
			par(mar = c(5.8, 4.5, .9, 8))
			enviMass:::plot_int_distrib(
				int_distrib,
				ion_mode = "positive",
				measurements = measurements,
				what = "boxplot"
			)
		},res = 100) 		
		# quantiles ######################################################################################
		output$int_quantiles_pos <- renderPlot({   
			par(mar = c(4.5, 4.5, .9, 8))
			enviMass:::plot_int_distrib(
				int_distrib,
				ion_mode = "positive",
				measurements = measurements,
				what = "quantiles_distrib"
			)
		},res = 100) 
		output$int_maxmed_pos <- renderPlot({   
			par(mar = c(4.5, 4.5, .9, 8))
			enviMass:::plot_int_distrib(
				int_distrib,
				ion_mode = "positive",
				measurements = measurements,
				what = "quantiles_out"
			)		
		},res = 100) 		
		##################################################################################################
	}else{
		output$int_box_pos <- renderPlot({ 
			plot.new(); plot.window(xlim = c(0, 10), ylim =c(0, 10)); text(5, 5, "Data not available")
		})					
		output$int_quantiles_pos <- renderPlot({ 
			plot.new(); plot.window(xlim = c(0, 10), ylim =c(0, 10)); text(5, 5, "Data not available")
		})					
		output$int_maxmed_pos <- renderPlot({ 
			plot.new(); plot.window(xlim = c(0, 10), ylim =c(0, 10)); text(5, 5, "Data not available")
		})							
	}
	######################################################################################################	


	######################################################################################################	
	# plot negative ######################################################################################
	if( any(measurements[(measurements$include == "TRUE"),]$Mode == "negative") ){
		# boxplots #######################################################################################
		output$int_box_neg <- renderPlot({   
			par(mar = c(5.8, 4.5, .9, 8))
			enviMass:::plot_int_distrib(
				int_distrib,
				ion_mode = "negative",
				measurements = measurements,
				what = "boxplot"
			)
		},res = 100) 		
		# quantiles ######################################################################################
		output$int_quantiles_neg <- renderPlot({   
			par(mar = c(4.5, 4.5, .9, 8))
			enviMass:::plot_int_distrib(
				int_distrib,
				ion_mode = "negative",
				measurements = measurements,
				what = "quantiles_distrib"
			)
		},res = 100) 
		output$int_maxmed_neg <- renderPlot({   
			par(mar = c(4.5, 4.5, .9, 8))
			enviMass:::plot_int_distrib(
				int_distrib,
				ion_mode = "negative",
				measurements = measurements,
				what = "quantiles_out"
			)		
		},res = 100) 		
		##################################################################################################
	}else{
		output$int_box_neg <- renderPlot({ 
			plot.new(); plot.window(xlim = c(0, 10), ylim =c(0, 10)); text(5, 5, "Data not available")
		})					
		output$int_quantiles_neg <- renderPlot({ 
			plot.new(); plot.window(xlim = c(0, 10), ylim =c(0, 10)); text(5, 5, "Data not available")
		})					
		output$int_maxmed_neg <- renderPlot({ 
			plot.new(); plot.window(xlim = c(0, 10), ylim =c(0, 10)); text(5, 5, "Data not available")
		})							
	}
	######################################################################################################	
	
	
	
	
	
	
	