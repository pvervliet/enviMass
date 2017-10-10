ranges_peaks_mz_RT <- reactiveValues(x = NULL, y = NULL, xchroma=FALSE, ychroma=FALSE) # used for several of the below observers
refresh_plot<-reactiveValues()
refresh_plot$a<-1
refresh_plot$b<-0
refresh_plot$c<-0




##############################################################################
# update results for individual measurements #################################
##############################################################################
observe({
    input$sel_meas
	input$blind_boxplot_log
	if(isolate(init$a)=="TRUE"){
	if(!is.na(isolate(input$sel_meas))){
    if(isolate(input$sel_meas)!=0){
		##########################################################################	
		measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
		if(any(measurements[,"ID"]==as.character(isolate(input$sel_meas)))){
			output$file_proc_name<-renderText(paste("File name: ",measurements[measurements[,"ID"]==as.character(isolate(input$sel_meas)),"Name"],sep=""))
			output$file_proc_type<-renderText(paste("File type: ",measurements[measurements[,"ID"]==as.character(isolate(input$sel_meas)),"Type"],sep=""))
			output$file_proc_mode<-renderText(paste("Ionization mode: ",measurements[measurements[,"ID"]==as.character(isolate(input$sel_meas)),"Mode"],sep=""))
			# peaklist info ##########################################################
			if(	file.exists(file.path(logfile$project_folder,"peaklist",as.character(isolate(input$sel_meas)))) &
				(measurements[measurements[,"ID"]==as.character(isolate(input$sel_meas)),"include"]=="TRUE")
			){
				load(file=file.path(logfile$project_folder,"peaklist",as.character(isolate(input$sel_meas))),envir=as.environment(".GlobalEnv"),verbose=FALSE);
				#load(file=file.path(logfile$project_folder,"peaklist",as.character(2165)),envir=as.environment(".GlobalEnv"),verbose=FALSE);			
				output$file_peak_number<-renderText(as.character(length(peaklist[,1])));	
				blind_aff<-round(
					(sum(peaklist[,colnames(peaklist)=="keep_2"]<Inf))/length(peaklist[,1])*100
				,digits=3)
				output$file_blind_aff<-renderText(as.character(blind_aff));
				output$file_blind_aff2<-renderText(as.character(blind_aff));
				blind_rem<-round(
					(sum(peaklist[,colnames(peaklist)=="keep_2"]<as.numeric(logfile$parameters$blind_threshold)))/length(peaklist[,1])*100
				,digits=3)
				output$file_blind_rem<-renderText(as.character(blind_rem));
				output$file_blind_rem2<-renderText(as.character(blind_rem));
				repl_rem<-round(
					(sum(peaklist[,colnames(peaklist)=="keep"]==0))/length(peaklist[,1])*100
				,digits=3)
				output$file_repl_rem<-renderText(as.character(repl_rem));
				######################################################################
				isolate(refresh_plot$c<-(refresh_plot$c+1))
	  			######################################################################
			}else{
				cat("\n no processed peaklist found for the selected file.")
				isolate(refresh_plot$b<-0)	
				isolate(refresh_plot$c<-0)	
				output$file_peak_number<-renderText("none");
				output$file_blind_aff<-renderText("none");
				output$file_blind_aff2<-renderText("none");
				output$file_blind_rem<-renderText("none");
				output$file_blind_rem2<-renderText("none");
				output$file_repl_rem<-renderText("none");
			}
			##########################################################################		
			# blind/blank peak tagging ###############################################
			if(
				(logfile$workflow[names(logfile$workflow)=="blind"]=="yes") & 
				(measurements[measurements[,"ID"]==as.character(isolate(input$sel_meas)),"Type"]!="blank") &
				(file.exists(file.path(logfile$project_folder,"peaklist",as.character(isolate(input$sel_meas))))) &
				(measurements[measurements[,"ID"]==as.character(isolate(input$sel_meas)),"include"]=="TRUE")
			){			
				output$showblank<-renderText("Blank/blind peak tagging (subtraction) results:")			
				load(file=file.path(logfile[[1]],"peaklist",as.character(isolate(input$sel_meas))),verbose=FALSE);
				if(any(peaklist[,"keep_2"]<Inf)){
					if(isolate(input$blind_boxplot_log)=="TRUE"){
						output$blind_boxplot <- renderPlot({				
							par(mar=c(4.2,4.2,0.8,0.8))
							boxplot(log10(peaklist[peaklist[,"keep_2"]<Inf,"keep_2"]),horizontal=TRUE)
							title(xlab=expression(log[10]*paste(" intensity ratio sample vs. blind",sep=" ")))
							abline(v=0,col="red")
						})
					}else{
						output$blind_boxplot <- renderPlot({		
							par(mar=c(4.2,4.2,0.8,0.8))					
							boxplot(peaklist[peaklist[,"keep_2"]<Inf,"keep_2"],horizontal=TRUE)
							title(xlab="Intensity ratio sample vs. blind")
							abline(v=1,col="red")
						})				
					}
				}else{
					output$blind_boxplot <- renderPlot({
						par(mar=c(4.2,4.2,0.8,0.8))	
						plot.new()
						plot.window(xlim=c(0,1),ylim=c(0,1))
						text(x=.5,y=.5,labels="No peaks matched to any blind/blank peaks",pos=1,cex=1.5)
					})
				}
				rm(peaklist)				
			}else{
				output$showblank<-renderText("No blank/blind subtraction results available.")
			}
			##########################################################################
			pics<-list.files(file.path(logfile[[1]],"pics"))			
			# recalibration ##########################################################
			if(
				any(pics==paste("recal_",as.character(isolate(input$sel_meas)),sep="")) & 
				(logfile$workflow[names(logfile$workflow)=="recal"]=="yes") &
				(measurements[measurements[,"ID"]==as.character(isolate(input$sel_meas)),"include"]=="TRUE")
			){
				output$showrecal<-renderText("Mass recalibration results:")
				expr1<-list(src=file.path(logfile[[1]],"pics",paste("recal_",as.character(isolate(input$sel_meas)),sep="")))
				output$recal_pic<-renderImage(expr1, deleteFile = FALSE)
				cat("\n Found recal_pic")
			}else{
				output$showrecal<-renderText("No mass recalibration results available.")
				cat("\n Not found recal_pic")			
			}
			##########################################################################
			# intensity distribution #################################################
			if(
				any(pics==paste("peakhist_",as.character(isolate(input$sel_meas)),sep="")) &
				(measurements[measurements[,"ID"]==as.character(isolate(input$sel_meas)),"include"]=="TRUE")
			){			
				output$showintensitydistrib<-renderText("Centroid & peak intensity distribution:")
				expr_peakhist<-list(src=file.path(logfile[[1]],"pics",paste("peakhist_",as.character(isolate(input$sel_meas)),sep="")))			
				output$peakhist_pic<-renderImage(expr_peakhist, deleteFile = FALSE)			
				cat("\n Found peakhist_ pic")				
			}else{
				output$showintensitydistrib<-renderText("No centroid/peak intensity distribution available.")
				cat("\n Not found peakhist_ pic")				
			}
			##########################################################################
			# LOD  ###################################################################
			if( 
				file.exists( file.path(logfile[[1]],"results","LOD",paste("plot_LOD_",as.character(isolate(input$sel_meas)),".png",sep="") ) ) &
				(measurements[measurements[,"ID"]==as.character(isolate(input$sel_meas)),"include"]=="TRUE")
			){
				output$showLOD<-renderText("LOD interpolation results:")
				expr_LOD<-list( src=file.path(logfile[[1]],"results","LOD",paste("plot_LOD_",as.character(isolate(input$sel_meas)),".png",sep="")) )
				output$LOD_pic<-renderImage(expr_LOD, deleteFile = FALSE)	
				cat("\n LOD pic file found")
			}else{
				output$showLOD<-renderText("No LOD interpolation available.")
				cat("\n LOD pic file not found")
			}
			##########################################################################			
			output$dowhat<-renderText("Processing per file viewed.");	
		}else{
			output$showblank<-renderText("FALSE")
			output$showrecal<-renderText("FALSE")
			output$showintensitydistrib<-renderText("FALSE")
			output$showLOD<-renderText("FALSE")
			output$dowhat<-renderText("Invalid ID chosen to view processing results.");	
			isolate(refresh_plot$b<-0)	
			isolate(refresh_plot$c<-0)	
		}
    }
	}
	}
})
##############################################################################
observe({
    refresh_plot$c
    input$peak_chromat_refresh
    if(isolate(init$a)=="TRUE" & isolate(refresh_plot$c>0)){				
				use_keep_2<-peaklist[,"keep_2"]
				use_keep_2[use_keep_2==Inf]<-0
	            output$exp_peaklist <- DT::renderDataTable(
	                DT::datatable(
						dat<-data.frame(
							I(as.character(peaklist[,"peak_ID"])),
							round(peaklist[,"m/z"],digits=5),
							round(log10(peaklist[,"var_m/z"]),digits=4),
							round(peaklist[,"m/z_corr"],digits=5),
							round(log10(peaklist[,"max_int"]),digits=5),
							round(log10(peaklist[,"sum_int"]),digits=5),
							round(log10(peaklist[,"int_corr"]),digits=5),
							round(peaklist[,"RT"],digits=1),
							round((peaklist[,"RT"]/60),digits=2),
							round(peaklist[,"minRT"],digits=1),
							round(peaklist[,"maxRT"],digits=1),
							round(peaklist[,"RT_corr"],digits=1),
							I(as.character(peaklist[,"keep"]==1)),
							round(use_keep_2,digits=2)
	                	),
	                	filter = list(position = 'top', clear = FALSE, plain = TRUE),
	                	colnames=c(
	                		"Peak ID","m/z",
	                		"log10 var(m/z)","m/z_recal",
	                		"log10 max int.","log10 sum int.","log10 norm int.",
	                		"RT [s]","RT [min]","min_RT [s]","max_RT [s]","RT_align [s]",
							"in replicates?","int. ratio sample/blind"
	                	),
	                	rownames = FALSE,
	                	selection = list(mode = 'multiple', target = 'row'),
	                	extensions = c('Buttons','FixedHeader','ColReorder'),
						options = list(
							lengthMenu = list(c(200, 500, 1000, -1), list('200', '500', '1000', 'All')),
							fixedHeader = FALSE,
							ordering = TRUE,
							dom = 'Blfrtip',
							buttons = c('excel', 'csv','colvis'),#buttons = c('excel', 'pdf', 'print', 'csv'),
							scrollX = TRUE,
							scrollY = "800px",
							colReorder = TRUE
						)
	                ),server = TRUE
	            )
	}else{
		output$exp_peaklist <- DT::renderDataTable(
			DT::datatable(cbind("no data available","")),
				colnames=c("","")
		)
	}
})
##############################################################################
observe({
    s5<-input$exp_peaklist_rows_selected
    if(isolate(init$a)=="TRUE"){
        if(length(s5)){
        	cat("\n Selected rows: ");print(s5)
        	these_peaks<<-peaklist[s5,"peak_ID"];print(these_peaks)
        	##################################################################
        	if(any(objects(envir=as.environment(".GlobalEnv"))=="MSlist")){
        		if(any(names(MSlist)=="File_ID")){
        			if(MSlist[["File_ID"]]!=as.character(isolate(input$sel_meas))){ # File_ID does not match
						load(file.path(logfile[[1]],"MSlist",as.character(isolate(input$sel_meas))),envir=as.environment(".GlobalEnv")) 
        			}
        		}else{ # available MSlist not with File_ID yet
					load(file.path(logfile[[1]],"MSlist",as.character(isolate(input$sel_meas))),envir=as.environment(".GlobalEnv"))  
					MSlist[["File_ID"]]<-as.character(isolate(input$sel_meas))
        		}
        	}else{ # no MSlist in GlobalEnv
				load(file.path(logfile[[1]],"MSlist",as.character(isolate(input$sel_meas))),envir=as.environment(".GlobalEnv"))  
        	}
        	##################################################################
	        isolate(refresh_plot$b<-(refresh_plot$b+1))
        	##################################################################
        }else{
        	cat("\n Selected nothing: ");print(s5)
        }
	}
})
##############################################################################
# observers on peak(s) chromatogram(s) plot ##################################
##############################################################################
observe({ # seconds <-> minutes switch when zoomed ###########################
	input$peak_chromat_time
	if(isolate(init$a)=="TRUE" & isolate(ranges_peaks_mz_RT$xchroma[1]!=FALSE)){
		if(isolate(input$peak_chromat_time)=="minutes"){
			isolate(ranges_peaks_mz_RT$xchroma<-(ranges_peaks_mz_RT$xchroma/60))
		}
		if(isolate(input$peak_chromat_time)=="seconds"){
			isolate(ranges_peaks_mz_RT$xchroma<-(ranges_peaks_mz_RT$xchroma*60))
		}
	}
})
##############################################################################
observe({ # nomralization switch when zoomed #################################
	input$peak_chromat_norm
	if(isolate(init$a)=="TRUE" & isolate(ranges_peaks_mz_RT$ychroma[1]!=FALSE)){
		isolate(ranges_peaks_mz_RT$ychroma<-FALSE)
	}
})
##############################################################################
observe({
	refresh_plot$b
	input$peak_chromat_norm
	input$peak_chromat_time
	input$peak_chromat_type
    if(isolate(init$a)=="TRUE" & isolate(refresh_plot$b>0)){
	    output$peak_chromat <- renderPlot({
	        par(mar=c(5,4,1,.8))
	        enviMass:::plotchromat(
	          	MSlist,
	           	peakIDs=these_peaks,
	           	RTlim=isolate(ranges_peaks_mz_RT$xchroma),
	           	Intlim=isolate(ranges_peaks_mz_RT$ychroma),
	           	normalize=as.logical(isolate(input$peak_chromat_norm)),
	           	n_col=dim(peaklist)[1],
	           	set_RT=isolate(input$peak_chromat_time),
	           	chromat_full=input$peak_chromat_type
	        );
	    },res=100) 
	}
})
##############################################################################
observeEvent(input$peak_chromat_dblclick, { 
	if(isolate(init$a)=="TRUE"){
          brush <- isolate(input$peak_chromat_brush)
          if (!is.null(brush)) {
            cat("\n Zoom in_1e")
            isolate(ranges_peaks_mz_RT$xchroma <- c(brush$xmin, brush$xmax))
            isolate(ranges_peaks_mz_RT$ychroma <- c(brush$ymin, brush$ymax))            
          } else {
            cat("\n Zoom out full_1e")
            isolate(ranges_peaks_mz_RT$xchroma <- FALSE)
            isolate(ranges_peaks_mz_RT$ychroma <- FALSE)
          }
          isolate(refresh_plot$b<-(refresh_plot$b+1)) # valid in both cases
    }
})
##############################################################################
observeEvent(input$peak_chromat_click, { # NOTE: brushing already triggers a click -> use brush with delay=0, which embeds the slower click
	if(isolate(init$a)=="TRUE"){
          cat("\n Zoom out part_1_ae")
          brush <- isolate(input$peak_chromat_brush)
          if (is.null(brush)) {
              cat("\n Zoom out part_1_be")
              if(isolate(ranges_peaks_mz_RT$xchroma[1])!=FALSE){
                old_range_dmass<-abs(isolate(ranges_peaks_mz_RT$xchroma[2]-ranges_peaks_mz_RT$xchroma[1]))
                isolate(ranges_peaks_mz_RT$xchroma[1]<-ranges_peaks_mz_RT$xchroma[1]-.3*old_range_dmass)
                isolate(ranges_peaks_mz_RT$xchroma[2]<-ranges_peaks_mz_RT$xchroma[2]+.3*old_range_dmass)
              }
              if(isolate(ranges_peaks_mz_RT$ychroma[1])!=FALSE){             
                old_range_dmass<-abs(isolate(ranges_peaks_mz_RT$ychroma[2]-ranges_peaks_mz_RT$ychroma[1]))
                isolate(ranges_peaks_mz_RT$ychroma[1]<-ranges_peaks_mz_RT$ychroma[1]-.3*old_range_dmass)
                isolate(if(ranges_peaks_mz_RT$ychroma[1]<0){ranges_peaks_mz_RT$ychroma[1]<-0})
                isolate(ranges_peaks_mz_RT$ychroma[2]<-ranges_peaks_mz_RT$ychroma[2]+.3*old_range_dmass)
              }  
              isolate(refresh_plot$b<-(refresh_plot$b+1))
          }else{
            cat("\n Doing hover_e - nothing")
          }   
    }
})     
##############################################################################
##############################################################################




##############################################################################
# retrieve peak information ##################################################
##############################################################################
observe({
	refresh_plot$a
    input$sel_meas_ID
	input$peaks_mz_RT_use_peaks
	input$peaks_mz_RT_use_raw
	input$peaks_mz_RT_use_IDs
	input$peaks_mz_RT_use_window
	input$peaks_mz_RT_use_window_mass
	input$peaks_mz_RT_use_window_RT
	input$peaks_mz_RT_use_bar	
	input$peaks_mz_RT_use_bar_value	
	input$peaks_mz_RT_use_window_RT_tol
	input$plot_filter_intensity
	input$plot_filter_blind
	input$plot_filter_replicates
	if(isolate(init$a)=="TRUE"){
		if(!is.na(isolate(input$sel_meas_ID))){ 
			if(!any(objects(envir=as.environment(".GlobalEnv"))=="atit")){ # atit -> dont load the same MSlist twice = too slow
				assign("atit",0,envir=as.environment(".GlobalEnv"))
			}	
			if(	
				file.exists(file.path(logfile[[1]],"MSlist",as.character(isolate(input$sel_meas_ID)))) &
				file.exists(file.path(logfile[[1]],"peaklist",as.character(isolate(input$sel_meas_ID)))) 			
			){
				if(verbose){cat("\n Plotting")}
				if((isolate(input$sel_meas_ID)!=atit)){
					if(any(objects(envir=as.environment(".GlobalEnv"))=="MSlist")){rm(MSlist,envir=as.environment(".GlobalEnv"))}
					if(any(objects()=="MSlist")){rm(MSlist)}				
					load(file.path(logfile[[1]],"MSlist",as.character(isolate(input$sel_meas_ID))), envir=as.environment(".GlobalEnv"))
					load(file.path(logfile[[1]],"peaklist",as.character(isolate(input$sel_meas_ID))), envir=as.environment(".GlobalEnv"))
					if(verbose){cat("\nafter:\n "); print(gc());cat("\n MSlist file loaded")}
					assign("atit",isolate(input$sel_meas_ID),envir=as.environment(".GlobalEnv"))
					max_int<-max(max(MSlist[["Peaklist"]][,"max_int"]),max(MSlist[["Scans"]][[2]][,"intensity"]))
					min_int<-min(min(MSlist[["Peaklist"]][,"max_int"]),min(MSlist[["Scans"]][[2]][,"intensity"]))
					updateSliderInput(session,"plot_filter_intensity", min=round(log10(min_int),digits=1), max=round(log10(max_int),digits=1), value=c(log10(min_int),log10(max_int)))
					########################################################################################
					measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
     				if(any(measurements[,"ID"]==as.character(isolate(input$sel_meas_ID)))){
    					output$file_viewer_name<-renderText(paste("File name: ",measurements[measurements[,"ID"]==as.character(isolate(input$sel_meas_ID)),"Name"],sep=""))
    					output$file_viewer_type<-renderText(paste("File type: ",measurements[measurements[,"ID"]==as.character(isolate(input$sel_meas_ID)),"Type"],sep=""))
     					output$file_viewer_mode<-renderText(paste("Ionization mode: ",measurements[measurements[,"ID"]==as.character(isolate(input$sel_meas_ID)),"Mode"],sep=""))
   					}else{
						output$file_viewer_name<-renderText("Invalid file ID")
						output$file_viewer_type<-renderText("")
						output$file_viewer_mode<-renderText("")
   					}
				}
				# prepare plotting information #############################################################
				if(
					(exists("MSlist")) &
					(exists("peaklist"))
				){
					if(!is.null(isolate(ranges_peaks_mz_RT$x))){
						x_lim<-isolate(ranges_peaks_mz_RT$x)
					}else{
						x_lim<-c(min(MSlist[["Scans"]][[2]][,"m/z"]),max(MSlist[["Scans"]][[2]][,"m/z"]))
					}
					if(!is.null(isolate(ranges_peaks_mz_RT$y))){
						y_lim<-isolate(ranges_peaks_mz_RT$y)
					}else{
						y_lim<-c(min(MSlist[["Scans"]][[2]][,"RT"]),max(MSlist[["Scans"]][[2]][,"RT"]))
					}		
					use_these<-which(
						(peaklist[,"m/z"]>=x_lim[1]) &
						(peaklist[,"m/z"]<=x_lim[2]) &
						(peaklist[,"RT"]>=y_lim[1]) &
						(peaklist[,"RT"]<=y_lim[2]) &
						(peaklist[,"max_int"]>=isolate(10^input$plot_filter_intensity[1])) &
						(peaklist[,"max_int"]<=isolate(10^input$plot_filter_intensity[2])) 
					)
					if((length(use_these)>0) & isolate(input$plot_filter_blind)){
						use_these<-use_these[peaklist[use_these,"keep_2"]>=as.numeric(logfile$parameters$blind_threshold)]
					}
					if((length(use_these)>0) & isolate(input$plot_filter_replicates)){
						use_these<-use_these[peaklist[use_these,"keep"]==1]
					}				
					if(length(use_these)>0){
						z_lim<-c(min(peaklist[use_these,"max_int"]),max(peaklist[use_these,"max_int"]))
					}else{
						z_lim<-c(min(peaklist[,"max_int"]),max(peaklist[,"max_int"]))
					}
					if(isolate(input$peaks_mz_RT_use_raw)){ # must be prerequisite to evaluate those, below
						if(is.null(isolate(ranges_peaks_mz_RT$x))){
							x_min<-min(MSlist[["Scans"]][[2]][,"m/z"])
							x_max<-max(MSlist[["Scans"]][[2]][,"m/z"])							
						}else{
							x_min<-min(MSlist[["Scans"]][[2]][
								MSlist[["Scans"]][[2]][,"m/z"]>=isolate(ranges_peaks_mz_RT$x[1])
							,"m/z"])
							x_max<-max(MSlist[["Scans"]][[2]][
								MSlist[["Scans"]][[2]][,"m/z"]<=isolate(ranges_peaks_mz_RT$x[2])
							,"m/z"])	
						}
						if(is.null(isolate(ranges_peaks_mz_RT$y))){
							y_min<-min(MSlist[["Scans"]][[2]][,"RT"])
							y_max<-max(MSlist[["Scans"]][[2]][,"RT"])							
						}else{
							y_min<-min(MSlist[["Scans"]][[2]][
								MSlist[["Scans"]][[2]][,"RT"]>=isolate(ranges_peaks_mz_RT$y[1])
							,"RT"])
							y_max<-max(MSlist[["Scans"]][[2]][
								MSlist[["Scans"]][[2]][,"RT"]<=isolate(ranges_peaks_mz_RT$y[2])
							,"RT"])	
						}					
						those<-which(
							(MSlist[["Scans"]][[2]][,"m/z"]>=x_min) &
							(MSlist[["Scans"]][[2]][,"m/z"]<=x_max) &
							(MSlist[["Scans"]][[2]][,"RT"]>=y_min) &
							(MSlist[["Scans"]][[2]][,"RT"]<=y_max) &
							(MSlist[["Scans"]][[2]][,"intensity"]>=isolate(10^input$plot_filter_intensity[1])) &
							(MSlist[["Scans"]][[2]][,"intensity"]<=isolate(10^input$plot_filter_intensity[2]))						
						)			
						if(length(those)<=1E6){
							colorit<-rep("gray",length(those))
							colorit[MSlist[["Scans"]][[2]][those,"peakID"]!=0]<-"red"
						}
					}else{
						those<-c()
					}
					####################################################################################			
					output$plot_peaks_mz_RT <- renderPlot({
						par(mar=c(4, 4, 3, .1))
						plot.new()
						plot.window(xlim=x_lim,ylim=y_lim)
						title(
							xlab="m/z [Th]", ylab="RT [s]",
							main="Draw rectangles and double-click into them to zoom in, double-click again to zoom fully out. Bottom plots adapt accordingly.",cex.main=.75
						)
						box();axis(1);axis(2);
						# add raw data ? #########################################################
						if(isolate(input$peaks_mz_RT_use_raw)){
							if(length(those)<=1E5 & length(those)>0){
								points(
									MSlist[["Scans"]][[2]][those,"m/z"],
									MSlist[["Scans"]][[2]][those,"RT"],
									pch=19,cex=.5,col=colorit
								)
							}else{
								smoothScatter(
									x=MSlist[["Scans"]][[2]][those,"m/z"], 
									y=MSlist[["Scans"]][[2]][those,"RT"],
									 colramp = colorRampPalette(c("white", "red")),
									nbin = 200, add = TRUE
								)
							}
						}
						# add picked peaks ? ######################################################
						if(isolate(input$peaks_mz_RT_use_peaks)){		
							if(length(use_these)>0){
								if(!isolate(input$peaks_mz_RT_use_IDs)){
									colorit<-"black"
								}else{
									colorit<-"darkgrey"
								}
								points(
									peaklist[use_these,"m/z"],
									peaklist[use_these,"RT"],						
									pch=21,cex=.8,col=colorit
								)
							}
						}
						# add peak IDs? ###########################################################
						if(isolate(input$peaks_mz_RT_use_IDs)){
							if(length(use_these)>0){
								text(						
									peaklist[use_these,"m/z"],
									peaklist[use_these,"RT"], 
									labels = as.character(peaklist[use_these,"peak_ID"]),
									pos = NULL, col	= "darkred", cex=.6
								)
							}					
						}					
						# add search window ? ######################################################
						if(isolate(input$peaks_mz_RT_use_window)){	
							del_mz<-((x_lim[2]-x_lim[1])/15)
							at_mz<-isolate(as.numeric(input$peaks_mz_RT_use_window_mass))
							at_RT<-isolate(as.numeric(input$peaks_mz_RT_use_window_RT))
							rect(
								xleft=(at_mz-del_mz), 
								ybottom=(at_RT-isolate(as.numeric(input$peaks_mz_RT_use_window_RT_tol))), 
								xright=(at_mz+del_mz), 
								ytop=(at_RT+isolate(as.numeric(input$peaks_mz_RT_use_window_RT_tol))),
								col=NULL, border="blue",lwd=2)
							if(isolate(input$peaks_mz_RT_use_bar)){		
								del_ppm<-(at_mz*isolate(as.numeric(input$peaks_mz_RT_use_bar_value))/1E6)
								lines(
									x=c((at_mz-del_ppm),(at_mz+del_ppm)),
									y=c(at_RT,at_RT),
									col="blue",lwd=2)
							}	
						}
						############################################################################	
					}, res = 100, execOnResize=TRUE)
					####################################################################################
					output$plot_peaks_mz_int <- renderPlot({		
							par(mar=c(4, 4, .8, .2))
							plot.new()
							plot.window(xlim=x_lim,ylim=c(0,z_lim[2]))
							title(xlab="m/z [Th]", ylab="Intensity")
							box();axis(1);axis(2);
							# add raw data ? #########################################################			
							if(isolate(input$peaks_mz_RT_use_raw)){			
								if(length(those)<=1E5 & length(those)>0){
									points(
										MSlist[["Scans"]][[2]][those,"m/z"],
										MSlist[["Scans"]][[2]][those,"intensity"],
										type="h",pch=19,cex=.5,col=colorit
									)
								}						
							}
							# add picked peaks ? #####################################################
							if(isolate(input$peaks_mz_RT_use_peaks)){						
								if(length(use_these)>0){
									if(!isolate(input$peaks_mz_RT_use_IDs)){
										colorit<-"black"
									}else{
										colorit<-"darkgrey"
									}
									points(
										peaklist[use_these,"m/z"],
										peaklist[use_these,"max_int"],						
										pch=21,cex=.8,col=colorit,type="h"
									)					
								}
							}
							# add ppm bar? #############################################################
							if(isolate(input$peaks_mz_RT_use_window)){	
								if(isolate(input$peaks_mz_RT_use_bar)){	
									at_mz<-isolate(as.numeric(input$peaks_mz_RT_use_window_mass))							
									del_ppm<-(at_mz*isolate(as.numeric(input$peaks_mz_RT_use_bar_value))/1E6)
									lines(
										x=c((at_mz-del_ppm),(at_mz+del_ppm)),
										y=c(0.5*z_lim[2],.5*z_lim[2]),
										col="blue",lwd=2)
								}	
							}						
					}, res = 100, execOnResize=TRUE)	
					####################################################################################
					output$plot_peaks_RT_int <- renderPlot({		
							par(mar=c(4, 4, .8, .2))
							plot.new()
							plot.window(xlim=y_lim,ylim=c(0,z_lim[2]))
							title(xlab="RT", ylab="Intensity")
							box();axis(1);axis(2);
							# add raw data ? #########################################################			
							if(isolate(input$peaks_mz_RT_use_raw)){			
								if(length(those)<=1E5 & length(those)>0){
									points(
										MSlist[["Scans"]][[2]][those,"RT"],
										MSlist[["Scans"]][[2]][those,"intensity"],
										type="h",pch=19,cex=.5,col=colorit
									)
								}						
							}
							# add picked peaks ? #####################################################
							if(isolate(input$peaks_mz_RT_use_peaks)){						
								if(length(use_these)>0){
									if(!isolate(input$peaks_mz_RT_use_IDs)){
										colorit<-"black"
									}else{
										colorit<-"darkgrey"
									}
									points(
										peaklist[use_these,"RT"],
										peaklist[use_these,"max_int"],						
										pch=21,cex=.8,col=colorit,type="h"
									)					
								}
							}
							# add RT window? ##########################################################
							if(isolate(input$peaks_mz_RT_use_window)){	
									at_RT<-isolate(as.numeric(input$peaks_mz_RT_use_window_RT))
									lines(
										x=c(
											(at_RT-isolate(as.numeric(input$peaks_mz_RT_use_window_RT_tol))),
											(at_RT+isolate(as.numeric(input$peaks_mz_RT_use_window_RT_tol)))
										),
										y=c(.5*z_lim[2],.5*z_lim[2]),
										col="blue",lwd=2)
								}	
					}, res = 100, execOnResize=TRUE)	
					####################################################################################		
					# plotly output ####################################################################
					# only peaks? ######################################################################
					if(
						( isolate(input$peaks_mz_RT_use_peaks) & !isolate(input$peaks_mz_RT_use_raw) & (length(use_these)>0) ) ||
						( isolate(input$peaks_mz_RT_use_peaks) & (length(those)==0) & (length(use_these)>0) )
					){
						if(length(use_these)<=1E5){
							if(verbose){cat("\n Plotting only peaks")}
							sub_peaks<-as.data.frame(peaklist[use_these,c("m/z","RT","max_int"),drop=FALSE])
							names(sub_peaks)<-c("m_z","RT","Intensity")				
							sub_peaks[,"Intensity"]<-(sub_peaks[,"Intensity"]/2)
							output$plot_peaks_3D <- renderPlotly({
								p <- plotly::plot_ly(type="scatter3d",mode="markers",showlegend=TRUE)%>%
								plotly::add_trace(p,
									x = ~m_z, y = ~RT, z = ~Intensity,
									data = sub_peaks,
									color=I("black"),
									size = I(1),
									name = "",
									error_z=list(
										color="black",
										thickness=0,
										symmetric = TRUE, 
										type = "data" ,
										array = sub_peaks$Intensity
									)
								)	
							})					
						}else{
							output$plot_peaks_3D <- renderPlotly({plotly::plot_ly(type="scatter3d",mode="markers",showlegend=TRUE)})
						}
					}
					# only raw data? ###########################################################						
					if(
						( !isolate(input$peaks_mz_RT_use_peaks) & isolate(input$peaks_mz_RT_use_raw) & (length(those)>0) ) ||
						( isolate(input$peaks_mz_RT_use_raw) & (length(use_these)==0) & (length(those)>0) )
					){
						if(length(those)<=1E5){  # implies number of peaks is lower, too
							if(verbose){cat("\n Plotting only raw data")}
							sub_MSlist<-as.data.frame(MSlist[["Scans"]][[2]][those,c("m/z","RT","intensity","peakID"),drop=FALSE])
							names(sub_MSlist)<-c("m_z","RT","Intensity","peakID")		
							sub_MSlist[,"Intensity"]<-(sub_MSlist[,"Intensity"]/2)
							if(any(sub_MSlist[,"peakID"]!=0)){ # raw data included in peaks available?
								output$plot_peaks_3D <- renderPlotly({
									p <- plotly::plot_ly(type="scatter3d",mode="markers",showlegend=TRUE)%>%
									plotly::add_trace(p,
										x = ~m_z, y = ~RT, z = ~Intensity,
										data = sub_MSlist[sub_MSlist[,"peakID"]==0,],
										color=I("gray"),
										size = I(1),
										name = "",
										error_z=list(
											color="gray",
											thickness=0,
											symmetric = TRUE, 
											type = "data" ,
											array = sub_MSlist[sub_MSlist[,"peakID"]==0,]$Intensity
										)
									)%>%	
									plotly::add_trace(p,
										x = ~m_z, y = ~RT, z = ~Intensity,
										data = sub_MSlist[sub_MSlist[,"peakID"]!=0,],
										color=I("red"),
										size = I(1),
										name = "",
										error_z=list(
											color="red",
											thickness=0,
											symmetric = TRUE, 
											type = "data" ,
											array = sub_MSlist[sub_MSlist[,"peakID"]!=0,]$Intensity
										)
									)
								})				
							}else{
								output$plot_peaks_3D <- renderPlotly({
									p <- plotly::plot_ly(type="scatter3d",mode="markers",showlegend=TRUE)%>%
									plotly::add_trace(p,
										x = ~m_z, y = ~RT, z = ~Intensity,
										data = sub_MSlist[sub_MSlist[,"peakID"]==0,],
										color=I("gray"),
										size = I(1),
										name = "",
										error_z=list(
											color="gray",
											thickness=0,
											symmetric = TRUE, 
											type = "data" ,
											array = sub_MSlist[sub_MSlist[,"peakID"]==0,]$Intensity
										)
									)
								})								
							}
						}else{
							output$plot_peaks_3D <- renderPlotly({plotly::plot_ly(type="scatter3d",mode="markers",showlegend=TRUE)})					
						}
					}					
					# peaks & raw data? ########################################################
					if( isolate(input$peaks_mz_RT_use_peaks) & isolate(input$peaks_mz_RT_use_raw) & (length(those)>0) & (length(use_these)>0) ){
						if(length(those)<=1E5){ # implies number of peaks is lower, too
							if(verbose){cat("\n Plotting peaks & raw data")}
							sub_MSlist<-as.data.frame(MSlist[["Scans"]][[2]][those,c("m/z","RT","intensity","peakID"),drop=FALSE])
							names(sub_MSlist)<-c("m_z","RT","Intensity","peakID")		
							sub_MSlist[,"Intensity"]<-(sub_MSlist[,"Intensity"]/2)
							sub_peaks<-as.data.frame(peaklist[use_these,c("m/z","RT","max_int"),drop=FALSE])
							names(sub_peaks)<-c("m_z","RT","Intensity")				
							sub_peaks[,"Intensity"]<-(sub_peaks[,"Intensity"]/2)
							output$plot_peaks_3D <- renderPlotly({
								p <- plotly::plot_ly(type="scatter3d",mode="markers",showlegend=TRUE)%>%
								plotly::add_trace(p,
									x = ~m_z, y = ~RT, z = ~Intensity,
									data = sub_MSlist[sub_MSlist[,"peakID"]==0,],
									color=I("gray"),
									size = I(1),
									name = "",
									error_z=list(
										color="gray",
										thickness=0,
										symmetric = TRUE, 
										type = "data" ,
										array = sub_MSlist[sub_MSlist[,"peakID"]==0,]$Intensity
									)
								)%>%	
								plotly::add_trace(p,
									x = ~m_z, y = ~RT, z = ~Intensity,
									data = sub_MSlist[sub_MSlist[,"peakID"]!=0,],
									color=I("red"),
									size = I(1),
									name = "",
									error_z=list(
										color="red",
										thickness=0,
										symmetric = TRUE, 
										type = "data" ,
										array = sub_MSlist[sub_MSlist[,"peakID"]!=0,]$Intensity
									)
								)%>%
								plotly::add_trace(p,
									x = ~m_z, y = ~RT, z = ~Intensity,
									data = sub_peaks,
									color=I("black"),
									size = I(1),
									name = "",
									error_z=list(
										color="black",
										thickness=0,
										symmetric = TRUE, 
										type = "data" ,
										array = sub_peaks$Intensity
									)
								)	
							})							
						}else{
							output$plot_peaks_3D <- renderPlotly({plotly::plot_ly(type="scatter3d",mode="markers",showlegend=TRUE)})					
						}
					}						
					# nothing? ##################################################################
					if( !isolate(input$peaks_mz_RT_use_peaks) & !isolate(input$peaks_mz_RT_use_raw)){
						output$plot_peaks_3D <- renderPlotly({plotly::plot_ly(type="scatter3d",mode="markers",showlegend=TRUE)})	
					}
					#############################################################################
				}else{
					output$file_viewer_name<-renderText("Invalid file ID")
					output$file_viewer_type<-renderText("")
					output$file_viewer_mode<-renderText("")					
				}	
			}else{
				output$plot_peaks_mz_RT <- renderPlot({})
				output$plot_peaks_mz_int <- renderPlot({})
				output$plot_peaks_RT_int <- renderPlot({})
				output$plot_peaks_3D <- renderPlotly({plotly::plot_ly(type="scatter3d",mode="markers",showlegend=TRUE)})	
				output$file_viewer_name<-renderText("Invalid file ID")
				output$file_viewer_type<-renderText("")
				output$file_viewer_mode<-renderText("")
			}
		}
	}
})

# When a double-click happens, check if there's a brush on the plot.
# If so, zoom to the brush bounds; if not, reset the zoom.
observeEvent(input$plot_peaks_mz_RT_dblclick, { 
    brush <- input$plot_peaks_mz_RT_brush
    if (!is.null(brush)) {
		isolate(ranges_peaks_mz_RT$x <- c(brush$xmin, brush$xmax))
		isolate(ranges_peaks_mz_RT$y <- c(brush$ymin, brush$ymax))
    } else {
		isolate(ranges_peaks_mz_RT$x <- NULL)
		isolate(ranges_peaks_mz_RT$y <- NULL)
    }
	refresh_plot$a<-(refresh_plot$a+1)
	cat("\n Zooming with brush")
})

observe({
	input$peaks_mz_RT_zoom_out
	if(isolate(init$a)=="TRUE"){
		if(!is.na(isolate(input$sel_meas_ID))){ 
			if(!is.null(isolate(ranges_peaks_mz_RT$x))){
				cat("\n Zooming out on X")
				old_range<-abs(isolate(ranges_peaks_mz_RT$x[2]-ranges_peaks_mz_RT$x[1]))
				isolate(ranges_peaks_mz_RT$x[1]<-ranges_peaks_mz_RT$x[1]-.5*old_range)
				isolate(ranges_peaks_mz_RT$x[2]<-ranges_peaks_mz_RT$x[2]+.5*old_range)
			}
			
			
		refresh_plot$a<-(refresh_plot$a+1)		
		}
	}
})

observe({
    input$sel_meas_ID
	input$sel_peak_ID
	if(!is.na(isolate(input$sel_meas_ID))){ # if user deletes entry!
		if(isolate(input$sel_meas_ID)!=0){
			if(atit==isolate(input$sel_meas_ID) & exists("MSlist")){ # above MSlist upload worked?
				if( !is.na(isolate(input$sel_peak_ID)) & 
					(isolate(input$sel_peak_ID)!=0) & 
					any(MSlist[[8]][,10]==isolate(input$sel_peak_ID)) &
					any(objects(envir=as.environment(".GlobalEnv"))=="MSlist")
				){
					EIC_ID<<-unique(MSlist[[8]][MSlist[[8]][,10]==isolate(input$sel_peak_ID),9]);
					peakit<<-MSlist[[4]][[2]][c(MSlist[[7]][as.numeric(isolate(input$sel_peak_ID)),1]:MSlist[[7]][
						as.numeric(isolate(input$sel_peak_ID)),2]),]			
					if(length(peakit)>7){
						EICit<<-MSlist[[4]][[2]][c(MSlist[[6]][EIC_ID,1]:MSlist[[6]][EIC_ID,2]),]
						output$EIC1 <- renderPlot({
							if(length(EICit)>7){
								plot(EICit[,3],EICit[,2],type="h",col="darkgrey",xlab="RT",ylab="Intensity",main="EIC (grey) & Peak (red)",xlim=c(min(MSlist[[4]][[1]]),max(MSlist[[4]][[1]])))
							}else{
								plot(EICit[3],EICit[2],type="h",col="darkgrey",xlab="RT",ylab="Intensity",main="EIC (grey) & Peak (red)")
							}
							if(length(peakit)>7){	
								points(peakit[,3],peakit[,2],type="h",col="red",lwd=2)
							}else{
								points(peakit[3],peakit[2],type="h",col="red",lwd=2)				
							}
						})	
						output$EIC2 <- renderPlot({
							if(length(EICit)>7){
								plot(EICit[,3],EICit[,2],type="h",col="darkgrey",xlab="RT",ylab="Intensity")
							}else{
								plot(EICit[3],EICit[2],type="h",col="darkgrey",xlab="RT",ylab="Intensity")
							}
							if(length(peakit)>7){	
								points(peakit[,3],peakit[,2],type="h",col="red",lwd=2)
							}else{
								points(peakit[3],peakit[2],type="h",col="red",lwd=2)				
							}
						})	
						output$EIC3 <- renderPlot({
							if(length(EICit)>7){
								plot(EICit[,3],EICit[,1],pch=19,col="darkgrey",xlab="RT",ylab="m/z")			
							}else{
								plot(EICit[3],EICit[1],pch=19,col="darkgrey",xlab="RT",ylab="m/z")
							}
							if(length(peakit)>7){	
								points(peakit[,3],peakit[,1],pch=19,col="red",cex=1.5)
							}else{
								points(peakit[,3],peakit[,1],pch=19,col="red",cex=1.5)				
							}
						})	
						cat("\n EIC & peak extracted")
					}else{
						cat("\n Peak based on single measurement - plotting skipped.")
					}
				}	
			}else{
				output$EIC1 <- renderPlot({plot.new()})
				output$EIC2 <- renderPlot({plot.new()})
				output$EIC3 <- renderPlot({plot.new()})
			}
		}
	}
})
##############################################################################

##############################################################################
# PROFILE FITERING - update results for changes in ion mode selection ########
##############################################################################
maincalc3<-reactive({
	input$Ion_mode
	if( (isolate(init$a)=="TRUE") & (isolate(input$Ion_mode)=="positive") ){
		exprprofnorm_pos<-list(src=file.path(logfile[[1]],"pics","profnorm_pos"))
		#output$profnorm<-renderImage(exprprofnorm_pos, deleteFile = FALSE)
		exprprofcount_pos<-list(src=file.path(logfile[[1]],"pics","profcount_pos"))
		output$profcount<-renderImage(exprprofcount_pos, deleteFile = FALSE)
		if(any(objects(envir=as.environment(".GlobalEnv"))=="profileList")){ rm(profileList, inherits = TRUE) }
		if(file.exists(file.path(logfile[[1]],"results","profileList_pos"))){
			if(any(objects(envir=as.environment(".GlobalEnv"))=="profileList_pos")){rm(profileList_pos,envir=as.environment(".GlobalEnv"))}
			if(any(objects()=="profileList_pos")){rm(profileList_pos)}				
			load(file=file.path(as.character(logfile[[1]]),"results","profileList_pos"),envir=as.environment(".GlobalEnv"), verbose=TRUE);
			assign("profileList",profileList_pos,envir=as.environment(".GlobalEnv"));
		}else{
			return("No profiles available for this ionization mode")
		}
		if(file.exists(file.path(logfile[[1]],"results","links_profiles_pos"))){
			if(any(objects(envir=as.environment(".GlobalEnv"))=="links_profiles_pos")){rm(links_profiles_pos,envir=as.environment(".GlobalEnv"))}
			if(any(objects()=="links_profiles_pos")){rm(links_profiles_pos)}				
			load(file=file.path(as.character(logfile[[1]]),"results","links_profiles_pos"),envir=as.environment(".GlobalEnv"), verbose=TRUE);
			assign("links_profiles",links_profiles_pos,envir=as.environment(".GlobalEnv"));
		}				
		expr4p<-list(src=file.path(logfile[[1]],"pics","boxprofile_pos"))
		output$boxprofile<-renderImage(expr4p, deleteFile = FALSE)		
		isolate(init$b<<-(init$b+1))
		if(any(objects()=="profileList")){stop("illegal profpeaks2 found, #1");}
		if(any(objects()=="profpeaks2")){stop("illegal profpeaks found, #1");}
		return("Select ionization (switch to negative):\n")
	}
	if( (isolate(init$a)=="TRUE") &  (isolate(input$Ion_mode)=="negative") ){
		exprprofnorm_neg<-list(src=file.path(logfile[[1]],"pics","profnorm_neg"))
		output$profnorm<-renderImage(exprprofnorm_neg, deleteFile = FALSE)
		exprprofcount_neg<-list(src=file.path(logfile[[1]],"pics","profcount_neg"))
		output$profcount<-renderImage(exprprofcount_neg, deleteFile = FALSE)
		if(any(objects(envir=as.environment(".GlobalEnv"))=="profileList")){ rm(profileList, inherits = TRUE) }
		if(file.exists(file.path(logfile[[1]],"results","profileList_neg"))){
			if(any(objects(envir=as.environment(".GlobalEnv"))=="profileList_neg")){rm(profileList_neg,envir=as.environment(".GlobalEnv"))}
			if(any(objects()=="profileList_neg")){rm(profileList_neg)}				
			load(file=file.path(as.character(logfile[[1]]),"results","profileList_neg"),envir=as.environment(".GlobalEnv"), verbose=TRUE);
			assign("profileList",profileList_neg,envir=as.environment(".GlobalEnv"));
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
		if(any(objects()=="profileList")){stop("illegal profpeaks2 found, #2");}
		if(any(objects()=="profpeaks2")){stop("illegal profpeaks found, #2");}
		return("Select ionization (switch to positive):\n")	
	}
})
output$had_ion<-renderText(paste(maincalc3())) 
##############################################################################

##############################################################################
# PROFILE FITERING - Sort and filter the profile list ########################
##############################################################################
maincalc6<-reactive({
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
	input$filterProf_count
	input$filterProf_medianblind
	input$filterProf_medianblind_value
	input$filterProf_minMD
	input$filterProf_maxMD
	input$filterProf_components
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
		if(any(objects(envir=as.environment(".GlobalEnv"))=="profpeaks2")){rm(profpeaks2,envir=as.environment(".GlobalEnv"))}
		if(any(objects()=="profpeaks2")){rm(profpeaks2)}		
		assign("profpeaks2",profileList[["index_prof"]],envir=as.environment(".GlobalEnv"));
		if(any(objects()=="profileList")){stop("illegal profileList found, #3");}
		if(any(objects()=="profpeaks2")){stop("illegal profileList found, #3");}
		###################################################################################################		
		profpeaks2<<-profpeaks2[profpeaks2[,"mean_mz"]>=isolate(input$filterProf_minmass),,drop = FALSE]
		profpeaks2<<-profpeaks2[profpeaks2[,"mean_mz"]<=isolate(input$filterProf_maxmass),,drop = FALSE]
		profpeaks2<<-profpeaks2[profpeaks2[,"mean_RT"]>=isolate(input$filterProf_minrt),,drop = FALSE]
		profpeaks2<<-profpeaks2[profpeaks2[,"mean_RT"]<=isolate(input$filterProf_maxrt),,drop = FALSE]
		profpeaks2<<-profpeaks2[profpeaks2[,"Mass defect"]>=isolate(input$filterProf_minMD),,drop = FALSE]
		profpeaks2<<-profpeaks2[profpeaks2[,"Mass defect"]<=isolate(input$filterProf_maxMD),,drop = FALSE]		
		###################################################################################################		
		if(isolate(input$filterProf_medianblind)=="yes"){
			profpeaks2<<-profpeaks2[(profpeaks2[,"above_blind?"]>=as.numeric(isolate(input$filterProf_medianblind_value))),,drop = FALSE] 
		}
		if(isolate(input$filterProf_notblind)=="yes"){
			profpeaks2<<-profpeaks2[profpeaks2[,"in_blind?"]==0,,drop = FALSE] # not in blind, = profileList[[7]][k,8]
		}
		###################################################################################################	
		if(isolate(input$filterProf_sort)=="ID (increasing)"){
			profpeaks2<<-profpeaks2[order(profpeaks2[,"profile_ID"], decreasing = FALSE),,drop = FALSE]
			sort_by<<-"profile_ID";sort_by_decreasing <<- "FALSE"
		}
		if(isolate(input$filterProf_sort)=="mean m/z (increasing)"){
			profpeaks2<<-profpeaks2[order(profpeaks2[,"mean_mz"], decreasing = FALSE),,drop = FALSE]
			sort_by<<-"mean_mz";sort_by_decreasing <<- "FALSE"
		}
		if(isolate(input$filterProf_sort)=="mean m/z (decreasing)"){
			profpeaks2<<-profpeaks2[order(profpeaks2[,"mean_mz"], decreasing = TRUE),,drop = FALSE]
			sort_by<<-"mean_mz";sort_by_decreasing <<- "TRUE"
		}		
		if(isolate(input$filterProf_sort)=="mean RT (increasing)"){
			profpeaks2<<-profpeaks2[order(profpeaks2[,"mean_RT"], decreasing = FALSE),,drop = FALSE]
			sort_by<<-"mean_RT";sort_by_decreasing <<- "FALSE"
		}	
		if(isolate(input$filterProf_sort)=="mean RT (decreasing)"){
			profpeaks2<<-profpeaks2[order(profpeaks2[,"mean_RT"], decreasing = TRUE),,drop = FALSE]
			sort_by<<-"mean_RT";sort_by_decreasing <<- "TRUE"
		}
		if(isolate(input$filterProf_sort)=="minimum RT (decreasing)"){		
			profpeaks2<<-profpeaks2[order(profpeaks2[,"min_RT"], decreasing = TRUE),,drop = FALSE]	
			sort_by<<-"min_RT";sort_by_decreasing <<- "TRUE"			
		}	
		if(isolate(input$filterProf_sort)=="maximum RT (decreasing)"){		
			profpeaks2<<-profpeaks2[order(profpeaks2[,"max_RT" ], decreasing = TRUE),,drop = FALSE]	
			sort_by<<-"max_RT";sort_by_decreasing <<- "TRUE"			
		}
		if(isolate(input$filterProf_sort)=="maximum overall intensity (decreasing)"){
			profpeaks2<<-profpeaks2[order(profpeaks2[,"max_int"], decreasing = TRUE),,drop = FALSE]
			sort_by<<-"max_int";sort_by_decreasing <<- "TRUE"
		}
		if(isolate(input$filterProf_sort)=="maximum intensity in samples (decreasing, zeros removed)"){		
			profpeaks2<<-profpeaks2[order(profpeaks2[,"max_int_sample"],decreasing=TRUE),,drop = FALSE]	
			profpeaks2<<-profpeaks2[profpeaks2[,"max_int_sample"] != 0,,drop = FALSE]			
			sort_by<<-"max_int_sample";sort_by_decreasing <<- "TRUE"
		}
		if(isolate(input$filterProf_sort)=="maximum intensity in blanks/blinds (decreasing, zeros removed)"){		
			profpeaks2<<-profpeaks2[order(profpeaks2[,"max_int_blind"], decreasing = TRUE),, drop = FALSE]	
			profpeaks2<<-profpeaks2[profpeaks2[,"max_int_blind"] != 0,, drop = FALSE]			
			sort_by<<-"max_int_blind";sort_by_decreasing <<- "TRUE"
		}
		if(isolate(input$filterProf_sort)=="mean intensity (decreasing)"){
			profpeaks2<<-profpeaks2[order(profpeaks2[,"mean_int"], decreasing = TRUE),, drop = FALSE]
			sort_by<<-"mean_int";sort_by_decreasing <<- "TRUE"
		}
		if(isolate(input$filterProf_sort)=="mean intensity in samples (decreasing, zeros removed)"){		
			profpeaks2<<-profpeaks2[order(profpeaks2[,"mean_int_sample"],decreasing=TRUE),,drop = FALSE]	
			profpeaks2<<-profpeaks2[profpeaks2[,"mean_int_sample"] != 0,,drop = FALSE]
			sort_by<<-"mean_int_sample";sort_by_decreasing <<- "TRUE"
		}
		if(isolate(input$filterProf_sort)=="mean intensity in blanks/blinds (decreasing, zeros removed)"){		
			profpeaks2<<-profpeaks2[order(profpeaks2[,"mean_int_blind"],decreasing=TRUE),, drop = FALSE]	
			profpeaks2<<-profpeaks2[profpeaks2[,"mean_int_blind"] != 0,,drop = FALSE]
			sort_by<<-"mean_int_blind";sort_by_decreasing <<- "TRUE"
		}
		if(isolate(input$filterProf_sort) == "past & current trend intensity (decreasing)"){		
			profpeaks2<<-profpeaks2[order(profpeaks2[,"deltaint_global"],decreasing=TRUE),, drop = FALSE]
			profpeaks2<<-profpeaks2[profpeaks2[,"deltaint_global"] != 0,, drop = FALSE]
			sort_by<<-"deltaint_global"; sort_by_decreasing <<- "TRUE"
		}
		if(isolate(input$filterProf_sort)=="current trend intensity (decreasing)"){
			profpeaks2 <<- profpeaks2[order(profpeaks2[,"deltaint_newest"], decreasing = TRUE),, drop = FALSE]
			profpeaks2 <<- profpeaks2[profpeaks2[,"deltaint_newest"] != 0,,drop = FALSE]
			sort_by <<- "deltaint_newest"; sort_by_decreasing <<- "TRUE"
		}
		if(isolate(input$filterProf_sort)=="total peak number (decreasing)"){
			profpeaks2<<-profpeaks2[order(profpeaks2[,"number_peaks_total"],profpeaks2[,"max_int"], decreasing = TRUE),, drop = FALSE]
			sort_by<<-"number_peaks_total";sort_by_decreasing <<- "TRUE"
		}			
		if(isolate(input$filterProf_sort)=="peak number in samples (decreasing, zeros removed)"){		
			profpeaks2<<-profpeaks2[order(profpeaks2[,"number_peaks_sample"], decreasing = TRUE),, drop = FALSE]	
			profpeaks2<<-profpeaks2[profpeaks2[,"number_peaks_sample"]!=0,,drop = FALSE]		
			sort_by<<-"number_peaks_sample";sort_by_decreasing <<- "TRUE"
		}
		if(isolate(input$filterProf_sort)=="peak number in blanks/blinds (decreasing, zeros removed)"){		
			profpeaks2<<-profpeaks2[order(profpeaks2[,"number_peaks_blind"],decreasing = TRUE),,drop = FALSE]	
			profpeaks2<<-profpeaks2[profpeaks2[,"number_peaks_blind"]!=0,,drop = FALSE]	
			sort_by<<-"number_peaks_blind";sort_by_decreasing <<- "TRUE"
		}
		if(isolate(input$filterProf_sort)=="mass defect (increasing)"){		
			profpeaks2<<-profpeaks2[order(profpeaks2[,"Mass defect"], decreasing = FALSE),,drop = FALSE]	
			sort_by<<-"Mass defect";sort_by_decreasing <<- "FALSE"				
		}
		if(isolate(input$filterProf_sort)=="mass defect (decreasing)"){		
			profpeaks2<<-profpeaks2[order(profpeaks2[,"Mass defect"], decreasing = TRUE),,drop = FALSE]		
			sort_by<<-"Mass defect";sort_by_decreasing <<- "TRUE"
		}		
		if(isolate(input$filterProf_sort)=="median sample above blind intensity (decreasing)"){		
			profpeaks2<<-profpeaks2[order(profpeaks2[,"in_blind?"], decreasing = TRUE),,drop = FALSE]		
			sort_by<<-"in_blind?";sort_by_decreasing <<- "TRUE"
		}
		###################################################################################################		
		if(
			(isolate(input$filterProf_components)=="TRUE") &
			(logfile$workflow[names(logfile$workflow)=="components_profiles"]=="yes") &
			(logfile$summary[logfile$summary[,1]=="components_profiles",2]=="TRUE") &
			(any(objects(envir=as.environment(".GlobalEnv"))=="links_profiles"))
		){
			keep_IDs<-enviMass::analyseE_links_profiles(
				profileList_index=profpeaks2, 
				links_profiles, 
				sort_what=sort_by, 
				sort_decreasing=sort_by_decreasing,
				use_profile = NULL, 
				with_bar = FALSE, 
				return_excl=FALSE		
			)
			if(length(keep_IDs) & all(!is.na(keep_IDs))){
				profpeaks2<<-profpeaks2[match(keep_IDs,profpeaks2[,"profile_ID"]),,drop=FALSE]
			}else{
				stop("\nDebug server_obs_res_meas.r @ 1")
			}
		}
		###################################################################################################			
		if(dim(profpeaks2)[1]>0){
			atit1<-sum(profpeaks2[,"number_peaks_total"]) 
			output$atprof1<-renderText({ atit1 })
			atit2<-dim(profpeaks2)[1]
			#output$atprof2<-renderText({ atit2 })
			atit3<-sum(profpeaks2[,"number_peaks_blind"]>0)
			output$atprof3<-renderText({ atit3 })
			atit4<-sum(profpeaks2[,"deltaint_global"]>0)
			output$atprof4<-renderText({ atit4 })
			atit5<-sum(profpeaks2[,"deltaint_newest"]>0)
			output$atprof5<-renderText({ atit5 })
			# intensity histogram
			path=file.path(logfile[[1]],"pics","profilehisto.png");
                png(filename = path, bg = "white", width = 600);
                plot_profiles_intensity_histograms(	mean_intensities=profpeaks2[,"mean_int"],
                                                    max_intensities=profpeaks2[,"max_int"],
                                                    past_incidents=profpeaks2[,"deltaint_global"],
                                                    current_incidents=profpeaks2[,"deltaint_newest"]);
            dev.off();
			expr6<-list(src=file.path(logfile[[1]],"pics","profilehisto.png"));
			output$profilehisto<-renderImage(expr6, deleteFile = FALSE);
			# generate output table ############################################
			if( (length(profpeaks2[,1])>isolate(input$filterProf_count))  &  !is.na(isolate(input$filterProf_count)) ){
				profpeaks2<<-profpeaks2[1:isolate(input$filterProf_count),,drop=FALSE]
			}
			profpeaks3<<-as.data.frame(profpeaks2)
			profpeaks3[,"mean_mz" ]<<-round(profpeaks3[,"mean_mz"],digits=6)
			profpeaks3[,"mean_RT" ]<<-round(profpeaks3[,"mean_RT"],digits=1)
			profpeaks3[,"min_RT" ]<<-round(profpeaks3[,"min_RT"],digits=1)
			profpeaks3[,"max_RT" ]<<-round(profpeaks3[,"max_RT"],digits=1)
			profpeaks3[,"mean_int"]<<-round(log10(profpeaks3[,"mean_int"]),digits=3)
			profpeaks3[,"max_int"]<<-round(log10(profpeaks3[,"max_int"]),digits=3)
			profpeaks3[,"in_blind?"]<<-as.integer(profpeaks3[,"in_blind?"])
			profpeaks3[,"above_blind?"]<<-round(profpeaks3[,"above_blind?"],digits=4)
			profpeaks3[,"var_mz"]<<-round(log10(profpeaks3[,"var_mz"]),digits=3)
			profpeaks3[,"profile_ID"]<<-as.character(profpeaks3[,"profile_ID"])
			profpeaks3[,"number_peaks_total"]<<-as.integer(profpeaks3[,"number_peaks_total"])
			profpeaks3[,"deltaint_global"]<<-round(log10(profpeaks3[,"deltaint_global"]),digits=2)
			profpeaks3[,"deltaint_newest"]<<-round(log10(profpeaks3[,"deltaint_newest"]),digits=2)
			profpeaks3<<-profpeaks3[,c("profile_ID","number_peaks_total","deltaint_global","deltaint_newest","mean_mz",
				"mean_int","mean_RT","mean_RT","max_int","in_blind?","above_blind?","var_mz","min_RT","max_RT")]
			profpeaks3[,8]<<-round((profpeaks3[,8]/60),digits=2)
			output$allproftable<-DT::renderDataTable(
				DT::datatable(profpeaks3,
					colnames=c("profile ID","number of peaks","log10 global trend intensity","log10 current trend intensity","mean m/z","log10 mean intensity", 
						"mean RT [s]","mean RT [min]","log10 maximum Intensity","in blind?","above blind?","log10 m/z variance","minimum RT [s]","maximum RT [s]"),
					rownames=FALSE,
					filter = 'top',
                    selection = list(mode = 'single', target = 'row'),
					extensions = c('Buttons','FixedHeader','ColReorder'),
					options = list(
						lengthMenu = list(c(25, 50, 100, 200, -1), list('25', '50', '100', '200', 'All')),
						fixedHeader = FALSE,
						ordering=T,
						dom = 'Blfrtip',
						buttons = c('excel', 'csv','colvis'),#buttons = c('excel', 'pdf', 'print', 'csv'),
						scrollX = TRUE,
						colReorder = TRUE
					)
				),
				server = TRUE
			)
			updateNumericInput(session,"profID",value = 0);
			updateNumericInput(session,"profentry",value = 0);
			return(as.character(atit2));
		}else{
			output$atprof1<-renderText({ "0" })
			#output$atprof2<-renderText({ "0" })
			output$atprof3<-renderText({ "0" })
			output$atprof4<-renderText({ "0" })
			output$atprof5<-renderText({ "0" })
			path=file.path(logfile[[1]],"pics","profilehisto.png");
			png(filename = path, bg = "white", width = 1100);
				plot.new()
				plot.window(xlim=c(0,1),ylim=c(0,1))
				text(0.5,0.5,labels="0 profiles for these filter settings",cex=1.8,col="red")
			dev.off()
			expr6<-list(src=file.path(logfile[[1]],"pics","profilehisto.png"))
			output$profilehisto<-renderImage(expr6, deleteFile = FALSE);
			output$allproftable<-renderText("No profiles left")
			updateNumericInput(session,"profID",value = 0)
			updateNumericInput(session,"profentry",value = 0)
			return("0")
		}
	}else{
		if( isolate(init$a)=="TRUE" ){
			cat("\n No profiles available\n")
			output$atprof1<-renderText({"0"}) # now used as reactive output
			#output$atprof2<-renderText({"0"})
			output$atprof3<-renderText({"0"})
			output$atprof4<-renderText({"0"})
			output$atprof5<-renderText({"0"})	
			output$allproftable<-renderText("No profiles available")
			path=file.path(logfile[[1]],"pics","profilehisto.png")
			png(filename = path, bg = "white", width = 1100)
				plot.new()
				plot.window(xlim=c(0,1),ylim=c(0,1))
				text(0.5,0.5,labels="0 profiles for these filter settings",cex=1.8,col="red")
			dev.off()
			expr6<-list(src=file.path(logfile[[1]],"pics","profilehisto.png"))
			output$profilehisto<-renderImage(expr6, deleteFile = FALSE)
			updateNumericInput(session,"profID",value = 0)
			updateNumericInput(session,"profentry",value = 0)
			return("0")
		}
	}
	if(any(objects()=="profileList")){stop("illegal profpeaks2 found, #4");}
	if(any(objects()=="profpeaks2")){stop("illegal profpeaks2 found, #4");}
	if(any(objects()=="profpeaks3")){stop("illegal profpeaks3 found, #4");}
})	
output$atprof2<-renderText(paste(maincalc6()))
output$prof_number<-renderText(paste(maincalc6())) 
##############################################################################

##############################################################################
# PROFILE FITERING - update results for individual profileIDs ################
##############################################################################
ranges_timeprofile <- reactiveValues(x = NULL, y = NULL)

observe({
    input$profID
	init$b
	cat("\n plotting profile _1")
    if( (isolate(init$a)=="TRUE") &  
		(!is.na(isolate(input$profID))) & 
		(isolate(input$profID)!=0) & 
		any(objects(envir=as.environment(".GlobalEnv"))=="profileList") & 
		any(objects(envir=as.environment(".GlobalEnv"))=="profpeaks2")
	){
		cat("\n plotting profile _2")
		if(any(objects()=="profileList")){stop("illegal profileList found, #5");}
		if(any(profileList[["index_prof"]][,"profile_ID"]==as.numeric(isolate(input$profID)))){
			cat("\n plotting profile with ID ");cat(as.numeric(isolate(input$profID)));
			##################################################################
			if(logfile$parameters$trend_blind=="yes"){
				blindsubtract<-TRUE
			}else{
				blindsubtract<-FALSE
			}
			lagit<-as.numeric(strsplit(logfile$parameters$trend_lags,",")[[1]])		
			if(isolate(input$prof_log)=="yes"){
				logscaled<-TRUE
			}else{
				logscaled<-FALSE
			}
			output$timeprofile <- renderPlot({			
				assign("peakTable",plotaprofile(
					profileList,
					profileID=as.numeric(isolate(input$profID)),
					logint=logscaled,
					blindsub=blindsubtract,
					blindfold=as.numeric(logfile$parameters$blind_threshold),
					lags=lagit,
					threshold=as.numeric(logfile$parameters$trend_threshold),
					ranges_x=ranges_timeprofile$x,
					ranges_y=ranges_timeprofile$y,
					),envir=as.environment(".GlobalEnv")
				)	
				if((!is.null(ranges_timeprofile$x))||(!is.null(ranges_timeprofile$y))){
					mtext("Now zoomed in",side=3,col="gray")
				}
			})			
			##################################################################
			prof_plot_IDs<<-as.numeric(isolate(input$profID))
			at_entry<-profileList[["index_prof"]][
				match(prof_plot_IDs,profileList[["index_prof"]][,"profile_ID"])
			,"links"]
			plot_similar_profiles<-FALSE
			if(at_entry!=0){
				if(length(links_profiles[[at_entry]]$group)>0){
					prof_plot_IDs<<-c(prof_plot_IDs,links_profiles[[at_entry]]$group)
					prof_plot_IDs<<-unique(prof_plot_IDs)
					plot_similar_profiles<-TRUE
				}
			}
			output$similar_profiles_plot <- renderPlot({	
				if(plot_similar_profiles){
						enviMass::plot_components(
							profileList=profileList,
							prof_IDs=prof_plot_IDs,
							links_profiles=links_profiles,
							what="profiles",
							xlim=FALSE,ylim=FALSE,await_input=FALSE,
							skipit=TRUE,
							min_peaks=NULL,
							norma=TRUE
						)	
				}else{
					plot.new()
					plot.window(xlim=c(0,1),ylim=c(0,1))
					text(0.5,0.5,labels="Nothing to plot - no other similar profiles available",cex=1.8,col="red")
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
			use_prof_IDs<-match(prof_plot_IDs,profileList[["index_prof"]][,"profile_ID"])
			use_prof_IDs<-use_prof_IDs[!is.na(use_prof_IDs)]
			similar_profiles_tab<-profileList[["index_prof"]][use_prof_IDs,c("profile_ID","mean_mz","mean_int","max_int_sample","mean_RT",
				"Mass defect","above_blind?","number_peaks_sample","number_peaks_blind","number_peaks_total"),drop=FALSE]
			similar_profiles_tab<-as.data.frame(similar_profiles_tab)
			similar_profiles_tab[,"mean_mz" ]<-format(similar_profiles_tab[,"mean_mz" ],digits=5)
			similar_profiles_tab[,"mean_int"]<-format(similar_profiles_tab[,"mean_int"],scientific=TRUE,digits=2)
			similar_profiles_tab[,"max_int_sample"]<-format(similar_profiles_tab[,"max_int_sample"],scientific=TRUE,digits=2)
			similar_profiles_tab[,"mean_RT"]<-format(similar_profiles_tab[,"mean_RT"],digits=1)
			similar_profiles_tab[,"Mass defect"]<-format(similar_profiles_tab[,"Mass defect"],digits=3)
			similar_profiles_tab[,"above_blind?"]<-format(similar_profiles_tab[,"above_blind?"],digits=3)
			similar_profiles_tab[,"number_peaks_sample"]<-as.integer(similar_profiles_tab[,"number_peaks_sample"])
			similar_profiles_tab[,"number_peaks_blind"]<-as.integer(similar_profiles_tab[,"number_peaks_blind"])
			similar_profiles_tab[,"number_peaks_total"]<-as.integer(similar_profiles_tab[,"number_peaks_total"])
			output$similar_profiles_table<-DT::renderDataTable({
				DT::datatable(similar_profiles_tab,
					colnames=c("Profile ID","Mean m/z","Mean int. overall","Max int. samples","Mean RT",
						"Mass defect","Median int. ratio blind vs. samples","Peak number samples","Peak number blinds","Peak number total"),
					rownames=FALSE,
					filter = 'top',
                    selection = list(mode = 'single', target = 'row'),
					extensions = c('Buttons','FixedHeader','ColReorder'),
					options = list(
						lengthMenu = c(15, 30, 50, 100),
						fixedHeader = FALSE,
						ordering=T,
						dom = 'Blfrtip',
						buttons = c('excel','colvis'),#buttons = c('excel', 'pdf', 'print', 'csv'),
						scrollX = TRUE,
						colReorder = TRUE
					)						
				)			
			});
			##################################################################
			output$oneproftable<-DT::renderDataTable(peakTable);
			updateNumericInput(session,"profpeakID",value = 0);
			path=file.path(logfile[[1]],"pics","massdens.png");
			png(filename = path, bg = "white", width = 550,height=200);			
				plot.new()
				plot.window(xlim=c(0,1),ylim=c(0,1))
				text(0.5,0.5,labels="Waiting...",cex=1.8,col="red")
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
# update results per profilepeak list entry index ############################
##############################################################################
observe({
    input$profentry
	init$b
    if(	(isolate(init$a)=="TRUE") &  
		!is.na(isolate(input$profentry)) & 
		(isolate(input$profentry)!=0) & 
		(any(objects(envir=as.environment(".GlobalEnv"))=="profileList")) & 
		any(objects(envir=as.environment(".GlobalEnv"))=="profpeaks2")
	){
		if( (isolate(input$profentry)<=length(profpeaks2[,1])) & 
			(isolate(input$profentry)>0) 
		){ 
				if( any( profileList[["index_prof"]][,"profile_ID"]==as.numeric(profpeaks2[isolate(input$profentry),"profile_ID"]) ) ){
					updateNumericInput(session,"profID",value = as.numeric(as.character(profpeaks2[isolate(input$profentry),"profile_ID"])))				
				}
		}else{
			output$timeprofile <- renderPlot({	
				plot.new()
				plot.window(xlim=c(0,1),ylim=c(0,1))
				text(0.5,0.5,labels="Invalid list entry",cex=1.8,col="red")
			})			
			output$oneproftable<-renderText("")
		}
	}else{
		if(isolate(init$a)=="TRUE"){
			output$timeprofile <- renderPlot({	
				plot.new()
				plot.window(xlim=c(0,1),ylim=c(0,1))
				text(0.5,0.5,labels="Nothing to plot - invalid ID",cex=1.8,col="red")
			})			
			output$oneproftable<-renderText("")
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
ranges_profiles <- reactiveValues(RTchrom=FALSE, intchrom=FALSE)

observe({ # seconds <-> minutes switch when zoomed ###########################
	input$profile_EIC_time
	if(isolate(init$a)=="TRUE" & isolate(ranges_profiles$RTchrom[1]!=FALSE)){
		if(isolate(input$profile_EIC_time)=="minutes"){
			isolate(ranges_profiles$RTchrom<-(ranges_profiles$RTchrom/60))
		}
		if(isolate(input$profile_EIC_time)=="seconds"){
			isolate(ranges_profiles$RTchrom<-(ranges_profiles$RTchrom*60))
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
maincalc5<-reactive({
	input$dens_mass
	if(	(isolate(input$dens_mass)) &
		(isolate(init$a)=="TRUE") & 
		(isolate(input$profID)!=0) & 
		(any(objects(envir=as.environment(".GlobalEnv"))=="profileList")) 
	){
		if(any(profileList[[7]][,4]==as.numeric(isolate(input$profID)))){
			######################################################################
			cat("\n kernel density ...")
			if(isolate(input$use_weight)=="yes"){use_weights<-TRUE}else{use_weights<-FALSE}
			path=file.path(logfile[[1]],"pics","massdens.png");
			png(filename = path, bg = "white", width = 550,height=400);	
				getmass<-mass_dens(
						profileList,
						profileID=as.numeric(isolate(input$profID)),
						bootstrap=TRUE,
						boot_size=as.numeric(isolate(input$boot_size)),
						use_weights
				)
			dev.off();
			expr_massdens<-list(src=file.path(logfile[[1]],"pics","massdens.png"));
			output$massdens<-renderImage(expr_massdens, deleteFile = FALSE);	
			path=file.path(logfile[[1]],"pics","mass_int.png");
			png(filename = path, bg = "white", width = 550,height=300);	
					mass_int(
						profileList,
						profileID=as.numeric(isolate(input$profID))
					)
			dev.off();
			expr_mass_int<-list(src=file.path(logfile[[1]],"pics","mass_int.png"));
			output$massint<-renderImage(expr_mass_int, deleteFile = FALSE);	
			return(getmass);
			######################################################################
		}else{
			return("Not available");
		}
	}else{
		######################################################################
		cat("\n kernel density ...")
		path=file.path(logfile[[1]],"pics","massdens.png");
		png(filename = path, bg = "white", width = 550,height=200);			
			plot.new()
			plot.window(xlim=c(0,1),ylim=c(0,1))
			text(0.5,0.5,labels="Not available",cex=1.8,col="red")
		dev.off();
		expr_massdens<-list(src=file.path(logfile[[1]],"pics","massdens.png"));
		output$massdens<-renderImage(expr_massdens, deleteFile = FALSE);		
		return("...");
		######################################################################
	}
})
output$prof_mass<-renderText(paste(maincalc5())) 
##############################################################################

##############################################################################
# Observe project reset buttons ##############################################
##############################################################################
observe({
    input$reset_1
    if( (isolate(init$a)=="TRUE") & isolate(input$reset_1) ){

		if(any(ls()=="logfile")){stop(paste("\n illegal logfile detected in server_obs_res_mean.r #1"))}
		logfile$Tasks_to_redo<<-replace(logfile$Tasks_to_redo,-1,TRUE)
		#logfile$Tasks_to_redo<<-replace(logfile$Tasks_to_redo,1,FALSE)
		measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
		measurements[,c("qc","recal","align","norm", "LOD","isotopologues","adducts","homologues","EIC_correlation","blind","components_files")]<-"FALSE"
		write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
		createAlert(session,anchorId = "reset", alertId="reset1", title = NULL, content="Project reset w/o peak picking",style = "warning",append=FALSE,dismiss=TRUE)
		save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
		cat("\nReset without peak picking \n")
	}
})
observe({
    input$reset_2
    if( (isolate(init$a)=="TRUE") & isolate(input$reset_2) ){
		if(any(ls()=="logfile")){stop(paste("\n illegal logfile detected in server_obs_res_mean.r #1"))}		
		logfile$Tasks_to_redo<<-replace(logfile$Tasks_to_redo,,TRUE)
		measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
		if(TRUE){
			measurements[,c("peakpicking","qc","recal","align","norm", "LOD","isotopologues","adducts","homologues","EIC_correlation","blind","components_files")]<-"FALSE"
		}
		write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
		rm(measurements)
		# delete all peaklists
		those<-list.files(file.path(logfile$project_folder,"peaklist"))
		if(length(those)>0){
			for(i in 1:length(those)){
				file.remove(file.path(logfile$project_folder,"peaklist",those[i]))
			}
		}
		createAlert(session,anchorId = "reset", alertId="reset2", title = NULL, content="Project reset",style = "warning",append=FALSE,dismiss=TRUE)
		save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
		cat("\nTotal reset \n")
	}
})
##############################################################################

##############################################################################
# Observe resolution data sets ###############################################
##############################################################################
observe({
    input$resolution
	#init$a
    if( (isolate(init$a)=="TRUE") & (isolate(input$resolution)!="none")  & (isolate(input$resolution)!="") ){
		#cat(input$resolution);cat("\n")
		path=file.path(logfile[[1]],"pics","resolution")
		png(filename = path, bg = "white")
			that<-resolution_list[names(resolution_list) == as.character(input$resolution)][[1]]
			plot(that[,1],that[,2],pch=19,cex=0.5,xlab="m/z",ylab="Resolution")		
		dev.off()
		exprres<-list(src=file.path(logfile[[1]],"pics","resolution"))
		output$plot_resolution<-renderImage(exprres, deleteFile = FALSE)	
	}
})
##############################################################################










