ranges_int_norm_IS_pos <- reactiveValues(xlim = FALSE)
ranges_int_norm_IS_neg <- reactiveValues(xlim = FALSE)



#####################################################################################
# -> POS
observeEvent(input$int_norm_ISTD_pos_median_dblclick, { 
          brush <- isolate(input$int_norm_ISTD_pos_median_brush)
          if (!is.null(brush)) {
            cat("\n Zoom in_1")
            isolate(ranges_int_norm_IS_pos$xlim <- c(brush$xmin, brush$xmax))
          } else {
            cat("\n Zoom out full_1")
            isolate(ranges_int_norm_IS_pos$xlim <- FALSE)
          }
          #refresh_homol$a<-(refresh_homol$a+1) # valid in both cases
})
observeEvent(input$int_norm_ISTD_pos_median_click, { # NOTE: brushing already triggers a click -> use brush with delay=0, which embeds the slower click
          cat("\n Zoom out part_1_a")
          brush <- isolate(input$int_norm_ISTD_pos_median_brush)
          if (is.null(brush)) {
            cat("\n Zoom out part_1_b")
            if(isolate(ranges_int_norm_IS_pos$xlim[1]) != FALSE){
              old_range_mass<-abs(isolate(ranges_int_norm_IS_pos$xlim[2] - ranges_int_norm_IS_pos$xlim[1]))
              isolate(ranges_int_norm_IS_pos$xlim[1] <- ranges_int_norm_IS_pos$xlim[1] - .3 * old_range_mass)
              isolate(ranges_int_norm_IS_pos$xlim[2] <- ranges_int_norm_IS_pos$xlim[2] + .3 * old_range_mass)
            }
            #refresh_homol$a<-(refresh_homol$a+1)
          }else{
            cat("\n Doing hover - and nothing")
          }   
})

observeEvent(input$int_norm_ISTD_pos_counts_dblclick, { 
          brush <- isolate(input$int_norm_ISTD_pos_counts_brush)
          if (!is.null(brush)) {
            cat("\n Zoom in_1")
            isolate(ranges_int_norm_IS_pos$xlim <- c(brush$xmin, brush$xmax))
          } else {
            cat("\n Zoom out full_1")
            isolate(ranges_int_norm_IS_pos$xlim <- FALSE)
          }
          #refresh_homol$a<-(refresh_homol$a+1) # valid in both cases
})
observeEvent(input$int_norm_ISTD_pos_counts_click, { # NOTE: brushing already triggers a click -> use brush with delay=0, which embeds the slower click
          cat("\n Zoom out part_1_a")
          brush <- isolate(input$int_norm_ISTD_pos_counts_brush)
          if (is.null(brush)) {
            cat("\n Zoom out part_1_b")
            if(isolate(ranges_int_norm_IS_pos$xlim[1]) != FALSE){
              old_range_mass<-abs(isolate(ranges_int_norm_IS_pos$xlim[2] - ranges_int_norm_IS_pos$xlim[1]))
              isolate(ranges_int_norm_IS_pos$xlim[1] <- ranges_int_norm_IS_pos$xlim[1] - .3 * old_range_mass)
              isolate(ranges_int_norm_IS_pos$xlim[2] <- ranges_int_norm_IS_pos$xlim[2] + .3* old_range_mass)
            }
            #refresh_homol$a<-(refresh_homol$a+1)
          }else{
            cat("\n Doing hover - and nothing")
          }   
})
# -> NEG
observeEvent(input$int_norm_ISTD_neg_median_dblclick, { 
          brush <- isolate(input$int_norm_ISTD_neg_median_brush)
          if (!is.null(brush)) {
            cat("\n Zoom in_1")
            isolate(ranges_int_norm_IS_neg$xlim <- c(brush$xmin, brush$xmax))
          } else {
            cat("\n Zoom out full_1")
            isolate(ranges_int_norm_IS_neg$xlim <- FALSE)
          }
          #refresh_homol$a<-(refresh_homol$a+1) # valid in both cases
})
observeEvent(input$int_norm_ISTD_neg_median_click, { # NOTE: brushing already triggers a click -> use brush with delay=0, which embeds the slower click
          cat("\n Zoom out part_1_a")
          brush <- isolate(input$int_norm_ISTD_neg_median_brush)
          if (is.null(brush)) {
            cat("\n Zoom out part_1_b")
            if(isolate(ranges_int_norm_IS_neg$xlim[1]) != FALSE){
              old_range_mass<-abs(isolate(ranges_int_norm_IS_neg$xlim[2] - ranges_int_norm_IS_neg$xlim[1]))
              isolate(ranges_int_norm_IS_neg$xlim[1] <- ranges_int_norm_IS_neg$xlim[1] - .3 * old_range_mass)
              isolate(ranges_int_norm_IS_neg$xlim[2] <- ranges_int_norm_IS_neg$xlim[2] + .3 * old_range_mass)
            }
            #refresh_homol$a<-(refresh_homol$a+1)
          }else{
            cat("\n Doing hover - and nothing")
          }   
})

observeEvent(input$int_norm_ISTD_neg_counts_dblclick, { 
          brush <- isolate(input$int_norm_ISTD_neg_counts_brush)
          if (!is.null(brush)) {
            cat("\n Zoom in_1")
            isolate(ranges_int_norm_IS_neg$xlim <- c(brush$xmin, brush$xmax))
          } else {
            cat("\n Zoom out full_1")
            isolate(ranges_int_norm_IS_neg$xlim <- FALSE)
          }
          #refresh_homol$a<-(refresh_homol$a+1) # valid in both cases
})
observeEvent(input$int_norm_ISTD_neg_counts_click, { # NOTE: brushing already triggers a click -> use brush with delay=0, which embeds the slower click
          cat("\n Zoom out part_1_a")
          brush <- isolate(input$int_norm_ISTD_neg_counts_brush)
          if (is.null(brush)) {
            cat("\n Zoom out part_1_b")
            if(isolate(ranges_int_norm_IS_neg$xlim[1]) != FALSE){
              old_range_mass<-abs(isolate(ranges_int_norm_IS_neg$xlim[2] - ranges_int_norm_IS_neg$xlim[1]))
              isolate(ranges_int_norm_IS_neg$xlim[1] <- ranges_int_norm_IS_neg$xlim[1] - .3 * old_range_mass)
              isolate(ranges_int_norm_IS_neg$xlim[2] <- ranges_int_norm_IS_neg$xlim[2] + .3 * old_range_mass)
            }
            #refresh_homol$a<-(refresh_homol$a+1)
          }else{
            cat("\n Doing hover - and nothing")
          }   
})



#####################################################################################

#####################################################################################
observe({ # -> POS
	ranges_int_norm_IS_pos$xlim
	if(isolate(init$a)=="TRUE"){
		#############################################################################
		if(file.exists(file.path(as.character(logfile[[1]]), "results", "int_norm_ISTD_pos"))){
			if(!any(objects(envir = as.environment(".GlobalEnv")) == "int_norm_ISTD_pos")){
				load(file.path(as.character(logfile[[1]]), "results", "int_norm_ISTD_pos"), envir = as.environment(".GlobalEnv"))
			}
			output$int_norm_ISTD_pos_median <- renderPlot({   
				par(mar = c(.2, 4.5, .9, 8))
				enviMass:::plot_ISTD_norm(
					int_norm_ISTD = int_norm_ISTD_pos,
					logfile = logfile,
					xlim = isolate(ranges_int_norm_IS_pos$xlim),
					what = "normalization"
				)
			},res = 100) 
			output$int_norm_ISTD_pos_counts <- renderPlot({   
				par(mar = c(4.5, 4.5, .9, 8))
				enviMass:::plot_ISTD_norm(
					int_norm_ISTD = int_norm_ISTD_pos,
					logfile = logfile,
					xlim = isolate(ranges_int_norm_IS_pos$xlim),
					what = "counts"
				)
			},res = 100) 					
		}else{
			output$int_norm_ISTD_pos_median <- renderPlot({ 
				plot.new()
			})
			output$int_norm_ISTD_pos_counts <- renderPlot({ 
				plot.new()
			})		
		}
		#############################################################################
	}
})
observe({ # -> NEG
	ranges_int_norm_IS_neg$xlim
	if(isolate(init$a)=="TRUE"){
		#############################################################################
		if(file.exists(file.path(as.character(logfile[[1]]), "results", "int_norm_ISTD_neg"))){
			if(!any(objects(envir = as.environment(".GlobalEnv")) == "int_norm_ISTD_neg")){
				load(file.path(as.character(logfile[[1]]), "results", "int_norm_ISTD_neg"), envir = as.environment(".GlobalEnv"))
			}
			output$int_norm_ISTD_neg_median <- renderPlot({   
				par(mar = c(.2, 4.5, .9, 8))
				enviMass:::plot_ISTD_norm(
					int_norm_ISTD = int_norm_ISTD_neg,
					logfile = logfile,
					xlim = isolate(ranges_int_norm_IS_neg$xlim),
					what = "normalization"
				)
			},res = 100) 
			output$int_norm_ISTD_neg_counts <- renderPlot({   
				par(mar = c(4.5, 4.5, .9, 8))
				enviMass:::plot_ISTD_norm(
					int_norm_ISTD = int_norm_ISTD_neg,
					logfile = logfile,
					xlim = isolate(ranges_int_norm_IS_neg$xlim),
					what = "counts"
				)
			},res = 100) 					
		}else{
			output$int_norm_ISTD_neg_median <- renderPlot({ 
				plot.new()
			})
			output$int_norm_ISTD_neg_counts <- renderPlot({ 
				plot.new()
			})		
		}
		#############################################################################
	}
})
#####################################################################################

