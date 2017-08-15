

  
  runApp(list(
    ##################################################################################
    ui = bootstrapPage(
			radioButtons("comp_hide", label="Show all component details?", choices = c("minimize"="minimize","maximize"="maximize"), selected = "maximize",
			  inline = TRUE, width = NULL, choiceNames = NULL, choiceValues = NULL),
			conditionalPanel(			
				#condition = "input.comp_hide == 'maximize'",
				condition = "1==2",	
				radioButtons("comp_hide2", label="Show all component details?", choices = c("minimize"="minimize","maximize"="maximize"), selected = "maximize",
				  inline = TRUE, width = NULL, choiceNames = NULL, choiceValues = NULL),				
				#conditionalPanel(			
				#	condition = "input.comp_hide2 == 'maximize'",					
					######################################################################
					plotOutput("homol_counts",
							click = "homol_counts_click",
							dblclick = "homol_counts_dblclick",
							#hover = "homol_counts_hover",
							brush = brushOpts(
							  id = "homol_counts_brush",
							  direction = c("x"),
							  resetOnNew = TRUE,
							  delay = 0
							),                
							height = "250px"
					),        
					######################################################################
					plotOutput("homol_RT",
							click = "homol_RT_click",
							dblclick = "homol_RT_dblclick",
							#hover = "homol_counts_hover",
							brush = brushOpts(
							  id = "homol_RT_brush",
							  direction = c("x"),
							  resetOnNew = TRUE,
							  delay = 0
							),                
							height = "250px"
					)
				#)
			),
                ######################################################################
                plotOutput("homol_plot",
                        click = "homol_plot_click",
                        dblclick = "homol_plot_dblclick",
                        #hover = "homol_plot_hover",
                        brush = brushOpts(
                          id = "homol_plot_brush",
                          resetOnNew = TRUE,
                          delay = 0
                        ),                
                        height = "700px"
                ), 
                ######################################################################
                bsCollapsePanel("Series table", 
                        tabsetPanel(  
                          tabPanel("Peaks in series ->",
                            HTML('<hr noshade="noshade" />'),
                            DT::dataTableOutput('homol_series_peaks')
                          ),
                          tabPanel("Series",
                            HTML('<hr noshade="noshade" />'),                
                            conditionalPanel(     
                              condition = "input.homol_series_table_rows_selected > 0",              
                                plotOutput("homol_chromat",
                                    height = "230px",
                                    click = "homol_chromat_click",
                                    dblclick = "homol_chromat_dblclick",
                                    #hover = "homol_chromat_hover",
                                    brush = brushOpts(
                                      id = "homol_chromat_brush",
                                      direction = c("x"),
                                      resetOnNew = TRUE,
                                      delay = 0
                                    )
                                ),
                                HTML('<hr noshade="noshade" />')
                            ), 
                            DT::dataTableOutput('homol_series_table')
                          )
                        ) 
                )
    ),
    ##################################################################################
    server = function(input, output) {  
	  at<-"1"
      load(file.path(logfile[[1]],"results","componentization","homologues",paste("full",at,sep="_")),envir=as.environment(".GlobalEnv"))
      load(file.path(logfile[[1]],"MSlist",at),envir=as.environment(".GlobalEnv"))      
    
      ranges_homol <- reactiveValues(mass = FALSE, RT = FALSE, massD = FALSE, dmass = FALSE, dRT = FALSE)
      refresh_homol <- reactiveValues()
      refresh_homol$a <- 0
      refresh_homol$b <- 0
      refresh_homol$c <- 0
      ################################################################################
      observe({
        refresh_homol$a
        cat("\n IN REFRESH")
                # filter segments to plot
                plot_those<<-enviMass:::filter_segments(
                  homol,
                  masslim = isolate(ranges_homol$mass),
                  RTlim = isolate(ranges_homol$RT),
                  massDlim = isolate(ranges_homol$massD),
                  dmasslim = isolate(ranges_homol$dmass),
                  dRTlim = isolate(ranges_homol$dRT)
                )
                sum_homol<-sum(plot_those)
                if(sum_homol>1000){
                  omit_theta<<-2
                  if(sum_homol>5000){omit_theta<<-10}
                  if(sum_homol>20000){omit_theta<<-15}    
                  if(sum_homol>100000){omit_theta<<-20}  
                  if(sum_homol>200000){omit_theta<<-30}           
                }else{
                  omit_theta<<-FALSE
                }
                # output homol. series plot ######################################
                output$homol_plot <- renderPlot({   
                  par(mar=c(5,4.5,.7,.5))
                  enviMass:::plothomol(homol,
                    xlim = isolate(ranges_homol$mass), ylim = isolate(ranges_homol$RT), 
                    dmasslim = isolate(ranges_homol$dmass), dRTlim = isolate(ranges_homol$dRT),
                    plot_what="mz_RT", omit_theta=omit_theta, plot_those = plot_those
                  );
                },res=100)  
                # output homol. mass difference counts plot #####################
                output$homol_counts <- renderPlot({   
                  par(mar=c(5,4.5,.7,.5))
                  enviMass:::plothomol(homol,
                    xlim = isolate(ranges_homol$mass), ylim = isolate(ranges_homol$RT), 
                    dmasslim = isolate(ranges_homol$dmass), dRTlim = isolate(ranges_homol$dRT),
                    plot_what="dmass", omit_theta=FALSE, plot_those = plot_those
                  );
                },res=100)          
                # output homol. RT difference vs mass scatter plot ##############
                output$homol_RT <- renderPlot({   
                  par(mar=c(5,4.5,.7,.5))
                  enviMass:::plothomol(homol,
                    xlim = isolate(ranges_homol$mass), ylim = isolate(ranges_homol$RT), 
                    dmasslim = isolate(ranges_homol$dmass), dRTlim = isolate(ranges_homol$dRT),
                    plot_what="dRT", omit_theta=FALSE, plot_those = plot_those
                  );
                },res=100)  
                # output tables ##################################################
                use_homol_peaks<-match(unique(c(homol[["homol_peaks_relat"]][plot_those,1],homol[["homol_peaks_relat"]][plot_those,2])),homol[["Peaks in homologue series"]][,"peak ID"])
                output$homol_series_peaks <- DT::renderDataTable(
                  DT::datatable(
                    dat<-data.frame(
                      I(as.character(homol[["Peaks in homologue series"]][use_homol_peaks,c("peak ID")])),
                      round(homol[["Peaks in homologue series"]][use_homol_peaks,c("mz")],digits=5),
                      round(log10(homol[["Peaks in homologue series"]][use_homol_peaks,c("intensity")]),digits=4),
                      round(homol[["Peaks in homologue series"]][use_homol_peaks,c("RT")],digits=1),
                      round((homol[["Peaks in homologue series"]][use_homol_peaks,c("RT")]/60),digits=1),        
                      I(homol[["Peaks in homologue series"]][use_homol_peaks,c("Targets")]),
                      I(homol[["Peaks in homologue series"]][use_homol_peaks,c("ISTDs")]),
                      I(homol[["Peaks in homologue series"]][use_homol_peaks,c("HS IDs")])             
                    ),
                    filter = 'top',
                    colnames=c("Peak ID","m/z","log10 intens.","RT [s]","RT [min]","Target matches","ISTD matches","Series ID(s)"),
                    rownames=FALSE,
                    selection = list(mode = 'single', target = 'row')
                  )
                )
                # output homol. series table #####################################
                output$homol_series_table <- DT::renderDataTable(
                  DT::datatable(
                    data.frame(
                      I(as.character(homol[["Homologue Series"]][,"HS IDs"])),
                      I(as.character(homol[["Homologue Series"]][,"peak IDs"])),
                      round(homol[["Homologue Series"]][,"m/z increment"],digits=5),
                      round(homol[["Homologue Series"]][,"RT increment"],digits=1),
                      round(log10(homol[["Homologue Series"]][,"max int."]),digits=4)
                    ),
                    filter = 'top',
                    colnames=c("Series ID","Peak IDs","m/z difference","RT difference [s]","Max log10 int."),
                    rownames=FALSE,
                    selection = list(mode = 'single', target = 'row')
                  )
                )
                ##################################################################
      })           
      ################################################################################ 
      observe({
          refresh_homol$b
          cat("\n IN REFRESH_2")
                # filter segments to plot
                plot_those<<-enviMass:::filter_segments(
                  homol,
                  masslim = isolate(ranges_homol$mass),
                  RTlim = isolate(ranges_homol$RT),
                  massDlim = isolate(ranges_homol$massD),
                  dmasslim = isolate(ranges_homol$dmass),
                  dRTlim = isolate(ranges_homol$dRT)
                )
                sum_homol<-sum(plot_those)
                if(sum_homol>1000){
                  omit_theta<<-2
                  if(sum_homol>5000){omit_theta<<-10}
                  if(sum_homol>20000){omit_theta<<-15}    
                  if(sum_homol>100000){omit_theta<<-20}  
                  if(sum_homol>200000){omit_theta<<-30}           
                }else{
                  omit_theta<<-FALSE
                }
                # output homol. series plot ######################################
                output$homol_plot <- renderPlot({   
                  par(mar=c(5,4.5,.7,.5))
                  enviMass:::plothomol(homol,
                    xlim = isolate(ranges_homol$mass), ylim = isolate(ranges_homol$RT), 
                    dmasslim = isolate(ranges_homol$dmass), dRTlim = isolate(ranges_homol$dRT),
                    plot_what="mz_RT", omit_theta=omit_theta, plot_those = plot_those
                  );
                },res=100) 
                # output tables ##################################################
                use_homol_peaks<-match(unique(c(homol[["homol_peaks_relat"]][plot_those,1],homol[["homol_peaks_relat"]][plot_those,2])),homol[["Peaks in homologue series"]][,"peak ID"])
                output$homol_series_peaks <- DT::renderDataTable(
                  DT::datatable(
                    dat<-data.frame(
                      I(as.character(homol[["Peaks in homologue series"]][use_homol_peaks,c("peak ID")])),
                      round(homol[["Peaks in homologue series"]][use_homol_peaks,c("mz")],digits=5),
                      round(log10(homol[["Peaks in homologue series"]][use_homol_peaks,c("intensity")]),digits=4),
                      round(homol[["Peaks in homologue series"]][use_homol_peaks,c("RT")],digits=1),
                      round((homol[["Peaks in homologue series"]][use_homol_peaks,c("RT")]/60),digits=1),        
                      I(homol[["Peaks in homologue series"]][use_homol_peaks,c("Targets")]),
                      I(homol[["Peaks in homologue series"]][use_homol_peaks,c("ISTDs")]),
                      I(homol[["Peaks in homologue series"]][use_homol_peaks,c("HS IDs")])             
                    ),
                    filter = 'top',
                    colnames=c("Peak ID","m/z","log10 intens.","RT [s]","RT [min]","Target matches","ISTD matches","Series ID(s)"),
                    rownames=FALSE,
                    selection = list(mode = 'single', target = 'row')
                  )
               )
      })
      ################################################################################ 
      observe({      
          s1<-input$homol_series_peaks_rows_selected
          if(length(s1)){
            if(s1>=1){
              print(s1);
              # output homol. series table ##################################### 
              use_homol_peaks<<-match(unique(c(homol[["homol_peaks_relat"]][plot_those,1],homol[["homol_peaks_relat"]][plot_those,2])),homol[["Peaks in homologue series"]][,"peak ID"])
              these_series<-as.numeric(strsplit(homol[["Peaks in homologue series"]][use_homol_peaks[s1],c("HS IDs")],"/")[[1]])
              print(these_series)
              use_these_series<-match(these_series,homol[["Homologue Series"]][,"HS IDs"])
              output$homol_series_table <- DT::renderDataTable(
                  DT::datatable(
                    data.frame(
                      I(as.character(homol[["Homologue Series"]][,"HS IDs"])),
                      I(as.character(homol[["Homologue Series"]][,"peak IDs"])),
                      round(homol[["Homologue Series"]][,"m/z increment"],digits=5),
                      round(homol[["Homologue Series"]][,"RT increment"],digits=1),
                      round(log10(homol[["Homologue Series"]][,"max int."]),digits=4)
                    ),
                    filter = 'top',
                    colnames=c("Series ID","Peak IDs","m/z difference","RT difference [s]","Max log10 int."),
                    rownames=FALSE,
                    selection = list(mode = 'single', target = 'row')
                  )
              )
              # output homol. series plot ######################################
              cat("\n for this peak: ");print(use_homol_peaks[s1]);
              output$homol_plot <- renderPlot({   
                par(mar=c(5,4.5,.7,.5))
                enviMass:::plothomol(homol,
                  xlim = isolate(ranges_homol$mass), ylim = isolate(ranges_homol$RT), 
                  dmasslim = isolate(ranges_homol$dmass), dRTlim = isolate(ranges_homol$dRT),
                  plot_what="mz_RT", omit_theta=omit_theta, plot_those = plot_those,
                  emph_point=homol[["Peaks in homologue series"]][use_homol_peaks[s1],"peak ID"],
                  emph_series=these_series
                );
              },res=100)
            }
          }else{
            cat("\n ALL DESELECTED_1")
            # output homol. series table #####################################    
            output$homol_series_table <- DT::renderDataTable(
                  DT::datatable(
                    data.frame(
                      I(as.character(homol[["Homologue Series"]][,"HS IDs"])),
                      I(as.character(homol[["Homologue Series"]][,"peak IDs"])),
                      round(homol[["Homologue Series"]][,"m/z increment"],digits=5),
                      round(homol[["Homologue Series"]][,"RT increment"],digits=1),
                      round(log10(homol[["Homologue Series"]][,"max int."]),digits=4)
                    ),
                    filter = 'top',
                    colnames=c("Series ID","Peak IDs","m/z difference","RT difference [s]","Max log10 int."),
                    rownames=FALSE,
                    selection = list(mode = 'single', target = 'row')
                  )
               )
            # output homol. series plot ######################################
            output$homol_plot <- renderPlot({
              par(mar=c(5,4.5,.7,.5))
              enviMass:::plothomol(homol,
                xlim = isolate(ranges_homol$mass), ylim = isolate(ranges_homol$RT), 
                dmasslim = isolate(ranges_homol$dmass), dRTlim = isolate(ranges_homol$dRT),
                plot_what="mz_RT", omit_theta=omit_theta, plot_those = plot_those
              );
            },res=100)  
          }
      })
      ################################################################################ 
      observe({      
          s2<-input$homol_series_table_rows_selected
          s1<-isolate(input$homol_series_peaks_rows_selected)
          if(length(s2)){
            if(s2>=1){
              cat("\n IN SELECT_2:")
              # output homol. series plot ######################################
              if(length(s1)){
                cat(" A ")
                use_homol_peaks<-match(unique(c(homol[["homol_peaks_relat"]][plot_those,1],homol[["homol_peaks_relat"]][plot_those,2])),homol[["Peaks in homologue series"]][,"peak ID"])
                these_series_1<-as.numeric(strsplit(homol[["Peaks in homologue series"]][use_homol_peaks[s1],c("HS IDs")],"/")[[1]][s2])
                use_emph_point<-homol[["Peaks in homologue series"]][use_homol_peaks[s1],"peak ID"]
              }else{
                cat(" B ")
                these_series_1<-homol[["Homologue Series"]][s2,"HS IDs"]
                use_emph_point<-FALSE
              }
              print(these_series_1)
              output$homol_plot <- renderPlot({
                par(mar=c(5,4.5,.7,.5))
                enviMass:::plothomol(homol,
                  xlim = isolate(ranges_homol$mass), ylim = isolate(ranges_homol$RT), 
                  dmasslim = isolate(ranges_homol$dmass), dRTlim = isolate(ranges_homol$dRT),
                  plot_what="mz_RT", omit_theta=omit_theta, plot_those = plot_those,
                  emph_point=use_emph_point,
                  emph_series=these_series_1
                );
              },res=100)   
              # output chromatograms #########################################
              these_peaks<-as.numeric(strsplit(homol[["Homologue Series"]][these_series_1,"peak IDs"],",")[[1]])
              isolate(refresh_homol$c<-these_peaks)
            }
          }else{
            cat("\n IN DESELECT_2:");print(s2);
            # output homol. series plot ######################################
            if(length(s1)){
              these_series<-as.numeric(strsplit(homol[["Peaks in homologue series"]][use_homol_peaks[s1],c("HS IDs")],"/")[[1]])
              output$homol_plot <- renderPlot({   
                par(mar=c(5,4.5,.7,.5))
                enviMass:::plothomol(homol,
                  xlim = isolate(ranges_homol$mass), ylim = isolate(ranges_homol$RT), 
                  dmasslim = isolate(ranges_homol$dmass), dRTlim = isolate(ranges_homol$dRT),
                  plot_what="mz_RT", omit_theta=omit_theta, plot_those = plot_those,
                  emph_point=homol[["Peaks in homologue series"]][use_homol_peaks[s1],"peak ID"],
                  emph_series=these_series
                );
              },res=100)
            }else{  
              output$homol_plot <- renderPlot({
                par(mar=c(5,4.5,.7,.5))
                enviMass:::plothomol(homol,
                  xlim = isolate(ranges_homol$mass), ylim = isolate(ranges_homol$RT), 
                  dmasslim = isolate(ranges_homol$dmass), dRTlim = isolate(ranges_homol$dRT),
                  plot_what="mz_RT", omit_theta=omit_theta, plot_those = plot_those
                );
              },res=100)  
            }
          }
      })
      ################################################################################ 
      observe({
        refresh_homol$c
        if(refresh_homol$c[1]>0){
          cat("\n IN CHROMAT");print(refresh_homol$c);
          # output chromatograms #########################################
          output$homol_chromat <- renderPlot({
            par(mar=c(4,4,.3,.2))
            enviMass:::plotchromat(
              MSlist,
              peakIDs=refresh_homol$c,
              RTlim=FALSE,
              normalize=FALSE
            );
          },res=100) 
        }      
      })
      ################################################################################
      observeEvent(input$homol_plot_dblclick, { 
          brush <- isolate(input$homol_plot_brush)
          if (!is.null(brush)) {
            cat("\n Zoom in_1")
            isolate(ranges_homol$mass <- c(brush$xmin, brush$xmax))
            isolate(ranges_homol$RT <- c(brush$ymin, brush$ymax))
          } else {
            cat("\n Zoom out full_1")
            isolate(ranges_homol$mass <- FALSE)
            isolate(ranges_homol$RT <- FALSE)
          }
          refresh_homol$a<-(refresh_homol$a+1) # valid in both cases
      })
      observeEvent(input$homol_plot_click, { # NOTE: brushing already triggers a click -> use brush with delay=0, which embeds the slower click
          cat("\n Zoom out part_1_a")
          brush <- isolate(input$homol_plot_brush)
          if (is.null(brush)) {
            cat("\n Zoom out part_1_b")
            if(isolate(ranges_homol$mass[1])!=FALSE){
              old_range_mass<-abs(isolate(ranges_homol$mass[2]-ranges_homol$mass[1]))
              isolate(ranges_homol$mass[1]<-ranges_homol$mass[1]-.3*old_range_mass)
              isolate(ranges_homol$mass[2]<-ranges_homol$mass[2]+.3*old_range_mass)
            }
            if(isolate(ranges_homol$RT[1])!=FALSE){
              old_range_RT<-abs(isolate(ranges_homol$RT[2]-ranges_homol$RT[1]))
              isolate(ranges_homol$RT[1]<-ranges_homol$RT[1]-.1*old_range_RT)
              isolate(ranges_homol$RT[2]<-ranges_homol$RT[2]+.1*old_range_RT)
            }
            refresh_homol$a<-(refresh_homol$a+1)
          }else{
            cat("\n Doing hover - and nothing")
          }   
      })
      ################################################################################
       observeEvent(input$homol_counts_dblclick, { 
          brush <- isolate(input$homol_counts_brush)
          if (!is.null(brush)) {
            cat("\n Zoom in_1d")
            isolate(ranges_homol$dmass <- c(brush$xmin, brush$xmax))
          } else {
            cat("\n Zoom out full_1d")
            isolate(ranges_homol$dmass <- FALSE)
          }
          refresh_homol$a<-(refresh_homol$a+1) # valid in both cases
      })
      observeEvent(input$homol_counts_click, { # NOTE: brushing already triggers a click -> use brush with delay=0, which embeds the slower click
          cat("\n Zoom out part_1_ad")
          brush <- isolate(input$homol_counts_brush)
          if (is.null(brush)) {
              cat("\n Zoom out part_1_bd")
              if(isolate(ranges_homol$dmass[1])!=FALSE){
                old_range_dmass<-abs(isolate(ranges_homol$dmass[2]-ranges_homol$dmass[1]))
                isolate(ranges_homol$dmass[1]<-ranges_homol$dmass[1]-.3*old_range_dmass)
                isolate(ranges_homol$dmass[2]<-ranges_homol$dmass[2]+.3*old_range_dmass)
              }  
              refresh_homol$a<-(refresh_homol$a+1)
          }else{
            cat("\n Doing hover_d")
            isolate(ranges_homol$dmass <- c(brush$xmin, brush$xmax))
            refresh_homol$b<-(refresh_homol$b+1)
          }   
      })     
      ################################################################################
       observeEvent(input$homol_RT_dblclick, { 
          brush <- isolate(input$homol_RT_brush)
          if (!is.null(brush)) {
            cat("\n Zoom in_1d")
            isolate(ranges_homol$dRT <- c(brush$xmin, brush$xmax))
          } else {
            cat("\n Zoom out full_1d")
            isolate(ranges_homol$dRT <- FALSE)
          }
          refresh_homol$a<-(refresh_homol$a+1) # valid in both cases
      })
      observeEvent(input$homol_RT_click, { # NOTE: brushing already triggers a click -> use brush with delay=0, which embeds the slower click
          cat("\n Zoom out part_1_ad")
          brush <- isolate(input$homol_RT_brush)
          if (is.null(brush)) {
              cat("\n Zoom out part_1_bd")
              if(isolate(ranges_homol$dRT[1])!=FALSE){
                old_range_dmass<-abs(isolate(ranges_homol$dRT[2]-ranges_homol$dRT[1]))
                isolate(ranges_homol$dRT[1]<-ranges_homol$dRT[1]-.3*old_range_dmass)
                isolate(ranges_homol$dRT[2]<-ranges_homol$dRT[2]+.3*old_range_dmass)
              }  
              refresh_homol$a<-(refresh_homol$a+1)
          }else{
            cat("\n Doing hover_d")
            isolate(ranges_homol$dRT <- c(brush$xmin, brush$xmax))
            refresh_homol$b<-(refresh_homol$b+1)
          }   
      })     
      ################################################################################

    }
  ))





