
    conditionalPanel( 
		condition = "output.textit != 'Waiting...'",  	
		tags$h5(""),
		bsAlert("alert_3"),
		tabsetPanel(
		#navbarPage("",
        ########################################################################
        # MEASUREMENTS #########################################################
        ########################################################################
		tabPanel("Measurements",
			tabsetPanel(	
				################################################################
				tabPanel("File manager",
					HTML('<p><a href="http://www.looscomputing.ch/eng/enviMass/inputs/files.htm" style="color:rgb(60, 100, 60); text-decoration: none"; target="_blank"><p align="left">&#8594; Check help for details</a></p>'),	
					bsCollapse(multiple = FALSE, open = "files_open", id = "files",
						# ADD FILE #################################################
						bsCollapsePanel("Add LC-MS file", 		
							helpText("To add a new file.mzXML to your project, first browse to its folder location and select the file. 
								Based on the filename, a guess for the file type, its ionization mode and date is made. 
								Then specify/check its properties and finally press the load button. "),		
							helpText("Permissible date formats in file names: 
								_YYYY_MM_DD_ / _YY_MM_DD_ / _YYYYMMDD_ / _YYMMDD_ / -YYYY-MM-DD- / -YY-MM-DD- / -YYYYMMDD- / -YYMMDD- 
								 (where - and _ at start and end of the date field can be empty spaces as well). First match in name used."),		
							HTML('<hr noshade="noshade" />'),
							fileInput("Measadd_path", "Select centroided .mzXML file:", multiple = FALSE, accept = c(".mzXML",".raw")),
							bsPopover("Measadd_path", 
									title = "WARNING",
									content = "Files must be centroided. If reload fails, press cancel in the file upload window first", 
									placement = "right", trigger = "hover"),
							textOutput("file_to_load"),
							HTML('<hr noshade="noshade" />'),
							HTML('<h1 align = "left"> &#x21e9; </h1> '),	
							fluidRow(
								column(width = 5, textInput("Measadd_name", "Name:", value = "File XY")),
								column(width = 5, selectInput("Measadd_type", "Type:", choices = c("sample", "blank", "calibration", "spiked"))),
								column(width = 5, selectInput("Measadd_incl", "Include?", choices = c("TRUE","FALSE"))),
								column(width = 5, selectInput("Measadd_mode", "Choose ionization mode:", choices = c("positive", "negative")))	
							),
							conditionalPanel(
								condition = "input.Measadd_type != 'calibration' & input.Measadd_type != 'spiked'",						
								HTML('<hr noshade="noshade" />'),
								fluidRow(
									column(width = 5, textInput("Measadd_place", "Place:", value = "Rhine")),		
									column(width = 5, dateInput("Measadd_date", "Date:", value = NULL, min = NULL,max = NULL, format = "yyyy-mm-dd", startview = "month",weekstart = 0, language = "en")),	
									column(width = 5, textInput("Measadd_time", "Time (HH:MM:SS):", value = "12:00:00")),							
									column(width = 5,
										conditionalPanel(
											condition = "input.Measadd_type == 'sample'",								
												textInput("Measadd_tag3", "Replicate group (tag3):", value = "FALSE")
										)
									),
									column(width = 5, textInput("Measadd_ID2", "Custom ID:", value = "FALSE")),
									column(width = 5, selectInput("Measadd_profiled", "Use for profiling (if Settings/Profiling/Omit adjusted to do so)?", choices = c("TRUE","FALSE"), selected = "TRUE"))
								)			
							),	
							conditionalPanel(
								condition = "input.Measadd_type == 'calibration'", 
								HTML('<hr noshade="noshade" />'),
								fluidRow(
									column(width = 5, textInput("Measadd_tag1", "Concentration (no units; tag1)", value = "FALSE")),
									column(width = 5, textInput("Measadd_tag2", "Name of calibration file set (required; tag2)", value = "Group A")),
									column(width = 5, dateInput("Measadd_cal_date1", "Date start", value = NULL, min = NULL,max = NULL, format = "yyyy-mm-dd", startview = "month",weekstart = 0, language = "en")),
									column(width = 5, textInput("Measadd_cal_time1", "Time start (HH:MM:SS)", value = "12:00:00")),
									column(width = 5, dateInput("Measadd_cal_date2", "Date end", value="2018-01-01", min = NULL,max = NULL, format = "yyyy-mm-dd", startview = "month",weekstart = 0, language = "en")),
									column(width = 5, textInput("Measadd_cal_time2", "Time end (HH:MM:SS)", value = "11:59:59")),
									column(width = 5, textInput("Measadd_cal_ID2", "Custom ID:", value = "FALSE"))					
								)
							),
							conditionalPanel(
								condition = "input.Measadd_type == 'spiked'", 
								HTML('<hr noshade="noshade" />'),
								fluidRow(
									column(width = 5, textInput("Measadd_spiked_tag2", "ID of reference file (tag2)", value = "FALSE")),
									column(width = 5, dateInput("Measadd_recov_date", "Date start", value = NULL, min = NULL,max = NULL, format = "yyyy-mm-dd", startview = "month",weekstart = 0, language = "en")),
									column(width = 5, textInput("Measadd_recov_time", "Time start (HH:MM:SS)", value = "12:00:00")),
									column(width = 5, textInput("Measadd_recov_ID2", "Custom ID:", value = "FALSE"))
								)
							),						
							HTML('<hr noshade="noshade" />'),
							HTML('<h1 align="left"> &#x21e9; </h1> '),
							bsButton("Load_file", "Load file into project", style = "success"),
							textOutput("had_meas_added")
						, style = "success"),
						# DELETE FILE ##################################################
						bsCollapsePanel("Delete LC-MS file", 		
							tags$h5("Delete file by its unique ID from the below file table"),
							textInput("Measdel_ID", "ID:", value = "123"),
							bsButton("Measdel","Remove",style="primary")		
						,style="success"),				
						# MODIFY FILE ##################################################
						bsCollapsePanel("Modify specifications for a single file", 
							fluidRow(
								column(width = 4, helpText("Load settings of a file into below mask by its ID, modify and then save the new settings into the main table. Modifications make a full recalculation default.")),
								column(width = 3, textInput("Modif_ID", "ID:", value = "123"), bsButton("Modif_load","Load",style="primary"))
							),
							HTML('<hr noshade="noshade" />'),
							fluidRow(
								column(width = 5, textInput("Modif_name", "Name:", value = "Sample 1")),
								column(width = 5, selectInput("Modif_type", "Type:", choices = c("sample", "blank", "calibration", "spiked"))),
								column(width = 5, selectInput("Modif_include","Include in workflow?",choices = c("TRUE","FALSE"),selected="TRUE")),
								column(width = 5, selectInput("Modif_mode", "Choose ionization mode:", choices = c("positive", "negative")))	
							),
							conditionalPanel(
								condition = "input.Modif_type == 'sample' | input.Modif_type == 'blank'", 						
								HTML('<hr noshade="noshade" />'),
								fluidRow(
									column(width = 5, textInput("Modif_place", "Place:", value = "Rhine")),										
									column(width = 5, dateInput("Modif_date", "Date", value = NULL, min = NULL,max = NULL, format = "yyyy-mm-dd", startview = "month",weekstart = 0, language = "en")),	
									column(width = 5, textInput("Modif_time", "Time:(HH:MM:SS)", value = "12:00:00")),
									column(width = 5, textInput("Modif_tag3", "Replicate group (tag3)", value = "FALSE")),
									column(width = 5, textInput("Modif_ID2", "Custom ID", value = "FALSE")),
									column(width = 5, selectInput("Modif_profiled","Use for profiling (if Settings/Profiling/Omit adjusted to do so)?",choices = c("TRUE","FALSE"),selected="TRUE"))								
								)
							),
							conditionalPanel(
								condition = "input.Modif_type == 'calibration'", 
								HTML('<hr noshade="noshade" />'),
								fluidRow(
									column(width = 5, textInput("Modif_tag1", "Concentration (no units; tag1)", value = "0")),
									column(width = 5, textInput("Modif_tag2", "Calibration group (tag2)", value = "FALSE")),
									column(width = 5, dateInput("Modif_cal_date1", "Date start", value = NULL, min = NULL,max = NULL, format = "yyyy-mm-dd", startview = "month",weekstart = 0, language = "en")),
									column(width = 5, textInput("Modif_cal_time1", "Time start (HH:MM:SS)", value = "12:00:00")),
									column(width = 5, dateInput("Modif_cal_date2", "Date end", value = NULL, min = NULL,max = NULL, format = "yyyy-mm-dd", startview = "month",weekstart = 0, language = "en")),
									column(width = 5, textInput("Modif_cal_time2", "Time end (HH:MM:SS)", value = "11:59:59")),
									column(width = 5, textInput("Modif_cal_ID2", "Custom ID", value = "FALSE"))							
								)
							),
							conditionalPanel(
								condition = "input.Modif_type == 'spiked'", 
								HTML('<hr noshade="noshade" />'),
								fluidRow(
									column(width = 5, textInput("Modif_spiked_tag2", "ID of file to subtract from (tag2)", value = "FALSE")),
									column(width = 5, dateInput("Modif_recov_date", "Date", value = NULL, min = NULL,max = NULL, format = "yyyy-mm-dd", startview = "month",weekstart = 0, language = "en")),	
									column(width = 5, textInput("Modif_recov_time", "Time:(HH:MM:SS)", value = "12:00:00")),
									column(width = 5, textInput("Modif_recov_ID2", "Custom ID", value = "FALSE"))
								)
							),						
							HTML('<hr noshade="noshade" />'),
							bsButton("Modif_export","Save",style="primary")
						,style="success"),
						# MODIFY CALIBRATION GROUP ######################################
						bsCollapsePanel("Modify, copy or delete a calibration group", 
							fluidRow(
								column(width = 4, helpText("Load a calibration group to be modified, deleted or copied below. To do so, select the ionization mode, insert the calibration group name (tag2 in the file table) and press Load.")),
								column(width = 3, selectInput("Modif_cal_mode", "Choose ionization mode:", choices = c("positive", "negative"))),
								column(width = 3, textInput("Modif_cal_group", "Calibration group (tag2):", value = "existinggroup"), bsButton("Load_cal","Load",style="primary"))
							),
							HTML('<hr noshade="noshade" />'),
							htmlOutput('Modif_cal_text_load'),
							HTML('<hr noshade="noshade" />'),
							helpText("Modify the specifications for all files of the above loaded calibration group and press Save to make the changes permanent"),
							fluidRow(
								column(width = 5, dateInput("Modif_calgroup_date1", "Date start", value = NULL, min = NULL,max = NULL, format = "yyyy-mm-dd", startview = "month",weekstart = 0, language = "en")),
								column(width = 5, textInput("Modif_calgroup_time1", "Time start (HH:MM:SS)", value = "12:00:00")),
								column(width = 5, dateInput("Modif_calgroup_date2", "Date end", value = NULL, min = NULL,max = NULL, format = "yyyy-mm-dd", startview = "month",weekstart = 0, language = "en")),
								column(width = 5, textInput("Modif_calgroup_time2", "Time end (HH:MM:SS)", value = "11:59:59"))
							),
							bsButton("Change_cal","Save",style="primary"),					
							HTML('<hr noshade="noshade" />'),
							helpText("Copy the calibration group to create a new calibration group with a different name. This will directly incorporate the above (changed) specifications for the new group."),
							textInput("Copy_cal_group", "New calibration group name (no underscores permitted!):", value = "newgroup"),
							bsButton("Copy_cal","Copy",style="primary"),	
							HTML('<hr noshade="noshade" />'),
							helpText("Delete all the files and any associated contents of the loaded calibration group."),
							bsButton("Del_cal","Delete",style="primary"),
							bsModal("Del_cal_confirm", "Sure about deleting the selected calibration file set?", "Del_cal", size = "small",
								bsButton("yes_delete_cal","Yes",style="warning")
							)
						,style="success"),		
						# BATCH UPLOAD ##################################################
						bsCollapsePanel("Batch upload from folder", 	
							tags$h5("Read in batches of files (.mzXML) from a folder; file specifications will be guessed and can later be change in the above file modification section. 
							File dates are set in the sequence of files provided in the folder."),	
							textInput("import_file_folder", "Folder path:", value = "C:\\...\\folder"),
							bsPopover("import_file_folder", 
								title = "Insert full path, including the project folder, but excluding the logfile.emp.",
								content = "Using your OS explorer, you may navigate into your project folder and copy/paste the full path.", 
								placement = "right", trigger = "hover"),
							checkboxInput("Import_file_folder_overwrite", "Import files with same name?", FALSE),
							bsPopover("Import_file_folder_overwrite", 
								title = "File import to existing files",
								content = "If a file in the folder with the same filename as one already existing in the project is found, should it be imported (new ID assigned)?", 
								placement = "right", trigger = "hover"),
							checkboxInput("Import_file_folder_namedate", "Guess date from file name?", FALSE),
							bsPopover("Import_file_folder_namedate", 
								title = "Set file date",
								content = "If checked, file dates are tried to be derived from the file name - otherwise, consecutive file dates are set (starting with the last unused file date, if available). Permissible date formats in file names: _YYYY_MM_DD_ / _YY_MM_DD_ / _YYYYMMDD_ / _YYMMDD_ / -YYYY-MM-DD- / -YY-MM-DD- / -YYYYMMDD- / -YYMMDD- (where - and _ at start and end of the date field can be empty spaces as well). First match in name used.", 
								placement = "right", trigger = "hover"),
							bsButton("Import_file_folder","Import",style="primary"),
							HTML('<hr noshade="noshade" />'),
							textOutput("had_import_folder")						
						,style="success"),
						# IMPORT from other project ######################################
						bsCollapsePanel("Import files from another project", 		
							tags$h5("Select project folder to import files from:"),
							textInput("import_pro_dir", "", value = "C:\\...\\old_project_name"),
							bsPopover("import_pro_dir", 
								title = "Insert full path, including the project folder, but excluding the logfile.emp.",
								content = "Using your OS explorer, you may navigate into your project folder and copy/paste the full path.", 
								placement = "right", trigger = "hover"),
							checkboxInput("Merge_project", "Omit duplicates?", FALSE),
							bsPopover("Merge_project", 
								title = "File duplicate handling",
								content = "A file with the same type, time, date, ionization and place as one which already exists will not be imported.", 
								placement = "right", trigger = "hover"),
							bsButton("Import_project","Import",style="primary"),		
							HTML('<hr noshade="noshade" />'),
							textOutput("had_import_project")	
						,style="info"),
						# File overview ##################################################
						bsCollapsePanel("File overview", 
							helpText("The below plot indicates available files as dots at their respective date and time, listed over the different file categories 
							and seperately for each of the two ion modes."),
							helpText("Select a time period via mouse brush to list the file IDs for it. Double-click into the area to zoom in; double-click again to zoom out."),
							plotOutput("file_overview", 
								dblclick = "file_overview_dblclick",
								brush = brushOpts(
									id = "file_overview_brush",
									resetOnNew = TRUE,
									direction = "x",
									fill="red"
								),
								height = "600px"
							),
							HTML('<hr noshade="noshade" />'),
							HTML('<font size="5"> + </font><font size="3"> Last selected file IDs, positive ionization:</font>'),
							htmlOutput("info_files_pos_samp"),htmlOutput("info_files_pos_blind"),htmlOutput("info_files_pos_cal"),htmlOutput("info_files_pos_calgroup"),htmlOutput("info_files_pos_spiked"),
							HTML('<hr noshade="noshade" />'),		
							HTML('<font size="5"> - </font><font size="3"> Last selected file IDs, negative ionization:</font>'),
							htmlOutput("info_files_neg_samp"),htmlOutput("info_files_neg_blind"),htmlOutput("info_files_neg_cal"),htmlOutput("info_files_neg_calgroup"),htmlOutput("info_files_neg_spiked")
						,style="warning")
					),	
					HTML('<hr noshade="noshade" />'),
					DT::dataTableOutput("measurements")	
				),
				################################################################
				tabPanel("Method setup",
					HTML('<hr noshade="noshade" />'),
					tags$p(align="justify","Select a file by its ID to show information on its method and its individual scans. Moreover, the below Scan type definition 
					section allows to select parameters by which the different scans are to be distinguished/combined. While enviMass uses the msLevel 1 scans 
					only at the time being, such Scan type definitions will mainly be used to process msLevel 2 scans in the future."),											
					div(style = widget_style3, numericInput("sel_scans_ID", "Type in file ID:", 0)),
					textOutput('scan_viewer_name'),
					conditionalPanel(			
						condition = "(input.sel_scans_ID != 0) & (output.scan_viewer_name != 'File name: Invalid file ID') & (output.scan_viewer_name != 'No scan type definition parameters selected - please select at least one!') & (output.scan_viewer_name != '.mzXML file not available')",						
						textOutput('scan_viewer_type'),
						textOutput('scan_viewer_mode'),
						HTML('<p style="background-color:darkgrey"; align="center"> <font color="#FFFFFF"> Method information </font></p> '),
						tableOutput("instrument_Info"),
						tableOutput("run_Info")
					),	
					conditionalPanel(			
						condition = "(input.sel_scans_ID != 0) & (output.scan_viewer_name != 'File name: Invalid file ID')",
						HTML('<h1 align="center"> &#x21e9; </h1> '),
						HTML('<p style="background-color:darkgrey"; align="center"> <font color="#FFFFFF"> Scan type definition </font></p> '),
						fluidRow(
							column(width = 6, 						
								selectInput("method_definition", label = "Parameters for defining scan types", 
									choices = c("polarity", "msLevel", "collisionEnergy", "precursorMZ", "ionisationEnergy", "basePeakMZ"), 
									selected = c("polarity", "msLevel", "collisionEnergy", "precursorMZ"), 
									multiple = TRUE, selectize = TRUE, width = NULL, size = NULL)
							),
							column(width = 4, 								
								checkboxInput("method_MS1_range_use", "Use MS1 mass seperator?", FALSE),
								numericInput("method_MS1_range", "-> m/z seperator:", 70, min = 0)
							)							
						)
					),
					conditionalPanel(			
						condition = "(input.sel_scans_ID != 0) & (output.scan_viewer_name != 'File name: Invalid file ID') & (output.scan_viewer_name != 'No scan type definition parameters selected - please select at least one!') & (output.scan_viewer_name != '.mzXML file not available')",
						tableOutput("heads_summary"),
						HTML('<h1 align="center"> &#x21e9; </h1> '),
						HTML('<p style="background-color:darkgrey"; align="center"> <font color="#FFFFFF"> Scan viewer </font></p> '),	
						div(style = widget_style3, numericInput("sel_scans_number", "Number of viewed scans", 50)),
						shinyTree(outputId = "scan_tree", checkbox = FALSE, search = FALSE, dragAndDrop = FALSE)
					)
				)
				################################################################								
			)
	    ),
        ########################################################################
        # Compounds ############################################################
        ########################################################################
        tabPanel("Compounds",
		 tabsetPanel(
            tabPanel("Internal standards",
				HTML('<p><a href="http://www.looscomputing.ch/eng/enviMass/inputs/IS.htm" style="color:rgb(60, 100, 60); text-decoration: none"; target="_blank"><p align="left">&#8594; Check help for details</a></p>'),	
				bsCollapse(multiple = FALSE, open = "IS_comp_open", id = "IS_comp",
					bsCollapsePanel("Add or modify an internal standard compound", 
						fluidRow(
							column(width = 5, helpText("To add a new internal standard, fill out the below form and press Add") ),
							column(width = 2, offset = 0.3, bsButton("AddIS","Add",style="primary"))
						),
						HTML('<hr noshade="noshade" />'),
						fluidRow(
							column(width = 5, helpText("To modify a compound insert its ID into the field to the right and press Load,
								which updates the below form. Modifications in the form can then be saved with the Modify button.
								If the compound ID in the form is changed, a new entry to the compound table with this ID will be made.
								Beware: if another compound with this new ID exists, it will be replaced in the table.") ),
							column(width = 4, 
								textInput("ISmodif_ID", "Compound ID:", value = "123_XYZ"),
								bsButton("LoadIS","Load",style="primary")),
							column(width = 3, offset = 0.3, bsButton("ModifIS","Modify",style="primary"))
						),
						HTML('<hr style="border: .8px solid " />'),
						fluidRow(
							column(width = 2, textInput("ISadd_ID", "Unique ID:", value = "123_XYZ")),  
							column(width = 8, helpText("Must be unique, numbers and characters permitted; no empty spaces or special characters or underscores permitted.") )							
						),
						HTML('<hr style="border: .8px solid " />'),
						fluidRow(         
							column(width = 3, textInput("ISadd_name", "Name:", value = "CompoundX")),
							column(width = 3, selectInput("ISadd_charge", label="Ionization mode", choices=c("positive","negative"), selected = "positive", multiple = FALSE))									
						),
						HTML('<hr style="border: .8px solid " />'),
						fluidRow(
							column(width = 3, textInput("ISadd_formula", "Formula:", value = "C10H12O10")),
							column(width = 8, helpText("Must consist of an upper case letter, possibly followed by lower case letters; to refer to individual isotopes (e.g., from isotope 
									labelling of a molecule, e.g., N5 vs. [15]N2N3), square brackets may precede the capital letter. Any other symbols which may be part of a chemical formula (e.g., charges (+), 
									dashes, asterisks, ...) are not permitted.") )										
						),
						HTML('<hr noshade="noshade" />'),
						fluidRow(
							column(width = 3, textInput("ISadd_RT", "Retention time (RT) [min]:", value = "7.52")),    
							column(width = 3, checkboxInput("ISadd_RTtol_use", "Use specific +/- RT tolerance [min]:", FALSE)),                    
							column(width = 3, textInput("ISadd_RTtol","", value = "2"))
						),
						HTML('<hr noshade="noshade" />'),
						fluidRow(						
							column(width = 4, 
								selectInput("ISadd_add", label="Main adduct:", choices= "FALSE", selected = "FALSE", multiple = FALSE),
								checkboxInput("ISadd_rest_adduct", "Restrict screening to main adduct?", FALSE)),
							column(width = 8, helpText("A compound-specific adduct can be defined here; general adducts to be considered for all compounds can be defined in the Settings/Screening/Adducts tab. 	
														Unless the below restriction, the compound-specific adduct is used alongside the general ones. 
														Seperate compound entries have to be made when including more than one compound-specific adduct."))							
						),
						HTML('<hr noshade="noshade" />'),
						fluidRow(	
							column(width = 4, checkboxInput("ISadd_use_recal", "To be used for m/z recalibration?", TRUE)),
							column(width = 4, checkboxInput("ISadd_use_screen", "To be used for screening?", TRUE))
						),
						HTML('<hr noshade="noshade" />'),
						fluidRow(
							column(width = 4, textInput("ISadd_remark", "remark", value = "none")),
							column(width = 4, textInput("ISadd_tag1", "tag1", value = "none")),
							column(width = 4, textInput("ISadd_tag2", "tag2", value = "none")),
							column(width = 4, textInput("ISadd_tag3", "tag3", value = "none"))
						),
						HTML('<hr noshade="noshade" />'),
						#div(style = widget_style,							
						#	checkboxInput("ISadd_date", "Restrict to temporal range?", FALSE),	
						#	dateRangeInput("ISadd_date_range", label="", start = NULL, end = NULL,
						#		min = NULL, max = NULL, format = "yyyy-mm-dd",
						#		startview = "month", weekstart = 0, language = "en",
						#		separator = " to ")
						#),
						fluidRow(
							column(width = 12, tags$h4("Quantification settings")),
							column(width = 4,
								helpText("Used adduct:"),
								selectInput("IS_quant_add", label=NULL, choices= "FALSE", selected = "FALSE", multiple = FALSE)
							),		
							column(width = 5,
								helpText("Used centroid peak (sorted by mass, 1=monoisotopic):"),
								numericInput("IS_quant_peak", label=NULL, 1)
							),
							column(width = 8,
								fluidRow(
									column(width = 5,
										HTML('<span class="help-block">Lower log<sub>10</sub> intensity bound of centroid peak used for quantification (.-separated):</span>'),
										textInput("Lower_intensity_bound", label=NULL, value = "0")),						
									column(width = 5,
										HTML('<span class="help-block">Upper log<sub>10</sub> intensity bound of centroid peak used for quantification (.-separated):</span>'),
										textInput("Upper_intensity_bound", label=NULL, value = "Inf"))
								)
							),
							column(width = 5,
								helpText("To rank several centroid peaks to quantify with, use:"),
								selectInput("IS_quant_rule", label=NULL, choices=c("most intense peak","closest RT","closest m/z"), selected = "most intense peak",width='400px'))							
						)
					,style="success"),
					bsCollapsePanel("Remove internal standard compound", 					
						helpText("To delete a compound from the list, type in its ID and press Delete"),
						textInput("ISdelete_ID", "ID for deletion:", value = "123_XYZ"),          
						bsButton("DeleteIS","Delete",style="primary")					
					,style="success"),
					bsCollapsePanel("Import / export internal standard compound list", 		
						helpText("Import IS compund list.txt file, e.g. from the dataframes folder of another project (where it can be found as IS.txt):"),					
						fileInput("ISlist_path", NULL, multiple = FALSE, accept = c(".txt")),
						checkboxInput("ISlist_save_copy", "Save a copy of the current IS compound table below if it is to be replaced by the import?", TRUE),
						bsPopover("ISlist_save_copy", 
							title = "Safety copy",
							content = "The copy can be found in the dataframes folder of your project, named as IS_date_time. You may reload this backup later to undo a compound import.", 
							placement = "right", trigger = "hover"),
						HTML('<hr noshade="noshade" />'),
						helpText("Export and save the below IS table as .txt file. The latter can again be reloaded after modifications, using the above import."),
						shinySaveButton(id="download_IS", label="Save", title="Save below IS table", filetype=list(text='txt'), buttonType = "default", class = NULL)	
					,style="info"),
					bsCollapsePanel("Modify in external editor", 					
						HTML('
							<p><font>
								The below compound table can be exported as .txt file and modified in external text editors, OpenOffice Calc or Excel and then again imported. 
								Export and import MUST be done with the above Import / export functions to check for correctness of the tables. The .txt files are tab-separated when exported and must be so for import.
								Mind the character encoding (e.g., ANSI) of .txt files when modifying and ensure your imported compound table does not contain empty rows, especially at the table end.
							</font></p>
							<p style="background-color:darkred"; align="center"> <font color="#FFFFFF"> 
								Any external modifications must strictly adhere to the following rules:
							</font></p> 
							<ol>
								<li>Input format: text file (.txt), tab delimited.</li>
								<li>Column headers: DO NOT MODIFY.</li>
								<li>No empty columns. If you are not sure what to fill in, simply stay with what is set in the exported file.</li>
								<li>Absolutely NO duplicated IDs or IDs containing underscores.</li>
								<li>Do not delete columns; their number, order and content types are all fixed for good reasons.</li>
								<li>No uncompleted entries per compounds,i.e., row.</li>
								<li>Numeric entries with decimal points: dot-separated.</li>	
							</ol>
							<p style="background-color:darkgreen"; align="center"> <font color="#FFFFFF"> 
								Help for column contents of the internal standard table:
							</font></p> 
							<ol>
								<li><b>ID</b>: must be unique, numbers and characters permitted; no empty spaces or special characters. Absolutely NO duplicate IDs and NO IDs with underscores permitted.</li>
								<li><b>Name</b>: (= compound name) no special signs permitted. Use big and small letters, numbers, underscores, hyphen, brackets and empty spaces - and absolutely nothing else.</li>
								<li><b>Formula</b>: molecular formula of compounds: must consist of an upper case letter, possibly followed by lower case letters; to refer to individual isotopes (e.g., from isotope 
									labelling of a molecule, e.g., N5 vs. [15]N2N3), square brackets may precede the capital letter. Any other symbols which may be part of a chemical formula (e.g., charges (+), 
									dashes, asterisks, ...) are not permitted. </li>
								<li><b>RT</b>: retention time, IN MINUTES. Must be specified. The value will be converted to seconds during workflow usage.</li>
								<li><b>RT_tolerance</b>: FALSE or a compound-specific retention time tolerance given IN MINUTES. Overwrites the one set as standard value in tab Settings/Screening/IS. 
									The RT_tolerance will be automatically converted to seconds during workflow usage.</li>
								<li><b>main_adduct</b>: FALSE or name of a special adduct to be used for this compound entry. Valid adduct names can be found in tab Settings/Adduct.</li>								
								<li><b>ion_mod</b>e: positive or negative. If a compound is to be screened in both modes, two entries (table rows) are required.</li>
								<li><b>use_for_recalibration</b>: TRUE or FALSE.</li>
								<li><b>use_for_screening</b>: TRUE or FALSE.</li>
								<li><b>restrict_adduct</b>: TRUE or FALSE. 
									TRUE: only use the <b>main_adduct</b> (must then be specified) for this compound and ignore the ones specified in tab Settings/Adduct? 
									FALSE: the adducts specified in tab Settings/Screening/Adducts and the one set in column <b>main_adduct</b> (if specified) are all considered for this compound.</li>
								<li><b>Remark</b>: Character string for your remark on this compound (no tabs.</li>
								<li><b>tag1</b>: Character string for specification of the compound (no tabs), e.g. pharmaceutical.</li>	
								<li><b>tag2</b>: Character string for further specifications (no tabs).</li>
								<li><b>tag3</b>: Character string for further specifications (no tabs).</li>
								<li>from: ignore.</li>
								<li>to: ignore.</li>
								<li><b>Lower_intensity_bound</b>: peaks with log<sub>10</sub> intensities below this bound are not used for quantification. Set to 0 to omit; can be interactively set for each IS compound after screening, in the Results tab.</li>	
								<li><b>Upper_intensity_bound</b>: peaks with log<sub>10</sub> intensities above this bound are not used for quantification. Set to Inf to omit; can be interactively set for each IS compound after screening, in the Results tab.</li>
								<li><b>Quant_adduct</b>: ESI adduct to be used for quantification purposes. Valid adduct names can be found among those selected in tab Settings/Adduct or must be in agreement with <b>main_adduct</b>. </li>
								<li><b>Quant_peak</b>: Integer number refering to the centroid peak (ordered by increasing mass; 1 = monoisotopic) of the adduct to be used for quantification purposes. </li>
								<li><b>Quant_rule</b>: Any ONE of exactly these choices "most intense peak","closest RT" or "closest m/z". Used for quantification if several peaks are available. 
									The one quantified concentration after applying the rule will be listed first in the quantification tables for targets linked to this internal standard.</li>
							</ol>
						')		
					,style="info")
				),
                HTML('<hr noshade="noshade" />'),
                DT::dataTableOutput("IS")
            ),
            tabPanel("Targets", 
				HTML('<p><a href="http://www.looscomputing.ch/eng/enviMass/inputs/targets.htm" style="color:rgb(60, 100, 60); text-decoration: none"; target="_blank"><p align="left">&#8594; Check help for details</a></p>'),	
				bsCollapse(multiple = FALSE, open = "target_comp_open", id = "target_comp",
					bsCollapsePanel("Add or modify a target compound", 
						fluidRow(
							column(width = 5, helpText("To add a new target, fill out the below form and press Add.") ),
							column(width = 2, offset = 0.3, bsButton("Addtargets","Add",style="primary"))
						),
						HTML('<hr noshade="noshade" />'),
						fluidRow(
							column(width = 5, helpText("To modify a compound insert its ID into the field to the right and press Load,
								which updates the below form. Modifications in the form can then be saved with the Modify button.
								If the compound ID in the form is changed, a new entry to the compound table with this ID will be made.
								Beware: if another compound with this new ID already exists, it will be replaced in the table.") ),
							column(width = 4, 
								textInput("targetmodif_ID", "Compound ID:", value = "123_XYZ"),
								bsButton("Loadtarget","Load",style="primary")),
							column(width = 3, offset = 0.3, bsButton("Modiftarget","Modify",style="primary"))
						),
						HTML('<hr style="border: .8px solid " />'),
						fluidRow(
							column(width = 2, textInput("targetsadd_ID", "Unique ID:", value = "123_XYZ")),          
							column(width = 8, helpText("Must be unique, numbers and characters permitted; no empty spaces or special characters or underscores permitted."))	
						),
						HTML('<hr style="border: .8px solid " />'),
						fluidRow(         
							column(width = 3, textInput("targetsadd_name", "Name:", value = "CompoundX")),
							column(width = 3, selectInput("targetsadd_charge", label="Ionization mode", choices=c("positive","negative"), selected = "positive", multiple = FALSE))									
						),						
						HTML('<hr style="border: .8px solid " />'),
						fluidRow(
							column(width = 3, textInput("targetsadd_formula", "Formula:", value = "C10H12O10")),
							column(width = 8, helpText("Must consist of an upper case letter, possibly followed by lower case letters; to refer to individual isotopes (e.g., from isotope 
									labelling of a molecule, e.g., N5 vs. [15]N2N3), square brackets may precede the capital letter. Any other symbols which may be part of a chemical formula (e.g., charges (+), 
									dashes, asterisks, ...) are not permitted."))		
						),						
						HTML('<hr noshade="noshade" />'),
						fluidRow(
							column(width = 3, textInput("targetsadd_RT", "Retention time (RT) [min]:", value = "7.52")),    
							column(width = 3, checkboxInput("targetsadd_RTtol_use", "Use specific +/- RT tolerance [min]:", FALSE)),                    
							column(width = 3, textInput("targetsadd_RTtol","", value = "2"))
						),
						HTML('<hr noshade="noshade" />'),
						fluidRow(						
							column(width = 4, 
								selectInput("targetsadd_add", label="Main adduct:", choices= "FALSE", selected = "FALSE", multiple = FALSE),
								checkboxInput("targetsadd_rest_adduct", "Restrict screening to main adduct?", FALSE)),
							column(width = 8, helpText("A compound-specific adduct can be defined here; general adducts to be considered for all compounds can be defined in the Settings/Screening/Adducts tab. 	
														Unless the below restriction, the compound-specific adduct is used alongside the general ones. 
														Seperate compound entries have to be made when including more than one compound-specific adduct."))							
						),
						HTML('<hr noshade="noshade" />'),
						fluidRow(	
							column(width = 4, checkboxInput("targetsadd_use_recal", "To be used for m/z recalibration?", TRUE)),
							column(width = 4, checkboxInput("targetsadd_use_screen", "To be used for screening?", TRUE))
						),
						HTML('<hr noshade="noshade" />'),
						fluidRow(
							column(width = 4, textInput("targetsadd_remark", "remark", value = "none")),
							column(width = 4, textInput("targetsadd_tag1", "tag1", value = "none")),
							column(width = 4, textInput("targetsadd_tag2", "tag2", value = "none")),
							column(width = 4, textInput("targetsadd_tag3", "tag3", value = "none"))
						),
						HTML('<hr noshade="noshade" />'),
						#div(style = widget_style,							
						#	checkboxInput("targetsadd_date", "Restrict to temporal range?", FALSE),	
						#	dateRangeInput("targetsadd_date_range", label="", start = NULL, end = NULL,
						#		min = NULL, max = NULL, format = "yyyy-mm-dd",
						#		startview = "month", weekstart = 0, language = "en",
						#		separator = " to ")
						#),
						fluidRow(
							column(width = 12, tags$h4("Quantification settings")),
							column(width = 4,
								helpText("Used adduct:"),
								selectInput("target_quant_add", label=NULL, choices= "FALSE", selected = "FALSE", multiple = FALSE)
							),		
							column(width = 5,
								helpText("Used centroid peak (sorted by mass, 1=monoisotopic):"),
								numericInput("target_quant_peak", label=NULL, 1)
							),
							column(width = 5,
								helpText("ID of internal standard used for calibration & quantification for this target:"),
								textInput("target_quant_ISID", label=NULL, value = "FALSE")
							),							
							column(width = 8,
								fluidRow(
									column(width = 5,
										helpText("First concentration warn level (.-separated):"),
										textInput("warn_1", label=NULL, value = "FALSE")),						
									column(width = 5,
										helpText("Second concentration warn level (.-separated):"),
										textInput("warn_2", label=NULL, value = "FALSE"))
								)
							),
							column(width = 5,
								helpText("To rank several centroid peaks to quantify with, use:"),
								selectInput("target_quant_rule", label=NULL, choices=c("most intense peak","closest RT","closest m/z"), selected = "most intense peak",width='400px'))							
						)
					,style="success"),
					bsCollapsePanel("Remove target compound", 
						helpText("To delete a compound from the list, type in its ID and press Delete"),
						textInput("targetsdelete_ID", "ID for deletion:", value = "123_XYZ"),          
						bsButton("Deletetargets","Delete",style="primary")					
					,style="success"),
					bsCollapsePanel("Import / export target compound list", 		
						helpText("Import target compound list.txt file, e.g. from the dataframes folder of another project (where it can be found as targets.txt):"),					
						fileInput("targetlist_path", NULL, multiple = FALSE, accept = c(".txt")),
						checkboxInput("targetlist_save_copy", "Save a copy of the current target compound table below if it is to be replaced by the import?", TRUE),
						bsPopover("targetlist_save_copy", 
							title = "Safety copy",
							content = "The copy can be found in the dataframes folder of your project, named as targets_date_time. You may reload this backup later to undo a compound import.", 
							placement = "right", trigger = "hover"),
						HTML('<hr noshade="noshade" />'),
						helpText("Export and save the below target compound table as .txt file. The latter can again be reloaded after modifications, using the above import."),
						shinySaveButton(id="download_target", label="Save", title="Save below target compound table", filetype=list(text='txt'), buttonType = "default", class = NULL)	
					,style="info"),
					bsCollapsePanel("Modify in external editor", 					
						HTML('
							<p><font>
								The below compound table can be exported as .txt file and modified in external text editors, OpenOffice Calc or Excel and then again imported. 
								Export and import MUST be done with the above Import / export functions to check for correctness of the tables. The .txt files are tab-separated when exported and must be so for import.
								Mind the character encoding (e.g., ANSI) of .txt files when modifying and ensure your imported compound table does not contain empty rows, especially at the table end.
							</font></p>
							<p style="background-color:darkred"; align="center"> <font color="#FFFFFF"> 
								Any external modifications must strictly adhere to the following rules:
							</font></p> 
							<ol>
								<li>Input format: text file (.txt), tab delimited.</li>
								<li>Column headers: DO NOT MODIFY.</li>
								<li>No empty columns. If you are not sure what to fill in, simply stay with what is set in the exported file.</li>
								<li>Absolutely NO duplicated IDs or IDs containing underscores.</li>
								<li>Do not delete columns; their number, order and content types are all fixed for good reasons.</li>
								<li>No uncompleted entries per compounds,i.e., row.</li>
								<li>Numeric entries with decimal points: dot-separated.</li>
								<li>If you want to make links to the internal standard table via column <b>ID_internal_standard</b>, the internal standard
								must exists in the according table. Thus, better start setting up the internal standard table first.</li>								
							</ol>
							<p style="background-color:darkgreen"; align="center"> <font color="#FFFFFF"> 
								Help for target table column contents:
							</font></p> 
							<ol>
								<li><b>ID</b>: must be unique, numbers and characters permitted; no empty spaces or special characters. Absolutely NO duplicate IDs and NO IDs with underscores permitted.</li>
								<li><b>Name</b>: (= compound name) no special signs permitted. Use big and small letters, numbers, underscores, hyphen, brackets and empty spaces - and absolutely nothing else.</li>
								<li><b>Formula</b>: molecular formula of compounds: must consist of an upper case letter, possibly followed by lower case letters; to refer to individual isotopes (e.g., from isotope 
									labelling of a molecule, e.g., N5 vs. [15]N2N3), square brackets may precede the capital letter. Any other symbols which may be part of a chemical formula (e.g., charges (+), 
									dashes, asterisks, ...) are not permitted. </li>
								<li><b>RT</b>: retention time, IN MINUTES. Must be specified. The value will be converted to seconds during workflow usage.</li>
								<li><b>ID_internal_standard</b>: FALSE or a valid ID of an internal standards to be found in the CURRENT internal standard table. Only targets linked with this ID will be considered for quantification.</li>
								<li><b>RT_tolerance</b>: FALSE or a compound-specific retention time tolerance given IN MINUTES. Overwrites the one set as standard value in tab Settings/Screening/IS. 
									The RT_tolerance will be automatically converted to seconds during workflow usage.</li>
								<li><b>main_adduct</b>: FALSE or name of a special adduct to be used for this compound entry. Valid adduct names can be found in tab Settings/Screening/Adducts.</li>								
								<li><b>ion_mod</b>e: positive or negative. If a compound is to be screened in both modes, two entries (table rows) are required.</li>
								<li><b>use_for_recalibration</b>: TRUE or FALSE.</li>
								<li><b>use_for_screening</b>: TRUE or FALSE.</li>
								<li><b>restrict_adduct</b>: TRUE or FALSE. 
									TRUE: only use the <b>main_adduct</b> (must then be specified) for this compound and ignore the ones specified in tab Settings/Screening/Adducts? 
									FALSE: the adducts specified in tab Settings/Adduct and the one set in column <b>main_adduct</b> (if specified) are all considered for this compound.</li>
								<li><b>Remark</b>: Character string for your remark on this compound (no tabs.</li>
								<li><b>tag1</b>: Character string for specification of the compound (no tabs), e.g. pharmaceutical.</li>	
								<li><b>tag2</b>: Character string for further specifications (no tabs).</li>
								<li><b>tag3</b>: Character string for further specifications (no tabs).</li>
								<li>from: ignore.</li>
								<li>to: ignore.</li>
								<li><b>warn_1</b>: First concentration warn level, without units. Quantified targets above this level are marked orange in the quantification table. </li>
								<li><b>warn_2</b>: Second concentration warn level, without units. Quantified targets above this level are marked red in the quantification table. </li>
								<li><b>Quant_adduct</b>: ESI adduct to be used for quantification purposes. Valid adduct names can be found among those selected in tab Settings/Adduct or must be in agreement with <b>main_adduct</b>. </li>
								<li><b>Quant_peak</b>: Integer number refering to the centroid peak (ordered by increasing mass; 1 = monoisotopic) of the adduct to be used for quantification purposes. </li>
								<li><b>Quant_rule</b>: Any ONE of exactly these choices "most intense peak","closest RT" or "closest m/z". Used for quantification if several peaks are available. 
									The one quantified concentration after applying the rule will be listed first in the quantification tables.</li>
							</ol>
						')		
					,style="info")
				),
				HTML('<hr noshade="noshade" />'),  
                DT::dataTableOutput("targets")
            )
          )  
        ),
        ########################################################################
        # WORKFLOW OPTIONS #####################################################
        ########################################################################
        tabPanel("Workflow options",
			tags$h5("Apply workflow settings to project?"), 
			bsButton("saveflow","Apply",style="warning"),
			bsAlert("alert_1"),
			tags$p(align="justify","Below is an illustrative network graph of dependencies between steps in the workflow. Enabled steps are shown in 
				dark blue; disabled ones in light blue; interactive and zoomable. When modifying workflow steps, parameters or inputs, enviMass 
				dynamically adjusts and minimizes all required recalculations via their relevant dependencies. These recalculations will ultimately be 
				enforced when pressing the left-sided Calculate button."),
				networkD3::forceNetworkOutput("force_workflow", width = 1500, height = 400),
				# block 1 ######################################################
				HTML('<hr noshade="noshade" />'),
				HTML('<p style="background-color:darkgrey"; align="center"> <font color="#FFFFFF"> File upload
					<a href="http://www.looscomputing.ch/eng/enviMass/topics/PW.htm" style="color:rgb(60, 100, 60); text-decoration: none"; target="_blank"><p align="right">&#8594; More info.</a>	
				</font></p> '),
				HTML('<p style="background-color:darkgrey"; align="center"> <font color="#FFFFFF"> Peak picking 
					<a href="http://www.looscomputing.ch/eng/enviMass/topics/peakpicking.htm" style="color:rgb(60, 100, 60); text-decoration: none"; target="_blank"><p align="right">&#8594; More info.</a>	
				</font></p> '),
				HTML('<p style="background-color:darkgrey"; align="center"> <font color="#FFFFFF"> Intensity control filter </font></p> '),
					fluidRow(
						#column(width = 2, radioButtons("qc", "Include?", c("yes"="yes","no"="no"), inline = TRUE)),          
						column(width = 10, offset = 0.3,
							tags$p(align="justify","The quantile distributions of peak intensities are compared between files to detect outliers.")#,
							#HTML('<p><a href="http://www.looscomputing.ch/eng/enviMass/topics/QC.htm" style="color:rgb(60, 100, 60); text-decoration: none"; target="_blank"><p align="right">&#8594; More info.</a></p>')								
						)
					),
				HTML('<hr noshade="noshade" />'),
				HTML('<h1 align="center"> &#x21e9; </h1> '),
				# block 2 ######################################################
				HTML('<p style="background-color:darkblue"; align="center"> <font color="#FFFFFF"> Mass recalibration </font></p> '),
					fluidRow(
						column(width = 2, radioButtons("recal", "Include?", c("yes"="yes","no"="no")) ),
						column(width = 10, offset = 0.3,
							tags$p(align="justify","Theoretical masses of internal standard compounds are used to correct systematic offsets in measured masses. This step can also be included separately for the different ionization modes in the Settings -> Mass recalibration tab."),
							HTML('<p><a href="http://www.looscomputing.ch/eng/enviMass/topics/recalibration.htm" style="color:rgb(60, 100, 60); text-decoration: none"; target="_blank"><p align="right">&#8594; More info.</a></p>')	
						)
					),
				HTML('<p style="background-color:darkblue"; align="center"> <font color="#FFFFFF"> Median intensity normalization </font></p> '),
					fluidRow(
						column(width = 2, radioButtons("norm", "Include?", c("yes"="yes","no"="no"))),          
						column(width = 10, offset = 0.3,
							tags$p(align="justify","Intensity normalization with the median peak intensity of each file. Only useful if overall intensities can be assumed to be constant among files.
							Should be avoided if this assumption is violated, e.g., when including blank/blind files into profiling. Better use the below normalization with internal standards."),
							HTML('<p><a href="http://www.looscomputing.ch/eng/enviMass/topics/median_normalization.htm" style="color:rgb(60, 100, 60); text-decoration: none"; target="_blank"><p align="right">&#8594; More info.</a></p>')	
						)
					),
				#HTML('<p style="background-color:darkgreen"; align="center"> <font color="#FFFFFF"> RT alignment </font></p> '),
				#radioButtons("RTalign", "Include?", c("yes"="yes","no"="no")),  
				HTML('<p style="background-color:darkblue"; align="center"> <font color="#FFFFFF"> Blank / blind peak detection </font></p> '),				
					fluidRow(
						column(width = 2, radioButtons("blind", "Detect?", c("yes" = "yes","no" = "no")) ),
						column(width = 2, radioButtons("blind_omit", "Remove?", c("yes" = "TRUE","no" = "FALSE")) ),
						column(width = 8, offset = 0.3,
							tags$p(align="justify","Tags sample peaks which are also present in blind/blank files. Check Settings -> Blind tab for selecting the latter. 
							By also enabling 'Remove?', affected peaks are fully excluded (not just tagged) from being assorted into profiles (but not from, e.g., nontargeted componentization), 
							which may be useful for e.g. its effect on compound screening (relies on profiling) - otherwise a later removal of such peaks is recommended (cp. red steps) ."),
							HTML('<p><a href="http://www.looscomputing.ch/eng/enviMass/topics/blind.htm" style="color:rgb(60, 100, 60); text-decoration: none"; target="_blank"><p align="right">&#8594; More info.</a></p>')	
						)
					),				
				HTML('<p style="background-color:darkblue"; align="center"> <font color="#FFFFFF"> Replicate filter </font></p> '),				
					fluidRow(
						column(width = 2, radioButtons("replicates", "Include?", c("yes"="yes","no"="no")) ),
						column(width = 2, radioButtons("replicates_prof", "Use in profiling?", c("yes"="yes","no"="no")) ),
						column(width = 8, offset = 0.3,
							tags$p(align="justify","Filters out picked peaks which are not ubiquitously present in a set of measurements, within tolerances specified in the Settings/Replicates tab.
							Typically, such a set would be composed of replicate measurements. 
							A set can be defined by a joint string other than FALSE (e.g. replicates_A) in the tag_3 column of the file table, as assembled in the Files tab.
							Optionally, profiles can be extracted in the replicates first, and later merged over all files."),
							HTML('<p><a href="http://www.looscomputing.ch/eng/enviMass/topics/replicates.htm" style="color:rgb(60, 100, 60); text-decoration: none"; target="_blank"><p align="right">&#8594; More info.</a></p>')	
						)
					),
				HTML('<hr noshade="noshade" />'),
				HTML('<h1 align="center"> &#x21e9; </h1> '),                     
				# block 3 ######################################################
				HTML('<p style="background-color:darkgreen"; align="center"> <font color="#FFFFFF"> Profile extraction </font></p>'),
					fluidRow(
						column(width = 2, radioButtons("profiling", "Include? ", c("yes"="yes","no"="no")) ),
						column(width = 10, offset = 0.3,
							tags$p(align="justify","An intensity descent assorts peaks into profiles. Uses a fixed retention time tolerance window and
							an adaptive mass tolerance window. Only peaks of files marked for profiling are used (column profiled in the files table / File Tab)."),
							HTML('<p><a href="http://www.looscomputing.ch/eng/enviMass/topics/profiling.htm" style="color:rgb(60, 100, 60); text-decoration: none"; target="_blank"><p align="right">&#8594; More info.</a></p>')	
						)
					),					
				HTML('<p style="background-color:darkgreen"; align="center"> <font color="#FFFFFF"> LOD interpolation </font></p>'),
					fluidRow(
						column(width = 2, radioButtons("LOD", "Include? ", c("yes"="yes","no"="no")) ),
						column(width = 10, offset = 0.3,
							tags$p(align="justify","For each measurement, estimates a RT-dependent intensity threshold below which peaks are not expected to get picked. 
							Can be used for the below compound screening."),
							HTML('<p><a href="http://www.looscomputing.ch/eng/enviMass/topics/lod.htm" style="color:rgb(60, 100, 60); text-decoration: none"; target="_blank"><p align="right">&#8594; More info.</a></p>')	
						)
					),
				HTML('<p style="background-color:darkgreen"; align="center"> <font color="#FFFFFF"> Compound screening </font></p> '),
					fluidRow(
						column(width = 3, 
							radioButtons("IS_screen", "Screen internal standards?", c("yes"="yes","no"="no"))
						),
						column(width = 3, 
							radioButtons("target_screen", "Screen target compounds?", c("yes"="yes","no"="no"))
						),						
						column(width = 6, offset = 0.3,
							tags$p(align="justify","Uses the LOD thresholds estimated in the above step. If the LOD interpolation is not run, a fixed intensity threshold as specified in the
							Settings/Screening tabs is used. All peaks of files marked for inclusion and profiling are used (column profiled in the files table / File Tab) - unless peaks are
							removed explicitly in the early blind detection (see blue workflow steps)."),
							HTML('<p><a href="http://www.looscomputing.ch/eng/enviMass/topics/screening.htm" style="color:rgb(60, 100, 60); text-decoration: none"; target="_blank"><p align="right">&#8594; More info.</a></p>')	
						)
					),	
				HTML('<p style="background-color:darkgreen"; align="center"> <font color="#FFFFFF"> Intensity normalization using IS-profiles </font></p> '),
					fluidRow(
						column(width = 2, radioButtons("IS_normaliz", "Include? ", c("yes"="yes","no"="no")) ),
						column(width = 10, offset = 0.3,
							tags$p(align="justify","Relies on the above profile extraction and internal standard screening (green steps). Intensities of picked peaks in each measurement 
							are normalized by the median deviation all internal standards have in the measurement from their individual median profile intensity.
							This approach can replace the less reliable Median intensity normalization above (blue steps). Internal standards must have been spiked at constant concentrations."),
							HTML('<p><a href="http://www.looscomputing.ch/eng/enviMass/topics/IS_normalization" style="color:rgb(60, 100, 60); text-decoration: none"; target="_blank"><p align="right">&#8594; More info.</a></p>')	
						)
					),	
				HTML('<hr noshade="noshade" />'),
				HTML('<h1 align="center"> &#x21e9; </h1> '),  					
				# block 4 ######################################################					
				HTML('<p style="background-color:orange"; align="center"> <font color="#FFFFFF"> Calibration </font></p> '),
					fluidRow(
						column(width = 2, radioButtons("calibration", "Include? ", c("yes"="yes","no"="no")) ),
						column(width = 10, offset = 0.3,
							tags$p(align="justify","This step screens for calibration sets of target and internal standard compound peaks, using
							the provided calibration files. The sets can be used in the Calibration tab to establish specific calibration models (curves)
							for quantification. If selected, the extraction of these calibration peaks will be affected by the above mass recalibration, 
							replicate intersection, blind subtraction and LOD interpolation steps. Once you have established the desired calibration models,
							you can deselect this step, until new calibration files have been added."),
							HTML('<p><a href="http://www.looscomputing.ch/eng/enviMass/topics/quantification.htm" style="color:rgb(60, 100, 60); text-decoration: none"; target="_blank"><p align="right">&#8594; More info.</a></p>')	
						)
					),
				HTML('<p style="background-color:orange"; align="center"> <font color="#FFFFFF"> Quantification </font></p> '),
					fluidRow(
						column(width = 2, radioButtons("quantification", "Include? ", c("yes"="yes","no"="no")) ),
						column(width = 10, offset = 0.3,
							tags$p(align="justify","Based on the calibration models from the above step, an estimation of target compound concentrations
							from their intensity ratios to their individual internal standard compounds is derived."),
							HTML('<a href="http://www.looscomputing.ch/eng/enviMass/topics/quantification.htm" style="color:rgb(60, 100, 60); text-decoration: none"; target="_blank"><p align="right">&#8594; More info.</a>')	
						)
					),
				HTML('<p style="background-color:orange"; align="center"> <font color="#FFFFFF"> Recovery </font></p> '),
					fluidRow(
						column(width = 2, radioButtons("recovery", "Include? ", c("yes"="yes","no"="no")) ),
						column(width = 10, offset = 0.3,
							tags$p(align="justify","Following the above calibration and quantification steps, a concentration recovery of spiked target compounds is calculated.
							Requires upload of spiked files."),
							HTML('<a href="http://www.looscomputing.ch/eng/enviMass/topics/quantification.htm" style="color:rgb(60, 100, 60); text-decoration: none"; target="_blank"><p align="right">&#8594; More info.</a>')	
						)
					),
				HTML('<hr noshade="noshade" />'),
				HTML('<h1 align="center"> &#x21e9; </h1> '),  		
				# block X ######################################################					
				HTML('<p style="background-color:black"; align="center"> <font color="#FFFFFF"> EIC correlation </font></p> '),
					fluidRow(
						column(width = 2, radioButtons("EIC_correlation", "Include?", c("yes"="yes","no"="no"))),
						column(width = 10, offset = 0.3, tags$p(align="justify","Checks for correlation between EICs (= chromatographic peak shapes) which share enough scans. 
						Also used to exclude non-correlated EIC peaks from the below isotopologue and adduct grouping.")),
						HTML('<p><a href="http://www.looscomputing.ch/eng/enviMass/topics/eic" style="color:rgb(60, 100, 60); text-decoration: none"; target="_blank"><p align="right">&#8594; More info.</a></p>')	
					),	
				HTML('<p style="background-color:black"; align="center"> <font color="#FFFFFF"> File-wise componentization </font></p> '),
					fluidRow(
						column(width = 2, offset = 0.3, tags$h5("Include...")),
						column(width = 3, 
							radioButtons("isotopologues", "... isotopologue grouping?", c("yes"="yes","no"="no"),inline=TRUE),
							tags$h5("(Requires availability of feasible isotopologue space.)")
						),
						column(width = 3, radioButtons("adducts", "... adduct grouping?", c("yes"="yes","no"="no"),inline=TRUE)),
						column(width = 3, radioButtons("homologues", "... homologue series detection?", c("yes"="yes","no"="no"),inline=TRUE))
					),		
					HTML('<p><a href="http://www.looscomputing.ch/eng/enviMass/topics/components" style="color:rgb(60, 100, 60); text-decoration: none"; 
						target="_blank"><p align="right">&#8594; More info.</a></p>'),	
				HTML('<hr noshade="noshade" />'),
				HTML('<h1 align="center"> &#x21e9; </h1> '),  				
				# block 5 ######################################################				
				HTML('<p style="background-color:darkred"; align="center"> <font color="#FFFFFF"> Profile filtering </font></p> '),
					fluidRow(
						column(width = 2, 
							radioButtons("subtr", "Include?", c("yes"="yes","no"="no"),inline=TRUE)
						),
						column(width = 10, offset = 0.3,
							tags$p(align="justify","Run a profile recalculation omitting..."),
							HTML('<p><a href="http://www.looscomputing.ch/eng/enviMass/topics/filtering" style="color:rgb(60, 100, 60); text-decoration: none"; target="_blank"><p align="right">&#8594; More info.</a></p>')				
						)
					),	
					HTML('<hr noshade="noshade" />'), 
					fluidRow(
						column(width = 6, offset = 0.3,
							tags$p(align="justify","... peaks belonging to matches >= the cutoff score defined in the Settings/Screening Tab.
							The original compound screening results remain visible, but affected peaks will not be part of the final profiles.")
						),
						column(width = 3, 
							radioButtons("subtr_IS", "Subtract internal standard peaks?", c("yes"="yes","no"="no"),inline=TRUE)
						),
						column(width = 3, 
							radioButtons("subtr_target", "Subtract target compound peaks?", c("yes"="yes","no"="no"),inline=TRUE)
						)
					),	
					HTML('<hr noshade="noshade" />') ,
					fluidRow(
						column(width = 10, offset = 0.3,
							tags$p(align="justify","... peaks which have also been detected in blind/blank samples but not removed yet (cp. blind detection of blue steps)?")
						),
						column(width = 2, 
							radioButtons("subtr_blind", "Subtract?", c("yes"="yes","no"="no"),inline=TRUE)
						)
					),	
					HTML('<hr noshade="noshade" />') ,
					fluidRow(
						column(width = 10, offset = 0.3,
							tags$p(align="justify","... peaks from spiked files?")
						),
						column(width = 2, 
							radioButtons("subtr_spiked", "Subtract?", c("yes"="yes","no"="no"),inline=TRUE)
						)
					),	
				HTML('<p style="background-color:darkred"; align="center"> <font color="#FFFFFF"> Profile blind detection </font></p> '),
					fluidRow(
						column(width = 2, radioButtons("profblind", "Include? ", c("yes"="yes","no"="no")) ),
						column(width = 10, offset = 0.3,
							tags$p(align="justify","Calculate median intensity ratio of sample vs. blank/blind peaks across each profile?"),
							HTML('<p><a href="http://www.looscomputing.ch/eng/enviMass/topics/blind.htm" style="color:rgb(60, 100, 60); text-decoration: none"; target="_blank"><p align="right">&#8594; More info.</a></p>')	
						)
					),						
				HTML('<p style="background-color:darkred"; align="center"> <font color="#FFFFFF"> Trend detection </font></p> '),
					fluidRow(
						column(width = 2, radioButtons("trendblind", "Include? ", c("yes"="yes","no"="no")) ),
						column(width = 10, offset = 0.3,
							tags$p(align="justify","Detects and ranks profile intensities which increase over time.
							Depending on settings, this can contain a separate blind detection step. Herein, intensities of blind/blank peaks are interpolated over the time series.
							This interpolation and subtraction is only applicable if the separate blind filter step is disabled (see above blue steps and the preceding red step)."),
							HTML('<p><a href="http://www.looscomputing.ch/eng/enviMass/topics/trends.htm" style="color:rgb(60, 100, 60); text-decoration: none"; target="_blank"><p align="right">&#8594; More info.</a></p>')	
						)
					),				
				HTML('<hr noshade="noshade" />'),
				HTML('<h1 align="center"> &#x21e9; </h1> '),					
				HTML('<p style="background-color:black"; align="center"> <font color="#FFFFFF"> Profile componentization </font></p> '),
					fluidRow(
						column(width = 2, radioButtons("components_profiles", "Include?", c("yes"="yes","no"="no"))),
						column(width = 10, offset = 0.3,
							tags$p(align="justify","Aggregate and filter filewise componentization results across profiles; group profiles with similar intensity patterns."),
							HTML('<p><a href="http://www.looscomputing.ch/eng/enviMass/topics/profile_components.htm" style="color:rgb(60, 100, 60); text-decoration: none"; target="_blank"><p align="right">&#8594; More info.</a></p>')	
						)
					),
				HTML('<hr noshade="noshade" />')
				################################################################
	
        ),
        ########################################################################
        # PARAMETER SETTINGS ###################################################
        ########################################################################
        tabPanel("Settings",     
		  tags$h5("Apply parameter settings to project?"), 
		  bsButton("savepar","Apply",style="warning"),
		  bsAlert("alert_2"),
		  HTML('<p><a href="http://www.looscomputing.ch/eng/enviMass/inputs/parameters.htm" style="color:rgb(60, 100, 60); text-decoration: none"; target="_blank"><p align="left">&#8594; Help on how to specify parameters.</a></p>'),	
          tabsetPanel(
            # PEAK PICKING #####################################################
            tabPanel("Peak picking",
  				HTML('<p><a href="http://www.looscomputing.ch/eng/enviMass/topics/peakpicking.htm" style="color:rgb(60, 100, 60); text-decoration: none"; target="_blank"><p align="left">&#8594; Check help for details & parameter descriptions.</a></p>'),			
				HTML('<hr noshade="noshade" />'),
					tags$h4("Data filtering"),
					fluidRow(
						column(width = 5, 
							selectInput("cut_RT", "Filter RT range?", choices = c("yes"="TRUE","no"="FALSE"), selected="FALSE"),
							numericInput("cut_RT_min", "Lower RT bound [minutes]", 0),
							numericInput("cut_RT_max", "Upper RT bound [minutes]", 25)
						),
						column(width = 5, 						
							selectInput("cut_mass", "Filter mass range?", choices = c("yes"="TRUE","no"="FALSE"), selected="FALSE"),
							numericInput("cut_mass_min", "Lower m/z bound", 0),
							numericInput("cut_mass_max", "Upper m/z bound", 2000)
						)
					),
				HTML('<hr noshade="noshade" />'),
				div(style = widget_style10,
					tags$h4("Parameter estimation"),
					tags$h5("This self-parametrization step estimates file-specific values for the below 'Maximum m/z deviation of a centroid data point from its EIC mean',
						the 'Minimum log10(intensity) threshold' and the 'Percentage of low-intense data points to discard'. It overwrites the provided values of these
						three parameters; if the estimation fails the provided values will be used instead. 
						The estimation is useful for peak-picking in files with large but variable amounts of low-intense noise."),
					selectInput("peak_estimate", "Include estimation?", choices = c("yes"="TRUE","no"="FALSE"), selected="TRUE")
				),
				div(style = widget_style2,
					tags$h4("Extraction of ion chromatogram (EIC)"),
					sliderInput("peak_drtgap", "Maximum retention time gap in an EIC [seconds]", min = 20, max = 1500, value = 300, step= 1),
					sliderInput("peak_dmzdens", "Maximum m/z deviation of a centroid data point from its EIC mean [ppm]", min = 1, max = 100, value = 3.5, step= 0.1)       				
				),
				div(style = widget_style2,
					tags$h4("Peak picking"), 
					numericInput("peak_minpeak", "Minimum number of centroid data points  per peak ...", 4),
					sliderInput("peak_drtsmall2", "... within a given RT window [s]", min = 1, max = 200, value = 20, step= 0.1),
					sliderInput("peak_drtfill", "Maximum RT gap length to be interpolated [s]", min = 0, max = 60, value = 10, step= 0.1),
					sliderInput("peak_drtdens2", "Maximum RT width of a single peak [s]", min = 10, max = 1500, value = 120, step= 0.1),
					HTML('<hr noshade="noshade" />') ,
					numericInput("peak_minint_log10", "Minimum log10(intensity) threshold", -10),
					numericInput("peak_SN", "Minimum Signal/Noise", 5),
					numericInput("peak_SB", "Minimum Signal/Base", 2),
					numericInput("peak_recurs", "Maximum possible number of peaks within a single EIC", 3),
					HTML('<hr noshade="noshade" />') ,
					selectInput("peak_which_intensity", "Peak intensity: use peak area or peak intensoid?", choices = c("intensoid (max int.)"="maximum","area (sum int.)"="area"), selected="intensoid"),
					selectInput("peak_get_mass", "Peak mass definition:", choices = c("Mean"="mean", "Weighted mean"="wmean"), selected="mean")					
				),
				div(style = widget_style11,
					tags$h4("Advanced settings"),
					numericInput("peak_maxint_log10", "Upper log10(intensity) safety threshold", 6.7),	
					numericInput("peak_ended", "How often can a peak detection fail to end the recursion? - peak picking", 1),
					numericInput("peak_weight", "Weight for assigning centroid data points to a peak - peak picking", 1),				
					sliderInput("peak_perc_cut", "Percentage of low-intense data points to discard", min = 0, max = 100, value = 0, step= .1)
				),
              tags$h4(""),
              tags$h4("")              
            ),
			# RESOLUTION #######################################################
            tabPanel("Instrument / Resolution",
                  div(style = widget_style3,selectInput("resolution", label="Instrument resolution:", choices="none", selected="none", width='600px')),
				  imageOutput("plot_resolution", height="auto")	  
            ),	
            # RECALIBRATION ####################################################
            tabPanel("Mass recalibration",
      		  HTML('<p><a href="http://www.looscomputing.ch/eng/enviMass/topics/recalibration.htm" style="color:rgb(60, 100, 60); text-decoration: none"; target="_blank"><p align="left">&#8594; Check help for details & parameter descriptions.</a></p>'),				
              HTML('<hr noshade="noshade" />'),
              div(style = widget_style2,
				tags$h4("Positive ionization:"), HTML('<hr noshade="noshade" />'),
                selectInput("recal_include_pos", "Include mass recalibration for positive ion. mode files?", c("Yes"="TRUE","No"="FALSE"),"TRUE",multiple=FALSE),                
                selectInput("recal_use_pos", "Reference compounds:", c("Internal standards","Target compounds","both"),"Internal standards",multiple=FALSE),                
                numericInput("recal_dmz_pos", "+/- m/z tolerance ...", 3),            
                numericInput("recal_maxdmz_pos", "Maximum allowable m/z correction ...", 30),  				
                selectInput("recal_ppm_pos", "... given in:", choices = c("ppm"="TRUE","absolute [mmu]"="FALSE"), "TRUE"),
                numericInput("recal_drt_pos", "RT tolerance [s]", 30)   
              ),
              div(style = widget_style2,
				tags$h4("Negative ionization:"), HTML('<hr noshade="noshade" />'),
                selectInput("recal_include_neg", "Include mass recalibration for negative ion. mode files?", c("Yes"="TRUE","No"="FALSE"),"TRUE",multiple=FALSE),                
                selectInput("recal_use_neg", "Reference compounds:", c("Internal standards","Target compounds","both"),"Internal standards",multiple=FALSE),                
                numericInput("recal_dmz_neg", "+/- m/z tolerance ...", 3),            
                numericInput("recal_maxdmz_neg", "Maximum allowable m/z correction ...", 30),  				
                selectInput("recal_ppm_neg", "... given in:", choices = c("ppm"="TRUE","absolute [mmu]"="FALSE"), "TRUE"),
                numericInput("recal_drt_neg", "RT tolerance [s]", 30)   
              )		  
            ),
			# REPLICATES #######################################################
            tabPanel("Replicates",
      		  	HTML('<p><a href="http://www.looscomputing.ch/eng/enviMass/topics/replicates.htm" style="color:rgb(60, 100, 60); text-decoration: none"; target="_blank"><p align="left">&#8594; Check help for details & parameter descriptions.</a></p>'),				
				HTML('<hr noshade="noshade" />'),
				tags$h5("Replicate files are defined (i.e., grouped) by the tag3 entry (not FALSE) in the measurements table"),
				numericInput("replicate_dmz", "+/- m/z tolerance ...", 3),                
				selectInput("replicate_ppm", "... given in:", choices = c("ppm"="TRUE","absolute [mmu]"="FALSE"), "TRUE"),	
				#selectInput("replicate_recalib", "... and corrected by recalibration results (if available)", choices = c("TRUE"="TRUE","FALSE"="FALSE"), "FALSE"),	
				numericInput("replicate_delRT", "RT tolerance window of peaks caused by the same analyte across replicate samples [s]", 30),
				numericInput("replicate_IS_dInt", "Absolute log intensity tolerance X (log 10 scale, 10^X):", 10)
			),	
            # ALLIGNMENT #######################################################
            #tabPanel("Alignment",
            #  tags$h5("RT alignment")
            #),
            # SCREENING ########################################################
            tabPanel("Screening",
              HTML('<p><a href="http://www.looscomputing.ch/eng/enviMass/topics/screening.htm" style="color:rgb(60, 100, 60); text-decoration: none"; target="_blank"><p align="left">&#8594; Check help for details & parameter descriptions.</a></p>'),				
              HTML('<hr noshade="noshade" />'),
              tabsetPanel(
                tabPanel("IS",
					div(style = widget_style2,
						tags$h5("Retention time"), HTML('<hr noshade="noshade" />'),
						numericInput("IS_drt1", "RT tolerance of peaks relative to their expected RT [s]", 30),
						numericInput("IS_drt2", "RT tolerance of peaks within an isotope pattern [s]", 50)
					),
				  	div(style = widget_style2,
						tags$h5("Mass"), HTML('<hr noshade="noshade" />'),
						numericInput("IS_dmz", "+/- m/z tolerance ...", 3),                
						selectInput("IS_ppm", "... given in:", choices = c("ppm"="TRUE","absolute [mmu]"="FALSE"), "TRUE")
 					),
					div(style = widget_style2,
						tags$h5("Intensity"), HTML('<hr noshade="noshade" />'),
						#sliderInput("IS_inttol", "Intensity tolerance %", min = 0, max = 100, value = 30, step= .2),
						knobInput("IS_inttol", label = "Intensity tolerance %", min = 0, max = 100, value = 30, step= .2, 
							thickness = .3, width = 160, height = 160, fgColor = "#7499AB", inputColor = "#428BCA", displayPrevious = TRUE),
						numericInput("IS_intcut", "Lower intensity threshold (if LOD interpolation disabled)", 5E4)                
					),
					div(style = widget_style2,
						tags$h5("Restrict screening to latest files (covered during profiling)?"), HTML('<hr noshade="noshade" />'),
						selectInput("screen_IS_restrict", "Include?", choices = c("TRUE"="TRUE","FALSE"="FALSE"), "FALSE"),
						numericInput("screen_IS_restrict_many", "Number of latest files to include", 10)
					),
					div(style = widget_style2,
						tags$h5("Scoring"), HTML('<hr noshade="noshade" />'),
						numericInput("IS_w1", "Cutoff score [0,1]", 0.8),         
						HTML('<h1 align="center"> &#x21f3; </h1> '),						
						selectInput("screen_IS_cutit", "Exclude matches below cutoff score?", choices = c("TRUE"="TRUE","FALSE"="FALSE"), "FALSE"),
						HTML('<h1 align="center"> &#x21f3; </h1> '),
						selectInput("screen_IS_maxonly", "Screen only most intense isotopologue peak?", choices = c("TRUE"="TRUE","FALSE"="FALSE"), "FALSE")
					)
                ),
                tabPanel("Targets & Suspects",
					div(style = widget_style2,
						tags$h5("Retention time"), HTML('<hr noshade="noshade" />'),
						numericInput("tar_drt1", "RT tolerance of peaks relative to their expected RT [s]", 30),
						numericInput("tar_drt2", "RT tolerance of peaks within an isotope pattern [s]", 50)
					),
				  	div(style = widget_style2,
						tags$h5("Mass"), HTML('<hr noshade="noshade" />'),
						numericInput("tar_dmz", "+/- m/z tolerance ...", 3),                
						selectInput("tar_ppm", "... given in:", choices = c("ppm"="TRUE","absolute [mmu]"="FALSE"), "TRUE")
 					),
					div(style = widget_style2,
						tags$h5("Intensity"), HTML('<hr noshade="noshade" />'),
						#sliderInput("tar_inttol", "Intensity tolerance %", min = 0, max = 100, value = 30, step= .2),
						knobInput("tar_inttol", label = "Intensity tolerance %", min = 0, max = 100, value = 30, step= .2, 
							thickness = .3, width = 160, height = 160, fgColor = "#7499AB", inputColor = "#428BCA", displayPrevious = TRUE),
						numericInput("tar_intcut", "Lower intensity threshold (if LOD interpolation disabled)", 5E4)                
					),
					div(style = widget_style2,
						tags$h5("Restrict screening to latest files (covered during profiling)?"), HTML('<hr noshade="noshade" />'),
						selectInput("screen_target_restrict", "Include?", choices = c("TRUE"="TRUE","FALSE"="FALSE"), "FALSE"),
						numericInput("screen_target_restrict_many", "Number of latest files to include", 10)
					),
					div(style = widget_style2,
						tags$h5("Scoring"), HTML('<hr noshade="noshade" />'),
						numericInput("tar_w1", "Cutoff score [0,1]", 0.8),
						HTML('<h1 align="center"> &#x21f3; </h1> '),
						selectInput("screen_target_cutit", "Exclude matches below cutoff score?", choices = c("TRUE"="TRUE","FALSE"="FALSE"), "FALSE"),
						HTML('<h1 align="center"> &#x21f3; </h1> '),
						selectInput("screen_target_maxonly", "Screen only most intense isotopologue peak?", choices = c("TRUE"="TRUE","FALSE"="FALSE"), "FALSE")					
					)
                ),
				tabPanel("Adducts",
					  helpText("Select adducts for the calculation of centroid masses of IS & targets isotope patterns.\n Used for recalibration and screening."),
					  div(style = widget_style3,checkboxGroupInput("adducts_pos", "Positive ions:", "none")),
					  div(style = widget_style3,checkboxGroupInput("adducts_neg", "Negative ions:", "none"))
				)
              )
            ),
			# QUANTIFICATION & RECOVERY ########################################
			tabPanel("Quantification",
				HTML('<p><a href="http://www.looscomputing.ch/eng/enviMass/topics/quantification.htm" style="color:rgb(60, 100, 60); text-decoration: none"; target="_blank"><p align="left">&#8594; Check help for details & parameter descriptions.</a></p>'),	
				HTML('<hr noshade="noshade" />'),
				numericInput("quant_files_included", "Number of latest file to include in the quantification:", 30),
				numericInput("recov_files_included", "Number of latest spiked files to include for recovery estimation:", 30),
				numericInput("quant_digits", "Decimal places / digits shown for concentrations:", value= 2, min = 0, max = 10, step = 1)
			),
            # INTENSITY NORMALIZATION ##########################################
            tabPanel("Normalization",
				HTML('<p><a href="http://www.looscomputing.ch/eng/enviMass/topics/quantification.htm" style="color:rgb(60, 100, 60); text-decoration: none"; target="_blank"><p align="left">&#8594; Check help for details & parameter descriptions.</a></p>'),	
				HTML('<hr noshade="noshade" />'),
				tags$h5("Intensity normalization based on IS-profiles"),
				div(style = widget_style2,
					tags$h4("Positive ionization"), HTML('<hr noshade="noshade" />'),
					selectInput("ISnorm_include_pos", "Include normalization for positive ion. mode files?", c("Yes"="TRUE","No"="FALSE"),"TRUE",multiple=FALSE), 
					#sliderInput("ISnorm_percfiles_pos", "Minimum of screened files covered by each IS profile (%)", min = 0, max = 100, value = 90, step= 1),
					knobInput("ISnorm_percfiles_pos", label = "Minimum of screened files covered by each IS profile (%)", min = 0, max = 100, value = 90, step= 1, 
							thickness = .3, width = 160, height = 160, fgColor = "#7499AB", inputColor = "#428BCA", displayPrevious = TRUE),
					#sliderInput("ISnorm_score_pos", "Screening threshold", min = 0, max = 1, value = 0.8, step= .01),
					knobInput("ISnorm_score_pos", label = "Screening threshold", min = 0, max = 1, value = 0.8, step= .01, 
							thickness = .3, width = 160, height = 160, fgColor = "#7499AB", inputColor = "#428BCA", displayPrevious = TRUE),
					numericInput("ISnorm_numbIS_pos", "Minimum number of IS profile peaks", 15),
					HTML('<hr noshade="noshade" />'),
					checkboxInput("ISnorm_medblank_pos", "Show median deviation of blank/blind profiles?", FALSE),
					checkboxInput("ISnorm_usesubblank_pos", "Use subsampling", FALSE),
					numericInput("ISnorm_numblank_pos", "Number of blank/blind profiles in subsample", 100),
					HTML('<hr noshade="noshade" />'),
					checkboxInput("ISnorm_medsam_pos", "Show median deviation of sample (i.e., non-blank) profiles?", FALSE),	
					checkboxInput("ISnorm_usesubsam_pos", "Use subsampling", FALSE),
					numericInput("ISnorm_numsam_pos", "Number of sample profiles in subsample", 100)
				),
				div(style = widget_style2,
					tags$h4("Negative ionization"), HTML('<hr noshade="noshade" />'),
					selectInput("ISnorm_include_neg", "Include normalization for positive ion. mode files?", c("Yes"="TRUE","No"="FALSE"),"TRUE",multiple=FALSE), 
					#sliderInput("ISnorm_percfiles_neg", "Minimum of screened files covered by each IS profile (%)", min = 0, max = 100, value = 90, step= 1),
					knobInput("ISnorm_percfiles_neg", label = "Minimum of screened files covered by each IS profile (%)", min = 0, max = 100, value = 90, step= 1, 
							thickness = .3, width = 160, height = 160, fgColor = "#7499AB", inputColor = "#428BCA", displayPrevious = TRUE),
					#sliderInput("ISnorm_score_neg", "Screening threshold", min = 0, max = 1, value = 0.8, step= .01),
					knobInput("ISnorm_score_neg", label = "Screening threshold", min = 0, max = 1, value = 0.8, step= .01, 
							thickness = .3, width = 160, height = 160, fgColor = "#7499AB", inputColor = "#428BCA", displayPrevious = TRUE),
					numericInput("ISnorm_numbIS_neg", "Minimum number of IS profile peaks", 15),
					HTML('<hr noshade="noshade" />'),
					checkboxInput("ISnorm_medblank_neg", "Show median deviation of blank/blind profiles?", FALSE),
					checkboxInput("ISnorm_usesubblank_neg", "Use subsampling", FALSE),
					numericInput("ISnorm_numblank_neg", "Number of blank/blind profiles in subsample", 100),
					HTML('<hr noshade="noshade" />'),
					checkboxInput("ISnorm_medsam_neg", "Show median deviation of sample (i.e., non-blank) profiles?", FALSE),	
					checkboxInput("ISnorm_usesubsam_neg", "Use subsampling", FALSE),
					numericInput("ISnorm_numsam_neg", "Number of sample profiles in subsample", 100)
				)
			),
            # PROFILING ########################################################
            tabPanel("Profiling",
				#tags$h5("Profile extraction"),
				HTML('<p><a href="http://www.looscomputing.ch/eng/enviMass/topics/profiling.htm" style="color:rgb(60, 100, 60); text-decoration: none"; target="_blank"><p align="left">&#8594; Check help for details & parameter descriptions.</a></p>'),	
				HTML('<hr noshade="noshade" />'),
				#sliderInput("prof_maxfiles", "Maximum number of newest (by date/time) files to be processed per ion mode:", min = 50, max = 5000, value = 100, step= 1),
				knobInput("prof_maxfiles", label = "Maximum number of newest (by date/time) files to be processed per ion mode:", min = 50, max = 5000, value = 100, step= 1, 
							thickness = .3, width = 160, height = 160, fgColor = "#7499AB", inputColor = "#428BCA", displayPrevious = TRUE),
				numericInput("prof_dmz", "Peak mass deviation within profiles: +/- m/z tolerance ...", 3),                
                selectInput("prof_ppm", "... given in:", choices = c("ppm"="TRUE","absolute [mmu]"="FALSE"), "TRUE"),
                numericInput("prof_drt", "Peak deviation within profiles: RT tolerance [s]", 60),     
				HTML('<hr noshade="noshade" />'),
				#sliderInput("prof_comp_maxfiles", "Maximum number of file-wise components to be considered for profile componentization (by date/time)", min = 1, max = 300, value = 15, step= 1),
				#HTML('<hr noshade="noshade" />'),				
				selectInput("prof_select", "Omit files with table entry profiled=FALSE from profiling?", choices = c("TRUE"="TRUE","FALSE"="FALSE"), selected="FALSE")			
			),
            # TREND ############################################################
            tabPanel("Trends",
				#tags$h5("Trend detection:"),
				HTML('<p><a href="http://www.looscomputing.ch/eng/enviMass/topics/trends.htm" style="color:rgb(60, 100, 60); text-decoration: none"; target="_blank"><p align="left">&#8594; Check help for details & parameter descriptions.</a></p>'),	
				HTML('<hr noshade="noshade" />'),
				textInput("trend_lags", "Time lags of trends [days], comma-separated:", value = "4,7,14"),
				numericInput("trend_threshold", "Trend vs. mean+variance intensity threshold:", 3),
				HTML('<hr noshade="noshade" />'),
				radioButtons("trend_blind", "Additional blind interpolation and subtraction per profile?", c("yes"="yes","no"="no")),
				HTML('<hr noshade="noshade" />'),
				selectInput("notrend", "Do not show global trend - instead, report it as maximum intensity above blind", choices = c("TRUE"="TRUE","FALSE"="FALSE"), "FALSE")
            ),
            # BLIND #############################################################
            tabPanel("Blind",
				#tags$h5("Blind subtraction:"),
				HTML('<p><a href="http://www.looscomputing.ch/eng/enviMass/topics/blind.htm" style="color:rgb(60, 100, 60); text-decoration: none"; target="_blank"><p align="left">&#8594; Check help for details & parameter descriptions.</a></p>'),	
				HTML('<hr noshade="noshade" />'),
				numericInput("blind_threshold", "Factor by which the sample peak intensity must exceed the blank/blind peak intensity to not be subtracted", 100),
				numericInput("blind_dmz", "+/- m/z tolerance ...", 3), 
                selectInput("blind_ppm", "... given in:", choices = c("ppm"="TRUE","absolute [mmu]"="FALSE"), "TRUE"),				
                numericInput("blind_drt", "RT tolerance [s]:", 60),       
				HTML('<hr noshade="noshade" />'),
				tags$h5("Positive ionization mode"),
				checkboxInput("subtract_pos_bydate", "Subtract with the next blank/blind file preceding a sample by its date & time?", FALSE),				
				checkboxInput("subtract_pos_byfile", "Additional non-sample files to subtract each sample file with (i.e. not preceding by date only), choose file ID:", FALSE),				
				checkboxGroupInput("files_pos_select_subtract", label="", choices=c("FALSE"), selected = NULL),
				HTML('<hr noshade="noshade" />'),
				tags$h5("Negative ionization mode"),
				checkboxInput("subtract_neg_bydate", "Subtract with the next blank/blind file preceding a sample by its date & time?", FALSE),				
				checkboxInput("subtract_neg_byfile", "Additional non-sample files to subtract each sample file with (i.e. not preceding by date only), choose file ID:", FALSE),							
				checkboxGroupInput("files_neg_select_subtract", label="", choices=c("FALSE"), selected = NULL)
            ),			
            # Componentization #################################################
            tabPanel("Componentization",		
				tabsetPanel(			
					tabPanel("File-wise componentization",
						HTML('<p><a href="http://www.looscomputing.ch/eng/enviMass/topics/file_components.htm" style="color:rgb(60, 100, 60); text-decoration: none"; target="_blank"><p align="left">&#8594; Check help on filewise componentization.</a></p>'),	
						HTML('<hr noshade="noshade" />'),
						tabsetPanel(
							tabPanel("Isotopologue grouping",	
								HTML('<p><a href="http://www.looscomputing.ch/eng/enviMass/topics/isotopologues.htm" style="color:rgb(60, 100, 60); text-decoration: none"; target="_blank"><p align="left">&#8594; Check help on isotopologue grouping parameters.</a></p>'),	
								HTML('<hr noshade="noshade" />'),
								numericInput("isotop_mztol", "+/- m/z tolerance ...", 3), 
								selectInput("isotop_ppm", "... given in:", choices = c("ppm"="TRUE","absolute [mmu]"="FALSE"), "TRUE"),				
								numericInput("isotop_rttol", "RT tolerance [s]:", 60),       
								#sliderInput("isotop_inttol", "Intensity tolerance %", min = 0, max = 100, value = 30, step= .2)
								knobInput("isotop_inttol", label = "Intensity tolerance %", min = 0, max = 100, value = 30, step= .2, 
											thickness = .3, width = 160, height = 160, fgColor = "#7499AB", inputColor = "#428BCA", displayPrevious = TRUE),
								HTML('<hr noshade="noshade" />'),
								helpText("The below switch allows to include an estimation of the maximum number of atoms for selected elements for each component. The results will be displayed in additional columns max_atoms_ in the table
								Results -> Grouping -> Components -> Full component list. Estimates are only based on the candidate monoisotopic peak of each component and the possible charges. The estimates are very conservative, i.e., likely too
								high. More fine-tuned estimates for selected peaks or components can be derived under Results -> Grouping -> atom bounds."),
								selectInput("do_atom_bounds_components", "Run atom bound estimation?", choices = c("TRUE" = "TRUE","FALSE" = "FALSE"), "FALSE"),				
								selectInput("atom_bounds_components", label = "Choose elements to include / exclude", choices = c("C","O","N","H","S","P","F","Cl","Br"), selected = c("Cl","Br"), multiple = TRUE)
							),
							tabPanel("Adduct grouping",	
								HTML('<p><a href="http://www.looscomputing.ch/eng/enviMass/topics/adducts.htm" style="color:rgb(60, 100, 60); text-decoration: none"; target="_blank"><p align="left">&#8594; Check help on adduct grouping parameters.</a></p>'),	
								HTML('<hr noshade="noshade" />'),
								numericInput("adducts_mztol", "+/- m/z tolerance ...", 3), 
								selectInput("adducts_ppm", "... given in:", choices = c("ppm"="TRUE","absolute [mmu]"="FALSE"), "TRUE"),				
								numericInput("adducts_rttol", "RT tolerance [s]:", 60),       
								HTML('<hr noshade="noshade" />'),
								div(style = widget_style3,checkboxGroupInput("adducts_pos_group", "Positive mode:", "none")),
								div(style = widget_style3,checkboxGroupInput("adducts_neg_group", "Negative mode:", "none"))						
							),
							tabPanel("Homologue series detection",	
								HTML('<p><a href="http://www.looscomputing.ch/eng/enviMass/topics/homol.htm" style="color:rgb(60, 100, 60); text-decoration: none"; target="_blank"><p align="left">&#8594; Check help on homologue detection parameters.</a></p>'),	
								HTML('<hr noshade="noshade" />'),
								textInput("homol_units", label="Homologue units (comma seperated, no empty spaces; insert FALSE to detect all series with m/z differences in between 10 and 120)", value = "CH2,CH4O", width = '100%', placeholder = "CH2,CH4O"),
								textInput("homol_charges", label="Charges z (comma seperated, no empty spaces)", value = "1,2", width = NULL, placeholder = "1,2"),
								HTML('<hr noshade="noshade" />'),
								numericInput("homol_minrt", "Minimum change in RT from one homologue to the next [s]", 10),       
								numericInput("homol_maxrt", "Maximum change in RT from one homologue to the next [s]", 60),   
								numericInput("homol_rttol", "Tolerance by which the RT differences between two adjacent homologue pairs are allowed to change [s]", 20), 
								HTML('<hr noshade="noshade" />'),
								numericInput("homol_mztol", "+/- m/z tolerance ...", 3), 
								selectInput("homol_ppm", "... given in:", choices = c("ppm"="TRUE","absolute [mmu]"="FALSE"), "TRUE"),				
								numericInput("homol_minlength", "Minimum number of homologues in a series:", 6),						
								HTML('<hr noshade="noshade" />'),						
								div(style = widget_style11,		
									tags$h5("Advanced settings"),
									numericInput("homol_vec_size", "Ignore unless a relevant error message is printed (then try to increase size): ", 1E8),
									HTML('<hr noshade="noshade" />'),
									selectInput("homol_blind", "Filter peaks by sample-vs-blind intensity ratio:", choices = c("yes"="TRUE","no"="FALSE"), "FALSE"),				
									numericInput("homol_blind_value", "Sample-vs-blind intensity ratio threshold:", 10)									
								)
							),
							tabPanel("EIC correlation",	
								HTML('<p><a href="http://www.looscomputing.ch/eng/enviMass/topics/eic.htm" style="color:rgb(60, 100, 60); text-decoration: none"; target="_blank"><p align="left">&#8594; Check help on EIC correlation parameter.</a></p>'),	
								HTML('<hr noshade="noshade" />'),
								numericInput("EICor_delRT", "RT tolerance window for candidate peak pairs [s]:", 5),
								numericInput("EICor_minpeaks", "Min. number of data points (scans) shared by EIC pairs to check for their correlation:", 15),
								numericInput("EICor_mincor", "Min. Pearson correlation [0,1] coefficient:", 0.95)
							)
						), style = "color: #123123;"
					),
					tabPanel("Profile componentization",		
						HTML('<p><a href="http://www.looscomputing.ch/eng/enviMass/topics/profile_components.htm" style="color:rgb(60, 100, 60); text-decoration: none"; target="_blank"><p align="left">&#8594; Check help on profile componentization.</a></p>'),						
						HTML('<hr noshade="noshade" />'), 
						fluidRow(
							column(5,numericInput("corr_min_peaks", "Minimum number of files over which peaks of different profiles have to co-occur to check their intensity correlation:", value = 5, min = 0, max = NA)),
							column(3,numericInput("comp_corr", "Minimum Pearson profile intensity correlation:", value = 0.9, min = 0, max = 1)),
							column(4,numericInput("corr_del_RT", "RT tolerance window for co-occuring peaks of different profiles [s]:", value = 5, min = 0, max = NA))
						),		
						HTML('<hr noshade="noshade" />'), 
						selectInput("prof_comp_link_only",label="Restrict profile componentization to isotopologue and selected adduct relations only?",choices=c("TRUE","FALSE"),selected = "FALSE", multiple = FALSE),
						HTML('<hr noshade="noshade" />'),
						HTML('<hr noshade="noshade" />'),
						tags$h4("File restrictions:"),
						fluidRow(
							column(5, selectInput("dofile_latest_profcomp",label="Restrict profile componentization to a set of latest files only?",choices=c("TRUE","FALSE"),selected = "TRUE", multiple = FALSE)),
							column(4, numericInput("numfile_latest_profcomp", "Number of latest files to include:", 50))
						),						
						HTML('<hr noshade="noshade" />'), 
						tags$h4("Filtering of outliers in profile component relations:"),
						fluidRow(
							column(5, 
								tags$h5("Positive ionization mode"),
								selectInput("filter_profcomp_pos",label="Filter positive mode components?",choices=c("TRUE","FALSE"),selected = "TRUE", multiple = FALSE),
								selectInput("for_which_profcomp_pos",label="Set standard profiles by:",
									choices=c("Internal standard profiles"="ISTD","Target compound profiles"="target","Internal standard and target profiles"="both","All profiles"="all"),selected = "all", multiple = FALSE)
							),			
							column(5, 
								tags$h5("Negative ionization mode"),
								selectInput("filter_profcomp_neg",label="Filter negative mode components?",choices=c("TRUE","FALSE"),selected = "TRUE", multiple = FALSE),
								selectInput("for_which_profcomp_neg",label="Set standard profiles by:",
									choices=c("Internal standard profiles"="ISTD","Target compound profiles"="target","Internal standard and target profiles"="both","All profiles"="all"),selected = "all", multiple = FALSE)							
							)
						),
						HTML('<hr noshade="noshade" />')
					)
				)
            ),			
            # GENERAL SETTINGS #################################################
            tabPanel("General",
				div(style = widget_style10,
					textInput("PWpath", "Path to Proteowizard (PW) MSConvert (use / in Windows, include msconvert.exe)", value = "C:/Program Files/ProteoWizard/ProteoWizard 3.0.5140/msconvert.exe"),
					checkboxInput("files_SIM", "Enable PW SIM sim as spectra?", FALSE)
				),
				div(style = widget_style10,
					tags$h5("Multi-core processing"),
					selectInput("parallel", "Enable parallelized workflow processing?", choices = c("TRUE", "FALSE"), selected = "TRUE"),
					HTML('<hr noshade="noshade" />'), 
					selectInput("parallel_restrict", "Restrict number of cores/threads?", choices = c("TRUE", "FALSE"), selected = "FALSE"),
					numericInput("parallel_cores", "... maximum number of cores/threads to use:", 4)
				),
				div(style = widget_style11,
					tags$h5("Debug tools"),
					selectInput("progressBar", "Show progress bars (Windows OS only)", choices = c("TRUE", "FALSE"), selected = "FALSE"),
					HTML('<hr noshade="noshade" />'), 
					selectInput("do_project_check", "Skip the project check before calculations?", choices = c("TRUE", "FALSE"), selected = "FALSE"),	
					selectInput("verbose", "Verbose console outputs?", choices = c("TRUE", "FALSE"), selected = "FALSE"),	
					selectInput("test", "Include testing during calculations?", choices = c("TRUE", "FALSE"), selected = "FALSE"),						
					selectInput("ignore_large_files", "Ignore .mzXML and MSlist files during check?", choices = c("TRUE", "FALSE"), selected="FALSE"),	
					HTML('<hr noshade="noshade" />'), 				
					textInput("upto_file", "Include latest files only up to file with ID:", value = "FALSE"),
					HTML('<hr noshade="noshade" />'), 
					tags$h6("Reset project excluding peak picking:"),
					bsButton("reset_1","Reset",style="danger"),
					tags$h6("Reset project including peak picking:"),
					bsButton("reset_2","Reset",style="danger")
				)
			),
			# IMPORT PARAMETERS FROM ANOTHER PROJECT ###########################
            tabPanel("Import",			
				tags$h5("Import all parameter settings from another project, excluding blind file selection. Insert full project path to import from, i.e., the path to the project folder in which the logfile.emp file of the project file is located. Does only work between projects made with the same enviMass version."),
						textInput("import_pro_dir_paras", "", value = "C:\\...\\other_project_name"),
						bsPopover("import_pro_dir_paras", 
							title = "Insert full project path to import from, i.e., the the project folder in which a logfile.emp file is located.",
							content = "Using your OS explorer, you may navigate into your project folder and copy/paste the full path.", 
							placement = "right", trigger = "hover"),
						actionButton("Import_project_para","Import")		
			)
          )
        ),
        ########################################################################
        # CALIBRATION ##########################################################
        ########################################################################
        tabPanel("Calibration", 	
			HTML('<a href="http://www.looscomputing.ch/eng/enviMass/topics/quantification.htm" style="color:rgb(60, 100, 60); text-decoration: none"; target="_blank"><p align="left">&#8594; More info.</a>'),
			tabsetPanel( 
				tabPanel("Create",
					helpText("Select the ionization mode to load the available calibration file sets into the below selection."),
					selectInput("Ion_mode_Cal", label="Ionization mode", c("none","positive","negative"), selected = ("none"), multiple = FALSE),	
					HTML('<hr noshade="noshade" />'),
					HTML('<h1 align="center"> &#x21e9; </h1> '),
					conditionalPanel(
						condition = "input.Ion_mode_Cal != 'none'", 	
						selectInput("Cal_file_set",label="Specify calibration file set",choices=c("none"),selected = "none", multiple = FALSE),
						conditionalPanel(
							condition = "input.Cal_file_set != 'none' & input.Ion_mode_Cal != 'none'", 						
							bsButton("Cal_file_set_delete", label="Delete all models for this calibration file set?", icon = icon("fire-extinguisher"))
						),
						HTML('<hr noshade="noshade" />'),
						HTML('<h1 align="center"> &#x21e9; </h1> '),
						conditionalPanel(
							condition = "input.Cal_file_set != 'none' & input.Ion_mode_Cal != 'none'", 
							helpText("Select Target compounds which are linked to an Internal standard (column 6 of the target compound table) to (resume) work on their individual calibration models below. 
							The screened compounds can also be viewed in the Results/Compound screening tab, choosing the calibration files. Target - Internal standard links can be modified either directly 
							in the target compound table (above tab Compound -> Targets) or with the below orange button Use for quantification"),
							fluidRow(
								column(3, selectInput("Cal_target_name",label="Target name",choices=c("none"),selected = "none", multiple = FALSE)),							
								column(3, selectInput("Cal_target_ID",label="Target ID",choices=c("none"),selected = "none", multiple = FALSE))
							),
							fluidRow(
								column(3, selectInput("Cal_IS_name",label="Internal standard name",choices=c("none"),selected = "none", multiple = FALSE)),
								column(3, selectInput("Cal_IS_ID",label="Internal standard ID",choices=c("none"),selected = "none", multiple = FALSE))
							),
							textOutput('number_missing_models'),
							HTML('<hr noshade="noshade" />'),			
							helpText("Proceed through entries in the target compound list:"),
							fluidRow(
								column(1, bsButton("Cal_first", label="", icon = icon("step-backward"))),
									bsTooltip("Cal_first", title="Jump to first Target", placement = "top", trigger = "hover", options = NULL),	
								column(1, bsButton("Cal_previous", label="", icon = icon("arrow-left"))),
									bsTooltip("Cal_previous", title="Jump to preceeding (or last, if beginning of table is reached) Target", placement = "top", trigger = "hover", options = NULL),
								column(1, bsButton("Cal_next", label="", icon = icon("arrow-right"))),		
									bsTooltip("Cal_next", title="Jump to next (or first, if end of table is reached) Target", placement = "top", trigger = "hover", options = NULL),	
								column(1, bsButton("Cal_next_missing", label="", icon = icon("arrow-right"), style="info")),		
									bsTooltip("Cal_next_missing", title="Jump to next Target with missing calibration model", placement = "top", trigger = "hover", options = NULL),	
								column(1, bsButton("Cal_missing", label="", icon = icon("map-pin"), style="info")),	
									bsTooltip("Cal_missing", title="Jump to first Target with missing calibration model", placement = "top", trigger = "hover", options = NULL),
								column(1, bsButton("Cal_last", label="", icon = icon("step-forward"))),
									bsTooltip("Cal_last", title="Jump to last Target", placement = "top", trigger = "hover", options = NULL)							
							),							
							HTML('<hr noshade="noshade" />'),
							HTML('<h1 align="center"> &#x21e9; </h1> '),
							bsButton("save_Cal","Save/replace model",style="success"),
							bsPopover("save_Cal", 
								title = "Save the selected calibration model for later modification or to quantify with? Overwrites any existing model.",
								content = "The model is only used for quantification if the Internal standard is linked to the selected Target compound in column 6 of the target compound table. You may modify this default link with the orange button on the right.", 
								placement = "bottom", trigger = "hover"),							
							bsButton("use_Cal","Use for quantification",style="warning"),
							bsPopover("use_Cal", 
								title = "Make entry to target compound list",
								content = "Save the selected Internal standard as the one to quantify with (= entry into column 6 of the target compound table)? NOTE: this does not save the model, but only the relevant relation between Target and Internal standard compound. Use the green button to save the model.", 
								placement = "bottom", trigger = "hover"),							
							bsButton("remove_Cal","Remove model",style="danger"),
							bsPopover("remove_Cal", 
								title = "Delete any existing model for the compound selection.",
								content = "Target compounds without a calibration model for their linked Internal standard  (= entry into column 6 of the target compound table) will not be quantified.", 
								placement = "bottom", trigger = "hover"),							
							bsButton("reload_Cal","Reload data",style="info"),	
							bsPopover("reload_Cal", 
								title = "Refresh below plot and table for the selected compound pair of Target and Internal standard ...",
								content = "... e.g., after you have changed the intensity bounds of the Internal standard in the tab <em> Results -> Screening& -> Ionization -> Internal standards (show all adducts) -> Characteristics for selected compound. </em>", 
								placement = "bottom", trigger = "hover"),							
							conditionalPanel(
								condition = "input.Cal_file_set != 'none' & input.Ion_mode_Cal != 'none' & input.Cal_target_ID != 'none' & input.Cal_target_name != 'none' & input.Cal_IS_ID != 'none' & input.Cal_IS_name != 'none'  ", 
								HTML('<hr noshade="noshade" />'),
								div(style = widget_style6,
									plotOutput("cal_plot", 
										dblclick = "cal_plot_dblclick",
										brush = brushOpts(
										  id = "cal_plot_brush",
										  resetOnNew = TRUE
										)
									)
								),
								HTML('<hr noshade="noshade" />'),
								fluidRow(
									column(4,selectInput("cal_model", "Select calibration model (red in above plot)", choices = c("linear","quadratic"), "linear")),
									column(4,checkboxInput("cal_model_0intercept", "Force 0-intercept?",  width = NULL)),
									column(4,checkboxInput("cal_model_weight", "Weight by inverse target intensity?",  width = NULL))
								),
								textOutput('cal_model_summary'),
								HTML('<hr noshade="noshade" />'),
								fluidRow(
									column(4,
										checkboxInput("cal_model_bound_low", "Set a lower bound for the intensity ratio?",  width = NULL),
										numericInput("cal_model_bound_low_value", "Lower ratio bound:", 0)
									),
									column(4,
										checkboxInput("cal_model_bound_up", "Set an upper bound for the intensity ratio?",  width = NULL),
										numericInput("cal_model_bound_up_value", "Upper ratio bound:", 20)
									)									
								),						
								HTML('<hr noshade="noshade" />'),
								helpText("Click into the below table rows to select and deselect data points for the above calibration model (red line):"),
								fluidRow(
									column(12,DT::dataTableOutput('cal_table'))
								),
								HTML('<hr noshade="noshade" />'),		
								numericInput("use_precision", "Numeric match precision", 0.01),
								bsPopover("use_precision", 
									title = "Precision with which intensity ratios at a given concentration from a calibration model are matched to those in the above table.",
									content = "Only applied if a calibration model exists. Rows in the table for which no matches exist are automatically deselected.", 
									placement = "right", trigger = "hover")
							)
						)
					)
				)#,	
				#tabPanel("Import",	
				#	HTML('<hr noshade="noshade" />'),
				#	helpText("Import calibration files from another project (WARNING - replaces calibration files in current project):"),
				#	textInput("import_pro_dir_Cal", "", value = "C:\\...\\other_project_name"),
				#	bsPopover("import_pro_dir_Cal", 
				#		title = "Insert full path, including the project folder, but excluding the logfile.emp.",
				#		content = "Using your OS explorer, you may navigate into your project folder and copy/paste the full path.", 
				#		placement = "right", trigger = "hover"),
				#	bsButton("Import_project_Cal","Import",style="info")	
				#)
			)
        ),
        ########################################################################
        # RESULTS ##############################################################
        ########################################################################
        tabPanel("Results", 	
			tabsetPanel( 
				######################################################################################################################s
				tabPanel("Intensity control",
					tabsetPanel(
						tabPanel("Positive ionization ",
							HTML('<hr noshade="noshade" />'),
							bsCollapse(multiple = TRUE, open = NULL, id = "collapse_intens_control_pos",
								bsCollapsePanel(title = "Boxplot peak intensity distributions", 

									plotOutput("int_box_pos", 
										dblclick = "int_box_pos_dblclick",
										click = "int_box_pos_click",
										brush = brushOpts(
											id = "int_box_pos_brush",
											resetOnNew = TRUE,
											direction = c("x"),
											delay = 0
										),
										height = "500px"
									),								
									HTML('<hr noshade="noshade" />'),
									tags$p(align = "justify", "The above plot compares the log intensity distribution of picked peaks among all the included positive ionization mode files, using boxplots. 
									Files are ordered by their Date/Time.
									For larger number of files, boxplot statistics are simplified to dots (median) and two lines connecting the first and third quartiles with the lower and upper
									whiskers, respectively. For smaller numbers of files depicted (e.g., after zooming), conventional boxplots are provided, i.e., 
									the distance between the first and third quartiles is additionally outlined by a box and the intensity median by a bold line in this box.
									Whiskers are set as the most extreme (lower and upper) data point which departs no more than 1.5 times the quartile range from the median. 
									Peak intensities departing further than that are given as gray dots."),
									tags$p(align = "justify", "In addition, the number of peaks picked for each file is indicated as a red cross, and associated with the right axis.")
								),
								bsCollapsePanel(title = "Quantile peak intensity distributions", 
									plotOutput("int_quantiles_pos", 
										height = "330px"
									),					
									HTML('<h1 align="center"> &#x21e7; </h1> '),
									plotOutput("int_maxmed_pos", 
										dblclick = "int_maxmed_pos_dblclick",
										click = "int_maxmed_pos_click",
										brush = brushOpts(
											id = "int_maxmed_pos_brush",
											resetOnNew = TRUE,
											delay = 0
										),
										height = "450px"
									),	
									HTML('<hr noshade="noshade" />'),
									tags$p(align = "justify", "The above top panel shows the log intensity quantile distribution of picked peaks for each of the included positive ionization 
									mode files (gray lines). In addition, the red line shows the mean intensity at a quantile over all these files. 
									In turn, the bottom plot details the maximum (x-axis) and median (y-axis) quantile intensity 
									deviation of each file from the mentioned quantile mean, with files stated by their IDs and connected by the gray line for their temporal order. 
									Please note that top panel quantile distributions become highlighted in black for depicted files upon zooming in the bottom plot.")
								),
								bsCollapsePanel(title = "IS-based intensity normalization", 
									plotOutput("int_norm_ISTD_pos_median", 
										dblclick = "int_norm_ISTD_pos_median_dblclick",
										click = "int_norm_ISTD_pos_median_click",
										brush = brushOpts(
											id = "int_norm_ISTD_pos_median_brush",
											resetOnNew = TRUE,
											direction = c("x"),
											delay = 0
										),
										height = "500px"
									),
									plotOutput("int_norm_ISTD_pos_counts", 
										dblclick = "int_norm_ISTD_pos_counts_dblclick",
										click = "int_norm_ISTD_pos_counts_click",
										brush = brushOpts(
											id = "int_norm_ISTD_pos_counts_brush",
											resetOnNew = TRUE,
											direction = c("x")
										),
										height = "300px"
									),
									HTML('<hr noshade="noshade" />'),
									tags$p(align = "justify", "The top panel above shows the deviation of a peak intensity from the median log-intensity of its profile (gray dots), 
										for each internal standard (ISTD) and over all the positive ionization files used during profiling. The deviation is expressed as an intensity 
										ratio between peak and median intensity; positive deviations hence indicate ISTD peak intensities above their longer-termed profile median. \n\n
										Files are ordered by their Date/Time, with gray and orange lines highlighting sample and blank measurements, respectively. 
										In turn, the red dots denote the median deviation over all internal standards for each such file; this value is then used for the intensity 
										normalization of all peaks derived from a file. The red line in the lower panel highlights the number of such ISTD peaks found per file. Please 
										note that the intensity normalization is only conducted if the preset number of ISTD peaks was reached per file (red dashed line; cp. Settings -> 
										Normalization -> Minimum number of IS profile peaks)."),
									tags$p(align = "justify", "In comparison, the green and blue lines of this lower panel show the overall peak count for randomly sampled non-ISTD 
										profiles which either contain blank peaks (blue line) or are not found by any picked peaks in blank measurements at all (green line). Their
										median intensity deviations from the individual profile medians are shown as blue and green dots in the top panel, too.")
								)	
							)
						),
						tabPanel("Negative ionization ",
							HTML('<hr noshade="noshade" />'),
							bsCollapse(multiple = TRUE, open = NULL, id = "collapse_intens_control_neg",
								bsCollapsePanel(title = "Boxplot peak intensity distributions", 
									plotOutput("int_box_neg", 
										dblclick = "int_box_neg_dblclick",
										click = "int_box_neg_click",
										brush = brushOpts(
											id = "int_box_neg_brush",
											resetOnNew = TRUE,
											direction = c("x"),
											delay = 0
										),
										height = "500px"
									),								
									HTML('<hr noshade="noshade" />'),
									tags$p(align = "justify", "The above plot compares the log intensity distribution of picked peaks among all the included negative ionization mode files, using boxplots. 
									For larger number of files, boxplot statistics are simplified to dots (median) and two lines connecting the first and third quartiles with the lower and upper
									whiskers, respectively. For smaller numbers of files depicted (e.g., after zooming), conventional boxplots are provided, i.e., 
									the distance between the first and third quartiles is additionally outlined by a box and the intensity median by a bold line in this box.
									Whiskers are set as the most extreme (lower and upper) data point which departs no more than 1.5 times the quartile range from the median. 
									Peak intensities departing further than that are given as gray dots."),
									tags$p(align = "justify", "In addition, the number of peaks picked for each file is indicated as a red cross, and associated with the right axis.")	
								),
								bsCollapsePanel(title = "Quantile peak intensity distributions", 
									plotOutput("int_quantiles_neg", 
										height = "330px"
									),		
									HTML('<h1 align="center"> &#x21e7; </h1> '),
									plotOutput("int_maxmed_neg", 
										dblclick = "int_maxmed_neg_dblclick",
										click = "int_maxmed_neg_click",
										brush = brushOpts(
											id = "int_maxmed_neg_brush",
											resetOnNew = TRUE,
											delay = 0
										),
										height = "450px"
									),	
									HTML('<hr noshade="noshade" />'),
									tags$p(align = "justify", "The above top panel shows the log intensity quantile distribution of picked peaks for each of the included negative ionization 
									mode files (gray lines). In addition, the red line shows the mean intensity at a quantile over all these files. 
									In turn, the bottom plot details the maximum (x-axis) and median (y-axis) quantile intensity 
									deviation of each file from the mentioned quantile mean, with files stated by their IDs and connected by the gray line for their temporal order. 
									Please note that top panel quantile distributions become highlighted in black for depicted files upon zooming in the bottom plot.")
								),
								bsCollapsePanel(title = "IS-based intensity normalization", 
									plotOutput("int_norm_ISTD_neg_median", 
										dblclick = "int_norm_ISTD_neg_median_dblclick",
										click = "int_norm_ISTD_neg_median_click",
										brush = brushOpts(
											id = "int_norm_ISTD_neg_median_brush",
											resetOnNew = TRUE,
											direction = c("x"),
											delay = 0
										),
										height = "500px"
									),
									plotOutput("int_norm_ISTD_neg_counts", 
										dblclick = "int_norm_ISTD_neg_counts_dblclick",
										click = "int_norm_ISTD_neg_counts_click",
										brush = brushOpts(
											id = "int_norm_ISTD_neg_counts_brush",
											resetOnNew = TRUE,
											direction = c("x"),
											delay = 0
										),
										height = "300px"
									),
									HTML('<hr noshade="noshade" />'),
									tags$p(align = "justify", "The top panel above shows the deviation of a peak intensity from the median log-intensity of its profile (gray dots), 
										for each internal standard (ISTD) and over all the negative ionization files used during profiling. The deviation is expressed as an intensity 
										ratio between peak and median intensity; positive deviations hence indicate ISTD peak intensities above their longer-termed profile median. \n\n
										Files are ordered by their Date/Time, with gray and orange lines highlighting sample and blank measurements, respectively. 
										In turn, the red dots denote the median deviation over all internal standards for each such file; this value is then used for the intensity 
										normalization of all peaks derived from a file. The red line in the lower panel highlights the number of such ISTD peaks found per file. Please 
										note that the intensity normalization is only conducted if the preset number of ISTD peaks was reached per file (red dashed line; cp. Settings -> 
										Normalization -> Minimum number of IS profile peaks)."),
									tags$p(align = "justify", "In comparison, the green and blue lines of this lower panel show the overall peak count for randomly sampled non-ISTD 
										profiles which either contain blank peaks (blue line) or are not found by any picked peaks in blank measurements at all (green line). Their
										median intensity deviations from the individual profile medians are shown as blue and green dots in the top panel, too.")
								)	
							)						
						)
					)	
				),
				###################################################################################################################### 
                tabPanel("Data viewer",
					HTML('<hr noshade="noshade" />'),
					tags$p(align="justify","This tab plots the centroided data points and picked peaks for an individual file - please first specify the ID of this file. 
					The ID associated with a file can be found in the second column of the files table in the Files tab."),											
					div(style = widget_style3,numericInput("sel_meas_ID", "Type in file ID:", 0)),
					textOutput('file_viewer_name'),
					textOutput('file_viewer_type'),
					textOutput('file_viewer_mode'),
					HTML('<hr noshade="noshade" />'),	
					conditionalPanel(			
						condition = "(input.sel_meas_ID != 0) & (output.file_viewer_name != 'Invalid file ID')",							
						bsCollapse(multiple = FALSE, open = NULL, id = "collapse_screen_pos_one",
							bsCollapsePanel(title="All picked peaks & raw data", 
								navbarPage("Settings:",
									tabPanel("Data",
										fluidRow(
											column(3,checkboxInput("peaks_mz_RT_use_peaks", "Plot peaks?", TRUE)),
											column(3,checkboxInput("peaks_mz_RT_use_IDs", "Display peak IDs?", FALSE))
										),
										HTML('<hr noshade="noshade" />'),
										fluidRow(											
											column(4,checkboxInput("peaks_mz_RT_use_raw", "Show raw data (as density grid >1E5 data points)?", FALSE)),
											column(3,checkboxInput("EIC_mz_RT_use_IDs", "Display EIC IDs?", FALSE))								
										),
										HTML('<hr noshade="noshade" />')
									),
									tabPanel("View",
										fluidRow(
											column(5,checkboxInput('showPanel1', 'Show marginal intensity distributions (for <1E5 raw data points)', FALSE)),
											column(5,checkboxInput('showPanel2', 'Show interactive 3D plot (for <1E5 raw data points)', FALSE))
										),
										HTML('<hr noshade="noshade" />')
									),
									tabPanel("Search & Filter",
										fluidRow(
											column(4,checkboxInput("peaks_mz_RT_use_window", "Add a blue search window centered at coordinates:", FALSE)),									
											column(2,	
												textInput("peaks_mz_RT_use_window_mass", "m/z",  "216.101"),
												helpText( a("Calculate a mass?", href="http://www.envipat.eawag.ch/index.php",target="_blank"))
											),
											column(2,textInput("peaks_mz_RT_use_window_RT", "RT [s]", "500")),
											column(3,textInput("peaks_mz_RT_use_window_RT_tol", "RT tolerance [s]", "60"))
										),
										fluidRow(
											column(4,checkboxInput("peaks_mz_RT_use_bar", "Add a blue ppm bar to search window?", FALSE)),									
											column(3,textInput("peaks_mz_RT_use_bar_value", "ppm", "10"))
										),
										HTML('<hr noshade="noshade" />'),
										fluidRow(
											column(6,sliderInput("plot_filter_intensity", "Filter intensity range (log10):", min=0, max=8, value=c(0,8), width = '100%', dragRange=TRUE))
										),
										HTML('<hr noshade="noshade" />'),
										helpText("If included in the workflow:"),							
										fluidRow(
											column(4,checkboxInput("plot_filter_blind", "Include blind subtraction?", FALSE)),									
											column(4,checkboxInput("plot_filter_replicates", "Include replicate filter?", FALSE))
										),
										HTML('<hr noshade="noshade" />')
									),
									tabPanel("Hide",HTML('<hr noshade="noshade" />')),
									inverse=FALSE
								),										
								plotOutput("plot_peaks_mz_RT", 
									width = "100%", height = "650px",
									dblclick = "plot_peaks_mz_RT_dblclick",
									#click = "plot_peaks_mz_RT_click",
									hover = "plot_peaks_mz_RT_hover",
									brush = brushOpts(
										id = "plot_peaks_mz_RT_brush",
										resetOnNew = TRUE
									)
								),
								HTML('<p><font>
									<span style="color:black"><b>&bull;</b></span> Picked peaks <span style="color:gray"><b>&bull;</b></span> Raw data centroids <span style="color:red"><b>&bull;</b></span> Peak centroids 							
								</p></font>'),
								conditionalPanel(condition = 'input.showPanel1',
									plotOutput("plot_peaks_mz_int", width = "100%", height = "320px"),
									plotOutput("plot_peaks_RT_int", width = "100%", height = "320px")
								),
								conditionalPanel(condition = 'input.showPanel2',
									plotly::plotlyOutput("plot_peaks_3D", width = "100%", height = "800px")
								)
							),
							bsCollapsePanel(title="Individual EICs and peaks", 
								div(style = widget_style3,numericInput("sel_peak_ID", "Specify peak ID:", 0)),
								HTML('<hr noshade="noshade" />'),								
								plotOutput("EIC1", 
									dblclick = "plot_EIC1_dblclick",
									brush = brushOpts(
										id = "plot_EIC1_brush",
										resetOnNew = TRUE
										)
								),
								plotOutput("EIC2", 
									dblclick = "plot_EIC2_dblclick",
									brush = brushOpts(
										id = "plot_EIC2_brush",
										resetOnNew = TRUE
										)
								),
								plotOutput("EIC3", 
									dblclick = "plot_EIC3_dblclick",
									brush = brushOpts(
										id = "plot_EIC3_brush",
										resetOnNew = TRUE
										)
								)						
							)
						)
					)
                ),
				######################################################################################################################
				tabPanel("Processing",      
					HTML('<hr noshade="noshade" />'),				
					tags$p(align="justify","Preprocessing and summary results of this tab are shown for each file individually; please first specify the ID of the file of interest. 
						The ID associated with a file can be found in the second column of the files table in the Files tab."),											
					numericInput("sel_meas", "Type in file ID:", 0),
					conditionalPanel(			
						condition = "(output.dowhat != 'Invalid ID chosen to view processing results.') && (input.sel_meas != 0)",
						textOutput('file_proc_name'),
						textOutput('file_proc_type'),
						textOutput('file_proc_mode'),
						HTML('<hr noshade="noshade" />'),
						tags$p(align="justify","Minimize the below Peak list & chromatograms section for faster result loading."),
						bsCollapse(multiple = TRUE, open = NULL, id = "collapse_processing",
							bsCollapsePanel(title="Peak inventory & intensity distributions",
								fluidRow(										
									column(3,tags$h5("Total number of peaks: "),textOutput('file_peak_number') ),
									column(3,tags$h5("% of peaks affected by blind annotation: "),textOutput('file_blind_aff')),
									column(3,tags$h5("% of peaks affected by blind annotation, below threshold: "),textOutput('file_blind_rem')),										
									column(3,tags$h5("% of peaks removed by replicate filter: "),textOutput('file_repl_rem'))										
								),	
								HTML('<hr noshade="noshade" />'),
								textOutput('showintensitydistrib'),
								conditionalPanel(			
									condition = "output.showintensitydistrib == 'Centroid & peak intensity distribution:'",	
									imageOutput("peakhist_pic", height="auto"),
									tags$p(align="justify","The white histogram bars in the above plot depict the absolute intensity distribution over all centroid data points available for the
											selected file. Normally, centroid data points are more frequent at lower intensities than at higher ones in this histogram.
											In contrast, the red bars show the same distribution for the subset of centroid data points incorporated into picked peaks only.
											Moreover, blue dots and their associated right axis indicate the fraction of centroids incorporated into peaks for each such histogram class.
											Usually, this fraction somewhat increases with intensity. However, a deviation from this increase occurs if highly intense baselines are present, 
											as these baselines can contain intense centroid peaks not commonly incorporated into picked peaks.")
								),
								HTML('<hr noshade="noshade" />'),
								HTML('<p><a href="http://www.looscomputing.ch/eng/enviMass/topics/peakpicking.htm" style="color:rgb(60, 100, 60); text-decoration: none"; target="_blank"><p align="left">&#8594; More information on EIC extraction and peak picking.</a></p>')
							),
							bsCollapsePanel(title="Peak list & chromatograms",
								tags$p(align="justify","The below sortable table lists all picked peaks and their characteristics for the above selected file. Click into table row(s) to inspect the chromatograms of selected peak(s)."),
								HTML('<hr noshade="noshade" />'),
						        conditionalPanel(     
						            condition = "(typeof input.exp_peaklist_rows_selected !== 'undefined') && (input.exp_peaklist_rows_selected.length > 0)",  
						            #####################################################      
						            tags$h5("Chromatograms of selected peaks:"),    				         
						            plotOutput("peak_chromat",
						                click = "peak_chromat_click",
						                dblclick = "peak_chromat_dblclick",
						                brush = brushOpts(
						                    id = "peak_chromat_brush",
						                    direction = c("xy"),
						                    resetOnNew = TRUE,
						                    delay = 0
						                ),
						                height = "360px"
						            ),
						            fluidRow(
						            	column(width = 3, radioButtons("peak_chromat_norm", "Normalize intensities?", c("No"="FALSE","Yes"="TRUE"), inline = TRUE)), 
						            	column(width = 3, radioButtons("peak_chromat_time", "Show RT in", c("Seconds"="seconds","Minutes"="minutes"), inline = FALSE)),
						            	column(width = 3, radioButtons("peak_chromat_type", "Show chromatograms for:", c("Full EICs with peaks"="TRUE","Peaks only"="FALSE"), inline = FALSE)),						            	
						            	column(width = 3.5, bsButton("peak_chromat_refresh"," Refresh and deselect peak table?", style="default",icon=NULL))
						            ),								         
						            #####################################################
						            HTML('<hr noshade="noshade" />')
						        ), 
						        HTML('<p><a href="http://www.looscomputing.ch/eng/enviMass/topics/peaklist_export.htm" style="color:rgb(60, 100, 60); text-decoration: none"; target="_blank"><p align="left">&#8594; Help on table navigation, filtering & table contents.</a></p>'),
								DT::dataTableOutput('exp_peaklist')
							),
							bsCollapsePanel(title="Blank/blind peak detection",
								fluidRow(										
									column(4,tags$h5("% of peaks affected by blind annotation: "),textOutput('file_blind_aff2')),
									column(7,tags$h5("% of peaks affected by blind annotation, below factor set in the Blind Settings: "),textOutput('file_blind_rem2'))										
								),
								HTML('<hr noshade="noshade" />'),
								textOutput('showblank'),
								conditionalPanel(			
									condition = "output.showblank == 'Blank/blind peak tagging (subtraction) results:'",	
									plotOutput("blind_boxplot", height = "200px"),
									radioButtons("blind_boxplot_log", "Plot intensity ratio on log scale?", c("no"="FALSE","yes"="TRUE"),inline=TRUE),
									tags$p(align="justify","The above boxplot shows the distribution of intensity ratios between picked peaks in the selected sample file and their counterpart peaks of maximum 
										intensity in the associated blank or blind files. 
										Sample peaks which could not be matched to any blank peaks do not contribute to this distribution; i.e., they are not part of the listed
										% of peaks affected by blind annotation.")							
								),
								HTML('<hr noshade="noshade" />'),
								HTML('<p><a href="http://www.looscomputing.ch/eng/enviMass/topics/blind.htm" style="color:rgb(60, 100, 60); text-decoration: none"; target="_blank"><p align="left">&#8594; More information on blind/blank peak subtraction and annotation.</a></p>')
							),
							bsCollapsePanel(title="Mass recalibration",
								textOutput('showrecal'),
								conditionalPanel(			
									condition = "output.showrecal != 'No mass recalibration results available.'",	
									imageOutput("recal_pic", height="auto"),
									tags$p(align="justify","The above plot shows the difference between calculated and thus expected masses of known compounds and their measured/observed peaks (ordinate) as a function 
										of their mass (abscissa). The considered masses are restricted to the most intense calculated centroid per compound, prior to any detailed screening in the workflow.
										Moreover, the red line is a smoothing spline model used to recalibrate this mass difference (i.e., the smoothed mass difference of the recalibrated masses would lie on the
										green zero-intercept). Any red dashed lines, if visible, denote the maximum allowable mass recalibration as set in the Mass recalibration section of the Settings tab.
										Grey dashed lines denote mass differences of 2, 5, and 10 ppm each.")
								),
								HTML('<hr noshade="noshade" />'),
								HTML('<p><a href="http://www.looscomputing.ch/eng/enviMass/topics/recalibration.htm" style="color:rgb(60, 100, 60); text-decoration: none"; target="_blank"><p align="left">&#8594; More information on mass recalibration.</a></p>')
							),
							bsCollapsePanel(title="LOD interpolation",						
								textOutput('showLOD'),
								conditionalPanel(			
									condition = "output.showLOD == 'LOD interpolation results:'",									
									imageOutput("LOD_pic", height = "280px"),
									tags$p(align="justify","The above plot outlines the RT-dependent LOD estimation used for workflow steps such as screening or the calculation of upper bounds on atom counts.
										Therein, each grey dot corresponds to a picked peak, whereas gray lines seperate the RT range into bins. 
										The lower 10th percentile peaks of each such bins are then used to smooth a conservative intensity threshold (red line), 
										below which peaks are not expected to be either present or picked.")
								),
								HTML('<hr noshade="noshade" />'),
								HTML('<p><a href="http://www.looscomputing.ch/eng/enviMass/topics/lod.htm" style="color:rgb(60, 100, 60); text-decoration: none"; target="_blank"><p align="left">&#8594; More information on LOD interpolation.</a></p>')
							)
						)		
					)
				),
				######################################################################################################################
				tabPanel("Screening & quantification",
					tabsetPanel(
						tabPanel("Positive ionization",
							fluidRow(
								column(3, selectInput("Pos_compound_select",label="",choices=c("Choose","Target compounds","Internal standards","Quantification","Recovery","File-wise counts"), 
									selected = "Choose", multiple = FALSE)),
								conditionalPanel(			
									condition = "input.Pos_compound_select == 'Internal standards' || input.Pos_compound_select == 'Target compounds'",										
										column(3, selectInput("screen_pos_summarize", label="", choices = c("Show all adducts"="yes","Collapse adducts"="no"), "yes")),
										column(4, selectInput("Pos_type_select",label="",choices=c("Sample/blind files","Calibration files"), 
											selected = "Non-calibration files", multiple = FALSE))			
								)							
							),
							conditionalPanel(
								condition = "input.Pos_compound_select == 'Internal standards' || input.Pos_compound_select == 'Target compounds'",	
									conditionalPanel(			
									condition = "input.screen_pos_summarize == 'yes'",	
									bsCollapse(multiple = FALSE, open = NULL, id = "collapse_screen_pos_one",
										bsCollapsePanel(title="Pattern match for selected compound", #style="info",
											textOutput('screening_details_comp_pos'),
											plotOutput("plot_pattern_pos", 
												dblclick = "plot_pattern_pos_dblclick",
												brush = brushOpts(
													id = "plot_pattern_pos_brush",
													resetOnNew = TRUE
												)
											),
											HTML('<hr noshade="noshade" />'),
											tags$p(align="justify","The above red bars show the theoretical centroid pattern for the selected compound. Screening matches of this pattern with measured peaks 
												are depicted as green circles. The grey lines connect peaks which have been jointly matched to this pattern, i.e., are present in the same peaklist of a file.
												The exact values are tabulated below."),
											HTML('<hr noshade="noshade" />'),
											fluidRow(
												column(4, tableOutput("table_pattern_pos")),
												column(5, textOutput("table_patternRT_pos"))
											)
										),
										bsCollapsePanel(title="Characteristics for selected compound",
											textOutput('screening_details_comp_pos2'),
											HTML('<hr noshade="noshade" />'),
											fluidRow(										
												column(4, selectInput("selec_pos_x",label="x axis",
													choices=c("m/z","RT","Intensity","Date&time","Type","Place"),selected = "m/z", multiple = FALSE)),
												column(4, selectInput("selec_pos_y",label="y axis",
													choices=c("m/z","RT","Intensity","Date&time","Type","Place"),selected = "Intensity", multiple = FALSE)),									
												column(4, radioButtons("selec_pos_log_rat", "Log intensity?", c("yes"="yes","no"="no"),inline=TRUE))
											),
											HTML('<hr noshade="noshade" />'),						
											plotOutput("plot_selec_dist_pos"),

											tags$p(align="justify","The above scatterplot shows peak characteristics for the selected compound, as set by the two drop-down selection for each axis.
												Peaks which have been matched to the theoretical pattern of this compound with a score below the Cutoff-score defined in the Settings -> Screening are shown
												in gray; those above as black dots."),
											HTML('<hr noshade="noshade" />'),
											textOutput('info_IS_bounds_pos'),
											conditionalPanel(				
												condition = "(input.Pos_compound_select == 'Internal standards') && (output.info_IS_bounds_pos != 'Compound/adduct not used for quantification')",					
													fluidRow(										
														column(3,numericInput("screen_int_pos_low", "Lower bound", 0,step=0.1)),
														column(3,numericInput("screen_int_pos_up", "Upper bound", 10,step=0.1)),
														column(2,bsButton("save_int_pos"," Adopt",style="warning",icon=icon("bookmark")))
													)	
											)
										),
										bsCollapsePanel(title="Screening table for selected compound", 
											textOutput('screening_details_comp_pos3'),
											conditionalPanel(     
						            			condition = "(typeof input.Table_screening_selected_pos_rows_selected !== 'undefined') && (input.Table_screening_selected_pos_rows_selected.length > 0)",  
												HTML('<hr noshade="noshade" />'),
												plotOutput("screening_chromat_pos",
								                    click = "screening_chromat_pos_click",
								                    dblclick = "screening_chromat_pos_dblclick",
								                    #hover = "homol_counts_hover",
								                    brush = brushOpts(
								                        id = "screening_chromat_pos_brush",
								                        direction = c("xy"),
								                        resetOnNew = TRUE,
								                        delay = 0
								                    ),       
													height = "330px"
												),
									            fluidRow(
									            	column(width = 3, radioButtons("screening_chromat_pos_norm", "Normalize intensities?", c("No"="FALSE","Yes"="TRUE"), inline = TRUE)), 
									            	column(width = 3, radioButtons("screening_chromat_pos_time", "Show RT in", c("Seconds"="seconds","Minutes"="minutes"), inline = FALSE)),
									            	column(width = 3, radioButtons("screening_chromat_pos_type", "Show chromatograms for:", c("Full EICs with peaks"="TRUE","Peaks only"="FALSE"), inline = FALSE))
									            )	
											),
											HTML('<hr noshade="noshade" />'),
											tags$p(align="justify","Click on a table row to show chromatograms."),
											HTML('<hr noshade="noshade" />'),
											DT::dataTableOutput('Table_screening_selected_pos')
										)
									)		
								),
								HTML('<hr noshade="noshade" />'),
								tags$p(align="justify","The below Sample and Blank matches show the number of files with screening matches above the cutoff score, 
								with multiple matches per file above this cutoff merged."),
###								
								DT::dataTableOutput('Table_screening_pos'),
								
								
								
								
###															
								HTML('<hr noshade="noshade" />'),
								bsCollapse(multiple = FALSE, open = NULL, id = "collapse_screen_pos_all",
									bsCollapsePanel(title="Summary plots",
										fluidRow(
											column(width = 6, offset = 0.6,
													tags$p(align="justify","Summary characteristics over all peaks matched to the screened compounds of the above table. Information based on peaks which have been 
														matched to the theoretical pattern of a compound with a score equal or above the Cutoff-score defined in the Settings -> Screening is shown
														in black or white; and otherwise in gray.")
											),
											column(3, selectInput("Summ_pos_x",label="x axis",choices=c("m/z","Measured RT","log Intensity","m/z deviation [ppm]","RT deviation within","Time sequence","Expected RT"),selected = "m/z", multiple = FALSE)),
											column(3, selectInput("Summ_pos_y",label="y axis",choices=c("m/z","Measured RT","log Intensity","m/z deviation [ppm]","RT deviation within","Time sequence","Expected RT"),selected = "RT", multiple = FALSE))									
										),
										plotOutput("plot_pattern_distrib_pos"),
										HTML('<hr noshade="noshade" />'),
										fluidRow(								
											column(7,textOutput('count_aboveBlank_pos')),
											column(3,radioButtons("screen_pos_log_rat", "Log scale?", c("yes"="yes","no"="no"),inline=TRUE))

										),
										plotOutput("plot_aboveBlank_pos",height = 250)
									)
								)	
							),		
							conditionalPanel(			
								condition = "input.Pos_compound_select == 'Quantification'",							
								HTML('<hr noshade="noshade" />'),
								DT::dataTableOutput('target_quant_table_pos')
							),		
							conditionalPanel(			
								condition = "input.Pos_compound_select == 'Recovery'",							
								HTML('<hr noshade="noshade" />'),
								DT::dataTableOutput('target_recov_table_pos')
							),												
							conditionalPanel(			
								condition = "input.Pos_compound_select == 'File-wise counts'",	
								tags$p(align="justify","The below table lists the number of compounds which have been positively screened above the cutoff score per file. Matches for different adducts of the same compound are counted separately."),
								HTML('<hr noshade="noshade" />'),
								DT::dataTableOutput('count_file_compound_pos'),
								HTML('<hr noshade="noshade" />')
							)						
						),
						tabPanel("Negative ionization",
							fluidRow(
								column(3, selectInput("Neg_compound_select",label="",choices=c("Choose","Target compounds","Internal standards","Quantification","Recovery","File-wise counts"), 
									selected = "Choose", multiple = FALSE)),
								conditionalPanel(			
									condition = "input.Neg_compound_select == 'Internal standards' || input.Neg_compound_select == 'Target compounds'",										
										column(3, selectInput("screen_neg_summarize", label="", choices = c("Show all adducts"="yes","Collapse adducts"="no"), "yes")),
										column(4, selectInput("Neg_type_select",label="",choices=c("Sample/blind files","Calibration files"), 
											selected = "Non-calibration files", multiple = FALSE))			
								)							
							),
							conditionalPanel(			
								condition = "input.Neg_compound_select == 'Internal standards' || input.Neg_compound_select == 'Target compounds'",	
								conditionalPanel(			
								condition = "input.screen_neg_summarize == 'yes'",	
									bsCollapse(multiple = TRUE, open = NULL, id = "collapse_screen_neg_one",
										bsCollapsePanel(title="Pattern match for selected compound", #style="info",
											textOutput('screening_details_comp_neg'),
											plotOutput("plot_pattern_neg", 
												dblclick = "plot_pattern_neg_dblclick",
												brush = brushOpts(
													id = "plot_pattern_neg_brush",
													resetOnNew = TRUE
												)
											),
											HTML('<hr noshade="noshade" />'),
											tags$p(align="justify","The above red bars show the theoretical centroid pattern for the selected compound. Screening matches of this pattern with measured peaks 
												are depicted as green circles. The grey lines connect peaks which have been jointly matched to this pattern, i.e., are present in the same peaklist of a file.
												The exact values are tabulated below."),
											HTML('<hr noshade="noshade" />'),
											fluidRow(
												column(3, tableOutput("table_pattern_neg")),
												column(3, textOutput("table_patternRT_neg"))
											)
										),
										bsCollapsePanel(title="Characteristics for selected compound",
											textOutput('screening_details_comp_neg2'),
											HTML('<hr noshade="noshade" />'),
											fluidRow(										
												column(4, selectInput("selec_neg_x",label="x axis",
													choices=c("m/z","RT","Intensity","Date&time","Type","Place"),selected = "m/z", multiple = FALSE)),
												column(4, selectInput("selec_neg_y",label="y axis",
													choices=c("m/z","RT","Intensity","Date&time","Type","Place"),selected = "Intensity", multiple = FALSE)),									
												column(4, radioButtons("selec_neg_log_rat", "Log intensity?", c("yes"="yes","no"="no"),inline=TRUE))
											),
											HTML('<hr noshade="noshade" />'),						
											plotOutput("plot_selec_dist_neg"),
											tags$p(align="justify","The above scatterplot shows peak characteristics for the selected compound, as set by the two drop-down selection for each axis.
												Peaks which have been matched to the theoretical pattern of this compound with a score below the Cutoff-score defined in the Settings -> Screening are shown
												in gray; those above as black dots."),
											textOutput('info_IS_bounds_neg'),
											HTML('<hr noshade="noshade" />'),
											conditionalPanel(				
												condition = "(input.Neg_compound_select == 'Internal standards') & (output.info_IS_bounds_neg != 'Compound/adduct not used for quantification')",					
													fluidRow(										
														column(3,numericInput("screen_int_neg_low", "Lower bound", 0,step=0.1)),
														column(3,numericInput("screen_int_neg_up", "Upper bound", 10,step=0.1)),
														column(2,bsButton("save_int_neg"," Adopt",style="warning",icon=icon("bookmark")))
													)	
											)
										),
										bsCollapsePanel(title="Screening table for selected compound", 
											textOutput('screening_details_comp_neg3'),
											conditionalPanel(     
						            			condition = "(typeof input.Table_screening_selected_neg_rows_selected !== 'undefined') && (input.Table_screening_selected_neg_rows_selected.length > 0)",  
												HTML('<hr noshade="noshade" />'),
												plotOutput("screening_chromat_neg",
								                    click = "screening_chromat_neg_click",
								                    dblclick = "screening_chromat_neg_dblclick",
								                    #hover = "homol_counts_hover",
								                    brush = brushOpts(
								                        id = "screening_chromat_neg_brush",
								                        direction = c("xy"),
								                        resetOnNew = TRUE,
								                        delay = 0
								                    ),       
													height = "330px"
												),
									            fluidRow(
									            	column(width = 3, radioButtons("screening_chromat_neg_norm", "Normalize intensities?", c("No"="FALSE","Yes"="TRUE"), inline = TRUE)), 
									            	column(width = 3, radioButtons("screening_chromat_neg_time", "Show RT in", c("Seconds"="seconds","Minutes"="minutes"), inline = FALSE)),
									            	column(width = 3, radioButtons("screening_chromat_neg_type", "Show chromatograms for:", c("Full EICs with peaks"="TRUE","Peaks only"="FALSE"), inline = FALSE))
									            )													
											),
											HTML('<hr noshade="noshade" />'),
											tags$p(align="justify","Click on a table row to show chromatograms."),										
											HTML('<hr noshade="noshade" />'),
											DT::dataTableOutput('Table_screening_selected_neg')
										)
									)		
								),
								HTML('<hr noshade="noshade" />'),
								tags$p(align="justify","The below sample and blank matches give the number of files with matches above the cutoff score, 
								with multiple matches per file above this cutoff merged."),
								DT::dataTableOutput('Table_screening_neg'),
								HTML('<hr noshade="noshade" />'),
								bsCollapse(multiple = FALSE, open = NULL, id = "collapse_screen_neg_all",
									bsCollapsePanel(title="Summary plots",
										fluidRow(
											column(width = 6, offset = 0.6,
													tags$p(align="justify","Summary characteristics over all peaks matched to the screened compounds of the above table. Information based on peaks which have been 
														matched to the theoretical pattern of a compound with a score equal or above the Cutoff-score defined in the Settings -> Screening is shown
														in black or white; and otherwise in gray.")
											),
											column(4, selectInput("Summ_neg_x",label="x axis",choices=c("m/z","Measured RT","log Intensity","m/z deviation [ppm]","RT deviation within","Time sequence","Expected RT"),selected = "m/z", multiple = FALSE)),
											column(4, selectInput("Summ_neg_y",label="y axis",choices=c("m/z","Measured RT","log Intensity","m/z deviation [ppm]","RT deviation within","Time sequence","Expected RT"),selected = "RT", multiple = FALSE))									
										),
										plotOutput("plot_pattern_distrib_neg"),
										HTML('<hr noshade="noshade" />'),
										fluidRow(											
											column(7,textOutput('count_aboveBlank_neg')),
											column(3,radioButtons("screen_neg_log_rat", "Log scale?", c("yes"="yes","no"="no"),inline=TRUE))
										),
										plotOutput("plot_aboveBlank_neg",height = 250)
									)
								)	
							),		
							conditionalPanel(			
								condition = "input.Neg_compound_select == 'Quantification'",							
								HTML('<hr noshade="noshade" />'),
								DT::dataTableOutput('target_quant_table_neg')
							),								
							conditionalPanel(			
								condition = "input.Neg_compound_select == 'Recovery'",							
								HTML('<hr noshade="noshade" />'),
								DT::dataTableOutput('target_recov_table_neg')
							),	
							conditionalPanel(			
								condition = "input.Neg_compound_select == 'File-wise counts'",	
								tags$p(align="justify","The below table lists the number of compounds which have been positively screened above the cutoff score per file. Matches for different adducts of the same compound are counted separately."),
								HTML('<hr noshade="noshade" />'),
								DT::dataTableOutput('count_file_compound_neg'),
								HTML('<hr noshade="noshade" />')
							)						
						)

					)	
				),
				######################################################################################################################
				tabPanel("Grouping",	
					HTML('<hr noshade="noshade" />'),
					tags$p(align="justify","To show nontarget grouping and homologue series detection results in the below tabs, please first specify the ID of the file of interest. 
						The ID associated with any file can be found in the second column of the files table in the Files tab."),											
					div(style = widget_style3,numericInput("sel_meas_comp", "Type in file ID:", 0, min=0)),
					textOutput('comp_file_name'),
					textOutput('comp_file_type'),
					textOutput('comp_file_mode'),
					HTML('<hr noshade="noshade" />'),
					textOutput('sel_meas_comp_state'),
					textOutput('sel_meas_comp_state1'),
					textOutput('sel_meas_comp_state2'),
					HTML('<hr noshade="noshade" />'),
					conditionalPanel(
						condition = "(output.sel_meas_comp_state != 'Invalid file ID') && (output.sel_meas_comp_state != '')",
						tabsetPanel( 					
							tabPanel("Components",	
							conditionalPanel(			
								condition = "(output.sel_meas_comp_state1 != 'no nontarget components available ') && (output.sel_meas_comp_state != 'Invalid file ID')",
									HTML('<hr noshade="noshade" />'),
									tags$p(align="justify","The below collapsible subsections give summarized and detailed information on the chemical component peaks (i.e., different adducts and isotopologues of the same analyte) 
									detected in the selected file. Any row selections made in the Full component list table subsection will update the indivdual component viewed in the Single component subsection. You can hide
									most component details in the latter subsection for easier navigation."),
									HTML('<hr noshade="noshade" />'),
									bsCollapse(multiple = TRUE, open = "col3", id = "collapse3",	
										bsCollapsePanel("Summary", 			
											textOutput('num_peaks_all'),
											textOutput('num_peaks_remain_replicate'),
											textOutput('num_peaks_remain_blind'),
											textOutput('num_peaks_remain'),
											HTML('<hr noshade="noshade" />'),
											textOutput('num_comp'),
											textOutput('reduc'),
											textOutput('min2_size_comp'),
											textOutput('median_size_comp'),
											textOutput('max_size_comp'),
											HTML('<hr noshade="noshade" />'),
											textOutput('num_comp_tar'),
											textOutput('num_comp_ISTD'),
											textOutput('num_comp_blind_any'),
											textOutput('num_comp_blind_all'),
											HTML('<hr noshade="noshade" />'),
											textOutput('num_comp_nontarget')										
										),
										bsCollapsePanel("Single component", 
											tags$p(align="justify","You can select components either via the ID of a single peak contained in the component or via the ID of the component itself. In the first case, the selection of the peak ID
												will automatically update the ID of the associated component."),											
											fluidRow(										
												column(4,numericInput("sel_meas_comp_peak", "Type in peak ID:", 0, min=0)),
												column(4,numericInput("sel_meas_comp_comp", "Type in component ID:", 0, min=0)),
												column(4,textOutput('found_compo'))
											),
											HTML('<hr noshade="noshade" />'),
											plotOutput("comp_plot_spec",
							                    click = "comp_plot_spec_click",
							                    dblclick = "comp_plot_spec_dblclick",
							                    #hover = "homol_counts_hover",
							                    brush = brushOpts(
							                        id = "comp_plot_spec_brush",
							                        direction = c("x"),
							                        resetOnNew = TRUE,
							                        delay = 0
							                    ),       
												height = "300px"
											),
											plotOutput("comp_plot_chromat",
							                    click = "comp_plot_chromat_click",
							                    dblclick = "comp_plot_chromat_dblclick",
							                    #hover = "homol_counts_hover",
							                    brush = brushOpts(
							                        id = "comp_plot_chromat_brush",
							                        direction = c("xy"),
							                        resetOnNew = TRUE,
							                        delay = 0
							                    ),       
												height = "330px"
											),
											fluidRow(
									           	column(width = 3, radioButtons("comp_chromat_norm", "Normalize intensities?", c("No"="FALSE","Yes"="TRUE"), inline = TRUE)), 
									            column(width = 3, radioButtons("comp_chromat_time", "Show RT in", c("Seconds"="seconds","Minutes"="minutes"), inline = FALSE)),
									            column(width = 3, radioButtons("comp_chromat_type", "Show chromatograms for:", c("Full EICs with peaks"="TRUE","Peaks only"="FALSE"), inline = FALSE))
									        ),
											conditionalPanel(			
												condition = "(output.found_compo != 'Invalid peak ID') && (input.sel_meas_comp_comp != 0)",
												tags$p(align="justify","The above mass spectrum shows all peaks grouped into the selected component (colored bars) and all other non-component peaks that fall into 
												the mass and RT range of the selected component (gray bars, if any). The chromatogram plot below that mass spectrum shows the individual chromatograms of these component 
												peaks; their selection if affected by any zooming in the mass spectrum."),											
												HTML('<hr noshade="noshade" />')
											),
											materialSwitch(inputId = "comp_hide", label = "Show more component details and tables below?", value = TRUE, status = "danger", 
												right = TRUE, inline = TRUE, width = NULL),
											HTML('<hr noshade="noshade" />'), 												
											conditionalPanel(	
												condition = "(output.found_compo != 'Invalid peak ID') && (input.sel_meas_comp_comp != 0) && (!input.comp_hide)",
												plotOutput("comp_plot_circ",height = "690px"),
												tags$p(align="justify","The above circular plot shows all links (relations) that exist between peaks in the selected component, with 
													individual peaks positioned at an outer circle by increasing m/z, starting at the gray arrow. The most intense peak is marked by a circle.
													Peaks are listed with their IDs and component-wise intensity ranking (in brackets); those also present in blank/blind files are marked with
													an asterisk. Adduct links are shown in red, links between different isotopologues of the same adduct in blue."),													
												HTML('<hr noshade="noshade" />'),
												textOutput('which_comp_tar'),
												HTML('<hr noshade="noshade" />'),
												textOutput('which_comp_ISTD'),												
												HTML('<hr noshade="noshade" />'),	
												conditionalPanel(
													condition = "(output.found_compo != 'single_peak') && (output.found_compo != 'Invalid peak ID') && (input.sel_meas_comp_comp != 0)",												
													HTML('<p style="background-color:darkblue"; align="center"> <font color="#FFFFFF"> Component linkages </font></p> '),
													tags$p(align="justify","The below table lists all pairwise relations between peaks being part of the selected component.  Each of the two peaks (i.e., 'first' & 'second' peak)
														forming a link are listed by their IDs in the first and second column of the table."), 
													DT::dataTableOutput('comp_table_a'),
													HTML('<hr noshade="noshade" />')
												),
												HTML('<p style="background-color:darkblue"; align="center"> <font color="#FFFFFF"> Component peaks </font></p> '),
												tags$p(align="justify","The below table lists all peaks which are part of the selected component."), 
												DT::dataTableOutput('comp_table_b'),
												HTML('<hr noshade="noshade" />'),
												HTML('<p style="background-color:darkblue"; align="center"> <font color="#FFFFFF"> Component homologues </font></p> '),
												tags$p(align="justify","If peaks of the selected component form part of a homologue series, the individual series are listed in the below table."), 
												DT::dataTableOutput('comp_table_c'),										
												HTML('<hr noshade="noshade" />'),
												HTML('<p style="background-color:darkgrey"; align="center"> <font color="#FFFFFF"> All peaks </font></p> '),
												tags$p(align="justify","The below table lists all peaks which range in the RT- and m/z-window of the selected component, i.e., also peaks which were not
													grouped into the selected component."),
												DT::dataTableOutput('comp_table_d')		
											)
										),
										bsCollapsePanel("Full component list", 
											tags$p(align="justify","Each row in the below table corresponds to a component and its different isotoplogue and adduct peaks. 
												Peaks which are also present in blank/blind files are marked with an asterisk in the below table.
												Click into a table row to have further details on the selected component highlighted in the above Single component subsection."), 
											HTML('<hr noshade="noshade" />'),
											DT::dataTableOutput('comp_table_full')					
										)
									)
								)
							),
							tabPanel("Homologues",	
							conditionalPanel(			
								condition = "(output.sel_meas_comp_state2 != ' no homologue series detection results available') && (output.sel_meas_comp_state != 'Invalid file ID')",
									HTML('<hr noshade="noshade" />'),
									tags$p(align="justify","The below two subsections visualize and tabulate peak series indicative of homologuous compounds. Any zooms and selections made in 
										any of the overview plots A to C will filter the listing of Peaks in series in the first subtab of the collapsible Series table subsection. 
										Any selection or deselection of rows made in this Peaks in series table will in turn list either all series of which the selected peak is part of or 
										simply all available series in the second Series subtab, respectively. 
										In addition, selecting rows in the latter Series table will show the chromatographic peaks of the concerned series."),
									HTML('<hr noshade="noshade" />'),
									bsCollapse(multiple = TRUE, open = "col4", id = "collapse4",
										bsCollapsePanel("Series plot overview", 	
											HTML('<hr noshade="noshade" />'),
											tags$p(align="justify","Selections by brushing or zooming made on the mass differences (plot A) or RT differences (plot B) between homologues adjacent in a series 
												will jointly filter the plotting of series segments between such homologues in plot C. In turn, zooming in plot C will filter plots A and C. 
												Please note that peaks incorporated in any series are plotted as bold darkgrey points, whereas all other peaks are set in lightgrey in plot C."),
											HTML('<hr noshade="noshade" />'),
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
											HTML('<h1 align="center"> &#x21D3; </h1> ')
										),
										bsCollapsePanel("Series table", 
											tabsetPanel( 	
												tabPanel("Peaks in series ->",
													HTML('<hr noshade="noshade" />'),
													tags$p(align="justify","The below table lists all series peaks for which series segments are visible in plot C of the above series plot overview. 
														Selection of a row in this table will filter the series in the subsequent Series subtab, i.e., all those series corresponding to the Series ID(s)
														shown in the selected row. In addition, the peak and all its corresponding series will be highlighted in the above plot C."),
													HTML('<hr noshade="noshade" />'),
													DT::dataTableOutput('homol_series_peaks')
												),
						                        tabPanel("Series",
						                        	HTML('<hr noshade="noshade" />'),
						                        	tags$p(align="justify","The below table lists either all available peak series or - if a row in the previous Peaks in series subtab is selected - all 
						                        		series in which the selected peak is part of. Click on a row in the below Series table to get a series highlighted in the above plot C and to see its 
						                        		chromatographic peaks."),
						                            HTML('<hr noshade="noshade" />'),                
						                            conditionalPanel(     
						                              condition = "(input.homol_series_table_rows_selected > 0) && (input.homol_series_table_rows_selected.length > 0)",  
						                                #####################################################            
						                                plotOutput("homol_chromat",
						                                    click = "homol_chromat_click",
						                                    dblclick = "homol_chromat_dblclick",
						                                    #hover = "homol_chromat_hover",
						                                    brush = brushOpts(
						                                      id = "homol_chromat_brush",
						                                      direction = c("xy"),
						                                      resetOnNew = TRUE,
						                                      delay = 0
						                                    ),
						                                    height = "250px"
						                                ),
						                                #####################################################
						                                HTML('<hr noshade="noshade" />')
						                            ), 
						                            DT::dataTableOutput('homol_series_table')
						                        )
											)	
										)
									)
								)
							),
							tabPanel("Atom bounds",	
							conditionalPanel(			
								condition = "output.sel_meas_comp_state != 'Invalid file ID'",
									HTML('<hr noshade="noshade" />'),
									tags$p(align="justify","Given a monoisotopic peak (via its ID) and any additional peaks of higher mass and similar RT (if any), this tool estimates upper bounds for the number of atoms of selected elements.
										This estimate can be based (a) either only on those additional peaks contained in the same isotopologue group of the monoisotopic peak (if any, and if the isotopologue grouping 
										is included in the workflow) or (b) all additional peaks with a RT similar to the monoisotopic one. 
										Parameters for estimation are taken from those set for the isotopologue grouping (e.g., the intensity and RT tolerance) and the compound screening (e.g., LOD estimation), except
										for the mass shift to be set below."), 
									HTML('<hr noshade="noshade" />'),
									fluidRow(										
										column(3,numericInput("atom_bound_peak", "Type in peak ID:", 0, width='200px',min=0)),
										column(5,selectInput("atom_bound_addpeaks", "Select additional peaks:", c("(b) all peaks with similar RT"),selected="(b) all peaks with similar RT",multiple=FALSE))
									),	
									plotOutput("atom_bound_plot_peak",height = "400px"),
									HTML('<hr noshade="noshade" />'),
									HTML('<h1 align="left"> &#x21e9; </h1> '),
									selectInput("atom_bounds_this", label="Choose elements to include / exclude", choices=c("C","H","N","O","Cl","Br","S"), selected = c("C","H","N","O","Cl","Br"), multiple = TRUE),
									tags$p(align="justify","The below numeric values set the maximum allowable shift in between centroided peaks and their exact underlying isotopologue masses for each selected element, in ppm.
									This maximum shift depends on the resolution of your instrument and can sometimes be quite large (even for high-res MS) - better use conservative values to be safe."), 
									uiOutput("atom_bounds_that"),
									HTML('<hr noshade="noshade" />'),
									HTML('<h1 align="left"> &#x21e9; </h1> '),
									tags$p(align="justify","Would you like the bounds of atom counts for elements other than carbon to be restricted by their maximum atom count ratio to that of carbon? Only functional if carbon is selected above."), 
									radioButtons("atom_bound_wcarbon", "Include?", c("yes"="TRUE","no"="FALSE"),inline=TRUE),
									HTML('<hr noshade="noshade" />'),
									bsButton("atom_bounds_calculate","Refresh",style="success"),
									HTML('<h1 align="left"> &#x21e9; </h1> '),
									tags$p(align="justify","Maximum possible number of atoms for each element at charges of z = 1 to 4 based on the above selections:"), 
									DT::dataTableOutput('atom_count_table')
								)
							)	
						)	
					)
				),				
				######################################################################################################################
				tabPanel("Profiles",
					div(style = widget_style5,
						textOutput("had_ion"),	
						selectInput("Ion_mode", label=NULL, c("positive","negative"), selected = ("positive"), multiple = FALSE)
					),
					conditionalPanel(			
					condition = "output.had_ion != 'No profiles available for this ionization mode'",
						HTML('<hr noshade="noshade" />'),  
						#navbarPage("", 
						tabsetPanel( 					
							tabPanel("Overview & Filtering",	
									HTML('<hr noshade="noshade" />'),							
									tags$p(align="justify","The below filtering options determine the set of profiles accessible in the rightmost Single Profiles tab."), 
									HTML('<hr noshade="noshade" />'),
									tags$h4("Profile inventory:"), 
									div(style = widget_style7,tags$h5("Total number of peaks:"), textOutput("atprof1")),
									div(style = widget_style8,tags$h5("Number of profiles:"), textOutput("atprof2")),
									div(style = widget_style7,tags$h5("...containing blind peaks:"), textOutput("atprof3")),
									div(style = widget_style7,tags$h5("...with a past trend:"), textOutput("atprof4")),
									div(style = widget_style9,tags$h5("...with a current trend:"), textOutput("atprof5")),
									HTML('<hr noshade="noshade" />'),						
									tags$h4("Filter profile list:"),
									HTML('<hr noshade="noshade" />'),
									div(style = widget_style3,numericInput("filterProf_minmass", "Minimum m/z:", 0, min=0)),
									div(style = widget_style3,numericInput("filterProf_maxmass", "Maximum m/z:", 3000, min=0)),
									div(style = widget_style3,numericInput("filterProf_minrt", "Minimum RT [s]:", 0, min=0)),
									div(style = widget_style3,numericInput("filterProf_maxrt", "Maximum RT [s]:", 100000, min=0)),
									div(style = widget_style3,numericInput("filterProf_minMD", "Minimum mass defect:", -0.5, min=-0.5, max=0.5)),
									div(style = widget_style3,numericInput("filterProf_maxMD", "Maximum mass defect:", 0.5, min=-0.5, max=0.5)),
									HTML('<hr noshade="noshade" />'),				
									div(style = widget_style3,
										fluidRow(									
											column(4,
												radioButtons("filterProf_medianblind", "Filter profiles by mean sample vs. blind intensity ratio ...", c("no"="no","yes"="yes"), selected="yes", inline = TRUE),
												bsPopover("filterProf_medianblind", 
													title = "Filter profiles by the specified value for their mean sample vs. blind/blank intensity ratio across all peaks?",
													content = "This filter is only available if blind/blank peaks have not been removed from the profiles yet. Indivual intensity ratios for peaks calculated according to Settings -> Blind.", 
													placement = "top", trigger = "hover")
											),
											column(4,
												numericInput("filterProf_medianblind_value", "Mean sample vs. blind intensity ratio:", 10, min = NA, max = NA)
											)
										),
										div(style = widget_style3,radioButtons("filterProf_notblind", "... or remove all profiles which contain any blind peaks at all?", c("no"="no","yes"="yes"), inline = TRUE))
									),
									HTML('<hr noshade="noshade" />'),
									selectInput("filterProf_sort", "Rank (and filter) profile list by:", 
										choices = c("ID (increasing)",
											"mean m/z (increasing)", 
											"mean m/z (decreasing)", 
											"mean RT (increasing)", 
											"mean RT (decreasing)",										
											"minimum RT (decreasing)", 
											"maximum RT (decreasing)", 
											"maximum overall intensity (decreasing)", 
											"maximum intensity in samples (decreasing, zeros removed)", 
											"maximum intensity in blanks/blinds (decreasing, zeros removed)", 
											"mean intensity (decreasing)", 
											"mean intensity in samples (decreasing, zeros removed)", 
											"mean intensity in blanks/blinds (decreasing, zeros removed)", 
											"current trend intensity (decreasing)",										
											"past & current trend intensity (decreasing)", 
											"total peak number (decreasing)", 
											"peak number in samples (decreasing, zeros removed)", 
											"peak number in blanks/blinds (decreasing, zeros removed)", 
											"mass defect (increasing)", 
											"mass defect (decreasing)", 
											"mean sample above blind intensity (decreasing)" 
											), 
										selected = "peak number in samples (decreasing, zeros removed)", width = '80%'),
									radioButtons("filterProf_components", "Omit lower-ranked profiles with redundant intensity patterns?", c("no" = "FALSE","yes" = "TRUE"), selected = "FALSE", inline = TRUE),
									HTML('<hr noshade="noshade" />'),
									searchInput(
									  inputId = "search_profile_compound", 
									  label = "Search profile(s) of a compound with its full name or its ID (as used in the compound tables):",
									  placeholder = "Type in compound name / ID or erase to clear filter", 
									  btnSearch = icon("search"), 
									  btnReset = icon("remove"), 
									  width = "100%"
									),
									HTML('<hr noshade="noshade" />'),
									bsButton("save_profile_filter","Save above filter options as project default?",style="success"),
									#HTML('<hr noshade="noshade" />'),
									#div(style = widget_style3,numericInput("filterProf_count", "Restrict size of filtered profile list table:", 500, min=1)),
									#conditionalPanel( # IS filter				
									#		condition = "input.screen_IS == 'yes'",
									#		tags$h5("IS compounds filter:"),										
									#		HTML('<hr noshade="noshade" />')
									#),
									#conditionalPanel( # target filter				
									#		condition = "input.screen_target == 'yes'",	
									#		tags$h5("Target compounds filter:"),	
									#		HTML('<hr noshade="noshade" />')
									#),															
									HTML('<hr noshade="noshade" />'),
									bsCollapse(multiple = FALSE, open = "col2", id = "collapse2",						
										bsCollapsePanel("Intensity histogram", 
											imageOutput("profilehisto", height="auto"),
											value="test5"),
										#bsCollapsePanel("Blind intensity ratio", 
										#	imageOutput("profileboxplot", height="auto"),
										#	value="test5"),
										bsCollapsePanel("Profile list", 		
											DT::dataTableOutput("allproftable"),
											value="test5")
									)				
							),
							tabPanel("Trends",
									tags$h5("Comparison of current vs. global trends by profile ID"),
									imageOutput("boxprofile", height="auto"),
									HTML(
										'<p>
										</br></br>
										The above boxplot (grey) shows the intensity distributions of all trends of concern, listing the IDs, mean masses (m/z) and mean retention time (RT)
										of the profiles with the most intense trends on the right.
										Colored points are used to elucidate the current trend intensities from the latest input file. Note that a file is only included if 
										surviving the quality check (QC). The red dots signify profiles with intensities in the outlier range of the global (past and latest) trend intensities; 
										green dots symbolize those below. 
										</br></br>
										For more information on the individual profiles, navigate to the Single Profile tab and use the profile ID to extract further profile information.
										</p>'
									)
							),
							#tabPanel("Profile grouping",
							#	tags$h5("To be completed")	
							#
							#
							#),
							tabPanel("Single Profiles",
									tags$h5("Number of filtered profiles:"), textOutput("prof_number"),		
									HTML('<hr noshade="noshade" />'),			
									HTML('<p> Enter the ID of a profile to extract relevant information. Profile IDs are listed both in the Summary tab and the Trends tab. 
										Alternatively, sort and filter the profile list in the Summary tab and choose an entry number to show a listed profile. </p>'),							
									div(style = widget_style3,numericInput("profentry", "Entry # in (filtered, sorted) profile list:", 0)),
									div(style = widget_style3,numericInput("profID", "(... or search by any profile ID - unfiltered:)", 0)),
									div(style = widget_style3,
										dropdownButton(
											radioButtons("prof_log", "Plot logarithmic intensity?", c("no" = "no", "yes" = "yes")),
											circle = TRUE, status = "info", icon = icon("gear"), width = "300px", tooltip = tooltipOptions(title = "Click to change plot settings")
										)
									),											
									plotOutput("timeprofile", 
										dblclick = "timeprofile_dblclick",
										brush = brushOpts(
										  id = "timeprofile_brush",
										  resetOnNew = TRUE
										)
									),	
									bsCollapse(multiple = FALSE, open = "col1", id = "collapse1",
										##################################################################
										bsCollapsePanel("Compound matches", 																			
											textOutput("prof_targets"),
											HTML('<hr noshade="noshade" />'),
											textOutput("prof_ISTD")							
										),
										##################################################################
										bsCollapsePanel("Profile EICs & Peak viewer", 
											div(style = widget_style3, numericInput("profpeakID", "Peak entry #:",  min = 0,  0)),
											bsPopover("profpeakID", title = "View extracted chromatograms & peaks of the selected profile (sample files only).",
												content = "Select peaks in the order listed in the profile peak table, i.e. from latest to oldest file.", placement = "right", trigger = "hover"),
											div(style = widget_style3, textOutput("prof_peak_text")),
											plotOutput("profile_position", height = "180px"),		
											HTML('<hr noshade="noshade" />'),													
											plotOutput("profile_EIC",
												click = "profile_EIC_click",
												dblclick = "profile_EIC_dblclick",
												brush = brushOpts(
													id = "profile_EIC_brush",
													direction = c("xy"),
													resetOnNew = TRUE,
													delay = 0
												)
											),
											fluidRow(
													#column(width = 3, radioButtons("profile_EIC_norm", "Normalize intensities?", c("No"="FALSE","Yes"="TRUE"), inline = TRUE)), 
													column(width = 3, radioButtons("profile_EIC_time",  "Show RT in",  c("Seconds" = "seconds", "Minutes" = "minutes"), inline = FALSE)),
													column(width = 3, radioButtons("profile_EIC_type",  "Show chromatograms for:",  c("Full EICs with peaks" = "TRUE", "Peaks only" = "FALSE"), inline = FALSE))
											),
											value="test1"),
										##################################################################
										bsCollapsePanel("Similar component profiles",
											plotOutput("similar_profiles_plot", 
												dblclick = "similar_profiles_plot_dblclick",
												brush = brushOpts(
											  		id = "similar_profiles_plot_brush",
											  		resetOnNew = TRUE
												)
											),
											HTML('<hr noshade="noshade" />'),
											#plotOutput("similar_profiles_relations", 
											#	dblclick = "similar_profiles_relations_dblclick",
											#	brush = brushOpts(
											#  		id = "similar_profiles_relations_brush",
											#  		resetOnNew = TRUE
											#	)
											#),
											#HTML('<hr noshade="noshade" />'),
											DT::dataTableOutput("similar_profiles_table")
										),
										##################################################################
										bsCollapsePanel("Profile mass estimation",
											div(style = widget_style3,
												bsButton("dens_mass","Get mass estimates",style="success"),
												numericInput("boot_size", "Size of bootstrap sample:",min=10, 200),
												radioButtons("use_weight", "Weight by intensity?", c("no"="no","yes"="yes"))
											),
											tags$h5("m/z estimate (dark blue):"),
											textOutput("prof_mass"),
											imageOutput("massdens", height="auto"),
											imageOutput("massint", height="auto"),
											value="test3"),
										##################################################################
										bsCollapsePanel("Concentration estimation", 
											tags$h5("To be completed")								
																			
										),			
										##################################################################
										bsCollapsePanel("Profile peak table", 
											DT::dataTableOutput("oneproftable")
										)
										##################################################################
									)
							)
							#id="navbar_prof",inverse=FALSE,collapsible=TRUE,fluid=TRUE
						)
					)
				)		
            )				
        ),
        ########################################################################
        # HELP & ABOUT #########################################################
        ########################################################################
        tabPanel("Help & About",
			tags$h4("Help"),
			helpText( a("For further help, instructions, topics, examples and requests please visit the enviMass website", href="http://www.envimass.ch/",target="_blank")),
			HTML('<hr noshade="noshade" />'),							
			tags$h4("Citing enviMass"),
			HTML('<p> Martin Loos, 2016, enviMass version 3.1. Zenodo. http://doi.org/10.5281/zenodo.48164 </p> '),
			HTML('<hr noshade="noshade" />'),	
			tags$h4("Contact, author, maintainer:"),
			helpText( a("Martin Loos, mloos@looscomputing.ch", href="http://looscomputing.ch/eng/contact.htm",target="_blank") ),
			tags$h4("Contributors:"),
			helpText("Rhine monitoring station CH-Basel, Eawag CH-Duebendorf, SCIEX"),
			HTML('<hr noshade="noshade" />'),	
			helpText( a("Licensing", href="http://www.looscomputing.ch/eng/enviMass/about.htm",target="_blank") )
		)
        ########################################################################
      ),
	  HTML('<font color="white">') # hide TRUE from sourcing
    )
