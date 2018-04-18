widget_style <-
  "display: inline-block;
  vertical-align: text-top;
  padding: 7px;
  border: solid;
  border-width: 2px;
  border-radius: 4px;
  border-color: #ACA;
  background-color:#cccccc;"

widget_style2 <-
  "display: inline-block;
  vertical-align: text-top;
  padding: 7px;
  border: solid;
  border-width: 2px;
  border-radius: 4px;
  border-color: #ACA;
  background-color:#FFFFFF;"

widget_style3 <-
  "display: inline-block;
  vertical-align: text-top;
  padding: 7px;
  border: solid;
  border-width: 2px;
  border-radius: 4px;
  border-color: #FFFFFF;
  background-color:#FFFFFF;"

widget_style4 <-
  "display: inline-block;
  vertical-align: text-top;
  align: center;
  padding: 7px;
  border: solid;
  border-width: 2px;
  border-radius: 4px;
  border-color: #FFFFFF;
  background-color:#FFFFFF;"

widget_style5 <-
  "display: inline-block;
  vertical-align: text-top;
  padding: 7px;
  border: solid;
  border-width: 3px;
  border-radius: 4px;
  border-color: darkgrey;
  background-color: #cccccc;" 

widget_style6 <-
  "max-width: 600px;"  
 
widget_style7 <-
  "display: inline-block;
  vertical-align: text-top;
  padding: 7px;
  border: solid;
  border-width: 2px;
  border-radius: 4px;
  border-color: darkgrey;
  background-color: #ddff99;" 

widget_style8 <-
  "display: inline-block;
  vertical-align: text-top;
  padding: 7px;
  border: solid;
  border-width: 2px;
  border-radius: 4px;
  border-color: #003300;
  background-color: #88cc00;" 

widget_style9 <-
  "display: inline-block;
  vertical-align: text-top;
  padding: 7px;
  border: solid;
  border-width: 2px;
  border-radius: 4px;
  border-color: #003300;
  background-color: #ff9933;"

widget_style10 <-
  "display: inline-block;
  vertical-align: text-top;
  padding: 7px;
  border: solid;
  border-width: 2px;
  border-radius: 4px;
  border-color: #003300;
  background-color: #b3e6b3;"

widget_style11 <-
  "display: inline-block;
  vertical-align: text-top;
  padding: 7px;
  border: solid;
  border-width: 2px;
  border-radius: 4px;
  border-color: #003300;
  background-color: #cce0ff"
  
widget_style11 <-
  "display: inline-block;
  vertical-align: text-top;
  padding: 7px;
  border: solid;
  border-width: 2px;
  border-radius: 4px;
  border-color: lightgrey;
  background-color: lightgrey"

widget_style12 <-
  "display: inline-block;
  vertical-align: text-top;
  padding: 7px;
  border: solid;
  border-width: 1px;
  border-radius: 1px;
  border-color: green;
  background-color: white"  

  widget_style13 <-
  "display: inline-block;
  vertical-align: text-top;
  padding: 7px;
  border: solid;
  border-width: 1px;
  border-radius: 1px;
  border-color: red;
  background-color: white" 
  
shinyUI(
	fluidPage(
	################################################################################
	################################################################################
		useShinyjs(),
		useToastr(),
		useShinyFeedback(), # include shinyFeedback
		##############################################################################
		conditionalPanel( 
			condition = "output.textit == 'Waiting...'", 
			titlePanel(
				HTML('
					<p style="background-color: darkgreen">
					<font color="#FFFFFF" size="6">
					&nbsp enviMass v3.5
					</font><br/></p>'),
					windowTitle="enviMass v3.5"
				)
		),
		sidebarLayout(
			##########################################################################
			#sidebarPanel( # included in sourced script - check why!
			source("ui_sidebar.R", local=TRUE)$value,
			#)  
			##########################################################################
			mainPanel(
				tags$head(	
					tags$style(HTML("
						li a{color: black; background-color: darkgrey}; 
						.tabs-above > .nav > li[class=active] > a {
							background-color: #870000;
							color: #FFFFFF;};	
					"))
				),
				#useShinyjs(),  # Set up shinyjs
				HTML('</font>'),
				source("ui_busy.r", local=TRUE)$value,  
				HTML('</font>'),
				source("ui_mainPanel_startup.r", local=TRUE)$value,
				HTML('</font>'),	
				source("ui_mainPanel.r", local=TRUE)$value 
			, style = "float:left")
			##########################################################################  
		) 
	################################################################################
	################################################################################
	)
)





