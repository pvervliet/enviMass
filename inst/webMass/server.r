options(shiny.maxRequestSize=2000*1024^2)
verbose<-FALSE
shinyServer(function(input, output, session){
cat("\n I run in ");print(environment())
in_envir <<- environment()
if(any(ls()=="logfile")){stop("\n illegal logfile detected in server.r #1")}
if(!any(objects(envir = as.environment(".GlobalEnv")) == "isotopes")){data(isotopes, package = "enviPat", envir = as.environment(".GlobalEnv"))}
updateSelectInput(session, inputId = "atom_bounds_this", choices = unique(isotopes[1:295,1]), selected = c("C","H","N","O","Cl","Br"))
if(!any(objects(envir = as.environment(".GlobalEnv")) == "adducts")){data(adducts, package = "enviPat", envir = as.environment(".GlobalEnv"))}
if(!any(objects(envir = as.environment(".GlobalEnv")) == "resolution_list")){data(resolution_list, package = "enviPat", envir = as.environment(".GlobalEnv"))}
if(any(names(resolution_list) == "Elite/R240000@400")){
shinyjs::info("library enviPat is not up to date - or you have loaded an old workspace containing old enviPat resolution data lists. Update enviPat and clean your workspace before continuing with enviMass!");
}
tried<-try(getVolumes()(),silent=FALSE)
if(!inherits(tried,"try-error")){
shinyFileChoose(input, "pro_dir3", updateFreq = 30000, roots = getVolumes(), filetypes = c("emp"))
shinyFileSave(input, "download_IS", updateFreq = 30000, roots = getVolumes() )
shinyFileSave(input, "download_target", updateFreq = 30000, roots = getVolumes() )
}else{
createAlert(session, anchorId = "alert_4", alertId="a4", title = NULL, content="logfile select via explorer disabled, please use folder path input instead",style = "alarm",append=FALSE,dismiss=TRUE)
}
output$textit <- renderText("Waiting...")
output$dowhat <- renderText("Open")
output$sel_meas_comp_state <- renderText("")
output$isotable <- renderTable(isotopes)
updateCheckboxGroupInput(session, "adducts_pos", "Positive ions:", choices =  as.character(adducts[adducts[,6]=="positive",1]), selected = as.character(adducts[adducts[,6]=="positive",1][1]))
updateCheckboxGroupInput(session, "adducts_neg", "Negative ions:", choices =  as.character(adducts[adducts[,6]=="negative",1]), selected = as.character(adducts[adducts[,6]=="negative",1][1]))
updateCheckboxGroupInput(session, "adducts_pos_group", "Positive mode:", choices =  as.character(adducts[adducts[,6]=="positive",1]), selected = as.character(adducts[adducts[,6]=="positive",1][1]))
updateCheckboxGroupInput(session, "adducts_neg_group", "Negative mode:", choices =  as.character(adducts[adducts[,6]=="negative",1]), selected = as.character(adducts[adducts[,6]=="negative",1][1]))
updateSelectInput(session, "resolution", "Instrument resolution:", choices =  names(resolution_list), selected= (names(resolution_list)[1]))
init <- reactiveValues() 	# reactive value to indicate ...
init$a <- "FALSE"
init$b <- 1
logobusy <- list(src = "circ.gif");
output$logobusy1 <- renderImage(logobusy, deleteFile = FALSE);
output$logobusy2 <- renderImage(logobusy, deleteFile = FALSE);
output$logobusy3 <- renderImage(logobusy, deleteFile = FALSE);
peakrot <- list(src = "peakrot.gif");
output$peakrot <- renderImage(peakrot, deleteFile = FALSE);
source("server_startup.r", local = TRUE)
source("server_obs_Add.r", local = TRUE)
source("server_obs_res_meas.r", local = TRUE)
source("server_variables_out.r", local = TRUE)
source("server_obs_screening.r", local = TRUE)
source("server_obs_plots.r", local = TRUE)
source("server_obs_calibration.r", local = TRUE)
source("server_obs_components.r", local = TRUE)
source("server_obs_profiles.r", local = TRUE)
source("server_force.r", local = TRUE)
source("server_warnings.r", local = TRUE)
source("server_comparison.r", local = TRUE)
source("server_calculation.r", local = TRUE)
observe({
input$Restart
if(input$Restart){
output$textit<<-renderText("Waiting...")
source("server_cleaner.R", local=TRUE);
cat("Restart\n")
}
})
observe({
input$Exit
if(input$Exit){
source("server_cleaner.R", local=TRUE);
stopApp(returnValue="Quit enviMass browser session")
}
})
})
