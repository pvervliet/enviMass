observeEvent(input$IS_drt2, {
feedback(
inputId = "IS_drt2",
condition = input$IS_drt2 > 10,
text = "For most applications, this tolerance would be too large",
color = "red"
)
})
observeEvent(input$tar_drt2, {
feedback(
inputId = "tar_drt2",
condition = input$tar_drt2 > 10,
text = "For most applications, this tolerance would be too large",
color = "red"
)
})
observeEvent(input$PWpath, {
PWpath_ok <- TRUE
if(
!file.exists(input$PWpath)
){
PWpath_ok <- FALSE
}
feedback(
inputId = "PWpath",
condition = !PWpath_ok,
text = "No msconvert.exe found for the specified path - please correct if you want to use msconvert for direct Thermo .raw file upload.",
color = "red"
)
})
observe({
input$isotop_license_ok
if(isolate(init$a) == "TRUE" ){
if(isolate(input$isotop_license_ok) >= 1){
removeModal(session = getDefaultReactiveDomain())
if(logfile$parameters$resolution == "OrbitrapXL,Velos,VelosPro_R60000@400"){
get_url <- "http://www.looscomputing.ch/eng/enviMass/inputs/quantiz/OrbitrapXL_Velos_VelosPro_R60K@400/quantiz"
}
if(logfile$parameters$resolution=="Q-Exactive,ExactivePlus_280K@200"){
get_url <- "http://www.looscomputing.ch/eng/enviMass/inputs/quantiz/Q-Exactive_ExactivePlus_280K@200/quantiz"
}
if(logfile$parameters$resolution=="Q-Exactive,ExactivePlus_R140000@200"){
get_url <- "http://www.looscomputing.ch/eng/enviMass/inputs/quantiz/Q-Exactive_ExactivePlus_R140K@200/quantiz"
}
if(logfile$parameters$resolution=="Q-Exactive,ExactivePlus_R70000@200"){
get_url <- "http://www.looscomputing.ch/eng/enviMass/inputs/quantiz/Q-Exactive_ExactivePlus_R70K@200/quantiz"
}
if(logfile$parameters$resolution=="Q-Exactive,ExactivePlus_R35000@200"){
get_url <- "http://www.looscomputing.ch/eng/enviMass/inputs/quantiz/Q-Exactive_ExactivePlus_R35K@200/quantiz"
}
if(logfile$parameters$resolution=="OTFusion,QExactiveHF_480000@200"){
get_url <- "http://www.looscomputing.ch/eng/enviMass/inputs/quantiz/OTFusion,QExactiveHF_480K@200/quantiz"
}
if(logfile$parameters$resolution=="OTFusion,QExactiveHF_240000@200"){
get_url <- "http://www.looscomputing.ch/eng/enviMass/inputs/quantiz/OTFusion,QExactiveHF_240K@200/quantiz"
}
if(logfile$parameters$resolution=="OTFusion,QExactiveHF_120000@200"){
get_url <- "http://www.looscomputing.ch/eng/enviMass/inputs/quantiz/OTFusion,QExactiveHF_120K@200/quantiz"
}
if(logfile$parameters$resolution=="Sciex_TripleTOF5600_R25000@200"){
get_url <- "http://www.looscomputing.ch/eng/enviMass/inputs/quantiz/Sciex_all/quantiz"
}
if(logfile$parameters$resolution=="Sciex_TripleTOF6600_R25000@200"){
get_url <- "http://www.looscomputing.ch/eng/enviMass/inputs/quantiz/Sciex_all/quantiz"
}
if(logfile$parameters$resolution=="Sciex_QTOFX500R_R25000@200"){
get_url <- "http://www.looscomputing.ch/eng/enviMass/inputs/quantiz/Sciex_all/quantiz"
}
if(logfile$parameters$resolution=="Agilent_QTOF6550_low_extended2GHz_highRes"){
get_url <- "http://www.looscomputing.ch/eng/enviMass/inputs/quantiz/Agilent_QTOF6550_low_extended2GHz_highRes/quantiz"
}
if(logfile$parameters$resolution=="Agilent_QTOF6550_low_highRes4GHz_highRes"){
get_url <- "http://www.looscomputing.ch/eng/enviMass/inputs/quantiz/Agilent_QTOF6550_low_highRes4GHz_highRes/quantiz"
}
dest_file <- file.path(logfile[[1]], "dataframes", "quantiz")
url_quantiz <- try(download.file(url = get_url, destfile = dest_file, mode = "wb"))
if(class(url_quantiz) == "try-error"){
cat("\n Download of missing isotopologue space failed.")
say <- "Unable to download missing isotopologue space. Please either check your internet connection and retry or
proceed manually as described on www.enviMass.ch -> Data input -> Download available isotopologue spaces."
shinytoastr::toastr_error(say, title = "Isotopologue space download:", closeButton = TRUE, position = c("top-center"), timeOut = 0);
}else{
cat("\n Download of missing isotopologue space completed.\n")
load_quantiz <- try(load(file.path(logfile[[1]], "dataframes", "quantiz")) )
if(class(load_quantiz) == "try-error"){
cat("\n Loading of missing isotopologue space failed.\n")
say <- "Loading failure of downloaded isotopologue space. Please proceed manually as described on www.enviMass.ch -> Data input -> Download available isotopologue spaces."
shinytoastr::toastr_error(say, title = "Isotopologue space download:", closeButton = TRUE, position = c("top-center"), timeOut = 0);
}
}
}
}
})
