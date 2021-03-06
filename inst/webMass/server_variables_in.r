if(any(ls()=="logfile")){stop("\n illegal logfile detected #1 in server_variabels_in.r!")}
path <- "ui_mainPanel.r"
enviMass::LIRYM3506Y(path, logfile, session)
updateCheckboxGroupInput(session, "adducts_pos", selected = as.character(logfile$adducts_pos))
updateCheckboxGroupInput(session, "adducts_neg", selected = as.character(logfile$adducts_neg))
updateCheckboxGroupInput(session, "adducts_pos_group", selected = as.character(logfile$adducts_pos_group))
updateCheckboxGroupInput(session, "adducts_neg_group", selected = as.character(logfile$adducts_neg_group))
updateTextInput(session, inputId = "PWpath", value = logfile$PW)
if(any( (measurements[,"ID"] != "-") & (measurements[,"Mode"] == "positive") & (measurements[,"Type"] != "sample"))){
IDs_pos <- measurements[
(measurements[,"Mode"] == "positive") & (measurements[,"Type"] != "sample")
,1]
names_pos <- measurements[
(measurements[,"Mode"] == "positive") & (measurements[,"Type"] != "sample")
,2]
IDs_pos <- paste(IDs_pos, names_pos, sep = " - ")
if(any(logfile[["Positive_subtraction_files"]] != "FALSE")){
select_pos <- logfile[["Positive_subtraction_files"]]
select_pos <- select_pos[select_pos!="FALSE"]
select_pos <- select_pos[!is.na(match(select_pos,IDs_pos))]
}else{
select_pos <- NULL
}
updateCheckboxGroupInput(session, inputId = "files_pos_select_subtract", label = "", choices = IDs_pos, selected = select_pos)
}
if(any( (measurements[,"ID"] != "-") & (measurements[,"Mode"] == "negative") & (measurements[,"Type"] != "sample"))){
IDs_neg <- measurements[
(measurements[,"Mode"]=="negative") & (measurements[,"Type"] != "sample")
,1]
names_neg <- measurements[
(measurements[,"Mode"]=="negative") & (measurements[,"Type"] != "sample")
,2]
IDs_neg <- paste(IDs_neg,names_neg,sep=" - ")
if(any(logfile[["Negative_subtraction_files"]] != "FALSE")){
select_neg <- logfile[["Negative_subtraction_files"]]
select_neg <- select_neg[select_neg != "FALSE"]
select_neg <- select_neg[!is.na(match(select_neg, IDs_neg))]
}else{
select_neg<-NULL
}
updateCheckboxGroupInput(session,inputId="files_neg_select_subtract", label="", choices = IDs_neg, selected = select_neg)
}
if(is.character(logfile$method_setup)){
output$heads_summary_existing <- renderTable(as.data.frame("No existing method available", optional = TRUE))
}else{
output$heads_summary_existing <- renderTable(logfile$method_setup)
}
if(any(ls()=="logfile")){stop("\n illegal logfile detected #2 in server_variabels_in.r!")}
