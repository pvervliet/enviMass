observe({
input$save_comparison
if(isolate(init$a)=="TRUE"){
if(logfile$parameters$verbose) cat("\n In comparison: saving")
if(length(logfile$comparisons)){
if(!any(names(logfile$comparisons) == isolate(input$comparison_name))){
found_empty <- FALSE
for(i in 1:length(logfile$comparisons)){
if(length(logfile$comparisons[[i]]) == 0){
found_empty <- TRUE
break;
}
}
}else{
found_empty <- TRUE
i <- which(names(logfile$comparisons) == isolate(input$comparison_name))
}
if(found_empty){
at_comparison <- i
}else{
at_comparison <- length(logfile$comparisons) + 1
}
}else{
at_comparison <- 1
}
logfile$comparisons[[at_comparison]] <<- list()
logfile$comparisons[[at_comparison]][[1]] <<- isolate(input$apply_comparison)
logfile$comparisons[[at_comparison]][[2]] <<- isolate(input$comparison_editor)
names(logfile$comparisons)[at_comparison] <<- isolate(input$comparison_name)
for_names <- names(logfile$comparisons)
for_names <- for_names[for_names != ""]
updateSelectInput(session, "load_comparison", choices = c("None", for_names), selected = isolate(input$comparison_name))
if(isolate(input$apply_comparison)){
enviMass::RNJJT9347D(logfile, down = "comparison")
output$summa_html <- renderText(enviMass::UCBRR9756L(logfile$summary, logfile$Tasks_to_redo));
}
save(logfile, file = file.path(as.character(logfile[["project_folder"]]), "logfile.emp"));
if(any(ls()=="logfile")){stop("\n illegal logfile detected #2 in server_updates.r!")}
}
})
observe({
input$check_comparison
if(isolate(init$a)=="TRUE"){
if(logfile$parameters$verbose) cat("\n In comparison: checking")
isolate(input$apply_comparison)
if(any(ls()=="logfile")){stop("\n illegal logfile detected #2 in server_updates.r!")}
}
})
observe({
input$load_comparison
if(isolate(init$a)=="TRUE"){
if(logfile$parameters$verbose) cat("\n In comparison: loading")
at_comparison <- which(names(logfile$comparisons) == isolate(input$load_comparison))
if(length(at_comparison)){
isolate(updateAceEditor(session, editorId = "comparison_editor", value = logfile$comparisons[[at_comparison]][[2]]))
updateTextInput(session, inputId = "comparison_name", value = names(logfile$comparisons)[at_comparison])
updateCheckboxInput(session, inputId = "apply_comparison", value = logfile$comparisons[[at_comparison]][[1]])
}else{
isolate(updateAceEditor(session, editorId = "comparison_editor", value = "Add new comparison here"))
}
if(any(ls()=="logfile")){stop("\n illegal logfile detected #2 in server_updates.r!")}
}
})
observe({
input$delete_comparison
if(isolate(init$a)=="TRUE"){
if(logfile$parameters$verbose) cat("\n In comparison: deleting")
at_comparison <- which(names(logfile$comparisons) == isolate(input$load_comparison))
if(length(at_comparison)){
if(logfile$comparisons[[at_comparison]][[1]] == TRUE){
enviMass::RNJJT9347D(logfile, down = "comparison")
output$summa_html <- renderText(enviMass::UCBRR9756L(logfile$summary, logfile$Tasks_to_redo));
}
logfile$comparisons[[at_comparison]] <<- list()
names(logfile$comparisons)[at_comparison] <<- ""
for_names <- names(logfile$comparisons)
for_names <- for_names[for_names != ""]
updateSelectInput(session, "load_comparison", choices = c("None", for_names), selected = "None")
}
save(logfile, file = file.path(as.character(logfile[["project_folder"]]), "logfile.emp"));
if(any(ls()=="logfile")){stop("\n illegal logfile detected #2 in server_updates.r!")}
}
})
