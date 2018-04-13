output$int_norm_ISTD_pos_median <- renderPlot({
plot.new()
})
output$int_norm_ISTD_pos_counts <- renderPlot({
plot.new()
})
output$int_norm_ISTD_neg_median <- renderPlot({
plot.new()
})
output$int_norm_ISTD_neg_counts <- renderPlot({
plot.new()
})
if(file.exists(file.path(as.character(logfile[[1]]), "results", "int_norm_ISTD_pos"))){
file.remove(file.path(as.character(logfile[[1]]), "results", "int_norm_ISTD_pos"))
}
if(file.exists(file.path(as.character(logfile[[1]]), "results", "int_norm_ISTD_neg"))){
file.remove(file.path(as.character(logfile[[1]]), "results", "int_norm_ISTD_neg"))
}
