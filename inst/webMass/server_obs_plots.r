ranges_int_norm_IS_pos <- reactiveValues(xlim = FALSE)
ranges_int_norm_IS_neg <- reactiveValues(xlim = FALSE)
ranges_int_box_pos <- reactiveValues(xlim = FALSE, redo = 0)
ranges_int_box_neg <- reactiveValues(xlim = FALSE, redo = 0)
ranges_int_quantiles_pos <- reactiveValues(xlim = FALSE, ylim = FALSE)
ranges_int_quantiles_neg <- reactiveValues(xlim = FALSE, ylim = FALSE)
observeEvent(input$int_box_pos_dblclick, {
brush <- isolate(input$int_box_pos_brush)
cat("\n Zoom in/out")
if (!is.null(brush)) {
cat("\n Zoom in_1")
isolate(ranges_int_box_pos$xlim <- c(brush$xmin, brush$xmax))
} else {
cat("\n Zoom out full_1")
isolate(ranges_int_box_pos$xlim <- FALSE)
}
isolate(ranges_int_box_pos$redo <- (ranges_int_box_pos$redo + 1))
})
observeEvent(input$int_box_pos_click, {
cat("\n Zoom out part_1_a")
brush <- isolate(input$int_box_pos_brush)
if (is.null(brush)) {
cat("\n Zoom out part_1_b")
if(isolate(ranges_int_box_pos$xlim[1]) != FALSE){
old_range_mass <- abs(isolate(ranges_int_box_pos$xlim[2] - ranges_int_box_pos$xlim[1]))
isolate(ranges_int_box_pos$xlim[1] <- ranges_int_box_pos$xlim[1] - .3 * old_range_mass)
isolate(ranges_int_box_pos$xlim[2] <- ranges_int_box_pos$xlim[2] + .3* old_range_mass)
}
isolate(ranges_int_box_pos$redo <- (ranges_int_box_pos$redo + 1))
}else{
cat("\n Doing hover - and nothing")
}
})
observeEvent(input$int_box_neg_dblclick, {
brush <- isolate(input$int_box_neg_brush)
cat("\n Zoom in/out")
if (!is.null(brush)) {
cat("\n Zoom in_1")
isolate(ranges_int_box_neg$xlim <- c(brush$xmin, brush$xmax))
} else {
cat("\n Zoom out full_1")
isolate(ranges_int_box_neg$xlim <- FALSE)
}
isolate(ranges_int_box_neg$redo <- (ranges_int_box_neg$redo + 1))
})
observeEvent(input$int_box_neg_click, {
cat("\n Zoom out part_1_a")
brush <- isolate(input$int_box_neg_brush)
if (is.null(brush)) {
cat("\n Zoom out part_1_b")
if(isolate(ranges_int_box_neg$xlim[1]) != FALSE){
old_range_mass <- abs(isolate(ranges_int_box_neg$xlim[2] - ranges_int_box_neg$xlim[1]))
isolate(ranges_int_box_neg$xlim[1] <- ranges_int_box_neg$xlim[1] - .3 * old_range_mass)
isolate(ranges_int_box_neg$xlim[2] <- ranges_int_box_neg$xlim[2] + .3* old_range_mass)
}
isolate(ranges_int_box_neg$redo <- (ranges_int_box_neg$redo + 1))
}else{
cat("\n Doing hover - and nothing")
}
})
observe({
ranges_int_box_pos$redo
if(isolate(init$a)=="TRUE"){
if(file.exists(file.path(as.character(logfile[[1]]), "results", "int_distrib"))){
if(!any(objects(envir = as.environment(".GlobalEnv")) == "int_distrib")){
load(file.path(as.character(logfile[[1]]), "results", "int_distrib"), envir = as.environment(".GlobalEnv"))
}
if(!any(objects(envir = as.environment(".GlobalEnv")) == "measurements")){
measurements <- read.csv(file = file.path(logfile[[1]], "dataframes", "measurements"), colClasses = "character");
}
output$int_box_pos <- renderPlot({
par(mar = c(5.8, 4.5, .9, 8))
enviMass:::MSBYD6677H(
int_distrib,
ion_mode = "positive",
measurements = measurements,
xlim = isolate(ranges_int_box_pos$xlim)
)
},res = 100)
}else{
output$int_box_pos <- renderPlot({
plot.new()
})
}
}
})
observe({
ranges_int_box_neg$redo
if(isolate(init$a)=="TRUE"){
if(file.exists(file.path(as.character(logfile[[1]]), "results", "int_distrib"))){
if(!any(objects(envir = as.environment(".GlobalEnv")) == "int_distrib")){
load(file.path(as.character(logfile[[1]]), "results", "int_distrib"), envir = as.environment(".GlobalEnv"))
}
if(!any(objects(envir = as.environment(".GlobalEnv")) == "measurements")){
measurements <- read.csv(file = file.path(logfile[[1]], "dataframes", "measurements"), colClasses = "character");
}
output$int_box_neg <- renderPlot({
par(mar = c(5.8, 4.5, .9, 8))
enviMass:::MSBYD6677H(
int_distrib,
ion_mode = "negative",
measurements = measurements,
xlim = isolate(ranges_int_box_neg$xlim)
)
},res = 100)
}else{
output$int_box_neg <- renderPlot({
plot.new()
})
}
}
})
observeEvent(input$int_maxmed_pos_dblclick, {
brush <- isolate(input$int_maxmed_pos_brush)
if (!is.null(brush)) {
cat("\n Zoom in_1")
isolate(ranges_int_quantiles_pos$xlim <- c(brush$xmin, brush$xmax))
isolate(ranges_int_quantiles_pos$ylim <- c(brush$ymin, brush$ymax))
} else {
cat("\n Zoom out full_1")
isolate(ranges_int_quantiles_pos$xlim <- FALSE)
isolate(ranges_int_quantiles_pos$ylim <- FALSE)
}
})
observeEvent(input$int_maxmed_pos_click, {
cat("\n Zoom out part_1_a")
brush <- isolate(input$int_maxmed_pos_brush)
if (is.null(brush)) {
cat("\n Zoom out part_1_b")
if(isolate(ranges_int_quantiles_pos$xlim[1]) != FALSE){
old_range_mass <- abs(isolate(ranges_int_quantiles_pos$xlim[2] - ranges_int_quantiles_pos$xlim[1]))
isolate(ranges_int_quantiles_pos$xlim[1] <- ranges_int_quantiles_pos$xlim[1] - .3 * old_range_mass)
isolate(ranges_int_quantiles_pos$xlim[2] <- ranges_int_quantiles_pos$xlim[2] + .3* old_range_mass)
}
if(isolate(ranges_int_quantiles_pos$ylim[1]) != FALSE){
old_range_mass <- abs(isolate(ranges_int_quantiles_pos$ylim[2] - ranges_int_quantiles_pos$ylim[1]))
isolate(ranges_int_quantiles_pos$ylim[1] <- ranges_int_quantiles_pos$ylim[1] - .3 * old_range_mass)
isolate(ranges_int_quantiles_pos$ylim[2] <- ranges_int_quantiles_pos$ylim[2] + .3* old_range_mass)
}
}else{
cat("\n Doing hover - and nothing")
}
})
observeEvent(input$int_maxmed_neg_dblclick, {
brush <- isolate(input$int_maxmed_neg_brush)
if (!is.null(brush)) {
cat("\n Zoom in_1")
isolate(ranges_int_quantiles_neg$xlim <- c(brush$xmin, brush$xmax))
isolate(ranges_int_quantiles_neg$ylim <- c(brush$ymin, brush$ymax))
} else {
cat("\n Zoom out full_1")
isolate(ranges_int_quantiles_neg$xlim <- FALSE)
isolate(ranges_int_quantiles_neg$ylim <- FALSE)
}
})
observeEvent(input$int_maxmed_neg_click, {
cat("\n Zoom out part_1_a")
brush <- isolate(input$int_maxmed_neg_brush)
if (is.null(brush)) {
cat("\n Zoom out part_1_b")
if(isolate(ranges_int_quantiles_neg$xlim[1]) != FALSE){
old_range_mass <- abs(isolate(ranges_int_quantiles_neg$xlim[2] - ranges_int_quantiles_neg$xlim[1]))
isolate(ranges_int_quantiles_neg$xlim[1] <- ranges_int_quantiles_neg$xlim[1] - .3 * old_range_mass)
isolate(ranges_int_quantiles_neg$xlim[2] <- ranges_int_quantiles_neg$xlim[2] + .3* old_range_mass)
}
if(isolate(ranges_int_quantiles_neg$ylim[1]) != FALSE){
old_range_mass <- abs(isolate(ranges_int_quantiles_neg$ylim[2] - ranges_int_quantiles_neg$ylim[1]))
isolate(ranges_int_quantiles_neg$ylim[1] <- ranges_int_quantiles_neg$ylim[1] - .3 * old_range_mass)
isolate(ranges_int_quantiles_neg$ylim[2] <- ranges_int_quantiles_neg$ylim[2] + .3* old_range_mass)
}
}else{
cat("\n Doing hover - and nothing")
}
})
observe({
ranges_int_quantiles_pos$xlim
ranges_int_quantiles_pos$ylim
if(isolate(init$a)=="TRUE"){
if(file.exists(file.path(as.character(logfile[[1]]), "results", "int_distrib"))){
if(!any(objects(envir = as.environment(".GlobalEnv")) == "int_distrib")){
load(file.path(as.character(logfile[[1]]), "results", "int_distrib"), envir = as.environment(".GlobalEnv"))
}
if(!any(objects(envir = as.environment(".GlobalEnv")) == "measurements")){
measurements <- read.csv(file = file.path(logfile[[1]], "dataframes", "measurements"), colClasses = "character");
}
output$int_quantiles_pos <- renderPlot({
par(mar = c(4.5, 4.5, .9, 8))
enviMass:::MSBYD6677H(
int_distrib,
ion_mode = "positive",
measurements = measurements,
what = "quantiles_distrib",
maxit_1_lim = isolate(ranges_int_quantiles_pos$xlim),
maxit_2_lim = isolate(ranges_int_quantiles_pos$ylim)
)
},res = 100)
output$int_maxmed_pos <- renderPlot({
par(mar = c(4.5, 4.5, .9, 8))
enviMass:::MSBYD6677H(
int_distrib,
ion_mode = "positive",
measurements = measurements,
what = "quantiles_out",
maxit_1_lim = isolate(ranges_int_quantiles_pos$xlim),
maxit_2_lim = isolate(ranges_int_quantiles_pos$ylim)
)
},res = 100)
}else{
output$int_quantiles_pos <- renderPlot({
plot.new(); plot.window(xlim = c(0, 10), ylim =c(0, 10)); text(5, 5, "Data not available")
})
output$int_maxmed_pos <- renderPlot({
plot.new(); plot.window(xlim = c(0, 10), ylim =c(0, 10)); text(5, 5, "Data not available")
})
}
}
})
observe({
ranges_int_quantiles_neg$xlim
ranges_int_quantiles_neg$ylim
if(isolate(init$a)=="TRUE"){
if(file.exists(file.path(as.character(logfile[[1]]), "results", "int_distrib"))){
if(!any(objects(envir = as.environment(".GlobalEnv")) == "int_distrib")){
load(file.path(as.character(logfile[[1]]), "results", "int_distrib"), envir = as.environment(".GlobalEnv"))
}
if(!any(objects(envir = as.environment(".GlobalEnv")) == "measurements")){
measurements <- read.csv(file = file.path(logfile[[1]], "dataframes", "measurements"), colClasses = "character");
}
output$int_quantiles_neg <- renderPlot({
par(mar = c(4.5, 4.5, .9, 8))
enviMass:::MSBYD6677H(
int_distrib,
ion_mode = "negative",
measurements = measurements,
what = "quantiles_distrib",
maxit_1_lim = isolate(ranges_int_quantiles_neg$xlim),
maxit_2_lim = isolate(ranges_int_quantiles_neg$ylim)
)
},res = 100)
output$int_maxmed_neg <- renderPlot({
par(mar = c(4.5, 4.5, .9, 8))
enviMass:::MSBYD6677H(
int_distrib,
ion_mode = "negative",
measurements = measurements,
what = "quantiles_out",
maxit_1_lim = isolate(ranges_int_quantiles_neg$xlim),
maxit_2_lim = isolate(ranges_int_quantiles_neg$ylim)
)
},res = 100)
}else{
output$int_quantiles_neg <- renderPlot({
plot.new(); plot.window(xlim = c(0, 10), ylim =c(0, 10)); text(5, 5, "Data not available")
})
output$int_maxmed_neg <- renderPlot({
plot.new(); plot.window(xlim = c(0, 10), ylim =c(0, 10)); text(5, 5, "Data not available")
})
}
}
})
observeEvent(input$int_norm_ISTD_pos_median_dblclick, {
brush <- isolate(input$int_norm_ISTD_pos_median_brush)
if (!is.null(brush)) {
cat("\n Zoom in_1")
isolate(ranges_int_norm_IS_pos$xlim <- c(brush$xmin, brush$xmax))
} else {
cat("\n Zoom out full_1")
isolate(ranges_int_norm_IS_pos$xlim <- FALSE)
}
})
observeEvent(input$int_norm_ISTD_pos_median_click, {
cat("\n Zoom out part_1_a")
brush <- isolate(input$int_norm_ISTD_pos_median_brush)
if (is.null(brush)) {
cat("\n Zoom out part_1_b")
if(isolate(ranges_int_norm_IS_pos$xlim[1]) != FALSE){
old_range_mass<-abs(isolate(ranges_int_norm_IS_pos$xlim[2] - ranges_int_norm_IS_pos$xlim[1]))
isolate(ranges_int_norm_IS_pos$xlim[1] <- ranges_int_norm_IS_pos$xlim[1] - .3 * old_range_mass)
isolate(ranges_int_norm_IS_pos$xlim[2] <- ranges_int_norm_IS_pos$xlim[2] + .3 * old_range_mass)
}
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
})
observeEvent(input$int_norm_ISTD_pos_counts_click, {
cat("\n Zoom out part_1_a")
brush <- isolate(input$int_norm_ISTD_pos_counts_brush)
if (is.null(brush)) {
cat("\n Zoom out part_1_b")
if(isolate(ranges_int_norm_IS_pos$xlim[1]) != FALSE){
old_range_mass<-abs(isolate(ranges_int_norm_IS_pos$xlim[2] - ranges_int_norm_IS_pos$xlim[1]))
isolate(ranges_int_norm_IS_pos$xlim[1] <- ranges_int_norm_IS_pos$xlim[1] - .3 * old_range_mass)
isolate(ranges_int_norm_IS_pos$xlim[2] <- ranges_int_norm_IS_pos$xlim[2] + .3* old_range_mass)
}
}else{
cat("\n Doing hover - and nothing")
}
})
observeEvent(input$int_norm_ISTD_neg_median_dblclick, {
brush <- isolate(input$int_norm_ISTD_neg_median_brush)
if (!is.null(brush)) {
cat("\n Zoom in_1")
isolate(ranges_int_norm_IS_neg$xlim <- c(brush$xmin, brush$xmax))
} else {
cat("\n Zoom out full_1")
isolate(ranges_int_norm_IS_neg$xlim <- FALSE)
}
})
observeEvent(input$int_norm_ISTD_neg_median_click, {
cat("\n Zoom out part_1_a")
brush <- isolate(input$int_norm_ISTD_neg_median_brush)
if (is.null(brush)) {
cat("\n Zoom out part_1_b")
if(isolate(ranges_int_norm_IS_neg$xlim[1]) != FALSE){
old_range_mass<-abs(isolate(ranges_int_norm_IS_neg$xlim[2] - ranges_int_norm_IS_neg$xlim[1]))
isolate(ranges_int_norm_IS_neg$xlim[1] <- ranges_int_norm_IS_neg$xlim[1] - .3 * old_range_mass)
isolate(ranges_int_norm_IS_neg$xlim[2] <- ranges_int_norm_IS_neg$xlim[2] + .3 * old_range_mass)
}
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
})
observeEvent(input$int_norm_ISTD_neg_counts_click, {
cat("\n Zoom out part_1_a")
brush <- isolate(input$int_norm_ISTD_neg_counts_brush)
if (is.null(brush)) {
cat("\n Zoom out part_1_b")
if(isolate(ranges_int_norm_IS_neg$xlim[1]) != FALSE){
old_range_mass<-abs(isolate(ranges_int_norm_IS_neg$xlim[2] - ranges_int_norm_IS_neg$xlim[1]))
isolate(ranges_int_norm_IS_neg$xlim[1] <- ranges_int_norm_IS_neg$xlim[1] - .3 * old_range_mass)
isolate(ranges_int_norm_IS_neg$xlim[2] <- ranges_int_norm_IS_neg$xlim[2] + .3 * old_range_mass)
}
}else{
cat("\n Doing hover - and nothing")
}
})
observe({
ranges_int_norm_IS_pos$xlim
if(isolate(init$a)=="TRUE"){
if(file.exists(file.path(as.character(logfile[[1]]), "results", "int_norm_ISTD_pos"))){
if(!any(objects(envir = as.environment(".GlobalEnv")) == "int_norm_ISTD_pos")){
load(file.path(as.character(logfile[[1]]), "results", "int_norm_ISTD_pos"), envir = as.environment(".GlobalEnv"))
}
output$int_norm_ISTD_pos_median <- renderPlot({
par(mar = c(.2, 4.5, .9, 8))
enviMass:::TICFG4083H(
int_norm_ISTD = int_norm_ISTD_pos,
logfile = logfile,
xlim = isolate(ranges_int_norm_IS_pos$xlim),
what = "normalization"
)
},res = 100)
output$int_norm_ISTD_pos_counts <- renderPlot({
par(mar = c(4.5, 4.5, .9, 8))
enviMass:::TICFG4083H(
int_norm_ISTD = int_norm_ISTD_pos,
logfile = logfile,
xlim = isolate(ranges_int_norm_IS_pos$xlim),
what = "counts"
)
},res = 100)
}else{
output$int_norm_ISTD_pos_median <- renderPlot({
plot.new(); plot.window(xlim = c(0, 10), ylim =c(0, 10)); text(5, 5, "Data not available")
})
output$int_norm_ISTD_pos_counts <- renderPlot({
plot.new(); plot.window(xlim = c(0, 10), ylim =c(0, 10)); text(5, 5, "Data not available")
})
}
}
})
observe({
ranges_int_norm_IS_neg$xlim
if(isolate(init$a)=="TRUE"){
if(file.exists(file.path(as.character(logfile[[1]]), "results", "int_norm_ISTD_neg"))){
if(!any(objects(envir = as.environment(".GlobalEnv")) == "int_norm_ISTD_neg")){
load(file.path(as.character(logfile[[1]]), "results", "int_norm_ISTD_neg"), envir = as.environment(".GlobalEnv"))
}
output$int_norm_ISTD_neg_median <- renderPlot({
par(mar = c(.2, 4.5, .9, 8))
enviMass:::TICFG4083H(
int_norm_ISTD = int_norm_ISTD_neg,
logfile = logfile,
xlim = isolate(ranges_int_norm_IS_neg$xlim),
what = "normalization"
)
},res = 100)
output$int_norm_ISTD_neg_counts <- renderPlot({
par(mar = c(4.5, 4.5, .9, 8))
enviMass:::TICFG4083H(
int_norm_ISTD = int_norm_ISTD_neg,
logfile = logfile,
xlim = isolate(ranges_int_norm_IS_neg$xlim),
what = "counts"
)
},res = 100)
}else{
output$int_norm_ISTD_neg_median <- renderPlot({
plot.new(); plot.window(xlim = c(0, 10), ylim =c(0, 10)); text(5, 5, "Data not available")
})
output$int_norm_ISTD_neg_counts <- renderPlot({
plot.new(); plot.window(xlim = c(0, 10), ylim =c(0, 10)); text(5, 5, "Data not available")
})
}
}
})
