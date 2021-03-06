observe({
input$newit
if(isolate(input$newit)){
if(any(ls()=="logfile")){stop("\n illegal logfile detected #1 in server_startup.r!")}
source("server_cleaner.R", local=TRUE);
createAlert(session, anchorId="alert_3", alertId="a3",
title = "Getting started?",
content="<p> Using the below panels, (1) add <span style='color: black'> Files</span>, (2) optionally add <span style='color: black'> Compounds</span>, (3) select the <span style='color: black'>
Workflow options </span> and (4) adjust the parameter <span style='color: black'> Settings </span> (make sure your data are suitable for enviPick peak picking, best by trying the seperate
<a href='http://cran.r-project.org/web/packages/enviPick/index.html' target='_blank'> package</a> beforehand). Then, press <span style='color: red'> Calculate </span> (a sidebar message will
tell you if your project is consistent) and wait for the <span style='color: black'> Results </span> to appear (depending on the number of files, this may take a while; especially if peak
picking has not been run yet).
<p>You can exit your project and reopen it later to add new files and compounds or to change settings. enviMass will then adjust recalculations to avoid unnecessary workflow
steps (e.g. to not redo peak picking for all files if only a few new have been added).</p>
",
style="info",append=FALSE,dismiss=TRUE)
if(!exists("IS",envir=as.environment(".GlobalEnv"))){data(IS,package="enviMass")}
if(!exists("targets",envir=as.environment(".GlobalEnv"))){data(targets,package="enviMass")}
say_path <- enviMass::JAYQN1284Q(isolate(input$pro_name),isolate(input$pro_dir))
if(say_path!="Project path ok"){
createAlert(session, anchorId="failed_new", alertId = "failed_new_id", title = "Invalid project path",
content = "Project already exists, the specified path is invalid or you lack permissions.",
style = "danger", dismiss = TRUE, append = FALSE)
cat("Invalid - project already exists or path invalid \n")
shinyjs::info(say_path)
}else{
logfile_path<-enviMass::KXUVE9576N(isolate(input$pro_name), isolate(input$pro_dir), IS, targets);
if(logfile_path != "FALSE"){
output$textit <- renderText(as.character(logfile_path));
load(logfile_path, envir = as.environment(".GlobalEnv"));
output$summa_html <- renderText(enviMass::UCBRR9756L(logfile$summary));
output$dowhat <- renderText("Started new project");
output$IS <- DT::renderDataTable(read.table(file=file.path(logfile[[1]],"dataframes","IS.txt"),header=TRUE,sep="\t",colClasses = "character"));
output$targets <- DT::renderDataTable(read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character"));
measurements <- read.csv(file=file.path(logfile$project_folder,"dataframes","measurements"),colClasses = "character")
measurements_tab <- measurements
measurements_tab[,1] <- as.numeric(measurements_tab[,1])
measurements_tab <- DT::datatable(
measurements_tab[,c("ID","Name","Type","Mode","Place","Date","Time","include","profiled","tag1",
"tag2","tag3","date_end","time_end","ID_2")],
extensions = c('Buttons','ColReorder','FixedHeader'),
rownames = FALSE,
options = list(
lengthMenu = list(c(25, 50, 100, -1), list('25', '50', '100', 'All')),
fixedHeader = TRUE,
ordering = T,
dom = 'Blfrtip',
buttons = c('excel', 'csv', 'colvis'),#buttons = c('excel', 'pdf', 'print', 'csv'),
scrollX = TRUE,
colReorder = TRUE
)
)
output$measurements <- DT::renderDataTable(measurements_tab);
output$sel_meas_comp_state<-renderText("")
enviMass::DODNL4951M(session)
output$int_box_pos <- renderPlot({
plot.new()
})
output$int_box_neg <- renderPlot({
plot.new()
})
path=file.path(logfile$project_folder,"pics","int_distr_pos")
png(filename = path, bg = "white")
plot.new();plot.window(xlim=c(0,1),ylim=c(0,1));text(0.5,0.5,"nothing selected \n or not available",cex=1)
dev.off()
expr3p<-list(src=path)
output$pic_int_distr_pos<-renderImage(expr3p, deleteFile = FALSE)
path=file.path(logfile$project_folder,"pics","int_distr_neg")
png(filename = path, bg = "white")
plot.new();plot.window(xlim=c(0,1),ylim=c(0,1));text(0.5,0.5,"nothing selected \n or not available",cex=1)
dev.off()
expr3n<-list(src=path)
output$pic_int_distr_neg<-renderImage(expr3n, deleteFile = FALSE)
path=file.path(logfile[[1]],"pics","recal_none")
png(filename = path, bg = "white")
plot.new();plot.window(xlim=c(0,1),ylim=c(0,1));text(0.5,0.5,"nothing selected \n or not available",cex=1)
dev.off()
exprrec<-list(src=path)
output$recal_pic<-renderImage(exprrec, deleteFile = FALSE);
output$peakhist_pic<-renderImage(exprrec, deleteFile = FALSE);
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
path=file.path(logfile$project_folder,"pics","boxprofile_pos")
png(filename = path, bg = "white")
plot.new();plot.window(xlim=c(0,1),ylim=c(0,1));text(0.5,0.5,"nothing selected \n or not available",cex=1)
dev.off()
expr4p<-list(src=path)
output$boxprofile<-renderImage(expr4p, deleteFile = FALSE)
path=file.path(logfile$project_folder,"pics","boxprofile_neg")
png(filename = path, bg = "white")
plot.new();plot.window(xlim=c(0,1),ylim=c(0,1));text(0.5,0.5,"nothing selected \n or not available",cex=1)
dev.off()
expr4n<-list(src=path)
output$boxprofile<-renderImage(expr4n, deleteFile = FALSE)
elements<-unique(as.character(isotopes[1:295,1]))
elements<-elements[order(elements)]
isotopos<-c()
for(k in 1:length(elements)){
if(length(isotopes[isotopes[,1]==elements[k],1])>1){
elem_isos<-isotopes[isotopes[,1]==elements[k],]
elem_isos<-as.character(elem_isos[-1,2])
isotopos<-c(isotopos,elem_isos)
}
}
if(any(ls()=="logfile")){stop("\n illegal logfile detected #1b in server_startup.r!")}
source("server_variables_in.R", local=TRUE)
isolate(init$a<-"TRUE")
cat("Started a new project\n");
}else{
createAlert(session, anchorId="failed_new", alertId = "failed_new_id", title = "Invalid project path",
content = "Project already exists or the specified path is invalid.",
style = "danger", dismiss = TRUE, append = FALSE)
cat("Invalid - project already exists or path invalid \n")
}
}
}
if(any(ls()=="logfile")){stop("\n illegal logfile detected #2 in server_startup.r!")}
})
maincalc2<-reactive({
input$openit
if(isolate(input$openit)){
if(any(ls()=="logfile")){stop("\n illegal logfile detected #3 in server_startup.r!")}
source("server_cleaner.R", local=TRUE);
cat("observed openit")
file_in <- NA
try({
file_in <- as.character(parseFilePaths(getVolumes()(), isolate(input$pro_dir3))[1,4])
})
if(is.na(file_in)){
#cat("\n is NA! \n");cat(file_in);
file_in <- as.character(isolate(input$pro_dir2))
if(grepl("\\", file_in, fixed = TRUE)){
file_in <- gsub("\\", .Platform$file.sep, file_in, fixed = TRUE)
}
}else{
#cat("\n is not NA! \n");cat(file_in);
file_in <- strsplit(file_in, .Platform$file.sep)[[1]]
file_in <- file_in[-length(file_in)]
file_in <- paste0(file_in, collapse = .Platform$file.sep)
}
if(file.exists(file.path(file_in, "logfile.emp"))){
load(file.path(file_in, "logfile.emp"), envir = as.environment(".GlobalEnv"))
logfile$project_folder <<- as.character(file_in);
source("server_updates.R", local=TRUE);
save(logfile,file = file.path(file_in,"logfile.emp"));
output$textit <- renderText(logfile$project_folder);
output$summa_html <- renderText(enviMass::UCBRR9756L(logfile$summary, logfile$Tasks_to_redo));
output$dowhat <- renderText("Opened existing project");
output$IS <- DT::renderDataTable(read.table(file = file.path(logfile$project_folder, "dataframes", "IS.txt"), header = TRUE, sep = "\t", colClasses = "character"));
output$targets <- DT::renderDataTable(read.table(file = file.path(logfile$project_folder, "dataframes", "targets.txt"), header = TRUE, sep = "\t", colClasses = "character"));
measurements <- read.csv(file = file.path(logfile$project_folder, "dataframes", "measurements"), colClasses = "character")
measurements_tab <- measurements
measurements_tab[,1] <- as.numeric(measurements_tab[,1])
measurements_tab <- DT::datatable(
measurements_tab[,c("ID","Name","Type","Mode","Place","Date","Time","include","profiled","tag1",
"tag2","tag3","date_end","time_end","ID_2")],
extensions = c('Buttons','ColReorder','FixedHeader'),
rownames = FALSE,
options = list(
lengthMenu = list(c(25, 50, 100, -1), list('25', '50', '100', 'All')),
fixedHeader = TRUE,
ordering = T,
dom = 'Blfrtip',
buttons = c('excel', 'csv', 'colvis'),#buttons = c('excel', 'pdf', 'print', 'csv'),
scrollX = TRUE,
colReorder = TRUE
)
)
output$measurements <- DT::renderDataTable(measurements_tab);
output$sel_meas_comp_state <- renderText("")
enviMass::DODNL4951M(session)
source("server_variables_in.R", local = TRUE)
if(
(file.exists(file.path(as.character(logfile[[1]]), "results", "int_distrib"))) &
!any(objects(envir = as.environment(".GlobalEnv")) == "no_load")
){
load(file.path(as.character(logfile[[1]]), "results", "int_distrib"), envir = as.environment(".GlobalEnv"))
output$int_box_pos <- renderPlot({
par(mar = c(5.8, 4.5, .9, 8))
enviMass:::MSBYD6677H(
int_distrib,
ion_mode = "positive",
measurements = measurements,
what = "boxplot"
)
},res = 100)
output$int_box_neg <- renderPlot({
par(mar = c(5.8, 4.5, .9, 8))
enviMass:::MSBYD6677H(
int_distrib,
ion_mode = "negative",
measurements = measurements,
what = "boxplot"
)
},res = 100)
output$int_quantiles_pos <- renderPlot({
par(mar = c(4.5, 4.5, .9, 8))
enviMass:::MSBYD6677H(
int_distrib,
ion_mode = "positive",
measurements = measurements,
what = "quantiles_distrib"
)
},res = 100)
output$int_quantiles_neg <- renderPlot({
par(mar = c(4.5, 4.5, .9, 8))
enviMass:::MSBYD6677H(
int_distrib,
ion_mode = "negative",
measurements = measurements,
what = "quantiles_distrib"
)
},res = 100)
output$int_maxmed_pos <- renderPlot({
par(mar = c(4.5, 4.5, .9, 8))
enviMass:::MSBYD6677H(
int_distrib,
ion_mode = "positive",
measurements = measurements,
what = "quantiles_out"
)
},res = 100)
output$int_maxmed_neg <- renderPlot({
par(mar = c(4.5, 4.5, .9, 8))
enviMass:::MSBYD6677H(
int_distrib,
ion_mode = "negative",
measurements = measurements,
what = "quantiles_out"
)
},res = 100)
}else{
output$int_box_pos <- renderPlot({
plot.new(); plot.window(xlim = c(0, 10), ylim =c(0, 10)); text(5, 5, "Data not available")
})
output$int_box_neg <- renderPlot({
plot.new(); plot.window(xlim = c(0, 10), ylim =c(0, 10)); text(5, 5, "Data not available")
})
output$int_quantiles_pos <- renderPlot({
plot.new(); plot.window(xlim = c(0, 10), ylim =c(0, 10)); text(5, 5, "Data not available")
})
output$int_quantiles_neg <- renderPlot({
plot.new(); plot.window(xlim = c(0, 10), ylim =c(0, 10)); text(5, 5, "Data not available")
})
output$int_maxmed_pos <- renderPlot({
plot.new(); plot.window(xlim = c(0, 10), ylim =c(0, 10)); text(5, 5, "Data not available")
})
output$int_maxmed_neg <- renderPlot({
plot.new(); plot.window(xlim = c(0, 10), ylim =c(0, 10)); text(5, 5, "Data not available")
})
}
if(file.exists(file.path(logfile$project_folder, "pics", "int_distr_pos"))){
expr3p <- list(src = file.path(logfile$project_folder, "pics", "int_distr_pos"))
output$pic_int_distr_pos <- renderImage(expr3p, deleteFile = FALSE)
}
if(file.exists(file.path(logfile$project_folder, "pics", "int_distr_neg"))){
expr3n <- list(src = file.path(logfile$project_folder, "pics", "int_distr_neg"))
output$pic_int_distr_neg <- renderImage(expr3n, deleteFile = FALSE)
}
path=file.path(logfile$project_folder,"pics","recal_none")
png(filename = path, bg = "white")
plot.new(); plot.window(xlim = c(0, 1), ylim = c(0, 1)); text(0.5,0.5,"nothing selected \n or not available",cex=1)
dev.off()
exprrec <- list(src = path)
output$recal_pic <- renderImage(exprrec, deleteFile = FALSE);
output$peakhist_pic <- renderImage(exprrec, deleteFile = FALSE);
if(file.exists(
file.path(as.character(logfile[[1]]), "results", "int_norm_ISTD_pos")) &
!any(objects(envir = as.environment(".GlobalEnv")) == "no_load")
){
load(file.path(as.character(logfile[[1]]), "results", "int_norm_ISTD_pos"), envir = as.environment(".GlobalEnv"))
output$int_norm_ISTD_pos_median <- renderPlot({
par(mar = c(.2, 4.5, .9, 8))
enviMass:::TICFG4083H(
int_norm_ISTD = int_norm_ISTD_pos,
logfile = logfile,
what = "normalization",
ion_mode = "positive"
)
},res = 100)
output$int_norm_ISTD_pos_counts <- renderPlot({
par(mar = c(4.5, 4.5, .9, 8))
enviMass:::TICFG4083H(
int_norm_ISTD = int_norm_ISTD_pos,
logfile = logfile,
what = "counts",
ion_mode = "positive"
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
if(
file.exists(file.path(as.character(logfile[[1]]), "results", "int_norm_ISTD_neg")) &
!any(objects(envir = as.environment(".GlobalEnv")) == "no_load")
){
load(file.path(as.character(logfile[[1]]), "results", "int_norm_ISTD_neg"), envir = as.environment(".GlobalEnv"))
output$int_norm_ISTD_neg_median <- renderPlot({
par(mar = c(.2, 4.5, .9, 8))
enviMass:::TICFG4083H(
int_norm_ISTD = int_norm_ISTD_neg,
logfile = logfile,
what = "normalization",
ion_mode = "negative"
)
},res = 100)
output$int_norm_ISTD_neg_counts <- renderPlot({
par(mar = c(4.5, 4.5, .9, 8))
enviMass:::TICFG4083H(
int_norm_ISTD = int_norm_ISTD_neg,
logfile = logfile,
what = "counts",
ion_mode = "negative"
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
if(file.exists(file.path(logfile$project_folder,"pics","boxprofile_pos"))){
if(isolate(input$Ion_mode)=="positive"){
expr4p<-list(src=file.path(logfile$project_folder,"pics","boxprofile_pos"))
output$boxprofile<-renderImage(expr4p, deleteFile = FALSE)
}
}
if(file.exists(file.path(logfile$project_folder,"pics","boxprofile_neg"))){
if(isolate(input$Ion_mode)=="negative"){
expr4n<-list(src=file.path(logfile$project_folder,"pics","boxprofile_neg"))
output$boxprofile<-renderImage(expr4n, deleteFile = FALSE)
}
}
enviMass:::VGMYK9509S(
logfile,
Ion_mode_profiles = isolate(input$Ion_mode)
)
if(file.exists(file.path(logfile$project_folder,"pics","profilehisto.png"))){
expr6<-list(src=file.path(logfile$project_folder,"pics","profilehisto.png"))
output$profilehisto<-renderImage(expr6, deleteFile = FALSE)
}
elements<-unique(as.character(isotopes[1:295,1]))
elements<-elements[order(elements)]
isotopos<-c()
for(k in 1:length(elements)){
if(length(isotopes[isotopes[,1]==elements[k],1])>1){
elem_isos<-isotopes[isotopes[,1]==elements[k],]
elem_isos<-as.character(elem_isos[-1,2])
isotopos<-c(isotopos,elem_isos)
}
}
if(length(logfile$comparisons)){
at_comparisons <- which(names(logfile$comparisons) != "")
if(length(at_comparisons)){
updateSelectInput(session, "load_comparison", choices = c("None", names(logfile$comparisons)[at_comparisons]), selected = "None")
}
}
cat(objects())
if(any(ls()=="logfile")){stop("\n illegal logfile detected #3b in server_startup.r!")}
if(isolate(init$a=="FALSE")){
isolate(init$a<-"TRUE")
}else{
isolate(init$b<-(init$b+1));cat(" - ")
}
cat("\nProject opened\n")
return("Project available\n")
}else{
createAlert(session, anchorId="failed_open", alertId = "failed_open_id", title = "Invalid project path",
content = "Project with specified path does not exist!",
style = "danger", dismiss = TRUE, append = FALSE)
cat("Invalid - project already exists or path invalid \n")
cat("Invalid project!\n")
return("Project invalid\n")
}
}
if(any(ls()=="logfile")){stop("\n illegal logfile detected #3 in server_startup.r!")}
})
output$had_opened<-renderText(paste(maincalc2()))
