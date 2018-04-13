ranges_peaks_mz_RT <- reactiveValues(x = NULL, y = NULL, xchroma=FALSE, ychroma=FALSE)
refresh_plot<-reactiveValues()
refresh_plot$a<-1
refresh_plot$b<-0
refresh_plot$c<-0
observe({
input$sel_meas
input$blind_boxplot_log
if(isolate(init$a)=="TRUE"){
if(!is.na(isolate(input$sel_meas))){
if(isolate(input$sel_meas)!=0){
measurements<-read.csv(file = file.path(logfile[[1]],"dataframes","measurements"), colClasses = "character");
if(any(measurements[,"ID"] == as.character(isolate(input$sel_meas)))){
output$file_proc_name <- renderText(paste("File name: ",measurements[measurements[,"ID"] == as.character(isolate(input$sel_meas)),"Name"],sep=""))
output$file_proc_type <- renderText(paste("File type: ",measurements[measurements[,"ID"] == as.character(isolate(input$sel_meas)),"Type"],sep=""))
output$file_proc_mode <- renderText(paste("Ionization mode: ",measurements[measurements[,"ID"] == as.character(isolate(input$sel_meas)),"Mode"],sep=""))
if(	file.exists(file.path(logfile$project_folder,"peaklist",as.character(isolate(input$sel_meas)))) &
(measurements[measurements[,"ID"] == as.character(isolate(input$sel_meas)),"include"] == "TRUE")
){
load(file = file.path(logfile$project_folder,"peaklist",as.character(isolate(input$sel_meas))),envir = as.environment(".GlobalEnv"),verbose=FALSE);
output$file_peak_number <- renderText(as.character(length(peaklist[,1])));
blind_aff <- round(
(sum(peaklist[,colnames(peaklist) == "keep_2"] < Inf))/length(peaklist[,1])*100
,digits=3)
output$file_blind_aff <- renderText(as.character(blind_aff));
output$file_blind_aff2 <- renderText(as.character(blind_aff));
blind_rem <- round(
(sum(peaklist[,colnames(peaklist) == "keep_2"] < as.numeric(logfile$parameters$blind_threshold))) / length(peaklist[,1])*100
,digits=3)
output$file_blind_rem <- renderText(as.character(blind_rem));
output$file_blind_rem2 <- renderText(as.character(blind_rem));
repl_rem<-round(
(sum(peaklist[,colnames(peaklist) == "keep"] == 0)) / length(peaklist[,1])*100
,digits=3)
output$file_repl_rem<-renderText(as.character(repl_rem));
isolate(refresh_plot$c<-(refresh_plot$c+1))
}else{
cat("\n no processed peaklist found for the selected file.")
isolate(refresh_plot$b <- 0)
isolate(refresh_plot$c <- 0)
output$file_peak_number <- renderText("none");
output$file_blind_aff <- renderText("none");
output$file_blind_aff2 <- renderText("none");
output$file_blind_rem <- renderText("none");
output$file_blind_rem2 <- renderText("none");
output$file_repl_rem <- renderText("none");
}
if(
(logfile$workflow[names(logfile$workflow)=="blind"]=="yes") &
(measurements[measurements[,"ID"]==as.character(isolate(input$sel_meas)),"Type"]!="blank") &
(file.exists(file.path(logfile$project_folder,"peaklist",as.character(isolate(input$sel_meas))))) &
(measurements[measurements[,"ID"]==as.character(isolate(input$sel_meas)),"include"]=="TRUE")
){
output$showblank<-renderText("Blank/blind peak tagging (subtraction) results:")
load(file=file.path(logfile[[1]],"peaklist",as.character(isolate(input$sel_meas))),verbose=FALSE);
if(any(peaklist[,"keep_2"]<Inf)){
if(isolate(input$blind_boxplot_log)=="TRUE"){
output$blind_boxplot <- renderPlot({
par(mar=c(4.2,4.2,0.8,0.8))
boxplot(log10(peaklist[peaklist[,"keep_2"]<Inf,"keep_2"]),horizontal=TRUE)
title(xlab=expression(log[10]*paste(" intensity ratio sample vs. blind",sep=" ")))
abline(v=0,col="red")
})
}else{
output$blind_boxplot <- renderPlot({
par(mar=c(4.2,4.2,0.8,0.8))
boxplot(peaklist[peaklist[,"keep_2"]<Inf,"keep_2"],horizontal=TRUE)
title(xlab="Intensity ratio sample vs. blind")
abline(v=1,col="red")
})
}
}else{
output$blind_boxplot <- renderPlot({
par(mar=c(4.2,4.2,0.8,0.8))
plot.new()
plot.window(xlim=c(0,1),ylim=c(0,1))
text(x=.5,y=.5,labels="No peaks matched to any blind/blank peaks",pos=1,cex=1.5)
})
}
rm(peaklist)
}else{
output$showblank<-renderText("No blank/blind subtraction results available.")
}
pics<-list.files(file.path(logfile[[1]],"pics"))
if(
any(pics==paste("recal_",as.character(isolate(input$sel_meas)),sep="")) &
(logfile$workflow[names(logfile$workflow)=="recal"]=="yes") &
(measurements[measurements[,"ID"]==as.character(isolate(input$sel_meas)),"include"]=="TRUE")
){
output$showrecal<-renderText("Mass recalibration results:")
expr1<-list(src=file.path(logfile[[1]],"pics",paste("recal_",as.character(isolate(input$sel_meas)),sep="")))
output$recal_pic<-renderImage(expr1, deleteFile = FALSE)
cat("\n Found recal_pic")
}else{
output$showrecal<-renderText("No mass recalibration results available.")
cat("\n Not found recal_pic")
}
if(
any(pics==paste("peakhist_",as.character(isolate(input$sel_meas)),sep="")) &
(measurements[measurements[,"ID"]==as.character(isolate(input$sel_meas)),"include"]=="TRUE")
){
output$showintensitydistrib<-renderText("Centroid & peak intensity distribution:")
expr_peakhist<-list(src=file.path(logfile[[1]],"pics",paste("peakhist_",as.character(isolate(input$sel_meas)),sep="")))
output$peakhist_pic<-renderImage(expr_peakhist, deleteFile = FALSE)
cat("\n Found peakhist_ pic")
}else{
output$showintensitydistrib<-renderText("No centroid/peak intensity distribution available.")
cat("\n Not found peakhist_ pic")
}
if(
file.exists( file.path(logfile[[1]],"results","LOD",paste("plot_LOD_",as.character(isolate(input$sel_meas)),".png",sep="") ) ) &
(measurements[measurements[,"ID"]==as.character(isolate(input$sel_meas)),"include"]=="TRUE")
){
output$showLOD<-renderText("LOD interpolation results:")
expr_LOD <- list( src=file.path(logfile[[1]],"results","LOD",paste("plot_LOD_",as.character(isolate(input$sel_meas)),".png",sep="")) )
output$LOD_pic <- renderImage(expr_LOD, deleteFile = FALSE)
cat("\n LOD pic file found")
}else{
output$showLOD <- renderText("No LOD interpolation available.")
cat("\n LOD pic file not found")
}
output$dowhat<-renderText("Processing per file viewed.");
}else{
output$showblank <- renderText("FALSE")
output$showrecal <- renderText("FALSE")
output$showintensitydistrib <- renderText("FALSE")
output$showLOD <- renderText("FALSE")
output$dowhat <- renderText("Invalid ID chosen to view processing results.");
isolate(refresh_plot$b <- 0)
isolate(refresh_plot$c <- 0)
}
}
}
}
})
observe({
refresh_plot$c
input$peak_chromat_refresh
if(isolate(init$a)=="TRUE" & isolate(refresh_plot$c>0)){
use_keep_2<-peaklist[,"keep_2"]
use_keep_2[use_keep_2==Inf]<-0
output$exp_peaklist <- DT::renderDataTable(
DT::datatable(
dat<-data.frame(
I(as.character(peaklist[,"peak_ID"])),
round(peaklist[,"m/z"],digits=5),
round(log10(peaklist[,"var_m/z"]),digits=4),
round(peaklist[,"m/z_corr"],digits=5),
round(log10(peaklist[,"max_int"]),digits=5),
round(log10(peaklist[,"sum_int"]),digits=5),
round(log10(peaklist[,"int_corr"]),digits=5),
round(peaklist[,"RT"],digits=1),
round((peaklist[,"RT"]/60),digits=2),
round(peaklist[,"minRT"],digits=1),
round(peaklist[,"maxRT"],digits=1),
round(peaklist[,"RT_corr"],digits=1),
I(as.character(peaklist[,"keep"]==1)),
round(use_keep_2,digits=2)
),
filter = list(position = 'top', clear = FALSE, plain = TRUE),
colnames=c(
"Peak ID","m/z",
"log10 var(m/z)","m/z_recal",
"log10 max int.","log10 sum int.","log10 norm int.",
"RT [s]","RT [min]","min_RT [s]","max_RT [s]","RT_align [s]",
"in replicates?","int. ratio sample/blind"
),
rownames = FALSE,
selection = list(mode = 'multiple', target = 'row'),
extensions = c('Buttons','FixedHeader','ColReorder'),
options = list(
lengthMenu = list(c(200, 500, 1000, -1), list('200', '500', '1000', 'All')),
fixedHeader = FALSE,
ordering = TRUE,
dom = 'Blfrtip',
buttons = c('excel', 'csv','colvis'),#buttons = c('excel', 'pdf', 'print', 'csv'),
scrollX = TRUE,
scrollY = "800px",
colReorder = TRUE
)
),server = TRUE
)
}else{
output$exp_peaklist <- DT::renderDataTable(
DT::datatable(cbind("no data available","")),
colnames=c("","")
)
}
})
observe({
s5<-input$exp_peaklist_rows_selected
if(isolate(init$a)=="TRUE"){
if(length(s5)){
if(logfile$parameters$verbose){ cat("\n Selected rows: ");print(s5)}
these_peaks<<-peaklist[s5,"peak_ID"];print(these_peaks)
if(any(objects(envir=as.environment(".GlobalEnv"))=="MSlist")){
if(any(names(MSlist)=="File_ID")){
if(MSlist[["File_ID"]]!=as.character(isolate(input$sel_meas))){
load(file.path(logfile[[1]],"MSlist",as.character(isolate(input$sel_meas))),envir=as.environment(".GlobalEnv"))
}
}else{
load(file.path(logfile[[1]],"MSlist",as.character(isolate(input$sel_meas))),envir=as.environment(".GlobalEnv"))
MSlist[["File_ID"]]<-as.character(isolate(input$sel_meas))
}
}else{
load(file.path(logfile[[1]],"MSlist",as.character(isolate(input$sel_meas))),envir=as.environment(".GlobalEnv"))
}
isolate(refresh_plot$b <- (refresh_plot$b+1))
}else{
if(logfile$parameters$verbose) cat("\n Selected nothing: ");print(s5)
}
}
})
observe({
input$peak_chromat_time
if(isolate(init$a)=="TRUE" & isolate(ranges_peaks_mz_RT$xchroma[1]!=FALSE)){
if(isolate(input$peak_chromat_time)=="minutes"){
isolate(ranges_peaks_mz_RT$xchroma<-(ranges_peaks_mz_RT$xchroma/60))
}
if(isolate(input$peak_chromat_time)=="seconds"){
isolate(ranges_peaks_mz_RT$xchroma<-(ranges_peaks_mz_RT$xchroma*60))
}
}
})
observe({
input$peak_chromat_norm
if(isolate(init$a)=="TRUE" & isolate(ranges_peaks_mz_RT$ychroma[1]!=FALSE)){
isolate(ranges_peaks_mz_RT$ychroma<-FALSE)
}
})
observe({
refresh_plot$b
input$peak_chromat_norm
input$peak_chromat_time
input$peak_chromat_type
if(isolate(init$a)=="TRUE" & isolate(refresh_plot$b>0)){
output$peak_chromat <- renderPlot({
par(mar=c(5,4,1,.8))
enviMass:::JJGET6467P(
MSlist,
peakIDs=these_peaks,
RTlim=isolate(ranges_peaks_mz_RT$xchroma),
Intlim=isolate(ranges_peaks_mz_RT$ychroma),
normalize=as.logical(isolate(input$peak_chromat_norm)),
n_col=dim(peaklist)[1],
set_RT=isolate(input$peak_chromat_time),
chromat_full=input$peak_chromat_type
);
},res=100)
}
})
observeEvent(input$peak_chromat_dblclick, {
if(isolate(init$a)=="TRUE"){
brush <- isolate(input$peak_chromat_brush)
if (!is.null(brush)) {
if(logfile$parameters$verbose) cat("\n Zoom in_1e")
isolate(ranges_peaks_mz_RT$xchroma <- c(brush$xmin, brush$xmax))
isolate(ranges_peaks_mz_RT$ychroma <- c(brush$ymin, brush$ymax))
} else {
if(logfile$parameters$verbose) cat("\n Zoom out full_1e")
isolate(ranges_peaks_mz_RT$xchroma <- FALSE)
isolate(ranges_peaks_mz_RT$ychroma <- FALSE)
}
isolate(refresh_plot$b<-(refresh_plot$b+1))
}
})
observeEvent(input$peak_chromat_click, {
if(isolate(init$a)=="TRUE"){
if(logfile$parameters$verbose) cat("\n Zoom out part_1_ae")
brush <- isolate(input$peak_chromat_brush)
if (is.null(brush)) {
if(logfile$parameters$verbose) cat("\n Zoom out part_1_be")
if(isolate(ranges_peaks_mz_RT$xchroma[1])!=FALSE){
old_range_dmass<-abs(isolate(ranges_peaks_mz_RT$xchroma[2]-ranges_peaks_mz_RT$xchroma[1]))
isolate(ranges_peaks_mz_RT$xchroma[1]<-ranges_peaks_mz_RT$xchroma[1]-.3*old_range_dmass)
isolate(ranges_peaks_mz_RT$xchroma[2]<-ranges_peaks_mz_RT$xchroma[2]+.3*old_range_dmass)
}
if(isolate(ranges_peaks_mz_RT$ychroma[1])!=FALSE){
old_range_dmass<-abs(isolate(ranges_peaks_mz_RT$ychroma[2]-ranges_peaks_mz_RT$ychroma[1]))
isolate(ranges_peaks_mz_RT$ychroma[1]<-ranges_peaks_mz_RT$ychroma[1]-.3*old_range_dmass)
isolate(if(ranges_peaks_mz_RT$ychroma[1]<0){ranges_peaks_mz_RT$ychroma[1]<-0})
isolate(ranges_peaks_mz_RT$ychroma[2]<-ranges_peaks_mz_RT$ychroma[2]+.3*old_range_dmass)
}
isolate(refresh_plot$b<-(refresh_plot$b+1))
}else{
if(logfile$parameters$verbose) cat("\n Doing hover_e - nothing")
}
}
})
observe({
refresh_plot$a
input$sel_meas_ID
input$peaks_mz_RT_use_peaks
input$peaks_mz_RT_use_raw
input$peaks_mz_RT_use_IDs
input$peaks_mz_RT_use_window
input$peaks_mz_RT_use_window_mass
input$peaks_mz_RT_use_window_RT
input$peaks_mz_RT_use_bar
input$peaks_mz_RT_use_bar_value
input$peaks_mz_RT_use_window_RT_tol
input$plot_filter_intensity
input$plot_filter_blind
input$plot_filter_replicates
input$EIC_mz_RT_use_IDs
if(isolate(init$a)=="TRUE"){
if(!is.na(isolate(input$sel_meas_ID))){
if(!any(objects(envir=as.environment(".GlobalEnv"))=="atit")){
assign("atit",0,envir=as.environment(".GlobalEnv"))
}
if(
file.exists(file.path(logfile[[1]],"MSlist",as.character(isolate(input$sel_meas_ID)))) &
file.exists(file.path(logfile[[1]],"peaklist",as.character(isolate(input$sel_meas_ID))))
){
if(verbose){cat("\n Plotting")}
if((isolate(input$sel_meas_ID)!=atit)){
if(any(objects(envir=as.environment(".GlobalEnv"))=="MSlist")){rm(MSlist,envir=as.environment(".GlobalEnv"))}
if(any(objects()=="MSlist")){rm(MSlist)}
load(file.path(logfile[[1]],"MSlist",as.character(isolate(input$sel_meas_ID))), envir=as.environment(".GlobalEnv"))
load(file.path(logfile[[1]],"peaklist",as.character(isolate(input$sel_meas_ID))), envir=as.environment(".GlobalEnv"))
if(verbose){cat("\nafter:\n "); print(gc());cat("\n MSlist file loaded")}
assign("atit",isolate(input$sel_meas_ID),envir=as.environment(".GlobalEnv"))
max_int<-max(max(MSlist[["Peaklist"]][,"max_int"]),max(MSlist[["Scans"]][[2]][,"intensity"]))
min_int<-min(min(MSlist[["Peaklist"]][,"max_int"]),min(MSlist[["Scans"]][[2]][,"intensity"]))
updateSliderInput(session,"plot_filter_intensity", min=round(log10(min_int),digits=1), max=round(log10(max_int),digits=1), value=c(log10(min_int),log10(max_int)))
measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
if(any(measurements[,"ID"]==as.character(isolate(input$sel_meas_ID)))){
output$file_viewer_name<-renderText(paste("File name: ",measurements[measurements[,"ID"]==as.character(isolate(input$sel_meas_ID)),"Name"],sep=""))
output$file_viewer_type<-renderText(paste("File type: ",measurements[measurements[,"ID"]==as.character(isolate(input$sel_meas_ID)),"Type"],sep=""))
output$file_viewer_mode<-renderText(paste("Ionization mode: ",measurements[measurements[,"ID"]==as.character(isolate(input$sel_meas_ID)),"Mode"],sep=""))
}else{
output$file_viewer_name<-renderText("Invalid file ID")
output$file_viewer_type<-renderText("")
output$file_viewer_mode<-renderText("")
}
}
if(
(exists("MSlist")) &
(exists("peaklist"))
){
if(!is.null(isolate(ranges_peaks_mz_RT$x))){
x_lim<-isolate(ranges_peaks_mz_RT$x)
}else{
x_lim<-c(min(MSlist[["Scans"]][[2]][,"m/z"]),max(MSlist[["Scans"]][[2]][,"m/z"]))
}
if(!is.null(isolate(ranges_peaks_mz_RT$y))){
y_lim<-isolate(ranges_peaks_mz_RT$y)
}else{
y_lim<-c(min(MSlist[["Scans"]][[2]][,"RT"]),max(MSlist[["Scans"]][[2]][,"RT"]))
}
use_these <- which(
(peaklist[,"m/z"]>=x_lim[1]) &
(peaklist[,"m/z"]<=x_lim[2]) &
(peaklist[,"RT"]>=y_lim[1]) &
(peaklist[,"RT"]<=y_lim[2]) &
(peaklist[,"max_int"]>=isolate(10^input$plot_filter_intensity[1])) &
(peaklist[,"max_int"]<=isolate(10^input$plot_filter_intensity[2]))
)
if((length(use_these)>0) & isolate(input$plot_filter_blind)){
use_these <- use_these[peaklist[use_these,"keep_2"]>=as.numeric(logfile$parameters$blind_threshold)]
}
if((length(use_these)>0) & isolate(input$plot_filter_replicates)){
use_these <- use_these[peaklist[use_these,"keep"]==1]
}
if(isolate(input$peaks_mz_RT_use_raw)){
if(is.null(isolate(ranges_peaks_mz_RT$x))){
x_min<-min(MSlist[["Scans"]][[2]][,"m/z"])
x_max<-max(MSlist[["Scans"]][[2]][,"m/z"])
}else{
x_min<-min(MSlist[["Scans"]][[2]][
MSlist[["Scans"]][[2]][,"m/z"]>=isolate(ranges_peaks_mz_RT$x[1])
,"m/z"])
x_max<-max(MSlist[["Scans"]][[2]][
MSlist[["Scans"]][[2]][,"m/z"]<=isolate(ranges_peaks_mz_RT$x[2])
,"m/z"])
}
if(is.null(isolate(ranges_peaks_mz_RT$y))){
y_min <- min(MSlist[["Scans"]][[2]][,"RT"])
y_max <- max(MSlist[["Scans"]][[2]][,"RT"])
}else{
y_min <- min(MSlist[["Scans"]][[2]][
MSlist[["Scans"]][[2]][,"RT"]>=isolate(ranges_peaks_mz_RT$y[1])
,"RT"])
y_max <- max(MSlist[["Scans"]][[2]][
MSlist[["Scans"]][[2]][,"RT"]<=isolate(ranges_peaks_mz_RT$y[2])
,"RT"])
}
those <- which(
(MSlist[["Scans"]][[2]][,"m/z"] >= x_min) &
(MSlist[["Scans"]][[2]][,"m/z"] <= x_max) &
(MSlist[["Scans"]][[2]][,"RT"] >= y_min) &
(MSlist[["Scans"]][[2]][,"RT"] <= y_max) &
(MSlist[["Scans"]][[2]][,"intensity"] >= isolate(10^input$plot_filter_intensity[1])) &
(MSlist[["Scans"]][[2]][,"intensity"] <= isolate(10^input$plot_filter_intensity[2]))
)
if(length(those) <= 1E6){
colorit <- rep("gray",length(those))
colorit[MSlist[["Scans"]][[2]][those,"peakID"]!=0] <- "red"
}
}else{
those <- c()
}
if(
(isolate(input$peaks_mz_RT_use_peaks)) &
(length(use_these))
){
z_max <- max(peaklist[use_these, "max_int"])
}else{
if(
(isolate(input$peaks_mz_RT_use_raw)) &
(length(those))
){
z_max <- max(MSlist[["Scans"]][[2]][those, "intensity"])
}else{
z_max <- max(peaklist[,"max_int"])
}
}
z_lim <- c(0, z_max)
output$plot_peaks_mz_RT <- renderPlot({
par(mar=c(4, 4, 3, .1))
plot.new()
plot.window(xlim=x_lim,ylim=y_lim)
title(
xlab="m/z [Th]", ylab="RT [s]",
main="Draw rectangles and double-click into them to zoom in, double-click again to zoom fully out. Bottom plots adapt accordingly.",cex.main=.75
)
box();axis(1);axis(2);
if(isolate(input$peaks_mz_RT_use_raw)){
if(length(those)<=1E5 & length(those)>0){
if(!isolate(input$EIC_mz_RT_use_IDs)){
points(
MSlist[["Scans"]][[2]][those,"m/z"],
MSlist[["Scans"]][[2]][those,"RT"],
pch=19,cex=.5,col=colorit
)
}else{
text(
MSlist[["Scans"]][[2]][those,"m/z"],
MSlist[["Scans"]][[2]][those,"RT"],
labels = as.character(MSlist[["Scans"]][[2]][those,"clustID"]),
pos = NULL, col	= "darkblue", cex=.6
)
}
}else{
smoothScatter(
x=MSlist[["Scans"]][[2]][those,"m/z"],
y=MSlist[["Scans"]][[2]][those,"RT"],
colramp = colorRampPalette(c("white", "red")),
nbin = 200, add = TRUE
)
}
}
if(isolate(input$peaks_mz_RT_use_peaks)){
if(length(use_these)>0){
if(!isolate(input$peaks_mz_RT_use_IDs)){
colorit<-"black"
}else{
colorit<-"darkgrey"
}
points(
peaklist[use_these,"m/z"],
peaklist[use_these,"RT"],
pch=21,cex=.8,col=colorit
)
}
}
if(isolate(input$peaks_mz_RT_use_IDs)){
if(length(use_these)>0){
text(
peaklist[use_these,"m/z"],
peaklist[use_these,"RT"],
labels = as.character(peaklist[use_these,"peak_ID"]),
pos = NULL, col	= "darkred", cex=.6
)
}
}
if(isolate(input$peaks_mz_RT_use_window)){
del_mz<-((x_lim[2]-x_lim[1])/15)
at_mz<-isolate(as.numeric(input$peaks_mz_RT_use_window_mass))
at_RT<-isolate(as.numeric(input$peaks_mz_RT_use_window_RT))
rect(
xleft=(at_mz-del_mz),
ybottom=(at_RT-isolate(as.numeric(input$peaks_mz_RT_use_window_RT_tol))),
xright=(at_mz+del_mz),
ytop=(at_RT+isolate(as.numeric(input$peaks_mz_RT_use_window_RT_tol))),
col=NULL, border="blue",lwd=2)
if(isolate(input$peaks_mz_RT_use_bar)){
del_ppm<-(at_mz*isolate(as.numeric(input$peaks_mz_RT_use_bar_value))/1E6)
lines(
x=c((at_mz-del_ppm),(at_mz+del_ppm)),
y=c(at_RT,at_RT),
col="blue",lwd=2)
}
}
}, res = 100, execOnResize=TRUE)
output$plot_peaks_mz_int <- renderPlot({
par(mar=c(4, 4, .8, .2))
plot.new()
plot.window(xlim=x_lim,ylim=c(0,z_lim[2]))
title(xlab="m/z [Th]", ylab="Intensity")
box();axis(1);axis(2);
if(isolate(input$peaks_mz_RT_use_raw)){
if(length(those)<=1E5 & length(those)>0){
points(
MSlist[["Scans"]][[2]][those,"m/z"],
MSlist[["Scans"]][[2]][those,"intensity"],
type="h",pch=19,cex=.5,col=colorit
)
}
}
if(isolate(input$peaks_mz_RT_use_peaks)){
if(length(use_these)>0){
if(!isolate(input$peaks_mz_RT_use_IDs)){
colorit<-"black"
}else{
colorit<-"darkgrey"
}
points(
peaklist[use_these,"m/z"],
peaklist[use_these,"max_int"],
pch=21,cex=.8,col=colorit,type="h"
)
}
}
if(isolate(input$peaks_mz_RT_use_window)){
if(isolate(input$peaks_mz_RT_use_bar)){
at_mz<-isolate(as.numeric(input$peaks_mz_RT_use_window_mass))
del_ppm<-(at_mz*isolate(as.numeric(input$peaks_mz_RT_use_bar_value))/1E6)
lines(
x=c((at_mz-del_ppm),(at_mz+del_ppm)),
y=c(0.5*z_lim[2],.5*z_lim[2]),
col="blue",lwd=2)
}
}
}, res = 100, execOnResize=TRUE)
output$plot_peaks_RT_int <- renderPlot({
par(mar=c(4, 4, .8, .2))
plot.new()
plot.window(xlim=y_lim,ylim=c(0,z_lim[2]))
title(xlab="RT", ylab="Intensity")
box();axis(1);axis(2);
if(isolate(input$peaks_mz_RT_use_raw)){
if(length(those)<=1E5 & length(those)>0){
points(
MSlist[["Scans"]][[2]][those,"RT"],
MSlist[["Scans"]][[2]][those,"intensity"],
type="h",pch=19,cex=.5,col=colorit
)
}
}
if(isolate(input$peaks_mz_RT_use_peaks)){
if(length(use_these)>0){
if(!isolate(input$peaks_mz_RT_use_IDs)){
colorit<-"black"
}else{
colorit<-"darkgrey"
}
points(
peaklist[use_these,"RT"],
peaklist[use_these,"max_int"],
pch=21,cex=.8,col=colorit,type="h"
)
}
}
if(isolate(input$peaks_mz_RT_use_window)){
at_RT<-isolate(as.numeric(input$peaks_mz_RT_use_window_RT))
lines(
x=c(
(at_RT-isolate(as.numeric(input$peaks_mz_RT_use_window_RT_tol))),
(at_RT+isolate(as.numeric(input$peaks_mz_RT_use_window_RT_tol)))
),
y=c(.5*z_lim[2],.5*z_lim[2]),
col="blue",lwd=2)
}
}, res = 100, execOnResize=TRUE)
if(
( isolate(input$peaks_mz_RT_use_peaks) & !isolate(input$peaks_mz_RT_use_raw) & (length(use_these)>0) ) ||
( isolate(input$peaks_mz_RT_use_peaks) & (length(those)==0) & (length(use_these)>0) )
){
if(length(use_these)<=1E5){
if(verbose){cat("\n Plotting only peaks")}
sub_peaks<-as.data.frame(peaklist[use_these,c("m/z","RT","max_int"),drop=FALSE])
names(sub_peaks)<-c("m_z","RT","Intensity")
sub_peaks[,"Intensity"]<-(sub_peaks[,"Intensity"]/2)
output$plot_peaks_3D <- renderPlotly({
p <- plotly::plot_ly(type="scatter3d",mode="markers",showlegend=TRUE)%>%
plotly::add_trace(p,
x = ~m_z, y = ~RT, z = ~Intensity,
data = sub_peaks,
color=I("black"),
size = I(1),
name = "",
error_z=list(
color="black",
thickness=0,
symmetric = TRUE,
type = "data" ,
array = sub_peaks$Intensity
)
)
})
}else{
output$plot_peaks_3D <- renderPlotly({plotly::plot_ly(type="scatter3d",mode="markers",showlegend=TRUE)})
}
}
if(
( !isolate(input$peaks_mz_RT_use_peaks) & isolate(input$peaks_mz_RT_use_raw) & (length(those)>0) ) ||
( isolate(input$peaks_mz_RT_use_raw) & (length(use_these)==0) & (length(those)>0) )
){
if(length(those)<=1E5){
if(verbose){cat("\n Plotting only raw data")}
sub_MSlist<-as.data.frame(MSlist[["Scans"]][[2]][those,c("m/z","RT","intensity","peakID"),drop=FALSE])
names(sub_MSlist)<-c("m_z","RT","Intensity","peakID")
sub_MSlist[,"Intensity"]<-(sub_MSlist[,"Intensity"]/2)
if(any(sub_MSlist[,"peakID"]!=0)){
output$plot_peaks_3D <- renderPlotly({
p <- plotly::plot_ly(type="scatter3d",mode="markers",showlegend=TRUE)%>%
plotly::add_trace(p,
x = ~m_z, y = ~RT, z = ~Intensity,
data = sub_MSlist[sub_MSlist[,"peakID"]==0,],
color=I("gray"),
size = I(1),
name = "",
error_z=list(
color="gray",
thickness=0,
symmetric = TRUE,
type = "data" ,
array = sub_MSlist[sub_MSlist[,"peakID"]==0,]$Intensity
)
)%>%
plotly::add_trace(p,
x = ~m_z, y = ~RT, z = ~Intensity,
data = sub_MSlist[sub_MSlist[,"peakID"]!=0,],
color=I("red"),
size = I(1),
name = "",
error_z=list(
color="red",
thickness=0,
symmetric = TRUE,
type = "data" ,
array = sub_MSlist[sub_MSlist[,"peakID"]!=0,]$Intensity
)
)
})
}else{
output$plot_peaks_3D <- renderPlotly({
p <- plotly::plot_ly(type="scatter3d",mode="markers",showlegend=TRUE)%>%
plotly::add_trace(p,
x = ~m_z, y = ~RT, z = ~Intensity,
data = sub_MSlist[sub_MSlist[,"peakID"]==0,],
color=I("gray"),
size = I(1),
name = "",
error_z=list(
color="gray",
thickness=0,
symmetric = TRUE,
type = "data" ,
array = sub_MSlist[sub_MSlist[,"peakID"]==0,]$Intensity
)
)
})
}
}else{
output$plot_peaks_3D <- renderPlotly({plotly::plot_ly(type="scatter3d",mode="markers",showlegend=TRUE)})
}
}
if( isolate(input$peaks_mz_RT_use_peaks) & isolate(input$peaks_mz_RT_use_raw) & (length(those)>0) & (length(use_these)>0) ){
if(length(those)<=1E5){
if(verbose){cat("\n Plotting peaks & raw data")}
sub_MSlist<-as.data.frame(MSlist[["Scans"]][[2]][those,c("m/z","RT","intensity","peakID"),drop=FALSE])
names(sub_MSlist)<-c("m_z","RT","Intensity","peakID")
sub_MSlist[,"Intensity"]<-(sub_MSlist[,"Intensity"]/2)
sub_peaks<-as.data.frame(peaklist[use_these,c("m/z","RT","max_int"),drop=FALSE])
names(sub_peaks)<-c("m_z","RT","Intensity")
sub_peaks[,"Intensity"]<-(sub_peaks[,"Intensity"]/2)
output$plot_peaks_3D <- renderPlotly({
p <- plotly::plot_ly(type="scatter3d",mode="markers",showlegend=TRUE)%>%
plotly::add_trace(p,
x = ~m_z, y = ~RT, z = ~Intensity,
data = sub_MSlist[sub_MSlist[,"peakID"]==0,],
color=I("gray"),
size = I(1),
name = "",
error_z=list(
color="gray",
thickness=0,
symmetric = TRUE,
type = "data" ,
array = sub_MSlist[sub_MSlist[,"peakID"]==0,]$Intensity
)
)%>%
plotly::add_trace(p,
x = ~m_z, y = ~RT, z = ~Intensity,
data = sub_MSlist[sub_MSlist[,"peakID"]!=0,],
color=I("red"),
size = I(1),
name = "",
error_z=list(
color="red",
thickness=0,
symmetric = TRUE,
type = "data" ,
array = sub_MSlist[sub_MSlist[,"peakID"]!=0,]$Intensity
)
)%>%
plotly::add_trace(p,
x = ~m_z, y = ~RT, z = ~Intensity,
data = sub_peaks,
color=I("black"),
size = I(1),
name = "",
error_z=list(
color="black",
thickness=0,
symmetric = TRUE,
type = "data" ,
array = sub_peaks$Intensity
)
)
})
}else{
output$plot_peaks_3D <- renderPlotly({plotly::plot_ly(type="scatter3d",mode="markers",showlegend=TRUE)})
}
}
if( !isolate(input$peaks_mz_RT_use_peaks) & !isolate(input$peaks_mz_RT_use_raw)){
output$plot_peaks_3D <- renderPlotly({plotly::plot_ly(type="scatter3d",mode="markers",showlegend=TRUE)})
}
}else{
output$file_viewer_name<-renderText("Invalid file ID")
output$file_viewer_type<-renderText("")
output$file_viewer_mode<-renderText("")
}
}else{
output$plot_peaks_mz_RT <- renderPlot({})
output$plot_peaks_mz_int <- renderPlot({})
output$plot_peaks_RT_int <- renderPlot({})
output$plot_peaks_3D <- renderPlotly({plotly::plot_ly(type="scatter3d",mode="markers",showlegend=TRUE)})
output$file_viewer_name<-renderText("Invalid file ID")
output$file_viewer_type<-renderText("")
output$file_viewer_mode<-renderText("")
}
}
}
})
observeEvent(input$plot_peaks_mz_RT_dblclick, {
brush <- input$plot_peaks_mz_RT_brush
if (!is.null(brush)) {
isolate(ranges_peaks_mz_RT$x <- c(brush$xmin, brush$xmax))
isolate(ranges_peaks_mz_RT$y <- c(brush$ymin, brush$ymax))
} else {
isolate(ranges_peaks_mz_RT$x <- NULL)
isolate(ranges_peaks_mz_RT$y <- NULL)
}
refresh_plot$a<-(refresh_plot$a+1)
if(logfile$parameters$verbose) cat("\n Zooming with brush")
})
observe({
input$peaks_mz_RT_zoom_out
if(isolate(init$a)=="TRUE"){
if(!is.na(isolate(input$sel_meas_ID))){
if(!is.null(isolate(ranges_peaks_mz_RT$x))){
if(logfile$parameters$verbose) cat("\n Zooming out on X")
old_range<-abs(isolate(ranges_peaks_mz_RT$x[2]-ranges_peaks_mz_RT$x[1]))
isolate(ranges_peaks_mz_RT$x[1]<-ranges_peaks_mz_RT$x[1]-.5*old_range)
isolate(ranges_peaks_mz_RT$x[2]<-ranges_peaks_mz_RT$x[2]+.5*old_range)
}
refresh_plot$a<-(refresh_plot$a+1)
}
}
})
observe({
input$sel_meas_ID
input$sel_peak_ID
if(!is.na(isolate(input$sel_meas_ID))){
if(isolate(input$sel_meas_ID)!=0){
if(atit==isolate(input$sel_meas_ID) & exists("MSlist")){
if( !is.na(isolate(input$sel_peak_ID)) &
(isolate(input$sel_peak_ID)!=0) &
any(MSlist[[8]][,10]==isolate(input$sel_peak_ID)) &
any(objects(envir=as.environment(".GlobalEnv"))=="MSlist")
){
EIC_ID<<-unique(MSlist[[8]][MSlist[[8]][,10]==isolate(input$sel_peak_ID),9]);
peakit<<-MSlist[[4]][[2]][c(MSlist[[7]][as.numeric(isolate(input$sel_peak_ID)),1]:MSlist[[7]][
as.numeric(isolate(input$sel_peak_ID)),2]),]
if(length(peakit)>7){
EICit<<-MSlist[[4]][[2]][c(MSlist[[6]][EIC_ID,1]:MSlist[[6]][EIC_ID,2]),]
output$EIC1 <- renderPlot({
if(length(EICit)>7){
plot(EICit[,3],EICit[,2],type="h",col="darkgrey",xlab="RT",ylab="Intensity",main="EIC (grey) & Peak (red)",xlim=c(min(MSlist[[4]][[1]]),max(MSlist[[4]][[1]])))
}else{
plot(EICit[3],EICit[2],type="h",col="darkgrey",xlab="RT",ylab="Intensity",main="EIC (grey) & Peak (red)")
}
if(length(peakit)>7){
points(peakit[,3],peakit[,2],type="h",col="red",lwd=2)
}else{
points(peakit[3],peakit[2],type="h",col="red",lwd=2)
}
})
output$EIC2 <- renderPlot({
if(length(EICit)>7){
plot(EICit[,3],EICit[,2],type="h",col="darkgrey",xlab="RT",ylab="Intensity")
}else{
plot(EICit[3],EICit[2],type="h",col="darkgrey",xlab="RT",ylab="Intensity")
}
if(length(peakit)>7){
points(peakit[,3],peakit[,2],type="h",col="red",lwd=2)
}else{
points(peakit[3],peakit[2],type="h",col="red",lwd=2)
}
})
output$EIC3 <- renderPlot({
if(length(EICit)>7){
plot(EICit[,3],EICit[,1],pch=19,col="darkgrey",xlab="RT",ylab="m/z")
}else{
plot(EICit[3],EICit[1],pch=19,col="darkgrey",xlab="RT",ylab="m/z")
}
if(length(peakit)>7){
points(peakit[,3],peakit[,1],pch=19,col="red",cex=1.5)
}else{
points(peakit[,3],peakit[,1],pch=19,col="red",cex=1.5)
}
})
if(logfile$parameters$verbose) cat("\n EIC & peak extracted")
}else{
if(logfile$parameters$verbose) cat("\n Peak based on single measurement - plotting skipped.")
}
}
}else{
output$EIC1 <- renderPlot({plot.new()})
output$EIC2 <- renderPlot({plot.new()})
output$EIC3 <- renderPlot({plot.new()})
}
}
}
})
observe({
input$reset_1
if( (isolate(init$a) == "TRUE") & isolate(input$reset_1) ){
if(any(ls() == "logfile")){stop(paste("\n illegal logfile detected in server_obs_res_mean.r #1"))}
logfile$Tasks_to_redo <<- replace(logfile$Tasks_to_redo,-1, TRUE)
measurements <- read.csv(file=file.path(logfile[[1]],"dataframes","measurements"), colClasses = "character");
measurements[,c("qc","recal","align","norm", "LOD","isotopologues","adducts","homologues","EIC_correlation","blind","components_files")]<-"FALSE"
write.csv(measurements, file = file.path(logfile[[1]],"dataframes","measurements"), row.names = FALSE);
createAlert(session, anchorId = "reset", alertId = "reset1", title = NULL, content = "Project reset w/o peak picking", style = "warning",append=FALSE,dismiss=TRUE)
output$summa_html <- renderText(enviMass::UCBRR9756L(logfile$summary, logfile$Tasks_to_redo));
save(logfile, file = file.path(as.character(logfile[[1]]),"logfile.emp"));
if(logfile$parameters$verbose) cat("\nReset without peak picking \n")
}
})
observe({
input$reset_2
if( (isolate(init$a)=="TRUE") & isolate(input$reset_2) ){
if(any(ls()=="logfile")){stop(paste("\n illegal logfile detected in server_obs_res_mean.r #1"))}
logfile$Tasks_to_redo <<- replace(logfile$Tasks_to_redo,,TRUE)
measurements <- read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
if(TRUE){
measurements[,c("peakpicking","qc","recal","align","norm", "LOD","isotopologues","adducts","homologues","EIC_correlation","blind","components_files")]<-"FALSE"
}
write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
rm(measurements)
those <- list.files(file.path(logfile$project_folder,"peaklist"))
if(length(those) > 0){
for(i in 1:length(those)){
file.remove(file.path(logfile$project_folder,"peaklist",those[i]))
}
}
createAlert(session, anchorId = "reset", alertId = "reset2", title = NULL, content = "Project reset", style = "warning", append = FALSE, dismiss = TRUE)
output$summa_html <- renderText(enviMass::UCBRR9756L(logfile$summary, logfile$Tasks_to_redo));
save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
if(logfile$parameters$verbose) cat("\nTotal reset \n")
}
})
observe({
input$resolution
if( (isolate(init$a)=="TRUE") & (isolate(input$resolution)!="none")  & (isolate(input$resolution)!="") ){
#cat(input$resolution);cat("\n")
path=file.path(logfile[[1]],"pics","resolution")
png(filename = path, bg = "white")
that<-resolution_list[names(resolution_list) == as.character(input$resolution)][[1]]
plot(that[,1],that[,2],pch=19,cex=0.5,xlab="m/z",ylab="Resolution")
dev.off()
exprres <- list(src=file.path(logfile[[1]],"pics","resolution"))
output$plot_resolution <- renderImage(exprres, deleteFile = FALSE)
}
})
observe({
input$sel_scans_ID
input$method_definition
input$sel_scans_number
input$method_MS1_separation
if( (isolate(init$a) == "TRUE") & (!is.na(as.character(isolate(input$sel_scans_ID))))){
if(logfile$parameters$verbose == "TRUE") cat("\n Retrieving scan information _A")
path = file.path(logfile[[1]], "files", paste0(as.character(isolate(input$sel_scans_ID)), ".mzXML"))
measurements <- read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
if((file.exists(path)) & any(measurements$ID == as.character(isolate(input$sel_scans_ID))) & length(isolate(input$method_definition))){
if(logfile$parameters$verbose == "TRUE") cat("\n Retrieving scan information _B")
mzXML_file <- mzR:::openMSfile(filename = path, backend = c("Ramp"), verbose = FALSE)
output$scan_viewer_name <- renderText({
paste0("File name: ", measurements[measurements$ID == as.character(isolate(input$sel_scans_ID)), "Name"])
})
output$scan_viewer_type <- renderText({
paste0("File type: ", measurements[measurements$ID == as.character(isolate(input$sel_scans_ID)), "Type"])
})
output$scan_viewer_mode <- renderText({
paste0("Tabulated ionization mode: ", measurements[measurements$ID == as.character(isolate(input$sel_scans_ID)), "Mode"])
})
file_name <- mzR::fileName(mzXML_file)
run_Info <- mzR::runInfo(mzXML_file)
run_Info$msLevels <- paste0(run_Info$msLevels, collapse = ", ")
run_Info <- as.data.frame(run_Info)
for_scans <- isolate(input$sel_scans_number)
if(for_scans < 1) for_scans <- 1
if(for_scans > run_Info$scanCount) for_scans <- run_Info$scanCount
instrument_Info <- as.data.frame(mzR::instrumentInfo(mzXML_file))
polar <- c("-", "+")
method_definition <- isolate(input$method_definition)
heads <- mzR::header(mzXML_file)
method_definition2 <- method_definition[!is.na(match(method_definition, names(heads)))]
if(length(method_definition2) != length(method_definition)){
stop("\n WARNING: scan header declaration conflict for your files - please report this problem!")
}
names(heads)[names(heads) == "peaksCount"] <- "CentroidCount"
peaksCount <- heads[, "CentroidCount"]
heads_short <- heads[, method_definition2, drop = FALSE]
heads_summary <- unique(heads_short)
if(any(names(heads_summary) == "msLevel")){
heads_summary <- heads_summary[order(heads_summary$msLevel, decreasing = FALSE),, drop = FALSE]
}
scanTypes2 <- seq(1:dim(heads_summary)[1])
scanTypes <- match(data.frame(t(heads_short)), data.frame(t(heads_summary)))
if(any(names(heads_summary) == "msLevel")){
if(
isolate(input$method_MS1_separation) &
any(heads_summary$msLevel == 1) &
any(heads_short[,"msLevel"] == 1)
){
separation_level <- rep(1, length(scanTypes2))
scan_level <- rep(1, length(scanTypes))
for(n in 2:length(scanTypes)){
if(heads_summary[scanTypes[n], "msLevel"] != 1) next
if(heads_summary[scanTypes[n - 1], "msLevel"] != 1) next
if(scanTypes[n] != scanTypes[n - 1]) next
scan_level[n] <- (scan_level[n - 1] + 1)
if( scan_level[n] > separation_level[scanTypes[n]] ) separation_level[scanTypes[n]] <- scan_level[n]
}
if(any(separation_level > 0)){
scanTypes2_new <- seq(1:(max(scanTypes2) + sum(separation_level - 1)))
scanTypes_new <- rep(0, length(scanTypes))
use_Type <- 0
MS1_consec <- c()
table_expand <- c()
for(n in 1:length(scanTypes2)){
for(m in 1:separation_level[n]){
use_Type <- (use_Type + 1)
scanTypes_new[
(scanTypes == scanTypes2[n]) &
(scan_level == m)
] <- use_Type
MS1_consec <- c(MS1_consec, (m - 1))
table_expand <- c(table_expand, n)
}
}
heads_summary <- heads_summary[table_expand,, drop = FALSE]
heads_summary <- cbind(heads_summary, MS1_consec)
scanTypes <- scanTypes_new
scanTypes2 <- scanTypes2_new
}
}
}
scanCounts <- rep(0, length(scanTypes2))
centroidCounts <- rep(0, length(scanTypes2))
for(n in 1:length(scanTypes2)){
these <- which(scanTypes == scanTypes2[n])
scanCounts[n] <- length(these)
centroidCounts[n] <- sum(peaksCount[these])
}
heads_summary <- cbind(scanTypes2, as.integer(scanCounts), as.integer(centroidCounts), heads_summary)
names(heads_summary)[1] <- "Scan type"
names(heads_summary)[2] <- "Scan counts"
names(heads_summary)[3] <- "Centroid counts"
if(any(names(heads_summary) == "polarity")){
heads_summary$polarity <- polar[(heads_summary$polarity)+1]
}
heads_summary <<- heads_summary
updateCheckboxGroupInput(session, "method_use_ScanTypes", choices = as.character(scanTypes2), inline = TRUE)
output$instrument_Info <- renderTable(instrument_Info)
output$run_Info <- renderTable(run_Info)
output$heads_summary_new <- renderTable(heads_summary)
use_m <- vector("list", for_scans)
which_col <- !names(heads) %in% c("msLevel", "seqNum", "acquisitionNum", "polarity")
head_dim <- dim(heads[, which_col])[2]
for(n in 1:for_scans){
sub_list <- as.list(heads[n, which_col])
sub_list[1] <- structure(sub_list[1], sticon = "signal")
names(sub_list) <- paste(names(sub_list), heads[n, which_col], sep = ": ")
if(heads[n, "msLevel"] == 1){
use_m[[n]] <- sub_list
}else{
use_m[[n]] <- structure(sub_list, sticon = "signal")
}
names(use_m)[n] <- paste0(
polar[heads[n, "polarity"] + 1],
" / MS level: ",
heads[n, "msLevel"],
" / Scan number: ",
heads[n, "seqNum"],
" - Scan type: ",
scanTypes[n]
)
for(m in 1:head_dim) attr(use_m[[n]][[m]], "sticon") <- "empty"
}
output$scan_tree <- renderTree(use_m, quoted = FALSE)
mzR:::close(mzXML_file)
}else{
if(logfile$parameters$verbose == "TRUE") cat("\n Retrieving scan information _C")
if(!length(isolate(input$method_definition))){
output$scan_viewer_name <- renderText({paste0("No scan type definition parameters selected - please select at least one!")})
}else{
if(!file.exists(path) & any(measurements$ID == as.character(isolate(input$sel_scans_ID)))){
output$scan_viewer_name <- renderText({paste0(".mzXML file not available")})
}else{
output$scan_viewer_name <- renderText({paste0("File name: Invalid file ID")})
}
}
}
}
})
observe({
input$save_method
if(
(isolate(init$a) == "TRUE") & isolate(input$save_method) &
(any(objects(envir = as.environment(".GlobalEnv")) == "heads_summary"))
){
if(logfile$parameters$verbose) cat("\n Saving method ...")
use_ScanTypes <- isolate(input$method_use_ScanTypes)
all_ok <- TRUE
if(is.null(use_ScanTypes)){
shinytoastr::toastr_error("No Scan type to include selected. Please use the green check box field next to the Save method button for this first.", title = "Method setup error:", closeButton = TRUE, position = c("top-center"), timeOut = 0);
output$dowhat <- renderText("Method setup failed");
if(logfile$parameters$verbose) cat(" failed.\n")
all_ok <- FALSE
}
if(!any(names(heads_summary) == "msLevel")){
shinytoastr::toastr_error("No msLevel included in method definition - should be included to differentiate for msLevel 1 scans. Please revise!", title = "Method setup error:", closeButton = TRUE, position = c("top-center"), timeOut = 0);
output$dowhat <- renderText("Method setup failed");
if(logfile$parameters$verbose) cat(" failed.\n")
all_ok <- FALSE
}
if(all_ok){
say <- "Existing method (if any) replaced."
used_scans <- rep(FALSE, dim(heads_summary)[1])
used_scans[as.numeric(use_ScanTypes)] <- TRUE
heads_summary_existing <- cbind(heads_summary, used_scans)
output$heads_summary_existing <- renderTable(heads_summary_existing)
if(any(names(heads_summary) == "msLevel")){
if(sum(heads_summary[used_scans,"msLevel"] == 1) > 1){
say <- paste(say, "BEWARE: more than one msLevel 1 Scan types included. These will be pooled during processing!")
}
}
if(any(names(heads_summary) == "polarity")){
say <- paste(say, "BEWARE: do you process switch mode files? The polarity filter may not make sense otherwise.")
}
logfile$method_setup <<- heads_summary_existing
save(logfile, file = file.path(as.character(logfile[["project_folder"]]), "logfile.emp"));
if( as.logical(logfile$parameters$method_use) ){
enviMass::RNJJT9347D(
down = "peakpicking",
check_node = FALSE,
single_file = FALSE
)
output$summa_html <- renderText(enviMass::UCBRR9756L(logfile$summary, logfile$Tasks_to_redo));
}
if(any(ls()=="logfile")){stop("\n illegal logfile detected during method saving!")}
shinytoastr::toastr_success(say, title = "Method setup saved", closeButton = TRUE);
output$dowhat <- renderText("Method setup saved");
if(logfile$parameters$verbose) cat(" done.\n")
}
}
})
