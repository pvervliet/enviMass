if(any(ls() == "logfile")){stop("\n illegal logfile detected #1 in server_obs_screening.r!")}
ee	<-	reactiveValues()
ee$entry <- 0
verbose <- TRUE
ranges_homol <- reactiveValues(mass = FALSE, RT = FALSE, massD = FALSE, dmass = FALSE, dRT = FALSE, RTchrom = FALSE, intchrom = FALSE)
refresh_homol <- reactiveValues()
refresh_homol$a <- 0
refresh_homol$b <- 0
refresh_homol$c <- 0
refresh_homol$d <- 0
ranges_compo <- reactiveValues(mass = FALSE, RTchrom = FALSE, intchrom = FALSE)
refresh_compo <- reactiveValues()
refresh_compo$a <- 0
observe({
input$sel_meas_comp
if(isolate(init$a)=="TRUE"){
if(logfile$parameters$verbose) cat("\n in Comp_A")
do_isot <- (logfile$workflow[names(logfile$workflow) == "isotopologues"] == "yes")
do_addu <- (logfile$workflow[names(logfile$workflow) == "adducts"] == "yes")
do_homol <- (logfile$workflow[names(logfile$workflow) == "homologues"] == "yes")
do_EIC <- (logfile$workflow[names(logfile$workflow) == "EIC_correlation"] == "yes")
measurements <- read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
if(
(!is.na(isolate(input$sel_meas_comp))) &
(isolate(input$sel_meas_comp) != "") &
(any(measurements$ID == isolate(input$sel_meas_comp))) &
file.exists(file.path(logfile[[1]], "MSlist", isolate(input$sel_meas_comp))) &
file.exists(file.path(logfile[[1]], "peaklist", isolate(input$sel_meas_comp)))
){
output$sel_meas_comp_state <- renderText("For this file:")
output$comp_file_name <- renderText(paste("File name: ", measurements[measurements[,"ID"] == as.character(isolate(input$sel_meas_comp)),"Name"],sep=""))
output$comp_file_type <- renderText(paste("File type: ", measurements[measurements[,"ID"] == as.character(isolate(input$sel_meas_comp)),"Type"],sep=""))
output$comp_file_mode <- renderText(paste("Ionization mode: ", measurements[measurements[,"ID"] == as.character(isolate(input$sel_meas_comp)),"Mode"],sep=""))
load(file.path(logfile[[1]], "MSlist", isolate(input$sel_meas_comp)), envir = as.environment(".GlobalEnv"))
load(file.path(logfile[[1]], "peaklist", isolate(input$sel_meas_comp)), envir = as.environment(".GlobalEnv"))
if(
file.exists(file.path(logfile[[1]], "results", "componentization", "components", isolate(input$sel_meas_comp))) &
(do_isot | do_addu) &
(measurements[measurements$ID==isolate(input$sel_meas_comp),"include"] == "TRUE")
){
if(do_isot){
updateSelectInput(session,inputId="atom_bound_addpeaks",
choices=c("(a) peaks in same isotopologue group","(b) all peaks with similar RT"),
selected="(a) peaks in same isotopologue group"
)
}else{
updateSelectInput(session,inputId="atom_bound_addpeaks",
choices=c("(b) all peaks with similar RT"),
selected="(b) all peaks with similar RT"
)
}
load(file.path(logfile[[1]],"results","componentization","components",isolate(input$sel_meas_comp)),envir=as.environment(".GlobalEnv"))
if(logfile$parameters$verbose) cat("\n in Comp_A_1: loaded file")
if(!is.null(dim(component[["pattern peak list"]]))){
num_peaks_remain <- dim(component[["pattern peak list"]])[1]
max_peak_ID <<- max(component[["pattern peak list"]][,"peak ID"])
}else{
num_peaks_remain <- dim(component[["adduct peak list"]])[1]
max_peak_ID <<- max(component[["adduct peak list"]][,"peak ID"])
}
num_comp<-dim(component[["Components"]])[1]
num_comp_tar<-sum(component[["Components"]][,"Target peaks"]!="-")
num_comp_ISTD<-sum(component[["Components"]][,"ISTD peaks"]!="-")
blind_num<-as.numeric(component[["Components"]][,"Blind peak number"])
tot_num<-as.numeric(component[["Components"]][,"Total peak number"])
num_comp_blind_any<-sum(blind_num>0)
num_comp_blind_all<-sum((blind_num/tot_num)==1)
comp_nontarget<-sum(
((blind_num/tot_num)<1) &
(component[["Components"]][,"Target peaks"]=="-") &
(component[["Components"]][,"ISTD peaks"]=="-")
)
reduc<-round((num_peaks_remain/num_comp),digits=2)
num_isot_peaks<-rep(0,num_comp)
num_adduc_peaks<-rep(0,num_comp)
for(i in 1:num_comp){
if(component[["Components"]][i,3]!="-"){
num_isot_peaks[i]<-length(strsplit(component[["Components"]][i,3],",")[[1]])
}
if(component[["Components"]][i,5]!="-"){
num_adduc_peaks[i]<-length(strsplit(component[["Components"]][i,5],",")[[1]])
}
}
min2_size_comp<-round((sum((num_isot_peaks+num_adduc_peaks)>1)/num_comp),digits=2)
median_size_comp<-round(mean(num_isot_peaks+num_adduc_peaks),digits=2)
max_size_comp<-max(num_isot_peaks+num_adduc_peaks)
output$num_peaks_all<-renderText(paste("Total number of picked peaks for selected file: ",as.character(length(peaklist[,1])),sep=""))
if(logfile$workflow[names(logfile$workflow)=="replicates"]=="yes"){
output$num_peaks_remain_replicate<-renderText(paste("Number of peaks removed by replicate filter: ",as.character((sum(peaklist[,colnames(peaklist)=="keep"]==0))),sep=""))
}else{
output$num_peaks_remain_replicate<-renderText(paste("Number of peaks removed by replicate filter: workflow step not enabled."))
}
if(logfile$workflow[names(logfile$workflow)=="blind"]=="yes"){
output$num_peaks_remain_blind<-renderText(paste("Number of peaks also found in blind files (peaks not removed): ",as.character((sum(peaklist[,colnames(peaklist)=="keep_2"]<Inf))),sep=""))
}else{
output$num_peaks_remain_blind<-renderText(paste("Number of peaks affected by blind peak detection: workflow step not enabled."))
}
output$num_peaks_remain<-renderText(paste("Remaining number of peaks: ",as.character(num_peaks_remain),sep=""))
output$num_comp<-renderText(paste("Total number of components build from remaining peaks: ",as.character(num_comp),sep=""))
output$reduc<-renderText(paste("Reduction factor: ",as.character(reduc),sep=""))
output$num_comp_tar<-renderText(paste("Components containing target or suspect compound peaks: ",as.character(num_comp_tar),sep=""))
output$num_comp_ISTD<-renderText(paste("Components containing ISTD peaks: ",as.character(num_comp_ISTD),sep=""))
output$num_comp_blind_any<-renderText(paste("Components containing any blank/blind peaks: ",as.character(num_comp_blind_any),sep=""))
output$num_comp_blind_all<-renderText(paste("Components containing only blank/blind peaks: ",as.character(num_comp_blind_all),sep=""))
output$min2_size_comp<-renderText(paste("Fraction of components with min. 2 peaks: ",as.character(min2_size_comp),sep=""))
output$median_size_comp<-renderText(paste("Mean number of peaks per component: ",as.character(median_size_comp),sep=""))
output$max_size_comp<-renderText(paste("Max number of peaks in a component: ",as.character(max_size_comp),sep=""))
output$num_comp_nontarget<-renderText(paste("Number of nontarget components with at least one non-blind peak: ",as.character(comp_nontarget),sep=""))
names_max_atom <- names(component[["Components"]])[grepl("max_atom_", names(component[["Components"]]))]
comp_table_full <<- component[["Components"]][,c(
"Component ID |",
"Monois. peak ID |",
"Monois. m/z |",
"Monois. RT |",
"Monois. RT |",
"Monois. int. |",
"Monois. sample/blind int. ratio",
"ID pattern peaks |",
"ID adduct peaks |",
"ID homologue series |",
"ID interfering peaks |",
"pattern group adduct|",
"adduct group adduct(s) |",
"Target peaks",
"ISTD peaks",
"Total peak number",
"Blind peak number",
"z",
names_max_atom
), drop = FALSE]
comp_table_full[,3]<-round(comp_table_full[,3],digits=5)
comp_table_full[,4]<-round(comp_table_full[,4],digits=2)
comp_table_full[,5]<-(comp_table_full[,5]/60)
comp_table_full[,5]<-round(comp_table_full[,5],digits=2)
comp_table_full[,6]<-round(comp_table_full[,6],digits=1)
comp_table_full[,7]<-round(comp_table_full[,7],digits=3)
output$comp_table_full <- DT::renderDataTable(
DT::datatable(
comp_table_full,
colnames=c(
"Component ID",
"Monois. peak ID",
"Monois. peak m/z",
"Monois. peak RT [s]","Monois. peak RT [min]",
"Monois. peak intens.",
"Monois. sample/blind int. ratio",
"ID(s) isot. peaks","ID(s) adduct peaks","ID(s) homol. series","ID(s) interfering peaks",
"Isot. peaks adducts","Adduct peak adducts",
"Target peaks","ISTD peaks",
"Total peak number","Blind peak number",
"Charge(s)", names_max_atom
),
rownames=FALSE,
extensions = c('Buttons','FixedHeader','ColReorder'),
options = list(
lengthMenu = list(c(25, 50, 200, -1), list('25', '50', '200', 'All')),
fixedHeader = FALSE,
ordering=T,
dom = 'Blfrtip',
buttons = c('excel', 'csv','colvis'),#buttons = c('excel', 'pdf', 'print', 'csv'),
scrollX = TRUE,
scrollY = "800px",
colReorder = TRUE
),
filter = 'top',
selection = list(mode = 'single', target = 'row')
),
server = TRUE
)
if((length(component[["pattern peak list"]])>1) & (do_isot)){found_isos<-TRUE}else{found_isos<-FALSE}
if((length(component[["adduct peak list"]])>1) & (do_addu)){found_addu<-TRUE}else{found_addu<-FALSE}
if( found_isos & !found_addu ){
output$sel_meas_comp_state1<-renderText("isotopologue grouping available ")
}
if( !found_isos & found_addu ){
output$sel_meas_comp_state1<-renderText("adduct grouping available ")
}
if( found_isos & found_addu ){
output$sel_meas_comp_state1<-renderText("isotopologue & adduct grouping available ")
}
}else{
output$sel_meas_comp_state1<-renderText("no nontarget components available ")
}
if(
file.exists(file.path(logfile[[1]],"results","componentization","homologues",isolate(input$sel_meas_comp))) &
do_homol &
measurements[measurements$ID==isolate(input$sel_meas_comp),"include"]=="TRUE"
){
if(logfile$parameters$verbose) cat("\n in Comp_A_3")
load(file.path(logfile[[1]],"results","componentization","homologues",paste("full",isolate(input$sel_meas_comp),sep="_")),envir=as.environment(".GlobalEnv"))
isolate(refresh_homol$a<-(refresh_homol$a+1))
output$sel_meas_comp_state2<-renderText(" homologue series detection results available")
}else{
output$sel_meas_comp_state2<-renderText(" no homologue series detection results available")
}
}else{
output$sel_meas_comp_state<-renderText("Invalid file ID")
output$comp_file_name<-renderText("")
output$comp_file_type<-renderText("")
output$comp_file_mode<-renderText("")
output$sel_meas_comp_state1<-renderText("")
output$sel_meas_comp_state2<-renderText("")
output$num_peaks_all<-renderText("")
output$num_peaks_remain_replicate<-renderText("")
output$num_peaks_remain_blind<-renderText("")
output$num_peaks_remain<-renderText("")
output$num_comp<-renderText("")
output$reduc<-renderText("")
output$min2_size_comp<-renderText("")
output$median_size_comp<-renderText("")
output$max_size_comp<-renderText("")
output$num_comp_tar<-renderText("")
output$num_comp_ISTD<-renderText("")
output$num_comp_blind_any<-renderText("")
}
}
})
observe({
refresh_homol$a
if(isolate(init$a)=="TRUE"){
cat("\n IN REFRESH_1")
plot_those<<-enviMass:::EPOZH3288F(
homol,
masslim = isolate(ranges_homol$mass),
RTlim = isolate(ranges_homol$RT),
massDlim = isolate(ranges_homol$massD),
dmasslim = isolate(ranges_homol$dmass),
dRTlim = isolate(ranges_homol$dRT)
)
sum_homol<-sum(plot_those)
if(sum_homol>1000){
omit_theta<<-2
if(sum_homol>5000){omit_theta<<-10}
if(sum_homol>20000){omit_theta<<-15}
if(sum_homol>100000){omit_theta<<-20}
if(sum_homol>200000){omit_theta<<-30}
}else{
omit_theta<<-FALSE
}
output$homol_plot <- renderPlot({
par(mar=c(4.5,4.5,.9,.8))
enviMass:::XTUOY4153L(homol,
xlim = isolate(ranges_homol$mass), ylim = isolate(ranges_homol$RT),
dmasslim = isolate(ranges_homol$dmass), dRTlim = isolate(ranges_homol$dRT),
plot_what="mz_RT", omit_theta=omit_theta, plot_those = plot_those
);
},res=100)
output$homol_counts <- renderPlot({
par(mar=c(4.5,4.5,.9,.8))
enviMass:::XTUOY4153L(homol,
xlim = isolate(ranges_homol$mass), ylim = isolate(ranges_homol$RT),
dmasslim = isolate(ranges_homol$dmass), dRTlim = isolate(ranges_homol$dRT),
plot_what="dmass", omit_theta=FALSE, plot_those = plot_those
);
},res=100)
# output homol. RT difference vs mass scatter plot ##############
output$homol_RT <- renderPlot({
par(mar=c(4.5,4.5,.9,.8))
enviMass:::XTUOY4153L(homol,
xlim = isolate(ranges_homol$mass), ylim = isolate(ranges_homol$RT),
dmasslim = isolate(ranges_homol$dmass), dRTlim = isolate(ranges_homol$dRT),
plot_what="dRT", omit_theta=FALSE, plot_those = plot_those
);
},res=100)
use_homol_peaks<-unique(c(homol[["homol_peaks_relat"]][plot_those,1],homol[["homol_peaks_relat"]][plot_those,2]))
output$homol_series_peaks <- DT::renderDataTable(
DT::datatable(
dat<-data.frame(
I(as.character(homol[["Peaks in homologue series"]][use_homol_peaks,c("peak ID")])),
round(homol[["Peaks in homologue series"]][use_homol_peaks,c("mz")],digits=5),
round(log10(homol[["Peaks in homologue series"]][use_homol_peaks,c("intensity")]),digits=4),
round(homol[["Peaks in homologue series"]][use_homol_peaks,c("RT")],digits=1),
round((homol[["Peaks in homologue series"]][use_homol_peaks,c("RT")]/60),digits=1),
I(homol[["Peaks in homologue series"]][use_homol_peaks,c("Targets")]),
I(homol[["Peaks in homologue series"]][use_homol_peaks,c("ISTDs")]),
I(homol[["Peaks in homologue series"]][use_homol_peaks,c("HS IDs")])
),
filter = 'top',
colnames=c("Peak ID","m/z","log10 intens.","RT [s]","RT [min]","Target matches","ISTD matches","Series ID(s)"),
rownames=FALSE,
extensions = c('Buttons','FixedHeader','ColReorder'),
options = list(
lengthMenu = list(c(25, 50, 200, -1), list('25', '50', '200', 'All')),
fixedHeader = FALSE,
ordering=T,
dom = 'Blfrtip',
buttons = c('excel', 'csv','colvis'),#buttons = c('excel', 'pdf', 'print', 'csv'),
scrollX = TRUE,
scrollY = "800px",
colReorder = TRUE
),
selection = list(mode = 'single', target = 'row')
),
server = TRUE
)
output$homol_series_table <- DT::renderDataTable({
DT::datatable(
data.frame(
I(as.character(homol[["Homologue Series"]][,"HS IDs"])),
I(as.character(homol[["Homologue Series"]][,"peak IDs"])),
round(homol[["Homologue Series"]][,"m/z increment"],digits=5),
round(homol[["Homologue Series"]][,"RT increment"],digits=1),
round(log10(homol[["Homologue Series"]][,"max int."]),digits=4)
),
filter = 'top',
colnames=c("Series ID","Peak IDs","m/z difference","RT difference [s]","Max log10 int."),
rownames=FALSE,
extensions = c('Buttons','FixedHeader','ColReorder'),
options = list(
lengthMenu = list(c(25, 50, 100, -1), list('25', '50', '100', 'All')),
fixedHeader = FALSE,
ordering=T,
dom = 'Blfrtip',
buttons = c('excel', 'csv','colvis'),#buttons = c('excel', 'pdf', 'print', 'csv'),
scrollX = TRUE,
scrollY = "800px",
colReorder = TRUE
),
selection = list(mode = 'single', target = 'row')
)
}, server = TRUE)
isolate(ranges_homol$RTchrom<-FALSE)
isolate(ranges_homol$intchrom<-FALSE)
}
})
observe({
refresh_homol$b
if(isolate(init$a)=="TRUE"){
cat("\n IN REFRESH_2")
plot_those<<-enviMass:::EPOZH3288F(
homol,
masslim = isolate(ranges_homol$mass),
RTlim = isolate(ranges_homol$RT),
massDlim = isolate(ranges_homol$massD),
dmasslim = isolate(ranges_homol$dmass),
dRTlim = isolate(ranges_homol$dRT)
)
sum_homol<-sum(plot_those)
if(sum_homol>1000){
omit_theta<<-2
if(sum_homol>5000){omit_theta<<-10}
if(sum_homol>20000){omit_theta<<-15}
if(sum_homol>100000){omit_theta<<-20}
if(sum_homol>200000){omit_theta<<-30}
}else{
omit_theta<<-FALSE
}
output$homol_plot <- renderPlot({
par(mar=c(4.5,4.5,.9,.8))
enviMass:::XTUOY4153L(homol,
xlim = isolate(ranges_homol$mass), ylim = isolate(ranges_homol$RT),
dmasslim = isolate(ranges_homol$dmass), dRTlim = isolate(ranges_homol$dRT),
plot_what="mz_RT", omit_theta=omit_theta, plot_those = plot_those
);
},res=100)
use_homol_peaks<-unique(c(homol[["homol_peaks_relat"]][plot_those,1],homol[["homol_peaks_relat"]][plot_those,2]))
output$homol_series_peaks <- DT::renderDataTable(
DT::datatable(
dat<-data.frame(
I(as.character(homol[["Peaks in homologue series"]][use_homol_peaks,c("peak ID")])),
round(homol[["Peaks in homologue series"]][use_homol_peaks,c("mz")],digits=5),
round(log10(homol[["Peaks in homologue series"]][use_homol_peaks,c("intensity")]),digits=4),
round(homol[["Peaks in homologue series"]][use_homol_peaks,c("RT")],digits=1),
round((homol[["Peaks in homologue series"]][use_homol_peaks,c("RT")]/60),digits=1),
I(homol[["Peaks in homologue series"]][use_homol_peaks,c("Targets")]),
I(homol[["Peaks in homologue series"]][use_homol_peaks,c("ISTDs")]),
I(homol[["Peaks in homologue series"]][use_homol_peaks,c("HS IDs")])
),
filter = 'top',
colnames=c("Peak ID","m/z","log10 intens.","RT [s]","RT [min]","Target matches","ISTD matches","Series ID(s)"),
rownames=FALSE,
extensions = c('Buttons','FixedHeader','ColReorder'),
options = list(
lengthMenu = list(c(25, 50, 100, -1), list('25', '50', '100', 'All')),
fixedHeader = FALSE,
ordering=T,
dom = 'Blfrtip',
buttons = c('excel', 'csv','colvis'),#buttons = c('excel', 'pdf', 'print', 'csv'),
scrollX = TRUE,
scrollY = "800px",
colReorder = TRUE
),
selection = list(mode = 'single', target = 'row')
),
server = TRUE
)
}
})
observe({
s1<-input$homol_series_peaks_rows_selected
if(isolate(init$a)=="TRUE"){
if(length(s1)){
if(s1>=1){
print(s1);
cat("\n IN SELECT_1:")
use_homol_peaks<<-unique(c(homol[["homol_peaks_relat"]][plot_those,1],homol[["homol_peaks_relat"]][plot_those,2]))
these_series<-as.numeric(strsplit(homol[["Peaks in homologue series"]][use_homol_peaks[s1],c("HS IDs")],"/")[[1]])
print(these_series)
use_these_series<-match(these_series,homol[["Homologue Series"]][,"HS IDs"])
output$homol_series_table <- DT::renderDataTable({
DT::datatable(
data.frame(
I(as.character(homol[["Homologue Series"]][use_these_series,"HS IDs"])),
I(as.character(homol[["Homologue Series"]][use_these_series,"peak IDs"])),
round(homol[["Homologue Series"]][use_these_series,"m/z increment"],digits=5),
round(homol[["Homologue Series"]][use_these_series,"RT increment"],digits=1),
round(log10(homol[["Homologue Series"]][use_these_series,"max int."]),digits=4)
),
filter = 'top',
colnames=c("Series ID","Peak IDs","m/z difference","RT difference [s]","Max log10 int."),
rownames=FALSE,
extensions = c('Buttons','FixedHeader','ColReorder'),
options = list(
lengthMenu = list(c(25, 50, 100, -1), list('25', '50', '100', 'All')),
fixedHeader = FALSE,
ordering=T,
dom = 'Blfrtip',
buttons = c('excel', 'csv','colvis'),#buttons = c('excel', 'pdf', 'print', 'csv'),
scrollX = TRUE,
scrollY = "800px",
colReorder = TRUE
),
selection = list(mode = 'single', target = 'row')
)
},server = TRUE)
isolate(ranges_homol$RTchrom<-FALSE)
isolate(ranges_homol$intchrom<-FALSE)
cat("\n for this peak: ");print(use_homol_peaks[s1]);
output$homol_plot <- renderPlot({
par(mar=c(4.5,4.5,.9,.8))
enviMass:::XTUOY4153L(homol,
xlim = isolate(ranges_homol$mass), ylim = isolate(ranges_homol$RT),
dmasslim = isolate(ranges_homol$dmass), dRTlim = isolate(ranges_homol$dRT),
plot_what="mz_RT", omit_theta=omit_theta, plot_those = plot_those,
emph_point=homol[["Peaks in homologue series"]][use_homol_peaks[s1],"peak ID"],
emph_series=these_series
);
},res=100)
}
}else{
cat("\n ALL DESELECTED_1")
output$homol_series_table <- DT::renderDataTable({
DT::datatable(
data.frame(
I(as.character(homol[["Homologue Series"]][,"HS IDs"])),
I(as.character(homol[["Homologue Series"]][,"peak IDs"])),
round(homol[["Homologue Series"]][,"m/z increment"],digits=5),
round(homol[["Homologue Series"]][,"RT increment"],digits=1),
round(log10(homol[["Homologue Series"]][,"max int."]),digits=4)
),
filter = 'top',
colnames=c("Series ID","Peak IDs","m/z difference","RT difference [s]","Max log10 int."),
rownames=FALSE,
extensions = c('Buttons', 'FixedHeader', 'ColReorder'),
options = list(
lengthMenu = list(c(25, 50, 100, -1), list('25', '50', '100', 'All')),
fixedHeader = FALSE,
ordering=T,
dom = 'Blfrtip',
buttons = c('excel', 'csv','colvis'),#buttons = c('excel', 'pdf', 'print', 'csv'),
scrollX = TRUE,
scrollY = "800px",
colReorder = TRUE
),
selection = list(mode = 'single', target = 'row')
)
},server = TRUE)
output$homol_plot <- renderPlot({
par(mar=c(4.5,4.5,.9,.8))
enviMass:::XTUOY4153L(homol,
xlim = isolate(ranges_homol$mass), ylim = isolate(ranges_homol$RT),
dmasslim = isolate(ranges_homol$dmass), dRTlim = isolate(ranges_homol$dRT),
plot_what="mz_RT", omit_theta=omit_theta, plot_those = plot_those
);
},res=100)
}
}
})
observe({
s2<-input$homol_series_table_rows_selected
if(isolate(init$a)=="TRUE"){
s1<-isolate(input$homol_series_peaks_rows_selected)
if(length(s2)){
if(s2>=1){
cat("\n IN SELECT_2:")
if(length(s1)){
cat(" A ")
use_homol_peaks <- unique(c(homol[["homol_peaks_relat"]][plot_those,1],homol[["homol_peaks_relat"]][plot_those,2]))
these_series_1 <- as.numeric(strsplit(homol[["Peaks in homologue series"]][use_homol_peaks[s1],c("HS IDs")],"/")[[1]][s2])
use_emph_point <- homol[["Peaks in homologue series"]][use_homol_peaks[s1],"peak ID"]
}else{
cat(" B ")
these_series_1 <- homol[["Homologue Series"]][s2,"HS IDs"]
use_emph_point <- FALSE
}
print(these_series_1)
output$homol_plot <- renderPlot({
par(mar=c(4.5,4.5,.9,.8))
enviMass:::XTUOY4153L(homol,
xlim = isolate(ranges_homol$mass), ylim = isolate(ranges_homol$RT),
dmasslim = isolate(ranges_homol$dmass), dRTlim = isolate(ranges_homol$dRT),
plot_what = "mz_RT", omit_theta = omit_theta, plot_those = plot_those,
emph_point = use_emph_point,
emph_series = these_series_1
);
}, res = 100)
these_peaks <- as.numeric(strsplit(homol[["Homologue Series"]][these_series_1,"peak IDs"],",")[[1]])
isolate(refresh_homol$c <- these_peaks)
isolate(refresh_homol$d <- (refresh_homol$d+1))
}
}else{
cat("\n IN DESELECT_2:");print(s2);
if(length(s1)){
these_series <- as.numeric(strsplit(homol[["Peaks in homologue series"]][use_homol_peaks[s1],c("HS IDs")],"/")[[1]])
if(length(these_series)){
output$homol_plot <- renderPlot({
par(mar=c(4.5,4.5,.9,.8))
enviMass:::XTUOY4153L(homol,
xlim = isolate(ranges_homol$mass), ylim = isolate(ranges_homol$RT),
dmasslim = isolate(ranges_homol$dmass), dRTlim = isolate(ranges_homol$dRT),
plot_what = "mz_RT", omit_theta=omit_theta, plot_those = plot_those,
emph_point = homol[["Peaks in homologue series"]][use_homol_peaks[s1],"peak ID"],
emph_series = these_series
);
},res=100)
}
}else{
output$homol_plot <- renderPlot({
par(mar=c(4.5,4.5,.9,.8))
enviMass:::XTUOY4153L(homol,
xlim = isolate(ranges_homol$mass), ylim = isolate(ranges_homol$RT),
dmasslim = isolate(ranges_homol$dmass), dRTlim = isolate(ranges_homol$dRT),
plot_what="mz_RT", omit_theta=omit_theta, plot_those = plot_those
);
},res=100)
}
}
}
})
observe({
refresh_homol$d
if(isolate(init$a)=="TRUE"){
if(isolate(refresh_homol$c[1]>0) & isolate(refresh_homol$d[1]>0)){
cat("\n IN CHROMAT");print(isolate(refresh_homol$c));
output$homol_chromat <- renderPlot({
par(mar=c(4.5,4,.8,.8))
enviMass:::JJGET6467P(
MSlist,
peakIDs=isolate(refresh_homol$c),
RTlim=isolate(ranges_homol$RTchrom),
Intlim=isolate(ranges_homol$intchrom),
normalize=FALSE,
chromat_full=TRUE
);
},res=100)
}
}
})
observeEvent(input$homol_plot_dblclick, {
brush <- isolate(input$homol_plot_brush)
if (!is.null(brush)) {
if(logfile$parameters$verbose) cat("\n Zoom in_1")
isolate(ranges_homol$mass <- c(brush$xmin, brush$xmax))
isolate(ranges_homol$RT <- c(brush$ymin, brush$ymax))
} else {
if(logfile$parameters$verbose) cat("\n Zoom out full_1")
isolate(ranges_homol$mass <- FALSE)
isolate(ranges_homol$RT <- FALSE)
}
refresh_homol$a<-(refresh_homol$a+1)
})
observeEvent(input$homol_plot_click, {
if(logfile$parameters$verbose) cat("\n Zoom out part_1_a")
brush <- isolate(input$homol_plot_brush)
if (is.null(brush)) {
if(logfile$parameters$verbose) cat("\n Zoom out part_1_b")
if(isolate(ranges_homol$mass[1])!=FALSE){
old_range_mass<-abs(isolate(ranges_homol$mass[2]-ranges_homol$mass[1]))
isolate(ranges_homol$mass[1]<-ranges_homol$mass[1]-.3*old_range_mass)
isolate(ranges_homol$mass[2]<-ranges_homol$mass[2]+.3*old_range_mass)
}
if(isolate(ranges_homol$RT[1])!=FALSE){
old_range_RT<-abs(isolate(ranges_homol$RT[2]-ranges_homol$RT[1]))
isolate(ranges_homol$RT[1]<-ranges_homol$RT[1]-.1*old_range_RT)
isolate(ranges_homol$RT[2]<-ranges_homol$RT[2]+.1*old_range_RT)
}
refresh_homol$a<-(refresh_homol$a+1)
}else{
if(logfile$parameters$verbose) cat("\n Doing hover - and nothing")
}
})
observeEvent(input$homol_counts_dblclick, {
brush <- isolate(input$homol_counts_brush)
if (!is.null(brush)) {
if(logfile$parameters$verbose) cat("\n Zoom in_1d")
isolate(ranges_homol$dmass <- c(brush$xmin, brush$xmax))
} else {
if(logfile$parameters$verbose) cat("\n Zoom out full_1d")
isolate(ranges_homol$dmass <- FALSE)
}
refresh_homol$a<-(refresh_homol$a+1)
})
observeEvent(input$homol_counts_click, {
if(logfile$parameters$verbose) cat("\n Zoom out part_1_ad")
brush <- isolate(input$homol_counts_brush)
if (is.null(brush)) {
if(logfile$parameters$verbose) cat("\n Zoom out part_1_bd")
if(isolate(ranges_homol$dmass[1])!=FALSE){
old_range_dmass<-abs(isolate(ranges_homol$dmass[2]-ranges_homol$dmass[1]))
isolate(ranges_homol$dmass[1]<-ranges_homol$dmass[1]-.3*old_range_dmass)
isolate(ranges_homol$dmass[2]<-ranges_homol$dmass[2]+.3*old_range_dmass)
}
refresh_homol$a<-(refresh_homol$a+1)
}else{
if(logfile$parameters$verbose) cat("\n Doing hover_d")
isolate(ranges_homol$dmass <- c(brush$xmin, brush$xmax))
refresh_homol$b<-(refresh_homol$b+1)
}
})
observeEvent(input$homol_RT_dblclick, {
brush <- isolate(input$homol_RT_brush)
if (!is.null(brush)) {
if(logfile$parameters$verbose) cat("\n Zoom in_1d")
isolate(ranges_homol$dRT <- c(brush$xmin, brush$xmax))
} else {
if(logfile$parameters$verbose) cat("\n Zoom out full_1d")
isolate(ranges_homol$dRT <- FALSE)
}
refresh_homol$a<-(refresh_homol$a+1)
})
observeEvent(input$homol_RT_click, {
if(logfile$parameters$verbose) cat("\n Zoom out part_1_ad")
brush <- isolate(input$homol_RT_brush)
if (is.null(brush)) {
if(logfile$parameters$verbose) cat("\n Zoom out part_1_bd")
if(isolate(ranges_homol$dRT[1])!=FALSE){
old_range_dmass<-abs(isolate(ranges_homol$dRT[2]-ranges_homol$dRT[1]))
isolate(ranges_homol$dRT[1]<-ranges_homol$dRT[1]-.3*old_range_dmass)
isolate(ranges_homol$dRT[2]<-ranges_homol$dRT[2]+.3*old_range_dmass)
}
refresh_homol$a<-(refresh_homol$a+1)
}else{
if(logfile$parameters$verbose) cat("\n Doing hover_d")
isolate(ranges_homol$dRT <- c(brush$xmin, brush$xmax))
refresh_homol$b<-(refresh_homol$b+1)
}
})
observeEvent(input$homol_chromat_dblclick, {
brush <- isolate(input$homol_chromat_brush)
if (!is.null(brush)) {
if(logfile$parameters$verbose) cat("\n Zoom in_1e")
isolate(ranges_homol$RTchrom <- c(brush$xmin, brush$xmax))
isolate(ranges_homol$intchrom <- c(brush$ymin, brush$ymax))
} else {
if(logfile$parameters$verbose) cat("\n Zoom out full_1e")
isolate(ranges_homol$RTchrom <- FALSE)
isolate(ranges_homol$intchrom <- FALSE)
}
refresh_homol$d<-(refresh_homol$d+1)
})
observeEvent(input$homol_chromat_click, {
if(logfile$parameters$verbose) cat("\n Zoom out part_1_ae")
brush <- isolate(input$homol_chromat_brush)
if (is.null(brush)) {
if(logfile$parameters$verbose)  cat("\n Zoom out part_1_be")
if(isolate(ranges_homol$RTchrom[1])!=FALSE){
old_range_dmass<-abs(isolate(ranges_homol$RTchrom[2]-ranges_homol$RTchrom[1]))
isolate(ranges_homol$RTchrom[1]<-ranges_homol$RTchrom[1]-.3*old_range_dmass)
isolate(ranges_homol$RTchrom[2]<-ranges_homol$RTchrom[2]+.3*old_range_dmass)
old_range_dmass<-abs(isolate(ranges_homol$intchrom[2]-ranges_homol$intchrom[1]))
isolate(ranges_homol$intchrom[1]<-ranges_homol$intchrom[1]-.3*old_range_dmass)
isolate(if(ranges_homol$intchrom[1]<0){ranges_homol$intchrom[1]<-0})
isolate(ranges_homol$intchrom[2]<-ranges_homol$intchrom[2]+.3*old_range_dmass)
}
refresh_homol$d<-(refresh_homol$d+1)
}else{
if(logfile$parameters$verbose) cat("\n Doing hover_e - nothing")
}
})
observe({
s1<-input$comp_table_full_rows_selected
if(isolate(init$a)=="TRUE"){
if(length(s1)){
if(logfile$parameters$verbose) print(s1)
if(s1>=1){
updateNumericInput(session, inputId = "sel_meas_comp_comp", value = as.numeric(comp_table_full[s1,"Component ID |"]))
isolate(ranges_compo$mass<-FALSE)
isolate(ranges_compo$RTchrom<-FALSE)
isolate(ranges_compo$intchrom<-FALSE)
}
}
}
})
observe({
input$sel_meas_comp_peak
if(isolate(init$a)=="TRUE"){ if(logfile$parameters$verbose) cat("\n in Comp_B") }
if(!is.na(as.numeric(isolate(input$sel_meas_comp_peak)))){
if(isolate(init$a)=="TRUE" & as.numeric(isolate(input$sel_meas_comp_peak))>0){
if(
file.exists(file.path(logfile[[1]],"results","componentization","components",isolate(input$sel_meas_comp))) &
(isolate(input$sel_meas_comp_peak)>0)
){
those<-gsub("*", "", component[["Components"]][,"ID pattern peaks |"], fixed=TRUE)
that<-which(!is.na(unlist(lapply(strsplit(those,","), match, x=as.numeric(isolate(input$sel_meas_comp_peak))))))
if(length(that)==0){
those<-gsub("*", "", component[["Components"]][,"ID adduct peaks |"], fixed=TRUE)
that<-which(!is.na(unlist(lapply(strsplit(those,","), match, x=as.numeric(isolate(input$sel_meas_comp_peak))))))
}
if(length(that)==0){
those<-gsub("*", "", component[["Components"]][,"ID interfering peaks |"], fixed=TRUE)
that<-which(!is.na(unlist(lapply(strsplit(those,","), match, x=as.numeric(isolate(input$sel_meas_comp_peak))))))
}
if(length(that)==1){
if(logfile$parameters$verbose) cat("\n in Comp_B_1")
updateNumericInput(session,"sel_meas_comp_comp",value=that)
}else{
if(verbose){cat("\n in Comp_B_2")}
if(logfile$parameters$verbose) cat("\n Invalid peak selected!")
updateNumericInput(session,"sel_meas_comp_comp",value=0)
output$comp_plot_spec <- renderPlot({
plot.new()
plot.window(xlim=c(0,1),ylim=c(0,1))
text(.5,.5,labels="No components for this peak ID. \n -> peak either removed during replicate intersection or peak ID invalid.")
},res=110)
output$comp_plot_chromat <- renderPlot({
plot.new()
},res=110)
output$comp_plot_circ <- renderPlot({
plot.new()
},res=110)
output$comp_table_a <- DT::renderDataTable(
datatable(data.frame("No results"),colnames="No results")
)
output$comp_table_b <- DT::renderDataTable(
datatable(data.frame("No results"),colnames="No results")
)
output$comp_table_c  <- DT::renderDataTable(
datatable(data.frame("No results"),colnames="No results")
)
output$comp_table_d <- DT::renderDataTable(
datatable(data.frame("No results"),colnames="No results")
)
}
}
}
}
})
observe({
input$sel_meas_comp_comp
if(isolate(init$a)=="TRUE"){ if(logfile$parameters$verbose) cat("\n in Comp_C") }
if(!is.na(as.numeric(isolate(input$sel_meas_comp_comp)))){
if(isolate(init$a) == "TRUE" & as.numeric(isolate(input$sel_meas_comp_comp))){
if(
file.exists(file.path(logfile[[1]],"results","componentization","components",isolate(input$sel_meas_comp))) &
(isolate(input$sel_meas_comp_comp) > 0)
){
if(
any(component[[1]][,1] == as.numeric(isolate(input$sel_meas_comp_comp)))
){
if(logfile$parameters$verbose) cat("\n in Comp_C_1")
ee$entry <- as.numeric(isolate(input$sel_meas_comp_comp))
if(logfile$parameters$verbose) print(ee$entry)
}else{
if(logfile$parameters$verbose) cat("\n in Comp_C_2")
if(logfile$parameters$verbose) cat("\n Invalid component selected!")
updateNumericInput(session,"sel_meas_comp_comp",value=0)
updateNumericInput(session,"sel_meas_comp_peak",value=0)
output$found_compo<-renderText("Invalid peak ID")
output$comp_plot_spec <- renderPlot({
plot.new()
plot.window(xlim=c(0,1),ylim=c(0,1))
text(.5,.5,labels="No components for this peak, \n peak removed during replicate intersection, blind subtraction or \n invalid peak ID?.")
},res=110)
output$comp_plot_chromat <- renderPlot({
plot.new()
},res=110)
output$comp_plot_circ <- renderPlot({
plot.new()
},res=110)
}
}else{
if(logfile$parameters$verbose) cat("\n in Comp_C_3")
if(logfile$parameters$verbose) cat("\n Invalid component selected!")
output$found_compo<-renderText("Invalid peak ID")
output$comp_plot_spec <- renderPlot({
plot.new()
plot.window(xlim=c(0,1),ylim=c(0,1))
text(.5,.5,labels="No components for this peak, \n peak removed during replicate intersection, blind subtraction or \n invalid peak ID?.")
},res=110)
output$comp_plot_chromat <- renderPlot({
plot.new()
},res=110)
output$comp_plot_circ <- renderPlot({
plot.new()
},res=110)
}
}else{
if(isolate(init$a)=="TRUE"){ if(logfile$parameters$verbose) cat("\n Invalid component selected!")}
output$found_compo<-renderText("Invalid peak ID")
output$comp_plot_spec <- renderPlot({
plot.new()
plot.window(xlim=c(0,1),ylim=c(0,1))
text(.5,.5,labels="No components for this peak, \n peak removed during replicate intersection, blind subtraction or \n invalid peak ID?.")
},res=110)
output$comp_plot_chromat <- renderPlot({
plot.new()
},res=110)
output$comp_plot_circ <- renderPlot({
plot.new()
},res=110)
}
}else{
if(logfile$parameters$verbose) cat("\n in Comp_C_5")
}
})
observe({
ee$entry
refresh_compo$a
if(isolate(init$a)=="TRUE"){ if(logfile$parameters$verbose) print(ee$entry) }
if( isolate(ee$entry)>0 || isolate(refresh_compo$a)>0 ){
if(logfile$parameters$verbose) cat("\n in Comp_D_1")
got_comp<-enviMass::IXEME0353N(component, compoID=as.numeric(isolate(ee$entry)), what="check")
if(got_comp=="available"){
output$found_compo<-renderText("")
output$comp_plot_spec <- renderPlot({
par(mar=c(4,4,.9,.8));
enviMass::IXEME0353N(
component,
compoID=as.numeric(isolate(ee$entry)),
what="spec",
n_col=max_peak_ID,
masslim=isolate(ranges_compo$mass)
)
},res=110)
at_comp<<-match(as.numeric(isolate(ee$entry)),comp_table_full[,"Component ID |"])
peaks1<-comp_table_full[at_comp,"ID pattern peaks |"]
peaks1<-gsub("*","",peaks1,fixed=TRUE)
peaks1<-as.numeric(strsplit(peaks1,",")[[1]])
if(comp_table_full[at_comp,"ID adduct peaks |"]!="-"){
peaks2<-comp_table_full[at_comp,"ID adduct peaks |"]
peaks2<-gsub("*","",peaks2,fixed=TRUE)
peaks2<-as.numeric(strsplit(peaks2,",")[[1]])
}else{
peaks2<-c()
}
if(comp_table_full[at_comp,"ID interfering peaks |"]!="-"){
peaks3<-c()
}else{
peaks3<-c()
}
all_peaks<<-c(peaks1,peaks2,peaks3)
if(any(objects(envir=as.environment(".GlobalEnv"))=="MSlist")){
if(any(names(MSlist)=="File_ID")){
if(MSlist[["File_ID"]]!=as.character(isolate(input$sel_meas_comp))){
load(file.path(logfile[[1]],"MSlist",as.character(isolate(input$sel_meas_comp))),envir=as.environment(".GlobalEnv"))
}
}else{
load(file.path(logfile[[1]],"MSlist",as.character(isolate(input$sel_meas_comp))),envir=as.environment(".GlobalEnv"))
MSlist[["File_ID"]]<-as.character(isolate(input$sel_meas_comp))
}
}else{
load(file.path(logfile[[1]],"MSlist",as.character(isolate(input$sel_meas_comp))),envir=as.environment(".GlobalEnv"))
}
output$comp_plot_chromat <- renderPlot({
par(mar=c(4.5,4,1,.8));
enviMass:::JJGET6467P(
MSlist,
peakIDs = all_peaks,
RTlim = isolate(ranges_compo$RTchrom),
Intlim = isolate(ranges_compo$intchrom),
masslim = isolate(ranges_compo$mass),
n_col = max_peak_ID,
set_RT = (input$comp_chromat_time),
normalize = (input$comp_chromat_norm),
chromat_full = (input$comp_chromat_type)
);
},res=110)
output$comp_plot_circ <- renderPlot({
enviMass::IXEME0353N(component, compoID=as.numeric(isolate(ee$entry)), what="circ")
},res=110)
comp_table<<-enviMass::IXEME0353N(component, compoID=as.numeric(isolate(ee$entry)), what="table")
inser<-rep("",length(comp_table$relations[,2]))
do_these<-which(grepl("<->",comp_table$relations[,2]))
if(length(do_these)>0){
inser[do_these]<-paste("Same isotopologues of different adducts:",comp_table$relations[do_these,2],sep=" ")
}
do_these<-which(grepl(",z=",comp_table$relations[,2]))
if(length(do_these)>0){
inser[do_these]<-"Different isotopologues of the same adduct"
}
if(component[[1]][as.numeric(isolate(ee$entry)),"Target peaks"]!="-"){
output$which_comp_tar<-renderText(paste0(
"Target/suspect peaks found for this component: ",
component[[1]][as.numeric(isolate(ee$entry)),"Target peaks"]
))
}else{
output$which_comp_tar<-renderText("No target or suspect peaks found for this component.")
}
if(component[[1]][as.numeric(isolate(ee$entry)),"ISTD peaks"]!="-"){
output$which_comp_ISTD<-renderText(paste0(
"ISTD peaks found for this component: ",
component[[1]][as.numeric(isolate(ee$entry)),"ISTD peaks"]
))
}else{
output$which_comp_ISTD<-renderText("No ISTD peaks found for this component.")
}
output$comp_table_a <- DT::renderDataTable(
datatable(
cbind(
unlist(lapply(strsplit(comp_table$relations[,1],"-"), `[[`,1)),
unlist(lapply(strsplit(comp_table$relations[,1],"-"), `[[`,2)),
inser,
round(comp_table$relations[,3],digits=3)
),
colnames=c("First peaks","Second peak","Links","Intensity ratio")
)
)
output$comp_table_b <- DT::renderDataTable(
datatable(
cbind(
round(comp_table$a[,1],digits=1),
round(comp_table$a[,2],digits=5),
format(comp_table$a[,3],scientific=TRUE,digits=2),
round(comp_table$a[,4],digits=2),
comp_table$a[,6]
),
colnames=c("Peak ID","m/z","Intensity","RT [s]","In blind?")
)
)
if(comp_table[[4]]=="Not part of a homologue series"){
output$comp_table_c  <- DT::renderDataTable(
datatable(
data.frame("Peaks in selected component are not part of any homologue series"),
colnames="No results"
)
)
}else{
these<-unique(strsplit(comp_table[[4]],"/")[[1]][-1])
output$comp_table_c  <- DT::renderDataTable(
datatable(
cbind(
homol[[3]][these,1],
round(homol[[3]][these,3],digits=4),
round(homol[[3]][these,4],digits=1)
),
colnames=c("Series ID","m/z difference","RT difference [s]")
)
)
}
output$comp_table_d <- DT::renderDataTable(
datatable(
cbind(
round(comp_table[[2]][,1],digits=1),
round(comp_table[[2]][,2],digits=5),
format(comp_table[[2]][,3],scientific=TRUE,digits=2),
round(comp_table[[2]][,4],digits=2)
),
colnames=c("Peak ID","m/z","Intensity","RT [s]")
)
)
}else{
if(got_comp=="single_peak"){
output$found_compo<-renderText("single_peak")
output$comp_plot_spec <- renderPlot({
par(mar=c(4,4,.9,.8));
enviMass::IXEME0353N(
component,
compoID=as.numeric(isolate(ee$entry)),
what="spec",
n_col=max_peak_ID,
masslim=isolate(ranges_compo$mass))
},res=110)
at_comp<<-match(as.numeric(isolate(ee$entry)),comp_table_full[,"Component ID |"])
peaks1<-comp_table_full[at_comp,"ID pattern peaks |"]
peaks1<-gsub("*","",peaks1,fixed=TRUE)
peaks1<-as.numeric(strsplit(peaks1,",")[[1]])
if(comp_table_full[at_comp,"ID adduct peaks |"]!="-"){
peaks2<-comp_table_full[at_comp,"ID adduct peaks |"]
peaks2<-gsub("*","",peaks2,fixed=TRUE)
peaks2<-as.numeric(strsplit(peaks2,",")[[1]])
}else{
peaks2<-c()
}
if(comp_table_full[at_comp,"ID interfering peaks |"]!="-"){
peaks3<-c()
}else{
peaks3<-c()
}
all_peaks<<-c(peaks1,peaks2,peaks3)
if(any(objects(envir=as.environment(".GlobalEnv"))=="MSlist")){
if(any(names(MSlist)=="File_ID")){
if(MSlist[["File_ID"]]!=as.character(isolate(input$sel_meas_comp))){
load(file.path(logfile[[1]],"MSlist",as.character(isolate(input$sel_meas_comp))),envir=as.environment(".GlobalEnv"))
}
}else{
load(file.path(logfile[[1]],"MSlist",as.character(isolate(input$sel_meas_comp))),envir=as.environment(".GlobalEnv"))
MSlist[["File_ID"]]<-as.character(isolate(input$sel_meas_comp))
}
}else{
load(file.path(logfile[[1]],"MSlist",as.character(isolate(input$sel_meas_comp))),envir=as.environment(".GlobalEnv"))
}
output$comp_plot_chromat <- renderPlot({
par(mar=c(5,4,1,.8))
enviMass:::JJGET6467P(
MSlist,
peakIDs=all_peaks,
RTlim=isolate(ranges_compo$RTchrom),
Intlim=isolate(ranges_compo$intchrom),
masslim=isolate(ranges_compo$mass),
n_col=max_peak_ID,
set_RT=(input$comp_chromat_time),
normalize=(input$comp_chromat_norm),
chromat_full=(input$comp_chromat_type)
);
},res=110)
comp_table<<-enviMass::IXEME0353N(component, compoID=as.numeric(isolate(ee$entry)), what="table")
output$comp_plot_circ <- renderPlot({
plot.new()
plot.window(xlim=c(0,1),ylim=c(0,1))
text(.5,.5,labels="The selected component contains only one peak.")
},res=110)
if(component[[1]][as.numeric(isolate(ee$entry)),"Target peaks"]!="-"){
output$which_comp_tar<-renderText(paste0(
"Target/suspect peaks found for this component: ",
component[[1]][as.numeric(isolate(ee$entry)),"Target peaks"]
))
}else{
output$which_comp_tar<-renderText("No target or suspect peaks found for this component.")
}
if(component[[1]][as.numeric(isolate(ee$entry)),"ISTD peaks"]!="-"){
output$which_comp_ISTD<-renderText(paste0(
"ISTD peaks found for this component: ",
component[[1]][as.numeric(isolate(ee$entry)),"ISTD peaks"]
))
}else{
output$which_comp_ISTD<-renderText("No ISTD peaks found for this component.")
}
output$comp_table_a <- DT::renderDataTable(
datatable(data.frame("No results"),colnames="No results")
)
output$comp_table_b <- DT::renderDataTable(
datatable(
cbind(
round(comp_table$a[,1],digits=1),
round(comp_table$a[,2],digits=5),
format(comp_table$a[,3],scientific=TRUE,digits=2),
round(comp_table$a[,4],digits=2),
comp_table$a[,6]
),
colnames=c("Peak ID","m/z","Intensity","RT [s]","In blind?")
)
)
if(comp_table[[4]]=="Not part of a homologue series"){
output$comp_table_c  <- DT::renderDataTable(
datatable(
data.frame("Peaks in selected component are not part of any homologue series"),
colnames="No results"
)
)
}else{
these<-unique(strsplit(comp_table[[4]],"/")[[1]][-1])
output$comp_table_c  <- DT::renderDataTable(
datatable(
cbind(
homol[[3]][these,1],
round(homol[[3]][these,3],digits=4),
round(homol[[3]][these,4],digits=1)
),
colnames=c("Series ID","m/z difference","RT difference [s]")
)
)
}
output$comp_table_d <- DT::renderDataTable(
datatable(
cbind(
round(comp_table[[2]][,1],digits=1),
round(comp_table[[2]][,2],digits=5),
format(comp_table[[2]][,3],scientific=TRUE,digits=2),
round(comp_table[[2]][,4],digits=2)
),
colnames=c("Peak ID","m/z","Intensity","RT [s]")
)
)
}else{
output$found_compo<-renderText("Invalid peak ID")
output$comp_plot_spec <- renderPlot({
plot.new()
plot.window(xlim=c(0,1),ylim=c(0,1))
text(.5,.5,labels="No components for this peak, \n peak removed during replicate intersection, blind subtraction or \n invalid peak ID?.")
},res=110)
output$comp_plot_chromat <- renderPlot({
plot.new()
},res=110)
output$comp_plot_circ <- renderPlot({
plot.new()
},res=110)
output$comp_table_a <- DT::renderDataTable(
datatable(data.frame("No results"),colnames="No results")
)
output$comp_table_b <- DT::renderDataTable(
datatable(data.frame("No results"),colnames="No results")
)
output$comp_table_c  <- DT::renderDataTable(
datatable(data.frame("No results"),colnames="No results")
)
output$comp_table_d <- DT::renderDataTable(
datatable(data.frame("No results"),colnames="No results")
)
}
}
}
})
observe({
input$comp_chromat_time
if(isolate(init$a)=="TRUE" & isolate(ranges_compo$RTchrom[1]!=FALSE)){
if(isolate(input$comp_chromat_time)=="minutes"){
isolate(ranges_compo$RTchrom<-(ranges_compo$RTchrom/60))
}
if(isolate(input$comp_chromat_time)=="seconds"){
isolate(ranges_compo$RTchrom<-(ranges_compo$RTchrom*60))
}
}
})
observe({
input$comp_chromat_norm
if(isolate(init$a)=="TRUE" & isolate(ranges_compo$intchrom[1]!=FALSE)){
isolate(ranges_compo$intchrom<-FALSE)
}
})
observeEvent(input$comp_plot_chromat_dblclick, {
brush <- isolate(input$comp_plot_chromat_brush)
if (!is.null(brush)) {
if(logfile$parameters$verbose) cat("\n Zoom comp in_1e")
isolate(ranges_compo$RTchrom <- c(brush$xmin, brush$xmax))
isolate(ranges_compo$intchrom <- c(brush$ymin, brush$ymax))
} else {
if(logfile$parameters$verbose) cat("\n Zoom out comp  full_1e")
isolate(ranges_compo$RTchrom <- FALSE)
isolate(ranges_compo$intchrom <- FALSE)
}
isolate(refresh_compo$a<-(refresh_compo$a+1))
})
observeEvent(input$comp_plot_chromat_click, {
cat("\n Zoom out comp  part_1_ae")
brush <- isolate(input$comp_plot_chromat_brush)
if (is.null(brush)) {
if(logfile$parameters$verbose) cat("\n Zoom out comp  part_1_be")
if(isolate(ranges_compo$RTchrom[1])!=FALSE){
old_range_dmass<-abs(isolate(ranges_compo$RTchrom[2]-ranges_compo$RTchrom[1]))
isolate(ranges_compo$RTchrom[1]<-ranges_compo$RTchrom[1]-.3*old_range_dmass)
isolate(ranges_compo$RTchrom[2]<-ranges_compo$RTchrom[2]+.3*old_range_dmass)
old_range_dmass<-abs(isolate(ranges_compo$intchrom[2]-ranges_compo$intchrom[1]))
isolate(ranges_compo$intchrom[1]<-ranges_compo$intchrom[1]-.3*old_range_dmass)
isolate(if(ranges_compo$intchrom[1]<0){ranges_compo$intchrom[1]<-0})
isolate(ranges_compo$intchrom[2]<-ranges_compo$intchrom[2]+.3*old_range_dmass)
}
isolate(refresh_compo$a<-(refresh_compo$a+1))
}else{
cat("\n Doing hover_e comp comp - nothing")
}
})
observeEvent(input$comp_plot_spec_dblclick, {
brush <- isolate(input$comp_plot_spec_brush)
if (!is.null(brush)) {
if(logfile$parameters$verbose) cat("\n Zoom comp in_1e")
isolate(ranges_compo$mass <- c(brush$xmin, brush$xmax))
} else {
if(logfile$parameters$verbose) cat("\n Zoom out comp  full_1e")
isolate(ranges_compo$mass <- FALSE)
}
isolate(refresh_compo$a<-(refresh_compo$a+1))
})
observeEvent(input$comp_plot_spec_click, {
if(logfile$parameters$verbose) cat("\n Zoom out comp  part_1_ae")
brush <- isolate(input$comp_plot_spec_brush)
if (is.null(brush)) {
if(logfile$parameters$verbose) cat("\n Zoom out comp  part_1_be")
if(isolate(ranges_compo$mass[1])!=FALSE){
old_range_dmass<-abs(isolate(ranges_compo$mass[2]-ranges_compo$mass[1]))
isolate(ranges_compo$mass[1]<-ranges_compo$mass[1]-.3*old_range_dmass)
isolate(ranges_compo$mass[2]<-ranges_compo$mass[2]+.3*old_range_dmass)
}
isolate(refresh_compo$a<-(refresh_compo$a+1))
}else{
if(logfile$parameters$verbose) cat("\n Doing hover_e comp comp - nothing")
}
})
output$atom_bounds_that <- renderUI({
if(length(input$atom_bounds_this)){
if(logfile$parameters$verbose) cat("\n Am observing")
those_elements <- input$atom_bounds_this
lapply(those_elements, function(i){
numericInput(inputId = paste0("ppm_", i),label = paste0("Maximum shifts for ", i),value = 30, width='200px')
})
}else{
if(logfile$parameters$verbose) cat("\n Am observing nothing")
}
})
observe({
input$atom_bound_peak
input$sel_meas_comp
input$atom_bound_addpeaks
input$atom_bounds_calculate
if(isolate(init$a)=="TRUE"){
if(file.exists(file.path(logfile[[1]],"peaklist",as.character(isolate(input$sel_meas_comp))))){
found_peaklist<-TRUE
}else{
found_peaklist<-FALSE
}
if(!is.na(isolate(input$atom_bound_peak))>0){
if(
(isolate(input$atom_bound_peak)>0) &
(isolate(input$sel_meas_comp)>0) &
found_peaklist
){
if(logfile$parameters$verbose) cat("\n in Atoms_1")
if(isolate(input$atom_bound_addpeaks)=="(b) all peaks with similar RT"){
if(logfile$parameters$verbose) cat("\n in Atoms_2")
load(file.path(logfile[[1]],"peaklist",as.character(isolate(input$sel_meas_comp))),envir=as.environment(".GlobalEnv"))
at_peak<<-which(peaklist[,"peak_ID"]==as.numeric(isolate(input$atom_bound_peak)))
if(length(at_peak)>0){
atom_peaks<<-peaklist[
(peaklist[,"m/z_corr"]>peaklist[at_peak,"m/z_corr"]) &
((abs(peaklist[,"RT_corr"]-peaklist[at_peak,"RT_corr"]))<=logfile$parameters$isotop_rttol) &
(peaklist[,"keep"]==1) & # omit replicate peaks
(abs(peaklist[,"m/z_corr"]-peaklist[at_peak,"m/z_corr"])<=10)
,,drop=FALSE]
atom_peaks<<-rbind(
peaklist[at_peak,],
atom_peaks)
}
}else{
if(logfile$parameters$verbose) cat("\n in Atoms_3")
at_peak<<-which(!is.na(unlist(lapply(strsplit(gsub("*","",component[["Components"]][,"ID pattern peaks |"],fixed=TRUE),","), match, x=as.numeric(isolate(input$atom_bound_peak))))))
if(length(at_peak)==0){
at_peak<<-which(!is.na(unlist(lapply(strsplit(gsub("*","",component[["Components"]][,"ID pattern peaks |"],fixed=TRUE),","), match, x=as.numeric(isolate(input$atom_bound_peak))))))
}
if(length(at_peak)==0){
at_peak<<-which(!is.na(unlist(lapply(strsplit(gsub("*","",component[["Components"]][,"ID pattern peaks |"],fixed=TRUE),","), match, x=as.numeric(isolate(input$atom_bound_peak))))))
}
if(length(at_peak)>0){
if(verbose){cat("\n in Atoms_4")}
get_peaks<<-c()
if(component[["Components"]][at_peak,"ID pattern peaks |"]!="-"){
get_peaks<<-c(get_peaks,strsplit(component[["Components"]][at_peak,"ID pattern peaks |"],",")[[1]])
}
if(component[["Components"]][at_peak,"ID adduct peaks |"]!="-"){
get_peaks<<-c(get_peaks,strsplit(component[["Components"]][at_peak,"ID adduct peaks |"],",")[[1]])
}
if(component[["Components"]][at_peak,"ID interfering peaks |"]!="-"){
get_peaks<<-c(get_peaks,strsplit(component[["Components"]][at_peak,"ID interfering peaks |"],",")[[1]])
}
get_peaks<<-sapply(get_peaks,gsub,pattern="*",replacement="",fixed = TRUE, USE.NAMES = FALSE)
get_peaks<<-as.numeric(get_peaks)
matched<-match(get_peaks,component[["pattern peak list"]][,"peak ID"])
peaklist<<-component[["pattern peak list"]][matched,1:4,drop=FALSE]
at_peak<<-which(peaklist[,"peak ID"]==isolate(input$atom_bound_peak))
atom_peaks<<-peaklist[
(peaklist[,"m/z_corr"]>peaklist[at_peak,"m/z_corr"]) &
((abs(peaklist[,"RT_corr"]-peaklist[at_peak,"RT_corr"]))<=logfile$parameters$isotop_rttol) &
(abs(peaklist[,"m/z_corr"]-peaklist[at_peak,"m/z_corr"])<=10)
,1:3,drop=FALSE]
atom_peaks<<-rbind(
peaklist[at_peak,1:3],
atom_peaks)
}else{
if(logfile$parameters$verbose) cat("\n Invalid peak selected - not found among components.")
}
}
if(length(at_peak)>0){
if(
(logfile$workflow[names(logfile$workflow)=="LOD"]=="yes") &
(file.exists(file=file.path(logfile$project_folder,"results","LOD","LOD_splined")))
){
if(logfile$parameters$verbose) cat("\n in Atoms_4")
load(file=file.path(logfile$project_folder,"results","LOD","LOD_splined"));
with_model<-which(names(LOD_splined)==paste("LOD_",as.character(isolate(input$sel_meas_comp)),sep=""))
if(length(with_model)>0){
use_LOD<<-10^(predict(LOD_splined[[with_model]],atom_peaks[1,"RT_corr"])$y[[1]])
}else{
cat("\n Shouldn`t there be a LOD spline? Could not find it! - using Lower intensity threshold set for target screening!")
use_LOD<<-as.numeric(logfile$parameters$tar_intcut)
}
}else{
if(logfile$parameters$verbose) cat("\n in Atoms_5")
if(logfile$parameters$verbose) cat("\n No LOD interpolation in workflow included - using Lower intensity threshold set for target screening!")
use_LOD<<-as.numeric(logfile$parameters$tar_intcut)
}
output$atom_bound_plot_peak <- renderPlot({
plot(
atom_peaks[,"m/z_corr"],atom_peaks[,"int_corr"],
xlab="m/z",ylab="Intensity",type="h",col="darkgreen",lwd=2,
ylim=c(0,max(atom_peaks[,"int_corr"])))
points(
atom_peaks[1,"m/z_corr"],atom_peaks[1,"int_corr"],
type="h",col="red",lwd=2)
abline(h=use_LOD,col="gray")
},res=110)
elements<<-isolate(input$atom_bounds_this)
charges<<-c(1,2,3,4)
if(length(elements)>0){
dmz<<-c()
lapply(elements,function(i){
x<-paste0("ppm_",i)
dmz<<-c(dmz,
input[[x]]
)
})
atom_counts<<-try({
enviMass::QZTGM3353W(
masses=atom_peaks[,"m/z_corr"],
intensities=atom_peaks[,"int_corr"],
elements,
dmz,
ppm=TRUE,
charges,
isotopes,
int_cut=use_LOD,
inttol=0.2,
use_C=as.logical(isolate(input$atom_bound_wcarbon)),
must_peak=FALSE
)
})
}
if(class(atom_counts)!="try-error"){
if(verbose){cat("\n in Atoms_6")}
output$atom_count_table <- DT::renderDataTable(
datatable(as.data.frame(cbind(charges,atom_counts),row.names = NULL),
colnames=c("Charge z",elements),rownames=FALSE)
)
}else{
output$atom_count_table <- DT::renderDataTable(
datatable(data.frame("No results- sth went wrong - debug?"),colnames="No results")
)
cat("\n Atom count estimation failed - debug!")
}
}else{
output$atom_bound_plot_peak <- renderPlot({
plot.new()
plot.window(xlim=c(0,1),ylim=c(0,1))
text(0.5,0.5,labels="Invalid peak ID")
},res=110)
}
}else{
output$atom_bound_plot_peak <- renderPlot({
plot.new()
plot.window(xlim=c(0,1),ylim=c(0,1))
text(0.5,0.5,labels="No peaks available")
},res=110)
}
}
}
})
