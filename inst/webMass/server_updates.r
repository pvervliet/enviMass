if(any(ls()=="logfile")){stop("\n illegal logfile detected _1 in server_updates.r!")}
#stop("\n\nMaintenance work; enviMass will be back in a couple of hours! Please update enviMass again later.")
if(!any(names(resolution_list)==logfile$parameters$resolution)){
shinyjs::info(paste0("Please specifiy your Instrument/Resolution for your instrument ",logfile$parameters$resolution," again (Settings tab): such specifications have changed and had to be reset."));
logfile$parameters$resolution<<-"Elite_R240000@400";
save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(as.numeric(installed.packages()[installed.packages()[,1]=="enviPat","Version"])<2.2){
shinyjs::info("Please first update enviPat (an enviMass package dependency)  to a version >=2.2. Check www.enviMass.ch installation section for how to update all dependencies. Aborting enviMass ...");
stop("\n package enviPat update required! Abort ...")
}
if(logfile$version<3.100){
cat("\n Updating to version 3.100 ...")
if(!file.exists(file.path(logfile$project_folder,"results","screening"))){
dir.create(file.path(logfile$project_folder,"results","screening"),recursive=TRUE)
}
if(!file.exists(file.path(logfile$project_folder,"quantification"))){
dir.create(file.path(logfile$project_folder,"quantification"),recursive=TRUE)   # subfolder
}
if(!file.exists(file.path(logfile$project_folder,"results","LOD"))){
dir.create(file.path(logfile$project_folder,"results","LOD"),recursive=TRUE)
}
if(!file.exists(file.path(logfile$project_folder,"results","recalibration"))){
dir.create(file.path(logfile$project_folder,"results","recalibration"),recursive=TRUE)
}
# another column in peaklists for the replicates!
IDs<-list.files(file.path(logfile[[1]],"peaklist"))
if(length(IDs)>0){
for(i in 1:length(IDs)){
load(file=file.path(logfile[[1]],"peaklist",as.character(IDs[i])),envir=as.environment(".GlobalEnv"),verbose=FALSE);
if(any(colnames(peaklist)=="keep")){next}
keep<-rep(1,length(peaklist[,1]))
peaklist<-cbind(peaklist,keep)
colnames(peaklist)[15]<-"keep";
save(peaklist,file=file.path(logfile[[1]],"peaklist",as.character(IDs[i])))
rm(peaklist)
}
}
first_ext<-FALSE;
if(!any(logfile$summary[,1]=="replicates")){
logfile$summary[,1]<<-as.character(logfile$summary[,1])
logfile$summary[,2]<<-as.character(logfile$summary[,2])
logfile$summary[11,1]<<-"replicates"
logfile$summary[11,2]<<-"FALSE"
first_ext<-TRUE;
save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(!any(logfile$summary[,1]=="IS_screen")){
logfile$summary[,1]<<-as.character(logfile$summary[,1])
logfile$summary[,2]<<-as.character(logfile$summary[,2])
logfile$summary[12,1]<<-"IS_screen"
logfile$summary[12,2]<<-"FALSE"
first_ext<-TRUE;
save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(!any(logfile$summary[,1]=="target_screen")){
logfile$summary[,1]<<-as.character(logfile$summary[,1])
logfile$summary[,2]<<-as.character(logfile$summary[,2])
logfile$summary[13,1]<<-"target_screen"
logfile$summary[13,2]<<-"FALSE"
first_ext<-TRUE;
save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(!any(logfile$summary[,1]=="LOD")){
logfile$summary[,1]<<-as.character(logfile$summary[,1])
logfile$summary[,2]<<-as.character(logfile$summary[,2])
logfile$summary[14,1]<<-"LOD"
logfile$summary[14,2]<<-"FALSE"
first_ext<-TRUE;
save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(!any(logfile$summary[,1]=="quantification")){
logfile$summary[,1]<<-as.character(logfile$summary[,1])
logfile$summary[,2]<<-as.character(logfile$summary[,2])
logfile$summary[15,1]<<-"quantification"
logfile$summary[15,2]<<-"FALSE"
first_ext<-TRUE;
save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(!any(logfile$summary[,1]=="blinds")){
logfile$summary[,1]<<-as.character(logfile$summary[,1])
logfile$summary[,2]<<-as.character(logfile$summary[,2])
logfile$summary[16,1]<<-"blinds"
logfile$summary[16,2]<<-"FALSE"
logfile$parameters$blind_dmz<<-"3";
logfile$parameters$blind_ppm<<-"TRUE";
logfile$parameters$blind_drt<<-"30";
first_ext<-TRUE;
save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(!any(logfile$summary[,1]=="IS_subtr")){
logfile$summary[17,1]<<-"IS_subtr"
logfile$summary[17,2]<<-"FALSE"
first_ext<-TRUE;
save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(!any(logfile$summary[,1]=="target_subtr")){
logfile$summary<<-rbind(logfile$summary,c("",""))
logfile$summary[18,1]<<-"target_subtr"
logfile$summary[18,2]<<-"FALSE"
first_ext<-TRUE;
save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(!any(logfile$summary[,1]=="blind_subtr")){
logfile$summary<<-rbind(logfile$summary,c("",""))
logfile$summary[19,1]<<-"blind_subtr"
logfile$summary[19,2]<<-"FALSE"
first_ext<-TRUE;
save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(first_ext){
logfile$summary[,1]<<-c(
"Data available?",
"peakpicking",
"qc",
"pattern",
"recal",
"align",
"norm",
"profiling",
"IS_normaliz",
"trendblind",
"replicates",
"IS_screen",
"target_screen",
"LOD",
"quantification",
"blinds",
"IS_subtr",
"target_subtr",
"blind_subtr"
)
}
if(!any(names(logfile$parameters)=="replicate_dmz")){
logfile$parameters$replicate_dmz<<-"3";
logfile$parameters$replicate_ppm<<-"TRUE";
logfile$parameters$replicate_recalib<<-"FALSE";
logfile$parameters$replicate_delRT<<-"30";
save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(!any(names(logfile$parameters)=="screen_IS_cutit")){
logfile$parameters$screen_IS_cutit<<-"FALSE";
logfile$parameters$screen_target_cutit<<-"FALSE";
save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(!any(names(logfile$parameters)=="peak_perc_cut")){
logfile$parameters$peak_perc_cut<<-"0";
save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(length(logfile[[2]])<23){
logfile[[2]]<<-rep("TRUE",23)
}
names(logfile[[2]])<<-c(
"peakpicking","qc","recal","norm","align","profiling","trendblind","pattern",
"replicates","IS_screen","target_screen","LOD","calibration","recovery","quantification","blinds","IS_normaliz","IS_subtr","target_subtr",
"blind_subtr","isotopologues","adducts","homologues"
)
names(logfile)[2]<<-c("Tasks_to_redo");
save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
if(!any(names(logfile$workflow)=="IS_screen")){
logfile$workflow[11]<<-"yes"; 	names(logfile$workflow)[11]<<-"IS_screen"
}
if(!any(names(logfile$workflow)=="target_screen")){
logfile$workflow[12]<<-"yes"; 	names(logfile$workflow)[12]<<-"target_screen"
}
if(!any(names(logfile$workflow)=="replicates")){
logfile$workflow[13]<<-"yes"; 	names(logfile$workflow)[13]<<-"replicates"
}
if(!any(names(logfile$workflow)=="LOD")){
logfile$workflow[6]<<-"yes"; 	names(logfile$workflow)[6]<<-"LOD"
}
if(!any(names(logfile$workflow)=="quantification")){
logfile$workflow[8]<<-"yes"; 	names(logfile$workflow)[8]<<-"quantification"
}
if(!any(names(logfile$workflow)=="blinds")){
logfile$workflow[14]<<-"yes"; 	names(logfile$workflow)[14]<<-"blinds"
}
if(any(names(logfile$workflow)=="profnorm")){
logfile$workflow[15]<<-"yes"; 	names(logfile$workflow)[15]<<-"IS_normaliz"
}
if(any(names(logfile$workflow)=="profiled")){
logfile$workflow[15]<<-"yes"; 	names(logfile$workflow)[9]<<-"profiling"
}
logfile$workflow[5]<<-"yes"; 	names(logfile$workflow)[5]<<-"pattern"
logfile$workflow[7]<<-"yes"; 	names(logfile$workflow)[7]<<-"peakpicking"
names(logfile$workflow)[10]<<-"trendblind"
if(!any(names(logfile$workflow)=="IS_subtr")){
logfile$workflow[16]<<-"yes"; 	names(logfile$workflow)[16]<<-"IS_subtr"
}
if(!any(names(logfile$workflow)=="target_subtr")){
logfile$workflow[17]<<-"yes"; 	names(logfile$workflow)[17]<<-"target_subtr"
}
if(!any(names(logfile$workflow)=="blind_subtr")){
logfile$workflow[18]<<-"yes"; 	names(logfile$workflow)[18]<<-"blind_subtr"
}
if(TRUE){
redo_pattern<-FALSE
if( file.exists(file.path(logfile[[1]],"results","pattern_pos_IS")) ){
load(file=file.path(logfile[[1]],"results","pattern_pos_IS"),envir=as.environment(".GlobalEnv"));
for(i in 1:length(pattern_pos_IS)){
if( any(pattern_pos_IS[[i]][,2]==0) ){
redo_pattern<-TRUE
}
}
rm(pattern_pos_IS,envir=as.environment(".GlobalEnv"))
}
if( file.exists(file.path(logfile[[1]],"results","pattern_neg_IS")) ){
load(file=file.path(logfile[[1]],"results","pattern_neg_IS"),envir=as.environment(".GlobalEnv"));
for(i in 1:length(pattern_neg_IS)){
if( any(pattern_neg_IS[[i]][,2]==0) ){
redo_pattern<-TRUE
}
}
rm(pattern_neg_IS,envir=as.environment(".GlobalEnv"))
}
if( file.exists(file.path(logfile[[1]],"results","pattern_pos_target")) ){
load(file=file.path(logfile[[1]],"results","pattern_pos_target"),envir=as.environment(".GlobalEnv"));
for(i in 1:length(pattern_pos_target)){
if( any(pattern_pos_target[[i]][,2]==0) ){
redo_pattern<-TRUE
}
}
rm(pattern_pos_target,envir=as.environment(".GlobalEnv"))
}
if( file.exists(file.path(logfile[[1]],"results","pattern_neg_target")) ){
load(file=file.path(logfile[[1]],"results","pattern_neg_target"),envir=as.environment(".GlobalEnv"));
for(i in 1:length(pattern_neg_target)){
if( any(pattern_neg_target[[i]][,2]==0) ){
redo_pattern<-TRUE
}
}
rm(pattern_neg_target,envir=as.environment(".GlobalEnv"))
}
if(redo_pattern){
logfile$workflow[names(logfile$workflow)=="pattern"]<<-"yes";
logfile$summary[logfile$summary[,1]=="pattern",2]<<-"TRUE";
logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)=="pattern"]<<-TRUE
}
}
work_names<-names(logfile$Tasks_to_redo)[1:23]
depend<-matrix(ncol=length(work_names),nrow=length(work_names),0)
colnames(depend)<-work_names
rownames(depend)<-work_names					# peakpicking	qc	recal	norm	align	profiling	trendblind	pattern		replicates	IS_screen	target_screeen	LOD		calibration	recovery	quantification	blinds IS_normaliz	IS_subtr	target_subtr	blind_subtr	isotopologues	adducts homologues
depend[,colnames(depend)=="peakpicking"]<-		c(0,			1,	1,		1,		1,		1,			1,			0,			1,			1,			1,				1,		1,			0,			1,				1,		1,			1,			1,				1,			0,				0,		0)
depend[,colnames(depend)=="qc"]<-				c(0,			0,	1,		1,		1,		1,			1,			0,			1,			1,			1,				1,		1,			0,			1,				1,		1,			1,			1,				1,			0,				0,		0)
depend[,colnames(depend)=="pattern"]<-			c(0,			0,	1,		0,		0,		0,			0,			0,			0,			1,			1,				0,		1,			0,			1,				0,		1,			1,			1,				0,			0,				0,		0)
depend[,colnames(depend)=="recal"]<-			c(0,			0,	0,		0,		0,		1,			1,			0,			1,			1,			1,				0,		1,			0,			1,				0,		1,			1,			1,				1,			0,				0,		0)
depend[,colnames(depend)=="align"]<-			c(0,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		1,			0,			1,				0,		0,			0,			0,				1,			0,				0,		0)
depend[,colnames(depend)=="norm"]<-				c(0,			0,	0,		0,		0,		1,			1,			0,			0,			1,			1,				1,		1,			0,			1,				1,		1,			0,			0,				1,			0,				0,		0)
depend[,colnames(depend)=="blinds"]<-			c(0,			0,	0,		0,		0,		1,			1,			0,			0,			1,			1,				0,		1,			0,			1,				0,		1,			1,			1,				1,			0,				0,		0)
depend[,colnames(depend)=="replicates"]<-		c(0,			0,	0,		0,		0,		1,			1,			0,			0,			1,			1,				1,		1,			0,			1,				0,		1,			1,			1,				1,			0,				0,		0)
depend[,colnames(depend)=="profiling"]<-		c(0,			0,	0,		0,		0,		0,			1,			0,			0,			1,			1,				0,		0,			0,			1,				0,		1,			1,			1,				1,			0,				0,		0)
depend[,colnames(depend)=="IS_screen"]<-		c(0,			0,	0,		0,		0,		1,			0,			0,			0,			0,			0,				0,		0,			0,			1,				0,		1,			1,			0,				0,			0,				0,		0)
depend[,colnames(depend)=="target_screen"]<-	c(0,			0,	0,		0,		0,		1,			0,			0,			0,			0,			0,				0,		0,			0,			1,				0,		1,			0,			1,				0,			0,				0,		0)
depend[,colnames(depend)=="IS_normaliz"]<-		c(0,			0,	0,		0,		0,		0,			1,			0,			0,			0,			0,				0,		0,			0,			0,				0,		1,			0,			0,				0,			0,				0,		0)
depend[,colnames(depend)=="trendblind"]<-		c(0,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
depend[,colnames(depend)=="LOD"]<-				c(0,			0,	0,		0,		0,		0,			0,			0,			0,			1,			1,				0,		1,			0,			0,				0,		1,			1,			1,				0,			0,				0,		0)
depend[,colnames(depend)=="calibration"]<-		c(0,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
depend[,colnames(depend)=="recovery"]<-			c(0,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,			1,				0,		0,			0,			0,				0,			0,				0,		0)
depend[,colnames(depend)=="quantification"]<-	c(0,			0,	0,		0,		0,		0,			0,			0,			0,			0,			1,				0,		0,			0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
depend[,colnames(depend)=="IS_subtr"]<-			c(0,			0,	0,		0,		0,		1,			1,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,			0,				1,			0,				0,		0)
depend[,colnames(depend)=="target_subtr"]<-		c(0,			0,	0,		0,		0,		1,			1,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,			0,				1,			0,				0,		0)
depend[,colnames(depend)=="blind_subtr"]<-		c(0,			0,	0,		0,		0,		1,			1,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
depend[,colnames(depend)=="isotopologues"]<-	c(0,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
depend[,colnames(depend)=="adducts"]<-			c(0,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
depend[,colnames(depend)=="homologues"]<-		c(0,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
logfile[[11]]<<-depend
names(logfile)[11]<<-"workflow_depend"
must<-matrix(ncol=length(work_names),nrow=length(work_names),0)
colnames(must)<-work_names
rownames(must)<-work_names					# peakpicking	qc	recal	norm	align	profiling	trendblind	pattern		replicates	IS_screen	target_screeen	LOD		calibration	recovery	quantification	blinds IS_normaliz	IS_subtr	target_subtr	blind_subtr	isotopologues	adducts homologues
must[,colnames(must)=="peakpicking"]<-		c(0,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
must[,colnames(must)=="qc"]<-				c(1,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
must[,colnames(must)=="pattern"]<-			c(0,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
must[,colnames(must)=="recal"]<-			c(1,			0,	0,		0,		0,		0,			0,			1,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
must[,colnames(must)=="align"]<-			c(1,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
must[,colnames(must)=="norm"]<-				c(1,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
must[,colnames(must)=="blinds"]<-			c(1,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
must[,colnames(must)=="replicates"]<-		c(1,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
must[,colnames(must)=="profiling"]<-		c(1,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
must[,colnames(must)=="IS_screen"]<-		c(1,			0,	0,		0,		0,		0,			0,			1,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
must[,colnames(must)=="target_screen"]<-	c(1,			0,	0,		0,		0,		0,			0,			1,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
must[,colnames(must)=="IS_normaliz"]<-		c(1,			0,	0,		0,		0,		0,			0,			1,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
must[,colnames(must)=="trendblind"]<-		c(1,			0,	0,		0,		0,		1,			0,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
must[,colnames(must)=="LOD"]<-				c(1,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
must[,colnames(must)=="calibration"]<-		c(1,			0,	0,		0,		0,		0,			0,			1,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
must[,colnames(must)=="recovery"]<-			c(1,			0,	0,		0,		0,		0,			0,			1,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
must[,colnames(must)=="quantification"]<-	c(1,			0,	0,		0,		0,		0,			0,			1,			0,			1,			1,				0,		0,			0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
must[,colnames(must)=="IS_subtr"]<-			c(1,			0,	0,		0,		0,		1,			0,			1,			0,			1,			1,				0,		0,			0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
must[,colnames(must)=="target_subtr"]<-		c(1,			0,	0,		0,		0,		1,			0,			1,			0,			1,			1,				0,		0,			0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
must[,colnames(must)=="blind_subtr"]<-		c(1,			0,	0,		0,		0,		1,			0,			0,			0,			0,			0,				0,		0,			0,			0,				1,		0,			0,			0,				0,			0,				0,		0)
must[,colnames(must)=="isotopologues"]<-	c(1,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
must[,colnames(must)=="adducts"]<-			c(1,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
must[,colnames(must)=="homologues"]<-		c(1,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,			0,				0,			0,				0,		0)
logfile[[12]]<<-must
names(logfile)[12]<<-"workflow_must"
logfile[[10]]<<-3.100
names(logfile)[10]<<-"version"
save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(logfile$version<3.101){
cat("\n Updating to version 3.101 ...")
if(!file.exists(file.path(logfile$project_folder,"quantification"))){
dir.create(file.path(logfile$project_folder,"quantification"),recursive=TRUE) # subfolder
}
intstand<-read.table(file=file.path(logfile[[1]],"dataframes","IS.txt"),header=TRUE,sep="\t",colClasses = "character");
if(any(names(intstand)=="Lower intensity bound")){
intstand<-intstand[,!(names(intstand)=="Lower intensity bound")]
intstand<-intstand[,!(names(intstand)=="Lower intensity bound")]
}
if(any(names(intstand)=="Lower.intensity.bound")){
intstand<-intstand[,!(names(intstand)=="Lower.intensity.bound")]
intstand<-intstand[,!(names(intstand)=="Upper.intensity.bound")]
}
if(!any(names(intstand)=="Lower_intensity_bound")){
names_1<-names(intstand)
intstand<-cbind(intstand,rep(0,length(intstand[,1])),rep(Inf,length(intstand[,1])))
names(intstand)<-c(names_1,"Lower_intensity_bound","Upper_intensity_bound")
}
if(any(names(intstand)=="Use_for_screening")){
names(intstand)[names(intstand)=="Use_for_screening"]<-"use_for_screening"
}
if(any(names(intstand)=="Use_for_recalibration")){
names(intstand)[names(intstand)=="Use_for_recalibration"]<-"use_for_recalibration"
}
write.table(intstand,file=file.path(logfile[[1]],"dataframes","IS.txt"),row.names=FALSE,sep="\t",quote=FALSE)
rm(intstand)
if(!any(names(logfile$parameters)=="subtract_pos_bydate")){
logfile$parameters$subtract_pos_bydate<<-"FALSE";
logfile$parameters$subtract_pos_byfile<<-"FALSE";
logfile$parameters$subtract_neg_bydate<<-"FALSE";
logfile$parameters$subtract_neg_byfile<<-"FALSE";
save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(!any(names(logfile$parameters)=="blind_omit")){
logfile$parameters$blind_omit<<-"FALSE";
save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(!any(names(logfile$parameters)=="prof_select")){
logfile$parameters$prof_select<<-"FALSE";
save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(!any(names(logfile$parameters)=="trend_blind")){
logfile$parameters$trend_blind<<-"yes";
save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(!any(names(logfile$parameters)=="replicate_IS_dInt")){
logfile$parameters$replicate_IS_dInt<<-"5";
save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(!any(names(logfile$parameters)=="replicates_prof")){
logfile$parameters$replicates_prof<<-"yes";
save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(!any(names(logfile$parameters)=="screen_IS_maxonly")){
logfile$parameters$screen_target_maxonly<<-"FALSE";
logfile$parameters$screen_IS_maxonly<<-"FALSE";
save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
targets<-read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character");
if(!any(names(targets)=="warn_1")){
names_1<-names(targets)
targets<-cbind(targets,rep("FALSE",length(targets[,1])),rep("FALSE",length(targets[,1])))
names(targets)<-c(names_1,"warn_1","warn_2")
}
if(any(names(targets)=="Use_for_screening")){
names(targets)[names(targets)=="Use_for_screening"]<-"use_for_screening"
}
if(any(names(targets)=="Use_for_recalibration")){
names(targets)[names(targets)=="Use_for_recalibration"]<-"use_for_recalibration"
}
if(any(names(targets)=="intensity_warn_1")){
targets<-targets[,!(names(targets)=="intensity_warn_1")]
targets<-targets[,!(names(targets)=="intensity_warn_2")]
}
if(any(names(targets)=="RT.tolerance")){
names(targets)[names(targets)=="RT.tolerance"]<-"RT_tolerance"
}
if(any(names(targets)=="intercept")){
targets<-targets[,names(targets)!="intercept"]
}
if(any(names(targets)=="slope")){
targets<-targets[,names(targets)!="slope"]
}
write.table(targets,file=file.path(logfile[[1]],"dataframes","targets.txt"),row.names=FALSE,sep="\t",quote=FALSE)
rm(targets)
IDs<-list.files(file.path(logfile[[1]],"peaklist"))
if(length(IDs)>0){
for(i in 1:length(IDs)){
load(file=file.path(logfile[[1]],"peaklist",as.character(IDs[i])),envir=as.environment(".GlobalEnv"),verbose=FALSE);
if(any(colnames(peaklist)=="keep_2")){next}
keep_2<-rep(1,length(peaklist[,1]))
peaklist<-cbind(peaklist,keep_2)
colnames(peaklist)[16]<-"keep_2";
save(peaklist,file=file.path(logfile[[1]],"peaklist",as.character(IDs[i])))
rm(peaklist)
}
}
if(!any(names(logfile)=="Positive_subtraction_files")){
logfile[[13]]<<-"FALSE"
names(logfile)[13]<<-"Positive_subtraction_files"
}
if(!any(names(logfile)=="Negative_subtraction_files")){
logfile[[14]]<<-"FALSE"
names(logfile)[14]<<-"Negative_subtraction_files"
}
logfile[[10]]<<-3.101
names(logfile)[10]<<-"version"
save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(logfile$version<3.102){
cat("\n Updating to version 3.102 ...")
intstand<-read.table(file=file.path(logfile[[1]],"dataframes","IS.txt"),header=TRUE,sep="\t",colClasses = "character");
if(!any(names(intstand)=="Quant_adduct")){
names_1<-names(intstand)
quant_add<-rep("",length(intstand[,1]))
quant_add[intstand[,names(intstand)=="ion_mode"]=="positive"]<-"M+H"
quant_add[intstand[,names(intstand)=="ion_mode"]=="negative"]<-"M-H"
intstand<-cbind(intstand,quant_add,rep("1",length(intstand[,1])))
names(intstand)<-c(names_1,"Quant_adduct","Quant_peak")
}
write.table(intstand,file=file.path(logfile[[1]],"dataframes","IS.txt"),row.names=FALSE,sep="\t",quote=FALSE)
rm(intstand)
targets<-read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character");
if(!any(names(targets)=="Quant_adduct")){
names_1<-names(targets)
quant_add<-rep("",length(targets[,1]))
quant_add[targets[,names(targets)=="ion_mode"]=="positive"]<-"M+H"
quant_add[targets[,names(targets)=="ion_mode"]=="negative"]<-"M-H"
targets<-cbind(targets,quant_add,rep("1",length(targets[,1])))
names(targets)<-c(names_1,"Quant_adduct","Quant_peak")
}
write.table(targets,file=file.path(logfile[[1]],"dataframes","targets.txt"),row.names=FALSE,sep="\t",quote=FALSE)
rm(targets)
if(!any(names(logfile$parameters)=="recal_maxdmz")){
logfile$parameters$recal_maxdmz<<-"30";
save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(!any(logfile$summary[,1]=="calibration")){
logfile$summary[,1]<<-as.character(logfile$summary[,1])
logfile$summary[,2]<<-as.character(logfile$summary[,2])
logfile$summary[20,1]<<-"calibration"
logfile$summary[20,2]<<-"FALSE"
save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(!any(logfile$summary[,1]=="isotopologues")){
logfile$summary[,1]<<-as.character(logfile$summary[,1])
logfile$summary[,2]<<-as.character(logfile$summary[,2])
logfile$summary[21,1]<<-"isotopologues"
logfile$summary[21,2]<<-"FALSE"
logfile$summary[22,1]<<-"adducts"
logfile$summary[22,2]<<-"FALSE"
logfile$summary[23,1]<<-"homologues"
logfile$summary[23,2]<<-"FALSE"
save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(!any(logfile$summary[,1]=="recovery",na.rm=TRUE)){
logfile$summary[,1]<<-as.character(logfile$summary[,1])
logfile$summary[,2]<<-as.character(logfile$summary[,2])
logfile$summary[24,1]<<-"recovery"
logfile$summary[24,2]<<-"FALSE"
first_ext<-TRUE;
save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(!any(names(logfile$workflow)=="calibration")){
logfile$workflow[19]<<-"yes"; 	names(logfile$workflow)[19]<<-"calibration"
}
if(!any(names(logfile$workflow)=="isotopologues")){
logfile$workflow[20]<<-"yes"; 	names(logfile$workflow)[20]<<-"isotopologues"
}
if(!any(names(logfile$workflow)=="adducts")){
logfile$workflow[21]<<-"yes"; 	names(logfile$workflow)[21]<<-"adducts"
}
if(!any(names(logfile$workflow)=="homologues")){
logfile$workflow[22]<<-"yes"; 	names(logfile$workflow)[22]<<-"homologues"
}
if(!any(names(logfile$workflow)=="recovery")){
logfile$workflow[23]<<-"yes"; 	names(logfile$workflow)[23]<<-"recovery"
}
measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
if( !any(names(measurements)=="date_end") ){
measurements<-cbind(
measurements,
rep("2019-06-08",length(measurements[,"ID"])),
rep("12:00:00",length(measurements[,"ID"]))
)
}
if( any(names(measurements)=="feat.") ){
names(measurements)[names(measurements)=="feat."]<-"profiled";
measurements[,names(measurements)=="profiled"]<-"TRUE";
}
names(measurements)<-c("ID","Name","Type","Mode","Place","Date","Time","include","copied","picked",
"checked","recal","align","norm","profiled","LOD","IS_screen","tar_screen","tag1","tag2","tag3","date_end","time_end")
write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
rm(measurements)
logfile[[10]]<<-3.102
names(logfile)[10]<<-"version"
save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(logfile$version<3.103){
cat("\n Updating to version 3.103 ...")
old_Tasks_to_redo<-logfile$Tasks_to_redo
logfile[[2]]<<-rep(FALSE,21);
names(logfile[[2]])<<-c(
"peakpicking","qc","recal","norm","align","profiled","trendblind","pattern",
"replicates","IS_screen","target_screen","LOD","calibration","recovery","quantification","blind",
"IS_normaliz","subtr","isotopologues","adducts","homologues"
)
names(logfile)[2]<<-c("Tasks_to_redo");
for(i in 1:length(old_Tasks_to_redo)){
if(any(names(logfile$Tasks_to_redo)==names(old_Tasks_to_redo)[i])){
logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)==names(old_Tasks_to_redo)[i]]<<-old_Tasks_to_redo[i]
}
}
old_workflow<-logfile$workflow
logfile$workflow<<-0
names(logfile)[6]<<-c("workflow")
for(i in 1:length(names(logfile[[2]]))){
logfile$workflow[i]<<-"yes";
names(logfile$workflow)[i]<<-names(logfile[[2]])[i]
}
for(i in 1:length(old_workflow)){
if(any(names(logfile$workflow)==names(old_workflow)[i])){
logfile$workflow[names(logfile$workflow)==names(old_workflow)[i]]<<-old_workflow[i]
}
}
old_summary<-logfile$summary
tasks<-names(logfile[[2]])
doneit<-rep(FALSE,length(tasks))
summar<-data.frame(tasks,doneit,stringsAsFactors = FALSE)
names(summar)<-c("Tasks","Done?")
logfile[[3]]<<-summar
names(logfile)[3]<<-c("summary")
for(i in 1:length(old_summary[,1])){
if(any(logfile$summary[,1]==as.character(old_summary[i,1]))){
logfile$summary[logfile$summary[,1]==as.character(old_summary[i,1]),2]<<-as.character(old_summary[i,2])
}
}
if(!any(names(logfile$parameters)=="subtr_IS")){
logfile$parameters$subtr_IS<<-"yes";
logfile$parameters$subtr_target<<-"yes";
logfile$parameters$subtr_blind<<-"yes";
logfile$parameters$subtr_spiked<<-"yes";
save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
# target_subtr or target screening result tables/lists for quantification ######################
work_names<-names(logfile$Tasks_to_redo)[1:21]
depend<-matrix(ncol=length(work_names),nrow=length(work_names),0)
colnames(depend)<-work_names
rownames(depend)<-work_names					# peakpicking	qc	recal	norm	align	profiled	trendblind	pattern		replicates	IS_screen	target_screeen	LOD		calibration	recovery	quantification	blind IS_normaliz	subtr	isotopologues	adducts homologues
depend[,colnames(depend)=="peakpicking"]<-		c(0,			1,	1,		1,		1,		1,			1,			0,			1,			1,			1,				1,		1,			1,			1,				1,		1,			1,		0,				0,		0)
depend[,colnames(depend)=="qc"]<-				c(0,			0,	1,		1,		1,		1,			1,			0,			1,			1,			1,				1,		1,			0,			1,				1,		1,			1,		0,				0,		0)
depend[,colnames(depend)=="pattern"]<-			c(0,			0,	1,		0,		0,		0,			0,			0,			0,			1,			1,				0,		0,			0,			0,				0,		0,			0,		0,				0,		0)
depend[,colnames(depend)=="recal"]<-			c(0,			0,	0,		0,		0,		1,			1,			0,			1,			1,			1,				0,		0,			0,			0,				0,		1,			0,		0,				1,		1)
depend[,colnames(depend)=="align"]<-			c(0,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,		0,				0,		0)
depend[,colnames(depend)=="norm"]<-				c(0,			0,	0,		0,		0,		1,			1,			0,			1,			1,			1,				1,		0,			0,			0,				1,		1,			0,		0,				0,		0)
depend[,colnames(depend)=="blind"]<-			c(0,			0,	0,		0,		0,		1,			1,			0,			0,			1,			1,				0,		0,			0,			0,				0,		1,			1,		1,				1,		1)
depend[,colnames(depend)=="replicates"]<-		c(0,			0,	0,		0,		0,		1,			1,			0,			0,			1,			1,				1,		0,			0,			0,				0,		1,			1,		1,				1,		1)
depend[,colnames(depend)=="profiled"]<-			c(0,			0,	0,		0,		0,		0,			1,			0,			0,			1,			1,				0,		0,			0,			0,				0,		1,			0,		0,				0,		0)
depend[,colnames(depend)=="IS_screen"]<-		c(0,			0,	0,		0,		0,		2,			0,			0,			0,			0,			0,				0,		0,			0,			1,				0,		1,			1,		0,				0,		0)
depend[,colnames(depend)=="target_screen"]<-	c(0,			0,	0,		0,		0,		2,			0,			0,			0,			0,			0,				0,		1,			0,			1,				0,		1,			1,		0,				0,		0)
depend[,colnames(depend)=="IS_normaliz"]<-		c(0,			0,	0,		0,		0,		0,			1,			0,			0,			0,			0,				0,		1,			0,			0,				0,		0,			0,		0,				0,		0)
depend[,colnames(depend)=="trendblind"]<-		c(0,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,		0,				0,		0)
depend[,colnames(depend)=="LOD"]<-				c(0,			0,	0,		0,		0,		0,			0,			0,			0,			1,			1,				0,		1,			0,			0,				0,		1,			1,		0,				0,		0)
depend[,colnames(depend)=="calibration"]<-		c(0,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			1,			1,				0,		0,			0,		0,				0,		0)
depend[,colnames(depend)=="recovery"]<-			c(0,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,		0,				0,		0)
depend[,colnames(depend)=="quantification"]<-	c(0,			0,	0,		0,		0,		0,			0,			0,			0,			0,			2,				0,		0,			1,			0,				0,		0,			0,		0,				0,		0)
depend[,colnames(depend)=="subtr"]<-			c(0,			0,	0,		0,		0,		2,			1,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,		0,				0,		0)
depend[,colnames(depend)=="isotopologues"]<-	c(0,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,		0,				0,		0)
depend[,colnames(depend)=="adducts"]<-			c(0,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,		0,				0,		0)
depend[,colnames(depend)=="homologues"]<-		c(0,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,		0,				0,		0)
logfile[[11]]<<-depend
names(logfile)[11]<<-"workflow_depend"
must<-matrix(ncol=length(work_names),nrow=length(work_names),0)
colnames(must)<-work_names
rownames(must)<-work_names					# peakpicking	qc	recal	norm	align	profiled	trendblind	pattern		replicates	IS_screen	target_screeen	LOD		calibration	recovery	quantification	blind IS_normaliz	subtr	isotopologues	adducts homologues
must[,colnames(must)=="peakpicking"]<-		c(0,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,		0,				0,		0)
must[,colnames(must)=="qc"]<-				c(1,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,		0,				0,		0)
must[,colnames(must)=="pattern"]<-			c(0,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,		0,				0,		0)
must[,colnames(must)=="recal"]<-			c(1,			0,	0,		0,		0,		0,			0,			1,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,		0,				0,		0)
must[,colnames(must)=="align"]<-			c(1,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,		0,				0,		0)
must[,colnames(must)=="norm"]<-				c(1,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,		0,				0,		0)
must[,colnames(must)=="blind"]<-			c(1,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,		0,				0,		0)
must[,colnames(must)=="replicates"]<-		c(1,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,		0,				0,		0)
must[,colnames(must)=="profiled"]<-			c(1,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,		0,				0,		0)
must[,colnames(must)=="IS_screen"]<-		c(1,			0,	0,		0,		0,		0,			0,			1,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,		0,				0,		0)
must[,colnames(must)=="target_screen"]<-	c(1,			0,	0,		0,		0,		0,			0,			1,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,		0,				0,		0)
must[,colnames(must)=="IS_normaliz"]<-		c(1,			0,	0,		0,		0,		0,			0,			1,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,		0,				0,		0)
must[,colnames(must)=="trendblind"]<-		c(1,			0,	0,		0,		0,		1,			0,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,		0,				0,		0)
must[,colnames(must)=="LOD"]<-				c(1,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,		0,				0,		0)
must[,colnames(must)=="calibration"]<-		c(1,			0,	0,		0,		0,		0,			0,			1,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,		0,				0,		0)
must[,colnames(must)=="recovery"]<-			c(1,			0,	0,		0,		0,		0,			0,			1,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,		0,				0,		0)
must[,colnames(must)=="quantification"]<-	c(1,			0,	0,		0,		0,		0,			0,			1,			0,			1,			1,				0,		0,			0,			0,				0,		0,			0,		0,				0,		0)
must[,colnames(must)=="subtr"]<-			c(1,			0,	0,		0,		0,		1,			0,			1,			0,			1,			1,				0,		0,			0,			0,				1,		0,			0,		0,				0,		0)
must[,colnames(must)=="isotopologues"]<-	c(1,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,		0,				0,		0)
must[,colnames(must)=="adducts"]<-			c(1,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,		0,				0,		0)
must[,colnames(must)=="homologues"]<-		c(1,			0,	0,		0,		0,		0,			0,			0,			0,			0,			0,				0,		0,			0,			0,				0,		0,			0,		0,				0,		0)
logfile[[12]]<<-must
names(logfile)[12]<<-"workflow_must"
schedule<-enviMass::PFMSW1067M(logfile$workflow_depend,logfile$workflow_must)
set_order<-match(schedule[,1],logfile$summary[,1])
logfile$summary<<-logfile$summary[set_order,]
logfile[[10]]<<-3.103
names(logfile)[10]<<-"version"
save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(logfile$version<3.104){
cat("\n Updating to version 3.104 ...")
measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
names(measurements)[names(measurements)=="picked"]<-"peakpicking"
write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
rm(measurements);
logfile[[10]]<<-3.104
names(logfile)[10]<<-"version"
save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(logfile$version<3.105){
cat("\n Updating to version 3.105 ...")
if(!file.exists(file.path(logfile$project_folder,"results","componentization"))){
dir.create(file.path(logfile$project_folder,"results","componentization"),recursive=TRUE)
}
if(!file.exists(file.path(logfile$project_folder,"results","componentization","adducts"))){
dir.create(file.path(logfile$project_folder,"results","componentization","adducts"),recursive=TRUE)
}
if(!file.exists(file.path(logfile$project_folder,"results","componentization","isotopologues"))){
dir.create(file.path(logfile$project_folder,"results","componentization","isotopologues"),recursive=TRUE)
}
if(!file.exists(file.path(logfile$project_folder,"results","componentization","EIC_corr"))){
dir.create(file.path(logfile$project_folder,"results","componentization","EIC_corr"),recursive=TRUE)
}
measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
if( !any(names(measurements)=="isotopologues") ){
measurements<-cbind(
measurements,
rep("FALSE",length(measurements[,"ID"])),
rep("FALSE",length(measurements[,"ID"])),
rep("FALSE",length(measurements[,"ID"])),
rep("FALSE",length(measurements[,"ID"]))
)
}
names(measurements)<-c("ID","Name","Type","Mode","Place","Date","Time","include","copied","peakpicking",
"checked","recal","align","norm","profiled","LOD","IS_screen","tar_screen","tag1","tag2","tag3","date_end","time_end",
"isotopologues","adducts","homologues","EIC_correlation");
write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
rm(measurements);
workflow_depend<-read.table(
file="workflow_depend"
)
workflow_depend<-as.matrix(workflow_depend)
workflow_must<-read.table(
file="workflow_must"
)
workflow_must<-as.matrix(workflow_must)
logfile[[11]]<<-workflow_depend
names(logfile)[11]<<-"workflow_depend"
logfile[[12]]<<-workflow_must
names(logfile)[12]<<-"workflow_must"
old_Tasks_to_redo<-logfile$Tasks_to_redo
logfile[[2]]<<-rep(FALSE,length(colnames(workflow_must)));
names(logfile[[2]])<<-colnames(workflow_must)
names(logfile)[2]<<-c("Tasks_to_redo");
for(i in 1:length(old_Tasks_to_redo)){
if(any(names(logfile$Tasks_to_redo)==names(old_Tasks_to_redo)[i])){
logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)==names(old_Tasks_to_redo)[i]]<<-old_Tasks_to_redo[i]
}
}
old_workflow<-logfile$workflow
logfile$workflow<<-0
names(logfile)[6]<<-c("workflow")
for(i in 1:length(names(logfile[[2]]))){
logfile$workflow[i]<<-"yes";
names(logfile$workflow)[i]<<-names(logfile[[2]])[i]
}
for(i in 1:length(old_workflow)){
if(any(names(logfile$workflow)==names(old_workflow)[i])){
logfile$workflow[names(logfile$workflow)==names(old_workflow)[i]]<<-old_workflow[i]
}
}
old_summary<-logfile$summary
tasks<-names(logfile[[2]])
doneit<-rep(FALSE,length(tasks))
summar<-data.frame(tasks,doneit,stringsAsFactors = FALSE)
names(summar)<-c("Tasks","Done?")
logfile[[3]]<<-summar
names(logfile)[3]<<-c("summary")
for(i in 1:length(old_summary[,1])){
if(any(logfile$summary[,1]==as.character(old_summary[i,1]))){
logfile$summary[logfile$summary[,1]==as.character(old_summary[i,1]),2]<<-as.character(old_summary[i,2])
}
}
schedule<-enviMass::PFMSW1067M(logfile$workflow_depend,logfile$workflow_must)
set_order<-match(schedule[,1],logfile$summary[,1])
logfile$summary<<-logfile$summary[set_order,]
if(!any(names(logfile$parameters)=="external")){
logfile$parameters$external<<-list()
}
logfile[[10]]<<-3.105
names(logfile)[10]<<-"version"
save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(logfile$version<3.106){
cat("\n Updating to version 3.106 ...")
workflow_depend<-read.table(
file="workflow_depend"
)
workflow_depend<-as.matrix(workflow_depend)
workflow_must<-read.table(
file="workflow_must"
)
workflow_must<-as.matrix(workflow_must)
logfile[[11]]<<-workflow_depend
names(logfile)[11]<<-"workflow_depend"
logfile[[12]]<<-workflow_must
names(logfile)[12]<<-"workflow_must"
old_Tasks_to_redo<-logfile$Tasks_to_redo
logfile[[2]]<<-rep(FALSE,length(colnames(workflow_must)));
names(logfile[[2]])<<-colnames(workflow_must)
names(logfile)[2]<<-c("Tasks_to_redo");
for(i in 1:length(old_Tasks_to_redo)){
if(any(names(logfile$Tasks_to_redo)==names(old_Tasks_to_redo)[i])){
logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)==names(old_Tasks_to_redo)[i]]<<-old_Tasks_to_redo[i]
}
}
old_workflow<-logfile$workflow
logfile$workflow<<-0
names(logfile)[6]<<-c("workflow")
for(i in 1:length(names(logfile[[2]]))){
logfile$workflow[i]<<-"yes";
names(logfile$workflow)[i]<<-names(logfile[[2]])[i]
}
for(i in 1:length(old_workflow)){
if(any(names(logfile$workflow)==names(old_workflow)[i])){
logfile$workflow[names(logfile$workflow)==names(old_workflow)[i]]<<-old_workflow[i]
}
}
old_summary<-logfile$summary
tasks<-names(logfile[[2]])
doneit<-rep(FALSE,length(tasks))
summar<-data.frame(tasks,doneit,stringsAsFactors = FALSE)
names(summar)<-c("Tasks","Done?")
logfile[[3]]<<-summar
names(logfile)[3]<<-c("summary")
for(i in 1:length(old_summary[,1])){
if(any(logfile$summary[,1]==as.character(old_summary[i,1]))){
logfile$summary[logfile$summary[,1]==as.character(old_summary[i,1]),2]<<-as.character(old_summary[i,2])
}
}
schedule<-enviMass::PFMSW1067M(logfile$workflow_depend,logfile$workflow_must)
set_order<-match(schedule[,1],logfile$summary[,1])
logfile$summary<<-logfile$summary[set_order,]
logfile[[10]]<<-3.106
names(logfile)[10]<<-"version"
save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(logfile$version<3.107){
cat("\n Updating to version 3.107 ...")
logfile$parameters$external$EICor_delRT<<-8
logfile$parameters$external$EICor_minpeaks<<-7
logfile$parameters$external$EICor_mincor<<-.9
logfile$version<<-3.107
save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(logfile$version<3.108){
cat("\n Updating to version 3.108 ...")
measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
if( !any(names(measurements)=="blind") ){
measurements<-cbind(
measurements,
rep("FALSE",length(measurements[,"ID"]))
)
names(measurements)<-c("ID","Name","Type","Mode","Place","Date","Time","include","copied","peakpicking",
"checked","recal","align","norm","profiled","LOD","IS_screen","tar_screen","tag1","tag2","tag3","date_end","time_end",
"isotopologues","adducts","homologues","EIC_correlation","blind")
}
write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
logfile$version<<-3.108
save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(logfile$version<3.109){
cat("\n Updating to version 3.109 ...")
if(!any(names(logfile$parameters$external)=="isotop_rttol")){
logfile$parameters$external$isotop_mztol<<-3
logfile$parameters$external$isotop_ppm<<-TRUE
logfile$parameters$external$isotop_inttol<<-0.5
logfile$parameters$external$isotop_rttol<<-15
logfile$parameters$external$isotop_use_charges<<-c(1,2)
}
logfile$version<<-3.109
save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(logfile$version<3.111){
cat("\n Updating to version 3.111 ...")
workflow_depend<-read.table(
file="workflow_depend"
)
workflow_depend<-as.matrix(workflow_depend)
workflow_must<-read.table(
file="workflow_must"
)
workflow_must<-as.matrix(workflow_must)
logfile[[11]]<<-workflow_depend
names(logfile)[11]<<-"workflow_depend"
logfile[[12]]<<-workflow_must
names(logfile)[12]<<-"workflow_must"
old_Tasks_to_redo<-logfile$Tasks_to_redo
logfile[[2]]<<-rep(FALSE,length(colnames(workflow_must)));
names(logfile[[2]])<<-colnames(workflow_must)
names(logfile)[2]<<-c("Tasks_to_redo");
for(i in 1:length(old_Tasks_to_redo)){
if(any(names(logfile$Tasks_to_redo)==names(old_Tasks_to_redo)[i])){
logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)==names(old_Tasks_to_redo)[i]]<<-old_Tasks_to_redo[i]
}
}
old_workflow<-logfile$workflow
logfile$workflow<<-0
names(logfile)[6]<<-c("workflow")
for(i in 1:length(names(logfile[[2]]))){
logfile$workflow[i]<<-"yes";
names(logfile$workflow)[i]<<-names(logfile[[2]])[i]
}
for(i in 1:length(old_workflow)){
if(any(names(logfile$workflow)==names(old_workflow)[i])){
logfile$workflow[names(logfile$workflow)==names(old_workflow)[i]]<<-old_workflow[i]
}
}
old_summary<-logfile$summary
tasks<-names(logfile[[2]])
doneit<-rep(FALSE,length(tasks))
summar<-data.frame(tasks,doneit,stringsAsFactors = FALSE)
names(summar)<-c("Tasks","Done?")
logfile[[3]]<<-summar
names(logfile)[3]<<-c("summary")
for(i in 1:length(old_summary[,1])){
if(any(logfile$summary[,1]==as.character(old_summary[i,1]))){
logfile$summary[logfile$summary[,1]==as.character(old_summary[i,1]),2]<<-as.character(old_summary[i,2])
}
}
schedule<-enviMass::PFMSW1067M(logfile$workflow_depend,logfile$workflow_must)
set_order<-match(schedule[,1],logfile$summary[,1])
logfile$summary<<-logfile$summary[set_order,]
logfile$parameters$external$adducts_rttol<<-10
logfile$parameters$external$adducts_mztol<<-3
logfile$parameters$external$adducts_ppm<<-TRUE
logfile$parameters$external$adducts_pos<<-c("M+H","M+Na","M+K","M+NH4")
logfile$parameters$external$adducts_neg<<-c("M-H","M-","2M-H")
logfile$version<<-3.111
save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(logfile$version<3.112){
cat("\n Updating to version 3.112 ...")
if(!file.exists(file.path(logfile$project_folder,"results","componentization","homologues"))){
dir.create(file.path(logfile$project_folder,"results","componentization","homologues"),recursive=TRUE)
}
logfile$version<<-3.112
save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(logfile$version<3.114){
cat("\n Updating to version 3.114 ...")
workflow_depend<-read.table(
file="workflow_depend"
)
workflow_depend<-as.matrix(workflow_depend)
workflow_must<-read.table(
file="workflow_must"
)
workflow_must<-as.matrix(workflow_must)
logfile[[11]]<<-workflow_depend
names(logfile)[11]<<-"workflow_depend"
logfile[[12]]<<-workflow_must
names(logfile)[12]<<-"workflow_must"
old_Tasks_to_redo<-logfile$Tasks_to_redo
logfile[[2]]<<-rep(FALSE,length(colnames(workflow_must)));
names(logfile[[2]])<<-colnames(workflow_must)
names(logfile)[2]<<-c("Tasks_to_redo");
for(i in 1:length(old_Tasks_to_redo)){
if(any(names(logfile$Tasks_to_redo)==names(old_Tasks_to_redo)[i])){
logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)==names(old_Tasks_to_redo)[i]]<<-old_Tasks_to_redo[i]
}
}
old_workflow<-logfile$workflow
logfile$workflow<<-0
names(logfile)[6]<<-c("workflow")
for(i in 1:length(names(logfile[[2]]))){
logfile$workflow[i]<<-"yes";
names(logfile$workflow)[i]<<-names(logfile[[2]])[i]
}
for(i in 1:length(old_workflow)){
if(any(names(logfile$workflow)==names(old_workflow)[i])){
logfile$workflow[names(logfile$workflow)==names(old_workflow)[i]]<<-old_workflow[i]
}
}
old_summary<-logfile$summary
tasks<-names(logfile[[2]])
doneit<-rep(FALSE,length(tasks))
summar<-data.frame(tasks,doneit,stringsAsFactors = FALSE)
names(summar)<-c("Tasks","Done?")
logfile[[3]]<<-summar
names(logfile)[3]<<-c("summary")
for(i in 1:length(old_summary[,1])){
if(any(logfile$summary[,1]==as.character(old_summary[i,1]))){
logfile$summary[logfile$summary[,1]==as.character(old_summary[i,1]),2]<<-as.character(old_summary[i,2])
}
}
schedule<-enviMass::PFMSW1067M(logfile$workflow_depend,logfile$workflow_must)
set_order<-match(schedule[,1],logfile$summary[,1])
logfile$summary<<-logfile$summary[set_order,]
logfile$parameters$external$homol_units<<-c("CH2","CH2O","CF2")
logfile$parameters$external$homol_charges<<-c(1,2,3)
logfile$parameters$external$homol_minmz<<-10
logfile$parameters$external$homol_maxmz<<-120
logfile$parameters$external$homol_minrt<<-60
logfile$parameters$external$homol_maxrt<<-60
logfile$parameters$external$homol_ppm<<-TRUE
logfile$parameters$external$homol_mztol<<-3
logfile$parameters$external$homol_rttol<<-15
logfile$parameters$external$homol_minlength<<-5
logfile$parameters$external$homol_vec_size<<-1E7
logfile$version<<-3.114
save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(logfile$version<3.115){
cat("\n Updating to version 3.115 ...")
measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
if( !any(names(measurements)=="ID_2") ){
measurements<-cbind(
measurements,
rep("FALSE",length(measurements[,"ID"]))
)
names(measurements)<-c("ID","Name","Type","Mode","Place","Date","Time","include","copied","peakpicking",
"checked","recal","align","norm","profiled","LOD","IS_screen","tar_screen","tag1","tag2","tag3","date_end","time_end",
"isotopologues","adducts","homologues","EIC_correlation","blind","ID_2")
}
write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
logfile$version<<-3.115
save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(logfile$version<3.116){
cat("\n Updating to version 3.116 ...")
if(!any(names(logfile$parameters)=="quant_files_included")){
logfile$parameters$quant_files_included<<-20
}
if(!any(names(logfile$parameters)=="recov_files_included")){
logfile$parameters$recov_files_included<<-20
}
targets<-read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character");
if(!any(names(targets)=="Quant_rule")){
names_1<-names(targets)
targets<-cbind(targets,rep("most intense peak",length(targets[,1])))
names(targets)<-c(names_1,"Quant_rule")
}
write.table(targets,file=file.path(logfile[[1]],"dataframes","targets.txt"),row.names=FALSE,sep="\t",quote=FALSE)
rm(targets)
intstand<-read.table(file=file.path(logfile[[1]],"dataframes","IS.txt"),header=TRUE,sep="\t",colClasses = "character");
if(!any(names(intstand)=="Quant_rule")){
names_1<-names(intstand)
intstand<-cbind(intstand,rep("most intense peak",length(intstand[,1])))
names(intstand)<-c(names_1,"Quant_rule")
}
write.table(intstand,file=file.path(logfile[[1]],"dataframes","IS.txt"),row.names=FALSE,sep="\t",quote=FALSE)
rm(intstand)
logfile$version<<-3.116
save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(logfile$version<3.117){
cat("\n Updating to version 3.117 ...")
if(!any(names(logfile$parameters)=="peak_which_intensity")){
logfile$parameters$peak_which_intensity<<-"maximum"
these_files<-list.files(file.path(logfile[[1]],"peaklist"))
for(i in 1:length(these_files)){
load(file=file.path(logfile[[1]],"peaklist",these_files[i]),envir=as.environment(".GlobalEnv"),verbose=FALSE);
colnames(peaklist)[13]<<-"int_corr";
save(peaklist,file=file.path(logfile[[1]],"peaklist",these_files[i]))
}
}
logfile$version<<-3.117
save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(logfile$version<3.118){
cat("\n Updating to version 3.118 ...")
if(!any(names(logfile$parameters)=="screen_IS_restrict")){
logfile$parameters$screen_IS_restrict<<-"FALSE";
logfile$parameters$screen_IS_restrict_many<<-"10";
logfile$parameters$screen_target_restrict<<-"FALSE";
logfile$parameters$screen_target_restrict_many<<-"10";
}
logfile$version<<-3.118
save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(logfile$version<3.119){
cat("\n Updating to version 3.119 ...")
if(!any(names(logfile$parameters)=="isotop_mztol")){
logfile$parameters$isotop_mztol<<-"2.5"
logfile$parameters$isotop_ppm<<-"TRUE"
logfile$parameters$isotop_inttol<<-"0.5"
logfile$parameters$isotop_rttol<<-"5"
logfile$parameters$isotop_use_charges<<-"FALSE"
}
logfile$version<<-3.119
save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(logfile$version<3.120){
cat("\n Updating to version 3.120 ...")
if(!any(names(logfile)=="adducts_pos_group")){
logfile[[15]]<<-c("M+H","M+Na","M+NH4","M+K");
names(logfile)[15]<<-c("adducts_pos_group")
logfile[[16]]<<-c("M-H","M-");
names(logfile)[16]<<-c("adducts_neg_group")
}
if(!any(names(logfile)=="isotopes")){
logfile[[9]]<<-"";
names(logfile)[9]<<-c("isotopes")
}
if(!any(names(logfile$parameters)=="adducts_rttol")){
logfile$parameters$adducts_rttol<<-"5"
logfile$parameters$adducts_mztol<<-"2.5"
logfile$parameters$adducts_ppm<<-TRUE
}
logfile$version<<-3.120
save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(logfile$version<3.121){
cat("\n Updating to version 3.121 ...")
if(!any(names(logfile$parameters)=="homol_units")){
logfile$parameters$homol_units<<-c("CH2,CH4O")
logfile$parameters$homol_charges<<-c("1,2")
logfile$parameters$homol_minmz<<-"10"
logfile$parameters$homol_maxmz<<-"120"
logfile$parameters$homol_minrt<<-"10"
logfile$parameters$homol_maxrt<<-"60"
logfile$parameters$homol_ppm<<-"TRUE"
logfile$parameters$homol_mztol<<-"2.5"
logfile$parameters$homol_rttol<<-"20"
logfile$parameters$homol_minlength<<-"6"
logfile$parameters$homol_vec_size<<-"1E8"
}
if(!any(names(logfile$parameters)=="EICor_delRT")){
logfile$parameters$EICor_delRT<<-"5"
logfile$parameters$EICor_minpeaks<<-"10"
logfile$parameters$EICor_mincor<<-".95"
}
if(!any(names(logfile$parameters)=="prof_comp_maxfiles")){
logfile$parameters$prof_comp_maxfiles<<-"15"
}
if(!file.exists(file.path(as.character(logfile$project_folder),"results","componentization","components"))){
dir.create(file.path(as.character(logfile$project_folder),"results","componentization","components"),recursive=TRUE)
}
logfile$version<<-3.121
save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(logfile$version<3.122){
cat("\n Updating to version 3.122 ...")
if(any(names(logfile$workflow)=="components")){
names(logfile$workflow)[names(logfile$workflow)=="components"]<<-"components_files"
logfile$workflow[names(logfile$workflow)=="components_files"]<<-"no"
logfile$workflow[24]<<-"no"
names(logfile$workflow)[24]<<-"components_profiles"
names(logfile$Tasks_to_redo)[names(logfile$Tasks_to_redo)=="components"]<<-"components_files"
logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)=="components_files"]<<-"FALSE"
logfile$Tasks_to_redo[24]<<-"FALSE"
names(logfile$Tasks_to_redo)[24]<<-"components_profiles"
logfile$summary[
logfile$summary[,1]=="components"
,1]<<-"components_files"
logfile$summary[
logfile$summary[,1]=="components"
,2]<<-"FALSE"
logfile$summary<<-rbind(
logfile$summary,c("components_profiles","FALSE")
)
}
logfile$version<<-3.122
save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(logfile$version<3.123){
cat("\n Updating to version 3.123 ...")
workflow_depend<-read.table(file="workflow_depend")
workflow_depend<-as.matrix(workflow_depend)
workflow_must<-read.table(file="workflow_must")
workflow_must<-as.matrix(workflow_must)
logfile[["workflow_depend"]]<<-workflow_depend
logfile[["workflow_must"]]<<-workflow_must
schedule<-enviMass::PFMSW1067M(logfile$workflow_depend,logfile$workflow_must)
if(!is.data.frame(schedule)){stop("\nschedule not a data frame")}
set_order<-match(schedule[,1],logfile$summary[,1])
logfile$summary<<-logfile$summary[set_order,]
logfile$version<<-3.123
save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(logfile$version<3.124){
cat("\n Updating to version 3.124 ...")
measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
if( !any(names(measurements)=="components_files") ){
measurements<-cbind(
measurements,
rep("FALSE",length(measurements[,"ID"]))
)
names(measurements)<-c(names(measurements)[-30],"components_files")
}
write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
rm(measurements)
logfile$version<<-3.124
save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(logfile$version<3.125){
cat("\n Updating to version 3.125 ...")
workflow_depend<-read.table(file="workflow_depend")
workflow_depend<-as.matrix(workflow_depend)
workflow_must<-read.table(file="workflow_must")
workflow_must<-as.matrix(workflow_must)
logfile[["workflow_depend"]]<<-workflow_depend
logfile[["workflow_must"]]<<-workflow_must
schedule<-enviMass::PFMSW1067M(logfile$workflow_depend,logfile$workflow_must)
if(!is.data.frame(schedule)){stop("\nschedule not a data frame")}
set_order<-match(schedule[,1],logfile$summary[,1])
logfile$summary<<-logfile$summary[set_order,]
logfile$version<<-3.125
save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(logfile$version<3.126){
cat("\n Updating to version 3.126 ...")
if(!any(names(logfile$parameters)=="cut_RT")){
logfile$parameters$cut_RT<<-"FALSE"
logfile$parameters$cut_RT_min<<-"0"
logfile$parameters$cut_RT_max<<-"25"
logfile$parameters$cut_mass<<-"FALSE"
logfile$parameters$cut_mass_min<<-"0"
logfile$parameters$cut_mass_max<<-"2000"
}
logfile$version<<-3.126
save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(logfile$version<3.127){
cat("\n Updating to version 3.127 ...")
workflow_depend<-read.table(
file="workflow_depend"
)
workflow_depend<-as.matrix(workflow_depend)
workflow_must<-read.table(
file="workflow_must"
)
workflow_must<-as.matrix(workflow_must)
logfile[[11]]<<-workflow_depend
names(logfile)[11]<<-"workflow_depend"
logfile[[12]]<<-workflow_must
names(logfile)[12]<<-"workflow_must"
old_Tasks_to_redo<-logfile$Tasks_to_redo
logfile[[2]]<<-rep(FALSE,length(colnames(workflow_must)));
names(logfile[[2]])<<-colnames(workflow_must)
names(logfile)[2]<<-c("Tasks_to_redo");
for(i in 1:length(old_Tasks_to_redo)){
if(any(names(logfile$Tasks_to_redo)==names(old_Tasks_to_redo)[i])){
logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)==names(old_Tasks_to_redo)[i]]<<-old_Tasks_to_redo[i]
}
}
old_workflow<-logfile$workflow
logfile$workflow<<-0
names(logfile)[6]<<-c("workflow")
for(i in 1:length(names(logfile[[2]]))){
logfile$workflow[i]<<-"yes";
names(logfile$workflow)[i]<<-names(logfile[[2]])[i]
}
for(i in 1:length(old_workflow)){
if(any(names(logfile$workflow)==names(old_workflow)[i])){
logfile$workflow[names(logfile$workflow)==names(old_workflow)[i]]<<-old_workflow[i]
}
}
old_summary<-logfile$summary
tasks<-names(logfile[[2]])
doneit<-rep(FALSE,length(tasks))
summar<-data.frame(tasks,doneit,stringsAsFactors = FALSE)
names(summar)<-c("Tasks","Done?")
logfile[[3]]<<-summar
names(logfile)[3]<<-c("summary")
for(i in 1:length(old_summary[,1])){
if(any(logfile$summary[,1]==as.character(old_summary[i,1]))){
logfile$summary[logfile$summary[,1]==as.character(old_summary[i,1]),2]<<-as.character(old_summary[i,2])
}
}
schedule<-enviMass::PFMSW1067M(logfile$workflow_depend,logfile$workflow_must)
set_order<-match(schedule[,1],logfile$summary[,1])
logfile$summary<<-logfile$summary[set_order,]
if(!file.exists(file.path(logfile$project_folder,"results","recalibration"))){
dir.create(file.path(logfile$project_folder,"results","recalibration"),recursive=TRUE)
}
if(!file.exists(file.path(logfile$project_folder,"quantification"))){
dir.create(file.path(logfile$project_folder,"quantification"),recursive=TRUE)   # subfolder
}
if(!file.exists(file.path(logfile$project_folder,"results","LOD"))){
dir.create(file.path(logfile$project_folder,"results","LOD"),recursive=TRUE)
}
if(!file.exists(file.path(logfile$project_folder,"results","recalibration"))){
dir.create(file.path(logfile$project_folder,"results","recalibration"),recursive=TRUE)
}
if(!file.exists(file.path(logfile$project_folder,"results","screening"))){
dir.create(file.path(logfile$project_folder,"results","screening"),recursive=TRUE)
}
if(!file.exists(file.path(logfile$project_folder,"results","componentization"))){
dir.create(file.path(logfile$project_folder,"results","componentization"),recursive=TRUE)
}
if(!file.exists(file.path(logfile$project_folder,"results","componentization","adducts"))){
dir.create(file.path(logfile$project_folder,"results","componentization","adducts"),recursive=TRUE)
}
if(!file.exists(file.path(logfile$project_folder,"results","componentization","isotopologues"))){
dir.create(file.path(logfile$project_folder,"results","componentization","isotopologues"),recursive=TRUE)
}
if(!file.exists(file.path(logfile$project_folder,"results","componentization","EIC_corr"))){
dir.create(file.path(logfile$project_folder,"results","componentization","EIC_corr"),recursive=TRUE)
}
if(!file.exists(file.path(logfile$project_folder,"results","componentization","homologues"))){
dir.create(file.path(logfile$project_folder,"results","componentization","homologues"),recursive=TRUE)
}
if(!file.exists(file.path(logfile$project_folder,"results","componentization","components"))){
dir.create(file.path(logfile$project_folder,"results","componentization","components"),recursive=TRUE)
}
logfile$version<<-3.127
save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(logfile$version<3.128){
cat("\n Updating to version 3.128 ...")
if(!any(names(logfile$parameters)=="peak_estimate")){
logfile$parameters$peak_estimate<<-"TRUE"
}
logfile$version<<-3.128
save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(logfile$version<3.2){
cat("\n Updating to version 3.2 ...")
if(!any(names(logfile$parameters)=="is_example")){
logfile$parameters$is_example<<-"FALSE"
}
logfile$version<<-3.2
save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(logfile$version<3.26){
cat("\n Updating to version 3.26 ...")
workflow_depend<-read.table(
file="workflow_depend"
)
workflow_depend<-as.matrix(workflow_depend)
workflow_must<-read.table(
file="workflow_must"
)
workflow_must<-as.matrix(workflow_must)
logfile[[11]]<<-workflow_depend
names(logfile)[11]<<-"workflow_depend"
logfile[[12]]<<-workflow_must
names(logfile)[12]<<-"workflow_must"
old_Tasks_to_redo<-logfile$Tasks_to_redo
logfile[[2]]<<-rep(FALSE,length(colnames(workflow_must)));
names(logfile[[2]])<<-colnames(workflow_must)
names(logfile)[2]<<-c("Tasks_to_redo");
for(i in 1:length(old_Tasks_to_redo)){
if(any(names(logfile$Tasks_to_redo)==names(old_Tasks_to_redo)[i])){
logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)==names(old_Tasks_to_redo)[i]]<<-old_Tasks_to_redo[i]
}
}
old_workflow<-logfile$workflow
logfile$workflow<<-0
names(logfile)[6]<<-c("workflow")
for(i in 1:length(names(logfile[[2]]))){
logfile$workflow[i]<<-"yes";
names(logfile$workflow)[i]<<-names(logfile[[2]])[i]
}
for(i in 1:length(old_workflow)){
if(any(names(logfile$workflow)==names(old_workflow)[i])){
logfile$workflow[names(logfile$workflow)==names(old_workflow)[i]]<<-old_workflow[i]
}
}
old_summary<-logfile$summary
tasks<-names(logfile[[2]])
doneit<-rep(FALSE,length(tasks))
summar<-data.frame(tasks,doneit,stringsAsFactors = FALSE)
names(summar)<-c("Tasks","Done?")
logfile[[3]]<<-summar
names(logfile)[3]<<-c("summary")
for(i in 1:length(old_summary[,1])){
if(any(logfile$summary[,1]==as.character(old_summary[i,1]))){
logfile$summary[logfile$summary[,1]==as.character(old_summary[i,1]),2]<<-as.character(old_summary[i,2])
}
}
schedule<-enviMass::PFMSW1067M(logfile$workflow_depend,logfile$workflow_must)
set_order<-match(schedule[,1],logfile$summary[,1])
logfile$summary<<-logfile$summary[set_order,]
if(!any(names(logfile$parameters)=="peak_which_intensity")){
logfile$parameters$peak_which_intensity<<-"maximum"
}
if(!any(names(logfile$parameters)=="quant_files_included")){
logfile$parameters$quant_files_included<<-"100"
}
if(!any(names(logfile$parameters)=="recov_files_included")){
logfile$parameters$recov_files_included<<-"20"
}
measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
if(any(names(measurements)=="checked")){
names(measurements)[names(measurements)=="checked"]<-"qc"
}
if(any(names(measurements)=="tar_screen")){
measurements<-measurements[,-(which(names(measurements)=="tar_screen"))]
}
if(any(names(measurements)=="IS_screen")){
measurements<-measurements[,-(which(names(measurements)=="IS_screen"))]
}
write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
rm(measurements)
logfile$version<<-3.26
save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(logfile$version<3.27){
cat("\n Updating to version 3.27 ...")
if(!any(names(logfile$parameters)=="recal_use_pos")){
logfile$parameters$recal_include_pos<<-"TRUE"
logfile$parameters$recal_use_pos<<-logfile$parameters$recal_use
logfile$parameters$recal_dmz_pos<<-logfile$parameters$recal_dmz
logfile$parameters$recal_ppm_pos<<-logfile$parameters$recal_ppm
logfile$parameters$recal_drt_pos<<-logfile$parameters$recal_drt
logfile$parameters$recal_maxdmz_pos<<-logfile$parameters$recal_maxdmz
logfile$parameters$recal_include_neg<<-"TRUE"
logfile$parameters$recal_use_neg<<-logfile$parameters$recal_use
logfile$parameters$recal_dmz_neg<<-logfile$parameters$recal_dmz
logfile$parameters$recal_ppm_neg<<-logfile$parameters$recal_ppm
logfile$parameters$recal_drt_neg<<-logfile$parameters$recal_drt
logfile$parameters$recal_maxdmz_neg<<-logfile$parameters$recal_maxdmz
logfile$parameters$recal_use<<-"NA"
logfile$parameters$recal_dmz<<-"NA"
logfile$parameters$recal_ppm<<-"NA"
logfile$parameters$recal_drt<<-"NA"
logfile$parameters$recal_maxdmz<<-"NA"
}
logfile$version<<-3.27
save(logfile,file=file.path(as.character(logfile[["project_folder"]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(logfile$version<3.28){
cat("\n Updating to version 3.28 ...")
IDs <- list.files(file.path(logfile[["project_folder"]],"peaklist"))
if(length(IDs) > 0){
for(i in 1:length(IDs)){
load(file=file.path(logfile[["project_folder"]],"peaklist",as.character(IDs[i])),envir=as.environment(".GlobalEnv"),verbose=FALSE);
peaklist[,"keep_2"] <- Inf;
save(peaklist,file=file.path(logfile[["project_folder"]],"peaklist",as.character(IDs[i])))
rm(peaklist)
}
}
enviMass::RNJJT9347D(
down="peakpicking",
except="peakpicking",
down_TF=c("TRUE","FALSE"),
check_node=TRUE,
single_file=FALSE
)
output$summa_html <- renderText(enviMass::UCBRR9756L(logfile$summary, logfile$Tasks_to_redo));
shinyjs::info(paste0("Sample peaklists now contain intensity ratios above blind - please press the Calculate button any time soon to make these changes permanent to your project results (entails a project recalculation except peakpicking)!"));
logfile$version<<-3.28
save(logfile,file=file.path(as.character(logfile[["project_folder"]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(logfile$version<3.29){
cat("\n Updating to version 3.29 ...")
if(
!any(names(logfile$Tasks_to_redo)=="profblind") ||
logfile$workflow_depend[rownames(logfile$workflow_depend)=="trendblind",colnames(logfile$workflow_depend)=="profblind"]==0
){
workflow_depend<-read.table(
file="workflow_depend"
)
workflow_depend<-as.matrix(workflow_depend)
workflow_must<-read.table(
file="workflow_must"
)
workflow_must<-as.matrix(workflow_must)
logfile[[11]]<<-workflow_depend
names(logfile)[11]<<-"workflow_depend"
logfile[[12]]<<-workflow_must
names(logfile)[12]<<-"workflow_must"
old_Tasks_to_redo<-logfile$Tasks_to_redo
logfile[[2]]<<-rep(FALSE,length(colnames(workflow_must)));
names(logfile[[2]])<<-colnames(workflow_must)
names(logfile)[2]<<-c("Tasks_to_redo");
for(i in 1:length(old_Tasks_to_redo)){
if(any(names(logfile$Tasks_to_redo)==names(old_Tasks_to_redo)[i])){
logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)==names(old_Tasks_to_redo)[i]]<<-old_Tasks_to_redo[i]
}
}
old_workflow<-logfile$workflow
logfile$workflow<<-0
names(logfile)[6]<<-c("workflow")
for(i in 1:length(names(logfile[[2]]))){
logfile$workflow[i]<<-"yes";
names(logfile$workflow)[i]<<-names(logfile[[2]])[i]
}
for(i in 1:length(old_workflow)){
if(any(names(logfile$workflow)==names(old_workflow)[i])){
logfile$workflow[names(logfile$workflow)==names(old_workflow)[i]]<<-old_workflow[i]
}
}
old_summary<-logfile$summary
tasks<-names(logfile[[2]])
doneit<-rep(FALSE,length(tasks))
summar<-data.frame(tasks,doneit,stringsAsFactors = FALSE)
names(summar)<-c("Tasks","Done?")
logfile[[3]]<<-summar
names(logfile)[3]<<-c("summary")
for(i in 1:length(old_summary[,1])){
if(any(logfile$summary[,1]==as.character(old_summary[i,1]))){
logfile$summary[logfile$summary[,1]==as.character(old_summary[i,1]),2]<<-as.character(old_summary[i,2])
}
}
schedule<-enviMass::PFMSW1067M(logfile$workflow_depend,logfile$workflow_must)
set_order<-match(schedule[,1],logfile$summary[,1])
logfile$summary<<-logfile$summary[set_order,]
}
if(logfile$parameters$blind_omit=="yes"){
logfile$parameters$blind_omit<<-"TRUE"
}
if(logfile$parameters$blind_omit=="no"){
logfile$parameters$blind_omit<<-"FALSE"
}
if(!any(names(logfile$parameters)=="dofile_latest_profcomp")){
logfile$parameters$dofile_latest_profcomp<<-"FALSE"
logfile$parameters$numfile_latest_profcomp<<-"100"
logfile$parameters$filter_profcomp_pos<<-"TRUE"
logfile$parameters$filter_profcomp_neg<<-"TRUE"
logfile$parameters$for_which_profcomp_pos<<-"all"
logfile$parameters$for_which_profcomp_neg<<-"all"
logfile$parameters$prof_comp_link_only<<-"FALSE"
logfile$parameters$corr_min_peaks<<-"5"
logfile$parameters$comp_corr<<-"0.9"
logfile$parameters$corr_del_RT<<-"5"
logfile$parameters$corr_skip_peaks<<-"TRUE"
}
logfile$version<<-3.29
save(logfile,file=file.path(as.character(logfile[["project_folder"]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(logfile$version<3.3){
cat("\n Updating to version 3.3 ...")
if(file.exists(file=file.path(as.character(logfile[[1]]),"results","profpeaks_pos"))){
file.remove(file=file.path(as.character(logfile[[1]]),"results","profpeaks_pos"))
}
if(file.exists(file=file.path(as.character(logfile[[1]]),"results","profpeaks_neg"))){
file.remove(file=file.path(as.character(logfile[[1]]),"results","profpeaks_neg"))
}
enviMass::RNJJT9347D(
down="profiling",
except=FALSE,
down_TF=c("TRUE","FALSE"),
check_node=TRUE,
single_file=FALSE
)
output$summa_html <- renderText(enviMass::UCBRR9756L(logfile$summary, logfile$Tasks_to_redo));
shinyjs::info(paste0("Profile list structure has been modified for profile componentization &  new filtering functionalities - please press the Calculate button any time soon to make these changes permanent to your project results (entails a project recalculation except peakpicking)!"));
logfile$version<<-3.3
save(logfile,file=file.path(as.character(logfile[["project_folder"]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(logfile$version<3.31){
cat("\n Updating to version 3.31 ...")
if(!any(names(logfile$parameters)=="homol_blind_value")){
logfile$parameters$homol_blind<<-"FALSE"
logfile$parameters$homol_blind_value<<-"10"
}
logfile$version<<-3.31
save(logfile,file=file.path(as.character(logfile[["project_folder"]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(logfile$version<3.311){
cat("\n Updating to version 3.311 ...")
if(
logfile$workflow_depend[rownames(logfile$workflow_depend)=="homologues",colnames(logfile$workflow_depend)=="target_screen"]==0
){
workflow_depend<-read.table(file="workflow_depend")
workflow_depend<-as.matrix(workflow_depend)
workflow_must<-read.table(file="workflow_must")
workflow_must<-as.matrix(workflow_must)
logfile[[11]]<<-workflow_depend
names(logfile)[11]<<-"workflow_depend"
logfile[[12]]<<-workflow_must
names(logfile)[12]<<-"workflow_must"
old_Tasks_to_redo<-logfile$Tasks_to_redo
logfile[[2]]<<-rep(FALSE,length(colnames(workflow_must)));
names(logfile[[2]])<<-colnames(workflow_must)
names(logfile)[2]<<-c("Tasks_to_redo");
for(i in 1:length(old_Tasks_to_redo)){
if(any(names(logfile$Tasks_to_redo)==names(old_Tasks_to_redo)[i])){
logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)==names(old_Tasks_to_redo)[i]]<<-old_Tasks_to_redo[i]
}
}
old_workflow<-logfile$workflow
logfile$workflow<<-0
names(logfile)[6]<<-c("workflow")
for(i in 1:length(names(logfile[[2]]))){
logfile$workflow[i]<<-"yes";
names(logfile$workflow)[i]<<-names(logfile[[2]])[i]
}
for(i in 1:length(old_workflow)){
if(any(names(logfile$workflow)==names(old_workflow)[i])){
logfile$workflow[names(logfile$workflow)==names(old_workflow)[i]]<<-old_workflow[i]
}
}
old_summary<-logfile$summary
tasks<-names(logfile[[2]])
doneit<-rep(FALSE,length(tasks))
summar<-data.frame(tasks,doneit,stringsAsFactors = FALSE)
names(summar)<-c("Tasks","Done?")
logfile[[3]]<<-summar
names(logfile)[3]<<-c("summary")
for(i in 1:length(old_summary[,1])){
if(any(logfile$summary[,1]==as.character(old_summary[i,1]))){
logfile$summary[logfile$summary[,1]==as.character(old_summary[i,1]),2]<<-as.character(old_summary[i,2])
}
}
schedule<-enviMass::PFMSW1067M(logfile$workflow_depend,logfile$workflow_must)
set_order<-match(schedule[,1],logfile$summary[,1])
logfile$summary<<-logfile$summary[set_order,]
shinyjs::info(paste0("Screening results and homologues series extraction are now intersected - please press Calculate a.s.a.p. to apply this change to your project!"));
}
if(logfile$workflow[names(logfile$workflow)=="blind"]=="yes"){
enviMass::RNJJT9347D(
down="blind",
except=FALSE,
down_TF=c("TRUE","FALSE"),
check_node=TRUE,
single_file=FALSE
)
output$summa_html <- renderText(enviMass::UCBRR9756L(logfile$summary, logfile$Tasks_to_redo));
shinyjs::info(paste0("Important patch in blind subtraction/annotation step - please recalculate your project!"));
}
logfile$version<<-3.311
save(logfile,file=file.path(as.character(logfile[["project_folder"]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(logfile$version < 3.4){
cat("\n Updating to version 3.4 ...")
if(
(logfile$workflow_depend[rownames(logfile$workflow_depend) == "components_profiles", colnames(logfile$workflow_depend) == "IS_normaliz"] == 0 ) ||
(logfile$workflow_depend[rownames(logfile$workflow_depend) == "IS_normaliz", colnames(logfile$workflow_depend) == "IS_screen"] == 0 ) ||
(logfile$workflow_depend[rownames(logfile$workflow_depend) == "calibration", colnames(logfile$workflow_depend) == "IS_normaliz"] == 0 )
){
workflow_depend<-read.table(
file="workflow_depend"
)
workflow_depend<-as.matrix(workflow_depend)
workflow_must<-read.table(
file="workflow_must"
)
workflow_must<-as.matrix(workflow_must)
logfile[[11]]<<-workflow_depend
names(logfile)[11]<<-"workflow_depend"
logfile[[12]]<<-workflow_must
names(logfile)[12]<<-"workflow_must"
old_Tasks_to_redo<-logfile$Tasks_to_redo
logfile[[2]]<<-rep(FALSE,length(colnames(workflow_must)));
names(logfile[[2]])<<-colnames(workflow_must)
names(logfile)[2]<<-c("Tasks_to_redo");
for(i in 1:length(old_Tasks_to_redo)){
if(any(names(logfile$Tasks_to_redo)==names(old_Tasks_to_redo)[i])){
logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)==names(old_Tasks_to_redo)[i]]<<-old_Tasks_to_redo[i]
}
}
old_workflow<-logfile$workflow
logfile$workflow<<-0
names(logfile)[6]<<-c("workflow")
for(i in 1:length(names(logfile[[2]]))){
logfile$workflow[i]<<-"yes";
names(logfile$workflow)[i]<<-names(logfile[[2]])[i]
}
for(i in 1:length(old_workflow)){
if(any(names(logfile$workflow)==names(old_workflow)[i])){
logfile$workflow[names(logfile$workflow)==names(old_workflow)[i]]<<-old_workflow[i]
}
}
old_summary<-logfile$summary
tasks<-names(logfile[[2]])
doneit<-rep(FALSE,length(tasks))
summar<-data.frame(tasks,doneit,stringsAsFactors = FALSE)
names(summar)<-c("Tasks","Done?")
logfile[[3]]<<-summar
names(logfile)[3]<<-c("summary")
for(i in 1:length(old_summary[,1])){
if(any(logfile$summary[,1]==as.character(old_summary[i,1]))){
logfile$summary[logfile$summary[,1]==as.character(old_summary[i,1]),2]<<-as.character(old_summary[i,2])
}
}
schedule<-enviMass::PFMSW1067M(logfile$workflow_depend,logfile$workflow_must)
set_order<-match(schedule[,1],logfile$summary[,1])
logfile$summary<<-logfile$summary[set_order,]
enviMass::RNJJT9347D(
down="pattern",
except=FALSE,
down_TF=c("TRUE","FALSE"),
check_node=TRUE,
single_file=FALSE
)
output$summa_html <- renderText(enviMass::UCBRR9756L(logfile$summary, logfile$Tasks_to_redo));
shinyjs::info(paste0("Peak intensity distribution checks now run as default - please recalculate your project!"));
}
if(!any(names(logfile$parameters)=="parallel")){
logfile$parameters$parallel<<-"TRUE"
logfile$parameters$parallel_restrict<<-"FALSE"
logfile$parameters$parallel_cores<<-"8"
logfile$parameters$parallel_x1<<-"FALSE"
logfile$parameters$parallel_x2<<-"FALSE"
logfile$parameters$parallel_x3<<-"FALSE"
}
if(!any(names(logfile$parameters)=="verbose")){
logfile$parameters$verbose<<-"TRUE"
}
if(!any(names(logfile) == "isotopes")){
logfile[[9]] <<- "";
names(logfile)[9] <<- c("isotopes")
}
logfile$version <<- 3.4
save(logfile,file=file.path(as.character(logfile[["project_folder"]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(logfile$version < 3.41){
cat("\n Updating to version 3.41 ...")
logfile$workflow[names(logfile$workflow)=="qc"] <<- "yes"
if(!any(names(logfile$parameters) == "test")){
logfile$parameters$test <<- "FALSE"
}
if(!any(names(logfile$parameters) == "ISnorm_include_pos")){
logfile$parameters$ISnorm_include_pos <<-"TRUE"
logfile$parameters$ISnorm_include_neg <<-"TRUE"
logfile$parameters$ISnorm_percfiles_pos <<- logfile$parameters$ISnorm_percfiles
logfile$parameters$ISnorm_numbIS_pos <<- logfile$parameters$ISnorm_numbIS
logfile$parameters$ISnorm_medblank_pos <<- logfile$parameters$ISnorm_medblank
logfile$parameters$ISnorm_usesubblank_pos <<- logfile$parameters$ISnorm_usesubblank
logfile$parameters$ISnorm_numblank_pos <<- logfile$parameters$ISnorm_numblank
logfile$parameters$ISnorm_medsam_pos <<- logfile$parameters$ISnorm_medsam
logfile$parameters$ISnorm_usesubsam_pos <<- logfile$parameters$ISnorm_usesubsam
logfile$parameters$ISnorm_numsam_pos <<- logfile$parameters$ISnorm_numsam
logfile$parameters$ISnorm_score_pos <<- logfile$parameters$ISnorm_score
logfile$parameters$ISnorm_percfiles_neg <<- logfile$parameters$ISnorm_percfiles
logfile$parameters$ISnorm_numbIS_neg <<- logfile$parameters$ISnorm_numbIS
logfile$parameters$ISnorm_medblank_neg <<- logfile$parameters$ISnorm_medblank
logfile$parameters$ISnorm_usesubblank_neg <<- logfile$parameters$ISnorm_usesubblank
logfile$parameters$ISnorm_numblank_neg <<- logfile$parameters$ISnorm_numblank
logfile$parameters$ISnorm_medsam_neg <<- logfile$parameters$ISnorm_medsam
logfile$parameters$ISnorm_usesubsam_neg <<- logfile$parameters$ISnorm_usesubsam
logfile$parameters$ISnorm_numsam_neg <<- logfile$parameters$ISnorm_numsam
logfile$parameters$ISnorm_score_neg <<- logfile$parameters$ISnorm_score
}
if(
logfile$workflow_depend[rownames(logfile$workflow_depend) == "components_files", colnames(logfile$workflow_depend) == "IS_screen"] == 1
){
workflow_depend <- read.table(
file="workflow_depend"
)
workflow_depend <- as.matrix(workflow_depend)
logfile[[11]] <<- workflow_depend
names(logfile)[11] <<- "workflow_depend"
}
logfile$version<<-3.41
save(logfile,file=file.path(as.character(logfile[["project_folder"]]),"logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"),envir=as.environment(".GlobalEnv"))
}
if(logfile$version < 3.411){
cat("\n Updating to version 3.411 ...")
if(!any(names(logfile$parameters) == "peak_get_mass")){
logfile$parameters$peak_get_mass <<- "mean"
}
logfile$version<<-3.411
save(logfile,file = file.path(as.character(logfile[["project_folder"]]), "logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"), envir = as.environment(".GlobalEnv"))
}
if(logfile$version < 3.412){
cat("\n Updating to version 3.412 ...")
enviMass::RNJJT9347D(
down="pattern",
except=FALSE,
down_TF=c("TRUE","FALSE"),
check_node=TRUE,
single_file=FALSE
)
output$summa_html <- renderText(enviMass::UCBRR9756L(logfile$summary, logfile$Tasks_to_redo));
shinyjs::info(paste0("Please recalculate your project by pressing the Calculate button."));
logfile$version <<- 3.412
save(logfile, file = file.path(as.character(logfile[["project_folder"]]), "logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"), envir = as.environment(".GlobalEnv"))
}
if(logfile$version < 3.413){
cat("\n Updating to version 3.413 ...")
logfile$version <<- 3.413
save(logfile, file = file.path(as.character(logfile[["project_folder"]]), "logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"), envir = as.environment(".GlobalEnv"))
}
if(logfile$version < 3.414){
cat("\n Updating to version 3.414 ...")
logfile$version <<- 3.414
save(logfile, file = file.path(as.character(logfile[["project_folder"]]), "logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"), envir = as.environment(".GlobalEnv"))
}
if(logfile$version < 3.4151){
cat("\n Updating to version 3.4151 ...")
if(!any(names(logfile) == "UI_options")){
logfile[[17]] <<- list(0)
names(logfile)[17] <<- c("UI_options")
logfile$UI_options$filterProf_minmass <<- "0"
logfile$UI_options$filterProf_maxmass <<- "3000"
logfile$UI_options$filterProf_minrt <<- "0"
logfile$UI_options$filterProf_maxrt <<- "100000"
logfile$UI_options$filterProf_minMD <<- "-0.5"
logfile$UI_options$filterProf_maxMD <<- "0.5"
logfile$UI_options$filterProf_medianblind <<- "yes"
logfile$UI_options$filterProf_medianblind_value <<- "10"
logfile$UI_options$filterProf_notblind <<- "no"
logfile$UI_options$filterProf_sort <<- "current trend intensity (decreasing)"
logfile$UI_options$filterProf_components <<- "TRUE"
}
logfile$version <<- 3.4151
save(logfile, file = file.path(as.character(logfile[["project_folder"]]), "logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"), envir = as.environment(".GlobalEnv"))
}
if(logfile$version < 3.43){
cat("\n Updating to version 3.43 ...")
if(!any(names(logfile$parameters) == "do_atom_bounds_components")){
logfile$parameters$do_atom_bounds_components <<- "FALSE"
logfile$parameters$atom_bounds_components <<- c("Cl","Br")
}
logfile$version <<- 3.43
save(logfile, file = file.path(as.character(logfile[["project_folder"]]), "logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"), envir = as.environment(".GlobalEnv"))
}
if(logfile$version < 3.431){
cat("\n Updating to version 3.431 ...")
if(
(logfile$workflow[names(logfile$workflow)=="components_profiles"]=="yes") &
(logfile$parameters$parallel == "FALSE")
){
logfile$Tasks_to_redo[["IS_screen"]] <<- "TRUE"
logfile$Tasks_to_redo[["target_screen"]] <<- "TRUE"
logfile$Tasks_to_redo[["components_profiles"]] <<- "TRUE"
shinyjs::info(paste0("Update requires a recalculation of some workflow nodes - please hit the Calculate button as soon as possible."));
}
logfile$version <<- 3.431
save(logfile, file = file.path(as.character(logfile[["project_folder"]]), "logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"), envir = as.environment(".GlobalEnv"))
}
if(logfile$version < 3.44){
cat("\n Updating to version 3.44 ...")
if(!any(names(logfile$parameters) == "quant_digits")){
logfile$parameters$quant_digits <<- "2"
}
if(!any(names(logfile$parameters) == "replicates_mean_prof")){
logfile$parameters$replicates_mean_prof <<- "TRUE";
}
logfile$version <<- 3.44
save(logfile, file = file.path(as.character(logfile[["project_folder"]]), "logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"), envir = as.environment(".GlobalEnv"))
}
if(logfile$version < 3.45){
cat("\n Updating to version 3.45 ...")
if(!file.exists(file.path(logfile$project_folder,"MSraw"))){
dir.create(file.path(logfile$project_folder,"MSraw"), recursive = TRUE)
}
logfile$version <<- 3.45
save(logfile, file = file.path(as.character(logfile[["project_folder"]]), "logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"), envir = as.environment(".GlobalEnv"))
}
if(logfile$version < 3.453){
cat("\n Updating to version 3.453 ...")
if(!any(names(logfile$parameters) == "files_SIM")){
logfile$parameters$files_SIM <<- "FALSE"
}
if(!any(names(logfile$parameters) == "method_use")){
logfile$parameters$method_use <<- "FALSE"
}
if(!any(names(logfile) == "method_setup")){
logfile[[18]] <<- "Not available"
names(logfile)[18] <<- c("method_setup")
}
if(!any(names(logfile$parameters) == "replicate_mean_profiles")){
logfile$parameters$replicate_mean_profiles <<- "FALSE"
}
logfile$version <<- 3.453
save(logfile, file = file.path(as.character(logfile[["project_folder"]]), "logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"), envir = as.environment(".GlobalEnv"))
}
if(logfile$version < 3.47){
cat("\n Updating to versions 3.46 / 3.47 ...")
if(!any(names(logfile) == "comparisons")){
logfile[[19]] <<- list()
names(logfile)[19] <<- "comparisons"
}
workflow_depend <- read.table(file = "workflow_depend")
workflow_depend <- as.matrix(workflow_depend)
workflow_must <- read.table(file = "workflow_must")
workflow_must <- as.matrix(workflow_must)
logfile[[11]] <<- workflow_depend
names(logfile)[11] <<- "workflow_depend"
logfile[[12]] <<- workflow_must
names(logfile)[12] <<- "workflow_must"
old_Tasks_to_redo <- logfile$Tasks_to_redo
logfile[[2]] <<- rep(FALSE, length(colnames(workflow_must)));
names(logfile[[2]]) <<- colnames(workflow_must)
names(logfile)[2] <<- c("Tasks_to_redo");
for(i in 1:length(old_Tasks_to_redo)){
if(any(names(logfile$Tasks_to_redo) == names(old_Tasks_to_redo)[i])){
logfile$Tasks_to_redo[names(logfile$Tasks_to_redo) == names(old_Tasks_to_redo)[i]] <<- old_Tasks_to_redo[i]
}
}
old_workflow <- logfile$workflow
logfile$workflow <<- 0
names(logfile)[6] <<- c("workflow")
for(i in 1:length(names(logfile[[2]]))){
logfile$workflow[i] <<- "no";
names(logfile$workflow)[i] <<- names(logfile[[2]])[i]
}
for(i in 1:length(old_workflow)){
if(any(names(logfile$workflow) == names(old_workflow)[i])){
logfile$workflow[names(logfile$workflow) == names(old_workflow)[i]] <<- old_workflow[i]
}
}
old_summary <- logfile$summary
tasks <- names(logfile[[2]])
doneit <- rep(FALSE,length(tasks))
summar <- data.frame(tasks, doneit, stringsAsFactors = FALSE)
names(summar) <- c("Tasks","Done?")
logfile[[3]] <<- summar
names(logfile)[3] <<- c("summary")
for(i in 1:length(old_summary[,1])){
if(any(logfile$summary[,1] == as.character(old_summary[i,1]))){
logfile$summary[logfile$summary[,1] == as.character(old_summary[i,1]),2] <<- as.character(old_summary[i,2])
}
}
schedule <- enviMass::PFMSW1067M(logfile$workflow_depend,logfile$workflow_must)
set_order <- match(schedule[,1],logfile$summary[,1])
logfile$summary <<- logfile$summary[set_order,]
logfile$version <<- 3.47
save(logfile, file = file.path(as.character(logfile[["project_folder"]]), "logfile.emp"));
load(file.path(logfile$project_folder,"logfile.emp"), envir = as.environment(".GlobalEnv"))
}
if(any(ls()=="logfile")){stop("\n illegal logfile detected _2 in server_updates.r!")}
