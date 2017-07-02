	##################################################################################
	# CLEAN!
	if(any(objects(envir=as.environment(".GlobalEnv"))=="profpeaks_pos")){rm(profpeaks_pos,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="profpeaks_pos")){rm(profpeaks_pos)}		
	if(any(objects(envir=as.environment(".GlobalEnv"))=="profpeaks_neg")){rm(profpeaks_neg,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="profpeaks_neg")){rm(profpeaks_neg)}		
	if(file.exists(file.path(as.character(logfile[[1]]),"results","profpeaks_pos"))){file.remove(file.path(as.character(logfile[[1]]),"results","profpeaks_pos"))}
	if(file.exists(file.path(as.character(logfile[[1]]),"results","profpeaks_neg"))){file.remove(file.path(as.character(logfile[[1]]),"results","profpeaks_neg"))}
	##################################################################################	

	##################################################################################		
	# POSITIVE #######################################################################
    if( file.exists(file.path(logfile[[1]],"results","profileList_pos")) ){
		if(any(objects(envir=as.environment(".GlobalEnv"))=="profileList_pos")){rm(profileList_pos,envir=as.environment(".GlobalEnv"))}
		if(any(objects()=="profileList_pos")){rm(profileList_pos)}				
		load(file.path(as.character(logfile[[1]]),"results","profileList_pos"),envir=as.environment(".GlobalEnv"))
		profileList_pos<-intensup(
				profileList=profileList_pos,
				from=FALSE,
				to=FALSE,
				progbar=logfile$parameters$progressBar,
				blindsub=FALSE,
				blindfold=as.numeric(logfile$parameters$blind_threshold), 	# dummy
				lags=10, 													# dummy
				threshold=3,												# dummy
				notrend=FALSE,												# dummy
				omit_trend=TRUE	# OMITTING TREND DETECTION!
		)
		profileList_pos<<-profileList_pos
		save(profileList_pos,file=file.path(as.character(logfile[[1]]),"results","profileList_pos"),compress=FALSE);
		profpeaks_pos<-enviMass:::profiletopeak(profileList_pos,progbar=logfile$parameters$progressBar)
		profpeaks_pos<-profpeaks_pos[order(profpeaks_pos[,13],decreasing=TRUE),]
		profpeaks_pos<<-profpeaks_pos;
		save(profpeaks_pos,file=file.path(as.character(logfile[[1]]),"results","profpeaks_pos"));
		if(isolate(input$Ion_mode)=="positive"){
			profileList<<-profileList_pos;
			profpeaks<<-profpeaks_pos;
		}
	}
	##################################################################################
		
	##################################################################################
	# NEGATIVE #######################################################################	
    if( file.exists(file.path(logfile[[1]],"results","profileList_neg")) ){
		if(any(objects(envir=as.environment(".GlobalEnv"))=="profileList_neg")){rm(profileList_neg,envir=as.environment(".GlobalEnv"))}
		if(any(objects()=="profileList_neg")){rm(profileList_neg)}				
		load(file.path(as.character(logfile[[1]]),"results","profileList_neg"),envir=as.environment(".GlobalEnv"))
		profileList_neg<-intensup(
				profileList=profileList_neg,
				from=FALSE,
				to=FALSE,
				progbar=logfile$parameters$progressBar,
				blindsub=FALSE,
				blindfold=as.numeric(logfile$parameters$blind_threshold), 	# dummy
				lags=10, 													# dummy
				threshold=3,												# dummy
				notrend=FALSE,												# dummy
				omit_trend=TRUE	# OMITTING TREND DETECTION!
		)
		profileList_neg<<-profileList_neg
		save(profileList_neg,file=file.path(as.character(logfile[[1]]),"results","profileList_neg"),compress=FALSE);
		profpeaks_neg<-enviMass:::profiletopeak(profileList_neg,progbar=logfile$parameters$progressBar)
		profpeaks_neg<-profpeaks_neg[order(profpeaks_neg[,13],decreasing=TRUE),]
		profpeaks_neg<<-profpeaks_neg;
		save(profpeaks_neg,file=file.path(as.character(logfile[[1]]),"results","profpeaks_neg"));
		if(isolate(input$Ion_mode)=="negative"){
			profileList<<-profileList_neg;
			profpeaks<<-profpeaks_neg;
		}
	}
	##################################################################################	

	##################################################################################
	path=file.path(logfile[[1]],"pics","boxprofile_pos")
		png(filename = path, bg = "white")
		plot.new();plot.window(xlim=c(0,1),ylim=c(0,1));text(0.5,0.5,"nothing selected \n or not available",cex=1)
	dev.off()
	expr4p<-list(src=path)
	output$boxprofile<-renderImage(expr4p, deleteFile = FALSE)
	path=file.path(logfile[[1]],"pics","boxprofile_neg")
		png(filename = path, bg = "white")
		plot.new();plot.window(xlim=c(0,1),ylim=c(0,1));text(0.5,0.5,"nothing selected \n or not available",cex=1)
	dev.off()
	expr4n<-list(src=path)
	output$boxprofile<-renderImage(expr4n, deleteFile = FALSE)			
	##################################################################################