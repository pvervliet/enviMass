	##################################################################################
	# CLEAN!
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
		if(isolate(input$Ion_mode)=="positive"){
			profileList<<-profileList_pos;
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
		if(isolate(input$Ion_mode)=="negative"){
			profileList<<-profileList_neg;
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