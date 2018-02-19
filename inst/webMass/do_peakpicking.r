# Peak picking ##################################################################

	output$dowhat<-renderText("Peak picking ... please wait");
	if(any(search()=="package:nlme")){detach(package:nlme,force=TRUE);addit<-TRUE}else{addit<-FALSE}
    measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
    leng<-dim(measurements)[1];         
	if(logfile$parameters$cut_RT=="TRUE"){
		use_minRT<-(as.numeric(logfile$parameters$cut_RT_min)*60)
		use_maxRT<-(as.numeric(logfile$parameters$cut_RT_max)*60)
		cat("(filter RT range)")
	}else{
		use_minRT<-FALSE
		use_maxRT<-FALSE				
	}		
	if(logfile$parameters$cut_mass=="TRUE"){
		use_minmass<-as.numeric(logfile$parameters$cut_mass_min)
		use_maxmass<-as.numeric(logfile$parameters$cut_mass_max)
		cat("(filter mass range)")
	}else{
		use_minmass<-FALSE
		use_maxmass<-FALSE				
	}				
				
	#i <- which(measurements[,"ID"] == "11") 
    for(i in 1:leng){ 
        # (measurement included & not yet picked) OR (peakpick forced) 
            if( (measurements[i,"include"]=="TRUE") & (measurements[i,"peakpicking"]=="FALSE") ){

				##################################################################
				cat(paste("\n    Peak picking sample ",as.character(i)," of ",as.character(leng),": "));    
				##################################################################
				if(#logfile$parameters$is_example=="FALSE"
					file.exists(file.path(logfile[[1]],"files",paste0(as.character(measurements[i,"ID"]),".mzXML")))
				){
					##############################################################					
					if(  as.logical(logfile$parameters$method_use) & !is.character(logfile$method_setup) ){	# method setup, mzR	
					
						##########################################################
						# based on method setup, retrieve scanTypes ##############
						# unused scans -> set index to NA ########################
						polar <- c("-", "+")
						method_setup <- logfile$method_setup[,, drop = FALSE]
						omit_col <- c("Scan type", "Scan counts", "Centroid counts", "used_scans")
						if(any(names(method_setup) == "MS1_consec")) omit_col <- c(omit_col, "MS1_consec")
						omit_col <- match(omit_col, names(method_setup))
						heads_summary <- method_setup[, -omit_col, drop = FALSE] 
						if(any(names(heads_summary) == "polarity")){
							heads_summary$polarity[heads_summary$polarity == "+"] <- "1"
							heads_summary$polarity[heads_summary$polarity == "-"] <- "0"
							heads_summary$polarity <- as.numeric(heads_summary$polarity)
						}
						path <- file.path(logfile[[1]], "files", paste0(as.character(measurements[i,1]), ".mzXML"))
						# path <- "D:/Projects/Sabine/new_project_name/files/2.mzXML"
						mzXML_file <- mzR:::openMSfile(filename = path, backend = c("Ramp"), verbose = FALSE)
						heads <- mzR::header(mzXML_file)
						get_heads <- match(names(heads_summary), names(heads)) # resort columns
						if(any(is.na(get_heads))){
							which_not_header <- names(heads_summary)[which(any(is.na(match(names(heads_summary), names(heads)))))][1]
							stop(paste("\n", which_not_header, " from the method definition not found in the file header - revise settings! Files from different experiments?"))
						}
						heads_short <- heads[,get_heads, drop = FALSE]
						scanTypes2 <- method_setup[,"Scan type"]			
						scanTypes <- match(data.frame(t(heads_short)), data.frame(t(heads_summary))) # convert to list of vectors
						if(any(is.na(scanTypes))) stop(
							paste("\n Problem with File ID", measurements[i,"ID"], "while applying Method setup - one or several scans could not be matched to the existing setup. Maybe the setup was defined for a different file type?")
						)
						if(
							any(names(method_setup) == "msLevel") & 
							any(names(method_setup) == "MS1_consec")
						){ # consecutive MS1 - scan separation
							if(
								any(heads_summary$msLevel == 1) &
								any(heads_short[,"msLevel"] == 1) &
								any(method_setup$MS1_consec > 0)
							){
								MS1_consec <- rep(0, length(scanTypes))
								for(n in 2:length(scanTypes)){ # collect additional levels from consecutive scans
									if(heads_summary[scanTypes[n], "msLevel"] != 1) next
									if(heads_summary[scanTypes[n - 1], "msLevel"] != 1) next
									if(scanTypes[n] != scanTypes[n - 1]) next
									MS1_consec[n] <- (MS1_consec[n - 1] + 1)
								}
								if(any(MS1_consec > 0)){ # consecutive MS1 scans of otherwise same definition found?
									scanTypes <- match(data.frame(t(cbind(heads_short, MS1_consec))), data.frame(t(cbind(heads_summary, method_setup$MS1_consec)))) # convert to list of vectors
								}
							}
						}
						if(any(is.na(logfile$method_setup$used_scans[scanTypes]))) stop("\n Unmatched Scan types detected in the methods setup - report this issue.")
						scanTypes[!logfile$method_setup$used_scans[scanTypes]] <- NA
						##########################################################
						# read data ##############################################								
						# read msLevel 1 data first 
						if(!any(names(method_setup) == "msLevel")) stop("\n No msLevel 1 scans remaining from method definition - ate least one such Scan type required to run the workflow. Please revise!")
						if(!any(method_setup$msLevel == 1)) stop("\n No msLevel 1 scans remaining from method definition - ate least one such Scan type required to run the workflow. Please revise!")						
						read_scanType <- unique(scanTypes)
						read_scanType <- read_scanType[!is.na(read_scanType)]
						read_scanType <- read_scanType[logfile$method_setup[read_scanType, "msLevel"] == 1]
						if(!length(read_scanType)) stop("No scanType remaining for msLevel 1 - is the method setup correct? Please revise.")
						MSlist <- enviMass:::convert_mzXML_MSlist(
							mzXML_file,
							scanTypes,
							read_scanType,						
							minRT = use_minRT,
							maxRT = use_maxRT,
							minmz = use_minmass,
							maxmz = use_maxmass,						
							get_acquisitionNum = FALSE			
						)
						##########################################################
					}else{ # default: all MS1 upload		
						MSlist <- enviPick::readMSdata(
							filepath.mzXML = file.path(logfile[[1]], "files", paste0(as.character(measurements[i,1]), ".mzXML")),
							MSlevel = logfile$parameters$peak_MSlevel,  # MSlevel
							progbar = logfile$parameters$progressBar, # progbar
							minRT = use_minRT,
							maxRT = use_maxRT,
							minmz = use_minmass,
							maxmz = use_maxmass,
							ion_mode = FALSE#measurements[i,"Mode"]
						);
					};					
					cat(" file read -"); 
					##############################################################
					if(logfile$parameters$peak_estimate=="TRUE"){
						use_peak_perc_cut<-0
						estim_values<-try({enviMass::dens_filter(MSlist,plotit=FALSE,n=2000,m=5)},silent=TRUE)
						if(class(estim_values)!="try-error"){
							use_peak_dmzdens<-estim_values[[1]]
							use_peak_minint_log10<-estim_values[[2]]
							if(as.numeric(logfile$parameters$peak_maxint_log10)<use_peak_minint_log10){
								use_peak_maxint_log10<-log10(max(MSlist[["Scans"]][[2]][,"intensity"])+1)
							}else{
								use_peak_maxint_log10<-as.numeric(logfile$parameters$peak_maxint_log10)
							} 
							len1<-dim(MSlist[["Scans"]][[2]])[1]
							MSlist[["Scans"]][[2]]<-MSlist[["Scans"]][[2]][
								MSlist[["Scans"]][[2]][,"intensity"]>=(10^use_peak_minint_log10)
							,]
							len2<-dim(MSlist[["Scans"]][[2]])[1]
							cat(paste(" ",as.character(len1-len2),"of",as.character(len1),"data points discarded by absolute threshold -"))
							cat(" data filtered -"); 
						}else{
							use_peak_dmzdens<-as.numeric(logfile$parameters$peak_dmzdens)
							use_peak_minint_log10<-as.numeric(logfile$parameters$peak_minint_log10)
							use_peak_maxint_log10<-as.numeric(logfile$parameters$peak_maxint_log10)					
							use_peak_perc_cut<-as.numeric(logfile$parameters$peak_perc_cut)
							cat(" no data filtering possible -"); 
						}	
					}else{
						use_peak_dmzdens <- as.numeric(logfile$parameters$peak_dmzdens)
						use_peak_minint_log10 <- as.numeric(logfile$parameters$peak_minint_log10)
						use_peak_maxint_log10 <- as.numeric(logfile$parameters$peak_maxint_log10)					
						use_peak_perc_cut <- as.numeric(logfile$parameters$peak_perc_cut)
						if(use_peak_perc_cut>0){ # not to be used with filtering estimates - that uses absolute threshold intensities
							len1<-dim(MSlist[["Scans"]][[2]])[1]
							MSlist[["Scans"]][[2]] <- MSlist[["Scans"]][[2]][
								MSlist[["Scans"]][[2]][,"intensity"] >= quantile(MSlist[["Scans"]][[2]][,"intensity"],(use_peak_perc_cut/100))
							,]
							len2<-dim(MSlist[["Scans"]][[2]])[1]
							cat(paste(" ",as.character(len1-len2),"of",as.character(len1),"data points discarded by fraction -"))
						}
					}	
					##############################################################
				}else{ # no mzXML files for example projects -> use MSlist
					load(file = file.path(logfile[[1]],"MSlist",as.character(measurements[i,"ID"]))); 
					MSlist[["Scans"]][[2]] <- MSlist[["Scans"]][[2]][ # re-set
						order(MSlist[["Scans"]][[2]][,"RT"], decreasing = FALSE)
					,]
					MSlist[["Scans"]][[2]][,"partID"] <- 0
					MSlist[["Scans"]][[2]][,"clustID"] <- 0					
					MSlist[["Scans"]][[2]][,"peakID"] <- 0
					use_peak_dmzdens <- as.numeric(logfile$parameters$peak_dmzdens)
					use_peak_minint_log10 <- 0
					use_peak_maxint_log10 <- as.numeric(logfile$parameters$peak_maxint_log10)					
					use_peak_perc_cut <- 0
					cat("example MSlist loaded -"); 
				}
				##################################################################
				if(any(MSlist[["Scans"]][[2]][,"intensity"] == 0)){
					cat("\n Note in do_peakpicking: zero intensities found and discarded.")
					MSlist[["Scans"]][[2]] <- MSlist[["Scans"]][[2]][MSlist[["Scans"]][[2]][,"intensity"]!=0,,drop=FALSE]
				}	
				MSlist <- enviPick::mzagglom(
					MSlist,
					((use_peak_dmzdens*2)+1),
					ppm=TRUE,  
					as.numeric(logfile$parameters$peak_drtgap),
					as.numeric(logfile$parameters$peak_minpeak),
					10^use_peak_maxint_log10,
					progbar=as.logical(logfile$parameters$progressBar)
				);
				cat(" partitioned -");
				##################################################################
				MSlist<-enviPick::mzclust(      
					MSlist,
					use_peak_dmzdens,
					ppm=TRUE,
					60,
					as.numeric(logfile$parameters$peak_minpeak),
					10^use_peak_maxint_log10,
					progbar=as.logical(logfile$parameters$progressBar),
					merged=TRUE,from=FALSE,to=FALSE
				);
				cat(" clustered -");
				##################################################################
				MSlist <- enviPick::mzpick(      
					MSlist = MSlist,
					minpeak = as.numeric(logfile$parameters$peak_minpeak), 
					drtsmall = as.numeric(logfile$parameters$peak_drtsmall2), 
					drtfill = as.numeric(logfile$parameters$peak_drtfill),        
					drttotal = as.numeric(logfile$parameters$peak_drtdens2),  
					recurs = as.numeric(logfile$parameters$peak_recurs),
					weight = as.numeric(logfile$parameters$peak_weight),
					SB = as.numeric(logfile$parameters$peak_SB),
					SN = as.numeric(logfile$parameters$peak_SN),               
					minint = (10^use_peak_minint_log10),
					maxint = (10^use_peak_maxint_log10),
					ended = as.numeric(logfile$parameters$peak_ended),
					get_mass = logfile$parameters$peak_get_mass,
					progbar = logfile$parameters$progressBar,
					from = FALSE,
					to = FALSE
				);
				if(any(MSlist[["Peaklist"]][,3] == 0)){stop("\n do_peakpicking: zero intensities found - resolve issue before proceding.")}
				MSlist[[9]] <- as.character(measurements[i,"ID"]);
				names(MSlist)[9] <- "File_ID";
				save(MSlist, file = file.path(logfile[[1]], "MSlist", as.character(measurements[i,"ID"])));   
				peaklist <- MSlist[["Peaklist"]];
				if(length(peaklist) == 0){
					stop("No peaks picked - wrong parameters (e.g., intensity thresholds too high)?")
				}
				##################################################################
# REMOVE ->		#peaklist<-cbind(peaklist,peaklist[,1],rep(0,length(peaklist[,4])),peaklist[,5])
				peaklist<-cbind(peaklist,
					rep(0,length(peaklist[,4])),
					rep(0,length(peaklist[,4])),
					peaklist[,5] # replace by rep(0) as soon as do_align.r is build!
				)	
				colnames(peaklist)[12] <- "m/z_corr";
				colnames(peaklist)[13] <- "int_corr";
				colnames(peaklist)[14] <- "RT_corr";      
				keep1 <- rep(1,length(peaklist[,1])) 		# replicates, 1 == TRUE
				keep2 <- rep(Inf,length(peaklist[,1])) 	# blind indicators, Inf == not affected
				peaklist <- cbind(peaklist,keep1,keep2) 	
				colnames(peaklist)[15] <- "keep";
				colnames(peaklist)[16] <- "keep_2";
				save(peaklist,file=file.path(logfile[[1]],"peaklist",as.character(measurements[i,"ID"])));   
				cat(" plotted -");  
				path=file.path(logfile[[1]],"pics",paste("peakhist_",as.character(measurements[i,"ID"]),sep=""))
				png(filename = path, bg = "white")    
				a<-hist(log10(MSlist[["Scans"]][[2]][,2]),breaks=200,plot=FALSE)
				hist(log10(MSlist[["Scans"]][[2]][,2]),breaks=a$breaks,
					xlim=c(min(a$breaks),max(a$breaks)),ylim=c(0,max(a$counts)),
					xlab="log10(Intensity)",main="All data points (white histrogram), those in peaks (red) and their count ratio (blue points)",cex.main=.8)
				b<-hist(log10(MSlist[["Scans"]][[2]][MSlist[["Scans"]][[2]][,7]!=0,2]),breaks=a$breaks,col="red",add=TRUE)        
				atfrac<-(b$counts/a$counts)
				plot.window(xlim=c(min(a$breaks),max(a$breaks)),ylim=c(0,1))
				points(a$mids[!is.na(atfrac)],atfrac[!is.na(atfrac)],col="blue",cex=.5,pch=19)
				axis(4)
				dev.off() 				
				path=file.path(logfile[[1]],"pics",paste("peakmzRT_",as.character(measurements[i,1]),sep=""))
				png(filename = path, bg = "white")    
				plot(MSlist[["Peaklist"]][,1],MSlist[["Peaklist"]][,5],xlab="m/z",ylab="RT",pch=19,cex=0.3,main="Picked peaks")
				dev.off()
				if(any(objects(envir=as.environment(".GlobalEnv"))=="MSlist")){rm(MSlist,envir=as.environment(".GlobalEnv"))}
				if(any(objects()=="MSlist")){rm(MSlist)}
				measurements[i,"peakpicking"]<-TRUE;
				write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
				cat(" picked."); 
				##################################################################

            }
    }
	if(addit){library(nlme)}
    cat("Peak picking completed \n"); 	  
    

