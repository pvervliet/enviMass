# Peak picking ##################################################################

	output$dowhat<-renderText("Peak picking ... please wait");
	if(any(search()=="package:nlme")){detach(package:nlme,force=TRUE);addit<-TRUE}else{addit<-FALSE}
    measurements <- read.csv(file = file.path(logfile[[1]], "dataframes", "measurements"), colClasses = "character");
    leng <- dim(measurements)[1];         
	
#measurements[,"peakpicking"]<-"FALSE"
#i<-1

    for(i in 1:leng){ 
            if( (measurements[i,"include"]=="TRUE") & (measurements[i,"peakpicking"]=="FALSE") ){

				##################################################################
				cat(paste("\n    Peak picking sample ",as.character(i)," of ",as.character(leng)," with ID ",measurements[i,"ID"],": "));    
				if(logfile$parameters$cut_RT=="TRUE"){
					use_minRT <- (as.numeric(logfile$parameters$cut_RT_min)*60)
					use_maxRT <- (as.numeric(logfile$parameters$cut_RT_max)*60)
					cat("(filter RT range)")
				}else{
					use_minRT <- FALSE
					use_maxRT <- FALSE				
				}
				if(logfile$parameters$cut_mass == "TRUE"){
					use_minmass <- as.numeric(logfile$parameters$cut_mass_min)
					use_maxmass <- as.numeric(logfile$parameters$cut_mass_max)
					cat("(filter mass range)")
				}else{
					use_minmass <- FALSE
					use_maxmass <- FALSE				
				}				
				##################################################################
				if(#logfile$parameters$is_example=="FALSE"
					file.exists(file.path(logfile[[1]], "files", paste0(as.character(measurements[i,"ID"]), ".mzXML")))
				){
					MSlist<-enviPick::readMSdata(
						filepath.mzXML = file.path(logfile[[1]], "files", paste0(as.character(measurements[i,1]), ".mzXML")),
						MSlevel = logfile$parameters$peak_MSlevel,  # MSlevel
						progbar = logfile$parameters$progressBar, # progbar
						minRT = use_minRT,
						maxRT = use_maxRT,
						minmz = use_minmass,
						maxmz = use_maxmass,
						ion_mode = FALSE#measurements[i,"Mode"]
					);
					cat(" file read -"); 
					##############################################################
					if( logfile$parameters$peak_estimate == "TRUE" ){
						use_peak_perc_cut <- 0
						estim_values <- try({enviMass::dens_filter(MSlist, plotit = FALSE, n = 2000, m = 5)}, silent = TRUE)
						if( (class(estim_values) != "try-error") & (!all(is.na(estim_values))) ){
							use_peak_dmzdens <- estim_values[[1]]
							use_peak_minint_log10 <- estim_values[[2]]
							if(as.numeric(logfile$parameters$peak_maxint_log10) < use_peak_minint_log10){
								use_peak_maxint_log10 <- log10(max(MSlist[["Scans"]][[2]][,"intensity"])+1)
							}else{
								use_peak_maxint_log10 <- as.numeric(logfile$parameters$peak_maxint_log10)
							} 
							len1<-dim(MSlist[["Scans"]][[2]])[1]
							MSlist[["Scans"]][[2]] <- MSlist[["Scans"]][[2]][
								MSlist[["Scans"]][[2]][,"intensity"] >= (10^use_peak_minint_log10)
							,]
							len2 <- dim(MSlist[["Scans"]][[2]])[1]
							cat(paste(" ", as.character(len1-len2), "of", as.character(len1), "data points discarded by absolute threshold -"))
							cat(" data filtered -"); 
						}else{
							use_peak_dmzdens <- as.numeric(logfile$parameters$peak_dmzdens)
							use_peak_minint_log10 <- as.numeric(logfile$parameters$peak_minint_log10)
							use_peak_maxint_log10 <- as.numeric(logfile$parameters$peak_maxint_log10)					
							use_peak_perc_cut <- as.numeric(logfile$parameters$peak_perc_cut)
							cat(" no data filtering possible -"); 
						}	
					}else{
						use_peak_dmzdens <- as.numeric(logfile$parameters$peak_dmzdens)
						use_peak_minint_log10 <- as.numeric(logfile$parameters$peak_minint_log10)
						use_peak_maxint_log10 <- as.numeric(logfile$parameters$peak_maxint_log10)					
						use_peak_perc_cut <- as.numeric(logfile$parameters$peak_perc_cut)
						if(use_peak_perc_cut > 0){ # not to be used with filtering estimates - that uses absolute threshold intensities
							len1<-dim(MSlist[["Scans"]][[2]])[1]
							MSlist[["Scans"]][[2]] <- MSlist[["Scans"]][[2]][
								MSlist[["Scans"]][[2]][,"intensity"] >= quantile(MSlist[["Scans"]][[2]][,"intensity"],(use_peak_perc_cut/100))
							,]
							len2 <- dim(MSlist[["Scans"]][[2]])[1]
							cat(paste(" ",as.character(len1-len2),"of",as.character(len1),"data points discarded by fraction -"))
						}
					}	
					##############################################################
				}else{ # no mzXML files for example projects -> use MSlist
					load(file=file.path(logfile[[1]],"MSlist",as.character(measurements[i,"ID"]))); 
					MSlist[["Scans"]][[2]]<-MSlist[["Scans"]][[2]][ # re-set
						order(MSlist[["Scans"]][[2]][,"RT"],decreasing=FALSE)
					,]
					MSlist[["Scans"]][[2]][,"partID"]<-0
					MSlist[["Scans"]][[2]][,"clustID"]<-0					
					MSlist[["Scans"]][[2]][,"peakID"]<-0
					use_peak_dmzdens<-as.numeric(logfile$parameters$peak_dmzdens)
					use_peak_minint_log10<-0
					use_peak_maxint_log10<-as.numeric(logfile$parameters$peak_maxint_log10)					
					use_peak_perc_cut<-0
					cat("example MSlist loaded -"); 
				}
				##################################################################
				if(any(MSlist[["Scans"]][[2]][,"intensity"]==0)){
					cat("\n Note in do_peakpicking: zero intensities found and discarded.")
					MSlist[["Scans"]][[2]]<-MSlist[["Scans"]][[2]][MSlist[["Scans"]][[2]][,"intensity"]!=0,,drop=FALSE]
				}	
				MSlist<-enviPick:::mzagglom(
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
				MSlist <- enviMass:::mzclust_pl( # no parallelization implemented yet
					MSlist,
					dmzdens = use_peak_dmzdens,
					ppm = TRUE,
					drtdens = 60,
					minpeak = as.numeric(logfile$parameters$peak_minpeak),
					maxint = (10^use_peak_maxint_log10),
					progbar = FALSE,
					merged = TRUE, 
					from = FALSE, 
					to = FALSE
				);
				if(mute(as.logical(logfile$parameters$test))){#if(TRUE){					
					num_EICs <- dim(MSlist[["EIC_index"]])[1]
					for(m in 1:num_EICs){
						if(
							anyDuplicated(MSlist[["Scans"]][[2]][
								(MSlist[["EIC_index"]][m, "start_ID"]):(MSlist[["EIC_index"]][m, "end_ID"]),
							"RT", drop = FALSE]) != 0
						) stop("Test result in do_peakpicking_pl: EIC centroids with non-unique RT found. Debug?")
					}
				}
				cat(" clustered -");				
				##################################################################
				minpeak <- as.numeric(logfile$parameters$peak_minpeak)
                drtsmall <- as.numeric(logfile$parameters$peak_drtsmall2)
                drtfill <- as.numeric(logfile$parameters$peak_drtdens2)
                drttotal <- as.numeric(logfile$parameters$peak_drtdens2)
                recurs <- as.numeric(logfile$parameters$peak_recurs)
                weight <- as.numeric(logfile$parameters$peak_weight)
                SB <- as.numeric(logfile$parameters$peak_SB)
                SN <- as.numeric(logfile$parameters$peak_SN)
                minint <- (10^use_peak_minint_log10)
                maxint <- (10^use_peak_maxint_log10)
                ended <- as.numeric(logfile$parameters$peak_ended)
				Retens <- MSlist[["Scans"]][[1]]
				MSlist[[7]] <- 0;
				MSlist[[8]] <- 0;
				MSlist[[4]][[2]][,7] <- rep(0,length(MSlist[[4]][[2]][,4]));
				####################################################################
				# check inputs #####################################################
				if(minpeak <= 0){stop("minpeak must be >0!")};
				if(drtsmall <= 0 || drtsmall<=0){stop("drt must be >0!")};
				if(!is.loaded("pickpeak")){stop(".Call to pickpeak failed!")};
				if(!is.loaded("gapfill")){stop(".Call to gapfill failed!")};
				if(!is.loaded("picklist")){stop(".Call to picklist failed!")};
				if(minint >= maxint){stop("Revise your minint & maxint settings!")}
				if(ended < 1){stop("Wrong ended argument!")}
				if(!is.integer(recurs)){recurs <- ceiling(recurs); recurs <- as.integer(recurs)}
				if(!is.integer(ended)){ended <- ceiling(ended); ended <- as.integer(ended)}
				if(!is.integer(minpeak)){minpeak <- ceiling(minpeak); minpeak <- as.integer(minpeak)}
				#####################################################################
				num_EICs <- dim(MSlist[["EIC_index"]])[1]
				clus_centroids <- vector("list", num_EICs)
				for(m in 1:num_EICs){
					clus_centroids[[m]] <- 
					MSlist[["Scans"]][[2]][
						(MSlist[["EIC_index"]][m, "start_ID"]):(MSlist[["EIC_index"]][m, "end_ID"]),
						c("m/z", "intensity", "RT", "measureID", "clustID"), drop = FALSE]
				}
				clusterEvalQ(cl = clus,{rm(list = ls()); gc(verbose = FALSE); NULL})
				clusterExport(cl = clus, 
					varlist = c("minpeak", "drtsmall", "drtfill", "drttotal", "recurs", "weight", "SB", "SN", "minint", "maxint", "ended" ,"Retens"), 
					envir = environment())	
				cluster_results <- clusterMap(
					cl = clus,
					fun = enviMass:::mzpick_pl_new,
					clus_centroids = clus_centroids,
					RECYCLE = TRUE,
					SIMPLIFY = FALSE, # no!
					USE.NAMES = FALSE,
					.scheduling = c("dynamic")
				)
				clusterEvalQ(cl = clus,{rm(list = ls()); gc(verbose = FALSE); NULL})
				startat <- 0
				for(m in 1:num_EICs){
					if(all(cluster_results[[m]]==0)) next
					those <- (cluster_results[[m]]!=0)
					cluster_results[[m]][those] <- (cluster_results[[m]][those] + startat)
					at_pos <- (MSlist[["EIC_index"]][m, "start_ID"] : MSlist[["EIC_index"]][m, "end_ID"])
					MSlist[["Scans"]][[2]][at_pos, "peakID"] <- cluster_results[[m]]
					MSlist[["Scans"]][[2]][at_pos,] <- 
						MSlist[["Scans"]][[2]][at_pos,][
							order(MSlist[["Scans"]][[2]][at_pos, "peakID"],decreasing = FALSE),, drop = FALSE
						]
					startat <- max(cluster_results[[m]][those])
				}
				if(as.logical(logfile$parameters$test)){
					IDpeak <- MSlist[["Scans"]][[2]][,"peakID"]
					IDpeak <- IDpeak[IDpeak!=0]
					IDpeak <- unique(IDpeak)	
					if(any(diff(IDpeak)>1)) stop("Test result in do_peakpicking_pl: non-continuous peak IDs. Debug?")	
				}
				# build index #############################################################	
				maxit <- max(MSlist[[4]][[2]][,7]);
				# generate peak ID table ##################################################
				if(maxit > 0){
					index <- .Call("indexed",
						as.integer(MSlist[[4]][[2]][,7]),
						as.numeric(MSlist[[4]][[2]][,2]),
						as.integer(minpeak),
						as.numeric(maxint),
						as.integer(maxit),
						PACKAGE = "enviPick"
					)
					if(any(index[,2] != 0)){
						index <- index[index[,2] != 0,, drop = FALSE];
						partID <- .Call("partID",
							as.integer(index),
							as.integer(length(MSlist[[4]][[2]][,7])),
							PACKAGE = "enviPick"  
						)
						MSlist[[4]][[2]][,7] <- partID
						colnames(index) <- c("start_ID", "end_ID", "number_peaks");
						MSlist[[7]] <- index
					}
				}
				############################################################################
				maxit <- max(MSlist[["Scans"]][[2]][,7]);
				# generate peaklist ########################################################
				if(maxit>0){
					peaklist<-matrix(0, ncol = 11, nrow = maxit)
					colnames(peaklist) <- c("m/z", "var_m/z", "max_int", "sum_int", "RT", "minRT", "maxRT", "part_ID", "EIC_ID", "peak_ID", "Score")
					for(k in 1:length(MSlist[["Peak_index"]][,1])){
						if(logfile$parameters$peak_get_mass == "mean"){
							peaklist[k,1] <- mean(MSlist[["Scans"]][[2]][MSlist[["Peak_index"]][k,1]:MSlist[["Peak_index"]][k,2],1])
						}
						if(logfile$parameters$peak_get_mass == "wmean"){
							peaklist[k,1] <- weighted.mean(
								x = (MSlist[[4]][[2]][MSlist[[7]][k,1]:MSlist[[7]][k,2],1]),
								w = (MSlist[[4]][[2]][MSlist[[7]][k,1]:MSlist[[7]][k,2],2])
							)						
						}
						peaklist[k,2] <- var(MSlist[["Scans"]][[2]][MSlist[["Peak_index"]][k,1]:MSlist[["Peak_index"]][k,2],1])
						peaklist[k,3] <- max(MSlist[["Scans"]][[2]][MSlist[["Peak_index"]][k,1]:MSlist[["Peak_index"]][k,2],2])
						peaklist[k,4] <- sum(MSlist[["Scans"]][[2]][MSlist[["Peak_index"]][k,1]:MSlist[["Peak_index"]][k,2],2])
						peaklist[k,5] <- (MSlist[["Scans"]][[2]][MSlist[["Peak_index"]][k,1]:MSlist[["Peak_index"]][k,2],3][
							MSlist[["Scans"]][[2]][MSlist[["Peak_index"]][k,1]:MSlist[["Peak_index"]][k,2],2] == max(MSlist[["Scans"]][[2]][MSlist[["Peak_index"]][k,1]:MSlist[["Peak_index"]][k,2],2])[1]
						])[1]
						peaklist[k,6] <- min(MSlist[["Scans"]][[2]][MSlist[["Peak_index"]][k,1]:MSlist[["Peak_index"]][k,2],3])
						peaklist[k,7] <- max(MSlist[["Scans"]][[2]][MSlist[["Peak_index"]][k,1]:MSlist[["Peak_index"]][k,2],3])        
						peaklist[k,8:10] <- MSlist[["Scans"]][[2]][MSlist[["Peak_index"]][k,1],5:7]
						peaklist[k,11] <- 0;
					}
					peaklist<-peaklist[order(peaklist[,3], decreasing = TRUE),];
					MSlist[["Peaklist"]] <- peaklist;
					MSlist[["Results"]][[6]] <- length(peaklist[,1])
				}else{
					MSlist[["Results"]][[6]] <- "No peaks picked"
				}
				############################################################################
				MSlist[["State"]][[5]]<-TRUE;
				MSlist[["Results"]][[7]]<-length(MSlist[["Scans"]][[2]][MSlist[["Scans"]][[2]][,7]!=0,2])  
				############################################################################
				MSlist[["Parameters"]][[2]][22] <- as.character(minpeak)
				MSlist[["Parameters"]][[2]][23] <- as.character(drtsmall) 	
				MSlist[["Parameters"]][[2]][24] <- as.character(drtfill) 	
				MSlist[["Parameters"]][[2]][25] <- as.character(drttotal) 	
				MSlist[["Parameters"]][[2]][26] <- as.character(recurs) 	
				MSlist[["Parameters"]][[2]][27] <- as.character(weight) 	
				MSlist[["Parameters"]][[2]][28] <- as.character(SB) 	
				MSlist[["Parameters"]][[2]][29] <- as.character(SN) 	
				MSlist[["Parameters"]][[2]][30] <- as.character(minint) 	
				MSlist[["Parameters"]][[2]][31] <- as.character(maxint) 	
				MSlist[["Parameters"]][[2]][32] <- as.character(ended) 	
				MSlist[["Parameters"]][[2]][33] <- as.character(1) 	
				MSlist[["Parameters"]][[2]][34] <- as.character(length(MSlist[[6]][,1]))			
				############################################################################	
				if(any(MSlist[["Peaklist"]][,3] == 0)){stop("\n do_peakpicking: zero intensities found - resolve issue before proceding.")}
				save(MSlist, file = file.path(logfile[[1]], "MSlist", as.character(measurements[i,"ID"])));   
				peaklist <- MSlist[["Peaklist"]];
				if(length(peaklist) == 0){
					stop("No peaks picked - wrong parameters (e.g., intensity thresholds too high)?")
				}
				##################################################################
				peaklist <- cbind(peaklist,
					rep(0, length(peaklist[,4])),
					rep(0, length(peaklist[,4])),
					peaklist[,5] # replace by rep(0) as soon as do_align.r is build!
				)	
				colnames(peaklist)[12] <- "m/z_corr";
				colnames(peaklist)[13] <- "int_corr";
				colnames(peaklist)[14] <- "RT_corr";      
				keep1 <- rep(1, length(peaklist[,1])) 		# replicates, 1 == TRUE
				keep2 <- rep(Inf, length(peaklist[,1])) 	# blind indicators, Inf == not affected
				peaklist <- cbind(peaklist, keep1, keep2) 	
				colnames(peaklist)[15] <- "keep";
				colnames(peaklist)[16] <- "keep_2";
				save(peaklist, file = file.path(logfile[[1]], "peaklist", as.character(measurements[i,"ID"])));   
				cat(" plotted -");  
				path = file.path(logfile[[1]], "pics", paste0("peakhist_", as.character(measurements[i,"ID"])))
				png(filename = path, bg = "white")    
				a <- hist(log10(MSlist[["Scans"]][[2]][,2]), breaks = 200, plot = FALSE)
				hist(log10(MSlist[["Scans"]][[2]][,2]), breaks = a$breaks,
					xlim = c(min(a$breaks), max(a$breaks)), ylim = c(0,max(a$counts)),
					xlab = "log10(Intensity)",main="All data points (white histrogram), those in peaks (red) and their count ratio (blue points)",cex.main=.8)
				b<-hist(log10(MSlist[["Scans"]][[2]][MSlist[["Scans"]][[2]][,7] != 0, 2]), breaks = a$breaks, col = "red", add = TRUE)        
				atfrac <- (b$counts/a$counts)
				plot.window(xlim = c(min(a$breaks), max(a$breaks)), ylim = c(0, 1))
				points(a$mids[!is.na(atfrac)], atfrac[!is.na(atfrac)], col = "blue", cex = .5, pch = 19)
				axis(4)
				dev.off() 				
				path = file.path(logfile[[1]], "pics", paste("peakmzRT_", as.character(measurements[i,1]), sep = ""))
				png(filename = path, bg = "white")    
				plot(MSlist[["Peaklist"]][,1], MSlist[["Peaklist"]][,5], xlab = "m/z", ylab = "RT", pch = 19, cex = 0.3, main = "Picked peaks")
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
    

