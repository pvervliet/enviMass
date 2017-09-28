

library(Rcpp); 

load("G:/PART_1/MS PROJECTS/Sabine/RUES_1411_to_1502/logfile.emp")



system.time({

		num_cores<-detectCores(all.tests = FALSE, logical = TRUE)
		#clus <- makeCluster(getOption("cl.cores", num_cores), useXDR = FALSE, methods = FALSE)
		clus <- makeCluster(getOption("cl.cores", num_cores))		
		clusterEvalQ(cl = clus,{library(enviMass, verbose=FALSE); NULL}) # time consuming!
		######################

})
		
		
		######################				
		#stopCluster(clus)





		#########################################################################################
		# set up profileList_pos ################################################################
	    profileList_pos <- list(0)
	    profileList_pos[[1]] <- data.frame(TRUE, FALSE, FALSE, FALSE)    # state
	    colnames(profileList_pos[[1]]) <- c("peaks?", "agglom?", "profiling", "trends?")
	    profileList_pos[[2]] <- 0  # peaks
	    profileList_pos[[3]] <- 0  # datetime
	    profileList_pos[[4]] <- 0  # time
	    profileList_pos[[5]] <- 0  # place
	    profileList_pos[[6]] <- 0  # index_agglom
	    profileList_pos[[7]] <- 0  # index_prof
	    profileList_pos[[8]] <- 0  # parameters
	    profileList_pos[[9]] <- 0  # sample type
	    names(profileList_pos) <- c("state","peaks","datetime","sampleID","place",
			"index_agglom","index_prof","parameters","type")
		for_files<-enviMass:::get_measurement_IDs(
				logfile,
				sets = as.numeric(logfile$parameters$prof_maxfiles),
				ion_mode = "positive",
				until = logfile$parameters$upto_file,
				selective = logfile$parameters$prof_select,
				types = c("sample", "blank", "spiked"),
				places = FALSE,
				check_exist = TRUE
			)
    	profileList_pos[["datetime"]] <- for_files[["datetime"]];
    	profileList_pos[["sampleID"]] <- for_files[["sampleID"]];
		profileList_pos[["place"]] <- for_files[["locus"]];
    	profileList_pos[["type"]] <- for_files[["typus"]];
		#########################################################################################
		# load peaks into profileList_pos #######################################################
		clusterEvalQ(cl = clus,{rm(list=ls()); NULL})
		cluster_results <- clusterApplyLB(cl = clus, 
			x = for_files[["sampleID"]], 
			fun = startprofiles_pl, 
			logfile = logfile,
			frac = FALSE,
			blind_omit=as.logical(logfile$parameters$blind_omit)
		)
		clusterEvalQ(cl = clus,{rm(list=ls()); NULL})
		profileList_pos[["peaks"]] <- do.call(rbind, cluster_results)
		rm(cluster_results)
		profileList_pos[["peaks"]] <- profileList_pos[["peaks"]][order(profileList_pos[["peaks"]][,"m/z"], decreasing = FALSE),]
		#########################################################################################
		if(any(profileList_pos[[2]][,2] == 0)){stop("\n issue in do_profiling: zero intensities detected. Try to rerun the workflow including the peakpicking, using -> Settings -> General -> Reset project including peak picking.")}
		profileList_pos <- agglomer(
			profileList_pos,
			dmass = (as.numeric(logfile$parameters$prof_dmz) + 1),
			ppm = as.logical(as.character(logfile$parameters$prof_ppm)),
			dret = (as.numeric(logfile$parameters$prof_drt) + 10)
		)
		#########################################################################################

		
sourceCpp("D:/Users/uchemadmin/Desktop/MS/R_packages/enviMass/enviMass/enviMass/inst/extractProfiles.cpp")


system.time({		
		#########################################################################################		
		profileList = profileList_pos
		dmass = as.numeric(logfile$parameters$prof_dmz)
		ppm = as.logical(as.character(logfile$parameters$prof_ppm))
		dret = as.numeric(logfile$parameters$prof_drt)
		from = FALSE
		to = FALSE
		progbar = FALSE
		plot_it = FALSE
		replicates = FALSE
		IDs = FALSE
		with_test = FALSE
		do_replicates<-FALSE
		m=1
		n=length(profileList[["index_agglom"]][,1])
		if(ppm){ppm2=1}else{ppm2=2};
	
		for(k in m:n){
			if(profileList[["index_agglom"]][k,3]>1){
				delmz<-(max(profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),1])-min(profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),1]))
				if(ppm){
					delmz<-(delmz/mean(profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),1])*1E6)
				}else{
					delmz<-(delmz/1000)
				}
				delRT<-(max(profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),3])-min(profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),3]))
				if( 
					(delmz>(dmass*2)) || (delRT>dret) || (any(duplicated(profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),6]))) 
				){	
				
					#################################################################################
					
					clusters <- 
						
						extractProfiles_new(
							peaks = profileList[["peaks"]][(profileList[["index_agglom"]][k,1] : 
								profileList[["index_agglom"]][k,2]), c("m/z", "intensity", "RT", "sampleIDs")],                   
							in_order = order(profileList[["peaks"]][(profileList[["index_agglom"]][k,1] :
								profileList[["index_agglom"]][k,2]),"intensity"], decreasing = TRUE),   	# intensity order 
							pregroup = 0, 		# no replicate pre-sorting
							dmass = dmass,
							ppm = ppm,
							dret = dret,
							run_pregroup = FALSE,
							verbose = FALSE
						)			
					
			stop()				
									
									
					#################################################################################					

					
					#if(any(clusters[,10]==0)){stop("missing assignment found_1!")} # check - all must be assigned

	
	

				}
			}
		}
		#########################################################################################	
})
		
		print(clusters)
		
		
		
		
		
		
		
		
		
		
		
		
		
		