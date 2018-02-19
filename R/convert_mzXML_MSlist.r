convert_mzXML_MSlist <-
	function( 
		mzXML_file,
		scanTypes,
		read_scanType,
		minRT = FALSE,
		maxRT = FALSE,
		minmz = FALSE,
		maxmz = FALSE,
		get_acquisitionNum = FALSE
	){
	
		##########################################################################
		# checks #################################################################
		if(mode(mzXML_file) != "S4") stop("\n mzXML_file not a S4 class from mzR - please revise.")
		if(is.na(match(read_scanType, scanTypes))) stop("\n read_scanType does not any of the scanTypes - please revise.")
		if(
			(minmz != FALSE) |  
			(maxmz != FALSE) |  			
			(minRT != FALSE) |  
			(maxRT != FALSE)
		){
			use_ranges <- TRUE
			if(minmz == FALSE) minmz <- 0  
			if(maxmz == FALSE) 	maxmz <- Inf
			if(minRT == FALSE)  minRT <- 0
			if(maxRT == FALSE)	maxRT <- Inf
		}else{
			use_ranges <- FALSE		
		}
		heads <- mzR::header(mzXML_file)
		##########################################################################
	  
		##########################################################################
		# local function definition ##############################################
		genMSlist <-  function(){
			# define MSlist list #################################################
			MSlist <- list(0);
			MSlist[[2]] <- 0; MSlist[[3]] <- 0; MSlist[[4]] <- 0; MSlist[[5]] <- 0; MSlist[[6]] <-0; MSlist[[7]] <- 0; MSlist[[8]] <- 0;
			names(MSlist) <- c("State", "Parameters", "Results", "Scans", "Partition_index", "EIC_index", "Peak_index", "Peaklist")
			stage <- data.frame(FALSE, FALSE, FALSE, FALSE, FALSE);
			names(stage) <- c("Raw?", "Partitioned?", "Clustered?", "Filtered?", "Picked?");
			MSlist[[1]] <- stage;
			p <- c(   
			"agglom_dmzgap",#1
			"agglom_ppm",
			"agglom_drtgap",
			"agglom_minpeak",
			"agglom_maxint",#5
			"part_dmzgap",#6
			"part_drtgap",
			"part_ppm",
			"part_minpeak",
			"part_peaklimit",
			"part_cutfrac",
			"part_drtsmall",
			"part_stoppoints",#13
			"clust_dmzdens",#14
			"clust_ppm",		  
			"clust_drtdens",
			"clust_minpeak",
			"clust_maxint",			
			"clust_merged",
			"clust_from",
			"clust_to",#21
			"pick_minpeak",#22
			"pick_drtsmall",
			"pick_drtfill",
			"pick_drttotal",
			"pick_recurs",
			"pick_weight",
			"pick_SB",
			"pick_SN",
			"pick_minint",
			"pick_maxint",
			"pick_ended",			
			"pick_from",
			"pick_to"#34	  
			);   
			parameters <- data.frame(p, rep("0", length(p)), stringsAsFactors = FALSE);
			names(parameters) <- c("parameters", "value")
			MSlist[[2]]	<- parameters;
			######################################################################
			return(MSlist)
		}
		##########################################################################
      	  
		##########################################################################
		peaknumb <- 0
		RT <- c()
		done_cen <- FALSE
		acquisitionNum <- c()
		use_scans <- rep(FALSE, length(scanTypes))
		for(i in 1:length(scanTypes)){
			if(is.na(scanTypes[i])) next;
			if(scanTypes[i] != read_scanType) next;
			if(heads$peaksCount[i] == 0) next;
			if(!use_ranges){
				peaknumb <- (peaknumb + heads$peaksCount[i])		
			}else{
				if(heads$retentionTime[i] < minRT) next;
				if(heads$retentionTime[i] > maxRT) next;			
				masses <- mzR::peaks(mzXML_file, scans = i)[,1]
				num_peaks <- sum((masses >= minmz) & (masses <= maxmz))
				if(num_peaks == 0) next;
				peaknumb <- (peaknumb + num_peaks)
			}
			RT <- c(RT, heads$retentionTime[i])
			acquisitionNum <- c(acquisitionNum, heads$acquisitionNum[i])
			use_scans[i] <- TRUE;
		}	
		if(get_acquisitionNum) return(acquisitionNum);
		if(peaknumb == 0) stop("\nWith this file & Method setup: no peaks available from file reading.\n")
			
		##########################################################################		
		scans <- list(0)
		scans[[1]] <- RT      
        getpeaks <- matrix(nrow = peaknumb, ncol = 7, 0)		
		from <- 1;
		for(i in 1:length(use_scans)){		
			if(!use_scans[i]) next;
			peaks <- mzR::peaks(mzXML_file, scans = i)
			if(!use_ranges){
				to <- (from + heads$peaksCount[i] - 1)
				getpeaks[from:to,1] <- peaks[,1]
				getpeaks[from:to,2] <- peaks[,2]				
				getpeaks[from:to,3] <- (heads$retentionTime[i])
				from<-(to+1)
			}else{
				peaks <- peaks[
					(peaks[,1] >= minmz) &
					(peaks[,1] <= maxmz) 			
				,, drop = FALSE]
				to <- (from + dim(peaks)[1] - 1)
				getpeaks[from:to,1] <- peaks[,1]
				getpeaks[from:to,2] <- peaks[,2]				
				getpeaks[from:to,3] <- (heads$retentionTime[i])	
				from<-(to+1)
			}
		}	
		##########################################################################
		getpeaks[,4] <- seq(1, dim(getpeaks)[1], 1)
		colnames(getpeaks) <- c("m/z", "intensity", "RT", "measureID", "partID", "clustID", "peakID")
		scans[[2]] <- getpeaks;
		rm(getpeaks);
		MSlist <- genMSlist();      
		MSlist[[4]] <- scans;
		rm(scans)
		MSlist[[1]][[1]] <- TRUE
		##########################################################################
		return(MSlist);
	      
}
