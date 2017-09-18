#' @title Extract profiles from peak set partitions, parallelized
#'
#' @export
#'
#' @description \code{partcluster} extract profiles from peak set partitions generated by \code{\link{agglomer}}.
#'
#' @param profileList A profile list.
#' @param dmass Numeric. m/z gap size
#' @param ppm Logical. \code{dmass} given in ppm?
#' @param dret Numeric. RT gap size; units equal to those of the input files
#' @param replicates. FALSE or vector of character strings with replicate labels for each LC-MS file, i.e., files of the same replicate have the same character string.
#' @param IDs. Integer vector. IDs of files, required if replicates is not set to FALSE.
#' @param clus. Socket cluster, to be set up by package parallel
#'
#' @return Updated profile list
#' 
#' @details enviMass workflow function. Works along decreasing intensities. The remaining peak of highest intensity not
#' yet part of a profile is either assigned to an existing profile (closest in mass) or initializes a new profile. With
#' addition of a peak to a new profile, profile mass tolerances are gradually adapted. If replicates are profiled, profiles
#' are first extracted in each replicate level and these profiles than further merged.
#' 
#' @seealso \code{\link{startprofiles}}, \code{\link{agglomer}}


partcluster_pl<-function(
	profileList,
	dmass = 3,
	ppm = TRUE,
	dret = 60,
	replicates = FALSE,
	IDs = FALSE,
	clus
){

	########################################################################################
	if(!profileList[[1]][[2]]){stop("run agglom first on that profileList; aborted.")}
	if(!is.numeric(dmass)){stop("dmass must be numeric; aborted.")}
	if(!is.numeric(dret)){stop("dret must be numeric; aborted.")}
	if(!is.logical(ppm)){stop("ppm must be logical; aborted.")}
	startat <- c(0);
	if(ppm){ppm2 = 1}else{ppm2 = 2};
	do_replicates <- FALSE; 
	replic <- FALSE
	if(any(replicates != "FALSE")){
		replic <- replicates[duplicated(replicates)]
		replic <- unique(replic)			  
		replic <- replic[replic != "FALSE"]
		if(length(replic) > 0){
			if(length(replicates) != length(IDs)){
				stop("\n replicates vector longer than file ID vector")
			}
			do_replicates <- TRUE
		}
	}
	########################################################################################
	# assign blocks & clean agglom index in each block
	along<-seq(1,dim(profileList[["index_agglom"]])[1],1)
	size_MB <- (as.numeric(object.size(profileList[["peaks"]])) / 1048576)
	for_split <- (size_MB / 2)
	for_split <- round( dim(profileList[["index_agglom"]])[1] / for_split )
	#for_split<-round( dim(profileList[["index_prof"]])[1] / (dim(summary(clus))[1]) )
	along <- split(along, ceiling(seq_along(along) / for_split) )
	clus_peaks <- list()
	clus_agglom <- list()
	for(i in 1:length(along)){
		clus_peaks[[i]] <- profileList[["peaks"]][
			(profileList[["index_agglom"]][along[[i]][1],"start_ID"]) : (profileList[["index_agglom"]][along[[i]][length(along[[i]])],"end_ID"])
		,,drop=FALSE]
		clus_agglom[[i]] <- profileList[["index_agglom"]][along[[i]][1] : along[[i]][length(along[[i]])],,drop=FALSE]
		start_at <- (clus_agglom[[i]][1,"start_ID"][[1]]-1)
		clus_agglom[[i]][,1] <- (clus_agglom[[i]][,1]-start_at)
		clus_agglom[[i]][,2] <- (clus_agglom[[i]][,2]-start_at)
	}
	clusterEvalQ(cl = clus,{rm(list = ls()); gc(verbose=FALSE); NULL})
	clusterExport(cl = clus, 
		varlist = c("dmass", "ppm2", "dret", "do_replicates", "replic", "IDs"), 
		envir = environment())	
	cluster_results <- clusterMap(
		cl = clus,
		fun = enviMass:::clust_func_partclust,
		peaks = clus_peaks,
		agglom = clus_agglom,
		RECYCLE = TRUE,
		SIMPLIFY = FALSE, # no!
		USE.NAMES = FALSE,
		.scheduling = c("dynamic")
	)
	clusterEvalQ(cl = clus,{rm(list = ls()); gc(verbose = FALSE); NULL})
	# assort results into original profileList
	if(length(cluster_results) > 1){
		for(i in 2:length(cluster_results)){ # insert profileIDs
			cluster_results[[i]][,"profileIDs"] <- (
				cluster_results[[i]][,"profileIDs"] + max(cluster_results[[i-1]][,"profileIDs"])
			)
		}
	}
	profileList[["peaks"]] <- do.call(rbind,cluster_results)
	rm(cluster_results, clus_agglom, clus_peaks)
	########################################################################################
	# assemble profile index matrix ########################################################
	index <- .Call("_enviMass_indexed",
		as.integer(profileList[["peaks"]][,"profileIDs"]),
		as.integer(max(profileList[["peaks"]][,"profileIDs"])),	
		as.integer(26),
		PACKAGE="enviMass"
	)
	index <- index[index[,1] != 0,]
	index[,4] <- seq(length(index[,4]))
	colnames(index) <- c(
		"start_ID",
		"end_ID",
		"number_peaks_total", #1
		"profile_ID",
		"deltaint_newest", #"current_incident"
		"deltaint_global", #4 #"past_incident"
		"absolute_mean_dev",
		"in_blind?",
		"above_blind?", #7
		"number_peaks_sample",
		"number_peaks_blind", #10
		"mean_int_sample",
		"mean_int_blind", #12
		"mean_mz",
		"mean_RT",
		"mean_int", #14
		"newest_intensity", #"newest_intensity"
		"links",
		"component", #17
		"max_int_sample",
		"max_int_blind",
		# new:
		"max_int",
		"var_mz",
		"min_RT",
		"max_RT",
		"Mass defect"
	)
	profileList[[7]] <- index
	########################################################################################
	# get characteristics of individual profiles ###########################################
	# -> could still be parallelized	
	m = 1
	n = dim(profileList[["index_prof"]])[1]
	mean_mz <- rep(0, n)
	mean_RT <- rep(0, n)
	mean_int <- rep(0, n)
	mean_int_blind <- rep(0, n)
	mean_int_sample <- rep(0, n)	
	max_int <- rep(0, n)
	max_int_blind <- rep(0, n)
	max_int_sample <- rep(0, n)		
	count_blind <- rep(0, n)
	count_sam <- rep(0, n)
	min_RT <- rep(0, n)
	max_RT <- rep(0, n)
	var_mz <- rep(0, n)
	Mass_defect <- rep(0, n)	
	sample_IDs <- as.numeric(profileList[["sampleID"]])
	sample_types <- profileList[["type"]]
    for(k in m:n){
		sub_mat <- profileList[["peaks"]][profileList[["index_prof"]][k,1]:  profileList[["index_prof"]][k,2], c(1:3, 6), drop=FALSE]
		len <- dim(sub_mat)[1]
		mean_mass <- (sum(sub_mat[,1])/len)
		mean_mz[k] <- mean_mass
		mean_RT[k] <- (sum(sub_mat[,3])/len)	  
		mean_int[k] <- (sum(sub_mat[,2])/len)	
		max_int[k] <- max(sub_mat[,2])
		min_RT[k] <- min(sub_mat[,3])	
		max_RT[k] <- max(sub_mat[,3])	
		var_mz <- (sum((sub_mat[,1]-mean_mass)^2)/(len-1))
		Mass_defect[k] <- (round(mean_mass)-mean_mass)
		are_blind <- sub_mat[,4] %in% (sample_IDs[sample_types == "blank"])
		count_blind[k] <- sum(are_blind)
		count_sam[k] <- (len - count_blind[k])
		if(count_blind[k] > 0){ 
			mean_int_blind[k] <- (sum(sub_mat[are_blind,"intensity"])/(len-1))
			max_int_blind[k] <- max(sub_mat[are_blind,"intensity"])
		}
		if(count_sam[k] > 0){
			mean_int_sample[k] <- (sum(sub_mat[!are_blind,"intensity"])/(len-1))
			max_int_sample[k] <- max(sub_mat[!are_blind,"intensity"])
		}	
	}
	profileList[["index_prof"]][,"mean_mz"] <- mean_mz
	profileList[["index_prof"]][,"mean_RT"] <- mean_RT
	profileList[["index_prof"]][,"mean_int"] <- mean_int
	profileList[["index_prof"]][,"mean_int_blind"] <- mean_int_blind 
	profileList[["index_prof"]][,"mean_int_sample"] <- mean_int_sample
	profileList[["index_prof"]][,"max_int"] <- max_int
	profileList[["index_prof"]][,"max_int_sample"] <- max_int_sample
	profileList[["index_prof"]][,"max_int_blind"] <- max_int_blind	
	profileList[["index_prof"]][,"in_blind?"] <- (count_blind > 0)
	profileList[["index_prof"]][,"number_peaks_blind"] <- count_blind 
	profileList[["index_prof"]][,"number_peaks_sample"] <- count_sam 
	profileList[["index_prof"]][,"min_RT"] <- min_RT
	profileList[["index_prof"]][,"max_RT"] <- max_RT
	profileList[["index_prof"]][,"Mass defect"] <- Mass_defect
	profileList[["index_prof"]][,"var_mz"] <- var_mz	
	profileList[["state"]][[3]] <- TRUE
	########################################################################################
	return(profileList)

}


