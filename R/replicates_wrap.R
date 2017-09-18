

replicates_wrap <-
function(
	x, 
	logfile,
	measurements,
	ppm,
	mz_tol,
	rt_tol,
	int_tol
){

	replic <- x
	report <- "nothing"
	IDs <- measurements$ID[measurements$tag3 == x]
	with_test <- TRUE
    ############################################################################			
	if(any(duplicated(IDs))){stop("replicates: non-unique IDs found!")} # should not happen anyway
	# initialize intersection rectangles with first peaklist ###################
	if(any(objects() == "peaklist")){rm(peaklist)}
	load(file = file.path(logfile$project_folder, "peaklist", as.character(IDs[1])), verbose = FALSE);
	peaklist <- peaklist[,c(12,13,14)]
	if(ppm){
		low_mass <- (peaklist[,1] - (mz_tol/1E6*peaklist[,1]))
		high_mass <- (peaklist[,1] + (mz_tol/1E6*peaklist[,1]))
	}else{
		low_mass <- (peaklist[,1] - (mz_tol/1000))
		high_mass <- (peaklist[,1] + (mz_tol/1000))
	}
	low_rt <- (peaklist[,3] - rt_tol)
	high_rt <- (peaklist[,3] + rt_tol)							
	low_int <- (peaklist[,2] - int_tol)
	low_int[low_int<0] <- 0
	high_int <- (peaklist[,2] + int_tol)				
	ID <- seq(1, length(peaklist[,1]),1)
	intersection <- cbind(low_mass, high_mass, low_rt, high_rt, low_int, high_int,ID) 
	if(with_test){ # test
		intersection <- rbind(
			intersection,
			c(20,20,5,15,1E4,1E5,-2),
			c(.3,.3001,39,41,1E4,1E5,-3),
			c(-800.001,-800.002,1400,1430,1E4,1E5,-4))
	}
	rm(peaklist)
	for(j in 2:length(IDs)){
		# load next peaklist
		load(file = file.path(logfile[[1]], "peaklist", as.character(IDs[j])), verbose=FALSE);
		peaklist <- peaklist[,c(12,13,14)]
		if(ppm){
			low_mass <- (peaklist[,1] - (mz_tol/1E6*peaklist[,1]))
			high_mass <- (peaklist[,1] + (mz_tol/1E6*peaklist[,1]))
		}else{
			low_mass <- (peaklist[,1] - (mz_tol/1000))
			high_mass <- (peaklist[,1] + (mz_tol/1000))	
		}
		low_rt <- (peaklist[,3] - rt_tol)
		high_rt <- (peaklist[,3] + rt_tol)		
		low_int <- (peaklist[,2] - int_tol)
		low_int[low_int<0] <- 0
		high_int <- (peaklist[,2] + int_tol)										
		query <- cbind(low_mass, high_mass, low_rt, high_rt, low_int, high_int)			
		if(with_test){
			query <- rbind(
				query,	
				c(20,20,5,15,1E4,1E5),
				c(.3,.3001,39,41,1E4,1E5),							
				c(-1000.001,-1000.002,1400,1430,1E4,1E5))								
		}
		rm(peaklist)
		# build a boxtree with next peaklist
		tree <- .Call("boxtree", 
			as.matrix(query),
			PACKAGE="nontarget"
		)	
		colnames(tree) <- c("LOSON", "MIDSON","HISON", "level", "disc")		
		# query
		search_bounds <- rep(0,6)					
		new_intersection <- matrix(nrow=10000,ncol=(6+j),-1)					
		len <- 10000
		at_new <- 1
		for(k in 1:length(intersection[,1])){
			search_bounds[1] <- intersection[k,1]
			search_bounds[2] <- intersection[k,2]				
			search_bounds[3] <- intersection[k,3]
			search_bounds[4] <- intersection[k,4]				
			search_bounds[5] <- intersection[k,5]
			search_bounds[6] <- intersection[k,6]	
			found <- .Call("search_boxtree", 
				query,
				tree,
				as.numeric(search_bounds),
				as.integer(1), # return full findings
				PACKAGE="nontarget"
			)
			if(length(found) > 0){
				for(m in 1:length(found)){
					new_intersection[at_new,1] <- max(c(intersection[k,1],query[found[m],1]))
					new_intersection[at_new,2] <- min(c(intersection[k,2],query[found[m],2]))							
					new_intersection[at_new,3] <- max(c(intersection[k,3],query[found[m],3]))							
					new_intersection[at_new,4] <- min(c(intersection[k,4],query[found[m],4]))
					new_intersection[at_new,5] <- max(c(intersection[k,5],query[found[m],5]))							
					new_intersection[at_new,6] <- min(c(intersection[k,6],query[found[m],6]))								
					new_intersection[at_new,(7:(7+j-2))] <- intersection[k,(7:(7+j-2))]
					new_intersection[at_new,(6+j)] <- found[m]
					at_new <- (at_new + 1)
					if(at_new > len){ # extend new_intersection
						new_intersection <- rbind(
							new_intersection,
							matrix(nrow=10000, ncol=(6+j), -1)
						)
						len <- (len + 10000)
					}
				}	
			}				
		}
		intersection <- new_intersection[new_intersection[,7] != -1,,drop=FALSE] # omit empty entries
		if(at_new == 1){break} # no further intersections detected - abort
	}
	if(with_test){					
		if(!any(intersection[,7] == -2)){stop("intersection test_1 failed")}
		if(!any(intersection[,7] == -3)){stop("intersection test_2 failed")}
		if(any(intersection[,7] == -4)){stop("intersection test_3 failed. Senseless mass tolerances in use?")}					
		intersection <- intersection[intersection[,7] > 0,]	# remove test cases	
	}
	# clean & fill peaklists ###################################################
	for(j in 1:length(IDs)){
		load(file = file.path(logfile[[1]], "peaklist", as.character(IDs[j])), verbose=FALSE);
		keep <- rep(0,length(peaklist[,1]))
		if(length(intersection[,1]) > 0){ # any intersections at all?
			keep_those <- unique(intersection[,(6+j)])						
			keep[keep_those] <- 1
		}
		peaklist[,colnames(peaklist) == "keep"] <- keep
		save(peaklist, file = file.path(logfile[[1]], "peaklist", as.character(IDs[j])))
		rm(peaklist)
	}
 	############################################################################
	report <- "done"
	return(report)
	
}




















