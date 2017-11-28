#' @title Pattern decomposition function 
#'
#' @description Pattern decomposition function 
#'
#' @param cent_peak_mat
#' @param pattern
#' @param profileList
#' @param LOD
#' @param RT_tol_inside
#' @param int_tol
#' @param use_score_cut Logical.
#' @param score_cut Numeric.
#' @param plot_it Logical.
#' @param verbose Logical.
#' @param RT_seperate Logical.
#' 
#' @details enviMass workflow function
#' 

recomb_score_pl<-function(
		cent_peak_mat,
		pattern_compound,
		peaks,
		LOD,
		RT_tol_inside,
		int_tol,
		use_score_cut = FALSE,
		score_cut = 0,		
		plot_it = FALSE,
		verbose = FALSE,
		RT_seperate = FALSE
	){
			

	#######################################################################		
	results <- list()
	at_results <- 1
	checked <- TRUE
	check_nodes <- list()
	check_nodes_index <- list()
	# Pre-decompose by RT gaps, especially for suspect screening 
	if((length(cent_peak_mat[,2])>1) & RT_seperate){
		cent_peak_mat <- cent_peak_mat[
			order(peaks[cent_peak_mat[,2], "RT"], decreasing = FALSE)
		,, drop = FALSE]	
		at_RT <- peaks[cent_peak_mat[,2], "RT"]
		in_node <- rep(1,length(cent_peak_mat[,1]))
		init_node <- 1
		for(i in 2:length(cent_peak_mat[,2])){
			if((at_RT[i] - at_RT[i-1]) > RT_tol_inside){
				init_node <- (init_node+1)
				in_node[i] <- init_node
			}else{
				in_node[i] <- init_node
			}
		}
		if(init_node>1){
			for(i in 1:max(init_node)){
				check_nodes[[i]] <- cent_peak_mat[in_node == i,, drop = FALSE]
				check_nodes[[i]] <- (check_nodes[[i]][order(check_nodes[[i]][,1],decreasing=FALSE),,drop=FALSE])
				check_nodes_index[[i]] <- sum(in_node == i)
			}
		}else{
			check_nodes[[1]] <- cent_peak_mat # initialize with full set
			check_nodes_index[[1]] <- length(cent_peak_mat[,1])
		}
	}else{
		check_nodes[[1]] <- cent_peak_mat # initialize with full set
		check_nodes_index[[1]] <- length(cent_peak_mat[,1])
	}
	#######################################################################	
	results <- while_checked( # Rcpp
		check_nodes,
		pattern_compound,
		peaks,
		RT_tol_inside = RT_tol_inside,
		int_tol = int_tol,
		verbose = verbose
	)
	#######################################################################	
	results2 <- list()
	at_results <- 1
	for(k in 1:length(results)){

		# order by increasing peak number #################################
		results[[k]] <- results[[k]][
			order(results[[k]][,1], decreasing = FALSE)
		,,drop = FALSE]
		# calculate score1 ################################################
		rescale <- weighted.mean(
			x = (pattern_compound[results[[k]][,1],2] / peaks[results[[k]][,2],2]),
			w = (peaks[results[[k]][,2],2] / (int_tol / 100 * peaks[results[[k]][,2],2]) )
		)
		above_LOD <- ((pattern_compound[,2] / rescale) > LOD)
		# ... measured "above LOD threshold:"
		if(any(above_LOD)){
			score1 <- (
				sum(pattern_compound[results[[k]][,1],2][above_LOD[results[[k]][,1]]]) / sum(pattern_compound[above_LOD,2])
			)
			score1 <- round(score1, digits = 4)
		}else{
			score1 <- NA
		}
		if(use_score_cut == "TRUE"){
			if(is.na(score1)){ 		# ... discarded if below LOD
				if(verbose){cat("- all below LOD.")}
				next;
			}
			if(score1 <= score_cut){ 	# 
				if(verbose){cat("- discarded by score.")}				
				next;
			}				
		}
		# calculate score2 ################################################
		if(any(!above_LOD)){
			score2 <- sum(pattern_compound[results[[k]][,1],2][!above_LOD[results[[k]][,1]]]) / sum(pattern_compound[,2])
		}else{
			score2 <- NA
		}
		###################################################################
		results2[[at_results]] <- list() 
		results2[[at_results]][[1]] <- results[[k]]
		results2[[at_results]][[2]] <- score1
		results2[[at_results]][[3]] <- score2
		results2[[at_results]][[4]] <- ((pattern_compound[results[[k]][,1],1] - peaks[results[[k]][,2],1]) / mean(pattern_compound[results[[k]][,1],1]) * 1E6)				
		results2[[at_results]][[5]] <- (mean(peaks[results[[k]][,2],3]) - peaks[results[[k]][,2],3])
		results2[[at_results]][[6]] <- rescale
		results2[[at_results]][[7]] <- peaks[results[[k]][,2],1]
		results2[[at_results]][[8]] <- peaks[results[[k]][,2],2]
		results2[[at_results]][[9]] <- peaks[results[[k]][,2],3]
		names(results2[[at_results]]) <- c("Peaks", "score_1", "score_2", "ppm deviation", "RT deviation from mean", "rescale factor", "m/z", "Intensity", "RT")
		###################################################################
		at_results <- (at_results + 1)

	}
	#######################################################################	
	# resort list - largest combinations go first #########################
	results3 <- list()
	len <- rep(0,length(results2))
	len2 <- rep(0,length(results)) 
	for(k in 1:length(results2)){
		len[k] <- length(results2[[k]][[1]][,1])
		len2[k] <- sum(results2[[k]]$Intensity)
	}
	ord <- order(len, len2, decreasing = TRUE)
	for(k in 1:length(results2)){	
		results3[[k]]<-results2[[ord[k]]]
	}
	#######################################################################
	return(results3)
	#######################################################################
	
}
