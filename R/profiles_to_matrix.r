#' @title Converts profileList to a matrix
#'
#'
#' @description Converts profileList peak list to a matrix of peak intensities. 
#'
#' @param profileList A profile list, i.e., profileList_pos or profileList_neg.
#' @param links_profiles List of links among profiles, i.e., links_profiles_pos or links_profiles_neg. Only used if reduce_comp=TRUE.
#' @param reduce_comp. Logical. Reduce profiles by their components. See details.
#' @param sort_by. Vector of character strings. One or several of colnames(profileList[["index_prof"]]) to prioritize profiles with.
#' @param sort_decreasing. Logical. Should the sort_by done in decreasing (TRUE) or increasing (FALSE) order? Applied to all strings of sort_by.
#' @param n_profiles Integer. How many of the most sort_by profiles to include? Set to NULL to include all.
#' @param only_sample_peaks. Logical. TRUE = should peaks of sample files (and not, e.g., blank files) be included only?
#' @param n_latest_peaks. NULL or Integer. If integer, number of latest file peaks to include.
#' @param median_above_blind. NULL or Integer. If integer, only profiles with a median sample-to-blank intensity >= median_above_blind are used. Set to Inf to only retain profiles entirely unaffected by blind/blank signals.
#' @param normalize. Logical. TRUE = should intensities of each profile be divided by their maximum to range in [0,1]?
#'
#' @return Matrix with intensities, cp. value section.
#' 
#' @details The function writes the (sample) peak intensities of the sorted n_profiles into a non-sparse matrix, cp. value section.
#'
#'
#'
#'
#' @value Matrix with peak intensities, with the following properties: 
#' Columns refer to profiles, rows to files.
#' Rows are ordered by date and time and named by file IDs.  
#' Columns are ordered by decreasing maximum profile intensity and named with profile IDs. 
#' Missing peak intensities (<LOD) are set to 0.
#'
#'
#'
#'
#'

profiles_to_matrix <- function(
	profileList,
	links_profiles = NULL,
	sort_by = c("number_peaks_sample", "mean_int"),
	sort_decreasing = TRUE, 
	reduce_comp = FALSE,
	n_profiles = NULL,
    only_sample_peaks = FALSE,
    n_latest_peaks = NULL,
	mean_above_blind = NULL,
    normalize = FALSE
){

	############################################################################
	if(any(is.na(match(sort_by, colnames(profileList[["index_prof"]]))))){stop("Argument sort_by not matching column names - abort.")}
    if(!profileList[["state"]]["profiling"][[1]]){stop("\nprofileList not profiled - abort.")}
    len <- dim(profileList[["index_prof"]])[1]
    if(!is.logical(sort_decreasing)){stop("\nArgument sort_decreasing must be logical - abort.")}
    if(!is.logical(reduce_comp)){stop("\nArgument reduce_comp must be logical - abort.")}	
    if(!is.logical(only_sample_peaks)){stop("\nArgument only_sample_peaks must be logical - abort.")}
    if(!is.logical(normalize)){stop("\nArgument normalize must be logical - abort.")}
	if((reduce_comp) & (!is.list(links_profiles))){stop("\n Either set reduce_comp=FALSE or provide a valid links_profiles list - abort.")}
	############################################################################

    ############################################################################    
    # get sample IDs ###########################################################
    ord <- order(as.POSIXct(profileList[["datetime"]]), decreasing = TRUE)
    type <- profileList[["type"]][ord]
    file_ID <- profileList[["sampleID"]][ord]
	if(only_sample_peaks){	
		file_ID <- file_ID[type == "sample"]
	}
    if(!is.null(n_latest_peaks)[1]){
        type <- type[1:n_latest_peaks]
		file_ID <- file_ID[1:n_latest_peaks]
    }
	keep <- rep(TRUE,len)
	profile_IDs <- profileList[["index_prof"]][,"profile_ID"]
    # filter out profiles which do not contain sample peaks ####################
	if(only_sample_peaks){	
		keep[profileList[["index_prof"]][,"number_peaks_sample"]==0] <- FALSE
	}
	# filter out profiles which range not above blind intensities ##############
	if(!is.null(mean_above_blind)[1]){
		keep[profileList[["index_prof"]][,"above_blind?"] < mean_above_blind] <- FALSE
	}
    ############################################################################ 

    ############################################################################    
	# based on sorting, remove correlated components files with a lower rank ###
	if(reduce_comp){
		# BEWARE: must include use_profiles argument to account for above filtering
		keep_IDs <- enviMass::analyseE_links_profiles(
				profileList_index = profileList[["index_prof"]], 
				links_profiles, 
				sort_what = sort_by, 	# internally sorted in function - again sorted below
				sort_decreasing = sort_decreasing, 
				use_profile = keep,	# omit/skip filtered profiles ab initio
				return_excl = TRUE) # returns IDs of excluded profiles 	
		if(length(keep_IDs)){
			keep[match(keep_IDs, profile_IDs)] <- FALSE
		}
	}
    # sort #####################################################################	
	max_ord <- profileList[["index_prof"]][,sort_by,drop=FALSE]
	profile_IDs <- profile_IDs[keep]
	max_ord <- max_ord[keep,,drop=FALSE]
	max_int_ord <- rev(do.call(order,as.data.frame(max_ord))) # ordering for multiple columns
	profile_IDs <- profile_IDs[max_int_ord]
    ############################################################################    

    ############################################################################    
    # write to matrix ##########################################################	
	if(is.null(n_profiles)){
		n_profiles <- length(profile_IDs)
	}else{
		if(n_profiles > length(profile_IDs)){
			n_profiles <- length(profile_IDs)
		}
	}
    mat<-matrix(nrow=length(file_ID),ncol=n_profiles,0)
    rownames(mat) <- as.character(file_ID)
	colnames(mat) <- as.character(profile_IDs[1:n_profiles])
    sub_peaks_ind <- (
        !is.na(
            match(
                profileList[["peaks"]][,"sampleIDs"],
                as.numeric(rownames(mat))
            )
        ) &
        !is.na(
            match(
                profileList[["peaks"]][,"profileIDs"],
                as.numeric(colnames(mat))
            )
        )
    )
	mat[
        cbind(
            match(as.character(profileList[["peaks"]][sub_peaks_ind,"sampleIDs"]),rownames(mat)),
            match(as.character(profileList[["peaks"]][sub_peaks_ind,"profileIDs"]),colnames(mat))
        )
    ] <- profileList[["peaks"]][sub_peaks_ind,"intensity"]
    ############################################################################
	
    ############################################################################
	# checks ###################################################################
	if(any(is.na(mat))){ 
		stop("Sth went wrong in function profiles_to_matrix - debug_1!")
	}	
	if(any(apply(mat,2,max)==0)){ 
		stop("Sth went wrong in function profiles_to_matrix - debug_2!")
	}
    ############################################################################ 	

    ############################################################################
    # (0,1)-normalize ##########################################################
    if(normalize) mat <- sweep(mat,2,apply(mat,2,max),"/")
    ############################################################################

    ############################################################################
    return(mat)
	
}
