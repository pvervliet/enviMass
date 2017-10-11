#' @title Removes isotopologue links from links_profiles
#'
#' @export
#'
#' @description Removes isotopologue links from links_profiles
#'
#' @param links_profiles links_profiles_pos or links_profiles_neg
#' @param profileList profileList_pos or profileList_neg
#' @param cut_delRT_isot Maximum RT deviation among same-sample isotopologue peaks in a profile
#' @param cut_cor_isot Minimum cross-sample intensity correlation between isotopologue peaks of any two profiles
#' @param cut_frac_iso Minimum co-occurence fraction of peaks within cut_delRT_isot
#' 
#' @details enviMass workflow function. 
#' 

cleanA_links_profiles<-function(
	links_profiles, 
	profileList, 
	cut_delRT_isot = 5, 
	cut_cor_isot = .9, 
	cut_frac_iso = .9
){

	################################################################
	for_profs <- (1:length(links_profiles))
	found <- 0
	removed <- 0
	for(n in 1:length(for_profs)){
		prof1 <- as.numeric(names(links_profiles)[for_profs[n]])
		# extract delRT_isot & int_cor
		if(dim(links_profiles[[for_profs[n]]]$isot)[1] > 0){
			keep <- rep(TRUE,dim(links_profiles[[for_profs[n]]]$isot)[1])
			found <- (found+length(keep))
			for(m in 1:dim(links_profiles[[for_profs[n]]]$isot)[1]){
			
				if(links_profiles[[for_profs[n]]]$isot[m,"use"] == 0) next; # will be dealt with by other partner profile
				prof2 <- links_profiles[[for_profs[n]]]$isot[m,"linked profile"]
				num_peaks <- links_profiles[[for_profs[n]]]$isot[m,"ref_1"]
				these <- profileList[["peaks"]][
					profileList[["index_prof"]][prof1,"start_ID"]:profileList[["index_prof"]][prof1,"end_ID"]
				,"sampleIDs"]
				those <- profileList[["peaks"]][
					profileList[["index_prof"]][prof2,"start_ID"]:profileList[["index_prof"]][prof2,"end_ID"]
				,"sampleIDs"]
				matched <- match(these,those)
				# first, check delRT ##############################
				RT_1 <- (profileList[["peaks"]][
						(profileList[["index_prof"]][prof1,"start_ID"]:profileList[["index_prof"]][prof1,"end_ID"])
						[!is.na(matched)]
					,"RT"])
				RT_2 <- (profileList[["peaks"]][
						(profileList[["index_prof"]][prof2,"start_ID"]:profileList[["index_prof"]][prof2,"end_ID"])
						[matched[!is.na(matched)]]
					,"RT"])
				dRT <- abs(RT_1-RT_2)
				got_frac <- (sum(dRT < cut_delRT_isot) / length(dRT))
				if(got_frac < cut_frac_iso){
					keep[m] <- FALSE
					next;
				}
				# second, check and insert correlation ############
				if(num_peaks > 2){
					int_1 <- (profileList[["peaks"]][
							(profileList[["index_prof"]][prof1,"start_ID"]:profileList[["index_prof"]][prof1,"end_ID"])
							[!is.na(matched)]
						,"intensity"])
					int_2 <- (profileList[["peaks"]][
							(profileList[["index_prof"]][prof2,"start_ID"]:profileList[["index_prof"]][prof2,"end_ID"])
							[matched[!is.na(matched)]]
						,"intensity"])
					int_cor <- cor(int_1, int_2)
					if(!is.na(int_cor)){
						if(int_cor < cut_cor_isot){
							keep[m] <- FALSE
						}
					}else{ # insert correlation
						int_cor <- as.integer(round(int_cor * 1000))
						links_profiles[[for_profs[n]]]$isot[m,"correl"] <- int_cor
						at_entry_2 <- profileList[["index_prof"]][prof2,"links"]
						at <- which(links_profiles[[at_entry_2]]$isot[,"linked profile"] == prof1)
						links_profiles[[at_entry_2]]$isot[at,"correl"] <- int_cor
					}
				}
			}
			
			# remove invalid isotopologue links ####################
			if(any(!keep)){
			
				removed <- (removed + sum(keep == FALSE))
				# mark the other partner profiles by "2"
				those<-links_profiles[[for_profs[n]]]$isot[!keep,"linked profile"]
				for(prof2 in those){
					# get link of prof2
					at_entry_2 <- profileList[["index_prof"]][prof2,"links"]
					at <- which(links_profiles[[at_entry_2]]$isot[,"linked profile"]==prof1)
					links_profiles[[at_entry_2]]$isot <- links_profiles[[at_entry_2]]$isot[-at,,drop=FALSE]
				}
				# delete invalid entries for prof1 
				links_profiles[[for_profs[n]]]$isot<-
					links_profiles[[for_profs[n]]]$isot[keep,,drop=FALSE]
					
			}
		}
	}
	perce <- as.character(round(removed/found,digits=3)*100)
	cat("\n");
	cat(paste0(perce, "% of ", found," isotopologue links filtered."));cat("\n");
	################################################################
	return(links_profiles)
	
}
