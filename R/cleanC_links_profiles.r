#' @title Removes/updates EIC links on links_profiles
#'
#' @export
#'
#' @description Removes/updates EIC links on links_profiles
#'
#' @param links_profiles links_profiles_pos or links_profiles_neg
#' @param profileList profileList_pos or profileList_neg
#' @param cut_EIC Minimum EIC correlation 
#' @param cut_frac_EIC Minimum co-occurence fraction of peaks within cut_delRT_adduc
#' @param cut_delRT_EIC Maximum RT deviation [s] among same-sample EIC peaks in a profile
#' 
#' @details enviMass workflow function. 
#' 

cleanC_links_profiles <- function(
	links_profiles, 
	profileList, 
	cut_EIC = .9, 
	cut_frac_EIC = .9, 
	cut_delRT_EIC = 5
){


	################################################################
	for_profs <- (1:length(links_profiles))
	found <- 0
	removed <- 0
	for(n in 1:length(for_profs)){
		if(dim(links_profiles[[for_profs[n]]]$EIC)[1]>0){		
			prof1 <- as.numeric(names(links_profiles)[for_profs[n]])
			keep <- rep(TRUE,dim(links_profiles[[for_profs[n]]]$EIC)[1])
			for(m in 1:dim(links_profiles[[for_profs[n]]]$EIC)[1]){
				if(links_profiles[[for_profs[n]]]$EIC[m,"use"] == 0) next	
				found <- (found + 1)				
				tot <- links_profiles[[for_profs[n]]]$EIC[m,"link counts"][[1]]
				links_profiles[[for_profs[n]]]$EIC[m,"link counts"] <- sum(links_profiles[[for_profs[n]]]$EIC_cor[[m]]>=cut_EIC)
				links_profiles[[for_profs[n]]]$EIC[m,"no-link counts"] <- (tot-links_profiles[[for_profs[n]]]$EIC[m,"link counts"])
				prof2 <- links_profiles[[for_profs[n]]]$EIC[m,"linked profile"]
				# which peaks co-occur in the same sample?
				these <- profileList[["peaks"]][
					profileList[["index_prof"]][prof1,"start_ID"]:profileList[["index_prof"]][prof1,"end_ID"]
				,"sampleIDs"]
				those <- profileList[["peaks"]][
					profileList[["index_prof"]][prof2,"start_ID"]:profileList[["index_prof"]][prof2,"end_ID"]
				,"sampleIDs"]
				matched <- match(these,those)
				RT_1 <- (profileList[["peaks"]][
						(profileList[["index_prof"]][prof1,"start_ID"]:profileList[["index_prof"]][prof1,"end_ID"])
						[!is.na(matched)]
					,"RT"])
				RT_2 <- (profileList[["peaks"]][
						(profileList[["index_prof"]][prof2,"start_ID"]:profileList[["index_prof"]][prof2,"end_ID"])
						[matched[!is.na(matched)]]
					,"RT"])
				del_RT <- abs(RT_1 - RT_2)
				links_profiles[[for_profs[n]]]$EIC[m,"ref_1"] <- length(del_RT)
				links_profiles[[for_profs[n]]]$EIC[m,"ref_2"] <- sum(del_RT <= cut_delRT_EIC)
				int_1 <- (profileList[["peaks"]][
						(profileList[["index_prof"]][prof1,"start_ID"]:profileList[["index_prof"]][prof1,"end_ID"])
						[!is.na(matched)]
					,"intensity"])
				int_2 <- (profileList[["peaks"]][
						(profileList[["index_prof"]][prof2,"start_ID"]:profileList[["index_prof"]][prof2,"end_ID"])
						[matched[!is.na(matched)]]
					,"intensity"])
				links_profiles[[for_profs[n]]]$EIC[m,"int_ratio"] <- mean(int_1 / int_2)
				if(links_profiles[[for_profs[n]]]$EIC[m,"ref_1"] > 2){
					int_cor <- cor(int_1, int_2)
					links_profiles[[for_profs[n]]]$EIC[m,"correl"] <- (round(int_cor * 1000))
				}
				at_entry_2 <- profileList[["index_prof"]][prof2,"links"]
				here2 <- which(links_profiles[[at_entry_2]][["EIC"]][,"linked profile"] == prof1)
				if(
					((links_profiles[[for_profs[n]]]$EIC[m,"ref_2"] / links_profiles[[for_profs[n]]]$EIC[m,"ref_1"]) < cut_frac_EIC) |
					(links_profiles[[for_profs[n]]]$EIC[m,"link counts"] / (links_profiles[[for_profs[n]]]$EIC[m,"link counts"] + links_profiles[[for_profs[n]]]$EIC[m,"no-link counts"]) < cut_frac_EIC)
					
				){ # delete this and partner profile ###############
					keep[m] <- FALSE					
					links_profiles[[at_entry_2]]$EIC <- links_profiles[[at_entry_2]]$EIC[-here2,, drop = FALSE]			
					links_profiles[[at_entry_2]]$EIC_cor <- links_profiles[[at_entry_2]]$EIC_cor[-here2]		
					removed<-(removed+1)					
				}else{ # complete other partner profile entry ######
					links_profiles[[at_entry_2]]$EIC[here2,"ref_1"] <- links_profiles[[for_profs[n]]]$EIC[m, "ref_1"]
					links_profiles[[at_entry_2]]$EIC[here2,"ref_2"] <- links_profiles[[for_profs[n]]]$EIC[m, "ref_2"]
					links_profiles[[at_entry_2]]$EIC[here2,"link counts"] <- links_profiles[[for_profs[n]]]$EIC[m, "link counts"]
					links_profiles[[at_entry_2]]$EIC[here2,"no-link counts"] <- links_profiles[[for_profs[n]]]$EIC[m, "no-link counts"]
					links_profiles[[at_entry_2]]$EIC[here2,"int_ratio"] <- (1/links_profiles[[for_profs[n]]]$EIC[m, "int_ratio"])		
					links_profiles[[at_entry_2]]$EIC[here2,"correl"] <- links_profiles[[for_profs[n]]]$EIC[m, "correl"]					
				}
			}# for m
			links_profiles[[for_profs[n]]]$EIC <- links_profiles[[for_profs[n]]]$EIC[keep,,drop=FALSE]
			links_profiles[[for_profs[n]]]$EIC_cor <- links_profiles[[for_profs[n]]]$EIC_cor[keep]
		}
	}
	perce <- as.character(round(removed / found, digits = 3)*100)
	cat("\n");cat(paste0(perce, "% of ", found," EIC links filtered."));cat("\n");
	################################################################
	return(links_profiles)
	
}
