#' @title Exlude other component profiles from sorted profileList
#'
#' @export
#'
#' @description Exlude other component profiles from sorted profileList
#'
#' @param profileList_index. profileList_pos[["index_prof"]] or profileList_neg[["index_prof"]]
#' @param links_profiles. links_profiles_pos or links_profiles_neg
#' @param sort_what. Vector with names of column in profileList[["index_prof"]] to sort with.
#' @param sort_decreasing. Logical. Sort decreasing (TRUE) or increasing (FALSE)?
#' @param use_profile. NULL or logical vector of length of profiles; FALSE excludes individual profiles.
#' @param with_bar. Show progress bar?
#' @param return_excl. Set to TRUE to return IDs of excluded instead or retained profiles.
#' 
#' @details enviMass workflow function. 
#' 

analyseE_links_profiles<-function(profileList_index, links_profiles, sort_what="deltaint_newest", sort_decreasing = TRUE, use_profile = NULL, with_bar = FALSE, return_excl=FALSE){

	if(any(is.na(match(sort_what,colnames(profileList_index))))){stop("\nFunction analyseE_links_profiles: wrong sort_what - debug!")}	
	################################################################	
	those_profiles<-profileList_index[,"profile_ID"]
	keep_out<-rep(TRUE,length(those_profiles))
	if(!is.null(use_profile)){
		if(length(use_profile)!=dim(profileList_index)[1]){
			stop("\nArgument use_profile must be of same length as number of profiles - abort.")
		}
		keep_out[!use_profile]<-FALSE
	}else{
		use_profile<-rep(TRUE,length(those_profiles))
	}
	if(sort_decreasing){
		along<-rev(do.call(order,as.data.frame(profileList_index[,sort_what,drop=FALSE]))) # for multiple columns
	}else{
		along<-do.call(order,as.data.frame(profileList_index[,sort_what,drop=FALSE])) # for multiple columns	
	}
	################################################################
	if(with_bar){pBar <- txtProgressBar(min = 0, max = length(along), style = 3)}
	for(i in 1:length(along)){
		if(with_bar){setTxtProgressBar(pBar, i, title = NULL, label = NULL)}
		if(keep_out[along[i]]==FALSE){next}
		if(use_profile[along[i]]==FALSE){next}		
		if(profileList_index[along[i],"links"]!=0){
			at_entry<-profileList_index[along[i],"links"]
			if(length(links_profiles[[at_entry]][["group"]])>0){
				keep_out[links_profiles[[at_entry]][["group"]]]<-FALSE
			}
		}
	}
	if(with_bar){close(pBar)}
	################################################################
	# print reduction factor #######################################
	cat(paste0("Reduction factor from profile grouping: ",as.character(round((length(keep_out)/sum(keep_out)),digits=3))))
	################################################################
	if(return_excl){
		return(those_profiles[!keep_out])
	}else{
		return(those_profiles[keep_out])
	}
	
}
