#' @title Exlude other component profiles from sorted profileList
#'
#' @export
#'
#' @description Exlude other component profiles from sorted profileList
#'
#' @param profileList. profileList_pos or profileList_neg
#' @param links_profiles. links_profiles_pos or links_profiles_neg
#' @param sort_what. Vector with names of column in profileList[["index_prof"]] to sort with (decreasing order).
#' @param use_profile. NULL or logical vector of length of profiles; FALSE excludes individual profiles.
#' @param with_bar. Show progress bar?
#' @param return_excl. Set to TRUE to return IDs of excluded instead or retained profiles.
#' 
#' @details enviMass workflow function. 
#' 

analyseE_links_profiles<-function(profileList, links_profiles, sort_what="deltaint_newest", use_profile = NULL, with_bar = FALSE, return_excl=FALSE){

	if(any(is.na(match(sort_what,colnames(profileList[["index_prof"]]))))){stop("\nFunction analyseE_links_profiles: wrong sort_what - debug!")}	
	################################################################	
	those_profiles<-profileList[["index_prof"]][,"profile_ID"]
	keep_out<-rep(TRUE,length(those_profiles))
	if(!is.null(use_profile)){
		if(length(use_profile)!=dim(profileList[["index_prof"]])[1]){
			stop("\nArgument use_profile must be of same length as number of profiles - abort.")
		}
		keep_out[!use_profile]<-FALSE
	}else{
		use_profile<-rep(TRUE,length(those_profiles))
	}
	along<-rev(do.call(order,as.data.frame(profileList[["index_prof"]][,sort_what,drop=FALSE]))) # for multiple columns
	################################################################
	if(with_bar){pBar <- txtProgressBar(min = 0, max = length(along), style = 3)}
	for(i in 1:length(along)){
		if(with_bar){setTxtProgressBar(pBar, i, title = NULL, label = NULL)}
		if(keep_out[along[i]]==FALSE){next}
		if(use_profile[along[i]]==FALSE){next}		
		if(profileList[["index_prof"]][along[i],"links"]!=0){
			at_entry<-profileList[["index_prof"]][along[i],"links"]
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
