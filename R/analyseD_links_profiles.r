#' @title Checks if profile link entry is empty
#'
#' @export
#'
#' @description Generate a new empty link list
#'
#' @param links_profiles links_profiles_pos or links_profiles_neg
#' @param for_which
#' 
#' @details enviMass workflow function. 
#' 

analyseC_links_profiles<-function(links_profiles, at_entry){

	################################################################
	is_empty<-TRUE
	if(length(links_profiles[[at_entry]]$targ)==0){
		if(length(links_profiles[[at_entry]]$IS)==0){
			if(length(links_profiles[[at_entry]]$EIC)==0){
				if(length(links_profiles[[at_entry]]$isot)==0){
					if(length(links_profiles[[at_entry]]$adduc)==0){
						if(length(links_profiles[[at_entry]]$homol)==0){
							is_empty<-FALSE
						}						
					}else{
						is_empty<-FALSE
					}				
				}else{
					is_empty<-FALSE
				}			
			}else{
				is_empty<-FALSE
			}		
		}else{
			is_empty<-FALSE
		}	
	}else{
		is_empty<-FALSE
	}
	################################################################
	return(is_empty)
	
}
