#' @title Extract non-isotopologue / adduct characteristics
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

analyseC_links_profiles<-function(links_profiles, profileList, for_which="all"){

	################################################################
	if(for_which=="ISTD"){
		for_profs<-which(unlist(lapply(lapply(links_profiles, `[[`, 2),length))!=0)
	}
	if(for_which=="target"){
		for_profs<-which(unlist(lapply(lapply(links_profiles, `[[`, 1),length))!=0)
	}	
	if(for_which=="both"){
		for_profs1<-which(unlist(lapply(lapply(links_profiles, `[[`, 1),length))!=0)
		for_profs2<-which(unlist(lapply(lapply(links_profiles, `[[`, 2),length))!=0)
		for_profs<-unique(c(for_profs1,for_profs2))
	}
	if(for_which=="all"){
		for_profs<-(1:length(links_profiles))
	}	
	################################################################
	del_mz<-c()
	at_mz1<-c()
	at_mz2<-c()
	at_RT<-c()
	at_cor<-c()
	at_rat<-c()
	EIC_cor_adduc<-c()
	for(n in 1:length(for_profs)){
		if(dim(links_profiles[[for_profs[n]]]$EIC)[1]==0) next
		for(m in 1:dim(links_profiles[[for_profs[n]]]$EIC)[1]){
			if(links_profiles[[for_profs[n]]]$EIC[m,"use"]==0) next	
			prof1<-as.numeric(names(links_profiles)[for_profs[n]])
			if(dim(links_profiles[[for_profs[n]]]$isot)[1]>0){	
				if( any(links_profiles[[for_profs[n]]]$isot[,"linked profile"] == prof1) ) next
			}
			if(dim(links_profiles[[for_profs[n]]]$adduc)[1]>0){
				if( any(links_profiles[[for_profs[n]]]$adduc[,"linked profile"] == prof1) ) next
			}
			prof2<-links_profiles[[for_profs[n]]]$EIC[m,"linked profile"]			
			del_mz<-c(del_mz,abs(
				(profileList[["index_prof"]][prof1,"mean_mz"])-(profileList[["index_prof"]][prof2,"mean_mz"])
			)[[1]])
			at_mz1<-c(at_mz1,profileList[["index_prof"]][prof1,"mean_mz"][[1]])
			at_mz2<-c(at_mz2,profileList[["index_prof"]][prof1,"mean_mz"][[1]])
			at_RT<-c(at_RT,profileList[["index_prof"]][prof1,"mean_RT"][[1]])
			at_cor<-c(at_cor,links_profiles[[for_profs[n]]]$EIC[m,"correl"][[1]])	
			at_rat<-c(at_rat,links_profiles[[for_profs[n]]]$EIC[m,"int_ratio"][[1]])				
		}
	}
	################################################################
	results<-as.matrix(cbind(at_mz1,at_mz2,at_RT,del_mz,(at_cor/1000),at_rat,deparse.level = 0))
	colnames(results)<-c("m/z1","m/z2","RT","del_m/z","cor","int_rat")
	return(results)
	
}
