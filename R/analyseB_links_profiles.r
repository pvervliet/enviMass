#' @title Analyse profiles - retrieve ISTD characteristics (2)
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

analyseB_links_profiles<-function(links_profiles, min_count=.4, for_which="ISTD"){

	################################################################
	num_sam<-(length(profileList[["sampleID"]])*min_count)
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
	EIC_cor_isot<-c()
	EIC_cor_adduc<-c()
	for(n in 1:length(for_profs)){
		if(links_profiles[[for_profs[n]]]$total<num_sam) next
		if(dim(links_profiles[[for_profs[n]]]$EIC)[1]==0) next
		if((dim(links_profiles[[for_profs[n]]]$adduc)[1]==0) & (dim(links_profiles[[for_profs[n]]]$isot)[1]==0)) next
		prof1<-as.numeric(names(links_profiles)[for_profs[n]])
		# isotopologue EIC correlation #############################
		if(dim(links_profiles[[for_profs[n]]]$isot)[1]>0){
			those<-match(links_profiles[[for_profs[n]]]$isot[,"linked profile"],links_profiles[[for_profs[n]]]$EIC[,"linked profile"])
			if(any(!is.na(those))){
				for(i in 1:length(those)){
					if(is.na(those[i])) next
					if(links_profiles[[for_profs[n]]]$isot[i,"use"]==0) next
					EIC_cor_isot<-c(EIC_cor_isot,links_profiles[[for_profs[n]]]$EIC_cor[[those[i]]])
				}
			}
		}		
		# adduct EIC correlation ###################################
		if(dim(links_profiles[[for_profs[n]]]$adduc)[1]>0){
			those<-match(links_profiles[[for_profs[n]]]$adduc[,"linked profile"],links_profiles[[for_profs[n]]]$EIC[,"linked profile"])
			if(any(!is.na(those))){
				for(i in 1:length(those)){
					if(is.na(those[i])) next
					if(links_profiles[[for_profs[n]]]$adduc[i,"use"]==0) next
					EIC_cor_adduc<-c(EIC_cor_adduc,links_profiles[[for_profs[n]]]$EIC_cor[[those[i]]])
				}
			}
		}		
	}
	################################################################
	listed<-list()
	if(length(EIC_cor_isot)) listed[[1]]<-EIC_cor_isot else listed[[1]]<-NA
	if(length(EIC_cor_adduc)) listed[[2]]<-EIC_cor_adduc else listed[[2]]<-NA
	names(listed)<-c("EIC_cor_isot","EIC_cor_adduc")
	return(listed)
	
}
