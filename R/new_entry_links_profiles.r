#' @title Generate a new empty link list
#'
#' @export
#'
#' @description Generate a new empty link list
#'
#' @param a Total peak number of the profile
#' 
#' @details enviMass workflow function. 
#' 

new_entry_links_profiles<-function(a){

	################################################################
	list_entry<-list()
	list_entry[[1]]<-list() # target
	list_entry[[2]]<-list()	# IS
	list_entry[[3]]<-list()	# EIC_correl
		list_entry[[3]]<-matrix(ncol=8,nrow=0) 
		# ID of linked profile / 
		# counts for peaks with correlated EICs /
		# counts for peaks with non-correlated EICs / 
		# ref_1: number of peaks co-occuring over samples
		# ref_2: number of peaks co-occuring, within delRT used for isotopologues
		# ref_3: number of peaks co-occuring, within delRT used for adducts		
		# use: use this profile to store EICs or to compare relations? avoid redundant and dispersed entries!
		colnames(list_entry[[3]])<-c("linked profile","link counts","no-link counts","ref_1","ref_2","ref_3","use","correl")
	list_entry[[4]]<-list()	# isotop
		list_entry[[4]]<-matrix(ncol=5,nrow=0)
		colnames(list_entry[[4]])<-c("linked profile","link counts","ref_1","use","correl")					
	list_entry[[5]]<-list()	# adducts
		list_entry[[5]]<-matrix(ncol=5,nrow=0)
		colnames(list_entry[[5]])<-c("linked profile","link counts","ref_1","use","correl")					
	list_entry[[6]]<-list()	# homol		
		list_entry[[6]]<-matrix(ncol=3,nrow=0)
		colnames(list_entry[[6]])<-c("linked profile","link counts","ref_1")
	list_entry[[7]]<-list()	# group
	list_entry[[8]]<-list()	# total number of peaks		
	list_entry[[8]][[1]]<-a
	list_entry[[9]]<-list() # EIC correlation
	names(list_entry)<-c("targ","IS","EIC","isot","adduc","homol","group","total","EIC_cor")
	################################################################
	return(list_entry)
	
}
