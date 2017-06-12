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
		list_entry[[3]]<-matrix(ncol=6,nrow=0)
		colnames(list_entry[[3]])<-c("linked profile","link counts","no-link counts","ref_1","ref_2","ref_3")
	list_entry[[4]]<-list()	# isotop
		list_entry[[4]]<-matrix(ncol=3,nrow=0)
		colnames(list_entry[[4]])<-c("linked profile","link counts","ref_1")					
	list_entry[[5]]<-list()	# adducts
		list_entry[[5]]<-matrix(ncol=3,nrow=0)
		colnames(list_entry[[5]])<-c("linked profile","link counts","ref_1")					
	list_entry[[6]]<-list()	# homol		
		list_entry[[6]]<-matrix(ncol=3,nrow=0)
		colnames(list_entry[[6]])<-c("linked profile","link counts","ref_1")
	list_entry[[7]]<-list()	# group
	list_entry[[8]]<-list()	# total number of peaks		
	list_entry[[8]][[1]]<-a
	names(list_entry)<-c("targ","IS","EIC","isot","adduc","homol","group","total")
	################################################################
	return(list_entry)
	
}
