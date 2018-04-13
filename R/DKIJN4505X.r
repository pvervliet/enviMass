DKIJN4505X <- function(
forID,
what,
...
){
if(what == "EICs"){
if(file.exists(file.path(logfile[[1]], "results", "componentization", "EIC_corr", as.character(forID)))){
load(file = file.path(logfile[[1]], "results", "componentization", "EIC_corr", as.character(forID)))
}else{
return("failed")
}
if(length(EIC_pairs[,1]) == 0) return("failed")
get1 <- cbind(
rep(forID, length(EIC_pairs[,1])), EIC_pairs[,1]
)
found1<-enviMass::GSZEL3799Y(get1,peaks[,c("sampleIDs","peakIDs")],row_order=FALSE,column_order_a=FALSE,column_order_b=FALSE,get_index=TRUE)
get2<-cbind(
rep(forID,length(EIC_pairs[,2])),EIC_pairs[,2]
)
found2<-enviMass::GSZEL3799Y(get2, peaks[,c("sampleIDs","peakIDs")],row_order=FALSE,column_order_a=TRUE,column_order_b=FALSE,get_index=TRUE)
found <- cbind(found1, found2, EIC_pairs[,4])
return(found)
}
if(what == "isotopologues"){
if(file.exists(file.path(logfile[[1]], "results", "componentization", "isotopologues", as.character(forID)))){
load(file = file.path(logfile[[1]], "results", "componentization", "isotopologues", as.character(forID)))
}else{
return("failed")
}
if(length(Isot_pairs[,1]) == 0) return("failed")
get1 <- cbind(
rep(forID, length(Isot_pairs[,1])), Isot_pairs[,1]
)
found1 <- enviMass::GSZEL3799Y(get1, peaks[,c("sampleIDs", "peakIDs")], row_order = FALSE, column_order_a = FALSE, column_order_b = FALSE, get_index = TRUE)
get2 <- cbind(
rep(forID, length(Isot_pairs[,2])), Isot_pairs[,2]
)
found2 <- enviMass::GSZEL3799Y(get2, peaks[,c("sampleIDs", "peakIDs")], row_order = FALSE, column_order_a = TRUE, column_order_b = FALSE, get_index = TRUE)
found <- cbind(found1, found2)
return(found)
}
if(what == "adducts"){
if(file.exists(file.path(logfile[[1]], "results", "componentization", "adducts", as.character(forID)))){
load(file = file.path(logfile[[1]], "results", "componentization", "adducts", as.character(forID)))
}else{
return("failed")
}
if(length(Adduct_pairs[,1]) == 0) return("failed")
get1<-cbind(
rep(forID,length(Adduct_pairs[,1])),Adduct_pairs[,1]
)
found1<-enviMass::GSZEL3799Y(get1,peaks[,c("sampleIDs","peakIDs")],row_order=FALSE,column_order_a=FALSE,column_order_b=FALSE,get_index=TRUE)
get2<-cbind(
rep(forID,length(Adduct_pairs[,2])),Adduct_pairs[,2]
)
found2<-enviMass::GSZEL3799Y(get2,peaks[,c("sampleIDs","peakIDs")],row_order=FALSE,column_order_a=TRUE,column_order_b=FALSE,get_index=TRUE)
found <- cbind(found1, found2)
return(found)
}
}
