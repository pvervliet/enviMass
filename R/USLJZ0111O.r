USLJZ0111O <- function(IDs){
if(!is.numeric(IDs)){stop("IDs not numeric")}
IDs<-IDs[order(IDs,decreasing=FALSE)]
lowestID<-(0)
for(i in 1:length(IDs)){
if(!any(IDs==i)){
lowestID<-i;
break;
}
}
if(lowestID==0){
lowestID<-(max(IDs)+1);
}
return(lowestID)
}
