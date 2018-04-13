VOLZN1972A <- function(
any_list
){
if(!is.list(any_list)){stop("any_list mut be a list")}
if(length(any_list)==0){return(c())}
positions<-c()
for(i in 1:length(any_list)){
if(length(any_list[[i]])==0){
positions<-c(positions,i)
}
}
return(positions)
}
