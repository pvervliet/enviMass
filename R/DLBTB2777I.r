DLBTB2777I <- function(
list1,
list2,
as_pairs=TRUE
){
if(length(list1)!=length(list2) & as_pairs){stop("comp_list: arguments not of same length!")}
if(as_pairs){
for(i in 1:length(list1)){
if(list1[[i]] != list2[[i]]){
return(FALSE)
}
}
}else{
for(i in 1:length(list1)){
if(!any(list1[[i]] == list2)){
return(FALSE)
}
}
for(i in 1:length(list2)){
if(!any(list2[[i]] == list1)){
return(FALSE)
}
}
}
return(TRUE)
}
