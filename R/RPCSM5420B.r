RPCSM5420B <- function(
profileList_index,
links_profiles,
sort_what = "deltaint_newest",
sort_decreasing = TRUE,
use_profile = NULL,
with_bar = FALSE,
return_excl = FALSE
){
if(any(is.na(match(sort_what, colnames(profileList_index))))){stop("\nFunction analyseE_links_profiles: wrong sort_what - debug!")}
those_profiles <- profileList_index[,"profile_ID"]
max_ID <- max(those_profiles)
keep_out <- rep(TRUE, max_ID)
if(!is.null(use_profile)){
if(length(use_profile) != dim(profileList_index)[1]){
stop("\nArgument use_profile must be of same length as number of profiles - abort.")
}
keep_out[those_profiles[!use_profile]] <- FALSE
}
if(sort_decreasing){
along <- rev(do.call(order, as.data.frame(profileList_index[, sort_what, drop = FALSE])))
}else{
along <- do.call(order, as.data.frame(profileList_index[, sort_what, drop = FALSE]))
}
not_done <- rep(TRUE, max_ID)
for(i in 1:length(along)){
not_done[those_profiles[along[i]]] <- FALSE
if(profileList_index[along[i],"links"] == 0) next
at_entry <- profileList_index[along[i], "links"]
if(!length(links_profiles[[at_entry]])) next
if(length(links_profiles[[at_entry]][["group"]]) > 0){
those <- links_profiles[[at_entry]][["group"]]
those <- those[those <= max_ID]
those <- those[not_done[those]]
if(length(those)){
keep_out[those] <- FALSE
}
}
}
keep_out <- keep_out[those_profiles]
cat(paste0("Reduction factor from profile grouping: ", as.character(round((length(keep_out) / sum(keep_out)), digits = 4))))
if(return_excl){
return(those_profiles[!keep_out])
}else{
return(those_profiles[keep_out])
}
}
