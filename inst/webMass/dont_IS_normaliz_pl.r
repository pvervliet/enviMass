if(file.exists(file.path(as.character(logfile[[1]]),"results","links_profiles_pos"))){file.remove(file.path(as.character(logfile[[1]]),"results","links_profiles_pos"))}
if(file.exists(file.path(as.character(logfile[[1]]),"results","links_profiles_neg"))){file.remove(file.path(as.character(logfile[[1]]),"results","links_profiles_neg"))}
if(file.exists(file.path(as.character(logfile[[1]]),"results","int_norm_ISTD_pos"))){file.remove(file.path(as.character(logfile[[1]]),"results","int_norm_ISTD_pos"))}
if(file.exists(file.path(as.character(logfile[[1]]),"results","int_norm_ISTD_neg"))){file.remove(file.path(as.character(logfile[[1]]),"results","int_norm_ISTD_neg"))}
if(any(objects(envir=as.environment(".GlobalEnv"))=="profileList_pos")){rm(profileList_pos,envir=as.environment(".GlobalEnv"))}
if(any(objects()=="profileList_pos")){rm(profileList_pos)}
if(any(objects(envir=as.environment(".GlobalEnv"))=="links_peaks_pos")){rm(links_peaks_pos,envir=as.environment(".GlobalEnv"))}
if(any(objects()=="links_peaks_pos")){rm(links_peaks_pos)}
if(any(objects(envir=as.environment(".GlobalEnv"))=="links_profiles_pos")){rm(links_profiles_pos,envir=as.environment(".GlobalEnv"))}
if(any(objects()=="links_profiles_pos")){rm(links_profiles_pos)}
if(any(objects(envir=as.environment(".GlobalEnv"))=="profileList_neg")){rm(profileList_neg,envir=as.environment(".GlobalEnv"))}
if(any(objects()=="profileList_neg")){rm(profileList_neg)}
if(any(objects(envir=as.environment(".GlobalEnv"))=="links_peaks_neg")){rm(links_peaks_neg,envir=as.environment(".GlobalEnv"))}
if(any(objects()=="links_peaks_neg")){rm(links_peaks_neg)}
if(any(objects(envir=as.environment(".GlobalEnv"))=="links_profiles_neg")){rm(links_profiles_neg,envir=as.environment(".GlobalEnv"))}
if(any(objects()=="links_profiles_neg")){rm(links_profiles_neg)}
output$int_norm_ISTD_pos_median <- renderPlot({
plot.new()
})
output$int_norm_ISTD_pos_counts <- renderPlot({
plot.new()
})
output$int_norm_ISTD_neg_median <- renderPlot({
plot.new()
})
output$int_norm_ISTD_neg_counts <- renderPlot({
plot.new()
})
if(file.exists(file.path(as.character(logfile[[1]]), "results", "int_norm_ISTD_pos"))){
file.remove(file.path(as.character(logfile[[1]]), "results", "int_norm_ISTD_pos"))
}
if(file.exists(file.path(as.character(logfile[[1]]), "results", "int_norm_ISTD_neg"))){
file.remove(file.path(as.character(logfile[[1]]), "results", "int_norm_ISTD_neg"))
}
if(
file.exists(file.path(logfile[[1]],"results","profileList_pos"))  &
logfile$workflow[names(logfile$workflow) == "components_profiles"] == "yes"
){
if(any(objects(envir=as.environment(".GlobalEnv"))=="profileList_pos")){rm(profileList_pos,envir=as.environment(".GlobalEnv"))}
if(any(objects()=="profileList_pos")){rm(profileList_pos)}
if(any(objects(envir=as.environment(".GlobalEnv"))=="links_peaks_pos")){rm(links_peaks_pos,envir=as.environment(".GlobalEnv"))}
if(any(objects()=="links_peaks_pos")){rm(links_peaks_pos)}
if(any(objects(envir=as.environment(".GlobalEnv"))=="links_profiles_pos")){rm(links_profiles_pos,envir=as.environment(".GlobalEnv"))}
if(any(objects()=="links_profiles_pos")){rm(links_profiles_pos)}
load(file.path(as.character(logfile[[1]]),"results","profileList_pos"),envir=as.environment(".GlobalEnv"));
load(file.path(as.character(logfile[[1]]),"results","links_peaks_pos"),envir=as.environment(".GlobalEnv"));
assign("links_profiles_pos", list(), envir = as.environment(".GlobalEnv"))
measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
use_entries_profiles <- enviMass::VOLZN1972A(links_profiles_pos)
profileList_pos[["index_prof"]][,"links"] <<- 0
with_bar<-FALSE
if(
(
(logfile$workflow[names(logfile$workflow)=="subtr"]=="no") |
(
((logfile$workflow[names(logfile$workflow)=="IS_screen"]=="yes") &  (logfile$parameters$subtr_IS!="yes")) |
((logfile$workflow[names(logfile$workflow)=="target_screen"]=="yes") & (logfile$parameters$subtr_target!="yes"))
)
) &
(length(links_peaks_pos)>0)
){
if(with_bar){pBar <- txtProgressBar(min = 0, max = dim(profileList_pos[["index_prof"]])[1], style = 3)}
for(i in 1:dim(profileList_pos[["index_prof"]])[1]){
if(with_bar){setTxtProgressBar(pBar, i, title = NULL, label = NULL)}
if(
any(profileList_pos[["peaks"]][
(profileList_pos[["index_prof"]][i,"start_ID"]:profileList_pos[["index_prof"]][i,"end_ID"]),"links"
]!=0)
){
if( profileList_pos[["index_prof"]][i,"links"]==0 ){
if(length(use_entries_profiles)>0){
at_entry<-use_entries_profiles[1]
use_entries<-use_entries_profiles[-1]
}else{
at_entry<-(length(links_profiles_pos)+1)
}
links_profiles_pos[[at_entry]] <<- enviMass::RHHPS2808T(profileList_pos[["index_prof"]][i,"number_peaks_total"][[1]])
names(links_profiles_pos)[at_entry] <<- as.character(i)
profileList_pos[["index_prof"]][i,"links"] <<- at_entry
}else{
at_entry <- profileList_pos[["index_prof"]][i,"links"]
}
for(j in (profileList_pos[["index_prof"]][i,"start_ID"]:profileList_pos[["index_prof"]][i,"end_ID"])){
if(profileList_pos[["peaks"]][j,"links"]!=0){
if(	length(links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[2]]) > 0 ){
for(k in 1:length(links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[2]])){
if( length(links_profiles_pos[[at_entry]][[2]])==0 ){
links_profiles_pos[[at_entry]][[2]] <<-
data.frame(
links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[2]][[k]][[1]],
1,
links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[2]][[k]][[2]],
stringsAsFactors = FALSE
)
names(links_profiles_pos[[at_entry]][[2]]) <<- c("Compound","Counts","max_score")
}else{
at<-which(links_profiles_pos[[at_entry]][[2]][,1] == links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[2]][[k]][[1]])
if(length(at)>0){
links_profiles_pos[[at_entry]][[2]][at,2] <<- (links_profiles_pos[[at_entry]][[2]][at,2]+1)
if(links_profiles_pos[[at_entry]][[2]][at,2] < links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[2]][[k]][[2]] ){
links_profiles_pos[[at_entry]][[2]][at,2] <<- links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[2]][[k]][[2]]
}
}else{
links_profiles_pos[[at_entry]][[2]] <<- data.frame(
c(
links_profiles_pos[[at_entry]][[2]][,1],
links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[2]][[k]][[1]]
),
c(links_profiles_pos[[at_entry]][[2]][,2],1),
c(
links_profiles_pos[[at_entry]][[2]][,3],
links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[2]][[k]][[2]]
),
stringsAsFactors = FALSE
)
names(links_profiles_pos[[at_entry]][[2]]) <<- c("Compound","Counts","max_score")
}
}
}
}
if(	length(links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[1]])>0 ){
for(k in 1:length(links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[1]])){
if( length(links_profiles_pos[[at_entry]][[1]])==0 ){
links_profiles_pos[[at_entry]][[1]] <<-
data.frame(
links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[1]][[k]][[1]],
1,
links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[1]][[k]][[2]],
stringsAsFactors = FALSE
)
names(links_profiles_pos[[at_entry]][[1]]) <<- c("Compound","Counts","max_score")
}else{
at<-which(links_profiles_pos[[at_entry]][[1]][,1] == links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[1]][[k]][[1]])
if(length(at)>0){
links_profiles_pos[[at_entry]][[1]][at,2] <<- (links_profiles_pos[[at_entry]][[1]][at,2]+1)
if(links_profiles_pos[[at_entry]][[1]][at,2] < links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[1]][[k]][[2]] ){
links_profiles_pos[[at_entry]][[1]][at,2] <<- links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[1]][[k]][[2]]
}
}else{
links_profiles_pos[[at_entry]][[1]] <<- data.frame(
c(
links_profiles_pos[[at_entry]][[1]][,1],
links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[1]][[k]][[1]]
),
c(links_profiles_pos[[at_entry]][[1]][,2],1),
c(
links_profiles_pos[[at_entry]][[1]][,3],
links_peaks_pos[[profileList_pos[["peaks"]][j,"links"]]][[1]][[k]][[2]]
),
stringsAsFactors = FALSE
)
names(links_profiles_pos[[at_entry]][[1]]) <<- c("Compound","Counts","max_score")
}
}
}
}
}
}
}
}
if(with_bar){close(pBar)}
}
output$int_norm_ISTD_pos_median <- renderPlot({
plot.new()
})
output$int_norm_ISTD_pos_counts <- renderPlot({
plot.new()
})
save(profileList_pos,file=file.path(as.character(logfile[[1]]),"results","profileList_pos"));
save(links_profiles_pos,file=file.path(as.character(logfile[[1]]),"results","links_profiles_pos"));
if(any(ls()=="profileList_pos")){stop("\n illegal profileList_pos detected _1 in do_IS_normaliz.r!")}
if(any(ls()=="links_profiles_pos")){stop("\n illegal links_profiles_pos detected _1 in do_IS_normaliz.r!")}
rm(links_profiles_pos, profileList_pos, envir=as.environment(".GlobalEnv"))
}
if(
file.exists(file.path(logfile[[1]],"results","profileList_neg"))  &
logfile$workflow[names(logfile$workflow) == "components_profiles"] == "yes"
){
if(any(objects(envir=as.environment(".GlobalEnv"))=="profileList_neg")){rm(profileList_neg,envir=as.environment(".GlobalEnv"))}
if(any(objects()=="profileList_neg")){rm(profileList_neg)}
if(any(objects(envir=as.environment(".GlobalEnv"))=="links_peaks_neg")){rm(links_peaks_neg,envir=as.environment(".GlobalEnv"))}
if(any(objects()=="links_peaks_neg")){rm(links_peaks_neg)}
if(any(objects(envir=as.environment(".GlobalEnv"))=="links_profiles_neg")){rm(links_profiles_neg,envir=as.environment(".GlobalEnv"))}
if(any(objects()=="links_profiles_neg")){rm(links_profiles_neg)}
load(file.path(as.character(logfile[[1]]),"results","profileList_neg"), envir = as.environment(".GlobalEnv"));
load(file.path(as.character(logfile[[1]]),"results","links_peaks_neg"), envir = as.environment(".GlobalEnv"));
assign("links_profiles_neg", list(), envir = as.environment(".GlobalEnv"))
measurements <- read.csv(file = file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
use_entries_profiles <- enviMass::VOLZN1972A(links_profiles_neg)
profileList_neg[["index_prof"]][,"links"] <<- 0
with_bar<-FALSE
if(
(
(logfile$workflow[names(logfile$workflow)=="subtr"]=="no") |
(
((logfile$workflow[names(logfile$workflow)=="IS_screen"]=="yes") &  (logfile$parameters$subtr_IS!="yes")) |
((logfile$workflow[names(logfile$workflow)=="target_screen"]=="yes") & (logfile$parameters$subtr_target!="yes"))
)
) &
(length(links_peaks_neg)>0)
){
if(with_bar){pBar <- txtProgressBar(min = 0, max = dim(profileList_neg[["index_prof"]])[1], style = 3)}
for(i in 1:dim(profileList_neg[["index_prof"]])[1]){
if(with_bar){setTxtProgressBar(pBar, i, title = NULL, label = NULL)}
if(
any(profileList_neg[["peaks"]][
(profileList_neg[["index_prof"]][i,"start_ID"]:profileList_neg[["index_prof"]][i,"end_ID"]),"links"
]!=0)
){
if( profileList_neg[["index_prof"]][i,"links"]==0 ){
if(length(use_entries_profiles)>0){
at_entry<-use_entries_profiles[1]
use_entries <- use_entries_profiles[-1]
}else{
at_entry <- (length(links_profiles_neg)+1)
}
links_profiles_neg[[at_entry]] <<- enviMass::RHHPS2808T(profileList_neg[["index_prof"]][i,"number_peaks_total"][[1]])
names(links_profiles_neg)[at_entry] <<- as.character(i)
profileList_neg[["index_prof"]][i,"links"] <<- at_entry
}else{
at_entry <- profileList_neg[["index_prof"]][i,"links"]
}
for(j in (profileList_neg[["index_prof"]][i,"start_ID"]:profileList_neg[["index_prof"]][i,"end_ID"])){
if(profileList_neg[["peaks"]][j,"links"]!=0){
if(	length(links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[2]]) > 0 ){
for(k in 1:length(links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[2]])){
if( length(links_profiles_neg[[at_entry]][[2]])==0 ){
links_profiles_neg[[at_entry]][[2]] <<-
data.frame(
links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[2]][[k]][[1]],
1,
links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[2]][[k]][[2]],
stringsAsFactors = FALSE
)
names(links_profiles_neg[[at_entry]][[2]]) <<- c("Compound","Counts","max_score")
}else{
at<-which(links_profiles_neg[[at_entry]][[2]][,1] == links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[2]][[k]][[1]])
if(length(at)>0){
links_profiles_neg[[at_entry]][[2]][at,2] <<- (links_profiles_neg[[at_entry]][[2]][at,2]+1)
if(links_profiles_neg[[at_entry]][[2]][at,2] < links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[2]][[k]][[2]] ){
links_profiles_neg[[at_entry]][[2]][at,2] <<- links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[2]][[k]][[2]]
}
}else{
links_profiles_neg[[at_entry]][[2]] <<- data.frame(
c(
links_profiles_neg[[at_entry]][[2]][,1],
links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[2]][[k]][[1]]
),
c(links_profiles_neg[[at_entry]][[2]][,2],1),
c(
links_profiles_neg[[at_entry]][[2]][,3],
links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[2]][[k]][[2]]
),
stringsAsFactors = FALSE
)
names(links_profiles_neg[[at_entry]][[2]]) <<- c("Compound","Counts","max_score")
}
}
}
}
if(	length(links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[1]])>0 ){
for(k in 1:length(links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[1]])){
if( length(links_profiles_neg[[at_entry]][[1]])==0 ){
links_profiles_neg[[at_entry]][[1]] <<-
data.frame(
links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[1]][[k]][[1]],
1,
links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[1]][[k]][[2]],
stringsAsFactors = FALSE
)
names(links_profiles_neg[[at_entry]][[1]]) <<- c("Compound","Counts","max_score")
}else{
at<-which(links_profiles_neg[[at_entry]][[1]][,1] == links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[1]][[k]][[1]])
if(length(at)>0){
links_profiles_neg[[at_entry]][[1]][at,2] <<- (links_profiles_neg[[at_entry]][[1]][at,2]+1)
if(links_profiles_neg[[at_entry]][[1]][at,2] < links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[1]][[k]][[2]] ){
links_profiles_neg[[at_entry]][[1]][at,2] <<- links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[1]][[k]][[2]]
}
}else{
links_profiles_neg[[at_entry]][[1]] <<- data.frame(
c(
links_profiles_neg[[at_entry]][[1]][,1],
links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[1]][[k]][[1]]
),
c(links_profiles_neg[[at_entry]][[1]][,2],1),
c(
links_profiles_neg[[at_entry]][[1]][,3],
links_peaks_neg[[profileList_neg[["peaks"]][j,"links"]]][[1]][[k]][[2]]
),
stringsAsFactors = FALSE
)
names(links_profiles_neg[[at_entry]][[1]]) <<- c("Compound","Counts","max_score")
}
}
}
}
}
}
}
}
if(with_bar){close(pBar)}
}
output$int_norm_ISTD_neg_median <- renderPlot({
plot.new()
})
output$int_norm_ISTD_neg_counts <- renderPlot({
plot.new()
})
save(profileList_neg,file=file.path(as.character(logfile[[1]]),"results","profileList_neg"));
save(links_profiles_neg,file=file.path(as.character(logfile[[1]]),"results","links_profiles_neg"));
if(any(ls()=="profileList_neg")){stop("\n illegal profileList_neg detected _1 in do_IS_normaliz.r!")}
if(any(ls()=="links_profiles_neg")){stop("\n illegal links_profiles_neg detected _1 in do_IS_normaliz.r!")}
rm(links_profiles_neg, profileList_neg, envir=as.environment(".GlobalEnv"))
}
