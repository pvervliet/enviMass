JGZUN3462U <- function(strings, check_mute = TRUE){
affected <- c()
affected_TF <- c()
files <- list.files()
for(i in 1:length(files)){
if(
(grepl("do_",files[i])) || (grepl("dont_",files[i]))
){
textit <- readLines(
file.path(files[i])
)
gotit <- FALSE
if(!length(textit)) next;
for(j in 1:length(textit)){
for(k in 1:length(strings)){
if(grepl(strings[k], textit[j], fixed = TRUE)){
if(check_mute){
all_muted <- TRUE
in_here <- gregexpr(strings[k],textit[j])[[1]]
for(n in 1:length(in_here)){
if(in_here[n] < 6){
all_muted <- FALSE;
break;
}
if(substr(textit[j],(in_here[n]-5),(in_here[n]-2))!="mute"){
all_muted <- FALSE
break;
}
}
if(all_muted){next}
}
if(grepl("do_", files[i])){
this <- strsplit(files[i],"do_")[[1]][2]
this <- strsplit(this, ".", fixed = TRUE)[[1]][1]
this <- strsplit(this,"_pl")[[1]][1]
affected <- c(affected, this)
affected_TF <- c(affected_TF, "TRUE")
}else{
this <- strsplit(files[i], "dont_")[[1]][2]
this <- strsplit(this, ".", fixed = TRUE)[[1]][1]
this <- strsplit(this,"_pl")[[1]][1]
affected <- c(affected, this)
affected_TF <- c(affected_TF,"FALSE")
}
gotit<-TRUE;
break;
}
}
if(gotit){
break
}
}
}
}
affected_table <- cbind(as.character(affected), as.character(affected_TF))
colnames(affected_table) <- c("node", "exec")
affected_table <- unique(affected_table)
return(affected_table)
}
