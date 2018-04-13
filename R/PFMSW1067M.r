PFMSW1067M <- function(depend,must){
for(i in 1:length(must[1,])){
if(any(must[,i])==1){
those<-rownames(must)[must[,i]==1]
for_that<-colnames(must)[i]
for(j in 1:length(those)){
if(depend[rownames(depend)==for_that,colnames(depend)==those[j]]==0){
cat(paste("\nVOID DEPENDENCY_1 detected between ",for_that," and ",those[j],sep=""))
}
depend[rownames(depend)==for_that,colnames(depend)==those[j]]<-1
}
}
if(any(must[,i])==2){
those<-rownames(must)[must[,i]==2]
for_that<-colnames(must)[i]
for(j in 1:length(those)){
if(depend[colnames(depend)==those[j],rownames(depend)==for_that]){
cat(paste("\nVOID DEPENDENCY_1 detected between ",those[j]," and ",for_that,sep=""))
}
depend[colnames(depend)==those[j],rownames(depend)==for_that]<-1
}
}
}
for(i in 1:length(depend[1,])){
if(any(depend[,i]==2)){
depend[i,
depend[,i]==2
]<-1
}
}
depend[depend==2]<-0
do_for<-length(depend[1,])
say<-"ok"
node<-c()
node_level<-c()
for(i in 1:do_for){
at<-apply(depend,1,sum)
do_these<-which(at==0)
if(length(node)>0){
do_these<-do_these[
is.na(match(
colnames(depend)[do_these],node
))
]
}
if(length(do_these)==0 & !all(at==0)){
say<-"workflow_scheduler issue_1"
break;
}
node<-c(node,colnames(depend)[do_these])
node_level<-c(node_level,rep(i,length(do_these)))
for(j in 1:length(do_these)){
depend[do_these[j],]<-0
depend[,do_these[j]]<-0
}
}
if(any(is.na(match(colnames(depend),node)))){say <- "workflow_scheduler issue _2"}
if(any(depend>0) & say!="ok"){say <- "workflow_scheduler issue_3"}
if(say == "ok"){
return(data.frame(node, node_level, stringsAsFactors = TRUE))
}else{
return(say)
}
}
