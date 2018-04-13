BYRIB9287R <- function(peaklist,mz,dmz=5,ppm=TRUE,RT,dRT,onlymax=FALSE,int_ratio=FALSE,int=FALSE,get_matches=TRUE,get_ratio=FALSE){
ord1<-order(peaklist[,1],decreasing=TRUE);
peaks<-peaklist[ord1,];
ord2<-order(mz,decreasing=TRUE);
mass<-mz[ord2];
if(length(dmz)==1){dmz<-rep(dmz,length(mz))}
if(length(dRT)==1){dRT<-rep(dRT,length(mz))}
if(!get_ratio){
result<-rep("FALSE",length(mz));
}else{
result<-rep(Inf,length(mz));
}
leng<-length(peaks[,1]);
k<-c(1);
options(digits=10)
for(i in 1:length(mass)){
deletes<-c();
if(ppm==TRUE){
target_low<-as.numeric(mass[i])-(as.numeric(dmz[i])*(as.numeric(mass[i])/1e6));
target_up<-as.numeric(mass[i])+(as.numeric(dmz[i])*(as.numeric(mass[i])/1e6));
}else{
target_low<-as.numeric(mass[i])-((as.numeric(dmz[i])/1000));
target_up<-as.numeric(mass[i])+((as.numeric(dmz[i])/1000));
}
while((k-1)>0 ){
if(target_up >= as.numeric(peaks[k-1,1])){
k<-c(k-1);
}else{
break;
}
};
while( (k+1)<=as.numeric(leng) ){
if(target_up < as.numeric(peaks[k+1,1])){
k<-c(k+1);
}else{
break
}
};
n<-k;
while( (n+1)<=as.numeric(leng) ){
if(target_low <= as.numeric(peaks[n+1,1])){
n<-c(n+1);
}else{
break;
}
};
for(f in k:n){
if(
(as.numeric(peaks[f,1]) >= as.numeric(target_low))  &
( as.numeric(peaks[f,1]) <= as.numeric(target_up))  &
( as.numeric(peaks[f,3]) >= (as.numeric(RT[ord2[i]])-as.numeric(dRT[ord2[i]]))) &
( as.numeric(peaks[f,3]) <= (as.numeric(RT[ord2[i]])+as.numeric(dRT[ord2[i]])))
){
deletes<-c(deletes,f);
}
}
if(is.numeric(int_ratio) & length(deletes)>0 & !get_ratio){
deletes<-deletes[
int_ratio>=(int[ord2[i]]/peaklist[ord1[deletes],2])
]
}
if(length(deletes)>0){
if(get_ratio){
if(length(deletes)>1){
deletes<-deletes[which.max(peaklist[ord1[deletes],2])]
}
result[ord2[i]]<-(int[ord2[i]]/peaklist[ord1[deletes],2])
}else{
if(!get_matches){
result[ord2[i]]<-"TRUE"
next;
}
if(onlymax){
if(length(deletes)==1){
result[ord2[i]]<-as.character(ord1[deletes]);
}else{
result[ord2[i]]<-as.character(
(ord1[deletes])[which.max(peaklist[ord1[deletes],2])]
)
}
}else{
result[ord2[i]]<-as.character(ord1[deletes[1]]);
if(length(deletes)>1){
for(j in 2:length(deletes)){
result[ord2[i]]<-paste(result[ord2[i]],"/",as.character(ord1[deletes[j]]))
}
}
}
}
}
}
return(result);
}
