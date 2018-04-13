KWICA7780L <- function(
peaklist,
mz,
tolmz=10,
ppm=TRUE,
ret,
tolret,
what="mass",
one=TRUE,
knot=5,
plotit=FALSE,
path_1=FALSE,
path_2=FALSE,
stopit=FALSE,
intermediate_results=FALSE,
plot_ppm=FALSE,
max_recal=FALSE
){
if(what!="mass" && what!="ret"){stop("what what?")};
peaks<-BYRIB9287R(peaklist,mz,tolmz,ppm,ret,tolret);
getit1<-c();
getit2<-c();
for(i in 1:length(peaks)){
if(peaks[i]!="FALSE"){
put<-as.numeric(strsplit(peaks[i],"/")[[1]]);
if(what=="mass"){
if(one==TRUE & length(put)==1){
getit1<-c(getit1,rep(mz[i],length(put)));
getit2<-c(getit2,peaklist[put,1]);
}else{
getit1<-c(getit1,rep(mz[i],length(put)));
getit2<-c(getit2,peaklist[put,1]);
}
}else{
if(one==TRUE & length(put)==1){
getit1<-c(getit1,rep(ret[i],length(put)));
getit2<-c(getit2,peaklist[put,3]);
}else{
getit1<-c(getit1,rep(ret[i],length(put)));
getit2<-c(getit2,peaklist[put,3]);
}
}
}
}
getit3<-c(getit1-getit2);
if(length(getit3)<15){
if(path_1!="FALSE"){png(filename = path_1, bg = "white")}
plot.new();
plot.window(xlim=c(1,1),ylim=c(1,1));
box();
text(1,1,label="not available",cex=1.5,col="darkred")
if(path_1!="FALSE"){dev.off()}
if(!stopit){
return("Too few data points for fit!\n");
}else{
stop("Too few data points for fit!\n")
}
}
that<-data.frame(getit2,getit3)
names(that)<-c("obs","delta")
model<-mgcv::gam(delta~s(obs,bs="ts",k=knot),data=that);
if(plotit==TRUE){
if(what=="mass"){
if(path_1!="FALSE"){png(filename = path_1, bg = "white")}
ylim<-c(min(getit3)*1000,max(getit3)*1000)
if(ylim[1]>0){ylim[1]<-0}
if(ylim[2]<0){ylim[2]<-0}
plot(getit2,getit3*1000,pch=19,cex=0.5,xlab="m/z",ylab="Expected m/z - observed m/z [mmu]",
main="Recalibration results",ylim=ylim);
abline(h=0,col="darkgreen");
points(getit2[order(getit2)],predict(model)[order(getit2)]*1000,col="red",type="l",lwd=2);
if(plot_ppm[1]!="FALSE"){
ppm_mass<-seq(0,max(getit2),10)
for(k in 1:length(plot_ppm)){
ppm_ppm<-(ppm_mass*plot_ppm[k]/1E6*1000)
lines(ppm_mass,ppm_ppm,lty=2,col="gray")
lines(ppm_mass,-ppm_ppm,lty=2,col="gray")
plotmass<-median(ppm_mass)
text(
plotmass,ppm_ppm[ppm_mass==plotmass],
labels=paste(as.character(plot_ppm[k]),"ppm"),
col="gray",cex=1.2
)
text(
plotmass,-ppm_ppm[ppm_mass==plotmass],
labels=paste("-",as.character(plot_ppm[k])," ppm",sep=""),
col="gray",cex=1.2
)
}
}
if(max_recal!="FALSE"){
if(ppm){
ppm_mass<-seq(0,max(getit2),10)
ppm_ppm<-(ppm_mass*max_recal/1E6*1000)
lines(ppm_mass,ppm_ppm,lty=2,lwd=1.5,col="red")
lines(ppm_mass,-ppm_ppm,lty=2,lwd=1.5,col="red")
}else{
abline(h=max_recal,lty=2,lwd=1.5,col="red")
abline(h=-max_recal,lty=2,lwd=1.5,col="red")
}
}
if(path_1!="FALSE"){dev.off()}
}else{
if(path_1!="FALSE"){png(filename = path_1, bg = "white")}
plot(getit2,getit3,pch=19,cex=0.5,xlab="Retention time",ylab="Expected RT - observed RT",
main="Recalibration results");
abline(h=0,col="red");
points(getit2[order(getit2)],predict(model)[order(getit2)],col="red",type="l",lwd=2);
if(path_1!="FALSE"){dev.off()}
}
}
if(what=="mass"){
that<-data.frame("obs"=peaklist[,1],"delta"=peaklist[,1]);
pred2<-mgcv::predict.gam(model,newdata=that);
newpeaks<-peaklist;
if(max_recal=="FALSE"){
newpeaks[,1]<-c(peaklist[,1]+pred2);
}else{
if(ppm){
if(!any((pred2/that[,1]*1E6)>max_recal)){
newpeaks[,1]<-c(peaklist[,1]+pred2);
}else{
cat("\n recalibration skipped - correction off limits!")
}
}else{
if(!any(abs(pred2)>max_recal)){
newpeaks[,1]<-c(peaklist[,1]+pred2);
}else{
cat("\n recalibration skipped - correction off limits!")
}
}
}
}else{
that<-data.frame("obs"=peaklist[,3],"delta"=peaklist[,3]);
pred2<-mgcv::predict.gam(model,newdata=that);
newpeaks<-peaklist;
newpeaks[,3]<-c(peaklist[,3]+pred2);
}
if(path_2!="FALSE"){
save(model,file=path_2)
}
if (intermediate_results) {
imr <- list()
imr[[1]] <- model;
imr[[2]] <- peaks;
imr[[3]] <- getit2;
imr[[4]] <- getit3;
names(imr) <- c("model", "matches", "x", "delta_y");
eval.parent(substitute(intermediate_results<-imr))
}
return(newpeaks)
}
