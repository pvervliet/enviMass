IXEME0353N <- function(
comp,
compoID,
n_col=FALSE,
masslim=FALSE,
what="circ"
){
if(!any(c("circ", "spec", "table", "check") == what)){stop("Wrong what in plotcomp_parts - debug me!")}
get1a <- strsplit(as.character(comp[[1]][compoID,"ID pattern peaks |"]),",")[[1]];
get1a_blind <- sapply(get1a, grepl, pattern = "*", fixed = TRUE, USE.NAMES = FALSE)
get1a <- as.numeric(sapply(get1a, gsub, pattern = "*", replacement = "",fixed = TRUE, USE.NAMES = FALSE))
if(as.character(comp[[1]][compoID, "ID adduct peaks |"])!="-"){
get1b <- strsplit(as.character(comp[[1]][compoID, "ID adduct peaks |"]),",")[[1]];
get1b_blind <- sapply(get1b, grepl, pattern = "*",fixed = TRUE, USE.NAMES = FALSE)
get1b <- as.numeric(sapply(get1b, gsub, pattern = "*", replacement = "", fixed = TRUE, USE.NAMES = FALSE))
}else{
get1b <- c()
get1b_blind <- c()
}
if(as.character(comp[[1]][compoID, "ID interfering peaks |"]) != "-"){
get1c <- strsplit(as.character(comp[[1]][compoID,"ID interfering peaks |"]),",")[[1]];
get1c_blind <- sapply(get1c, grepl, pattern = "*", fixed = TRUE, USE.NAMES = FALSE)
get1c <- as.numeric(sapply(get1c, gsub, pattern = "*", replacement = "",fixed = TRUE, USE.NAMES = FALSE))
}else{
get1c <- c()
get1c_blind <- c()
}
get1 <- c(get1a, get1b);
get1 <- as.numeric(unique(get1));
get3 <- c(get1a,get1b,get1c);
get3_blind <- c(get1a_blind, get1b_blind, get1c_blind);
dub <- duplicated(get3)
get3 <- get3[!dub]
get3_blind <- get3_blind[!dub]
if(what == "check"){
if(length(get3) < 2){
if(length(get3) < 1){
return("no_peak")
}else{
return("single_peak")
}
}else{
return("available")
}
}
if(length(comp[["pattern peak list"]])>1){
matched <- match(get3,comp[["pattern peak list"]][,"peak ID"])
dat1 <- comp[["pattern peak list"]][matched,];
ord <- rank(1/dat1[,2]);
get3 <- get3[order(dat1[,1],decreasing=FALSE)];
get3_blind <- get3_blind[order(dat1[,1], decreasing=FALSE)]
ord <- ord[order(dat1[,1],decreasing=FALSE)];
}else{
dat1 <- FALSE;
}
if(length(comp[["adduct peak list"]])>1){
matched <- match(get3, comp[["adduct peak list"]][,"peak ID"])
dat2 <- comp[["adduct peak list"]][matched,];
ord <- rank(1 / dat2[,2]);
get3 <- get3[order(dat2[,1], decreasing=FALSE)];
get3_blind <- get3_blind[order(dat2[,1], decreasing=FALSE)]
ord <- ord[order(dat2[,1], decreasing=FALSE)];
}else{
dat2 <- FALSE;
}
if(length(comp[["pattern peak list"]]) > 1){
relat1 <- matrix(ncol=length(get3), nrow = length(get3),"");
rownames(relat1) <- get3;
colnames(relat1) <- get3;
for(i in 1:length(get3)){
that1 <- as.numeric(strsplit(as.character(dat1[dat1[,4]==get3[i],7]),"/")[[1]]);
that2 <- strsplit(as.character(dat1[dat1[,4]==get3[i],8]),"/")[[1]];
that3 <- strsplit(as.character(dat1[dat1[,4]==get3[i],10]),"/")[[1]];
that4 <- strsplit(as.character(dat1[dat1[,4]==get3[i],9]),"/")[[1]];
if(that1[1] != 0){
for(j in 1:length(that1)){
relat1[get3==get3[i],get3==that1[j]]<-paste(relat1[get3==get3[i],get3==that1[j]],"/",that2[j],"(",that4[j],")",",z=",that3[j],sep="");
}
}
}
}
if(length(comp[["adduct peak list"]])>1){
relat2<-matrix(ncol=length(get3),nrow=length(get3),"");
rownames(relat2)<-get3;colnames(relat2)<-get3;
for(i in 1:length(get3)){
that1<-as.numeric(strsplit(as.character(dat2[dat2[,4]==get3[i],6]),"/")[[1]]);
that2<-strsplit(as.character(dat2[dat2[,4]==get3[i],7]),"//")[[1]];
if(that1[1]!=0){
for(j in 1:length(that1)){
relat2[get3==get3[i],get3==that1[j]]<-paste(relat2[get3==get3[i],get3==that1[j]],that2[j],sep="/");
}
}
}
rm(i)
}
if(what=="circ"){
par(mar=c(1,1,1,1));
plot.new();
plot.window(xlim=c(-1.1,1.1),ylim=c(-1.1,1.1));
coordx<-rep(0,length(get3));
coordy<-rep(0,length(get3));
b=0;
a=720/(length(get3)+3)/100;
for(i in 1:length(get3)){
coordx[i]=sin(b);
coordy[i]=cos(b);
b=b+a;
}
rm(i);
for(i in 1:length(get3)){
if(any(get3[i]==get1)){
if(!get3_blind[i]){
text(coordx[i],coordy[i],labels=paste0(get3[i]," (",ord[i],")"),col="darkgreen",cex=.8)
}else{
text(coordx[i],coordy[i],labels=paste0(get3[i],"* (",ord[i],")"),col="darkgreen",cex=.8)
}
}else{
if(!get3_blind[i]){
text(coordx[i],coordy[i],labels=paste0(get3[i]," (",ord[i],")"),col="darkgrey",cex=.8)
}else{
text(coordx[i],coordy[i],labels=paste0(get3[i],"* (",ord[i],")"),col="darkgrey",cex=.8)
}
}
}
rm(i);
coordx<-coordx*0.8;
coordy<-coordy*0.8;
if(length(comp[["pattern peak list"]])>1){
for(i in 1:length(get3)){
for(j in 1:length(get3)){
if(relat1[i,j]!=""){
if(any(get1==get3[i]) & any(get1==get3[j])){
arrows(coordx[i],coordy[i],coordx[j],coordy[j],col="blue",length=.1,lwd=2);
}else{
arrows(coordx[i],coordy[i],coordx[j],coordy[j],col="blue",length=.1,lwd=1);
}
}
}
}
}
if(length(comp[["adduct peak list"]])>1){
for(i in 1:length(get3)){
for(j in 1:length(get3)){
if(relat2[i,j]!=""){
if(any(get1==get3[i]) & any(get1==get3[j])){
lines(c(coordx[i],coordx[j]),c(coordy[i],coordy[j]),col="red",lwd=2);
}else{
lines(c(coordx[i],coordx[j]),c(coordy[i],coordy[j]),col="red",lwd=1);
}
}
}
}
}
if(length(comp[["pattern peak list"]])>1){
matched<-match(get1,comp[["pattern peak list"]][,"peak ID"])
dat3<-comp[["pattern peak list"]][matched,];
that<-dat3[dat3[,2]==max(dat3[,2]),][,4];
points(coordx[match(that,get3)],coordy[match(that,get3)],pch=21,cex=3);
}
lines(c(-0.23,-0.21),c(1.15,1.03),col="darkgrey");
arrows(-0.23,1.15,0.1,1.15,length=0.1,col="darkgrey");
text(0.2,1.15,labels="m/z",col="darkgrey",cex=.8);
text(-1,1.1,labels="Isotopologue links",col="blue",cex=.8,pos=4);
text(-1,1,labels="Adduct links",col="red",cex=.8,pos=4);
}
if(what=="spec"){
if(length(comp[["pattern peak list"]])>1 & length(comp[[3]])>1){
mintol<-c(min(dat1[,3],dat2[,3])-comp[["Parameters"]][1]);
maxtol<-c(max(dat1[,3],dat2[,3])+comp[["Parameters"]][1]);
if(comp[["Parameters"]][4]==TRUE){
minmz<-c(min(dat1[,1],dat2[,1])-(comp[["Parameters"]][2]*min(dat1[,1],dat2[,1])/1e6));
maxmz<-c(max(dat1[,1],dat2[,1])+(comp[["Parameters"]][2]*max(dat1[,1],dat2[,1])/1e6));
}else{
minmz<-c(min(dat1[,1],dat2[,1])-comp[["Parameters"]][2]);
maxmz<-c(max(dat1[,1],dat2[,1])+comp[["Parameters"]][2]);
}
}else{
if(length(comp[["pattern peak list"]])>1){
mintol<-c(min(dat1[,3])-comp[["Parameters"]][1]);
maxtol<-c(max(dat1[,3])+comp[["Parameters"]][1]);
if(comp[[6]][4]==TRUE){
minmz<-c(min(dat1[,1])-(comp[["Parameters"]][2]*min(dat1[,1])/1e6));
maxmz<-c(max(dat1[,1])+(comp[["Parameters"]][2]*max(dat1[,1])/1e6));
}else{
minmz<-c(min(dat1[,1])-comp[["Parameters"]][2]);
maxmz<-c(max(dat1[,1])+comp[["Parameters"]][2]);
}
}
if(length(comp[[3]])>1){
mintol<-c(min(dat2[,3])-comp[["Parameters"]][1]);
maxtol<-c(max(dat2[,3])+comp[["Parameters"]][1]);
if(comp[[6]][4]==TRUE){
minmz<-c(min(dat2[,1])-(comp[["Parameters"]][2]*min(dat2[,1])/1e6));
maxmz<-c(max(dat2[,1])+(comp[["Parameters"]][2]*max(dat2[,1])/1e6));
}else{
minmz<-c(min(dat2[,1])-comp[["Parameters"]][2]);
maxmz<-c(max(dat2[,1])+comp[["Parameters"]][2]);
}
}
}
if(length(comp[["pattern peak list"]])>1){
dat4<-comp[["pattern peak list"]][
comp[["pattern peak list"]][,3]>=mintol &
comp[["pattern peak list"]][,3]<=maxtol &
comp[["pattern peak list"]][,1]>=minmz &
comp[["pattern peak list"]][,1]<=maxmz
,]
}else{
dat4<-comp[["adduct peak list"]][
comp[["adduct peak list"]][,3]>=mintol &
comp[["adduct peak list"]][,3]<=maxtol &
comp[["adduct peak list"]][,1]>=minmz &
comp[["adduct peak list"]][,1]<=maxmz
,]
}
if(n_col[1]==FALSE){
colo<-"darkgreen"
}else{
colo<-rainbow(n_col[1],s = .5, v = .7)
set.seed(1)
colo<-sample(colo,size=n_col[1], replace = FALSE)
}
if(masslim[1]==FALSE){
use_xlim<-c(min(dat4[,1]),max(dat4[,1]))
}else{
use_xlim<-masslim
}
plot(dat4[,1],dat4[,2],type="h",xlab="",ylab="",lwd=1,col="lightgrey",cex.axis=.9,xlim=use_xlim);
title(xlab="m/z",ylab="Intensity",cex.lab=.9,line=2.6)
if(length(comp[["pattern peak list"]])>2){
points(dat1[,1],dat1[,2],type="h",lwd=2,col=colo[as.numeric(dat1[,"peak ID"])]);
}
if(length(comp[["adduct peak list"]])>2){
points(dat2[,1],dat2[,2],type="h",lwd=2,col=colo[as.numeric(dat2[,"peak ID"])]);
}
if(masslim[1]!=FALSE){
mtext("Zoomed in - click to zoom out partly or double-click to zoom out fully.", side = 3, line=0.1, cex=.8, col="darkgrey", at=use_xlim[1], adj = 0)
}else{
mtext("Brush and doubleclick to zoom.", side = 3, line=0.1, cex=.8, col="darkgrey", at=use_xlim[1], adj = 0)
}
}
if(what=="table"){
if(length(comp[["pattern peak list"]])>1 & length(comp[["adduct peak list"]])>1){
mintol<-c(min(dat1[,3],dat2[,3])-comp[["Parameters"]][1]);
maxtol<-c(max(dat1[,3],dat2[,3])+comp[["Parameters"]][1]);
if(comp[["Parameters"]][4]==TRUE){
minmz<-c(min(dat1[,1],dat2[,1])-(comp[["Parameters"]][2]*min(dat1[,1],dat2[,1])/1e6));
maxmz<-c(max(dat1[,1],dat2[,1])+(comp[["Parameters"]][2]*max(dat1[,1],dat2[,1])/1e6));
}else{
minmz<-c(min(dat1[,1],dat2[,1])-comp[["Parameters"]][2]);
maxmz<-c(max(dat1[,1],dat2[,1])+comp[["Parameters"]][2]);
}
}else{
if(length(comp[["pattern peak list"]])>1){
mintol<-c(min(dat1[,3])-comp[["Parameters"]][1]);
maxtol<-c(max(dat1[,3])+comp[["Parameters"]][1]);
if(comp[[6]][4]==TRUE){
minmz<-c(min(dat1[,1])-(comp[["Parameters"]][2]*min(dat1[,1])/1e6));
maxmz<-c(max(dat1[,1])+(comp[["Parameters"]][2]*max(dat1[,1])/1e6));
}else{
minmz<-c(min(dat1[,1])-comp[["Parameters"]][2]);
maxmz<-c(max(dat1[,1])+comp[["Parameters"]][2]);
}
}
if(length(comp[["adduct peak list"]])>1){
mintol<-c(min(dat2[,3])-comp[["Parameters"]][1]);
maxtol<-c(max(dat2[,3])+comp[["Parameters"]][1]);
if(comp[[6]][4]==TRUE){
minmz<-c(min(dat2[,1])-(comp[["Parameters"]][2]*min(dat2[,1])/1e6));
maxmz<-c(max(dat2[,1])+(comp[["Parameters"]][2]*max(dat2[,1])/1e6));
}else{
minmz<-c(min(dat2[,1])-comp[["Parameters"]][2]);
maxmz<-c(max(dat2[,1])+comp[["Parameters"]][2]);
}
}
}
if(length(comp[["pattern peak list"]])>1){
dat4<-comp[["pattern peak list"]][
comp[["pattern peak list"]][,3]>=mintol &
comp[["pattern peak list"]][,3]<=maxtol &
comp[["pattern peak list"]][,1]>=minmz &
comp[["pattern peak list"]][,1]<=maxmz
,]
}else{
dat4<-comp[["adduct peak list"]][
comp[["adduct peak list"]][,3]>=mintol &
comp[["adduct peak list"]][,3]<=maxtol &
comp[["adduct peak list"]][,1]>=minmz &
comp[["adduct peak list"]][,1]<=maxmz
,]
}
these1<-c();
these2<-c();
these3<-c();
for(i in 1:length(get3)){
for(j in 1:length(get3)){
if(length(comp[["pattern peak list"]])>1){
if(relat1[i,j]!=""){
these1<-c(these1,paste(get3[i],"-",get3[j],sep=""));
these2<-c(these2,substr(relat1[i,j],2,nchar(relat1[i,j])));
these3<-c(these3,(dat1[dat1[,4]==get3[i],2]/dat1[dat1[,4]==get3[j],2]));
}
}
}
}
for(i in 1:length(get3)){
for(j in 1:length(get3)){
if(length(comp[["adduct peak list"]])>1){
if(relat2[i,j]!=""){
these1<-c(these1,paste(get3[i],"-",get3[j],sep=""));
these2<-c(these2,substr(relat2[i,j],2,nchar(relat2[i,j])));
these3<-c(these3,(dat2[dat2[,4]==get3[i],2]/dat2[dat2[,4]==get3[j],2]));
}
}
}
}
if(length(these1)){
relat1<-data.frame(these1,these2,these3,stringsAsFactors=FALSE);
names(relat1)<-c("peaks","relation","intensity ratio");
}else{
relat1<-NULL
}
if(comp[[1]][compoID,6]!="-"){
if(length(dat1)>1){
relat<-list(dat1[order(dat1[,1],decreasing=FALSE),c(4,1,2,3,6)],dat4[order(dat4[,1],decreasing=FALSE),
c(4,1,2,3)],relat1,as.character(comp[[1]][compoID,6]))
}else{
relat<-list(dat2[order(dat2[,1],decreasing=FALSE),c(4,1,2,3,6)],dat4[order(dat4[,1],decreasing=FALSE),
c(4,1,2,3)],relat1,as.character(comp[[1]][compoID,6]))
}
}else{
if(length(dat1)>1){
relat<-list(dat1[order(dat1[,1],decreasing=FALSE),c(4,1,2,3,6)],dat4[order(dat4[,1],decreasing=FALSE),
c(4,1,2,3)],relat1,"Not part of a homologue series")
}else{
relat<-list(dat2[order(dat2[,1],decreasing=FALSE),c(4,1,2,3,6)],dat4[order(dat4[,1],decreasing=FALSE),
c(4,1,2,3)],relat1,"Not part of a homologue series")
}
}
names(relat)<-c("a","all peaks within range","relations","Part of homologue series:");
relat$a<-cbind(relat$a,get3_blind[match(relat$a[,"peak ID"],get3)])
names(relat$a)[6]<-"in blind?"
return(relat);
}
}
