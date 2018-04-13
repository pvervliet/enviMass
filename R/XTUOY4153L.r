XTUOY4153L <- function(
homol,
xlim=FALSE,
ylim=FALSE,
dmasslim = FALSE,
dRTlim=FALSE,
plot_what="mz_RT",
plotdefect=FALSE,
omit_theta=FALSE,
deldel=.001,
plot_those=FALSE,
emph_point=FALSE,
emph_series=FALSE,
with_num=TRUE
){
if(xlim[1]!=FALSE){if(length(xlim)>2){stop("xlim not correct!")}}
if(ylim[1]!=FALSE){if(length(ylim)>2){stop("xlim not correct!")}}
if(length(homol[[5]])<1){stop("no homologue series found!")}
if(!is.logical(plotdefect)){stop("plotdefect must be TRUE or FALSE")}
if(plotdefect){
mass_def<-c(homol[[1]][,1]-round(homol[[1]][,1]))
homol[[1]][,3]<-mass_def
}
if(length(plot_those)>1){
if(length(plot_those) != dim(homol[["homol_peaks_relat"]])[1]){stop("plot_those must match homol")}
}else{
plot_those<-rep(TRUE,dim(homol[["homol_peaks_relat"]])[1])
}
that<-round(homol[["homol_peaks_relat"]][,3],digits=2);
this<-unique(that)
colo<-rainbow(length(this),s = .4, v = .8)
coloring<-colo[match(that,this)]
if(omit_theta[1]!=FALSE){
use_segm<-((homol[["homol_peaks_relat"]][,5]>=omit_theta) & plot_those)
}else{
use_segm<-plot_those
}
if(plot_what=="mz_RT"){
if(xlim[1]!=FALSE){use_xlim<-xlim}else{use_xlim<-c(min(homol[["Peaks in homologue series"]][,1]),max(homol[["Peaks in homologue series"]][,1]))}
if(ylim[1]!=FALSE){use_ylim<-ylim}else{use_ylim<-c(min(homol[["Peaks in homologue series"]][,3]),max(homol[["Peaks in homologue series"]][,3]))}
plot.new();
plot.window(xlim=use_xlim,ylim=use_ylim);
axis(1);axis(2);
if(!plotdefect){
title(ylab="Retention time [s]",xlab="m/z");
}else{
title(ylab="mass defect",xlab="m/z");
}
if(!plotdefect){
points(
homol[["Peaks in homologue series"]][homol[["Peaks in homologue series"]][,"HS IDs"]=="0","mz"],
homol[["Peaks in homologue series"]][homol[["Peaks in homologue series"]][,"HS IDs"]=="0","RT"],
cex=0.2,pch=19,col="lightgrey");
points(
homol[["Peaks in homologue series"]][homol[["Peaks in homologue series"]][,"HS IDs"]!="0","mz"],
homol[["Peaks in homologue series"]][homol[["Peaks in homologue series"]][,"HS IDs"]!="0","RT"],
cex=0.5,pch=19,col="darkgrey");
segments(
x0=homol[["Peaks in homologue series"]][homol[["homol_peaks_relat"]][use_segm,1],"mz"],
y0=homol[["Peaks in homologue series"]][homol[["homol_peaks_relat"]][use_segm,1],"RT"],
x1=homol[["Peaks in homologue series"]][homol[["homol_peaks_relat"]][use_segm,2],"mz"],
y1=homol[["Peaks in homologue series"]][homol[["homol_peaks_relat"]][use_segm,2],"RT"],
col=coloring[use_segm]
)
if(emph_series[1]!=FALSE){
these_series<-match(emph_series,homol[["Homologue Series"]][,"HS IDs"])
for(i in these_series){
these_peaks<-as.numeric(strsplit(homol[["Homologue Series"]][i,"peak IDs"],",")[[1]])
these_peaks<-match(these_peaks,homol[["Peaks in homologue series"]][,"peak ID"])
for(j in 2:length(these_peaks)){
segments(
x0=homol[["Peaks in homologue series"]][these_peaks[j],"mz"],
y0=homol[["Peaks in homologue series"]][these_peaks[j],"RT"],
x1=homol[["Peaks in homologue series"]][these_peaks[j-1],"mz"],
y1=homol[["Peaks in homologue series"]][these_peaks[j-1],"RT"],
col="black",lwd=1
)
}
for(j in 1:length(these_peaks)){
points(
homol[["Peaks in homologue series"]][these_peaks[j],"mz"],homol[["Peaks in homologue series"]][these_peaks[j],"RT"],
cex=0.5,pch=19,col="black");
}
}
}
if(emph_point[1]!=FALSE){
that_peak<-match(emph_point[1],homol[["Peaks in homologue series"]][,"peak ID"])
points(
homol[["Peaks in homologue series"]][that_peak,"mz"],homol[["Peaks in homologue series"]][that_peak,"RT"],
cex=0.5,pch=19,col="black");
points(
homol[["Peaks in homologue series"]][that_peak,"mz"],homol[["Peaks in homologue series"]][that_peak,"RT"],
cex=2,pch=1,col="black");
}
if(xlim[1]!=FALSE | ylim[1]!=FALSE){
mtext("Zoomed in - click to zoom out partly or double-click to zoom out fully.", side = 3, line=0.1, cex=.8, col="darkgrey", at=use_xlim[1], adj = 0)
}else{
mtext("Brush and doubleclick to zoom.", side = 3, line=0.1, cex=.8, col="darkgrey", at=use_xlim[1], adj = 0)
}
if(any(!plot_those)){
mtext("Series subsets selected.", side = 3, line=0.1, cex=.8, col="darkgrey", at=use_xlim[2], adj = 1)
}
}
box();
if(with_num){
plot.window(xlim=c(0,1),ylim=c(0,1))
text(0,1,labels="C",cex=1.3)
plot.window(xlim=use_xlim,ylim=use_ylim);
}
}
if( plot_what=="dmass"){
if(dmasslim[1]!=FALSE){
use_dmasslim<-dmasslim
}else{
use_dmasslim<-c(min(homol[["homol_peaks_relat"]][,3]),max(homol[["homol_peaks_relat"]][,3]))
}
if(any(plot_those)){
counts<-.Call("_enviMass_moving_count",
homol[["homol_peaks_relat"]][plot_those,,drop=FALSE],
deldel,
PACKAGE="enviMass"
)
plot(
homol[["homol_peaks_relat"]][plot_those,3],counts,
type="h",col=coloring[plot_those],lwd=1.5,
xlab="Series m/z difference",ylab="Moving window count",
xlim=use_dmasslim,ylim=c(0,max(counts)))
if(dmasslim[1]!=FALSE){
mtext("Series mass differences zoomed in - click to zoom out partly or double-click to zoom out fully.", side = 3, line=0.1, cex=.8, col="darkgrey", at=use_dmasslim[1], adj = 0)
}else{
mtext("Brush to select subrange. Doubleclick into brush to zoom in.", side = 3, line=0.1, cex=.8, col="darkgrey", at=use_dmasslim[1], adj = 0)
}
}else{
plot.new();
plot.window(xlim=use_dmasslim,ylim=c(0,1));
title(xlab="Series m/z difference");
axis(1);
box();
}
if(with_num){
plot.window(xlim=c(0,1),ylim=c(0,1))
text(0,.85,labels="A",cex=1.3)
plot.window(xlim=use_dmasslim,ylim=c(0,max(counts)))
}
}
if( plot_what=="dRT"){
if(dRTlim[1]!=FALSE){
use_dRTlim<-dRTlim
}else{
use_dRTlim<-c(min(homol[["homol_peaks_relat"]][,4]),max(homol[["homol_peaks_relat"]][,4]))
}
if(dmasslim[1]!=FALSE){
use_dmasslim<-dmasslim
}else{
use_dmasslim<-c(min(homol[["homol_peaks_relat"]][,3]),max(homol[["homol_peaks_relat"]][,3]))
}
if(any(plot_those)){
plot(
homol[["homol_peaks_relat"]][plot_those,4],homol[["homol_peaks_relat"]][plot_those,3],
col=coloring[plot_those],pch=15,cex=.8,
xlab="Series RT difference [s]",ylab="m/z difference",
xlim=use_dRTlim,ylim=use_dmasslim)
if(dRTlim[1]!=FALSE){
mtext("Series RT differences zoomed in - click to zoom out partly or double-click to zoom out fully.", side = 3, line=0.1, cex=.8, col="darkgrey", at=use_dRTlim[1], adj = 0)
}else{
mtext("Brush to select subrange. Doubleclick into brush to zoom in.", side = 3, line=0.1, cex=.8, col="darkgrey", at=use_dRTlim[1], adj = 0)
}
}else{
plot.new();
plot.window(xlim=use_dRTlim,ylim=c(0,1));
title(xlab="Series RT difference");
axis(1);
box();
}
if(with_num){
plot.window(xlim=c(0,1),ylim=c(0,1))
text(0,.85,labels="B",cex=1.3)
plot.window(xlim=use_dRTlim,ylim=use_dmasslim)
}
}
}
