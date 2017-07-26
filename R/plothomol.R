plothomol<-
function(
	homol,
	xlim=FALSE,
	ylim=FALSE,
	plotlegend=TRUE,
	plotdefect=FALSE
){

  ############################################################################
  # check inputs #############################################################
  if(xlim[1]!=FALSE){if(length(xlim)>2){stop("xlim not correct!")}}
  if(ylim[1]!=FALSE){if(length(ylim)>2){stop("xlim not correct!")}}
  if(length(homol[[5]])<1){stop("no homologue series found!")}
	if(!is.logical(plotlegend)){stop("plotlegend must be TRUE or FALSE")}
	if(!is.logical(plotdefect)){stop("plotdefect must be TRUE or FALSE")}	
	if(plotdefect){
    mass_def<-c(homol[[1]][,1]-round(homol[[1]][,1]))		# calculate mass defect
		homol[[1]][,3]<-mass_def # local function definition
	}
  ############################################################################
    
  ############################################################################
  # plot #####################################################################
  sc<-close.screen();if(sc[1]!=FALSE){for(m in 1:length(sc)){close.screen(sc[m])}};
  plot.new();
  if(xlim[1]!=FALSE & ylim[1]==FALSE){plot.window(ylim=ylim,xlim=c(min(homol[[1]][,1]),max(homol[[1]][,1])));}
  if(ylim[1]!=FALSE & xlim[1]==FALSE){
    if(plotlegend==TRUE){
      plot.window(ylim=c(min(homol[[1]][,3]),max(homol[[1]][,3])*1.2),xlim=xlim);
    }else{
      plot.window(ylim=c(min(homol[[1]][,3]),max(homol[[1]][,3])),xlim=xlim);
    }
  }
  if(xlim[1]==FALSE & ylim[1]==FALSE){
    if(plotlegend==TRUE){
      plot.window(ylim=c(min(homol[[1]][,3]),max(homol[[1]][,3])*1.2),xlim=c(min(homol[[1]][,1]),max(homol[[1]][,1])));
    }else{
      plot.window(ylim=c(min(homol[[1]][,3]),max(homol[[1]][,3])),xlim=c(min(homol[[1]][,1]),max(homol[[1]][,1])));    
    }
  }
  if(xlim[1]!=FALSE & ylim[1]!=FALSE){plot.window(xlim=xlim,ylim=ylim);}
  box();axis(1);axis(2);
	if(!plotdefect){
		title(ylab="Retention time [s]",xlab="m/z");
  }else{
		title(ylab="mass defect",xlab="m/z");	
	}
	points(homol[[1]][,1],homol[[1]][,3],cex=0.3,pch=19,col="lightgrey");
  ############################################################################
    
  ############################################################################
  that<-round(homol[["homol_peaks_relat"]][,3],digits=1);
  this<-unique(that)
  colo<-rainbow(length(this),s = .5, v = .9)
  coloring<-colo[match(that,this)]
  segments(
    x0=homol[["Peaks in homologue series"]][homol[["homol_peaks_relat"]][,1],"mz"],
    y0=homol[["Peaks in homologue series"]][homol[["homol_peaks_relat"]][,1],"RT"],
    x1=homol[["Peaks in homologue series"]][homol[["homol_peaks_relat"]][,2],"mz"],
    y1=homol[["Peaks in homologue series"]][homol[["homol_peaks_relat"]][,2],"RT"],
    col=coloring
  )
  ############################################################################

  ############################################################################
  # add a legend
  if(plotlegend==TRUE){
    plot.window(xlim=c(0,1),ylim=c(min(as.numeric(that)),max(as.numeric(that))));
    lines(c(0.95,0.95),c(min(as.numeric(that)),max(as.numeric(that))),col="lightgrey",lwd=6)
    it<-2
    for(i in 1:length(that)){
      points(0.95,as.numeric(that[i]),pch=19,col=coloring[i])
      text(0.95,as.numeric(that[i]),labels=that[i],col=coloring[i],cex=0.65,pos=it)
      if(it==2){it<-4}else{it<-2}
    }
  }
  ############################################################################
    
}
