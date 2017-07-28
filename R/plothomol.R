plothomol<-
function(
	homol,
	xlim=FALSE,
	ylim=FALSE,
	plotlegend=FALSE,
	plotdefect=FALSE,
	omit_theta=FALSE,
	deldel=.001
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
  that<-round(homol[["homol_peaks_relat"]][,3],digits=2);
  this<-unique(that)
  colo<-rainbow(length(this),s = .5, v = .9)
  coloring<-colo[match(that,this)]
  if(omit_theta[1]!=FALSE){
    use_segm<-(homol[["homol_peaks_relat"]][,5]>=omit_theta)
  }else{
    use_segm<-rep(TRUE,dim(homol[["homol_peaks_relat"]])[1])
  }
  ############################################################################
    
  ############################################################################
  # plot #####################################################################
  if(plotlegend==FALSE){
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
    segments(
      x0=homol[["Peaks in homologue series"]][homol[["homol_peaks_relat"]][use_segm,1],"mz"],
      y0=homol[["Peaks in homologue series"]][homol[["homol_peaks_relat"]][use_segm,1],"RT"],
      x1=homol[["Peaks in homologue series"]][homol[["homol_peaks_relat"]][use_segm,2],"mz"],
      y1=homol[["Peaks in homologue series"]][homol[["homol_peaks_relat"]][use_segm,2],"RT"],
      col=coloring[use_segm]
    )
}
############################################################################

############################################################################
# add a legend?
if(plotlegend==TRUE){
  ##########################################################################
	counts<-.Call("moving_count", 
		homol[["homol_peaks_relat"]], 
		deldel,
		PACKAGE="enviMass"
	)
  plot(homol[["homol_peaks_relat"]][,3],counts,type="h",col=coloring,lwd=1.5,xlab="m/z difference",ylab="Moving window count",xlim=c(min(homol[["homol_peaks_relat"]][,3]),max(homol[["homol_peaks_relat"]][,3])))
  ##########################################################################


}
############################################################################
    
}
