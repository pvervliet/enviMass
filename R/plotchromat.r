#' @title Plot chromatograms for peak(s)
#'
#' @export
#'
#' @description \code{plotprofiles} plots EICs for peak(s)
#'
#' @details enviMass workflow function
#' 

plotchromat <-
function(
  MSlist,
  peakIDs,
  RTlim=FALSE,
  normalize=FALSE
	){

    ############################################################################
    these_peaks<-peakIDs
    #if(any(is.na(these_peaks))){stop("\n Invalid peakIDs in plotchromat detected!")}
    if(RTlim[1]==FALSE){
      min_RT<-min(MSlist[["Scans"]][[2]][,"RT"])
      max_RT<-max(MSlist[["Scans"]][[2]][,"RT"])
      x_lim<-c(min_RT,max_RT)
    }else{
      x_lim<-RTlim
    }
    if(normalize){
      y_lim<-c(0,1)
    }else{
      y_lim<-c(0,0)
      for(i in these_peaks){
        if(
          max(MSlist[["Scans"]][[2]][
            MSlist[["Peak_index" ]][i,"start_ID"]:MSlist[["Peak_index"]][i,"end_ID"]
          ,"intensity"])
          > y_lim[2]
        ){
          y_lim[2]<-max(MSlist[["Scans"]][[2]][
            MSlist[["Peak_index"]][i,"start_ID"]:MSlist[["Peak_index"]][i,"end_ID"]
          ,"intensity"])
        }
      }
    }
    all_RTs<-unique(MSlist[["Scans"]][[2]][,"RT"])
    all_RTs<-all_RTs[order(all_RTs,decreasing=FALSE)]
    ############################################################################ 

    ############################################################################    
    colo<-rainbow(length(these_peaks),s = .5, v = .7)
    plot.new()
    plot.window(xlim=x_lim,ylim=y_lim)
    for(i in 1:length(these_peaks)){    
        RTs<-MSlist[["Scans"]][[2]][
            MSlist[["Peak_index"]][these_peaks[i],"start_ID"]:MSlist[["Peak_index"]][these_peaks[i],"end_ID"]
          ,"RT"]
        intens<-MSlist[["Scans"]][[2]][
            MSlist[["Peak_index"]][these_peaks[i],"start_ID"]:MSlist[["Peak_index"]][these_peaks[i],"end_ID"]
          ,"intensity"]
        if(normalize){
          intens<-(intens/max(intens))
        }
        all_int<-rep(0,length(all_RTs))
        all_int[match(RTs,all_RTs)]<-intens
        points(all_RTs,all_int,type="l",col=colo[i])
    }
    if(!normalize){
      title(xlab="RT [s]",ylab="Intensity")
    }else{
      title(xlab="RT [s]",ylab="Norm. intensity")
    }
    box();axis(1);axis(2)
    ############################################################################    
    return("done")
    ############################################################################

}
