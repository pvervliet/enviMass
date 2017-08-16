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
  Intlim=FALSE,
  masslim=FALSE,
  normalize=FALSE,
  n_col=FALSE,
  set_RT="seconds",
  chromat_full=FALSE
){

    ############################################################################
    normalize<-as.logical(normalize)
    these_peaks<-peakIDs
    #if(any(is.na(these_peaks))){stop("\n Invalid peakIDs in plotchromat detected!")}
    if(RTlim[1]==FALSE){
      min_RT<-min(MSlist[["Scans"]][[2]][,"RT"])
      max_RT<-max(MSlist[["Scans"]][[2]][,"RT"])
      if(set_RT=="minutes"){
        min_RT<-(min_RT/60)
        max_RT<-(max_RT/60)
      }
      x_lim<-c(min_RT,max_RT)
    }else{
      x_lim<-RTlim
    }
    if(normalize){
      y_lim<-c(0,1)
    }else{
      if(Intlim[1]==FALSE){  
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
      }else{
        y_lim<-Intlim
      }
    }
    all_RTs<-unique(MSlist[["Scans"]][[2]][,"RT"])
    all_RTs<-all_RTs[order(all_RTs,decreasing=FALSE)]
    if(set_RT=="minutes"){all_RTs<-(all_RTs/60)}
    ############################################################################ 

    ############################################################################    
    if(n_col[1]==FALSE){
      colo<-rainbow(length(these_peaks),s = .5, v = .7)
    }else{
      colo<-rainbow(n_col[1],s = .5, v = .7)
      set.seed(1)
      colo<-sample(colo,size=n_col[1], replace = FALSE)
      colo<-colo[these_peaks]
    }
    plot.new()
    plot.window(xlim=x_lim,ylim=y_lim)
    for(i in 1:length(these_peaks)){    
        RTs<-MSlist[["Scans"]][[2]][
            MSlist[["Peak_index"]][these_peaks[i],"start_ID"]:MSlist[["Peak_index"]][these_peaks[i],"end_ID"]
          ,"RT"]
        if(set_RT=="minutes"){RTs<-(RTs/60)}
        intens<-MSlist[["Scans"]][[2]][
            MSlist[["Peak_index"]][these_peaks[i],"start_ID"]:MSlist[["Peak_index"]][these_peaks[i],"end_ID"]
          ,"intensity"]
        all_int<-rep(0,length(all_RTs))
        all_int[match(RTs,all_RTs)]<-intens
        if(chromat_full){  
          #####################################################################
          this_EIC<-MSlist[["Peaklist"]][match(these_peaks[i],MSlist[["Peaklist"]][,"peak_ID"]),"EIC_ID"][[1]]
          RTsb<-MSlist[["Scans"]][[2]][
              MSlist[["EIC_index"]][this_EIC,"start_ID"]:MSlist[["EIC_index"]][this_EIC,"end_ID"]
            ,"RT"]          
          if(set_RT=="minutes"){RTsb<-(RTsb/60)}
          intensb<-MSlist[["Scans"]][[2]][
              MSlist[["EIC_index"]][this_EIC,"start_ID"]:MSlist[["EIC_index"]][this_EIC,"end_ID"]
            ,"intensity"]
          if(normalize){
            intensc<-(intensb/max(intensb))
          }else{
            intensc<-intensb
          }
          all_intb<-rep(0,length(all_RTs))
          all_intb[match(RTsb,all_RTs)]<-intensc
          points(all_RTs,all_intb,type="l",col="darkgrey",lwd=.5,lty="dashed")       
          #####################################################################
          if(normalize){
            all_int<-(all_int/max(intensb))
          }
          points(all_RTs[all_int>0],all_int[all_int>0],type="l",col=colo[i],lwd=2)
        }else{
          if(normalize){
            all_int<-(all_int/max(all_int))
          }
          points(all_RTs,all_int,type="l",col=colo[i],lwd=1.5)
        }
    }
    if(!normalize & set_RT=="seconds"){title(xlab="RT [s]",ylab="Intensity",cex.lab=.9,line=3)}
    if(!normalize & set_RT=="minutes"){title(xlab="RT [min]",ylab="Intensity",cex.lab=.9,line=3)}
    if(normalize & set_RT=="seconds"){title(xlab="RT [s]",ylab="Norm. intens.",cex.lab=.9,line=3)}
    if(normalize & set_RT=="minutes"){title(xlab="RT [min]",ylab="Norm. intens.",cex.lab=.9,line=3)}
    box();axis(1,cex.axis=.9);axis(2,cex.axis=.9)
    if(RTlim[1]!=FALSE | Intlim[1]!=FALSE){
      mtext("Zoomed in - click to zoom out partly or double-click to zoom out fully.", side = 3, line=0.1, cex=.8, col="darkgrey", at=x_lim[1], adj = 0)
    }else{
      mtext("Brush and doubleclick to zoom.", side = 3, line=0.1, cex=.8, col="darkgrey", at=x_lim[1], adj = 0)
    }
    ############################################################################    
    return("done")
    ############################################################################

}
