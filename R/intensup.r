#' @title Trend detection and blind subtraction for time profiles.
#'
#' @export
#'
#' @description \code{intensup} runs a trend detection and blind subtraction for the list of time profiles.
#'
#' @param profileList A profile list.
#' @param from Logical or integer of index.
#' @param to Logical or integer of index.
#' @param progbar Logical. Should a progress bar be shown? Only for Windows.
#' @param blindsub Logical. Run blind subtraction?
#' @param blindfold Numerical. Blind definition; above blind, if \code{blindfold} larger in intensity.
#' @param lags Vector of numericals.
#' @param threshold Numerical. A trend is reported if its intensity is \code{threshold} above the mean intensity plus the intensity deviation of other trends.
#' @param notrend Logical. Report global trend intensity as maximum intensity after blind subtraction.
#' @param omit_trend Logical. Omit trend detection altogether? Used in dont_trendblind.r
#'
#' @return Updated \code{profileList[[7]]}.
#' 
#' @details enviMass workflow function for trend detection and blind subtraction.
#' 

intensup<-function(
	profileList,
	from=FALSE,
	to=FALSE,
	progbar=FALSE,
	blindsub=TRUE, # do a blind subtraction?
	blindfold=100, # how much higher in intensity than blind?
	lags=c(5,14),  # time lags
	threshold=3,   # trend threshold: 
	notrend=FALSE,  # no global threshold, but global maximum above blind	
	omit_trend=FALSE
){

    ############################################################################
    if(!profileList[[ "state"]][[3]]){stop("profileList not profiled; aborted.")}
    if(!from){m=1}else{m=from}
    if(!to){n=length(profileList[["index_prof"]][,1])}else{n=to}
	if(blindsub!=FALSE){if(!is.numeric(blindfold) || (blindfold<0)){stop("Invalid blindfold argument; aborted.")}}
    if(blindsub!=FALSE){subit=1;subrat=blindfold;}else{subit=2;subrat=0;}
	if(!is.numeric(lags)){stop("lags argument must be numeric; aborted.")}
	if(!is.logical(notrend)){stop("notrend must be logical.")}
	if(!is.logical(omit_trend)){stop("notrend must be logical.")}
	############################################################################
    # set matrix to sort & store data from a profile ###########################
    atPOSIX<-profileList[["datetime"]];
    sampletype<-profileList[["type"]];
    sampleID<-profileList[["sampleID"]];
	# filter out other file types such as spiked ones
	keep<-((sampletype=="sample")|(sampletype=="blank"))
	atPOSIX<-atPOSIX[keep]
	sampletype<-sampletype[keep]
	sampleID<-sampleID[keep]
	#
    atdate<-c();
    attime<-c();
    for(i in 1:length(atPOSIX)){
        atdate<-c(atdate, strsplit(atPOSIX[i]," ")[[1]][1]);
        attime<-c(attime, strsplit(atPOSIX[i]," ")[[1]][2]);
    }
    attime<-as.difftime(attime);
    atdate<-as.Date(atdate, tz="GMT");
    ord<-order(as.numeric(atdate),as.numeric(attime),sampleID);
    atPOSIXsort<-atPOSIX[ord];
    atdate<-atdate[ord];
    attime<-attime[ord];
	sampleID<-sampleID[ord];
	sampletype<-sampletype[ord];
    timeset<-matrix(nrow=length(atPOSIX),ncol=(5+(4*length(lags))),0);
    for(i in 1:length(sampleID)){
      if(sampletype[i]=="sample"){
        timeset[i,2]<-as.numeric(sampleID[i]);
      }
      if(sampletype[i]=="blank"){
        timeset[i,3]<-as.numeric(sampleID[i]);
      }
    }
    numtime<-(as.numeric(atdate)+as.numeric(attime/(24*60*60)))
	colnames(timeset)<-c("above blank?","sampleID","blankID","sample_int","blank_int",rep("lag_int",length(lags)),rep("del_int",length(lags)),rep("max_time",length(lags)),rep("blind_int",length(lags)))			
	leng<-max(seq(1,length(timeset[,1]),1)[timeset[,2]!=0])	
	latestID<-timeset[leng,2][[1]]
	############################################################################
	# check & adjust lags ######################################################
	if(any(lags>(max(numtime)-min(numtime)+1))){
		lags<-lags[lags<=(max(numtime)-min(numtime)+1)]
		cat("WARNING: at least one lag longer than covered time period - omitted!\n")
		if(length(lags)==0){
			stop("...no lags left; aborted.")
		}		
	}
	# INSERT ... also check any of the inter-sample distances ##################
	if(max(diff(numtime[timeset[,2]!=0],lag=1))>max(lags)){
		warning("\n At least one gap in sample time series larger than largest lag!")
	}
	############################################################################
    profileList[["index_prof"]][,5:7]<-0;
    if(progbar==TRUE){  prog<-winProgressBar("Extract intensity differences...",min=m,max=n);
						setWinProgressBar(prog, 0, title = "Extract intensity differences...", label = NULL);}
    for(k in m:n){
      if(progbar==TRUE){setWinProgressBar(prog, k, title = "Extract intensity differences...", label = NULL)}
      if(profileList[["index_prof"]][k,"number_peaks_total"]>1){
        ########################################################################
        # fill timeset #########################################################
        timeset[,4:length(timeset[1,])]<-0;
        timeset[,c(4,5)] <-.Call("fill_timeset",
                                as.numeric(timeset),
                                as.numeric(profileList[["peaks"]][(profileList[["index_prof"]][k,"start_ID"]:profileList[[7]][k,"end_ID"]),"sampleIDs"]), 
                                as.numeric(profileList[["peaks"]][(profileList[["index_prof"]][k,"start_ID"]:profileList[[7]][k,"end_ID"]),"intensity"]), 
                                as.integer(length(timeset[,1])),
                                PACKAGE="enviMass"
                            )	
		########################################################################
        # interpolate & subtract blind,  #######################################
        ########################################################################
        # subtract blank #######################################################
        ########################################################################
        if(any(timeset[,4]>0)){ # any non-blind peak present?
			what<-1 # !=1 -> get raw output, i.e., peak series
			if(!omit_trend){
				that<-.Call("meandel",
					as.numeric(timeset),
					as.integer(subit),
					as.numeric(subrat),
					as.numeric(numtime),
					as.integer(what),
					as.numeric(lags),
					as.numeric(threshold),
					as.integer(notrend),
					PACKAGE="enviMass"
				)
			}
			if( (what!=1) & (!omit_trend)){ 
				that<-matrix(nrow=length(atPOSIX),ncol=(5+(4*length(lags))),that);
				colnames(that)<-c("above blank?","sampleID","blankID","sample_int","blank_int",rep("lag_int",length(lags)),rep("del_int",length(lags)),rep("max_time",length(lags)),rep("blind_int",length(lags)))			
				# plot smoothed series ... and abort ###########################
				plot.new();
				plot.window(xlim=c(min(numtime),max(numtime)),ylim=c(min(that[that[,2]!=0,4]),max(that[that[,2]!=0,4])));
				box();axis(1);axis(2);
				title(xlab="Time",ylab="Intensity")
				points(numtime[that[,2]!=0],that[that[,2]!=0,4],col="red",type="l");
				for(i in 1:length(lags)){
					points(numtime[that[,2]!=0],that[that[,2]!=0,(5+i)],col="darkgrey",type="l");
					points(that[that[,(5+i+(length(lags)*2))]!=0,(5+i+(length(lags)*2))],that[that[,(5+i+(length(lags)*2))]!=0,(5+i+length(lags))],col="blue",pch=19);					
				}
				stop(" YOU wanted the smoothed series...\n")
			}else{		
				if(!omit_trend){
					profileList[["index_prof"]][k,"deltaint_newest"]<-max(that[4,]); # current>abs.dev
					profileList[["index_prof"]][k,"deltaint_global"]<-max(that[5,]); # global>abs.dev
					profileList[["index_prof"]][k,"absolute_mean_dev"]<-max(that[3,]); # abs.dev
				}
				if(any(timeset[,5]>0)){ 			  # in blind?
					profileList[["index_prof"]][k,"blind?"]<-1 # in blind
					profileList[["index_prof"]][k,"number_peaks_blind"]<-length(timeset[timeset[,5]!=0,5]) 	# number_peaks_blind
					profileList[["index_prof"]][k,"mean_int_blind"]<-mean(timeset[timeset[,5]!=0,5]) 	# mean_int_blind
					profileList[["index_prof"]][k,"max_int_blind"]<-max(timeset[timeset[,5]!=0,5]) 	# mean_int_blind					
				}else{
					profileList[["index_prof"]][k,"blind?"]<-0 # in blind
					profileList[["index_prof"]][k,"number_peaks_blind"]<-0 # number_peaks_blind
					profileList[["index_prof"]][k,"mean_int_blind"]<-0 # mean_int_blind		
					profileList[["index_prof"]][k,"max_int_blind"]<-0 # mean_int_blind						
				}
				#############################################
				# Replicates: mean sample above mean blind ? <- TO BE DELETED -> now in script do_profblind.r
				#if(any(timeset[,5]>0)){
				#	if( mean(timeset[,4][timeset[,2]!=0])>=(mean(timeset[,5][timeset[,3]!=0])*blindfold) ){
				#		profileList[["index_prof"]][k,"above_blind?"]<-1 
				#	}else{
				#		profileList[["index_prof"]][k,"above_blind?"]<-0 
				#	}
				#}else{
				#	profileList[["index_prof"]][k,"above_blind?"]<-1
				#}				
				#############################################
				profileList[["index_prof"]][k,"number_peaks_sample"]<-length(timeset[timeset[,4]!=0,4]) 	# number_peaks_sample
				profileList[["index_prof"]][k,"mean_int_sample"]<-mean(timeset[timeset[,4]!=0,4]) 	# mean_int_sample	
				profileList[["index_prof"]][k,"max_int_sample"]<-max(timeset[timeset[,4]!=0,4]) 						
			}
		}else{
			# profileList[[7]][k,7]<-() # abs.dev = not of interest for blind
			profileList[["index_prof"]][k,"blind?"]<-1 # only in blind = not above (single peak)
			#profileList[["index_prof"]][k,"above_blind?"]<-0	# mean sample above mean blind ? <- TO BE DELETED -> now in script do_profblind.r
			profileList[["index_prof"]][k,"number_peaks_sample"]<-0 # number_peaks_sample
			profileList[["index_prof"]][k,"number_peaks_blind"]<-length(timeset[timeset[,5]!=0,5])# number_peaks_blind
			profileList[["index_prof"]][k,"mean_int_sample"]<-0 # mean_int_sample
			profileList[["index_prof"]][k,"max_int_sample"]<-0 
			profileList[["index_prof"]][k,"mean_int_blind"]<-mean(timeset[timeset[,5]!=0,5]) # mean_int_blind			
			profileList[["index_prof"]][k,"max_int_blind"]<-max(timeset[timeset[,5]!=0,5])			
		}
		profileList[["index_prof"]][k,"newest_intensity"]<-timeset[leng,4]
		########################################################################
      }else{ # single-peaked profile
		profileList[["index_prof"]][k,"absolute_mean_dev"]<-0;
		if( any( timeset[,3]== (profileList[[2]][profileList[["index_prof"]][k,1],6]) ) ){ # only in blind
			profileList[["index_prof"]][k,"blind?"]<-1 # in blind = not above (single peak)
			#profileList[["index_prof"]][k,"above_blind?"]<-0	# mean sample above mean blind ? <- TO BE DELETED -> now in script do_profblind.r
			profileList[["index_prof"]][k,"number_peaks_sample"]<-0 # number_peaks_sample
			profileList[["index_prof"]][k,"number_peaks_blind"]<-1 # number_peaks_blind
			profileList[["index_prof"]][k,"mean_int_sample"]<-0 # mean_int_sample
			profileList[["index_prof"]][k,"max_int_sample"]<-0 
			profileList[["index_prof"]][k,"mean_int_blind"]<-(profileList[[2]][profileList[[7]][k,1],2]) # mean_int_blind	
			profileList[["index_prof"]][k,"max_int_blind"]<-(profileList[[2]][profileList[[7]][k,1],2])			
		}else{ # only in sample
			profileList[["index_prof"]][k,"blind?"]<-0 
			#profileList[["index_prof"]][k,"above_blind?"]<-1 # not in blind = above (single peak) <- TO BE DELETED -> now in script do_profblind.r	
			profileList[["index_prof"]][k,"deltaint_global"]<-(profileList[[2]][profileList[["index_prof"]][k,1],2])
			if(any(profileList[[2]][profileList[["index_prof"]][k,1],6]==latestID)){
			  profileList[["index_prof"]][k,"deltaint_newest"]<-(profileList[["peaks"]][profileList[["index_prof"]][k,1],2])
			  profileList[["index_prof"]][k,"newest_intensity"]<-(profileList[["peaks"]][profileList[["index_prof"]][k,1],2])
			}
			profileList[["index_prof"]][k,"number_peaks_sample"]<-1 # number_peaks_sample
			profileList[["index_prof"]][k,"number_peaks_blind"]<-0 # number_peaks_blind
			profileList[["index_prof"]][k,"mean_int_sample"]<-(profileList[["peaks"]][profileList[["index_prof"]][k,1],2]) # mean_int_sample
			profileList[["index_prof"]][k,"max_int_sample"]<-(profileList[["peaks"]][profileList[["index_prof"]][k,1],2]) 
			profileList[["index_prof"]][k,"mean_int_blind"]<-0 # mean_int_blind		
			profileList[["index_prof"]][k,"max_int_blind"]<-0				
		}
      }
    }
    if(progbar==TRUE){close(prog);}
	profileList[[1]][[4]]<-TRUE;
    ############################################################################
    return(profileList)

}
