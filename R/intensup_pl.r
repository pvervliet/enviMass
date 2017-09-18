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

# below outcommented prof_index features moved to partlust_pl.r -> required earlier in workflow

intensup_pl<-function(
	peaks,
	index_prof,
	...
){

    ############################################################################
	if(blindsub!=FALSE){subit=1;subrat=blindfold;}else{subit=2;subrat=0;}
	if(blindsub!=FALSE){if(!is.numeric(blindfold) || (blindfold<0)){stop("Invalid blindfold argument; aborted.")}}
	m=1
	n=dim(index_prof)[1]
    ############################################################################
    for(k in m:n){
      if(index_prof[k,"number_peaks_total"]>1){
        ########################################################################
        # fill timeset #########################################################
        timeset[,4:length(timeset[1,])]<-0;
        timeset[,c(4,5)] <-.Call("_enviMass_fill_timeset",
                                as.numeric(timeset),
                                as.numeric(peaks[(index_prof[k,"start_ID"]:index_prof[k,"end_ID"]),"sampleIDs"]), 
                                as.numeric(peaks[(index_prof[k,"start_ID"]:index_prof[k,"end_ID"]),"intensity"]), 
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
				that<-.Call("_enviMass_meandel",
					as.numeric(timeset), 	# ok
					as.integer(subit),		# ok
					as.numeric(subrat),		# ok
					as.numeric(numtime),	# ok
					as.integer(what),		# ok
					as.numeric(lags),		# ok
					as.numeric(threshold),	# ok
					as.integer(notrend),	# 
					PACKAGE="enviMass"
				)
			}
			if(!omit_trend){
				index_prof[k,"deltaint_newest"]<-max(that[4,]); # current>abs.dev
				index_prof[k,"deltaint_global"]<-max(that[5,]); # global>abs.dev
				index_prof[k,"absolute_mean_dev"]<-max(that[3,]); # abs.dev
			}
			#if(any(timeset[,5]>0)){ 			  # in blind?
			#	index_prof[k,"in_blind?"]<-1 # in blind
			#	index_prof[k,"number_peaks_blind"]<-length(timeset[timeset[,5]!=0,5]) 	# number_peaks_blind
			#	index_prof[k,"mean_int_blind"]<-mean(timeset[timeset[,5]!=0,5]) 	# mean_int_blind
			#	index_prof[k,"max_int_blind"]<-max(timeset[timeset[,5]!=0,5]) 	# mean_int_blind					
			#}else{
			#	index_prof[k,"in_blind?"]<-0 # in blind
			#	index_prof[k,"number_peaks_blind"]<-0 # number_peaks_blind
			#	index_prof[k,"mean_int_blind"]<-0 # mean_int_blind		
			#	index_prof[k,"max_int_blind"]<-0 # mean_int_blind						
			#}
			#############################################
			#index_prof[k,"number_peaks_sample"]<-length(timeset[timeset[,4]!=0,4]) 	# number_peaks_sample
			#index_prof[k,"mean_int_sample"]<-mean(timeset[timeset[,4]!=0,4]) 	# mean_int_sample	
			#index_prof[k,"max_int_sample"]<-max(timeset[timeset[,4]!=0,4]) 						
		}#else{
			#index_prof[k,"in_blind?"]<-1 # only in blind = not above (single peak)
			#index_prof[k,"number_peaks_sample"]<-0 # number_peaks_sample
			#index_prof[k,"number_peaks_blind"]<-length(timeset[timeset[,5]!=0,5])# number_peaks_blind
			#index_prof[k,"mean_int_sample"]<-0 # mean_int_sample
			#index_prof[k,"max_int_sample"]<-0 
			#index_prof[k,"mean_int_blind"]<-mean(timeset[timeset[,5]!=0,5]) # mean_int_blind			
			#index_prof[k,"max_int_blind"]<-max(timeset[timeset[,5]!=0,5])			
		#}
		index_prof[k,"newest_intensity"]<-timeset[leng,4]
		########################################################################
      }else{ # single-peaked profile
		index_prof[k,"absolute_mean_dev"]<-0;
		if( any( timeset[,3]== (peaks[index_prof[k,1],6]) ) ){ # only in blind
			#index_prof[k,"in_blind?"]<-1 # in blind = not above (single peak)
			#index_prof[k,"number_peaks_sample"]<-0 # number_peaks_sample
			#index_prof[k,"number_peaks_blind"]<-1 # number_peaks_blind
			#index_prof[k,"mean_int_sample"]<-0 # mean_int_sample
			#index_prof[k,"max_int_sample"]<-0 
			#index_prof[k,"mean_int_blind"]<-(peaks[index_prof[k,1],2]) # mean_int_blind	
			#index_prof[k,"max_int_blind"]<-(peaks[index_prof[k,1],2])			
		}else{ # only in sample
			#index_prof[k,"in_blind?"]<-0 
			index_prof[k,"deltaint_global"]<-(peaks[index_prof[k,1],2])
			if(any(peaks[index_prof[k,1],6]==latestID)){
			  index_prof[k,"deltaint_newest"]<-(peaks[index_prof[k,1],2])
			  index_prof[k,"newest_intensity"]<-(peaks[index_prof[k,1],2])
			}
			#index_prof[k,"number_peaks_sample"]<-1 # number_peaks_sample
			#index_prof[k,"number_peaks_blind"]<-0 # number_peaks_blind
			#index_prof[k,"mean_int_sample"]<-(peaks[index_prof[k,1],2]) # mean_int_sample
			#index_prof[k,"max_int_sample"]<-(peaks[index_prof[k,1],2]) 
			#index_prof[k,"mean_int_blind"]<-0 # mean_int_blind		
			#index_prof[k,"max_int_blind"]<-0				
		}
      }
    }
    ############################################################################
    return(index_prof)

}
