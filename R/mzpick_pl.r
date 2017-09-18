mzpick_pl <- function(           
	clus_EICs, 		# MSlist[[6]]
	clus_centroids,	# MSlist[[4]][[2]]
	...
){

	############################################################################
    startat<-c(0);
    clus_centroids[,"measureID"]<-seq(1,length(clus_centroids[,"measureID"]),1); 
	m <- 1
	n <- length(clus_EICs[,1])
    for(i in m:n){
		if(clus_EICs[i,3]>=minpeak){
		  # same RT, interpolate ###################################################
		  out1 <- .Call("gapfill",
			as.numeric(clus_centroids[clus_EICs[i,"start_ID"]:clus_EICs[i,"end_ID"],"RT"]),
			as.numeric(clus_centroids[clus_EICs[i,"start_ID"]:clus_EICs[i,"end_ID"],"intensity"]),
			as.integer(order(clus_centroids[clus_EICs[i,"start_ID"]:clus_EICs[i,"end_ID"],"RT"],decreasing=FALSE)),
			as.numeric(clus_centroids[clus_EICs[i,"start_ID"]:clus_EICs[i,"end_ID"],"m/z"]),
			as.numeric(clus_centroids[clus_EICs[i,"start_ID"]:clus_EICs[i,"end_ID"],"measureID"]),
			as.numeric(Retens),
			as.numeric(drtfill),
			PACKAGE="enviPick"
			)
		  out1<-matrix(out1,ncol=10);
		  colnames(out1)<-c("m/z","intens","RT","index","intens_filt","1pick","pickcrit","baseline","intens_corr","2pick");                       
		  # Filter step ############################################################
		  out1[,5]<-out1[,2];   # yet to be implemented!
		  # 1st peak detection & baseline subtraction & 2nd peak detection ##########  
		  out2 <- .Call("pickpeak",          
			as.numeric(out1),
			as.numeric(drtsmall),
			as.numeric(drttotal),
			as.integer(minpeak),
			as.integer(recurs),
			as.numeric(weight),   # weight
			as.numeric(SB),       # SB 
			as.numeric(SN),       # SN                                                                           
			as.numeric(minint),   # minimum intensity
			as.numeric(maxint),   # maximum intensity threshold
			as.integer(ended),
			as.integer(2),
			PACKAGE="enviPick"
		  )     
		  out2<-matrix(out2,ncol=10);
		  colnames(out2)<-c("m/z","intens","RT","index","intens_filt","1pick","pickcrit","baseline","intens_corr","2pick");
		  # assign final entries ####################################################
		  if(!all(out2[,"2pick"]==0)){      
			out2[,"2pick"]<-out2[,"2pick"]+startat;         
			for(k in 1:length(out2[,"2pick"])){
			  if(out2[k,"2pick"]!=startat){
				clus_centroids[out2[k,4],"peakID"]<-out2[k,"2pick"]
			  }
			}        
			startat<-c(max(out2[,"2pick"]));
			clus_centroids[clus_EICs[i,"start_ID"]:clus_EICs[i,"end_ID"],]<-
			  clus_centroids[clus_EICs[i,"start_ID"]:clus_EICs[i,"end_ID"],][order(
				clus_centroids[clus_EICs[i,"start_ID"]:clus_EICs[i,"end_ID"],"peakID"],decreasing=FALSE),];   
		  }
		}
    }              
    ############################################################################	
	return(clus_centroids)
	
}



             