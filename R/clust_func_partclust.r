

clust_func_partclust<-function(
		peaks, 
		agglom, 
		...
	){ 
		###################################################################################
		m = 1
		n = dim(agglom)[1]
		startat<-c(0);
		for(k in m:n){
			if(agglom[k,"number_peaks"]>1){
				delmz<-(max(peaks[agglom[k,"start_ID"]:agglom[k,"end_ID"],"m/z"])-min(peaks[agglom[k,"start_ID"]:agglom[k,"end_ID"],"m/z"]))
				if(ppm2 == 1){
					delmz<-(delmz/mean(peaks[agglom[k,"start_ID"]:agglom[k,"end_ID"],"m/z"])*1E6)
				}else{
					delmz<-(delmz/1000)
				}
				delRT<-(max(peaks[agglom[k,"start_ID"]:agglom[k,"end_ID"],"RT"])-min(peaks[agglom[k,"start_ID"]:agglom[k,"end_ID"],"RT"]))
				if( 
					(delmz>(dmass*2)) || 
					(delRT>dret) || 
					(any(duplicated(peaks[agglom[k,"start_ID"]:agglom[k,"end_ID"],"sampleIDs"]))) 	
				){  # check dmass & dret & uniqueness & replicates
					#######################################################################
					if(!do_replicates){
						peaks[agglom[k,"start_ID"]:agglom[k,"end_ID"],"profileIDs"]<-0
						clusters <-.Call("_enviMass_extractProfiles",
										  as.numeric(peaks[agglom[k,"start_ID"]:agglom[k,"end_ID"],"m/z"]),       	# mz
										  as.numeric(peaks[agglom[k,"start_ID"]:agglom[k,"end_ID"],"RT"]),       	# RT
										  as.numeric(peaks[agglom[k,"start_ID"]:agglom[k,"end_ID"],"intensity"]), 	# intens
										  as.integer(peaks[agglom[k,"start_ID"]:agglom[k,"end_ID"],"sampleIDs"]), 	# sampleID                          
										  as.integer(order(peaks[agglom[k,"start_ID"]:agglom[k,"end_ID"],"intensity" ],decreasing=TRUE)),   	# intensity order 
										  as.integer(0), 		# no replicate pre-sorting
										  as.numeric(dmass),
										  as.integer(ppm2),
										  as.numeric(dret),
										  as.integer(0),		# replicates
										  PACKAGE="enviMass"
										)				 
						if(any(clusters[,10]==0)){stop("missing assignment found_1!")} # check - all must be assigned
						clusters[clusters[,10]!=0,10]<-(clusters[clusters[,10]!=0,10]+startat); 
						peaks[agglom[k,"start_ID"]:agglom[k,"end_ID"],"profileIDs"]<-clusters[,10];
						peaks[agglom[k,"start_ID"]:agglom[k,"end_ID"],]<-(peaks[agglom[k,"start_ID"]:agglom[k,"end_ID"],][order(clusters[,10],decreasing=FALSE),]);
						startat<-c(max(clusters[,10]))
						###################################################################
					}else{
						###################################################################
						# first run profile extraction in replicates only #################
						startat_inner<-0
						peaks[agglom[k,"start_ID"]:agglom[k,"end_ID"],"profileIDs"]<-0
						for(i in 1:length(replic)){
							IDs2<-as.numeric(IDs[replicates==replic[i]])
							those<-!is.na(match(peaks[agglom[k,"start_ID"]:agglom[k,"end_ID"],"sampleIDs"],IDs2))
							if(sum(those)>1){
								often<-c(often+1)					
								clusters_rep <-.Call("_enviMass_extractProfiles",
												  as.numeric(peaks[agglom[k,"start_ID"]:agglom[k,"end_ID"],"m/z"][those]),       	# mz
												  as.numeric(peaks[agglom[k,"start_ID"]:agglom[k,"end_ID"],"RT"][those]),       	# RT
												  as.numeric(peaks[agglom[k,"start_ID"]:agglom[k,"end_ID"],"intensity"][those]),  	# intens
												  as.integer(peaks[agglom[k,"start_ID"]:agglom[k,"end_ID"],"sampleIDs" ][those]), 	# sampleID                          
												  as.integer(order(peaks[agglom[k,"start_ID"]:agglom[k,"end_ID"],"intensity"][those],decreasing=TRUE)),   # intensity order                        
												  as.integer(0),	# pre-ordering
												  as.numeric(dmass),
												  as.integer(ppm2),
												  as.numeric(dret),					  
												  as.integer(0),
												  PACKAGE="enviMass"
												)
								clusters_rep[clusters_rep[,10]!=0,10]<-(clusters_rep[clusters_rep[,10]!=0,10]+startat_inner); 
								peaks[agglom[k,"start_ID"]:agglom[k,"end_ID"],"profileIDs"][those]<-clusters_rep[,10];	
								startat_inner<-c(max(clusters_rep[,10]))
								if(any(clusters_rep[,10]==0)){stop("missing assignment found_2!")}							
							}						
						}
						###################################################################
						# pre-order peaks in the partition by their replicate clusters = "profileIDs": required by "getProfiles" ####
						peaks[agglom[k,"start_ID"]:agglom[k,"end_ID"],]<-(
							peaks[agglom[k,"start_ID"]:agglom[k,"end_ID"],][
								order(peaks[agglom[k,"start_ID"]:agglom[k,"end_ID"],"profileIDs"],decreasing=FALSE)
							,]
						);
						###################################################################
						# extract replicates with replicate-preordering ###################
						clusters <-.Call("_enviMass_extractProfiles",
										as.numeric(peaks[agglom[k,"start_ID"]:agglom[k,"end_ID"],"m/z"][those]),       	# mz
										as.numeric(peaks[agglom[k,"start_ID"]:agglom[k,"end_ID"],"RT"][those]),       	# RT
										as.numeric(peaks[agglom[k,"start_ID"]:agglom[k,"end_ID"],"intensity"][those]),  # intens
										as.integer(peaks[agglom[k,"start_ID"]:agglom[k,"end_ID"],"sampleIDs" ][those]), # sampleID                          
										as.integer(order(peaks[agglom[k,"start_ID"]:agglom[k,"end_ID"],"intensity"][those],decreasing=TRUE)),  # intensity order 
										as.integer(peaks[agglom[k,"start_ID"]:agglom[k,"end_ID"],"profileIDs"]), 				# with replicate pre-sorting
										as.numeric(dmass),
										as.integer(ppm2),
										as.numeric(dret),
										as.integer(1),		# replicates
										PACKAGE="enviMass"
									)			
						clusters[clusters[,10]!=0,10]<-(clusters[clusters[,10]!=0,10]+startat); 
						peaks[agglom[k,"start_ID"]:agglom[k,"end_ID"],"profileIDs"]<-clusters[,10];		
						peaks[agglom[k,"start_ID"]:agglom[k,"end_ID"],]<-(peaks[agglom[k,"start_ID"]:agglom[k,"end_ID"],][order(clusters[,10],decreasing=FALSE),]);
						startat<-c(max(clusters[,10]))
						###################################################################
					}
					#######################################################################
				}else{
					startat<-(startat+1);
					peaks[agglom[k,"start_ID"]:agglom[k,"end_ID"],"profileIDs"]<-startat					
				}
			}else{
				startat<-(startat+1);
				peaks[agglom[k,"start_ID"]:agglom[k,"end_ID"],"profileIDs"]<-startat
			}
		}
		###################################################################################
		return(peaks)
	}









