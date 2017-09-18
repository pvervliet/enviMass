#' @title Extract profiles from peak set partitions
#'
#' @export
#'
#' @description \code{partcluster} extract profiles from peak set partitions generated by \code{\link{agglomer}}.
#'
#' @param profileList A profile list.
#' @param dmass Numeric. m/z gap size
#' @param ppm Logical. \code{dmass} given in ppm?
#' @param dret Numeric. RT gap size; units equal to those of the input files
#' @param from Logical FALSE or integer. Restrict to certain partitons.
#' @param to Logical FALSE or integer. Restrict to certain partitions.
#' @param progbar Logical. Should a progress bar be shown? Only for Windows.
#' @param plot_it Logical. Plot profile extraction? For debugging.
#' @param replicates. FALSE or vector of character strings with replicate labels for each LC-MS file, i.e., files of the same replicate have the same character string.
#' @param IDs. Integer vector. IDs of files, required if replicates is not set to FALSE.
#' @param with_test. Logical. Do some sporadic internal testing of results?
#'
#' @return Updated profile list
#' 
#' @details enviMass workflow function. Works along decreasing intensities. The remaining peak of highest intensity not
#' yet part of a profile is either assigned to an existing profile (closest in mass) or initializes a new profile. With
#' addition of a peak to a new profile, profile mass tolerances are gradually adapted. If replicates are profiled, profiles
#' are first extracted in each replicate level and these profiles than further merged.
#' 
#' @seealso \code{\link{startprofiles}}, \code{\link{agglomer}}


partcluster<-function(
	profileList,
	dmass=3,
	ppm=TRUE,
	dret=60,
	from=FALSE,
	to=FALSE,
	progbar=FALSE,
	plot_it=FALSE,
	replicates=FALSE,
	IDs=FALSE,
	with_test=FALSE
){

	########################################################################################
	if(!profileList[[1]][[2]]){stop("run agglom first on that profileList; aborted.")}
	if(!is.numeric(dmass)){stop("dmass must be numeric; aborted.")}
	if(!is.numeric(dret)){stop("dret must be numeric; aborted.")}
	if(!is.logical(ppm)){stop("ppm must be logical; aborted.")}
	if(!from){m=1}else{m=from}
	if(!to){n=length(profileList[["index_agglom"]][,1])}else{n=to}
	if( (from!=FALSE) || (to!=FALSE) ){
		startat<-c(1);
		profileList[[2]][,8]<-1;
	}else{
		startat<-c(0);
	}
	if(ppm){ppm2=1}else{ppm2=2};
	do_replicates<-FALSE
	if(any(replicates!="FALSE")){
		replic<-replicates[duplicated(replicates)]
		replic<-unique(replic)			  
		replic<-replic[replic!="FALSE"]
		if(length(replic)>0){
			if(length(replicates)!=length(IDs)){
				stop("\n replicates vector longer than file ID vector")
			}
			do_replicates<-TRUE
		}
	}
	########################################################################################
	if(progbar==TRUE){prog<-winProgressBar("Extract time profiles...",min=m,max=n);setWinProgressBar(prog, 0, title = "Extract time profiles...", label = NULL);}
	often<-c(0);
	for(k in m:n){
		#if(progbar==TRUE){setWinProgressBar(prog, k, title = paste("Extract time profiles for ",(profileList[["index_agglom"]][k,2]-profileList[["index_agglom"]][k,1]+1)," peaks",sep=""), label = NULL)}
		if(progbar==TRUE){setWinProgressBar(prog, k, title = paste("Extract time profiles for partition ",k,sep=""), label = NULL)}
		if(profileList[["index_agglom"]][k,3]>1){
			delmz<-(max(profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),1])-min(profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),1]))
			if(ppm){
				delmz<-(delmz/mean(profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),1])*1E6)
			}else{
				delmz<-(delmz/1000)
			}
			delRT<-(max(profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),3])-min(profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),3]))
			if( (delmz>(dmass*2)) || (delRT>dret) || (any(duplicated(profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),6]))) 
#|| (do_replicates & (...))			
			){  # check dmass & dret & uniqueness & replicates
				if(!do_replicates){
					if(TRUE){	# REMOVE STATEMENT!
					often<-c(often+1)
					clusters <-.Call("_enviMass_extractProfiles",
									  as.numeric(profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),"m/z"]),       # mz
									  as.numeric(profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),"RT"]),       # RT
									  as.numeric(profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),"intensity"]),       # intens
									  as.integer(profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),"sampleIDs"]),       # sampleID                          
									  as.integer(order(profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),"intensity"],decreasing=TRUE)),   	# intensity order 
									  as.integer(0), 		# no replicate pre-sorting
									  as.numeric(dmass),
									  as.integer(ppm2),
									  as.numeric(dret),
									  as.integer(0),		# replicates
									  PACKAGE="enviMass"
									)				 
					if(any(clusters[,10]==0)){stop("missing assignment found_1!")} # check - all must be assigned
					clusters[clusters[,10]!=0,10]<-(clusters[clusters[,10]!=0,10]+startat); 
					profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),8]<-clusters[,10];
					profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),]<-(profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),][order(clusters[,10],decreasing=FALSE),]);
					###
					if(with_test){
						cat("\nTEST\n");
						these<-match( # all clusters matched in peaks?
							clusters[,10],
							profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),"profileIDs"]
						)
						if(any(is.na(these))){ stop("\n partcluster issue_A")}
						these<-match( # clusters continuous?
							clusters[,10],
							seq(max(clusters[,10]))
						)
						if(any(is.na(these))){ stop("\n partcluster issue_B")}
						these<-match( # sequential profileIDs?
							seq(
								min(profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),"profileIDs"]),
								max(profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),"profileIDs"])
							),
							profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),"profileIDs"]
						)
						if(any(is.na(these))){ stop("\n partcluster issue_C")}		
						if(sum(clusters[,5])!=length(clusters[,5])){stop("\n partcluster issue: wrong peak sum")}						
					}
					###
					startat<-c(max(clusters[,10]))
					########################################################################
					if(plot_it){
							profileList[[2]][(profileList[[6]][k,1]:profileList[[6]][k,2]),]<-
							profileList[[2]][(profileList[[6]][k,1]:profileList[[6]][k,2]),][
							order(profileList[[2]][(profileList[[6]][k,1]:profileList[[6]][k,2]),6]),]
							plot(
							  profileList[[2]][(profileList[[6]][k,1]:profileList[[6]][k,2]),1],
							  profileList[[2]][(profileList[[6]][k,1]:profileList[[6]][k,2]),3],
							  pch=19,xlab="m/z",ylab="RT",cex=0.5,col="lightgrey"
							)
							atID<-unique(profileList[[2]][(profileList[[6]][k,1]:profileList[[6]][k,2]),6])
							seqID<-c(1:length(profileList[[2]][(profileList[[6]][k,1]:profileList[[6]][k,2]),6]))
							for(i in 1:length(atID)){
							  if(length(seqID[profileList[[2]][(profileList[[6]][k,1]:profileList[[6]][k,2]),6]==atID[i]])>1){
								subseqID<-seqID[profileList[[2]][(profileList[[6]][k,1]:profileList[[6]][k,2]),6]==atID[i]]
								for(m in 1:(length(subseqID)-1)){
								  for(n in (m+1):length(subseqID)){
								   lines(
									  c(profileList[[2]][(profileList[[6]][k,1]:profileList[[6]][k,2]),1][subseqID[m]],profileList[[2]][(profileList[[6]][k,1]:profileList[[6]][k,2]),1][subseqID[n]]),
									  c(profileList[[2]][(profileList[[6]][k,1]:profileList[[6]][k,2]),3][subseqID[m]],profileList[[2]][(profileList[[6]][k,1]:profileList[[6]][k,2]),3][subseqID[n]]),
									  col="lightgrey",lwd=1
									)
								  }
								}
							  }
							}
							meanmass<-mean(profileList[[2]][(profileList[[6]][k,1]:profileList[[6]][k,2]),1])
							meanRT<-mean(profileList[[2]][(profileList[[6]][k,1]:profileList[[6]][k,2]),3])
							lengppm<-(as.numeric(dmass)*meanmass/1e6)
							lines(
							  c((meanmass-0.5*lengppm),(meanmass+0.5*lengppm)),
							  c(meanRT,meanRT),
							  col="blue",lwd=3
							)
							lines(
							  c(meanmass,meanmass),
							  c((meanRT-0.5*dret),(meanRT+0.5*dret)),
							  col="blue",lwd=3
							)
							clust<-profileList[[2]][(profileList[[6]][k,1]:profileList[[6]][k,2]),8]
							clust<-(clust-min(clust)+1)
							colorit<-sample(colors(),max(clust))
							points(
							  profileList[[2]][(profileList[[6]][k,1]:profileList[[6]][k,2]),1],
							  profileList[[2]][(profileList[[6]][k,1]:profileList[[6]][k,2]),3],
							  pch=19,cex=1,col=colorit[clust]
							)
							text(
							  profileList[[2]][(profileList[[6]][k,1]:profileList[[6]][k,2]),1],
							  profileList[[2]][(profileList[[6]][k,1]:profileList[[6]][k,2]),3],
							  labels=as.character(clust),
							  pch=19,cex=1,col="darkgrey"
							)
						}
					########################################################################
					}
				}else{	
					# first run profile extraction in replicates only ######################
					startat_inner<-0
					profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),"profileIDs"]<-0
					for(i in 1:length(replic)){
						IDs2<-as.numeric(IDs[replicates==replic[i]])
						those<-!is.na(match(profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),6],IDs2))
						if(sum(those)>1){
							often<-c(often+1)					
							clusters_rep <-.Call("_enviMass_extractProfiles",
											  as.numeric(profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),"m/z"][those]),       		# mz
											  as.numeric(profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),"RT"][those]),       		# RT
											  as.numeric(profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),"intensity"][those]),       # intens
											  as.integer(profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),"sampleIDs"][those]),       # sampleID                          
											  as.integer(order(profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),"intensity"][those],decreasing=TRUE)),   # intensity order                        
											  as.integer(0),	# pre-ordering
											  as.numeric(dmass),
											  as.integer(ppm2),
											  as.numeric(dret),					  
											  as.integer(0),
											  PACKAGE="enviMass"
											)
							clusters_rep[clusters_rep[,10]!=0,10]<-(clusters_rep[clusters_rep[,10]!=0,10]+startat_inner); 
							if(with_test){
								if(any(profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),"profileIDs"][those]!=0)){stop("\n partcluster issue: overwriting replicates?!")}
							}
							profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),"profileIDs"][those]<-clusters_rep[,10];	
							startat_inner<-c(max(clusters_rep[,10]))
							if(any(clusters_rep[,10]==0)){stop("missing assignment found_2!")}							
						}						
					}
					if(with_test){
						cat("\nTEST\n");
						# no duplicated sample IDs in a cluster!
						clstr<-unique(profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),"profileIDs"])
						for(w in 1:length(clstr)){
							if(clstr[w]==0){next}
							these<-which(profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),"profileIDs"]==clstr[w])
							if(any(duplicated(profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),"sampleIDs"][these]))){stop("\n partcluster issue in duplication")}
							min_RT<-min(profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),"RT"][these])
							max_RT<-max(profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),"RT"][these])
							if((max_RT-min_RT)>dret){stop("\n partcluster issue: too large RT windows detected in replicates!")}
							min_mass<-min(profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),"m/z"][these])
							max_mass<-max(profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),"m/z"][these])	
							mean_mass<-mean(profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),"m/z"][these])								
							if(ppm){
								test_delmz<-((max_mass-min_mass)/mean_mass*1E6)
							}else{
								test_delmz<-(max_mass-min_mass)
							}
							if(test_delmz>(dmass*2)){stop("\n partcluster issue: too large mass windows detected in replicates!")}
						}	
					}
					# pre-order peaks in the partition by their replicate clusters = "profileIDs": required by "getProfiles" ####
					profileList[["peaks"]][(profileList[["index_agglom"]][k,"start_ID"]:profileList[["index_agglom"]][k,"end_ID"]),]<-(
						profileList[["peaks"]][(profileList[["index_agglom"]][k,"start_ID"]:profileList[["index_agglom"]][k,"end_ID"]),][
							order(profileList[["peaks"]][(profileList[["index_agglom"]][k,"start_ID"]:profileList[["index_agglom"]][k,"end_ID"]),"profileIDs"],decreasing=FALSE)
						,]
					);
					# extract replicates with replicate-preordering ########################
					clusters <-.Call("_enviMass_extractProfiles",
									  as.numeric(profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),"m/z"]),       # mz
									  as.numeric(profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),"RT"]),       # RT
									  as.numeric(profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),"intensity"]),       # intens
									  as.integer(profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),"sampleIDs"]),       # sampleID                          
									  as.integer(order(profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),"intensity"],decreasing=TRUE)),   	# intensity order 
									  as.integer(profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),"profileIDs"]), 				# with replicate pre-sorting
									  as.numeric(dmass),
									  as.integer(ppm2),
									  as.numeric(dret),
									  as.integer(1),		# replicates
									  PACKAGE="enviMass"
									)			
					clusters[clusters[,10]!=0,10]<-(clusters[clusters[,10]!=0,10]+startat); 
					profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),"profileIDs"]<-clusters[,10];		
					###
					if(with_test){
						cat("\nTEST\n");
						these<-match(
							clusters[,10],
							profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),"profileIDs"]
						)
						if(any(is.na(these))){ stop("\n partcluster issue_A")}
						these<-match(
							clusters[,10],
							seq(max(clusters[,10]))
						)
						if(any(is.na(these))){ stop("\n partcluster issue_B")}
						these<-match( # sequential profileIDs?
							seq(
								min(profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),"profileIDs"]),
								max(profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),"profileIDs"])
							),
							profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),"profileIDs"]
						)
						if(any(is.na(these))){ stop("\n partcluster issue_C")}		
						if(sum(clusters[,5])!=length(clusters[,5])){stop("\n partcluster issue: wrong peak sum at replicates")}
						# unique sample IDs in a cluster?
						clstr<-profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),"profileIDs"]
						for(w in 1:length(clstr)){
							these<-which(profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),"profileIDs"]==clstr[w])
							if(any(duplicated(profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),"sampleIDs"][these]))){stop("\n partcluster issue in duplication")}
							min_RT<-min(profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),"RT"][these])
							max_RT<-max(profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),"RT"][these])
							if((max_RT-min_RT)>dret){stop("\n partcluster issue: too large RT windows detected!")}
							min_mass<-min(profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),"m/z"][these])
							max_mass<-max(profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),"m/z"][these])	
							mean_mass<-mean(profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),"m/z"][these])								
							if(ppm){
								test_delmz<-((max_mass-min_mass)/mean_mass*1E6)
							}else{
								test_delmz<-(max_mass-min_mass)
							}
							if(test_delmz>(dmass*2)){stop("\n partcluster issue: too large mass windows detected!")}
						}	
					}
					###
					profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),]<-(profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[[6]][k,2]),][order(clusters[,10],decreasing=FALSE),]);
					startat<-c(max(clusters[,10]))
					########################################################################
				}
			}else{
				startat<-(startat+1);
				profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),"profileIDs"]<-startat;
			}
		}else{
			startat<-(startat+1);
			profileList[["peaks"]][(profileList[["index_agglom"]][k,1]:profileList[["index_agglom"]][k,2]),"profileIDs"]<-startat;
		}
	}
	if(progbar==TRUE){close(prog);}
	########################################################################################
	# assemble index matrix ################################################################
	if(with_test){
		those<-match(
			as.integer(profileList[["peaks"]][,"profileIDs"]),
			seq(startat)
		)	
		if(any(is.na(those))){stop("\n partcluster issue_D")}
		those<-match(
			seq(startat),
			as.integer(profileList[["peaks"]][,"profileIDs"])
		)	
		if(any(is.na(those))){stop("\n partcluster issue_E")}
		for(j in 2:length(profileList[["peaks"]][,1])){
			if(
				(profileList[["peaks"]][j,"profileIDs"]-profileList[["peaks"]][(j-1),"profileIDs"])>1
			){
				stop()
			}
		}
	}
	index <- .Call("_enviMass_indexed",
		as.integer(profileList[[2]][,"profileIDs"]),
		as.integer(startat),
		as.integer(26),
		PACKAGE="enviMass"
	)
	index<-index[index[,1]!=0,]
	index[,4]<-seq(length(index[,4]))
	colnames(index)<-c(
		"start_ID",
		"end_ID",
		"number_peaks_total", #1
		"profile_ID",
		"deltaint_newest", #"current_incident"
		"deltaint_global", #4 #"past_incident"
		"absolute_mean_dev",
		"in_blind?",
		"above_blind?", #7
		"number_peaks_sample",
		"number_peaks_blind", #10
		"mean_int_sample",
		"mean_int_blind", #12
		"mean_mz",
		"mean_RT",
		"mean_int", #14
		"newest_intensity", #"newest_intensity"
		"links",
		"component", #17
		"max_int_sample",
		"max_int_blind",
		# new:
		"max_int",
		"var_mz",
		"min_RT",
		"max_RT",
		"Mass defect"
	)
	profileList[[7]]<-index
	if(with_test){
		those<-match(
			as.integer(profileList[["peaks"]][,"profileIDs"]),
			index[,4]
		)	
		if(any(is.na(those))){stop("\n partcluster issue_F")}
		those<-match(
			index[,4],
			as.integer(profileList[["peaks"]][,"profileIDs"])
		)	
		if(any(is.na(those))){stop("\n partcluster issue_G")}
	}
	########################################################################################
	# get characteristics of individual profiles ###########################################
	m=1
	n=length(profileList[["index_prof"]][,1])
    if(progbar==TRUE){  prog<-winProgressBar("Extract profile data...",min=m,max=n);
						setWinProgressBar(prog, 0, title = "Extract profile data...", label = NULL);}
    for(k in m:n){
		if(progbar==TRUE){setWinProgressBar(prog, k, title = "Extract profile data...", label = NULL)}
			profileList[["index_prof"]][k,"mean_mz"]<-mean(profileList[["peaks"]][(profileList[["index_prof"]][k,1]:profileList[["index_prof"]][k,2]),"m/z"])
			profileList[["index_prof"]][k,"mean_RT"]<-mean(profileList[["peaks"]][(profileList[["index_prof"]][k,1]:profileList[["index_prof"]][k,2]),"RT"])	  
			profileList[["index_prof"]][k,"mean_int"]<-mean(profileList[["peaks"]][(profileList[["index_prof"]][k,1]:profileList[["index_prof"]][k,2]),"intensity"])	
			profileList[["index_prof"]][k,"max_int"]<-max(profileList[["peaks"]][(profileList[["index_prof"]][k,"start_ID"]:profileList[["index_prof"]][k,"end_ID"]),"intensity"])
			profileList[["index_prof"]][k,"var_mz"]<-var(profileList[["peaks"]][(profileList[["index_prof"]][k,"start_ID"]:profileList[["index_prof"]][k,"end_ID"]),"m/z"])
			profileList[["index_prof"]][k,"min_RT"]<-min(profileList[["peaks"]][(profileList[["index_prof"]][k,"start_ID"]:profileList[["index_prof"]][k,"end_ID"]),"RT"])	
			profileList[["index_prof"]][k,"max_RT"]<-max(profileList[["peaks"]][(profileList[["index_prof"]][k,"start_ID"]:profileList[["index_prof"]][k,"end_ID"]),"RT"])	
			profileList[["index_prof"]][k,"Mass defect"]<-(round(profileList[["index_prof"]][k,"mean_mz"])-profileList[["index_prof"]][k,"mean_mz"])
	}
	if(progbar==TRUE){close(prog);}
	profileList[[1]][[3]]<-TRUE
	########################################################################################
	return(profileList)

}




