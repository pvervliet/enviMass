# Run homologue series detection
	
	####################################################################################	
    measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
	cat("Homologue series detection: ")
	if(logfile$parameters$homol_units[1]!="FALSE"){
		these<-enviPat::check_chemform(isotopes,strsplit(logfile$parameters$homol_units,",")[[1]])[,3]
		mzfilter<-c(enviPat::check_chemform(isotopes,strsplit(logfile$parameters$homol_units,",")[[1]])[,3] %*% t(1/as.numeric(strsplit(logfile$parameters$homol_charges,",")[[1]])))
		mzfilter<-unique(mzfilter);
		elements<-unique(unlist(sapply(enviMass::check_chemform(isotopes,strsplit(logfile$parameters$homol_units,",")[[1]],get_list=TRUE),names)))
		use_minmz<-(min(mzfilter)-.1)
		use_maxmz<-(max(mzfilter)+.1)
	}else{
		mzfilter<-FALSE
		elements<-unique(as.character(isotopes[,1])[1:295]) #then use all available elements
		use_minmz<-as.numeric(logfile$parameters$homol_minmz)
		use_maxmz<-as.numeric(logfile$parameters$homol_maxmz)		
	}
	# intersect with compounds 
	intstand<-read.table(file=file.path(logfile[[1]],"dataframes","IS.txt"),header=TRUE,sep="\t",colClasses = "character");
	targets<-read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character");
	####################################################################################		

# >	
#measurements[,names(measurements)=="homologues"]<-"FALSE"
# <

	####################################################################################		
	for(b in 1:length(measurements[,"ID"])){
		if( 
			(measurements[b,names(measurements)=="include"]=="TRUE") & 			# included?
			(measurements[b,names(measurements)=="homologues"]=="FALSE")  	# not yet done
		){ 
		
			##########################################################################
			# exclude files that do not end up in profiles ###########################
			if( (mute(logfile$parameters$prof_select=="TRUE")) & (measurements$profiled[b]=="FALSE") ){
				cat("\n Skip file - not included in profile building.");next;
			}
			##########################################################################		
			cat(paste("\n Doing ",as.character(b)," of ",as.character(length(measurements[,"ID"]))," files: ",sep=""))			
			##########################################################################
			# LOAD FILES & REMOVE OLD RESULTS ########################################
			cat("loading - ")
			for_file<-measurements[b,1]
			if( file.exists(file.path(logfile[[1]],"results","componentization","homologues",for_file) ) ){
				file.remove(file.path(logfile[[1]],"results","componentization","homologues",for_file) )
			}			
			# Peaklist 
			load(file=file.path(logfile[[1]],"peaklist",as.character(for_file))); 
			peaklist4<-peaklist[order(peaklist[,"peak_ID"],decreasing=FALSE),] # match with IDs 
			if(logfile$parameters$homol_blind=="TRUE"){ # remove blind peaks
				peaklist4<-peaklist4[peaklist4[,"keep_2"]>=as.numeric(logfile$parameters$homol_blind_value),,drop=FALSE]
				cat("blind peaks removed -")
			}
			peaklist4<-as.data.frame(peaklist4[(peaklist4[,"keep"]==1),c("m/z_corr","int_corr","RT_corr","peak_ID"),drop=FALSE])
			##########################################################################
			cat("series extraction - ")		
			if(logfile$parameters$homol_ppm=="TRUE"){
				use_mztol<-as.numeric(logfile$parameters$homol_mztol)
			}else{ # mmu
				use_mztol<-(as.numeric(logfile$parameters$homol_mztol)/1000)
			}			
			homol<-try(
				enviMass::homol_search2(
					peaklist=peaklist4[,c("m/z_corr","int_corr","RT_corr","peak_ID")],
					isotopes,
					elements=elements,
					use_C=FALSE,
					minmz=use_minmz,
					maxmz=use_maxmz,
					minrt=as.numeric(logfile$parameters$homol_minrt),
					maxrt=as.numeric(logfile$parameters$homol_maxrt),
					ppm=as.logical(logfile$parameters$homol_ppm),
					mztol=use_mztol,
					rttol=as.numeric(logfile$parameters$homol_rttol),
					minlength=as.numeric(logfile$parameters$homol_minlength),
					mzfilter,
					vec_size=as.numeric(logfile$parameters$homol_vec_size),
					mat_size=3,
					R2=.98,
					spar=.45,
					plotit=FALSE,
					deb=0
				)	
			)
			if(class(homol)=="try-error"){
				cat("\n Homologue series detection failed - adapt parameters?");
				next;
			}	
			if(length(homol[["Homologue Series"]][,1])==0){
				cat("\n No homologue series detected for this file. Continue to next ...");
				next;			
			}
			atby<-1000
			Homol_groups<-matrix(nrow=atby,ncol=3,0)
			at<-atby
			from<-0
			for(i in 1:length(homol[["Homologue Series"]][,1])){
				those<-as.numeric(strsplit(homol[["Homologue Series"]][i,2],",")[[1]])
				#these<-match(those,peaklist2[,"peak_ID"]) # nonsense? remove?
				those<-those[order(peaklist4[those,1],decreasing=FALSE)] # by increasing mass!
				for(j in 2:length(those)){
					from<-(from+1)
					if(from>at){
						Homol_groups<-rbind(
							Homol_groups,
							matrix(nrow=atby,ncol=3,0)
						)
						at<-dim(Homol_groups)[1]
					}
					Homol_groups[from,1]<-those[j-1]
					Homol_groups[from,2]<-those[j]					
					Homol_groups[from,3]<-i	
				}
			}
			Homol_groups<-Homol_groups[1:from,]
			those<-(Homol_groups[,1]>Homol_groups[,2])
			if(any(those)){ # should not be envoked - peaklists are ordered by mass
				Homol_groups[those,c(1,2)]<-Homol_groups[those,c(2,1)]
			}
			Homol_groups<-Homol_groups[order(Homol_groups[,1],Homol_groups[,2],decreasing=FALSE),]
			save(Homol_groups,file=(file.path(logfile[[1]],"results","componentization","homologues",paste(for_file,sep="_"))))
			##########################################################################	
			# remove blind peaks - impute removed peaks ##############################			
			if(logfile$parameters$homol_blind=="TRUE"){ 
				those<-is.na(match(peaklist[,"peak_ID"],peaklist4[,"peak_ID"]))		
				if(any(those)){
					# impute (1) - "Peaks in homologue series"
					homol_left<-cbind(
						as.data.frame(peaklist[those,c("m/z_corr","int_corr","RT_corr","peak_ID")]),
						rep(0,sum(those)), 		# HS IDs
						rep(0,sum(those)), 		# series level
						rep(0,sum(those)), 		# to ID
						rep("none",sum(those)),	# m/z increment				
						rep("none",sum(those)),	# RT increment					
						rep(0,sum(those)),		# HS cluster	
						rep("",sum(those)),		# Targets
						rep("",sum(those))		# ISTDs						
					)
					names(homol_left)<-names(homol[["Peaks in homologue series"]])
					homol[["Peaks in homologue series"]]<-rbind(homol[["Peaks in homologue series"]],homol_left)
					ord_homol<-order(homol[["Peaks in homologue series"]][,"peak ID"])
					homol[["Peaks in homologue series"]]<-homol[["Peaks in homologue series"]][ord_homol,]
					# impute (2) - "Peaks per level"
					find_peak<-match(seq(1,dim(homol[["Peaks in homologue series"]])[1],1),ord_homol)
					for(n in 1:length(homol[["Peaks per level"]])){
						for(m in 1:length(homol[["Peaks per level"]][[n]])){	
							homol[["Peaks per level"]][[n]][[m]]<-
								find_peak[homol[["Peaks per level"]][[n]][[m]]]
						}
					}
				}
			}
			##########################################################################	
			# intersect with target screening results ################################
			if(
				(logfile$workflow[names(logfile$workflow)=="target_screen"]=="yes")	||
				(logfile$workflow[names(logfile$workflow)=="IS_screen"]=="yes")
			){
				cat("Annotating compounds ... ")
				for_file<-as.numeric(measurements[b,"ID"])
				for_mode<-measurements[b,"Mode"]
				######################################################################
				if(for_mode=="positive" & (any(intstand[,"ion_mode"] =="positive") || any(targets[,"ion_mode"] =="positive"))){
					if(exists("profileList_pos_copy",envir=as.environment(".GlobalEnv"))){rm("profileList_pos_copy",envir=as.environment(".GlobalEnv"))}	
					if(exists("profileList_pos_copy")){rm("profileList_pos_copy")}	
					if(exists("links_peaks_pos",envir=as.environment(".GlobalEnv"))){rm("links_peaks_pos",envir=as.environment(".GlobalEnv"))}	
					if(exists("links_peaks_pos")){rm("links_peaks_pos")}	
					load(file=file.path(as.character(logfile[[1]]),"results","profileList_pos_copy"),envir=as.environment(".GlobalEnv"))
					load(file=file.path(as.character(logfile[[1]]),"results","links_peaks_pos"),envir=as.environment(".GlobalEnv"))
					peaks<-profileList_pos_copy[["peaks"]][which(profileList_pos_copy[["peaks"]][,"sampleIDs"]==for_file),]		
					if(dim(peaks)[1]>0){
						peaks<-peaks[order(peaks[,"peakIDs"],decreasing=FALSE),]
						list_delmass<-list()
						for(i in 1:dim(homol[["Peaks in homologue series"]])[1]){
							if(homol[["Peaks in homologue series"]][i,"HS IDs"]!=0){
								# any screening matches ? #############################
								that<-match(homol[["Peaks in homologue series"]][i,"peak ID"],peaks[,"peakIDs"])
								if(is.na(that)){next}
								if(peaks[that,"links"]==0){next}
								got_tar<-unlist(links_peaks_pos[[peaks[that,"links"]]][[1]])
								got_IS<-unlist(links_peaks_pos[[peaks[that,"links"]]][[2]])
								if(!length(got_tar)&!length(got_IS)){next} # compounds in links entry?
								# check if matches belong to homologuous compounds ####
								if(length(got_tar)){
									keep_tar<-rep(TRUE,length(got_tar))
									for(j in 1:length(got_tar)){
# > BAUSTELLE
									}
									got_tar<-got_tar[keep_tar];
								}
								if(length(got_IS)){
									keep_IS<-rep(TRUE,length(got_IS))
									for(j in 1:length(got_IS)){
# > BAUSTELLE
									}
									got_IS<-got_IS[keep_IS]
								}
								if(!length(got_tar)&!length(got_IS)){next} # any homologuous compounds?					
								# collect series mass differences for each peak ########
								series<-as.numeric(strsplit(homol[["Peaks in homologue series"]][i,"HS IDs"],"/")[[1]])
								list_delmass[[i]]<-homol[["Homologue Series"]][series[1],"m/z increment"]
								for(j in series[-1]){
									if(homol[["Homologue Series"]][j,"HS IDs"]!=j){stop("\n Debug do_homologues at_1.")} # just a check
									list_delmass[[i]]<-c(list_delmass[[i]],homol[["Homologue Series"]][j,"m/z increment"])
								}
								list_delmass[[i]]<-unique(list_delmass[[i]])
								# check & add target compounds ##########################
								if(length(got_tar)){
									got_tar<-unique(got_tar) # wtf?
									for(j in 1:length(got_tar)){
										this_target<-strsplit(got_tar[j],"_")[[1]]
										at_target<-match(this_target[[1]],targets[,"ID"])
										if(!length(at_target)){stop("\n Debug do_homologues at_2.")} # just a check
										# -> calculate mass differences for chains ######
# > BAUSTELLE										
										# -> check if all chain masses were found #######
# > BAUSTELLE										
										# -> make entry #################################
										use_label<-paste0(c(this_target[[1]],targets[at_target,"Name"],this_target[[2]],this_target[[3]]),collapse="_")
										if(homol[["Peaks in homologue series"]][i,"Targets"]==""){ # first enty?
											homol[["Peaks in homologue series"]][i,"Targets"]<-use_label
										}else{
											homol[["Peaks in homologue series"]][i,"Targets"]<-paste0(c(homol[["Peaks in homologue series"]][i,"Targets"],use_label),collapse=", ")
										}
									}
								}
								# check & add ISTD compounds ###########################
								if(length(got_IS)){
									got_IS<-unique(got_IS)
									for(j in 1:length(got_IS)){
										this_IS<-strsplit(got_IS[j],"_")[[1]]
										at_IS<-match(this_IS[[1]],intstand[,"ID"])
										if(!length(at_IS)){stop("\n Debug do_homologues at_2.")} # just a check
										# -> calculate mass differences for chains ######
# > BAUSTELLE										
										# -> check if all chain masses were found #######
# > BAUSTELLE										
										# -> make entry #################################
										use_label<-paste0(c(this_IS[[1]],intstand[at_IS,"Name"],this_IS[[2]],this_IS[[3]]),collapse="_")
										if(homol[["Peaks in homologue series"]][i,"ISTDs"]==""){ # first enty?
											homol[["Peaks in homologue series"]][i,"ISTDs"]<-use_label
										}else{
											homol[["Peaks in homologue series"]][i,"ISTDs"]<-paste0(c(homol[["Peaks in homologue series"]][i,"Targets"],use_label),collapse=", ")
										}
									}
								}
								#########################################################
							}
						}
					}
				}
				######################################################################
				if(for_mode=="negative" & (any(intstand[,"ion_mode"] =="negative") || any(targets[,"ion_mode"] =="negative"))){
					if(exists("profileList_neg_copy",envir=as.environment(".GlobalEnv"))){rm("profileList_neg_copy",envir=as.environment(".GlobalEnv"))}	
					if(exists("profileList_neg_copy")){rm("profileList_neg_copy")}	
					if(exists("links_peaks_neg",envir=as.environment(".GlobalEnv"))){rm("links_peaks_neg",envir=as.environment(".GlobalEnv"))}	
					if(exists("links_peaks_neg")){rm("links_peaks_neg")}	
					load(file=file.path(as.character(logfile[[1]]),"results","profileList_neg_copy"),envir=as.environment(".GlobalEnv"))
					load(file=file.path(as.character(logfile[[1]]),"results","links_peaks_neg"),envir=as.environment(".GlobalEnv"))
					peaks<-profileList_neg_copy[["peaks"]][which(profileList_neg_copy[["peaks"]][,"sampleIDs"]==for_file),]		
					if(dim(peaks)[1]>0){
						peaks<-peaks[order(peaks[,"peakIDs"],decreasing=FALSE),]
						list_delmass<-list()
						for(i in 1:dim(homol[["Peaks in homologue series"]])[1]){
							if(homol[["Peaks in homologue series"]][i,"HS IDs"]!=0){
								# any screening matches ? #############################
								that<-match(homol[["Peaks in homologue series"]][i,"peak ID"],peaks[,"peakIDs"])
								if(is.na(that)){next}
								if(peaks[that,"links"]==0){next}
								got_tar<-unlist(links_peaks_neg[[peaks[that,"links"]]][[1]])
								got_IS<-unlist(links_peaks_neg[[peaks[that,"links"]]][[2]])
								if(!length(got_tar)&!length(got_IS)){next} # compounds in links entry?
								# check if matches belong to homologuous compounds ####
								if(length(got_tar)){
									keep_tar<-rep(TRUE,length(got_tar))
									for(j in 1:length(got_tar)){
# > BAUSTELLE
									}
									got_tar<-got_tar[keep_tar];
								}
								if(length(got_IS)){
									keep_IS<-rep(TRUE,length(got_IS))
									for(j in 1:length(got_IS)){
# > BAUSTELLE
									}
									got_IS<-got_IS[keep_IS]
								}
								if(!length(got_tar)&!length(got_IS)){next} # any homologuous compounds?					
								# collect series mass differences for each peak ########
								series<-as.numeric(strsplit(homol[["Peaks in homologue series"]][i,"HS IDs"],"/")[[1]])
								list_delmass[[i]]<-homol[["Homologue Series"]][series[1],"m/z increment"]
								for(j in series[-1]){
									if(homol[["Homologue Series"]][j,"HS IDs"]!=j){stop("\n Debug do_homologues at_1.")} # just a check
									list_delmass[[i]]<-c(list_delmass[[i]],homol[["Homologue Series"]][j,"m/z increment"])
								}
								list_delmass[[i]]<-unique(list_delmass[[i]])
								# check & add target compounds ##########################
								if(length(got_tar)){
									got_tar<-unique(got_tar) # wtf?
									for(j in 1:length(got_tar)){
										this_target<-strsplit(got_tar[j],"_")[[1]]
										at_target<-match(this_target[[1]],targets[,"ID"])
										if(!length(at_target)){stop("\n Debug do_homologues at_2.")} # just a check
										# -> calculate mass differences for chains ######
# > BAUSTELLE										
										# -> check if all chain masses were found #######
# > BAUSTELLE										
										# -> make entry #################################
#stop()
										use_label<-paste0(c(this_target[[1]],targets[at_target,"Name"],this_target[[2]],this_target[[3]]),collapse="_")
										if(homol[["Peaks in homologue series"]][i,"Targets"]==""){ # first enty?
											homol[["Peaks in homologue series"]][i,"Targets"]<-use_label
										}else{
											homol[["Peaks in homologue series"]][i,"Targets"]<-paste0(c(homol[["Peaks in homologue series"]][i,"Targets"],use_label),collapse=", ")
										}
									}
								}
								# check & add ISTD compounds ###########################
								if(length(got_IS)){
									got_IS<-unique(got_IS)
									for(j in 1:length(got_IS)){
										this_IS<-strsplit(got_IS[j],"_")[[1]]
										at_IS<-match(this_IS[[1]],intstand[,"ID"])
										if(!length(at_IS)){stop("\n Debug do_homologues at_2.")} # just a check
										# -> calculate mass differences for chains ######
# > BAUSTELLE										
										# -> check if all chain masses were found #######
# > BAUSTELLE										
										# -> make entry #################################
										use_label<-paste0(c(this_IS[[1]],intstand[at_IS,"Name"],this_IS[[2]],this_IS[[3]]),collapse="_")
										if(homol[["Peaks in homologue series"]][i,"ISTDs"]==""){ # first enty?
											homol[["Peaks in homologue series"]][i,"ISTDs"]<-use_label
										}else{
											homol[["Peaks in homologue series"]][i,"ISTDs"]<-paste0(c(homol[["Peaks in homologue series"]][i,"Targets"],use_label),collapse=", ")
										}
									}
								}
								#########################################################
							}
						}
					}
				}
				######################################################################	
			}
			##########################################################################				
			save(homol,file=(file.path(logfile[["project_folder"]],"results","componentization","homologues",paste("full",for_file,sep="_"))))			
			rm(peaklist,peaklist4,homol,Homol_groups)
			measurements[b,"homologues"]<-"TRUE"
			write.csv(measurements,file=file.path(logfile[["project_folder"]],"dataframes","measurements"),row.names=FALSE);
			cat("done.")
			##########################################################################		

		}else{
			cat("\n Homologues extracted before or file not included.")
		}
	}
	####################################################################################	
	
	####################################################################################	
	rm(mzfilter,elements,use_minmz,use_maxmz,measurements,use_mztol,intstand,targets)
	####################################################################################	

	
	
	
	