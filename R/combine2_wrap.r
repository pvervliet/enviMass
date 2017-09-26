combine2_wrap <-
function(
	x,
	logfile,
	measurements,
	...
){
	
	for_file<-x
	for_mode<-measurements[measurements$ID==for_file,"Mode"]
	##########################################################################
	# get isotopologue grouping results ######################################
	if(
		do_isot & file.exists(file.path(logfile[[1]],"results","componentization","isotopologues",paste("full",for_file,sep="_")))
	){
		load(file.path(logfile[[1]],"results","componentization","isotopologues",paste("full",for_file,sep="_")))
	}else{
		pattern<-FALSE
	}
	##########################################################################
	# get adduct grouping results ############################################
	if(
		do_addu & file.exists(file.path(logfile[[1]],"results","componentization","adducts",paste("full",for_file,sep="_")))
	){
		load(file.path(logfile[[1]],"results","componentization","adducts",paste("full",for_file,sep="_")))
	}else{
		adduct<-FALSE
	}
	##########################################################################
	# get adduct grouping results ############################################
	if(
		do_homol & file.exists(file.path(logfile[[1]],"results","componentization","homologues",paste("full",for_file,sep="_")))
	){
		load(file.path(logfile[[1]],"results","componentization","homologues",paste("full",for_file,sep="_")))
	}else{
		homol<-FALSE
	}
	##########################################################################
	if(
		pattern[1]=="FALSE" & adduct[1]=="FALSE" & homol[1]=="FALSE" 
	){
		return(paste0("No adduct, isotopologue or HS results available for file with ID ", for_file))
	}
	##########################################################################
	# build components #######################################################
	component<-enviMass::combine2(
		pattern, 
		adduct, 
		homol, 
		rules = c(FALSE, FALSE, FALSE), 
		dont = FALSE
	)		
	component[[1]][,17]<-as.character(component[[1]][,17]) # please debug in nontarget!
	component[[1]][,18]<-as.character(component[[1]][,18]) # please debug in nontarget!
	named<-colnames(component[[1]]) # add 2 columns for target / ISTD screening intersection
	component[[1]]<-cbind(
		component[[1]],
		rep("-",dim(component[[1]])[1]),
		rep("-",dim(component[[1]])[1]),
		rep(0,dim(component[[1]])[1]),
		rep(0,dim(component[[1]])[1]),
		rep(0,dim(component[[1]])[1]),
		rep(0,dim(component[[1]])[1]),
		rep(0,dim(component[[1]])[1]),
		rep(0,dim(component[[1]])[1]),
		rep(0,dim(component[[1]])[1])
	)
	component[[1]][,19]<-as.character(component[[1]][,19]) # please debug in nontarget!
	component[[1]][,20]<-as.character(component[[1]][,20]) # please debug in nontarget!
	#component[[1]][,21]<-as.character(component[[1]][,21]) # please debug in nontarget!
	#component[[1]][,22]<-as.character(component[[1]][,22]) # please debug in nontarget!
	colnames(component[[1]])<-c(named,
		"Target peaks","ISTD peaks","Total peak number","Blind peak number",
		"Monois. peak ID |","Monois. m/z |","Monois. RT |","Monois. int. |","Monois. sample/blind int. ratio"
	)
	##########################################################################
	# extend component table (1) - mark TARGET peaks #########################
	# feasible via profileList_XXX_copy and links_peaks_pos
	if(logfile$workflow[names(logfile$workflow)=="target_screen"]=="yes"){
		######################################################################
		if(for_mode == "positive"){
			if(exists("profileList_pos_copy")){rm("profileList_pos_copy")}	
			if(exists("links_peaks_pos")){rm("links_peaks_pos")}	
			load(file = file.path(as.character(logfile[[1]]),"results","profileList_pos_copy"), envir = environment())
			load(file = file.path(as.character(logfile[[1]]),"results","links_peaks_pos"), envir = environment())
			for_peaks <- which(
				profileList_pos_copy[["peaks"]][,"sampleIDs"] == for_file
			)
			peaks <- profileList_pos_copy[["peaks"]][for_peaks,]
			rm(for_peaks)
			if(dim(peaks)[1] > 0){
				peaks<-peaks[order(peaks[,"peakIDs"],decreasing=FALSE),]
				for(i in 1:dim(component[["Components"]])[1]){
					collect_comp<-c()
					collect_peak<-c()
					##########################################################
					if(component[["Components"]][i,"ID pattern peaks |"] != "-"){
						those <- as.numeric(strsplit(component[["Components"]][i, "ID pattern peaks |"], ",")[[1]])
						these <- as.numeric(those)
						these2 <- these
						these <- match(these, peaks[,"peakIDs"])
						for(n in 1:length(these) ){
							if(peaks[these[n], "links"] != 0){
								if(these2[n] != peaks[these[n], "peakIDs"]){stop("Debug do_components.files.r - issue 2")} # TEST
								got <- unlist(links_peaks_pos[[peaks[these[n], "links"]]]$target)
								if(!length(got)) next
								got <- got[seq(1, length(got), 2)]
								collect_comp <- c(collect_comp, got)
								collect_peak <- c(collect_peak, rep(these[n], length(got)))	
							}
						}
					}	
					##########################################################
					if(component[["Components"]][i,"ID adduct peaks |"] != "-"){
						those <- as.numeric(strsplit(component[["Components"]][i,"ID adduct peaks |"],",")[[1]])
						these <- as.numeric(those)
						these2 <- these
						these <- match(these, peaks[,"peakIDs"])
						for(n in 1:length(these) ){
							if(peaks[these[n], "links"] != 0){
								if(these2[n] != peaks[these[n],"peakIDs"]){stop("Debug do_components.files.r - issue 3")} # TEST
								got <- unlist(links_peaks_pos[[peaks[these[n], "links"]]]$target)
								if(!length(got)) next
								got <- got[seq(1, length(got), 2)]
								collect_comp <- c(collect_comp, got)
								collect_peak <- c(collect_peak, rep(these[n], length(got)))	
							}
						}
					}	
					##########################################################
					if(component[["Components"]][i,"ID interfering peaks |"] != "-"){
						those <- as.numeric(strsplit(component[["Components"]][i, "ID interfering peaks |"], ",")[[1]])
						these <- as.numeric(those)
						these2 <- these
						these <- match(these,peaks[,"peakIDs"])
						for(n in 1:length(these) ){
							if(peaks[these[n],"links"] != 0){
								if(these2[n] != peaks[these[n],"peakIDs"]){stop("Debug do_components.files.r - issue 4")} # TEST
								got <- unlist(links_peaks_pos[[peaks[these[n], "links"]]]$target)
								if(!length(got)) next
								got <- got[seq(1, length(got), 2)]
								collect_comp <- c(collect_comp, got)
								collect_peak <- c(collect_peak, rep(these[n], length(got)))	
							}
						}
					}	
					##########################################################
					if(length(collect_comp)>0){ # found sth?
						got_comp<-unique(collect_comp)
						collect_all<-""
						for(n in 1:length(got_comp)){
							tar<-paste0(strsplit(got_comp[n],"_")[[1]][1:2],collapse="_")
							collect_all<-paste0(c(collect_all,
								paste(tar,
									paste0(unique(collect_peak[collect_comp==got_comp[n]]),collapse=","),
									sep=": "
								)),
								collapse=" / "
							)
						}
						if(nchar(collect_all)>3){
							collect_all<-substr(collect_all,4,nchar(collect_all))
						}
						component[[1]][i,"Target peaks"]<-collect_all
					}
				}
				rm(collect_comp,collect_peak)
			}	
			rm(profileList_pos_copy,links_peaks_pos,envir=as.environment(".GlobalEnv"))
		}	
		######################################################################
		if(for_mode=="negative"){
			if(exists("profileList_neg_copy",envir=as.environment(".GlobalEnv"))){rm("profileList_neg_copy",envir=as.environment(".GlobalEnv"))}	
			if(exists("profileList_neg_copy")){rm("profileList_neg_copy")}	
			if(exists("links_peaks_neg",envir=as.environment(".GlobalEnv"))){rm("links_peaks_neg",envir=as.environment(".GlobalEnv"))}	
			if(exists("links_peaks_neg")){rm("links_peaks_neg")}	
			load(file=file.path(as.character(logfile[[1]]),"results","profileList_neg_copy"), envir = environment())
			load(file=file.path(as.character(logfile[[1]]),"results","links_peaks_neg"), envir = environment())
			for_peaks<-which(
				profileList_neg_copy[["peaks"]][,"sampleIDs"]==for_file
			)
			peaks<-profileList_neg_copy[["peaks"]][for_peaks,]
			rm(for_peaks)
			if(dim(peaks)[1]>0){
				peaks<-peaks[order(peaks[,"peakIDs"],decreasing=FALSE),]
				for(i in 1:dim(component[[1]])[1]){
					collect_comp<-c()
					collect_peak<-c()
					##########################################################
					if(component[[1]][i,3]!="-"){
						those<-as.numeric(strsplit(component[[1]][i,3],",")[[1]])
						these<-as.numeric(those)
						these2<-these
						these<-match(these,peaks[,"peakIDs"])
						for(n in 1:length(these) ){
							if(peaks[these[n],"links"]!=0){
								if(these2[n]!=peaks[these[n],"peakIDs"]){stop("Debug do_components.files.r - issue 2")} # TEST
								got<-unlist(links_peaks_neg[[peaks[these[n],"links"]]]$target)
								if(!length(got)) next
								got <- got[seq(1,length(got),2)]
								collect_comp<-c(collect_comp,got)
								collect_peak<-c(collect_peak,rep(these[n],length(got)))	
							}
						}
					}	
					##########################################################
					if(component[[1]][i,5]!="-"){
						those<-as.numeric(strsplit(component[[1]][i,5],",")[[1]])
						these<-as.numeric(those)
						these2<-these
						these<-match(these,peaks[,"peakIDs"])
						for(n in 1:length(these) ){
							if(peaks[these[n],"links"]!=0){
								if(these2[n]!=peaks[these[n],"peakIDs"]){stop("Debug do_components.files.r - issue 3")} # TEST
								got<-unlist(links_peaks_neg[[peaks[these[n],"links"]]]$target)
								if(!length(got)) next
								got <- got[seq(1,length(got),2)]
								collect_comp<-c(collect_comp,got)
								collect_peak<-c(collect_peak,rep(these[n],length(got)))	
							}
						}
					}	
					##########################################################
					if(component[[1]][i,7]!="-"){
						those<-as.numeric(strsplit(component[[1]][i,7],",")[[1]])
						these<-as.numeric(those)
						these2<-these
						these<-match(these,peaks[,"peakIDs"])
						for(n in 1:length(these) ){
							if(peaks[these[n],"links"]!=0){
								if(these2[n]!=peaks[these[n],"peakIDs"]){stop("Debug do_components.files.r - issue 3")} # TEST
								got<-unlist(links_peaks_neg[[peaks[these[n],"links"]]]$target)
								if(!length(got)) next
								got <- got[seq(1,length(got),2)]
								collect_comp<-c(collect_comp,got)
								collect_peak<-c(collect_peak,rep(these[n],length(got)))	
							}
						}
					}	
					##########################################################
					if(length(collect_comp)>0){ # found sth?
						got_comp<-unique(collect_comp)
						collect_all<-""
						for(n in 1:length(got_comp)){
							tar<-paste0(strsplit(got_comp[n],"_")[[1]][1:2],collapse="_")
							collect_all<-paste0(c(collect_all,
								paste(tar,
									paste0(unique(collect_peak[collect_comp==got_comp[n]]),collapse=","),
								sep=": "
							)),
							collapse=" / "
							)
						}
						if(nchar(collect_all)>3){
							collect_all<-substr(collect_all,4,nchar(collect_all))
						}
						component[[1]][i,"Target peaks"]<-collect_all
					}
				}
				rm(collect_comp,collect_peak)
			}	
			rm(profileList_neg_copy,links_peaks_neg,envir=as.environment(".GlobalEnv"))
		}	
		######################################################################										
	}
	##########################################################################
	# extend component table (2) - mark ISTD peaks ###########################
	if(logfile$workflow[names(logfile$workflow)=="IS_screen"]=="yes"){
		######################################################################
		if(for_mode=="positive"){
			if(exists("profileList_pos_copy")){rm("profileList_pos_copy")}	
			if(exists("links_peaks_pos")){rm("links_peaks_pos")}	
			load(file=file.path(as.character(logfile[[1]]),"results","profileList_pos_copy"), envir= environment())
			load(file=file.path(as.character(logfile[[1]]),"results","links_peaks_pos"), envir= environment())
			for_peaks<-which(
				profileList_pos_copy[["peaks"]][,"sampleIDs"]==for_file
			)
			peaks<-profileList_pos_copy[["peaks"]][for_peaks,]
			rm(for_peaks)
			if(dim(peaks)[1]>0){
				peaks<-peaks[order(peaks[,"peakIDs"],decreasing=FALSE),]
				for(i in 1:dim(component[[1]])[1]){
					collect_comp<-c()
					collect_peak<-c()
					##########################################################
					if(component[[1]][i,3]!="-"){
						those<-as.numeric(strsplit(component[[1]][i,3],",")[[1]])
						these<-as.numeric(those)
						these2<-these
						these<-match(these,peaks[,"peakIDs"])
						for(n in 1:length(these) ){
							if(peaks[these[n],"links"]!=0){
								if(these2[n]!=peaks[these[n],"peakIDs"]){stop("Debug do_components.files.r - issue 2")} # TEST
								got<-unlist(links_peaks_pos[[peaks[these[n],"links"]]]$IS)
								if(!length(got)) next
								got <- got[seq(1,length(got),2)]
								collect_comp<-c(collect_comp,got)
								collect_peak<-c(collect_peak,rep(these[n],length(got)))	
							}
						}
					}	
					##########################################################
					if(component[[1]][i,5]!="-"){
						those<-as.numeric(strsplit(component[[1]][i,5],",")[[1]])
						these<-as.numeric(those)
						these2<-these
						these<-match(these,peaks[,"peakIDs"])
						for(n in 1:length(these) ){
							if(peaks[these[n],"links"]!=0){
								if(these2[n]!=peaks[these[n],"peakIDs"]){stop("Debug do_components.files.r - issue 3")} # TEST
								got<-unlist(links_peaks_pos[[peaks[these[n],"links"]]]$IS)
								if(!length(got)) next
								got <- got[seq(1,length(got),2)]
								collect_comp<-c(collect_comp,got)
								collect_peak<-c(collect_peak,rep(these[n],length(got)))	
							}
						}
					}	
					##########################################################
					if(component[[1]][i,7]!="-"){
						those<-as.numeric(strsplit(component[[1]][i,7],",")[[1]])
						these<-as.numeric(those)
						these2<-these
						these<-match(these,peaks[,"peakIDs"])
						for(n in 1:length(these) ){
							if(peaks[these[n],"links"]!=0){
								if(these2[n]!=peaks[these[n],"peakIDs"]){stop("Debug do_components.files.r - issue 3")} # TEST
								got<-unlist(links_peaks_pos[[peaks[these[n],"links"]]]$IS)
								if(!length(got)) next
								got <- got[seq(1,length(got),2)]
								collect_comp<-c(collect_comp,got)
								collect_peak<-c(collect_peak,rep(these[n],length(got)))	
							}
						}
					}	
					##########################################################
					if(length(collect_comp)>0){ # found sth?
						got_comp<-unique(collect_comp)
						collect_all<-""
						for(n in 1:length(got_comp)){
							tar<-paste0(strsplit(got_comp[n],"_")[[1]][1:2],collapse="_")
							collect_all<-paste0(c(collect_all,
								paste(tar,
									paste0(unique(collect_peak[collect_comp==got_comp[n]]),collapse=","),
									sep=": "
								)),
								collapse=" / "
							)
						}
						if(nchar(collect_all)>3){
							collect_all<-substr(collect_all,4,nchar(collect_all))
						}
						component[[1]][i,"ISTD peaks"]<-collect_all
					}
				}
				rm(collect_comp,collect_peak)
			}	
			rm(profileList_pos_copy,links_peaks_pos,envir=as.environment(".GlobalEnv"))
		}
		######################################################################
		if(for_mode=="negative"){
			if(exists("profileList_neg_copy")){rm(profileList_neg_copy)}		
			if(exists("links_peaks_neg")){rm(links_peaks_neg)}	
			load(file=file.path(as.character(logfile[[1]]),"results","profileList_neg_copy"), envir = environment())
			load(file=file.path(as.character(logfile[[1]]),"results","links_peaks_neg"), envir = environment())
			for_peaks<-which(
				profileList_neg_copy[["peaks"]][,"sampleIDs"]==for_file
			)
			peaks<-profileList_neg_copy[["peaks"]][for_peaks,]
			rm(for_peaks)
			if(dim(peaks)[1]>0){
				peaks<-peaks[order(peaks[,"peakIDs"],decreasing=FALSE),]
				for(i in 1:dim(component[[1]])[1]){
					collect_comp<-c()
					collect_peak<-c()
					##########################################################
					if(component[[1]][i,3]!="-"){
						those<-as.numeric(strsplit(component[[1]][i,3],",")[[1]])
						these<-as.numeric(those)
						these2<-these
						these<-match(these,peaks[,"peakIDs"])
						for(n in 1:length(these) ){
							if(peaks[these[n],"links"]!=0){
							if(these2[n]!=peaks[these[n],"peakIDs"]){stop("Debug do_components.files.r - issue 2")} # TEST
							got<-unlist(links_peaks_neg[[peaks[these[n],"links"]]]$IS)
							if(!length(got)) next
							got <- got[seq(1,length(got),2)]
							collect_comp<-c(collect_comp,got)
							collect_peak<-c(collect_peak,rep(these[n],length(got)))	
							}
						}
					}	
					##########################################################
					if(component[[1]][i,5]!="-"){
						those<-as.numeric(strsplit(component[[1]][i,5],",")[[1]])
						these<-as.numeric(those)
						these2<-these
						these<-match(these,peaks[,"peakIDs"])
						for(n in 1:length(these) ){
							if(peaks[these[n],"links"]!=0){
								if(these2[n]!=peaks[these[n],"peakIDs"]){stop("Debug do_components.files.r - issue 3")} # TEST
								got<-unlist(links_peaks_neg[[peaks[these[n],"links"]]]$IS)
								if(!length(got)) next
								got <- got[seq(1,length(got),2)]
								collect_comp<-c(collect_comp,got)
								collect_peak<-c(collect_peak,rep(these[n],length(got)))	
							}
						}
					}	
					##########################################################
					if(component[[1]][i,7]!="-"){
						those<-as.numeric(strsplit(component[[1]][i,7],",")[[1]])
						these<-as.numeric(those)
						these2<-these
						these<-match(these,peaks[,"peakIDs"])
						for(n in 1:length(these) ){
							if(peaks[these[n],"links"]!=0){
								if(these2[n]!=peaks[these[n],"peakIDs"]){stop("Debug do_components.files.r - issue 3")} # TEST
								got<-unlist(links_peaks_neg[[peaks[these[n],"links"]]]$IS)
								if(!length(got)) next
								got <- got[seq(1,length(got),2)]
								collect_comp<-c(collect_comp,got)
								collect_peak<-c(collect_peak,rep(these[n],length(got)))	
							}
						}
					}	
					##########################################################
					if(length(collect_comp)>0){ # found sth?
						got_comp<-unique(collect_comp)
						collect_all<-""
						for(n in 1:length(got_comp)){
							tar<-paste0(strsplit(got_comp[n],"_")[[1]][1:2],collapse="_")
							collect_all<-paste0(c(collect_all,
								paste(tar,
									paste0(unique(collect_peak[collect_comp==got_comp[n]]),collapse=","),
									sep=": "
								)),
								collapse=" / "
							)
						}
						if(nchar(collect_all)>3){
							collect_all<-substr(collect_all,4,nchar(collect_all))
						}
						component[[1]][i,"ISTD peaks"]<-collect_all
					}
				}
				rm(collect_comp,collect_peak)
			}	
			rm(profileList_neg_copy,links_peaks_neg,envir=as.environment(".GlobalEnv"))	
		}	
		######################################################################		
	}
	##########################################################################
	# extend component table (3) - mark blind peaks by asterisk ##############				
	if(logfile$workflow[names(logfile$workflow) == "blind"] == "yes"){	
		if(exists("peaklist")){rm("peaklist")}	
		load(file = file.path(logfile[[1]], "peaklist", as.character(for_file)), envir = environment()); # Peaklist  
		peaklist <- peaklist[order(peaklist[,10], decreasing = FALSE),] # match with IDs - for saving pattern; IDs are retrieved for pairs seperately
		for(i in 1:dim(component[[1]])[1]){
			len_tot <- 0
			len_blind <- 0
			if(component[["Components"]][i, "ID pattern peaks |"] != "-"){
				those <- as.numeric(strsplit(component[[1]][i, "ID pattern peaks |"], ",")[[1]])
				these <- as.numeric(those)
				these2 <- these
				these <- match(these, peaklist[,"peak_ID"])
				len_tot <- (len_tot + length(these))
				for(j in 1:length(these)){
					if(these2[j] != peaklist[these[j], "peak_ID"]){stop("Debug do_components.files.r - issue 5")} # TEST
					if(peaklist[these[j],"keep_2"] < as.numeric(logfile$parameters$blind_threshold)){
						those[j] <- paste0(those[j],"*")
						len_blind <- (len_blind+1)
					}
				}
				component[["Components"]][i,"ID pattern peaks |"] <- paste0(those,collapse=",")
			}
			if(component[["Components"]][i,"ID adduct peaks |"] != "-"){
				those <- as.numeric(strsplit(component[[1]][i,"ID adduct peaks |"], ",")[[1]])
				these <- as.numeric(those)
				these2 <- these
				these <- match(these, peaklist[,"peak_ID"])
				len_tot <- (len_tot + length(these))
				for(j in 1:length(these)){
					if(these2[j] != peaklist[these[j], "peak_ID"]){stop("Debug do_components.files.r - issue 6")} # TEST
					if(peaklist[these[j],"keep_2"] < as.numeric(logfile$parameters$blind_threshold)){
						those[j] <- paste0(those[j],"*")
						len_blind <- (len_blind+1)
					}
				}
				component[["Components"]][i, "ID adduct peaks |"] <- paste0(those,collapse=",")
			}
			if(component[["Components"]][i,"ID interfering peaks |"]!="-"){
				those<-as.numeric(strsplit(component[[1]][i,"ID interfering peaks |"],",")[[1]])
				these<-as.numeric(those)
				these2<-these
				these<-match(these,peaklist[,"peak_ID"])
				len_tot<-(len_tot+length(these))
				for(j in 1:length(these)){
					if(these2[j]!=peaklist[these[j],"peak_ID"]){stop("Debug do_components.files.r - issue 7")} # TEST
					if(peaklist[these[j],"keep_2"]<as.numeric(logfile$parameters$blind_threshold)){
						those[j]<-paste0(those[j],"*")
						len_blind<-(len_blind+1)
					}
				}
				component[["Components"]][i,"ID interfering peaks |"]<-paste0(those,collapse=",")
			}
			component[["Components"]][i,"Total peak number"]<-as.character(len_tot)
			component[["Components"]][i,"Blind peak number"]<-as.character(len_blind)
		}
		rm(peaklist,envir=as.environment(".GlobalEnv"))
	}
	##########################################################################
	save(component,file=(file.path(logfile[[1]],"results","componentization","components",paste(for_file))))
	##########################################################################	
	return(NULL);

}
