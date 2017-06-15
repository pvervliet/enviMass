

measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
IDs<-(measurements[,"ID"])
incl<-(measurements[,"include"])
filetypus<-(measurements[,"Type"])
ionmode<-(measurements[,"Mode"])
atdate<-(measurements[,"Date"])
attime<-(measurements[,"Time"])
attime2<-as.difftime(attime);
atdate<-as.Date(atdate);
sampleID<-(measurements[,"ID"])
old_samplewise<-(measurements[,"blind"])
new_samplewise<-old_samplewise
ord<-order(as.numeric(atdate),as.numeric(attime2),sampleID);

ppm<-logfile$parameters$blind_ppm
dmz<-as.numeric(logfile$parameters$blind_dmz)
dRT<-as.numeric(logfile$parameters$blind_drt)
int_ratio<-as.numeric(logfile$parameters$blind_threshold)

if(FALSE){ # debug parameters - ignore
	ppm<-TRUE
	dmz<-3
	dRT<-30
	int_ratio<-10
}

# clean old entries ####################################################################################
for(i in 1:length(IDs)){
	if(incl[i]=="FALSE"){next}	
	if(filetypus[i]=="blank"){next}
	if(old_samplewise[i]=="TRUE"){next}
	load(file=file.path(logfile[[1]],"peaklist",as.character(IDs[i])),envir=as.environment(".GlobalEnv"),verbose=FALSE);
	#keep_2<-rep(1,length(peaklist[,1])) # 1 == TRUE
	peaklist[,"keep_2"]<-1
	save(peaklist,file=file.path(logfile[[1]],"peaklist",as.character(IDs[i])))
	rm(peaklist)
}


#######################################################################################################
# run last blank by date & time subtraction ###########################################################
if((logfile$parameters$subtract_pos_bydate=="TRUE") || (logfile$parameters$subtract_neg_bydate=="TRUE")){
	blank_ID_last<-"FALSE"
	for(i in 2:length(ord)){ # can skip first file
		if((logfile$parameters$subtract_pos_bydate=="FALSE") & (ionmode[ord[i]]=="positive")){next}
		if((logfile$parameters$subtract_neg_bydate=="FALSE") & (ionmode[ord[i]]=="negative")){next}	
		if(old_samplewise[ord[i]]=="TRUE"){next}
		if(incl[ord[i]]=="FALSE"){next}		
		if(filetypus[ord[i]]=="sample"){
			sam_ID<-sampleID[ord[i]]
			found_blank<-FALSE
			for(j in (i-1):1){ # backward
				if((filetypus[ord[j]]=="blank") & (ionmode[ord[i]]==ionmode[ord[j]])){
					blank_ID<-sampleID[ord[j]]
					found_blank<-TRUE
					break;
				}
			}
			if(!found_blank){
				next;
			}
			if(blank_ID!=blank_ID_last){ # only load, not reload
				load(file=file.path(logfile[[1]],"peaklist",as.character(blank_ID)),verbose=FALSE);
				peaks_blank<-peaklist[,c("m/z_corr","int_corr","RT_corr")]
				rm(peaklist)
			}
			load(file=file.path(logfile[[1]],"peaklist",as.character(sam_ID)),verbose=FALSE);
			peaks_sample<-peaklist[,c("m/z_corr","int_corr","RT_corr")]		
			getit <- search_peak( 
				peaklist=peaks_blank, 
				mz=peaks_sample[,"m/z_corr"], 
				dmz=(dmz*2), 
				ppm=ppm, 
				RT=peaks_sample[,"RT_corr"], 
				dRT=dRT,
				onlymax=TRUE,
				int_ratio=int_ratio,
				int=peaks_sample[,"int_corr"],
				get_matches=FALSE
			)	
			peaklist[getit=="TRUE","keep_2"]<-0
			save(peaklist,file=file.path(logfile[[1]],"peaklist",as.character(sam_ID)))
			blank_ID_last<-blank_ID
			cat(paste("\n",
				round(sum(peaklist[,"keep_2"]==0)/length(peaklist[,1])*100,digits=1),
				" % of ",
				length(peaklist[,1]),
				" peaks blind filtered (files ",
				sam_ID," vs. ",blank_ID,", ",ionmode[ord[i]],", by date & time)."
			,sep=""))
			rm(peaklist);
			new_samplewise[ord[i]]<-"TRUE"
		}
	}
}
#######################################################################################################

#######################################################################################################
# run the selective subtraction #######################################################################
# POSITIVE ############################################################################################
if( (logfile$parameters$subtract_pos_byfile=="TRUE") & any(logfile$Positive_subtraction_files!="FALSE") ){
	selec_pos<-logfile$Positive_subtraction_files
	selec_pos<-selec_pos[selec_pos!="FALSE"]
	for(i in 1:length(IDs)){
		if(any(measurements[,"ID"]==IDs[i])){ # how not though?
			if( filetypus[measurements[,"ID"]==IDs[i]]=="sample" &  ionmode[measurements[,"ID"]==IDs[i]]=="positive" ){
				if(old_samplewise[measurements[,"ID"]==IDs[i]]=="TRUE"){next} # done before, samplewise
				if(incl[measurements[,"ID"]==IDs[i]]=="FALSE"){next}	
				load(file=file.path(logfile[[1]],"peaklist",as.character(IDs[i])),verbose=FALSE);
				sam_peaklist<-peaklist;rm(peaklist);
				for(j in 1:length(selec_pos)){
					subID<-strsplit(selec_pos[j]," - ")[[1]][1]
					load(file=file.path(logfile[[1]],"peaklist",as.character(subID)),verbose=FALSE)
					peaks_blank<-peaklist[,c("m/z_corr","int_corr","RT_corr")];rm(peaklist);
					getit <- search_peak(
						peaklist=peaks_blank, 
						mz=sam_peaklist[,"m/z_corr"], 
						dmz=(dmz*2), # precheck for profiles
						ppm=ppm, 
						RT=sam_peaklist[,"RT_corr"], 
						dRT=dRT,
						onlymax=TRUE,
						int_ratio=int_ratio,
						int=sam_peaklist[,"int_corr"],
						get_matches=FALSE
					)	
					sam_peaklist[getit=="TRUE","keep_2"]<-0
					rm(peaks_blank)
				}
				peaklist<-sam_peaklist
				save(peaklist,file=file.path(logfile[[1]],"peaklist",as.character(IDs[i])))
				cat(paste("\n",
					round(sum(peaklist[,colnames(peaklist)=="keep_2"]==0)/length(peaklist[,1])*100,digits=1),
					" % of ",
					length(peaklist[,1]),
					" peaks blind filtered (selective, file ",
					as.character(IDs[i]),"). "
				,sep=""))
				rm(peaklist,sam_peaklist);
				new_samplewise[measurements[,"ID"]==IDs[i]]<-"TRUE"
			}
		}else{
			cat("\n Orphaned peaklist detected - from an older workflow run?")
		}
	}
}
#######################################################################################################

#######################################################################################################
# run the selective subtraction #######################################################################
# NEGATIVE ############################################################################################
if( (logfile$parameters$subtract_neg_byfile=="TRUE") & any(logfile$Negative_subtraction_files!="FALSE") ){
	selec_neg<-logfile$Negative_subtraction_files
	selec_neg<-selec_neg[selec_neg!="FALSE"]
	for(i in 1:length(IDs)){
		if(any(measurements[,"ID"]==IDs[i])){	
			if(filetypus[measurements[,"ID"]==IDs[i]]=="sample" &  ionmode[measurements[,"ID"]==IDs[i]]=="negative"){
				if(old_samplewise[measurements[,"ID"]==IDs[i]]=="TRUE"){next} # done before, samplewise
				if(incl[measurements[,"ID"]==IDs[i]]=="FALSE"){next}
				load(file=file.path(logfile[[1]],"peaklist",as.character(IDs[i])),verbose=FALSE);
				sam_peaklist<-peaklist;rm(peaklist);
				for(j in 1:length(selec_neg)){
					subID<-strsplit(selec_neg[j]," - ")[[1]][1]
					load(file=file.path(logfile[[1]],"peaklist",as.character(subID)),verbose=FALSE)
					peaks_blank<-peaklist[,c("m/z_corr","int_corr","RT_corr")];rm(peaklist);
					getit <- search_peak( 
						peaklist=peaks_blank, 
						mz=sam_peaklist[,"m/z_corr"], 
						dmz=(dmz*2), # precheck for profiles
						ppm=ppm, 
						RT=sam_peaklist[,"RT_corr"], 
						dRT=dRT,
						onlymax=TRUE,
						int_ratio=int_ratio,
						int=sam_peaklist[,"int_corr"],
						get_matches=FALSE
					)	
					sam_peaklist[getit=="TRUE",colnames(sam_peaklist)=="keep_2"]<-0
					rm(peaks_blank)
				}
				peaklist<-sam_peaklist
				save(peaklist,file=file.path(logfile[[1]],"peaklist",as.character(IDs[i])))
				cat(paste("\n",
					round(sum(peaklist[,colnames(peaklist)=="keep_2"]==0)/length(peaklist[,1])*100,digits=1),
					" % of ",
					length(peaklist[,1]),
					" peaks blind filtered (selective, file ",
					as.character(IDs[i]),"). "
				,sep=""))
				rm(peaklist,sam_peaklist);
				new_samplewise[measurements[,"ID"]==IDs[i]]<-"TRUE"
			}
		}else{
			cat("\n Orphaned peaklist detected - from an older workflow run?")
		}		
	}
}
#######################################################################################################

#######################################################################################################
# update filewise switches for blind detection ########################################################
measurements[,"blind"]<-new_samplewise
write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);	
rm(measurements)
#######################################################################################################

