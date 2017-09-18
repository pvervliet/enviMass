#' @title Initialize profile list object
#'
#' @export
#'
#' @description \code{startprofiles} initializes a list object containing all peaks from the files in an enviMass project, associated metadata and placeholder.
#'
#' @param logfile logfile object of an enviMass project.
#' @param sets Integer. Number of latest files to include.
#' @param ion_mode Character string, either "positive" or "negative".
#' @param until Integer, ID of file. All peaks of files up to the date of this file will be included.
#' @param selective Logical. Should only peaklist with measurements$profiled==TRUE be inluded?
#' @param types. File types to include in profiling, e.g., "sample", "blind", "calibration" or "spiked". For "spiked", all related files to subtract from are also included.
#' @param places. For which locations? Otherwise set to FALSE.
#' @param check_exist. Check (TRUE) if peaklist file exists for a ID - otherwise set to FALSE.
#'
#' @return profile list
#' 
#' @details  enviMass workflow function, parallel. Used sequentially: filter by type -> places -> check_exists -> until ID -> latest counts
#' 

get_measurement_IDs<-function(
	logfile,
	sets=FALSE,
	ion_mode="positive",
	until=FALSE,
	selective=FALSE,
	types=FALSE,
	places = FALSE,
	check_exist = TRUE
){

    ############################################################################
    # (0) check inputs #########################################################
	if((sets!=FALSE) & (!is.numeric(sets))){stop("\n sets must be FALSE or numeric; aborted.")}
	if(!is.logical(selective) & (selective!="FALSE" & selective!="TRUE")){stop("\n Argument selective must be logical")}
    ############################################################################
    # (1) read in data #########################################################
    measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
    measurements<-measurements[measurements[,names(measurements)=="include"]=="TRUE",]
	if(selective=="TRUE"){
		measurements<-measurements[measurements[,names(measurements)=="profiled"]=="TRUE",]
	}
	measurements<-measurements[measurements[,names(measurements)=="Mode"]==ion_mode,]	
    ############################################################################	
	# (2) adjust time period, sort #############################################
    dated<-measurements[,"Date"]
    timed<-measurements[,"Time"]
    datetime<-c()
    for(i in 1:length(timed)){
      datetime<-c(datetime,paste(dated[i],timed[i],"CET",sep=" "))
    }
	atPOSIX<-as.POSIXct(datetime);
	sampleID<-measurements[,"ID"];
	locus<-measurements[,"Place"];
	typus<-measurements[,"Type"];
	ord<-order(atPOSIX,decreasing=TRUE);
	atPOSIX<-atPOSIX[ord];
	sampleID<-sampleID[ord];
	locus<-locus[ord];
	typus<-typus[ord];
	datetime<-datetime[ord];
    ############################################################################	
	# (3) filter by types ######################################################
	if(types[1]!="FALSE"){
		remain<-rep(TRUE,length(sampleID))
		remain[is.na(match(typus,types))]<-FALSE
		if(any(types=="spiked") & any(measurements[,"Type"]=="spiked")){ # include subtraction files, too
			subtr_files<-measurements[
				match(sampleID[typus=="spiked"],measurements[,"ID"])
			,]$tag2
			remain[match(subtr_files,sampleID)]<-TRUE	
		}
		if(!any(remain)){stop("\n get_measurement_IDs: No remaining files to process/profile with these settings!")}
		sampleID<-sampleID[remain]
		locus<-locus[remain]		
		typus<-typus[remain]	
		atPOSIX<-atPOSIX[remain]
		datetime<-datetime[remain]
	}
    ############################################################################
    # (4) check if peaklist exists -> filter ##################################
    if(check_exist){
		remain<-rep(TRUE,length(sampleID))
		remain[!file.exists(file.path(logfile[[1]],"peaklist",as.character(sampleID)))]<-FALSE
		if(all(!remain)){stop("\n get_measurement_IDs: No remaining files to process/profile with these settings!")}
		sampleID<-sampleID[remain]
		locus<-locus[remain]		
		typus<-typus[remain]	
		atPOSIX<-atPOSIX[remain]	
		datetime<-datetime[remain]	
    }
    ############################################################################
    # (5) filter by places #####################################################
    if(places[1] != FALSE){
 		remain<-rep(TRUE,length(sampleID))
		remain[is.na(match(locus,places))]<-FALSE   	
		if(!any(remain)){stop("\n get_measurement_IDs: No remaining files to process/profile with these settings!")}
		sampleID<-sampleID[remain]	
		typus<-typus[remain]	
		atPOSIX<-atPOSIX[remain]
		datetime<-datetime[remain]
		locus<-locus[remain]
    }
    ############################################################################
	# (6) filter by until ID ###################################################
	if(until!="FALSE" & any(sampleID==until)){
		untilPOSIX<-atPOSIX[sampleID==until]
		from<-(1:length(sampleID))[atPOSIX<=untilPOSIX][1]
	}else{
		if(until!="FALSE"){warning("\n problem in startprofiles: until used in conjunction with an invalid sampleID. revise or debug?")}
		from<-1
	}
    ############################################################################
	# (7) filter by latest counts ##############################################
	if(sets!=FALSE){
		if((from+sets-1)>length(sampleID)){
			to<-length(sampleID)
		}else{
			to<-(from+sets-1)
		}
	}else{
		to<-length(sampleID)
	}
	sampleID<-sampleID[from:to];
	datetime<-datetime[from:to];
	locus<-locus[from:to];
	typus<-typus[from:to];
    ############################################################################
    for_samples<-list()
    for_samples[[1]]<-sampleID;
    for_samples[[2]]<-datetime;
    for_samples[[3]]<-locus;
	for_samples[[4]]<-typus;
	names(for_samples)<-c("sampleID","datetime","locus","typus")
    return(for_samples)

}

