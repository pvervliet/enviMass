
recalib_wrap<-function(
	x,
	logfile,
	mz,
	tolmz,
	ppm,
	ret,
	tolret,
	max_recal
){

	##############################################################################
	for_file <- x
	if(any(objects()=="peaklist")){rm(peaklist)}
	load( file = file.path(logfile[[1]], "peaklist", for_file), envir=environment());   
	if( length(mz)>0 ){
		peak_recal<-recalib(
			peaklist = peaklist[,c("m/z","max_int","RT")],
			mz = mz,
			tolmz = tolmz,
			ppm = ppm,
			ret = ret,
			tolret = tolret,
			what = "mass",
			one = TRUE,
			knot = 5,
			plotit = TRUE,
			path_1 = file.path(logfile[[1]],"pics",paste("recal_", for_file, sep="")),
			path_2 = file.path(logfile[[1]],"results","recalibration",paste("recal_gam_", for_file, sep="")),
			plot_ppm = c(2,5,10),
			max_recal = max_recal
		)
		if(length(peak_recal)>1){
			peaklist[,c(12,13,14)]<-peak_recal
			save(peaklist,file=file.path(logfile[[1]], "peaklist", for_file));         
		}else{
			peaklist[,"m/z_corr"] <- peaklist[,"m/z"]; # use dummy value!
			save(peaklist,file=file.path(logfile[[1]], "peaklist", for_file));         
			png(filename = file.path(logfile[[1]],"pics",paste("recal_", for_file, sep="")), bg = "white", width = 1100, height= 300)
			plot.new()
			plot.window(xlim=c(0,1),ylim=c(0,1))
			text(0.5,0.5,"No recalibration \n feasible.",cex=1)
			dev.off();
			return(paste("  Mass recalibration for file ID ", for_file, " - infeasible.\n", sep=""))        
		}
	}else{
		png(filename = file.path(logfile[[1]],"pics",paste("recal_", for_file, sep="")), bg = "white", width = 1100, height= 300)
		plot.new()
		plot.window(xlim=c(0,1),ylim=c(0,1))
		text(0.5,0.5,"No recalibration \n feasible.",cex=1)
		dev.off();	
		return(paste("  Mass recalibration for file ID ", for_file, " - infeasible: not enough positive recalibration masses available!\n", sep="")) 
	}
	##############################################################################
	return("done")
  
}
