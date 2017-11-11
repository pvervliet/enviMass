#' @title Compare two lists for matching elements in the first list with at least one in the second list
#'
#' @export
#'
#' @description Compare two lists for equality of entries
#'
#' @param list1. First list
#' @param list2. Second list
#' @param as_pairs. Strictly pairwise matching?
#'
#' @return Logical whether all pairwise entries are equal
#' 
#' @details  enviMass workflow function.
#' 

LOD_wrap<-function(
	x,
	logfile
){

	for_ID <- x
	report <- "nothing"
	##############################################################################
	load(file.path(logfile$project_folder, "peaklist", for_ID), envir=as.environment(".GlobalEnv"))
	#peaklist<-peaklist[peaklist[,colnames(peaklist)=="keep"]==1,,drop=FALSE]
	if(length(peaklist[,1]) == 0){ return(report) }
	# LOD ###################################################################
	his <- hist(peaklist[,"RT"], breaks=100, plot=FALSE)
	get_int <- c()
	get_ret <- c()
	get_w <- c()
	for(j in 2:length(his$breaks)){
		ret <- peaklist[(peaklist[,"RT"] >= his$breaks[j-1] & peaklist[,"RT"] < his$breaks[j]),"RT"]
		int <- log10(peaklist[(peaklist[,"RT"] >= his$breaks[j-1] & peaklist[,"RT"] < his$breaks[j]),"int_corr"])
		ret <- ret[order(int,decreasing=FALSE)]
		int <- int[order(int,decreasing=FALSE)]
		getit <- ceiling(length(int)*0.1)
		get_int <- c(get_int,int[getit])
		get_ret <- c(get_ret,ret[getit])
		if(length(ret) > 0){get_w <- c(get_w,length(int))}
	}
	model <- smooth.spline(x = get_ret, y = get_int)	
	if(TRUE){
		png(file = file.path(logfile$project_folder, "results", "LOD", paste("plot_LOD_", for_ID, ".png", sep="")),
			width = 720, height = 250)
		par(mar=c(4.2, 4.2, 0.8, 0.8))
		plot(peaklist[,5], log10(peaklist[,13]), pch=19, cex=0.6, col="darkgrey", xlab="RT [s]",
			ylab = expression(log[10]*paste(" Intensity",sep=" "))
		)
		his <- hist(peaklist[,5], breaks=100, plot=FALSE)
		abline(v = his$breaks, col="grey")
		lines(get_ret, predict(model)$y, col="red", lwd=2)
		points(get_ret, get_int, col="black", pch=19, cex=0.5)
		box()
		dev.off();
	}
	##############################################################################
	return(model)
  
}

