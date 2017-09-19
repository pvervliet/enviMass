#' @title 
#'
#' @description 
#'
#' @param 
#' 
#' @details enviMass workflow function
#' 

plot_ISTD_norm <-
function(
	int_norm_ISTD,
	xlim = FALSE,
	logfile,
	what
){

	if( is.na( match(what, c("normalization", "counts")) ) ) stop("\n Wrong what argument in plot_ISTD_norm!")
    ############################################################################
 	ylimit_del = c(0,0)
	count_IS = c(10000,0)
	count_b = c(10000,0)
	count_nb = c(10000,0)
	for( k in 1:length(int_norm_ISTD$lis_delint_IS) ){
		if(length(int_norm_ISTD$lis_delint_IS[[k]]) > 0){ # on IS
			if( min(int_norm_ISTD$lis_delint_IS[[k]]) < ylimit_del[1] ){
				ylimit_del[1] <- min(int_norm_ISTD$lis_delint_IS[[k]])
			}
			if( max(int_norm_ISTD$lis_delint_IS[[k]]) > ylimit_del[2] ){
				ylimit_del[2] <- max(int_norm_ISTD$lis_delint_IS[[k]])
			}
		}
		if(length(int_norm_ISTD$lis_delint_IS[[k]]) < count_IS[1]){
			count_IS[1] <- length(int_norm_ISTD$lis_delint_IS[[k]])
		}
		if(length(int_norm_ISTD$lis_delint_IS[[k]]) > count_IS[2]){
			count_IS[2] <- length(int_norm_ISTD$lis_delint_IS[[k]])
		}
		if( logfile$parameters$ISnorm_medblank == "TRUE" ){ # on blank
			if(length(int_norm_ISTD$lis_delint_b[[k]]) > 0){ 
				if( median(int_norm_ISTD$lis_delint_b[[k]]) < ylimit_del[1] ){
					ylimit_del[1] <- median(int_norm_ISTD$lis_delint_b[[k]])
				}
				if( median(int_norm_ISTD$lis_delint_b[[k]]) > ylimit_del[2] ){
					ylimit_del[2] <- median(int_norm_ISTD$lis_delint_b[[k]])
				}
			}
			if(length(int_norm_ISTD$lis_delint_b[[k]]) < count_b[1]){
				count_b[1] <- length(int_norm_ISTD$lis_delint_b[[k]])
			}
			if(length(int_norm_ISTD$lis_delint_b[[k]]) > count_b[2]){
				count_b[2] <- length(int_norm_ISTD$lis_delint_b[[k]])
			}
		}
		if( logfile$parameters$ISnorm_medsam == "TRUE" ){ # on non-blank
			if(length(int_norm_ISTD$lis_delint_nb[[k]]) > 0){ 
				if( median(int_norm_ISTD$lis_delint_nb[[k]]) < ylimit_del[1] ){
					ylimit_del[1] <- median(int_norm_ISTD$lis_delint_nb[[k]])
				}
				if( median(int_norm_ISTD$lis_delint_nb[[k]]) > ylimit_del[2] ){
					ylimit_del[2] <- median(int_norm_ISTD$lis_delint_nb[[k]])
				}
			}
			if(length(int_norm_ISTD$lis_delint_nb[[k]]) < count_nb[1]){
				count_nb[1] <- length(int_norm_ISTD$lis_delint_nb[[k]])
			}
			if(length(int_norm_ISTD$lis_delint_nb[[k]]) > count_nb[2]){
				count_nb[2] <- length(int_norm_ISTD$lis_delint_nb[[k]])
			}
		}					
	}				
    ############################################################################
	atPOSIX <- int_norm_ISTD$atPOSIX
	keep <- ((int_norm_ISTD$sampletype=="sample") | (int_norm_ISTD$sampletype=="blank"))
	atPOSIX <- atPOSIX[keep]
	sampleID <- int_norm_ISTD$sampleID[keep]
	timeset <- matrix(nrow=length(atPOSIX),ncol=5,0);
	for(i in 1:length(sampleID)){
		if(int_norm_ISTD$sampletype[i] == "sample"){
			timeset[i,2] <- as.numeric(int_norm_ISTD$sampleID[i]);
		}
		if(int_norm_ISTD$sampletype[i] == "blank"){
			timeset[i,3] <- as.numeric(int_norm_ISTD$sampleID[i]);
		}
	}	
	if(xlim[1] == FALSE){
		use_xlim <- c(-1, length(timeset[,1])+1)
	}else{
		use_xlim <- xlim
	}
    ############################################################################	
	if(what == "normalization"){

		########################################################################	
		plot.new()
		plot.window(xlim = use_xlim, ylim = c(ylimit_del[1] - 0.3, ylimit_del[2]))
		#axis(side = 1, at = seq(1, length(int_norm_ISTD$lis_delint_IS), 1), labels = sampleID, las = 2, cex.axis = 1)
		axis(side = 2, cex.axis = 1);
		#title(xlab="Temporal sequence of file IDs",ylab="Deviation from median log10 intensity",cex.lab=1.5)
		title(ylab = "Deviation from median log10 intensity", cex.lab = .9, cex = 1.2)
		abline(h = 0, col = "red")	
		for(k in 1:length(int_norm_ISTD$lis_delint_IS)){
			if(timeset[k,3] != 0){
				abline(v = k, col = "orange", lwd= 1)
			}else{
				abline(v = k, col = "lightgrey", lwd= .5)				
			}
		}
		plot.window(xlim = c(0, 10), ylim = c(0, 10))
		rect( # for legend; box can still be overplotted by points; but not the above gray lines
			xleft = -1, ybottom = -10, xright = 3.5, ytop = 2.,
			col = "white", border = NA
		)	
		plot.window(xlim = use_xlim, ylim = c(ylimit_del[1] - 0.3, ylimit_del[2]))
		for(k in 1:length(int_norm_ISTD$lis_delint_IS)){
			points( rep(k, length(int_norm_ISTD$lis_delint_IS[[k]])), int_norm_ISTD$lis_delint_IS[[k]], pch = 19, cex=  0.6, col = "darkgrey" )
		}		
		if( (logfile$parameters$ISnorm_medblank == "TRUE") ){	
			for(k in 1:length(int_norm_ISTD$lis_delint_b)){
				if(length(int_norm_ISTD$lis_delint_b[[k]]) > 0){
					points(k, median(int_norm_ISTD$lis_delint_b[[k]]), pch = 21, cex = .9, bg = "blue")
				}	
			}
		}
		if( (logfile$parameters$ISnorm_medsam=="TRUE") ){	
			for(k in 1:length(int_norm_ISTD$lis_delint_nb)){
				if(length(int_norm_ISTD$lis_delint_nb[[k]]) > 0){ 
					points(k, median(int_norm_ISTD$lis_delint_nb[[k]]), pch = 21, cex = .9, bg="green3")
				}
			}
		}
		for(k in 1:length(int_norm_ISTD$lis_delint_IS)){
			points( k, median(int_norm_ISTD$lis_delint_IS[[k]]), pch = 21, cex = .9, bg = "red")
		}			
		########################################################################
		plot.window(xlim = c(0,10), ylim = c(0,10))
		legend(-0.2, 2,
			pch = c(19, 21, 21, 21, 19),
			pt.cex = .9, cex = .9, bty = "n",
			legend = c(	"IS single profile deviation",
						"IS median deviation",
						"Blank median deviation",
						"Non-blank median deviation",
						"Blank file"),
			pt.bg = c("lightgrey", "red", "blue", "green3", "white"),
			col = c("lightgrey", "black", "black", "black", "white")
		)
		lines(x = c(-0.1, 0), y = c(0, 0), col = "orange")
	    if(xlim[1] != FALSE){
	      mtext("Zoomed in - click to zoom out partly or double-click to zoom out fully.", side = 3, line = 0.1, cex = .8, col = "darkgrey", at = 0, adj = 0)
	    }else{
	      mtext("Brush and doubleclick to zoom.", side = 3, line = 0.1, cex = .8, col = "darkgrey", at = 0, adj = 0)
	    }
		########################################################################
		plot.window(xlim = use_xlim, ylim = c(ylimit_del[1] - 0.3, ylimit_del[2]))	
		box();

	}
	############################################################################
	if(what == "counts"){	

		########################################################################
		plot.new()
		plot.window(xlim = use_xlim, ylim = c(count_IS[1]-1, max(c(count_IS[2]+1, as.numeric(logfile$parameters$ISnorm_numbIS)))))
		for(k in 1:length(int_norm_ISTD$lis_delint_IS)){
			if(timeset[k,3] != 0){
				abline(v = k, col = "orange", lwd = 1)
			}else{
				abline(v = k, col = "lightgrey", lwd = .5)				
			}
		}
		use_seq <- seq(1, length(int_norm_ISTD$lis_delint_IS), 1)
		use_seq <- use_seq[ (use_seq >= use_xlim[1]) & (use_seq <= use_xlim[2]) ]
		use_seq <- use_seq[seq(1, length(use_seq), ceiling(length(use_seq) / 35))]
		axis(side = 1 , at = use_seq, labels = sampleID[use_seq], las = 2, cex.axis = 1)
		axis(2, col = "blue", col.ticks = "red", col.axis = "red", cex.axis = .9);
		box()
		title(xlab = "Temporal sequence of file IDs", ylab = "", cex.lab = .9, cex = 1.2)
		countit <- c()	
		for(k in 1:length(int_norm_ISTD$lis_delint_IS)){
			countit <- c(countit, length(int_norm_ISTD$lis_delint_IS[[k]]))	
		}
		lines(countit, col = "red", lwd = 1.5)
		abline(h = as.numeric(logfile$parameters$ISnorm_numbIS), col = "red", lwd = 1, lty = 2)
		mtext("Number of IS peaks", side = 2, line = 2.5, col = "red", cex = 1)			
		if( logfile$parameters$ISnorm_medblank=="TRUE" ){	
		plot.window( xlim = use_xlim, ylim = c(count_b[1] - 1, count_b[2] + 1) )	
		countit <- c()	
		for( k in 1:length(int_norm_ISTD$lis_delint_IS) ){
				countit <- c(countit, length(int_norm_ISTD$lis_delint_b[[k]]))	
			}
			lines(countit, col = "blue", lwd = 1.5)
			axis(4, col="blue", col.ticks = "blue", col.axis = "blue", cex.axis = .9)
			mtext("Number of blank peaks", side = 4, line = 2.3, col = "blue", cex = 1)
		}
		if( logfile$parameters$ISnorm_medsam == "TRUE" ){	
			plot.window( xlim = use_xlim, ylim = c(count_nb[1] - 1, count_nb[2] + 1) )	
			countit <- c()	
			for( k in 1:length(int_norm_ISTD$lis_delint_IS) ){
				countit <- c(countit, length(int_norm_ISTD$lis_delint_nb[[k]]))	
			}
			lines(countit, col = "green3", lwd = 1.5)
			axis(4, col = "green3", col.ticks = "green3", col.axis = "green3", line = 4.5 , cex.axis = .9)
			mtext("Number of non-blank peaks", side = 4, line = 6.8, col = "green3", cex = 1)	
		}
		########################################################################
		box();

	}
	############################################################################
	return("done")
	
}




















