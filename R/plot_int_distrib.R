
plot_int_distrib <-
function(
	int_distrib,
	ion_mode = "positive",
	types = FALSE,
	measurements,
	xlim = FALSE,
	what = "boxplot",
	maxit_1_lim = FALSE,
	maxit_2_lim = FALSE
){

	those_files <- which(measurements[, "Mode"] == ion_mode)
	if(!length(those_files)) return("nothing")
    ############################################################################
	# available in int_distrib ? ###############################################
	keep_IDs <- rep(TRUE, length(those_files))
	for(i in 1:length(those_files)){
		if(length(int_distrib[[as.numeric(measurements[those_files[i],]$ID)]]) < 2){
			keep_IDs[i] <- FALSE
		}
	}
	those_files <- those_files[keep_IDs]
	if(!length(those_files)) return("nothing")
    ############################################################################	
	# only specific types? #####################################################
	if(types[1] != FALSE){
		those_files <- those_files[!is.na(match(measurements[those_files,]$Type, types))]
	}
	if(!length(those_files)) return("nothing")
    ############################################################################			
    dated <- measurements[those_files, "Date"]
    timed <- measurements[those_files, "Time"]
	type <- measurements[those_files, "Type"]
    datetime <- c()
    for(i in 1:length(timed)){
      datetime <- c(datetime, paste(dated[i], timed[i], "CET", sep = " "))
    }
	atPOSIX <- as.POSIXct(datetime)
	ord <- order(atPOSIX, type,decreasing=TRUE)
	those_files <- those_files[ord]
	type <- type[ord]
	atPOSIX <- atPOSIX[ord]
	datetime <- datetime[ord]
    ############################################################################	
	if(what == "boxplot"){
		
		########################################################################	
		# restrict to xlim #####################################################
		file_seq <- seq(1,length(those_files),1)
		if(xlim[1] != FALSE){	
			file_seq  <- file_seq[
				(file_seq >= xlim[1]) & 
				(file_seq <= xlim[2]) 		
			] 
			use_xlim <- c(xlim[1], xlim[2])
		}else{
			use_xlim <- c(-1, length(those_files)+1)	
		}
		if(!length(those_files)) return("nothing")	
		########################################################################	
		ylim <- c(Inf, -Inf)
		peak_counts <- rep(0, length(file_seq))
		for(i in 1:length(file_seq)){
			at_entry <- as.numeric(measurements[those_files[file_seq[i]],]$ID)
			range_value <- range( c(int_distrib[[at_entry]]$boxplot$stats, int_distrib[[at_entry]]$boxplot$out))
			if(range_value[1] < ylim[1]) {ylim[1] <- range_value[1]}
			if(range_value[2] > ylim[2]) {ylim[2] <- range_value[2]}	
			if(range_value[1] < ylim[1]) {ylim[1] <- range_value[1]}
			if(range_value[2] > ylim[2]) {ylim[2] <- range_value[2]}	
			peak_counts[i] <- int_distrib[[at_entry]]$peak_count  
		}		
		plot.new()
		# n > 20 ->  Overview plot 
		plot.window(xlim = use_xlim, ylim = ylim)
		axis(2)
		for(i in 1:length(file_seq)){		
			at_entry <- as.numeric(measurements[those_files[file_seq[i]],]$ID)
			use_col <- "black" # default
			if( type[file_seq[i]] == "blank"){ use_col <- "green" }			
			if( type[file_seq[i]] == "calibration"){ use_col <- "orange" }				
			if( type[file_seq[i]] == "spiked"){ use_col <- "darkblue" }				
			if(length(file_seq) <= 20){	
				rect(
					xleft = (file_seq[i] - .3),
					ybottom = int_distrib[[at_entry]]$boxplot$stats[[2]],
					xright = (file_seq[i] + .3),
					ytop = int_distrib[[at_entry]]$boxplot$stats[[4]],
					border = use_col
				)
				lines(
					x = c((file_seq[i] - .3), (file_seq[i] + .3)),
					y = c(int_distrib[[at_entry]]$boxplot$stats[[3]], int_distrib[[at_entry]]$boxplot$stats[[3]]),
					lwd = 2,
					col = use_col
				)		
			}else{ # just a dot
				points(
					file_seq[i],
					int_distrib[[at_entry]]$boxplot$stats[[3]],
					pch = 19, cex = .7, col = use_col
				)
			}			
			points(
				rep(file_seq[i], length(int_distrib[[at_entry]]$boxplot$out)),
				int_distrib[[at_entry]]$boxplot$out,
				pch = 19, cex = .3, col = "lightgray"
			)
			lines(
				x = c(file_seq[i],file_seq[i]),
				y = c(int_distrib[[at_entry]]$boxplot$stats[[1]], int_distrib[[at_entry]]$boxplot$stats[[2]]),
				lwd = 2,
				col = use_col
			)
			lines(
				x = c(file_seq[i],file_seq[i]),
				y = c(int_distrib[[at_entry]]$boxplot$stats[[4]], int_distrib[[at_entry]]$boxplot$stats[[5]]),
				lwd = 2,
				col = use_col
			)
		}
		########################################################################
		use_seq <- seq(1, length(file_seq), 1)
		use_seq <- use_seq[seq(1, length(use_seq), ceiling(length(use_seq) / 35))]
		axis(side = 1 , at = file_seq[use_seq], labels = measurements[those_files[file_seq],]$ID[use_seq], las = 2, cex.axis = 1)	
		title(xlab = "Temporal sequence of file IDs", ylab = "log10 intensity", cex.lab = .9, cex = 1.2)
		mtext("Number of peaks", side = 4, line = 2.5, col = "red", cex = 1)		
		########################################################################	
		plot.window(xlim = use_xlim, ylim = c(0, max(peak_counts)))
		axis(4, col = "red", col.ticks = "red", col.axis = "red")
		points(
			file_seq,
			peak_counts,
			col = "red",
			pch = 4, cex = .8
		)
		########################################################################
		plot.window(xlim = c(0,10), ylim = c(0,10))
	    if(xlim[1] != FALSE){
			mtext("Zoomed in - click to zoom out partly or double-click to zoom out fully.", side = 3, line = 0.1, cex = .8, col = "darkgrey", at = 0, adj = 0)
	    }else{
			mtext("Brush and doubleclick to zoom.", side = 3, line = 0.1, cex = .8, col = "darkgrey", at = 0, adj = 0)
	    }
		mtext("Files types:", side = 1, line = 4.6, col = "black", cex = .8, at = 0, adj = 0)
		mtext("sample", side = 1, line = 4.6, col = "black", cex = .8, at = 1.5, adj = 0)
		mtext("blank/blind", side = 1, line = 4.6, col = "green", cex = .8, at = 2.5, adj = 0)		
		mtext("calibration", side = 1, line = 4.6, col = "orange", cex = .8, at = 4, adj = 0)		
		mtext("spiked", side = 1, line = 4.6, col = "darkblue", cex = .8, at = 5.5, adj = 0)	
		plot.window(xlim = use_xlim, ylim = c(0, max(peak_counts)))		# reset -> for proper zooming
		########################################################################
	
	}
    ############################################################################
	
    ############################################################################	
	if(what == "quantiles_distrib"){
		
		########################################################################	
		iles <- seq(0, 1, .01)
		mean_quantiles <- rep(0, length(iles))
		file_seq <- seq(1, length(those_files), 1)
		quantiles_val <- c(Inf, -Inf)
		for(i in 1:length(those_files)){
			at_entry <- as.numeric(measurements[those_files[i],]$ID)
			range_iles <- range(int_distrib[[at_entry]]$quantiles)
			if(range_iles[1] < quantiles_val[1]) {quantiles_val[1] <- range_iles[1]}
			if(range_iles[2] > quantiles_val[2]) {quantiles_val[2] <- range_iles[2]}			
			mean_quantiles <- (mean_quantiles + int_distrib[[at_entry]]$quantiles)
		}
		mean_quantiles <- (mean_quantiles / (length(those_files)) )
		########################################################################
		if( (maxit_1_lim[1] != FALSE) || (maxit_2_lim[1] != FALSE) ){
			####################################################################			
			maxit_1 <- rep(0, length(those_files))
			maxit_2 <- rep(0, length(those_files))
			for(i in 1:length(those_files)){		
				at_entry <- as.numeric(measurements[those_files[i],]$ID)
				maxit_1[i] <- max(abs(int_distrib[[at_entry]]$quantiles - mean_quantiles))
				maxit_2[i] <- median(int_distrib[[at_entry]]$quantiles - mean_quantiles)	
			}
			use_col <- rep("black", length(those_files))
			use_lwd <- rep(2, length(those_files))
			####################################################################
			if(maxit_1_lim[1] != FALSE){
				use_col[
					(maxit_1 < maxit_1_lim[1]) |
					(maxit_1 > maxit_1_lim[2])
				] <- "lightgray"
				use_lwd[
					(maxit_1 < maxit_1_lim[1]) |
					(maxit_1 > maxit_1_lim[2])
				] <- 1
			}
			####################################################################
			if(maxit_2_lim[1] != FALSE){
				use_col[
					(maxit_2 < maxit_2_lim[1]) |
					(maxit_2 > maxit_2_lim[2])
				] <- "lightgray"
				use_lwd[
					(maxit_2 < maxit_2_lim[1]) |
					(maxit_2 > maxit_2_lim[2])
				] <- 1
			}
			####################################################################					
		}else{
			use_col <- rep("lightgray", length(those_files))
			use_lwd <- rep(1, length(those_files))
		}
		########################################################################
		plot.new()
		plot.window(xlim = quantiles_val, ylim = c(0, 1))
		title(ylab = "Quantiles", xlab = "log10 intensity")
		for(i in 1:length(those_files)){			
			at_entry <- as.numeric(measurements[those_files[i],]$ID)
			lines(
				int_distrib[[at_entry]]$quantiles,
				iles,
				col = use_col[i],
				lwd = use_lwd[i]
			)
		}
		lines(
			mean_quantiles,
			iles,
			col = "red", lwd = 2
		)
		box();axis(1);axis(2);
		########################################################################	
	
	}
    ############################################################################

    ############################################################################	
	if(what == "quantiles_out"){

		########################################################################
		iles <- seq(0, 1, .01)
		mean_quantiles <- rep(0, length(iles))
		maxit_1 <- rep(0, length(those_files))
		maxit_2 <- rep(0, length(those_files))
		maxit_ID <- rep(0, length(those_files))
		maxit_col <- rep("", length(those_files))
		for(i in 1:length(those_files)){
			at_entry <- as.numeric(measurements[those_files[i],]$ID)			
			mean_quantiles <- (mean_quantiles + int_distrib[[at_entry]]$quantiles)
		}
		mean_quantiles <- (mean_quantiles / (length(those_files)) )
		for(i in 1:length(those_files)){		
			maxit_ID[i] <- measurements[those_files[i],]$ID
			at_entry <- as.numeric(maxit_ID[i])
			maxit_1[i] <- max(abs(int_distrib[[at_entry]]$quantiles - mean_quantiles))
			maxit_2[i] <- median(int_distrib[[at_entry]]$quantiles - mean_quantiles)	
			if(measurements[those_files[i], "Type"] == "sample"){ maxit_col[i] <- "black" }
			if(measurements[those_files[i], "Type"] == "blank"){ maxit_col[i] <- "green" }				
			if(measurements[those_files[i], "Type"] == "calibration"){ maxit_col[i] <- "orange" }				
			if(measurements[those_files[i], "Type"] == "spiked"){ maxit_col[i] <- "darkblue" }	
		}
		########################################################################	
		if( maxit_1_lim[1] != FALSE ){
			use_maxit1 <- maxit_1_lim
		}else{
			use_maxit1 <- c(min(maxit_1), max(maxit_1))
		}
		if( maxit_2_lim[1] != FALSE ){
			use_maxit2 <- maxit_2_lim 
		}else{
			use_maxit2 <- c(min(maxit_2), max(maxit_2))	
		}
		plot.new()
		plot.window(xlim = c(0, 1), ylim = c(0 ,1))
		use_types <- c("sample", "blank", "calibration", "spiked")
		use_types_col <- c("black", "green", "orange", "darkblue")
		use_types_col <- use_types_col[!is.na(match( use_types, measurements[those_files, "Type"]))]
		use_types <- use_types[!is.na(match( use_types, measurements[those_files, "Type"]))]			
		legend(x = -.02, y = 1,
			title = "File types:",
			legend = use_types,
			text.col = use_types_col,
			box.col = "white",
			cex = .7
		)
		plot.window(xlim = use_maxit1, ylim = use_maxit2)
		title(xlab = "Maximum intensity deviation", ylab = "Median intensity deviation")
		abline(h = 0, col = "lightgray", lty = 2)
		abline(v = 0, col = "lightgray", lty = 2)
		lines(maxit_1, maxit_2, col = "lightgray")
		text(maxit_1, maxit_2, labels = maxit_ID, col = maxit_col, cex = .8)
		box(); axis(1); axis(2);
		########################################################################
		plot.window(xlim = c(0,10), ylim = c(0,10))
	    if( (maxit_1_lim[1] != FALSE) | (maxit_2_lim[1] != FALSE) ){
			mtext("Zoomed in - click to zoom out partly or double-click to zoom out fully.", side = 3, line = 0.1, cex = .8, col = "darkgrey", at = 0, adj = 0)
	    }else{
			mtext("Brush and doubleclick to zoom.", side = 3, line = 0.1, cex = .8, col = "darkgrey", at = 0, adj = 0)
	    }
		plot.window(xlim = use_maxit1, ylim = use_maxit2)	# reset -> for proper zooming
		########################################################################

	}
    ############################################################################	
	
    ############################################################################
	return("done")
	
}




















