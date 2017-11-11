#' @title 
#'
#' @description 
#'
#' @param masses Masses
#' @param intensities 
#' @param elements 
#' @param dmz 
#' @param ppm 
#' @param charges
#' @param isotopes
#' @param int_cut
#' @param inttol
#' @param use_C
#' @param must_peak Must a non-monoisotopic peak be present to have an atom count >0 for an element? Non applicable to single-isotopic elements.
#' 
#' @details enviMass workflow function
#' 

file_guess <-
function(
	file_name,
	propose = TRUE
){

	
    ############################################################################
  	# defaults
	guesses <- list()
	if(propose){
		guesses[[1]] <- "positive" 
		guesses[[2]] <- "sample" 	
		guesses[[3]] <- "2018-01-01" 
	}else{
		guesses[[1]] <- "FALSE" 
		guesses[[2]] <- "FALSE" 	
		guesses[[3]] <- "FALSE" 
	}
	names(guesses) <- c("Mode", "Type", "Date")
	############################################################################	
	# guess ionization mode #################################################### 
	if( grepl("pos", file_name, fixed = TRUE) ) guesses["Mode"] <- "positive" 
	if( grepl("Pos", file_name, fixed = TRUE) ) guesses["Mode"] <- "positive" 
	if( grepl("POS", file_name, fixed = TRUE) ) guesses["Mode"] <- "positive" 
	if( grepl("+", file_name, fixed = TRUE) ) guesses["Mode"] <- "positive" 
	if( grepl("neg", file_name, fixed = TRUE) ) guesses["Mode"] <- "negative" 
	if( grepl("Neg", file_name, fixed = TRUE) ) guesses["Mode"] <- "negative" 	
	if( grepl("NEG", file_name, fixed = TRUE) ) guesses["Mode"] <- "negative" 		
	############################################################################
	# guess file type ##########################################################
	if( grepl("sample", file_name, fixed = TRUE) ) guesses["Type"] <- "sample" 		
	if( grepl("Sample", file_name, fixed = TRUE) ) guesses["Type"] <- "sample" 
	if( grepl("SAMPLE", file_name, fixed = TRUE) ) guesses["Type"] <- "sample" 
	if( grepl("SAM", file_name, fixed = TRUE) ) guesses["Type"] <- "sample" 
	if( grepl("bl", file_name, fixed = TRUE) ) guesses["Type"] <- "blank" 
	if( grepl("Bl", file_name, fixed = TRUE) ) guesses["Type"] <- "blank" 
	if( grepl("BL", file_name, fixed = TRUE) ) guesses["Type"] <- "blank" 
	if( grepl("blank", file_name, fixed = TRUE) ) guesses["Type"] <- "blank" 	
	if( grepl("Blank", file_name, fixed = TRUE) ) guesses["Type"] <- "blank" 		
	if( grepl("BLANK", file_name, fixed = TRUE) ) guesses["Type"] <- "blank" 	
	if( grepl("blind", file_name, fixed = TRUE) ) guesses["Type"] <- "blank" 	
	if( grepl("Blind", file_name, fixed = TRUE) ) guesses["Type"] <- "blank" 		
	if( grepl("BLIND", file_name, fixed = TRUE) ) guesses["Type"] <- "blank" 	
	if( grepl("dot", file_name, fixed = TRUE) ) guesses["Type"] <- "spiked" 	
	if( grepl("Dot", file_name, fixed = TRUE) ) guesses["Type"] <- "spiked" 		
	if( grepl("DOT", file_name, fixed = TRUE) ) guesses["Type"] <- "spiked"  		
	if( grepl("cal", file_name, fixed = TRUE) ) guesses["Type"] <- "calibration" 		
	if( grepl("Cal", file_name, fixed = TRUE) ) guesses["Type"] <- "calibration" 
	if( grepl("CAL", file_name, fixed = TRUE) ) guesses["Type"] <- "calibration" 	
	############################################################################
	# guess date ###############################################################
	strsplits <- function(
		x, 
		splits,
		...
	){
		for (split in splits) x <- unlist(strsplit(x, split, ...))
		return(x[!x == ""]) # Remove empty values
	}
	that <- strsplits( file_name,  splits = c("_", "-", " ") )
	len <- length(that)
	found <- FALSE
	if(len >= 3 & !found){
		for(i in 1:(length(that)-2)){
			if(nchar(that[i]) != 4) next
			if(is.na(as.numeric(that[i]))) next
			if(nchar(that[i+1]) != 2) next
			if(is.na(as.numeric(that[i+1]))) next
			if(nchar(that[i+2]) != 2) next
			if(is.na(as.numeric(that[i+2]))) next			
			use_date <- paste(that[i], that[i+1], that[i+2], sep="-")
			use_date <- try(as.Date(use_date), silent = TRUE)
			if(class(use_date) != "try-error"){
				guesses["Date"] <- as.character(use_date)
				found <- TRUE
				break
			}
		}
	}
	if(len >= 3 & !found){
		for(i in 1:(length(that)-2)){
			if(nchar(that[i]) != 2) next
			if(is.na(as.numeric(that[i]))) next
			if(nchar(that[i+1]) != 2) next
			if(is.na(as.numeric(that[i+1]))) next
			if(nchar(that[i+2]) != 2) next
			if(is.na(as.numeric(that[i+2]))) next
			use_date <- paste(paste0("20",that[i]), that[i+1], that[i+2], sep="-")
			use_date <- try(as.Date(use_date), silent = TRUE)
			if(class(use_date) != "try-error"){
				guesses["Date"] <- as.character(use_date)
				found <- TRUE
				break
			}
		}
	}
	if(len >= 3 & !found){
		for(i in 1:length(that)){
			if(nchar(that[i]) != 6) next
			if(is.na(as.numeric(that[i]))) next
			use_date <- paste(paste0("20", substr(that[i], 1, 2)), substr(that[i], 3, 4), substr(that[i], 5, 6), sep="-")
			use_date <- try(as.Date(use_date), silent = TRUE)
			if(class(use_date) != "try-error"){
				guesses["Date"] <- as.character(use_date)
				found <- TRUE
				break
			}	
		}
	}
	if(len >= 3 & !found){
		for(i in 1:length(that)){
			if(nchar(that[i]) != 8) next
			if(is.na(as.numeric(that[i]))) next
			use_date <- paste(substr(that[i], 1, 4), substr(that[i], 5, 6), substr(that[i], 7, 8), sep="-")
			use_date <- try(as.Date(use_date), silent = TRUE)
			if(class(use_date) != "try-error"){
				guesses["Date"] <- as.character(use_date)
				found <- TRUE
				break
			}	
		}
	}
	############################################################################
	return(guesses)
	
}

