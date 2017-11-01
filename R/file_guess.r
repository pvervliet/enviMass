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
	file_name
){

	
    ############################################################################
  	# defaults
	guesses <- list()
  	guesses[[1]] <- "positive" 
  	guesses[[2]] <- "sample" 	
	names(guesses) <- c("Mode", "Type")
	############################################################################	
	if( grepl("pos", file_name, fixed = TRUE) ) guesses["Mode"] <- "positive" 
	if( grepl("Pos", file_name, fixed = TRUE) ) guesses["Mode"] <- "positive" 
	if( grepl("POS", file_name, fixed = TRUE) ) guesses["Mode"] <- "positive" 
	if( grepl("+", file_name, fixed = TRUE) ) guesses["Mode"] <- "positive" 
	if( grepl("neg", file_name, fixed = TRUE) ) guesses["Mode"] <- "negative" 
	if( grepl("Neg", file_name, fixed = TRUE) ) guesses["Mode"] <- "negative" 	
	if( grepl("NEG", file_name, fixed = TRUE) ) guesses["Mode"] <- "negative" 		
	############################################################################
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
	return(guesses)
	
}




















