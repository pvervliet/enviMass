#' @title Reset selection in the enviMass UI
#'
#' @export
#'
#' @description Reset selection in the enviMass UI
#'
#' @param session shiny session ID
#' 
#' @details enviMass workflow UI function. 
#' 

reset_selections<-function(session){

	################################################################
	updateSelectInput(session,"Ion_mode_Cal",selected = "none")
	updateSelectInput(session,inputId="Cal_file_set",selected = "none")
	updateNumericInput(session,"sel_meas_ID", value = 0)
	updateNumericInput(session,"sel_meas", value = 0)
	updateSelectInput(session,"Pos_compound_select",selected = "Choose")	
	updateSelectInput(session,"Neg_compound_select",selected = "Choose")	
	updateNumericInput(session,"sel_meas_comp", value = 0)	
	updateNumericInput(session,"sel_meas_comp_peak", value = 0)
	updateNumericInput(session,"sel_meas_comp_comp", value = 0)
	updateNumericInput(session,"atom_bound_peak", value = 0)
	updateNumericInput(session,"sel_scans_ID", value = 0)		
	################################################################
	#output$recal_pic<-renderImage({
	#	outfile <- tempfile(fileext='.png')
	#	png(outfile, width=400, height=400)
	#		plot.new();
	#		plot.window(xlim=c(1,1),ylim=c(1,1));
	#		box();
	#	text(1,1,label="not available",cex=1.5,col="darkred")
	#	dev.off()
	#	list(src = outfile,
	#		 alt = "")
	#}, deleteFile = TRUE)


	################################################################
	
}
