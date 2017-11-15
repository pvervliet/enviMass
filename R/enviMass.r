#' @title Trend detection in LC-MS data
#'
#' @description enviMass version 3 provides basic workflow functions and a webbrowser user interface for the proccessing of high resolution
#' LC-MS data, with emphasis on trend detection in time series of LC-MS files. 
#'
#' @details 
#' Incorporates functions from several other packages for peak picking, adduct grouping, isotope pattern calculation and grouping, etc. 
#' Has only been tested for Thermo Orbitrap data, centroided and converted to .mzXML via Proteo Wizard`s msconvert tool.
#' The workflow user interface (UI) has only been tested with google chrome and firefox, which must be set as default browsers. 
#'
#'
#' Exported functions:
#' \describe{
#' \item{\code{\link{webMass}}}{Start enviMass Workflow GUI.}
#' \item{\code{\link{recalib}}}{Function for m/z and RT recalibration.}
#' \item{\code{\link{PWfile}}}{Thermo .raw file conversion and centroidization with ProteoWizard.}
#' \item{\code{\link{profiletopeak}}}{Summarise profile information.}
#' \item{\code{\link{plotprofiles}}}{Interactive profile viewer.}
#' \item{\code{\link{plotglobal}}}{Plot all time series from a list of profiles.}
#' \item{\code{\link{plotaprofile}}}{Plot profiles.}
#' \item{\code{\link{plotagroup}}}{Plot time series of component profiles.}
#' \item{\code{\link{intensup}}}{Trend detection and blind subtraction for time profiles.}
#' \item{\code{\link{getID}}}{Find smallest available ID.}
#' \item{\code{\link{search_peak}}}{Find matches with peaklist.}
#' \item{\code{\link{screening}}}{Screens a sample and blank/blind peaklist for compound matches.}
#' \item{\code{\link{startprofiles}}}{Initializes a profile list object.}
#' \item{\code{\link{agglomer}}}{Partition the pooled peak lists generated by \code{\link{startprofiles}}.}
#' \item{\code{\link{partcluster}}}{Extract profiles from peak set partitions generated by \code{\link{agglomer}}.}
#' \item{\code{\link{mass_dens}}}{Kernel-density based mass estimation with confidence intervals.}
#' \item{\code{\link{mass_int}}}{Plot mass against intensity.}
#' }
#'
#'
#'
#' @docType package
#' @name enviMass
NULL

#' Exemplary data set of internal standard (IS) compounds.
#'
#' \itemize{
#'   \item ID. Unique compound IDs - character. 
#'   \item Name. Commpound name - character.
#'   \item Formula. Compound molecular formula - character.
#'   \item RT. Retention time - numeric or character.
#'   \item RT_tolerance. Character of RT Tolerance or FALSE if default used.
#'   \item main_adduct. Name of restrictive main adduct or FALSE if default used.
#'   \item ion_mode. "positive" or "negative".
#'   \item Use_for_recalibration. Logical, character.
#'   \item Use_for_screening. Logical, character.
#'   \item restrict_adduct. Logical, character.
#'   \item Remark. Character.
#'   \item tag1. Character.
#'   \item tag2. Character.
#'   \item tag3. Character.
#'   \item from. POSIXct date, time - character.
#'   \item to. POSIXct date, time - character.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name IS
#' @usage data(IS)
#' @format A data frame with 19 rows and 16 variables
NULL

#' Exemplary data set of target compounds.
#'
#' \itemize{
#'   \item ID. Unique compound IDs - character. 
#'   \item Name. Commpound name - character.
#'   \item Formula. Compound molecular formula - character.
#'   \item RT. Retention time - numeric or character.
#'   \item RT_tolerance. Character of RT Tolerance or FALSE if default used.
#'   \item ID_internal_standard. ID of internal standard - first column of IS data.frame.
#'   \item main_adduct. Name of restrictive main adduct or FALSE if default used.
#'   \item ion_mode. "positive" or "negative".
#'   \item Use_for_recalibration. Logical, character.
#'   \item Use_for_screening. Logical, character.
#'   \item restrict_adduct. Logical, character.
#'   \item Remark. Character.
#'   \item tag1. Character.
#'   \item tag2. Character.
#'   \item tag3. Character.
#'   \item from. POSIXct date, time - character.
#'   \item to. POSIXct date, time - character.
#'   \item intercept. Quantization parameter - character.
#'   \item slope. Quantization parameter - character.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name targets
#' @usage data(targets)
#' @format A data frame with 68 rows and 19 variables
NULL

#' @import enviPick
#' @import enviPat
#' @import nontarget
#' @import mgcv
#' @import shiny
#' @import shinyFiles
#' @import shinyBS
#' @import tcltk
#' @useDynLib enviMass extractProfiles neighbour agglom indexed fill_timeset meandel intdiff plot_prof binRT_prof binmz_prof
#'
#'


.onAttach <- function(lib, pkg)
{

	packageStartupMessage("\n \n Welcome to enviMass version 3.413 \n Run webMass() to start the enviMass browser UI \n\n");
	
	# add menus for enviMass & enviPick 
	if(interactive() && .Platform$OS.type == "windows" && .Platform$GUI == "Rgui" ){
		winMenuAdd("enviMass");
		winMenuAddItem("enviMass", "enviMass", "webMass()");
		#winMenuAddItem("enviMass", "enviPick", "webpick()");		
	}
	# copy the enviMass favicon to shiny`s localhost folder
	file.copy(	
			from=paste(path.package("enviMass", quiet = FALSE),"/favicon.ico",sep=""), 
			to=paste(path.package("shiny", quiet = FALSE),"/www/favicon.ico",sep=""), 
			overwrite = TRUE, 
			recursive = FALSE,
			copy.mode = TRUE
	)	
	
	
}

