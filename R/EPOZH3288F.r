EPOZH3288F <- function(
homol,
masslim = FALSE,
RTlim = FALSE,
massDlim = FALSE,
dmasslim = FALSE,
dRTlim = FALSE
){
if(!any(names(homol) == "homol_peaks_relat")){stop("Wrong homol argument in filter_segments!")}
if(any(dim(homol[["homol_peaks_relat"]]) == 0)){stop("Incomplete homol argument in filter_segments!")}
plot_those <- rep(TRUE, dim(homol[["homol_peaks_relat"]])[1])
if(masslim[1] != FALSE){
plot_those[
((homol[["homol_peaks_relat"]][,6] < masslim[1]) & (homol[["homol_peaks_relat"]][,8] < masslim[1]))
] <- FALSE
plot_those[
((homol[["homol_peaks_relat"]][,6] > masslim[2]) & (homol[["homol_peaks_relat"]][,8] > masslim[2]))
] <- FALSE
}
if(RTlim[1] != FALSE){
plot_those[
((homol[["homol_peaks_relat"]][,7] < RTlim[1]) & (homol[["homol_peaks_relat"]][,9] < RTlim[1]))
] <- FALSE
plot_those[
((homol[["homol_peaks_relat"]][,7] > RTlim[2]) & (homol[["homol_peaks_relat"]][,9] > RTlim[2]))
] <- FALSE
}
if(dmasslim[1] != FALSE){
plot_those[
(homol[["homol_peaks_relat"]][,3] < dmasslim[1]) | (homol[["homol_peaks_relat"]][,3] > dmasslim[2])
] <- FALSE
}
if(dRTlim[1] != FALSE){
plot_those[
(homol[["homol_peaks_relat"]][,4] < dRTlim[1]) | (homol[["homol_peaks_relat"]][,4] > dRTlim[2])
] <- FALSE
}
return(plot_those)
}
