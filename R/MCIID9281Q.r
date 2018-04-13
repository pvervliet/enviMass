MCIID9281Q <- function(
cent_peak_combi,
pattern_centro,
peaks,
RT_tol_inside,
int_tol
){
if(any(duplicated(cent_peak_combi[,1]))){return(FALSE)}
if(any(duplicated(cent_peak_combi[,2]))){return(FALSE)}
rangeRT<-range(peaks[cent_peak_combi[,2],3])
if((rangeRT[2]-rangeRT[1])>RT_tol_inside){return(FALSE)}
if(int_tol<100){
if(length(cent_peak_combi[,1])>1){
for(n in 1:(length(cent_peak_combi[,1])-1)){
for(m in (n+1):length(cent_peak_combi[,1])){
ratio_int<-(peaks[cent_peak_combi[n,2],2]/peaks[cent_peak_combi[m,2],2])
ratio_int_theo_high<-(
(pattern_centro[cent_peak_combi[n,1],2]+(pattern_centro[cent_peak_combi[n,1],2]*int_tol/100))/
(pattern_centro[cent_peak_combi[m,1],2]-(pattern_centro[cent_peak_combi[m,1],2]*int_tol/100))
)
if(ratio_int_theo_high<ratio_int){return(FALSE)}
ratio_int_theo_low<-(
(pattern_centro[cent_peak_combi[n,1],2]-(pattern_centro[cent_peak_combi[n,1],2]*int_tol/100))/
(pattern_centro[cent_peak_combi[m,1],2]+(pattern_centro[cent_peak_combi[m,1],2]*int_tol/100))
)
if(ratio_int_theo_low>ratio_int){return(FALSE)}
}
}
}
}
return(TRUE)
}
