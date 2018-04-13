YEJGA7087N <- function(profileList,progbar){
if(!profileList[[1]][[3]]){stop("profileList not profiled; aborted.")}
peaklist<-matrix(ncol=17,nrow=length(profileList[["index_prof"]][,1]),0)
colnames(peaklist)<-c(
"var_mz",
"profileID",
"number_peaks_total",
"past_incident",
"current_incident",
"current_intensity",
"Component",
"Homologue",
"Mass defect"
)
if(progbar==TRUE){  prog<-winProgressBar("Convert profiles to peaklist",min=1,max=length(profileList[["index_prof"]][,8]));
setWinProgressBar(prog, 0, title = "Convert profiles to peaklist...", label = NULL);}
for(k in 1:length(profileList[["index_prof"]][,1])){
if(progbar==TRUE){setWinProgressBar(prog, k, title = "Convert profiles to peaklist...", label = NULL)}
peaklist[k,"mean_m/z"]<-profileList[["index_prof"]][k,"mean_mz"]
peaklist[k,"mean_intensity"]<-mean(profileList[["peaks"]][(profileList[["index_prof"]][k,"start_ID"]:profileList[["index_prof"]][k,"end_ID"]),"intensity"])
peaklist[k,"mean_RT"]<-profileList[["index_prof"]][k,"mean_RT"]
peaklist[k,"max_intensity"]<-max(profileList[["peaks"]][(profileList[["index_prof"]][k,"start_ID"]:profileList[["index_prof"]][k,"end_ID"]),"intensity"])
peaklist[k,"in_blind?"]<-profileList[["index_prof"]][k,"in_blind?"]
peaklist[k,"above_blind?"]<-profileList[["index_prof"]][k,"above_blind?"]
peaklist[k,"var_mz"]<-var(profileList[["peaks"]][(profileList[["index_prof"]][k,"start_ID"]:profileList[["index_prof"]][k,"end_ID"]),"m/z"])
peaklist[k,"min_RT"]<-min(profileList[["peaks"]][(profileList[["index_prof"]][k,"start_ID"]:profileList[["index_prof"]][k,"end_ID"]),"RT"])
peaklist[k,"max_RT"]<-max(profileList[["peaks"]][(profileList[["index_prof"]][k,"start_ID"]:profileList[["index_prof"]][k,"end_ID"]),"RT"])
peaklist[k,"profileID"]<-profileList[["index_prof"]][k,"profile_ID"]
peaklist[k,"number_peaks_total"]<-profileList[["index_prof"]][k,"number_peaks_total"]
peaklist[k,"past_incident"]<-profileList[["index_prof"]][k,"deltaint_global"]
peaklist[k,"current_incident"]<-profileList[["index_prof"]][k,"deltaint_newest"]
peaklist[k,"current_intensity"]<-profileList[["index_prof"]][k,"newest_intensity"]
}
if(progbar==TRUE){ close(prog); }
return(peaklist);
}
