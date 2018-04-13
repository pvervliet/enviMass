BTYMI6379F <- function(
x,
logfile,
measurements,
ppm,
dmz,
dRT,
int_ratio
){
for_file <- x
for_mode <- measurements[measurements$ID == for_file, "Mode"]
if(any(objects()=="peaklist")){rm(peaklist)}
if(any(objects(envir=as.environment(".GlobalEnv"))=="peaklist")){rm(peaklist,envir=as.environment(".GlobalEnv"))}
load(file = file.path(logfile[[1]], "peaklist", for_file), envir = environment(), verbose=FALSE);
peaklist[,"keep_2"] <- Inf
sample_peaklist <- peaklist;
rm(peaklist)
report <- "nothing"
if(
(logfile$parameters$subtract_pos_bydate == "TRUE" & for_mode == "positive") ||
(logfile$parameters$subtract_neg_bydate == "TRUE" & for_mode == "negative")
){
atdate <- (measurements[,"Date"])
atdate <- as.Date(atdate, tz="GMT");
attime <- (measurements[,"Time"])
attime2 <- as.difftime(attime);
filetypus <- (measurements[,"Type"])
sampleID <- (measurements[,"ID"])
ord <- order(as.numeric(atdate), as.numeric(attime2), filetypus, sampleID);
i <- which(measurements$ID[ord] == for_file)
ionmode<-(measurements[,"Mode"])
found_blank <- FALSE
if(i > 1){
for(j in (i-1):1){
if(
(filetypus[ord[j]] == "blank") &
(ionmode[ord[j]] == for_mode)
){
blank_ID <- sampleID[ord[j]]
found_blank <- TRUE
break;
}
}
}
if(found_blank){
load(file = file.path(logfile[[1]],"peaklist",as.character(blank_ID)), verbose=FALSE, envir = environment());
peaks_blank <- peaklist[,c("m/z_corr","int_corr","RT_corr")]
rm(peaklist)
getit <- BYRIB9287R(
peaklist = peaks_blank,
mz = sample_peaklist[,"m/z_corr"],
dmz = (dmz * 2),
ppm = ppm,
RT = sample_peaklist[,"RT_corr"],
dRT = dRT,
onlymax = TRUE,
int_ratio = int_ratio,
int = sample_peaklist[,"int_corr"],
get_matches = FALSE,
get_ratio = TRUE
)
report <- "done"
which_affected <- (sample_peaklist[, "keep_2"] > getit)
sample_peaklist[which_affected, "keep_2"] <- getit[which_affected]
}
}
if(
(for_mode == "positive") &
(logfile$parameters$subtract_pos_byfile == "TRUE") &
any(logfile$Positive_subtraction_files != "FALSE")
){
selec_pos <- logfile$Positive_subtraction_files
selec_pos <- selec_pos[selec_pos != "FALSE"]
for(j in 1:length(selec_pos)){
subID <- strsplit(selec_pos[j]," - ")[[1]][1]
load(file = file.path(logfile[[1]], "peaklist", as.character(subID)), verbose = FALSE)
peaks_blank <- peaklist[,c("m/z_corr", "int_corr", "RT_corr")]; rm(peaklist);
getit <- BYRIB9287R(
peaklist = peaks_blank,
mz = sample_peaklist[,"m/z_corr"],
dmz = (dmz * 2),
ppm = ppm,
RT = sample_peaklist[,"RT_corr"],
dRT = dRT,
onlymax = TRUE,
int_ratio = int_ratio,
int = sample_peaklist[,"int_corr"],
get_matches = FALSE,
get_ratio = TRUE
)
which_affected <- (sample_peaklist[,"keep_2"] > getit)
sample_peaklist[which_affected,"keep_2"] <- getit[which_affected]
report <- "done"
rm(peaks_blank)
}
}
if(
(for_mode == "negative") &
(logfile$parameters$subtract_neg_byfile == "TRUE") &
any(logfile$Negative_subtraction_files != "FALSE")
){
selec_neg <- logfile$Negative_subtraction_files
selec_neg <- selec_neg[selec_neg != "FALSE"]
for(j in 1:length(selec_neg)){
subID <- strsplit(selec_neg[j]," - ")[[1]][1]
load(file = file.path(logfile[[1]],"peaklist",as.character(subID)),verbose=FALSE)
peaks_blank <- peaklist[,c("m/z_corr","int_corr","RT_corr")];
rm(peaklist);
getit <- BYRIB9287R(
peaklist = peaks_blank,
mz = sample_peaklist[,"m/z_corr"],
dmz = (dmz*2),
ppm = ppm,
RT = sample_peaklist[,"RT_corr"],
dRT = dRT,
onlymax = TRUE,
int_ratio = int_ratio,
int = sample_peaklist[,"int_corr"],
get_matches = FALSE,
get_ratio = TRUE
)
which_affected <- (sample_peaklist[,"keep_2"] > getit)
sample_peaklist[which_affected,"keep_2"] <- getit[which_affected]
report <- "done"
rm(peaks_blank)
}
}
if(report == "done"){
peaklist <- sample_peaklist;
save(peaklist, file = file.path(logfile[[1]], "peaklist", for_file))
rm(peaklist)
}
rm(sample_peaklist)
return(report)
}
