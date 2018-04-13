KXUVE9576N <- function(
pro_name,
pro_dir,
IS,
targets,
workflow_depend = FALSE,
workflow_must = FALSE,
skript_check = TRUE
){
if(any(ls() == "logfile")){stop("\n illegal logfile detected _1 in newproject.r!")}
if(grepl("\\", pro_dir, fixed = TRUE)){
pro_dir <- gsub("\\",.Platform$file.sep, pro_dir, fixed = TRUE)
}
if(!file.exists(file.path(pro_dir,pro_name))){
dir.create(file.path(pro_dir,pro_name), recursive = TRUE)
}
dir.create(file.path(pro_dir, pro_name, "files"), recursive = TRUE)
dir.create(file.path(pro_dir, pro_name, "MSlist"), recursive = TRUE)
dir.create(file.path(pro_dir, pro_name, "MSraw"), recursive = TRUE)
dir.create(file.path(pro_dir, pro_name, "peaklist"), recursive = TRUE)
dir.create(file.path(pro_dir, pro_name, "features"), recursive = TRUE)
dir.create(file.path(pro_dir, pro_name, "results"), recursive = TRUE)
dir.create(file.path(pro_dir, pro_name, "results", "screening"),recursive = TRUE)
dir.create(file.path(pro_dir, pro_name, "results", "LOD"), recursive = TRUE)
dir.create(file.path(pro_dir, pro_name, "results", "recalibration"), recursive = TRUE)
dir.create(file.path(pro_dir, pro_name, "results", "componentization"), recursive = TRUE)
dir.create(file.path(pro_dir, pro_name, "results", "componentization", "adducts"), recursive = TRUE)
dir.create(file.path(pro_dir, pro_name, "results", "componentization", "isotopologues"), recursive = TRUE)
dir.create(file.path(pro_dir, pro_name, "results", "componentization", "EIC_corr"), recursive = TRUE)
dir.create(file.path(pro_dir, pro_name, "results", "componentization", "homologues"), recursive = TRUE)
dir.create(file.path(pro_dir, pro_name, "results", "componentization", "components"), recursive = TRUE)
dir.create(file.path(pro_dir, pro_name, "dataframes"), recursive = TRUE)
dir.create(file.path(pro_dir, pro_name, "pics"), recursive = TRUE)
dir.create(file.path(pro_dir, pro_name, "exports"), recursive = TRUE)
dir.create(file.path(pro_dir, pro_name, "quantification"), recursive = TRUE)
write.table(IS, file = file.path(pro_dir, pro_name, "dataframes", "IS.txt"), row.names = FALSE, sep = "\t", quote = FALSE)
write.table(targets, file = file.path(pro_dir,pro_name, "dataframes", "targets.txt"), row.names = FALSE, sep = "\t", quote = FALSE)
measurements <- data.frame(c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),
c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),
c("-"),c("-"),c("-"),c("-"),c("-"),c("-"),c("-"));
names(measurements) <- c("ID", "Name", "Type", "Mode", "Place", "Date", "Time", "include", "copied", "peakpicking",
"qc", "recal", "align", "norm", "profiled", "LOD", "tag1", "tag2", "tag3", "date_end", "time_end",
"isotopologues", "adducts", "homologues", "EIC_correlation", "blind", "ID_2", "components_files")
write.csv(measurements, file = file.path(pro_dir, pro_name, "dataframes", "measurements"), row.names = FALSE)
if(!is.matrix(workflow_depend)){
workflow_depend <- read.table(file = "workflow_depend")
workflow_depend <- as.matrix(workflow_depend)
}
if(!is.matrix(workflow_must)){
workflow_must <- read.table(file = "workflow_must")
workflow_must <- as.matrix(workflow_must)
}
if(length(colnames(workflow_depend)) != length(rownames(workflow_depend))){
stop("\n Invalid matrix workflow_depend: matrix not quadratic; please revise!")
}
if(length(colnames(workflow_must)) != length(rownames(workflow_must))){
stop("\n Invalid matrix workflow_must: matrix not quadratic; please revise!")
}
if(any(is.na(match(colnames(workflow_depend),rownames(workflow_depend))))){
stop("\n Invalid matrix workflow_depend: unequal row vs. column names!")
}
if(any(is.na(match(colnames(workflow_must),rownames(workflow_must))))){
stop("\n Invalid matrix workflow_must: unequal row vs. column names!")
}
if(any(duplicated(colnames(workflow_must)))){
stop("\n Invalid matrix workflow_must: duplicated node names!")
}
if(any(duplicated(colnames(workflow_depend)))){
stop("\n Invalid matrix workflow_depend: duplicated node names!")
}
if(any(is.na(match(colnames(workflow_depend),colnames(workflow_must))))){
stop("\n Different nodes in workflow_depend vs. workflow_must: revise!")
}
if(any(is.na(match(colnames(workflow_must),colnames(workflow_depend))))){
stop("\n Different nodes in workflow_depend vs. workflow_must: revise!")
}
if(skript_check){
files <- list.files()
for(i in 1:length(colnames(workflow_depend))){
if(!any(files == paste("do_", colnames(workflow_depend)[i],".r", sep = ""))){
stop(paste("Missing do_ scripts.r for node ", colnames(workflow_depend)[i],sep = ""))
}
if(!any(files == paste("dont_", colnames(workflow_depend)[i],".r", sep = ""))){
stop(paste("Missing dont_ scripts.r for node ", colnames(workflow_depend)[i], sep = ""))
}
}
}
logfile <- list(0);
logfile[[1]] <- file.path(pro_dir,pro_name);
names(logfile)[1] <- c("project_folder")
logfile[[2]] <- rep(FALSE,length(colnames(workflow_must)));
names(logfile[[2]]) <- colnames(workflow_must)
names(logfile)[2] <- c("Tasks_to_redo");
tasks <- names(logfile[[2]])
doneit <- rep(FALSE,length(tasks))
summar <- data.frame(tasks,doneit,stringsAsFactors = FALSE)
names(summar) <- c("Tasks","Done?")
logfile[[3]] <- summar
names(logfile)[3] <- c("summary")
logfile[[4]] <- "C:/Program Files/ProteoWizard/ProteoWizard 3.0.11028/msconvert.exe"
names(logfile)[4] <- c("PW MSconvert path")
logfile[[5]] <- list(0)
names(logfile)[5] <- c("parameters")
logfile$parameters$verbose  <-  "TRUE"
logfile$parameters$test  <-  "FALSE"
logfile$parameters$method_use <-  "FALSE"
logfile$parameters$files_SIM <-  "FALSE"
logfile$parameters$peak_MSlevel  <-  "1";
logfile$parameters$peak_drtgap  <-  "300";
logfile$parameters$peak_dmzdens  <-  "3.5";
logfile$parameters$peak_minpeak  <-  "4";
logfile$parameters$peak_drtsmall2  <-  "8";
logfile$parameters$peak_drtfill  <-  "10";
logfile$parameters$peak_drtdens2  <-  "120";
logfile$parameters$peak_minint_log10  <-  "4";
logfile$parameters$peak_SN  <-  "5";
logfile$parameters$peak_SB <- "2";
logfile$parameters$peak_recurs <- "3";
logfile$parameters$peak_ended <- "1";
logfile$parameters$peak_weight <- "1";
logfile$parameters$peak_maxint_log10 <- "6.5";
logfile$parameters$peak_perc_cut <- "0";
logfile$parameters$peak_which_intensity <- "maximum"
logfile$parameters$cut_RT <- "FALSE"
logfile$parameters$cut_RT_min <- "0"
logfile$parameters$cut_RT_max <- "25"
logfile$parameters$cut_mass <- "FALSE"
logfile$parameters$cut_mass_min <- "0"
logfile$parameters$cut_mass_max <- "2000"
logfile$parameters$peak_estimate <- "TRUE"
logfile$parameters$peak_get_mass <- "mean"
logfile$parameters$progressBar <- "FALSE";
logfile$parameters$resolution <- "Elite_R240000@400";
logfile$parameters$recal_include_pos <- "TRUE"
logfile$parameters$recal_use_pos <- "Internal standards";
logfile$parameters$recal_dmz_pos <- "2";
logfile$parameters$recal_ppm_pos <- "FALSE";
logfile$parameters$recal_drt_pos <- "30";
logfile$parameters$recal_maxdmz_pos <- "2";
logfile$parameters$recal_include_neg <- "TRUE"
logfile$parameters$recal_use_neg <- "Internal standards";
logfile$parameters$recal_dmz_neg <- "2";
logfile$parameters$recal_ppm_neg <- "FALSE";
logfile$parameters$recal_drt_neg <- "30";
logfile$parameters$recal_maxdmz_neg <- "2";
logfile$parameters$replicate_dmz <- "12";
logfile$parameters$replicate_ppm <- "TRUE";
logfile$parameters$replicate_recalib <- "FALSE";
logfile$parameters$replicate_delRT <- "30";
logfile$parameters$replicate_IS_dInt <- "10";
logfile$parameters$replicate_mean_profiles <- "TRUE";
logfile$parameters$notrend <- "TRUE";
logfile$parameters$trend_lags <- "4,7,14";
logfile$parameters$trend_threshold <- "3";
logfile$parameters$trend_blind <- "yes";
logfile$parameters$blind_threshold <- "100";
logfile$parameters$blind_dmz <- "12";
logfile$parameters$blind_ppm <- "TRUE";
logfile$parameters$blind_drt <- "30";
logfile$parameters$subtract_pos_bydate <- "FALSE";
logfile$parameters$subtract_pos_byfile <- "FALSE";
logfile$parameters$subtract_neg_bydate <- "FALSE";
logfile$parameters$subtract_neg_byfile <- "FALSE";
logfile$parameters$blind_omit <- "no";
logfile$parameters$prof_maxfiles <- "100";
logfile$parameters$upto_file <- "FALSE";
logfile$parameters$prof_dmz <- "12";
logfile$parameters$prof_ppm <- "TRUE";
logfile$parameters$prof_drt <- "60";
logfile$parameters$prof_comp_maxfiles <- "15"
logfile$parameters$prof_select <- "TRUE";
logfile$parameters$replicates_prof <- "yes";
logfile$parameters$replicates_mean_prof <- "TRUE";
logfile$parameters$IS_drt1 <- "30";
logfile$parameters$IS_drt2 <- "2";
logfile$parameters$IS_dmz <- "8";
logfile$parameters$IS_ppm <- "TRUE";
logfile$parameters$IS_inttol <- "30";
logfile$parameters$IS_intcut <- "50000";
logfile$parameters$IS_w1 <- "0.8";
logfile$parameters$screen_IS_cutit <- "FALSE";
logfile$parameters$screen_IS_maxonly <- "FALSE";
logfile$parameters$screen_IS_restrict <- "FALSE";
logfile$parameters$screen_IS_restrict_many <- "10";
logfile$parameters$tar_drt1 <- "30";
logfile$parameters$tar_drt2 <- "2";
logfile$parameters$tar_dmz <- "8";
logfile$parameters$tar_ppm <- "TRUE";
logfile$parameters$tar_inttol <- "30";
logfile$parameters$tar_intcut <- "50000";
logfile$parameters$tar_w1 <- "0.8";
logfile$parameters$screen_target_cutit <- "FALSE";
logfile$parameters$screen_target_maxonly <- "FALSE";
logfile$parameters$screen_target_restrict <- "FALSE";
logfile$parameters$screen_target_restrict_many <- "10";
logfile$parameters$ISnorm_include_pos <-"TRUE"
logfile$parameters$ISnorm_percfiles_pos <- "90";
logfile$parameters$ISnorm_numbIS_pos <- "15";
logfile$parameters$ISnorm_medblank_pos <- "FALSE";
logfile$parameters$ISnorm_usesubblank_pos <- "TRUE";
logfile$parameters$ISnorm_numblank_pos <- "100";
logfile$parameters$ISnorm_medsam_pos <- "FALSE";
logfile$parameters$ISnorm_usesubsam_pos <- "TRUE";
logfile$parameters$ISnorm_numsam_pos <- "100";
logfile$parameters$ISnorm_score_pos <- "0.8";
logfile$parameters$ISnorm_include_neg <-"TRUE"
logfile$parameters$ISnorm_percfiles_neg <- "90";
logfile$parameters$ISnorm_numbIS_neg <- "15";
logfile$parameters$ISnorm_medblank_neg <- "FALSE";
logfile$parameters$ISnorm_usesubblank_neg <- "TRUE";
logfile$parameters$ISnorm_numblank_neg <- "100";
logfile$parameters$ISnorm_medsam_neg <- "FALSE";
logfile$parameters$ISnorm_usesubsam_neg <- "TRUE";
logfile$parameters$ISnorm_numsam_neg <- "100";
logfile$parameters$ISnorm_score_neg <- "0.8";
logfile$parameters$subtr_IS <- "yes";
logfile$parameters$subtr_target <- "yes";
logfile$parameters$subtr_blind <- "yes";
logfile$parameters$subtr_spiked <- "yes";
logfile$parameters$quant_files_included <- "100"
logfile$parameters$recov_files_included <- "20"
logfile$parameters$quant_digits <- "2"
logfile$parameters$isotop_mztol <- "8"
logfile$parameters$isotop_ppm <- "TRUE"
logfile$parameters$isotop_inttol <- "50"
logfile$parameters$isotop_rttol <- "2"
logfile$parameters$isotop_use_charges <- "FALSE"
logfile$parameters$adducts_rttol <- "2"
logfile$parameters$adducts_mztol <- "8"
logfile$parameters$adducts_ppm <- "TRUE"
logfile$parameters$do_atom_bounds_components <- "FALSE"
logfile$parameters$atom_bounds_components <- c("Cl","Br")
logfile$parameters$homol_units <- c("CH2,C2H4O")
logfile$parameters$homol_charges <- c("1,2")
logfile$parameters$homol_minmz <- "10"
logfile$parameters$homol_maxmz <- "120"
logfile$parameters$homol_minrt <- "10"
logfile$parameters$homol_maxrt <- "60"
logfile$parameters$homol_ppm <- "TRUE"
logfile$parameters$homol_mztol <- "8"
logfile$parameters$homol_rttol <- "20"
logfile$parameters$homol_minlength <- "6"
logfile$parameters$homol_vec_size <- "1E8"
logfile$parameters$homol_blind <- "FALSE"
logfile$parameters$homol_blind_value <- "10"
logfile$parameters$EICor_delRT <- "2"
logfile$parameters$EICor_minpeaks <- "15"
logfile$parameters$EICor_mincor <- ".95"
logfile$parameters$is_example <- "FALSE"
logfile$parameters$dofile_latest_profcomp <- "FALSE"
logfile$parameters$numfile_latest_profcomp <- "100"
logfile$parameters$filter_profcomp_pos <- "TRUE"
logfile$parameters$filter_profcomp_neg <- "TRUE"
logfile$parameters$for_which_profcomp_pos <- "all"
logfile$parameters$for_which_profcomp_neg <- "all"
logfile$parameters$prof_comp_link_only <- "FALSE"
logfile$parameters$corr_min_peaks <- "5"
logfile$parameters$comp_corr <- "0.9"
logfile$parameters$corr_del_RT <- "5"
logfile$parameters$corr_skip_peaks <- "TRUE"
logfile$parameters$parallel <- "TRUE"
logfile$parameters$parallel_restrict <- "FALSE"
logfile$parameters$parallel_cores <- "8"
logfile$parameters$parallel_x1 <- "FALSE"
logfile$parameters$parallel_x2 <- "FALSE"
logfile$parameters$parallel_x3 <- "FALSE"
if(skript_check){
source(file = "workflow_parameters.r", local = TRUE)
if(any(duplicated(names(logfile$parameters)))){stop("Duplicated parameter names found - revise!")}
}
logfile[[6]] <- list()
names(logfile)[6] <- c("workflow")
logfile$workflow <- 0
for(i in 1:length(names(logfile[[2]]))){
if(any(names(logfile[[2]])[i] == c("peakpicking", "qc") )){
logfile$workflow[i] <- "yes";
}else{
logfile$workflow[i] <- "no";
}
names(logfile$workflow)[i] <- names(logfile[[2]])[i]
}
logfile[[11]] <- workflow_depend
names(logfile)[11] <- "workflow_depend"
logfile[[12]] <- workflow_must
names(logfile)[12] <- "workflow_must"
schedule <- enviMass::PFMSW1067M(logfile$workflow_depend, logfile$workflow_must)
if(!is.data.frame(schedule)){stop("\nschedule not a data frame")}
set_order <- match(schedule[,1], logfile$summary[,1])
logfile$summary <- logfile$summary[set_order,]
logfile[[7]] <- 0
names(logfile)[7] <- c("adducts_pos")
logfile[[7]] <- c("M+H","M+Na","M+K","M+NH4");
logfile[[8]] <- 0
names(logfile)[8] <- c("adducts_neg")
logfile[[8]] <- "M-H";
logfile[[9]] <- "";
names(logfile)[9] <- c("isotopes")
logfile[[10]] <- as.character(packageVersion("enviMass"))
names(logfile)[10] <- c("version")
logfile[[13]] <- "FALSE"
names(logfile)[13] <- "Positive_subtraction_files"
logfile[[14]] <- "FALSE"
names(logfile)[14] <- "Negative_subtraction_files"
logfile[[15]] <- 0
names(logfile)[15] <- c("adducts_pos_group")
logfile[[15]] <- c("M+H","M+Na","M+K","M+NH4")
logfile[[16]] <- 0
names(logfile)[16] <- c("adducts_neg_group")
logfile[[16]] <- c("M-H","M-","2M-H")
logfile[[17]] <- list(0)
names(logfile)[17] <- c("UI_options")
logfile$UI_options$filterProf_minmass <- "0"
logfile$UI_options$filterProf_maxmass <- "3000"
logfile$UI_options$filterProf_minrt <- "0"
logfile$UI_options$filterProf_maxrt <- "100000"
logfile$UI_options$filterProf_minMD <- "-0.5"
logfile$UI_options$filterProf_maxMD <- "0.5"
logfile$UI_options$filterProf_medianblind <- "yes"
logfile$UI_options$filterProf_medianblind_value <- "10"
logfile$UI_options$filterProf_notblind <- "no"
logfile$UI_options$filterProf_sort <- "current trend intensity (decreasing)"
logfile$UI_options$filterProf_components <- "TRUE"
logfile[[18]] <- "Not available"
names(logfile)[18] <- "method_setup"
logfile[[19]] <- list()
names(logfile)[19] <- "comparisons"
cal_models_pos <- list()
save(cal_models_pos, file = file.path(logfile$project_folder, "quantification", "cal_models_pos"));
cal_models_neg <- list()
save(cal_models_neg, file= file.path(logfile$project_folder, "quantification", "cal_models_neg"));
save(logfile,file = file.path(pro_dir, pro_name, "logfile.emp"));
rm(logfile)
return(file.path(pro_dir, pro_name, "logfile.emp"));
}
