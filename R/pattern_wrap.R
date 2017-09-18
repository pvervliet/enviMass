
pattern_wrap <-
function(
	x,
	compound_table,
	logfile,
	...
){

	i <- x
	############################################################################	
	if( (compound_table[i,"restrict_adduct"] == "TRUE") & (compound_table[i,"main_adduct"] == "FALSE") ){
		return("nothing")
	}
	############################################################################
	if(compound_table[i,"ion_mode"] == "positive"){

		########################################################################
		takeall<-TRUE # use all adducts for recal?
		if(compound_table[i,"restrict_adduct"]=="FALSE"){ # use main adduct & all adducts - main adduct MUST come first -> for mz_ and RT_ below
			if(compound_table[i,"main_adduct"]=="FALSE"){ # no main adduct chosen?
				with_add<-c(logfile[["adducts_pos"]])
				takeall<-TRUE
			}else{ # with main adduct?
				with_add<-c(compound_table[i,"main_adduct"],logfile[["adducts_pos"]])
				with_add<-unique(with_add)
				takeall<-FALSE
			}
		}else{ # use main adduct only
			with_add<-c(compound_table[i,"main_adduct"])
			takeall<-TRUE # there is just one entry anyway
		}
		########################################################################

	}
	############################################################################
	if(compound_table[i,"ion_mode"] == "negative"){

		########################################################################
		takeall<-TRUE # use all adducts for recal?
		if(compound_table[i,"restrict_adduct"]=="FALSE"){ # use main adduct & all adducts - main adduct MUST come first -> for mz_ and RT_ below
			if(compound_table[i,"main_adduct"]=="FALSE"){ # no main adduct chosen?
				with_add<-c(logfile[["adducts_neg"]])
				takeall<-TRUE
			}else{ # with main adduct?
				with_add<-c(compound_table[i,"main_adduct"],logfile[["adducts_neg"]])
				with_add<-unique(with_add)
				takeall<-FALSE
			}
		}else{ # use main adduct only
			with_add<-c(compound_table[i,"main_adduct"])
			takeall<-TRUE # there is just one entry anyway
		}
		########################################################################

	}
	############################################################################
	if(!length(with_add)){return("nothing")}
	############################################################################		
	pattern_all <- list()	
	counter <- 1
	for(j in 1:length(with_add)){

		formelone<-compound_table[i,"Formula"];
		formelone<-enviPat::multiform(formelone,as.numeric(adducts[with_add[j] == adducts[,1],9]))
		if(as.character(adducts[with_add[j] == adducts[,1],7]) != "FALSE"){
			formelone<-enviPat::mergeform(formelone,as.character(adducts[with_add[j] == adducts[,1],7]))
		}
		if(as.character(adducts[with_add[j] == adducts[,1],8]) != "FALSE"){
			if(
				(enviPat::check_ded(formelone, as.character(adducts[with_add[j] == adducts[,1],8]))) == FALSE
			){
				formelone <- enviPat::subform(formelone,as.character(adducts[with_add[j] == adducts[,1],8]))
			}else{
				next; # just skip that deduct!
			}
		}
		pattern<-enviPat::isopattern(
			isotopes,
			formelone,
			threshold = 0.5,
			plotit = FALSE,
			charge = as.numeric(strsplit(as.character(adducts[with_add[j] == adducts[,1],3]),"+",fixed=TRUE)[[1]][1]),
			emass = 0.00054858,
			algo = 1
		)
		checked<-enviPat::check_chemform(isotopes, formelone)
		res<-try(enviPat::getR(checked, resolution_list[names(resolution_list) == logfile$parameters$resolution][[1]], nknots = 7, spar = 0.1, plotit = FALSE), silent=TRUE)
		if(grepl("Error",res)){ # no notification!
			if(checked[[3]] <= min(resolution_list[names(resolution_list) == logfile$parameters$resolution][[1]][,1])){
				use_this <- resolution_list[names(resolution_list) == logfile$parameters$resolution][[1]][,1] == min(resolution_list[names(resolution_list) == logfile$parameters$resolution][[1]][,1])
				res <- resolution_list[names(resolution_list) == logfile$parameters$resolution][[1]][use_this,2][1]
			}else{
				use_this <- resolution_list[names(resolution_list) == logfile$parameters$resolution][[1]][,1] == max(resolution_list[names(resolution_list) == logfile$parameters$resolution][[1]][,1])
				res <- resolution_list[names(resolution_list) == logfile$parameters$resolution][[1]][use_this,2][1]									
			}
		}
		pattern <- enviPat::envelope(
			pattern,
			ppm = FALSE,
			dmz = 0.0001,
			frac = 1/4,
			env = "Gaussian",
			resolution = res,
			plotit = FALSE
		)
		pattern <- enviPat::vdetect(
			pattern,
			detect = "centroid",
			plotit = FALSE
		)
		pattern_all[[counter]] <- pattern[[1]][order(pattern[[1]][, 2], decreasing = TRUE), ,drop=FALSE]
		names(pattern_all)[counter] <- paste0(
			compound_table$ID[i], "_", 
			with_add[j], "_", 
			compound_table$tag1[i], "_", 
			compound_table$tag2[i], "_", 
			compound_table$tag3[i]
		);
		counter <- (counter +1)
	}
	############################################################################
	pattern_all[[length(pattern_all) + 1]] <- takeall 
	return(pattern_all)
	
}



