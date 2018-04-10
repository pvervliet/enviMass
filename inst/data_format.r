
	#######################################################################	
	
	system.time({
		load("D:/Users/uchemadmin/Desktop/MS/Projects/Aurea/Greifensee/Greifensee/MSlist/10")
	})
	path <-"D:/Users/uchemadmin/Desktop/MS/R_packages/data_files/data_format/output.msraw"
	
	system.time({
		load("D:/PART_4/MS/R_packages/data_files/data_format/2000")
	})
	path <-"D:/PART_4/MS/R_packages/data_files/data_format/output.bin"

	
	verbose <- TRUE
	#######################################################################

	#######################################################################	
	init_MSraw( # Rcpp
		file_path = path,
		sep_dim = 1,
		type = 1
	)	
	#######################################################################
	
	#######################################################################		
	read_header_MSraw_num(
		file_path = path
	)
	#######################################################################		
	MSraw_size()
	#######################################################################		

	#######################################################################		
	check_sum <- tools:::md5sum(path)[[1]]; print(check_sum);	
	write_checksum_MSraw(
		file_path = path,
		checksum = check_sum
	)
	#
	read_checksum_MSraw(
		file_path = path
	)
	#
	erase_checksum_MSraw(
		file_path = path
	)
	#######################################################################	
	
	#######################################################################	
	centroid_scan_numbers <- match(MSlist[["Scans"]][[2]][,"RT"], MSlist[["Scans"]][[1]]) # use package fastmatch?
	write_centroids_MSraw( # inits MSraw file, too
		file_path = path, 						# file path -> check for compatibility!!
		MSlist[["Scans"]][[1]], 				# scan_RTs
		seq(1:length(MSlist[["Scans"]][[1]])),	# scan_IDs
		MSlist[["Scans"]][[2]], 				# centroids matrix
		centroid_scan_numbers, 					# centroid -> scan_RTs index
		skip_prof_index = TRUE
	)
	#######################################################################	

	#######################################################################		
	read_scanRTs_MSraw(
		file_path = path
	)
	#######################################################################		
	
	#######################################################################		
	read_scanIDs_MSraw(
		file_path = path
	)
	#######################################################################		
	
	#######################################################################		
	system.time({
		centroids <- read_centroids_MSraw(
			file_path = path,
			insert_RT = TRUE,
			use_ID = FALSE,
			read_all = FALSE,
			get_index = TRUE
		)
	})
	head(centroids)
	head(MSlist[["Scans"]][[2]])
	tail(centroids)
	tail(MSlist[["Scans"]][[2]])
	identical(centroids[,"partID"], MSlist[["Scans"]][[2]][,"partID"])	
	identical(centroids[,"clustID"], MSlist[["Scans"]][[2]][,"clustID"])	
	identical(centroids[,"peakID"], MSlist[["Scans"]][[2]][,"peakID"])
	dim(centroids)	
	dim(MSlist[["Scans"]][[2]])
	#######################################################################		
	
	#######################################################################
	write_partition_MSraw(
		file_path = path,
		partition_index = MSlist[[5]]
	)
	#######################################################################	
	
	#######################################################################	
	partition_index <- read_partition_MSraw(
		file_path = path,
		insert_RT = TRUE,
		index_convert = TRUE
	)
	head(partition_index)
	head(MSlist[[5]])
	tail(partition_index)
	tail(MSlist[[5]])
	dim(partition_index)
	dim(MSlist[[5]])
	#######################################################################		
	
	#######################################################################
	write_EIC_MSraw(
		file_path = path,
		EIC_index = MSlist[[6]]
	)
	#######################################################################	
	
	#######################################################################	
	EIC_index <- read_EIC_MSraw(
		file_path = path,
		index_convert = TRUE
	)
	head(EIC_index)
	head(MSlist[[6]])
	tail(EIC_index)
	tail(MSlist[[6]])
	dim(EIC_index)
	dim(MSlist[[6]])
	#######################################################################		

	#######################################################################
	write_peaks_MSraw(
		file_path = path,
		peak_index = MSlist[[7]],
		peak_info = MSlist[[8]]
	)
	#######################################################################	
	
	#######################################################################	
	peak_index <- read_peaks_MSraw(
		file_path = path,
		index_convert = TRUE
	)
	
	head(peak_index[,1:3])
	head(MSlist[[7]])
	tail(peak_index[,1:3])
	tail(MSlist[[7]])
	dim(peak_index[,1:3])
	dim(MSlist[[7]])
	
	head(peak_index[,4:14])
	head(MSlist[[8]])
	#######################################################################		
		
	#######################################################################	
	EIC_indices <- c(1,1000,2000,3000)
	read_single_EIC_MSraw(
		file_path = path,
		EIC_indices,
		insert_RT = TRUE,
		use_ID = FALSE
	)
	#######################################################################		
	


	
	
	#######################################################################		
	system.time({
		for(i in 1:10){
			centroids <- filter_centroids_MSraw(
				file_path = path,
				min_RT = 300,
				max_RT = 500,
				min_mass = 400,
				max_mass = 552,
				insert_RT = TRUE,
				read_all = FALSE,
				use_partitions = FALSE
			)	
		}
	})			
	dim(centroids)
	#######################################################################		
	
	
	
	
	
	
	#######################################################################	
	m <- 1000
	run_check <- FALSE
	EIC_numbers <- sample(1:dim(MSlist[[6]])[1], m)
	system.time({
		for(i in 1:m){
			EIC_number <- EIC_numbers[i]
			EIC_list <- read_single_EIC_MSraw(
				file_path = path,
				EIC_indices = EIC_number,
				insert_RT = TRUE,
				use_ID = FALSE
			)
			if(run_check){
				EIC_list_orig <- 
					MSlist[["Scans"]][[2]][
						MSlist[[6]][EIC_number,1]:MSlist[[6]][EIC_number,2]
					,]
				if(!identical(EIC_list_orig, EIC_list)) stop("FUCKED")	
			}
		}
	})
	#######################################################################	
	
		
	
	
	
	
	
	


	