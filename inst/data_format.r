
	#######################################################################
	system.time({
		load("D:/Projects/LUBW/LUBW_Rhine_enviMass_Januar18/MSlist/10")
	})
	path <-"D:/PART_4/MS/R_packages/data_files/data_format/output.bin"
	
	
	verbose <- TRUE
	#######################################################################

	#######################################################################	
	init_MSraw( # Rcpp
		file_path = path
	)	
	#######################################################################
	
	#######################################################################	
	centroid_scan_numbers <- match(MSlist[["Scans"]][[2]][,"RT"], MSlist[["Scans"]][[1]]) # use package fastmatch?
	write_centroids_MSraw( # inits MSraw file, too
		file_path = path, 		# file path -> check for compatibility!!
		MSlist[["Scans"]][[1]], # scan_RTs
		MSlist[["Scans"]][[2]], # matrix
		centroid_scan_numbers
	)
	#######################################################################	

	#######################################################################
	write_partition_MSraw(
		file_path = path,
		partition_index = MSlist[[5]],
		add_ranges = TRUE
	)
	#######################################################################	
	
	#######################################################################		
	read_header_MSraw(
		file_path = path
	)
	#######################################################################		
			
	#######################################################################		
	read_scans_MSraw(
		file_path = path
	)
	#######################################################################		
	
	#######################################################################	
	partition_index <- read_partition_MSraw(
		file_path = path,
		insert_RT = TRUE,
		index_convert = FALSE
	)
	partition_index[1:10,]
	MSlist[[5]][1:10,]
	dim(partition_index)
	dim(MSlist[[5]])
	#######################################################################		
	
	
	#######################################################################		
	system.time({
		centroids <- read_centroids_MSraw(
			file_path = path,
			insert_RT = TRUE
		)
	})
	centroids[1:10,]
	MSlist[["Scans"]][[2]][1:10,]
	dim(centroids)	
	dim(MSlist[["Scans"]][[2]])
	#######################################################################		

	#######################################################################		
	system.time({		
		centroids <- filter_centroids_MSraw(
			file_path = path,
			min_RT = 410,
			max_RT = 450,
			min_mass = 200,
			max_mass = 201,
			insert_RT = TRUE
		)	
	})			
	#######################################################################		
	
	




	at_entry = (header_entries_start + 0 * centroid_size);	
	infile.seekg (at_entry, ios::beg);
	skip_1 = (centroid_size - double_size); 
	skip_2 = (centroid_size - (2 * double_size + int_size));
	for(n = 0; n < MSlist_centroids_nrow; n++){				
		infile.read (reinterpret_cast<char *> (&store_data_double), double_size);
		if((store_data_double < min_mass) || (store_data_double > max_mass)){
			infile.seekg (skip_1, ios::cur);
			continue;
		}
		infile.read (reinterpret_cast<char *> (&store_data_double), double_size);
		infile.read (reinterpret_cast<char *> (&store_data_int), int_size);		
		if((store_data_int < min_RT_ind) || (store_data_int > max_RT_ind)){ 
			infile.seekg (skip_2, ios::cur);
			continue;
		}
		keep_centroids(n) = 1;
		number_keep_centroids ++;
		infile.seekg (skip_2, ios::cur);	
	}



	
	#######################################################################
	write_MSlist( # Rcpp
		file_path = path, 		# file path -> check for compatibility!!
		MSlist[["Scans"]][[1]], # scan_RTs
		MSlist[["Scans"]][[2]], # matrix
		verbose = verbose
	)
	#######################################################################


	#######################################################################	
	from <- 1
	to <- 100
	
	got <- read_MSlist( # Rcpp
		file_path = path, 		# file path -> check for compatibility!!
		MSlist_scans_size = length(MSlist[["Scans"]][[1]])+1,
		from = from, # R indices, starting at 1
		to = to,
		verbose = verbose
	)
	print(got)
	print(MSlist[["Scans"]][[1]][from:to])
	#######################################################################	
	






	large2 <- rnorm(n = 1E9)

	#######################################################################
	write_MSlist( # Rcpp
		file_path = path, 		# file path -> check for compatibility!!
		large2, # scan_RTs
		MSlist[["Scans"]][[2]], # matrix
		verbose = verbose
	)
	#######################################################################

	#######################################################################	
	from <- 60000000
	to <- from + 10000
	system.time({
		for(i in 1:1000){
		got <- read_MSlist( # Rcpp
			file_path = path, 		# file path -> check for compatibility!!
			MSlist_scans_size = length(large2)+1,
			from = from, # R indices, starting at 1
			to = to,
			verbose = verbose
		)
		}
	})
	#print(got)
	#print(large2[from:to])
	#######################################################################	


	
	
	names(MSlist)


	