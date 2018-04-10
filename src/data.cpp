# include <iostream>
# include <Rcpp.h>
# include <fstream>

using namespace std;
using namespace Rcpp;

// -> Initialize a MSraw file: to contain file checks?
// [[Rcpp::export]]
bool init_MSraw(
	StringVector file_path,
	int sep_dim,
	int type
){
	int f = 0; 
	size_t n;	
	unsigned long int int_size = sizeof (int); 
	unsigned long int double_size = sizeof (double);
	unsigned long int char_size = sizeof (char);
	unsigned long int unsignedlongint_size = sizeof (unsigned long int);	
	unsigned long int header_entry_1, sep_dim_long = (unsigned long int) sep_dim, type_long = (unsigned long int) type;
	
	ofstream outfile;
    outfile.open(file_path(0), ios::binary | ios::out | ios::in | ios::trunc); // truncate -> init new file
	if(!outfile) return false;
	for(n = 0; n < 100; n++){ // header consists of bytes for 100 unsigned long int entries + ...
		outfile.write (reinterpret_cast<char *> (&f), unsignedlongint_size);
	}
	for(n = 0; n < 2000; n++){ // .... bytes for 2000 char entries 
		outfile.write (reinterpret_cast<char *> (&f), char_size);
	}
	outfile.seekp (0 * unsignedlongint_size, ios::beg); 						// write type
	outfile.write (reinterpret_cast<char *> (&type_long), unsignedlongint_size);
	outfile.seekp (1 * unsignedlongint_size, ios::beg); 						// write index for start of data blocks
	header_entry_1 = (100 * unsignedlongint_size + 2000 * char_size); 			
	outfile.write (reinterpret_cast<char *> (&header_entry_1), unsignedlongint_size);	
	outfile.seekp (89 * unsignedlongint_size, ios::beg);	
	outfile.write (reinterpret_cast<char *> (&int_size), unsignedlongint_size);
	outfile.seekp (90 * unsignedlongint_size, ios::beg);	
	outfile.write (reinterpret_cast<char *> (&double_size), unsignedlongint_size);
	outfile.seekp (91 * unsignedlongint_size, ios::beg);	
	outfile.write (reinterpret_cast<char *> (&char_size), unsignedlongint_size);
	outfile.seekp (92 * unsignedlongint_size, ios::beg);	
	outfile.write (reinterpret_cast<char *> (&unsignedlongint_size), unsignedlongint_size);
	outfile.seekp (93 * unsignedlongint_size, ios::beg);	// got_profile_data?	
	outfile.seekp (94 * unsignedlongint_size, ios::beg);	// got_centroids?	
	outfile.seekp (95 * unsignedlongint_size, ios::beg);	// got_filtered_centroids?	
	outfile.seekp (96 * unsignedlongint_size, ios::beg);	// got_part?	
	outfile.seekp (97 * unsignedlongint_size, ios::beg);	// got_EICs?
	outfile.seekp (98 * unsignedlongint_size, ios::beg);	// got_picked?	
	outfile.seekp (99 * unsignedlongint_size, ios::beg);	// write sep_dim
	outfile.write (reinterpret_cast<char *> (&sep_dim_long), unsignedlongint_size);	
	outfile.close ();
	
	return true;
}

// -> Report type sizes
// [[Rcpp::export]]
IntegerVector MSraw_size(){
	
	IntegerVector reading(4);
	reading(0) = sizeof (int);
	reading(1) = sizeof (double);	
	reading(2) = sizeof (char);	
	reading(3) = sizeof (unsigned long int);
	
	return reading;
	
}

// -> Read header of MSraw file
// WARNING: conversion of (unsigned long int) to R IntegerVector (int)
// [[Rcpp::export]]
IntegerVector read_header_MSraw_num(
	StringVector file_path
){
	
	int unsignedlongint_size = sizeof (unsigned long int);
	size_t n;
	unsigned long int header_entries;
	IntegerVector reading(100);
	
	ifstream infile;
    infile.open(file_path(0), ios::binary | ios::in); 	
	if(!infile){
		reading(0) = 0;
		return reading;
	}
	infile.seekg (0, ios::beg);
	for(n = 0; n < 100; n++){ // header consists of 500 bytes / unsignedlongint_size=4 bytes
		infile.read (reinterpret_cast<char *> (&header_entries), unsignedlongint_size);	
		reading(n) = (int) header_entries;
	}
	infile.close();	
	
	return reading;
}

// -> Write a MD5 checksum to MSraw file
// [[Rcpp::export]]
bool write_checksum_MSraw(
	StringVector file_path,
	StringVector checksum
){

	unsigned int checksum_size = 32; 
	int char_size = sizeof (char), unsignedlongint_size = sizeof (unsigned long int);
	size_t n;
	vector<char> writesum(32);
	
	ofstream outfile;
    outfile.open(file_path(0), ios::binary | ios::out | ios::in);
	if(!outfile) return false;	
	outfile.seekp (100 * unsignedlongint_size, ios::beg);	// write sep_dim
	for(n = 0; n < checksum_size; n++){
		writesum[n] = checksum[0][n];
		outfile.write (reinterpret_cast<char *> (&writesum[n]), char_size);
	}
	outfile.close ();
	
	return true;
}

// -> Read the MD5 checksum from a MSraw file
// [[Rcpp::export]]
CharacterVector read_checksum_MSraw(
	StringVector file_path
){

	unsigned int checksum_size = 32; 
	int char_size = sizeof (char), unsignedlongint_size = sizeof (unsigned long int);
	size_t n;
	std::string writesum;
	writesum.resize(checksum_size);
	
	ifstream infile;
    infile.open(file_path(0), ios::binary | ios::out | ios::in);	
	if(!infile) return false;	
	infile.seekg (100 * unsignedlongint_size, ios::beg);
	for(n = 0; n < checksum_size; n++){
		infile.read (reinterpret_cast<char *> (&writesum[n]), char_size);
	}
	infile.close ();

	return wrap(writesum);
}

// -> Erase the MD5 checksum from a MSraw file (set to 0 for a valid recalculation of the checksum)
// [[Rcpp::export]]
bool erase_checksum_MSraw(
	StringVector file_path
){
	unsigned int checksum_size = 32;
	int unsignedlongint_size = sizeof (unsigned long int), f = 0;
	size_t n;
	ofstream outfile;
    outfile.open(file_path(0), ios::binary | ios::out | ios::in);
	if(!outfile) return false;	
	outfile.seekp (100 * unsignedlongint_size, ios::beg);
	for(n = 0; n < checksum_size; n++){
		outfile.write (reinterpret_cast<char *> (&f), unsignedlongint_size);
	}
	outfile.close ();
	
	return true;
}

// -> Write centroids into MSraw file
// [[Rcpp::export]]
bool write_centroids_MSraw(
	StringVector file_path,
	NumericVector scan_RTs,
	IntegerVector scan_IDs,
	NumericMatrix centroids,
	IntegerVector match_scans_RT,
	bool skip_prof_index
){
	
	int measureID;
	size_t n;
	int char_size = sizeof (char), int_size = sizeof (int), unsignedlongint_size = sizeof (unsigned long int), double_size = sizeof (double);
	unsigned int scans_size = scan_RTs.size();
	unsigned int centroids_nrow = centroids.nrow();
	unsigned long int centroid_size = (2 * double_size + 2 * int_size);
	unsigned long int header_entry_1, header_entry_2, header_entry_3, header_entry_4, header_entry_5, header_entry_6, header_entry_7;
	
	ofstream outfile;
    outfile.open(file_path(0), ios::binary | ios::out | ios::in);
	if(!outfile) return false;
	// first three entries to MSraw header:
	// index start of scan RT data block - written in init_MSraw, fixed size
	header_entry_1 = (100 * unsignedlongint_size + 2000 * char_size); 
	outfile.seekp (1 * unsignedlongint_size, ios::beg); 
	outfile.write (reinterpret_cast<char *> (&header_entry_1), unsignedlongint_size);
	// index start of scan ID data block, variable size
	header_entry_2 = (header_entry_1 + (scans_size * double_size));
	outfile.seekp (2 * unsignedlongint_size, ios::beg);
	outfile.write (reinterpret_cast<char *> (&header_entry_2), unsignedlongint_size);	
	// index start of profile data block, variable size - inserted, but empty for the time being
	header_entry_3 = (header_entry_2 + (scans_size * int_size));
	outfile.seekp (3 * unsignedlongint_size, ios::beg);
	outfile.write (reinterpret_cast<char *> (&header_entry_3), unsignedlongint_size);	
	// index start of filtered centroid data block, variable size - inserted, but empty for the time being
	header_entry_4 = (header_entry_3 + 0);	
	outfile.seekp (4 * unsignedlongint_size, ios::beg);
	outfile.write (reinterpret_cast<char *> (&header_entry_4), unsignedlongint_size);	
	// index start of non-filtered centroid data block, variable size
	header_entry_5 = (header_entry_4 + 0);	
	outfile.seekp (5 * unsignedlongint_size, ios::beg);
	outfile.write (reinterpret_cast<char *> (&header_entry_5), unsignedlongint_size);	
	// index start of profile data index block, variable size
	header_entry_6 = (header_entry_5 + (centroids_nrow * centroid_size));
	outfile.seekp (6 * unsignedlongint_size, ios::beg); 
	outfile.write (reinterpret_cast<char *> (&header_entry_6), unsignedlongint_size);	
	// if skip_prof_index: no profile data index for centroids - set start of centroid index, else omit
	if(skip_prof_index){
		header_entry_7 = header_entry_6;
		outfile.seekp (7 * unsignedlongint_size, ios::beg); 
		outfile.write (reinterpret_cast<char *> (&header_entry_7), unsignedlongint_size);		
	}
	// entries to scan RT data block:
	outfile.seekp (header_entry_1, ios::beg); 
	for(n = 0; n < scans_size; n++){
		outfile.write (reinterpret_cast<char *> (&scan_RTs(n)), double_size);	
	}
	// entries to scan RT ID block:
	outfile.seekp (header_entry_2, ios::beg); 
	for(n = 0; n < scans_size; n++){
		outfile.write (reinterpret_cast<char *> (&scan_IDs(n)), int_size);	
	}
	// entries to centroid data block:
	outfile.seekp (header_entry_5, ios::beg);
	for(n = 0; n < centroids_nrow; n++){	
		outfile.write (reinterpret_cast<char *> (&centroids(n, 0)), double_size);	// m/z
		outfile.write (reinterpret_cast<char *> (&centroids(n, 1)), double_size);	// intensity
		outfile.write (reinterpret_cast<char *> (&match_scans_RT(n)), int_size); 				// RT index
		measureID = (int) centroids(n, 3);
		outfile.write (reinterpret_cast<char *> (&measureID), int_size);	// measureID		
	}	
	outfile.seekp (95 * unsignedlongint_size, ios::beg);	// got_filtered_centroids?
	n = 1;
	outfile.write (reinterpret_cast<char *> (&n), unsignedlongint_size);
	
	outfile.close ();

	return true;
}

// -> Read scan RTs from MSraw file
// [[Rcpp::export]]
NumericVector read_scanRTs_MSraw(
	StringVector file_path
){

	int double_size = sizeof (double), unsignedlongint_size = sizeof (unsigned long int);
	unsigned long int header_entries_start, header_entries_end, scans_size;
	size_t n;
	
	ifstream infile;
    infile.open(file_path(0), ios::binary | ios::in); 	
	if(!infile) return 0;

	infile.seekg (1 * unsignedlongint_size, ios::beg);
	infile.read (reinterpret_cast<char *> (&header_entries_start), unsignedlongint_size);
	infile.read (reinterpret_cast<char *> (&header_entries_end), unsignedlongint_size);
	if(header_entries_end <= header_entries_start) return 0;
	scans_size = ((header_entries_end - header_entries_start) / double_size);
	NumericVector reading(scans_size);
	infile.seekg (header_entries_start, ios::beg);
	for(n = 0; n < scans_size; n++){		
		infile.read (reinterpret_cast<char *> (&reading[n]), double_size);
	}
	infile.close();	

	return reading;
}

// -> Read scan IDs from MSraw file
// [[Rcpp::export]]
IntegerVector read_scanIDs_MSraw(
	StringVector file_path
){

	int int_size = sizeof (int), unsignedlongint_size = sizeof (unsigned long int);
	unsigned long int header_entries_start, header_entries_end, scans_size;
	size_t n;
	
	ifstream infile;
    infile.open(file_path(0), ios::binary | ios::in); 	
	if(!infile) return 0;

	infile.seekg (2 * unsignedlongint_size, ios::beg);
	infile.read (reinterpret_cast<char *> (&header_entries_start), unsignedlongint_size);
	infile.read (reinterpret_cast<char *> (&header_entries_end), unsignedlongint_size);
	if(header_entries_end <= header_entries_start) return 0;
	scans_size = ((header_entries_end - header_entries_start) / int_size);
	IntegerVector reading(scans_size);
	infile.seekg (header_entries_start, ios::beg);
	for(n = 0; n < scans_size; n++){		
		infile.read (reinterpret_cast<char *> (&reading[n]), int_size);
	}
	infile.close();	

	return reading;
}

// -> Read all centroids from MSraw file
// [[Rcpp::export]]
NumericMatrix read_centroids_MSraw(
	StringVector file_path,
	bool insert_RT,
	bool use_ID,
	bool read_all,
	bool get_index
){

	size_t n, m;
	int int_size = sizeof (int), double_size = sizeof (double), unsignedlongint_size = sizeof (unsigned long int);
	int RT, measureID;
	unsigned long int centroid_size = (unsigned long int) (2 * double_size + 2 * int_size);
	unsigned long int header_entries_start, header_entries_end, centroids_nrow, scans_size;
	unsigned long int f1, f2, header_entry_5, header_entry_7, header_entry_8, header_entry_9, header_entry_10, partition_index_size;
	unsigned int partition_size = (3 * unsignedlongint_size + 2 * int_size + 2 * double_size);
	unsigned long int got_partition_index = 0, got_EIC_index = 0, got_peak_index = 0;
	int skip_1 = (unsignedlongint_size + 2 * double_size + 2 * int_size);
	int skip_2 = (unsignedlongint_size + 11 * double_size);
	unsigned long int EIC_index_size, peak_index_size;
	unsigned int EIC_size = (3 * unsignedlongint_size); 
	unsigned int peak_size = (3 * unsignedlongint_size + 10 * double_size);
	
	ifstream infile;
    infile.open(file_path(0), ios::binary | ios::in); 
	if(!infile) return 0;
	infile.seekg (95 * unsignedlongint_size, ios::beg);	// got_filtered_centroids?
	infile.read (reinterpret_cast<char *> (&n), unsignedlongint_size);	
	if(n == 0){ // no centroids added yet
		infile.close();	
		return 0;
	}
	infile.seekg (96 * unsignedlongint_size, ios::beg);	// got_partition?
	infile.read (reinterpret_cast<char *> (&got_partition_index), unsignedlongint_size);	
	infile.seekg (97 * unsignedlongint_size, ios::beg);	// got_partition?
	infile.read (reinterpret_cast<char *> (&got_EIC_index), unsignedlongint_size);		
	infile.seekg (98 * unsignedlongint_size, ios::beg);	// got_partition?
	infile.read (reinterpret_cast<char *> (&got_peak_index), unsignedlongint_size);		
	if(read_all){ 	// include filtered centroids
		infile.seekg (4 * unsignedlongint_size, ios::beg);
		infile.read (reinterpret_cast<char *> (&header_entries_start), unsignedlongint_size);
	}else{			// omit filtered centroids
		infile.seekg (5 * unsignedlongint_size, ios::beg);
		infile.read (reinterpret_cast<char *> (&header_entries_start), unsignedlongint_size);	
	}
	infile.seekg (6 * unsignedlongint_size, ios::beg);
	infile.read (reinterpret_cast<char *> (&header_entries_end), unsignedlongint_size);	
	centroids_nrow = ((header_entries_end - header_entries_start) / centroid_size);
	NumericMatrix reading(centroids_nrow , 7);
	infile.seekg (header_entries_start, ios::beg);
	for(n = 0; n < centroids_nrow; n++){		
		infile.read (reinterpret_cast<char *> (&reading(n, 0)), double_size);
		infile.read (reinterpret_cast<char *> (&reading(n, 1)), double_size);
		infile.read (reinterpret_cast<char *> (&RT), int_size);		
		infile.read (reinterpret_cast<char *> (&measureID), int_size);	
		reading(n, 2) = (double) RT;
		reading(n, 3) = (double) measureID;
		if(!get_index || (got_partition_index == 0)) reading(n, 4) = -1;
		if(!get_index || (got_EIC_index == 0)) reading(n, 5) = -1;
		if(!get_index || (got_peak_index == 0)) reading(n, 6) = -1;		
	}
	// insert scan RTs
	if(insert_RT & !use_ID){
		infile.seekg (1 * unsignedlongint_size, ios::beg);
		infile.read (reinterpret_cast<char *> (&header_entries_start), unsignedlongint_size);
		infile.read (reinterpret_cast<char *> (&header_entries_end), unsignedlongint_size);
		scans_size = ((header_entries_end - header_entries_start) / double_size);
		NumericVector readScans(scans_size);
		infile.seekg (header_entries_start, ios::beg);
		for(n = 0; n < scans_size; n++){		
			infile.read (reinterpret_cast<char *> (&readScans[n]), double_size);
		}
		for(n = 0; n < centroids_nrow; n++){		
			reading(n, 2) = readScans(((int) reading(n, 2)) -1);
		}
	}
	// insert scan IDs
	if(insert_RT & use_ID){
		infile.seekg (2 * unsignedlongint_size, ios::beg);
		infile.read (reinterpret_cast<char *> (&header_entries_start), unsignedlongint_size);
		infile.read (reinterpret_cast<char *> (&header_entries_end), unsignedlongint_size);
		scans_size = ((header_entries_end - header_entries_start) / int_size);
		IntegerVector readScans(scans_size);
		infile.seekg (header_entries_start, ios::beg);
		for(n = 0; n < scans_size; n++){		
			infile.read (reinterpret_cast<char *> (&readScans[n]), int_size);
		}
		for(n = 0; n < centroids_nrow; n++){		
			reading(n, 2) = (double) readScans(((int) reading(n, 2)) -1);
		}
	}
	// insert indices
	if(get_index){
		// insert partition IDs
		if(got_partition_index == 1){ 
			infile.seekg (5 * unsignedlongint_size, ios::beg);
			infile.read (reinterpret_cast<char *> (&header_entry_5), unsignedlongint_size);
			infile.seekg (7 * unsignedlongint_size, ios::beg);
			infile.read (reinterpret_cast<char *> (&header_entry_7), unsignedlongint_size);
			infile.read (reinterpret_cast<char *> (&header_entry_8), unsignedlongint_size);
			partition_index_size = ((header_entry_8 - header_entry_7) / (unsigned long int) partition_size);
			infile.seekg (header_entry_7, ios::beg);
			for(n = 0; n < partition_index_size; n++){
				infile.read (reinterpret_cast<char *> (&f1), unsignedlongint_size);	
				f1 = (((f1 - header_entry_5) / centroid_size));
				infile.read (reinterpret_cast<char *> (&f2), unsignedlongint_size);	
				f2 = (((f2 - header_entry_5) / centroid_size));	
				for(m = f1; m <= f2; m++) reading(m, 4) = (n + 1);
				infile.ignore(skip_1);  // skip remaining partition index infos
			}
		}
		// insert EIC IDs
		if(got_EIC_index == 1){ 
			infile.seekg (8 * unsignedlongint_size, ios::beg);
			infile.read (reinterpret_cast<char *> (&header_entry_8), unsignedlongint_size);
			infile.read (reinterpret_cast<char *> (&header_entry_9), unsignedlongint_size);
			EIC_index_size = ((header_entry_9 - header_entry_8) / EIC_size);
			infile.seekg (header_entry_8, ios::beg);		
			for(n = 0; n < EIC_index_size; n++){
				infile.read (reinterpret_cast<char *> (&f1), unsignedlongint_size);	
				f1 = (((f1 - header_entry_5) / centroid_size));
				infile.read (reinterpret_cast<char *> (&f2), unsignedlongint_size);	
				f2 = (((f2 - header_entry_5) / centroid_size));
				for(m = f1; m <= f2; m++) reading(m, 5) = (n + 1);
				infile.ignore(unsignedlongint_size);  // skip remaining partition index infos				
			}
		}
		// insert peak IDs
		if(got_peak_index == 1){ 
			infile.seekg (5 * unsignedlongint_size, ios::beg);
			infile.read (reinterpret_cast<char *> (&header_entry_5), unsignedlongint_size);
			infile.seekg (9 * unsignedlongint_size, ios::beg);
			infile.read (reinterpret_cast<char *> (&header_entry_9), unsignedlongint_size);
			infile.read (reinterpret_cast<char *> (&header_entry_10), unsignedlongint_size);
			peak_index_size = ((header_entry_10 - header_entry_9) / peak_size);
			infile.seekg (header_entry_9, ios::beg);		
			for(n = 0; n < peak_index_size; n++){	
				infile.read (reinterpret_cast<char *> (&f1), unsignedlongint_size);	
				f1 = (((f1 - header_entry_5) / centroid_size));
				infile.read (reinterpret_cast<char *> (&f2), unsignedlongint_size);	
				f2 = (((f2 - header_entry_5) / centroid_size));		
				for(m = f1; m <= f2; m++) reading(m, 6) = (n + 1);
				infile.ignore(skip_2);  // skip remaining peak index infos	
			}
		}		
	}
	
	colnames(reading) = CharacterVector::create("m/z", "intensity", "RT", "measureID", "partID", "clustID", "peakID");
	infile.close();	

	return reading;
}

// -> Write partition index to MSraw file
// [[Rcpp::export]]
bool write_partition_MSraw(
	StringVector file_path,
	NumericMatrix partition_index
){

	unsigned long int d, partition_index_size = (unsigned long int) partition_index.nrow(); 
	size_t n, m;
	unsigned long int header_entry_5, header_entry_7, header_entry_8;
	int int_size = sizeof (int), unsignedlongint_size = sizeof (unsigned long int), double_size = sizeof (double);
	unsigned int dif;
	int partition_size = (3 * unsignedlongint_size + 2 * int_size + 2 * double_size); // RT range set as int -> cp. centroid block -> cp. scans block
	unsigned long int start_centroid, end_centroid, centroid_size = (unsigned long int) (2 * double_size + 2 * int_size);
	double minmass, maxmass, get_mass;
	int minRT, maxRT, get_RT; 
	
	fstream inoutfile;
    inoutfile.open(file_path(0), ios::binary | ios::out | ios::in); 
	if(!inoutfile) return false;
	// index start of non-filtered centroid data block	
	inoutfile.seekg (5 * unsignedlongint_size, ios::beg);
	inoutfile.read (reinterpret_cast<char *> (&header_entry_5), unsignedlongint_size);
	// index start of partition index block - abort if not set yet (0) during eic or profile index writing
	inoutfile.seekg (7 * unsignedlongint_size, ios::beg);
	inoutfile.read (reinterpret_cast<char *> (&header_entry_7), unsignedlongint_size);
	if(header_entry_7 == 0){
		inoutfile.close ();
		return false;
	}
	inoutfile.seekg (0 * unsignedlongint_size, ios::beg);
	// index start of EIC index block
	header_entry_8 = (header_entry_7 + partition_index_size * partition_size);
	inoutfile.seekp (8 * unsignedlongint_size, ios::beg);
	inoutfile.write (reinterpret_cast<char *> (&header_entry_8), unsignedlongint_size);	
	inoutfile.seekp (0 * unsignedlongint_size, ios::beg);
	NumericMatrix mass_ranges(partition_index_size, 2);		
	IntegerMatrix RT_ranges(partition_index_size, 2);	
	inoutfile.seekg (header_entry_5, ios::beg);	
	for(n = 0; n < partition_index_size; n++){	
		start_centroid = (header_entry_5 + ((((unsigned long int) partition_index(n, 0)) - 1) * centroid_size));
		//end_centroid = (header_entry_5 + ((((unsigned long int) partition_index(n, 1)) - 1) * centroid_size));
		dif = (int) (partition_index(n, 1) - partition_index(n, 0));
		inoutfile.seekg (start_centroid, ios::beg);
		minRT = 1E6, maxRT = 0, minmass = 1E6, maxmass = 0;
		for(m = 0; m < (dif + 1); m++){ // centroids must all be in consecutive blocks
				inoutfile.read (reinterpret_cast<char *> (&get_mass), double_size);
				inoutfile.ignore(double_size); // skip intensity
				inoutfile.read (reinterpret_cast<char *> (&get_RT), int_size);	
				if(get_mass < minmass) minmass = get_mass;
				if(get_mass > maxmass) maxmass = get_mass;
				if(get_RT < minRT) minRT = get_RT;
				if(get_RT > maxRT) maxRT = get_RT;		
				inoutfile.ignore(int_size);  // skip centroid membership indices		
		}
		mass_ranges(n, 0) = minmass;
		mass_ranges(n, 1) = maxmass;		
		RT_ranges(n, 0) = minRT;
		RT_ranges(n, 1) = maxRT;		
	}
	inoutfile.seekg (0, ios::beg); // -> otherwise ends at header_entry3 = next seekp -> data race problem
	inoutfile.seekp (header_entry_7, ios::beg);	// write partition index including RT and mass ranges
	for(n = 0; n < partition_index_size; n++){
		start_centroid = (header_entry_5 + ((((unsigned long int) partition_index(n, 0)) - 1) * centroid_size));
		inoutfile.write (reinterpret_cast<char *> (&start_centroid), unsignedlongint_size);	// start_ID of first centroid in partition
		end_centroid = (header_entry_5 + ((((unsigned long int) partition_index(n, 1)) - 1) * centroid_size));
		inoutfile.write (reinterpret_cast<char *> (&end_centroid), unsignedlongint_size);	// start_ID of last centroid in partition
		d = (unsigned long int) partition_index(n, 2);
		inoutfile.write (reinterpret_cast<char *> (&d), unsignedlongint_size);	// number_peaks
		inoutfile.write (reinterpret_cast<char *> (&RT_ranges(n, 0)), int_size);	//		
		inoutfile.write (reinterpret_cast<char *> (&RT_ranges(n, 1)), int_size);	//	
		inoutfile.write (reinterpret_cast<char *> (&mass_ranges(n, 0)), double_size);	//	
		inoutfile.write (reinterpret_cast<char *> (&mass_ranges(n, 1)), double_size);	//				
	}
	inoutfile.seekp (96 * unsignedlongint_size, ios::beg);	// got_partition?
	n = 1;
	inoutfile.write (reinterpret_cast<char *> (&n), unsignedlongint_size);
	inoutfile.close ();

	return true;
}

// -> Read partition matrix from MSraw file
// [[Rcpp::export]]
NumericMatrix read_partition_MSraw(
	StringVector file_path,
	bool insert_RT,
	bool index_convert
){
	
	unsigned long int fu, partition_index_size, scans_size;
	size_t n;
	unsigned long int header_entry_1, header_entry_2, header_entry_5, header_entry_7, header_entry_8;
	int fi, int_size = sizeof (int), unsignedlongint_size = sizeof (unsigned long int), double_size = sizeof (double);	
	unsigned int partition_size = (3 * unsignedlongint_size + 2 * int_size + 2 * double_size); 
	double fd;
	unsigned long int centroid_size = (unsigned long int) (2 * double_size + 2 * int_size);
	
	ifstream infile;
    infile.open(file_path(0), ios::binary | ios::in); 	
	if(!infile) return 0;
	infile.seekg (96 * unsignedlongint_size, ios::beg);	// got_partitions?
	infile.read (reinterpret_cast<char *> (&n), unsignedlongint_size);	
	if(n == 0){ // no centroids added yet
		infile.close();	
		return 0;
	}
	infile.seekg (1 * unsignedlongint_size, ios::beg);
	infile.read (reinterpret_cast<char *> (&header_entry_1), unsignedlongint_size);
	infile.read (reinterpret_cast<char *> (&header_entry_2), unsignedlongint_size);
	infile.seekg (5 * unsignedlongint_size, ios::beg);
	infile.read (reinterpret_cast<char *> (&header_entry_5), unsignedlongint_size);
	infile.seekg (7 * unsignedlongint_size, ios::beg);
	infile.read (reinterpret_cast<char *> (&header_entry_7), unsignedlongint_size);
	infile.read (reinterpret_cast<char *> (&header_entry_8), unsignedlongint_size);
	partition_index_size = ((header_entry_8 - header_entry_7) / (unsigned long int) partition_size);
	NumericMatrix reading(partition_index_size, 7);	
	infile.seekg (header_entry_7, ios::beg);		
	for(n = 0; n < partition_index_size; n++){		
		infile.read (reinterpret_cast<char *> (&fu), unsignedlongint_size);	
		if(index_convert) fu = (((fu - header_entry_5) / centroid_size) + 1);
		reading(n, 0) = (double) fu;
		infile.read (reinterpret_cast<char *> (&fu), unsignedlongint_size);	
		if(index_convert) fu = (((fu - header_entry_5) / centroid_size) + 1);
		reading(n, 1) = (double) fu;		
		infile.read (reinterpret_cast<char *> (&fu), unsignedlongint_size);	
		reading(n, 2) = (double) fu;		
		infile.read (reinterpret_cast<char *> (&fi), int_size);		// RT_min
		reading(n, 3) = (double) fi;	
		infile.read (reinterpret_cast<char *> (&fi), int_size);		// RT_max
		reading(n, 4) = (double) fi;		
		infile.read (reinterpret_cast<char *> (&fd), double_size);	// mass_min
		reading(n, 5) = fd; // R numeric = double
		infile.read (reinterpret_cast<char *> (&fd), double_size);	// mass_max	
		reading(n, 6) = fd; // R numeric = double		
	}
	if(insert_RT){
		scans_size = ((header_entry_2 - header_entry_1) / double_size);
		NumericVector readScans(scans_size);
		infile.seekg (header_entry_1, ios::beg);
		for(n = 0; n < scans_size; n++){		
			infile.read (reinterpret_cast<char *> (&readScans[n]), double_size);
		}		
		for(n = 0; n < partition_index_size; n++){				
			reading(n, 3) = readScans(((int) reading(n, 3)) -1);
			reading(n, 4) = readScans(((int) reading(n, 4)) -1);			
		}
	}
	colnames(reading) = CharacterVector::create("start_ID", "end_ID", "number_peaks", "min_RT", "max_RT", "min_mass", "max_mass");
	infile.close();
	
	return reading;	

}

// -> Write EIC index to MSraw file
// [[Rcpp::export]]
bool write_EIC_MSraw(
	StringVector file_path,
	NumericMatrix EIC_index
){

	unsigned long int d, EIC_index_size = (unsigned long int) EIC_index.nrow(); 
	size_t n;
	unsigned long int header_entry_5, header_entry_8, header_entry_9;
	
	int int_size = sizeof (int), unsignedlongint_size = sizeof (unsigned long int), double_size = sizeof (double);
	int EIC_size = (3 * unsignedlongint_size);
	unsigned long int start_centroid, end_centroid, centroid_size = (unsigned long int) (2 * double_size + 2 * int_size);
	fstream inoutfile;
    inoutfile.open(file_path(0), ios::binary | ios::out | ios::in); 
	if(!inoutfile) return 0;
	// index start of centroid data block
	inoutfile.seekg (5 * unsignedlongint_size, ios::beg);
	inoutfile.read (reinterpret_cast<char *> (&header_entry_5), unsignedlongint_size);
	// index start of EIC index block 
	inoutfile.seekg (8 * unsignedlongint_size, ios::beg);
	inoutfile.read (reinterpret_cast<char *> (&header_entry_8), unsignedlongint_size);	
	if(header_entry_8 == 0){
		inoutfile.close ();
		return false;
	}	
	header_entry_9 = (header_entry_8 + EIC_index_size * EIC_size);
	inoutfile.seekp (9 * unsignedlongint_size, ios::beg);
	inoutfile.write (reinterpret_cast<char *> (&header_entry_9), unsignedlongint_size);	
	inoutfile.seekg (0, ios::beg); // -> otherwise ends at header_entry4 = next seekp -> data race problem
	inoutfile.seekp (header_entry_8, ios::beg);
	for(n = 0; n < EIC_index_size; n++){
		start_centroid = (header_entry_5 + ((((unsigned long int) EIC_index(n, 0)) - 1) * centroid_size));
		inoutfile.write (reinterpret_cast<char *> (&start_centroid), unsignedlongint_size);	// start_ID of first centroid in partition
		end_centroid = (header_entry_5 + ((((unsigned long int) EIC_index(n, 1)) - 1) * centroid_size));
		inoutfile.write (reinterpret_cast<char *> (&end_centroid), unsignedlongint_size);	// start_ID of last centroid in partition
		d = (unsigned long int) EIC_index(n, 2);
		inoutfile.write (reinterpret_cast<char *> (&d), unsignedlongint_size);	// number_peaks				
	}
	inoutfile.seekp (97 * unsignedlongint_size, ios::beg);	// got_partition?
	n = 1;
	inoutfile.write (reinterpret_cast<char *> (&n), unsignedlongint_size);
	inoutfile.close ();

	return true;
}

// -> Read EIC matrix from MSraw file
// [[Rcpp::export]]
NumericMatrix read_EIC_MSraw(
	StringVector file_path,
	bool index_convert
){
	
	unsigned long int fu, EIC_index_size;
	size_t n;
	unsigned long int header_entry_5, header_entry_8, header_entry_9;	
	int int_size = sizeof (int), unsignedlongint_size = sizeof (unsigned long int), double_size = sizeof (double);	
	int EIC_size = (3 * unsignedlongint_size); 
	unsigned long int centroid_size = (unsigned long int) (2 * double_size + 2 * int_size);
	
	ifstream infile;
    infile.open(file_path(0), ios::binary | ios::in); 	
	if(!infile) return 0;
	infile.seekg (97 * unsignedlongint_size, ios::beg);	// got_EICs?
	infile.read (reinterpret_cast<char *> (&n), unsignedlongint_size);	
	if(n == 0){ // no centroids added yet
		infile.close();	
		return 0;
	}
	infile.seekg (5 * unsignedlongint_size, ios::beg);
	infile.read (reinterpret_cast<char *> (&header_entry_5), unsignedlongint_size);
	infile.seekg (8 * unsignedlongint_size, ios::beg);
	infile.read (reinterpret_cast<char *> (&header_entry_8), unsignedlongint_size);
	infile.read (reinterpret_cast<char *> (&header_entry_9), unsignedlongint_size);
	EIC_index_size = ((header_entry_9 - header_entry_8) / EIC_size);
	NumericMatrix reading(EIC_index_size, 3);	
	infile.seekg (header_entry_8, ios::beg);		
	for(n = 0; n < EIC_index_size; n++){		
		infile.read (reinterpret_cast<char *> (&fu), unsignedlongint_size);	
		if(index_convert) fu = (((fu - header_entry_5) / centroid_size) + 1);
		reading(n, 0) = (double) fu;
		infile.read (reinterpret_cast<char *> (&fu), unsignedlongint_size);	
		if(index_convert) fu = (((fu - header_entry_5) / centroid_size) + 1);
		reading(n, 1) = (double) fu;		
		infile.read (reinterpret_cast<char *> (&fu), unsignedlongint_size);	
		reading(n, 2) = (double) fu;				
	}
	colnames(reading) = CharacterVector::create("start_ID", "end_ID", "number_peaks");
	infile.close();
	
	return reading;	

}

// -> Write peak index to MSraw file
// [[Rcpp::export]]
bool write_peaks_MSraw(
	StringVector file_path,
	NumericMatrix peak_index,
	NumericMatrix peak_info
){

	unsigned long int d, peak_index_size = (unsigned long int) peak_index.nrow(), peak_info_size = (unsigned long int) peak_info.nrow(); 
	size_t n;
	if(peak_index_size != peak_info_size) return 0; 
	unsigned long int header_entry_5, header_entry_9, header_entry_10;
	int int_size = sizeof (int), unsignedlongint_size = sizeof (unsigned long int), double_size = sizeof (double);
	int peak_size = (3 * unsignedlongint_size + 10 * double_size);
	
	unsigned long int start_centroid, end_centroid, centroid_size = (unsigned long int) (2 * double_size + 2 * int_size);
	fstream inoutfile;
    inoutfile.open(file_path(0), ios::binary | ios::out | ios::in); 
	if(!inoutfile) return 0;
	// index start of centroid data block
	inoutfile.seekg (5 * unsignedlongint_size, ios::beg);
	inoutfile.read (reinterpret_cast<char *> (&header_entry_5), unsignedlongint_size);
	// index start of peak index block 
	inoutfile.seekg (9 * unsignedlongint_size, ios::beg);
	inoutfile.read (reinterpret_cast<char *> (&header_entry_9), unsignedlongint_size);	
	if(header_entry_9 == 0){
		inoutfile.close ();
		return false;
	}	
	header_entry_10 = (header_entry_9 + peak_index_size * peak_size);
	inoutfile.seekp (10 * unsignedlongint_size, ios::beg);
	inoutfile.write (reinterpret_cast<char *> (&header_entry_10), unsignedlongint_size);	
	inoutfile.seekg (0, ios::beg); // -> otherwise data race problem
	inoutfile.seekp (header_entry_9, ios::beg);
	for(n = 0; n < peak_index_size; n++){
		start_centroid = (header_entry_5 + ((((unsigned long int) peak_index(n, 0)) - 1) * centroid_size));
		inoutfile.write (reinterpret_cast<char *> (&start_centroid), unsignedlongint_size);	// start_ID of first centroid in peak
		end_centroid = (header_entry_5 + ((((unsigned long int) peak_index(n, 1)) - 1) * centroid_size));
		inoutfile.write (reinterpret_cast<char *> (&end_centroid), unsignedlongint_size);	// start_ID of last centroid in peak
		d = (unsigned long int) peak_index(n, 2);
		inoutfile.write (reinterpret_cast<char *> (&d), unsignedlongint_size);	// number_centroids			
		inoutfile.write (reinterpret_cast<char *> (&peak_info(n, 0)), double_size); // m/z
		inoutfile.write (reinterpret_cast<char *> (&peak_info(n, 1)), double_size); // var m/z
		inoutfile.write (reinterpret_cast<char *> (&peak_info(n, 2)), double_size); // max_int
		inoutfile.write (reinterpret_cast<char *> (&peak_info(n, 3)), double_size); // sum_int
		inoutfile.write (reinterpret_cast<char *> (&peak_info(n, 4)), double_size); // RT
		inoutfile.write (reinterpret_cast<char *> (&peak_info(n, 5)), double_size); // minRT
		inoutfile.write (reinterpret_cast<char *> (&peak_info(n, 6)), double_size); // maxRT
		inoutfile.write (reinterpret_cast<char *> (&peak_info(n, 7)), double_size); // part_ID
		inoutfile.write (reinterpret_cast<char *> (&peak_info(n, 8)), double_size); // EIC_ID
		inoutfile.write (reinterpret_cast<char *> (&peak_info(n, 9)), double_size); // peak_ID
		inoutfile.write (reinterpret_cast<char *> (&peak_info(n, 10)), double_size); // score
	}
	inoutfile.seekp (98 * unsignedlongint_size, ios::beg);	// got_peaks?
	n = 1;
	inoutfile.write (reinterpret_cast<char *> (&n), unsignedlongint_size);
	inoutfile.close ();

	return true;
}

// -> Read peak matrix from MSraw file
// [[Rcpp::export]]
NumericMatrix read_peaks_MSraw(
	StringVector file_path,
	bool index_convert
){
	
	unsigned long int fu, peak_index_size;
	size_t n;
	unsigned long int header_entry_5, header_entry_9, header_entry_10;
	int int_size = sizeof (int), unsignedlongint_size = sizeof (unsigned long int), double_size = sizeof (double);	
	int peak_size = (3 * unsignedlongint_size + 10 * double_size);
	unsigned long int centroid_size = (unsigned long int) (2 * double_size + 2 * int_size);
	
	ifstream infile;
    infile.open(file_path(0), ios::binary | ios::in); 	
	if(!infile) return 0;
	infile.seekg (98 * unsignedlongint_size, ios::beg);	// got_peaks?
	infile.read (reinterpret_cast<char *> (&n), unsignedlongint_size);	
	if(n == 0){ // no centroids added yet
		infile.close();	
		return 0;
	}
	infile.seekg (5 * unsignedlongint_size, ios::beg);
	infile.read (reinterpret_cast<char *> (&header_entry_5), unsignedlongint_size);
	infile.seekg (9 * unsignedlongint_size, ios::beg);
	infile.read (reinterpret_cast<char *> (&header_entry_9), unsignedlongint_size);
	infile.read (reinterpret_cast<char *> (&header_entry_10), unsignedlongint_size);
	peak_index_size = ((header_entry_10 - header_entry_9) / peak_size);
	NumericMatrix reading(peak_index_size, 14);	
	infile.seekg (header_entry_9, ios::beg);		
	for(n = 0; n < peak_index_size; n++){		
		infile.read (reinterpret_cast<char *> (&fu), unsignedlongint_size);	
		if(index_convert) fu = (((fu - header_entry_5) / centroid_size) + 1);
		reading(n, 0) = (double) fu;
		infile.read (reinterpret_cast<char *> (&fu), unsignedlongint_size);	
		if(index_convert) fu = (((fu - header_entry_5) / centroid_size) + 1);
		reading(n, 1) = (double) fu;		
		infile.read (reinterpret_cast<char *> (&fu), unsignedlongint_size);	
		reading(n, 2) = (double) fu;		
		infile.read (reinterpret_cast<char *> (&reading(n, 3)), double_size);	// m/z
		infile.read (reinterpret_cast<char *> (&reading(n, 4)), double_size);	// var m/z
		infile.read (reinterpret_cast<char *> (&reading(n, 5)), double_size);	// max_int
		infile.read (reinterpret_cast<char *> (&reading(n, 6)), double_size);	// sum_int
		infile.read (reinterpret_cast<char *> (&reading(n, 7)), double_size);	// RT
		infile.read (reinterpret_cast<char *> (&reading(n, 8)), double_size);	// minRT
		infile.read (reinterpret_cast<char *> (&reading(n, 9)), double_size);	// maxRT
		infile.read (reinterpret_cast<char *> (&reading(n, 10)), double_size);	// part_ID
		infile.read (reinterpret_cast<char *> (&reading(n, 11)), double_size);	// EIC_ID		
		infile.read (reinterpret_cast<char *> (&reading(n, 12)), double_size);	// peak_ID
		infile.read (reinterpret_cast<char *> (&reading(n, 13)), double_size);	// score
	}
	colnames(reading) = CharacterVector::create("start_ID", "end_ID", "number_peaks", "m/z", "var_m/z", "max_int", 
		"sum_int", "RT", "minRT", "maxRT", "part_ID", "EIC_ID", "peak_ID", "Score");
	infile.close();
	
	return reading;	

}

// -> Read single EICs from MSraw file
// [[Rcpp::export]]
List read_single_EIC_MSraw(
	StringVector file_path,
	IntegerVector EIC_indices,
	bool insert_RT,
	bool use_ID
){
	
	size_t n, m;
	unsigned long int centroid_start, centroid_end, centroids_nrow, at_EIC_index;
	unsigned long int header_entry_5, header_entry_8, header_entry_9;	
	int max_EIC_index, int_size = sizeof (int), unsignedlongint_size = sizeof (unsigned long int), double_size = sizeof (double);	
	int k, EIC_size = (3 * unsignedlongint_size), EIC_indices_size = EIC_indices.size(); 
	int RT, measureID;
	unsigned long int header_entries_start, header_entries_end, scans_size;
	
	ifstream infile;
    infile.open(file_path(0), ios::binary | ios::in); 	
	if(!infile) return 0;
	infile.seekg (97 * unsignedlongint_size, ios::beg);	// got_EICs?
	infile.read (reinterpret_cast<char *> (&n), unsignedlongint_size);	
	if(n == 0){ // no centroids added yet
		infile.close();	
		return 0;
	}
	infile.seekg (5 * unsignedlongint_size, ios::beg);
	infile.read (reinterpret_cast<char *> (&header_entry_5), unsignedlongint_size);
	infile.seekg (8 * unsignedlongint_size, ios::beg);
	infile.read (reinterpret_cast<char *> (&header_entry_8), unsignedlongint_size);
	infile.read (reinterpret_cast<char *> (&header_entry_9), unsignedlongint_size);
	max_EIC_index = (int) ((header_entry_9 - header_entry_8) / EIC_size);
	List reading_list(EIC_indices_size);

	if(!use_ID){
		infile.seekg (1 * unsignedlongint_size, ios::beg);
		infile.read (reinterpret_cast<char *> (&header_entries_start), unsignedlongint_size);
		infile.read (reinterpret_cast<char *> (&header_entries_end), unsignedlongint_size);
		scans_size = ((header_entries_end - header_entries_start) / double_size);
	}else{
		infile.seekg (2 * unsignedlongint_size, ios::beg);
		infile.read (reinterpret_cast<char *> (&header_entries_start), unsignedlongint_size);
		infile.read (reinterpret_cast<char *> (&header_entries_end), unsignedlongint_size);	
		scans_size = ((header_entries_end - header_entries_start) / int_size);
	}
	NumericVector readScans(scans_size);
	infile.seekg (header_entries_start, ios::beg);
	if(!use_ID){	
		for(n = 0; n < scans_size; n++){		
			infile.read (reinterpret_cast<char *> (&readScans[n]), double_size);
		}
	}else{	
		for(n = 0; n < scans_size; n++){		
			infile.read (reinterpret_cast<char *> (&k), int_size);
			readScans[n] = (double) k;
		}
	}
	for(k = 0; k < EIC_indices_size; k++){
		if((EIC_indices[k] - 1) > max_EIC_index) continue;
		// find position of EIC index entry & find position of EIC
		at_EIC_index = header_entry_8 + ((EIC_indices[k] - 1) * EIC_size);
		infile.seekg (at_EIC_index, ios::beg);
		infile.read (reinterpret_cast<char *> (&centroid_start), unsignedlongint_size);
		infile.read (reinterpret_cast<char *> (&centroid_end), unsignedlongint_size);	
		infile.read (reinterpret_cast<char *> (&centroids_nrow), unsignedlongint_size);				
		// allocate EIC matrix
		NumericMatrix reading(centroids_nrow, 4);	
		// read centroids
		infile.seekg (centroid_start, ios::beg);
		for(m = 0; m < centroids_nrow; m++){
			infile.read (reinterpret_cast<char *> (&reading(m, 0)), double_size);
			infile.read (reinterpret_cast<char *> (&reading(m, 1)), double_size);
			infile.read (reinterpret_cast<char *> (&RT), int_size);		
			infile.read (reinterpret_cast<char *> (&measureID), int_size);	
			reading(m, 2) = (double) RT;
			if(insert_RT) reading(m, 2) = readScans(((int) reading(m, 2)) -1);
			reading(m, 3) = (double) measureID;
		}
		colnames(reading) = CharacterVector::create("m/z", "intensity", "RT", "measureID");
		reading_list[k] = reading;
	
	}
	infile.close();	
	
	return reading_list;	

}

// -> Filter centroids from MSraw file
// [[Rcpp::export]]
NumericMatrix filter_centroids_MSraw(
	StringVector file_path,
	int min_RT,
	int max_RT,
	int min_mass,
	int max_mass,
	bool insert_RT,
	bool read_all,
	bool use_partitions
){

	int int_size = sizeof (int), double_size = sizeof (double), unsignedlongint_size = sizeof (unsigned long int);
	int centroid_size = (2 * double_size + 2 * int_size), skip_1, skip_2, skip_3, skip_4, fi;
	int min_RT_ind, max_RT_ind, store_data_int;
	int RT, measureID;
	size_t m, n;
	unsigned long int f1, f2, header_entries_start, header_entries_end, header_entry_7, header_entry_8, at_entry, got_partition_index = 0;
	unsigned long int MSlist_centroids_nrow, MSlist_scans_size, number_keep_centroids = 0, partition_index_size;
	double store_data_double;
	unsigned int partition_size = (3 * unsignedlongint_size + 2 * int_size + 2 * double_size);
	ifstream infile;
    infile.open(file_path(0), ios::binary | ios::in); 
	if(!infile) return 0;
	// get RT index range
	infile.seekg (1 * unsignedlongint_size, ios::beg);
	infile.read (reinterpret_cast<char *> (&header_entries_start), unsignedlongint_size);
	infile.read (reinterpret_cast<char *> (&header_entries_end), unsignedlongint_size);
	MSlist_scans_size = ((header_entries_end - header_entries_start) / double_size);
	NumericVector readScans(MSlist_scans_size);
	infile.seekg (header_entries_start, ios::beg);
	for(n = 0; n < MSlist_scans_size; n++) infile.read (reinterpret_cast<char *> (&readScans[n]), double_size);
	min_RT_ind = MSlist_scans_size;
	for(n = 0; n < MSlist_scans_size; n++){			
		if(readScans[n] >= min_RT){ 
			min_RT_ind = (n + 1);
			break;
		}
	}
	max_RT_ind = 1;
	for(n = (MSlist_scans_size - 1); n >= 0; n--){	
		if(readScans[n] <= max_RT){ 	
			max_RT_ind = (n + 1);
			break;
		}
	}
	if(max_RT_ind < min_RT_ind) return 0;
	infile.seekg (96 * unsignedlongint_size, ios::beg);	// got_partition?
	infile.read (reinterpret_cast<char *> (&got_partition_index), unsignedlongint_size);	
	if(read_all){ 	// include filtered centroids
		infile.seekg (4 * unsignedlongint_size, ios::beg);
		infile.read (reinterpret_cast<char *> (&header_entries_start), unsignedlongint_size);
	}else{			// omit filtered centroids
		infile.seekg (5 * unsignedlongint_size, ios::beg);
		infile.read (reinterpret_cast<char *> (&header_entries_start), unsignedlongint_size);	
	}
	infile.seekg (6 * unsignedlongint_size, ios::beg);
	infile.read (reinterpret_cast<char *> (&header_entries_end), unsignedlongint_size);	
	MSlist_centroids_nrow = ((header_entries_end - header_entries_start) / centroid_size);
	IntegerVector keep_centroids(MSlist_centroids_nrow);	
	// -> prefilter via partitions ----------------------------------- //
	if(use_partitions && (got_partition_index == 1) && !read_all){
		for(n = 0; n < MSlist_centroids_nrow; n++) keep_centroids(n) = 0;
		skip_1 = unsignedlongint_size;
		skip_2 = (int_size + 2 * double_size);
		skip_3 = (2 * double_size);		
		skip_4 = double_size;			
		infile.seekg (7 * unsignedlongint_size, ios::beg);
		infile.read (reinterpret_cast<char *> (&header_entry_7), unsignedlongint_size);
		infile.read (reinterpret_cast<char *> (&header_entry_8), unsignedlongint_size);
		partition_index_size = ((header_entry_8 - header_entry_7) / (unsigned long int) partition_size);
		infile.seekg (header_entry_7, ios::beg);
		for(n = 0; n < partition_index_size; n++){
			infile.read (reinterpret_cast<char *> (&f1), unsignedlongint_size);	
			infile.read (reinterpret_cast<char *> (&f2), unsignedlongint_size);
			infile.ignore(skip_1);
			infile.read (reinterpret_cast<char *> (&fi), int_size);		// RT_min
			if(fi > max_RT_ind){
				infile.ignore(skip_2);
				continue;
			}
			infile.read (reinterpret_cast<char *> (&fi), int_size);		// RT_max
			if(fi < min_RT_ind){
				infile.ignore(skip_3);
				continue;
			}	
			infile.read (reinterpret_cast<char *> (&store_data_double), double_size);		// mass_min	
			if(store_data_double > max_mass){
				infile.ignore(skip_4);
				continue;
			}	
			infile.read (reinterpret_cast<char *> (&store_data_double), double_size);		// mass_max	
			if(store_data_double < min_mass){
				continue;
			}		
			// if not continued, mark the centroids in this partition
			f1 = (((f1 - header_entries_start) / centroid_size));
			f2 = (((f2 - header_entries_start) / centroid_size));
			for(m = f1; m <= f2; m++) keep_centroids(m) = 2; // mark centroids as candidates
		}	
	}else{
		for(n = 0; n < MSlist_centroids_nrow; n++) keep_centroids(n) = 2;		
	}
	// search centroids
	infile.seekg (header_entries_start, ios::beg);
	skip_1 = (2 * double_size + 2 * int_size); 
	skip_2 = (double_size + 2 * int_size); 
	for(n = 0; n < MSlist_centroids_nrow; n++){	
		if(keep_centroids(n) != 2){
			infile.ignore(skip_1);
			continue;			
		}
		infile.read (reinterpret_cast<char *> (&store_data_double), double_size);
		if((store_data_double < min_mass) || (store_data_double > max_mass)){
			infile.ignore(skip_2);
			continue;
		}
		infile.ignore(double_size); // skip intensity
		infile.read (reinterpret_cast<char *> (&store_data_int), int_size);		
		if((store_data_int < min_RT_ind) || (store_data_int > max_RT_ind)){ 
			infile.ignore(int_size);
			continue;
		}
		keep_centroids(n) = 1;
		number_keep_centroids ++;
		infile.ignore(int_size);	
	}	
	// write
	NumericMatrix reading(number_keep_centroids, 4);
	at_entry = 0;
	m = 0;	
	infile.seekg (header_entries_start, ios::beg);
	for(n = 0; n < MSlist_centroids_nrow; n++){
		if(keep_centroids(n) == 1){
			infile.ignore(centroid_size * at_entry); // shift by previously omitted positions
			infile.read (reinterpret_cast<char *> (&reading(m, 0)), double_size);
			infile.read (reinterpret_cast<char *> (&reading(m, 1)), double_size);
			infile.read (reinterpret_cast<char *> (&RT), int_size);		
			infile.read (reinterpret_cast<char *> (&measureID), int_size);	
			reading(m, 2) = (double) RT;
			reading(m, 3) = (double) measureID;
			m++;
			at_entry = 0;
		}else{
			at_entry++;
		}
	}
	// insert RTs ---------------------------------------------------- //
	if(insert_RT){
		for(n = 0; n < number_keep_centroids; n++){		
			reading(n, 2) = readScans(((int) reading(n, 2)) -1);
		}
	}
	colnames(reading) = CharacterVector::create("m/z", "intensity", "RT", "measureID");
	infile.close();	
	
	return reading;
}







