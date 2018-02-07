# include <iostream>
# include <Rcpp.h>
# include <fstream>

using namespace std;
using namespace Rcpp;


// -> Initialize a MSraw file: to contain file checks?
// [[Rcpp::export]]
bool init_MSraw( 
	StringVector file_path
){
	int unsignedlongint_size = sizeof (unsigned long int), f = 0, n;
	
	ofstream outfile;
    outfile.open(file_path(0), ios::binary | ios::out | ios::trunc); // truncate -> init new file
	if(!outfile) return false;
	for(n = 0; n < 10; n++){
		outfile.write (reinterpret_cast<char *> (&f), unsignedlongint_size);
	}
	outfile.close ();
	
	return true;
}


// -> Write centroids into MSraw file
// [[Rcpp::export]]
bool write_centroids_MSraw(
	StringVector file_path,
	NumericVector MSlist_scans,
	NumericMatrix MSlist_centroids,
	IntegerVector match_scans_RT
){
	
	int n, measureID, partID, clustID, peakID, profDataID = 0;
	int int_size = sizeof (int), unsignedlongint_size = sizeof (unsigned long int), double_size = sizeof (double);
	int MSlist_scans_size = MSlist_scans.size();
	int MSlist_centroids_nrow = MSlist_centroids.nrow();
	unsigned long int header_entry_1, header_entry_2, header_entry_3;
	
	ofstream outfile;
    outfile.open(file_path(0), ios::binary | ios::out);
	if(!outfile) return false;
	// first three entries to MSraw header:
	// index start of scan RT data block = (number of bytes - 1), fixed sized
	header_entry_1 = (10 * unsignedlongint_size); 
	outfile.seekp (0 * unsignedlongint_size, ios::beg); 
	outfile.write (reinterpret_cast<char *> (&header_entry_1), unsignedlongint_size);
	// index start of centroid data block = (number of bytes - 1), variably sized
	header_entry_2 = (header_entry_1 + (MSlist_scans_size * double_size));
	outfile.seekp (1 * unsignedlongint_size, ios::beg);
	outfile.write (reinterpret_cast<char *> (&header_entry_2), unsignedlongint_size);	
	// index start of partition data block = (number of bytes - 1), variably sized
	header_entry_3 = (header_entry_2 + (MSlist_centroids_nrow * 2 * double_size) + (MSlist_centroids_nrow * 6 * int_size));
	outfile.seekp (2 * unsignedlongint_size, ios::beg); 
	outfile.write (reinterpret_cast<char *> (&header_entry_3), unsignedlongint_size);	
	// entries to scan RT data block:
	outfile.seekp (header_entry_1, ios::beg); 
	for(n = 0; n < MSlist_scans_size; n++){
		outfile.write (reinterpret_cast<char *> (&MSlist_scans(n)), double_size);	
	}
	// entries to centroid data block:
	outfile.seekp (header_entry_2, ios::beg);
	for(n = 0; n < MSlist_centroids_nrow; n++){	
		outfile.write (reinterpret_cast<char *> (&MSlist_centroids(n, 0)), double_size);	// m/z
		outfile.write (reinterpret_cast<char *> (&MSlist_centroids(n, 1)), double_size);	// intensity
		outfile.write (reinterpret_cast<char *> (&match_scans_RT(n)), int_size); 				// RT index
		measureID = (int) MSlist_centroids(n, 3);
		outfile.write (reinterpret_cast<char *> (&measureID), int_size);	// measureID
		partID = (int) MSlist_centroids(n, 4); 
		outfile.write (reinterpret_cast<char *> (&partID), int_size);		// partID
		clustID = (int) MSlist_centroids(n, 5); 
		outfile.write (reinterpret_cast<char *> (&clustID), int_size);		// clustID
		peakID = (int) MSlist_centroids(n, 6); 
		outfile.write (reinterpret_cast<char *> (&peakID), int_size);		// peakID
		//profDataID = placeholder for future linking to profile raw data
		outfile.write (reinterpret_cast<char *> (&profDataID), int_size);		// peakID		
	}	
	outfile.close ();

	return true;
}


// -> Write partition index to MSraw file
// [[Rcpp::export]]
bool write_partition_MSraw(
	StringVector file_path,
	NumericMatrix partition_index,
	bool add_ranges
){

	unsigned long int n, m, header_entry_2, header_entry_3, header_entry_4, d, partition_index_size = (unsigned long int) partition_index.nrow(); 
	int int_size = sizeof (int), unsignedlongint_size = sizeof (unsigned long int), double_size = sizeof (double);	
	int partition_size = (3 * unsignedlongint_size + 2 * int_size + 2 * double_size); // RT range set as int -> cp. centroid block -> cp. scans block
	unsigned long int start_centroid, end_centroid, centroid_size = (unsigned long int) (2 * double_size + 6 * int_size);
	double fi = 0, minRT, maxRT, minmass, maxmass, get_mass, get_RT;
	int fd = 0, skip_1 = (5 * int_size);
	
	fstream inoutfile;
    inoutfile.open(file_path(0), ios::binary | ios::out | ios::in); 
	if(!inoutfile) return 0;
	// index start of centroid data block = (number of bytes - 1), variably sized	
	inoutfile.seekg (1 * unsignedlongint_size, ios::beg);
	inoutfile.read (reinterpret_cast<char *> (&header_entry_2), unsignedlongint_size);
	// index start of partition data block = (number of bytes - 1), variably sized
	inoutfile.seekg (2 * unsignedlongint_size, ios::beg);
	inoutfile.read (reinterpret_cast<char *> (&header_entry_3), unsignedlongint_size);
	// index start of EIC data block = (number of bytes - 1), variably sized
	header_entry_4 = (header_entry_3 + partition_index_size * partition_size);
	inoutfile.seekp (3 * unsignedlongint_size, ios::beg);
	inoutfile.write (reinterpret_cast<char *> (&header_entry_4), unsignedlongint_size);	
	inoutfile.seekp (header_entry_3, ios::beg);	
	for(n = 0; n < partition_index_size; n++){
		start_centroid = (header_entry_2 + ((((unsigned long int) partition_index(n, 0)) - 1) * centroid_size));
		inoutfile.write (reinterpret_cast<char *> (&start_centroid), unsignedlongint_size);	// start_ID of first centroid in partition
		end_centroid = (header_entry_2 + ((((unsigned long int) partition_index(n, 1)) - 1) * centroid_size));
		inoutfile.write (reinterpret_cast<char *> (&end_centroid), unsignedlongint_size);	// start_ID of last centroid in partition
		d = (unsigned long int) partition_index(n, 2);
		inoutfile.write (reinterpret_cast<char *> (&d), unsignedlongint_size);	// number_peaks
		if(add_ranges){
			minRT = 1E10, maxRT = 0, minmass = 1E10, maxmass = 0;
			inoutfile.seekg (start_centroid, ios::beg);
			for(m = start_centroid; m < (end_centroid + 1); m++){ // centroids must all be in one block
				inoutfile.read (reinterpret_cast<char *> (&get_mass), double_size);
				inoutfile.ignore(double_size); // skip intensity
				inoutfile.read (reinterpret_cast<char *> (&get_RT), int_size);	
				if(get_mass < minmass) minmass = get_mass;
				if(get_mass > maxmass) maxmass = get_mass;
				if(get_RT < minRT) minRT = get_RT;
				if(get_RT > maxRT) maxRT = get_RT;		
				inoutfile.ignore(skip_1); // skip centroid membership indices			
			}
			inoutfile.seekp ((header_entry_3 + (n * partition_size) + (3 * unsignedlongint_size)), ios::beg);	
			inoutfile.write (reinterpret_cast<char *> (&minRT), int_size);	//		
			inoutfile.write (reinterpret_cast<char *> (&maxRT), int_size);	//	
			inoutfile.write (reinterpret_cast<char *> (&minmass), double_size);	//	
			inoutfile.write (reinterpret_cast<char *> (&maxmass), double_size);	//				
		}else{	// set to zero
			inoutfile.write (reinterpret_cast<char *> (&fi), int_size);	//		
			inoutfile.write (reinterpret_cast<char *> (&fi), int_size);	//	
			inoutfile.write (reinterpret_cast<char *> (&fd), double_size);	//	
			inoutfile.write (reinterpret_cast<char *> (&fd), double_size);	//				
		}
	}
	inoutfile.close ();

	return true;
}


// -> Read partition matrix from MSraw file
// [[Rcpp::export]]
NumericVector read_partition_MSraw(
	StringVector file_path,
	bool insert_RT,
	bool index_convert
){
	
	unsigned long int n, fu, header_entry_2, header_entry_3, header_entry_4, partition_index_size, header_entries_start, header_entries_end, MSlist_scans_size;
	int fi, int_size = sizeof (int), unsignedlongint_size = sizeof (unsigned long int), double_size = sizeof (double);	
	int partition_size = (3 * unsignedlongint_size + 2 * int_size + 2 * double_size); 
	double fd;
	unsigned long int centroid_size = (unsigned long int) (2 * double_size + 6 * int_size);
	
	ifstream infile;
    infile.open(file_path(0), ios::binary | ios::in); 	
	if(!infile) return 0;
	infile.seekg (1 * unsignedlongint_size, ios::beg);
	infile.read (reinterpret_cast<char *> (&header_entry_2), unsignedlongint_size);
	infile.read (reinterpret_cast<char *> (&header_entry_3), unsignedlongint_size);
	infile.read (reinterpret_cast<char *> (&header_entry_4), unsignedlongint_size);
	partition_index_size = ((header_entry_4 - header_entry_3) / partition_size);
	NumericMatrix reading(partition_index_size, 7);	
	infile.seekg (header_entry_3, ios::beg);		
	for(n = 0; n < partition_index_size; n++){		
		infile.read (reinterpret_cast<char *> (&fu), unsignedlongint_size);	
		if(index_convert) fu = (((fu - header_entry_2) / centroid_size) + 1);
		reading(n, 0) = (double) fu;
		infile.read (reinterpret_cast<char *> (&fu), unsignedlongint_size);	
		if(index_convert) fu = (((fu - header_entry_2) / centroid_size) + 1);
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
	//if(insert_RT){	
		infile.seekg (0 * unsignedlongint_size, ios::beg);
		infile.read (reinterpret_cast<char *> (&header_entries_start), unsignedlongint_size);
		infile.read (reinterpret_cast<char *> (&header_entries_end), unsignedlongint_size);
		MSlist_scans_size = ((header_entries_end - header_entries_start) / double_size);
		NumericVector readScans(MSlist_scans_size);
		infile.seekg (10 * unsignedlongint_size, ios::beg);
		for(n = 0; n < MSlist_scans_size; n++){		
			infile.read (reinterpret_cast<char *> (&readScans[n]), double_size);
		}
		for(n = 0; n < partition_index_size; n++){				
			//reading(n, 3) = readScans(((int) reading(n, 3)) -1);
			//reading(n, 4) = readScans(((int) reading(n, 4)) -1);			
		}
	//}	
	infile.close();
	
	//return reading;	
	return readScans;
}

// -> Read header of MSraw file
// WARNING: conversion of (unsigned long int) to R IntegerVector (int)
// [[Rcpp::export]]
IntegerVector read_header_MSraw(
	StringVector file_path
){
	
	int n, unsignedlongint_size = sizeof (unsigned long int);
	unsigned long int header_entries;
	IntegerVector reading(10);
	
	ifstream infile;
    infile.open(file_path(0), ios::binary | ios::in); 	
	if(!infile){
		reading(0) = 0;
		return reading;
	}
	infile.seekg (0, ios::beg);
	for(n = 0; n < 10; n++){	
		infile.read (reinterpret_cast<char *> (&header_entries), unsignedlongint_size);	
		reading(n) = (int) header_entries;
	}
	infile.close();	
	
	return reading;
}


// -> Read scan RTs from MSraw file
// [[Rcpp::export]]
NumericVector read_scans_MSraw(
	StringVector file_path
){

	int double_size = sizeof (double), unsignedlongint_size = sizeof (unsigned long int);
	unsigned long int n, header_entries_start, header_entries_end, MSlist_scans_size;
	
	ifstream infile;
    infile.open(file_path(0), ios::binary | ios::in); 	
	if(!infile) return 0;

	infile.seekg (0 * unsignedlongint_size, ios::beg);
	infile.read (reinterpret_cast<char *> (&header_entries_start), unsignedlongint_size);
	infile.read (reinterpret_cast<char *> (&header_entries_end), unsignedlongint_size);
	MSlist_scans_size = ((header_entries_end - header_entries_start) / double_size);
	NumericVector reading(MSlist_scans_size);
	infile.seekg (10 * unsignedlongint_size, ios::beg);
	for(n = 0; n < MSlist_scans_size; n++){		
		infile.read (reinterpret_cast<char *> (&reading[n]), double_size);
	}
	infile.close();	

	return reading;
}


// -> Read all centroids from MSraw file
// [[Rcpp::export]]
NumericMatrix read_centroids_MSraw(
	StringVector file_path,
	bool insert_RT
){

	int int_size = sizeof (int), double_size = sizeof (double), unsignedlongint_size = sizeof (unsigned long int);
	int RT, measureID, partID, clustID, peakID, profDataID;
	unsigned long int n, header_entries_start, header_entries_end, MSlist_centroids_nrow, MSlist_scans_size;
	
	ifstream infile;
    infile.open(file_path(0), ios::binary | ios::in); 
	if(!infile) return 0;
	infile.seekg (1 * unsignedlongint_size, ios::beg);
	infile.read (reinterpret_cast<char *> (&header_entries_start), unsignedlongint_size);
	infile.read (reinterpret_cast<char *> (&header_entries_end), unsignedlongint_size);	
	MSlist_centroids_nrow = ((header_entries_end - header_entries_start) / (2 * double_size + 6 * int_size));
	NumericMatrix reading(MSlist_centroids_nrow , 7);
	
	infile.seekg (header_entries_start, ios::beg);
	for(n = 0; n < MSlist_centroids_nrow; n++){		
		infile.read (reinterpret_cast<char *> (&reading(n, 0)), double_size);
		infile.read (reinterpret_cast<char *> (&reading(n, 1)), double_size);
		infile.read (reinterpret_cast<char *> (&RT), int_size);		
		infile.read (reinterpret_cast<char *> (&measureID), int_size);	
		infile.read (reinterpret_cast<char *> (&partID), int_size);	
		infile.read (reinterpret_cast<char *> (&clustID), int_size);	
		infile.read (reinterpret_cast<char *> (&peakID), int_size);	
		infile.read (reinterpret_cast<char *> (&profDataID), int_size);	// dummy
		reading(n, 2) = (double) RT;
		reading(n, 3) = (double) measureID;
		reading(n, 4) = (double) partID;
		reading(n, 5) = (double) clustID;
		reading(n, 6) = (double) peakID;
	}
	// insert RTs
	if(insert_RT){
		infile.seekg (0 * unsignedlongint_size, ios::beg);
		infile.read (reinterpret_cast<char *> (&header_entries_start), unsignedlongint_size);
		infile.read (reinterpret_cast<char *> (&header_entries_end), unsignedlongint_size);
		MSlist_scans_size = ((header_entries_end - header_entries_start) / double_size);
		NumericVector readScans(MSlist_scans_size);
		infile.seekg (10 * unsignedlongint_size, ios::beg);
		for(n = 0; n < MSlist_scans_size; n++){		
			infile.read (reinterpret_cast<char *> (&readScans[n]), double_size);
		}
		for(n = 0; n < MSlist_centroids_nrow; n++){		
			reading(n, 2) = readScans(((int) reading(n, 2)) -1);
		}
	}
	infile.close();	

	return reading;
}


// -> Filter centroids from MSraw file
// [[Rcpp::export]]
NumericMatrix filter_centroids_MSraw(
	StringVector file_path,
	int min_RT,
	int max_RT,
	int min_mass,
	int max_mass,
	bool insert_RT
){

	int int_size = sizeof (int), double_size = sizeof (double), unsignedlongint_size = sizeof (unsigned long int);
	int centroid_size = (2 * double_size + 6 * int_size), skip_1, skip_2;
	int min_RT_ind, max_RT_ind, store_data_int;
	int RT, measureID, partID, clustID, peakID, profDataID;
	unsigned long int m, n, header_entries_start, header_entries_end, at_entry, MSlist_centroids_nrow, MSlist_scans_size, number_keep_centroids = 0;
	double store_data_double;
	
	ifstream infile;
    infile.open(file_path(0), ios::binary | ios::in); 
	if(!infile) return 0;

	// -> via direct search ------------------------------------------ //
	// get RT index range
	infile.seekg (0 * unsignedlongint_size, ios::beg);
	infile.read (reinterpret_cast<char *> (&header_entries_start), unsignedlongint_size);
	infile.read (reinterpret_cast<char *> (&header_entries_end), unsignedlongint_size);
	MSlist_scans_size = ((header_entries_end - header_entries_start) / double_size);
	NumericVector readScans(MSlist_scans_size);
	infile.seekg (10 * unsignedlongint_size, ios::beg);
	for(n = 0; n < MSlist_scans_size; n++){		
		infile.read (reinterpret_cast<char *> (&readScans[n]), double_size);
	}
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
	// search
	infile.seekg (1 * unsignedlongint_size, ios::beg);
	infile.read (reinterpret_cast<char *> (&header_entries_start), unsignedlongint_size);
	infile.read (reinterpret_cast<char *> (&header_entries_end), unsignedlongint_size);	
	MSlist_centroids_nrow = ((header_entries_end - header_entries_start) / (2 * double_size + 6 * int_size));
	IntegerVector keep_centroids(MSlist_centroids_nrow);	
	for(n = 0; n < MSlist_centroids_nrow; n++){		
		keep_centroids(n) = 0;
	}
	at_entry = (header_entries_start + 0 * centroid_size);	
	infile.seekg (at_entry, ios::beg);
	skip_1 = (centroid_size - double_size); 
	skip_2 = (centroid_size - (2 * double_size + int_size));
	for(n = 0; n < MSlist_centroids_nrow; n++){				
		infile.read (reinterpret_cast<char *> (&store_data_double), double_size);
		if((store_data_double < min_mass) || (store_data_double > max_mass)){
			infile.ignore(skip_1);
			continue;
		}
		infile.ignore(double_size); // skip intensity
		infile.read (reinterpret_cast<char *> (&store_data_int), int_size);		
		if((store_data_int < min_RT_ind) || (store_data_int > max_RT_ind)){ 
			infile.ignore(skip_2);
			continue;
		}
		keep_centroids(n) = 1;
		number_keep_centroids ++;
		infile.ignore(skip_2);	
	}
	// write
	NumericMatrix reading(number_keep_centroids, 7);
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
			infile.read (reinterpret_cast<char *> (&partID), int_size);	
			infile.read (reinterpret_cast<char *> (&clustID), int_size);	
			infile.read (reinterpret_cast<char *> (&peakID), int_size);	
			infile.read (reinterpret_cast<char *> (&profDataID), int_size);	// dummy
			reading(m, 2) = (double) RT;
			reading(m, 3) = (double) measureID;
			reading(m, 4) = (double) partID;
			reading(m, 5) = (double) clustID;
			reading(m, 6) = (double) peakID;
			m++;
			at_entry = 0;
		}else{
			at_entry++;
		}
	}
	// insert RTs
	if(insert_RT){
		for(n = 0; n < number_keep_centroids; n++){		
			reading(n, 2) = readScans(((int) reading(n, 2)) -1);
		}
	}

	
	// -> via search in partitions ----------------------------------- //
	// if ...
	
	
	
	
	
	// insert RTs ---------------------------------------------------- //
	
	infile.close();	
	
	return reading;
}
