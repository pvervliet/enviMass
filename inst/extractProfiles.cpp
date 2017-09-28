# include <Rcpp.h>

using namespace Rcpp;

// resize a List to a smaller or larger List .................................//
List resize( 
	const List& a,
	int n
){
    int formersize = a.size(), use_size, i;
	List y(n);
	if(n < formersize){
		use_size = n;
	}else{
		use_size = formersize;
	}
    for(i = 0; i < use_size; i++) y[i] = a[i];
    return y;
}

// set RT windows for a cluster ...............................................//
NumericVector set_RT_window(
	NumericMatrix peaks,
	IntegerVector for_peaks,
	double dret
){
	NumericVector RT_window(2);
	RT_window[0] = R_PosInf;
	RT_window[1] = R_NegInf;
	int n, size = for_peaks.size();
	for(n = 0; n < size; n++){
		if( peaks(for_peaks[n], 2) < RT_window[0] ){RT_window[0] = peaks(for_peaks[n], 2);};
		if( peaks(for_peaks[n], 2) > RT_window[1] ){RT_window[1] = peaks(for_peaks[n], 2);};
	};
	RT_window[0] = (RT_window[0] - dret);
	RT_window[1] = (RT_window[1] + dret);	
	return RT_window;
}

// get tolerance windows for a mass ..........................................//
NumericVector set_mass_window(
	NumericMatrix peaks,
	IntegerVector for_peaks,
	bool ppm,
	double dmass
){
	NumericVector mass_lim(2);
	double mass, delmass;
	int n, size = for_peaks.size();	
	
	mass = peaks(for_peaks[0], 0);
	if(ppm){
		delmass = (dmass * mass / 1000000);
	}else{
		delmass = (dmass / 1000);
	};	
	mass_lim[0] = (mass - delmass);
	mass_lim[1] = (mass + delmass);	
	
	if(size > 1){
		for(n = 1; n < size; n++){
			mass = peaks(for_peaks[n], 0);
			if(ppm){
				delmass = (dmass * mass / 1000000);
			}else{
				delmass = (dmass / 1000);
			};	
			if( (mass + delmass) < mass_lim[1] ){mass_lim[1] = (mass + delmass);};
			if( (mass - delmass) > mass_lim[0] ){mass_lim[0] = (mass - delmass);};	
		}
	}
	
	return mass_lim;
}

// expand existing vector by one new scalar .................................//
IntegerVector expand_vector_scalar(
	IntegerVector cluster_with_peaks_vector,
	int at_peak
){
	int size = cluster_with_peaks_vector.size(), n;
	IntegerVector new_peak_vector(size + 1); 
	for(n = 0; n < size; n++){
		new_peak_vector[n] = cluster_with_peaks_vector[n];
	}
	new_peak_vector[size] = at_peak;
	return new_peak_vector;
}

// does a remaining peak fit to an existing cluster? .........................//
bool check_fits_cluster(
	NumericMatrix peaks,
	List cluster_with_peaks,
	List cluster_with_mass_lim,
	List cluster_with_RT_lim,
	int which_cluster,
	int new_peak,
	bool ppm,
	double dmass	
){
	// new_peak ranges within RT limits?
	NumericVector RT_lim = as<NumericVector>(cluster_with_RT_lim[which_cluster]);
	if( peaks(new_peak, 2) < RT_lim[0] ) return false;
	if( peaks(new_peak, 2) > RT_lim[1] ) return false;
	// new_peak overlaps with mass limits?
	double delmass = 0;
	if(ppm){
		delmass = (dmass * peaks(new_peak, 0) / 1000000);
	}else{
		delmass = (dmass / 1000);
	};	
	NumericVector mass_lim = as<NumericVector>(cluster_with_mass_lim[which_cluster]);
	if( (peaks(new_peak, 0) + delmass) < mass_lim[0] ) return false;
	if( (peaks(new_peak, 0) - delmass) > mass_lim[1] ) return false;	
	// already a peak with same sampleID in the cluster?
	IntegerVector peak_vector = as<IntegerVector>(cluster_with_peaks[which_cluster]);
	int n, size = peak_vector.size();
	for(n = 0; n < size; n++){
		if(peaks(peak_vector(n), 3) == peaks(new_peak, 3)){
			return false;
		}
	}
	// otherwise - all fits:
	return true;
}


// [[Rcpp::export]]
NumericMatrix extractProfiles_new(
	NumericMatrix peaks,
	IntegerVector in_order,
	int pregroup,
	double dmass,
	bool ppm,
	double dret,
	bool run_pregroup, 
	bool verbose
){

	int n_rows = peaks.nrow(), cluster_with_peaks_size = 10, cluster_with_peaks_at = 0, n = 0, m = 0;
	int peak_vector_size = 0, at_peak = 0, which_cluster = 0;
	bool fits_to_cluster = false;
	NumericMatrix clusters( n_rows, 15);
	std::fill(clusters.begin(), clusters.end(), 0);
	if(verbose & (n_rows == 1)) { Rprintf("Single rowed input?");}; 	
	List cluster_with_peaks(cluster_with_peaks_size); 	// contains the peaks of each cluster - variable size
	List cluster_with_mass_lim(cluster_with_peaks_size); 	// contains the peaks of each cluster - variable size
	List cluster_with_RT_lim(cluster_with_peaks_size); 		// contains the peaks of each cluster - variable size	
	IntegerVector peak_vector(1);
	IntegerVector one_peak_vector(1);
	// > DELETE	
	int checkpeak = 196;
	// < DELETE	

	// initiate first cluster .................................................... //
	one_peak_vector[0] = (in_order[0] - 1);
	cluster_with_peaks[cluster_with_peaks_at] = one_peak_vector;
	cluster_with_mass_lim[cluster_with_peaks_at] = set_mass_window(peaks, one_peak_vector, ppm, dmass);
	cluster_with_RT_lim[cluster_with_peaks_at] = set_RT_window(peaks, one_peak_vector, dret);
	cluster_with_peaks_at++;
	
	// cluster over remaining peaks .............................................. //
	if(n_rows > 1){
		for(n = 1; n < n_rows; n++){
			at_peak = (in_order[n] - 1);
			// > DELETE			
			Rprintf("\n next");
			if(at_peak == checkpeak) Rprintf("\n\n FOUND it!");
			// < DELETE			
			fits_to_cluster = false; // check: fits to any existing cluster? ......... //
			for(m = 0; m < cluster_with_peaks_at; m++){
				if( check_fits_cluster(peaks, cluster_with_peaks, cluster_with_mass_lim, cluster_with_RT_lim, m, at_peak, ppm, dmass) ) {
					fits_to_cluster = true;
					break;
				}
			}
			if(fits_to_cluster){ // add to existing cluster 
				// > DELETE			
				if(at_peak == checkpeak) Rprintf(" -> into existing! ");
				// < DELETE				
				peak_vector = as<IntegerVector>(cluster_with_peaks[m]);
				peak_vector = expand_vector_scalar(peak_vector, at_peak);
				cluster_with_peaks[m] = peak_vector;
				cluster_with_mass_lim[m] = set_mass_window(peaks, peak_vector, ppm, dmass);
				cluster_with_RT_lim[m] = set_RT_window(peaks, peak_vector, dret);			
			}else{ // create new cluster
				// > DELETE			
				if(at_peak == checkpeak) Rprintf(" -> into new! ");
				// < DELETE				
				// expand cluster list?
				if(cluster_with_peaks_at >= cluster_with_peaks_size){
					if(at_peak == checkpeak) Rprintf(" -> expand ");
					cluster_with_peaks_size = (cluster_with_peaks_size + 20);
					cluster_with_peaks = resize(cluster_with_peaks, cluster_with_peaks_size);
					cluster_with_mass_lim = resize(cluster_with_mass_lim, cluster_with_peaks_size);
					cluster_with_RT_lim = resize(cluster_with_RT_lim, cluster_with_peaks_size);		
				}		
				// add new cluster & get values
				one_peak_vector[0] = at_peak;
				cluster_with_peaks[cluster_with_peaks_at] = one_peak_vector;
				cluster_with_mass_lim[cluster_with_peaks_at] = set_mass_window(peaks, one_peak_vector, ppm, dmass);
				cluster_with_RT_lim[cluster_with_peaks_at] = set_RT_window(peaks, one_peak_vector, dret);		
				cluster_with_peaks_at++;
			}
		}
	}
	
	
	// write results to clusters output .......................................... //	
	for(n = 0; n < cluster_with_peaks_at; n++){
		peak_vector = as<IntegerVector>(cluster_with_peaks[n]);
		peak_vector_size = peak_vector.size();
		for(m = 0; m < peak_vector_size; m++){	
			clusters(peak_vector(m), 0) = (n+1); // insert cluster membership
			// > DELETE			
			clusters(peak_vector(m), 1) = peak_vector_size; // insert cluster membership
			// < DELETE
		}
	}
	
	
	return clusters;

}








