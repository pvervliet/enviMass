# include <Rcpp.h>

using namespace Rcpp;

List resize( // to smaller or larger Lists
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

bool check_plaus(
	NumericMatrix check_nodes_sub,
	NumericMatrix pattern_compound,
	NumericMatrix peaks,
	double RT_tol_inside,
	double int_tol
){

	int n_rows = check_nodes_sub.nrow();
	if(n_rows == 1) return true;
	if(is_true(any(duplicated(check_nodes_sub( _ , 0))))) return false; // unique centroids combined?
	if(is_true(any(duplicated(check_nodes_sub( _ , 1))))) return false; // unique peaks combined - should always be the case, not?
	// all within small RT window ?
	int n, m;
	double RT, RT_min, RT_max, ratio_int, ratio_int_theo_low, ratio_int_theo_high;
	RT_min = peaks( (check_nodes_sub(0, 1) -1) , 2);
	RT_max = RT_min;
	if(n_rows > 1){
		for(m = 1; m < n_rows; m++){
			RT = peaks( (check_nodes_sub(m, 1) -1) , 2);
			if( RT	> RT_max ){
				RT_max = RT;
			}else{
				if( RT < RT_min ){
					RT_min = RT;
				}
			}
		}
	}
	if((RT_max - RT_min) > RT_tol_inside) return false;
	// does intensity pattern match?
	if((int_tol < 100) && (n_rows > 1)){
		for(n = 0; n < (n_rows-1); n++){
			for(m = (n+1); m < n_rows; m++){
				ratio_int = ( peaks( (check_nodes_sub(m, 1) -1) , 1) / peaks( (check_nodes_sub(n, 1) -1) , 1) );
				ratio_int_theo_high = (
					( pattern_compound( (check_nodes_sub(m, 0) -1) , 1) + ( pattern_compound( (check_nodes_sub(m, 0) -1) , 1) * int_tol / 100)) /
					( pattern_compound( (check_nodes_sub(n, 0) -1) , 1) - ( pattern_compound( (check_nodes_sub(n, 0) -1) , 1) * int_tol / 100))
				);
				if(ratio_int_theo_high < ratio_int) return false;
				ratio_int_theo_low = (
					( pattern_compound( (check_nodes_sub(m, 0) -1) , 1) - ( pattern_compound( (check_nodes_sub(m, 0) -1) , 1) * int_tol / 100)) /
					( pattern_compound( (check_nodes_sub(n, 0) -1) , 1) + ( pattern_compound( (check_nodes_sub(n, 0) -1) , 1) * int_tol / 100))
				);
				if(ratio_int_theo_low > ratio_int) return false;
			}
		}
	}
	// else -> all fine
	return true;
}

int get_row_size(
	NumericMatrix check_nodes_sub
){

	int n = check_nodes_sub.nrow();
	return n;
}

NumericMatrix drop_row(
	NumericMatrix check_nodes_sub,
	int skip_row
){

	int n_rows = check_nodes_sub.nrow(), n, m_cols = check_nodes_sub.ncol(), m, at_row = 0;
	NumericMatrix new_matrix( (n_rows-1) , 2 );
	for(n = 0;  n < n_rows; n++){
		if(n == skip_row) continue; // skip that row
		for(m = 0; m < m_cols; m++){
			new_matrix(at_row, m) = check_nodes_sub(n, m);
		}
		at_row++;
	}
	return new_matrix;
}

NumericVector return_matrix_column(
	NumericMatrix x,
	int col_num
){

	return x(_,col_num);
}

bool result_exists(
	NumericMatrix check_nodes_sub,
	List results_peaks,
	int at_size
){

	if(at_size == 0) return false;
	int n;
	NumericVector sub_vector, check_vector;
	check_vector = check_nodes_sub(_,1);
	for(n = 0; n < at_size; n++){ // requires above return false statement!
		sub_vector = return_matrix_column(results_peaks[n], 1);
		if(is_false(any(is_na(match(check_vector, sub_vector))))) return true;
	}
	return false;
}

// [[Rcpp::export]]
List while_checked(
	List check_nodes,
	NumericMatrix pattern_compound,
	NumericMatrix peaks,
	double RT_tol_inside,
	double int_tol
){

	int list_size_results = 10, list_size_new_nodes = 10, at_size = 0, at_new_nodes = 0, m, n, old_size, check_nodes_size;
	bool checked = true, get_plaus = false, verbose = false;
	List results_peaks(list_size_results); // variable size
	List results(list_size_results); // variable size
	List new_nodes(list_size_new_nodes); // variable size

	while(checked){
        Rcpp::checkUserInterrupt();
        if(verbose) Rprintf( "\n NEW - " );
		at_new_nodes = 0;
		checked = false;
		check_nodes_size = check_nodes.size();
		for(m = 0; m < check_nodes_size; m++){
            if(verbose) Rprintf( " * " );
			// contained in results?
			if(at_size > 0){
                if(result_exists(check_nodes[m], results_peaks, at_size)) continue;
			}
            // contained in previous combinations?
            if(m > 0){
                if(result_exists(check_nodes[m], check_nodes, m - 1)) continue;
            }
			// < complete contained section!
            if(verbose) Rprintf( " . " );
			// check: peak combination is plausible?
			get_plaus = check_plaus(
				check_nodes[m],
				pattern_compound,
				peaks,
				RT_tol_inside,
				int_tol
			);
			if(get_plaus){ // write to results
                if(verbose) Rprintf( " + " );
				// save to results, eventually increase its size
				results[at_size] = check_nodes[m];
				results_peaks[at_size] = check_nodes[m];
				at_size++;
				if(at_size >= list_size_results){
					list_size_results = (list_size_results + 10);
					results = resize(results, list_size_results);
					results_peaks = resize(results_peaks, list_size_results);
				}
			}else{	// re-arrange: drop a line at a time
                if(verbose) Rprintf( " -" );
				old_size = get_row_size(check_nodes[m]);
				for(n = 0; n < old_size; n++){
                    if(verbose) Rprintf( "," );
					new_nodes[at_new_nodes] = drop_row(check_nodes[m], n);
					at_new_nodes++;
					if(at_new_nodes >= list_size_new_nodes){
						list_size_new_nodes = (list_size_new_nodes + 10);
						new_nodes = resize(new_nodes, list_size_new_nodes);
					}
				}
                if(verbose) Rprintf( " " );
				checked = true;
			}
		}
		if(checked){
			check_nodes = resize(new_nodes, at_new_nodes);
		}
	}

    if(verbose) Rprintf( " \n\n " );
	results = resize(results, at_size);
    return results;

}


