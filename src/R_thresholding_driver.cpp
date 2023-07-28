/*
* Wrapped code for the package
*/


#include "R_thresholding_driver.h"

/* C Api Route --> Will wait and see if this is needed
SEXP R_thresh_analysis (SEXP outfile_prefix, SEXP G, SEXP l, 
                        SEXP u, SEXP increment, SEXP windowSize, 
                        SEXP minimumpartitionsize, SEXP minimum_cliquesize, 
                        SEXP min_alpha, SEXP max_alpha, SEXP alpha_increment, 
                        SEXP num_samples, SEXP significance_alpha, 
                        SEXP bonferroni_corrected, SEXP methods 
                        ){
    SEXP R_data;

    R_data = 1;
    
    return R_data;
}
*/

// [[Rcpp::export]]
int Rcpp_thresh_analysis(std::string& outfile_prefix,
                        igraph_t &G,
                        double l,
                        double u,
                        double increment,
                        int windowsize,
                        int minimumpartitionsize,
                        int minimum_cliquesize,
                        double min_alpha,
                        double max_alpha,
                        double alpha_increment,
                        int num_samples,
                        double significance_alpha,
                        double bonferroni_corrected,
                        std::set<int> methods){

  Rcpp::Rcout << "Function Test" << "\n";
}

