/*
* Wrapped code for the package
*/


#include "R_thresholding_driver.h"

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
