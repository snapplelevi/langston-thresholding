/* 
* R wrapper of Thresholding Code and Analysis Tools
* 
* Author: Levi Hochstetler (lhochste@vols.utk.edu)
* Created: July 2023
* 
Last Modified: July 2023
*/

#ifndef __R_THRESHOLDING_DRIVER__
#define __R_THRESHOLDING_DRIVER__

#include <R.h>
#include <Rinternals.h>
#include <stdio.h>
#include <stdlib.h>

/* Utility functions used through programs */
// #include "utils.h"

SEXP R_thresh_analysis(SEXP outfile_prefix, SEXP G, SEXP l, 
                        SEXP u, SEXP increment, SEXP windowSize, 
                        SEXP minimumpartitionsize, SEXP minimum_cliquesize, 
                        SEXP min_alpha, SEXP max_alpha, SEXP alpha_increment, 
                        SEXP num_samples, SEXP significance_alpha, 
                        SEXP bonferroni_corrected, SEXP methods 
                        );


#endif