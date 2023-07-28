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

//#include <R.h>
//#include <Rinternals.h>
#include <Rcpp.h>

#include <igraph/igraph.h>

#include <vector>     // std::vector
#include <iostream>   // std::cout, std::cerr, std::endl
#include <fstream>    // fopen, fclose (to read igraph), ofstream
#include <algorithm>  // std::nth_element, std::min_element, std::max_element
#include <math.h>     // pow, sqrt, fabs, M_PI
#include <getopt.h>   // commandline argument parsing
#include <stdlib.h>   // atoi, atof
#include <sstream>    // stringstream
#include <cmath>      // std::copysign
#include <math.h>     // isnormal
#include <set>        // sets
#include <string>     // getline
#include <csignal>    // signal

//   *****************************************************
// ***TEST THIS LATER --->> #include <igraph/igraph.h>   ***
//  *********************************************************

/* Utility functions used through programs */
// #include "utils.h"

/* C API Method - Going to try Rcpp as of 7/27/23
SEXP R_thresh_analysis(SEXP outfile_prefix, SEXP G, SEXP l, 
                        SEXP u, SEXP increment, SEXP windowSize, 
                        SEXP minimumpartitionsize, SEXP minimum_cliquesize, 
                        SEXP min_alpha, SEXP max_alpha, SEXP alpha_increment, 
                        SEXP num_samples, SEXP significance_alpha, 
                        SEXP bonferroni_corrected, SEXP methods 
                        ); 
*/

/* Trying out Rcpp 7/27/23, struggling to check package in R Studio*/
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
                        std::set<int> methods);
#endif