#include <Rcpp.h>

#include "math_ext.h"
#include "spectral_methods.h"
#include <igraph.h>

//' Driver code for primary graph analysis
//' ADD INFO FOR PARAMS
// [[Rcpp::export]]
int thresholdAnalysis(double thing, int other_thing, int optional = 0){
    Rcpp::Rcout << "THIS IS THE OPTIONAL PARAM: " << optional << '\n';

    return 0;

}