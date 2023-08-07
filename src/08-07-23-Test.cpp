#include <Rcpp.h>

// [[Rcpp::export]]
int function_test(int c){
  Rcpp::Rcout << "This is a test of Rcpp" << '\n';
  return c;
}