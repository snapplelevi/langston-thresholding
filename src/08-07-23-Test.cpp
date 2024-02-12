#include <Rcpp.h>

#include <igraph.h>
#include <vector>


//' Manually exported in NAMESPACE
int function_test(int c){
  Rcpp::Rcout << "This is a test of Rcpp" << '\n';
  igraph_t test;
  return c;
}

//' Prints out 10 numbers
//' @param i - an integer that gets printed right back out
//' Manually exported in NAMESPACE
void print_vector(int i){
  std::vector<int> tmp;
  Rcpp::Rcout << "your i value: " << i << '\n';
  for(int i = 0; i < 10; i++){
    tmp.push_back(i);
    Rcpp::Rcout << tmp[i] << '\n';
  }
  
}
