

#include <Rcpp.h>


#include "fortran_bmpop.h"


using namespace Rcpp;

//' @export
// [[Rcpp::export]]
SEXP cpp_bmpop(
        NumericMatrix data, NumericMatrix sig, 
        NumericMatrix rytpar,
        double alpha, int k, double perd
    ){
    
    bmpopsub_(data.begin(), data.nrow(), data.ncol(), sig.begin(), rytpar.begin(),
             &alpha, k, perd);
    
    return List::create(Named("sig") = sig,
                        Named("rytpar") = rytpar);
}

//' @export
// [[Rcpp::export]]
SEXP cpp_test(NumericMatrix X){
    
    int nr = X.nrow();
    int nc = X.ncol();
    
    test_(X.begin(), nr, nc);
    
    Rcpp::Rcout << X << std::endl;
    
    return X;
}



