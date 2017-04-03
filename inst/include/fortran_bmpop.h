
#include <Rcpp.h>

extern "C"
{
    void bmpopsub_(double *, int, int, double *,double *,double *, int,
                   double);
    
    void test_(double *, int, int);
}