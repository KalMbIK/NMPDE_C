#include <math.h>
#include <stdio.h>
#include "mygrid.h"
#define TYPE int
#include "myarrays.h"

double exactSolution(double x, double t){
    return sin(M_PI*x)*exp(-M_PI*M_PI*t);
}

int main() {
    int N = 11, M = 11;
    double x0=0, x1 = 1;
    double T = 0.2;
    double a = 1;
    double h = (x1-x0)/(N-1);
    double t = T/(M-1);
    double la = a*t/h/h;
    Array * ar = newArray(N);
    printArray(ar);
    ar->g(0) = 1;
    printf("%lg",ar->g(0));
//    set(ar,1,1);
//    double* gr = getGrid(x0,x1,N);
//    printArray(gr,N+1);
    double * f = &exactSolution;
    return 0;
}