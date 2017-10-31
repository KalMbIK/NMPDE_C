#define TYPE int
#include <math.h>
#include <stdio.h>
#include "mygrid.h"

int N = 11, M = 11;
double x0=0, x1 = 1;
double T = 0.2;
double a = 1;

double exactSolution(double x, double t){
    return sin(M_PI*x)*exp(-M_PI*M_PI*t);
}

double endSolution(double x){
    return exactSolution(x,T);
}

int main() {
    double h = (x1-x0)/(N-1);
    double t = T/(M-1);
    double la = a*t/h/h;
    int t1 = 0, t2 = 5, t3 = 6;
    TYPE* g = linspace(t1,t2,t3);
    printArray(g,N);
    printf("%d", sizeof(TYPE));
//    double * f =int &exactSolution;
    return 0;
}