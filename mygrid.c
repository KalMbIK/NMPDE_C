//
// Created by pusheen on 30.10.17.
//

#include "mygrid.h"
#include <stdlib.h>

Array * getGrid(double a, double b, int N){
    double step = (b-a)/N;
    Array * gr = newArray(N);
    for (int i = 1; i < N; i++){
        gr->ar[i]=gr->ar[i-1]+step;
    }
    return gr;
}

double* vectorizeDoubleFun(double (*f)(double, double), double* gr, int N, double t){
    double * result = (double*)calloc(N, sizeof(double));
    for (int i = 0; i <= N; i++){
        result[i]=f(gr[i],t);
    }
    return result;
}
