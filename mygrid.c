//
// Created by pusheen on 30.10.17.
//

#include "mygrid.h"
#include <stdlib.h>

TYPE* linspace(TYPE a, TYPE b, int N) {
    TYPE step = (b-a)/(N-1);
    TYPE* gr = newArray(N);
    gr[0] = a;
    for (int i = 1; i < N; i++){
        gr[i]=gr[i-1]+step;
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
