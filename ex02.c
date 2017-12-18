//
// Created by pusheen on 30.10.17.
//
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "myNumPy.h"
#include "myTDMAlgorithms.h"


//IMPORTANT: I used notation from the Wikipedia to declare the TDM

int N = 200;
double x0=0, x1 = 1;
double T = 0.2;
double a_diff = 1;

double u0(double x){
    return sin(M_PI*x);
}

double f2(double t){
    return exp(-M_PI*M_PI*t);
}

double exactSolution(double x){
    return u0(x)*f2(T);
}

void genCoeffs(double lambda, double* a, double* b, double* c, int size){
    double A = -lambda/2;
    double C = A;
    double B = 1 + lambda;
    a[0] = 0;
    b[0] = B;
    c[0] = C;
    a[size-1] = A;
    b[size-1] = B;
    c[size-1] = 0;
    for (int i = 1; i < size-1; i++){
        b[i] = B;
        c[i] = C;
        a[i] = A;
    }
}

void solver(double lambda, double t, double* u0, double * res, int size){
    double * a = newArray(size);
    double * b = newArray(size);
    double * c = newArray(size);
    double * aRight = newArray(size);
    double * bRight = newArray(size);
    double * cRight = newArray(size);
    double * d = newArray(size);
    double * u = res;
    double * v = newArray(size);
    double * dPrime = newArray(size);
    double * cPrime = newArray(size);
    size_t nbytes = size*sizeof(double);

    genCoeffs(lambda, a, b, c, size);
    genCoeffs(-lambda, aRight, bRight, cRight, size);
    memcpy(u, u0, nbytes);

    double tt = 0;
    while (tt < T){
        TDM_dot(aRight, bRight, cRight, u, d, size);
        TDMA_solver(a, b, c, cPrime, d, dPrime, v, size);
        tt += t;
        double* temp = u;
        u = v;
        v = temp;
    }

    res = u;

    free(a);
    free(b);
    free(c);
    free(aRight);
    free(bRight);
    free(cRight);
    free(d);
    free(v);
    free(dPrime);
    free(cPrime);
}

int main(void){
    double h = (x1-x0)/(N);
    double t = h;
    double la = a_diff*t/h/h;
    double * gr = linspace(x0,x1,N+1);
    double * u = vectorize(&u0,gr,N+1);
    double * v = newArray(N+1);
    double * exact = vectorize(&exactSolution,gr,N+1);

    solver(la,t,&u[1],&v[1],N-1);

    printf("Error = %g", getError(v,exact,N+1));

    free(gr);
    free(u);
    free(v);
    free(exact);
    return 0;
}