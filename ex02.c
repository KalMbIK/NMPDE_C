//#define TYPE int
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "myNumPy.h"

int N = 81;
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

void TDA_solver(const double *a, const double *b, const double *c, double *cPrime, const double *d, double *dPrime, double *result, int size) {
    cPrime[0] = c[0]/b[0];
    dPrime[0] = d[0]/b[0];
    for (int i = 1; i < size; i++){
        if (i != size-1)
            cPrime[i] = c[i]/(b[i]-a[i]*cPrime[i-1]);
        dPrime[i] = (d[i] - a[i]*dPrime[i-1])/(b[i]-a[i]*cPrime[i-1]);
    }
    result[size-1] = dPrime[size-1];
    for (int i = size-2; i >= 0; i--)
        result[i] = dPrime[i] - cPrime[i]*result[i+1];
}

void generateRightHandSide(const double* a, const double* b, const double* c, const double* u, double* d, int size){
    d[0] = b[0]*u[0] + c[0] * u[1];
    d[size-1] = a[size-1]*u[size-2]+b[size-1]*u[size-1];
    for (int i = 1; i < size-2; i++)
        d[i] = a[i]*u[i-1] + b[i]*u[i] + c[i]*u[i+1];
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
        generateRightHandSide(aRight, bRight, cRight, u, d, size);
        TDA_solver(a, b, c, cPrime, d, dPrime, v, size);
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

double getError(double* a, double* b, int size){
    double * delta = subtr(a,b,size);
    double * absDelta = vectorize(&fabs,delta,size);
    double max = getMaxElement(absDelta,size);
    free(absDelta);
    free(delta);
    return max;
}

int main(void){
    double h = (x1-x0)/(N-1);
    double t = h;
    double la = a_diff*t/h/h;
    double * gr = linspace(x0,x1,N);
    double * u = vectorize(&u0,gr,N);
    double * v = newArray(N);
    double * exact = vectorize(&exactSolution,gr,N);

    solver(la,t,&u[1],&v[1],N-2);

    printf("Error = %g", getError(v,exact,N));

    free(gr);
    free(u);
    free(v);
    free(exact);
    return 0;
}