//
// Created by pusheen on 18.01.2018.
//

#include <stdlib.h>
#include <math.h>
#include "myNumPyV2.h"
#include "myTDMAlgorithms.h"

void nonLinspace(double alpha, double* gr, size_t size){
    double div = exp(alpha)-1;
    int N = size-1;
    double cc = alpha/N;
    for (int i = 0; i < size; i++)
        gr[i]=(exp(cc*i)-1)/div;
}

void createStepsArray(double* gr, double* steps, size_t size){
    for (int i = 0; i < size-1; i++)
        steps[i]=gr[i+1]-gr[i];
}

void createTDM(double* invertSteps, size_t size, double* a, double* b, double* c){
    b[0] = invertSteps[0]+invertSteps[1];
    c[0] = -invertSteps[1];
    for (int i = 1; i < size-1; i++){
        b[i] = invertSteps[i]+invertSteps[i+1];
        c[i] = -invertSteps[i+1];
        a[i] = -invertSteps[i];
    }
}

double f(double x){
    return M_PI*M_PI*sin(M_PI*x);
}

double ex(double x){
    return sin(M_PI*x);
}

void solver(double *g, double * steps, double*u, size_t size){
    double *invertSteps = newArray(size-1);

    double *a = newArray(size-2);
    double *b = newArray(size-2);
    double *c = newArray(size-2);
    double *d = newArray(size-2);
    double *cPrime = newArray(size-2);
    double *dPrime = newArray(size-2);

    for (int i = 0; i < size-1; i++)
        invertSteps[i]=1/steps[i];

    for (int i = 0; i < size-2; i++)
        d[i] = .5*(g[i+1]*steps[i]+g[i+2]*steps[i+1]);

    createTDM(invertSteps,size-1,a,b,c);

    TDMA_solver(a,b,c,cPrime,d,dPrime,u,size-2);

    free(a);
    free(b);
    free(c);
    free(d);
    free(cPrime);
    free(dPrime);
    free(invertSteps);
}

double calculateError(double* steps, double* u, double* exact, size_t size){
    double * temp = newArray(size);
    double res = 0;
    subtr(u,exact,temp,size);
    for (int i = 0; i < size; i++)
        temp[i] *= temp[i];

    for (int i = 0; i < size-1; i++)
        res += .5*steps[i]*(temp[i+1]+temp[i]);

    free(temp);
    return res;
}

int main(void){

    size_t N = 1000;
    double al = 10;
    double *gr = newArray(N+1);
    double *steps = newArray(N);
    double *g = newArray(N+1);
    double *u = newArray(N+1);
    double *exact = newArray(N+1);
    nonLinspace(al,gr,N+1);
    createStepsArray(gr,steps,N+1);
    vectorize(&f,gr,g,N+1);
    vectorize(&ex,gr,exact,N+1);
    solver(g,steps,&u[1],N+1);

    printf("STEP=%lg\n",getMaxElement(steps,N));
    printf("ERROR=%lg\n",calculateError(steps,u,exact,N+1));
    free(gr);
    free(g);
    free(steps);
    free(exact);
    free(u);

    return 0;
}