//
// Created by pusheen on 29.11.17.
//
//#define TYPE int
#include "myNumPyV2.h"


//IMPORTANT: I used notation from the Wikipedia to declare the TDM

int N = 100;
double x0=0, x1 = 1;
double T = 1;

double u0(double x){
    return 1+sin(2*M_PI*x)/2;
}

double getError(double* a, double* b, int size){
    double * delta = newArray(size);
    double * absDelta = newArray(size);
    subtr(a,b,delta,size);
    vectorize(&fabs,delta,absDelta,size);
    double max = getMaxElement(absDelta,size);
    free(absDelta);
    free(delta);
    return max;
}

int main(void){
    double h = (x1-x0)/(N);
    double * gr = newArray(N+1);
    double * U_abar = newArray(N+1);
    double * U_2abar = newArray(N+1);
    double * v = newArray(N+1);
    linspace(x0,x1,gr,N+1);
    vectorize(&u0,gr,U_abar,N+1);

    //filename, size of grid, number of vectors to print, ... - pointers to vectors
    vectorsToCsv("test.csv",N+1,1,gr);
    printf("Error=%lg\n",getError(gr,gr,N+1));
//    solver(la,t,&U_abar[1],&v[1],N-1);

    free(gr);
    free(U_abar);
    free(U_2abar);
    free(v);
    return 0;
}