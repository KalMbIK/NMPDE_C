//
// Created by pusheen on 30.10.17.
//

#include "myTDMAlgorithms.h"

void TDMA_solver(const double *a, const double *b, const double *c, double *cPrime, const double *d, double *dPrime,
                 double *result, int size) {
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

void TDM_dot(const double *a, const double *b, const double *c, const double *u, double *d, int size){
    d[0] = b[0]*u[0] + c[0] * u[1];
    d[size-1] = a[size-1]*u[size-2]+b[size-1]*u[size-1];
    for (int i = 1; i < size-2; i++)
        d[i] = a[i]*u[i-1] + b[i]*u[i] + c[i]*u[i+1];
}
