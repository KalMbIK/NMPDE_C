//
// Created by pusheen on 30.10.17.
//

#ifndef NMPDE_C_MYTDMALGORITHMS_H
#define NMPDE_C_MYTDMALGORITHMS_H

void TDMA_solver(const double *a, const double *b, const double *c, double *cPrime, const double *d, double *dPrime,
                 double *result, int size);

void TDM_dot(const double *a, const double *b, const double *c, const double *u, double *d, int size);


#endif //NMPDE_C_MYTDMALGORITHMS_H
