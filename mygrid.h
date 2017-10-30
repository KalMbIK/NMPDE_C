//
// Created by pusheen on 30.10.17.
//

#ifndef NMPDE_C_MYGRID_H
#define NMPDE_C_MYGRID_H
#include "myarrays.h"
Array * gethGrid(double a, double b, int N);
double* vectorizeDoubleFun(double (*f)(double, double), double* gr, int N, double t);
#endif //NMPDE_C_MYGRID_H
