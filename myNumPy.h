//
// Created by pusheen on 30.10.17.
//

#ifndef NMPDE_C_MYGRID_H
#define NMPDE_C_MYGRID_H

#include "myArrays.h"

double dot(TYPE* a, TYPE* b, int size);
double getNorm(TYPE* a, int size);
TYPE* add(TYPE *a, TYPE *b, int size);
TYPE* subtr(TYPE *a, TYPE *b, int size);

TYPE* linspace(TYPE a, TYPE b, int size);
TYPE* vectorize(TYPE (*f)(TYPE), TYPE *gr, int size);
TYPE* vectorizeOfArrayInPoint(TYPE (*f)(TYPE, TYPE), TYPE *gr, int size, TYPE t);

#endif //NMPDE_C_MYGRID_H
