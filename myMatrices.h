//
// Created by Hell on 04.01.2018.
//

#ifndef NMPDE_C_MYMATRICES_H
#define NMPDE_C_MYMATRICES_H

#include "myNumPyV2.h"

TYPE** newMatrix(size_t hSize, size_t vSize);
void freeMatrix(TYPE **a);
void printMatrix(TYPE **a, size_t hSize, size_t vSize);

TYPE** newMatrixTransposed(TYPE **a, size_t hSize, size_t vSize);
void transposeMatrix(TYPE **a, TYPE **ta, size_t hSize, size_t vSize);

// dim(a) = [hSize, vSize], dim(b) = [hSize], dim(res) = [vSize]
void matrixVector(TYPE **a, TYPE *b, TYPE *res, size_t hSize, size_t vSize);

void matrixReduce(TYPE (*vFun)(TYPE *, size_t), TYPE **vectors, size_t size, TYPE *res, int numElements);

#endif //NMPDE_C_MYMATRICES_H
