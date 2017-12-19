//
// Updated by pusheen on 17.12.17.
//

#ifndef NMPDE_C_MYGRID_H
#define NMPDE_C_MYGRID_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <stdarg.h>

#ifndef TYPE
#define TYPE double
#endif //TYPE

TYPE* newArray(size_t size);
void printArray(const TYPE *a, int size);

//void arrayOfVectorsToCsv(char *filename, TYPE **vectors, int size, unsigned num);
//void vectorsToCsv(char *filename, int size, unsigned num, ...);

void vectorToFile(FILE *filePointer, TYPE *vector, int size);
void vectorsToCsv(char *filename, int size, unsigned num, ...);

double dot(const TYPE *a, const TYPE *b, int size);
double getNorm(const TYPE *a, int size);
//TYPE reduce(TYPE (*f)(TYPE), TYPE *gr, int size);
TYPE getMaxElement(const TYPE *a, int size);
TYPE getMinElement(const TYPE *a, int size);

void add(const TYPE *a, const TYPE *b, TYPE *res, int size);
void subtr(const TYPE *a, const TYPE *b, TYPE *res, int size);

void linspace(TYPE a, TYPE b, TYPE* gr, int size);
void vectorize(TYPE (*f)(TYPE), const TYPE *gr, TYPE *res, int size);
void vectorizeOfArrayInPoint(TYPE (*f)(TYPE, TYPE), const TYPE *gr, TYPE *res, int size, TYPE t);

#endif //NMPDE_C_MYGRID_H
