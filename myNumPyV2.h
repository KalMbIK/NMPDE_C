//
// Created by pusheen on 30.10.17.
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
void printArray(TYPE *a, int size);

//void arrayOfVectorsToCsv(char *filename, TYPE **vectors, int size, unsigned num);
//void vectorsToCsv(char *filename, int size, unsigned num, ...);

void vectorToFile(FILE *filePointer, TYPE *vector, int size);
void vectorsToCsv(char *filename, int size, unsigned num, ...);

double dot(TYPE* a, TYPE* b, int size);
double getNorm(TYPE* a, int size);
TYPE getMaxElement(TYPE *a, int size);

void add(TYPE *a, TYPE *b, TYPE* res, int size);
void subtr(TYPE *a, TYPE *b, TYPE* res, int size);

void linspace(TYPE a, TYPE b, TYPE* gr, int size);
void vectorize(TYPE (*f)(TYPE), TYPE *gr, TYPE* res, int size);
void vectorizeOfArrayInPoint(TYPE (*f)(TYPE, TYPE), TYPE *gr, TYPE* res, int size, TYPE t);

#endif //NMPDE_C_MYGRID_H
