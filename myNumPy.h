//
// Created by pusheen on 30.10.17.
//

#ifndef NMPDE_C_MYGRID_H
#define NMPDE_C_MYGRID_H

#ifndef TYPE
#define TYPE double
#endif //TYPE

TYPE* newArray(int size);
void printArray(TYPE *a, int size);

double dot(TYPE* a, TYPE* b, int size);
double getNorm(TYPE* a, int size);
double getMaxElement(TYPE* a, int size);

TYPE* add(TYPE *a, TYPE *b, int size);
TYPE* subtr(TYPE *a, TYPE *b, int size);

TYPE* linspace(TYPE a, TYPE b, int size);
TYPE* vectorize(TYPE (*f)(TYPE), TYPE *gr, int size);
TYPE* vectorizeOfArrayInPoint(TYPE (*f)(TYPE, TYPE), TYPE *gr, int size, TYPE t);

#endif //NMPDE_C_MYGRID_H
