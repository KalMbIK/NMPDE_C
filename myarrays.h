//
// Created by pusheen on 30.10.17.
//

#ifndef NMPDE_C_ARRAYS_H
#define NMPDE_C_ARRAYS_H

#ifndef TYPE
#define TYPE double
#endif //TYPE

#define g(i) ar[(i)]

struct array {
    int size;
    TYPE* ar;
};
typedef struct array Array;

Array* newArray(int N);
void freeArray(Array* a);

//TYPE get(Array* ar, int idx);

void set(Array *ar, int idx, TYPE value);

void printArray(Array *ar);

#endif //NMPDE_C_ARRAYS_H
