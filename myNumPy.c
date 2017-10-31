//
// Created by pusheen on 30.10.17.
//

#include "myNumPy.h"
#include <math.h>
#include <stdlib.h>

double dot(TYPE* a, TYPE* b, int size){
    double sum = 0;
    for (int i = 0; i < size; i++)
        sum += a[i]*b[i];
    return sum;
}
double getNorm(TYPE* a, int size){
    return sqrt(dot(a,a,size));
}
TYPE* add(TYPE *a, TYPE *b, int size){
    TYPE* res = newArray(size);
    for (int i = 0; i < size; i++)
        res[i] = a[i]+b[i];
    return res;
}
TYPE* subtr(TYPE *a, TYPE *b, int size){
    TYPE* res = newArray(size);
    for (int i = 0; i < size; i++)
        res[i] = a[i]-b[i];
    return res;
}

TYPE* linspace(TYPE a, TYPE b, int size) {
    TYPE step = (b-a)/(size-1);
    TYPE* gr = newArray(size);
    gr[0] = a;
    for (int i = 1; i < size; i++){
        gr[i]=gr[i-1]+step;
    }
    return gr;
}

TYPE* vectorize(TYPE (*f)(TYPE), TYPE *gr, int size){
    TYPE * result = newArray(size);
    for (int i = 0; i < size; i++){
        result[i]=f(gr[i]);
    }
    return result;
}

TYPE* vectorizeOfArrayInPoint(TYPE (*f)(TYPE, TYPE), TYPE *gr, int size, TYPE t){
    TYPE * result = newArray(size);
    for (int i = 0; i < size; i++){
        result[i]=f(gr[i],t);
    }
    return result;
}
