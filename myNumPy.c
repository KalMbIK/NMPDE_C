//
// Created by pusheen on 30.10.17.
//

#include "myNumPy.h"
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <stdio.h>

TYPE* newArray(int size){
    return (TYPE*)calloc(size,sizeof(TYPE));
}
void printArray(TYPE *a, int size){
    for (int i = 0; i < size; i++){
        printf("%lg ", a[i]);
    }
    printf("\n");
}


double dot(const TYPE *a, const TYPE *b, int size){
    double sum = 0;
    for (int i = 0; i < size; i++)
        sum += a[i]*b[i];
    return sum;
}
double getNorm(TYPE* a, int size){
    return sqrt(dot(a,a,size));
}

TYPE getMaxElement(const TYPE *a, int size){
    TYPE max = INT_MIN;
    for (int i = 0; i < size; i++){
        if (a[i] > max)
            max = a[i];
    }
    return max;
}

TYPE getError(TYPE* a, TYPE* b, int size){
    TYPE * delta = subtr(a,b,size);
    TYPE * absDelta = vectorize(&fabs,delta,size);
    TYPE max = getMaxElement(absDelta,size);
    free(absDelta);
    free(delta);
    return max;
}

TYPE* add(const TYPE *a, const TYPE *b, int size){
    TYPE* res = newArray(size);
    for (int i = 0; i < size; i++)
        res[i] = a[i]+b[i];
    return res;
}
TYPE* subtr(const TYPE *a, const TYPE *b, int size){
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

TYPE* vectorize(TYPE (*f)(TYPE), const TYPE *gr, int size){
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
