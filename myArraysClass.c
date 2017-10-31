//
// Created by pusheen on 30.10.17.
//
#include "myArraysClass.h"
#include <stdio.h>
#include "stdlib.h"

void printArray(Array *ar) {
    for (int i = 0; i < ar->size; i++){
        printf("%lg ", ar->ar[i]);
    }
    printf("\n");
}

Array* newArray(int N){
    Array* newOne = (Array*)malloc(sizeof(int)+N*sizeof(TYPE)+sizeof(TYPE*));
    newOne->size = N;
    newOne->ar = (TYPE*)(&newOne->size+sizeof(int));
    for (int i = 0; i < N; i++)
        newOne->ar[i] = 0;
    return newOne;
}
void freeArray(Array* a){
    free((void *)a);
}

//TYPE get(Array* ar, int idx) {
//    if (ar != NULL)
//        if ((idx >= 0) && (idx < ar->size))
//            return ar->ar[idx];
//    return -100;
//}

void set(Array *ar, int idx, TYPE value) {
    if (ar != NULL)
        if ((idx >= 0) && (idx < ar->size))
            ar->ar[idx] = value;
}