//
// Created by pusheen on 31.10.17.
//
#include <stdlib.h>
#include <stdio.h>
#include "myArrays.h"


TYPE* newArray(int size){
    return (TYPE*)calloc(size,sizeof(TYPE));
}
void printArray(TYPE *a, int size){
    for (int i = 0; i < size; i++){
        printf("%lg ", a[i]);
    }
    printf("\n");
}
