//
// Created by pusheen on 31.10.17.
//

#ifndef NMPDE_C_MYARRAYS_H
#define NMPDE_C_MYARRAYS_H

#ifndef TYPE
#define TYPE int
#endif //TYPE

TYPE* newArray(int size);
void freeArray(TYPE* a);
void printArray(TYPE *a, int size);


#endif //NMPDE_C_MYARRAYS_H
