cmake_minimum_required(VERSION 3.8)
project(NMPDE_C)

set(CMAKE_C_STANDARD 99)

set(SOURCE_FILES myMatrices.h myMatrices.c myNumPyV2.h myNumPyV2.c ex06.c)
link_libraries(-lm)
add_executable(NMPDE_C ${SOURCE_FILES})