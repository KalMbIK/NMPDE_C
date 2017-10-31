cmake_minimum_required(VERSION 3.8)
project(NMPDE_C)

set(CMAKE_C_STANDARD 99)

set(SOURCE_FILES ex02.c mygrid.h mygrid.c myArrays.h myArrays.c)
link_libraries(-lm)
add_executable(NMPDE_C ${SOURCE_FILES})