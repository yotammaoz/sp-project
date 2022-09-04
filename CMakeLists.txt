cmake_minimum_required(VERSION 3.21)
project(final C)

set(CMAKE_C_STANDARD 90)

find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

add_executable(final matrix_calculations.c spkmeansmodule.c spkmeans.h useful.h useful.c kmeans.h kmeans.c jacobi.c jacobi.h matrix_calculations.h)
