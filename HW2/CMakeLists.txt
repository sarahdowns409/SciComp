cmake_minimum_required(VERSION 3.23)
project(MATH233-hw2)

set(CMAKE_CXX_STANDARD 14)

add_executable(hw2 main.cpp Grid2d.cpp Grid2d.h SL_method.cpp SL_method.h Godunov.h Godunov.cpp)
