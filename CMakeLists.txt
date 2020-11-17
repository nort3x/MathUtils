cmake_minimum_required(VERSION 3.17)
project(MathUtils)

set(CMAKE_CXX_STANDARD 11)

add_executable(MathUtils main.cpp MathUtils.cpp MathUtils.h DataType.cpp DataType.h gnuplot.h Utils.h)