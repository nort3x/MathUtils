cmake_minimum_required(VERSION 3.17)
project(MathUtils)

set(CMAKE_CXX_STANDARD 11)

add_executable(MathUtils main.cpp library.cpp library.h DataType.cpp DataType.h)