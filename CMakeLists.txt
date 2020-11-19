cmake_minimum_required(VERSION 3.17)
project(MathUtils)

set(CMAKE_CXX_STANDARD 11)

add_executable(MathUtils main.cpp MathUtils.cpp MathUtils.h DataType.cpp DataType.h gnuplot.h Utils.h)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
target_link_libraries(MathUtils PRIVATE Threads::Threads)