cmake_minimum_required(VERSION 3.0)

project(c_functions)


add_library(cfunc SHARED
            src/interface.c
            src/pivi.c
            src/pivi.h
            src/erfinv.h src/erfinv.c
            src/boost_wrapper.cpp
            src/boost_wrapper.hpp
            )

add_executable(test src/test.c src/pivi.h src/pivi.c src/erfinv.h src/erfinv.c src/boost_wrapper.cpp src/boost_wrapper.hpp)

find_package(GSL REQUIRED)
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
target_link_libraries(cfunc GSL::gsl GSL::gslcblas m Threads::Threads dpivi)
target_link_libraries(test GSL::gsl GSL::gslcblas m Threads::Threads dpivi)

