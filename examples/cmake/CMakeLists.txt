# Sample CMake file illustrating how to reference the NUMFORT library

cmake_minimum_required(VERSION 3.1)

project(numfort_consumer Fortran)

find_package(numfort REQUIRED)

set(SOURCE_FILES numfort_consumer.f90)

add_executable(consumer ${SOURCE_FILES})
target_link_libraries(consumer numfort::numfort)
