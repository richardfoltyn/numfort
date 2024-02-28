# Sample CMake file illustrating how to reference the NUMFORT library

cmake_minimum_required(VERSION 3.12)

project(numfort_client Fortran)

find_package(numfort REQUIRED)

add_executable(client client.f90)
target_link_libraries(client numfort::numfort)
