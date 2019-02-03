# Sample build file to test linking against numfort

cmake_minimum_required(VERSION 3.0)
project(numfort_client Fortran)

find_package(numfort REQUIRED)

add_executable(hellonumfort hellonumfort.f90)
target_link_libraries(hellonumfort numfort::numfort)
