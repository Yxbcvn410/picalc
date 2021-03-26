cmake_minimum_required(VERSION 3.17)
project(picalc)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "-O3 -Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -Wextra")

add_executable(picalc main.cpp)

find_package(MPI)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})
target_link_libraries(picalc ${MPI_C_LIBRARIES})