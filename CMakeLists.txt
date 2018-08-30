cmake_minimum_required(VERSION 3.0.0)
project(photo VERSION 0.1.0)

include(CTest)
enable_testing()

set(CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pthread -I/usr/local/include/eigen3")

set (source_dir "${PROJECT_SOURCE_DIR}/sources/")
file (GLOB source_files "${source_dir}/*.cpp")

add_executable(photo ${source_files})

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})

include(CPack)
