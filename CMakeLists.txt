cmake_minimum_required(VERSION 3.0.0)
project(photo VERSION 0.1.0)

include(CTest)
enable_testing()

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

set(CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pipe -Wall -pthread -I/usr/local/include/eigen3")

set (source_dir "${PROJECT_SOURCE_DIR}/sources/")
file (GLOB source_files "${source_dir}/*.cpp")

add_executable(photo ${source_files})

target_link_libraries (photo Eigen3::Eigen)
target_link_libraries (photo lapacke)
target_link_libraries (photo lapack)
target_link_libraries (photo blas)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})

include(CPack)
