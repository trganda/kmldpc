# CMakeList.txt : CMake project for kmldpc, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

# include src
include_directories(src, include)

aux_source_directory(src SRC)

# Add source to this project's executable.
add_executable (kmldpc "kmldpc.cpp" ${SRC})
# TODO: Add tests and install targets if needed.
