# CMakeList.txt : CMake project for kmldpc, include source and define

# include src
include_directories(src, include)

aux_source_directory(src SRC)

link_directories(/usr/local/lib)
# Add source to this project's executable.
add_executable (kmldpc "kmldpc.cpp" ${SRC})
target_link_libraries(kmldpc matio)
# TODO: Add tests and install targets if needed.
