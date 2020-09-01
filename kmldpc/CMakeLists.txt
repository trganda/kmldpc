# CMakeList.txt : CMake project for kmldpc, include source and define

# Will add __FILENAME__ macros for all source files, which is the filename without full find_path
# Courtesy of SO
function(define_file_basename_for_sources targetname)
    get_target_property(source_files "${targetname}" SOURCES)
    foreach(sourcefile ${source_files})
        # Get source file's current list of compile definitions.
        get_property(defs SOURCE "${sourcefile}"
                PROPERTY COMPILE_DEFINITIONS)
        # Add the FILE_BASENAME=filename compile definition to the list.
        get_filename_component(basename "${sourcefile}" NAME)
        list(APPEND defs "__FILENAME__=\"${basename}\"")
        # Set the updated compile definitions on the source file.
        set_property(
                SOURCE "${sourcefile}"
                PROPERTY COMPILE_DEFINITIONS ${defs})
    endforeach()
endfunction()

# include src
include_directories(src, include)
aux_source_directory(src SRC)

#link_directories(/usr/local/lib)
#target_link_libraries(kmldpc matio)

find_package(matio REQUIRED)
if(MATIO_FOUND)
    include_directories(${MATIO_INCLUDE_DIR})
else()
    message("\nmatio directory not found.")
endif()

# Add source to this project's executable.
add_executable (kmldpc "kmldpc.cpp" ${SRC})
target_link_libraries(kmldpc ${MATIO_LIBRARY})

define_file_basename_for_sources(kmldpc)
# TODO: Add tests and install targets if needed.
