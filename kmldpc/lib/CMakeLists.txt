# include src
set(SUR_DIR lab)
include_directories(${CMAKE_CURRENT_LIST_DIR}/${SUR_DIR})

# using set to include the source code since the aux_source_directory
# cannot find *.hpp files.
set(LAB_SRC
        ${SUR_DIR}/binary5gldpccodec.hpp
        ${SUR_DIR}/binaryldpccodec.hpp
        ${SUR_DIR}/linearsystem.hpp
        ${SUR_DIR}/log.hpp
        ${SUR_DIR}/mat.hpp
        ${SUR_DIR}/modem.hpp
        ${SUR_DIR}/modemlinearsystem.hpp
        ${SUR_DIR}/randnum.hpp
        ${SUR_DIR}/sourcesink.hpp
        ${SUR_DIR}/utility.hpp
        ${SUR_DIR}/xorsegcodec.hpp)

#aux_source_directory(${SUR_DIR} LAB_SRC)

# build as library
add_library(lab ${LAB_SRC})

# set using c++
set_target_properties(lab PROPERTIES LINKER_LANGUAGE CXX)

target_include_directories(lab
        INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
        )