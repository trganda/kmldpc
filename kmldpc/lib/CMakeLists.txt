# include src
set(LAB_DIR lab)
set(TOML11_DIR toml11)
include_directories(${CMAKE_CURRENT_LIST_DIR}/${LAB_DIR})

# using set to include the source code since the aux_source_directory
# cannot find *.hpp files.
set(LAB_SRC
        ${LAB_DIR}/binary5gldpccodec.hpp
        ${LAB_DIR}/binaryldpccodec.hpp
        ${LAB_DIR}/linearsystem.hpp
        ${LAB_DIR}/log.hpp
        ${LAB_DIR}/mat.hpp
        ${LAB_DIR}/modem.hpp
        ${LAB_DIR}/modemlinearsystem.hpp
        ${LAB_DIR}/randnum.hpp
        ${LAB_DIR}/sourcesink.hpp
        ${LAB_DIR}/utility.hpp
        ${LAB_DIR}/xorsegcodec.hpp)

# add toml11 https://github.com/ToruNiina/toml11/releases/tag/v3.6.0
set(TOML11_SRC
        ${TOML11_DIR}/color.hpp
        ${TOML11_DIR}/combinator.hpp
        ${TOML11_DIR}/comments.hpp
        ${TOML11_DIR}/datetime.hpp
        ${TOML11_DIR}/exception.hpp
        ${TOML11_DIR}/from.hpp
        ${TOML11_DIR}/get.hpp
        ${TOML11_DIR}/into.hpp
        ${TOML11_DIR}/lexer.hpp
        ${TOML11_DIR}/literal.hpp
        ${TOML11_DIR}/parser.hpp
        ${TOML11_DIR}/region.hpp
        ${TOML11_DIR}/result.hpp
        ${TOML11_DIR}/serializer.hpp
        ${TOML11_DIR}/source_location.hpp
        ${TOML11_DIR}/storage.hpp
        ${TOML11_DIR}/string.hpp
        ${TOML11_DIR}/toml.hpp
        ${TOML11_DIR}/traits.hpp
        ${TOML11_DIR}/types.hpp
        ${TOML11_DIR}/utility.hpp
        ${TOML11_DIR}/value.hpp
        )

# build as library
add_library(lab ${LAB_SRC})
add_library(toml11 ${TOML11_SRC})

# set using c++
set_target_properties(lab PROPERTIES LINKER_LANGUAGE CXX)
set_target_properties(toml11 PROPERTIES LINKER_LANGUAGE CXX)