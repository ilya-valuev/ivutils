project(ivutils)
include(sources.cmake)
set(IVUTILS_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include" CACHE STRING "IvUtils include dir")
message("IvUtils include dir:" ${IVUTILS_INCLUDE_DIR})

set(ERF_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/src/erf_namd.c" "${CMAKE_CURRENT_SOURCE_DIR}/src/erf.c" CACHE STRINGS "Source for erf implementation")
set(${IVUTILS_SOURCES} ${IVUTILS_SOURCES} ${ERF_SOURCE})


include_directories(${IVUTILS_INCLUDE_DIR})
add_library(ivutils ${BUILDCONF} ${IVUTILS_SOURCES} ${IVUTILS_HEADERS})
target_link_libraries(ivutils)
