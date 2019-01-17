project(ivutils)
include(sources.cmake)
get_filename_component(IVUTILS_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include ABSOLUTE)
message("IvUtils include dir:" ${IVUTILS_INCLUDE_DIR})

if(NOT ${HAVE_ERF})   
  set(ERF_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/src/erf_namd.c" CACHE FILEPATH "Source for erf implementation")
  set(IVUTILS_SOURCES ${IVUTILS_SOURCES} ${ERF_SOURCE})
endif()

include_directories(${IVUTILS_INCLUDE_DIR})
add_library(ivutils ${BUILDCONF} ${IVUTILS_SOURCES} ${IVUTILS_HEADERS})
target_link_libraries(ivutils)
