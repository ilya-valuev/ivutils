project(ivutils)
include(sources.cmake)
set(IVUTILS_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include" CACHE STRING "IvUtils include dir")
message("IvUtils include dir:" ${IVUTILS_INCLUDE_DIR})

include_directories(${IVUTILS_INCLUDE_DIR})
add_library(ivutils ${BUILDCONF} ${IVUTILS_SOURCES} ${IVUTILS_HEADERS})
target_link_libraries(ivutils)
