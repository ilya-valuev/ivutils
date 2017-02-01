file(GLOB_RECURSE IVUTILS_HEADERS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "include/*.h")

set(IVUTILS_SOURCES "src/graph/basis_3.cpp" "src/graph/cluster.cpp" "src/graph/vecor_3.cpp" "src/graph/vector_n.cpp"  
                    "src/common.cpp" "src/four/four.c"
                    "src/rparam/rparam.c" "src/logexc.cpp" "src/mcarlo.cpp" "src/record.cpp" "src/seqpack.cpp"
                    "src/statist.cpp" "src/timers.cpp")
