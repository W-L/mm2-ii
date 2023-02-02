cmake_minimum_required(VERSION 3.24)
set(CMAKE_C_STANDARD 17)

project(mm2ii C)


add_executable(mm2ii
        main.c
        index_mm2ii.c
        ../minimap2/index.c
        ../minimap2/sketch.c
        ../minimap2/kalloc.c
        ../minimap2/misc.c
        ../minimap2/kthread.c
        ../minimap2/bseq.c
        ../minimap2/options.c)


#find_library(HASHMAP_LIB NAMES libhashmap.a hashmap.h PATHS "../hashmap/include")
#find_path(HASHMAP_DIR hashmap.h "../hashmap/include")

#include_directories(${HASHMAP_DIR})
#target_link_libraries(mm2ii pthread z ${HASHMAP_LIB})
target_link_libraries(mm2ii pthread z)

target_include_directories(mm2ii PUBLIC ../minimap2)
target_include_directories(mm2ii PUBLIC ../)


