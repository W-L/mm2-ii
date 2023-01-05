cmake_minimum_required(VERSION 3.24)
set(CMAKE_C_STANDARD 17)

project(mm2ii C)


add_executable(mm2ii
        src/main.c
        src/index_mm2ii.c
        minimap2/index.c
        minimap2/sketch.c
        minimap2/kalloc.c
        minimap2/misc.c
        minimap2/kthread.c
        minimap2/bseq.c
        minimap2/options.c
        )

target_link_libraries(mm2ii z hashmap)

target_include_directories(mm2ii PUBLIC minimap2)



