cmake_minimum_required(VERSION 3.10.0)

project(dsp_lab VERSION 0.1.0 LANGUAGES C CXX)

add_executable(dsp_lab
    main.cpp
    src/t.cpp
)

target_include_directories(dsp_lab 
    PUBLIC 
    include
) 

find_package(SndFile REQUIRED)
find_package(SDL2 REQUIRED CONFIG REQUIRED COMPONENTS SDL2)
# find_package(Gnuplot REQUIRED)

target_link_libraries(dsp_lab PRIVATE SndFile::sndfile)
target_link_libraries(dsp_lab PRIVATE SDL2::SDL2)
