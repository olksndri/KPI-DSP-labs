cmake_minimum_required(VERSION 3.10.0)

project(dsp_lab VERSION 0.1.0 LANGUAGES C CXX)

file(GLOB_RECURSE MATH_SRC "math/*.cpp")
file(GLOB_RECURSE DSP_SRC "dsp/*.cpp")
file(GLOB_RECURSE UTILS_SRC "utils/*.cpp")

add_executable(dsp_lab
    main.cpp
    ${MATH_SRC}
    ${DSP_SRC}
    ${UTILS_SRC}
)


target_include_directories(dsp_lab 
    PUBLIC 
    dsp
    libs
    math
    utils
) 

find_package(SndFile REQUIRED)
# find_package(SDL2 REQUIRED CONFIG REQUIRED COMPONENTS SDL2)
find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)

target_link_libraries(dsp_lab PRIVATE SndFile::sndfile)
# target_link_libraries(dsp_lab PRIVATE SDL2::SDL2)
target_link_libraries(dsp_lab PRIVATE Python3::Python Python3::NumPy)
