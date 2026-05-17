cmake_minimum_required(VERSION 3.10)
project(dsp_lab VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(raylib REQUIRED)
find_package(SDL2 REQUIRED)
find_package(PkgConfig REQUIRED)
find_package(TIFF REQUIRED)
pkg_check_modules(SNDFILE REQUIRED IMPORTED_TARGET sndfile)

set(Python3_EXECUTABLE "${CMAKE_SOURCE_DIR}/.venv/bin/python")
set(CMAKE_PREFIX_PATH "${CMAKE_SOURCE_DIR}/.venv")
find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)

file(GLOB_RECURSE MATH_SRC "math/*.cpp")
file(GLOB_RECURSE DSP_SRC "dsp/*.cpp")
file(GLOB_RECURSE UTILS_SRC "utils/*.cpp")
file(GLOB_RECURSE IMAGE_COMPRESSION_SRC "image_compression/*.cpp")

add_executable(dsp_lab
				main.cpp
				${MATH_SRC}
				${DSP_SRC}
				${UTILS_SRC}
				${IMAGE_COMPRESSION_SRC}
)

target_include_directories(dsp_lab PUBLIC
				${CMAKE_SOURCE_DIR}
				dsp
				libs
				math
				utils
				${Python3_NumPy_INCLUDE_DIRS}
)

target_link_libraries(dsp_lab PRIVATE
				raylib
				TIFF::TIFF
				SDL2::SDL2
				SDL2::SDL2main
				PkgConfig::SNDFILE
				Python3::Python
				Python3::NumPy
				m
				pthread
				dl
				rt
				X11
)
