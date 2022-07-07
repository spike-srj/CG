cmake_minimum_required(VERSION 3.0)
project(mywork)

include_directories(include "src/")
set(CMAKE_BUILD_TYPE Debug)
#FILE(GLOB_RECURSE app_sources ${CMAKE_SOURCE_DIR}/src/*.*)
add_executable(mywork src/main.cpp src/material.cpp src/vec3.cpp src/hittable.cpp src/sphere.cpp src/hittable_list.cpp)