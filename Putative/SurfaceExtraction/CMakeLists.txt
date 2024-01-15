cmake_minimum_required(VERSION 2.8)
 
project(Surface)
 
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})
 
add_executable(Surface main.cpp)
target_link_libraries(Surface ${ITK_LIBRARIES})
