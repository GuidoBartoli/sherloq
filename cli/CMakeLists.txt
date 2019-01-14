cmake_minimum_required (VERSION 3.5)
project (Sherlock)

find_package (OpenCV REQUIRED)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -O3")

add_executable (Sherlock src/sherlock.cpp)
target_link_libraries (Sherlock ${OpenCV_LIBS} boost_filesystem boost_system)
