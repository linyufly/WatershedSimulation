cmake_minimum_required(VERSION 2.8)
project(VincentSoille)

find_package(VTK)
if(VTK_FOUND)
	include(${VTK_USE_FILE})
	message("VTK found")
endif()

add_executable(watershed main.cpp)

target_link_libraries(watershed ${VTK_LIBRARIES})
