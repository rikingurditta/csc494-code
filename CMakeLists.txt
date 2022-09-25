cmake_minimum_required(VERSION 3.16)
project(cages)

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Libigl
include(libigl)

# Enable the target igl::glfw
igl_include(glfw)

# Add your project files

include_directories(${PROJECT_SOURCE_DIR}/TinyAD/include/TinyAD)
include_directories(${PROJECT_SOURCE_DIR}/include)

file(GLOB TINYAD_FILES ${PROJECT_SOURCE_DIR}/TinyAD/include/TinyAD/*.hh)
file(GLOB H_FILES ${PROJECT_SOURCE_DIR}/include/*.h)
file(GLOB SRC_FILES ${PROJECT_SOURCE_DIR}/src/*.cpp)

add_executable(${PROJECT_NAME} main.cpp ${TINYAD_FILES} ${H_FILES} ${SRC_FILES})
target_link_libraries(${PROJECT_NAME} PUBLIC igl::glfw)
