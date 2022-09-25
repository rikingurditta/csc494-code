cmake_minimum_required(VERSION 3.16)
project(cages)

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Libigl
include(libigl)

# Enable the target igl::glfw
igl_include(glfw)

# Add your project files
file(GLOB H_FILES ${PROJECT_SOURCE_DIR}/include/*.h)
file(GLOB TINYAD ${PROJECT_SOURCE_DIR}/TinyAD/include/TinyAD/*.hh)
file(GLOB SRC_FILES ${PROJECT_SOURCE_DIR}/src/*.cpp)
add_executable(${PROJECT_NAME} main.cpp ${TINYAD} ${H_FILES} ${SRC_FILES})
target_link_libraries(${PROJECT_NAME} PUBLIC igl::glfw)
