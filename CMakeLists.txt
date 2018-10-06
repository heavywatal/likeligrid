cmake_minimum_required(VERSION 3.1)
project(likeligrid
  VERSION 0.1.0
  LANGUAGES CXX)

include(GNUInstallDirs)
set(CMAKE_VERBOSE_MAKEFILE ON)

set(Boost_NO_BOOST_CMAKE ON)
find_package(Boost REQUIRED COMPONENTS program_options context filesystem system)
message(STATUS "Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}")

find_package(ZLIB REQUIRED)
find_package(Threads REQUIRED)
find_package(sfmt 0.7.0)
find_package(wtl 0.7.2)

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
set(CMAKE_CXX_FLAGS_DEV "-O2 -g")
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
add_compile_options(-march=native -Wall -Wextra -Wpedantic)

message(STATUS "CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")
if(${CMAKE_INSTALL_PREFIX} MATCHES linuxbrew)
  set(CMAKE_INSTALL_LIBDIR lib)
endif()

add_subdirectory(src)

add_executable(${PROJECT_NAME}-exe main.cpp $<TARGET_OBJECTS:obj>)
target_link_libraries(${PROJECT_NAME}-exe PRIVATE obj)
set_target_properties(${PROJECT_NAME}-exe PROPERTIES
  OUTPUT_NAME ${PROJECT_NAME}
)
install(TARGETS ${PROJECT_NAME}-exe
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

add_library(${PROJECT_NAME} SHARED $<TARGET_OBJECTS:obj>)
target_include_directories(${PROJECT_NAME} BEFORE
  INTERFACE $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
)
target_link_libraries(${PROJECT_NAME} PRIVATE obj)

include(CTest)
if(BUILD_TESTING)
  add_subdirectory(test)
endif()
