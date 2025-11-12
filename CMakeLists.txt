cmake_minimum_required(VERSION 3.30)
project(likeligrid
  VERSION 0.1.0
  LANGUAGES CXX)

set(CMAKE_VERBOSE_MAKEFILE ON)
include(CMakePrintHelpers)
include(FetchContent)
include(GNUInstallDirs)
cmake_print_variables(CMAKE_INSTALL_PREFIX)

add_compile_options(
  -Wall -Wextra -Wpedantic -Wconversion -Wsign-conversion -Wshadow
  $<IF:$<BOOL:$<CONFIG>>,,-Werror>
  $<IF:$<BOOL:$<CONFIG>>,,-O2>
  $<IF:$<BOOL:$<CONFIG>>,,-g>
  $<$<STREQUAL:${CMAKE_SYSTEM_PROCESSOR},x86_64>:-march=native>
  $<$<STREQUAL:${CMAKE_SYSTEM_PROCESSOR},arm64>:-march=armv8.4-a>
)

macro(Pkg_MakeAvailable)
  cmake_parse_arguments(arg "" "" "" ${ARGN})
  FetchContent_MakeAvailable(${arg_UNPARSED_ARGUMENTS})
  foreach(pkg IN LISTS arg_UNPARSED_ARGUMENTS)
    cmake_print_variables(${pkg}_FOUND ${pkg}_DIR)
    if(NOT ${pkg}_FOUND)
      cmake_print_variables(${pkg}_SOURCE_DIR ${pkg}_BINARY_DIR)
    endif()
  endforeach()
endmacro()

FetchContent_Declare(clippson
  GIT_REPOSITORY https://github.com/heavywatal/clippson.git
  GIT_TAG fb915f6f3b00e4854ab159aea38f620c4934da5c
  EXCLUDE_FROM_ALL
  FIND_PACKAGE_ARGS 0.8.11
)

FetchContent_Declare(wtl
  GIT_REPOSITORY https://github.com/heavywatal/cxxwtl.git
  GIT_TAG fdaa1420abdb25a54ae7262aeaddf622c21a5935
  EXCLUDE_FROM_ALL
  FIND_PACKAGE_ARGS 0.11.0
)

Pkg_MakeAvailable(clippson wtl)

find_package(Boost CONFIG REQUIRED COMPONENTS context)
cmake_print_variables(Boost_INCLUDE_DIRS)

add_library(${PROJECT_NAME} STATIC)
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_17)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_EXTENSIONS OFF)
target_include_directories(${PROJECT_NAME} PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
)
target_link_libraries(${PROJECT_NAME} PUBLIC
  wtl::wtl wtl::zlib wtl::threads
  PRIVATE clippson::clippson Boost::context
)
add_subdirectory(src)

add_executable(${PROJECT_NAME}-exe src/main.cpp)
target_link_libraries(${PROJECT_NAME}-exe PRIVATE ${PROJECT_NAME})
set_target_properties(${PROJECT_NAME}-exe PROPERTIES
  OUTPUT_NAME ${PROJECT_NAME}
)
install(TARGETS ${PROJECT_NAME}-exe
  EXPORT ${PROJECT_NAME}-config
)

set(config_destination ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
install(EXPORT ${PROJECT_NAME}-config
  DESTINATION ${config_destination}
  NAMESPACE ${PROJECT_NAME}::
)

set(version_file ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake)
include(CMakePackageConfigHelpers)
write_basic_package_version_file(${version_file}
  COMPATIBILITY AnyNewerVersion
)
install(FILES ${version_file}
  DESTINATION ${config_destination}
)

include(CTest)
if(BUILD_TESTING AND PROJECT_IS_TOP_LEVEL)
  add_subdirectory(test)
endif()
