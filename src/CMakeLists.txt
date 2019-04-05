execute_process(
  COMMAND git describe --always --dirty
  OUTPUT_VARIABLE GIT_DESCRIBE
  ERROR_QUIET
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
if(GIT_DESCRIBE)
  set(PROJECT_VERSION ${GIT_DESCRIBE})
endif()

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/version.hpp.in
  ${CMAKE_CURRENT_SOURCE_DIR}/version.hpp @ONLY
)

# Be patient until 3.13 is popularized
add_library(objlib OBJECT
  ${CMAKE_CURRENT_SOURCE_DIR}/genotype.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/gradient_descent.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/gridsearch.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/pathtype.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/program.cpp
)
target_include_directories(objlib INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
)
target_link_libraries(objlib PUBLIC
  wtl::wtl
  wtl::sfmt
  clippson::clippson
  Threads::Threads
  ZLIB::ZLIB
  Boost::context
)
