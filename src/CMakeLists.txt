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

file(GLOB srcfiles *.cpp)
add_library(objlib OBJECT ${srcfiles})
target_include_directories(objlib INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
)
target_link_libraries(objlib PUBLIC wtl::wtl wtl::sfmt Threads::Threads ZLIB::ZLIB)
target_link_libraries(objlib PUBLIC Boost::program_options Boost::context Boost::filesystem)
