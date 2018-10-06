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
add_library(obj OBJECT ${srcfiles})
target_include_directories(obj BEFORE
  PUBLIC ${Boost_INCLUDE_DIRS}
)
set_property(TARGET obj PROPERTY POSITION_INDEPENDENT_CODE ON)
target_link_libraries(obj PUBLIC wtl::wtl wtl::sfmt ${Boost_LIBRARIES} Threads::Threads ZLIB::ZLIB)
