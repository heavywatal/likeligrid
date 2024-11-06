execute_process(
  COMMAND git describe --always --dirty
  OUTPUT_VARIABLE GIT_DESCRIBE
  ERROR_QUIET
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
if(GIT_DESCRIBE)
  set(PROJECT_VERSION ${GIT_DESCRIBE})
endif()
cmake_print_variables(PROJECT_VERSION)

configure_file(
  version.hpp.in
  version.hpp @ONLY
)

target_sources(${PROJECT_NAME} PRIVATE
  genotype.cpp
  gradient_descent.cpp
  gridsearch.cpp
  pathtype.cpp
  program.cpp
)
