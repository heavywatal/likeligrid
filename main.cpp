/*! @file main.cpp
    @brief Only defines tiny main()
*/
#include "src/program.hpp"

//! Just instantiate and run Program
int main(int argc, char* argv[]) {
    std::vector<std::string> arguments(argv, argv + argc);
    likeligrid::Program program(arguments);
    program.run();
    return 0;
}
