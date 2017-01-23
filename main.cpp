// -*- mode: c++; coding: utf-8 -*-
/*! @file main.cpp
    @brief Only defines tiny main()
*/
#include <likeligrid/program.hpp>

//! Just instantiate and run Program
int main(int argc, char* argv[]) {
    std::vector<std::string> arguments(argv + 1, argv + argc);
    try {
        likeligrid::Program program(arguments);
        program.run();
    } catch (const wtl::ExitSuccess& e) {
        std::cerr << e.what() << std::endl;
    }
    return 0;
}
