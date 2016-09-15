// -*- mode: c++; coding: utf-8 -*-
/*! @file main.cpp
    @brief Only defines tiny main()
*/
#include <lmpp/program.hpp>

//! Just instantiate and run Simulation
int main(int argc, char* argv[]) {
    std::vector<std::string> arguments(argv + 1, argv + argc);
    try {
        lmpp::Program program(arguments);
        program.run();
        program.write();
    } catch (wtl::ExitSuccess) {}
    return 0;
}
