// -*- mode: c++; coding: utf-8 -*-
/*! @file main.cpp
    @brief Only defines tiny main()
*/
#include <lmpp/model.hpp>

//! Just instantiate and run Simulation
int main(int argc, char* argv[]) {
    std::vector<std::string> arguments(argv + 1, argv + argc);
    try {
        lmpp::Model instance(arguments);
        instance.run();
        instance.write();
    } catch (wtl::ExitSuccess) {}
    return 0;
}
