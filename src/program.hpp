/*! @file program.hpp
    @brief Interface of Program class
*/
#pragma once
#ifndef LIKELIGRID_PROGRAM_HPP_
#define LIKELIGRID_PROGRAM_HPP_

#include <string>
#include <vector>

namespace boost {
    namespace program_options {
        class options_description;
    }
}

namespace likeligrid {

/*! @brief Represents single run
*/
class Program {
  public:
    //! Parse command arguments
    Program(const std::vector<std::string>& args);

    //! Top level function that should be called once from main()
    void run();

    /////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
  private:
    boost::program_options::options_description options_desc();
    boost::program_options::options_description positional_desc();
    [[noreturn]] void help_and_exit();
    void test(const int flag);
    std::string make_outdir() const;

    unsigned int CONCURRENCY = 1u;
    size_t MAX_SITES = 3u;
    std::string INFILE = "-";
    bool GRADIENT_MODE = false;
    std::vector<size_t> EPISTASIS_PAIR = {0u, 0u};
    bool PLEIOTROPY = false;
};

} // namespace likeligrid

#endif /* LIKELIGRID_PROGRAM_HPP_ */
