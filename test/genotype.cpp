#include "genotype.hpp"

#include <iostream>
#include <sstream>

int main() {
    std::stringstream sst;
    sst <<
R"({
  "pathway": ["A", "B"],
  "annotation": ["0011", "1100"],
  "sample": ["0011", "0101", "1001", "0110", "1010", "1100"]
})";
    likeligrid::GenotypeModel model(sst, 4u);
    std::cerr << model.calc_loglik({1.0, 1.0}) << std::endl;
    return 0;
}
