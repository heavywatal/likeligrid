#include "pathtype.hpp"

#include <iostream>
#include <sstream>

int main() {
    std::stringstream sst;
    sst <<
R"(A B
0 1
1 0
1 1
0 2
)";
    likeligrid::PathtypeModel model(std::move(sst), 3u);
    std::cerr << model.calc_loglik({0.8, 1.2}) << std::endl;
    return 0;
}
