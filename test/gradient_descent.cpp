#include "gradient_descent.hpp"

#include "wtl/iostr.hpp"

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
    likeligrid::GradientDescent searcher(sst, 4, {0, 1}, false);
    searcher.run(std::cout);
    std::cout << *searcher.const_max_iterator() << std::endl;
    return 0;
}
