#include "gridsearch.hpp"

#include <sstream>

int main() {
    std::stringstream sst;
    sst <<
R"({
  "pathway": ["A", "B"],
  "annotation": ["0011", "1100"],
  "sample": ["0011", "0101", "1001", "0110", "1010", "1100"]
})";
    likeligrid::GridSearch searcher(sst, 4u, {0, 1});
    searcher.run_cout();
    return 0;
}
