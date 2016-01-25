#include "discretization.hh"

int main() {
  using namespace dealii;
  ElasticityEoc::Parameter param({{"lambda", {1.,}}, {"mu", {1.,}}});
  ElasticityEoc eoc(3, 5, param);
  auto res = eoc.run();
}
