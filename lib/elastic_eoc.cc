#include "example.hh"

int main() {
  using namespace dealii;
  ElasticityExample::Parameter param({{"lambda",
                                       {
                                           1.,
                                       }},
                                      {"mu",
                                       {
                                           1.,
                                       }}});
  ElasticityEoc eoc(3, 5, param);
  auto res = eoc.run();
}
