/*
GRIT (https://github.com/GRIT-RBSim/GRIT)
Licensed under the Apache-2.0 License.
Copyright (c) 2021, [Renyi Chen, Gongjie Li, Molei Tao]. All rights reserved.
*/

#include "math_func.hpp"
#include "matrix.hpp"
#include "system.hpp"
#include <cmath>

std::vector<ld> CustomizeData(const rb_sim::System &sys) {
  std::vector<ld> ret;

  // Add the first customized variable.
  rb_sim::Vec3<ld> pos = sys.GetPos(0);
  ld val1 = pos[0] * pos[1];
  ret.push_back(val1);

  // Add the second customized variable.
  auto mass0 = sys.GetMass(0);
  auto mass1 = sys.GetMass(1);
  ld val2 = mass0 * mass1;
  ret.push_back(val2);

  // ...

  return ret;
}
