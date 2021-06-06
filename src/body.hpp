/*
GRIT (https://github.com/GRIT-RBSim/GRIT)
Licensed under the Apache-2.0 License.
Copyright (c) 2021, [Renyi Chen, Gongjie Li, Molei Tao]. All rights reserved.
*/

#ifndef BODY_H_
#define BODY_H_

#include "math_func.hpp"
#include <string>
#include <vector>

namespace rb_sim {

class Body {
  friend class System;

private:
  // constants
  const std::string name;
  const ld mass;
  const ld radius;
  const bool rigid;
  std::string body_type;
  bool symmetric;
  ld I1, I2, I3; // Inertia
  // 1/I2-1/I1
  ld delta;
  Matrix<ld> I, Iinv;
  ld traceI;

  // variables: depending on time
  Matrix<ld> R;
  Vec3<ld> q;  // position
  Vec3<ld> p;  // linear momentum
  Vec3<ld> pi; // angularMomentum

  // tide constants
  ld time_lag, Q;

  // compensated summation
  Vec3<ld> qc, pc;

  std::string key;

public:
  Body(const json &j);
  void InitInertiaTensorFromSemiAxes(ld ra, ld rb, ld rc);
  // See whether the tidal force will be applied to this body.
  bool CheckTidalParameters() const { return (time_lag >= 0) && (Q >= 0); }
  Vec3<ld> GetPos() const;
  Vec3<ld> GetVel() const;
  Vec3<ld> GetAngMomentum() const;
  void GetSemiAxes(ld &ra, ld &rb, ld &rc) const;
  Matrix<ld> GetRotMat() const;
  void GetInertia(int &, int &, int &) const;
  friend std::ostream &operator<<(std::ostream &os, const Body &s);
  friend std::istream &operator>>(std::istream &is, Body &s);
  std::ostream &Dump(std::ostream &os) const {
    if (body_type == "rigid_body") {
      os << R << '\n' << p << '\n' << q << '\n' << pi << '\n';
    } else {
      os << p << '\n' << q << '\n';
    }
    return os;
  }

  friend void to_json(json &j, const Body &s);
};
} // namespace rb_sim

#endif
