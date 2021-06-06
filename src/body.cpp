/*
GRIT (https://github.com/GRIT-RBSim/GRIT)
Licensed under the Apache-2.0 License.
Copyright (c) 2021, [Renyi Chen, Gongjie Li, Molei Tao]. All rights reserved.
*/

#include "body.hpp"
#include <cassert>

namespace rb_sim {

Body::Body(const json &j)
    : name(j.value<std::string>("name", "")), mass(j.value<ld>("mass", -1)),
      radius(j.value<ld>("radius", -1) / kAU * 1000),
      rigid(j.value<bool>("rigid", false)) {
  qc = Vec3<ld>(0.0, 0.0, 0.0);
  pc = Vec3<ld>(0.0, 0.0, 0.0);
  assert(name != "" && "Body parameters: need body name.");
  assert(mass >= 0 && "Body parameters: need mass.");
  assert(radius >= 0 && "Body parameters: need radius.");
  // First search for orbital elements, if there is no input orbital elements,
  // search for position and velocity.
  if (j.contains("orbital_elements")) {
    assert(j.contains("center_mass") &&
           "Body parameters: need center mass to convert the orbital elements "
           "to position and velocity.");
    // Convert orbital elements to position and velocity.
    ld CM = j["center_mass"];
    std::vector<ld> orb = j["orbital_elements"];
    ld T = 2 * kPI * sqrt(orb[0] * orb[0] * orb[0] / kG / (CM + mass));
    OrbElements2PosVel(j["orbital_elements"], T, q, p);
    auto tmp = PosVel2OrbElements(q, p, mass, CM);
    p = p * mass;
  } else {
    assert((j.contains("position") &&
            (j.contains("velocity") || j.contains("_p"))) &&
           "Body parameters: need initial (position, velocity) or initial "
           "orbital elements.");
    q = j["position"].get<Vec3<ld>>();
    if (j.contains("_p"))
      p = j["_p"];
    else
      p = j["velocity"].get<Vec3<ld>>() * kDaysOfYear * mass;
  }

  // For rigid body, more properties.
  if (rigid) {
    if (j.contains("symmetric"))
      symmetric = j["symmetric"];
    else if (j.contains("oblateness")) {
      symmetric = true;
    } else
      symmetric = false;
    if (j.contains("oblateness")) {
      ld oblateness = j["oblateness"];
      ld axis_a = radius / pow(1 - oblateness, ld(1.0) / ld(3.0));
      ld axis_b = axis_a;
      ld axis_c = axis_a * (ld(1.0) - oblateness);
      InitInertiaTensorFromSemiAxes(axis_a, axis_b, axis_c);
    } else if (j.contains("inertia_tensor")) {
      Vec3<ld> I = j["inertia_tensor"];
      I1 = I[0];
      I2 = I[1];
      I3 = I[2];
    } else if (j.contains("semi_axes")) {
      Vec3<ld> axes = j["semi_axes"];
      ld axis_a = axes[0] / kAU * ld(1000);
      ld axis_b = axes[1] / kAU * ld(1000);
      ld axis_c = axes[2] / kAU * ld(1000);
      InitInertiaTensorFromSemiAxes(axis_a, axis_b, axis_c);
    } else {
      assert(false && "Body parameters: need inertia variables.");
    }
    if (symmetric)
      delta = 1 / I2 - 1 / I1;
    else
      delta = 0;
    I = Matrix<ld>(std::vector<ld>{I1, 0, 0, 0, I2, 0, 0, 0, I3}, 3, 3);
    Iinv = Matrix<ld>(
        std::vector<ld>{ld(1) / I1, 0, 0, 0, ld(1) / I2, 0, 0, 0, ld(1) / I3},
        3, 3);
    traceI = I1 + I2 + I3;

    assert((j.contains("euler_angles") || j.contains("_R")) &&
           "Body parameters: need euler_angles.");
    auto euler_angles = j["euler_angles"].get<Vec3<ld>>();
    if (j.contains("_R")) {
      R = j["_R"].get<Matrix<ld>>();
    } else
      EulerAngles2RotMat(R, euler_angles[0], euler_angles[1], euler_angles[2]);
    assert((j.contains("angular_velocity") || j.contains("_pi")) &&
           "Body parameters: need rotational angular velocity.");
    if (j.contains("_pi"))
      pi = j["_pi"];
    else
      pi = I * j["angular_velocity"].get<Vec3<ld>>();
  }

  if (j.contains("time_lag"))
    time_lag = j["time_lag"];
  else
    time_lag = -1;
  if (j.contains("Q"))
    Q = j["Q"];
  else
    Q = -1;
}

void Body::InitInertiaTensorFromSemiAxes(ld ra, ld rb, ld rc) {
  I1 = mass / ld(5.0) * (rb * rb + rc * rc);
  I2 = mass / ld(5.0) * (ra * ra + rc * rc);
  I3 = mass / ld(5.0) * (ra * ra + rb * rb);
}

void Body::GetSemiAxes(ld &ra, ld &rb, ld &rc) const {
  ld rsqr = (I1 + I2 + I3) / mass * ld(5.0) / 2;
  ra = sqrt(rsqr - I1 / mass * ld(5.0));
  rb = sqrt(rsqr - I2 / mass * ld(5.0));
  rc = sqrt(rsqr - I3 / mass * ld(5.0));
}

Vec3<ld> Body::GetPos() const { return q; }

Vec3<ld> Body::GetVel() const { return p / mass; }

Vec3<ld> Body::GetAngMomentum() const { return pi; }

Matrix<ld> Body::GetRotMat() const { return R; }

void Body::GetInertia(int &retI1, int &retI2, int &retI3) const {
  retI1 = I1;
  retI2 = I2;
  retI3 = I3;
}

std::ostream &operator<<(std::ostream &os, const Body &s) {
  auto velocity = s.p / s.mass / kDaysOfYear;
  auto radius = s.radius * kAU / 1000;

  if (s.rigid) {
    auto euler_angle = Vec3<ld>(0);
    RotMat2EulerAngles(s.R, euler_angle);
    auto ang_vel = s.I.inverse() * s.pi;
    os << s.name << '\t' << s.mass << '\t' << radius << '\t' << s.q << '\t'
       << velocity << '\t' << ang_vel << '\t' << euler_angle << '\t'
       << s.symmetric << '\t' << s.I1 << '\t' << s.I2 << '\t' << s.I3 << '\n';
  } else {
    os << s.name << '\t' << s.mass << '\t' << s.q << '\t' << velocity << '\t'
       << '\n';
  }
  os << s.time_lag << '\t' << s.Q << '\n';
  return os;
}

std::istream &operator>>(std::istream &is, Body &s) {
  auto velocity = s.p / s.mass / kDaysOfYear;
  auto radius = s.radius * kAU / 1000;

  auto name = s.name;
  auto mass = s.mass;

  if (s.rigid) {
    Vec3<ld> ang_vel;
    Vec3<ld> euler_angle;
    is >> name >> mass >> radius >> s.q >> velocity >> ang_vel >> euler_angle >>
        s.symmetric >> s.I1 >> s.I2 >> s.I3;
    s.pi = s.I * ang_vel;
    EulerAngles2RotMat(s.R, euler_angle);
  } else {
    is >> name >> mass >> s.q >> velocity;
  }

  is >> s.time_lag >> s.Q;

  if (name != s.name)
    throw std::runtime_error(std::string("name is not ") + s.name + "!");

  s.p = velocity * s.mass * kDaysOfYear;

  return is;
}

void to_json(json &j, const Body &s) {
  j = json::object();
  j["name"] = s.name;
  j["mass"] = s.mass;
  j["radius"] = s.radius * kAU / 100;
  j["position"] = s.q;
  j["velocity"] = s.p / s.mass / kDaysOfYear;
  j["_p"] = s.p;
  j["rigid"] = s.rigid;
  j["time_lag"] = s.time_lag;
  j["Q"] = s.Q;
  if (s.rigid) {
    j["symmetric"] = s.symmetric;
    j["inertia_tensor"] = Vec3<ld>(s.I1, s.I2, s.I3);
    j["_R"] = s.R;
    j["_pi"] = s.pi;
    auto euler_angle = Vec3<ld>(0);
    RotMat2EulerAngles(s.R, euler_angle);
    auto angVel = s.I.inverse() * s.pi;
    j["angular_velocity"] = angVel;
    j["euler_angles"] = euler_angle;
  }
}

} // namespace rb_sim
