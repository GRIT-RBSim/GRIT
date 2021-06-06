/*
GRIT (https://github.com/GRIT-RBSim/GRIT)
Licensed under the Apache-2.0 License.
Copyright (c) 2021, [Renyi Chen, Gongjie Li, Molei Tao]. All rights reserved.
*/

#ifndef MATHFUNC_H_
#define MATHFUNC_H_

#include "matrix.hpp"
#include <array>
#include <cmath>
#include <vector>

using std::sin;
using std::cos;
using std::acos;

namespace rb_sim {
// pi
const ld kPI = acos(-1.0L);
const ld k2PI = 2.0*kPI;
// The gravational constant 4*pi^2/(365.2422)^2 [AU^3 days^{-2} M_sun^{-1}]
const ld kG = ld(4.0) * kPI * kPI;
// Days of a year: for converting the time unit
const ld kDaysOfYear = ld(365.25);
//  const ld kDaysOfYear = ld(365.25636);
// Astronomical Unit in meters: 1 AU = 1.495978707E+11 meters
const ld kAU = ld(1.495978707E+11);
// The threshold: used as small values in calculating orbital elements
const ld kEps = ld(1e-15);
// The Earth Tilt Angle
const ld kEarthTilt = ld(23.44) * kPI / ld(180.0);
// Three standard directions in Cartesian coordinates
const Vec3<ld> kE1(ld(1.0), ld(0.0), ld(0.0));
const Vec3<ld> kE2(ld(0.0), ld(1.0), ld(0.0));
const Vec3<ld> kE3(ld(0.0), ld(0.0), ld(1.0));
// Identity matrix
const Matrix<ld> kZero =
    Matrix<ld>(std::vector<ld>{0, 0, 0, 0, 0, 0, 0, 0, 0}, 3, 3);
const Matrix<ld> kIdentity =
    Matrix<ld>(std::vector<ld>{1, 0, 0, 0, 1, 0, 0, 0, 1}, 3, 3);

// Return the rotation matrix of rotating around the axis [`x`, `y`, `z`] with
// angle `theta`.
Matrix<ld> RotAngle2Mat(ld x, ld y, ld z, ld theta);
// Return the rotation matrix of rotating around `axis` with angle `theta`.
Matrix<ld> RotAngle2Mat(const Vec3<ld> &axis, ld theta);

// Return the rotation matrix of rotating around X,Y,Z axes with angle `theta`.
Matrix<ld> RotX(ld theta);
Matrix<ld> RotY(ld theta);
Matrix<ld> RotZ(ld theta);

// Update `axis` and `theta` according to rotation matrix mat such that rotating
// around `axis` with angle `theta` can be represented by the rotation matrix
// `mat`.
void RotMat2Angle(const Matrix<ld> &mat, Matrix<ld> &axis, ld &theta);

Vec3<ld> Rot(const Matrix<ld> &mat);

// Return the orbital elements from the position `pos` and the velocity `vel`;
// orbital elements order: (a, e, i, Omega, omega, M) with angle in degree.
// `m1` and `m2` are masses of two bodies.
// Reference:
// https://en.wikibooks.org/wiki/Astrodynamics/Classical_Orbit_Elements Time
// unit: yr Angle unit: rad
std::vector<ld> PosVel2OrbElements(const Vec3<ld> &pos, const Vec3<ld> &vel,
                                   ld m1, ld m2);
// Return E: Solve kepler equation $M = E - e \sin E$ with M && e constants
ld SlvKeplerEq(ld M, ld e);
// Update the position `pos` and the velocity `vel` from orbital elements `orbs`
void OrbElements2PosVel(const std::vector<ld> &orbs, ld T, Vec3<ld> &pos,
                        Vec3<ld> &vel);
// Conversion among True Anomaly, EccentricAnomaly and Mean Anomaly
ld TrueAnomaly2MeanAnomaly(ld nu, ld e);
ld TrueAnomaly2EccentricAnomaly(ld nu, ld e);
ld MeanAnomaly2EccentricAnomaly(ld M, ld e);
ld MeanAnomaly2TrueAnomaly(ld M, ld e);
ld EccentricAnomaly2TrueAnomaly(ld E, ld e);
ld EccentricAnomaly2MeanAnomaly(ld E, ld e);

// Update the Eular angles (`phi`, `theta`, `psi`) in order X Z X from the
// rotation matrix `R`
void RotMat2EulerAngles(Matrix<ld> R, ld &phi, ld &theta, ld &psi);
// Update the Eular angles `euler_angles` in order X Z X from the rotation
// matrix `R`
void RotMat2EulerAngles(Matrix<ld> R, Vec3<ld> &euler_angle);
// Update the rotation matrix `R` from Eular angles (`phi`, `theta`, `psi`) in
// order X Z X
void EulerAngles2RotMat(Matrix<ld> &R, ld phi, ld theta, ld psi);
// Update the rotation matrix `R` from Eular angles `euler_angles` in order X Z
// X
void EulerAngles2RotMat(Matrix<ld> &R, const Vec3<ld> &euler_angles);

// Danby "Fundamentals of celestial mechanics"
// TODO: SOTA kepler orbit solvers
void KeplerSlvDanby(ld dt, ld mu, Vec3<ld>& p, Vec3<ld>& v);
void Stumpff(ld x, ld& c0, ld& c1, ld& c2, ld& c3);

template <typename T>
inline std::ostream &operator<<(std::ostream &os, const std::vector<T> &c) {
  for (const auto &e : c)
    os << e << '\t';
  return os;
}
} // namespace rb_sim

#endif
