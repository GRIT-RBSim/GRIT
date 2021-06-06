/*
GRIT (https://github.com/GRIT-RBSim/GRIT)
Licensed under the Apache-2.0 License.
Copyright (c) 2021, [Renyi Chen, Gongjie Li, Molei Tao]. All rights reserved.
*/

#include "math_func.hpp"

namespace rb_sim {
// calc the rotation matrix of rotation about $axis=(x, y, z) with angle $theta
// with right hand rule (unit: rad)
Matrix<ld> RotAngle2Mat(ld x, ld y, ld z, ld theta) {
  std::vector<ld> ret(9);

  ld len = sqrt(x * x + y * y + z * z);
  if (std::abs(len) != 0) {
    x = x / len;
    y = y / len;
    z = z / len;
  }
  ld C = cos(theta);
  ld S = sin(theta);
  ld t = ld(1.0) - C;

  ret[0] = t * x * x + C;
  ret[1] = t * x * y - S * z;
  ret[2] = t * x * z + S * y;
  ret[3] = t * x * y + S * z;
  ret[4] = t * y * y + C;
  ret[5] = t * y * z - S * x;
  ret[6] = t * x * z - S * y;
  ret[7] = t * y * z + S * x;
  ret[8] = t * z * z + C;

  return Matrix<ld>(ret, 3, 3);
}

Matrix<ld> RotAngle2Mat(const Vec3<ld> &axis, ld theta) {
  return RotAngle2Mat(axis.GetElement(1), axis.GetElement(2),
                      axis.GetElement(3), theta);
}

Matrix<ld> RotX(ld theta) {
  ld S = sin(theta);
  ld C = cos(theta);
  return Matrix<ld>(std::vector<ld>{ld(1.0), 0, 0, 0, C, -S, 0, S, C}, 3, 3);
}

Matrix<ld> RotY(ld theta) {
  ld S = sin(theta);
  ld C = cos(theta);
  return Matrix<ld>(std::vector<ld>{C, 0, S, 0, ld(1.0), 0, -S, 0, C}, 3, 3);
}

Matrix<ld> RotZ(ld theta) {
  ld S = sin(theta);
  ld C = cos(theta);
  return Matrix<ld>(std::vector<ld>{C, -S, 0, S, C, 0, 0, 0, ld(1.0)}, 3, 3);
}

void RotMat2Angle(const Matrix<ld> &mat, Matrix<ld> &axis, ld &theta) {
  ld a, b, c, d, e, f, g, h, i;
  a = mat.GetElement(1, 1);
  b = mat.GetElement(1, 2);
  c = mat.GetElement(1, 3);
  d = mat.GetElement(2, 1);
  e = mat.GetElement(2, 2);
  f = mat.GetElement(2, 3);
  g = mat.GetElement(3, 1);
  h = mat.GetElement(3, 2);
  i = mat.GetElement(3, 3);

  ld trace = a + e + i;
  theta = acos((trace - 1.0) / 2.0);

  ld x, y, z;

  x = (h - f);
  y = (c - g);
  z = (d - b);
  ld norm = sqrt(x * x + y * y + z * z);
  if (norm != 0) {
    x = x / norm;
    y = y / norm;
    z = z / norm;
  }
  axis = Matrix<ld>(std::vector<ld>{x, y, z}, 3, 1);
}

Vec3<ld> Rot(const Matrix<ld> &mat) {
  auto ret = (mat - mat.Transpose()).SkewInv();
  return Vec3<ld>(ret.GetElement(1, 1), ret.GetElement(2, 1),
                  ret.GetElement(3, 1));
}

std::vector<ld> PosVel2OrbElements(const Vec3<ld> &pos, const Vec3<ld> &vel,
                                   ld m1, ld m2) {

  auto r = pos;
  auto v = vel;
  ld mu = kG * (m1 + m2);

  Vec3<ld> K = {0, 0, 1};

  Vec3<ld> h = r.CrossProduct(v);
  Vec3<ld> n = K.CrossProduct(h);

  ld r_norm = r.Norm();
  ld v_norm = v.Norm();
  ld h_norm = h.Norm();
  ld n_norm = n.Norm();

  ld epsilon = v_norm * v_norm / 2.0 - mu / r_norm;
  ld a = -mu / 2.0 / epsilon;
  ld val = 1.0 + 2.0 * epsilon * h_norm * h_norm / mu / mu;
  if (val < 0)
    val = 0;
  ld e = sqrt(val);
  Vec3<ld> vec_e =
      r * (v_norm * v_norm / mu - 1.0 / r_norm) - v * ((r.DotProduct(v)) / mu);

  val = h.GetElement<3>() / h_norm;
  ld i;
  if (val >= -ld(1.0) && val <= ld(1.0))
    i = acos(val);
  else {
    if (val > ld(0.0))
      i = ld(0.0);
    else
      i = kPI;
  }

  val = n.GetElement<1>() / n_norm;
  ld Omega;
  if (val > -ld(1.0) && val < ld(1.0))
    Omega = acos(val);
  else {
    if (val > ld(0.0))
      Omega = ld(0.0);
    else
      Omega = kPI;
  }
  if (n.GetElement(2) < 0)
    Omega = 2 * kPI - Omega;

  val = n.DotProduct(vec_e) / n_norm / e;
  ld omega;
  if (val > -ld(1.0) && val < ld(1.0))
    omega = acos(val);
  else {
    if (val > ld(0.0))
      omega = ld(0.0);
    else
      omega = kPI;
  }
  if (vec_e.GetElement(3) < 0)
    omega = 2 * kPI - omega;

  val = vec_e.DotProduct(r) / e / r_norm;
  ld true_anomaly;
  if (val > -ld(1.0) && val < ld(1.0))
    true_anomaly = acos(val);
  else {
    if (val > ld(0.0))
      true_anomaly = ld(0.0);
    else
      true_anomaly = kPI;
  }
  if (r.DotProduct(v) < 0.0)
    true_anomaly = 2 * kPI - true_anomaly;

  return {a, e, i, Omega, omega, true_anomaly};
}

ld SlvKeplerEq(ld M, ld e) {
  // Newton-Raphson iteration with a simple starting estimate by Danby, 1987
  ld E = M + 0.85 * e;
  while (std::abs(M - E + e * sin(E)) > 1e-15) {
    E = E - (E - e * sin(E) - M) / (1 - e * cos(E));
  }
  return E;
}

void OrbElements2PosVel(const std::vector<ld> &orbs, ld T, Vec3<ld> &pos,
                        Vec3<ld> &vel) {

  ld a = orbs.at(0);
  ld e = orbs.at(1);
  ld i = orbs.at(2);
  ld Omega = orbs.at(3);
  ld omega = orbs.at(4);
  ld M = orbs.at(5);

  ld n = ld(2.0) * kPI / T;

  ld E = SlvKeplerEq(M, e);
  ld r = a * (ld(1.0) - e * cos(E));
  ld V = atan(sqrt((ld(1.0) + e) / (ld(1.0) - e)) * tan(E / ld(2.0))) * ld(2.0);

  Matrix<ld> R = RotZ(Omega) * RotX(i) * RotZ(omega);

  Matrix<ld> P = R.getCol(1);
  Matrix<ld> Q = R.getCol(2);

  auto tmp_pos =
      R * Matrix<ld>(std::vector<ld>{r * cos(V), r * sin(V), ld(0.0)}, 3, 1);
  auto tmp_vel = n * a * a / r * (Q * sqrt(1 - e * e) * cos(E) - P * sin(E));
  pos = Vec3<ld>{tmp_pos.GetElement(1, 1), tmp_pos.GetElement(2, 1),
                 tmp_pos.GetElement(3, 1)};
  vel = Vec3<ld>{tmp_vel.GetElement(1, 1), tmp_vel.GetElement(2, 1),
                 tmp_vel.GetElement(3, 1)};
}

ld TrueAnomaly2MeanAnomaly(ld nu, ld e) {
  ld E = TrueAnomaly2EccentricAnomaly(nu, e); 
  ld M = EccentricAnomaly2MeanAnomaly(E, e);
  return M;
}

ld TrueAnomaly2EccentricAnomaly(ld nu, ld e) {
  ld E;
  E = atan2(sqrt(1 - e * e) * sin(nu), e + cos(nu));
  return E;
}

ld MeanAnomaly2EccentricAnomaly(ld M, ld e) {
  ld E = SlvKeplerEq(M, e);
  return E;
}

ld MeanAnomaly2TrueAnomaly(ld M, ld e) {
  ld E = MeanAnomaly2EccentricAnomaly(M, e);
  ld nu = EccentricAnomaly2TrueAnomaly(E, e);
  return nu;
}

ld EccentricAnomaly2TrueAnomaly(ld E, ld e) {
  ld nu;
  nu = atan2(sqrt(1 - e * e) * sin(E), cos(E) - e);
  return nu;
}

ld EccentricAnomaly2MeanAnomaly(ld E, ld e) {
  ld M;
  M = E - e * sin(E);
  return M;
}

void RotMat2EulerAngles(Matrix<ld> R, ld &phi, ld &theta, ld &psi) {
  ld a, b, c, d, e, f, g, h, i;
  a = R.GetElement(1, 1);
  b = R.GetElement(1, 2);
  c = R.GetElement(1, 3);
  d = R.GetElement(2, 1);
  e = R.GetElement(2, 2);
  f = R.GetElement(2, 3);
  g = R.GetElement(3, 1);
  h = R.GetElement(3, 2);
  i = R.GetElement(3, 3);

  if (std::abs(std::abs(a) - 1.0) > kEps) {
    theta = acos(a);
    psi = atan2(c / sin(theta), -b / sin(theta));
    phi = atan2(g / sin(theta), d / sin(theta));
  } else {
    psi = 0.0;
    if (std::abs(a - 1.0) < kEps) {
      theta = 0;
      phi = atan2(h, e) - psi;
    } else {
      theta = kPI;
      phi = psi + atan2(h, -e);
    }
  }

}

void RotMat2EulerAngles(Matrix<ld> R, Vec3<ld> &euler_angle) {
  ld phi, psi, theta;
  RotMat2EulerAngles(R, phi, theta, psi);
  euler_angle = Vec3<ld>(phi, theta, psi);
}

void EulerAngles2RotMat(Matrix<ld> &R, ld phi, ld theta, ld psi) {
  R = RotX(phi) * RotZ(theta) * RotX(psi);
}

void EulerAngles2RotMat(Matrix<ld> &R, const Vec3<ld> &euler_angle) {
  EulerAngles2RotMat(R, euler_angle.GetElement(1), euler_angle.GetElement(2),
                     euler_angle.GetElement(3));
}

void KeplerSlvDanby(ld dt, ld mu, Vec3<ld>& p, Vec3<ld>& v) {

    ld r_norm = p.Norm();
    ld v_norm_sqr = v.NormSquared();
    ld u = p.DotProduct(v);
    ld alpha = 2.0 * mu / r_norm - v_norm_sqr;

    // initial guess `s`
    ld s = 0.0;
    if (alpha > 0) {
      if (dt / r_norm < 0.2) {
        s = dt / r_norm - (dt * dt * u) / (2.0 * r_norm * r_norm * r_norm);
      } else {
        ld a = mu / alpha;
        ld en = sqrt(mu / (a * a * a));
        ld ec = 1.0 - r_norm / a;
        ld es = u / (en * a * a);
        ld e = sqrt(ec * ec + es * es);
        ld y = en * dt - es;
        ld sigma = 1.0;
        if (es * cos(y) + ec * sin(y) < 0)
          sigma = -1.0;
        ld x = y + sigma * 0.85 * e;
        s = x / sqrt(alpha);
      }
    } else {
      ld div = (mu - alpha * r_norm) / 6.0;
      ld a2 = 0.5 * u / div;
      ld a1 = r_norm / div;
      ld a0 = -dt / div;

      ld q = (a1 - a2 * a2 / 3.0) / 3.0;
      ld r = (a1 * a2 - 3.0 * a0) / 6.0 - (a2 * 3.0) / 27.0;
      ld tmp = q * q * q + r * r;

      if (tmp > 0) {
        ld sq = sqrt(tmp);

        ld p1, p2;
        if ((r + sq) < 0) {
          p1 = -pow(-(r + sq), 1.0 / 3.0);
        } else {
          p1 = pow(r + sq, 1.0 / 3.0);
        }
        if ((r - sq) < 0) {
          p2 = -pow(-(r - sq), 1.0 / 3.0);
        } else {
          p2 = pow(r - sq, 1.0 / 3.0);
        }

        s = p1 + p2 - a2 / 3.0;
      }
    }
    // drift using Laguerre's method
    int n_max = 400;
    ld ln = 5.0;
    int iter = 0;
    ld f, fp, fpp, c0, c1, c2, c3;
    while (iter < n_max) {
      ld x = s * s * alpha;
      Stumpff(x, c0, c1, c2, c3);
      c1 = c1 * s;
      c2 = c2 * s * s;
      c3 = c3 * s * s * s;
      f = r_norm * c1 + u * c2 + mu * c3 - dt;
      fp = r_norm * c0 + u * c1 + mu * c2;
      fpp = (-40.0 * alpha + mu) * c1 + u * c0;
      ld ds = -ln * f /
              (fp + copysign(1.0, fp) *
                        sqrt(std::abs((ln - 1.0) * (ln - 1.0) * fp * fp -
                                      (ln - 1.0) * ln * f * fpp)));
      s = s + ds;
      ld fdt = f / dt;
      if (std::abs(fdt) < kEps)
        break;
      iter++;
    }

    f = 1.0 - (mu / r_norm) * c2;
    ld g = dt - mu * c3;
    ld fdot = -(mu / (fp * r_norm)) * c1;
    ld gdot = 1.0 - (mu / fp) * c2;

    auto p1 = p * f + v * g;
    auto v1 = p * fdot + v * gdot;
    p = p1;
    v = v1;
}


void Stumpff(ld x, ld& c0, ld& c1, ld& c2, ld& c3) {
    ld xm = 0.1;
    int k = 0;
    while (std::abs(x)>xm) { x*=0.25; k++; }

    c2 = (1 -
          x *
              (1 -
               x *
                   (1 -
                    x *
                        (1 - x * (1 - x * (1 - x * (1.0 - x / 182.0)) / 132.0) /
                                 90.0) /
                        56.0) /
                   30.0) /
              12.0) /
         2.0;

    c3 = (1 -
          x *
              (1 -
               x *
                   (1 -
                    x *
                        (1 - x * (1 - x * (1 - x * (1.0 - x / 210.0)) / 156.0) /
                                 110.0) /
                        72.0) /
                   42.0) /
              20.0) /
         6.0;

    c1 = 1.0 - x * c3;
    c0 = 1.0 - x * c2;

    for (int i = 0; i < k; i++) {
      c3 = (c2 + c0 * c3) * 0.25;
      c2 = c1 * c1 * 0.5;
      c1 = c0 * c1;
      c0 = 2.0 * c0 * c0 - 1.0;
      x = x * 4.0;
    }
}

} // namespace rb_sim
