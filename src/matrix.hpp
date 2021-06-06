#ifndef MATRIX_H_
#define MATRIX_H_

#include "json.hpp"
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <functional>
#include <initializer_list>
#include <iostream>
#include <vector>

typedef double ld;

namespace rb_sim {

using json = nlohmann::basic_json<
    // template <typename U, typename V, typename... Args> class ObjectType =
    std::map,
    // template <typename U, typename... Args> class ArrayType =
    std::vector,
    // class StringType =
    std::string,
    // class BooleanType =
    bool,
    // class NumberIntegerType =
    std::int64_t,
    // class NumberUnsignedType =
    std::uint64_t,
    // class NumberFloatType =
    ld,
    // template <typename U> class AllocatorType =
    std::allocator,
    // template <typename T, typename SFINAE = void> class JSONSerializer =
    nlohmann::adl_serializer,
    // class BinaryType =
    std::vector<std::uint8_t>>;

template <typename T> class Matrix;
template <typename T> class Vec3;
template <typename T> class Vec3 {
  friend class Matrix<T>;

public:
  Vec3() : m{{0, 0, 0}} {}
  Vec3(T t) : m{{t, t, t}} {}
  Vec3(T t1, T t2, T t3) : m{{t1, t2, t3}} {}
  Vec3(std::array<T, 3> m0) : m(m0) {}
  std::array<T, 3> Get() const { return m; }
  Vec3 operator+(const Vec3 &rhs) const {
    return {m[0] + rhs.m[0], m[1] + rhs.m[1], m[2] + rhs.m[2]};
  }
  Vec3 operator-(const Vec3 &rhs) const {
    return {m[0] - rhs.m[0], m[1] - rhs.m[1], m[2] - rhs.m[2]};
  }
  Vec3 operator+(const Matrix<T> &rhs) const {
    assert(rhs.nr == 3 && rhs.nc == 1);
    return {m[0] + rhs.m[0], m[1] + rhs.m[1], m[2] + rhs.m[2]};
  }
  Vec3 operator-(const Matrix<T> &rhs) const {
    assert(rhs.nr == 3 && rhs.nc == 1);
    return {m[0] - rhs.m[0], m[1] - rhs.m[1], m[2] - rhs.m[2]};
  }

  Vec3 operator+(const T &scalar) const {
    return {m[0] + scalar, m[1] + scalar, m[2] + scalar};
  }
  Vec3 operator-(const T &scalar) const {
    return {m[0] - scalar, m[1] - scalar, m[2] - scalar};
  }
  Vec3 operator*(const T &scalar) const {
    return {m[0] * scalar, m[1] * scalar, m[2] * scalar};
  }
  Vec3 operator/(const T &scalar) const {
    return {m[0] / scalar, m[1] / scalar, m[2] / scalar};
  }
  template <typename TT> Vec3 &operator+=(const TT &rhs) {
    return *this = (this->operator+(rhs));
  }
  template <typename TT> Vec3 &operator-=(const TT &rhs) {
    return *this = (this->operator-(rhs));
  }
  template <typename TT> Vec3 &operator*=(const TT &rhs) {
    return *this = (this->operator*(rhs));
  }

  Vec3 CrossProduct(const Vec3 &v) const {
    return {m[1] * v.m[2] - m[2] * v.m[1], m[2] * v.m[0] - m[0] * v.m[2],
            m[0] * v.m[1] - m[1] * v.m[0]};
  }

  Matrix<T> Skew() {
    return Matrix<T>(
        std::array<T, 9>{0, -m[2], m[1], m[2], 0, -m[0], -m[1], m[0], 0}, 3, 3);
  }

  T &operator[](int index) { return m[index]; }
  template <int n> T GetElement() const;
  T GetElement(int n) const { return m[n - 1]; }
  std::array<T, 3> GetM() const { return m; }
  T DotProduct(const Vec3 &rhs) const {
    return m[0] * rhs.m[0] + m[1] * rhs.m[1] + m[2] * rhs.m[2];
  }
  T Norm() const { return sqrt(m[0] * m[0] + m[1] * m[1] + m[2] * m[2]); }
  void Normalize() { *this = this->operator/(Norm()); }
  T NormSquared() const { return (m[0] * m[0] + m[1] * m[1] + m[2] * m[2]); }
  std::ostream &operator>>(std::ostream &os) const {
    return os << m[0] << '\t' << m[1] << '\t' << m[2];
  }
  std::istream &operator<<(std::istream &is) {
    return is >> m[0] >> m[1] >> m[2];
  }

private:
  std::array<T, 3> m;
};

template <typename T> void to_json(json &j, const Vec3<T> &m) { j = m.Get(); }
template <typename T> void from_json(const json &j, Vec3<T> &m) {
  m = Vec3<T>(j.get<std::array<T, 3>>());
}
template <>
template <>
inline long double Vec3<long double>::GetElement<3>() const {
  return m[2];
}
template <> template <> inline double Vec3<double>::GetElement<3>() const {
  return m[2];
}
template <> template <> inline float Vec3<float>::GetElement<3>() const {
  return m[2];
}
template <>
template <>
inline long double Vec3<long double>::GetElement<2>() const {
  return m[1];
}
template <> template <> inline double Vec3<double>::GetElement<2>() const {
  return m[1];
}
template <> template <> inline float Vec3<float>::GetElement<2>() const {
  return m[1];
}
template <>
template <>
inline long double Vec3<long double>::GetElement<1>() const {
  return m[0];
}
template <> template <> inline double Vec3<double>::GetElement<1>() const {
  return m[0];
}
template <> template <> inline float Vec3<float>::GetElement<1>() const {
  return m[0];
}

template <typename T>
inline std::ostream &operator<<(std::ostream &os, const Vec3<T> &v) {
  return v >> os;
}

template <typename T>
inline std::istream &operator>>(std::istream &is, Vec3<T> &v) {
  return v << is;
}

template <typename T> class Matrix {
  friend class Vec3<T>;
  template <typename U> friend void to_json(json &j, const Matrix<U> &);
  template <typename U> friend void from_json(const json &j, Matrix<U> &);

public:
  Matrix() {}
  Matrix(T t) : m{{t, t, t, t, t, t, t, t, t}} {}
  Matrix(const Matrix &rhs) = default;
  Matrix(Matrix &&rhs) : m(std::move(rhs.m)), nc(rhs.nc), nr(rhs.nr) {}
  Matrix &operator=(const Matrix &rhs) = default;
  Matrix &operator=(Matrix &&rhs) {
    m = std::move(rhs.m);
    nc = rhs.nc;
    nr = rhs.nr;
    return *this;
  }
  Matrix(T t, size_t nr_, size_t nc_) : nc(nc_), nr(nr_) {
    for (int i = 0; i < nc * nr; i++)
      m[i] = t;
  }
  Matrix(const std::vector<T> &rhs, size_t nr_, size_t nc_) : nc(nc_), nr(nr_) {
    for (unsigned int i = 0; i < nc * nr; i++)
      m[i] = rhs[i];
  }
  Matrix(const std::array<T, 9> &rhs, size_t nr_, size_t nc_)
      : m(rhs), nc(nc_), nr(nr_) {}

  Matrix operator+(const Matrix &rhs) const {
    Matrix ret(nr, nc);
    std::transform(m.begin(), m.end(), rhs.m.begin(), ret.m.begin(),
                   std::plus<T>());
    return ret;
  }

  Matrix operator-(const Matrix &rhs) const {
    Matrix ret(nr, nc);
    std::transform(m.begin(), m.end(), rhs.m.begin(), ret.m.begin(),
                   std::minus<T>());
    return ret;
  }

  Matrix operator+(const T &scalar) const {
    Matrix ret(nr, nc);
    std::transform(m.begin(), m.end(), ret.m.begin(), [scalar](const T &elem) {
      return std::plus<T>()(elem, scalar);
    });
    return ret;
  }

  Matrix operator-(const T &scalar) const {
    Matrix ret(nr, nc);
    std::transform(m.begin(), m.end(), ret.m.begin(), [scalar](const T &elem) {
      return std::minus<T>()(elem, scalar);
    });
    return ret;
  }

  Matrix operator*(const T &scalar) const {
    Matrix ret(nr, nc);
    std::transform(m.begin(), m.end(), ret.m.begin(), [scalar](const T &elem) {
      return std::multiplies<T>()(elem, scalar);
    });
    return ret;
  }

  Matrix operator*(const Matrix &rhs) const {
    if (nr == 1 && nc == 1) {
      return rhs * m[0];
    }
    if (rhs.nr == 1 && rhs.nc == 1) {
      return *this * rhs.m[0];
    }
    std::array<T, 9> ret = {{0}};
    assert(nc == rhs.nr);
    for (auto i = 0u; i < nr; i++)
      for (auto j = 0u; j < rhs.nc; j++)
        for (auto k = 0u; k < nc; k++)
          ret[i * rhs.nc + j] += m[i * nc + k] * rhs.m[k * rhs.nc + j];
    return Matrix(ret, nr, rhs.nc);
  }

  Vec3<T> operator*(const Vec3<T> &rhs) const {
    assert(nc == 3 && nr == 3);
    return Vec3<T>{m[0] * rhs.m[0] + m[1] * rhs.m[1] + m[2] * rhs.m[2],
                   m[3] * rhs.m[0] + m[4] * rhs.m[1] + m[5] * rhs.m[2],
                   m[6] * rhs.m[0] + m[7] * rhs.m[1] + m[8] * rhs.m[2]};
  }

  Matrix operator/(const T &scalar) const {
    Matrix ret(nr, nc);
    std::transform(m.begin(), m.end(), ret.m.begin(), [scalar](const T &elem) {
      return std::divides<T>()(elem, scalar);
    });
    return ret;
  }

  T GetElement(size_t i, size_t j) const { return m[(i - 1) * nc + j - 1]; }

  Matrix<T> getCol(size_t col) const {
    Matrix<T> ret;
    for (auto i = 0u; i < nr; i++) {
      ret.m[i] = m[i * nc + col - 1];
    }
    ret.nc = 1;
    ret.nr = nr;
    return ret;
  }

  explicit operator T() const {
    assert(nr == 1 && nc == 1);
    return m[0];
  }

  std::ostream &operator>>(std::ostream &os) const {
    for (auto i = 0u; i < nr; i++) {
      for (auto j = 0u; j < nc; j++) {
        os << m[i * nc + j] << '\t';
      }
      // os << '\n';
    }

    return os;
  }

  Matrix Transpose() const {
    Matrix ret(nc, nr);
    for (auto i = 0u; i < nr; i++) {
      for (auto j = 0u; j < nc; j++) {
        ret.m[j * nr + i] = m[i * nc + j];
      }
    }
    return ret;
  }

  T Trace() const {
    T ret = m[0] + m[4] + m[8];
    return ret;
  }

  Matrix inverse() const {
    assert(nr == 3 && nc == 3);
    T a1 = m[0];
    T a2 = m[1];
    T a3 = m[2];
    T b1 = m[3];
    T b2 = m[4];
    T b3 = m[5];
    T c1 = m[6];
    T c2 = m[7];
    T c3 = m[8];

    T det = a1 * b2 * c3 - a1 * b3 * c2 - a2 * b1 * c3 + a2 * b3 * c1 +
            a3 * b1 * c2 - a3 * b2 * c1;

    Matrix ret(std::array<T, 9>{b2 * c3 - c2 * b3, a3 * c2 - a2 * c3,
                                a2 * b3 - a3 * b2, b3 * c1 - b1 * c3,
                                a1 * c3 - a3 * c1, a3 * b1 - a1 * b3,
                                b1 * c2 - b2 * c1, a2 * c1 - a1 * c2,
                                a1 * b2 - a2 * b1},
               3, 3);

    return ret * (T(1) / det);
  }

  T Norm() const {
    T ret = 0;
    for (auto i = 0; i < nc * nr; i++) {
      ret += m[i] * m[i];
    }

    return sqrt(ret);
  }

  Matrix Skew() {
    assert(nr == 3 && nc == 1);
    // return Matrix(std::array<T,9>{0, m[2], -m[1], -m[2], 0, m[0], m[1],
    // -m[0], 0}, 3, 3);
    return Matrix(
        std::array<T, 9>{0, -m[2], m[1], m[2], 0, -m[0], -m[1], m[0], 0}, 3, 3);
  }

  Matrix SkewInv() {
    assert(nr == 3 && nc == 3);
    // return Matrix(std::array<T,9>{m[5], m[6], m[1]}, 3, 1);
    return Matrix(std::array<T, 9>{-m[5], -m[6], -m[1]}, 3, 1);
  }

  Matrix<T> CrossProduct(const Matrix<T> &v) const {
    assert(nr == 3 && v.nr == 3);
    std::array<T, 9> ret = {0};
    ret[0] = m[1] * v.m[2] - m[2] * v.m[1];
    ret[1] = m[2] * v.m[0] - m[0] * v.m[2];
    ret[2] = m[0] * v.m[1] - m[1] * v.m[0];
    return Matrix(ret, 3, 1);
  }

  template <typename TT> Matrix &operator=(const TT &rhs) {
    return *this = (rhs);
  }
  template <typename TT> Matrix &operator+=(const TT &rhs) {
    return *this = (this->operator+(rhs));
  }
  template <typename TT> Matrix &operator-=(const TT &rhs) {
    return *this = (this->operator-(rhs));
  }
  template <typename TT> Matrix &operator*=(const TT &rhs) {
    return *this = (this->operator*(rhs));
  }

  const std::array<T, 9> &mat() const { return m; }

  std::istream &operator<<(std::istream &is) const {
    for (auto &i : m) {
      is >> i;
    }
    return is;
  }

private:
  Matrix(size_t nr_, size_t nc_) : nc(nc_), nr(nr_) {}
  std::array<T, 9> m;
  size_t nc, nr;
};

template <typename T> void to_json(json &j, const Matrix<T> &m) {
  j = std::remove_reference_t<decltype(j)>::object();
  j["nc"] = m.nc;
  j["nr"] = m.nr;
  j["m"] = m.m;
}
template <typename T> void from_json(const json &j, Matrix<T> &m) {
  m.m = j["m"].get<std::array<T, 9>>();
  m.nc = j["nc"];
  m.nr = j["nr"];
}

template <typename T>
inline std::ostream &operator<<(std::ostream &os, const Matrix<T> &mat) {
  return mat >> os;
}

template <typename T>
inline Matrix<T> operator*(const T &scalar, const Matrix<T> &mat) {
  return mat * scalar;
}

template <typename T>
inline Matrix<T> operator+(const T &scalar, const Matrix<T> &mat) {
  return mat + scalar;
}

template <typename T>
inline Matrix<T> operator-(const T &scalar, const Matrix<T> &mat) {
  return scalar + mat * (T)(-1.0);
}

} // namespace rb_sim

#endif
