/**
\file latitude.h
\brief Functions to compute moments of the latitude distribution.

*/

#ifndef _STARRY_GP_LATITUDE_H_
#define _STARRY_GP_LATITUDE_H_

#include "../utils.h"
#define STARRY_GP_2F1_MAXITER 200
#define STARRY_GP_2F1_MAXTOL 1e-15
#define STARRY_GP_2F1_MINTOL 1e-12

namespace starry {
namespace gp {

using namespace utils;

/**
  The Gauss hypergeometric function 2F1.

*/
template <typename T>
inline T hyp2f1(const T &a_, const T &b_, const T &c_, const T &z) {

  // Compute the value
  T a = a_;
  T b = b_;
  T c = c_;
  T term = a * b * z / c;
  T value = 1.0 + term;
  int n = 1;
  while ((abs(term) > STARRY_GP_2F1_MAXTOL) && (n < STARRY_GP_2F1_MAXITER)) {
    a += 1;
    b += 1;
    c += 1;
    n += 1;
    term *= a * b * z / c / n;
    value += term;
  }
  if ((n == STARRY_GP_2F1_MAXITER) && (abs(term) > STARRY_GP_2F1_MINTOL)) {
    std::cout << abs(term) << std::endl;
    std::stringstream args;
    args << "a_ = " << a_ << ", "
         << "b_ = " << b_ << ", "
         << "c_ = " << c_ << ", "
         << "z = " << z;
    throw StarryException("Series for 2F1 did not converge.", "gp/latitude.h",
                          "hyp2f1", args.str());
  }
  return value;
}

template <class Scalar> class LatitudeGP {

protected:
  int ydeg;
  int N;
  int n;

  Vector<Scalar> B;
  Vector<Scalar> F;
  Matrix<Scalar> term;

public:
  Vector<Scalar> q;
  Matrix<Scalar> Q;

  explicit LatitudeGP(int ydeg)
      : ydeg(ydeg), N((ydeg + 1) * (ydeg + 1)), n(4 * ydeg + 1) {
    B.setZero(n);
    F.setZero(n);
    term.setZero(n, n);
    q.setZero(N);
    Q.setZero(N, N);
  };

  inline void compute(const Scalar &alpha, const Scalar &beta) {

    int n1, n2, j1, i1, j2, i2;

    // B functions
    B(0) = 1.0;
    for (int k = 1; k < n; ++k) {
      B(k) = (alpha - 1.0 + k) / (alpha + beta - 1.0 + k) * B(k - 1);
    }

    // F functions
    Scalar ab = alpha + beta;
    F(0) = sqrt(2.0) * hyp2f1(-0.5, beta, ab, 0.5);
    F(1) = sqrt(2.0) * hyp2f1(-0.5, beta, ab + 1.0, 0.5);
    for (int k = 2; k < n; ++k) {
      F(k) = ((ab + k - 1.0) / ((alpha + k - 1.0) * (ab + k - 0.5)) *
              ((ab + k - 2.0) * F(k - 2) + (1.5 - beta) * F(k - 1)));
    }
    F.array() = F.array().cwiseProduct(B.array()).eval();

    // Terms
    Eigen::Map<Vector<Scalar>> func(NULL, n);
    Scalar fac1, fac2;
    term.setZero();
    for (int i = 0; i < n; ++i) {
      if (is_even(i)) {
        new (&func) Eigen::Map<Vector<Scalar>>(B.data(), n);
        i2 = i / 2;
      } else {
        new (&func) Eigen::Map<Vector<Scalar>>(F.data(), n);
        i2 = (i - 1) / 2;
      }
      for (int j = 0; j < n; j += 2) {
        j2 = j / 2;
        fac1 = 1.0;
        for (int k1 = 0; k1 < i2 + 1; ++k1) {
          fac2 = fac1;
          for (int k2 = 0; k2 < j2 + 1; ++k2) {
            term(i, j) += fac2 * func(k1 + k2);
            fac2 *= (k2 - j2) / (k2 + 1.0);
          }
          fac1 *= (i2 - k1) / (k1 + 1.0);
        }
      }
    }

    // Beta normalization
    term /= B(0);

    // Moment integrals
    n1 = 0;
    Scalar inv_two_l1 = 1.0;
    Scalar inv_two_l1l2;
    for (int l1 = 0; l1 < ydeg + 1; ++l1) {
      for (int m1 = -l1; m1 < l1 + 1; ++m1) {
        j1 = m1 + l1;
        i1 = l1 - m1;
        q(n1) = term(j1, i1) * inv_two_l1;
        n2 = 0;
        inv_two_l1l2 = inv_two_l1;
        for (int l2 = 0; l2 < ydeg + 1; ++l2) {
          for (int m2 = -l2; m2 < l2 + 1; ++m2) {
            j2 = m2 + l2;
            i2 = l2 - m2;
            Q(n1, n2) = term(j1 + j2, i1 + i2) * inv_two_l1l2;
            n2 += 1;
          }
          inv_two_l1l2 *= 0.5;
        }
        n1 += 1;
      }
      inv_two_l1 *= 0.5;
    }
  }

}; // class Latitude

} // namespace gp
} // namespace starry

#endif