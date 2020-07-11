/**
\file spot.h
\brief Functions to compute moments of the spot size & amplitude distribution.

*/

#ifndef _STARRY_GP_SPOT_H_
#define _STARRY_GP_SPOT_H_

#include "../utils.h"

namespace starry {
namespace gp {
namespace integrals {

using namespace utils;

/*
  Shifts elements in a row vector.

*/
template <typename Derived>
RowVector<typename Derived::Scalar> roll(const Eigen::MatrixBase<Derived> &v,
                                         int shift) {
  if (!shift)
    return v;
  RowVector<typename Derived::Scalar> w(v.size());
  if (shift > 0)
    shift = shift % v.size();
  else
    shift = v.size() - (-shift % v.size());
  int rest = v.size() - shift;
  w.head(shift) = v.tail(shift);
  w.tail(rest) = v.head(rest);
  return w;
}

/**
  Spot size and amplitude integrals.

*/
template <class Scalar> class Spot {
protected:
  // Sizes
  const int ydeg; /**< Degree of the spherical harmonic vector */
  const int Ny;   /**< Number of spherical harmonic `(l, m)` coefficients */

  Matrix<Scalar, RowMajor> term;
  Matrix<Scalar, RowMajor> IP;
  Matrix<Scalar, RowMajor> ID;

  Matrix<Scalar, RowMajor> sigma_basis_mat;
  Vector<Scalar> sigma_basis_vec;

  Scalar amp_factor1;
  Scalar amp_factor2;

  Scalar sigma_mu;
  Scalar sigma_sigma;
  Scalar amp_mu;
  Scalar amp_sigma;
  int sign;

  inline void precompute() {

    Scalar A, B, C;

    term.resize(ydeg + 1, ydeg + 1);
    IP.setZero(ydeg + 1, ydeg + 1);
    ID.setZero(ydeg + 1, ydeg + 1);

    /*
    The Ylm coefficient at index n(n + 1) is

        amp * self._term[n].dot(basis)

    where `basis` is `[1, sigma, sigma ** 2, ...]`
    */

    // Seeding values
    IP(0, 0) = 1.0;
    IP(1, 0) = 1.0;
    IP(1, 1) = -sqrt(Scalar(2.0) / pi<Scalar>());
    ID(1, 0) = 1.0;
    term.row(0) = IP.row(0);
    term.row(1) = sqrt(Scalar(3.0)) * IP.row(1);

    // Recurse
    for (int n = 2; n < ydeg + 1; ++n) {
      C = 2.0 * n - 1.0;
      A = C / n;
      B = A - 1;
      IP.row(n) =
          A * roll(ID.row(n - 1), 2) + A * IP.row(n - 1) - B * IP.row(n - 2);
      IP(n, 1) += A * IP(1, 1);
      ID.row(n) = C * IP.row(n - 1) + ID.row(n - 2);
      term.row(n) = sqrt(Scalar(2.0 * n + 1)) * IP.row(n);
    }

    // Reshape stuff
    sigma_basis_vec.resize(ydeg + 1);
    sigma_basis_mat.resize(ydeg + 1, ydeg + 1);
  }

public:
  /*

  */
  inline void set_params(const Scalar &sigma_mu_, const Scalar &sigma_sigma_,
                         const Scalar &amp_mu_, const Scalar &amp_sigma_,
                         const int sign_) {
    // Store
    sigma_mu = sigma_mu_;
    sigma_sigma = sigma_sigma_;
    amp_mu = amp_mu_;
    amp_sigma = amp_sigma_;
    sign = sign_;

    // Sigma Taylor basis
    RowVector<Scalar> tmp(2 * ydeg + 1);
    for (int n = 0; n < 2 * ydeg + 1; ++n)
      tmp(n) = exp(n * (sigma_mu + 0.5 * n * sigma_sigma * sigma_sigma));
    sigma_basis_vec = tmp.head(ydeg + 1);
    for (int n = 0; n < ydeg + 1; ++n) {
      for (int m = 0; m < ydeg + 1; ++m) {
        sigma_basis_mat(n, m) = tmp(n + m);
      }
    }

    // Amplitude integrals
    amp_factor1 = sign * exp(amp_mu + 0.5 * amp_sigma * amp_sigma);
    amp_factor2 = exp(2 * amp_mu + 2 * amp_sigma * amp_sigma);
  }

  /*
    Return the first moment vector.

  */
  inline RowVector<Scalar> mom1() {
    Vector<Scalar> S0;
    RowVector<Scalar> S;
    S0 = amp_factor1 * term * sigma_basis_vec;
    S.setZero(Ny);
    for (int l = 0; l < ydeg + 1; ++l) {
      S[l * (l + 1)] = S0[l];
    }
    S(0) = 1.0;
    return S;
  }

  /*
    Return the second moment matrix.

  */
  inline Matrix<Scalar, RowMajor> mom2() {
    Matrix<Scalar, RowMajor> S(Ny, Ny);
    S.setZero();
    RowVector<Scalar> uS;
    for (int l1 = 0; l1 < ydeg + 1; ++l1) {
      int i = l1 * (l1 + 1);
      uS = term.row(l1) * sigma_basis_mat;
      for (int l2 = 0; l2 < ydeg + 1; ++l2) {
        int j = l2 * (l2 + 1);
        S(i, j) = uS.dot(term.row(l2).transpose());
      }
    }
    S *= amp_factor2;
    S.row(0) = mom1();
    S.col(0) = S.row(0).transpose().eval();
    return S;
  }

  Spot(int ydeg) : ydeg(ydeg), Ny((ydeg + 1) * (ydeg + 1)) { precompute(); }
};

} // namespace integrals
} // namespace gp
} // namespace starry

#endif