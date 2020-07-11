/**
\file latitude.h
\brief Functions to compute moments of the longitude distribution.

*/

#ifndef _STARRY_GP_LONGITUDE_H_
#define _STARRY_GP_LONGITUDE_H_

#include "../utils.h"
#include "wigner.h"

namespace starry {
namespace gp {

using namespace wigner;

namespace integrals {

using namespace utils;

/**
  Spot longitude integrals.

*/
template <class Scalar> class Longitude {
protected:
  // Sizes
  const int ydeg; /**< Degree of the spherical harmonic vector */
  const int Ny;   /**< Number of spherical harmonic `(l, m)` coefficients */

  Wigner<Scalar, 2> W;

  Matrix<Scalar, RowMajor> term;

  Matrix<Scalar, RowMajor> q;
  Matrix<Scalar, RowMajor> Q;

  Matrix<Scalar, RowMajor> matrix;
  RowVector<Scalar> vector;

  inline void precompute() {
    term.resize(4 * ydeg + 1, 4 * ydeg + 1);
    term.setZero();
    term(0, 0) = 2 * pi<Scalar>();
    term(1, 0) = 4.0;
    for (int i = 0; i < 4 * ydeg + 1; ++i) {
      if (i >= 2) {
        term(i, 0) = (i - 1.0) / i * term(i - 2, 0);
      }
      for (int j = 2; j < 4 * ydeg + 1; j += 2) {
        term(i, j) = (j - 1.0) / (i + j) * term(i, j - 2);
      }
    }
  }

public:
  /*

  */
  inline void set_vector(const RowVector<Scalar> &vector_) {
    vector = vector_;
    q = W.mom1(vector);
  }

  /*

  */
  inline void set_matrix(const Matrix<Scalar, RowMajor> &matrix_) {
    matrix = matrix_;
    Q = W.mom2(matrix);
  }

  /*
    Return the first moment vector.

  */
  inline RowVector<Scalar> mom1() {
    RowVector<Scalar> m(Ny);
    for (int l = 0; l < ydeg + 1; ++l) {
      auto t =
          term.block(0, 0, 2 * l + 1, 2 * l + 1).rowwise().reverse().diagonal();
      int i0 = l * l;
      int ilen = 2 * l + 1;
      m.segment(i0, ilen) = q.block(i0, 0, ilen, 2 * l + 1) * t;
    }
    return m / (2 * pi<Scalar>());
  }

  Longitude(int ydeg) : ydeg(ydeg), Ny((ydeg + 1) * (ydeg + 1)), W(ydeg) {
    precompute();
  }
};

} // namespace integrals
} // namespace gp
} // namespace starry

#endif