/**
\file wigner.h
\brief Spherical harmonic rotation utilities.

*/

#ifndef _STARRY_GP_WIGNER_H_
#define _STARRY_GP_WIGNER_H_

#include "basis.h"
#include "utils.h"

namespace starry {
namespace gp {
namespace wigner {

using namespace utils;

/**
Rotation matrix class for the spherical harmonics.

*/
template <class Scalar> class Wigner {
protected:
  using RowMatrixMap = Eigen::Map<Matrix<Scalar, RowMajor>>;

  // Sizes
  const int ydeg; /**< Degree of the spherical harmonic vector */
  const int Ny;   /**< Number of spherical harmonic `(l, m)` coefficients */

  // Matrices
  std::vector<Matrix<Scalar, RowMajor>> D; /**< The complex Wigner matrix */
  std::vector<Matrix<Scalar, RowMajor>> R; /**< The real Wigner matrix */

  std::vector<Matrix<Scalar, RowMajor>> Rzi;
  std::vector<Matrix<Scalar, RowMajor>> Rzj;
  std::vector<Matrix<Scalar, RowMajor>> Rxi;
  std::vector<Matrix<Scalar, RowMajor>> Rxj;

  // Basis terms
  RowVector<Scalar> b001, b010, b100, bXXX, b10201;

  // Misc
  Scalar tol;
  Scalar root_two;

  // -*- Functions -*-

  // Multiplication by tangent (a shift-left operation on the basis).
  inline RowVector<Scalar> shift_left(const RowVector<Scalar> &x) {
    RowVector<Scalar> xs(x.size());
    xs.segment(0, xs.size() - 1) = x.segment(1, x.size() - 1);
    xs(xs.size() - 1) = 0;
    return xs;
  }

  // Product of two basis vectors.
  inline RowVector<Scalar> prod(const RowVector<Scalar> &x1,
                                const RowVector<Scalar> &x2) {

    int l1 = (x1.size() - 1) / 2;
    int l2 = (x2.size() - 1) / 2;
    int l = l1 + l2;
    RowVector<Scalar> x1x2(2 * l + 1);
    x1x2.setZero();
    for (int m = 0; m < x1.size(); ++m) {
      for (int n = 0; n < x2.size(); ++n) {
        x1x2(m + n) += x1(m) * x2(n);
      }
    }
    return x1x2;
  }

  /**
  Compute the Wigner d matrices.

  */
  inline void dlmn(int l, const Scalar &c1, const Scalar &s1, const Scalar &c3,
                   const Scalar &s3) {
    int iinf = 1 - l;
    int isup = -iinf;
    int m, mp;
    int laux, lbux, lauz, lbuz;
    int sign;
    Scalar a, b;
    Scalar auz, aux, cux, fact, cuz;
    Scalar cosmal, sinmal, cosag, sinag, cosagm, sinagm, cosmga, sinmga;
    RowVector<Scalar> d1, d2;
    RowVector<Scalar> vaux;

    // Compute the D[l;m',m) matrix.
    // First row by recurrence (Eq. 19 and 20 in Alvarez Collado et al.)
    D.at(l).row(4 * l * (l + 1)) =
        prod(D.at(l - 1).row(2 * l * (isup + l - 1)), b001);
    D.at(l).row(2 * l * (2 * l + 1)) = prod(
        D.at(l - 1).row((2 * l - 1) * (isup + l - 1) - isup + l - 1), b100);
    for (m = isup; m > iinf - 1; --m) {
      // Multiplication by s/c
      D.at(l).row(2 * l * (2 * l + 1) + m + l) =
          shift_left(-sqrt(Scalar(l + m + 1) / Scalar(l - m)) *
                     D.at(l).row(2 * l * (2 * l + 1) + l + m + 1));
    }

    // The rows of the upper quarter triangle of the D[l;m',m) matrix
    // (Eq. 21 in Alvarez Collado et al.)
    for (mp = l - 1; mp > -1; --mp) {
      laux = l + mp;
      lbux = l - mp;
      aux = 1.0 / (Scalar(l - 1) * sqrt(Scalar(laux * lbux)));
      cux = sqrt(Scalar((laux - 1) * (lbux - 1))) * l;
      for (m = isup; m > iinf - 1; --m) {
        lauz = l + m;
        lbuz = l - m;
        auz = Scalar(1.0) / sqrt(Scalar(lauz * lbuz));
        fact = aux * auz;
        a = l * (l - 1);
        b = -(m * mp) / a;
        bXXX << b - 1, 0, b + 1;
        D.at(l).row((2 * l + 1) * (mp + l) + m + l) =
            prod(fact * (2 * l - 1) * a *
                     D.at(l - 1).row((2 * l - 1) * (mp + l - 1) + m + l - 1),
                 bXXX);
        if ((lbuz != 1) && (lbux != 1)) {
          cuz = sqrt(Scalar((lauz - 1) * (lbuz - 1)));
          D.at(l).row((2 * l + 1) * (mp + l) + m + l) -=
              (fact * cux * cuz) *
              prod(D.at(l - 2).row((2 * l - 3) * (mp + l - 2) + m + l - 2),
                   b10201);
        }
      }
      ++iinf;
      --isup;
    }

    // The remaining elements of the D[l;m',m) matrix are calculated
    // using the corresponding symmetry relations:
    // reflection ---> ((-1)**(m-m')) D[l;m,m') = D[l;m',m), m'<=m
    // inversion ---> ((-1)**(m-m')) D[l;-m',-m) = D[l;m',m)

    // Reflection
    sign = 1;
    iinf = -l;
    isup = l - 1;
    for (m = l; m > 0; --m) {
      for (mp = iinf; mp < isup + 1; ++mp) {
        D.at(l).row((2 * l + 1) * (mp + l) + m + l) =
            sign * D.at(l).row((2 * l + 1) * (m + l) + mp + l);
        sign *= -1;
      }
      ++iinf;
      --isup;
    }

    // Inversion
    iinf = -l;
    isup = iinf;
    for (m = l - 1; m > -(l + 1); --m) {
      sign = -1;
      for (mp = isup; mp > iinf - 1; --mp) {
        D.at(l).row((2 * l + 1) * (mp + l) + m + l) =
            sign * D.at(l).row((2 * l + 1) * (-mp + l) - m + l);
        sign *= -1;
      }
      ++isup;
    }

    // Compute the real rotation matrices R from the complex ones D
    R.at(l).row((2 * l + 1) * l + l) = D.at(l).row((2 * l + 1) * l + l);
    cosmal = c1;
    sinmal = s1;
    sign = -1;
    for (mp = 1; mp < l + 1; ++mp) {
      cosmga = c3;
      sinmga = s3;
      vaux = root_two * D.at(l).row((2 * l + 1) * l + mp + l);
      R.at(l).row((2 * l + 1) * (mp + l) + l) = vaux * cosmal;
      R.at(l).row((2 * l + 1) * (-mp + l) + l) = vaux * sinmal;
      for (m = 1; m < l + 1; ++m) {
        vaux = root_two * D.at(l).row((2 * l + 1) * (m + l) + l);
        R.at(l).row((2 * l + 1) * l + m + l) = vaux * cosmga;
        R.at(l).row((2 * l + 1) * l - m + l) = -vaux * sinmga;
        d1 = D.at(l).row((2 * l + 1) * (-mp + l) - m + l);
        d2 = sign * D.at(l).row((2 * l + 1) * (mp + l) - m + l);
        cosag = cosmal * cosmga - sinmal * sinmga;
        cosagm = cosmal * cosmga + sinmal * sinmga;
        sinag = sinmal * cosmga + cosmal * sinmga;
        sinagm = sinmal * cosmga - cosmal * sinmga;
        R.at(l).row((2 * l + 1) * (mp + l) + m + l) = d1 * cosag + d2 * cosagm;
        R.at(l).row((2 * l + 1) * (mp + l) - m + l) = -d1 * sinag + d2 * sinagm;
        R.at(l).row((2 * l + 1) * (-mp + l) + m + l) = d1 * sinag + d2 * sinagm;
        R.at(l).row((2 * l + 1) * (-mp + l) - m + l) = d1 * cosag - d2 * cosagm;
        aux = cosmga * c3 - sinmga * s3;
        sinmga = sinmga * c3 + cosmga * s3;
        cosmga = aux;
      }
      sign *= -1;
      aux = cosmal * c1 - sinmal * s1;
      sinmal = sinmal * c1 + cosmal * s1;
      cosmal = aux;
    }
  }

  /**
  Compute the Wigner D and R matrices.

  */
  inline void rotar(const Scalar &c1, const Scalar &s1, const Scalar &c3,
                    const Scalar &s3) {
    Scalar cosag, cosamg, sinag, sinamg;

    // Compute the initial matrices D0, R0, D1 and R1
    D.at(0).row(0) << 1.0;
    R.at(0).row(0) << 1.0;
    D.at(1).row(8) << b001;
    D.at(1).row(7) = -root_two * b010;
    D.at(1).row(6) = b100;
    D.at(1).row(5) = -D.at(1).row(7);
    D.at(1).row(4) = D.at(1).row(8) - D.at(1).row(6);
    D.at(1).row(3) = D.at(1).row(7);
    D.at(1).row(2) = D.at(1).row(6);
    D.at(1).row(1) = D.at(1).row(5);
    D.at(1).row(0) = D.at(1).row(8);
    cosag = c1 * c3 - s1 * s3;
    cosamg = c1 * c3 + s1 * s3;
    sinag = s1 * c3 + c1 * s3;
    sinamg = s1 * c3 - c1 * s3;
    R.at(1).row(4) = D.at(1).row(4);
    R.at(1).row(7) = root_two * D.at(1).row(5) * c1;
    R.at(1).row(1) = root_two * D.at(1).row(5) * s1;
    R.at(1).row(5) = root_two * D.at(1).row(7) * c3;
    R.at(1).row(3) = -root_two * D.at(1).row(7) * s3;
    R.at(1).row(8) = D.at(1).row(8) * cosag - D.at(1).row(6) * cosamg;
    R.at(1).row(6) = -D.at(1).row(8) * sinag - D.at(1).row(6) * sinamg;
    R.at(1).row(2) = D.at(1).row(8) * sinag - D.at(1).row(6) * sinamg;
    R.at(1).row(0) = D.at(1).row(8) * cosag + D.at(1).row(6) * cosamg;

    // The remaining matrices are calculated using
    // symmetry and and recurrence relations
    for (int l = 2; l < ydeg + 1; ++l) {
      dlmn(l, c1, s1, c3, s3);
    }

    return;
  }

public:
  /*
    Return the second moment matrix `Rz * S * Rz^T`.

  */
  inline Matrix<Scalar, RowMajor> RzMom2(const Matrix<Scalar> &S) {
    Matrix<Scalar, RowMajor> M(Ny * Ny, 4 * ydeg + 1);
    M.setZero();
    for (int l1 = 0; l1 < ydeg + 1; ++l1) {
      int I = 2 * l1 + 1;
      int l1_2 = l1 * l1;
      for (int l2 = 0; l2 < l1 + 1; ++l2) {
        int J = 2 * l2 + 1;
        int K = 2 * (l1 + l2) + 1;
        int l2_2 = l2 * l2;

        // R . S
        auto Sij = S.block(l1_2, l2_2, I, J);
        Matrix<Scalar, RowMajor> tmp = (Rzi.at(l1) * Sij).transpose();
        RowMatrixMap RziS(tmp.data(), J * I, I);

        // (R . S) . R^T
        for (int ii = 0; ii < I; ++ii) {
          for (int jj = 0; jj < J; ++jj) {
            for (int kk = 0; kk < J; ++kk) {
              int ij = Ny * (l1_2 + ii) + l2_2 + jj;
              M.row(ij).segment(0, K) += prod(RziS.col(ii).segment(I * kk, I),
                                              Rzj[l2].row(J * jj + kk));
            }
          }
        }
      }
    }

    return M;
  }

  Wigner(int ydeg) : ydeg(ydeg), Ny((ydeg + 1) * (ydeg + 1)) {

    // Allocate the Wigner matrices
    D.resize(ydeg + 1);
    R.resize(ydeg + 1);
    for (int l = 0; l < ydeg + 1; ++l) {
      int sz = 2 * l + 1;
      D.at(l).resize(sz * sz, sz);
      R.at(l).resize(sz * sz, sz);
    }

    // Basis terms
    b001.resize(3);
    b010.resize(3);
    b100.resize(3);
    bXXX.resize(3);
    b10201.resize(5);
    b001 << 0, 0, 1;
    b010 << 0, 1, 0;
    b100 << 1, 0, 0;
    b10201 << 1, 0, 2, 0, 1;

    // Misc
    tol = 10 * mach_eps<Scalar>();
    root_two = sqrt(Scalar(2.0));

    // Pre-compute the lat/lon rotation matrices
    rotar(0.0, 1.0, 0.0, -1.0);
    Rzi.resize(ydeg + 1);
    Rzj.resize(ydeg + 1);
    for (int l = 0; l < ydeg + 1; ++l) {

      int I = 2 * l + 1;

      // Left operator
      Matrix<Scalar, RowMajor> tmp = R.at(l).transpose();
      RowMatrixMap Rzil(tmp.data(), I * I, I);
      Rzi.at(l) = Rzil;

      // Right operator
      Rzj.at(l) = R.at(l);
    }
  }
};

} // namespace wigner
} // namespace gp
} // namespace starry

#endif