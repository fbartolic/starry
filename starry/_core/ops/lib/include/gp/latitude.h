/**
\file latitude.h
\brief Functions to compute moments of the latitude distribution.

*/

#ifndef _STARRY_GP_LATITUDE_H_
#define _STARRY_GP_LATITUDE_H_

#include "../utils.h"

namespace starry {
namespace gp {
namespace integrals {

using namespace utils;

/**
  Spot latitude integrals.

*/
template <class Scalar> class Latitude {
protected:
  // Sizes
  const int ydeg; /**< Degree of the spherical harmonic vector */
  const int Ny;   /**< Number of spherical harmonic `(l, m)` coefficients */

public:
  Latitude(int ydeg) : ydeg(ydeg), Ny((ydeg + 1) * (ydeg + 1)) {}
};

} // namespace integrals
} // namespace gp
} // namespace starry

#endif