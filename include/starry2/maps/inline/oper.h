/**
Generate a random isotropic map with a given power spectrum.

*/
template<typename V, typename U=S, typename=IsSingleColumn<U>>
inline void random (
    const Vector<Scalar>& power,
    const V& seed
) {
    random_(power, seed, 0);
}

/**
Generate a random isotropic map with a given power spectrum.
NOTE: If `col = -1`, sets all columns to the same map.

*/
template<typename V, typename U=S, typename=IsMultiColumn<U>>
inline void random (
    const Vector<Scalar>& power,
    const V& seed,
    int col=-1
) {
    random_(power, seed, col);
}

/**
Compute the Taylor expansion basis at a point in time.
Static specialization (does nothing).

*/
template<typename U=S>
inline IsStatic<U, void> computeTaylor (
    const Scalar & t
) {
}

/**
Compute the Taylor expansion basis at a point in time.
Temporal specialization.

*/
template<typename U=S>
inline IsTemporal<U, void> computeTaylor (
    const Scalar & t
) {
    if (t != cache.taylort) {
        for (int n = 1; n < ncoly; ++n)
            taylor(n) = taylor(n - 1) * t / n;
        cache.taylort = t;
    }
}

/**
Temporal contraction operation for static maps: 
does nothing, and returns a reference to the 
original map.

*/
template<typename U=S, typename T1>
inline IsStatic<U, MatrixBase<T1>&> contract (
    MatrixBase<T1> const & mat
) {
    return MBCAST(mat, T1);
}

/**
Contracts a temporal map by dotting the map matrix with the
Taylor expansion basis.

*/
template<typename U=S>
inline IsTemporal<U, Vector<Scalar>> contract (
    const Matrix<Scalar> & mat
) {
    return mat * taylor;
}

/**
Set the zeroth order limb darkening coefficient.
This is a **constant** whose value ensures that
I(mu = 0) / I0 = 1.

*/
template<typename U=S>
inline IsStatic<U, void> setU0 () {
    u.row(0).setConstant(-1.0);
}

/**
Set the zeroth order limb darkening coefficient
for a temporal map. All derivatives are set to
zero.

*/
template<typename U=S>
inline IsTemporal<U, void> setU0 () {
    u.row(0).setZero();
    u(0, 0) = -1.0;
}

/**
Rotate the map by an angle `theta` and store
the result in `cache.Ry`. Optionally compute
and cache the Wigner rotation matrices and
their derivatives. Static map specialization.

*/
template<typename U=S>
inline IsStatic<U, void> rotateIntoCache (
    const Scalar& theta,
    bool compute_matrices=false
) 
{
    Scalar theta_rad = theta * radian;
    computeWigner();
    if ((!compute_matrices) && (cache.theta != theta)) {
        W.rotate(cos(theta_rad), sin(theta_rad), cache.Ry);
        cache.theta = theta;
    } else if (compute_matrices && (cache.theta_with_grad != theta)) {
        W.compute(cos(theta_rad), sin(theta_rad));
        for (int l = 0; l < lmax + 1; ++l) {
            cache.Ry.block(l * l, 0, 2 * l + 1, ncoly) =
                W.R[l] * y.block(l * l, 0, 2 * l + 1, ncoly);
            cache.DRDthetay.block(l * l, 0, 2 * l + 1, ncoly) =
                W.DRDtheta[l] * y.block(l * l, 0, 2 * l + 1, ncoly);
        }
        cache.theta_with_grad = theta;
    }
}

/**
Rotate the map by an angle `theta` and store
the result in `cache.Ry`. Optionally compute
and cache the Wigner rotation matrices and
their derivatives. Temporal map specialization.

*/
template<typename U=S>
inline IsTemporal<U, void> rotateIntoCache (
    const Scalar& theta,
    bool compute_matrices=false
) 
{    
    Scalar theta_rad = theta * radian;
    computeWigner();
    if (!compute_matrices) {
        if (cache.theta != theta) {
            W.rotate(cos(theta_rad), sin(theta_rad), cache.RyUncontracted);
            cache.theta = theta;
        }
        cache.Ry = contract(cache.RyUncontracted);
    } else { 
        if (cache.theta_with_grad != theta) {
            W.compute(cos(theta_rad), sin(theta_rad));
            for (int l = 0; l < lmax + 1; ++l) {
                cache.RyUncontracted.block(l * l, 0, 2 * l + 1, ncoly) =
                    W.R[l] * y.block(l * l, 0, 2 * l + 1, ncoly);
            }
            cache.theta_with_grad = theta;
        }
        cache.Ry = contract(cache.RyUncontracted);
        for (int l = 0; l < lmax + 1; ++l) {
            cache.DRDthetay.block(l * l, 0, 2 * l + 1, nflx) =
                contract(W.DRDtheta[l] * y.block(l * l, 0, 2 * l + 1, ncoly));
        }
    }
}

/**
Normalize the Agol g basis and its derivatives. 
Default map specialization.

*/
template<typename U=S, typename T1, typename T2>
inline IsDefault<U, void> normalizeAgolGBasis (
    MatrixBase<T1> const & g,
    MatrixBase<T2> const & DgDu
) {
    // The total flux is given by `y00 * (s . g)`
    Scalar norm = Scalar(1.0) / (pi<Scalar>() * (g(0) + 2.0 * g(1) / 3.0));
    MBCAST(g, T1) = g * norm;
    MBCAST(DgDu, T2) = DgDu * norm;
}

/**
Normalize the Agol g basis and its derivatives. 
Spectral map specialization.

*/
template<typename U=S, typename T1, typename T2>
inline IsSpectral<U, void> normalizeAgolGBasis (
    MatrixBase<T1> const & g,
    MatrixBase<T2> const & DgDu
) {
    // The total flux is given by `y00 * (s . g)`
    for (int n = 0; n < ncoly; ++n) {
        Scalar norm = Scalar(1.0) / (pi<Scalar>() * (g(0, n) + 2.0 * g(1, n) / 3.0));
        MBCAST(g, T1).col(n) = g.col(n) * norm;
        MBCAST(DgDu, T2).block(n * lmax, 0, lmax, lmax + 1) = 
            DgDu.block(n * lmax, 0, lmax, lmax + 1) * norm;
    }
}

/**
Normalize the Agol g basis and its derivatives. 
Temporal map specialization.

*/
template<typename U=S, typename T1, typename T2>
inline IsTemporal<U, void> normalizeAgolGBasis (
    MatrixBase<T1> const & g,
    MatrixBase<T2> const & DgDu
) {
    // The total flux is given by `y00 * (s . c)`
    Scalar norm = Scalar(1.0) / (pi<Scalar>() * (g(0) + 2.0 * g(1) / 3.0));
    MBCAST(g, T1) = g * norm;
    MBCAST(DgDu, T2) = DgDu * norm;
}

/** 
Compute the limb darkening polynomial `p` and optionally
propagate gradients. Default specialization.

*/
template<bool GRADIENT=false, typename U=S>
inline IsDefault<U, void> computeLDPolynomial () {
    if (!GRADIENT) {
        if (cache.compute_p) {
            UType tmp = B.U1 * u;
            Scalar norm = pi<Scalar>() * y(0) / B.rT.dot(tmp);
            cache.p = tmp * norm;
            cache.compute_p = false;
        }
    } else {
        if (cache.compute_p_grad) {
            UType tmp = B.U1 * u;
            Scalar rTU1u = B.rT.dot(tmp);
            Scalar norm0 = pi<Scalar>() / rTU1u;
            Scalar norm = norm0 * y(0);
            RowVector<Scalar> dnormdu = -(norm / rTU1u) * B.rTU1;
            cache.p = tmp * norm;
            cache.DpuDu = B.U1 * norm + tmp * dnormdu;
            cache.DpuDy0 = tmp * norm0;
            cache.compute_p_grad = false;
        }
    }
}

/** 
Compute the limb darkening polynomial `p` and optionally
propagate gradients. Spectral specialization.

*/
template<bool GRADIENT=false, typename U=S>
inline IsSpectral<U, void> computeLDPolynomial () {
    if (!GRADIENT) {
        if (cache.compute_p) {
            UType tmp = B.U1 * u;
            UCoeffType norm = pi<Scalar>() * 
                              y.row(0).cwiseQuotient(B.rT * tmp);
            cache.p = tmp.array().rowwise() * norm.array();
            cache.compute_p = false;
        }
    } else {
        if (cache.compute_p_grad) {
            UType tmp = B.U1 * u;
            UCoeffType rTU1u = B.rT * tmp;
            for (int i = 0; i < ncolu; ++i) {
                Scalar norm0 = pi<Scalar>() / rTU1u(i);
                Scalar norm = norm0 * y(0, i);
                RowVector<Scalar> dnormdu = -(norm / rTU1u(i)) * B.rTU1;
                cache.p.col(i) = tmp.col(i) * norm;
                cache.DpuDu[i] = B.U1 * norm + tmp.col(i) * dnormdu;
                cache.DpuDy0.col(i) = tmp.col(i) * norm0;
            }
            cache.compute_p_grad = false;
        }
    }
}

/** 
Compute the limb darkening polynomial `p` and optionally
propagate gradients. Temporal specialization.

*/
template<bool GRADIENT=false, typename U=S>
inline IsTemporal<U, void> computeLDPolynomial () {
    if (!GRADIENT) {
        if (cache.compute_p) {
            UType tmp = B.U1 * u;
            Scalar norm = pi<Scalar>() * contract(y.row(0))(0) / B.rT.dot(tmp);
            cache.p = tmp * norm;
            cache.compute_p = false;
        }
    } else {
        if (cache.compute_p_grad) {
            UType tmp = B.U1 * u;
            Scalar rTU1u = B.rT.dot(tmp);
            Scalar norm0 = pi<Scalar>() / rTU1u;
            Scalar norm = norm0 * contract(y.row(0))(0);
            RowVector<Scalar> dnormdu = -(norm / rTU1u) * B.rTU1;
            cache.p = tmp * norm;
            cache.DpuDu = B.U1 * norm + tmp * dnormdu;
            cache.DpuDy0 = tmp * norm0;
            cache.compute_p_grad = false;
        }
    }
}

/**
Limb-darken a polynomial map, and optionally compute the
gradient of the resulting map with respect to the input
polynomial map and the input limb-darkening map.
Default specialization.

*/
template<bool GRADIENT=false, typename U=S>
inline IsDefault<U, void> limbDarken (
    const CtrYType& poly, 
    CtrYType& poly_ld
) {
    // Compute the limb darkening polynomial
    computeLDPolynomial<GRADIENT>();

    // Compute the gradient
    if (GRADIENT) {
        // Multiply the polynomials
        basis::polymul(y_deg, poly, u_deg, cache.p, lmax, poly_ld, B.rT,
                       cache.rTDpupyDpy, cache.rTDpupyDpu);
        // Propagate the gradient to d(polynomial) / du
        // and d(polynomial) / dy 
        cache.rTDpupyDu = cache.rTDpupyDpu * cache.DpuDu;
        cache.rTDpupyDpyA1R = cache.rTDpupyDpy * B.A1;
        for (int l = 0; l < lmax + 1; ++l)
            cache.rTDpupyDpyA1R.segment(l * l, 2 * l + 1) *= W.R[l];
        cache.rTDpupyDy = cache.rTDpupyDpyA1R.transpose();
        cache.rTDpupyDy(0) += (cache.rTDpupyDpu * cache.DpuDy0)(0);
    } else {
        // Multiply the polynomials
        basis::polymul(y_deg, poly, u_deg, cache.p, lmax, poly_ld);
    }
}

/**
Limb-darken a polynomial map, and optionally compute the
gradient of the resulting map with respect to the input
polynomial map and the input limb-darkening map.
Spectral specialization.

*/
template<bool GRADIENT=false, typename U=S>
inline IsSpectral<U, void> limbDarken (
    const CtrYType& poly, 
    CtrYType& poly_ld
) {
    // Compute the limb darkening polynomial
    computeLDPolynomial<GRADIENT>();

    // Compute the gradient
    if (GRADIENT) {
        // Multiply the polynomials
        basis::polymul(y_deg, poly, u_deg, cache.p, lmax, poly_ld, B.rT,
                       cache.rTDpupyDpy, cache.rTDpupyDpu);
        // Propagate the gradient to d(polynomial) / du
        // and d(polynomial) / dy 
        for (int i = 0; i < ncolu; ++i)
            cache.rTDpupyDu.col(i) = cache.rTDpupyDpu.row(i) * cache.DpuDu[i];
        cache.rTDpupyDpyA1R = cache.rTDpupyDpy * B.A1;
        for (int l = 0; l < lmax + 1; ++l)
            cache.rTDpupyDpyA1R.block(0, l * l, nflx, 2 * l + 1) *= W.R[l];
        for (int i = 0; i < ncoly; ++i) {
            cache.rTDpupyDy.col(i) = cache.rTDpupyDpyA1R.row(i);
            cache.rTDpupyDy(0, i) += (cache.rTDpupyDpu.row(i) * 
                                      cache.DpuDy0.col(i))(0);
        }
    } else {
        // Multiply the polynomials
        basis::polymul(y_deg, poly, u_deg, cache.p, lmax, poly_ld);
    }
}

/**
Limb-darken a polynomial map, and optionally compute the
gradient of the resulting map with respect to the input
polynomial map and the input limb-darkening map.
Temporal specialization.

*/
template<bool GRADIENT=false, typename U=S>
inline IsTemporal<U, void> limbDarken (
    const CtrYType& poly, 
    CtrYType& poly_ld
) {
    // Compute the limb darkening polynomial
    computeLDPolynomial<GRADIENT>();

    // Compute the gradient
    if (GRADIENT) {
        // Multiply the polynomials
        basis::polymul(y_deg, poly, u_deg, cache.p, lmax, poly_ld, B.rT,
                       cache.rTDpupyDpy, cache.rTDpupyDpu);
        // Propagate the gradient to d(polynomial) / du
        // and d(polynomial) / dy 
        cache.rTDpupyDu = cache.rTDpupyDpu * cache.DpuDu;
        cache.rTDpupyDpyA1R = cache.rTDpupyDpy * B.A1;
        for (int l = 0; l < lmax + 1; ++l)
            cache.rTDpupyDpyA1R.segment(l * l, 2 * l + 1) *= W.R[l];
        cache.rTDpupyDy = cache.rTDpupyDpyA1R.replicate(ncoly, 1).transpose();
        cache.rTDpupyDy.row(0) += (cache.rTDpupyDpu * cache.DpuDy0)
                                  .replicate(ncoly, 1).transpose();
    } else {
        // Multiply the polynomials
        basis::polymul(y_deg, poly, u_deg, cache.p, lmax, poly_ld);
    }
}