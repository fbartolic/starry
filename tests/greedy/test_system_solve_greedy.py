# -*- coding: utf-8 -*-
"""
System linear solve tests.

"""
import starry
import numpy as np
from scipy.linalg import cho_solve
from scipy.stats import multivariate_normal


def test_solve():
    # Instantiate a star with a dipole map
    A = starry.Primary(starry.Map(ydeg=1), prot=0.0)
    y_true = [0.1, 0.2, 0.3]
    A.map[1, :] = y_true

    # Instantiate two transiting planets with different longitudes of
    # ascending node. This ensures there's no null space!
    b = starry.Secondary(starry.Map(), porb=1.0, r=0.1, t0=-0.05, Omega=30.0)
    c = starry.Secondary(starry.Map(), porb=1.0, r=0.1, t0=0.05, Omega=-30.0)
    sys = starry.System(A, b, c)

    # Generate a synthetic light curve with just a little noise
    t = np.linspace(-0.1, 0.1, 100)
    flux = sys.flux(t)
    sigma = 1e-5
    np.random.seed(0)
    flux += np.random.randn(len(t)) * sigma

    # Place a generous prior on the map coefficients
    A.map.set_prior(L=1)
    sys.set_data(flux, C=sigma ** 2)

    # Solve the linear problem
    mu, cho_cov = sys.solve(t=t)

    # Ensure the likelihood of the true value is close to that of
    # the MAP solution
    mean = mu[0]
    cov = cho_cov[0].dot(cho_cov[0].T)
    LnL0 = multivariate_normal.logpdf(mean, mean=mean, cov=cov)
    LnL = multivariate_normal.logpdf([0.1, 0.2, 0.3], mean=mean, cov=cov)
    assert LnL0 - LnL < 5.00


def test_lnlike():
    # Instantiate a star with a dipole map
    A = starry.Primary(starry.Map(ydeg=1), prot=0.0)
    y_true = [0.1, 0.2, 0.3]
    A.map[1, :] = y_true

    # Instantiate a transiting planet
    b = starry.Secondary(starry.Map(), porb=1.0, r=0.1, Omega=30.0)
    sys = starry.System(A, b)

    # Generate a synthetic light curve with just a little noise
    t = np.linspace(-0.1, 0.1, 100)
    flux = sys.flux(t)
    sigma = 1e-5
    np.random.seed(0)
    flux += np.random.randn(len(t)) * sigma

    # Place a generous prior on the map coefficients
    A.map.set_prior(L=1)

    # Provide the dataset
    sys.set_data(flux, C=sigma ** 2)

    # Compute the marginal log likelihood for different secondari radii
    rs = [0.05, 0.075, 0.1, 0.125, 0.15]
    ll = np.zeros_like(rs)
    for i, r in enumerate(rs):
        b.r = r
        ll[i] = sys.lnlike(t=t)

    # Verify that we get the correct radius
    assert rs[np.argmax(ll)] == 0.1