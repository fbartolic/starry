import numpy as np
import theano
import theano.tensor as tt
import theano.sparse as ts
import pytest

tt.config.floatX = "float64"
sizes = [10, 100, 1000, 10000, 100000]


@pytest.mark.parametrize("M", sizes)
def test_dot_product(M, N=300, L=1):
    def dot_product():
        U = tt.ones((M, N))
        V = tt.ones((N, L))
        dot = tt.dot(U, V)
        return tt.sum(dot) / (M * N * L)

    func = theano.function([], dot_product())
    print(func())
