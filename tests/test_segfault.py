import numpy as np
import theano
import theano.tensor as tt
import theano.sparse as ts
import pytest
import itertools

floatX = ["float32", "float64"]
sizes = [10, 100, 1000, 10000, 100000, 1000000]


@pytest.mark.parametrize("floatX,M", itertools.product(floatX, sizes))
def test_dot_product(floatX, M, N=300, L=10):

    # The theano function we'll compile
    def dot_sum(x, y):
        dot = tt.dot(x, y)
        return tt.sum(dot)

    # Compile it
    x = tt.matrix(dtype=floatX)
    y = tt.matrix(dtype=floatX)
    func = theano.function([x, y], dot_sum(x, y))

    # Evaluate it
    np.random.seed(0)
    u = np.array(np.random.randn(M, N), dtype=floatX)
    v = np.array(np.random.randn(N, L), dtype=floatX)
    print(func(u, v))
