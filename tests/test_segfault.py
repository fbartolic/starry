import numpy as np
import theano
import theano.tensor as tt
import theano.sparse as ts
import pytest

tt.config.floatX = "float64"
sizes = [10, 100, 1000, 10000, 100000]


def dot_sum(x, y):
    dot = tt.dot(x, y)
    return tt.sum(dot)


x = tt.dmatrix()
y = tt.dmatrix()
func = theano.function([x, y], dot_sum(x, y))


@pytest.mark.parametrize("M", sizes)
def test_dot_product(M, N=300, L=10):
    u = np.random.randn(M, N)
    v = np.random.randn(N, L)
    print(func(u, v))
