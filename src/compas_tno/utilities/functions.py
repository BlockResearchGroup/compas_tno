import math


__all__ = [
    'paraboloid',
    'dome'
]


def paraboloid(x, y, coef=-0.1):

    return coef * (x*x + y*y)


def dome(x, y, r):

    return math.sqrt(r*r - x*x - y*y)
