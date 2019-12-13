import math

__author__    = ['Ricardo Maia Avelino <mricardo@ethz.ch>']
__copyright__ = 'Copyright 2019, BLOCK Research Group - ETH Zurich'
__license__   = 'MIT License'
__email__     = 'mricardo@ethz.ch'


__all__ = [
    'paraboloid',
    'dome'
]

def paraboloid(x, y, coef= -0.1):

    return coef * (x*x + y*y)

def dome(x, y, r):

    return math.sqrt(r*r - x*x - y*y)