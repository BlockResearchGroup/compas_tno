import math


__all__ = [
    'paraboloid',
    'dome'
]


def paraboloid(x, y, coef=-0.1):
    """Function of a paraboloid

    Parameters
    ----------
    x : float
        x-coordinate
    y : float
        y-coordinate
    coef : float, optional
        Coefficient of the equation, by default -0.1

    Returns
    -------
    z : float
        Value at (x, y)
    """

    return coef * (x*x + y*y)


def dome(x, y, r):
    """Function of a dome

    Parameters
    ----------
    x : float
        x-coordinate
    y : float
        y-coordinate
    r : float
        Radius

    Returns
    -------
    z : float
        Value at (x, y)
    """

    return math.sqrt(r*r - x*x - y*y)
