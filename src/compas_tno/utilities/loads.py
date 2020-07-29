

__all__ = [
    'not_sym_load',
]


def not_sym_load(form, x0=0, x1=5.0, magnitude=2.0):

    tol = 0.01

    for key in form.vertices():
        x, _, _ = form.vertex_coordinates(key)
        if x > x0 - tol and x < x1 + tol:
            pz0 = form.vertex_attribute(key, 'pz')
            if x > x1 - tol:
                form.vertex_attribute(key, 'pz', value=((magnitude-1)/2 + 1) * pz0)
            else:
                form.vertex_attribute(key, 'pz', value=magnitude * pz0)

    return form
