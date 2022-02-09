def rectangular_topology(x, y):
    """Draw the vertices position and dictionary for a rectangular base besh

    Parameters
    ----------
    x : list
        List with the x-coordinates required for the mesh
    y : list
        List with the y-coordinates required for the mesh

    Returns
    -------
    xi: list
        List with the x-coordinate of the vertices.
    yi: list
        List with the y-coordinate of the vertices.
    faces_i: list of lists
        Lists with the connectivity required to create the mesh.
    """

    index = 0
    uv_i = {}
    faces = []
    faces_i = []

    xi = []
    yi = []

    for i in range(len(x)):
        for j in range(len(y)):
            uv_i[(i, j)] = index
            xi.append(x[i])
            yi.append(y[j])

            if i < len(x) - 1 and j < len(y) - 1:
                p1 = (i, j)
                p2 = (i, j+1)
                p3 = (i+1, j)
                p4 = (i+1, j+1)
                if i != j and i + j != len(x) - 2:
                    face = [p1, p2, p4, p3]
                    faces.append(face)
                else:
                    if i == j:
                        faces.append([p1, p2, p4])
                        faces.append([p1, p4, p3])
                    else:
                        faces.append([p2, p3, p1])
                        faces.append([p2, p4, p3])
            index = index + 1

    for face in faces:
        face_i = []
        for uv in face:
            u, v = uv
            i = uv_i[(u, v)]
            face_i.append(i)
        faces_i.append(face_i)

    return xi, yi, faces_i
