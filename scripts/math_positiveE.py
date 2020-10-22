from compas_tno.diagrams import FormDiagram
from compas_tno.plotters import plot_form
from compas_tno.problems import initialise_problem
from numpy.linalg import matrix_rank
from compas.numerical import equilibrium_matrix
from compas.numerical import connectivity_matrix

discr = 2
span = 1.0
sag = span/10

lines_A = []
lines_B = []
lines_C = []

pts = []
pts_sag = []

for j in range(discr+1):
    for i in range(discr+1):
        xi = float(span * i / discr)
        yi = float(span * j / discr)
        pts.append([xi, yi, 0.0])


lines_A = [[pts[0], pts[1]],
            [pts[1], pts[2]],
            [pts[0], pts[3]],
            [pts[1], pts[4]],
            [pts[2], pts[5]],
            [pts[3], pts[4]],
            [pts[4], pts[5]],
            [pts[3], pts[6]],
            [pts[4], pts[7]],
            [pts[5], pts[8]],
            [pts[6], pts[7]],
            [pts[7], pts[8]],
            ]
lines_A_i = [
    [0,1],
    [1,2],
    [0,3],
    [1,4],
    [2,5],
    [3,4],
    [4,5],
    [3,6],
    [4,7],
    [5,8],
    [6,7],
    [7,8],
]

formA = FormDiagram.from_lines(lines_A, 'list')
fixed = [0, 2, 6, 8]
for key in fixed:
    formA.vertex_attribute(key, 'is_fixed', True)
# plot_form(formA).show()

args = initialise_problem(formA)
ind = args[1]
EA = args[3]
print('\n------- EA -------')
print(EA.shape)
print(matrix_rank(EA))
print(ind)
print(EA)

print(lines_A)

C = connectivity_matrix(lines_A_i)
free = list(set([key for key in formA.vertices()]) - set(fixed))
xyz = formA.vertices_attributes('xyz')
E = equilibrium_matrix(C, xyz, free)

print(E)
print(e)

lines_B = lines_A[:]
lines_B.append([pts[0], pts[4]])
lines_B.append([pts[2], pts[4]])
lines_B.append([pts[6], pts[4]])
lines_B.append([pts[8], pts[4]])

formB = FormDiagram.from_lines(lines_B)
for key in fixed:
    formB.vertex_attribute(key, 'is_fixed', True)
# plot_form(formB).show()

args = initialise_problem(formB)
ind = args[1]
EB = args[3]
print('\n------- EB -------')
print(EB.shape)
print(matrix_rank(EB))
print(ind)
print(EB)


pts[1][1] += sag
pts[3][0] += sag
pts[5][0] -= sag
pts[7][1] -= sag

formC = FormDiagram.from_lines(lines_A)
for key in fixed:
    formC.vertex_attribute(key, 'is_fixed', True)
# plot_form(formC).show()

args = initialise_problem(formC)
ind = args[1]
EC = args[3]
print('\n------- EC -------')
print(EC.shape)
print(matrix_rank(EC))
print(ind)
print(EC)


lines_D = [
    lines_A[0],
    lines_A[1],
    lines_A[2],
    lines_A[4],
    lines_A[7],
    lines_A[9],
    lines_A[10],
    lines_A[11],
]

print(lines_D)

formD = FormDiagram.from_lines(lines_D)
fixed = [0, 2, 5, 6]
for key in fixed:
    formD.vertex_attribute(key, 'is_fixed', True)

from compas_plotters import MeshPlotter
plotter = MeshPlotter(formD)
plotter.draw_vertices(text={key: key for key in formD.vertices()}, facecolor={key: '#ff0000' for key in fixed})
plotter.draw_faces()
plotter.draw_edges()
plotter.show()

# plot_form(formD).show()

args = initialise_problem(formD)
ind = args[1]
ED = args[3]
print('\n------- ED -------')
print(ED.shape)
print(matrix_rank(ED))
print(ind)
print(EC)



# # print:

# #------- EA -------
# #(10, 12)
# EA=[[ 0.   0.5 -0.5  0.   0.   0.   0.   0.   0.   0.   0.   0. ]   # 1 x
#     [ 0.   0.   0.   0.   0.   0.  -0.5  0.   0.   0.   0.   0. ]   # 3 x
#     [ 0.   0.   0.   0.   0.   0.   0.5 -0.5  0.   0.   0.   0. ]   # 4 x
#     [ 0.   0.   0.   0.   0.   0.   0.   0.5  0.   0.   0.   0. ]   # 5 x
#     [ 0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.5 -0.5]   # 7 x
#     [ 0.   0.   0.  -0.5  0.   0.   0.   0.   0.   0.   0.   0. ]   # 1 y
#     [ 0.5  0.   0.   0.   0.  -0.5  0.   0.   0.   0.   0.   0. ]   # 3 y
#     [ 0.   0.   0.   0.5  0.   0.   0.   0.  -0.5  0.   0.   0. ]   # 4 y
#     [ 0.   0.   0.   0.   0.5  0.   0.   0.   0.  -0.5  0.   0. ]   # 5 y
#     [ 0.   0.   0.   0.   0.   0.   0.   0.   0.5  0.   0.   0. ]]  # 7 y

# [ 1 _ _ _ _ _ 1 _ _ _]
# # 1 3 4 5 7 1 3 4 5 7
# # x x x x x y y y y y


# #------- EB -------
# #(10, 16)
# EB=[[ 0.   0.5  0.  -0.5  0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0. ]   # 1
#     [ 0.   0.   0.   0.   0.   0.   0.   0.  -0.5  0.   0.   0.   0.   0.   0.   0. ]   # 3
#     [ 0.   0.   0.5  0.   0.   0.  -0.5  0.   0.5 -0.5  0.5 -0.5  0.   0.   0.   0. ]   # 4
#     [ 0.   0.   0.   0.   0.   0.   0.   0.   0.   0.5  0.   0.   0.   0.   0.   0. ]   # 5
#     [ 0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.5 -0.5]   # 7
#     [ 0.   0.   0.   0.  -0.5  0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0. ]   # 1
#     [ 0.5  0.   0.   0.   0.   0.   0.  -0.5  0.   0.   0.   0.   0.   0.   0.   0. ]   # 3
#     [ 0.   0.   0.5  0.   0.5  0.   0.5  0.   0.   0.  -0.5 -0.5 -0.5  0.   0.   0. ]   # 4
#     [ 0.   0.   0.   0.   0.   0.5  0.   0.   0.   0.   0.   0.   0.  -0.5   0.   0. ]  # 5
#     [ 0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.5  0.   0.   0. ]]  # 7


# # ------- EC -------
# #(10, 12)
# EC=[[ 0.   0.5 -0.5  0.   0.   0.   0.   0.   0.   0.   0.   0. ]   # 1
#     [ 0.1  0.   0.   0.   0.   0.1 -0.4  0.   0.   0.   0.   0. ]   # 3
#     [ 0.   0.   0.   0.   0.   0.   0.4 -0.4  0.   0.   0.   0. ]   # 4
#     [ 0.   0.   0.   0.  -0.1  0.   0.   0.4  0.  -0.1  0.   0. ]   # 5
#     [ 0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.5 -0.5]   # 7
#     [ 0.   0.1  0.1 -0.4  0.   0.   0.   0.   0.   0.   0.   0. ]   # 1
#     [ 0.5  0.   0.   0.   0.  -0.5  0.   0.   0.   0.   0.   0. ]   # 3
#     [ 0.   0.   0.   0.4  0.   0.   0.   0.  -0.4  0.   0.   0. ]   # 4
#     [ 0.   0.   0.   0.   0.5  0.   0.   0.   0.  -0.5  0.   0. ]   # 5
#     [ 0.   0.   0.   0.   0.   0.   0.   0.   0.4  0.  -0.1 -0.1]]  # 7


# # ------- ED -------
# #(10, 12)
# EC=[[ 0.   0.5 -0.5  0.   0.   0.   0.   0.   0.   0.   0.   0. ]   # 1
#     [ 0.1  0.   0.   0.   0.   0.1 -0.4  0.   0.   0.   0.   0. ]   # 3
#     [ 0.   0.   0.   0.   0.   0.   0.4 -0.4  0.   0.   0.   0. ]   # 4
#     [ 0.   0.   0.   0.  -0.1  0.   0.   0.4  0.  -0.1  0.   0. ]   # 5
#     [ 0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.5 -0.5]   # 7
#     [ 0.   0.1  0.1 -0.4  0.   0.   0.   0.   0.   0.   0.   0. ]   # 1
#     [ 0.5  0.   0.   0.   0.  -0.5  0.   0.   0.   0.   0.   0. ]   # 3
#     [ 0.   0.   0.   0.4  0.   0.   0.   0.  -0.4  0.   0.   0. ]   # 4
#     [ 0.   0.   0.   0.   0.5  0.   0.   0.   0.  -0.5  0.   0. ]   # 5
#     [ 0.   0.   0.   0.   0.   0.   0.   0.   0.4  0.  -0.1 -0.1]]  # 7
