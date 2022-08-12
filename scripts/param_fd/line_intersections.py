from compas_tno.utilities import split_intersection_lines

pt0, pt1, pt2, pt3 = [[0, 0], [0, 1], [1, 1], [1, 0]]

lines = [[pt0, pt2], [pt1, pt3]]
ints = split_intersection_lines(lines)
print(len(ints))
print(ints)

lines = [[pt0, pt1], [pt1, pt2], [pt2, pt3], [pt3, pt0]]
ints = split_intersection_lines(lines)
print(len(ints))
print(ints)

## ----

lines = [[pt0, pt2], [pt1, pt3], [pt0, pt2]]
ints = split_intersection_lines(lines)
print(len(ints))
print(ints)


lines = [[pt0, pt1], [pt0, pt1], [pt1, pt2], [pt2, pt3], [pt3, pt0]]
ints = split_intersection_lines(lines)
print(len(ints))
print(ints)
