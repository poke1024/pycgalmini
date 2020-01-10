import pycgalmini
import numpy

polygon = pycgalmini.Polygon2D(numpy.array([
	[-1, -1],
	[0, -12],
	[1, -1],
	[12, 0],
	[1, 1],
	[0, 12],
	[-1, 1],
	[-12, 0]
], dtype=numpy.float64))

skeleton = polygon.create_interior_straight_skeleton()

for a, b in skeleton.bisectors:
	print(a, "<->", b)
