"""
Polygon meshing
===============
"""

# %%
import copy
import math
import openturns as ot
import openturns.viewer as otv
import otmeshing

# %%
# Define a non-convex regular polygon from a set of points
polyline = []
n = 20
for i in range(n):
    r = 2.0 + i % 2
    theta = i * 2.0 * math.pi / n
    x = r * math.cos(theta)
    y = r * math.sin(theta)
    polyline.append([x, y])

# %%
# Repeat the first point to draw the closed polygon
polyline_closed = copy.copy(polyline)
polyline_closed.append(polyline[0])
graph = ot.Graph("Polygon", "X1", "X2", True, '')
curve = ot.Curve(polyline_closed)
graph.add(curve)
view = otv.View(graph)

# %%
# Triangulate the polygon
mesher = otmeshing.PolygonMesher()
triangulation = mesher.build(polyline)

# %%
# Plot triangulation
graph = triangulation.draw()
graph.setLegendPosition("upper left")
graph.setTitle("Polygon triangulation")
view = otv.View(graph)

# %%
otv.View.ShowAll()
