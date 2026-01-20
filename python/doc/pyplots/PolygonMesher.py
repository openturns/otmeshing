import math
import otmeshing
import openturns.viewer as otv

mesher = otmeshing.PolygonMesher()
polyline = []
n = 20
for i in range(n):
    r = 2.0 + i % 2
    theta = i * 2.0 * math.pi / n
    x = r * math.cos(theta)
    y = r * math.sin(theta)
    polyline.append([x, y])
triangulation = mesher.build(polyline)

graph = triangulation.draw()
graph.setTitle("Polygon triangulation")
view = otv.View(graph)
otv.View.ShowAll()
