# import math
import openturns as ot
import openturns.viewer as otv
import otmeshing as otm

dim = 2
mesh = ot.IntervalMesher([1] * dim).build(ot.Interval(dim))
domain = otm.MeshDomain2(mesh)

# sample points and sort points inside/outside from the distance sign
n = 200
points = ot.Normal(dim).getSample(n)
distances = domain.computeDistance(points)
points_inside = []
points_outside = []
for i in range(n):
    if distances[i, 0] < 0.0:
        points_inside.append(points[i])
    else:
        points_outside.append(points[i])

graph = mesh.draw()
cloud_inside = ot.Cloud(points_inside)
cloud_inside.setColor("orange")
cloud_inside.setLegend("inside")
graph.add(cloud_inside)
cloud_outside = ot.Cloud(points_outside)
cloud_outside.setColor("grey")
cloud_outside.setLegend("outside")
graph.add(cloud_outside)

graph.setTitle("Mesh domain")
view = otv.View(graph)
otv.View.ShowAll()
