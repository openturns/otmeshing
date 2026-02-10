import math
import otmeshing
import openturns as ot
import openturns.viewer as otv

polyline = []
n = 8
for i in range(n):
    theta = i * 2.0 * math.pi / n
    polyline.append([math.cos(theta), math.sin(theta)])
base = otmeshing.PolygonMesher().build(polyline)
extension = ot.Interval([0.0], [5.0])
injection = [2]
cylinder = otmeshing.Cylinder(base, extension, injection, 4)
mesh = otmeshing.CloudMesher().build(cylinder.getVertices())

shading = True
thetaX, thetaY, thetaZ = 6.1, 2.7, 3.3
drawEdge = True
graph = mesh.draw3D(drawEdge, thetaX, thetaY, thetaZ, shading)
graph.setTitle("Cylinder mesh")
view = otv.View(graph)
otv.View.ShowAll()
