import openturns as ot
import openturns.viewer as otv
import otmeshing

tetra_pts = ot.Sample([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
surface = otmeshing.ConvexHullMesher().build(tetra_pts)

mesher = otmeshing.VolumeMesher()
mesher.setApexStrategy(otmeshing.VolumeMesher.CENTROID)
volume = mesher.build(surface)

graph = volume.draw3D(True, 6.1, 3.7, 4.3, True)
V = volume.getVerticesNumber()
T = volume.getSimplicesNumber()
graph.setTitle(f"VolumeMesher centroid V={V} T={T}")

view = otv.View(graph)
otv.View.ShowAll()
