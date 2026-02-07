import otmeshing
import openturns as ot
import openturns.viewer as otv

mesher = otmeshing.ConvexIntersectionMesher()
dim = 2
mesh1 = ot.IntervalMesher([1] * dim).build(ot.Interval([0.0] * dim, [3.0] * dim))
mesh2 = ot.IntervalMesher([1] * dim).build(ot.Interval([1.0] * dim, [4.0] * dim))
try:
    intersection = mesher.build(mesh1, mesh2)
except Exception:
    # if cddlib is not available
    intersection = mesh1

graph = intersection.draw()
graph.setTitle("Polygon intersection")
view = otv.View(graph)
otv.View.ShowAll()
