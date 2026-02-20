import otmeshing
import openturns as ot
import openturns.viewer as otv

mesher = otmeshing.UnionMesher()
dim = 2
mesh1 = ot.IntervalMesher([1] * dim).build(ot.Interval([0.0] * dim, [1.0] * dim))
mesh2 = ot.IntervalMesher([1] * dim).build(ot.Interval([3.0] * dim, [4.0] * dim))
union = mesher.build([mesh1, mesh2])

graph = union.draw()
graph.setTitle("Polygon union")
view = otv.View(graph)
otv.View.ShowAll()
