"""
Convex hull meshing
===================
"""

# %%
import openturns as ot
import openturns.viewer as otv
import otmeshing

# %%
# Generate a sample of Gaussian points
points = ot.Normal(2).getSample(1000)
graph = ot.Graph("sample", "X1", "X2", True, '')
cloud = ot.Cloud(points)
graph.add(cloud)
view = otv.View(graph)

# %%
# Compute the convex hull
mesher = otmeshing.ConvexHullMesher()
hull = mesher.build(points)

# %%
# Plot convex hull
graph = hull.draw()
graph.add(cloud)
graph.setLegendPosition("upper left")
graph.setTitle("Convex hull of sample")
view = otv.View(graph)

# %%
otv.View.ShowAll()
