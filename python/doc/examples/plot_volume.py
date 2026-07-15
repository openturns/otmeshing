"""
Volume meshing
==============
"""

# %%
import openturns as ot
import otmeshing

# %%
# Generate a tetrahedron point cloud and build its convex hull surface
tetra_pts = ot.Sample(
    [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
)
hull = otmeshing.ConvexHullMesher()
surface = hull.build(tetra_pts)
print("surface: vertices=", surface.getVerticesNumber(),
      "facets=", surface.getSimplicesNumber())

# %%
# Tetrahedralize the surface with the centroid strategy
mesher = otmeshing.VolumeMesher()
mesher.setApexStrategy(otmeshing.VolumeMesher.CENTROID)
volume = mesher.build(surface)
print("volume (centroid):", volume.getVerticesNumber(),
      "vertices,", volume.getSimplicesNumber(),
      "tets, volume=", volume.getVolume())

# %%
# First-vertex strategy produces the same volume but fewer tetrahedra
mesher.setApexStrategy(otmeshing.VolumeMesher.FIRST_VERTEX)
volume_first = mesher.build(surface)
print("volume (first):", volume_first.getVerticesNumber(),
      "vertices,", volume_first.getSimplicesNumber(),
      "tets, volume=", volume_first.getVolume())

# %%
# Cube example
cube_corners = ot.IntervalMesher([1] * 3).build(ot.Interval(3)).getVertices()
cube_surface = hull.build(cube_corners)
mesher.setApexStrategy(otmeshing.VolumeMesher.CENTROID)
cube_volume = mesher.build(cube_surface)
print("cube volume:", cube_volume.getSimplicesNumber(), "tets, volume=", cube_volume.getVolume())
