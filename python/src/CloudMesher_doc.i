%feature("docstring") OTMESHING::CloudMesher
"Mesher from a set of points.

Triangulates a set of points.

This can also be seen as the volumetric mesh (full cells) of the convex hull.

Parameters
----------
triangulationMethod : int
    Triangulation method to use, either:

    - CloudMesher.BASIC (default)
    - CloudMesher.DELAUNAY triangulation with the empty ball property (slower)

Examples
--------
Triangulate a set of points:

>>> import openturns as ot
>>> import otmeshing
>>> points = ot.JointDistribution([ot.Uniform()] * 2).getSample(15)
>>> mesher = otmeshing.CloudMesher(otmeshing.CloudMesher.DELAUNAY)
>>> triangulation = mesher.build(points)"

// ---------------------------------------------------------------------

%feature("docstring") OTMESHING::CloudMesher::build
"Triangulate a set of points.

Parameters
----------
points : :class:`~openturns.Sample`
    A set of points.

Returns
-------
mesh : :class:`~openturns.Mesh`
    The mesh built."
