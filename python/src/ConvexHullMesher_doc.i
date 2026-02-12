%feature("docstring") OTMESHING::ConvexHullMesher
"Meshing of the convex hull of a set of points.

Yields a surface mesh of the convex hull (of intrinsic dimension d-1).

Examples
--------
>>> import openturns as ot
>>> import otmeshing
>>> points = ot.Normal(2).getSample(100)
>>> mesher = otmeshing.ConvexHullMesher()
>>> hull = mesher.build(points)"

// ---------------------------------------------------------------------

%feature("docstring") OTMESHING::ConvexHullMesher::build
"Buld the convex hull a set of points.

Parameters
----------
points : :class:`~openturns.Sample`
    A set of points.

Returns
-------
mesh : :class:`~openturns.Mesh`
    The convex hull."
