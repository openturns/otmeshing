%feature("docstring") OTMESHING::ConvexIntersectionMesher
"Intersection meshing algorithm for convex meshes.

Examples
--------
Mesh the intersection of two meshed intervals:

>>> import otmeshing
>>> import openturns as ot
>>> mesher = otmeshing.ConvexIntersectionMesher()
>>> dim = 2
>>> mesh1 = ot.IntervalMesher([1] * dim).build(ot.Interval([0.0] * dim, [3.0] * dim))
>>> mesh2 = ot.IntervalMesher([1] * dim).build(ot.Interval([1.0] * dim, [4.0] * dim))
>>> intersection = mesher.build(mesh1, mesh2)  # doctest: +SKIP
"

// ---------------------------------------------------------------------

%feature("docstring") OTMESHING::ConvexIntersectionMesher::build
"Generate the mesh of the intersection.

Parameters
----------
mesh1, mesh2 : :py:class:`openturns.Mesh`
   Input meshes, supposed to be convex (not checked).

Returns
-------
mesh : :py:class:`openturns.Mesh`
    The mesh of the intersection."
