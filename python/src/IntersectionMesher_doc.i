%feature("docstring") OTMESHING::IntersectionMesher
"Intersection meshing algorithm.

Examples
--------
Triangulate a parallelogram:

>>> import otmeshing
>>> import openturns as ot
>>> mesher = otmeshing.IntersectionMesher()
>>> dim = 2
>>> mesh1 = ot.IntervalMesher([1] * dim).build(ot.Interval([0.0] * dim, [3.0] * dim))
>>> mesh2 = ot.IntervalMesher([1] * dim).build(ot.Interval([1.0] * dim, [4.0] * dim))
>>> intersection = mesher.build([mesh1, mesh2])  # doctest: +SKIP
"

// ---------------------------------------------------------------------

%feature("docstring") OTMESHING::IntersectionMesher::build
"Generate the mesh of the intersection.

Parameters
----------
coll : sequence of :py:class:`openturns.Mesh`
    Input meshes.

Returns
-------
mesh : :py:class:`openturns.Mesh`
    The mesh of the intersection."

// ---------------------------------------------------------------------

%feature("docstring") OTMESHING::IntersectionMesher::buildConvex
"Generate the mesh of the intersection of convexes.

Parameters
----------
coll : sequence of :py:class:`openturns.Mesh`
    Input convex meshes.

Returns
-------
mesh : :py:class:`openturns.Mesh`
    The mesh of the intersection."

// ---------------------------------------------------------------------

%feature("docstring") OTMESHING::IntersectionMesher::setRecompress
"Recompression flag accessor.

Parameters
----------
recompress : bool
    Whether to eliminate duplicate vertices.
"

// ---------------------------------------------------------------------

%feature("docstring") OTMESHING::IntersectionMesher::getRecompress
"Recompression flag accessor.

Returns
-------
recompress : bool
    Whether to eliminate duplicate vertices.
"
