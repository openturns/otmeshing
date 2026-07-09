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

%feature("docstring") OTMESHING::IntersectionMesher::buildWithConvexParts
"Generate the mesh of the intersection of a mesh with pre-decomposed convex pieces.

Parameters
----------
mesh : :py:class:`openturns.Mesh`
    The first mesh, to be internally decomposed into convex parts.
convexPieces : sequence of :py:class:`openturns.Sample`
    Pre-decomposed convex vertex sets for the second operand, e.g. from :meth:`buildCylinderConvex`.

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

%feature("docstring") OTMESHING::IntersectionMesher::buildConvexSample
"Generate the vertices of the intersection of convexes.

Parameters
----------
coll : sequence of :py:class:`openturns.Sample`
    Input convex vertices.

Returns
-------
vertices : :py:class:`openturns.Sample`
    The vertices of the intersection."

// ---------------------------------------------------------------------

%feature("docstring") OTMESHING::IntersectionMesher::buildCylinder
"Generate the mesh of the intersection of cylinders.

Parameters
----------
coll : sequence of :class:`~otmeshing.Cylinder`
    Input cylinders.

Returns
-------
mesh : :py:class:`openturns.Mesh`
    The mesh of the intersection."

// ---------------------------------------------------------------------

%feature("docstring") OTMESHING::IntersectionMesher::buildCylinderConvex
"Generate the convex decomposition of the intersection of cylinders.

Parameters
----------
coll : sequence of :class:`~otmeshing.Cylinder`
    Input cylinders.

Returns
-------
pieces : sequence of :py:class:`openturns.Sample`
    Convex vertex sets representing the intersection, suitable for use with :meth:`build`."


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

// ---------------------------------------------------------------------

%feature("docstring") OTMESHING::IntersectionMesher::setUseSimplicesDecomposition
"Simplicial decomposition flag accessor.

Parameters
----------
useSimplicesDecomposition : bool
    Whether to decompose the mesh by its simplices.
"

// ---------------------------------------------------------------------

%feature("docstring") OTMESHING::IntersectionMesher::getUseSimplicesDecomposition
"Simplicial decomposition flag accessor.

Returns
-------
useSimplicesDecomposition : bool
    Whether to decompose the mesh by its simplices.
"
