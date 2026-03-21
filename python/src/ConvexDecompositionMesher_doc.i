%feature("docstring") OTMESHING::ConvexDecompositionMesher
"Build a convex decomposition."

// ---------------------------------------------------------------------

%feature("docstring") OTMESHING::ConvexDecompositionMesher::build
"Build a convex decomposition.

Parameters
----------
mesh : :class:`~openturns.Mesh`
    A mesh.

Returns
-------
decomposition : sequence of :class:`~openturns.Mesh`
    A sequence of non-convex polyhedra."

// ---------------------------------------------------------------------

%feature("docstring") OTMESHING::ConvexDecompositionMesher::IsConvex
"Test whether a mesh is convex.

We test the volume of the mesh versus the volume of its convex hull.

Parameters
----------
mesh : :py:class:`openturns.Mesh`
    A mesh.

Returns
-------
isConvex : bool
    Whether the mesh is convex.
"

// ---------------------------------------------------------------------

%feature("docstring") OTMESHING::ConvexDecompositionMesher::setUseSimplicesDecomposition
"Simplicial decomposition flag accessor.

Parameters
----------
useSimplicesDecomposition : bool
    Whether to decompose the mesh by its simplices.
"

// ---------------------------------------------------------------------

%feature("docstring") OTMESHING::ConvexDecompositionMesher::getUseSimplicesDecomposition
"Simplicial decomposition flag accessor.

Returns
-------
useSimplicesDecomposition : bool
    Whether to decompose the mesh by its simplices.
"
