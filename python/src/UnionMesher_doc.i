%feature("docstring") OTMESHING::UnionMesher
"Union of disjoint meshes.

Examples
--------
>>> import openturns as ot
>>> import otmeshing
>>> mesher = otmeshing.UnionMesher()
>>> dim = 2
>>> mesh1 = ot.IntervalMesher([1] * dim).build(ot.Interval(dim))
>>> mesh2 = ot.IntervalMesher([1] * dim).build(ot.Interval([2.0] * dim, [3.0] * dim))
>>> union = mesher.build([mesh1, mesh2])"

// ---------------------------------------------------------------------

%feature("docstring") OTMESHING::UnionMesher::build
"Generate a mesh from the union of disjoint meshes.

Parameters
----------
coll : sequence of :py:class:`openturns.Mesh`
    Non-overlapping meshes.

Returns
-------
mesh : :py:class:`openturns.Mesh`
    The mesh of the union."

// ---------------------------------------------------------------------

%feature("docstring") OTMESHING::UnionMesher::CompressMesh
"Deduplicate mesh vertices.

A k-D tree radius search is used to filter out duplicate vertices.

Parameters
----------
mesh : :py:class:`openturns.Mesh`
    A mesh.

Returns
-------
compressedMesh : :py:class:`openturns.Mesh`
    The crompressed mesh.
"
