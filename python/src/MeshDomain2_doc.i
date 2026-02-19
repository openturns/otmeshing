%feature("docstring") OTMESHING::MeshDomain2
"Adaptor to convert a Mesh to a Domain, with signed distance.

Parameters
----------
mesh : :class:`~openturns.Mesh`
    Underlying mesh.

Examples
--------
Compute the distance from a point to a domain defined from
the mesh of the unit square

>>> import otmeshing
>>> import openturns as ot
>>> dim = 2
>>> mesh = ot.IntervalMesher([1] * dim).build(ot.Interval(dim))
>>> domain = otmeshing.MeshDomain2(mesh)
>>> p = [5.0] * dim
>>> distance = domain.computeDistance(p)
"
