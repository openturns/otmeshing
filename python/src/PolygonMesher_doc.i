%feature("docstring") OTMESHING::PolygonMesher
"2-d Polygon meshing algorithm.

Examples
--------
Triangulate a simple polygon:

>>> import otmeshing
>>> mesher = otmeshing.PolygonMesher()
>>> polyline = [[3.0, 0.0], [2.0, 0.0], [2.0, 0.75], [2.5, 0.75], [3.0, 0.2]]
>>> triangulation = mesher.build(polyline)"

// ---------------------------------------------------------------------

%feature("docstring") OTMESHING::PolygonMesher::build
"Generate a mesh from polygon coordinates.

Parameters
----------
polyline : :py:class:`openturns.Sample`
    An ordered set of points defining a 2-d polygon (possibly non-convex).

Returns
-------
mesh : :py:class:`openturns.Mesh`
    The triangulation generated."
