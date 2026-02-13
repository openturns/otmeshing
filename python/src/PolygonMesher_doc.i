%feature("docstring") OTMESHING::PolygonMesher
"2-d Polygon meshing algorithm.

Examples
--------
Triangulate a parallelogram:

>>> import otmeshing
>>> mesher = otmeshing.PolygonMesher()
>>> polyline = [[0, 0], [3, 0], [4, 2], [1, 2]]
>>> triangulation = mesher.build(polyline)"

// ---------------------------------------------------------------------

%feature("docstring") OTMESHING::PolygonMesher::build
"Generate a mesh from polygon coordinates.

Parameters
----------
polyline : :py:class:`openturns.Sample`
    An ordered set of vertices defining a 2-d polygon, possibly non-convex.
    The polygon must be simple (with no redundant vertex).
    The vertices can be of dimension greater than 2,
    but the polygon itsef should be within a single plane.

Returns
-------
mesh : :py:class:`openturns.Mesh`
    The triangulation generated."
