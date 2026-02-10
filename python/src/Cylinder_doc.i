%feature("docstring") OTMESHING::Cylinder
"Generalized cylinder.

Parameters
----------
base : :py:class:`openturns.Mesh`
    Cylinder base
extension : :py:class:`openturns.Interval`
    Extension range 
injection : sequence of int
    Dimension indices of the extension
discretization : int
    Discretization number along dimensions of the extension

Examples
--------
Create a cylinder with an octogonal base, extended in the z-axis over [0, 5]:

>>> import otmeshing
>>> import openturns as ot
>>> import math
>>> polyline = []
>>> n = 8
>>> for i in range(n):
...     theta = i * 2.0 * math.pi / n
...     polyline.append([math.cos(theta), math.sin(theta)])
>>> base = otmeshing.PolygonMesher().build(polyline)
>>> extension = ot.Interval([0.0], [5.0])  # z-range
>>> injection = [2]  # extension in the z-axis
>>> M = 4  # number of cells in the z-axis
>>> cylinder = otmeshing.Cylinder(base, extension, injection, M)"

// ---------------------------------------------------------------------

%feature("docstring") OTMESHING::Cylinder::getVertices
"Vertices accessor.

Returns
-------
vertices : :py:class:`openturns.Sample`
    Cylinder vertices."

// ---------------------------------------------------------------------

%feature("docstring") OTMESHING::Cylinder::getBoundingBox
"Bounding box accessor.

Returns
-------
bbox : :py:class:`openturns.Interval`
    Cylinder bounding box."
    
// ---------------------------------------------------------------------

%feature("docstring") OTMESHING::Cylinder::getVolume
"Volume accessor.

Returns
-------
volume : float
    Cylinder volume."