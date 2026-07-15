%feature("docstring") OTMESHING::VolumeMesher
"Build a volume mesh from a surface mesh.

Tetrahedralizes a closed triangular surface mesh by creating a fan of
tetrahedra from an apex vertex to each boundary triangle.
"

// ---------------------------------------------------------------------

%feature("docstring") OTMESHING::VolumeMesher::build
"Build a volume mesh from a surface mesh.

Parameters
----------
surface : :class:`~openturns.Mesh`
    A surface mesh (triangles in dimension 3).

Returns
-------
volume : :class:`~openturns.Mesh`
    A volume mesh (tetrahedra in dimension 3).
"

// ---------------------------------------------------------------------

%feature("docstring") OTMESHING::VolumeMesher::setApexStrategy
"Apex strategy accessor.

Parameters
----------
strategy : int 
    VolumeMesher.CENTROID or VolumeMesher.FIRSTVERTEX.
"

// ---------------------------------------------------------------------

%feature("docstring") OTMESHING::VolumeMesher::getApexStrategy
"Apex strategy accessor.

Returns
-------
strategy : int
    VolumeMesher.CENTROID or VolumeMesher.FIRSTVERTEX.
"
