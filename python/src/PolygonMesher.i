// SWIG file PolygonMesher.i

%{
#include "otmeshing/PolygonMesher.hxx"
%}

%include PolygonMesher_doc.i

%include otmeshing/PolygonMesher.hxx

%copyctor OTMESHING::PolygonMesher;
