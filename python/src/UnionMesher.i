// SWIG file UnionMesher.i

%{
#include "otmeshing/UnionMesher.hxx"
%}

%apply const MeshCollection & { const OTMESHING::UnionMesher::MeshCollection & };

%include UnionMesher_doc.i

%include otmeshing/UnionMesher.hxx

%copyctor OTMESHING::UnionMesher;
