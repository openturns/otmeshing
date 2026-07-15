// SWIG file VolumeMesher.i

%{
#include "otmeshing/VolumeMesher.hxx"
%}

%include VolumeMesher_doc.i

%copyctor OTMESHING::VolumeMesher;

%include otmeshing/VolumeMesher.hxx
