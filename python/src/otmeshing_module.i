// SWIG file otmeshing_module.i

%module(docstring="otmeshing module") otmeshing

%{
#include <openturns/OT.hxx>
#include <openturns/PythonWrappingFunctions.hxx>
%}

// Prerequisites needed
%include typemaps.i
%include exception.i
%ignore *::load(OT::Advocate & adv);
%ignore *::save(OT::Advocate & adv) const;

%import base_module.i
%import uncertainty_module.i

// The new classes
%include otmeshing/otmeshingprivate.hxx
%include CloudMesher.i
%include ConvexHullMesher.i
%include ConvexDecompositionMesher.i
%include Cylinder.i
%include IntersectionMesher.i
%include PolygonMesher.i
