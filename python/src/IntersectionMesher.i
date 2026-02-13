// SWIG file IntersectionMesher.i

%{
#include "otmeshing/IntersectionMesher.hxx"
%}

%typemap(in) const MeshCollection & {
  if (SWIG_IsOK(SWIG_ConvertPtr($input, (void **) &$1, $1_descriptor, 0))) {
    // From interface class, ok
  } else {
    try {
      $1 = OT::buildCollectionFromPySequence< OT::Mesh >( $input );
    } catch (const OT::InvalidArgumentException &) {
      SWIG_exception(SWIG_TypeError, "Object passed as argument is not convertible to a collection of Mesh");
    }
  }
}

%typemap(typecheck,precedence=SWIG_TYPECHECK_POINTER) const MeshCollection & {
  $1 = SWIG_IsOK(SWIG_ConvertPtr($input, NULL, $1_descriptor, 0))
    || OT::canConvertCollectionObjectFromPySequence< OT::Mesh >( $input );
}

%apply const MeshCollection & { const OTMESHING::IntersectionMesher::MeshCollection & };

%include IntersectionMesher_doc.i

%include otmeshing/IntersectionMesher.hxx

%copyctor OTMESHING::IntersectionMesher;
