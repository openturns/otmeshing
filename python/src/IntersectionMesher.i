// SWIG file IntersectionMesher.i

%{
#include "otmeshing/IntersectionMesher.hxx"

namespace OT {
  template <>
  struct traitsPythonType<OTMESHING::Cylinder>
  {
    typedef _PyObject_ Type;
  };

  template <>
  inline
  OTMESHING::Cylinder
  convert< _PyObject_, OTMESHING::Cylinder >(PyObject * pyObj)
  {
    void * ptr = 0;
    if (SWIG_IsOK(SWIG_ConvertPtr(pyObj, &ptr, SWIG_TypeQuery("OTMESHING::Cylinder *"), SWIG_POINTER_NO_NULL))) {
      OTMESHING::Cylinder * p_it = reinterpret_cast< OTMESHING::Cylinder * >(ptr);
      return *p_it;
    }
    else {
      throw OT::InvalidArgumentException(HERE) << "Object passed as argument is not convertible to a Cylinder";
    }
    return OTMESHING::Cylinder();
  }

  // Needed by canConvertCollectionObjectFromPySequence<Sample> which is
  // called from the SampleCollection typecheck typemap (ProcessSample.i
  // is imported but its percent-brace block is not included by SWIG for imports).
  template <>
  inline
  bool canConvert<OT::_PySequence_, OT::Sample>(_object * pyObj)
  {
    return SWIG_IsOK(SWIG_ConvertPtr(pyObj, NULL, SWIG_TypeQuery("OT::Sample *"), SWIG_POINTER_NO_NULL))
             || OT::isAPythonBufferOf<OT::Scalar, 2>(pyObj) || OT::isAPythonSequenceOf<OT::_PySequence_>(pyObj);
  }
}
%}

// MeshCollection typemaps
%typemap(in) const MeshCollection & {
  if (SWIG_IsOK(SWIG_ConvertPtr($input, (void **) &$1, $1_descriptor, 0))) {
    // From interface class, ok
  } else {
    try {
      $1 = OT::buildCollectionFromPySequence< OT::Mesh >($input);
    } catch (const OT::InvalidArgumentException &) {
      SWIG_exception(SWIG_TypeError, "Object passed as argument is not convertible to a collection of Mesh");
    }
  }
}

%typemap(typecheck,precedence=SWIG_TYPECHECK_POINTER) const MeshCollection & {
  $1 = SWIG_IsOK(SWIG_ConvertPtr($input, NULL, $1_descriptor, 0))
    || OT::canConvertCollectionObjectFromPySequence< OT::Mesh >($input);
}

%apply const MeshCollection & { const OTMESHING::IntersectionMesher::MeshCollection & };

// CylinderCollection typemaps
%typemap(in) const CylinderCollection & {
  if (SWIG_IsOK(SWIG_ConvertPtr($input, (void **) &$1, $1_descriptor, 0))) {
    // From interface class, ok
  } else {
    try {
      $1 = OT::buildCollectionFromPySequence< OTMESHING::Cylinder >($input);
    } catch (const OT::InvalidArgumentException &) {
      SWIG_exception(SWIG_TypeError, "Object passed as argument is not convertible to a collection of Cylinder");
    }
  }
}

%typemap(typecheck,precedence=SWIG_TYPECHECK_POINTER) const CylinderCollection & {
  $1 = SWIG_IsOK(SWIG_ConvertPtr($input, NULL, $1_descriptor, 0))
    || OT::canConvertCollectionObjectFromPySequence< OTMESHING::Cylinder >($input);
}

%apply const CylinderCollection & { const OTMESHING::IntersectionMesher::CylinderCollection & };

// SampleCollection typemaps
%typemap(in) const SampleCollection & {
  if (SWIG_IsOK(SWIG_ConvertPtr($input, (void **) &$1, $1_descriptor, 0))) {
    // From interface class, ok
  } else {
    try {
      $1 = OT::buildCollectionFromPySequence< OT::Sample >($input);
    } catch (const OT::InvalidArgumentException &) {
      SWIG_exception(SWIG_TypeError, "Object passed as argument is not convertible to a collection of Sample");
    }
  }
}

%typemap(typecheck,precedence=SWIG_TYPECHECK_POINTER) const SampleCollection & {
  $1 = SWIG_IsOK(SWIG_ConvertPtr($input, NULL, $1_descriptor, 0))
    || OT::canConvertCollectionObjectFromPySequence< OT::Sample >($input);
}

%apply const SampleCollection & { const OTMESHING::IntersectionMesher::SampleCollection & };

%include IntersectionMesher_doc.i

%include otmeshing/IntersectionMesher.hxx

%copyctor OTMESHING::IntersectionMesher;
