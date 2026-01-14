// SWIG file ConvexDecompositionMesher.i

%{
#include "otmeshing/ConvexDecompositionMesher.hxx"
%}

%include ConvexDecompositionMesher_doc.i

%copyctor OTMESHING::ConvexDecompositionMesher;

%include otmeshing/ConvexDecompositionMesher.hxx

%{
#include "openturns/Mesh.hxx"
namespace OT {
  template <>
  struct traitsPythonType<OT::Mesh>
  {
    typedef _PyObject_ Type;
  };

  template <>
  inline
  OT::Mesh
  convert< _PyObject_, OT::Mesh >(PyObject * pyObj)
  {
    void * ptr = 0;
    if (SWIG_IsOK(SWIG_ConvertPtr(pyObj, &ptr, SWIG_TypeQuery("OT::Mesh *"), SWIG_POINTER_NO_NULL))) {
      OT::Mesh * p_it = reinterpret_cast< OT::Mesh * >(ptr);
      return *p_it;
    }
    else {
      throw OT::InvalidArgumentException(HERE) << "Object passed as argument is not convertible to a Mesh";
    }
    return OT::Mesh();
  }
}
%}

%template(_MeshCollection) OT::Collection<OT::Mesh>;
