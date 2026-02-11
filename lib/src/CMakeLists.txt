

ot_add_current_dir_to_include_dirs ()

ot_add_source_file (CloudMesher.cxx)
ot_add_source_file (ConvexDecompositionMesher.cxx)
ot_add_source_file (ConvexHullMesher.cxx)
ot_add_source_file (Cylinder.cxx)
ot_add_source_file (PolygonMesher.cxx)
ot_add_source_file (IntersectionMesher.cxx)

ot_install_header_file (CloudMesher.hxx)
ot_install_header_file (ConvexDecompositionMesher.hxx)
ot_install_header_file (ConvexHullMesher.hxx)
ot_install_header_file (Cylinder.hxx)
ot_install_header_file (PolygonMesher.hxx)
ot_install_header_file (IntersectionMesher.hxx)

include_directories (${INTERNAL_INCLUDE_DIRS})

add_library (otmeshing ${SOURCEFILES})
set_target_properties (otmeshing PROPERTIES POSITION_INDEPENDENT_CODE ON)
if (BUILD_SHARED_LIBS)
  set_target_properties (otmeshing PROPERTIES COMPILE_DEFINITIONS "OTMESHING_DLL_EXPORTS")
endif ()

if (NOT DEFINED LIB_VERSION) 
  set (LIB_VERSION 0.0.0)
endif ()
if (NOT DEFINED LIB_SOVERSION) 
  set (LIB_SOVERSION 0 )
endif ()
set_target_properties (otmeshing PROPERTIES VERSION ${LIB_VERSION})
set_target_properties (otmeshing PROPERTIES SOVERSION ${LIB_SOVERSION})
target_link_libraries (otmeshing PUBLIC ${OPENTURNS_LIBRARY})
target_link_libraries (otmeshing PRIVATE Eigen3::Eigen)
target_link_libraries (otmeshing PRIVATE CGAL::CGAL)

if (Qhull_FOUND)
  target_compile_definitions (otmeshing PRIVATE OPENTURNS_HAVE_QHULL)
  target_include_directories(otmeshing PRIVATE ${QHULL_INCLUDE_DIRS})
  target_link_libraries(otmeshing PRIVATE ${QHULL_LIBRARIES})
endif ()

if (cddlib_FOUND)
  target_compile_definitions (otmeshing PRIVATE OPENTURNS_HAVE_CDDLIB)
  target_include_directories(otmeshing PRIVATE ${CDDLIB_INCLUDE_DIRS})
  target_link_libraries(otmeshing PRIVATE ${CDDLIB_LIBRARIES})
endif ()

# Add targets to the build-tree export set
export (TARGETS otmeshing FILE ${PROJECT_BINARY_DIR}/otmeshing-Targets.cmake)

# Install the export set for use with the install-tree
install(EXPORT otmeshing-Targets
        DESTINATION "${OTMESHING_CONFIG_CMAKE_PATH}"
        COMPONENT Development)

install(TARGETS otmeshing
        EXPORT  otmeshing-Targets 
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION ${OTMESHING_LIBRARY_PATH}
        ARCHIVE DESTINATION ${OTMESHING_LIBRARY_PATH}
)

