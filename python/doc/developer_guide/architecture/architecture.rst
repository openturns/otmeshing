Architecture considerations
===========================

Dependencies
------------

Several dependencies are needed in order to build the module:

 - OpenTURNS >=1.26
 - CGAL >=5.2
 - Qhull >=2020.2 (optional)
 - Sphinx-doc (optional for this doc)

Compilation
-----------

.. code-block:: bash

    cd otmeshing
    cmake \
      -DCMAKE_INSTALL_PREFIX=$PWD/install \
      -DOpenTURNS_DIR=$PWD/../../openturns/build/install/lib/cmake/openturns \
      -B build .
    cmake --build build --target install

Source code structure
---------------------

Here is the global class diagram for each layer:

.. image:: class_diagram.png
    :align: center

