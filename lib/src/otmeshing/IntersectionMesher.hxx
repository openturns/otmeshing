//                                               -*- C++ -*-
/**
 *  @brief Intersection meshing algorithm
 *
 *  Copyright 2005-2026 Airbus-EDF-IMACS-ONERA-Phimeca
 *
 *  This library is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this library.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#ifndef OTMESHING_INTERSECTIONMESHER_HXX
#define OTMESHING_INTERSECTIONMESHER_HXX

#include <openturns/Mesh.hxx>
#include "otmeshing/otmeshingprivate.hxx"

namespace OTMESHING
{

/**
 * @class IntersectionMesher
 */
class OTMESHING_API IntersectionMesher
  : public OT::PersistentObject
{
  CLASSNAME
public:

  /** Default constructor */
  IntersectionMesher();

   /** Virtual constructor */
  IntersectionMesher * clone() const override;

  /** String converter */
  OT::String __repr__() const override;

  /** String converter */
  OT::String __str__(const OT::String & offset = "") const override;

  /* Here is the interface that all derived class must implement */
  virtual OT::Mesh build(const OT::Mesh & mesh1, const OT::Mesh & mesh2) const;

  /** Recompression flag accessor */
  void setRecompress(const OT::Bool recompress);
  OT::Bool getRecompress() const;

protected:
  OT::Bool recompress_ = true;
private:

}; /* class IntersectionMesher */

}

#endif /* OPENTURNS_INTERVALMESHER_HXX */
