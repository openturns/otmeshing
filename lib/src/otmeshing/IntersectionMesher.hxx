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
  typedef OT::Collection<OT::Mesh> MeshCollection;

  /** Default constructor */
  IntersectionMesher();

   /** Virtual constructor */
  IntersectionMesher * clone() const override;

  /** String converter */
  OT::String __repr__() const override;

  /** String converter */
  OT::String __str__(const OT::String & offset = "") const override;

  /** intersection */
  virtual OT::Mesh build(const MeshCollection & coll) const;

  /** intersection of convexes */
  virtual OT::Mesh buildConvex(const MeshCollection & coll) const;

  /** Recompression flag accessor */
  void setRecompress(const OT::Bool recompress);
  OT::Bool getRecompress() const;

  /** Method save() stores the object through the StorageManager */
  void save(OT::Advocate & adv) const override;

  /** Method load() reloads the object from the StorageManager */
  void load(OT::Advocate & adv) override;

protected:
  OT::Mesh build2(const OT::Mesh & mesh1, const OT::Mesh & mesh2) const;
  OT::Mesh build2Convex(const OT::Mesh & mesh1, const OT::Mesh & mesh2) const;

  OT::Bool recompress_ = true;
private:

}; /* class IntersectionMesher */

}

#endif /* OPENTURNS_INTERVALMESHER_HXX */
