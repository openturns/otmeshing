//                                               -*- C++ -*-
/**
 *  @brief Union meshing
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
#ifndef OTMESHING_UNIONMESHER_HXX
#define OTMESHING_UNIONMESHER_HXX

#include <openturns/Mesh.hxx>
#include "otmeshing/otmeshingprivate.hxx"

namespace OTMESHING
{

/**
 * @class UnionMesher
 */
class OTMESHING_API UnionMesher
  : public OT::PersistentObject
{
  CLASSNAME
public:
  typedef OT::Collection<OT::Mesh> MeshCollection;

  /** Default constructor */
  UnionMesher();

   /** Virtual constructor */
  UnionMesher * clone() const override;

  /** String converter */
  OT::String __repr__() const override;

  /** String converter */
  OT::String __str__(const OT::String & offset = "") const override;

  /** Generate mesh */
  virtual OT::Mesh build(const MeshCollection & coll) const;

  /** Deduplicate vertices */
  static OT::Mesh CompressMesh(const OT::Mesh & mesh);

  /** Method save() stores the object through the StorageManager */
  void save(OT::Advocate & adv) const override;

  /** Method load() reloads the object from the StorageManager */
  void load(OT::Advocate & adv) override;

protected:
private:

}; /* class UnionMesher */

}

#endif /* OTMESHING_UNIONMESHER_HXX */
