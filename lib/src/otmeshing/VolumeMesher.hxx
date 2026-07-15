//                                               -*- C++ -*-
/**
 *  @brief Volume mesh generation from a surface mesh
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
#ifndef OTMESHING_VOLUMEMESHER_HXX
#define OTMESHING_VOLUMEMESHER_HXX

#include <openturns/PersistentObject.hxx>
#include <openturns/StorageManager.hxx>
#include <openturns/Mesh.hxx>
#include "otmeshing/otmeshingprivate.hxx"

namespace OTMESHING
{

/**
 * @class VolumeMesher
 *
 * Tetrahedralize a surface mesh by fan triangulation from an apex.
 */
class OTMESHING_API VolumeMesher
  : public OT::PersistentObject
{
  CLASSNAME

public:

  /** Strategy for apex choice */
  enum ApexStrategy { CENTROID = 0, FIRST_VERTEX = 1 };

  /** Default constructor */
  VolumeMesher();

  /** Virtual constructor method */
  VolumeMesher * clone() const override;

  /** Build a volume mesh from a surface mesh */
  OT::Mesh build(const OT::Mesh & surface) const;

  /** Apex strategy accessor */
  void setApexStrategy(const ApexStrategy strategy);
  ApexStrategy getApexStrategy() const;

  /** String converter */
  OT::String __repr__() const override;

  /** Method save() stores the object through the StorageManager */
  void save(OT::Advocate & adv) const override;

  /** Method load() reloads the object from the StorageManager */
  void load(OT::Advocate & adv) override;

protected:
  ApexStrategy apexStrategy_ = CENTROID;

}; /* class VolumeMesher */

} /* namespace OTMESHING */

#endif /* OTMESHING_VOLUMEMESHER_HXX */
