//                                               -*- C++ -*-
/**
 *  @brief Meshing algorithm for points
 *
 *  Copyright 2005-2026 Airbus-EDF-IMACS-ONERA-Phimeca
 *
 *  This library is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this library.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#ifndef OTMESHING_CONVEXDECOMPOSITIONMESHER_HXX
#define OTMESHING_CONVEXDECOMPOSITIONMESHER_HXX

#include <openturns/PersistentObject.hxx>
#include <openturns/StorageManager.hxx>
#include <openturns/Mesh.hxx>
#include "otmeshing/otmeshingprivate.hxx"

namespace OTMESHING
{

/**
 * @class ConvexDecompositionMesher
 *
 * ConvexDecompositionMesher is some cloudmesher type to illustrate how to add some classes in OpenTURNS
 */
class OTMESHING_API ConvexDecompositionMesher
  : public OT::PersistentObject
{
  CLASSNAME

public:
 
  /** Default constructor */
  ConvexDecompositionMesher();

  /** Virtual constructor method */
  ConvexDecompositionMesher * clone() const override;

  /** example of a func that return a point squared. **/
  OT::Collection<OT::Mesh> build(const OT::Mesh & mesh) const;

  /** String converter */
  OT::String __repr__() const override;

  /** Method save() stores the object through the StorageManager */
  void save(OT::Advocate & adv) const override;

  /** Method load() reloads the object from the StorageManager */
  void load(OT::Advocate & adv) override;

private:

}; /* class ConvexDecompositionMesher */

} /* namespace OTMESHING */

#endif /* OTMESHING_CONVEXDECOMPOSITIONMESHER_HXX */
