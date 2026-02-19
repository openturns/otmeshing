//                                               -*- C++ -*-
/**
 *  @brief Mesh domain with distance
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
#ifndef OTMESHING_MESHDOMAIN2_HXX
#define OTMESHING_MESHDOMAIN2_HXX

#include <openturns/MeshDomain.hxx>
#include "otmeshing/otmeshingprivate.hxx"

namespace OTMESHING
{

/**
 * @class MeshDomain2
 */
class OTMESHING_API MeshDomain2
  : public OT::MeshDomain
{
  CLASSNAME
public:

  /** Default constructor */
  MeshDomain2();

  /** Default constructor */
  explicit MeshDomain2(const OT::Mesh & mesh);

   /** Virtual constructor */
  MeshDomain2 * clone() const override;

  /** Compute the Euclidean distance from a given point to the domain */
  OT::Scalar computeDistance(const OT::Point & point) const override;
  OT::Sample computeDistance(const OT::Sample & point) const override;

protected:
  
private:

}; /* class MeshDomain2 */

}

#endif /* OTMESHING_MESHDOMAIN2_HXX */
