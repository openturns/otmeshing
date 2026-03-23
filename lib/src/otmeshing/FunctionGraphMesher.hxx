//                                               -*- C++ -*-
/**
 *  @brief Function meshing algorithm
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
#ifndef OTMESHING_FUNCTIONGRAPHMESHER_HXX
#define OTMESHING_FUNCTIONGRAPHMESHER_HXX

#include <openturns/Mesh.hxx>
#include <openturns/KDTree.hxx>
#include <openturns/Function.hxx>

#include "otmeshing/otmeshingprivate.hxx"

namespace OTMESHING
{

/**
 * @class FunctionGraphMesher
 */
class OTMESHING_API FunctionGraphMesher
  : public OT::PersistentObject
{
  CLASSNAME
public:

  /** Default constructor */
  FunctionGraphMesher();

  /** Parameters constructor */
  FunctionGraphMesher(const OT::Interval & inputInterval,
                      const OT::Indices & inputDiscretization);

   /** Virtual constructor */
  FunctionGraphMesher * clone() const override;

  /** String converter */
  OT::String __repr__() const override;

  /** String converter */
  OT::String __str__(const OT::String & offset = "") const override;

  /** Generate mesh */
  virtual OT::Mesh build(const OT::Function & function,
                         const OT::UnsignedInteger outputIndex,
                         const OT::Scalar minOutput,
                         const OT::Scalar maxOutput,
                         const OT::UnsignedInteger outputDiscretization = 1,
                         const OT::Bool subGraph = true) const;
                         
  /** Method save() stores the object through the StorageManager */
  void save(OT::Advocate & adv) const override;

  /** Method load() reloads the object from the StorageManager */
  void load(OT::Advocate & adv) override;

protected:
  void initialize();

  OT::Interval inputInterval_;
  OT::Indices inputDiscretization_;

  OT::Point minInput_;
  OT::Point maxInput_;
  OT::Sample inputVertices_;
  OT::KDTree kdTree_;
  
private:

}; /* class FunctionGraphMesher */

}

#endif /* OTMESHING_FUNCTIONGRAPHMESHER_HXX */
