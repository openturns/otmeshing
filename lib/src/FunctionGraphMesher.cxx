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
#include <openturns/PersistentObjectFactory.hxx>
#include <openturns/IntervalMesher.hxx>

#include "otmeshing/FunctionGraphMesher.hxx"

using namespace OT;

namespace OTMESHING
{

CLASSNAMEINIT(FunctionGraphMesher)
static const Factory<FunctionGraphMesher> Factory_FunctionGraphMesher;


/* Default constructor */
FunctionGraphMesher::FunctionGraphMesher()
  : PersistentObject()
{
  // Nothing to do
}

/* Parameters constructor */
FunctionGraphMesher::FunctionGraphMesher(const Interval & inputInterval,
                      const Indices & inputDiscretization)
  : PersistentObject()
  , inputInterval_(inputInterval)
  , inputDiscretization_(inputDiscretization)
{
  initialize();
}

void FunctionGraphMesher::initialize()
{
  if (inputDiscretization_.getSize() != inputInterval_.getDimension())
    throw InvalidArgumentException(HERE) << "FunctionGraphMesher expected a discretization of dimension " <<  inputInterval_.getDimension() << " got " << inputDiscretization_.getSize();
  minInput_ = inputInterval_.getLowerBound();
  maxInput_ = inputInterval_.getUpperBound();
  if (inputInterval_.getDimension())
  {
    inputVertices_ = IntervalMesher(inputDiscretization_).build(inputInterval_).getVertices();
    kdTree_ = KDTree(inputVertices_);
  }
}

/* Virtual constructor */
FunctionGraphMesher * FunctionGraphMesher::clone() const
{
  return new FunctionGraphMesher(*this);
}

/* String converter */
String FunctionGraphMesher::__repr__() const
{
  OSS oss(true);
  oss << "class=" << FunctionGraphMesher::GetClassName();
  return oss;
}

/* String converter */
String FunctionGraphMesher::__str__(const String & ) const
{
  return __repr__();
}

/* Here is the interface that all derived class must implement */
Mesh FunctionGraphMesher::build(const OT::Function & function,
                         const OT::UnsignedInteger outputIndex,
                         const OT::Scalar minOutput,
                         const OT::Scalar maxOutput,
                         const OT::UnsignedInteger outputDiscretization,
                         const OT::Bool subGraph) const
{
  const UnsignedInteger inputDimension = inputInterval_.getDimension();
  if (function.getInputDimension() != inputDimension)
    throw InvalidArgumentException(HERE) << "FunctionGraphMesher expected a function of input dimension " << inputDimension << " got " << function.getInputDimension();
  if (function.getOutputDimension() != 1)
    throw InvalidArgumentException(HERE) << "FunctionGraphMesher expected a function of output dimension 1 got " << function.getOutputDimension();
  
  Point minBound;
  Point maxBound;
  Indices discretization;
  Description  description;
  UnsignedInteger index = 0;
  for (UnsignedInteger i = 0; i <= inputDimension; ++ i)
  {
    if (i == outputIndex)
    {
      minBound.add(minOutput);
      maxBound.add(maxOutput);
      discretization.add(outputDiscretization);
      description.add(function.getOutputDescription()[0]);
    }
    else
    {
      minBound.add(minInput_[index]);
      maxBound.add(maxInput_[index]);
      discretization.add(inputDiscretization_[index]);
      description.add(function.getInputDescription()[index]);
      ++ index;
    }
  }
  const Mesh initialMesh(IntervalMesher(discretization).build(Interval(minBound, maxBound)));
  // Evaluate the function over its input
  Point outputValues = function(inputVertices_).asPoint();
  // Move the vertices of the initial mesh according to these values
  Sample vertices(initialMesh.getVertices());
  Point point(inputDiscretization_.getSize());
  const Scalar a = (subGraph ? minOutput : maxOutput);
  const Scalar b = (subGraph ? maxOutput : minOutput);
  for (UnsignedInteger i = 0; i < vertices.getSize(); ++ i)
  {
    index = 0;
    for (UnsignedInteger j = 0; j < vertices.getDimension(); ++ j)
    {
      if (j != outputIndex)
      {
        point[index] = vertices(i, j);
        ++ index;
      }
      const UnsignedInteger baseIndex = kdTree_.query(point);
      const Scalar value = std::clamp(outputValues[baseIndex], minOutput, maxOutput);
      vertices(i, outputIndex) = a + (value - a) * (vertices(i, outputIndex) - a) / (b - a);
    } // for j
  } // for i
  Mesh result(vertices, initialMesh.getSimplices());
  result.setName(function.getName());
  result.setDescription(description);
  return result;
}

/* Method save() stores the object through the StorageManager */
void FunctionGraphMesher::save(Advocate & adv) const
{
  PersistentObject::save(adv);
  adv.saveAttribute("inputInterval_", inputInterval_);
  adv.saveAttribute("inputDiscretization_", inputDiscretization_);
}

/* Method load() reloads the object from the StorageManager */
void FunctionGraphMesher::load(Advocate & adv)
{
  PersistentObject::load(adv);
  adv.loadAttribute("inputInterval_", inputInterval_);
  adv.loadAttribute("inputDiscretization_", inputDiscretization_);
  initialize();
}

}
