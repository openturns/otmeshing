//                                               -*- C++ -*-
/**
 *  @brief Cylinder
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

#include "otmeshing/Cylinder.hxx"

using namespace OT;

namespace OTMESHING
{

CLASSNAMEINIT(Cylinder)
static const Factory<Cylinder> Factory_Cylinder;


/* Default constructor */
Cylinder::Cylinder()
  : PersistentObject()
{
  // Nothing to do
}

/* Parameters constructor */
Cylinder::Cylinder(const Mesh & base,
           const Interval & extension,
           const Indices & injection,
           const UnsignedInteger discretization)
: PersistentObject()
, base_(base)
, extension_(extension)
, injection_(injection)
, discretization_(discretization)
{
  initialize();
}

void Cylinder::initialize()
{
  baseDimension_ = base_.getDimension();
  extensionDimension_ = extension_.getDimension();
  dimension_ = baseDimension_ + extensionDimension_;
  complement_ = injection_.complement(dimension_);
  if (injection_.getSize() != extension_.getDimension())
    throw InvalidArgumentException(HERE) << "The injection indices size must be equal to the extension dimension.";
}

/* Virtual constructor */
Cylinder * Cylinder::clone() const
{
  return new Cylinder(*this);
}

/* String converter */
String Cylinder::__repr__() const
{
  OSS oss(true);
  oss << "class=" << Cylinder::GetClassName();
  return oss;
}

/* String converter */
String Cylinder::__str__(const String & ) const
{
  return __repr__();
}

Point Cylinder::combine(const Point & pBase, const Point & pExtension) const
{
  Point p(dimension_);
  for (UnsignedInteger j = 0; j < baseDimension_; ++ j)
    p[complement_[j]] = pBase[j];
  for (UnsignedInteger j = 0; j < extensionDimension_; ++ j)
    p[injection_[j]] = pExtension[j];
  return p;
}

/* Vertices accessor */
Sample Cylinder::getVertices() const
{
  const Indices discretization(extension_.getDimension(), discretization_);
  const Sample baseVertices(base_.getVertices());
  const Sample extensionVertices(IntervalMesher(discretization).build(extension_).getVertices());
  Sample vertices(baseVertices.getSize() * extensionVertices.getSize(), dimension_);
  UnsignedInteger index = 0;
  for (UnsignedInteger i1 = 0; i1 < baseVertices.getSize(); ++ i1)
    for (UnsignedInteger i2 = 0; i2 < extensionVertices.getSize(); ++ i2)
    {
      vertices[index] = combine(baseVertices[i1], extensionVertices[i2]);
      ++ index;
    }
  return vertices;
}

/* BBox accessor */
Interval Cylinder::getBoundingBox() const
{
  return Interval(combine(base_.getVertices().getMin(), extension_.getLowerBound()),
                  combine(base_.getVertices().getMax(), extension_.getUpperBound()));
}

/* Volume accessor */
Scalar Cylinder::getVolume() const
{
  return base_.getVolume() * extension_.getVolume();
}

/* Method save() stores the object through the StorageManager */
void Cylinder::save(Advocate & adv) const
{
  PersistentObject::save(adv);
  adv.saveAttribute("base_", base_);
  adv.saveAttribute("extension_", extension_);
  adv.saveAttribute("injection_", injection_);
  adv.saveAttribute("discretization_", discretization_);
}

/* Method load() reloads the object from the StorageManager */
void Cylinder::load(Advocate & adv)
{
  PersistentObject::load(adv);
  adv.loadAttribute("base_", base_);
  adv.loadAttribute("extension_", extension_);
  adv.loadAttribute("injection_", injection_);
  adv.loadAttribute("discretization_", discretization_);
  initialize();
}

}
