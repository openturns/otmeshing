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
#ifndef OTMESHING_CYLINDER_HXX
#define OTMESHING_CYLINDER_HXX

#include <openturns/Mesh.hxx>
#include "otmeshing/otmeshingprivate.hxx"

namespace OTMESHING
{

/**
 * @class Cylinder
 */
class OTMESHING_API Cylinder
  : public OT::PersistentObject
{
  CLASSNAME
public:

  /** Default constructor */
  Cylinder();

  /** Parameters constructor */
  Cylinder(const OT::Mesh & base,
           const OT::Interval & extension,
           const OT::Indices & injection,
           const OT::UnsignedInteger discretization);

   /** Virtual constructor */
  Cylinder * clone() const override;

  /** String converter */
  OT::String __repr__() const override;

  /** String converter */
  OT::String __str__(const OT::String & offset = "") const override;

  /** Vertices accessor */
  OT::Sample getVertices() const;

  /** BBox accessor */
  OT::Interval getBoundingBox() const;

  /** Volume accessor */
  OT::Scalar getVolume() const;

  /** Method save() stores the object through the StorageManager */
  void save(OT::Advocate & adv) const override;

  /** Method load() reloads the object from the StorageManager */
  void load(OT::Advocate & adv) override;

protected:
  void initialize();

  OT::Point combine(const OT::Point & pBase, const OT::Point & pExtension) const;
  
  OT::Mesh base_;
  OT::Interval extension_;
  OT::Indices injection_;
  OT::UnsignedInteger discretization_;

  OT::UnsignedInteger baseDimension_ = 0;
  OT::UnsignedInteger extensionDimension_ = 0;
  OT::UnsignedInteger dimension_ = 0;
  OT::Indices complement_;
private:

}; /* class Cylinder */

}

#endif /* OTMESHING_CYLINDER_HXX */
