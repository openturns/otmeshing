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
#include "otmeshing/VolumeMesher.hxx"

#include <openturns/PersistentObjectFactory.hxx>
#include <openturns/SpecFunc.hxx>

using namespace OT;

namespace OTMESHING
{

CLASSNAMEINIT(VolumeMesher)

static const Factory<VolumeMesher> Factory_VolumeMesher;

/* Default constructor */
VolumeMesher::VolumeMesher()
  : PersistentObject()
{
  // Nothing to do
}

/* Virtual constructor method */
VolumeMesher * VolumeMesher::clone() const
{
  return new VolumeMesher(*this);
}

/* String converter */
String VolumeMesher::__repr__() const
{
  OSS oss(true);
  oss << "class=" << VolumeMesher::GetClassName()
      << " apexStrategy=" << (apexStrategy_ == CENTROID ? "centroid" : "firstVertex");
  return oss;
}

/* Apex strategy accessor */
void VolumeMesher::setApexStrategy(const ApexStrategy strategy)
{
  apexStrategy_ = strategy;
}

VolumeMesher::ApexStrategy VolumeMesher::getApexStrategy() const
{
  return apexStrategy_;
}

/* Volume meshing from a surface mesh */
Mesh VolumeMesher::build(const Mesh & surface) const
{
  const UnsignedInteger dimension = surface.getDimension();
  if (surface.getIntrinsicDimension() + 1 != dimension)
    throw InvalidArgumentException(HERE) << "VolumeMesher requires a surface mesh (intrinsic dimension = "
                                         << (dimension - 1) << "), got intrinsic dimension "
                                         << surface.getIntrinsicDimension();

  const UnsignedInteger nbVertices = surface.getVerticesNumber();
  if (nbVertices < dimension + 1)
    throw InvalidArgumentException(HERE) << "VolumeMesher requires at least " << (dimension + 1)
                                         << " vertices, got " << nbVertices;

  const Sample vertices(surface.getVertices());
  const IndicesCollection simplices(surface.getSimplices());
  const UnsignedInteger nbFacets = simplices.getSize();

  // Compute apex
  UnsignedInteger apexIndex = 0;
  Sample verticesWithApex(nbVertices + 1, dimension);

  if (apexStrategy_ == CENTROID)
  {
    Point centroid(dimension);
    for (UnsignedInteger i = 0; i < nbVertices; ++ i)
      for (UnsignedInteger j = 0; j < dimension; ++ j)
        centroid[j] += vertices(i, j);
    centroid /= static_cast<Scalar>(nbVertices);

    apexIndex = nbVertices;
    for (UnsignedInteger i = 0; i < nbVertices; ++ i)
      for (UnsignedInteger j = 0; j < dimension; ++ j)
        verticesWithApex(i, j) = vertices(i, j);
    for (UnsignedInteger j = 0; j < dimension; ++ j)
      verticesWithApex(apexIndex, j) = centroid[j];
  }
  else
  {
    apexIndex = 0;
    for (UnsignedInteger i = 0; i < nbVertices; ++ i)
      for (UnsignedInteger j = 0; j < dimension; ++ j)
        verticesWithApex(i, j) = vertices(i, j);
  }

  // Build simplices of dimension d+1 from facets of dimension d-1
  Collection<Indices> simplexColl;
  for (UnsignedInteger i = 0; i < nbFacets; ++ i)
  {
    // Surface facet has dimension vertices in first dimension entries (last is repeated)
    // Build volume simplex: [apex, f[0], ..., f[dimension-1]]
    Indices simplex(dimension + 1);
    simplex[0] = apexIndex;
    for (UnsignedInteger j = 0; j < dimension; ++ j)
      simplex[j + 1] = simplices(i, j);

    // In FIRST_VERTEX mode, skip facets incident to the apex
    if (apexStrategy_ == FIRST_VERTEX)
    {
      Bool incident = false;
      for (UnsignedInteger j = 1; j <= dimension; ++ j)
        if (simplex[j] == apexIndex)
        {
          incident = true;
          break;
        }
      if (incident) continue;
    }

    simplexColl.add(simplex);
  }

  if (simplexColl.isEmpty())
    return Mesh();

  // Compact unused vertices
  Indices usedVertices(verticesWithApex.getSize());
  for (const auto & simplex : simplexColl)
    for (const UnsignedInteger vi : simplex)
      usedVertices[vi] = 1;

  Indices oldToNew(verticesWithApex.getSize());
  Sample verticesCompact(0, dimension);
  for (UnsignedInteger i = 0; i < verticesWithApex.getSize(); ++ i)
  {
    if (usedVertices[i])
    {
      Point p(dimension);
      for (UnsignedInteger j = 0; j < dimension; ++ j)
        p[j] = verticesWithApex(i, j);
      verticesCompact.add(p);
      oldToNew[i] = verticesCompact.getSize() - 1;
    }
  }

  Collection<Indices> simplexCollCompact;
  for (const auto & simplex : simplexColl)
  {
    Indices simplexNew(simplex.getSize());
    for (UnsignedInteger j = 0; j < simplex.getSize(); ++ j)
      simplexNew[j] = oldToNew[simplex[j]];
    simplexCollCompact.add(simplexNew);
  }

  Mesh result(verticesCompact, IndicesCollection(simplexCollCompact));
  result.fixOrientation();
  return result;
}

void VolumeMesher::save(OT::Advocate & adv) const
{
  PersistentObject::save(adv);
  const UnsignedInteger apex = static_cast<UnsignedInteger>(apexStrategy_);
  adv.saveAttribute("apexStrategy_", apex);
}

void VolumeMesher::load(OT::Advocate & adv)
{
  PersistentObject::load(adv);
  UnsignedInteger apex = 0;
  adv.loadAttribute("apexStrategy_", apex);
  apexStrategy_ = static_cast<ApexStrategy>(apex);
}

} /* namespace OTMESHING */
