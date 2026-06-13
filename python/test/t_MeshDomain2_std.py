#! /usr/bin/env python

import openturns as ot
import openturns.testing as ott
import otmeshing as otm

ot.TESTPREAMBLE()

# single-point API and sample API
for dim in range(2, 6):
    mesh = ot.IntervalMesher([1] * dim).build(ot.Interval([0.0] * dim, [1.0] * dim))
    domain = otm.MeshDomain2(mesh)

    # inside via single point
    p = [0.2] * dim
    distance = domain.computeDistance(p)
    ott.assert_almost_equal(distance, -0.2)

    # inside via sample
    ps = ot.Sample(1, p)
    distances = domain.computeDistance(ps)
    ott.assert_almost_equal(distances[0, 0], -0.2)

    # outside
    p = [1.2] * dim
    distance = domain.computeDistance(p)
    ott.assert_almost_equal(distance, 0.2 * dim**0.5)

    # far outside
    p = [10.0] * dim
    distance = domain.computeDistance(p)
    ott.assert_almost_equal(distance, 9.0 * dim**0.5)

    # point on the boundary
    p = [0.0] * dim
    distance = domain.computeDistance(p)
    ott.assert_almost_equal(distance, 0.0)

    # wrong dimension
    try:
        domain.computeDistance([0.0] * (dim + 1))
        assert False
    except RuntimeError:
        pass

# disjoint domain [0,1] U [2,3]
for dim in range(2, 6):
    mesh1 = ot.IntervalMesher([1] * dim).build(ot.Interval([0.0] * dim, [1.0] * dim))
    mesh2 = ot.IntervalMesher([1] * dim).build(ot.Interval([2.0] * dim, [3.0] * dim))
    mesh = otm.UnionMesher().build([mesh1, mesh2])
    domain = otm.MeshDomain2(mesh)

    # inside block 1
    p = [0.1] * dim
    distance = domain.computeDistance(p)
    ott.assert_almost_equal(distance, -0.1)

    # near block 1 (outside)
    p = [1.1] * dim
    distance = domain.computeDistance(p)
    ott.assert_almost_equal(distance, 0.1 * dim**0.5)

    # inside block 2
    p = [2.1] * dim
    distance = domain.computeDistance(p)
    ott.assert_almost_equal(distance, -0.1)

    # near block 2 (outside)
    p = [1.9] * dim
    distance = domain.computeDistance(p)
    ott.assert_almost_equal(distance, 0.1 * dim**0.5)

# 3D unit cube with known distances
mesh = ot.IntervalMesher([2] * 3).build(ot.Interval([0.0] * 3, [1.0] * 3))
domain = otm.MeshDomain2(mesh)

# center
distance = domain.computeDistance([0.5] * 3)
ott.assert_almost_equal(distance, -0.5)

# outside
distance = domain.computeDistance([1.5] * 3)
ott.assert_almost_equal(distance, 0.5 * 3.0**0.5)

# edge
distance = domain.computeDistance([0.5, 0.5, 0.0])
ott.assert_almost_equal(distance, 0.0)
