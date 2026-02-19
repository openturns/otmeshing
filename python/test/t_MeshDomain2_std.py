#! /usr/bin/env python

import openturns as ot
import openturns.testing as ott
import otmeshing as otm

ot.TESTPREAMBLE()

for dim in [2, 3]:
    # disjoint domain [0,1] U [2,3]
    mesh1 = ot.IntervalMesher([1] * dim).build(ot.Interval([0.0] * dim, [1.0] * dim))
    mesh2 = ot.IntervalMesher([1] * dim).build(ot.Interval([2.0] * dim, [3.0] * dim))
    mesh = otm.UnionMesher().build([mesh1, mesh2])
    domain = otm.MeshDomain2(mesh)

    # inside block 1
    p = [0.1] * dim
    distance = domain.computeDistance(p)
    print(f"inside {distance=:.6g}")
    ott.assert_almost_equal(distance, -0.1)

    # near block 1
    p = [1.1] * dim
    distance = domain.computeDistance(p)
    print(f"b1 {distance=:.6g}")
    ott.assert_almost_equal(distance, 0.1 * dim**0.5)

    # inside block 2
    p = [2.1] * dim
    distance = domain.computeDistance(p)
    print(f"inside {distance=:.6g}")
    ott.assert_almost_equal(distance, -0.1)

    # near block 2
    p = [1.9] * dim
    distance = domain.computeDistance(p)
    print(f"b2 {distance=:.6g}")
    ott.assert_almost_equal(distance, 0.1 * dim**0.5)
