#! /usr/bin/env python

import openturns as ot
import openturns.testing as ott
import otmeshing

ot.TESTPREAMBLE()

mesher = otmeshing.UnionMesher()
print(mesher)
print(repr(mesher))
for dim in range(2, 6):
    mesh1 = ot.IntervalMesher([1] * dim).build(ot.Interval(dim))
    mesh2 = ot.IntervalMesher([1] * dim).build(ot.Interval([2.0] * dim, [3.0] * dim))
    union = mesher.build([mesh1, mesh2])
    print(union)
    assert union.isValid()
    volume = union.getVolume()
    print(f"{dim=} {volume=}")
    ott.assert_almost_equal(union.getVolume(), 2.0)
