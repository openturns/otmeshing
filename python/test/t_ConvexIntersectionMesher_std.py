#! /usr/bin/env python

import openturns as ot
import openturns.testing as ott
import otmeshing

ot.TESTPREAMBLE()

mesher = otmeshing.ConvexIntersectionMesher()
print("mesher=", mesher)

# intersection of two cubes
for dim in range(2, 6):
    for N in range(1, 25, 5):
        # mesher.setRecompress(compression)
        discr = [N, N] + [1] * (dim - 2)
        mesh1 = ot.IntervalMesher(discr).build(ot.Interval([0.0] * dim, [3.0] * dim))
        mesh2 = ot.IntervalMesher(discr).build(ot.Interval([1.0] * dim, [4.0] * dim))
        intersection = mesher.build(mesh1, mesh2)
        volume = intersection.getVolume()
        print(f"{dim=} {N=} intersection={intersection} {volume=:.3g}")
        ott.assert_almost_equal(volume, 2.0**dim)
