#! /usr/bin/env python

import openturns as ot
import openturns.testing as ott
import otmeshing

ot.TESTPREAMBLE()

mesher = otmeshing.IntersectionMesher()
print("mesher=", mesher)

# intersection of two cubes
for dim in range(2, 6):
    mesh1 = ot.IntervalMesher([1] * dim).build(ot.Interval([0.0] * dim, [3.0] * dim))
    mesh2 = ot.IntervalMesher([1] * dim).build(ot.Interval([1.0] * dim, [4.0] * dim))
    intersection = mesher.build(mesh1, mesh2)
    volume = intersection.getVolume()
    print(f"dim={dim} intersection={intersection} volume={volume}")
    ott.assert_almost_equal(volume, 2.0**dim)
