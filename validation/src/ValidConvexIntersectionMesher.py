import openturns as ot
import otmeshing as otm
from time import time

useConvex = False
resName = "res_convex.csv" if useConvex else "res_general.csv"
dimMax = 11
NMax = 100
tMax = 60.0
res = ot.Sample(0, 3)
res.setDescription(["dimension", "N", "t"])

for dim in range(2, dimMax + 1):
    N = 1
    t = 0.0
    while (N <= NMax) and (t <= tMax):
        disc = [N, N] + [1]*(dim-2)
        I1 = ot.Interval([-5.0]*dim, [5.0]*dim)
        mesh1 = ot.IntervalMesher(disc).build(I1)
        I2 = ot.Interval([-2.5]*dim, [7.5]*dim)
        mesh2 = ot.IntervalMesher(disc).build(I2)
        print("Intersect, N=", N, "dim=", dim)
        if useConvex:
            algo = otm.ConvexIntersectionMesher()
        else:
            algo = otm.IntersectionMesher()
            #algo.setRecompress(False) # There is a bug here
        t0 = time()
        inter12 = algo.build(mesh1, mesh2)
        t1 = time()
        t = t1 - t0
        print("t=", t1 - t0, "s")
        vol1 = inter12.getVolume()
        print("inter volume=", vol1)
        vol2 = (I1.intersect(I2)).getVolume()
        print("inter volume=", vol2)
        err = abs(1 - vol1 / vol2)
        print("err=%.3e" % err)
        if err > 1e-10:
            raise Exception(f"Error for {N=} {dim=}, got {vol1=} and {vol2=}, {err=}")
        print("#"*50)
        res.add([dim, N, t])
        res.exportToCSVFile(resName)
        N += 1
    if (N == 1) and (t > tMax):
        break
