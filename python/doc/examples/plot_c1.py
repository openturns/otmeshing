import openturns as ot
import openturns.viewer as otv
import otmeshing as otm
import pyvista as pv

ds1 = pv.read("C_1.vtk")
ds2 = pv.read("C_2.vtk")
ds3 = pv.read("C_3.vtk")
plotter = pv.Plotter(off_screen=False)
plotter.set_background("gray")
plotter.add_mesh(ds1, show_edges=True, color="red", opacity=0.1, line_width=0.5, lighting=True, smooth_shading=False)
plotter.add_mesh(ds2, show_edges=True, color="green", opacity=0.1, line_width=0.5, lighting=True, smooth_shading=False)
plotter.add_mesh(ds3, show_edges=True, color="cyan", opacity=0.02, line_width=0.5, lighting=True, smooth_shading=False)
plotter.view_isometric()
plotter.enable_parallel_projection()
plotter.show_grid(color="white", location="outer")
plotter.show()
