import os
import matplotlib.pyplot as plt
from kuibit.simdir import SimDir
import kuibit.visualize_matplotlib as viz

datadir = os.environ["HOME"]+"/simulations/tov_ET"
datadir2 = os.environ["HOME"]+"/simulations/tov_ET2"

sim = SimDir(datadir)
sim2 = SimDir(datadir2)

iteration_number = 0
rho = sim.timeseries.maximum['rho']
rho2 = sim2.timeseries.maximum['rho']
rho_xy = sim.gridfunctions.xy["rho"][iteration_number]
rho_xy2 = sim2.gridfunctions.xy["rho"][iteration_number]

# We cannot plot rho_xy directly because it contains all
# the information for the various refinement levels. 
# We need to resample the data onto a uniform grid.

# shape is the resolution at which we resample
# x0, x1 are the bottom left and top right coordiantes
# that we want to consider

# Here we choose x0=[0,0] because we have reflection 
# symmetry

# resample=True activates multilinear resampling

rho_xy_unif = rho_xy.to_UniformGridData(shape=[100, 100], 
                                        x0=[0,0],
                                        x1=[10, 10],
                                        resample=True)

# Undo reflection symmetry on the x axis
rho_xy_unif.reflection_symmetry_undo(dimension=0)
# Undo reflection symmetry on the y axis
rho_xy_unif.reflection_symmetry_undo(dimension=1)

rho_xy_unif2 = rho_xy2.to_UniformGridData(shape=[100, 100], 
                                        x0=[0,0],
                                        x1=[10, 10],
                                        resample=True)

# Undo reflection symmetry on the x axis
rho_xy_unif2.reflection_symmetry_undo(dimension=0)
# Undo reflection symmetry on the y axis
rho_xy_unif2.reflection_symmetry_undo(dimension=1)

viz.plot_color(rho_xy_unif - rho_xy_unif2,
               logscale=True,
               colorbar=True,
               label="rho",
               xlabel="x",
               ylabel="y",
              )
