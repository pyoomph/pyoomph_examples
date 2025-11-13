from problem_def import *
from problem_def.plotters import *

# Geometry / problem parameters
b = 2.1                              # dimensionless distance between rivulet centres
theta = 40*degree                    # contact angle (convert degrees to radians)
base_radius = 1                      # liquid base radius (in mm)
Ma = 1e4                             # Marangoni number

# Mesh resolutions (can be tuned for accuracy / speed)
resolution_factor = 0.2
resolution_interface = 0.1

# Create and configure a quasi-stationary rivulet problem.
with QuasiStationaryRivulet() as pr:
    # Setup plotter
    extension = "png"  # file extension for output plots
    pr.plotter = QuasiStationaryRivuletPlotter(pr, fileext=extension)

    # Set geometric parameter on the problem object
    pr.base_radius = base_radius
    pr.theta.value = theta
    pr.b.value = b
    pr.Ma.value = Ma

    # Toggle whether to use flat rivulet approximation for evaporation model
    pr.use_flat_rivulet_limit_for_evap = False

    # Set mesh parameters
    pr.resolution_factor = resolution_factor
    pr.resolution_interface = resolution_interface

    # Run the quasi-stationary simulation:
    pr.solve()
    pr.output_at_increased_time()