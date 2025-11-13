from problem_def import *
from problem_def.plotters import *

# Geometry / problem parameters
b = 2.1                              # dimensionless distance between droplet centres
theta = 40*degree                    # contact angle (convert degrees to radians)
base_radius = 1                      # liquid base radius (in mm)
Ma = 1e4                             # Marangoni number

# Mesh properties (can be tuned for accuracy / speed)
mesh_properties = {
            "layer_thickness": 0.2, # thickness of refined layer near contact line
            "mesh_size_center": 0.5, # mesh size at the center of the drop
            "boundsize": 0.15, # size of boundary elements
            "max_refinement_level": 2, # maximum refinement level (applied at contact line)
            "drop_refinement_level": 1, # general refinement level for the drop
            "with_macro_elements": True # use macro elements for the drop domain
        }

# Create and configure a quasi-stationary droplet problem.
with QuasiStationaryLubricationDroplet() as pr:
    # Setup plotter
    extension = "png"  # file extension for output plots
    pr.plotter = QuasiStationaryLubricationDropletPlotter(pr, fileext=extension)

    # Set geometric parameter on the problem object
    pr.base_radius = base_radius
    pr.theta.value = theta
    pr.b.value = b
    pr.Ma.value = Ma

    # Set mesh parameters
    pr.mesh_properties = mesh_properties

    # Run the quasi-stationary simulation:
    pr.solve(max_newton_iterations=20)
    pr.output_at_increased_time()

