from problem_def import *
from problem_def.plotters import *

# Geometry / problem parameters
b = 2.1                              # dimensionless distance between rivulet centres
theta = 50*degree                    # contact angle (convert degrees to radians)
base_radius = 1.0*milli*meter        # liquid base radius (in mm)

# Environment / composition
relative_humidity = 0.7              # ambient relative humidity (0-1)
mass_fraction_solute = 0.05          # mass fraction of solute in the liquid mixture
solute = "12hexanediol"              # solute name for mixture definition

# Mesh resolutions (can be tuned for accuracy / speed)
resolution_factor = 0.3
resolutions = {
    "axis": 0.5,          # mesh size along axis lines
    "far": 5,             # coarse far-field mesh size
    "contact_line": 0.05, # fine mesh near contact lines
    "apex": 0.2,          # mesh size near apex
    "center": 0.2         # mesh size near rivulet center
}

# Create and configure a transient rivulet problem.
with TransientRivulet() as pr:
    # Setup plotter
    extension = "png"  # file extension for output plots
    pr.plotter = RivuletTransientPlotter(pr, fileext=extension)

    # Set geometric parameter on the problem object
    pr.b.value = b
    pr.liquid.contact_angle = theta
    pr.liquid.base_radius = base_radius

    # Define liquid mixture: pure water (volatile liquid) plus a small fraction of 1,2-hexanediol
    pr.liquid.mixture = Mixture(get_pure_liquid("water") + mass_fraction_solute*get_pure_liquid(solute))

    # Define gas phase mixture using relative humidity: air + water vapor at specified RH.
    # The 'temperature' and 'quantity' arguments ensure the mixture is initialized using RH.
    pr.gas.mixture = Mixture(get_pure_gas("air") + relative_humidity*get_pure_gas("water"), temperature=pr.temperature, quantity="RH")

    # Set mesh parameters
    pr.resolution_factor = resolution_factor
    pr.resolutions = resolutions

    # Run the transient simulation:
    # - total time: 300 minutes
    # - startstep: initial time step (1 second)
    # - temporal_error: allowed relative temporal error
    # - outstep: interval between output saves (50 seconds)
    # - out_initially: do not output at time zero
    pr.run(
        300*minute,
        startstep=1*second,
        temporal_error=1,
        outstep=50*second,
        out_initially=False
    )
 