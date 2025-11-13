"""
This module defines a transient rivulet problem consisting of binary rivulet
and gas phases, their coupling, mesh handling, and boundary conditions.
The droplet is pinned on a substrate with a tiny puddle to avoid numerical
issues at the contact line. 
"""

from pyoomph import *
from pyoomph.expressions import *
from pyoomph.equations.poisson import *
from pyoomph.equations.advection_diffusion import *
from pyoomph.equations.navier_stokes import *
from pyoomph.equations.multi_component import *
from pyoomph.equations.ALE import *
from pyoomph.utils.dropgeom import *
from pyoomph.materials.default_materials import *
from .properties import _LiquidSpecs, _GasSpecs
from .meshes import *
from .equations import *


class TransientRivulet(Problem):
    """
    Problem that models a transient rivulet composed of a liquid on a substrate
    coupled to a surrounding gas phase. The problem sets up material properties,
    scaling, meshes and coupled equations (composition diffusion + Navier-Stokes).
    """

    def __init__(self):
        super(TransientRivulet, self).__init__()

        # Physical/specification objects for the liquid and gas
        self.liquid = _LiquidSpecs()  # liquid properties 
        self.gas = _GasSpecs()       # Gas properties 

        # Environmental / initial conditions
        self.temperature = 24 * celsius             # Ambient temperature
        self.absolute_pressure = 1 * atm            # Ambient pressure

        # Global (tunable) parameters exposed to the solver or user
        self.Ma = self.define_global_parameter(Ma=0)         # Marangoni number
        self.b = self.define_global_parameter(b=2.1)         # Centre-to-centre distance / liquid diameter
        self.theta = self.define_global_parameter(theta=50 * degree)  # Initial contact angle

        # Mesh / numerical controls
        self.resolution_factor = 0.4                # Overall mesh resolution factor
        self.gas_radius = 3 * self.b/2              # Gas domain radius, in units of base_radius
        self.far_field_length = 100                 # Far-field extension, in units of base_radius, to impose infinity BCs of logarithmic nature
        self.hpre = 0.005                           # Pre-wetting film height to avoid numerical singularities
        # Local mesh resolutions (can be used by the mesh generator)
        self.resolutions = {"axis": 0.5, "far": 5, "contact_line": 0.05, "apex": 0.2, "center": 0.2}

    def define_problem(self):
        """
        Assemble the full problem: finalize material properties, compute
        reference scales, add meshes and equations for gas and liquid, and
        connect them through interfaces and boundary conditions.
        """

        # Finalize specifications (compute derived params, set up mixture models)
        self.liquid._finalize()
        self.gas._finalize()

        # Get name of a required advected species in the binary liquid phase. Raise an error if not binary mixture.
        w = list(self.liquid.mixture.required_adv_diff_fields)[0]
        if len(self.liquid.mixture.components) != 2:
            raise RuntimeError("The transient rivulet problem only works for binary liquid mixtures.") 

        # Evaluate material properties at the initial/reference condition.
        # These are used to define non-dimensional/scaling factors for the solver.
        rho0 = self.liquid.mixture.evaluate_at_condition("mass_density", cond=self.liquid.mixture.initial_condition, temperature=self.temperature)
        mu0 = self.liquid.mixture.evaluate_at_condition("dynamic_viscosity", cond=self.liquid.mixture.initial_condition, temperature=self.temperature)
        D0 = self.liquid.mixture.evaluate_at_condition(self.liquid.mixture.get_diffusion_coefficient(w), cond=self.liquid.mixture.initial_condition, temperature=self.temperature)
        sigma0 = self.liquid.mixture.evaluate_at_condition(self.liquid.mixture.default_surface_tension["gas"], cond=self.liquid.mixture.initial_condition, temperature=self.temperature)

        # Set characteristic time and scaling based on liquid geometry and material properties.
        # TC chosen from lubrication-like scaling: TC ~ (R^4 * mu) / (H^3 * sigma)
        TC = 3 * self.liquid.base_radius ** 4 * mu0 / (self.liquid.apex_height ** 3 * sigma0)
        self.set_scaling(spatial=self.liquid.base_radius, temporal=TC, pressure=2 * sigma0 / self.liquid.curv_radius, velocity=milli * meter / second, mass_density=rho0)
        # Define named constants for use inside expressions and BCs
        self.define_named_var(temperature=self.temperature, absolute_pressure=self.absolute_pressure)

        # --- Geometry and mesh -------------------------------------------------
        # Add the rivulet mesh to the problem. The mesh object defines domains,
        # markers and (possibly) local refinement strategies.
        self += RivuletMesh(base_radius=self.liquid.base_radius, apex_height=self.liquid.apex_height, b=self.b, gas_radius=self.gas_radius, hpre=self.hpre, resolution_factor=self.resolution_factor, resolutions=self.resolutions)

        # ---------------- Gas-phase equations ---------------------------------
        # Create a container for gas-phase related equations and mesh handling
        geqs = MeshFileOutput()                      # mesh output for gas
        geqs += TextFileOutput()                     # generic text output for gas
        geqs += PseudoElasticMesh()                  # moving mesh for gas domain (ALE)
        geqs += RemeshWhen(RemeshingOptions())       # remeshing criteria (default options)

        # Add composition diffusion (isothermal) for the gas phase.
        # This sets up advection-diffusion equations for gas species mass fractions.
        geqs += CompositionDiffusionEquations(fluid_props=self.gas.mixture, isothermal=True, initial_temperature=self.temperature)

        # Far-field (infinity) boundary conditions for gas composition.
        # req_adv_diff lists the advected/diffusive species required by the gas mixture.
        req_adv_diff = assert_gas_properties(self.gas.mixture).required_adv_diff_fields
        ic = self.gas.mixture.initial_condition
        # Build dictionary of far-field mass fractions for the required species
        farfield = {n: ic["massfrac_" + n] for n in req_adv_diff if "massfrac_" + n in ic.keys()}
        # Composition diffusion to infinity (imposes far-field values at a large x distance)
        geqs += CompositionDiffusionInfinityEquations(farfield_length=self.far_field_length * self.liquid.base_radius, **farfield) @ "gas_infinity"
        geqs += InitialCondition(**{"massfrac_" + n: ic["massfrac_" + n] for n in req_adv_diff if "massfrac_" + n in ic.keys()})

        # Mesh boundary conditions for the gas infinite region and symmetry/ substrate
        geqs += DirichletBC(mesh_x=True, mesh_y=True) @ "gas_infinity"   # fix mesh at infinity
        geqs += DirichletBC(mesh_x=True) @ "gas_axis"                     # axisymmetric axis (if present)
        geqs += DirichletBC(mesh_y=True) @ "gas_substrate"                # substrate line for gas mesh

        # ---------------- Liquid-phase equations ------------------------------
        # Container for liquid-phase equations and outputs
        leqs = MeshFileOutput()
        leqs += TextFileOutput()                     # generic liquid output
        leqs += TextFileOutput() @ "liquid_gas"     # output at liquid-gas interface
        leqs += TextFileOutput() @ "liquid_substrate"  # output at substrate beneath liquid

        # Mesh handling for liquid: moving mesh + remeshing
        leqs += PseudoElasticMesh()
        leqs += RemeshWhen(RemeshingOptions())

        # Substrate BCs for mesh and velocity: fix mesh in y at substrate, prevent slip normal motion
        leqs += DirichletBC(mesh_y=True) @ "liquid_substrate"
        leqs += NoSlipBC() @ "liquid_substrate"

        # Walls bounding the puddle region (fix mesh in x and prevent tangential motion)
        leqs += DirichletBC(mesh_x=True) @ "puddle_wall"
        leqs += DirichletBC(velocity_x=0) @ "puddle_wall"

        # Full Navier-Stokes + composition coupling for the liquid (isothermal)
        leqs += CompositionFlowEquations(fluid_props=self.liquid.mixture, isothermal=True, initial_temperature=self.temperature)

        # Add the multi-component Navier-Stokes interface operator connecting the two fluids
        leqs += MultiComponentNavierStokesInterface(interface_props=self.liquid.mixture | self.gas.mixture) @ "liquid_gas"
        # Ensure mesh connectivity (matching meshes) at the liquid-gas interface
        leqs += ConnectMeshAtInterface() @ "liquid_gas"

        # Initial condition for the advected species inside liquid
        leqs += InitialCondition(**{"massfrac_" + w: self.liquid.mixture.initial_condition["massfrac_" + w]})

        # Enforce kinematic boundary condition at the contact between liquid and puddle wall:
        # the vertical fluid velocity equals the mesh vertical velocity (mesh moves with interface).
        leqs += EnforcedBC(velocity_y=partial_t(var("mesh_y"))-0) @ "liquid_gas/puddle_wall"

        # ---------------- Output functions ------------------------------------
        leqs += IntegralObservables(eqvl_volume=pi*absolute(var("coordinate_x")-self.b*self.liquid.base_radius/2)) # Volume of a corresponding droplet
        leqs += IntegralObservables(_avg_massfrac_w=var("massfrac_"+w),_volume=1)+IntegralObservables(avg_massfrac_w=lambda _avg_massfrac_w,_volume: _avg_massfrac_w/_volume) # Average mass fraction of species w in liquid
        leqs += IntegralObservableOutput(first_column=["time",self.b,self.Ma])

        # ---------------- Final assembly --------------------------------------
        # Add liquid and gas equation groups to the problem. The @-syntax labels
        # the groups/domains for output organization. Both containers are combined
        # under different top-level names to keep outputs separated.
        self.add_equations(leqs @ "liquid" + geqs @ "gas")



class QuasiStationaryRivulet(Problem):
    """
    Problem that models a quasi-stationary rivulet composed of a liquid on a substrate
    coupled to a surrounding gas phase. The problem sets up material properties,
    scaling, meshes and coupled equations (composition diffusion + Stokes/Navier-Stokes).
    """
    def __init__(self):
        super(QuasiStationaryRivulet, self).__init__()
        
        # Toggle to use a simplified flat-droplet limit (does not require gas domain and uses analytical evaporation profile)
        self.use_flat_rivulet_limit_for_evap=False
        
        # Characteristic radius used in analytical/construction parts
        self.base_radius=1

        # Global parameters
        self.theta=self.define_global_parameter(theta=20*degree)   # Contact angle
        self.b=self.define_global_parameter(b=2.1)                 # Centre-to-centre distance / liquid diameter
        self.Ma=self.get_global_parameter("Ma")                    # Marangoni number
        self.cavg=self.define_global_parameter(cavg=0)             # Average concentration constraint for the advected field w
        # Linear surface tension
        self.surface_tension=self.Ma*var("w")

        # Mesh / numerical controls
        self.resolution_factor = 0.4                # Overall mesh resolution factor
        self.gas_radius = 5 * self.b/2              # Gas domain radius, in units of base_radius
        self.far_field_length = 100                 # Far-field extension, in units of base_radius, to impose infinity BCs of logarithmic nature
        self.resolution_interface = 0.05            # Local mesh resolution at the liquid-gas interface

    def asinh(self, x):
        return log(maximum(x  + square_root(x**2 + 1), 1e-10))  # avoid log(0) for small x

    # Compute the base evaporation rate profile (J0) for an isolated rivulet,
    # optionally modified by a neighbor rivulet at distance `b` into J1.
    def get_analytical_evaporation_rate(self):
        rsqr = var("coordinate_x")**2
        # Use expressions from https://doi.org/10.1007/s10665-019-10033-7
        J0 = 1 / (square_root(self.base_radius**2-rsqr) * self.asinh(self.far_field_length/self.base_radius))
        if self.b is None:
            return J0
        else:
            a,b=self.base_radius,self.b
            if float(2*a/b)>=1:
                raise RuntimeError("b must be larger than two times the rivulet radius")
            x = var("coordinate_x")
            eta = b/2+x
            Omega = b/2+a
            Ivar = b/2-a
            L = self.far_field_length * self.base_radius
            J1 = 1 / self.asinh(L / square_root(Omega**2 - Ivar**2)) * eta / (square_root(Omega**2 - eta**2) * square_root(eta**2 - Ivar**2))
            return J1

    def define_problem(self):
        # --- Mesh / Gas-phase treatment -------------------------------------
        # Two alternative setups:
        #  - flat rivulet limit: uses a simplified RivuletMesh and an analytically supplied evaporation
        #  - full rivulet mesh: solves a Poisson problem in the gas to obtain the vapour field
        if self.use_flat_rivulet_limit_for_evap:
            # Simplified geometry/mesh for the flat rivulet limit
            self+=RivuletMeshNoGas(base_radius=self.base_radius,theta=self.theta,resolution_factor=self.resolution_factor,resolution_interface=self.resolution_interface)
            evap=self.get_analytical_evaporation_rate()
        else:   
            # Detailed rivulet mesh, suitable for coupling liquid and gas domains
            geom = DropletGeometry(base_radius=self.base_radius, contact_angle=self.theta, evalf=False)
            mesh = RivuletMesh(base_radius=self.base_radius,apex_height=geom.apex_height,b=self.b,gas_radius=self.gas_radius,hpre=0,resolution_factor=self.resolution_factor,resolutions={"contact_line":self.resolution_interface, "apex":self.resolution_interface, "center": 0.2, "axis": 0.5, "far": 5})
            self.add_mesh(mesh)

            # Container for gas-phase equations / outputs
            geqs = MeshFileOutput()
            geqs+=TextFileOutput()

            # Poisson equation for the vapor concentration (cvap) in the gas domain
            geqs+=PoissonEquation(name="cvap",space="C2")
            # Axisymmetry BCs for gas axis and far-field
            geqs+=AxisymmetryBC()@"gas_axis"
            geqs+=AxisymmetryBC() @ "gas_axis/gas_infinity"
            # Impose far-field Dirichlet value (cvap -> 0) at a large distance
            geqs+=PoissonFarFieldMonopoleCondition(name="cvap", far_value=0, farfield_length=self.far_field_length)@"gas_infinity"
            # Set cvap = 1 on the liquid-gas interface (driving evaporation)
            geqs+=DirichletBC(cvap=1)@"liquid_gas"
            geqs+=InitialCondition(cvap=0)
            # Evaporation flux is the normal derivative of the vapor field at the interface
            evap=-dot(grad(var("cvap",domain="gas")),var("normal"))

        # --- Liquid-phase equations / boundary conditions --------------------
        deqs = MeshFileOutput()
        deqs+=TextFileOutput()
        deqs+=TextFileOutput()@"liquid_gas"

        # Use Stokes (or weakly inertial) solver for the quasi-stationary flow;
        # "TH" denotes a particular Taylor-Hood element choice.
        deqs+=StokesEquations(mode="TH").with_pressure_integral_constraint(self)
        # Slip-length model on the substrate to regularize contact-line motion
        deqs+=NavierStokesSlipLength(sliplength=1e-4)@"liquid_substrate"
        # No-penetration at the substrate (vertical velocity = 0)
        deqs+=DirichletBC(velocity_y=0)@"liquid_substrate"
        # Advection-diffusion for the composition field 'w' with an integral constraint to fix average concentration
        deqs+=AdvectionDiffusionEquations(fieldnames="w",space="C2").with_integral_constraint(self,average=self.cavg)
        # Neumann BC on the interface: mass loss due to evaporation
        deqs+=NeumannBC(w=evap)@"liquid_gas"
        # Start from a quiescent initial state
        deqs+=InitialCondition(velocity_x=0,velocity_y=0)
        # Project the computed evaporation flux to output files for inspection
        deqs+=ProjectExpression(J=evap)@"liquid_gas"  # Output file
        # Marangoni stresses and Laplace pressure at the liquid-gas interface
        deqs+=NavierStokesFreeSurface(surface_tension=self.surface_tension)@"liquid_gas"

        # --- Output and assembly --------------------------------------------
        # Global time/parameter outputs (ODE-style globals) and integral observables
        self += ODEFileOutput()@"globals"
        deqs += IntegralObservables(volume=1)

        # Add liquid and (if present) gas equation groups to the problem
        self.add_equations(deqs@"liquid")
        if not self.use_flat_rivulet_limit_for_evap:
            self.add_equations(geqs@"gas")


class QuasiStationaryLubricationDroplet(Problem):
    """
    Quasi-stationary lubrication-model droplet problem.
    This class sets up a 2D lubrication-type mesh and equations for a single
    droplet in the right position (the neighbor droplet, if any, lies at -x).
    The model uses an analytical evaporation profile (possibly modified by a
    neighbor droplet) and optionally includes Taylor dispersion.
    """

    def __init__(self):
        super().__init__()
        # Characteristic base radius of the droplet (used for scaling/geometry)
        self.base_radius = 1
        # Contact angle of the droplet (global parameter exposed to solver/user)
        self.theta = self.define_global_parameter(theta=20 * degree)
        # Centre-to-centre distance to a neighboring droplet (global parameter)
        self.b = self.define_global_parameter(b=2.1)
        # Marangoni number obtained from the global parameter namespace
        self.Ma = self.get_global_parameter("Ma")
        # Mesh refinement and sizing properties used by the mesh/refinement operators
        self.mesh_properties = {
            "layer_thickness": 0.2,
            "mesh_size_center": 0.5,
            "boundsize": 0.08,
            "max_refinement_level": 2,
            "drop_refinement_level": 1,
            "with_macro_elements": True
        }
        # Height at which to compute velocity for output/diagnostics
        self.velocity_height_for_output = 0.001

    # Analytical evaporation profile for the left droplet. If self.b is set,
    # the expression is modified to account for the neighbour droplet located at +b.
    def get_analytical_evaporation_rate(self):
        # Polar/angle at each point (useful when constructing neighbour-modified factor)
        phi = atan2(var("coordinate_y"), var("coordinate_x"))
        # Radial coordinate squared and radius
        rsqr = dot(var("coordinate"), var("coordinate"))
        r = square_root(rsqr)

        # Base Popov evaporation profile J0 for an isolated droplet
        J0 = get_analytical_popov_evaporation_rate(contact_angle=self.theta, base_radius=self.base_radius)
        # Substitute radial coordinate and film height variable into the analytical expression
        J0 = substitute_in_expression(J0, {"coordinate_x": r, "coordinate_y": var("h")})

        # If a neighbor droplet distance b is specified, apply a correction factor
        if self.b:
            a, b = self.base_radius, self.b
            # Geometric prefactor (from approximate analytical correction)
            F = 4 * a / (1 + (2 / pi) * asin(a / b))
            # Correction factor that depends on the local position (r, phi)
            factor = 1 - F * square_root(b**2 - a**2) / (2 * pi * (rsqr + b**2 - 2 * r * b * cos(phi)))
        else:
            # No neighbour: no modification
            factor = 1

        # Return the (possibly modified) evaporation flux
        return J0 * factor

    def define_problem(self):
        # --- Mesh setup -------------------------------------------------
        # Create and add a lubrication-specific mesh for the droplet domain.
        mesh = DropletLubricationMesh(base_radius=self.base_radius, mesh_properties=self.mesh_properties)
        self.add_mesh(mesh)

        # --- Equation group / outputs container ------------------------
        # Create an equation/output container.
        eqs = MeshFileOutput()
        eqs += MapNodesOnCircle()@"contact_line"  # Move interface nodes onto the unit circle exactly for improved mesh quality

        # Apply mesh refinement operators. We refine up to the configured level
        # generally and also specifically at the contact line region.
        eqs += RefineToLevel(self.mesh_properties["drop_refinement_level"])  # general refinement
        eqs += RefineToLevel(self.mesh_properties["max_refinement_level"]) @ "contact_line"

        # --- Evaporation ------------------------------------------------
        # Obtain the analytical (or neighbour-modified) evaporation flux expression.
        J = self.get_analytical_evaporation_rate()

        # Add the lubrication equations coupled with evaporation and optional
        # Taylor dispersion correction. 
        eqs += QuasiStationaryLubricationEquationsWithEvaporation(
            base_radius=self.base_radius,
            total_evap=J,
            theta=self.theta,
            Ma=self.Ma,
            velocity_height_for_output=0.001
        )

        # --- Global constraints via Lagrange multipliers ----------------
        # Add global Lagrange multipliers for enforcing integral constraints
        # (average concentration and pressure).
        self += GlobalLagrangeMultiplier(avg_c=0, avg_p=0) @ "globals"
        avg_c, avg_p = var("avg_c", domain="globals"), var("avg_p", domain="globals")
        h = var("h")

        # The lubrication formulation couples local fields with global multipliers.
        # Add weak contributions to enforce the integral constraints:
        eqs += WeakContribution(h * var("c"), testfunction(avg_c)) + WeakContribution(avg_c * h, testfunction("c"))
        eqs += WeakContribution(h * var("pressure"), testfunction(avg_p)) + WeakContribution(avg_p * h, testfunction("pressure"))

        # --- Outputs / diagnostics -------------------------------------
        # Project the evaporation flux and flow rate to output files for inspection
        eqs += ProjectExpression(J=J) + ProjectExpression(Q=- h**3 / 3 * grad(var("pressure")) + h**2 / 2 * grad(self.Ma * self.theta * var("c")), space="C1", field_type="vector")

        # Add text outputs for general domain, contact-line region and symmetry axis
        eqs += TextFileOutput() + TextFileOutput() @ "contact_line" + TextFileOutput() @ "symmetry"

        # Project additional diagnostic fields: pressure-gradient-based Marangoni indicator,
        # and a Marangoni indicator computed from M and Q vectors.
        eqs += ProjectPressureGradientForMarangoniIndicator()
        eqs += ProjectExpression(Mindicator=dot(var("M"), var("Q")))

        # --- Final assembly --------------------------------------------
        # Add the liquid equation group labeled as "liquid" and also add the global ODE
        # output group for time/parameter histories.
        self += eqs @ "liquid" + ODEFileOutput() @ "globals"