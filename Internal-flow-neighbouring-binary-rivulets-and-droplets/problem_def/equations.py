from pyoomph import *
from pyoomph.expressions import *
from pyoomph.utils.dropgeom import DropletGeometry

# Equations for the lubrication droplet
class QuasiStationaryLubricationEquationsWithEvaporation(Equations):
    """
    Lubrication-model equations for a quasi-stationary evaporating droplet.

    Parameters:
    - base_radius: ExpressionOrNum, droplet base radius 
    - total_evap: ExpressionOrNum, evaporation/source term for the solute
    - theta: Expression, contact angle (small parameter used in scaling)
    - Ma: GlobalParameter, Marangoni number (controls surface-tension-driven flow)
    - taylor_dispersion: bool, include Taylor-dispersion corrections if True

    Fields defined:
    - pressure: pressure field used in lubrication flux
    - h: film height (spherical-cap profile set as Dirichlet BC)
    - c: solute concentration
    """
    def __init__(self,base_radius:ExpressionOrNum,total_evap:ExpressionOrNum,theta:Expression,Ma:GlobalParameter,velocity_height_for_output:float=0.001):
        super().__init__()
        # Model parameters and switches
        self.base_radius = base_radius
        self.total_evap = total_evap
        self.theta = theta
        self.Ma = Ma
        self.velocity_height_for_output=velocity_height_for_output

        # Surface-tension variation (Marangoni contribution) assumed proportional to concentration
        # surface_tension here represents the gradient of surface tension driving Marangoni flows
        self.surface_tension = self.Ma * self.theta * var("c")

    def height(self):
        """
        Compute the spherical-cap height profile h(r) as a subexpression.

        Uses DropletGeometry to obtain the curvature radius of the spherical cap,
        then computes the vertical height above the substrate. A maximum() is used
        to avoid taking sqrt of a negative number for r beyond the cap (numerical safety).
        """
        rsqr = dot(var("coordinate"), var("coordinate"))  # r^2 = x^2 + y^2 in 2D/3D coordinates
        curvrad = DropletGeometry(base_radius=self.base_radius, contact_angle=self.theta, evalf=False).curv_radius
        # spherical-cap height: (const - curvrad + sqrt(curvrad^2 - r^2)) / theta
        h = ((1 - cos(self.theta)) / sin(self.theta) * self.base_radius - curvrad + square_root(maximum(curvrad**2 - rsqr, 0))) / self.theta
        return h

    def define_fields(self):
        # Scalar fields used by the PDE system
        self.define_scalar_field("pressure", space="C2")
        self.define_scalar_field("h", space="C2")      # film height (Dirichlet from height())
        self.define_scalar_field("c", space="C2")      # solute concentration
        self.define_vector_field("u", space="C1")       # height-averaged velocity field
        self.define_vector_field("uh", space="C1")      # u at the height of the droplet

    def define_residuals(self):
        # Pressure test and trial functions
        p, ptest = var_and_test("pressure")

        # Height is a known profile -> treat as a subexpression and apply Dirichlet BC
        h = subexpression(self.height())
        self.set_Dirichlet_condition("h", h)

        # Lubrication flux Q = - (h^3 / 3) * grad(p) + (h^2 / 2) * grad(surface_tension)
        # First term: pressure-driven Poiseuille-like flow
        # Second term: Marangoni-driven surface-tension flow
        Q = -h**3 / 3 * grad(p) + h**2 / 2 * grad(self.surface_tension)

        # Pressure equation from lubrication: div(Q) = 0  (weak form)
        self.add_weak(Q, grad(ptest))

        # Concentration field and its test function
        c, ctest = var_and_test("c")

        # Base molecular diffusivity (nondimensionalized)
        D = 1

        # Taylor-dispersion correction (effective longitudinal dispersion
        # due to velocity profile in thin film), see 10.1103/PhysRevFluids.7.L022001
        delta_m = h**2 / 2 * grad(self.surface_tension)
        delta_c = -h**3 / 3 * grad(p)
        D += self.theta**2 * (2/105 * dot(delta_c, delta_c) + 1/20 * dot(delta_m, delta_c) + 1/30 * dot(delta_m, delta_m))

        # Weak form contributions for the solute transport equation:
        #   time-derivative: h * partial_t(c)
        #   advection: div(Q * c) -> implemented as weak form add_weak(dot(Q, grad(c)), ctest)
        #   diffusion: div(D * h * grad(c))
        #   source/sink: total_evap / theta (evaporation contribution scaled by theta)
        self.add_weak(var("h") * partial_t(c), ctest)
        self.add_weak(dot(Q, grad(c)), ctest)
        self.add_weak(D * h * grad(c), grad(ctest))
        self.add_weak(self.total_evap / self.theta, ctest)

        # u is velocity at height z and uh is velocity at interface height h
        # This is for output/diagnostics only. We compute u and uh based on parabolic profile of velocity in lubrication flow.
        u,utest=var_and_test("u")
        uh, uhtest=var_and_test("uh")
        u1=(grad(self.surface_tension)-h*grad(p))
        u2=grad(p)
        z=self.velocity_height_for_output
        udef=u1*z+u2/2*z**2        
        uhdef=u1+u2/2    
        self.add_weak(u-udef,utest)
        self.add_weak(uh-uhdef,uhtest)


# Helper equations to understand Marangoni vs pressure-gradient contributions
class ProjectPressureGradientForMarangoniIndicator(Equations):
    """
    Project the pressure gradient into a vector field M for diagnostics.

    M stores grad(p) (the pressure-gradient contribution). In the lubrication
    flux Q the pressure-driven contribution appears as -h^3/3 * grad(p), so
    comparing M with grad(surface_tension) lets one determine whether the height-averaged
    flow is dominated by Marangoni (surface-tension) forcing or by the
    pressure-gradient forcing.
    """
    def define_fields(self):
        # Define a C1 vector field to hold the projected pressure gradient.
        self.define_vector_field("M", "C1")

    def define_residuals(self):
        # Project grad(p) into M via the weak form:
        #   ∫ (M - grad(p)) · Mtest = 0
        M, Mtest = var_and_test("M")
        p = var("pressure")

        # Add weak residual to enforce M ≈ grad(p). This serves only as a
        # diagnostic/indicator field, not a dynamic variable.
        self.add_weak(M - grad(p), Mtest)

# Helper equations to map nodes exactly onto a circle (to improve mesh quality)
class MapNodesOnCircle(InterfaceEquations):
    """
    Move mesh interface nodes onto the unit circle.
    """
    def after_mapping_on_macro_elements(self):
        # Obtain the mesh object handled by this InterfaceEquations instance
        mesh = self.get_mesh()

        # Iterate over nodes and project each one radially onto the unit circle.
        # Uses numpy.sqrt for the radial distance; ensure numpy is imported in the module.
        for n in mesh.nodes():
            x, y = n.x(0), n.x(1)

            # Compute radius r = sqrt(x^2 + y^2)
            r = numpy.sqrt(x**2 + y**2)

            # Guard against the degenerate center node (r == 0) to avoid division by zero.
            # If such a node exists, skip it (it cannot be projected radially).
            if r == 0:
                continue

            # Set node coordinates to (x/r, y/r) so the node lies exactly on the unit circle.
            n.set_x(0, x / r)
            n.set_x(1, y / r)
