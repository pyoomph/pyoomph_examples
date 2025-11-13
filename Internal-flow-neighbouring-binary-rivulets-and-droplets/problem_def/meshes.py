from pyoomph import *
from pyoomph.expressions import *
from pyoomph.meshes.remesher import Remesher2d
from pyoomph.utils.dropgeom import DropletGeometry

'''
Mesh class for a liquid rivulet on a substrate surrounded by gas.

Constructs a 2D Gmsh template describing:
 - a liquid rivulet (with an apex and contact lines to the substrate)
 - the surrounding gas domain up to a far-field radius

The mesh is built from points, circular arcs and line segments, then
grouped into plane surfaces for "liquid" and "gas". A Remesher2d is
attached for later adaptive refinement.
'''

# Simple mesh class of a liquid rivulet on a substrate surrounded by gas
class RivuletMesh(GmshTemplate):
    def __init__(
        self,
        base_radius:ExpressionOrNum,
        apex_height:ExpressionOrNum,
        b:GlobalParameter,
        gas_radius:ExpressionOrNum=5,
        hpre:ExpressionOrNum=0.005,
        resolution_factor:ExpressionOrNum=1.0,
        resolutions:dict=None
    ):
        """
        Initialize the rivulet mesh template.

        Parameters
        - base_radius: physical base radius of the rivulet (used to scale geometry)
        - apex_height: height of the rivulet apex above substrate
        - b: global parameter controlling neighbor distance (dimensionless based on base_radius)
        - gas_radius: far-field gas radius (scaled by base_radius)
        - hpre: pre-wet film thickness (scaled by base_radius)
        - resolution_factor: global multiplier for point sizes (mesh sizing)
        - resolutions: optional dict of per-region sizing factors
        """
        super().__init__()
        # Attach a 2D remesher to support future adaptive remeshing operations
        self.remesher = Remesher2d(self)

        # store geometric and meshing parameters
        self.base_radius = base_radius
        self.apex_height = apex_height
        self.b = b
        self.gas_radius = gas_radius
        self.hpre = hpre
        self.resolution_factor = resolution_factor
        self.resolutions = resolutions

        # default resolution values for key regions if none provided
        if resolutions is None:
            self.resolutions = {
                "axis": 0.5,          # mesh size along axis lines
                "far": 5,             # coarse far-field mesh size
                "contact_line": 0.05, # fine mesh near contact lines
                "apex": 0.2,          # mesh size near apex
                "center": 0.2         # mesh size near rivulet center
            }

    def define_geometry(self) -> None:
        """
        Define the Gmsh geometry: points, arcs, lines and surfaces.
        """

        # Use triangular elements
        self.mesh_mode = "tris"

        # scale common geometric quantities by base_radius where appropriate
        rc = self.base_radius                # base radius (scaling)
        ha = self.apex_height                # apex height (scaled later)
        b = self.b.value * rc  # neighbor distance offset
        R0 = self.gas_radius * rc            # far-field radius
        hpre = self.hpre * rc                # pre-wet film thickness
        res = self.resolution_factor         # global resolution multiplier

        # --- Points for the liquid puddle / rivulet geometry ---
        # right and left contact-line points on substrate (where liquid meets substrate)
        pbr0 = self.point(b/2 + rc, 0, size=res * self.resolutions["contact_line"])
        pbmr0 = self.point(b/2 - rc, 0, size=res * self.resolutions["contact_line"])
        # apex point above substrate (highest point of the rivulet)
        pb0h = self.point(b/2, ha, size=res * self.resolutions["apex"])
        # points at the substrate offset by prewet film thickness (below substrate level)
        pbmr0h = self.point(b/2 - rc, -hpre, size=res * self.resolutions["contact_line"])
        pbr0h = self.point(b/2 + rc, -hpre, size=res * self.resolutions["contact_line"])

        # --- Liquid-gas interface approximated by two circular arcs ---
        # arc from right contact to apex through left contact (forming the droplet cap)
        self.circle_arc(pbr0, pb0h, through_point=pbmr0, name="liquid_gas")
        # symmetric arc from left contact to apex through right contact
        self.circle_arc(pbmr0, pb0h, through_point=pbr0, name="liquid_gas")

        # --- Substrate and puddle wall lines ---
        # straight substrate segment between the two lower substrate-offset points
        self.create_lines(pbmr0h, "liquid_substrate", pbr0h)

        # --- Define all lines ---
        lines = {"liquid_gas", "liquid_substrate"}

        # vertical (or slanted) wall segments from contact line down to substrate-offset
        if float(hpre/rc) > 0:
            self.create_lines(pbmr0, "puddle_wall", pbmr0h)
            self.create_lines(pbr0, "puddle_wall", pbr0h)
            lines.add("puddle_wall")

        # Create a closed plane surface for the liquid region (bounded by the above)
        self.plane_surface(*lines, name="liquid")

        # --- Points for the gas / far-field domain ---
        p00 = self.point(0, 0, size=res * self.resolutions["axis"])    # origin / symmetry axis
        pgr0 = self.point(R0, 0, size=res * self.resolutions["far"])   # far point on +x
        pg0H = self.point(0, R0, size=res * self.resolutions["far"])  # far point on +y (axis)
        pgmr0 = self.point(-R0, 0, size=res * self.resolutions["far"]) # far point on -x

        # gas-substrate boundaries: connect origin and right contact to far-field
        self.create_lines(p00, "gas_substrate", pbmr0)
        self.create_lines(pbr0, "gas_substrate", pgr0)
        # axis boundary from origin up to far-field axis point
        self.create_lines(p00, "gas_axis", pg0H)
        # far-field circular arc closing the outer gas boundary
        self.circle_arc(pg0H, pgr0, through_point=pgmr0, name="gas_infinity")

        # Assemble the gas plane surface, subtracting the liquid region
        self.plane_surface("gas_infinity", "gas_axis", "gas_substrate", "liquid_gas", name="gas")


'''
Mesh class for a liquid rivulet on a substrate (no surrounding gas domain).

Constructs a simple 2D Gmsh template describing just the liquid rivulet
bounded by the substrate. The free surface (liquid-gas interface) is
approximated by two symmetric circular arcs computed from DropletGeometry.

This simplified mesh is useful for tests or cases where the external
gas domain is not required.
'''
class RivuletMeshNoGas(GmshTemplate):
    def __init__(self, 
                 base_radius:ExpressionOrNum, 
                 theta:ExpressionOrNum,
                 resolution_factor:ExpressionOrNum=0.02,
                 resolution_interface:ExpressionOrNum=0.25
    ):
        """
        Initialize the rivulet mesh with a given radius and contact angle.

        Parameters
        - R: base radius of the rivulet (used to position contact points)
        - theta: contact angle of the rivulet on the substrate
        - resolution: base mesh size used for most points (smaller => finer mesh)
        - resolution_interface: multiplier applied to resolution for points on the interface
        """
        super().__init__()
        # geometric parameters
        self.base_radius = base_radius
        self.theta = theta

        # meshing parameters
        # - resolution: typical element size in the domain
        # - resolution_interface: relative sizing for interface/contact points
        self.resolution_factor = resolution_factor
        self.resolution_interface = resolution_interface

    def define_geometry(self):
        # set a default element size for the template (used if points don't specify size)
        self.default_resolution = self.resolution_factor

        # use triangular elements for the planar mesh
        self.mesh_mode = "tris"

        # Create the three base points on the substrate/symmetry axis:
        # pl: left contact point at x = -R, y = 0
        # pr: right contact point at x = +R, y = 0
        # p0: origin / substrate midpoint (used to close substrate segments)
        # The interface/contact points use a slightly scaled size to allow finer
        # meshing along the liquid-gas interface if desired.
        pl = self.point(-self.base_radius, 0, size=self.resolution_interface * self.resolution_factor)
        pr = self.point(self.base_radius, 0, size=self.resolution_interface * self.resolution_factor)
        p0 = self.point(0, 0, size=self.resolution_factor)  # keep origin default size (falls back to default_resolution)

        # Compute droplet geometry (apex height, etc.) using the helper class.
        # DropletGeometry encapsulates the cap shape given base radius and contact angle.
        geom = DropletGeometry(base_radius=self.base_radius, contact_angle=self.theta, evalf=False)

        # Apex point: highest point of the rivulet (on symmetry axis at x=0)
        pc = self.point(0, geom.apex_height)

        # Create the liquid-gas interface as two symmetric circular arcs.
        # Each circle_arc is defined by its start, end and a 'through_point'
        # which forces the arc to pass via the opposite contact to form the cap.
        # Naming both arcs "liquid_gas" groups them for surface construction.
        self.circle_arc(pl, pc, through_point=pr, name="liquid_gas")
        self.circle_arc(pr, pc, through_point=pl, name="liquid_gas")

        # Create substrate segments: straight lines along the substrate that connect
        # the contacts via the origin. These are grouped under the name "substrate".
        # Note: orientation/order ensures a closed loop when combined with the liquid_gas.
        self.line(pr, p0, name="liquid_substrate")
        self.line(p0, pl, name="liquid_substrate")

        # Define the liquid plane surface bounded by the liquid_gas (top) and liquid_substrate (bottom).
        # The named curve groups "liquid_gas" and "liquid_substrate" form the boundary of this surface.
        self.plane_surface("liquid_gas", "liquid_substrate", name="liquid")



'''
Mesh class for a small contact angle liquid droplet on a substrate, to be solved in the lubrication limit.
A semi-circular droplet shape is used for efficiency.

Constructs a simple 2D Gmsh template describing just the liquid droplet. 
'''

class DropletLubricationMesh(GmshTemplate):
    def __init__(self,
                base_radius:ExpressionOrNum=1,
                mesh_properties:dict={"layer_thickness": 0.2, "mesh_size_center": 0.5, "boundsize": 0.08, "with_macro_elements": True}
               ):

        """
        Initialize the droplet lubrication mesh with given mesh properties.
        
        Parameters
        - base_radius: physical base radius of the droplet (default is 1)
        - mesh_properties: dictionary containing mesh parameters such as layer_thickness, mesh_size_center, and boundsize.
        """
        super().__init__()
        # store geometric and meshing parameters
        self.base_radius=base_radius
        self.mesh_properties=mesh_properties

    def define_geometry(self):

        # set a default element size for the template (used if points don't specify size)
        self.default_resolution = self.mesh_properties["mesh_size_center"]

        # use only quad elements for the planar mesh
        self.mesh_mode = "only_quads"

        sc=self.mesh_properties["mesh_size_center"]
        lt=self.mesh_properties["layer_thickness"]
        Rl=self.base_radius - lt
        bs=self.mesh_properties["boundsize"]

        # Create key points defining the droplet geometry
        p00=self.point(0,0,size=sc)
        pE=self.point(self.base_radius,0,size=bs)
        pW=self.point(-self.base_radius,0,size=bs)
        pN=self.point(0,self.base_radius,size=bs)
        pEm=self.point(Rl,0,size=bs)
        pWm=self.point(-Rl,0,size=bs)
        pNm=self.point(0,Rl,size=bs)

        # Inner part of the droplet mesh
        E20m=self.line(pEm,p00,name="symmetry")
        W20m=self.line(pWm,p00,name="symmetry")
        N20m=self.line(pNm,p00)
        E2Nm=self.circle_arc(pEm,pNm,center=p00)
        W2Nm=self.circle_arc(pWm,pNm,center=p00)
        self.plane_surface(N20m,E20m,E2Nm,name="liquid")
        self.plane_surface(W20m,N20m,W2Nm,name="liquid")

        # Outer part of the droplet mesh
        E20o=self.line(pEm,pE,name="symmetry")
        W20o=self.line(pWm,pW,name="symmetry")
        N20o=self.line(pNm,pN)
        E2No=self.circle_arc(pE,pN,center=p00,name="contact_line")
        W2No=self.circle_arc(pW,pN,center=p00,name="contact_line")
        ps1=self.plane_surface(E20o,E2Nm,E2No,N20o,name="liquid")
        ps2=self.plane_surface(W2Nm,W20o,W2No,N20o,name="liquid")

        # Apply transfinite meshing to key lines and surfaces for structured mesh
        self.make_lines_transfinite(E2Nm,W2Nm,E2No,W2No,numnodes=math.ceil(float(pi/2/bs)))
        self.make_lines_transfinite(E20o,W20o,N20o,numnodes=4*math.ceil(float(lt/bs)))
        self.make_surface_transfinite(ps1,corners=[pN,pNm,pEm,pE])
        self.make_surface_transfinite(ps2,corners=[pN,pNm,pWm,pW])
        self.use_macro_elements=self.mesh_properties.get("with_macro_elements",True)