from pyoomph import *
from pyoomph.expressions import *
from pyoomph.materials import *
from pyoomph.utils.dropgeom import *

# Class storing the Liquid properties
class _LiquidSpecs:
    def __init__(self) -> None:
        # Geometric properites: Set exactly two of these
        self.base_radius:ExpressionNumOrNone=None
        self.contact_angle:ExpressionNumOrNone=None
        self.volume:ExpressionNumOrNone=None
        self.apex_height:ExpressionNumOrNone=None
        # Mixture. Must be set
        self.mixture:AnyLiquidProperties=None
        self.pinned:bool=True # Pinned or unpinned contact angle
        # Radius of curvature, read only
        self.curv_radius:ExpressionNumOrNone=None
        
    # Calculate missing properties and perform checks        
    def _finalize(self):
        geom=DropletGeometry(base_radius=self.base_radius,apex_height=self.apex_height,contact_angle=self.contact_angle,volume=self.volume,evalf=False)
        self.base_radius,self.contact_angle,self.volume,self.apex_height,self.curv_radius=geom.base_radius,geom.contact_angle,geom.volume,geom.apex_height,geom.curv_radius
        if self.mixture is None:
            raise RuntimeError("Liquid mixture must be defined")
        if float(self.contact_angle-90*degree)>0:
            raise RuntimeError("Contact angle must be less than 90 degrees")
        
# Class storing all gas properties
class _GasSpecs:
    def __init__(self) -> None:
        # Gas mixture: Must be set
        self.mixture:AnyGasProperties=None
        
    # Perform checks
    def _finalize(self):
        if self.mixture is None:
            raise RuntimeError("Gas mixture must be defined")