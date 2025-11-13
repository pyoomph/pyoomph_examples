from pyoomph.output.plotting import *
from matplotlib.patches import ArrowStyle
from pyoomph.utils.dropgeom import DropletGeometry
from .problems import *

# Plotter class specialised for liquid rivulet problem
class RivuletTransientPlotter(MatplotlibPlotter):
    def __init__(self, problem, filetrunk="plot_{:05d}", fileext="pdf", eigenvector=None, eigenmode="abs"):
        # Initialize base class with same signature
        super(RivuletTransientPlotter, self).__init__(problem, filetrunk, fileext, eigenvector, eigenmode)
        # Aspect ratio for the plotted region (since contact angle is small, use wider aspect ratio)
        self.aspect_ratio = 2.5

    def define_plot(self):
        # Retrieve problem instance for convenience
        p = cast(TransientRivulet, self.get_problem())

        # General styling and defaults
        self.background_color = "gray"       # overall figure background
        self.useLaTeXFont()              # use LaTeX fonts for all text
        self.defaults("colorbar").hide_some_ticks = False

        # --- View / axes limits ---
        # Base geometry parameters used to set an appropriate view window
        br = p.liquid.base_radius
        h = p.liquid.apex_height
        b = p.b.value * p.liquid.base_radius
        w = list(p.liquid.mixture.required_adv_diff_fields)[0]
        # Keep some margin around the liquid using base radius fractions
        self.set_view(-(b + 0.35 * br), -0.2 * br, (b + 0.35 * br), 0.75 * br)

        # --- Colorbars definitions ---
        # Gas-phase water mass fraction (% wt) at top-left
        cb_gas = self.add_colorbar("water vapour [$\\%$wt]", cmap="Blues", position="top left", factor=100, length=0.38, thickness=0.06)
        cb_gas.ymargin += 0.025

        # Liquid-phase solute concentration colorbar at bottom-right
        cb_T = self.add_colorbar(w+" [$\\%$wt]", cmap="coolwarm", position="bottom right", factor=100, length=0.38, thickness=0.06)

        # Prepare velocity colorbar with log scaling to handle wide dynamic range
        domain_data = self._get_mesh_data("liquid")
        udata = 1000 * numpy.sqrt(domain_data.get_data("velocity_x")**2 + domain_data.get_data("velocity_y")**2)  # mm/s
        umin, umax = 1e-9, 1e-8 # default min/max to avoid errors
        if udata[udata > 0].size > 0:
            umin = numpy.amin(udata[udata > 0])  # avoid log(0) by taking smallest positive value
            umax = numpy.amax(udata)
        cb_v = self.add_colorbar("velocity [mm/s]", position="bottom left", cmap="viridis", norm=matplotlib.colors.LogNorm(vmin=umin, vmax=umax), vmin=umin, vmax=umax, length=0.38, thickness=0.06, factor=1000)  # factor converts to desired units
        # Small layout tweaks for the colorbars
        cb_v.ymargin += 0.02
        cb_v.xmargin += 0.02
        cb_T.ymargin += 0.02
        cb_T.xmargin += 0.02

        # --- Plots / overlays ---
        # Add concentration field inside liquid (use first required advective-diffusive field)
        self.add_plot("liquid/massfrac_" + w, colorbar=cb_T, transform=None)

        # Streamlines for velocity inside the liquid (white streamlines)
        streams = self.add_plot("liquid/velocity", mode="streamlines", linewidths=6, linecolor="white", transform=[None, "mirror_x"])
        # Velocity field plotted with mirror transformation to show symmetric half
        self.add_plot("liquid/velocity", colorbar=cb_v, transform="mirror_x")

        # Adjust arrow style and density for each streamline set
        for stream in streams:
            stream.arrowstyle = ArrowStyle("-|>", head_length=2.0, head_width=1.8)
            stream.density *= 0.5  # reduce density for clarity

        # --- Arrow key (legend) for evaporation rate vectors ---
        arrk = self.add_arrow_key(title="evap. rate [g/(m$^2$s)]", position="top right", factor=1000, format="{:.3f}")
        arrk.maxlength = 0.55
        arrk.minlength_key = 0.2
        arrk.maxlength_key = 0.6
        # Styling for the arrow key (color and thickness)
        arrk.arrowprops = dict(color="lawngreen", fc="lawngreen", ec="lawngreen", arrowstyle="simple", linewidth=3, mutation_scale=80)
        arrk.ymargin += 0.06
        arrk.xmargin += -0.15

        # Add evaporation rate arrows as a separate plot, mirrored for symmetry
        self.add_plot("liquid/liquid_gas/masstrans_water", arrowkey=arrk, arrowdensity=30,
                      transform=[None, "mirror_x"])

        # Gas-phase water concentration plot, mirrored
        self.add_plot("gas/massfrac_water", colorbar=cb_gas, transform=[None, "mirror_x"])

        # --- Time label ---
        time = self.add_time_label("top center")
        time.format = "{:.0f}"   # show integer time
        time.ymargin += 0.2

        # --- Volume & contact angle annotation ---
        R = float(br / (milli * meter))  # convert base radius to desired units (if needed elsewhere)
        # Compute actual liquid volume (exclude puddle volume)
        V = p.get_mesh("liquid").evaluate_observable("eqvl_volume") - p.hpre * p.liquid.base_radius**3 * pi
        theta = round(float(DropletGeometry(base_radius=br, volume=V).contact_angle / numpy.pi * 180), 1)
        # Add contact angle text at top center inside a white rounded box
        text = self.add_text("$\\theta={}^\circ$".format(float(theta)), position="top center", bbox=dict(facecolor='white', boxstyle='round'))

        # Scale bar at bottom center for reference
        scale_bar = self.add_scale_bar("bottom center")

        # If the liquid is very flat (small contact angle), increase streamline spacing for readability
        if theta < 15:
            for stream in streams:
                stream.density *= 0.7

        # Add a light grey box across the top portion for aesthetic/background of labels
        ybox = self.ymax - 0.2 * (self.ymax - self.ymin)
        self.add_polygon([(self.xmin, self.ymax), (self.xmax, self.ymax), (self.xmax, ybox), (self.xmin, ybox)], edgecolor="black", facecolor="lightgrey", alpha=0.75).zindex = 5

        # --- Text and element sizes for publication-quality figures ---
        # Increase text sizes and tick sizes for all colorbars and keys
        cb_T.textsize = 100
        cb_T.ticsize = 100
        cb_gas.textsize = 100
        cb_gas.ticsize = 100
        cb_v.textsize = 100
        cb_v.ticsize = 100
        arrk.textsize = 100

        # Scale bar visual adjustments
        scale_bar.textsize = 100
        scale_bar.linewidths = 6
        scale_bar.text_yoffset = 0.02

        # Text label and time label sizes
        text.textsize = 100
        time.textsize = 100

# Plotter class for the 2D side view of the droplet
class QuasiStationaryRivuletPlotter(MatplotlibPlotter):
    def __init__(self, problem, filetrunk="plot_{:05d}", fileext="png", eigenvector=None, eigenmode="abs"):
        # Initialize base plotter with same signature
        super(QuasiStationaryRivuletPlotter, self).__init__(problem, filetrunk, fileext, eigenvector, eigenmode)
        # Aspect ratio for the plotted region (side-view is tall and narrow)
        self.aspect_ratio = 2.5

    def define_plot(self):
        # Retrieve problem instance for convenience
        p = cast(QuasiStationaryRivulet, self.get_problem())

        # --- General styling and defaults ---
        # Background and colorbar defaults (transparent figure background for overlays)
        self.background_color = "gray"
        self.useLaTeXFont()              # use LaTeX fonts for all text
        self.defaults("colorbar").hide_some_ticks = False
        offset_x = 0 if not p.use_flat_rivulet_limit_for_evap else p.b.value/2

        # --- View / axes limits ---
        # Base geometry used to set an appropriate view window
        br = p.base_radius
        geom = DropletGeometry(base_radius=br, contact_angle=p.theta)
        h = geom.apex_height
        b = p.b.value
        # Keep small margins around droplet using fractions of geometry
        self.set_view(-(b+0.35*br), -0.25*br, (b+0.35*br), 0.75*br)

        # Draw substrate line (thick black baseline)
        plt.axhline(0, color="black", lw=20)

        # --- Concentration colorbar & plot ---
        # solute mass fraction with fixed range tuned for this dataset
        cb_T = self.add_colorbar("solute [$\\%$wt]", cmap="coolwarm", position="bottom left", length=0.35, thickness=0.1)
        cb_T.ymargin += 0.03
        cb_T.xmargin += 0.05
        # Plot concentration field shifted so droplet sits centered in view
        self.add_plot("liquid/w", colorbar=cb_T, transform=PlotTransformMirror(x=True, offset_x=-offset_x, offset_y=0))

        # --- Velocity colorbar & streamlines ---
        # Prepare velocity data and use log scaling to handle wide dynamic range
        domain_data = self._get_mesh_data("liquid")
        udata = numpy.sqrt(domain_data.get_data("velocity_x")**2 + domain_data.get_data("velocity_y")**2)
        umin = max(numpy.amin(udata), 1e-7)  # avoid log(0) by ensuring positive min
        umax = min(numpy.amax(udata), 1e2)
        cb_v = self.add_colorbar("velocity", position="bottom right", cmap="viridis", norm=matplotlib.colors.LogNorm(vmin=umin, vmax=umax),vmin=umin, vmax=umax, length=0.35, thickness=0.1)
        cb_v.ymargin += 0.03

        # Streamlines plotted on the mirrored domain (show symmetric halves)
        streams=self.add_plot("liquid/velocity", mode="streamlines", linewidths=8, linecolor="white", transform=[PlotTransformMirror(x=False, offset_x=offset_x,offset_y=0),PlotTransformMirror(x=True, offset_x=-offset_x,offset_y=0)])
        # Also plot velocity magnitude with the same mirror transform so colorbar matches
        self.add_plot("liquid/velocity", colorbar=cb_v, transform=PlotTransformMirror(x=False, offset_x=offset_x, offset_y=0))
        # Tweak streamline arrow style and density for clarity
        for stream in streams:
            stream.arrowstyle = ArrowStyle("-|>", head_length=3.5, head_width=3.3)
            stream.density *= 0.5
            if float(p.theta.value) < 20:
                stream.density *= 0.7

        # If not using the flat rivulet limit, plot the fas phase water vapour concentration
        if not p.use_flat_rivulet_limit_for_evap:
            cb_gas=self.add_colorbar("c", position="top left", cmap="Blues", length=0.35, thickness=0.06)
            cb_gas.ymargin += 0.03
            cb_gas.xmargin += 0.05
            cb_gas.ticsize = 80
            cb_gas.textsize = 80
            self.add_plot("gas/cvap", colorbar=cb_gas, transform=[PlotTransformMirror(x=False, offset_x=offset_x,offset_y=0),PlotTransformMirror(x=True, offset_x=-offset_x,offset_y=0)])

        # --- Arrow key (legend) for evaporation rate vectors ---
        arrk = self.add_arrow_key(title="evap. rate", position="top right")
        arrk.maxlength=0.4
        arrk.minlength_key=0.2
        arrk.maxlength_key=0.4
        arrk.arrowprops=dict(color="lawngreen", fc="lawngreen", ec="lawngreen", arrowstyle="simple",linewidth=4,mutation_scale=80)

        # Nudges for the arrow key placement
        arrk.ymargin+=0.05

        # Add evaporation arrows on both mirrored sides of the droplet interface
        arrows=self.add_plot("liquid/liquid_gas/J", arrowkey=arrk, arrowdensity=30, transform=[PlotTransformMirror(x=False, offset_x=offset_x,offset_y=0),PlotTransformMirror(x=True, offset_x=-offset_x,offset_y=0)])
        
        # Ensure arrows render on top
        for arrow in arrows:
            arrow.zindex = 100

        # --- Volume / parameter label ---
        # Format Marangoni number for display (special-case zero)
        text= self.add_text("Ma={}\nb={}\n$\\theta={:.1f}^\circ$".format(p.Ma.value,float(p.b.value),float(p.theta.value/degree)), "top center", bbox=dict(facecolor='wheat', boxstyle='round', ec="black", lw=2))
        text.ymargin+=0.2

        # --- Background box for labels ---
        ybox=self.ymax-0.2*(self.ymax-self.ymin)
        self.add_polygon([(self.xmin,self.ymax),(self.xmax,self.ymax),(self.xmax,ybox),(self.xmin,ybox)],edgecolor="black",facecolor="lightgrey",alpha=0.75).zindex=5

        # --- Text and element sizes for publication-quality figures ---
        cb_T.textsize = 80
        cb_T.ticsize = 80
        cb_v.textsize = 80
        cb_v.ticsize = 80
        arrk.textsize = 80
        text.textsize = 80



class QuasiStationaryLubricationDropletPlotter(MatplotlibPlotter):
    def __init__(self, problem = None, filetrunk = "plot_{:05d}", fileext = "png", eigenvector = None, eigenmode = "abs", add_eigen_to_mesh_positions = True, position_eigen_scale = 1, eigenscale = 1):
        # Initialize base MatplotlibPlotter with optional eigenvector-related args
        super().__init__(problem, filetrunk, fileext, eigenvector, eigenmode, add_eigen_to_mesh_positions, position_eigen_scale, eigenscale)

    def define_plot(self):
        # --- General plotting defaults and fonts ---
        self.defaults("colorbar").hide_some_ticks=False
        self.useLaTeXFont()

        # Retrieve typed problem instance for convenience
        pr = cast(QuasiStationaryLubricationDroplet, self.get_problem())

        # Geometry / view parameters (base radius and liquid half-width b)
        Rf = pr.base_radius
        b = pr.b.value

        # Set view window centered on droplet with some vertical margin
        self.set_view(-1.1*b, -Rf*2, 1.1*b, Rf*2)

        # --- Flow (Q) magnitude colorbar and plot ---
        cb_Q = self.add_colorbar(r"$\mathbf{Q}$ magn.", cmap="seismic", position="top left")
        # Plot Q field on a mirrored transform so droplet sits centered
        self.add_plot("liquid/Q", colorbar=cb_Q, transform=PlotTransformMirror(x=False, y=False, offset_x=-b/2, offset_y=-0.1*Rf), levels=64)

        # --- Prepare data for quiver plots (vector fields) ---
        liquid_data = self._get_mesh_data("liquid")
        Mindicator = liquid_data.get_data("Mindicator")
        Q_x, Q_y = liquid_data.get_data("Q_x"), liquid_data.get_data("Q_y")

        # Coordinates of the mesh nodes used for triangulation/interpolation
        coords = liquid_data.get_coordinates()
        X, Y = coords[0, :], coords[1, :]

        # Create a regular grid for plotting arrows (coarse mesh for clarity)
        xl = numpy.linspace(-1, 1, 12, endpoint=True)
        yl = numpy.linspace(-1, 1, 12, endpoint=True)
        xs, ys = numpy.meshgrid(xl, yl)

        # Triangulate and linearly interpolate Q components and indicator onto regular grid
        triang = matplotlib.tri.Triangulation(X, Y)
        Qxinter = tri.LinearTriInterpolator(triang, Q_x)(xs, ys)
        Qyinter = tri.LinearTriInterpolator(triang, Q_y)(xs, ys)
        Minter = tri.LinearTriInterpolator(triang, Mindicator)(xs, ys)

        # Mask regions where Mindicator > 0 to only show relevant Q
        Qxconti = numpy.ma.masked_where(Minter > 0, Qxinter)
        Qyconti = numpy.ma.masked_where(Minter > 0, Qyinter)
        Qmagnitude = numpy.sqrt(Qxconti**2 + Qyconti**2)
        Qxconti = Qxconti / numpy.where(Qmagnitude / numpy.max(Qmagnitude) > 0.01, Qmagnitude, 1)
        Qyconti = Qyconti / numpy.where(Qmagnitude / numpy.max(Qmagnitude) > 0.01, Qmagnitude, 1)

        # Plot green arrows for the primary Q field on both symmetric halves (mirror y)
        plt.gca().quiver(xs - b/2, ys - 0.1*Rf, Qxconti, Qyconti, pivot="mid", width=0.005, fc="lightgreen", linewidth=1, ec="black", zorder=10, scale=30)
        plt.gca().quiver(xs - b/2, -ys - 0.1*Rf, Qxconti, -Qyconti, pivot="mid", width=0.005, fc="lightgreen", linewidth=1, ec="black", zorder=10, scale=30)

        # Now mask complementary region (Mindicator < 0) to plot a different arrow style/color (yellow)
        Qxmara = numpy.ma.masked_where(Minter < 0, Qxinter)
        Qymara = numpy.ma.masked_where(Minter < 0, Qyinter)
        Qmagnitude = numpy.sqrt(Qxmara**2 + Qymara**2)
        Qxmara = Qxmara / numpy.where(Qmagnitude / numpy.max(Qmagnitude) > 0.01, Qmagnitude, 1)
        Qymara = Qymara / numpy.where(Qmagnitude / numpy.max(Qmagnitude) > 0.01, Qmagnitude, 1)

        # Plot masked yellow arrows on both mirrored halves
        plt.gca().quiver(xs - b/2, ys - 0.1*Rf, Qxmara, Qymara, pivot="mid", width=0.005, fc="yellow", linewidth=1, ec="black", zorder=10, scale=30)
        plt.gca().quiver(xs - b/2, -ys - 0.1*Rf, Qxmara, -Qymara, pivot="mid", width=0.005, fc="yellow", linewidth=1, ec="black", zorder=10, scale=30)

        # --- Concentration xi colorbar and plot ---
        cb_c = self.add_colorbar("$\\xi$", cmap="coolwarm", position="bottom left")
        # Plot xi on mirrored vertical axis so top/bottom symmetry is shown
        self.add_plot("liquid/c", colorbar=cb_c, transform=PlotTransformMirror(x=False, y=True, offset_x=-b/2, offset_y=-0.1*Rf), levels=64)

        # --- Velocity magnitude colorbars (at two z-levels) ---
        udata = numpy.sqrt(liquid_data.get_data("u_x")**2 + liquid_data.get_data("u_y")**2)
        umin = numpy.amin(udata)
        umax = numpy.amax(udata)

        uhdata = numpy.sqrt(liquid_data.get_data("uh_x")**2 + liquid_data.get_data("uh_y")**2)
        uhmin = numpy.amin(uhdata)
        uhmax = numpy.amax(uhdata)

        # Colorbar for velocity at specified z (pr.velocity_height_for_output)
        cb_u = self.add_colorbar(r"$\|\mathbf{u}\|$ at $z=$" + str(round(float(pr.velocity_height_for_output), 4)), cmap="viridis", position="bottom right", norm=matplotlib.colors.LogNorm(vmin=umin, vmax=umax), vmin=umin, vmax=umax)
        # Colorbar for velocity at the free surface z=h
        cb_uh = self.add_colorbar(r"$\|\mathbf{u}\|$ at $z=h$", cmap="viridis", position="top right", norm=matplotlib.colors.LogNorm(vmin=uhmin, vmax=uhmax), vmin=uhmin, vmax=uhmax)

        # Plot velocity magnitude fields and arrow overlays using mirrored transforms so both halves are visible
        self.add_plot("liquid/u", colorbar=cb_u, transform=[PlotTransformMirror(x=True, y=True, offset_x=b/2, offset_y=-0.1*Rf)], levels=64)
        self.add_plot("liquid/uh", colorbar=cb_uh, transform=[PlotTransformMirror(x=True, y=False, offset_x=b/2, offset_y=-0.1*Rf)], levels=64)

        # Add vector arrows for u and uh overlaid (combined into a single list for styling)
        arrows = self.add_plot("liquid/u", mode="arrows", transform=[PlotTransformMirror(x=True, y=True, offset_x=b/2, offset_y=-0.1*Rf)], arrowdensity=20)
        arrows += self.add_plot("liquid/uh", mode="arrows", transform=[PlotTransformMirror(x=True, y=False, offset_x=b/2, offset_y=-0.1*Rf)], arrowdensity=20)

        # Tweak arrow styles for visibility
        for arrow in arrows:
            arrow.arrowstyle = ArrowStyle("->", head_length=1, head_width=0.7)
            arrow.linewidths = 4
            arrow.arrowlength = 0.15

        # --- Final colorbar/text sizing and layout tweaks ---
        for cbi in [cb_u, cb_uh, cb_c, cb_Q]:
            cbi.textsize = 40
            cbi.ticsize = 40
            cbi.length = 0.4 * Rf
            cbi.xmargin += 0.02 * Rf
            cbi.labelpad = 10 * Rf
            cbi.ymargin += 0.02 * Rf

        # --- Parameter text annotation (Ma, b, theta) ---
        text = self.add_text("Ma={}\nb={}\n$\\theta={:.1f}^\circ$".format(pr.Ma.value, float(pr.b.value), float(pr.theta.value / degree)), "top center", bbox=dict(facecolor='wheat', boxstyle='round', ec="black", lw=2))
        text.ymargin += 0.15 * Rf
        text.textsize = 25