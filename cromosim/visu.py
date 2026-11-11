"""visualization utilities for Cromosim microscopic simulations.

Two classes: 
- pyvistaVisualization
- matplotlibVisualization

Old style visualization functions:
- plot_people()
- plot_sensors() still used at the end of the simulation
    
The module contains a 2-D renderers able to display the simulation
layers : pedestrians, contacts, actual and desired velocities, 
trajectories, and sensor lines.
"""

# Authors:
#     Sylvain Faure <sylvain.faure@universite-paris-saclay.fr>
#     Bertrand Maury <bertrand.maury@universite-paris-saclay.fr>
# License: GPL

from __future__ import annotations

import logging

import numpy as np
from PIL import Image

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Ellipse, Rectangle, Polygon
from matplotlib.lines import Line2D
from matplotlib.collections import LineCollection, EllipseCollection

try:
    import pyvista as pv
except ImportError:  # PyVista is optional when Matplotlib is selected.
    pv = None

LOGGER = logging.getLogger(__name__)


class _VisualizationDataMixin:
    """Common data preparation helpers for visualization backends.

    This private mixin contains only backend-independent operations.  Rendering
    remains implemented separately by the Matplotlib and PyVista classes.
    """

    @staticmethod
    def _as_xyrv(people):
        """Return a normalized ``(N, 4)`` pedestrian array.

        ``None``, a missing ``"xyrv"`` key, ``"xyrv": None`` and empty
        arrays are all treated as an empty population.
        """
        if people is None:
            return np.empty((0, 4), dtype=float)

        values = people.get("xyrv") if isinstance(people, dict) else people
        if values is None:
            return np.empty((0, 4), dtype=float)

        xyrv = np.asarray(values, dtype=float)
        if xyrv.size == 0:
            return np.empty((0, 4), dtype=float)

        if xyrv.ndim == 1:
            if xyrv.shape[0] != 4:
                raise ValueError(
                    "A one-dimensional xyrv array must contain exactly 4 values"
                )
            xyrv = xyrv.reshape((1, 4))

        if xyrv.ndim != 2 or xyrv.shape[1] < 4:
            raise ValueError(
                "xyrv must be a two-dimensional array with at least 4 columns"
            )

        return xyrv

    @staticmethod
    def _combine_field(people, virtual_people, key, n_components=1):
        """Combine real and virtual pedestrian fields when available."""
        empty_shape = (0,) if n_components == 1 else (0, n_components)

        if not isinstance(people, dict):
            return np.empty(empty_shape, dtype=float)

        real_values = people.get(key)
        if real_values is None:
            real = np.empty(empty_shape, dtype=float)
        else:
            real = np.asarray(real_values, dtype=float)
            if real.size == 0:
                real = np.empty(empty_shape, dtype=float)

        if not isinstance(virtual_people, dict):
            return real

        virtual_values = virtual_people.get(key)
        if virtual_values is None:
            return real

        virtual = np.asarray(virtual_values, dtype=float)
        if virtual.size == 0:
            return real

        if real.size == 0:
            return virtual
        return np.concatenate((real, virtual), axis=0)

    @staticmethod
    def _contact_segments(xyrv, contacts):
        """Return person-person and person-wall contact segments."""
        if contacts is None or len(contacts) == 0:
            return []

        xyrv = np.asarray(xyrv, dtype=float)
        contacts = np.asarray(contacts)
        if xyrv.size == 0 or contacts.size == 0:
            return []

        segments = []
        for contact in contacts:
            i = int(contact[0])
            j = int(contact[1])
            if i < 0 or i >= xyrv.shape[0]:
                continue

            p0 = xyrv[i, :2]
            if j >= 0:
                if j >= xyrv.shape[0]:
                    continue
                p1 = xyrv[j, :2]
            else:
                distance = contact[2]
                normal = contact[3:5]
                p1 = p0 - (xyrv[i, 2] + distance) * normal
            segments.append(np.vstack((p0, p1)))
        return segments

    @staticmethod
    def _path_segments(people):
        """Return one polyline array for each stored pedestrian path."""
        if not isinstance(people, dict):
            return []

        segments = []
        for path in people.get("paths", {}).values():
            coordinates = np.asarray(path, dtype=float)
            if coordinates.size < 4:
                continue
            segments.append(coordinates.reshape((-1, 2)))
        return segments

    @staticmethod
    def _sensor_segments(sensors):
        """Return sensor line segments."""
        segments = []
        for sensor in sensors or []:
            line = np.asarray(sensor.get("line", []), dtype=float)
            if line.size == 4:
                segments.append(line.reshape((2, 2)))
        return segments

    @staticmethod
    def _compute_velocity_scale(velocity_reference, velocity_arrow_length):
        """Return the fixed conversion from velocity to displayed arrow length.

        ``velocity_reference`` is expressed in simulation velocity units
        (typically m/s), while ``velocity_arrow_length`` is expressed in the
        physical coordinates of the domain (typically m).  Consequently, a
        vector whose norm is ``velocity_reference`` is always drawn with
        length ``velocity_arrow_length``, independently of the time step.
        """
        velocity_reference = float(velocity_reference)
        velocity_arrow_length = float(velocity_arrow_length)
        if velocity_reference <= 0:
            raise ValueError("velocity_reference must be strictly positive")
        if velocity_arrow_length <= 0:
            raise ValueError("velocity_arrow_length must be strictly positive")
        return velocity_arrow_length / velocity_reference


class matplotlibVisualization(_VisualizationDataMixin):
    """2-D Matplotlib renderer for a Cromosim domain.

    Parameters
    ----------
    domain : cromosim.domain.Domain
        Computational domain to display.
    figure_id : int, optional
        Matplotlib figure identifier.
    axis : sequence of float, optional
        Display limits ``[xmin, xmax, ymin, ymax]``.  By default the domain
        limits are used.
    cmap : str, default="winter"
        Colormap used for real pedestrians.
    figsize : tuple of float, optional
        Matplotlib figure size in inches.
    velocity_reference : float, default=1.5
        Reference velocity used to define a fixed arrow scale. A vector with
        this norm is displayed with length ``velocity_arrow_length``.
    velocity_arrow_length : float, default=0.5
        Displayed length, in domain coordinates, of a velocity vector whose
        norm is ``velocity_reference``. The same scale is used for ``U`` and
        ``Vd`` and remains constant throughout the simulation.

    Notes
    -----
    The background image and axes are created only once.  Pedestrian positions,
    contacts, trajectories, velocity vectors and sensor lines are subsequently
    updated in place whenever possible.
    """

    def __init__(
        self,
        domain,
        figure_id=None,
        axis=None,
        cmap="winter",
        figsize=None,
        velocity_reference=1.5,
        velocity_arrow_length=0.5,
    ):
        self.domain = domain
        self.cmap = cmap
        self.axis = axis
        self.velocity_reference = float(velocity_reference)
        self.velocity_arrow_length = float(velocity_arrow_length)
        self.velocity_scale = self._compute_velocity_scale(
            self.velocity_reference, self.velocity_arrow_length
        )

        self.fig = plt.figure(num=figure_id, figsize=figsize)
        self.fig.clf()
        self.ax = self.fig.add_subplot(111)

        self.background = self.ax.imshow(
            domain.image,
            interpolation="nearest",
            extent=[domain.xmin, domain.xmax, domain.ymin, domain.ymax],
            origin="lower",
            zorder=0,
        )

        if axis is None:
            self.ax.set_xlim(domain.xmin, domain.xmax)
            self.ax.set_ylim(domain.ymin, domain.ymax)
        else:
            self.ax.set_xlim(axis[0], axis[1])
            self.ax.set_ylim(axis[2], axis[3])

        self.ax.set_xticks([])
        self.ax.set_yticks([])
        self.ax.set_axis_off()
        self.time_text = self.ax.set_title("")

        self.people_collection = None
        self.virtual_people_collection = None
        self.contacts_collection = None
        self.paths_collection = None
        self.sensors_collection = None
        self.velocity_quiver = None
        self.desired_velocity_quiver = None
        self.colorbar = None

        self._people_radii = None
        self._virtual_people_radii = None
        self._velocity_count = 0
        self._desired_velocity_count = 0

        self.fig.canvas.draw_idle()

    @staticmethod
    def _remove_artist(artist):
        """Remove a Matplotlib artist if it currently exists."""
        if artist is not None:
            artist.remove()
        return None

    def _create_people_collection(
        self,
        xyrv,
        *,
        colors=None,
        cmap=None,
        vmin=None,
        vmax=None,
        virtual=False,
    ):
        """Create a pedestrian ellipse collection."""
        xyrv = np.asarray(xyrv, dtype=float)
        if xyrv.size == 0:
            return None

        kwargs = {
            "widths": 2.0 * xyrv[:, 2],
            "heights": 2.0 * xyrv[:, 2],
            "angles": 0,
            "units": "xy",
            "offsets": xyrv[:, :2],
            "transOffset": self.ax.transData,
            "zorder": 3 if not virtual else 3.1,
        }

        if virtual:
            kwargs.update(facecolors="black", edgecolors="none", alpha=0.5)
        else:
            kwargs.update(cmap=plt.get_cmap(cmap or self.cmap))

        collection = EllipseCollection(**kwargs)

        if not virtual:
            if colors is None:
                colors = xyrv[:, 2]
            collection.set_array(np.asarray(colors))
            if vmin is not None and vmax is not None:
                collection.set_clim(vmin, vmax)

        self.ax.add_collection(collection)
        return collection

    def _update_people(
        self,
        people,
        virtual_people,
        *,
        enabled,
        color_scalar_values,
        color_minmax,
        color_cmap,
        show_scalar_bar,
        color_scalar_name,
    ):
        """Update pedestrian ellipse collections."""
        if not enabled:
            # The colorbar must be removed before its mappable collection.
            # Otherwise Matplotlib may find ``mappable.axes is None`` while
            # restoring the subplot geometry in ``Colorbar.remove()``.
            self._remove_colorbar()
            self.people_collection = self._remove_artist(self.people_collection)
            self.virtual_people_collection = self._remove_artist(
                self.virtual_people_collection
            )
            self._people_radii = None
            self._virtual_people_radii = None
            return

        xyrv = self._as_xyrv(people)
        virtual_xyrv = self._as_xyrv(virtual_people)

        if color_scalar_values is None:
            colors = xyrv[:, 2] if xyrv.size else np.empty((0,), dtype=float)
        else:
            colors = np.asarray(color_scalar_values)
            if colors.shape[0] != xyrv.shape[0]:
                raise ValueError(
                    "color_scalar_values must contain one value per real pedestrian"
                )

        radii = xyrv[:, 2].copy() if xyrv.size else np.empty((0,))
        geometry_changed = (
            self.people_collection is None
            or self._people_radii is None
            or radii.shape != self._people_radii.shape
            or not np.array_equal(radii, self._people_radii)
        )

        if geometry_changed:
            # A colorbar references ``self.people_collection`` as its mappable.
            # Remove the colorbar first whenever the collection is replaced or
            # removed (for instance when the population becomes empty).
            if self.colorbar is not None:
                self._remove_colorbar()
            self.people_collection = self._remove_artist(self.people_collection)
            self.people_collection = self._create_people_collection(
                xyrv,
                colors=colors,
                cmap=color_cmap,
                vmin=None if color_minmax is None else color_minmax[0],
                vmax=None if color_minmax is None else color_minmax[1],
            )
            self._people_radii = radii
        elif self.people_collection is not None:
            self.people_collection.set_offsets(xyrv[:, :2])
            self.people_collection.set_array(colors)
            if color_minmax is not None:
                self.people_collection.set_clim(*color_minmax)

        virtual_radii = (
            virtual_xyrv[:, 2].copy()
            if virtual_xyrv.size
            else np.empty((0,), dtype=float)
        )
        virtual_geometry_changed = (
            self.virtual_people_collection is None
            or self._virtual_people_radii is None
            or virtual_radii.shape != self._virtual_people_radii.shape
            or not np.array_equal(virtual_radii, self._virtual_people_radii)
        )

        if virtual_geometry_changed:
            self.virtual_people_collection = self._remove_artist(
                self.virtual_people_collection
            )
            self.virtual_people_collection = self._create_people_collection(
                virtual_xyrv,
                virtual=True,
            )
            self._virtual_people_radii = virtual_radii
        elif self.virtual_people_collection is not None:
            self.virtual_people_collection.set_offsets(virtual_xyrv[:, :2])

        self._update_colorbar(show_scalar_bar, color_scalar_name)

    def _remove_colorbar(self):
        """Remove the pedestrian colorbar safely, if present.

        Matplotlib's ``Colorbar.remove()`` expects the mappable to still be
        attached to an Axes.  Normally this method is called before removing
        the pedestrian collection.  The fallback below also handles a stale
        colorbar left behind by an earlier update.
        """
        if self.colorbar is None:
            return

        colorbar = self.colorbar
        self.colorbar = None

        try:
            colorbar.remove()
        except (AttributeError, ValueError):
            # Fallback for a colorbar whose mappable has already been removed.
            colorbar_ax = getattr(colorbar, "ax", None)
            if colorbar_ax is not None:
                try:
                    colorbar_ax.remove()
                except (AttributeError, KeyError, ValueError):
                    pass
            self.fig.canvas.draw_idle()

    def _update_colorbar(self, show_scalar_bar, scalar_name):
        """Create or update the pedestrian scalar bar."""
        if not show_scalar_bar or self.people_collection is None:
            self._remove_colorbar()
            return

        if self.colorbar is None:
            self.colorbar = self.fig.colorbar(
                self.people_collection,
                ax=self.ax,
                orientation="horizontal",
                fraction=0.04,
                pad=0.02,
            )
        else:
            self.colorbar.update_normal(self.people_collection)
        self.colorbar.set_label(scalar_name)

    def _update_line_collection(
        self,
        artist,
        segments,
        *,
        enabled,
        color,
        linewidth,
        alpha,
        zorder,
    ):
        """Create or update a line collection."""
        if not enabled:
            return self._remove_artist(artist)

        if artist is None:
            artist = LineCollection(
                segments,
                colors=color,
                linewidths=linewidth,
                alpha=alpha,
                zorder=zorder,
            )
            self.ax.add_collection(artist)
        else:
            artist.set_segments(segments)
            artist.set_color(color)
            artist.set_linewidth(linewidth)
            artist.set_alpha(alpha)
        return artist

    def _update_quiver(
        self,
        artist,
        centers,
        vectors,
        *,
        enabled,
        color,
        zorder,
        count_attribute,
    ):
        """Create or update a velocity quiver."""
        if not enabled:
            setattr(self, count_attribute, 0)
            return self._remove_artist(artist)

        centers = np.asarray(centers, dtype=float)
        vectors = np.asarray(vectors, dtype=float)
        count = min(len(centers), len(vectors))
        centers = centers[:count]
        vectors = vectors[:count]
        scaled_vectors = self.velocity_scale * vectors

        previous_count = getattr(self, count_attribute)
        if artist is None or previous_count != count:
            artist = self._remove_artist(artist)
            if count:
                artist = self.ax.quiver(
                    centers[:, 0],
                    centers[:, 1],
                    scaled_vectors[:, 0],
                    scaled_vectors[:, 1],
                    angles="xy",
                    scale_units="xy",
                    scale=1.0,
                    color=color,
                    width=0.003,
                    zorder=zorder,
                )
        elif count:
            artist.set_offsets(centers)
            artist.set_UVC(scaled_vectors[:, 0], scaled_vectors[:, 1])

        setattr(self, count_attribute, count)
        return artist

    def update(
        self,
        *,
        people,
        contacts=[],
        virtual_people=None,
        sensors=None,
        time=None,
        plot_people=True,
        plot_contacts=True,
        plot_velocities=False,
        plot_desired_velocities=False,
        plot_paths=False,
        plot_sensors=False,
        color_scalar_values=None,
        color_scalar_name="r",
        color_minmax=None,
        color_cmap="winter",
        show_scalar_bar=False,
        contacts_line_width=1.0,
        paths_line_width=0.5,
        sensors_line_width=1.0,
        savefig=False,
        filename=None,
        dpi=150,
    ):
        """Update one Matplotlib frame without rebuilding the figure.

        Parameters mirror the PyVista visualization layer flags.  Only artists
        whose topology changes are recreated; otherwise their underlying data
        are updated in place.
        """
        real_xyrv = self._as_xyrv(people)
        virtual_xyrv = self._as_xyrv(virtual_people)

        self._update_people(
            people,
            virtual_people,
            enabled=plot_people,
            color_scalar_values=color_scalar_values,
            color_minmax=color_minmax,
            color_cmap=color_cmap,
            show_scalar_bar=show_scalar_bar,
            color_scalar_name=color_scalar_name,
        )

        if virtual_xyrv.size:
            combined_xyrv = np.concatenate((real_xyrv, virtual_xyrv), axis=0)
        else:
            combined_xyrv = real_xyrv

        self.contacts_collection = self._update_line_collection(
            self.contacts_collection,
            self._contact_segments(combined_xyrv, contacts),
            enabled=plot_contacts,
            color="black",
            linewidth=contacts_line_width,
            alpha=0.6,
            zorder=4,
        )

        self.paths_collection = self._update_line_collection(
            self.paths_collection,
            self._path_segments(people),
            enabled=plot_paths,
            color="black",
            linewidth=paths_line_width,
            alpha=0.5,
            zorder=2,
        )

        self.sensors_collection = self._update_line_collection(
            self.sensors_collection,
            self._sensor_segments(sensors),
            enabled=plot_sensors,
            color="green",
            linewidth=sensors_line_width,
            alpha=0.8,
            zorder=6,
        )

        centers = combined_xyrv[:, :2] if combined_xyrv.size else np.empty((0, 2))

        velocities = self._combine_field(people, virtual_people, "U", n_components=2)
        self.velocity_quiver = self._update_quiver(
            self.velocity_quiver,
            centers,
            velocities,
            enabled=plot_velocities,
            color="red",
            zorder=5,
            count_attribute="_velocity_count",
        )

        desired_velocities = self._combine_field(
            people, virtual_people, "Vd", n_components=2
        )
        self.desired_velocity_quiver = self._update_quiver(
            self.desired_velocity_quiver,
            centers,
            desired_velocities,
            enabled=plot_desired_velocities,
            color="blue",
            zorder=5,
            count_attribute="_desired_velocity_count",
        )

        self.time_text.set_text("" if time is None else f"time = {time:.2f} s, domain = {self.domain.name}")

        self.fig.canvas.draw_idle()
        self.fig.canvas.flush_events()

        if savefig:
            if filename is None:
                raise ValueError("filename must be provided when savefig=True")
            self.fig.savefig(
                filename,
                dpi=dpi,
                bbox_inches="tight",
                pad_inches=0,
            )

    def show(self):
        """Show the figure without blocking the simulation loop."""
        self.fig.show()
        self.fig.canvas.draw_idle()
        self.fig.canvas.flush_events()

    def process_events(self):
        """Process GUI events without recreating the figure."""
        self.fig.canvas.flush_events()

    def close(self):
        """Close the Matplotlib figure."""
        plt.close(self.fig)


class pyvistaVisualization(_VisualizationDataMixin):
    """Display a Cromosim microscopic simulation with PyVista.

    Parameters
    ----------
    image_path : str or path-like
        Path to the background image representing the domain.
    pixel_size : float
        Physical size represented by one image pixel.
    movie_path : str or path-like, optional
        Output movie path. When provided, rendered frames can be written with
        :meth:`update` by setting ``write_frame=True``.
    fps : int, default=20
        Movie frame rate.
    window_size : tuple of int, default=(1000, 800)
        PyVista window size in pixels.
    velocity_reference : float, default=1.5
        Reference velocity used to define a fixed arrow scale. A vector with
        this norm is displayed with length ``velocity_arrow_length``.
    velocity_arrow_length : float, default=0.5
        Displayed length, in domain coordinates, of a velocity vector whose
        norm is ``velocity_reference``. The same scale is used for ``U`` and
        ``Vd`` and remains constant throughout the simulation.

    Notes
    -----
    The public class name is intentionally kept in French for backward
    compatibility with existing Cromosim scripts.
    """

    _LAYER_BACKGROUND = 0.000
    _LAYER_PATHS      = 0.010
    _LAYER_PEOPLE     = 0.020
    # Contacts must be above pedestrian disks. Since person-person contact
    # segments run from center to center, drawing them below the disks makes
    # almost the entire segment invisible.
    _LAYER_CONTACTS   = 0.030
    _LAYER_VECTORS    = 0.040
    _LAYER_SENSORS    = 0.050

    def __init__(
        self,
        domain,
        image_path,
        pixel_size,
        movie_path=None,
        fps=20,
        window_size=(1000, 800),
        movie_macro_block_size=1,
        velocity_reference=1.5,
        velocity_arrow_length=0.5,
    ):
        if pv is None:
            raise ImportError(
                "PyVista visualization was requested, but the 'pyvista' "
                "package is not installed. Install PyVista or select the "
                "Matplotlib visualization backend."
            )

        self.domain = domain
        self.pixel_size = float(pixel_size)
        self.movie_path = movie_path
        self.movie_macro_block_size = movie_macro_block_size
        self.velocity_reference = float(velocity_reference)
        self.velocity_arrow_length = float(velocity_arrow_length)
        self.velocity_scale = self._compute_velocity_scale(
            self.velocity_reference, self.velocity_arrow_length
        )

        image = Image.open(image_path).convert("RGB")
        self.image_array = np.asarray(image)
        self.height, self.width = self.image_array.shape[:2]
        self.W = self.width * self.pixel_size
        self.H = self.height * self.pixel_size

        self.plotter = pv.Plotter(window_size=window_size)
        
        # if self.plotter.render_window is not None:
        #     self.plotter.render_window.SetCoverable(False)
        
        self._add_background_plane()
        self._configure_camera()

        if movie_path is not None:
            self.plotter.open_movie(
                movie_path, 
                framerate=fps,
                macro_block_size=movie_macro_block_size)

    def _add_background_plane(self):
        """Create the textured plane used as the domain background."""
        plane = pv.Plane(
            center=(self.W / 2, self.H / 2, self._LAYER_BACKGROUND),
            direction=(0, 0, 1),
            i_size=self.W,
            j_size=self.H,
            i_resolution=1,
            j_resolution=1,
        )
        plane.texture_map_to_plane(inplace=True)
        texture = pv.Texture(self.image_array)
        self.plotter.add_mesh(
            plane,
            texture=texture,
            lighting=False,
            name="background",
            render=False,
        )

    def _configure_camera(self):
        """Configure a fixed orthographic camera for an XY-plane view."""
        self.plotter.set_background("white")
        self.plotter.view_xy()
        self.plotter.camera.parallel_projection = True
        self.plotter.enable_2d_style()
        self.plotter.camera.focal_point = (self.W / 2, self.H / 2, 0)
        self.plotter.camera.position = (self.W / 2, self.H / 2, 100)
        self.plotter.camera.up = (0, 1, 0)

    def show(self):
        """Open the interactive PyVista window without blocking execution."""
        self.plotter.show(interactive_update=True, auto_close=False)

    def process_events(self):
        """Process pending events for the interactive PyVista window."""
        self.plotter.update()

    def _remove_actor(self, name):
        """Remove a named actor without triggering an intermediate render."""
        self.plotter.remove_actor(name, reset_camera=False, render=False)

    def _set_mesh(self, name, mesh, **kwargs):
        """Add or replace a named mesh actor without resetting the camera."""
        if mesh is None or mesh.n_points == 0:
            self._remove_actor(name)
            return

        self.plotter.add_mesh(
            mesh,
            name=name,
            reset_camera=False,
            render=False,
            **kwargs,
        )

    @classmethod
    def _build_particle_glyphs(
        cls,
        xyrv,
        color_scalar_values=None,
        color_scalar_name="r",
        z_layer=None,
    ):
        """Build one disk glyph per pedestrian."""
        xyrv = np.asarray(xyrv)
        if xyrv.size == 0:
            return pv.PolyData(np.empty((0, 3), dtype=float))

        if z_layer is None:
            z_layer = cls._LAYER_PEOPLE

        centers = np.zeros((len(xyrv), 3), dtype=float)
        centers[:, :2] = xyrv[:, :2]
        centers[:, 2] = z_layer

        points = pv.PolyData(centers)
        points["r"] = xyrv[:, 2]

        use_custom_scalar = (
            color_scalar_values is not None and color_scalar_name != "r"
        )
        if use_custom_scalar:
            color_scalar_values = np.asarray(color_scalar_values)
            if color_scalar_values.shape[0] != xyrv.shape[0]:
                raise ValueError(
                    "color_scalar_values must contain one value per pedestrian"
                )
            points[color_scalar_name] = color_scalar_values

        disk = pv.Disc(
            center=(0, 0, 0),
            inner=0,
            outer=1,
            r_res=1,
            c_res=24,
        )
        glyphs = points.glyph(geom=disk, scale="r", orient=False)

        if use_custom_scalar:
            glyphs[color_scalar_name] = np.repeat(
                points[color_scalar_name], disk.n_points
            )

        return glyphs

    @staticmethod
    def _build_segments_mesh(segments, z_layer):
        """Build a PolyData mesh from independent 2-D line segments."""
        if not segments:
            return pv.PolyData(np.empty((0, 3), dtype=float))

        points = np.empty((2 * len(segments), 3), dtype=float)
        lines = np.empty((len(segments), 3), dtype=np.int64)

        for index, (p0, p1) in enumerate(segments):
            points[2 * index, :] = (p0[0], p0[1], z_layer)
            points[2 * index + 1, :] = (p1[0], p1[1], z_layer)
            lines[index, :] = (2, 2 * index, 2 * index + 1)

        mesh = pv.PolyData(points)
        mesh.lines = lines.ravel()
        return mesh

    @classmethod
    def _build_contact_mesh(cls, xyrv, contacts):
        """Build line segments representing person-person and wall contacts."""
        return cls._build_segments_mesh(
            cls._contact_segments(xyrv, contacts),
            cls._LAYER_CONTACTS,
        )

    @classmethod
    def _build_paths_mesh(cls, people):
        """Build polyline trajectories stored in ``people['paths']``."""
        paths = cls._path_segments(people)
        if not paths:
            return pv.PolyData(np.empty((0, 3), dtype=float))

        points_blocks = []
        lines = []
        offset = 0

        for coordinates in paths:
            xyz = np.zeros((coordinates.shape[0], 3), dtype=float)
            xyz[:, :2] = coordinates
            xyz[:, 2] = cls._LAYER_PATHS
            points_blocks.append(xyz)
            lines.append(
                np.concatenate(
                    ([coordinates.shape[0]], np.arange(offset, offset + coordinates.shape[0]))
                )
            )
            offset += coordinates.shape[0]

        mesh = pv.PolyData(np.vstack(points_blocks))
        mesh.lines = np.concatenate(lines).astype(np.int64)
        return mesh

    @classmethod
    def _build_sensor_mesh(cls, sensors):
        """Build line segments representing sensor counting lines."""
        return cls._build_segments_mesh(
            cls._sensor_segments(sensors),
            cls._LAYER_SENSORS,
        )

    @classmethod
    def _build_vector_glyphs(cls, centers_xy, vectors_xy, scale=1.0):
        """Build arrow glyphs for a 2-D vector field attached to pedestrians."""
        centers_xy = np.asarray(centers_xy, dtype=float)
        vectors_xy = np.asarray(vectors_xy, dtype=float)
        if centers_xy.size == 0 or vectors_xy.size == 0:
            return pv.PolyData(np.empty((0, 3), dtype=float))

        count = min(centers_xy.shape[0], vectors_xy.shape[0])
        centers_xy = centers_xy[:count]
        vectors_xy = vectors_xy[:count]
        magnitude = np.linalg.norm(vectors_xy, axis=1)
        active = magnitude > 0
        if not np.any(active):
            return pv.PolyData(np.empty((0, 3), dtype=float))

        centers = np.zeros((active.sum(), 3), dtype=float)
        centers[:, :2] = centers_xy[active]
        centers[:, 2] = cls._LAYER_VECTORS

        vectors = np.zeros((active.sum(), 3), dtype=float)
        vectors[:, :2] = vectors_xy[active]

        points = pv.PolyData(centers)
        points["vectors"] = vectors
        points["magnitude"] = magnitude[active]

        arrow = pv.Arrow(
            start=(0, 0, 0),
            direction=(1, 0, 0),
            tip_length=0.25,
            tip_radius=0.08,
            shaft_radius=0.025,
            tip_resolution=8,
            shaft_resolution=8,
        )
        return points.glyph(
            orient="vectors",
            scale="magnitude",
            factor=float(scale),
            geom=arrow,
        )

    def update(
        self,
        xyrv=None,
        *,
        people=None,
        contacts=[],
        virtual_people=None,
        sensors=None,
        time=None,
        plot_people=True,
        plot_contacts=True,
        plot_velocities=False,
        plot_desired_velocities=False,
        plot_paths=False,
        plot_sensors=False,
        color_scalar_values=None,
        color_scalar_name="r",
        color_minmax=None,
        color_cmap="viridis",
        show_scalar_bar=True,
        write_frame=True,
    ):
        """Update all visible layers of the PyVista scene.

        Parameters
        ----------
        xyrv : array-like, optional
            Backward-compatible particle array of shape ``(N, 4)``. Prefer
            passing ``people`` when velocities or paths must be displayed.
        people : dict, optional
            Cromosim people dictionary. Expected fields include ``"xyrv"``,
            ``"U"``, ``"Vd"`` and ``"paths"`` as required by the selected
            plotting options.
        contacts : numpy.ndarray, optional
            Contact array ``i, j, dij, eij_x, eij_y``.
        virtual_people : dict, optional
            Virtual population used at domain interfaces.
        sensors : sequence of dict, optional
            Sensor definitions containing a ``"line"`` entry.
        time : float, optional
            Simulation time displayed in the upper-left corner.
        plot_people, plot_contacts, plot_velocities, plot_desired_velocities,
        plot_paths, plot_sensors : bool
            Enable the corresponding simulation layer. These flags mirror the
            Matplotlib options used by :func:`cromosim.micro.plot_people`.
        color_scalar_values : array-like, optional
            Per-person scalar values used to color real pedestrians.
        color_scalar_name : str, default="r"
            Name of the scalar field used for coloring.
        color_minmax : tuple of float, optional
            Minimum and maximum values of the color scale.
        color_cmap : str, default="viridis"
            Colormap name.
        show_scalar_bar : bool, default=True
            Display a scalar bar for real pedestrians.
        write_frame : bool, default=True
            Write the rendered frame when movie recording is enabled.
        """
        if people is None:
            people = {"xyrv": self._as_xyrv(xyrv)}
        real_xyrv = self._as_xyrv(people)
          
        virtual_xyrv = self._as_xyrv(virtual_people)

        if plot_people:
            real_mesh = self._build_particle_glyphs(
                real_xyrv,
                color_scalar_values=color_scalar_values,
                color_scalar_name=color_scalar_name,
            )
            self._set_mesh(
                "people",
                real_mesh,
                scalars=color_scalar_name,
                cmap=color_cmap,
                clim=color_minmax,
                show_scalar_bar=show_scalar_bar,
                lighting=False,
                scalar_bar_args={
                    "title": color_scalar_name,
                    "position_x": 0.3,
                    "position_y": 0.05,
                    "width": 0.4,
                    "height": 0.08,
                    "vertical": False,
                },
            )
            
            virtual_mesh = self._build_particle_glyphs(
                virtual_xyrv,
                z_layer=self._LAYER_PEOPLE + 0.001,
            )
            self._set_mesh(
                "virtual_people",
                virtual_mesh,
                color="black",
                opacity=0.5,
                show_scalar_bar=False,
                lighting=False,
            )
        else:
            self._remove_actor("people")
            self._remove_actor("virtual_people")

        combined_xyrv = real_xyrv
        if virtual_xyrv.size:
            combined_xyrv = np.concatenate((real_xyrv, virtual_xyrv), axis=0)

        if plot_contacts:
            self._set_mesh(
                "contacts",
                self._build_contact_mesh(combined_xyrv, contacts),
                color="black",
                line_width=2,
                show_scalar_bar=False,
                lighting=False,
            )
        else:
            self._remove_actor("contacts")

        if plot_paths:
            self._set_mesh(
                "paths",
                self._build_paths_mesh(people),
                color="gray",
                opacity=0.2,
                line_width=1,
                show_scalar_bar=False,
                lighting=False,
            )
        else:
            self._remove_actor("paths")

        if plot_sensors:
            self._set_mesh(
                "sensors",
                self._build_sensor_mesh(sensors),
                color="green",
                line_width=2.0,
                show_scalar_bar=False,
                lighting=False,
            )
        else:
            self._remove_actor("sensors")

        combined_centers = combined_xyrv[:, :2] if combined_xyrv.size else np.empty((0, 2))

        if plot_velocities:
            velocities = self._combine_field(people, virtual_people, "U", n_components=2)
            self._set_mesh(
                "velocities",
                self._build_vector_glyphs(
                    combined_centers, velocities, scale=self.velocity_scale
                ),
                color="red",
                show_scalar_bar=False,
                lighting=False,
            )
        else:
            self._remove_actor("velocities")

        if plot_desired_velocities:
            desired_velocities = self._combine_field(
                people, virtual_people, "Vd", n_components=2
            )
            self._set_mesh(
                "desired_velocities",
                self._build_vector_glyphs(
                    combined_centers, desired_velocities, scale=self.velocity_scale
                ),
                color="blue",
                show_scalar_bar=False,
                lighting=False,
            )
        else:
            self._remove_actor("desired_velocities")

        if time is not None:
            self.plotter.add_text(
                f"time = {time:.2f} s, domain = {self.domain.name}",
                position="upper_left",
                font_size=10,
                name="simulation_time",
                render=False,
            )

        self.plotter.render()
        if write_frame and self.movie_path is not None:
            self.plotter.write_frame()

    def close(self):
        """Close the PyVista plotter and its movie writer, if any."""
        self.plotter.close()



def plot_people(ifig, dom, people, contacts, colors, time=-1, axis=None,
                virtual_people=None, plot_people=True, plot_contacts=True,
                plot_velocities=False, plot_desired_velocities=False,
                plot_paths=False, plot_sensors=False, vmin=None, vmax=None,
                sensors=None, savefig=False, filename="fig.png", dpi=150,
                cmap="winter"):
    """Plot pedestrians, contacts, velocities, paths and sensors.

    Parameters
    ----------
    ifig : int
        Matplotlib figure number.
    dom : Domain
        Domain whose background image and geometry are displayed.
    people : dict
        Pedestrian data.
    contacts : numpy.ndarray
        Contact rows ``(i, j, dij, eij_x, eij_y)``.
    colors : numpy.ndarray
        Scalar values used to color pedestrians.
    time : float, optional
        Simulation time displayed in the figure title when non-negative.
    axis : sequence of float or None, optional
        Plot limits ``[xmin, xmax, ymin, ymax]``.
    virtual_people : dict or None, optional
        Optional virtual pedestrians drawn with partial transparency.
    plot_people : bool, optional
        Draw pedestrians.
    plot_contacts : bool, optional
        Draw active contacts.
    plot_velocities : bool, optional
        Draw actual velocities.
    plot_desired_velocities : bool, optional
        Draw desired velocities.
    plot_paths : bool, optional
        Draw stored pedestrian trajectories.
    plot_sensors : bool, optional
        Draw sensor line segments.
    vmin, vmax : float or None, optional
        Color limits for pedestrian scalar values.
    sensors : sequence of dict or None, optional
        Sensor definitions.
    savefig : bool, optional
        Save the figure when ``True``.
    filename : str, optional
        Output image filename.
    dpi : int, optional
        Output resolution in dots per inch.
    cmap : str, optional
        Matplotlib colormap name.

    Returns
    -------
    None
        The function updates the selected Matplotlib figure in place.
    """
    if sensors is None:
        sensors = []

    fig = plt.figure(ifig)
    plt.clf()
    ax = fig.add_subplot(111)
    ax.imshow(
        dom.image,
        interpolation="nearest",
        extent=[dom.xmin, dom.xmax, dom.ymin, dom.ymax],
        origin="lower",
    )

    has_virtual_people = (
        virtual_people is not None
        and virtual_people.get("xyrv") is not None
    )

    if plot_people:
        offsets = people["xyrv"][:, :2]
        collection = EllipseCollection(
            widths=2 * people["xyrv"][:, 2],
            heights=2 * people["xyrv"][:, 2],
            angles=0,
            units="xy",
            cmap=plt.get_cmap(cmap),
            offsets=offsets,
            transOffset=ax.transData,
        )
        collection.set_array(colors)
        if vmin is not None and vmax is not None:
            collection.set_clim([vmin, vmax])
        ax.add_collection(collection)

        if has_virtual_people:
            virtual_offsets = virtual_people["xyrv"][:, :2]
            virtual_collection = EllipseCollection(
                widths=2 * virtual_people["xyrv"][:, 2],
                heights=2 * virtual_people["xyrv"][:, 2],
                angles=0,
                units="xy",
                facecolors="black",
                alpha=0.5,
                offsets=virtual_offsets,
                transOffset=ax.transData,
            )
            ax.add_collection(virtual_collection)

    if plot_contacts:
        if has_virtual_people:
            xyrv = np.concatenate((people["xyrv"], virtual_people["xyrv"]))
        else:
            xyrv = people["xyrv"]

        for contact in contacts:
            i = int(contact[0])
            j = int(contact[1])
            if j != -1:
                line = Line2D(
                    [xyrv[i, 0], xyrv[j, 0]],
                    [xyrv[i, 1], xyrv[j, 1]],
                    lw=1.0,
                    alpha=0.6,
                    color="k",
                )
            else:
                line = Line2D(
                    [xyrv[i, 0], xyrv[i, 0] - (xyrv[i, 2] + contact[2]) * contact[3]],
                    [xyrv[i, 1], xyrv[i, 1] - (xyrv[i, 2] + contact[2]) * contact[4]],
                    lw=1.0,
                    alpha=0.6,
                    color="k",
                )
            ax.add_line(line)

    if plot_velocities:
        for position, velocity in zip(people["xyrv"], people["U"]):
            ax.add_patch(
                Arrow(
                    position[0],
                    position[1],
                    velocity[0],
                    velocity[1],
                    width=0.3,
                    color="red",
                )
            )
        if has_virtual_people:
            for position, velocity in zip(virtual_people["xyrv"], virtual_people["U"]):
                ax.add_patch(
                    Arrow(
                        position[0],
                        position[1],
                        velocity[0],
                        velocity[1],
                        width=0.3,
                        color="red",
                        alpha=0.6,
                    )
                )

    if plot_desired_velocities:
        for position, velocity in zip(people["xyrv"], people["Vd"]):
            ax.add_patch(
                Arrow(
                    position[0],
                    position[1],
                    velocity[0],
                    velocity[1],
                    width=0.3,
                    color="blue",
                    alpha=0.8,
                )
            )
        if has_virtual_people:
            for position, velocity in zip(virtual_people["xyrv"], virtual_people["Vd"]):
                ax.add_patch(
                    Arrow(
                        position[0],
                        position[1],
                        velocity[0],
                        velocity[1],
                        width=0.3,
                        color="blue",
                        alpha=0.4,
                    )
                )

    if plot_paths:
        pathlines = [
            np.reshape(people["paths"][person_id], (-1, 2))
            for person_id in people["paths"]
        ]
        path_collection = LineCollection(
            pathlines,
            color="black",
            linewidths=0.5,
            linestyle="solid",
            alpha=0.5,
            cmap=plt.get_cmap(cmap),
        )
        ax.add_collection(path_collection)

    if plot_sensors:
        for sensor_data in sensors:
            line = Line2D(
                sensor_data["line"][::2],
                sensor_data["line"][1::2],
                lw=1.0,
                alpha=0.6,
                color="g",
            )
            ax.add_line(line)

    if axis is not None:
        ax.set_xlim(axis[0], axis[1])
        ax.set_ylim(axis[2], axis[3])
    ax.set_xticks([])
    ax.set_yticks([])
    ax.axis("off")
    if time >= 0:
        ax.set_title(f"time = {time:.2f} s")
    fig.canvas.draw()

    if savefig:
        fig.savefig(filename, dpi=dpi, bbox_inches="tight", pad_inches=0)

    return fig

    
def plot_sensors(ifig, sensors, time, flux_timestep=1, savefig=False,
                 filename='fig.png', dpi=150, cmap='winter'):
    """Plot sensor crossing times and entry/exit fluxes.

    Parameters
    ----------
    ifig : int
        Matplotlib figure number.
    sensors : sequence of dict
        Sensor records containing at least ``"name"``, ``"times"`` and
        ``"dir"``.
    time : float
        Current simulation time in seconds.
    flux_timestep : float, optional
        Width, in seconds, of bins used for entry and exit counts.
    savefig : bool, optional
        Save the figure when ``True``.
    filename : str, optional
        Output image filename.
    dpi : int, optional
        Output resolution in dots per inch.
    cmap : str, optional
        Matplotlib colormap name. Kept for API compatibility.

    Returns
    -------
    None
        The function updates the selected Matplotlib figure in place.
    """
    tmin = 0
    tmax = time

    Ns = len(sensors)
    if (Ns > 0):
        fig = plt.figure(ifig)
        # ncol = np.floor(np.sqrt(Ns))
        # nrow = np.ceil(Ns / ncol)
        ncol = 2
        nrow = Ns
        for i, s in enumerate(sensors):
            name = s["name"]  # sensor name
            tt = np.array(s["times"])
            dir = np.array(s["dir"])
            Np = tt.shape[0]
            ax1 = fig.add_subplot(nrow, ncol, 2*i+1)
            ax1.plot(np.arange(Np), tt, 'b+')
            ax1.set_title(name+": crossing times (s)", fontsize="small")

            tgrid = np.arange(tmin, tmax, step=flux_timestep)
            tgrid = np.append(tgrid, tgrid[-1]+flux_timestep)
            flux_exits = np.zeros(tgrid.shape)
            flux_entries = np.zeros(tgrid.shape)
            exits = np.where(dir == 1)[0]
            entries = np.where(dir == -1)[0]
            t_exits = np.ceil((tt[exits]-tmin)/flux_timestep)
            t_entries = np.ceil((tt[entries]-tmin)/flux_timestep)
            # t_exits = np.floor((tt[exits]-tmin)/flux_timestep)
            # t_entries = np.floor((tt[entries]-tmin)/flux_timestep)
            unique_exits, counts_exits = np.unique(t_exits, return_counts=True)
            unique_entries, counts_entries = np.unique(t_entries, return_counts=True)
            flux_exits[unique_exits.astype(int)] = counts_exits
            flux_entries[unique_entries.astype(int)] = counts_entries
            ax2 = fig.add_subplot(nrow, ncol, 2*i+2)
            ax2.plot(tgrid, flux_entries, ':og',
                     tgrid, flux_exits, ':or')
            ax2.set_title(name + ": entries (green) and exits (red) per " +
                          str(flux_timestep)+"s", fontsize="small")
        fig.set_tight_layout(True)
        fig.canvas.draw()
        if (savefig):
            fig.savefig(filename, dpi=dpi)
