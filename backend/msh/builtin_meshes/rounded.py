from dataclasses import dataclass

import numpy as np

from usermesh import UserGeometryBase
from ._common import nmtoum


def _corner_side_distance(apex_x, apex_y, side_x, side_y, normal_x, normal_y, rad_corner, bowlen_l):
    """Return the distance from a side's apex to the tangent point of its outer corner fillet.

    The fillet has radius rad_corner, its centre lies rad_corner inside the side along the
    unit normal, and it touches the circle of radius bowlen_l from inside. Returns None if
    no such fillet exists.
    """
    offset_x = apex_x + rad_corner * normal_x
    offset_y = apex_y + rad_corner * normal_y
    offset_along_side = offset_x * side_x + offset_y * side_y
    discriminant = (
        offset_along_side**2 + (bowlen_l - rad_corner) ** 2 - offset_x**2 - offset_y**2
    )
    if discriminant <= 0:
        return None
    return -offset_along_side + np.sqrt(discriminant)


@dataclass(frozen=True)
class _RoundedBowtieGeometry:
    tip_center_x: float
    tip_tangent_x: float
    tip_tangent_y: float
    corner_center_x: float
    corner_center_y: float
    side_end_x: float
    side_end_y: float
    outer_end_x: float
    outer_end_y: float
    side_distance: float
    tip_tangent_distance: float


def _make_rounded_bowtie_geometry(rad_tip, alpha_tip, gap_g, bowlen_l, rad_corner):
    """Return the tangent construction for the right rounded bowtie inclusion."""
    tip_x = gap_g / 2
    virtual_tip_x = tip_x - rad_tip * (1 / np.sin(alpha_tip) - 1)
    tip_center_x = tip_x + rad_tip
    tip_tangent_x = tip_center_x - rad_tip * np.sin(alpha_tip)
    tip_tangent_y = rad_tip * np.cos(alpha_tip)

    side_x = np.cos(alpha_tip)
    side_y = np.sin(alpha_tip)
    inner_normal_x = np.sin(alpha_tip)
    inner_normal_y = -np.cos(alpha_tip)
    side_distance = _corner_side_distance(
        virtual_tip_x, 0.0, side_x, side_y, inner_normal_x, inner_normal_y, rad_corner, bowlen_l
    )
    if side_distance is None:
        return None

    corner_center_x = virtual_tip_x + rad_corner * inner_normal_x + side_distance * side_x
    corner_center_y = rad_corner * inner_normal_y + side_distance * side_y
    side_end_x = corner_center_x - rad_corner * inner_normal_x
    side_end_y = corner_center_y - rad_corner * inner_normal_y
    outer_scale = bowlen_l / np.hypot(corner_center_x, corner_center_y)

    return _RoundedBowtieGeometry(
        tip_center_x=tip_center_x,
        tip_tangent_x=tip_tangent_x,
        tip_tangent_y=tip_tangent_y,
        corner_center_x=corner_center_x,
        corner_center_y=corner_center_y,
        side_end_x=side_end_x,
        side_end_y=side_end_y,
        outer_end_x=outer_scale * corner_center_x,
        outer_end_y=outer_scale * corner_center_y,
        side_distance=side_distance,
        tip_tangent_distance=rad_tip / np.tan(alpha_tip),
    )


@dataclass(frozen=True)
class _RoundedWedgeGeometry:
    tip_center_y: float
    tip_tangent_x: float
    tip_tangent_y: float
    corner_center_x: float
    corner_center_y: float
    side_end_x: float
    side_end_y: float
    outer_end_x: float
    outer_end_y: float
    side_distance: float
    tip_tangent_distance: float


def _make_rounded_wedge_geometry(rad_tip, alpha_tip, bowlen_l, rad_corner):
    """Return the tangent construction for the right side of the rounded wedge.

    The sides of the wedge meet at a sharp virtual tip on the circle of radius bowlen_l.
    The actual tip is rounded with radius rad_tip, lying just inside that point.
    """
    tip_center_y = bowlen_l - rad_tip / np.sin(alpha_tip)
    tip_tangent_x = rad_tip * np.cos(alpha_tip)
    tip_tangent_y = tip_center_y + rad_tip * np.sin(alpha_tip)

    side_x = np.sin(alpha_tip)
    side_y = -np.cos(alpha_tip)
    inner_normal_x = -np.cos(alpha_tip)
    inner_normal_y = -np.sin(alpha_tip)
    side_distance = _corner_side_distance(
        0.0, bowlen_l, side_x, side_y, inner_normal_x, inner_normal_y, rad_corner, bowlen_l
    )
    if side_distance is None:
        return None

    corner_center_x = rad_corner * inner_normal_x + side_distance * side_x
    corner_center_y = bowlen_l + rad_corner * inner_normal_y + side_distance * side_y
    side_end_x = corner_center_x - rad_corner * inner_normal_x
    side_end_y = corner_center_y - rad_corner * inner_normal_y
    outer_scale = bowlen_l / np.hypot(corner_center_x, corner_center_y)

    return _RoundedWedgeGeometry(
        tip_center_y=tip_center_y,
        tip_tangent_x=tip_tangent_x,
        tip_tangent_y=tip_tangent_y,
        corner_center_x=corner_center_x,
        corner_center_y=corner_center_y,
        side_end_x=side_end_x,
        side_end_y=side_end_y,
        outer_end_x=outer_scale * corner_center_x,
        outer_end_y=outer_scale * corner_center_y,
        side_distance=side_distance,
        tip_tangent_distance=rad_tip / np.tan(alpha_tip),
    )


def _rounded_wedge_shape_problems(rad_tip, alpha_tip, bowlen_l, rad_corner):
    """Return a message describing invalid wedge shape parameters, or "" if they are valid."""
    if rad_tip <= 0:
        return "Tip radius (rad_tip) must be positive.\n"
    if not 0 < alpha_tip < np.pi / 2:
        return "Tip half-angle (alpha_tip) must be between 0 and pi/2 radians.\n"
    if bowlen_l <= 0:
        return "Wedge length (bowlen_l) must be positive.\n"
    if rad_corner <= 0:
        return "Corner radius (rad_corner) must be positive.\n"

    corner_limit = 2 * bowlen_l * np.sin(alpha_tip) / (1 + np.sin(alpha_tip))
    if rad_corner >= corner_limit:
        return (
            "Corner radius (rad_corner) must be less than "
            "2 * bowlen_l * sin(alpha_tip) / (1 + sin(alpha_tip)) so that the corner "
            "transitions do not meet on the vertical axis.\n"
            f"For the current parameters, rad_corner < {corner_limit:.3f} is required.\n"
        )

    geometry = _make_rounded_wedge_geometry(rad_tip, alpha_tip, bowlen_l, rad_corner)
    if geometry is None or geometry.side_distance <= geometry.tip_tangent_distance:
        return "Corner radius (rad_corner) cannot form a corner transition that clears the rounded tip.\n"
    return ""


def _ccw_arc_points(center_x, center_y, start_x, start_y, end_x, end_y, num_points=41):
    """Return points on the counter-clockwise arc from start to end about the centre."""
    start_angle = np.arctan2(start_y - center_y, start_x - center_x)
    end_angle = np.arctan2(end_y - center_y, end_x - center_x)
    sweep = (end_angle - start_angle) % (2 * np.pi)
    radius = np.hypot(start_x - center_x, start_y - center_y)
    angles = start_angle + np.linspace(0, sweep, num_points)
    return np.column_stack(
        (center_x + radius * np.cos(angles), center_y + radius * np.sin(angles))
    )


def _rotate(x, y, angle):
    """Rotate points (x, y) about the origin by angle in radians. Accepts scalars or arrays."""
    cos_a, sin_a = np.cos(angle), np.sin(angle)
    return cos_a * x - sin_a * y, sin_a * x + cos_a * y


def _valid_wedge_geometry(rad_tip, alpha_tip, gap_g, bowlen_l, rad_corner):
    """Return the tangent construction for a valid wedge, or raise ValueError."""
    geometry = _make_rounded_bowtie_geometry(rad_tip, alpha_tip, gap_g, bowlen_l, rad_corner)
    if (
        geometry is None
        or geometry.corner_center_y <= 0
        or geometry.side_distance <= geometry.tip_tangent_distance
    ):
        raise ValueError("Invalid rounded wedge geometry; call validate_dimensions() for details.")
    return geometry


def _outer_diameter_messages(rad_wg, dom_x, dom_y):
    msg = ""
    diam_outer = rad_wg * 2
    if diam_outer >= dom_x:
        msg += "Outer cylinder has total diameter larger than domain width (domain_x).\n"
    if diam_outer >= dom_y:
        msg += "Outer cylinder has total diameter larger than domain height (domain_y).\n"
    return msg


def _wedge_geometry_messages(rad_tip, alpha_tip, gap_g, bowlen_l, rad_corner):
    """Return error text for one rounded wedge with these parameters, or an empty string."""
    msg = ""

    if rad_tip <= 0:
        msg += "Tip radius (rad_tip) must be positive.\n"
    if not 0 < alpha_tip < np.pi / 2:
        msg += "Tip half-angle (alpha_tip) must be between 0 and pi/2 radians.\n"
    if gap_g <= 0:
        msg += "Gap (gap_g) must be positive.\n"
    if bowlen_l <= gap_g / 2:
        msg += "Bowtie length (bowlen_l) must exceed half the gap (gap_g).\n"
    if not 0 < rad_corner < bowlen_l:
        msg += "Corner radius (rad_corner) must be positive and less than bowlen_l.\n"

    if not rad_corner < bowlen_l * np.sin(alpha_tip) / (1 + np.sin(alpha_tip)):
        msg+= "Corner radius (rad_corner) must be less than bowlen_l * sin(alpha_tip) / (1 + sin(alpha_tip)) to avoid overlap with the tip.\n"
        msg+= "For the current parameters, rad_corner < %.3f is required.\n" % (bowlen_l * np.sin(alpha_tip) / (1 + np.sin(alpha_tip)))

    if (
        0 < alpha_tip < np.pi / 2
        and rad_tip > 0
        and gap_g > 0
        and 0 < rad_corner < bowlen_l
    ):
        geometry = _make_rounded_bowtie_geometry(
            rad_tip, alpha_tip, gap_g, bowlen_l, rad_corner
        )

        if geometry is None:
            msg += (
                "Corner radius (rad_corner) cannot form a tangent transition between "
                "the bowtie side and outer boundary.\n"
            )
        else:
            if geometry.corner_center_y <= 0:
                msg += (
                    "Corner radius (rad_corner) makes the outer corner transitions "
                    "cross at the horizontal axis.\n"
                )
            if geometry.side_distance <= geometry.tip_tangent_distance:
                msg += (
                    "Corner radius (rad_corner) overlaps the rounded tip transition.\n"
                )
    return msg


def _arc_points(center_x, center_y, start_x, start_y, end_x, end_y, short=True):
    """Return points along the circle about (center) from (start) to (end), for plotting."""
    start_angle = np.arctan2(start_y - center_y, start_x - center_x)
    end_angle = np.arctan2(end_y - center_y, end_x - center_x)
    if short:
        end_angle = start_angle + (end_angle - start_angle + np.pi) % (2 * np.pi) - np.pi
    radius = np.hypot(start_x - center_x, start_y - center_y)
    angles = np.linspace(start_angle, end_angle, 41)
    return np.column_stack(
        (center_x + radius * np.cos(angles), center_y + radius * np.sin(angles))
    )


def _right_wedge_outline(geometry):
    """Return the closed outline of the right wedge (unrotated, in nm) as an (N, 2) array."""
    tip_arc = _arc_points(
        geometry.tip_center_x,
        0,
        geometry.tip_tangent_x,
        -geometry.tip_tangent_y,
        geometry.tip_tangent_x,
        geometry.tip_tangent_y,
    )
    upper_corner_arc = _arc_points(
        geometry.corner_center_x,
        geometry.corner_center_y,
        geometry.side_end_x,
        geometry.side_end_y,
        geometry.outer_end_x,
        geometry.outer_end_y,
    )
    outer_arc = _arc_points(
        0,
        0,
        geometry.outer_end_x,
        geometry.outer_end_y,
        geometry.outer_end_x,
        -geometry.outer_end_y,
        short=False,
    )
    lower_corner_arc = _arc_points(
        geometry.corner_center_x,
        -geometry.corner_center_y,
        geometry.outer_end_x,
        -geometry.outer_end_y,
        geometry.side_end_x,
        -geometry.side_end_y,
    )
    return np.vstack(
        (
            tip_arc,
            [geometry.side_end_x, geometry.side_end_y],
            upper_corner_arc,
            outer_arc,
            lower_corner_arc,
            [geometry.tip_tangent_x, -geometry.tip_tangent_y],
        )
    )


class RoundedBowtie(UserGeometryBase):
    """NumBAT geometry template for a Jesper Mork bowtie resonator."""

    def init_geometry(self):

        self.set_properties("rounded_bowtie", is_curvi=True)

        self.set_required_parameters(["rad_tip", "rad_wg", "alpha_tip", "gap_g", "bowlen_l", "rad_corner"], num_req_mats=3)
        self.set_allowed_parameters(["lc_bkg", "lc_refine_1", "lc_refine_2"])

        self.set_parameter_help(
            {
                "rad_tip": "radius of curvature at the tip of the bowtie",
                "rad_wg": "radius of the whole waveguide",
                "alpha_tip": "angle of the tip of the bowtie",
                "gap_g": "gap between the two bowtie inclusions",
                "bowlen_l": "length of the bowtie",
                "rad_corner": "radius of curvature at the corners of the bowtie",
                "material_a": "material of left bowtie inclusion",
                "material_b": "material of right bowtie inclusion",
                "lc_bkg": "mesh spacing on outer boundary",
                "lc_refine_1": "mesh refinement at the tip",
                "lc_refine_2": "mesh refinement at the outer boundary of the bowtie inclusion",

            }
        )

    def apply_parameters(self):
        subs = [
           ("dx_in_nm = 1000;", "dx_in_nm = %f;", "domain_x"),
           ("dy_in_nm = 1000;", "dy_in_nm = %f;", "domain_y"),
           ("rad_wg = 1000;", "rad_wg = %f;", "rad_wg"),
           ("rad_tip = 10;", "rad_tip = %f;", "rad_tip"),
           ("alpha_tip = Pi/6;", "alpha_tip = %f;", "alpha_tip"),
           ("gap_g = 20;", "gap_g = %f;", "gap_g"),
           ("bowlen_l = 100;", "bowlen_l = %f;", "bowlen_l"),
           ("rad_corner = 5;", "rad_corner = %f;", "rad_corner"),
           ("lc_bkg = 0.1;", "lc_bkg = %f;", "lc_bkg"),
           ("lc_refine_1 = lc_bkg/1;", "lc_refine_1 = lc_bkg/%f;", "lc_refine_1"),
           ("lc_refine_2 = lc_bkg/2;", "lc_refine_2 = lc_bkg/%f;", "lc_refine_2"),
          ]
        return subs

    def check_dimensions(self):
        dom_x = self.get_param("domain_x")
        dom_y = self.get_param("domain_y")
        rad_tip = self.get_param("rad_tip")
        rad_wg = self.get_param("rad_wg")
        gap_g = self.get_param("gap_g")
        bowlen_l = self.get_param("bowlen_l")
        rad_corner = self.get_param("rad_corner")
        alpha_tip = self.get_param("alpha_tip")

        msg = _outer_diameter_messages(rad_wg, dom_x, dom_y)
        msg += _wedge_geometry_messages(rad_tip, alpha_tip, gap_g, bowlen_l, rad_corner)

        dims_ok = not len(msg)
        return dims_ok, msg

    def draw_mpl_frame(self, ax, styles):
        geometry = _valid_wedge_geometry(
            self.get_param("rad_tip"),
            self.get_param("alpha_tip"),
            self.get_param("gap_g"),
            self.get_param("bowlen_l"),
            self.get_param("rad_corner"),
        )
        right_outline = _right_wedge_outline(geometry)

        for outline in (right_outline, right_outline * np.array([-1, 1])):
            ax.plot(
                outline[:, 0] * nmtoum,
                outline[:, 1] * nmtoum,
                color=styles["edgecolor"],
                linewidth=styles["linewidth"],
            )


class RoundedWedge(UserGeometryBase):
    """NumBAT geometry template for a single rounded wedge inside a uniform fibre.

    The wedge points upwards, is symmetric about the vertical axis and is centred in the fibre.
    """

    def init_geometry(self):

        self.set_properties("rounded_wedge", is_curvi=True)

        self.set_required_parameters(["rad_tip", "rad_wg", "alpha_tip", "bowlen_l", "rad_corner"], num_req_mats=2)
        self.set_allowed_parameters(["lc_bkg", "lc_refine_1", "lc_refine_2"])

        self.set_parameter_help(
            {
                "rad_tip": "radius of curvature at the tip of the wedge",
                "rad_wg": "radius of the whole waveguide",
                "alpha_tip": "half-angle of the tip of the wedge, measured from the vertical axis",
                "bowlen_l": "distance from the centre of the waveguide to the sharp (unrounded) tip of the wedge, which is also the radius of its base",
                "rad_corner": "radius of curvature where the sides of the wedge meet its base",
                "material_a": "material of the wedge",
                "lc_bkg": "mesh spacing on outer boundary",
                "lc_refine_1": "mesh refinement at the tip",
                "lc_refine_2": "mesh refinement at the outer boundary of the wedge",
            }
        )

    def apply_parameters(self):
        subs = [
           ("dx_in_nm = 1000;", "dx_in_nm = %f;", "domain_x"),
           ("dy_in_nm = 1000;", "dy_in_nm = %f;", "domain_y"),
           ("rad_wg = 1000;", "rad_wg = %f;", "rad_wg"),
           ("rad_tip = 10;", "rad_tip = %f;", "rad_tip"),
           ("alpha_tip = Pi/6;", "alpha_tip = %f;", "alpha_tip"),
           ("bowlen_l = 100;", "bowlen_l = %f;", "bowlen_l"),
           ("rad_corner = 5;", "rad_corner = %f;", "rad_corner"),
           ("lc_bkg = 0.1;", "lc_bkg = %f;", "lc_bkg"),
           ("lc_refine_1 = lc_bkg/1;", "lc_refine_1 = lc_bkg/%f;", "lc_refine_1"),
           ("lc_refine_2 = lc_bkg/2;", "lc_refine_2 = lc_bkg/%f;", "lc_refine_2"),
          ]
        return subs

    def check_dimensions(self):
        dom_x = self.get_param("domain_x")
        dom_y = self.get_param("domain_y")
        rad_tip = self.get_param("rad_tip")
        rad_wg = self.get_param("rad_wg")
        bowlen_l = self.get_param("bowlen_l")
        rad_corner = self.get_param("rad_corner")
        alpha_tip = self.get_param("alpha_tip")

        msg = ""

        diam_outer = rad_wg * 2

        if diam_outer >= dom_x:
            msg += f"Outer cylinder has total diameter={diam_outer} exceeding the domain width (domain_x)={dom_x}.\n"
        if diam_outer >= dom_y:
            msg += f"Outer cylinder has total diameter={diam_outer} exceeding the domain width (domain_y)={dom_y}.\n"
        if rad_wg <= bowlen_l:
            msg += "Wedge length (bowlen_l) must be less than the waveguide radius (rad_wg).\n"
        msg += _rounded_wedge_shape_problems(rad_tip, alpha_tip, bowlen_l, rad_corner)

        dims_ok = not len(msg)
        return dims_ok, msg

    def draw_mpl_frame(self, ax, styles):
        alpha_tip = self.get_param("alpha_tip")
        bowlen_l = self.get_param("bowlen_l")
        rad_tip = self.get_param("rad_tip")
        rad_corner = self.get_param("rad_corner")

        problems = _rounded_wedge_shape_problems(rad_tip, alpha_tip, bowlen_l, rad_corner)
        if problems:
            raise ValueError("Invalid rounded wedge geometry:\n" + problems)

        g = _make_rounded_wedge_geometry(rad_tip, alpha_tip, bowlen_l, rad_corner)
        tx, ty = g.tip_tangent_x, g.tip_tangent_y
        sx, sy = g.side_end_x, g.side_end_y
        cx, cy = g.corner_center_x, g.corner_center_y
        ox, oy = g.outer_end_x, g.outer_end_y

        outline = np.vstack(
            (
                _ccw_arc_points(0.0, g.tip_center_y, tx, ty, -tx, ty),
                _ccw_arc_points(-cx, cy, -sx, sy, -ox, oy),
                _ccw_arc_points(0.0, 0.0, -ox, oy, ox, oy),
                _ccw_arc_points(cx, cy, ox, oy, sx, sy),
                [[tx, ty]],
            )
        )

        ax.plot(
            outline[:, 0] * nmtoum,
            outline[:, 1] * nmtoum,
            color=styles["edgecolor"],
            linewidth=styles["linewidth"],
        )


class _RoundedWedgeN(UserGeometryBase):
    """Base class for N identical rounded wedges arranged rotationally symmetrically about the origin.

    Each wedge is the right wedge of RoundedBowtie, rotated about the origin, so the
    tangent construction is shared with RoundedBowtie rather than repeated per wedge.
    """

    _num_wedges = 0
    _first_wedge_angle_deg = 0.0
    _alpha_tip_limit_deg = 0.0
    _inc_shape = ""

    def init_geometry(self):
        self.set_properties(self._inc_shape, is_curvi=True)

        self.set_required_parameters(
            ["rad_tip", "rad_wg", "alpha_tip", "gap_g", "bowlen_l", "rad_corner"], num_req_mats=2
        )
        self.set_allowed_parameters(["lc_bkg", "lc_refine_1", "lc_refine_2", "rot_deg"])

        self.set_parameter_help(
            {
                "rad_tip": "radius of curvature at the tip of each wedge",
                "rad_wg": "radius of the whole waveguide",
                "alpha_tip": f"half-angle of the tip of each wedge, must be less than {self._alpha_tip_limit_deg:g} degrees",
                "gap_g": "gap between the centre and the tip of each wedge is gap_g/2",
                "bowlen_l": "length of each wedge",
                "rad_corner": "radius of curvature at the corners of each wedge",
                "material_a": "material of all wedges",
                "lc_bkg": "mesh spacing on outer boundary",
                "lc_refine_1": "mesh refinement at the tip (divides lc_bkg)",
                "lc_refine_2": "mesh refinement at the outer boundary of each wedge (divides lc_bkg)",
                "rot_deg": "rotation of the whole structure about the origin, in degrees",
            }
        )

    def check_dimensions(self):
        alpha_tip = self.get_param("alpha_tip")

        msg = _outer_diameter_messages(
            self.get_param("rad_wg"), self.get_param("domain_x"), self.get_param("domain_y")
        )
        if not 0 < alpha_tip < np.radians(self._alpha_tip_limit_deg):
            msg += (
                f"Tip half-angle (alpha_tip) must be less than {self._alpha_tip_limit_deg:g} degrees "
                f"so that the {self._num_wedges} wedges do not overlap.\n"
            )
        msg += _wedge_geometry_messages(
            self.get_param("rad_tip"),
            alpha_tip,
            self.get_param("gap_g"),
            self.get_param("bowlen_l"),
            self.get_param("rad_corner"),
        )

        dims_ok = not len(msg)
        return dims_ok, msg

    def _single_wedge_geometry(self):
        return _valid_wedge_geometry(
            self.get_param("rad_tip"),
            self.get_param("alpha_tip"),
            self.get_param("gap_g"),
            self.get_param("bowlen_l"),
            self.get_param("rad_corner"),
        )

    def _wedge_angle_deg(self, k):
        """Return the angle of the k-th wedge axis about the origin, including the global rotation."""
        return (
            self._first_wedge_angle_deg
            + self.get_param("rot_deg", 0.0)
            + k * 360.0 / self._num_wedges
        )

    def draw_mpl_frame(self, ax, styles):
        outline = _right_wedge_outline(self._single_wedge_geometry())

        for k in range(self._num_wedges):
            x, y = _rotate(outline[:, 0], outline[:, 1], np.radians(self._wedge_angle_deg(k)))
            ax.plot(
                x * nmtoum,
                y * nmtoum,
                color=styles["edgecolor"],
                linewidth=styles["linewidth"],
            )

    def make_geometry(self, p_dir_templates):
        """Return the Gmsh script for all wedges, built by rotating the single-wedge construction.

        Coordinates are scaled by 1/domain_x, as in the RoundedBowtie template.
        """
        geometry = self._single_wedge_geometry()
        scale = 1.0 / self.get_param("domain_x")
        lc_bkg = self.get_param("lc_bkg")
        lc_refine_1 = lc_bkg / self.get_param("lc_refine_1", 1)
        lc_refine_2 = lc_bkg / self.get_param("lc_refine_2", 2)
        outer_radius = self.get_param("rad_wg") * scale

        def fmt(v):
            return f"{v:.12g}"

        lines = [
            "Mesh.MshFileVersion = 2.2;",
            "",
            f"Point(1) = {{0, 0, 0, {fmt(lc_bkg)}}};",
            f"Point(2) = {{0, {fmt(outer_radius)}, 0, {fmt(lc_bkg)}}};",
            f"Point(3) = {{{fmt(-outer_radius)}, 0, 0, {fmt(lc_bkg)}}};",
            f"Point(4) = {{0, {fmt(-outer_radius)}, 0, {fmt(lc_bkg)}}};",
            f"Point(5) = {{{fmt(outer_radius)}, 0, 0, {fmt(lc_bkg)}}};",
            "Circle(100) = {2, 1, 3};",
            "Circle(101) = {3, 1, 4};",
            "Circle(102) = {4, 1, 5};",
            "Circle(103) = {5, 1, 2};",
            "Line Loop(110) = {100, 101, 102, 103};",
        ]

        # Single-wedge points of the upper half, then the mirrored lower half, as in the bowtie template.
        g = geometry
        local_pts = [
            (g.tip_tangent_x, g.tip_tangent_y),
            (g.side_end_x, g.side_end_y),
            (g.outer_end_x, g.outer_end_y),
            (g.outer_end_x, -g.outer_end_y),
            (g.side_end_x, -g.side_end_y),
            (g.tip_tangent_x, -g.tip_tangent_y),
            (g.tip_center_x, 0.0),
            (g.corner_center_x, g.corner_center_y),
            (g.corner_center_x, -g.corner_center_y),
        ]
        tip_point_offsets = (0, 5, 6)

        wedge_loops = []
        wedge_surfaces = []
        for k in range(self._num_wedges):
            p0 = 10 + 9 * k
            c0 = 10 + 6 * k
            loop = 20 + 2 * k
            theta = np.radians(self._wedge_angle_deg(k))

            for j, (x, y) in enumerate(local_pts):
                rx, ry = _rotate(x * scale, y * scale, theta)
                lc = lc_refine_1 if j in tip_point_offsets else lc_refine_2
                lines.append(f"Point({p0 + j}) = {{{fmt(rx)}, {fmt(ry)}, 0, {fmt(lc)}}};")

            p = [p0 + j for j in range(9)]
            lines += [
                f"Circle({c0}) = {{{p[5]}, {p[6]}, {p[0]}}};",
                f"Line({c0 + 1}) = {{{p[0]}, {p[1]}}};",
                f"Circle({c0 + 2}) = {{{p[1]}, {p[7]}, {p[2]}}};",
                f"Circle({c0 + 3}) = {{{p[2]}, 1, {p[3]}}};",
                f"Circle({c0 + 4}) = {{{p[3]}, {p[8]}, {p[4]}}};",
                f"Line({c0 + 5}) = {{{p[4]}, {p[5]}}};",
                f"Line Loop({loop}) = {{{', '.join(str(c0 + j) for j in range(6))}}};",
                f"Plane Surface({loop + 1}) = {{{loop}}};",
            ]
            wedge_loops.append(str(loop))
            wedge_surfaces.append(str(loop + 1))

        lines += [
            f"Plane Surface(111) = {{110, {', '.join(wedge_loops)}}};",
            "",
            "Physical Surface(1) = {111};  // Background material.",
            f"Physical Surface(2) = {{{', '.join(wedge_surfaces)}}};  // All wedges, material_a.",
        ]
        return "\n".join(lines) + "\n"


class RoundedWedge2(_RoundedWedgeN):
    """NumBAT geometry template for three identical rounded wedges, one pointing up."""

    _num_wedges = 2
    _first_wedge_angle_deg = 0.0
    _alpha_tip_limit_deg = 90.0
    _inc_shape = "rounded_wedge2"

class RoundedWedge3(_RoundedWedgeN):
    """NumBAT geometry template for three identical rounded wedges, one pointing up."""

    _num_wedges = 3
    _first_wedge_angle_deg = 90.0
    _alpha_tip_limit_deg = 60.0
    _inc_shape = "rounded_wedge3"


class RoundedWedge4(_RoundedWedgeN):
    """NumBAT geometry template for four identical rounded wedges, centred on NE, NW, SW and SE."""

    _num_wedges = 4
    _first_wedge_angle_deg = 45.0
    _alpha_tip_limit_deg = 45.0
    _inc_shape = "rounded_wedge4"


class RoundedWedge6(_RoundedWedgeN):
    """NumBAT geometry template for six identical rounded wedges, two of them aligned with E and W."""

    _num_wedges = 6
    _first_wedge_angle_deg = 0.0
    _alpha_tip_limit_deg = 30.0
    _inc_shape = "rounded_wedge6"

