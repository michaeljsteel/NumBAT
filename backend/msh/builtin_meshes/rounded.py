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

        msg = ""

        diam_outer = rad_wg * 2

        if diam_outer >= dom_x:
            msg += "Outer cylinder has total diameter larger than domain width (domain_x).\n"
        if diam_outer >= dom_y:
            msg += "Outer cylinder has total diameter larger than domain height (domain_y).\n"
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
            alpha_tip = self.get_param("alpha_tip")
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

        dims_ok = not len(msg)
        return dims_ok, msg

    def draw_mpl_frame(self, ax, styles):
        alpha_tip = self.get_param("alpha_tip")
        bowlen_l = self.get_param("bowlen_l")
        gap_g = self.get_param("gap_g")
        rad_tip = self.get_param("rad_tip")
        rad_corner = self.get_param("rad_corner")

        geometry = _make_rounded_bowtie_geometry(
            rad_tip, alpha_tip, gap_g, bowlen_l, rad_corner
        )
        if (
            geometry is None
            or geometry.corner_center_y <= 0
            or geometry.side_distance <= geometry.tip_tangent_distance
        ):
            raise ValueError("Invalid rounded bowtie geometry; call validate_dimensions() for details.")

        def arc_points(center_x, center_y, start_x, start_y, end_x, end_y, short=True):
            start_angle = np.arctan2(start_y - center_y, start_x - center_x)
            end_angle = np.arctan2(end_y - center_y, end_x - center_x)
            if short:
                end_angle = start_angle + (end_angle - start_angle + np.pi) % (2 * np.pi) - np.pi
            radius = np.hypot(start_x - center_x, start_y - center_y)
            angles = np.linspace(start_angle, end_angle, 41)
            return np.column_stack(
                (center_x + radius * np.cos(angles), center_y + radius * np.sin(angles))
            )

        tip_arc = arc_points(
            geometry.tip_center_x,
            0,
            geometry.tip_tangent_x,
            -geometry.tip_tangent_y,
            geometry.tip_tangent_x,
            geometry.tip_tangent_y,
        )
        upper_corner_arc = arc_points(
            geometry.corner_center_x,
            geometry.corner_center_y,
            geometry.side_end_x,
            geometry.side_end_y,
            geometry.outer_end_x,
            geometry.outer_end_y,
        )
        outer_arc = arc_points(
            0,
            0,
            geometry.outer_end_x,
            geometry.outer_end_y,
            geometry.outer_end_x,
            -geometry.outer_end_y,
            short=False,
        )
        lower_corner_arc = arc_points(
            geometry.corner_center_x,
            -geometry.corner_center_y,
            geometry.outer_end_x,
            -geometry.outer_end_y,
            geometry.side_end_x,
            -geometry.side_end_y,
        )
        right_outline = np.vstack(
            (
                tip_arc,
                [geometry.side_end_x, geometry.side_end_y],
                upper_corner_arc,
                outer_arc,
                lower_corner_arc,
                [geometry.tip_tangent_x, -geometry.tip_tangent_y],
            )
        )

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
