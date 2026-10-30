// Template mesh geometry file for two rounded bowtie inclusions.

Mesh.MshFileVersion = 2.2;

dx_in_nm = 1000;
dy_in_nm = 1000;
rad_tip = 10;
alpha_tip = Pi/6;
gap_g = 20;
bowlen_l = 100;
rad_wg = 1000;
rad_corner = 5;
lc_bkg = 0.1;
lc_refine_1 = lc_bkg/1;
lc_refine_2 = lc_bkg/2;

scale = 1/dx_in_nm;
outer_radius = rad_wg * scale;
inclusion_radius = bowlen_l * scale;
tip_radius = rad_tip * scale;
corner_radius = rad_corner * scale;
gap = gap_g * scale;

// Tangent construction for the upper side of the right inclusion.
tip_x = gap/2;
virtual_tip_x = tip_x - tip_radius * (1/Sin(alpha_tip) - 1);
tip_center_x = tip_x + tip_radius;
tip_tangent_x = tip_center_x - tip_radius * Sin(alpha_tip);
tip_tangent_y = tip_radius * Cos(alpha_tip);

side_x = Cos(alpha_tip);
side_y = Sin(alpha_tip);
inner_normal_x = Sin(alpha_tip);
inner_normal_y = -Cos(alpha_tip);
offset_x = virtual_tip_x + corner_radius * inner_normal_x;
offset_y = corner_radius * inner_normal_y;
offset_along_side = offset_x * side_x + offset_y * side_y;
corner_center_distance = inclusion_radius - corner_radius;
side_distance = -offset_along_side + Sqrt(
    offset_along_side^2 + corner_center_distance^2 - offset_x^2 - offset_y^2
);
corner_center_x = offset_x + side_distance * side_x;
corner_center_y = offset_y + side_distance * side_y;
side_end_x = corner_center_x - corner_radius * inner_normal_x;
side_end_y = corner_center_y - corner_radius * inner_normal_y;
outer_scale = inclusion_radius / Sqrt(corner_center_x^2 + corner_center_y^2);
outer_end_x = outer_scale * corner_center_x;
outer_end_y = outer_scale * corner_center_y;

// Origin and outer circular simulation boundary.
Point(1) = {0, 0, 0, lc_bkg};
Point(2) = {0, outer_radius, 0, lc_bkg};
Point(3) = {-outer_radius, 0, 0, lc_bkg};
Point(4) = {0, -outer_radius, 0, lc_bkg};
Point(5) = {outer_radius, 0, 0, lc_bkg};

// Right bowtie inclusion.
Point(10) = {tip_tangent_x, tip_tangent_y, 0, lc_refine_1};
Point(11) = {side_end_x, side_end_y, 0, lc_refine_2};
Point(12) = {outer_end_x, outer_end_y, 0, lc_refine_2};
Point(13) = {outer_end_x, -outer_end_y, 0, lc_refine_2};
Point(14) = {side_end_x, -side_end_y, 0, lc_refine_2};
Point(15) = {tip_tangent_x, -tip_tangent_y, 0, lc_refine_1};
Point(16) = {tip_center_x, 0, 0, lc_refine_1};
Point(17) = {corner_center_x, corner_center_y, 0, lc_refine_2};
Point(18) = {corner_center_x, -corner_center_y, 0, lc_refine_2};

Circle(10) = {15, 16, 10};
Line(11) = {10, 11};
Circle(12) = {11, 17, 12};
Circle(13) = {12, 1, 13};
Circle(14) = {13, 18, 14};
Line(15) = {14, 15};
Line Loop(20) = {10, 11, 12, 13, 14, 15};
Plane Surface(21) = {20};

// Left bowtie inclusion, mirrored about the y axis.
Point(30) = {-tip_tangent_x, tip_tangent_y, 0, lc_refine_1};
Point(31) = {-side_end_x, side_end_y, 0, lc_refine_2};
Point(32) = {-outer_end_x, outer_end_y, 0, lc_refine_2};
Point(33) = {-outer_end_x, -outer_end_y, 0, lc_refine_2};
Point(34) = {-side_end_x, -side_end_y, 0, lc_refine_2};
Point(35) = {-tip_tangent_x, -tip_tangent_y, 0, lc_refine_1};
Point(36) = {-tip_center_x, 0, 0, lc_refine_1};
Point(37) = {-corner_center_x, corner_center_y, 0, lc_refine_2};
Point(38) = {-corner_center_x, -corner_center_y, 0, lc_refine_2};

Circle(30) = {35, 36, 30};
Line(31) = {30, 31};
Circle(32) = {31, 37, 32};
Circle(33) = {32, 1, 33};
Circle(34) = {33, 38, 34};
Line(35) = {34, 35};
Line Loop(40) = {30, 31, 32, 33, 34, 35};
Plane Surface(41) = {40};

// Background: the outer disk with the two inclusion loops removed.
Circle(100) = {2, 1, 3};
Circle(101) = {3, 1, 4};
Circle(102) = {4, 1, 5};
Circle(103) = {5, 1, 2};
Line Loop(110) = {100, 101, 102, 103};
Plane Surface(111) = {110, 20, 40};

Physical Surface(1) = {111};  // Background material.
Physical Surface(2) = {21};   // Right bowtie inclusion, material_a.
Physical Surface(3) = {41};   // Left bowtie inclusion, material_b.
