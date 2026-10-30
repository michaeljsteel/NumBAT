// Template mesh geometry file for a single rounded wedge inclusion.

Mesh.MshFileVersion = 2.2;

dx_in_nm = 1000;
dy_in_nm = 1000;
rad_tip = 10;
alpha_tip = Pi/6;
bowlen_l = 100;
rad_wg = 1000;
rad_corner = 5;
lc_bkg = 0.1;
lc_refine_1 = lc_bkg/1;
lc_refine_2 = lc_bkg/2;

scale = 1/dx_in_nm;
outer_radius = rad_wg * scale;
wedge_radius = bowlen_l * scale;
tip_radius = rad_tip * scale;
corner_radius = rad_corner * scale;

// Tangent construction for the right side of the wedge, which has its tip at (0, wedge_radius).
tip_center_y = wedge_radius - tip_radius/Sin(alpha_tip);
tip_tangent_x = tip_radius * Cos(alpha_tip);
tip_tangent_y = tip_center_y + tip_radius * Sin(alpha_tip);

side_x = Sin(alpha_tip);
side_y = -Cos(alpha_tip);
inner_normal_x = -Cos(alpha_tip);
inner_normal_y = -Sin(alpha_tip);
offset_x = corner_radius * inner_normal_x;
offset_y = wedge_radius + corner_radius * inner_normal_y;
offset_along_side = offset_x * side_x + offset_y * side_y;
corner_center_distance = wedge_radius - corner_radius;
side_distance = -offset_along_side + Sqrt(
    offset_along_side^2 + corner_center_distance^2 - offset_x^2 - offset_y^2
);
corner_center_x = offset_x + side_distance * side_x;
corner_center_y = offset_y + side_distance * side_y;
side_end_x = corner_center_x - corner_radius * inner_normal_x;
side_end_y = corner_center_y - corner_radius * inner_normal_y;
outer_scale = wedge_radius / Sqrt(corner_center_x^2 + corner_center_y^2);
outer_end_x = outer_scale * corner_center_x;
outer_end_y = outer_scale * corner_center_y;

// Origin and outer circular simulation boundary.
Point(1) = {0, 0, 0, lc_bkg};
Point(2) = {0, outer_radius, 0, lc_bkg};
Point(3) = {-outer_radius, 0, 0, lc_bkg};
Point(4) = {0, -outer_radius, 0, lc_bkg};
Point(5) = {outer_radius, 0, 0, lc_bkg};

// Single wedge outline, traversed counter-clockwise from the right tip tangent point.
Point(10) = {tip_tangent_x, tip_tangent_y, 0, lc_refine_1};
Point(11) = {-tip_tangent_x, tip_tangent_y, 0, lc_refine_1};
Point(12) = {-side_end_x, side_end_y, 0, lc_refine_2};
Point(13) = {-outer_end_x, outer_end_y, 0, lc_refine_2};
Point(14) = {0, -wedge_radius, 0, lc_refine_2};
Point(15) = {outer_end_x, outer_end_y, 0, lc_refine_2};
Point(16) = {0, tip_center_y, 0, lc_refine_1};
Point(17) = {-corner_center_x, corner_center_y, 0, lc_refine_2};
Point(18) = {corner_center_x, corner_center_y, 0, lc_refine_2};
Point(19) = {side_end_x, side_end_y, 0, lc_refine_2};

Circle(20) = {10, 16, 11};
Line(21) = {11, 12};
Circle(22) = {12, 17, 13};
Circle(23) = {13, 1, 14};
Circle(24) = {14, 1, 15};
Circle(25) = {15, 18, 19};
Line(26) = {19, 10};
Line Loop(30) = {20, 21, 22, 23, 24, 25, 26};
Plane Surface(31) = {30};

// Background: the outer disk with the wedge loop removed.
Circle(100) = {2, 1, 3};
Circle(101) = {3, 1, 4};
Circle(102) = {4, 1, 5};
Circle(103) = {5, 1, 2};
Line Loop(110) = {100, 101, 102, 103};
Plane Surface(111) = {110, 30};

Physical Surface(1) = {111};  // Background material.
Physical Surface(2) = {31};   // Rounded wedge inclusion, material_a.
