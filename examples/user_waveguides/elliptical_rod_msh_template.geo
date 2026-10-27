// A single elliptical core surrounded by one background region.
Mesh.MshFileVersion = 2.2;

d = 1;
dx_in_nm = 100;
dy_in_nm = 50;
dy = dy_in_nm / dx_in_nm;
a1 = 20;
a1y = 10;
radius_x = a1 / (2 * dx_in_nm) * d;
radius_y = a1y / (2 * dx_in_nm) * d;
lc = 0.1;
lc_refine_1 = lc/1;
lc_refine_2 = lc/1;

Point(1) = {-d / 2, dy / 2, 0, lc};
Point(2) = {-d / 2, -dy / 2, 0, lc};
Point(3) = {d / 2, -dy / 2, 0, lc};
Point(4) = {d / 2, dy / 2, 0, lc};

Point(5) = {0, 0, 0, lc_refine_2};
Point(6) = {0, radius_y, 0, lc_refine_1};
Point(7) = {-radius_x, 0, 0, lc_refine_1};
Point(8) = {0, -radius_y, 0, lc_refine_1};
Point(9) = {radius_x, 0, 0, lc_refine_1};

Line(1) = {1, 4};
Line(2) = {4, 3};
Line(3) = {3, 2};
Line(4) = {2, 1};
Ellipsis(5) = {9, 5, 6, 6};
Ellipsis(6) = {6, 5, 7, 7};
Ellipsis(7) = {7, 5, 8, 8};
Ellipsis(8) = {8, 5, 9, 9};

Line Loop(9) = {1, 2, 3, 4};
Line Loop(10) = {5, 6, 7, 8};
Plane Surface(1) = {9, 10};
Plane Surface(2) = {10};

Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3};
Physical Line(4) = {4};
Physical Surface(1) = {1};
Physical Surface(2) = {2};
