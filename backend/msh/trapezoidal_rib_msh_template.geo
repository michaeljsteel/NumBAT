// Template mesh geometry file for a slot waveguide.

// Force Gmsh to use legacy msh file format v2
Mesh.MshFileVersion = 2.2;

//Unscaled dimensions
dx_in_nm = 4000.0; // total width
dy_in_nm = 2000.0; // total height


rib_top_width = 600.0;
rib_base_width = 900.0;
rib_height = 500.0;

slab_thickness = 300.0;
slab_width = 1800.0;

//Scaled dimensions
dy = dy_in_nm/dx_in_nm;  // height (fraction or width)
hy2 = -.65*dy; // top of the substrate
hy2 = -.5*dy; // top of the substrate

x0 = -.5;
y0 = .5*dy;

//Mesh
lc = 0.020000; // background and unitcell edge
lc_refine_1 = lc/10.0; // rib
lc_refine_2 = lc/5.0; // slab

//////////////////////////////////////////


Point(1) = {0+x0, 0+y0, 0, lc};
Point(2) = {0+x0, -dy+y0, 0, lc};
Point(3) = {1+x0, -dy+y0, 0, lc};
Point(4) = {1+x0, 0+y0, 0, lc};

// slab (buried part of the rib)
Point(5) = {.5+x0-slab_width/dx_in_nm/2., hy2+y0, 0, lc_refine_2};
Point(6) = {.5+x0+slab_width/dx_in_nm/2., hy2+y0, 0, lc_refine_2};
Point(13) = {.5+x0-slab_width/dx_in_nm/2., hy2+y0-slab_thickness/dx_in_nm, 0, lc_refine_2};
Point(14) = {.5+x0+slab_width/dx_in_nm/2., hy2+y0-slab_thickness/dx_in_nm, 0, lc_refine_2};

// sticking out part of the rib
Point(20) = {.5+x0-rib_base_width/dx_in_nm/2., hy2+y0, 0, lc_refine_1};
Point(21) = {.5+x0+rib_base_width/dx_in_nm/2., hy2+y0, 0, lc_refine_1};
Point(22) = {.5+x0-rib_top_width/dx_in_nm/2., hy2+y0+rib_height/dx_in_nm, 0, lc_refine_1};
Point(23) = {.5+x0+rib_top_width/dx_in_nm/2., hy2+y0+rib_height/dx_in_nm, 0, lc_refine_1};

Point(11) = {0+x0, hy2+y0, 0, lc};
Point(12) = {1+x0, hy2+y0, 0, lc};

Point(24) = {.5+x0+rib_top_width/dx_in_nm/2., y0, 0, lc};
Point(25) = {.5+x0-rib_top_width/dx_in_nm/2., y0, 0, lc};


Line(1) = {2, 3};
Line(2) = {24, 25};
Line(3) = {22, 23};
Line(4) = {20, 21};
Line(11) = {6, 14};
Line(12) = {14, 13};
Line(13) = {13, 5};
Line(14) = {1, 11};
Line(15) = {11, 2};
Line(16) = {4, 12};
Line(17) = {12, 3};
Line(19) = {11, 5};
Line(20) = {6, 12};
Line(47) = {1, 25};
Line(49) = {25, 22};
Line(51) = {22, 20};
Line(53) = {20, 5};
Line(55) = {21, 6};
Line(56) = {21, 23};
Line(57) = {23, 24};
Line(59) = {24, 4};


// elevated part of the rib
Line Loop(1) = {56, -3, 51, 4};
Plane Surface(1) = {1};

// buried part of the rib
Line Loop(2) = {-12, -11, -55, -4, 53, -13};
Plane Surface(2) = {2};

// entire bottom substrate
Line Loop(3) = {-13, -12, -11, 20, 17, -1, -15, 19};
Plane Surface(3) = {3};

// surrounding - left
Line Loop(4) = {19, -53, -51, -49, -47, 14};
Plane Surface(4) = {4};

// surrounding - top
Line Loop(5) = {3, 57, 2, 49};
Plane Surface(5) = {5};

// surrounding - right
Line Loop(6) = {55, 20, -16, -59, -57, -56};
Plane Surface(6) = {6};

// surroundings
Physical Surface(1) = {4, 5, 6};

// rib elevated
Physical Surface(2) = {1};

// rib buried
Physical Surface(3) = {2};

// substrate
Physical Surface(4) = {3};
