// Template mesh geometry file for multiple concentric circular inclusions.

// Force Gmsh to use legacy msh file format v2
Mesh.MshFileVersion = 2.2;

d = 1; // grating period
d_in_nm = 2000;
dy_in_nm = 2000;
dy = dy_in_nm/d_in_nm;
a1 = 100;
a2 = 100;
a3 = 100;
a4 = 100;
a5 = 100;
a6 = 100;
a7 = 100;
a8 = 100;
a9 = 100;
a10 = 100;
a11 = 100;
a12 = 100;
a13 = 100;
a14 = 100;
a15 = 100;
a1r = a1/2; // Radius of central rod
scale = d/d_in_nm;

rad1 = a1r * scale;
rad2 = (a1r+a2)*scale;
rad3 = (a1r+a2+a3)*scale;
rad4 = (a1r+a2+a3+a4)*scale;
rad5 = (a1r+a2+a3+a4+a5)*scale;
rad6 = (a1r+a2+a3+a4+a5+a6)*scale;
rad7 = (a1r+a2+a3+a4+a5+a6+a7)*scale;
rad8 = (a1r+a2+a3+a4+a5+a6+a7+a8)*scale;
rad9 = (a1r+a2+a3+a4+a5+a6+a7+a8+a9)*scale;
rad10 = (a1r+a2+a3+a4+a5+a6+a7+a8+a9+a10)*scale;
rad11 = (a1r+a2+a3+a4+a5+a6+a7+a8+a9+a10+a11)*scale;
rad12 = (a1r+a2+a3+a4+a5+a6+a7+a8+a9+a10+a11+a12)*scale;
rad13 = (a1r+a2+a3+a4+a5+a6+a7+a8+a9+a10+a11+a12+a13)*scale;
rad14 = (a1r+a2+a3+a4+a5+a6+a7+a8+a9+a10+a11+a12+a13+a14)*scale;
rad15 = (a1r+a2+a3+a4+a5+a6+a7+a8+a9+a10+a11+a12+a13+a14+a15)*scale;

lc = 0.1;
lc_refine_1 = lc/1; // on cylinder surfaces
lc_refine_2 = lc/1; // cylinder1 centres

//hy = dy; // Thickness: square profile => hy=d
//hx = d;
radb = d_in_nm * scale;  // Outer boundary radius

// Origin coords
x0 = 0.0;
y0 = 0.0;

// Origin
Point(1) = {x0,      y0,      0, lc_refine_2};   // Center

//Boundary vertices
Point(6) = {x0,      y0+radb, 0, lc};   // N
Point(7) = {x0-radb, y0,      0, lc};   // W
Point(8) = {x0,      y0-radb, 0, lc};   // S
Point(9) = {x0+radb, y0,      0, lc};   // E


// Circle 1 vertices
Point(16) = {x0,      y0+rad1, 0, lc_refine_2};   // N
Point(17) = {x0-rad1, y0,      0, lc_refine_2};   // W
Point(18) = {x0,      y0-rad1, 0, lc_refine_2};   // S
Point(19) = {x0+rad1, y0,      0, lc_refine_2};   // E

//
//Radial lines

// Ring 1 radii
Line(10) = {17,1};   // west
Line(11) = {1,19};   // east
Line(14) = {16,1};   // north
Line(15) = {1,18};   // south



// Final radii out to Boundary
Line(90) = {6,16};
Line(91) = {7,17};
Line(92) = {8,18};
Line(93) = {9,19};

//
//Ellipse arcs
//

// Outer boundary
Ellipsis(101) = {6,1,7,7};
Ellipsis(102) = {8,1,7,7};
Ellipsis(103) = {8,1,9,9};
Ellipsis(104) = {9,1,6,6};

// Circle 1
Ellipsis(201) = {19,1,16,16};
Ellipsis(202) = {16,1,17,17};
Ellipsis(203) = {17,1,18,18};
Ellipsis(204) = {18,1,19,19};

//
//Surfaces
//

// Background outer region
Line Loop(901) = {101, 91, -202, -90};
Plane Surface(902) = {901};

Line Loop(903) = {-102, 92, -203, -91};
Plane Surface(904) = {903};
Line Loop(905) = {103, 93, -204, -92};
Plane Surface(906) = {905};
Line Loop(907) = {104, 90, -201, -93};
Plane Surface(908) = {907};

Physical Surface(1) = {902,904,906,908};      // Outer space, mat_bkg


// Inner rings

//Ring 1
Line Loop(947) = {202, 10, -14};
Plane Surface(948) = {947};
Line Loop(953) = {15, -203, 10};
Plane Surface(954) = {953};

Line Loop(949) = {14, 11, 201};
Plane Surface(950) = {949};
Line Loop(951) = {11, -204, -15};
Plane Surface(952) = {951};

Physical Surface(2) = {948, 950, 952, 954};   // Center, mat_a, r<rad1



