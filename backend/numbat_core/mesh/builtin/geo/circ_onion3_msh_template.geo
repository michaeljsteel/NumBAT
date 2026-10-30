// Template mesh geometry file for multiple concentric circular inclusions.

// Force Gmsh to use legacy msh file format v2
Mesh.MshFileVersion = 2.2;

d = 1; // grating period
dx_in_nm = 2000;
dy_in_nm = 2000;
dy = dy_in_nm/dx_in_nm;
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
scale = d/dx_in_nm;

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
radb = dx_in_nm/2 * scale;  // Outer boundary radius

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

// Circle 2 vertices
Point(26) = {x0,      y0+rad2, 0, lc_refine_2};   // N
Point(27) = {x0-rad2, y0,      0, lc_refine_2};   // W
Point(28) = {x0,      y0-rad2, 0, lc_refine_2};   // S
Point(29) = {x0+rad2, y0,      0, lc_refine_2};   // E

// Circle 3 vertices
Point(36) = {x0,      y0+rad3, 0, lc_refine_2};   // N
Point(37) = {x0-rad3, y0,      0, lc_refine_2};   // W
Point(38) = {x0,      y0-rad3, 0, lc_refine_2};   // S
Point(39) = {x0+rad3, y0,      0, lc_refine_2};   // E

//
//Radial lines

// Ring 1 radii
Line(10) = {17,1};   // west
Line(11) = {1,19};   // east
Line(14) = {16,1};   // north
Line(15) = {1,18};   // south


// Ring 2 radii
Line(30) = {17,27}; //west
Line(76) = {19,29}; //east
Line(36) = {16,26}; //north
Line(56) = {18,28}; //south

// Ring 3 radii
Line(29) = {37,27};
Line(89) = {39,29};
Line(49) = {36,26};
Line(69) = {38,28};


// Final radii out to Boundary
Line(90) = {6,36};
Line(91) = {7,37};
Line(92) = {8,38};
Line(93) = {9,39};

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

// Circle 2
Ellipsis(205) = {29,1,26,26};
Ellipsis(206) = {26,1,27,27};
Ellipsis(207) = {27,1,28,28};
Ellipsis(208) = {28,1,29,29};

// Circle 3
Ellipsis(301) = {39,1,36,36};
Ellipsis(302) = {36,1,37,37};
Ellipsis(303) = {37,1,38,38};
Ellipsis(304) = {38,1,39,39};


//
//Surfaces
//

// Background outer region
Line Loop(901) = {101, 91, -302, -90};
Plane Surface(902) = {901};

Line Loop(903) = {-102, 92, -303, -91};
Plane Surface(904) = {903};
Line Loop(905) = {103, 93, -304, -92};
Plane Surface(906) = {905};
Line Loop(907) = {104, 90, -301, -93};
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


//Ring 2
Line Loop(945) = {206, -30, -202, 36};
Plane Surface(946) = {945};
Line Loop(990) = {56, -207, -30, 203};
Plane Surface(1008) = {990};

Line Loop(955) = {201, 36, -205, -76};
Plane Surface(956) = {955};
Line Loop(988) = {204, 76, -208, -56};
Plane Surface(989) = {988};

Physical Surface(3) = {946, 956, 989, 1008};  // Ring 1, mat_b, rad1<r<rad2

//Ring 3
Line Loop(943) = {302, 29, -206, -49};
Plane Surface(944) = {943};
Line Loop(1009) = {69, -207, -29, 303};
Plane Surface(1010) = {1009};

Line Loop(957) = {89, 205, -49, -301};
Plane Surface(958) = {957};
Line Loop(1006) = {89, -208, -69, 304};
Plane Surface(1007) = {1006};

Physical Surface(4) = {944, 958, 1007, 1010}; // Ring 2, mat_c

