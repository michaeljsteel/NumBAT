// Template mesh geometry file for a single suspended inclusion.
// Inclusion can be circular/elliptical (default), or square/rectangular.

// Force Gmsh to use legacy msh file format v2
Mesh.MshFileVersion = 2.2;

d = 1; // grating period
dx_in_nm = 1000;
dy_in_nm = 1000;
dy = dy_in_nm/dx_in_nm;
a1 = 100;
a2 = 100;
a3 = 100;
rad1 = (a1/(2*dx_in_nm))*d;
rad2 = ((a1+a2)/(2*dx_in_nm))*d;
rad3 = ((a1+a2+a3)/(2*dx_in_nm))*d;
lc = 0.1; 
lc_refine_1 = lc/1; // on cylinder surfaces
lc_refine_2 = lc/1; // cylinder1 surfaces

hy = dy; // Thickness: square profile => hy=d
hx = d;

x0 = -dy/2;
y0 = dy/2;


//Outer box corners
Point(1) = {x0,    y0,    0, lc};         // NW
Point(2) = {x0,    y0-hy, 0, lc};         // SW
Point(3) = {x0+hx, y0-hy, 0, lc};         // SE
Point(4) = {x0+hx, y0,    0,lc};          // NE

//Outer box midpoints
Point(6) = {x0+hx/2, y0,      0, lc};     // N
Point(7) = {x0,      y0-hy/2, 0, lc};     // W
Point(8) = {x0+hx/2, y0-hy,   0, lc};     // S
Point(9) = {x0+hx,   y0-hy/2, 0, lc};     // E 

// Circle 1 vertices
Point(10) = {x0+hx/2,      y0-hy/2,      0, lc_refine_2};   // Center
Point(11) = {x0+hx/2,      y0-hy/2+rad1, 0, lc_refine_2};   // N
Point(12) = {x0+hx/2-rad1, y0-hy/2,      0, lc_refine_2};   // W
Point(13) = {x0+hx/2,      y0-hy/2-rad1, 0, lc_refine_2};   // S
Point(14) = {x0+hx/2+rad1, y0-hy/2,      0, lc_refine_2};   // E

// Outer Box
Line(1) = {1,6};  // top left
Line(2) = {6,4};  // top right
Line(3) = {2,8};  // bottom left
Line(4) = {8,3};  // bottom right
Line(5) = {1,7};  // left upper
Line(6) = {7,2};  // left lower
Line(7) = {4,9};  // right upper
Line(8) = {9,3};  // right lower

//Circle 1 radii 
Line(11) = {11,10};  // N
Line(12) = {12,10};  // W
Line(13) = {10,13};  // S
Line(14) = {10,14};  // E

// Lines from outer box to outer circle (1)
Line(2001) = {6,11}; //N
Line(2002) = {7,12}; //W
Line(2003) = {13,8}; //S
Line(2004) = {14,9}; //E

//Circle 1 circumference
Ellipsis(211) = {14,10,11,11};  //NE
Ellipsis(212) = {11,10,12,12};  //NW
Ellipsis(213) = {12,10,13,13};  //SW
Ellipsis(214) = {13,10,14,14};  //SE

//
//Surfaces 
//

// Sectors circle 1
Line Loop(511) = {11, 14, 211};  //NE
Plane Surface(516) = {511};
Line Loop(512) = {212, 12, -11}; //NW 
Plane Surface(517) = {512};
Line Loop(513) = {13, -213, 12}; //SW
Plane Surface(518) = {513};
Line Loop(514) = {14, -214, -13}; //SE
Plane Surface(519) = {514};

// Outer circle to boundary surfaces (Squares with quarter circles removed at the origin)
Line Loop(501) = {2001, -211, 2004, -7, -2}; //NE
Plane Surface(506) = {501};
Line Loop(502) = {212, -2002, -5, 1, 2001};  //NW
Plane Surface(507) = {502};
Line Loop(503) = {2003, -3, -6, 2002, 213};   //SW
Plane Surface(508) = {503};
Line Loop(504) = {2004, 8, -4, -2003, 214};  //SE
Plane Surface(509) = {504};

//
//Physical features
//

// Outer boundary
Physical Line(1) = {1,2};  // Top
Physical Line(2) = {3,4};  // Bottom
Physical Line(3) = {5,6};  // Left
Physical Line(4) = {7,8};  // Right

Physical Surface(1) = {507,506,509,508};     //Outer region, mat_bkg, r>rad_1
Physical Surface(2) = {517, 516, 519, 518};  //Circle 1, mat_a, r<rad_1



