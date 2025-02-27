// Template mesh geometry file for a slot waveguide.

// Force Gmsh to use legacy msh file format v2
Mesh.MshFileVersion = 2.2;

//d = 1; // grating period
dx_in_nm = 1000.0; // total width
dy_in_nm = 1000.0; // total height

dy = dy_in_nm/dx_in_nm;  // height (fraction or width)

scale = 1.0/dx_in_nm;

//hy2 = -.65*dy; // top of the substrate

base_width = 600.0;   // width of horizontal base
peak_xoff = 200.0;    // horizontal offset of peak from left end of base
peak_height = 400.0;  // perpendicular height of peak from base

bw = base_width*scale;
pxo = peak_xoff*scale;
ph = peak_height*scale;

lc = 0.1;         // background and unitcell edge
lc_refine_1 = lc/1.0;  // triangle boundaries
lc_refine_2 = lc/1.0;  // peak corner (which should be the sharp one if there is one)

//////////////////////////////////////////

boxx=1;
boxy=dy_in_nm*scale;

//
// Points
// 

// Outer box 
Point(1) = {-boxx/2.,  boxy/2.,  0, lc};    //NW
Point(2) = { boxx/2.,  boxy/2.,  0, lc};    //NE
Point(3) = {-boxx/2., -boxy/2.,  0, lc};    //SW
Point(4) = { boxx/2., -boxy/2.,  0, lc};    //SE

Point(5) = { -bw/2+pxo,  boxy/2., 0, lc};     // Boundary above peak
Point(6) = {-boxx/2,     -ph/2., 0, lc};           // Boundary left of base
Point(7) = { boxx/2,     -ph/2., 0, lc};           // Boundary right of base

// Triangle 
Point(10) = {-bw/2.,       -ph/2., 0, lc_refine_1}; // Left base
Point(11) = {bw/2.,        -ph/2., 0, lc_refine_1}; // Right base
Point(12) = {-bw/2.+pxo,    ph/2., 0, lc_refine_2}; // Peak


//
// Regions
// 

// Triangle 
Line(11) = {10, 11};
Line(12) = {11, 12};
Line(13) = {12, 10};
Line Loop(14) = {11, 12, 13};
Plane Surface(15) = {14};

// Top left background
Line(21) = {1, 5};
Line(22) = {5, 12};
//Line(23) = {12, 10};
Line(24) = {10, 6};
Line(25) = {6, 1};
Line Loop(26) = {-25, -24, -13, -22, -21};
Plane Surface(27) = {26};



// Top right background
Line(31) = {5, 2};
Line(32) = {2, 7};
Line(33) = {7, 11};
//Line(34) = {11, 12};
//Line(35) = {12, 5};
Line Loop(36) = {22, -12, -33, -32, -31};
Plane Surface(37) = {36};

// Bottom background
//Line(41) = {6, 10};
Line(42) = {7, 4};
Line(43) = {4, 3};
Line(44) = {3, 6};
Line Loop(46) = {24, -44, -43, -42, 33, -11};
Plane Surface(47) = {46};



// Final physical structures

// Outer boundary
Physical Line(1) = {21,31};  //Top 
Physical Line(2) = {32,42};  //Right
Physical Line(3) = {43};     //Bottom
Physical Line(4) = {44,25};  //Left 

Physical Surface(1) = {27,37,47};  //Background
Physical Surface(2) = {15};        //Triangle

