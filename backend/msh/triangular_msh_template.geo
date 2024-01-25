// Template mesh geometry file for a slot waveguide.

// Force Gmsh to use legacy msh file format v2
Mesh.MshFileVersion = 2.2;

d = 1; // grating period
//d_in_nm = 8200.0; // total width
//dy_in_nm = 2500.0; // total height

d_in_nm = 1000.0; // total width
dy_in_nm = 1000.0; // total height

dy = dy_in_nm/d_in_nm;  // height (fraction or width)

scale = d/d_in_nm;

//hy2 = -.65*dy; // top of the substrate

base_width = 600.0;   // width of horizontal base
peak_xoff = 200.0;    // horizontal offset of peak from left end of base
peak_height = 400.0;  // perpendicular height of peak from base

bw = base_width*scale;
pxo = peak_xoff*scale;
ph = peak_height*scale;

lc = 0.020000;         // background and unitcell edge
lc_refine_1 = lc/1.0; // triangle boundaries
//lc_refine_2 = lc/1.0;  // 

//////////////////////////////////////////

boxside=1.0;

//
// Points
// 

// Outer box 
Point(1) = {-boxside/2.,  boxside/2.,  0, lc};    //NW
Point(2) = { boxside/2.,  boxside/2.,  0, lc};    //NE
Point(3) = {-boxside/2., -boxside/2.,  0, lc};    //SW
Point(4) = { boxside/2., -boxside/2.,  0, lc};    //SE

Point(5) = { -bw/2+pxo,  boxside/2., 0, lc};     // Boundary above peak
Point(6) = {-boxside/2, -ph/2., 0, lc};           // Boundary left of base
Point(7) = { boxside/2, -ph/2., 0, lc};           // Boundary right of base

// Triangle 
Point(10) = {-bw/2.,       -ph/2., 0, lc_refine_1}; // Left base
Point(11) = {bw/2.,        -ph/2., 0, lc_refine_1}; // Right base
Point(12) = {-bw/2.+pxo,    ph/2., 0, lc_refine_1}; // Peak


//
// Regions
// 

// Triangle 
Line(11) = {10, 11};
Line(12) = {11, 12};
Line(13) = {12, 10};
Line Loop(11) = {11, 12, 13};
Plane Surface(12) = {11};

// Top left background
Line(21) = {1, 5};
Line(22) = {5, 12};
//Line(23) = {12, 10};
Line(24) = {10, 6};
Line(25) = {6, 1};
Line Loop(26) = {21, 22, 13, 24, 25};
Plane Surface(27) = {26};



// Top right background
Line(31) = {5, 2};
Line(32) = {2, 7};
Line(33) = {7, 11};
//Line(34) = {11, 12};
//Line(35) = {12, 5};
Line Loop(36) = {31, 32, 33, 12, -22};
Plane Surface(37) = {36};

// Bottom background
Line(41) = {6, 10};
Line(42) = {7, 4};
Line(43) = {4, 3};
Line(44) = {3, 6};
Line Loop(46) = {41, 11, -33, 42, 43, 44};
Plane Surface(47) = {46};



// Final physical structures

// Outer boundary
Physical Line(1) = {21,31};  //Top 
Physical Line(2) = {32,42};  //Right
Physical Line(3) = {43};     //Bottom
Physical Line(4) = {44,25};  //Left 

Physical Surface(1) = {27,37,47};  //Background
Physical Surface(2) = {12};        //Triangle

