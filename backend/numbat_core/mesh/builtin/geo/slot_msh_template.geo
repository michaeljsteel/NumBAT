// Template mesh geometry file for a slot waveguide.

// Force Gmsh to use legacy msh file format v2
Mesh.MshFileVersion = 2.2;

d = 1; // grating period
dx_in_nm = 1000;
dy_in_nm = 500;
dy = dy_in_nm/dx_in_nm;

slot_w = 200;
rib_h = 100;
slot_hwidth = (slot_w/(2*dx_in_nm))*d;
rib_height = (rib_h/(dx_in_nm))*d;

rib_w = 200;
rib_width = (rib_w/(dx_in_nm))*d;

slab_w = 800;
slab_h = 100;
slab_width = slab_w/dx_in_nm;
slab_height = slab_h/dx_in_nm;

lc_bkg = 0.1; // background and unitcell edge
lc_refine_1 = lc_bkg/1; // rib
lc_refine_2 = lc_bkg/1; // slab

hy = dy/2 + (slab_height/2) + rib_height/2; // 
x0 = -d/2;
x1 = x0+d;
xmid = x0+d/2;

y0 = hy-slab_height;


Point(1) = {x0, y0, 0, lc_bkg};
Point(2) = {x0, y0-dy, 0, lc_bkg};
Point(3) = {x0+d, y0-dy, 0, lc_bkg};
Point(4) = {x0+d, y0, 0,lc_bkg};

// Slab
Point(5) = {x0+d/2-slab_width/2, y0-hy+slab_height, 0, lc_refine_2};
Point(6) = {x0+d/2+slab_width/2, y0-hy+slab_height, 0, lc_refine_2};
Point(13) = {x0+d/2-slab_width/2,y0 -hy, 0, lc_refine_2};
Point(14) = {x0+d/2+slab_width/2, y0-hy, 0, lc_refine_2};

// Slot 
Point(7) = {x0+d/2-slot_hwidth, y0-hy+slab_height, 0, lc_refine_1};
Point(8) = {x0+d/2+slot_hwidth, y0-hy+slab_height, 0, lc_refine_1};
Point(9) = {x0+d/2-slot_hwidth, y0-hy+rib_height+slab_height, 0, lc_refine_1};
Point(10) = {x0+d/2+slot_hwidth, y0-hy+rib_height+slab_height, 0, lc_refine_1};

// Left rib left wall
Point(20) = {x0+d/2-slot_hwidth-rib_width, y0-hy+slab_height, 0, lc_refine_1};
Point(22) = {x0+d/2-slot_hwidth-rib_width, y0-hy+rib_height+slab_height, 0, lc_refine_1};

// Right rib right wall
Point(21) = {x0+d/2+slot_hwidth+rib_width, y0-hy+slab_height, 0, lc_refine_1};
Point(23) = {x0+d/2+slot_hwidth+rib_width, y0-hy+rib_height+slab_height, 0, lc_refine_1};

Point(11) = {x0+0, y0-hy+slab_height, 0, lc_bkg};
Point(12) = {x0+d, y0-hy+slab_height, 0, lc_bkg};
Point(15) = {x0+d/2+slot_hwidth, y0, 0, lc_bkg};
Point(16) = {x0+d/2-slot_hwidth, y0, 0, lc_bkg};
Point(17) = {x0+d/2+slot_hwidth, y0-dy, 0, lc_bkg};
Point(18) = {x0+d/2-slot_hwidth, y0-dy, 0, lc_bkg};
Point(24) = {x0+d/2+slot_hwidth+rib_width, y0, 0, lc_bkg};
Point(25) = {x0+d/2-slot_hwidth-rib_width, y0, 0, lc_bkg};
Point(26) = {x0+d/2+slot_hwidth+rib_width, y0-dy, 0, lc_bkg};
Point(27) = {x0+d/2-slot_hwidth-rib_width, y0-dy, 0, lc_bkg};


Line(6) = {7, 9};
Line(7) = {9, 10};
Line(8) = {10, 8};
Line(9) = {8, 7};
Line(11) = {6, 14};
Line(12) = {14, 13};
Line(13) = {13, 5};
Line(14) = {1, 11};
Line(15) = {11, 2};
Line(16) = {4, 12};
Line(17) = {12, 3};
Line(19) = {11, 5};
Line(20) = {6, 12};
Line(33) = {16, 9};
Line(34) = {10, 15};
Line(35) = {15, 16};
Line(45) = {18, 17};
Line(47) = {1, 25};
Line(48) = {25, 16};
Line(49) = {25, 22};
Line(50) = {22, 9};
Line(51) = {22, 20};
Line(52) = {20, 7};
Line(53) = {20, 5};
Line(54) = {8, 21};
Line(55) = {21, 6};
Line(56) = {21, 23};
Line(57) = {23, 24};
Line(58) = {15, 24};
Line(59) = {24, 4};
Line(60) = {10, 23};
Line(85) = {2, 27};
Line(86) = {27, 18};
Line(87) = {17, 26};
Line(88) = {26, 3};
Line Loop(61) = {47, 49, 51, 53, -19, -14};
Plane Surface(62) = {61};
Line Loop(89) = {19, -13, -12, -11, 20, 17, -88, -87, -45, -86, -85, -15};
Plane Surface(90) = {89};
Line Loop(65) = {53, -13, -12, -11, -55, -54, 9, -52};
Plane Surface(66) = {65};
Line Loop(67) = {52, 6, -50, 51};
Plane Surface(68) = {67};
Line Loop(69) = {50, -33, -48, 49};
Plane Surface(70) = {69};
Line Loop(71) = {33, 7, 34, 35};
Plane Surface(72) = {71};
Line Loop(73) = {7, 8, 9, 6};
Plane Surface(74) = {73};
Line Loop(75) = {8, 54, 56, -60};
Plane Surface(76) = {75};
Line Loop(77) = {55, 20, -16, -59, -57, -56};
Plane Surface(78) = {77};
Line Loop(79) = {57, -58, -34, 60};
Plane Surface(80) = {79};
Physical Line(81) = {47, 48, 35, 58, 59};
Physical Line(82) = {16, 17};
Physical Line(84) = {14, 15};
Physical Line(91) = {85, 86, 45, 87, 88};
Physical Surface(1) = {62, 70, 72, 80, 78, 64, 90};
Physical Surface(2) = {74};
Physical Surface(3) = {66};
Physical Surface(4) = {68, 76};
