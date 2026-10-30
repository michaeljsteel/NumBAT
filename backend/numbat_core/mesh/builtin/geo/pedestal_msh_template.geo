// Template mesh geometry file for a pedestal waveguide.

// Force Gmsh to use legacy msh file format v2
Mesh.MshFileVersion = 2.2;

d = 1; // grating period
dx_in_nm = 100;
dy_in_nm = 50;
dy = dy_in_nm/dx_in_nm;

a1 = 20;    // ped base width
a1top = 15; //ped top width
a1y = 10;   //ped height

pedestal_whalf = (a1/(2*dx_in_nm))*d;
pedestal_hhalf = (a1y/(2*dx_in_nm))*d;
pedestal_wtophalf = (a1top/(2*dx_in_nm))*d;

slabx = 80;
slaby = 10;
slab_w = slabx/dx_in_nm;
slab_h = slaby/dx_in_nm;

px = 2;  //pillar width
py = 5;  //pillar height
p_w = px/dx_in_nm;
p_h = py/dx_in_nm;

lc = 0.1; // background and unitcell edge
lc_refine_1 = lc/1; // pedestal
lc_refine_2 = lc/1; // slab

hy = dy/2 + (slab_h/2) + pedestal_hhalf; // 
hx = 0.;

//Border
Point(1) = {-d/2,    dy/2, 0, lc};     // NW
Point(2) = {-hx-d/2, -dy/2, 0, lc};    // SW
Point(3) = {-hx+d/2, -dy/2, 0, lc};    // SE
Point(4) = {d/2,     dy/2, 0,lc};      // NE

// Slab
Point(5) = {-slab_w/2, -hy+slab_h+dy/2, 0, lc_refine_2};
Point(6) = {+slab_w/2, -hy+slab_h+dy/2, 0, lc_refine_2};
Point(13) = {-slab_w/2, -hy+dy/2,       0, lc_refine_2};
Point(14) = {+slab_w/2, -hy+dy/2,       0, lc_refine_2};

// Supported pedestal
Point(7) = {-hx-pedestal_whalf, -hy+slab_h+p_h+dy/2, 0, lc_refine_1};
Point(8) = {-hx+pedestal_whalf, -hy+slab_h+p_h+dy/2, 0, lc_refine_1};
Point(9) = {-hx-pedestal_wtophalf, -hy+2*pedestal_hhalf+slab_h+p_h+dy/2, 0, lc_refine_1};
Point(10) = {-hx+pedestal_wtophalf, -hy+2*pedestal_hhalf+slab_h+p_h+dy/2, 0, lc_refine_1};

// Pillar
Point(19) = {-hx-p_w, -hy+slab_h+dy/2, 0, lc_refine_2};
Point(20) = {-hx+p_w, -hy+slab_h+dy/2, 0, lc_refine_2};
Point(21) = {-hx-p_w, -hy+slab_h+p_h+dy/2, 0, lc_refine_1};
Point(22) = {-hx+p_w, -hy+slab_h+p_h+dy/2, 0, lc_refine_1};

// Pillar surrounds (air)
Point(23) = {-hx-pedestal_whalf, -hy+slab_h+dy/2, 0, lc_refine_1};
Point(24) = {-hx+pedestal_whalf, -hy+slab_h+dy/2, 0, lc_refine_1};

// X-axis
Point(11) = {-d/2,       -hy+slab_h+dy/2, 0, lc};
Point(12) = {d/2,        -hy+slab_h+dy/2, 0, lc};

Point(15) = {-hx+pedestal_whalf, dy/2, 0, lc};
Point(16) = {-hx-pedestal_whalf, dy/2, 0, lc};
Point(17) = {-hx+pedestal_whalf, -dy/2, 0, lc};
Point(18) = {-hx-pedestal_whalf, -dy/2, 0, lc};

Line(6) = {7, 9};
Line(8) = {10, 8};
Line(11) = {6, 14};
Line(12) = {14, 13};
Line(13) = {13, 5};
Line(14) = {1, 11};
Line(15) = {11, 2};
Line(16) = {4, 12};
Line(17) = {12, 3};
Line(19) = {11, 5};
Line(20) = {6, 12};
Line(32) = {1, 16};
Line(33) = {16, 9};
Line(34) = {10, 15};
Line(35) = {15, 16};
Line(36) = {15, 4};
Line(44) = {2, 18};
Line(45) = {18, 17};
Line(46) = {17, 3};


Physical Line(29) = {14, 15};
Physical Line(31) = {17, 16};
Physical Line(43) = {32, 35, 36};
Physical Line(47) = {44, 45, 46};
Line(48) = {9, 10};
Line(49) = {8, 22};
Line(50) = {21, 22};
Line(51) = {22, 20};
Line(52) = {20, 19};
Line(53) = {19, 21};
Line(54) = {21, 7};
Line(55) = {23, 19};
Line(56) = {7, 23};
Line(57) = {23, 5};
Line(58) = {20, 24};
Line(59) = {24, 8};
Line(60) = {24, 6};
Line Loop(61) = {32, 33, -6, 56, 57, -19, -14};
Plane Surface(62) = {61};
Line Loop(63) = {33, 48, 34, 35};
Plane Surface(64) = {63};
Line Loop(65) = {48, 8, 49, -50, 54, 6};
Plane Surface(66) = {65};
Line Loop(67) = {49, 51, 58, 59};
Plane Surface(68) = {67};
Line Loop(69) = {50, 51, 52, 53};
Plane Surface(70) = {69};
Line Loop(71) = {53, 54, 56, 55};
Plane Surface(72) = {71};
Line Loop(73) = {55, -52, 58, 60, 11, 12, 13, -57};
Plane Surface(74) = {73};
Line Loop(75) = {19, -13, -12, -11, 20, 17, -46, -45, -44, -15};
Plane Surface(76) = {75};
Line Loop(77) = {16, -20, -60, 59, -8, 34, 36};
Plane Surface(78) = {77};
Physical Surface(1) = {62, 64, 78, 76, 68, 72};  // bkg
Physical Surface(2) = {66};  // pedestal
Physical Surface(3) = {74};  // slab
Physical Surface(4) = {70};  // pillar
