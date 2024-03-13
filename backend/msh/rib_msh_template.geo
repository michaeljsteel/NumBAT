// Template mesh geometry file for a rib waveguide.

// Force Gmsh to use legacy msh file format v2
Mesh.MshFileVersion = 2.2;

d = 1; // grating period
dx_in_nm = 100;
dy_in_nm = 50;
dy = dy_in_nm/dx_in_nm;
a1 = 20;
a1y = 10;
rib_halfw = (a1/(2*dx_in_nm))*d;
rib_halfh = (a1y/(2*dx_in_nm))*d;

slabx = 80;
slaby = 10;
slab_w = slabx/dx_in_nm;
slab_h = slaby/dx_in_nm;

lc = 0.1; // background and unitcell edge
lc_refine_1 = lc/1; // rib
lc_refine_2 = lc/1; // slab

hy = dy/2 + (slab_h/2) + rib_halfh; // 
y0 = hy-slab_h;
ytop = y0;
ybot = y0-dy;
ymid = y0-hy;


// Outer box
Point(1) = {-d/2, ytop,    0, lc};     // NW
Point(2) = {-d/2, ybot, 0, lc};     // SW
Point(3) = {d/2,  ybot, 0, lc};     // SE
Point(4) = {d/2,  ytop,    0,lc};      // NE

// Slab
Point(5)  = {-slab_w/2, ymid+slab_h, 0, lc_refine_2}; // NW
Point(6)  = {+slab_w/2, ymid+slab_h, 0, lc_refine_2}; // NE
Point(13) = {-slab_w/2, ymid,        0, lc_refine_2}; // SW
Point(14) = {+slab_w/2, ymid,        0, lc_refine_2}; // SE

// Rib
Point(7)  = {-rib_halfw, ymid+slab_h, 0, lc_refine_1};             // SW
Point(8)  = {+rib_halfw, ymid+slab_h, 0, lc_refine_1};             // SE
Point(9)  = {-rib_halfw, ymid+2*rib_halfh+slab_h, 0, lc_refine_1}; // NW
Point(10) = {+rib_halfw, ymid+2*rib_halfh+slab_h, 0, lc_refine_1}; // NE


Point(11) = {-d/2,       ymid+slab_h, 0, lc}; // Left end of surface
Point(12) = {d/2,        ymid+slab_h, 0, lc}; // Right end of surface
Point(15) = {+rib_halfw, ytop, 0, lc};           // Intersec of top and right rib vertical
Point(16) = {-rib_halfw, ytop, 0, lc};           // Intersec of top and left rib vertical
Point(17) = {+rib_halfw, ybot, 0, lc};           // Intersec of bot and right rib vertical (no purpose)
Point(18) = {-rib_halfw, ybot, 0, lc};           // Intersec of bot and right rib vertical (no purpose)

Line(2) = {2,3};
Line(4) = {4,1};
Line(5) = {5, 7};
Line(6) = {7, 9};
Line(7) = {9, 10};
Line(8) = {10, 8};
Line(9) = {8, 7};
Line(10) = {8, 6};
Line(11) = {6, 14};
Line(12) = {14, 13};
Line(13) = {13, 5};
Line(14) = {1, 11};
Line(15) = {11, 2};
Line(16) = {4, 12};
Line(17) = {12, 3};
Line Loop(18) = {14, 15, 2, -17, -16, 4}; //Lines 2 and 4 don't exist?!
Line(19) = {11, 5};
Line(20) = {6, 12};
Line Loop(25) = {6, 7, 8, 9};
Plane Surface(26) = {25};
Line Loop(27) = {9, -5, -13, -12, -11, -10};
Plane Surface(28) = {27};
Line(32) = {1, 16};
Line(33) = {16, 9};
Line(34) = {10, 15};
Line(35) = {15, 16};
Line(36) = {15, 4};
Line Loop(37) = {5, 6, -33, -32, 14, 19};
Plane Surface(38) = {37};
Line Loop(39) = {7, 34, 35, 33};
Plane Surface(40) = {39};
Line Loop(41) = {34, 36, 16, -20, -10, -8};
Plane Surface(42) = {41};
Line(44) = {2, 18};
Line(45) = {18, 17};
Line(46) = {17, 3};
Line Loop(48) = {19, -13, -12, -11, 20, 17, -46, -45, -44, -15};
Plane Surface(49) = {48};


Physical Line(29) = {14, 15};
Physical Line(31) = {17, 16};
Physical Line(43) = {32, 35, 36};
Physical Line(47) = {44, 45, 46};

Physical Surface(1) = {38, 40, 42, 49};
Physical Surface(2) = {26};
Physical Surface(3) = {28};

