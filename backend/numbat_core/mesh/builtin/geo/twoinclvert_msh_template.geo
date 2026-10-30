// Template mesh geometry file for two suspended inclusions.
// Inclusions can be circular/elliptical (default), or square/rectangular.

// Force Gmsh to use legacy msh file format v2
Mesh.MshFileVersion = 2.2;

dx_in_nm = 100.0;
dy_in_nm = 50.0;
dy = dy_in_nm/dx_in_nm;

inc_a_w = 20.0;
inc_a_h = 10.0;
inc_b_w = 20.0;
inc_b_h = 10.0;
inc_sep_x = 5.0;
inc_sep_y = 15.0;

halfwid1 = (inc_a_w/(2*dx_in_nm));
halfhgt1 = (inc_a_h/(2*dx_in_nm));

halfwid2 = (inc_b_w/(2*dx_in_nm));
halfhgt2 = (inc_b_h/(2*dx_in_nm));

yoff = inc_sep_y/(2*dx_in_nm);
xoff = inc_sep_x/(2*dx_in_nm);

//rect = 1;

lc_bkg = 0.05;           // boundary
lc_refine_1 = lc_bkg/2.0; // on cylinder surfaces
lc_refine_2 = lc_bkg/2.0; // cylinder centres

hy = dy; // Thickness: square profile => hy=d
x0=-.5;
x1=.5;
y0=hy/2;
y1=-hy/2;
xm=(x0+x1)/2;
ym=(y0+y1)/2;

// Boundary box

Point(1) = {x0, y0, 0, lc_bkg};    //NW
Point(2) = {x0, y1, 0, lc_bkg};    //SW
Point(3) = {x1, y1, 0, lc_bkg};    //SE
Point(4) = {x1, y0, 0, lc_bkg};    //NE

box1x=xm-xoff;
box1y=ym+yoff;
box1e=box1x+halfwid1;
box1w=box1x-halfwid1;
box1n=box1y+halfhgt1;
box1s=box1y-halfhgt1;

box2x=xm+xoff;
box2y=ym-yoff;
box2e=box2x+halfwid1;
box2w=box2x-halfwid2;
box2n=box2y+halfhgt2;
box2s=box2y-halfhgt2;

// Upper box
Point(101) = {x0, box1y, 0, lc_bkg};          // left wall
Point(102) = {x1, box1y, 0, lc_bkg};          // right wall
Point(103) = {box1x, y0, 0, lc_bkg};          // top wall

Point(104) = {box1x, box1y, 0, lc_refine_2};  // centre
Point(105) = {box1x, box1n, 0, lc_refine_2};  // N
Point(106) = {box1w, box1n, 0, lc_refine_2};  // NW
Point(107) = {box1w, box1y, 0, lc_refine_2};  // W
Point(108) = {box1w, box1s, 0, lc_refine_2};  // SW
Point(109) = {box1x, box1s, 0, lc_refine_2};  // S
Point(110) = {box1e, box1s, 0, lc_refine_2};  // SE
Point(111) = {box1e, box1y, 0, lc_refine_2};  // E
Point(112) = {box1e, box1n, 0, lc_refine_2};  // NE

// Lower box
Point(201) = {x0, box2y, 0, lc_bkg};          // left wall
Point(202) = {x1, box2y, 0, lc_bkg};          // right wall
Point(203) = {box2x, y1, 0, lc_bkg};          // bottom wall

Point(204) = {box2x, box2y, 0, lc_refine_2};  // centre
Point(205) = {box2x, box2n, 0, lc_refine_2};  // N
Point(206) = {box2w, box2n, 0, lc_refine_2};  // NW
Point(207) = {box2w, box2y, 0, lc_refine_2};  // W
Point(208) = {box2w, box2s, 0, lc_refine_2};  // SW
Point(209) = {box2x, box2s, 0, lc_refine_2};  // S
Point(210) = {box2e, box2s, 0, lc_refine_2};  // SE
Point(211) = {box2e, box2y, 0, lc_refine_2};  // E
Point(212) = {box2e, box2n, 0, lc_refine_2};  // NE


// Outer box
Line(301) = {103,1};   // N border
Line(302) = {4,103};   // ...
Line(303) = {1,101};   // W border
Line(304) = {101,201}; // ...
Line(305) = {201,2};   // ...
Line(306) = {2,203};   // S border
Line(307) = {203,3};   // ...
Line(308) = {3,202};   // E border
Line(309) = {202,102}; // ...
Line(310) = {102,4};   // ...


// Upper box lines
Line(311) = {101,107};  //W border
Line(312) = {111,102};  //E border
Line(313) = {103,105};  //N border

Line(314) = {104,105};  //0->N
Line(315) = {104,107};  //0->W
Line(316) = {104,109};  //0->S
Line(317) = {104,111};  //0->E

Line(318) = {105,106};  //outer going anticlockwise
Line(319) = {106,107};
Line(320) = {107,108};
Line(321) = {108,109};
Line(322) = {109,110};
Line(323) = {110,111};
Line(324) = {111,112};
Line(325) = {112,105};


// Lower box lines
Line(341) = {201,207};  //W border
Line(342) = {211,202};  //E border
Line(343) = {203,209};  //S border

Line(344) = {204,205};  //0->N
Line(345) = {204,207};  //0->W
Line(346) = {204,209};  //0->S
Line(347) = {204,211};  //0->E

Line(348) = {205,206};  //outer going anticlockwise
Line(349) = {206,207};
Line(350) = {207,208};
Line(351) = {208,209};
Line(352) = {209,210};
Line(353) = {210,211};
Line(354) = {211,212};
Line(355) = {212,205};


// All loops defined counterclockwise

//Background regions
Line Loop(40) = {303,311,-319,-318,-313,301};  //top left
Plane Surface(41) = {40};
Line Loop(42) = {313,-325,-324,312,310,302};  //top right
Plane Surface(43) = {42};
Line Loop(44) = {304,341,-349,-348,-355,-354,342,309,-312,-323,-322,-321,-320,-311};  //middle
Plane Surface(45) = {44};
Line Loop(46) = {305,306,343,-351,-350,-341};  //lower left
Plane Surface(47) = {46};
Line Loop(48) = {307,308,-342,-353,-352,-343};  //lower right
Plane Surface(49) = {48};

//Box 1
Line Loop(50) = {319,320,321,322,323,324,325,318};
Plane Surface(51) = {50};

//Box 2
Line Loop(52) = {349,350,351,352,353,354,355,348};
Plane Surface(53) = {52};

// Outer physical boundary
Physical Line(80) = {303,304,305};
Physical Line(81) = {306,307};
Physical Line(82) = {308,309,310};
Physical Line(83) = {302,301};

Physical Surface(1) = {45,41,43,47,49};  // background regions
Physical Surface(2) = {51};              // upper wg
Physical Surface(3) = {53};              // lower wg
