// Template mesh geometry file for two suspended inclusions.
// Inclusions can be circular/elliptical (default), or square/rectangular.

// Force Gmsh to use legacy msh file format v2
Mesh.MshFileVersion = 2.2;

dx_in_nm = 1500.0;
dy_in_nm = 1500.0;
dy = dy_in_nm/dx_in_nm;

nslabs = 4;

// Unscaled vars
un_rib_w = 300.0;
un_rib_h = 200.0;
un_slab1_h = 100.0;
un_slab2_h = 100.0;
un_slab3_h = 100.0;
un_slab4_h = 100.0;
un_slab5_h = 100.0;

// Scaled vars
rib_w = un_rib_w/dx_in_nm;
rib_h = un_rib_h/dx_in_nm;
slab1_h = un_slab1_h/dx_in_nm;
slab2_h = un_slab2_h/dx_in_nm;
slab3_h = un_slab3_h/dx_in_nm;
slab4_h = un_slab4_h/dx_in_nm;
slab5_h = un_slab5_h/dx_in_nm;

x0=-.5;
x1=.5;
y0=dy/2;
y1=-dy/2;
xm=(x0+x1)/2;
ym=(y0+y1)/2;

slab1_y = ym - slab1_h;       // lower edge of slab 1
slab2_y = slab1_y - slab2_h;  // lower edge of slab 2
slab3_y = slab2_y - slab3_h;  // lower edge of slab 3
slab4_y = slab3_y - slab4_h;  // lower edge of slab 4
slab5_y = slab4_y - slab5_h;  // lower edge of slab 5


lc_bkg = 0.05;           // boundary
lc_refine_1 = lc_bkg/2.0; // on cylinder surfaces
lc_refine_2 = lc_bkg/2.0; // cylinder centres

// Boundary box

Point(1) = {x0, y0, 0, lc_bkg};    //NW
Point(2) = {x1, y0, 0, lc_bkg};    //NE
Point(3) = {x0, y1, 0, lc_bkg};    //SW
Point(4) = {x1, y1, 0, lc_bkg};    //SE

// Rib points
Point(5) = {xm-rib_w/2, ym, 0, lc_refine_1};
Point(6) = {xm+rib_w/2, ym, 0, lc_refine_1};
Point(7) = {xm-rib_w/2, ym+rib_h, 0, lc_refine_1};
Point(8) = {xm+rib_w/2, ym+rib_h, 0, lc_refine_1};

// Slab ends
Point(11) = {x0, ym, 0, lc_bkg};
Point(12) = {x1, ym, 0, lc_bkg};
Point(13) = {x0, slab1_y, 0, lc_bkg};
Point(14) = {x1, slab1_y, 0, lc_bkg};


// Top vacuum
Line(100) = {1,11};
Line(101) = {11,5};
Line(102) = {5,6};
Line(103) = {6,8};
Line(104) = {8,7};
Line(105) = {7,5};
Line(106) = {6,12};
Line(107) = {12,2};
Line(108) = {2,1};

Line Loop(40) = {100,101,-105,-104,-103,106,107,108};
Plane Surface(40) = {40};

//Rib
Line Loop(41) = {102,103,104,105};
Plane Surface(41) = {41};

//Bottom
Line(115) = {3,4};

//Slab 1
Line(110) = {11,13};
Line(111) = {13,14};
Line(112) = {14,12};

Line Loop(51) = {110,111,112,-106,-102,-101};
Plane Surface(51) = {51};


// Extra slabs

If (nslabs>=2)
    Point(15) = {x0, slab2_y, 0, lc_bkg};
    Point(16) = {x1, slab2_y, 0, lc_bkg};

    Line(120) = {13,15};
    Line(121) = {15,16};
    Line(122) = {16,14};
    Line Loop(52) = {120,121,122,-111};
    Plane Surface(52) = {52};


    If (nslabs>=3)
        Point(17) = {x0, slab3_y, 0, lc_bkg};
        Point(18) = {x1, slab3_y, 0, lc_bkg};
        Line(130) = {15,17};
        Line(131) = {17,18};
        Line(132) = {18,16};
        Line Loop(53) = {130,131,132,-121};
        Plane Surface(53) = {53};

        If (nslabs>=4)
            Point(19) = {x0, slab4_y, 0, lc_bkg};
            Point(20) = {x1, slab4_y, 0, lc_bkg};
            Line(140) = {17,19};
            Line(141) = {19,20};
            Line(142) = {20,18};
            Line Loop(54) = {140,141,142,-131};
            Plane Surface(54) = {54};

            If (nslabs>=5)
                Point(21) = {x0, slab5_y, 0, lc_bkg};
                Point(22) = {x1, slab5_y, 0, lc_bkg};
                Line(150) = {19,21};
                Line(151) = {21,22};
                Line(152) = {22,20};
                Line Loop(55) = {150,151,152,-141};
                Plane Surface(55) = {55};

            EndIf
        EndIf
    EndIf
EndIf

// Bottom background section

If (nslabs == 1)
    Line(200) = {13,3};
    Line(201) = {4,14};
    Line Loop(200) = {200,115,201,-111};
    Plane Surface(200) = {200};
EndIf

If (nslabs == 2)
    Line(200) = {15,3};
    Line(201) = {4,16};
    Line Loop(200) = {200,115,201,-121};
    Plane Surface(200) = {200};
EndIf

If (nslabs == 3)
    Line(200) = {17,3};
    Line(201) = {4,18};
    Line Loop(200) = {200,115,201,-131};
    Plane Surface(200) = {200};
EndIf

If (nslabs == 4)
    Line(200) = {19,3};
    Line(201) = {4,20};
    Line Loop(200) = {200,115,201,-141};
    Plane Surface(200) = {200};
EndIf

If (nslabs == 5)
    Line(200) = {21,3};
    Line(201) = {4,22};
    Line Loop(200) = {200,115,201,-151};
    Plane Surface(200) = {200};
EndIf


// Physical surfaces

Physical Surface(1) = {40,200};   // background regions
Physical Surface(2) = {41};   // rib
Physical Surface(3) = {51};   // slab 1

If (nslabs>=2)
    Physical Surface(4) = {52};   // slab 2
EndIf

If (nslabs>=3)
    Physical Surface(5) = {53};   // slab 3
EndIf

If (nslabs>=4)
    Physical Surface(6) = {54};   // slab 3
EndIf

If (nslabs>=5)
    Physical Surface(7) = {55};   // slab 3
EndIf


If (nslabs==1)
    Physical Line(80) = {100,110,200};
    Physical Line(81) = {115};
    Physical Line(82) = {201,112,107};
    Physical Line(83) = {108};
EndIf

If (nslabs==2)
    Physical Line(80) = {100,110,120,200};
    Physical Line(81) = {115};
    Physical Line(82) = {201,122,112,107};
    Physical Line(83) = {108};
EndIf

If (nslabs==3)
    Physical Line(80) = {100,110,120,130,200};
    Physical Line(81) = {115};
    Physical Line(82) = {201,132,122,112,107};
    Physical Line(83) = {108};
EndIf

If (nslabs==4)
    Physical Line(80) = {100,110,120,130,140,200};
    Physical Line(81) = {115};
    Physical Line(82) = {201,142,132,122,112,107};
    Physical Line(83) = {108};
EndIf


If (nslabs==5)
    Physical Line(80) = {100,110,120,130,140,150,200};
    Physical Line(81) = {115};
    Physical Line(82) = {201,152,142,132,122,112,107};
    Physical Line(83) = {108};
EndIf
