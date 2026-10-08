side = 0.5;
n = 11;  // 11 points per edge -> 10x10 quads, 121 vertices (matches current mesh)

Point(1) = {0, 0, 0};
Point(2) = {side, 0, 0};
Point(3) = {side, side, 0};
Point(4) = {0, side, 0};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Transfinite Line{1,2,3,4} = n;
Transfinite Surface{1};
Recombine Surface{1};

Mesh 2;
