
wall = 0.12;

L = 1.0;


/* 
 *        4                4L                  3
 *         o -------------------------------- o         
 *         |                                  |       
 *         |                                  |       
 *        1|      5                           |       Y         
 *      L  ------                             | L     ^
 *         |     |                            |       |
 *         |     |  6                         |       |
 *         o -------------------------------- o       o -----> X
 *                          4L                 2
 * */

Point(1)  = {0.0, 0.5*L, 0.0,  wall}; // p0
Point(2)  = {8*L, 0.0, 0.0,  wall}; // p3
Point(3)  = {8*L, 1*L, 0.0,  wall}; // p4
Point(4)  = {0.0, 1*L, 0.0,  wall}; // p5
Point(5)  = {1.0*L, 0.5*L, 0.0,  wall}; // p5
Point(6)  = {1.0*L, 0.0, 0.0,  wall}; // p5

Line(1) = {1, 5};
Line(2) = {5, 6};
Line(3) = {6, 2};
Line(4) = {2, 3};
Line(5) = {3, 4};
Line(6) = {4, 1};

//+
Line Loop(1) = {1, 2, 3, 4, 5, 6};
Plane Surface(1) = {1};

//+ boundary conditions for stream function
Physical Line('inlet') = {6};
Physical Line('outlet') = {4};
Physical Line('paredeInf') = {1, 2, 3};
Physical Line('paredeSup') = {5};

Physical Surface('surface') = {1};
