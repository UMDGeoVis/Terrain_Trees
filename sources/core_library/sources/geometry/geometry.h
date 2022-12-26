/*
    This file is part of the Terrain Trees library.

    Author(s): Riccardo Fellegara (riccardo.fellegara@gmail.com)

    This project has been supported by the Italian Ministry of Education and
    Research under the PRIN 2009 program, and by the National Science Foundation
    under grant number IIS-1116747.

    The Terrain Trees library is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    The Terrain Trees library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with the Terrain Trees library.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _GEOMETRY_H
#define	_GEOMETRY_H

#include <cmath>
#include <iostream>
#include <stdio.h>
using namespace std;

// ------------ defines and macros for tetra_in_box --------------------------
#ifndef PI
#define PI 3.14159265358979323846
#endif
/* ------------------------------------------------------------------------ */
/*                    Tolerance used in computations                        */
/* ------------------------------------------------------------------------ */
#define ZERO (10E-14)
#define Coincide(a,b) (fabs((a)-(b))<=ZERO)
/* ------------------------------------------------------------------------ */
/*                     Basic arithmetic functions                           */
/* ------------------------------------------------------------------------ */
#define Square(x) ((x)*(x))
/*
Calculate determinant  |a b|
                       |c d|
*/
#define Det2D(a,b,c,d)  ( (a*d)-(b*c) )
/*
Calculate determinant  |a1 a2 a3|
                       |b1 b2 b3|
                       |c1 c2 c3|
*/
#define Det3D(a1,a2,a3,b1,b2,b3,c1,c2,c3)	\
    ( a1*Det2D(b2,b3,c2,c3) - a2*Det2D(b1,b3,c1,c3) + a3*Det2D(b1,b2,c1,c2) )

/*
 *  Turns
 */
#define LEFT_TURN -1
#define UP_TURN -1
#define NO_TURN 0
#define RIGHT_TURN 1
#define DOWN_TURN 1

class Geometry
{
protected:
    // ------------ Main functions -------------
    // this version of triangle_in_box_strict execute the following tests:
    // - if all vertices are on the same side of the box, no intersection
    // - consider all the faces of the box open (if a vertex is adjacent to one of these face is considered external)
    // - check if one of the vertex is inside the triangle (PointInTriangle_strict)
    // - check if the box center is inside the triangle (PointInTriangle_strict)
    // - check if an edge is at least partially inside the box, then the triangle intersects the box (ClipLine2D_strict)
    // - check if one triangle edge is alinged with, and contains one box edge, then if the triangle is extended on the interior side of the box w.r.t. to that box edge, then there is intersection (OverlapXSegment)
    //This version is suited for box query, while for building we have to do external tests if we use different assumption on the closeness of box faces      
    static int triangle_in_box_strict(double * minF, double * maxF, double **c);


    // --------------- Auxiliaries for tetra_in_box -----------------
    static int PointInTriangle2D (double x, double y,
                           double x1, double y1,
                           double x2, double y2,
                           double x3, double y3);    

    /*
    Here the triangle is given as two arrays x[3], y[3] containing the coordinates of its three vertices
    */
    static int PointInTriangle (double xp, double yp, double * v1, double * v2, double * v3);
    static int PointInTriangle_old (double xp, double yp, double * v1, double * v2, double * v3);
    static int PointInTriangle_strict (double xp, double yp, double * v1, double * v2, double * v3);


    /*
    Calculate sign of determinant  |a b|
                                   |c d|
    */
    static int DetSign2D (double a, double b, double c, double d);
    /*
    Calculate sign of determinant |a1 a2 a3|
                                  |b1 b2 b3|
                                  |c1 c2 c3|
    */
    static int DetSign3D(double a1, double a2, double a3,
                  double b1, double b2, double b3,
                  double c1, double c2, double c3);

    /* ------------------------------------------------------------------------ */
    /*                                 Turns                                    */
    /* ------------------------------------------------------------------------ */

    static int PointTurn2D(double x,double y,double x1,double y1,double x2,double y2)	{ return Geometry::DetSign2D((x)-(x1), (y)-(y1), (x2)-(x1), (y2)-(y1)); }

    /* ------------------------------------------------------------------------ */
    /*                       Intersection test w.r.t. a box (2D)                */
    /* ------------------------------------------------------------------------ */
    /*
    Restrict the admissible interval [u1,u2] by intersecting it with the
    half-line, solution of inequality u*p <= q.
    Return 1 if the resulting interval is not empty, 0 otherwise.
     */
    static int ClipTest2D_strict(double p, double q, double * u1, double * u2);
    /*
    Return 1 if the segment is at least partially inside the box.
     */
    static int ClipLine2D_strict(double minX, double minY, double maxX, double maxY, /* box */
                          double x1, double y1, double x2, double y2); /* line */
    /*
    Return 1 iff edge (x1,y1) (x2,y2) overlaps edge having x=x0 and y01<=y<y02
     */
    static int OverlapXSegment(double x1, double y1, double x2, double y2,
                        double x0, double y01, double y02);
    /*
    Return 1 if the triangle is at least partially inside the box.
     */
    static int ClipTriangle2D_strict(double minX, double minY, double maxX, double maxY, /* box */
                              double x[3], double y[3]); /* triangle */
};



#endif	/* _GEOMETRY_H */

