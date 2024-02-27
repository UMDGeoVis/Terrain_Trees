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

#include "geometry.h"

/******************************************************************************/

int Geometry::DetSign2D(double a, double b, double c, double d)
{
    double t1, t2;
    t1 = (a * d);
    t2 = (b * c);
    if (t1 > (t2 + ZERO)) return 1;
    if (t2 > (t1 + ZERO)) return -1;
    return 0;
}

/*
Calculate sign of determinant |a1 a2 a3|
                              |b1 b2 b3|
                              |c1 c2 c3|
 */
int Geometry::DetSign3D(double a1, double a2, double a3,
              double b1, double b2, double b3,
              double c1, double c2, double c3)
{
    double d = Det3D(a1, a2, a3, b1, b2, b3, c1, c2, c3);
    if (fabs(d) <= ZERO) return 0;
    return ( (d > 0.0) ? 1 : -1);
}


/******************************************/

/******************************************************************************/

int Geometry::PointInTriangle2D(double x, double y,
                      double x1, double y1,
                      double x2, double y2,
                      double x3, double y3)
{
    if ((Geometry::PointTurn2D(x, y, x1, y1, x2, y2) == LEFT_TURN) &&
            (Geometry::PointTurn2D(x, y, x2, y2, x3, y3) == LEFT_TURN) &&
            (Geometry::PointTurn2D(x, y, x3, y3, x1, y1) == LEFT_TURN))
        return 1;
    if ((PointTurn2D(x, y, x1, y1, x2, y2) == RIGHT_TURN) &&
            (PointTurn2D(x, y, x2, y2, x3, y3) == RIGHT_TURN) &&
            (PointTurn2D(x, y, x3, y3, x1, y1) == RIGHT_TURN))
        return 1;
    return 0;
}

/* ------------------------------------------------------------------------ */
/*                       Intersection test w.r.t. a box                     */
/* ------------------------------------------------------------------------ */

/*
WITH A TWO-DIMENSIONAL BOX
Implementation based on the line clipping algorithm by Liang-Barsky.
The basic idea is the following:
Consider the parametric equation of the segment to be clipped:

x(u) = x1 + u*(x2-x1)             with  0<=u<=1
y(u) = y1 + u*(y2-y1)
The segment intersects the box if there exists u, 0<=u<=1, such that
minX <= x(u) <= maxX   and   minY <= y(u) <= maxY
or, equivalently
minX-x1 <= u*(x2-x1) <= maxX-x1   and   minY-y1 <= u*(y2-y1) <= maxY-y1

We have four inequalities that give four conditions on u.
Check if the four conditions are mutually consistent and
consistent with condition 0<=u<=1.
 */

/*
Restrict the admissible interval [u1,u2] by intersecting it with the
half-line, solution of inequality u*p <= q.
Return 1 if the resulting interval is not empty, 0 otherwise.
 */
int Geometry::ClipTest2D_strict(double p, double q, double * u1, double * u2)
{
    double r;

    if (p < 0.0) {
        r = q / p;
        if (r >= (*u2)) return 0;
        else if (r > (*u1)) (*u1) = r;
    } else {
        if (p > 0.0) {
            r = q / p;
            if (r <= (*u1)) return 0;
            else if (r < (*u2)) (*u2) = r;
        } else {
            /* p==0.0 line parallel to clipping edge */
            if (q <= 0.0) return 0;
        }
    }
    return 1;
}

/*
Return 1 if the segment is at least partially inside the box.
 */
int Geometry::ClipLine2D_strict(double minX, double minY, double maxX, double maxY, /* box */
                      double x1, double y1, double x2, double y2) /* line */
{
    double u1 = 0.0, u2 = 1.0; /* admissible interval, initially all [0,1] */
    double dx = x2 - x1, dy = y2 - y1;
    if ( ClipTest2D_strict(-dx, x1 - minX, &u1, &u2) &&
         ClipTest2D_strict(dx, maxX - x1, &u1, &u2) &&
         ClipTest2D_strict(-dy, y1 - minY, &u1, &u2) &&
         ClipTest2D_strict(dy, maxY - y1, &u1, &u2) )
    {
        return 1;
    }
    return 0;
}

/*
Return 1 iff edge (x1,y1) (x2,y2) overlaps edge having x=x0 and y01<=y<y02
 */
int Geometry::OverlapXSegment(double x1, double y1, double x2, double y2,
                    double x0, double y01, double y02)
{
    if ((x1 != x0) || (x2 != x0)) return 0; /* edge does not lie on x=x0 */
    if ((y1 <= y01) && (y2 <= y01)) return 0; /* edge <= y01 */
    if ((y1 >= y02) && (y2 >= y02)) return 0; /* edge >= y02 */
    return 1;
}

/*
Return 1 if the triangle is at least partially inside the box.
 */
int Geometry::ClipTriangle2D_strict(double minX, double minY, double maxX, double maxY, /* box */
                          double x[3], double y[3]) /* triangle */
{
    int i;
    /* if all vertices are on the same side of the box, no intersection */
    //  cout<<"controllo che tutti i vertici siano da uno stesso lato del triangolo"<<endl;
    if (x[0] <= minX && x[1] <= minX && x[2] <= minX) return 0;
    if (x[0] >= maxX && x[1] >= maxX && x[2] >= maxX) return 0;
    if (y[0] <= minY && y[1] <= minY && y[2] <= minY) return 0;
    if (y[0] >= maxY && y[1] >= maxY && y[2] >= maxY) return 0;


    /* if a vertex is inside the box, then the triangle intersects the box */
    for (i = 0; i < 3; i++) {
        if ((x[i] < maxX) && (x[i] > minX) && (y[i] < maxY) && (y[i] > minY)) {
            return 1;
        }
    }
    /* if an edge is at least partially inside the box, then the triangle
     intersects the box */
    for (i = 0; i < 3; i++) {
        if (ClipLine2D_strict(minX, minY, maxX, maxY,
                              x[i], y[i], x[(i + 1) % 3], y[(i + 1) % 3])) {
            return 1;
        }
    }

    /* none of the three triangle edges intersects the box */
    /* check if the triangle completely contains the box
     by applying the point-in-triangle test on the center of the box */
    if (PointInTriangle2D(0.5 * (minX + maxX), 0.5 * (minY + maxY),
                          x[0], y[0], x[1], y[1], x[2], y[2])) {
        return 1;
    }

    /* if one triangle edge is alinged with, and contains one box edge,
     then if the triangle is extended on the interior side of the box
     w.r.t. to that box edge, then there is intersection */
    for (i = 0; i < 3; i++)
    {
        if (OverlapXSegment(x[i], y[i], x[(i + 1) % 3], y[(i + 1) % 3], minX, minY, maxY)
                && (x[(i + 2) % 3] > minX)) /* third triangle vertex on the right */
            return 1;
        if (OverlapXSegment(x[i], y[i], x[(i + 1) % 3], y[(i + 1) % 3], maxX, minY, maxY)
                && (x[(i + 2) % 3] < maxX)) /* third triangle vertex on the left */
            return 1;
        if (OverlapXSegment(y[i], x[i], y[(i + 1) % 3], x[(i + 1) % 3], minY, minX, maxX)
                && (y[(i + 2) % 3] > minY)) /* third triangle vertex below */
            return 1;
        if (OverlapXSegment(y[i], x[i], y[(i + 1) % 3], x[(i + 1) % 3], maxY, minX, maxX)
                && (y[(i + 2) % 3] < maxY)) /* third triangle vertex below */
            return 1;
    }
    return 0;
}

/******************************************************************************/

/**
  **
  **
  **
  **/
int Geometry::PointInTriangle(double xp, double yp, double * v1, double * v2, double * v3)
{
    // we check if the point is a vertex of the triangle
    if (xp == v1[0] && yp == v1[1]) return 1;
    if (xp == v2[0] && yp == v2[1]) return 1;
    if (xp == v3[0] && yp == v3[1]) return 1;
    v2[0] -= v1[0]; v3[0] -= v1[0];
    v2[1] -= v1[1]; v3[1] -= v1[1];
    xp -= v1[0]; yp -= v1[1];
    v1[0] = 0; v1[1] = 0;
    //  static double * v[4] = {v1, v2, v3, v4};
    int d1, d2, d3, orientation;
    orientation = DetSign3D(v1[0], v1[1], 1,
                            v2[0], v2[1], 1,
                            v3[0], v3[1], 1);

    d1 = DetSign3D(xp, yp, 1,
                   v2[0], v2[1], 1,
                   v3[0], v3[1], 1);

    if (d1 != orientation && d1 != 0) return 0;

    d2 = DetSign3D(v1[0], v1[1], 1,
                   xp, yp, 1,
                   v3[0], v3[1], 1);

    if (d2 != orientation && d2 != 0) return 0;

    d3 = DetSign3D(v1[0], v1[1], 1,
                   v2[0], v2[1], 1,
                   xp, yp, 1);

    if (d3 != orientation && d3 != 0) return 0;

    return 1;
}

int Geometry::PointInTriangle_old(double xp, double yp, double * v1, double * v2, double * v3)
{
    // we check if the point is a vertex of the triangle
    if (xp == v1[0] && yp == v1[1]) return 1;
    if (xp == v2[0] && yp == v2[1]) return 1;
    if (xp == v3[0] && yp == v3[1]) return 1;

    //  static double * v[4] = {v1, v2, v3, v4};
    int d1, d2, d3, orientation;
    orientation = DetSign3D(v1[0], v1[1], 1,
                            v2[0], v2[1], 1,
                            v3[0], v3[1], 1);

    d1 = DetSign3D(xp, yp, 1,
                   v2[0], v2[1], 1,
                   v3[0], v3[1], 1);

    if (d1 != orientation && d1 != 0) return 0;

    d2 = DetSign3D(v1[0], v1[1], 1,
                   xp, yp, 1,
                   v3[0], v3[1], 1);

    if (d2 != orientation && d2 != 0) return 0;

    d3 = DetSign3D(v1[0], v1[1], 1,
                   v2[0], v2[1], 1,
                   xp, yp, 1);

    if (d3 != orientation && d3 != 0) return 0;

    return 1;
}

int Geometry::PointInTriangle_strict(double xp, double yp, double * v1, double * v2, double * v3)
{
    double d1, d2, d3, d;
    d = DetSign3D(v1[0], v1[1], 1,
                  v2[0], v2[1], 1,
                  v3[0], v3[1], 1);

    d1 = DetSign3D(xp, yp, 1,
                   v2[0], v2[1], 1,
                   v3[0], v3[1], 1);

    if (d1 != d) return 0;

    d2 = DetSign3D(v1[0], v1[1], 1,
                   xp, yp, 1,
                   v3[0], v3[1], 1);

    if (d2 != d) return 0;

    d3 = DetSign3D(v1[0], v1[1], 1,
                   v2[0], v2[1], 1,
                   xp, yp, 1);

    if (d3 != d) return 0;
    return 1;
}

/******************************************************************************/

/******************************************************************************/

int Geometry::triangle_in_box_strict(double * minF, double * maxF, double **c)
{
    /* if all vertices are on the same side of the box, no intersection */
    for (int j = 0; j < 2; j++)
    {
        if ((c[0][j] <= minF[j]) && (c[1][j] <= minF[j]) && (c[2][j] <= minF[j]) )
        {
            return 0;
        }
        if ((c[0][j] >= maxF[j]) && (c[1][j] >= maxF[j]) && (c[2][j] >= maxF[j]))
        {
            return 0;
        }
    }

    int i;

    /* check if one triangle vertex is inside box */
    for (i = 0; i < 3; i++)
    {
        if ((minF[0] < c[i][0]) && (c[i][0] < maxF[0]) &&
                (minF[1] < c[i][1]) && (c[i][1] < maxF[1]))
        {
            return 1;
        }
    }

    /* check if one box vertex is inside the triangle */
    if (PointInTriangle_strict(minF[0], minF[1], c[0], c[1], c[2]) || PointInTriangle_strict(minF[0], maxF[1], c[0], c[1], c[2]) ||
            PointInTriangle_strict(maxF[0], minF[1], c[0], c[1], c[2]) || PointInTriangle_strict(maxF[0], maxF[1], c[0], c[1], c[2]))
    {
        return 1;
    }

    /* check if box center is inside triangle */
    if (PointInTriangle_strict(0.5 * (minF[0] + maxF[0]),
                            0.5 * (minF[1] + maxF[1]), c[0], c[1], c[2]))
    {
        return 1;
    }

    /* if an edge is at least partially inside the box, then the triangle
     intersects the box */
    for (i = 0; i < 3; i++) {
        if (ClipLine2D_strict(minF[0], minF[1], maxF[0], maxF[1],
                              c[i][0], c[i][1], c[(i + 1) % 3][0], c[(i + 1) % 3][1])) {
            return 1;
        }
    }

    /* if one triangle edge is alinged with, and contains one box edge,
     then if the triangle is extended on the interior side of the box
     w.r.t. to that box edge, then there is intersection */
    for (i = 0; i < 3; i++)
    {
        if (OverlapXSegment(c[i][0], c[i][1], c[(i + 1) % 3][0], c[(i + 1) % 3][1], minF[0], minF[1], maxF[1])
                && (c[(i + 2) % 3][0] > minF[0])) /* third triangle vertex on the right */
            return 1;
        if (OverlapXSegment(c[i][0], c[i][1], c[(i + 1) % 3][0], c[(i + 1) % 3][1], maxF[0], minF[1], maxF[1])
                && (c[(i + 2) % 3][0] < maxF[0])) /* third triangle vertex on the left */
            return 1;
        if (OverlapXSegment(c[i][1], c[i][0], c[(i + 1) % 3][1], c[(i + 1) % 3][0], minF[1], minF[0], maxF[0])
                && (c[(i + 2) % 3][1] > minF[1])) /* third triangle vertex below */
            return 1;
        if (OverlapXSegment(c[i][1], c[i][0], c[(i + 1) % 3][1], c[(i + 1) % 3][0], maxF[1], minF[0], maxF[0])
                && (c[(i + 2) % 3][1] < maxF[1])) /* third triangle vertex below */
            return 1;
    }
    /* box and triangle do not intersect each other */
    return 0;
}
