/*
    This file is part of the Terrain Trees library.

    Author(s): Riccardo Fellegara (riccardo.fellegara@gmail.com)

    This project has been supported by the Italian Ministry of Education and
    Research under the PRIN 2009 program, and by the National Science Foundation
    under grant number IIS-1116747.

    The triangle Trees library is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    The triangle Trees library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with the triangle Trees library.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef GEOMETRY_WRAPPER_H
#define GEOMETRY_WRAPPER_H

#include "basic_types/mesh.h"
#include "geometry.h"
/**
 * @brief The Geometry_Wrapper class provides an interface for executing geometric tests for generating trees and answering queries
 */
class Geometry_Wrapper : public Geometry
{
public:
    /**
     * @brief A public static method that computes the centroid of a triangle
     *
     * @param t_id an integer representing the position index of the triangle
     * @param p a Point& that will contains the centroid coordinates
     * @param mesh a Mesh&, the triangle mesh
     */
    static void get_triangle_centroid(int t_id, Point& p, Mesh &mesh);
    /**
     * @brief A public static method that computes the point-in-triangle geometric test
     *
     * @param t_id an integer representing the position index of the triangle
     * @param p a Point& representing the point to test
     * @param mesh a Mesh&, the triangle mesh
     * @return true if the point is contained in the triangle, false otherwise
     */
    static bool point_in_triangle(int t_id, Point& point, Mesh &mesh);
    /**
     * @brief A public static method that computes the triangle-in-box geometric test
     * NOTA: this procedure is used during the generation process of a tree.
     * It considers a triangle internal even if it has only a vertex inside the box
     *
     * @param t_id an integer representing the position index of the triangle
     * @param box a Box& representing the domain of a block
     * @param mesh a Mesh&, the triangle mesh
     * @return true if exists an intersection between t_id and box, false otherwise
     */
    static bool triangle_in_box_build(int t_id, Box& box, Mesh& mesh);
    /**
     * @brief A public static method that computes the triangle-in-box geometric test
     * NOTA: this procedure is used in box queries.
     * It considers all the faces of the box open as we want to output the tetrahedra
     * with a real intersection with the box and not only a vertex adjacent to one of box faces.
     *
     * @param t_id an integer representing the position index of the triangle
     * @param box a Box& representing the domain of a block
     * @param mesh a Mesh&, the triangle mesh
     * @return true if exists a real intersection between t_id and box, false otherwise
     */
    static bool triangle_in_box(int t_id, Box& box, Mesh& mesh);
    static bool same_side_of_edge(int v1,int v2, int e1, int e2, Mesh& mesh);
    static bool point_in_circle(int t_id, int v_id, Mesh& mesh);
    static bool point_in_circle_range(int t_id, int v_start, int v_end, Mesh& mesh);
    static bool point_in_circle_all(int t_id,  Mesh& mesh);

    static bool point_in_circle_range(Point& center, coord_type radius, int v_start, int v_end, Mesh& mesh);
    static double triangle_compactness(int t_id, Mesh& mesh);

};

#endif // GEOMETRY_WRAPPER_H
