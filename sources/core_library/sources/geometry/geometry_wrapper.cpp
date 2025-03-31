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

#include "geometry_wrapper.h"
#include <boost/dynamic_bitset.hpp>

void Geometry_Wrapper::get_triangle_centroid(int t_id, Point& p, Mesh &mesh)
{
    Triangle &tri = mesh.get_triangle(t_id);

    Vertex& v0 = mesh.get_vertex(tri.TV(0));
    Vertex& v1 = mesh.get_vertex(tri.TV(1));
    Vertex& v2 = mesh.get_vertex(tri.TV(2));

    for(int i=0; i<v0.get_dimension(); i++)
        p.set_c(i,(v0.get_c(i) + v1.get_c(i) + v2.get_c(i)) / 3.0);
}

bool Geometry_Wrapper::point_in_triangle(int t_id, Point& point, Mesh &mesh)
{
    Triangle &tri = mesh.get_triangle(t_id);
    coord_type **c;
    c = new coord_type* [3];
    for (int i = 0; i < 3; i++) {
        c[i] = new coord_type[2];
        Vertex &v = mesh.get_vertex(tri.TV(i));
        c[i][0] = v.get_x();
        c[i][1] = v.get_y();
    }

    bool ret = false;

    if(PointInTriangle(point.get_x(),point.get_y(),c[0],c[1],c[2]))
        ret = true;

    for (int i = 0; i < 3; i++)
        delete c[i];
    delete c;

    return ret;
}

coord_type Geometry_Wrapper::elevation_interpolate_tri(int t_id, Point& point, Mesh &mesh)
{
    Triangle &tri = mesh.get_triangle(t_id);
    Vertex& v1 = mesh.get_vertex(tri.TV(0));
    Vertex& v2 = mesh.get_vertex(tri.TV(1));
    Vertex& v3 = mesh.get_vertex(tri.TV(2));
    if(v1 == v2 || v1 == v3 || v2 == v3){
        cerr << "Two vertices of a triangle overlapped, this is an error." <<endl;
        return 0;
    }
    coord_type w1 = ((v2.get_y() - v3.get_y()) * (point.get_x() - v3.get_x()) + (v3.get_x() - v2.get_x()) * (point.get_y() - v3.get_y()))/
                    ((v2.get_y() - v3.get_y()) * (v1.get_x() - v3.get_x()) + (v3.get_x() - v2.get_x()) * (v1.get_y() - v3.get_y()));
    coord_type w2 = ((v3.get_y() - v1.get_y()) * (point.get_x() - v3.get_x()) + (v1.get_x() - v3.get_x()) * (point.get_y() - v3.get_y()))/
                    ((v2.get_y() - v3.get_y()) * (v1.get_x() - v3.get_x()) + (v3.get_x() - v2.get_x()) * (v1.get_y() - v3.get_y()));
    coord_type w3 = 1 - w1 - w2;
    return w1 * v1.get_z() + w2 * v2.get_z() + w3 * v3.get_z();

}

bool Geometry_Wrapper::triangle_in_box_build(int t_id, Box& box, Mesh& mesh)
{
    Triangle &t = mesh.get_triangle(t_id);

    //I consider internal a triangle that at least has a vertex inside the box (also only a vertex)
    for(int v=0; v<t.vertices_num(); v++)
    {
        if(box.contains(mesh.get_vertex(t.TV(v)),mesh.get_domain().get_max()))
            return true;
    }

    return Geometry_Wrapper::triangle_in_box(t_id,box,mesh);
}

bool Geometry_Wrapper::triangle_in_box(int t_id, Box& box, Mesh& mesh)
{
    Triangle &t = mesh.get_triangle(t_id);

    coord_type* minf = new coord_type[2];
    minf[0] = box.get_min().get_x(); minf[1] = box.get_min().get_y();
    coord_type* maxf = new coord_type[2];
    maxf[0] = box.get_max().get_x(); maxf[1] = box.get_max().get_y();

    coord_type **c;
    c = new coord_type* [t.vertices_num()];
    for (int i = 0; i < t.vertices_num(); i++)
    {        
        Vertex& v = mesh.get_vertex(t.TV(i));
        c[i] = new coord_type[v.get_dimension()];

        c[i][0] = v.get_x();
        c[i][1] = v.get_y();
    }

    bool ret = false;

    if(triangle_in_box_strict(minf, maxf, c))
        ret = true;

    // free memory
    delete[] minf;
    delete[] maxf;
    for (int i = 0; i < t.vertices_num(); i++) {
        delete c[i];
    }
    delete[] c;
    //

    return ret;
}
