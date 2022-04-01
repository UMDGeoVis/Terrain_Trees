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
        delete[] c[i];
    }
    delete[] c;
    //

    return ret;
}
bool Geometry_Wrapper::same_side_of_edge(int v1,int v2, int e1, int e2, Mesh&mesh){

    Vertex va= mesh.get_vertex(e1);
    Vertex vb= mesh.get_vertex(e2);   
    int turn_v1=Geometry::PointTurn2D(mesh.get_vertex(v1).get_x(),mesh.get_vertex(v1).get_y(),va.get_x(),va.get_y(),vb.get_x(),vb.get_y());
    int turn_v2=Geometry::PointTurn2D(mesh.get_vertex(v2).get_x(),mesh.get_vertex(v2).get_y(),va.get_x(),va.get_y(),vb.get_x(),vb.get_y());
    if(turn_v1*turn_v2==1)
     return true;  // Do not return true if v1/v2 is on the edge. 
    else
     {
         return false;
     }

}

bool Geometry_Wrapper::point_in_circle(int t_id, int v_id, Mesh& mesh){

    Triangle tri = mesh.get_triangle(t_id); 
    Vertex v0=mesh.get_vertex(tri.TV(0));
    Vertex v1=mesh.get_vertex(tri.TV(1));
    Vertex v2=mesh.get_vertex(tri.TV(2));
    coord_type x_circle;
    coord_type y_circle;
    
   double A = v1.get_x() - v0.get_x();
   double B = v1.get_y() - v0.get_y();
   double C =v2.get_x()- v0.get_x();
   double D = v2.get_y() - v0.get_y();
   
   double E= A*(v0.get_x()+v1.get_x())+B*(v0.get_y()+v1.get_y());
   double F= C*(v0.get_x()+v2.get_x())+D*(v0.get_y()+v2.get_y());
   double G=2.0*((A*(v2.get_y()-v1.get_y()))-(B*(v2.get_x()-v1.get_x())));
   x_circle = ( (D*E) - (B*F) ) / G;
   y_circle = ( (A*F) - (C*E) ) / G;
    
   coord_type radius=sqrt((v0.get_x()-x_circle)*(v0.get_x()-x_circle)+(v0.get_y()-y_circle)*(v0.get_y()-y_circle));
   
    Vertex v = mesh.get_vertex(v_id);

   coord_type dist=sqrt((v.get_x()-x_circle)*(v.get_x()-x_circle)+(v.get_y()-y_circle)*(v.get_y()-y_circle));
   double EQ_TOLL=1E-8;
   if (dist>radius+EQ_TOLL)
       return false;
   else //if(dist<radius)
       return true;  //ignore the situation of four points on the same circle
//   else 
//       cout<<"four points are on the same circle"<<endl;
   

}


bool Geometry_Wrapper::point_in_circle_all(int t_id, Mesh& mesh){
    return point_in_circle_range(t_id, 1, mesh.get_vertices_num(),mesh);
}

bool Geometry_Wrapper::point_in_circle_range(int t_id, int v_start, int v_end, Mesh& mesh){

    Triangle tri = mesh.get_triangle(t_id); 
    Vertex v0=mesh.get_vertex(tri.TV(0));
    Vertex v1=mesh.get_vertex(tri.TV(1));
    Vertex v2=mesh.get_vertex(tri.TV(2));
    coord_type x_circle;
    coord_type y_circle;
    
   double A = v1.get_x() - v0.get_x();
   double B = v1.get_y() - v0.get_y();
   double C =v2.get_x()- v0.get_x();
   double D = v2.get_y() - v0.get_y();
   
   double E= A*(v0.get_x()+v1.get_x())+B*(v0.get_y()+v1.get_y());
   double F= C*(v0.get_x()+v2.get_x())+D*(v0.get_y()+v2.get_y());
   double G=2.0*((A*(v2.get_y()-v1.get_y()))-(B*(v2.get_x()-v1.get_x())));
   x_circle = ( (D*E) - (B*F) ) / G;
   y_circle = ( (A*F) - (C*E) ) / G;
    
   coord_type radius=sqrt((v0.get_x()-x_circle)*(v0.get_x()-x_circle)+(v0.get_y()-y_circle)*(v0.get_y()-y_circle));
   
    for(int v_id = v_start; v_id < v_end; v_id++){
        Vertex v = mesh.get_vertex(v_id);
        coord_type dist=sqrt((v.get_x()-x_circle)*(v.get_x()-x_circle)+(v.get_y()-y_circle)*(v.get_y()-y_circle));
        double EQ_TOLL=1E-8;
        if (dist+EQ_TOLL<radius)
            return true;

    }
    return false;


}


bool Geometry_Wrapper::point_in_circle_range(Point& center, coord_type radius, int v_start, int v_end, Mesh& mesh){

 
    coord_type x_circle = center.get_x();
    coord_type y_circle = center.get_y();
    

    for(int v_id = v_start; v_id < v_end; v_id++){
        Vertex v = mesh.get_vertex(v_id);
        coord_type dist=sqrt((v.get_x()-x_circle)*(v.get_x()-x_circle)+(v.get_y()-y_circle)*(v.get_y()-y_circle));
        double EQ_TOLL=1E-8;
        if (dist+EQ_TOLL<radius)
            return true;

    }
    return false;


}



// bool Geometry_Wrapper::box_circumcircle_disjoint(int t_id, Box& box, Mesh& mesh){

//     Triangle tri = mesh.get_triangle(t_id); 
//     Vertex v0=mesh.get_vertex(tri.TV(0));
//     Vertex v1=mesh.get_vertex(tri.TV(1));
//     Vertex v2=mesh.get_vertex(tri.TV(2));
//     coord_type x_circle;
//     coord_type y_circle;
    
//    double A = v1.get_x() - v0.get_x();
//    double B = v1.get_y() - v0.get_y();
//    double C =v2.get_x()- v0.get_x();
//    double D = v2.get_y() - v0.get_y();
   
//    double E= A*(v0.get_x()+v1.get_x())+B*(v0.get_y()+v1.get_y());
//    double F= C*(v0.get_x()+v2.get_x())+D*(v0.get_y()+v2.get_y());
//    double G=2.0*((A*(v2.get_y()-v1.get_y()))-(B*(v2.get_x()-v1.get_x())));
//    x_circle = ( (D*E) - (B*F) ) / G;
//    y_circle = ( (A*F) - (C*E) ) / G;
    
//    coord_type radius=sqrt((v0.get_x()-x_circle)*(v0.get_x()-x_circle)+(v0.get_y()-y_circle)*(v0.get_y()-y_circle));
   
//     Point max(x_circle+radius, y_circle+radius);

//     Point min(x_circle-radius, y_circle-radius);
//     Box bounding_box(min,max);
//     return !box.intersects(bounding_box);

// }

double Geometry_Wrapper::triangle_compactness(int t_id, Mesh& mesh){

    Triangle& t = mesh.get_triangle(t_id);

    coord_type area = 0;

    Vertex &v1 = mesh.get_vertex(t.TV(0));
    Vertex &v2 = mesh.get_vertex(t.TV(1));
    Vertex &v3 = mesh.get_vertex(t.TV(2));

    dvect e1 = {v1.get_x()-v2.get_x(),v1.get_y()-v2.get_y(),v1.get_z()-v2.get_z()};
    dvect e2 = {v3.get_x()-v2.get_x(),v3.get_y()-v2.get_y(),v3.get_z()-v2.get_z()};

    area = (e1[1]*e2[2]-e1[2]*e2[1])*(e1[1]*e2[2]-e1[2]*e2[1]);
    area += (e1[0]*e2[2]-e1[2]*e2[0])*(e1[0]*e2[2]-e1[2]*e2[0]);
    area += (e1[0]*e2[1]-e1[1]*e2[0])*(e1[0]*e2[1]-e1[1]*e2[0]);
    // |x1 y1 z1|
    // |x2 y2 z2|
    // |x3 y3 z3|

    // area += v1.get_x()*(v2.get_y()*v3.get_z()-v3.get_y()*v2.get_z());
    // area -= v2.get_x()*(v1.get_y()*v3.get_z()-v3.get_y()*v1.get_z());
    // area += v3.get_x()*(v1.get_y()*v2.get_z()-v2.get_y()*v1.get_z());    
    // //  area = abs(Det3D(v1.get_x(),v1.get_y(),v1.get_z(),v2.get_x(),v2.get_y(),v2.get_z(),v3.get_x(),v3.get_y(),v3.get_z()));
    
    area = sqrt(area) / 2;
    coord_type edgeSum = 0;

    for(int vid =0; vid<3;vid++){
    
        Vertex va = mesh.get_vertex(t.TV((vid+1)%3));
        Vertex vb = mesh.get_vertex(t.TV((vid+2)%3));
        coord_type xdist = va.get_x() - vb.get_x();
        coord_type ydist = va.get_y() - vb.get_y();
        coord_type zdist = va.get_z() - vb.get_z();
        edgeSum += xdist*xdist + ydist*ydist + zdist*zdist;
       // cout<<edgeSum<<endl;
    }

    return 4*sqrt(3)*area/edgeSum;



}