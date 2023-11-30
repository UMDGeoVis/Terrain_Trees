#include "geometry_slope.h"

// coord_type Geometry_Slope::compute_triangle_slope(Triangle &t, Mesh &mesh)
// {
//     Vertex &v1 = mesh.get_vertex(t.TV(0));
//     Vertex &v2 = mesh.get_vertex(t.TV(1));
//     Vertex &v3 = mesh.get_vertex(t.TV(2));

//     // get the two vectors in the triangle
//     dvect u = { v2.get_x()-v1.get_x() , v2.get_y()-v1.get_y() , v2.get_z()-v1.get_z() };
//     dvect v = { v3.get_x()-v1.get_x() , v3.get_y()-v1.get_y() , v3.get_z()-v1.get_z() };

//     // get the cross product to get the normal
//     dvect n = { u[1]*v[2] - u[2]*v[1] , u[2]*v[0] - u[0]*v[2], u[0]*v[1] - u[1]*v[0] };
//     // get the magnitude of the normal
//     coord_type magnitude_n = sqrt( n[0]*n[0] + n[1]*n[1] + n[2]*n[2] );

//     // the dot product between n and the ground normal is only formed by the normalized y component of n
//     n[1] /= magnitude_n;
//     // thus, we do not have to compute the real dot product but we use the y part of the norm
//     coord_type angle = acos(n[1])*180/3.14159265358979323846;
//     coord_type slope = tan(angle);
//     return angle;
// }

//// Updated version of slope computation
coord_type Geometry_Slope::compute_triangle_slope(Triangle &t, Mesh &mesh)
{
    Vertex &v1 = mesh.get_vertex(t.TV(0));
    Vertex &v2 = mesh.get_vertex(t.TV(1));
    Vertex &v3 = mesh.get_vertex(t.TV(2));

    // get the two vectors in the triangle
    dvect u = { v2.get_x()-v1.get_x() , v2.get_y()-v1.get_y() , v2.get_z()-v1.get_z() };
    dvect v = { v3.get_x()-v1.get_x() , v3.get_y()-v1.get_y() , v3.get_z()-v1.get_z() };

    // get the cross product to get the normal
    dvect n = { u[1]*v[2] - u[2]*v[1] , u[2]*v[0] - u[0]*v[2], u[0]*v[1] - u[1]*v[0] };
    // get the magnitude of the normal
    //coord_type slope = sqrt( n[0]*n[0] + n[1]*n[1])/n[2];
    coord_type magnitude_n = sqrt( n[0]*n[0] + n[1]*n[1] + n[2]*n[2] );

    // the dot product between n and the ground normal is only formed by the normalized y component of n
    n[2] /= magnitude_n;
    // thus, we do not have to compute the real dot product but we use the y part of the norm
    coord_type angle = acos(abs(n[2]))*180/3.14159265358979323846;

   // coord_type slope2= tan(angle);
    //if (abs(slope-slope2)>0.001)
      //  cerr<<"difference is  "<<abs(slope-slope2)<<endl;
    
            

    return angle;
}

coord_type Geometry_Slope::compute_edge_slope(ivect &e, Mesh &mesh)
{
    Vertex &v1 = mesh.get_vertex(e[0]);
    Vertex &v2 = mesh.get_vertex(e[1]);
    Vertex v3 = Vertex(v2.get_x(), v2.get_y(), v1.get_z()); // temp vertex to get a ground

    // get the two vectors in the triangle
    dvect u = { v2.get_x() - v1.get_x() , v2.get_y() - v1.get_y() , v2.get_z() - v1.get_z() };
    dvect v = { v3.get_x() - v1.get_x() , v3.get_y() - v1.get_y() , v3.get_z() - v1.get_z() };

    coord_type u_mag = sqrt( u[0]*u[0] + u[1]*u[1] + u[2]*u[2] );
    dvect u_norm = { u[0]/u_mag , u[1]/u_mag , u[2]/u_mag };

    coord_type v_mag = sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] );
    dvect v_norm = { v[0]/v_mag , v[1]/v_mag , v[2]/v_mag };

    coord_type dot_prod = u_norm[0]*v_norm[0] + u_norm[1]*v_norm[1] + u_norm[2]*v_norm[2];

    coord_type angle = acos(dot_prod);

    return angle;
}
