#include "Roughness_circle.h"


void Roughness_circle::compute(Node_V& root, Mesh& mesh, const coord_type radius, Spatial_Subdivision &division)
{
    this->radius = radius;
    int v_num = mesh.get_vertices_num();
    // v_roughness.assign(v_num, 0);
    #pragma omp paralllel for
    for(unsigned i = 1; i <= v_num; i++)
    {
        vector<coord_type> elevations;
        roughness_vertex(mesh.get_vertex(i), elevations, root, mesh, mesh.get_domain(), 0, division);
        // first field is the roughness, second field is the count of vertices in that range
        mesh.get_vertex(i).add_field(calc_roughness(elevations));
        mesh.get_vertex(i).add_field(elevations.size());
    
    }

    this->print_roughness_stats(mesh, mesh.get_vertex(1).get_fields_num() - 1);


}

void Roughness_circle::roughness_vertex(Point& center, vector<coord_type>& elev, Node_V& n, Mesh& mesh, Box& n_dom, int level, Spatial_Subdivision& division)
{
    Point max(center.get_x() + radius, center.get_y() + radius);
    Point min(center.get_x() - radius, center.get_y() - radius);
    Box bounding_box(min, max);
    if(!n_dom.intersects(bounding_box))
        return;
    coord_type x_circle = center.get_x();
    coord_type y_circle = center.get_y();
    coord_type radius_sq = radius*radius;
    if(n.is_leaf())
    {
        if(!n.indexes_vertices()) return;
        itype v_start = n.get_v_start();
        itype v_end = n.get_v_end();
        itype v_range = v_end - v_start;
        for(unsigned i=0;i<v_range;i++)
        {
            itype real_v_id=v_start+i;
            Vertex v = mesh.get_vertex(real_v_id);
            coord_type dist_sq = (v.get_x()-x_circle)*(v.get_x()-x_circle)+(v.get_y()-y_circle)*(v.get_y()-y_circle);
            double EQ_TOLL = 1E-10;
            if(dist_sq + EQ_TOLL <= radius_sq)
            {
                elev.push_back(v.get_z());
            }
        }
    }
    else
    {
        for (int i = 0; i < division.son_number(); i++)
        {
            Box son_dom = division.compute_domain(n_dom,level,i);
            int son_level = level +1;
            if (n.get_son(i) != NULL)
            {
                roughness_vertex(center, elev, *n.get_son(i),  mesh, son_dom, son_level, division);
            }
        }        
    }
}

coord_type Roughness_circle::calc_roughness(vector<coord_type>& elevations)
{
    coord_type zDistSum = 0;
    coord_type zSum=0;

    for(auto z:elevations)
    { 
        zSum += z;
        
    }
    coord_type zAVG=zSum/elevations.size();

    for(auto z:elevations)
    {
       
        zDistSum+= (z - zAVG)*(z - zAVG);
                
    }
    
    coord_type roughness=sqrt(zDistSum/(elevations.size()));
    return roughness;
}

