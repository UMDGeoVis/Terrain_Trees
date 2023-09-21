#ifndef ROUGHNESS_CIRCLE_H
#define ROUGHNESS_CIRCLE_H

#pragma once

#include "terrain_trees/node_v.h"
#include "terrain_trees/node_t.h"
#include "terrain_trees/spatial_subdivision.h"
#include <limits.h>
#include <float.h>

using namespace std;

class Roughness_circle
{
public:
    Roughness_circle() {radius = 0;};
    ~Roughness_circle() {};
    void compute(Node_V& root, Mesh &mesh, const coord_type radius, Spatial_Subdivision &division);

    inline void print_roughness_stats(Mesh &mesh, int f_pos){
        coord_type r_min=INFINITY, r_max=-INFINITY, r_sum=0;
        for(itype v=1; v<=mesh.get_vertices_num(); v++)
        {
            coord_type roughness = mesh.get_vertex(v).get_field(f_pos);
            r_max = max(r_max, roughness);
            r_min = min(r_min, roughness);
            r_sum += roughness;
        }
        cerr<<"[STATS] roughness min: "<< r_min<<" avg: "<<r_sum/(coord_type)mesh.get_vertices_num()<<" max: "<< r_max<<endl;

    }

private:
    // void roughness_vertex_leaf(Point& center, Node_V &n, Mesh &mesh,  PRT_Tree &tree);

    void roughness_vertex(Point& center, vector<coord_type>& elev, Node_V &n, Mesh &mesh,  Box &n_dom, int level, Spatial_Subdivision &division);
    coord_type calc_roughness(vector<coord_type>& elev);

    // vector<coord_type> v_roughness;
    coord_type radius;
};

#endif