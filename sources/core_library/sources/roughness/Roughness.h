/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Roughness.h
 * Author: ytsong
 *
 * Created on January 6, 2019, 2:33 PM
 */

#ifndef ROUGHNESS_H
#define ROUGHNESS_H
#include "terrain_trees/node_v.h"
#include "terrain_trees/node_t.h"
#include "terrain_trees/spatial_subdivision.h"
using namespace std;

class Roughness 
{
public:
    Roughness() {}
    Roughness(Mesh &mesh){}
    virtual ~Roughness(){}
  void compute(Node_V &n, Mesh &mesh, Spatial_Subdivision &division);

  void compute(Node_T &n, Box &domain, int level, Mesh &mesh, Spatial_Subdivision &division);
  inline void print_roughness_stats(Mesh &mesh, int c_pos)
    {
        coord_type min=INFINITY, max=-INFINITY, avg=0;
        for(itype v=1; v<=mesh.get_vertices_num(); v++)
        {
            coord_type c = mesh.get_vertex(v).get_field(c_pos);

            if(c < min)
                min = c;
            if(c > max)
                max = c;
            avg += c;
        }
        cerr<<"[STATS] roughness min: "<<min<<" avg: "<<avg/(coord_type)mesh.get_vertices_num()<<" max: "<<max<<endl;
    }
  
// void store_result(Mesh &mesh){

//     for(itype v=1;v<=mesh.get_vertices_num();v++)
//     {
//         mesh.get_vertex(v).add_field(roughness[v]);
       
//     }

// };
//  void print_time(){cout<<"Overall time in leaves"<<time_in_leaves<<endl;};
private:
 
 //   dvect roughness;
   
   leaf_VV_vec  vvs;
   void roughness_leaf(Node_V& n, Mesh& mesh);
   void roughness_leaf(Node_T &n, Box &n_dom, Mesh &mesh);

    /**
     * @brief roughness_leaf
     * @param n
     * @param mesh
     */
 //   void roughness_leaf(Node_T& n, Box &n_dom, Mesh& mesh)=0;

        //void finalize_local_roughness(itype v_start, leaf_VT &vts, boost::dynamic_bitset<> &is_v_border, dvect &local_curvatures, Mesh &mesh);
};

#endif /* ROUGHNESS_H */

