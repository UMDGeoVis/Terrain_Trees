/*
    This file is part of the Terrain Trees library.

    Author(s): Riccardo Fellegara (riccardo.fellegara@gmail.com)
               Federico Iuricich (federico.iuricich@gmail.com)

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

#ifndef FORMAN_GRADIENT_H
#define FORMAN_GRADIENT_H

#include <sys/time.h>
#include <sys/resource.h>

#include <list>
#include <sstream>
#include <queue>
#include <cstring>

#include "basic_types/mesh.h"
#include "utilities/timer.h"
//#include "terrain_trees/node_v.h"
//#include "terrain_trees/node_t.h"
//#include "terrain_trees/run_iterator.h"
//#include "terrain_trees/spatial_subdivision.h"

#include "ig.h"

#include "forman_gradient_aux_structure.h"
#include "forman_gradient_caching_structure.h"
#include "forman_gradient_stats.h"
using namespace forman_aux_structures;

typedef std::vector<Triangle_Arrows>			    CompressedToExpandedLUT;
typedef std::map<Triangle_Arrows, unsigned short>	ExpandedToCompressedLUT;

class Forman_Gradient
{
public:
    Forman_Gradient(itype num_t);
    //
    /// new uncompress gradient function
    inline TriGradient convert_compressed_to_expand(itype t_id)
    {
        if(forman_gradient[t_id-1] == OUTSIDECASES) return TriGradient();
        return TriGradient(compressed_to_expanded[forman_gradient[t_id-1]]);
    }
    inline ushort convert_expand_to_compressed(Triangle_Arrows ga) { return expanded_to_compressed[ga]; }
    bool is_valid_case(Triangle_Arrows arrow);

    inline bool is_triangle_critical(itype t) { return convert_compressed_to_expand(t).is_triangle_unpaired(); }

    inline bool is_edge_critical(const ivect &e, const pair<itype,itype> &et, Mesh &mesh)
    {
        Triangle& t = mesh.get_triangle(et.first);
        if(t.vertex_index(e[0])==-1)
        cout<<et.first<<endl;
        if(!convert_compressed_to_expand(et.first).is_edge_unpaired(t.vertex_index(e[0]),t.vertex_index(e[1])))
            return false;
        if(et.second!=-1)
        {
            Triangle& t2 = mesh.get_triangle(et.second);
            if(!convert_compressed_to_expand(et.second).is_edge_unpaired(t2.vertex_index(e[0]),t2.vertex_index(e[1])))
                return false;
        }
        return true;
    }

    inline bool is_edge_critical(const ivect &e, itype etstar, Mesh &mesh)
    {
        Triangle& t = mesh.get_triangle(etstar);
        return (convert_compressed_to_expand(etstar).is_edge_unpaired(t.vertex_index(e[0]),t.vertex_index(e[1])));
    }

    inline bool is_vertex_critical(itype v, itype t_id, Mesh &mesh)
    {
        //A TRIANGLE PAIRED WITH THE VERTEX (IF EXISTS) IS THE FIRST ENTRY OF VT(v)
        //THUS, IT IS SUFFICIENT TO CHECK IT TO UNDERSTAND IF v IS CRITICAL
        Triangle &t = mesh.get_triangle(t_id);
        TriGradient grad = convert_compressed_to_expand(t_id);
        return grad.is_vertex_unpaired(t.vertex_index(v));
    }

    void set_VE(itype v, itype v2, pair<itype,itype> &et, Mesh &mesh);
    void set_VE(itype v, itype v2, leaf_ET &et, Mesh &mesh);
    void free_VE(itype v1, itype v2, pair<itype,itype> &et, Mesh &mesh);

    void set_ET(itype t, const ivect &edge, Mesh &mesh);
    inline void set_ET(int v_pos, itype t)
    {
        TriGradient ga = convert_compressed_to_expand(t);
        ga.setEF(v_pos);
       // #pragma omp critical
        forman_gradient[t-1] = convert_expand_to_compressed(ga.getArrow()); //TESTING
    }

    inline void free_ET(int v_pos, itype t)
    {
        TriGradient ga = convert_compressed_to_expand(t);
        ga.clearEF(v_pos);
     //   #pragma omp critical
        forman_gradient[t-1] = convert_expand_to_compressed(ga.getArrow());//TESTING
    }

    inline void update_VE_adj_T(itype tid,itype v1, itype v2, Mesh &mesh, Forman_Gradient &gradient){

    TriGradient grad = gradient.convert_compressed_to_expand(tid);
    // cout<<"[DEBUG]tid:"<<tid<<"v1 and v2:"<<v1<<", "<<v2<<endl;
    grad.erase_edge_relation(mesh.get_triangle(tid).vertex_index(v1),mesh.get_triangle(tid).vertex_index(v2));
    grad.setVE(mesh.get_triangle(tid).vertex_index(v1),mesh.get_triangle(tid).vertex_index(v2));
    
    forman_gradient[tid-1]=convert_expand_to_compressed(grad.getArrow());
   
    //short v1i = gradient.convert_compressed_to_expand(tid).get_vertex_pair(mesh.get_triangle(tid).vertex_index(v1));
     //cout<<"Now vertex "<<v1<<" is paired with "<<mesh.get_triangle(tid).TV(v1i)<<endl;
    }
    
    /// for debug only
    vector<ushort>& get_gradient() { return forman_gradient; }

    inline void reorder_forman_gradient(Mesh &mesh, ivect& new_t_pos){
        // vector<ushort> new_forman_gradient;
        // int t=0;
        for(int i=0;i<new_t_pos.size();i++){
        if(new_t_pos[i]!=-1)
        {
            swap(forman_gradient[new_t_pos[i]-1], forman_gradient[i]);
            // new_forman_gradient.push_back(forman_gradient[t]);
        }
 
        }
        
        forman_gradient.resize(mesh.get_triangles_num());

        // cout<<"Old number of triangles:"<<forman_gradient.size()<<endl;
        // forman_gradient = new_forman_gradient;
        cout<<"Forman gradient updated"<<endl;
        // cout<<"number of triangles:"<<forman_gradient.size()<<endl;
    }

private: /// --- PRIVATE VARIABLES --- ///        
    vector<ushort> forman_gradient;
    CompressedToExpandedLUT compressed_to_expanded;
    ExpandedToCompressedLUT expanded_to_compressed;
};

#endif // FORMAN_GRADIENT_H
