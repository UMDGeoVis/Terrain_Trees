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

#ifndef FORMAN_GRADIENT_COMPUTATION_H
#define FORMAN_GRADIENT_COMPUTATION_H

#include <boost/bind.hpp>
#include "forman_gradient.h"
#include "terrain_trees/node_v.h"
#include "terrain_trees/node_t.h"
#include "terrain_trees/spatial_subdivision.h"

class Forman_Gradient_Computation
{
public:
    Forman_Gradient_Computation() {}

    void compute_gradient_vector(Forman_Gradient &gradient, Node_V &n, Mesh &mesh, Spatial_Subdivision &division);
    void compute_gradient_vector(Forman_Gradient &gradient, Node_T &n, Box &n_dom, Mesh &mesh, Spatial_Subdivision &division);
    

    
    void initial_filtering(Mesh &mesh);
    void initial_filtering_IA(Mesh &mesh);
    inline uvect get_filtration(){return this->filtration;}
    
    // update the filtration after building the tree since now the vertices are reindexed. 
    void reset_filtering(Mesh &mesh,ivect& original_vertex_indices, bool new_to_old);
    inline map<short,ivect_set>& get_critical_simplices() { return this->critical_simplices; }
    //TODO: add a function to count all critical simplices in the current mesh. 
    //To be used after simplification. 
    inline void count_critical_simplices(Forman_Gradient &gradient,Mesh &mesh){

        

    };
protected:
    uvect filtration; //for each vertex its filtration value
    map<short,ivect_set> critical_simplices;
    itype extract_missing_id(ivect &smaller, ivect &bigger);

    bool cmp_filtered_simplices(const ivect& lhs, const ivect& rhs);
    utype simplex_filtration(const ivect &simpl);
        
    void push_in_coboundary(simplex *sface, simplex *co_sface);
    

    void add_VT_entry(itype v_id, itype t_id, Node_V &n, leaf_VT &vts);
    void add_ET_entry(ivect &e, itype max_field_id, itype t_id, Node_V &n, leaf_ET &ets);

private:
    int saddle, _max, _min;
    inline bool filtration_cmp(const pair<coord_type,itype>& v1, const pair<coord_type,itype>& v2) const
    {
        if(v1.first == v2.first){
            return v2.second < v1.second;
        }

        return v1.first < v2.first;
    }
    inline bool sort_filtered_vertices(const itype& v1,const itype& v2) { return filtration[v1-1] > filtration[v2-1]; }

    vector<uvect> componentBasedFiltration; //injective function for each component [number of fields x Vertices ]

    coord_type star_time;
    coord_type gradient_time;
    void compute_local_gradient_vector(Forman_Gradient &gradient, Node_V &n, Mesh &mesh, Spatial_Subdivision &division);
    void compute_local_gradient_vector(Forman_Gradient &gradient, Node_T &n, Box &n_dom, int level, Mesh &mesh, Spatial_Subdivision &division);
    //functions used to compute the gradient vector
    void compute_local_gradient_vector_leaf(Forman_Gradient &fg, Node_V &n, Mesh& mesh);
    void compute_local_gradient_vector_leaf(Forman_Gradient &fg, Node_T &n, Box &n_dom, Mesh& mesh);
    //compute the gradient in a leaf block
    void compute_gradient(Forman_Gradient &fg, itype v_start, leaf_VT &vts, leaf_ET &ets,  vector< lower_star_map > &lower_stars, Mesh &mesh);
    //compute gradient vector field for a given vertex
    void compute_gradient_in_a_vertex(ivect &v, lower_star_map &lower_star, VT &vts,
                                                              simplices_map &critical_points, gradient &local_vectors, Mesh &mesh);
   // void compute_gradient_in_a_vertex_SF(lower_star_map &lower_star_group, Forman_Gradient &fg, VT &vt, leaf_ET &ets, Mesh &mesh);

    

    void extract_topological_relations(Node_V &n, Mesh &mesh, vector<lower_star_map> &lower_stars, leaf_VT &vts, leaf_ET &ets);
    void extract_topological_relation(Node_V &n, itype t_id, Triangle &t, vector<lower_star_map> &lower_stars, leaf_VT &vts, leaf_ET &ets,Mesh &mesh);
    void extract_topological_relations(Node_T &n, itype v_start, itype v_end, Mesh &mesh, vector<lower_star_map> &lower_stars,
                                       leaf_VT &vts, leaf_ET &ets);
    void extract_topological_relation(Node_T &n, itype v_start, itype v_end, itype t_id, Triangle &t, vector<lower_star_map> &lower_stars,
                                      leaf_VT &vts, leaf_ET &ets,Mesh &mesh);
    
//    bool cmp_filtered_simplices_SF(const double& lhs, const double& rhs);
 
   

    
    

    void init_lower_stars_groups(itype v_id, vector<lower_star_map> &lw_groups, lower_star_map &lower_star);
   
    void add_lower_star_entry(itype v_id, Triangle &t,  lower_star_map &lower_stars, Mesh &mesh);
  
    int numPairableLowerStar(const ivect &next, const lower_star_map &sset, ivect& pair);

    //int numPairableLowerStar_SF(const ivect &next, const lower_star_map_SF &sset, ivect& pair);
   // void setPair_SF(ivect &next, ivect &pair, Forman_Gradient &fg, VT &vt, leaf_ET &ets, Mesh &mesh);

    void initialize_forman_gradient(Forman_Gradient &fg, gradient &grad,leaf_ET &ets, Mesh &mesh);

    //void initialize_forman_gradient_SF(gradient &grad,leaf_ET &ets, Mesh &mesh);
//    void get_lower_stars_VT_and_EF_SF(Node_V &n,lower_star_map &lower_star,leaf_VT &vts, leaf_ET &ets);
//    void get_single_lower_stars_VT_and_EF_SF(itype t_id,Triangle &t,Node_V &n,lower_star_map &lower_star,leaf_VT &vts, leaf_ET &ets);

// Functions for single field

  
    inline void push_ordered(list<vector<int> >& coda, vector<int>& simplesso){

        bool same=false;

        list<vector<int> >::iterator it = coda.begin();
        for(; it != coda.end(); it++)
        {
            if(!is_lower(*it, simplesso, &same)) break;
        }

        if(!same)
        {
            if(it == coda.begin()) coda.push_front(simplesso);
            else if(it == coda.end()) coda.push_back(simplesso);
            else coda.insert(it, simplesso);
        }
    }
    
    
    int num_unpared_faces(vector<int> &co_face, gradient &local_vectors, simplices_map &critical_points, Mesh &mesh);
    bool is_lower(const vector<int> &simplex1, const vector<int> &simplex2, bool* same);
    vector<int> unique_pairable_face(vector<int> &s_face, gradient &local_vectors, simplices_map &critical_points, Mesh &mesh);
    int get_t_id(const vector<int> &popped, vector<int> &vt, Mesh &mesh);
    void sort_simplex(vector<int> &simplex, Mesh &mesh);

    
};


#endif // FORMAN_GRADIENT_COMPUTATION_H
