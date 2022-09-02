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

#ifndef FORMAN_GRADIENT_FEATURES_EXTRACTOR_H
#define FORMAN_GRADIENT_FEATURES_EXTRACTOR_H

#include "forman_gradient.h"
#include "dangling_paths_handler.h"
#include "forman_gradient_topological_relations.h"
#include "forman_gradient_computation.h"
#include "boundary_matrix.h"

class Forman_Gradient_Features_Extractor
{
public:
    Forman_Gradient_Features_Extractor()
    {
        stats = FGV_cache_dangling_paths_stats();
        new_paths = 0;
        u_paths = 0;
        cache_t = 0;
        rels_gathering = 0;
        triangles_2celle = simplices_map();
        edges_1celle = simplices_map();

        found_max = /*found_2selle =*/ found_1selle = found_min = 0;
    }

    void extract_descending_2cells(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, Spatial_Subdivision &division, Node_V &root,
                                   OpType operation, int cache_size);
    void extract_descending_1cells(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, Spatial_Subdivision &division, Node_V &root,
                                   OpType operation, int cache_size);

    void extract_ascending_2cells(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, Spatial_Subdivision &division, Node_V &root,
                                  OpType operation, int cache_size);
    void extract_ascending_1cells(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, Spatial_Subdivision &division, Node_V &root,
                                  OpType operation, int cache_size);

    void extract_incidence_graph(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, Spatial_Subdivision &division,  OpType operation, int cache_size);
    void extract_incidence_graph(Node_T &n, Box &n_dom, Mesh &mesh, Forman_Gradient &gradient, Spatial_Subdivision &division,  OpType operation, int cache_size);

    void extract_critical_clusters(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, map<short,ivect_set> &critical_simplices,
                                   Spatial_Subdivision &division,  OpType operation, int cache_size, string mesh_name);
    void extract_critical_clusters(Node_T &n, Mesh &mesh, Forman_Gradient &gradient, map<short, ivect_set> &critical_simplices,
                                   Spatial_Subdivision &division, OpType operation, int cache_size, string mesh_name);

    // for saving to vtk file
    inline void init_segmentation_vector(Mesh &mesh) { segmentation  = ivect(mesh.get_triangles_num(), -1); }
    inline ivect& get_segmentation_vector() { return segmentation; }
    
    inline void reset_output_structures(Mesh &mesh)
    {
        segmentation.clear();
        segmentation  = ivect(mesh.get_triangles_num(), -1);

        triangles_2celle.clear();
        edges_1celle.clear();
        manifold_2celle_asc.clear();

        forman_ig.clear();
    }
    
    //for extracting ascending manifolds
    inline void init_ascending_segmentation_vector(Mesh &mesh) { manifold_2celle_asc  = ivect(mesh.get_vertices_num(), -1); }
    inline ivect& get_ascending_segmentation() { return manifold_2celle_asc; }

    inline simplices_map& get_extracted_cells(TopType cell_type)
    {
        if(cell_type == TRIANGLE)
            return triangles_2celle;
        else if(cell_type == EDGE)
            return edges_1celle;
        else
        {
            cerr<<"[get_extracted_cells] not a valid cell type"<<endl;
            return edges_1celle;
        }
    }

    inline IG& get_incidence_graph() { return forman_ig; }

    //reset di tutte le variabili di classe
    inline void reset_timer_variables()
    {
        new_paths = 0;
        u_paths = 0;
        cache_t = 0;
        rels_gathering = 0;
    }

    inline void print_feature_extraction_time()
    {
        cerr<<rels_gathering<<" "<<new_paths<<" "<<u_paths<<" "<<cache_t<<endl;
    }

    inline void print_stats() { stats.print_stats(); }
    inline void reset_stats() { stats.reset_stats_counter(); }
    inline void set_filtration_vec(uvect input){this->filtration=input;}

    inline void reset_extraction_critica_counters()
    {
        found_max = found_1selle = found_min = 0;
        sad_min=sad_max=level2=0;
    }
    inline void print_arc_nums(){cout<<"saddle and minimum: "<<sad_min<<endl;
    cout<<"saddle and maximum: "<<sad_max<<endl;
    cout<<"Level 2 arcs:"<<level2<<endl;}
    
    void getBoundarySimplices(ivect simpl_id, int dim, vector<vector<int>>& boundary);

protected:    
    /// variables neede for output purposes
    ivect segmentation;
    simplices_map triangles_2celle;
    simplices_map edges_1celle;
    ivect manifold_2celle_asc;
    IG forman_ig;
    uvect filtration; //for each vertex its filtration value
    
    //contains the net timing obtained for a single manifold extraction
    coord_type new_paths, u_paths, cache_t, rels_gathering;

    //various stats from feature extraction
    FGV_cache_dangling_paths_stats stats;

    // for debug
    int found_max, found_1selle, found_min;

    int sad_min,sad_max,level2;
//    // for debug only
//    map<ivect,ET> saddles_around_maxima;
    
private:
    /// ----- DESCENDING MANIFOLD EXTRACTION FUNCTIONS ----- ///
    //functions used to compute the descending_2cells (implemented in forman_gradient_descending2cells.h)
    void descending_2cells_extraction(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, Spatial_Subdivision &division,
                                      Node_V &root, OpType operation,
                                      leaves_2_desc_map &dangling_paths, et_cache &cache);
    void descending_2cells_extraction_leaf(Node_V &n, Mesh &mesh, Forman_Gradient &gradient,
                                           Node_V &root, Spatial_Subdivision &division, OpType operation,
                                           leaves_2_desc_map &dangling_paths, et_cache &cache);
    void get_new_descending_2cells(Node_V &n, Mesh &mesh, Forman_Gradient &gradient,
                                   leaf_ET &local_ef, Node_V &root, Spatial_Subdivision &division,
                                   OpType operation, leaves_2_desc_map &dangling_paths, et_cache &cache);
    void get_dangling_descending_2cells(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, leaf_ET &local_ef,
                                        Node_V &root, Spatial_Subdivision &division, OpType operation,
                                        leaves_2_desc_map &dangling_paths, et_cache &cache);
    void get_one_descending_2cells(Node_V &n, itype t_id, Mesh &mesh,  Forman_Gradient &gradient, itype label, leaf_ET &local_ef,
                                   Node_V &root, Spatial_Subdivision &division, OpType operation,
                                   leaves_2_desc_map &dangling_paths, et_cache &cache);

    //functions used to compute the descending_1cells (implemented in forman_gradient_descending1cells.h)
    void descending_1cells_extraction(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, Spatial_Subdivision &division,  Node_V &root,
                                      OpType operation, leaves_1_desc_map &dangling_paths,
                                      vtstar_cache &cache);
    void descending_1cells_extraction_leaf(Node_V &n, Mesh &mesh, Forman_Gradient &gradient,
                                           Node_V &root, Spatial_Subdivision &division, OpType operation,
                                           leaves_1_desc_map &dangling_paths, vtstar_cache &cache);
    void get_new_descending_1cells(Node_V &n, Mesh &mesh, Forman_Gradient &gradient,  desc1rels &all_rels,
                                   Node_V &root, Spatial_Subdivision &division,
                                   OpType operation, leaves_1_desc_map &dangling_paths, vtstar_cache &cache);
    void get_dangling_descending_1cells(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, desc1rels &all_rels,
                                        Node_V &root, Spatial_Subdivision &division, OpType operation, leaves_1_desc_map &dangling_paths,
                                        vtstar_cache &cache);
    void get_one_descending_1cells(Node_V &n, Edge* edge, Mesh &mesh, Forman_Gradient &gradient, itype label, desc1rels &all_rels,
                                   Node_V &root, Spatial_Subdivision &division, OpType operation,
                                   leaves_1_desc_map &dangling_paths, vtstar_cache &cache);

    /// ----- ASCENDING MANIFOLD EXTRACTION FUNCTIONS ----- ///
    void ascending_2cells_extraction(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, Spatial_Subdivision &division,  Node_V &root,
                                     OpType operation, leaves_2_asc_map &dangling_paths, asc3_cache &cache);
    void ascending_2cells_extraction_leaf(Node_V &n, Mesh &mesh, Forman_Gradient &gradient,
                                          Node_V &root, Spatial_Subdivision &division, OpType operation,
                                          leaves_2_asc_map &dangling_paths, asc3_cache &cache);
    void get_new_ascending_2cells(Node_V &n, Mesh &mesh, Forman_Gradient &gradient,  asc2rels &all_rels, Node_V &root, Spatial_Subdivision &division,
                                  OpType operation, leaves_2_asc_map &dangling_paths, asc3_cache &cache);
    void get_dangling_ascending_2cells(Node_V &n, Mesh &mesh, Forman_Gradient &gradient,  asc2rels &all_rels,
                                       Node_V &root, Spatial_Subdivision &division, OpType operation,
                                       leaves_2_asc_map &dangling_paths, asc3_cache &cache);
    void get_one_ascending_2cells(Node_V &n, itype crit_v, Mesh &mesh, Forman_Gradient &gradient, itype label, asc2rels &all_rels,
                                  Node_V &root, Spatial_Subdivision &division, OpType operation, leaves_2_asc_map &dangling_paths, asc3_cache &cache);

    void ascending_1cells_extraction(Node_V &n, Mesh &mesh,  Forman_Gradient &gradient, Spatial_Subdivision &division, Node_V &root,
                                     OpType operation, leaves_asc_map &dangling_paths, et_cache &cache);
    void ascending_1cells_extraction_leaf(Node_V &n, Mesh &mesh, Forman_Gradient &gradient,  Node_V &root, Spatial_Subdivision &division, OpType operation,
                                          leaves_asc_map &dangling_paths, et_cache &cache);
    void get_new_ascending_1cells(Node_V &n, Mesh &mesh, Forman_Gradient &gradient,  leaf_ET &local_ef,
                                  Node_V &root, Spatial_Subdivision &division, OpType operation,
                                  leaves_asc_map &dangling_paths, et_cache &cache);
    void get_dangling_ascending_1cells(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, leaf_ET &local_ef,
                                       Node_V &root, Spatial_Subdivision &division, OpType operation,
                                       leaves_asc_map &dangling_paths, et_cache &cache);
    void get_one_ascending_1cells(Node_V &n, const pair<itype,short> &tf, Mesh &mesh, Forman_Gradient &gradient,
                                  itype label, leaf_ET &local_ef, Node_V &root, Spatial_Subdivision &division, OpType operation,
                                  leaves_asc_map &dangling_paths, et_cache &cache);

    /// ----- CRITICAL CLUSTERS EXTRACTION FUNCTIONS ----- ///
    void critical_clusters_extraction(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, ivect_set &saddles, Spatial_Subdivision &division,
                                      Node_V &root, vtstar_cache &cache, critical_clusters &cc);
    void critical_clusters_extraction(Node_T &n, Box &n_dom, int n_level, Mesh &mesh, Forman_Gradient &gradient, ivect_set &saddles,
                                      Spatial_Subdivision &division, Node_T &root, vtstar_cache &cache, critical_clusters &cc);
    void critical_clusters_extraction_leaf(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, ivect_set &saddles, Spatial_Subdivision &division,
                                           Node_V &root, vtstar_cache &cache, critical_clusters &cc);
    void critical_clusters_extraction_leaf(Node_T &n, Box &n_dom, Mesh &mesh, Forman_Gradient &gradient, ivect_set &saddles,
                                           Spatial_Subdivision &division, Node_T &root, vtstar_cache &cache, critical_clusters &cc);
    void init_critical_t_labels(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, ivect_set &maxima, Spatial_Subdivision &division,
                                critical_clusters &cc);
    void init_critical_t_labels(Node_T &n, Box &n_dom, int level, Mesh &mesh, Forman_Gradient &gradient, ivect_set &maxima,
                                Spatial_Subdivision &division, critical_clusters &cc);
    void clusterize_minima(const ivect &saddle, leaf_VTstar &vtstars, vtstar_cache &cache, Node_V &n, Mesh &mesh, Forman_Gradient &gradient,
                           Spatial_Subdivision &division, Node_V &root, critical_clusters &cc);
    void clusterize_minima(const ivect &saddle, leaf_VTstar &vtstars, vtstar_cache &cache, Node_T &n, itype v_start, itype v_end, Mesh &mesh,
                           Forman_Gradient &gradient, Spatial_Subdivision &division, Node_T &root, critical_clusters &cc);
    void clusterize_maxima(const ivect &saddle, leaf_ET &ets, Forman_Gradient &gradient, critical_clusters &cc);

protected:
    /// ----- MORSE INCIDENCE GRAPH FUNCTIONS ----- ///
    void incidence_graph_extraction(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, Spatial_Subdivision &division,  Node_V &root,
                                    OpType operation, mig_cache &cache, ig_paths &paths);
    void incidence_graph_extraction_leaf(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, Node_V &root, Spatial_Subdivision &division,
                                         OpType operation, mig_cache &cache, ig_paths &paths);
    void get_new_IG_paths(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, IG &ig, local_VTstar_ET &leaf_struct, Node_V &root,
                          Spatial_Subdivision &division, OpType operation, mig_cache &cache, ig_paths &paths);
    void get_dangling_IG_paths(Node_V &n, Mesh &mesh, Forman_Gradient &gradient, IG &ig, local_VTstar_ET &leaf_struct, Node_V &root,
                               Spatial_Subdivision &division, OpType operation, mig_cache &cache, ig_paths &paths);
    void explore_desc1cell_mig(Node_V &n, Edge *edge, itype last_v, short label, Mesh &mesh, Forman_Gradient &gradient, IG &ig,
                               local_VTstar_ET &local_rels, iNode *saddle_node, Node_V &root, Spatial_Subdivision &division,
                               OpType operation, mig_cache &cache, ig_paths &paths);
    void explore_asc1cell_mig(Node_V &n, const pair<itype, short> &tf, iNode *saddle_node, itype first_t_id, leaf_ET &local_ef, Mesh &mesh,
                              Forman_Gradient &gradient, IG &ig, Node_V &root, Spatial_Subdivision &division, OpType operation,
                              leaves_1_asc_mig_map &paths, et_cache &cache);
    /// return true if the edge has been found and thus the vertex is not a minimum, false otherwise
    static bool get_next_edge(Node_V &n, itype v, Edge *&e, leaf_VTstar &all_vtstar, vtstar_cache &cache,
                              Node_V &root, Spatial_Subdivision &division, Mesh &mesh, Forman_Gradient &gradient);
    /// PMR-T tree procedures for extracting the MIG
    void incidence_graph_extraction(Node_T &n, Box &n_dom, int n_level, Mesh &mesh, Forman_Gradient &gradient, Spatial_Subdivision &division,  Node_T &root,
                                    OpType operation, mig_cache &cache, ig_paths &paths);
    void incidence_graph_extraction_leaf(Node_T &n, Box &n_dom, Mesh &mesh, Forman_Gradient &gradient, Node_T &root, Spatial_Subdivision &division,
                                         OpType operation, mig_cache &cache, ig_paths &paths);
    void get_new_IG_paths(Node_T &n, itype v_start, itype v_end, Mesh &mesh, Forman_Gradient &gradient, IG &ig, local_VTstar_ET &leaf_struct, Node_T &root,
                          Spatial_Subdivision &division, OpType operation, mig_cache &cache, ig_paths &paths);
    void get_dangling_IG_paths(Node_T &n, itype v_start, itype v_end, Mesh &mesh, Forman_Gradient &gradient, IG &ig, local_VTstar_ET &leaf_struct, Node_T &root,
                               Spatial_Subdivision &division, OpType operation, mig_cache &cache, ig_paths &paths);
    void explore_desc1cell_mig(Node_T &n, itype v_start, itype v_end, Edge *edge, itype last_v, short label, Mesh &mesh, Forman_Gradient &gradient, IG &ig,
                               local_VTstar_ET &local_rels, iNode *saddle_node, Node_T &root, Spatial_Subdivision &division,
                               OpType operation, mig_cache &cache, ig_paths &paths);
    void explore_asc1cell_mig(Node_T &n, itype v_start, itype v_end, const pair<itype, short> &tf, iNode *saddle_node, itype first_t_id, leaf_ET &local_ef, Mesh &mesh,
                              Forman_Gradient &gradient, IG &ig, Node_T &root, Spatial_Subdivision &division, OpType operation,
                              leaves_1_asc_mig_map &paths, et_cache &cache);
    /// return true if the edge has been found and thus the vertex is not a minimum, false otherwise
    static bool get_next_edge(Node_T &n, itype v_start, itype v_end, itype v, Edge *&e, leaf_VTstar &all_vtstar, vtstar_cache &cache, Node_T &root,
                              Spatial_Subdivision &division, Mesh &mesh, Forman_Gradient &gradient);

    inline void pair_minimum_saddle(iNode* saddle_node, itype label, itype v, itype last_v, IG &ig)
    {
        nNode* minimum = ig.find_minimum(v);
        if(minimum == NULL)
        {
            minimum = new nNode(v);
                            sad_min++;
            ig.add_minimum(v,minimum);
            ig.addArc(minimum,last_v,saddle_node,label,0);
        }
        else
        {
            Arc *arc = ig.already_connected(minimum,saddle_node);
            if(arc==NULL)
            {
                sad_min++;
                ig.addArc(minimum,last_v,saddle_node,label,0);
            }
            else
               { arc->setLabel(2); ///this edge cannot be simplified
                 level2++;
               }
        }
    }

    inline void pair_saddle_maximum(iNode* saddle_node, itype first_t_id, itype t, itype last_t, IG &ig)
    {
        nNode* maximum = ig.find_maximum(t);
        if(maximum == NULL)
        {
            maximum = new nNode(t);
            ig.add_maximum(t,maximum);
                            sad_max++;
            ig.addArc(saddle_node,first_t_id,maximum,last_t,1);
        }
        else
        {
            Arc *arc = ig.already_connected(maximum,saddle_node);
            if(arc==NULL)
            {
                sad_max++;
                ig.addArc(saddle_node,first_t_id,maximum,last_t,1);
            }
            else
                {arc->setLabel(2); ///this edge cannot be simplified
                level2++;
                }
        }
    }
    inline itype get_max_elevation_vertex(const ivect &vect)
    {
        if(vect.size() == 0)
        {
            cerr << "[ERROR] get_max_field_vertex -> empty vector" << endl;
            return -1;
        }

        itype max = vect[0];
        coord_type max_field = filtration[vect[0]-1];
        for(unsigned v=1; v<vect.size(); v++)
        {
            if(filtration[vect[v]-1] > max_field)
            {
                max = vect[v];
                max_field = filtration[vect[v]-1];
            }
        }
        return max;
    }
    
    inline itype get_max_elevation_vertex(Triangle &t)
    {
        itype max = t.TV(0);
        coord_type max_field = filtration[t.TV(0)-1];
        for(itype v=1; v<t.vertices_num(); v++)
        {
            if(filtration[t.TV(v)-1] > max_field)
            {
                max = t.TV(v);
                max_field = filtration[t.TV(v)-1];
            }
        }
        return max;
    }
    
};



#endif // FORMAN_GRADIENT_FEATURES_EXTRACTOR_H
