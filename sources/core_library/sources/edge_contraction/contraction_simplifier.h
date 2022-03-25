#ifndef CONTRACTION_SIMPLIFIER_H
#define CONTRACTION_SIMPLIFIER_H

#include "terrain_trees/node_v.h"
#include "terrain_trees/node_t.h"
#include "basic_types/lru_cache.h"
#include "terrain_trees/prt_tree.h"
#include "statistics/statistics.h"
#include "utilities/container_utilities.h"
#include "terrain_trees/mesh_updater.h"
#include "utilities/usage.h"
#include "utilities/cli_parameters.h"
#include "queries/topological_queries.h"
#include "simplification_aux_structures.h"

/// Contraction_Simplifier class ///
class Contraction_Simplifier
{
public:
    Contraction_Simplifier() {}
    /// this procedure simplify the simplicial complex without considering any weight for the edges
    void simplify(PRT_Tree &tree, Mesh &mesh, cli_parameters &cli);
    void simplify_parallel(PRT_Tree &tree, Mesh &mesh, cli_parameters &cli);
    //Generate conflict leaf list and vertex in leaf list
    void preprocess(PRT_Tree &tree, Mesh &mesh, cli_parameters &cli);
    void error_range(PRT_Tree &tree, Mesh &mesh, cli_parameters &cli, itype num_bin);
    
protected:
    
    void update_mesh_and_tree(PRT_Tree &tree, Mesh &mesh, contraction_parameters &params);
    void simplify_compute(Node_V &n, Mesh &mesh, LRU_Cache<int, leaf_VT> &cache, Spatial_Subdivision &division, contraction_parameters &params, PRT_Tree &tree);
    void simplify_compute_parallel(Mesh &mesh,  Spatial_Subdivision &division, contraction_parameters &params, PRT_Tree &tree);
    void simplify_leaf(Node_V &n, Mesh &mesh, LRU_Cache<int, leaf_VT> &cache, contraction_parameters &params, PRT_Tree &tree);
    void simplify_leaf_cross(Node_V &n,int n_id, Mesh &mesh, contraction_parameters &params, PRT_Tree &tree);
    void simplify_leaf_cross_QEM(Node_V &n,int n_id, Mesh &mesh, contraction_parameters &params, PRT_Tree &tree);
    ///edge contraction based on QEM criterion
    void simplify_leaf_QEM(Node_V &n, Mesh &mesh, LRU_Cache<int, leaf_VT> &cache, contraction_parameters &params, PRT_Tree &tree);
    /// skips an edge that has one of the two extreme deleted
    inline bool skip_edge(const ivect &e, Mesh &mesh) /// the procedure checks if one of the two extreme is already deleted
    {
        return (mesh.is_vertex_removed(e[0]) || mesh.is_vertex_removed(e[1]));
    }

    ///
    void contract_edge(ivect &e, ET &et, VT &vt0, VT &vt1, Node_V &outer_v_block, edge_queue &edges,
                       Node_V &n, Mesh &mesh, LRU_Cache<int, leaf_VT> &cache, contraction_parameters &params, map<vector<int>, double>& updated_edges);
    void contract_edge(ivect &e, ET &et, VT &vt0, VT &vt1, Node_V &outer_v_block, edge_queue &edges,
                       Node_V &n, Mesh &mesh, contraction_parameters &params, map<vector<int>, double>& updated_edges);
  // void contract_edge(ivect &e, ET &et, VT &vt0, VT &vt1, Node_V &outer_v_block, edge_queue &edges,
  //                     Node_V &n, Mesh &mesh, contraction_parameters &params);
    /// initializes the VTop and ETop of the edge and its two vertices
    /// plus saves the second leaf block if we are processing a cross-edge
    void get_edge_relations(ivect &e, ET &et, VT *&vt0, VT *&vt1, Node_V *&outer_v_block,
                            Node_V &n, Mesh &mesh, leaf_VT &vts, LRU_Cache<int, leaf_VT> &cache, contraction_parameters &params, PRT_Tree &tree);
    
    void get_edge_relations(ivect &e, ET &et, VT *&vt0, VT *&vt1,bool& v1_is_border, bool& v2_is_border, Node_V *&outer_v_block,
                            Node_V &n, Mesh &mesh, leaf_VT &vts, boost::dynamic_bitset<>is_border_edge, map<int, leaf_VT>&cache, contraction_parameters &params, PRT_Tree &tree);
  
    void get_edge_relations(ivect &e, ET &et, VT *&vt0, VT *&vt1,bool& v1_is_border, bool& v2_is_border, Node_V *&outer_v_block,
                            Node_V &n, Mesh &mesh, leaf_VT &vts,boost::dynamic_bitset<>is_border_edge, LRU_Cache<int, leaf_VT> &cache, contraction_parameters &params,PRT_Tree &tree);

    void update_cached_VT(int v_id, LRU_Cache<int, leaf_VT> &cache);
    /// the VTop is always without removed top-simplices
    VT *get_VT(int v_id, Node_V &n, Mesh &mesh, leaf_VT &vts, LRU_Cache<int, leaf_VT> &cache,
               PRT_Tree &tree, Node_V *&v_block, contraction_parameters &params);

    VT * get_VT(int v_id, Node_V &n, Mesh &mesh, leaf_VT &vts, map<int, leaf_VT>&cache,
               PRT_Tree &tree, Node_V *&v_block, contraction_parameters &params);


    leaf_VT &get_VTS(Node_V &n, Mesh &mesh, LRU_Cache<int, leaf_VT> &cache,
                     PRT_Tree &tree, contraction_parameters &params,boost::dynamic_bitset<>& is_v_border);

    // Find two adjacent triangles of edge e.
    void get_ET(ivect &e, ET &et, Node_V &n, Mesh &mesh, leaf_VT &vts);
    ///

    /// the procedure removes from the VT relation the deleted top d-simplices
    void clean_coboundary(VT &cob, Mesh &mesh);

    void find_candidate_edges_QEM(Node_V &n, Mesh &mesh, leaf_VT &vts, edge_queue &edges, contraction_parameters &params);
    void find_candidate_edges(Node_V &n, Mesh &mesh, leaf_VT &vts, edge_queue &edges, contraction_parameters &params);

    void find_candidate_edges_parallel(Node_V &n, Mesh &mesh, leaf_VT &vts, edge_queue &edges, contraction_parameters &params, bool is_cross);
    //void find_candidate_edges_parallel_QEM(Node_V &n, Mesh &mesh, leaf_VT &vts, edge_queue &edges, contraction_parameters &params, bool is_cross);
    void error_range_leaf(Node_V &n, Mesh &mesh, dvect& edge_costs, coord_type& min, coord_type& max);
  
    /// the procedure updatess
    /// (1) the VTop relation of the surviving vertex
    /// (2) the top simplices for which we change the v2 with v1
    /// (3) checks and updates the top list in the target leaf block n (if a top simplex is now incident in it)
    /// (4) and add the new edges to the edge queue
    void update(const ivect &e, VT &vt, VT &difference, Node_V &n, Node_V &v_block, edge_queue &edges,
                Mesh &mesh, contraction_parameters &params,map<vector<int>, double>& updated_edges);
    void remove_from_mesh(int to_delete_v, ET &et, Mesh &mesh, contraction_parameters &params);
    bool link_condition(int v0, int v1, VT &vt0, VT &vt1, ET &et, Mesh &mesh);
   // bool link_condition(int v0, int v1, VT &vt0, VT &vt1, ET &et,Node_V &n,VV& vv_locks, Mesh &mesh);
    bool link_condition(int v0, int v1, VT &vt0, VT &vt1, ET &et,Node_V &n, Node_V &v_block, VV& vv_locks, Mesh &mesh);

    void update_parallel(const ivect &e, VT &vt, VT &difference, Node_V &n, Node_V &v_block, edge_queue &edges,
                Mesh &mesh, contraction_parameters &params,map<vector<int>, double>& updated_edges);
    void compute_initial_QEM(Mesh &mesh, vector<dvect> &planes);
    void compute_plane_and_QEM(Node_V &n, Mesh &mesh,Spatial_Subdivision &division, PRT_Tree &tree);
    void QEM_leaf(Node_V &n, Mesh &mesh, PRT_Tree &tree);
    void compute_initial_QEM_parallel(PRT_Tree &tree, Mesh &mesh, vector<dvect> &planes);
    void compute_initial_plane_and_QEM_parallel(PRT_Tree &tree, Mesh &mesh);
    void compute_triangle_plane(Mesh &mesh, vector<dvect> &trPl);
    // new_vertex is the pos of the remaining vertex in the edge to be contracted. 0 refers to e[0], 1 refers to e[1]
    double compute_error(int v1, int v2, Mesh &mesh, int &new_vertex_pos);

    bool valid_boundary_condition(int v1,int v2,VT &vt1, VT &vt2,ET& et, bool v1_is_border, bool v2_is_border, Mesh &mesh);
    bool not_fold_over(int v1,int v2,VT &vt1, VT &vt2,ET& et, Mesh &mesh);
    void compute_leaf_QEM(Node_V &n, Mesh &mesh, leaf_VT &vts );
    void update_QEM(ivect& surviving_vertices, Mesh &mesh);
    inline double vertex_error(Matrix q, double x, double y, double z)
    {

        double result = q[0] * x * x;
        result += 2 * q[1] * x * y;
        result += 2 * q[2] * x * z;
        result += 2 * q[3] * x;
        result += q[5] * y * y;
        result += 2 * q[6] * y * z;
        result += 2 * q[7] * y;
        result += q[10] * z * z;
        result += 2 * q[11] * z;
        result += q[15];

        //   cout<<"calculated result: "<<result<<endl;
        //result = round(result * 100000) / 100000.0;
        return result;
    }

    inline double median(vector<double>::const_iterator begin, vector<double>::const_iterator end) {
    int len = end - begin;
    auto it = begin + len / 2;
    double m = *it;
    if ((len % 2) == 0) m = (m + *(--it)) / 2;
    return m;
}
    // update conflict_leafs based on the vertices in the VV(v0) & VV(v1) that are not contained by n or outer_v_block
    void update_conflict_nodes(VV & vv_locks,int n_id, PRT_Tree &tree);

    void generate_v_in_leaf(PRT_Tree &tree, Mesh &mesh, cli_parameters &cli);

    vector<Matrix> initialQuadric;
    vector<dvect> trianglePlane;

    vector<coord_type> contracted_costs;
    vector<coord_type> skipped_costs;

    vector<omp_lock_t> t_locks;
    vector<omp_lock_t> v_locks;
    vector<omp_lock_t> l_locks;
    lists_leafs conflict_leafs;
    vector<int> v_in_leaf;
    //int count_round=0;
};

#endif // CONTRACTION_SIMPLIFIER_H
