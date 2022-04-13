#ifndef GRAD_AWARE_CONTRACTION_H
#define GRAD_AWARE_CONTRACTION_H

#include "contraction_simplifier.h"
#include "morse/forman_gradient.h"
#include "morse/forman_gradient_topological_relations.h"

class Gradient_Aware_Simplifier: public Contraction_Simplifier{


public:
    Gradient_Aware_Simplifier(){}
    void gradient_aware_simplify(PRT_Tree &tree, Mesh &mesh, cli_parameters &cli,Forman_Gradient &gradient);
    void gradient_aware_simplify_parallel(PRT_Tree &tree, Mesh &mesh, cli_parameters &cli,Forman_Gradient &gradient);
protected:
    void simplify_compute(Node_V &n,  Mesh &mesh, LRU_Cache<int, leaf_VT> &cache,Spatial_Subdivision &division,  contraction_parameters &params,PRT_Tree &tree,Forman_Gradient &gradient);
    void simplify_compute_parallel( Mesh &mesh, Spatial_Subdivision &division,  contraction_parameters &params,PRT_Tree &tree,Forman_Gradient &gradient);
    void simplify_leaf(Node_V &n, Mesh &mesh, LRU_Cache<int, leaf_VT> &cache, contraction_parameters &params, PRT_Tree& tree,Forman_Gradient &gradient); 
    void simplify_leaf_QEM(Node_V &n, Mesh &mesh, LRU_Cache<int, leaf_VT> &cache, contraction_parameters &params, PRT_Tree& tree,Forman_Gradient &gradient);
    void simplify_leaf_cross(Node_V &n,int n_id, Mesh &mesh, contraction_parameters &params, PRT_Tree &tree,Forman_Gradient &gradient); 
    void simplify_leaf_cross_QEM(Node_V &n,int n_id, Mesh &mesh, contraction_parameters &params,PRT_Tree &tree, Forman_Gradient &gradient);
   
    void contract_edge(ivect &e, ET &et, VT &vt0, VT &vt1, Node_V &outer_v_block, edge_queue &edges,
                              Node_V &n, Mesh &mesh, LRU_Cache<int, leaf_VT> &cache, contraction_parameters &params,Forman_Gradient &gradient, map<vector<int>, double>& updated_edges);

    void contract_edge(ivect &e, ET &et, VT &vt0, VT &vt1, Node_V &outer_v_block, edge_queue &edges,
                              Node_V &n, Mesh &mesh,  contraction_parameters &params,Forman_Gradient &gradient, map<vector<int>, double>& updated_edges);
    // void get_edge_relations(ivect &e, ET &et, VT *&vt0, VT *&vt1, Node_V *& outer_v_block, 
    //                                Node_V &n, Mesh &mesh, leaf_VT &vts, LRU_Cache<int,leaf_VT> &cache, contraction_parameters &params, PRT_Tree &tree);
    
    // valid_gradient_configuration should always be the last condition to be checked for a candidate edge
    // because in this function, if an edge satisfies the gradient condition, the gradient field will be modified. 
    // if later another condition delines the edge contraction, then the modification on the gradient field will be wrong. 
 
    bool valid_gradient_configuration(int v1,int v2,VT &vt1, VT &vt2,ET& et, bool v1_is_border, bool v2_is_border, Forman_Gradient &gradient, Mesh &mesh);
    // void update(const ivect &e, VT& vt, VT& difference, Node_V &n, Node_V &v_block, edge_queue &edges,
    //                                       Mesh &mesh, contraction_parameters &params);

    //void get_edge_relations(ivect &e, ET &et, VT *&vt0, VT *&vt1,bool& v1_is_border, bool& v2_is_border, Node_V *&outer_v_block,
    //                                            Node_V &n, Mesh &mesh, leaf_VT &vts,boost::dynamic_bitset<>is_border_edge, LRU_Cache<int, leaf_VT> &cache, contraction_parameters &params,PRT_Tree &tree);
    void update_mesh_and_tree(PRT_Tree &tree, Mesh &mesh, contraction_parameters &params, Forman_Gradient &gradient,  cli_parameters &cli);

    
};

#endif // GRAD_AWARE_CONTRACTION_H