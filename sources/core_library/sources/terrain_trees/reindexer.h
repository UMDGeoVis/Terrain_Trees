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

#ifndef REINDEXER_H
#define REINDEXER_H

#include <vector>
#include <set>
#include <iostream>
#include <map>

#include "basic_types/mesh.h"
#include "pmrt_tree.h"
#include "node_v.h"

/**
 * @brief A class that exploits the spatial coherence of the indexed enties of a trianglesl tree and resort accordingly these in the indexed mesh and trianglesl tree representation
 * This class also compresses the tree representation using the Compressed encoding.
 */
class Reindexer
{
public:
    Reindexer() { this->indices_counter=1; }
    /**
     * @brief A public method that exploit the spatial coherence of vertices and triangles and compresses their representation within the tree
     *
     * @param tree a T& argument representing the tree to reindex and compress
     *
     * NOTA: the other three parameters have been added for backward compatibility
     */
    void reindex_tree_and_mesh(PMRT_Tree& tree, bool save_v_indices, ivect &original_vertex_indices, bool save_t_indices, ivect &original_triangle_indices);
    //for PMR-T tree
    /**
     * @brief A public method that exploit the spatial coherence of vertices and triangles and compresses their representation within the tree
     *
     * @param tree a P_Tree& argument representing the tree to reindex and compress
     */
    template<class T> void reindex_tree_and_mesh(T& tree, bool save_v_indices, ivect &original_vertex_indices, bool save_t_indices, ivect &original_triangle_indices);
    // for PR-T and PM-T trees

    template<class T> void reindex_triangle_array(T& tree, bool save_t_indices, ivect &original_triangle_indices);

private:
    ///A private vector containing the coherent position indices of vertices/triangles
    ivect coherent_indices;
    ///A private variable representing the counter of position indices
    int indices_counter;
    /// A private variable that encodes for each triangles type an associative array, in which the key represents a leaves tuple
    /// and the value is the number of triangles indexed by the leaves tuple and the starting position index for the triangles group.
    map<ivect,pair<utype,utype> > leaf_tuples_array;

    // FOR VERTICES
    /**
     * @brief A private method that exploits the spatial coherence of the vertices, compresses their representation
     * in the tree and resort the corresponding array of the mesh
     * The procedure is compatible for PM-T trees and PMR-T trees
     *
     * @param n a Node_V& argument, represents the node
     * @param domain a Box& argument, represents the node domain
     * @param level an integer argument representing the level of n in the hierarchy
     * @param decomposition a Spatial_Decomposition& argument, representing the tree subdivision
     * @param mesh a Mesh& argument, the triangle mesh
     */
    void reindex_vertices(Node_T& n, Box& domain, int level, Spatial_Subdivision& decomposition, Mesh& mesh,
                          bool save_v_indices, ivect &original_vertex_indices);
    /**
     * @brief A private method that exploits the spatial coherence of the vertices, compresses their representation in the tree and resort the
     * corresponding array of the mesh
     * The procedure is compatible for PR-T trees and PM-T trees
     *
     * @param n a Node_V& argument, represents the node
     * @param decomposition a Spatial_Decomposition& argument, representing the tree subdivision
     */
    void reindex_vertices(Node_V& n, Spatial_Subdivision& decomposition, bool save_v_indices, ivect &original_vertex_indices);
    /**
     * @brief A private method that resort the vertices array of the mesh
     *
     * @param mesh a Mesh& argument, the triangle mesh
     */
    void update_mesh_vertices(Mesh& mesh);
    /**
     * @brief A private method that updates the Triangle-Vertex relation for each triangle of the mesh
     *
     * @param mesh a Mesh& argument, the triangle mesh
     */
    void update_triangles_boundaries(Mesh &mesh);

    // FOR TRIANGLES
    /**
     * @brief A private method that compresses the triangles array of a leaf block
     *
     * @param n a N& argument, represents the node
     * @param decomposition a Spatial_Decomposition& argument, representing the tree subdivision
     * @param mesh a Mesh& argument, the triangle mesh
     */
    template<class N> void compress_t_array(N& n,ivect &new_t_list);
    /**
     * @brief A private method that resort the triangles array of the mesh
     *
     * @param mesh a Mesh& argument, the triangle mesh
     */
    void update_mesh_triangles(Mesh& mesh);
    /**
     * @brief A private method that resets the auxiliary variables of the class
     *
     */
    inline void reset()
    {
        this->coherent_indices.clear();
        this->indices_counter=1;
    }

    /**
     * @brief A private method that extract for each triangle indexed in the tree the leaves tuple indexing it (visiting procedure).
     * The procedure visits recursively the nodes of the tree, and for each leaf block calls extract_leaf_top_association_leaf
     *
     * @param n a N& argument representing the current node
     * @param mesh a Mesh& argument representing the indexed mesh to reindex
     * @param root a N& argument representing the root of the tree
     */
    template<class N> void extract_leaf_tri_association(N& n, Mesh &mesh, N &root);
    void extract_leaf_tri_association(Node_T& n, Box& dom, int level, Spatial_Subdivision& decomposition, Mesh &mesh, Node_T &root);
    /**
     * @brief A private method that extract for each triangle indexed in the tree the leaves tuple indexing it.
     *
     * @param n a N& argument representing the current node
     * @param mesh a Mesh& argument representing the indexed mesh to reindex
     * @param root a N& argument representing the root of the tree
     */
    template<class N> void extract_leaf_tri_association_leaf(N& n, Mesh &mesh, N& root);
    void extract_leaf_tri_association_leaf(Node_T& n, Box& dom, Spatial_Subdivision& decomposition, Mesh &mesh, Node_T& root);
    /**
     * @brief A private method that visit the tree in order to extract the leaves tuple for a triangle indexed in more than a leaf block
     * The procedure visits recursively the nodes of the tree.
     *
     * @param n a N& argument representing the current node
     * @param v_ids an integer vector containing the vertices of the triangle indexed outside the target node
     * @param key an integer vector encoding the leaves tuple
     */
    template<class N> void get_tri_tuple_key(N& n, ivect &v_ids, ivect &key);
    void get_tri_tuple_key(Node_T& n, Box& dom, int level, Spatial_Subdivision& decomposition, Mesh &mesh, ivect &v_ids, ivect &key);
    /**
     * @brief A private method that extracts the spatially coherent position indexes for the triangles of the mesh
     *
     * @param save_original_positions a boolean that enables the backup of the triangles original position indexes in a auxiliary array
     */
    void get_triangles_reordered_indexes(bool save_t_indices, ivect &original_triangle_indices);
    /**
     * @brief A private method that compress the triangle arrays in the tree (visiting procedure).
     *
     * The procedure visits recursively the nodes of the tree and call for each leaf block the compress_tree_representation_top_leaf procedure.
     *
     * @param n a N& argument representing the current node
     */
    template<class N> void compress_tree_representation_top(N& n);
    /**
     * @brief A private method that compress the triangle arrays in the tree.
     *
     * The procedure extracts the new list with the updated triangles position indexes and then compress it using the Sequential-range encoding (SRE)
     * (internally calls compress_t_list).
     *
     * @param n a N& argument representing the current node
     */
    template<class N> void compress_tree_representation_top_leaf(N& n);

};

template<class T> void Reindexer::reindex_tree_and_mesh(T &tree, bool save_v_indices, ivect &original_vertex_indices, bool save_t_indices, ivect &original_triangle_indices)
{
    coherent_indices.assign(tree.get_mesh().get_vertices_num(),-1);
    reindex_vertices(tree.get_root(),tree.get_subdivision(),save_v_indices,original_vertex_indices);
    update_triangles_boundaries(tree.get_mesh());
    update_mesh_vertices(tree.get_mesh());
    reset();

    coherent_indices.assign(tree.get_mesh().get_triangles_num(),-1);
    extract_leaf_tri_association(tree.get_root(),tree.get_mesh(),tree.get_root());
    get_triangles_reordered_indexes(save_t_indices,original_triangle_indices);
    compress_tree_representation_top(tree.get_root());
    update_mesh_triangles(tree.get_mesh());

    reset();
    return;
}

template<class T> void Reindexer::reindex_triangle_array(T &tree, bool save_t_indices, ivect &original_triangle_indices)
{
    coherent_indices.assign(tree.get_mesh().get_triangles_num(),-1);

    extract_leaf_tri_association(tree.get_root(),tree.get_mesh(),tree.get_root());

    get_triangles_reordered_indexes(save_t_indices,original_triangle_indices);

    compress_tree_representation_top(tree.get_root());

    update_mesh_triangles(tree.get_mesh());

    reset();
    return;
}



template<class N> void Reindexer::compress_t_array(N& n, ivect &new_t_list)
{
    sort(new_t_list.begin(),new_t_list.end());

    int count=0;
    itype start_t_id = new_t_list[0];

    if(new_t_list.size()==1)
    {
       n.add_triangle(start_t_id);
       return;
    }
    //otherwise visit the t_list

    //now obtain the new encoding
    for(unsigned i=0; i<new_t_list.size(); i++)
    {
        if((i+1<new_t_list.size()) && (new_t_list[i]+1) == new_t_list[i+1]) //I have a consecutive range in the t_list
        {
            count++;
        }
        else //found a possible new range
        {
            if(count > 1) //more than two consecutive triangles
            {
                //the range should be from start_t_id to start_t_id + count
                //example: 12 - 13 - 14 - 15
                //is encoded: -12 - 3
                //the first index is implicit
                n.add_triangle(-start_t_id);
                n.add_triangle(count);
            }
            else //less or equal to two
            {
                n.add_triangle(start_t_id);
                if(count==1)
                    n.add_triangle(start_t_id+count);
            }
            //re-init
            count = 0;
            start_t_id = new_t_list[i+1];
        }
    }
}

///////////////////////***************************/////////////////////////////
template<class N> void Reindexer::extract_leaf_tri_association(N &n, Mesh &mesh, N &root)
{
    if (n.is_leaf())
    {        
        extract_leaf_tri_association_leaf(n,mesh,root);
    }
    else
    {
        for(typename N::child_iterator it=n.begin(); it!=n.end(); ++it)
        {
            if(*it != NULL)
                extract_leaf_tri_association(**it,mesh,root);
        }
    }
}

template<class N>  void Reindexer::extract_leaf_tri_association_leaf(N &n, Mesh &mesh, N &root)
{
    /// if there are no vertices in the leaf we have nothing to do..
    if(!n.indexes_vertices())
        return;
    ivect internal_t_key; internal_t_key.push_back(n.get_v_start());

    for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        Triangle &t = mesh.get_triangle(*t_id);
        itype min_v = t.minindex();
        if(n.indexes_vertex(min_v))
        {
            ivect key = internal_t_key;
            pair<utype,utype> value = make_pair(this->leaf_tuples_array.size(),1);

            if(!n.completely_indexes_triangle_vertices(t))
            {
                ivect outer_v_ids;
                for(int i=0; i<t.vertices_num(); i++)
                {
                    if(!n.indexes_vertex(t.TV(i)))
                        outer_v_ids.push_back(t.TV(i));
                }
                get_tri_tuple_key(root,outer_v_ids,key);
            }

            pair<map<ivect,pair<utype,utype> >::iterator,bool> ret = this->leaf_tuples_array.insert(make_pair(key,value));
            if(ret.second) // inserted
            {
                this->coherent_indices[*t_id-1] = ret.first->second.first;
            }
            else
            {
                this->coherent_indices[*t_id-1] = ret.first->second.first;
                ret.first->second.second++;
            }
        }
    }
}

template<class N>  void Reindexer::compress_tree_representation_top(N &n)
{
    if (n.is_leaf())
    {
        compress_tree_representation_top_leaf(n);
    }
    else
    {
        for(typename N::child_iterator it=n.begin(); it!=n.end(); ++it)
        {
            if(*it != NULL)
                compress_tree_representation_top(**it);
        }
    }
}

template<class N>  void Reindexer::compress_tree_representation_top_leaf(N& n)
{
    ivect new_t_list;
    for(RunIteratorPair itPair = n.make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        new_t_list.push_back(coherent_indices[*t_id-1]);
    }
    n.clear_t_array();

    if(new_t_list.size()>0)
    {
        compress_t_array(n,new_t_list);
    }
}

template<class N>  void Reindexer::get_tri_tuple_key(N &n, ivect &v_ids, ivect &key)
{
    if (n.is_leaf())
    {
        // when I'm entering a leaf the v of a top are for sure indexed..
        key.push_back(n.get_v_start());
        // clean-up phase for v_ids array
        for(ivect_iter it=v_ids.begin(); it!=v_ids.end();)
        {
            if(n.indexes_vertex(*it))
            {
                v_ids.erase(it);
            }
            else
                ++it;
        }
    }
    else
    {
        for(typename N::child_iterator it=n.begin(); it!=n.end(); ++it)
        {
            if(*it != NULL)
            {
                for(ivect_iter itv=v_ids.begin(); itv!=v_ids.end(); ++itv)
                {
                    if((*it)->indexes_vertex(*itv))
                    {
                        get_tri_tuple_key(**it,v_ids,key);
                        break;
                    }
                }
            }
        }
    }
}

#endif // REINDEXER_H
