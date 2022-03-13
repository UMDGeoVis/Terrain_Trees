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

#ifndef NODE_V_H
#define NODE_V_H

#include "terrain_trees/node.h"

/**
 * @brief The Node_V class extends the Node class, encoding also a list of vertices.
 * This class is used by the PR-T and PM-T trees
 */
class Node_V : public Node<Node_V>
{
public:
    ///A constructor
    Node_V() : Node<Node_V>() { }
    ///A copy-constructor
    Node_V(const Node_V& orig) : Node<Node_V>(orig) { this->vertices = orig.vertices; }
    ///A destructor
    virtual ~Node_V() {}
    ///A public method that add a vertex index to the corresponding node list
    /*!
     * \param ind an itype argument, representing the vertex index
     */
    inline void add_vertex(itype ind) { this->vertices.push_back(ind); }
    ///A public method that free space used by the two lists
    inline void clear_v_array() { this->vertices.clear(); }
    inline void set_v_array(ivect new_v_list){this->vertices=new_v_list;}
    /**
     * @brief A public method that set the vertices range after exploiting the vertices spatial coherence
     *
     * @param start an itype with the first vertex position index of the node
     * @param end an itype with the first vertex position index outside the node
     */
    inline void set_v_range(itype start, itype end) { vertices.push_back(-start); vertices.push_back(end-start-1); }
    ///
    inline itype get_v_start() const { return abs(vertices[0]); }
    ///
    inline itype get_v_end() const { return (abs(vertices[0])+vertices[1])+1; }    
    /**
     * @brief A public method that checks if a vertex is indexed by the node
     * NOTA: to use only after the index is built/loaded from file
     *
     * @param v_id an itype representing the position index of the vertex
     * @return true if v_id is indexed, false otherwise
     */
    inline bool indexes_vertex(itype v_id)
    {
        return (indexes_vertices() && v_id >= get_v_start() && v_id < get_v_end());
    }

    /**
     * @brief A public method that checks if a triangle is indexed by the node
     * The method checks if at least one of the vertices of the triangle is indexed by the node
     * NOTA: only range comparisons are executed
     *
     * @param t a Triangle& representing the triangle
     * @return true if at least one vertex is indexed, false otherwise
     */
    inline bool indexes_triangle_vertices(Triangle &t)
    {
        if(!this->indexes_vertices())
            return false;
        for(int i=0; i<t.vertices_num(); i++)
            if(this->indexes_vertex(t.TV(i)))
                return true;
        return false;
    }
    /**
     * @brief A public method that checks if a triangle is completely indexed by a node
     * The method checks if all the vertices of the triangle are indexed by the node
     * NOTA: only range comparisons are executed
     *
     * @param t a Triangle& representing the triangle
     * @return true if all the vertices are indexed, false otherwise
     */
    inline bool completely_indexes_triangle_vertices(Triangle &t)
    {
        if(!this->indexes_vertices())
            return false;
        for(int i=0; i<t.vertices_num(); i++)
            if(!this->indexes_vertex(t.TV(i)))
                return false;
        return true;
    }

    inline bool partial_indexes_triangle_vertices(Triangle &t)
    {
        if(!this->indexes_vertices())
            return false;
        // We can use the counter here instead.    
        bool first=this->indexes_vertex(t.TV(0));
        for(int i=1; i<t.vertices_num(); i++)
            if(this->indexes_vertex(t.TV(i))!=first)
                return true;
        return false;
    }


    // inline bool completely_indexes_edge_vertices(Triangle &t)
    // {
    //     if(!this->indexes_vertices())
    //         return false;
    //     for(int i=0; i<t.vertices_num(); i++)
    //         if(!this->indexes_vertex(t.TV(i)))
    //             return false;
    //     return true;
    // }

   inline bool completely_indexes_simplex(const ivect &s)
   {
       if(!this->indexes_vertices())
           return false;
       for(ivect_citer it=s.begin(); it!=s.end(); ++it)
           if(!this->indexes_vertex(*it))
               return false;
       return true;
   }

    /**
     * @brief operator <<
     * @param out
     * @param p
     * @return
     */
    friend std::ostream& operator<<(std::ostream& out, const Node_V& p)
    {
        if(p.is_leaf())
        {
            if(p.is_v_array_run_compressed()) //it has been reindexed
               { out <<"Leaf["<< p.get_v_start() << " " << p.get_v_end()<<"]";
                 out << "Leaf[v->"<< p.get_real_v_array_size() << " t->"<< p.get_real_t_array_size() <<"]";
               }
            else   //should be reindexed
            {
                 
                //   out <<"Leaf["<< p.get_v_start() << " " << p.get_v_end()<<"]";
                out << "Leaf[v->"<< p.get_real_v_array_size() << " t->"<< p.get_real_t_array_size() <<"]";
   
            }
        }
        else
            out <<"Node["<< p.get_v_start() << " " << p.get_v_end()<<"]";
        return out;
    }
    /**
     * @brief A public method that returns a pair of iterators to the vertices array
     *
     * @return RunIteratorPair
     */
    inline RunIteratorPair make_v_array_iterator_pair() { return run_iterator<itype>::make_run_iterator_pair(vertices); }
    /**
     * @brief A public method that returns the iterator to the begin of the vertices array
     *
     * @return RunIterator
     */
    inline RunIterator v_array_begin_iterator() { return run_iterator<itype>(vertices.begin(),vertices.end()); }
    /**
     * @brief A public method that returns the iterator to the end of the vertices array
     *
     * @return RunIterator
     */
    inline RunIterator v_array_end_iterator() { return run_iterator<itype>(vertices.end()); }
    /**
     * @brief A public method that returns the number of indexed vertices
     *
     * @return utype
     */
    inline utype get_real_v_array_size() const { return run_iterator<itype>(vertices.begin(),vertices.end()).elementCountFast(vertices); }
    /**
     * @brief A public method that returns the size of the vertices array
     *
     * @return itype
     */
    inline utype get_v_array_size() const { return this->vertices.size(); }

    /**
     * @brief A public method that checks if a simplex has all the vertices indexed in already processed leaf blocks
     *
     * @param s a templated variable representing the simplex to process
     *
     * @return true if the condition is verified, false otherwise
     */
    template<class S> inline bool has_all_vertices_visited(S &s) { return (s.maxindex() < this->get_v_end()); }
    /**
     * @brief A public method that checks if a vertex is indexed in already processed leaf blocks
     *
     * @param v an integer variable representing the vertex position index
     *
     * @return true if the condition is verified, false otherwise
     */
    inline bool visited_vertex(itype v) { return (v < this->get_v_end()); }
    /**
     * @brief A public method that checks if the current leaf block indexes vertices
     *
     * @return true if the condition is verified, false otherwise
     */
    inline bool indexes_vertices() { return (this->get_v_array_size() > 0); }

    /**
     * @brief A public method that extracts the Vertex-Triangle (VT) relations for the vertices indexed in the current block
     * @param all_vt a leaf_VT variable, that encodes the VT relations
     * @param mesh a Mesh& variable representing the triangle mesh
     */
    void get_VT(leaf_VT &all_vt, Mesh &mesh);

    /**
     * @brief A public method that extracts the Vertex-Triangle (VT) relations for the vertices indexed in the current block
     * @param all_vt a leaf_VT variable, that encodes the VT relations
     * @param mesh a Mesh& variable representing the triangle mesh
     */
    void get_VT_and_border(leaf_VT &all_vt, boost::dynamic_bitset<> &is_v_border,Mesh &mesh);

    /**
     * @brief A public method that extracts the Vertex-Vertex (VV) relations for the vertices indexed in the current block
     * @param all_vv a leaf_VV variable, that encodes the VV relations
     * @param mesh a Mesh& variable representing the triangle mesh
     */
    void get_VV(leaf_VV &all_vv, Mesh& mesh);
    
    void get_VV_vector(leaf_VV_vec &all_vv, Mesh& mesh);
    
    /**
     * @brief A public method that extracts the Vertex-Triangle (VT) and Vertex-Vertex (VV) relations for the vertices indexed in the current block
     * @param all_vv a leaf_VV variable, that encodes the VV relations
     * @param all_vt a leaf_VT variable, that encodes the VT relations
     * @param mesh a Mesh& variable representing the triangle mesh
     */
    void get_VV_VT(leaf_VV &all_vv, leaf_VT &all_vt, Mesh& mesh);
    /**
     * @brief A public method that extracts the Edge-Triangle (ET) relations for the edges indexed in the current block
     * Notice that an edge is considered processed only if the current leaf block indexes the extreme with the higher position index
     * in this way each edge is processed once during each traversal of the tree
     *
     * @param ets a leaf_ET variable, that encodes the ET relations
     * @param mesh a Mesh& variable representing the triangle mesh
     */
    void get_ET(leaf_ET &ets,  Mesh &mesh);
    /**
     * @brief A public method that removes the vertices that have been deleted during a simplification process
     * The procedure saves in an auxiliary vector the position indexes of those not deleted
     *
     * @param mesh is the mesh indexed by the tree
     * @param surviving_vertices, a vector of integer containing the position indexes of the vertices not deleted
     */
    void compact_vertices_array(Mesh &mesh, ivect &surviving_vertices);
     /**
     * @brief A public method that updates the position indexes of the vertices
     *
     * @param new_v_indices, a vector of integers containing the updated position indexes
     *
     * NOTA: the SRE encoding is not valid any more. A reindexing operation is needed after calling this procedure on all the leaf blocks
     */
  
    void update_vertex_indices(ivect &new_v_indices, itype &index_counter);
    /**
     * @brief A public method that updates the top cells list of the leaf block after a simplification process.

     * @param mesh a Mesh& argument representing the mesh indexed by the Stellar tree
     */
    inline void compact_triangle_array(Mesh &mesh);// Change compact_triangle_arrays into compact_triangle_array



 /**
     * @brief A public method that updates the top cells encoded in the leaf block by updating the position indexes
     *
     * @param new_top_positions, a vector of vectors of integers that contains the new position indexes of the top cells
     * @param all_deleted, a bit-vector stating if a specific top d-cell type has been completely erased by a simplification procedure
     *
     * NOTA: this update removes the run compression!
     */
    void update_and_compress_triangles_arrays(ivect  &new_t_positions, bool all_deleted);
    /**
     * @brief A public method that updates the top d-cells encoded in the leaf block by updating the position indexes
     *
     * @param d an integer representing the list to update
     * @param new_indices, a vector of integers that contains the new position indexes of the top d-cells
     *
     * NOTA: this update removes the run compression!
     */
    void update_and_compress_triangles_array(ivect &new_indices);
    /**
     * @brief A public procedure that compress a top d-cells array of the leaf block
     *
     * The procedure compresses the top d-cells array using the Sequential-Range Encoding (SRE).
     *
     * @param pos an integer identifying the position index of the array to update/compress
     * @param new_t_list an integer vector containing the updated array
     */
    void compress_triangle_array(ivect &new_t_list);

    /////NOTE: this functions is supposed to be used for the simplification only. 

    inline bool index_triangle(itype tid, Mesh &mesh, itype v_id=-1){
        if(v_id !=-1) {
            /// then extract the other (up to) two vertex ids in the boundary of tid
            /// if one of such ids is in this (the leaf)
            /// then return true
            itype v1=-1,v2=-1;
            Triangle t=mesh.get_triangle(tid);
            itype v_pos=t.vertex_index(v_id);
            v1=t.TV((v_pos+1)%3);
            v2=t.TV((v_pos+2)%3);
            if(this->indexes_vertex(v1)||this->indexes_vertex(v2))
            return true;
    
        }
        ///otherwise iterate over the triangle list
        for(RunIteratorPair itPair = make_t_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        if(*t_id==tid){
            return true;
        }
    }
        return false;
    }

    /**
     * @brief A public method that check if the vertices list is encoded using the run compression encoding
     *
     * @return true if the run compression is used, false otherwise
     */
    inline bool is_v_array_run_compressed() const { return (this->vertices.size()==2 && this->vertices[0]<0); }


////ONLY FOR DEBUG, SHOULD REMOVE BEFORE MERGING
// inline void print_vertices(){
//     cout<<"V list:";
//     for(auto it = vertices.begin();it!=vertices.end();it++)
//     cout<<*it<<", ";
//     cout<<endl;
// }

protected:    
   ivect vertices;
};

#endif // NODE_V_H
