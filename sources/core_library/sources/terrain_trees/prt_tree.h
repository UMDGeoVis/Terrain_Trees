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

#ifndef P_TREE_H
#define	P_TREE_H

#include "tree.h"
#include "node_v.h"
#include <geometry/geometry_wrapper.h>

///An inner-class, implementing Tree, that represents a tree which uses the PR-T criterion
class PRT_Tree : public Tree<Node_V>
{
public:
    ///A constructor method
    PRT_Tree() { }
    ///A constructor method
    /*!
     * \param vertices_per_leaf an integer argument, represents the maximum number of vertices admitted for a tree node
     */
    PRT_Tree(int vertices_per_leaf, int sons_num);
    ///A constructor method
    PRT_Tree(const PRT_Tree& orig) : Tree<Node_V>(orig) { this->vertices_threshold = orig.vertices_threshold; }
    ///A destructor method
    virtual ~PRT_Tree() {}
    ///A public method that sets the maximum number of vertices admitted for a node
    /*!
     * \param maxV an integer argument, represents the maximum number of vertices
     */
    inline void set_vertices_threshold(int maxV) { this->vertices_threshold = maxV; }
    ///A public method that returns the maximum number of vertices admitted for a node
    /*!
     * \return an integer value representing the maximum number of vertices admitted
     */
    inline int get_vertices_threshold() { return this->vertices_threshold; }
    ///A public method that builds the tree
    /*!
     * This method before inserts all the mesh vertices and then all the tetrahedra.
     * The insertion of tetrahedra does not change the hierarchy.
     */
    void build_tree();
    /*!
     * \brief A public method that builds the tree from a soup of triangles
     * \param soup
     */
    void build_tree(Soup &soup);
    ///A public method that builds the tree from a points cloud
    void build_tree_from_cloud(vertex_multifield &multifield);


    void compact_vertices_lists(Node_V &n, Mesh &mesh, ivect &surviving_vertices);
    void reinsert_triangles();

    void update_tree(Node_V &n, ivect &new_v_positions,ivect &new_t_positions, bool all_deleted, itype& index_counter);
    void update_vertex_index(Node_V &n, ivect &new_v_positions, itype& index_counter);
    void get_leaf_indexing_vertex(Node_V &n, int v_id, Node_V *&res);


    void init_leaves_list(Node_V &n);
        /**
     * @brief A public method that returns the leaf block at i-th position in the leaves array
     *
     * @param i an integer argument representing the position index of the leaf
     * @return Node_V* a pointer reference to the leaf block
     */
    inline Node_V* get_leaf(unsigned i) { return this->leaves[i]; }


    /**
     * @brief A public method that returns the number of leaf blocks in the tree
     *
     * @return unsigned
     */
    inline unsigned get_leaves_number() { return this->leaves.size(); }

    /**
     * @brief A public method that recursively visits the Terrain tree and unify the tree branches that index less than kv vertices.
     *
     * The procedure returns the number of vertices indexed in a tree branch.
     *
     * @param n a Node_V& argument representing the current node
     * @param mesh a Mesh& argument representing the mesh indexed by the Terrain tree
     * @return int an integer containing the number of vertices indexed in a tree branch
     */
    int visit_and_unify(Node_V &n, Mesh &mesh);

private:
    ///A private variable representing the maximum number of vertices admitted for a node
    int vertices_threshold;
    ///A private method that adds a vertex to the tree structure
    /*!
     * This method checks if a vertex is contained by a node, and then insert the vertex into
     * the node list if the node is a leaf, otherwise if the node is internal a recursive call is activated.
     * If the maximum number of vertices for a node is exceeded, a splitting operation is activated.
     *
     * \param n a Node_V& argument, represents the node in which we try to insert the vertex
     * \param domain a Box& argument, represents the node domain
     * \param level an integer argument representing the level of n in the hierarchy
     * \param v an integer argument, represents the vertex to insert
     */
    void add_vertex(Node_V& n, Box& domain, int level, itype v);
    ///A private method that adds a triangle to the tree structure
    /*!
     * This method checks if a triangle has a proper intersection with the node, and then insert the triangle into
     * the node list if the node is a leaf, otherwise if the node is internal a recursive call is activated.
     *
     * \param n a Node_V& argument, represents the node in which we try to insert the triangle
     * \param domain a Box& argument, represents the node domain
     * \param level an integer argument representing the level of n in the hierarchy
     * \param t an integer argument, represents the triangle to insert
     */
    void add_triangle(Node_V& n, Box& domain, int level, itype t);
    ///A protected method that split a node, creating the sons node, following the current division type
    /*!
     * This method also reinsert all the vertices and the tetrahedra, which are in the splitted node, into the sons node
     * \param n a Node_V& argument, represents the node to split
     * \param domain a Box& argument, represents the node domain
     * \param level an integer argument representing the level of n in the hierarchy
     */
    void split(Node_V& n, Box& domain, int level);
    ///A public method that checks if a node has raised the maximum block capacity
    /*!
     * \param n a Node_V& argument, represents the node to check
     * \return a boolean value, true if the limit is exceeded, false otherwise
     */
    inline bool is_full(Node_V &n) { return (n.get_v_array_size() > this->vertices_threshold); }

    /**
     * @brief A private method that insert a vertex from a cells soup in the leaf block that contains it
     *
     * @param n a Node_V& argument representing the current node
     * @param domain a Box& argument, represents the node domain
     * @param level an integer argument representing the level of n in the hierarchy
     * @param v a Vertex& argument representing the vertex to insert
     * @param vertex_index an integer argument that after the insertion of v saves the position index of v in the indexed mesh representation
     * @param indexed_triangle a vector of integers representing the indexed representation of the cell containing the vertex v
     *
     * This method checks if a vertex is contained by a node, and then insert the vertex into
     * the vertices array if the node is a leaf, otherwise if the node is internal a recursive call is activated.
     * If the maximum number of vertices for a node is exceeded, a splitting operation is activated.
     */
    void add_vertex_from_soup(Node_V& n, Box& domain, int level, Vertex& v, itype vertex_index, ivect &indexed_triangle, bool &first_time);
    /**
     * @brief A private method that checks if a vertex is already into the vertices array of leaf n.
     *
     * The method returns the position index of the vertex if a vertex is already into the vertices array, -1 otherwise
     *
     * @param n a Node_V& argument representing the current node
     * @param v a Vertex& argument representing the target vertex
     * @return int the position index of the vertex if v is already into the array, -1 otherwise
     */
    itype is_already_inserted(Node_V& n, Vertex &v);

    ///A private method that adds a vertex from a points cloud to the tree structure
    /*!
     * This method checks if a vertex is contained by a node, and then insert the vertex into
     * the node list if the node is a leaf only if it is not already present, otherwise if the node is internal a recursive call is activated.
     * If the maximum number of vertices for a node is exceeded, a splitting operation is activated.
     *
     * If the vertex is already indexed by the corresponding leaf block, we flag the vertex as removed and we save its field value
     * to the other vertex indexed by the leaf.
     *
     * \param n a Node_V& argument, represents the node in which we try to insert the vertex
     * \param domain a Box& argument, represents the node domain
     * \param level an integer argument representing the level of n in the hierarchy
     * \param v a Vertex& argument representing the vertex to insert
     * \param vertex_index an itype representing the position index of v in the vertex array
     * \param multifield a map representing the multifield assiciated to a vertex
    */
    void add_vertex_from_cloud(Node_V& n, Box& domain, int level, Vertex& v, itype vertex_index, vertex_multifield &multifield);

    vector<Node_V*> leaves;
};

#endif	/* P_TREE_H */

