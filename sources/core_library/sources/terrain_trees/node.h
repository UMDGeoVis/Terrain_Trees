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

#ifndef _NODE_H
#define	_NODE_H

#include <vector>
#include <cstddef>
#include <set>
#include "basic_types/basic_wrappers.h"
#include "basic_types/box.h"
#include "basic_types/mesh.h"
#include "run_iterator.h"

/**
 * @brief A super-class, not instantiable, that represents a generic node of the tree with associated an array of triangles
 */
template<class N> class Node
{
public:
    ///A destructor method
    virtual ~Node() {}
    ///A public method that sets a node son
    /*!
     * \param n a N* argument, represents the son to set
     * \param pos an integer, represents the son position in the array
     */
    inline void set_son(N* n, int pos) { (*this->sons)[pos] = n; }
    ///A public method that checks if the node is a leaf node
    /*!
     * \return a boolean, true if the node is a leaf, false otherwise
     */
    inline bool is_leaf() const { return (this->sons==NULL); }
    ///A public method that returns a node son
    /*!
     * \param i an integer argument, represents the son position into the list
     * \return a Node*, representing the son at the i-th position
     */
    inline N* get_son(int i) { return (*this->sons)[i]; }
    /**
     * @brief A public method that initializes the sons array
     * @param son_number an integer containing the sons number
     */
    inline void init_sons(int son_number) { this->sons = new vector<N*>(); this->sons->assign(son_number,NULL); }
    /**
     * @brief A public method that clear the son array
     *
     */
    inline void delete_sons() { this->sons->clear(); delete this->sons; this->sons=NULL; }

    ///A public method that adds a triangle index to the array
    /*!
     * \param ind an itype argument, representing the triangle index
     */
    inline void add_triangle(itype ind) { this->triangles.push_back(ind); }

    ///A public method that returns the run_iterator pair to navigate the triangles array
    inline RunIteratorPair make_t_array_iterator_pair() { return run_iterator<itype>::make_run_iterator_pair(triangles); }
    ///A public method that returns the begin run_iterator to navigate the triangles array
    inline RunIterator t_array_begin_iterator() { return run_iterator<itype>(triangles.begin(),triangles.end()); }
    ///A public method that returns the end run_iterator to navigate the triangles array
    inline RunIterator t_array_end_iterator() { return run_iterator<itype>(triangles.end()); }

    /**
     * @brief A public method that returns the number of indexed triangles
     * NOTA: this method expand the runs and returns the real number of top d-cells indexed in the node
     *
     * @return int
     */
    inline int get_real_t_array_size() const { return run_iterator<itype>(triangles.begin(),triangles.end()).elementCountFast(triangles); }
    /**
     * @brief A public method that returns the size of the trianglesl array
     *
     * @return utype
     */
    inline utype get_t_array_size() { return this->triangles.size(); }
    /**
     * @brief A public method returning the triangles array
     *
     * @return
     */
    inline ivect get_t_array() const { return this->triangles; }
    /**
     * @brief A public method that clears the space used by the triangles array
     */
    inline void clear_t_array() { triangles.clear(); }
    ///A public method that return the begin iterator of the triangles array for explicitly unroll the runs of triangles
    inline ivect_iter get_t_array_begin() { return this->triangles.begin(); }
    ///A public method that return the end iterator of the triangles array for explicitly unroll the runs of triangles
    inline ivect_iter get_t_array_end() { return this->triangles.end(); }

    inline void set_t_array(ivect new_t_list){this->triangles=new_t_list;}
    // geometric procedures //
    /**
     * @brief A public method that checks if all the four vertices of a triangle are indexed by the node
     *
     * @param t a Triangle& argument
     * @param domain a Box& representing the node domain
     * @param mesh a Mesh& representing the trianglesl mesh
     * @return true if all the four vertices are indexed, false otherwise
     */
    inline bool completely_indexes_triangle_vertices_dom(Triangle &t, Box& domain, Mesh& mesh)
    {
        for(int i=0; i<t.vertices_num(); i++)
            if(!domain.contains(mesh.get_vertex(t.TV(i)),mesh.get_domain().get_max()))
                return false;
        return true;
    }
    /**
     * @brief A public method that checks if at least one vertex of a triangle is indexed by the node
     *
     * @param t a Triangle& argument
     * @param domain a Box& representing the node domain
     * @param mesh a Mesh& representing the trianglesl mesh
     * @return true if at least one vertex is indexed, false otherwise
     */
    inline bool indexes_triangle_vertices_dom(Triangle &t, Box& domain, Mesh& mesh)
    {
        for(int i=0; i<t.vertices_num(); i++)
            if(domain.contains(mesh.get_vertex(t.TV(i)),mesh.get_domain().get_max()))
                return true;
        return false;
    }
    /**
     * @brief A public method that computes the bounding box of a run
     * @param id an iterator to the current array entry
     * @param bb a Box& argument, that is set with the run bounding box (if a run is found)
     * @param mesh a Mesh& representing the trianglesl mesh
     * @param run a pair that will contains the run, if any
     * @return true if a run has been encounter, false otherwise
     */
    bool get_run_bounding_box(ivect_iter &id, Box& bb, Mesh &mesh, pair<itype,itype> &run);

    /**
     * @brief A typename alias for accessing the iterator to the children nodes vector
     *
     */
    typedef typename std::vector<N*>::iterator child_iterator;
    /**
     * @brief A public method that returns the begin iterator to the son array
     *
     * @return child_iterator an iterator to the first position of the array
     */
    inline child_iterator begin() { return sons->begin(); }
    /**
     * @brief A public method that returns the end iterator to the son array
     *
     * @return child_iterator an iterator to the last position of the array
     */
    inline child_iterator end() { return sons->end(); }
    /**
     * @brief set_parent
     * @param p
     */
    inline void set_parent(N* p) { this->parent = p; }
    /**
     * @brief get_parent
     * @return
     */
    inline N* get_parent() { return this->parent; }
    


protected:    
    ///A constructor method
    Node() { this->sons = NULL; this->parent = NULL; }
    ///A copy-constructor method
    Node(const Node& orig)
    {
        this->sons = orig.sons;
        this->triangles = orig.triangles;
    }
    ///A private variable representing the list of node sons
    vector<N*> *sons;
    ///A private variable representing the list containing the triangles indexed by the node
    ivect triangles;
    ///A private variable representing the parent node
    N* parent;
};

template<class N> bool Node<N>::get_run_bounding_box(ivect_iter &id, Box& bb, Mesh &mesh, pair<itype, itype> &run)
{
    if(*id<0) //I have a run
    {
        run.first = abs(*id);
        ++id;
        run.second = run.first + *id;

        coord_type min_p[2]={0,0}, max_p[2]={0,0};

        // init phase
        itype t_id=run.first;
        Triangle &t_first = mesh.get_triangle(t_id);
        min_p[0] = mesh.get_vertex(t_first.TV(0)).get_x();
        min_p[1] = mesh.get_vertex(t_first.TV(0)).get_y();
        max_p[0] = mesh.get_vertex(t_first.TV(0)).get_x();
        max_p[1] = mesh.get_vertex(t_first.TV(0)).get_y();

        for(int i=1; i<t_first.vertices_num(); i++)
        {
            Vertex &v = mesh.get_vertex(t_first.TV(i));
            for(int j=0;j<v.get_dimension();j++)
            {
                if(v.get_c(j) < min_p[j])
                    min_p[j] = v.get_c(j);
                else if(v.get_c(j) > max_p[j])
                    max_p[j] = v.get_c(j);
            }
        }

        t_id++;
        for(; t_id<=run.second; t_id++)
        {
            Triangle &tri = mesh.get_triangle(t_id);

            for(int i=0; i<tri.vertices_num(); i++)
            {
                Vertex &v = mesh.get_vertex(tri.TV(i));
                for(int j=0;j<v.get_dimension();j++)
                {
                    if(v.get_c(j) < min_p[j])
                        min_p[j] = v.get_c(j);
                    else if(v.get_c(j) > max_p[j])
                        max_p[j] = v.get_c(j);
                }
            }
        }
        //save the computed bounding box
        bb.set_min(min_p[0],min_p[1]);
        bb.set_max(max_p[0],max_p[1]);
        return true;
    }
    else
    {
        return false;
    }
}

#endif	/* _NODE_H */

