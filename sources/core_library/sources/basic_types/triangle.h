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

#ifndef _TRIANGLE_H
#define	_TRIANGLE_H

using namespace std;

#include "basic_wrappers.h"
#include "edge.h"
#include "utilities/sorting_structure.h"

#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include <boost/dynamic_bitset.hpp>

///A class representing a triangle of the mesh
class Triangle {
public:
    ///A constructor method
    Triangle() {}
    /**
     * @brief A copy-constructor method
     * @param orig
     */
    Triangle(const Triangle& orig) { this->vertices = orig.vertices; }
    ///
    Triangle(ivect &v) { this->vertices = v; }
    ///A constructor method
    /*!
     * \param v1 a itype argument, represents the first triangle vertex
     * \param v2 a itype argument, represents the second triangle vertex
     * \param v3 a itype argument, represents the third triangle vertex
     */
    Triangle(itype v1, itype v2, itype v3) { this->set(v1,v2,v3); }
    ///A destructor method
    virtual ~Triangle() {}
    /**
     * @brief A public method that sets the current triangle
     *
     * \param v1 a itype argument, represents the first triangle vertex
     * \param v2 a itype argument, represents the second triangle vertex
     * \param v3 a itype argument, represents the third triangle vertex
     */
    inline void set(itype v1, itype v2, itype v3) { this->vertices = { v1, v2, v3 }; }
    ///A public method that returns the vertex at the pos-th position in the boundary array
    /*!
     * we return absolute (positive) index value stored in the triangle.
     * NOTA: a vertex index is make negative if the opposite edge is a an edge on the mesh borders
     *       (see Border_Checker class for additional informations)
     *
     * \param pos an integer argument, represents the vertex position into the list
     * \return an integer value, representing the position index of the vertex
     */
    inline itype TV(int pos) const {  return abs(this->vertices[pos]); }
    /**
     * @brief A public procedure that updates the index of a vertex in the boundary array
     * @param pos an integer argument, represents the vertex position into the list
     * @param newId an integer representing the new vertex in the boundary
     */
    inline void setTV(int pos, itype newId) { this->vertices[pos] = newId; }

    inline void setTV_keep_border(int pos, itype newId) { this->vertices[pos] = (this->vertices[pos]>0)? newId:-newId; }
    /**
     * @brief A public procedure that returns an edge in the boundary of the triangle
     * NOTA: as a vertex index can be negative we return the positive indices during an edge initialization
     *
     * @param pos an integer representing the edge position in the boundary
     * @param e an integer vector that it is set with the sorted edge extrema
     */
    void TE(int pos, ivect& e)
    {        
        e = { abs(this->TV((pos+1)%this->vertices_num())), abs(this->TV((pos+2)%this->vertices_num())) };
        sort(e.begin(),e.end());
    }
    /**
     * @brief A public procedure that returns an edge in the boundary of the triangle
     * NOTA: as a vertex index can be negative we return the positive indices during an edge initialization
     *
     * @param pos an integer representing the edge position in the boundary
     * @param e an Edge& variable that it is set with the sorted edge extrema
     */
    inline void TE(int pos, Edge &e) { e.init_edge(abs(this->TV((pos+1)%this->vertices_num())),abs(this->TV((pos+2)%this->vertices_num()))); }
    /**
     * @brief A public procedure that returns the tuple composed by the two indices forming the edge and the triangle id
     *
     * @param pos an integer representing the edge position in the boundary
     * @param f a triangle_triangle_tuple& that represents the tuple
     * @param t_id an integer representing the triangle index
     */
    inline void edge_tuple(int pos, edge_triangle_tuple &f, itype t_id) const { f.sort_and_set(this->TV((pos+1)%3),this->TV((pos+2)%3),t_id,pos); }
    /**
     * @brief A public method that returns the position index of a vertex (if exists) in the local boundary array of the triangle
     *
     * @param v_id an integer representing the position index of the vertex in the mesh
     * @return itype an integer containing the position index of the vertex in the boundary array or -1 if the vertex is not in the boundary of the cell
     */
    inline int vertex_index(itype v_id)
    {
        for(int j=0; j<vertices_num(); j++)
        {
            if(this->TV(j) == v_id)
                return j;
        }
        return -1;
    }
    /**
     * @brief A public method that checks if a vertex is on the boundary of the triangle
     *
     * @param v_id an itype representing the position index of the vertex in the mesh
     * @return true if the vertex is on the boundary, false otherwise
     */
    inline bool has_vertex(itype v_id) { return (this->vertex_index(v_id)!=-1); }
    /**
     * @brief A public method that checks if an edge is on the mesh borders
     *
     * NOTA: a vertex index is make negative if the opposite edge is a an edge on the mesh borders
     *       (see Border_Checker class for additional informations)
     *
     * @param pos an integer representing the edge position in the boundary
     * @return true if the edge is on the mesh borders, false otherwise
     */
    inline bool is_border_edge(int pos) { return (this->vertices[pos] < 0); }

    /// Smaller-index vertex getter
    inline itype minindex() const { return min(min(abs(vertices[0]),abs(vertices[1])),abs(vertices[2])); }
    /// Bigger-index vertex getter
    inline itype maxindex() const { return max(max(abs(vertices[0]),abs(vertices[1])),abs(vertices[2])); }

    inline bool is_removed() {return this->vertices[0]==0;}
    /**
     * @brief operator ==
     * @param p
     * @param q
     * @return
     */
    inline friend bool operator== (const Triangle &p, const Triangle &q)
    {
        bool b[3];
        b[0] = false; b[1] = false; b[2] = false;
        for(int i=0;i<p.vertices_num();i++)
        {
            for(int j=0;j<q.vertices_num();j++)
            {
                if(!b[j] && p.TV(i)==q.TV(j))
                {
                    b[j] = true;
                    break;
                }
            }
        }
        return b[0] && b[1] && b[2];
    }
    /**
     * @brief operator !=
     * @param p
     * @param q
     * @return
     */
    inline friend bool operator!= (const Triangle &p, const Triangle &q) { return !(p==q); }
    /**
     * @brief operator <<
     * @param out
     * @param p
     * @return
     */
    inline friend std::ostream& operator<<(std::ostream& out, const Triangle& p)
    {
        out <<"T[" << p.vertices[0] << " " << p.vertices[1] << " " << p.vertices[2] << "]";
        return out;
    }
    /**
     * @brief A public procedure returning the number of vertices in the triangle (constant value)
     * @return an integer
     */
    inline int vertices_num() const { return vertices.size(); }

    /**
     * @brief A public procedure checking if an edge, represented as an array, is in the boundary of the triangle
     * @param e, a vector variable representing the edge
     * @return true if e is on the boundary of the triangle, false otherwise
     */
    inline bool has_edge(const ivect &e)
    {
        int count = 0;
        for(int i=0;i<this->vertices_num();i++)
        {
            if(e[0]==abs(vertices[i]) || e[1]==abs(vertices[i]))
                count++;
            if(count==2)
                return true;
        }
        return false;
    }

    /**
     * @brief A public procedure checking if a simplex (i.e., a vertex, an edge or a triangle), represented as an array, is in the boundary of the triangle
     * @param s, a vector variable representing the simplex
     * @return true if s is on the boundary of the triangle, false otherwise
     */
    inline bool has_simplex(ivect &s)
    {
        boost::dynamic_bitset<> b(s.size());

        for(unsigned i=0;i<s.size();i++)
        {
            for(int j=0;j<this->vertices_num();j++)
            {
                if(s[i]==abs(vertices[j]))
                {
                    b[i] = 1;
                    break;
                }
            }
        }
        return b.count() == b.size();
    }

    /**
     * @brief A public method that returns the position index of an edge (if exists) in the local boundary array of the triangle
     *
     * @param e is a vector representing the edge
     * @return short an integer containing the position index of the edge in the boundary array or -1 if the edge is not in the boundary of the cell
     */
    inline short edge_index(const ivect &e)
    {
        for(int i=0; i<this->vertices_num(); i++)
        {
            if(e[0] != abs(vertices[i]) && e[1] != abs(vertices[i]))
                return i;
        }
        return -1;
    }

    /**
     * @brief A public procedure converting the triangle to a positive array of vertices indices
     * @param t, a vector variable containing the three vertices indices
     */
    inline void convert_to_vec(ivect &t) const
    {
        t.clear();
        for(int i=0;i<vertices_num();i++)
            t.push_back(abs(vertices[i]));
    }

private:
    ///A private variable representing the array of vertices in the boundary
    ivect vertices;
};

#endif	/* _TRIANGLE_H */

