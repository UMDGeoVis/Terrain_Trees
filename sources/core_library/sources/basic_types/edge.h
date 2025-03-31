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

#ifndef EDGE_H
#define EDGE_H

#include "basic_wrappers.h"
#include <sstream>
#include <cmath>
#include <vector>

/**
 * @brief The Edge class defines a representation for an edge in the triangle mesh
 * NOTA: an edge is represented with a pair of positive integers
 */
using namespace std;
class Edge
{
public:
    /// A constructor procedure
    Edge()
    {
        this->vertices = {0,0};
    }
    /// A constructor procedure
    Edge(itype v1, itype v2)
    {
        this->vertices = {min(v1,v2),max(v1,v2)};
    }
    /// A constructor procedure
    Edge(const Edge &e)
    {
        this->vertices = e.vertices;
    }
    /**
     * @brief A public procedure initializing an edge with two integer values
     * NOTA: the two integers are sorted
     * @param v1
     * @param v2
     */
    inline void init_edge(itype v1, itype v2)
    {
        this->vertices = {min(v1,v2),max(v1,v2)};
    }
    /**
     * @brief A public procedure returning the extreme index at a given position
     * @param pos an integer representing the extreme position in the array
     * @return an integer
     */
    inline itype EV(itype pos) const { return this->vertices[pos]; }
    /**
     * @brief A public procedure returning the number of vertices in the edge (constant value)
     * @return an integer
     */
    inline int vertices_num() const { return vertices.size(); }
//    /**
//     * @brief A public procedure returning the opposite extreme w.r.t. to a given one
//     * @param v
//     * @return
//     */
//    inline itype get_other_extreme(itype v)
//    {
//        if(v == vertices[0])
//            return vertices[1];
//        else if(v == vertices[1])
//            return vertices[0];
//        else
//        {
//            cerr<<"[get_other_extreme] wrong call"<<endl;
//            return -1;
//        }
//    }
    /**
     * @brief A public method that checks if a vertex is on the boundary of the edge
     *
     * @param v_id an itype representing the position index of the vertex in the mesh
     * @return true if the vertex is on the boundary, false otherwise
     */
    inline bool has_vertex(itype v) { return (v == vertices[0] || v == vertices[1]); }

    /// Smaller-index endpoint getter
    inline itype minindex() const {return vertices[0];}

    /// Bigger-index endpoint getter
    inline itype maxindex() const {return vertices[1];}

//    /**
//     * @brief A public procedure converting the edge to a positive array of vertices indices
//     * @param t, a vector variable containing the three vertices indices
//     */
//    inline void convert_to_vec(ivect &e) const
//    {
//        for(int i=0;i<vertices_num();i++)
//            e.push_back(abs(vertices[i]));
//    }

//    inline std::pair<itype,itype> make_edge_pair() const
//    {
//        return make_pair(abs(vertices[0]),abs(vertices[1]));
//    }
    /**
     * @brief operator ==
     * @param p
     * @param q
     * @return
     */
    inline friend bool operator== (const Edge &p, const Edge &q)
    {
        if(p.vertices[0]!=q.vertices[0]) return false;
        if(p.vertices[1]!=q.vertices[1]) return false;
        return true;
    }
    /**
     * @brief operator !=
     * @param p
     * @param q
     * @return
     */
    inline friend bool operator!= (const Edge &p, const Edge &q) { return !(p==q); }
    /**
     * @brief operator <
     * @param s
     * @return
     */
    inline bool operator<(const Edge& s) const
    {
        return ((this->vertices[0] < s.vertices[0]) || (this->vertices[0] == s.vertices[0] && this->vertices[1] < s.vertices[1]));
    }
    /**
     * @brief operator >
     * @param s
     * @return
     */
    inline bool operator>(const Edge& s) const
    {
        return ((this->vertices[0] > s.vertices[0]) || (this->vertices[0] == s.vertices[0] && this->vertices[1] > s.vertices[1]));
    }
    /**
     * @brief operator <<
     * @param out
     * @param q
     * @return
     */
    inline friend std::ostream& operator<< (std::ostream &out, Edge &q)
    {
        out << "["<<q.vertices[0]<< " " << q.vertices[1] << "]";
        return out;
    }

private:
    /// A private variable containing the vertex indicies in the boundary of the edge
    ivect vertices;
};

#endif // EDGE_H
