/*
    This file is part of the Terrain Trees library.

    Author(s): Riccardo Fellegara (riccardo.fellegara@gmail.com)

    This project has been supported by the Italian Ministry of Education and
    Research under the PRIN 2009 program, and by the National Science Foundation
    under grant number IIS-1116747.

    The triangle Trees library is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    The triangle Trees library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with the Terrain Trees library.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef CRITICAL_POINTS_EXTRACTOR_H
#define CRITICAL_POINTS_EXTRACTOR_H

#include "terrain_trees/node_v.h"
#include "terrain_trees/node_t.h"
#include "terrain_trees/spatial_subdivision.h"

// the key corresponds to the field value of the flat area
// the value correspond to the vertices and triangles belonging to the flat area
typedef map<coord_type,iset> flat_areas;

enum class Point_Type: short {REGULAR=0, MINIMUM=1, SADDLE=2, MULTIPLE_SADDLE=3, MAXIMUM=4};

class Critical_Points_Extractor
{
public:
    //
    Critical_Points_Extractor() { }
    //
    void compute_critical_points(Node_V &n, Mesh &mesh, Spatial_Subdivision &division);
    void compute_critical_points(Node_T &n, Box &dom, Mesh &mesh, Spatial_Subdivision &division);
    //
    void print_stats();
    inline vector<Point_Type>& get_critical_points() { return this->critical_points; }

private:
    vector<Point_Type> critical_points;

    void extract_local_critical_points(Node_V &n, Mesh &mesh, Spatial_Subdivision &division, flat_areas &fa); ///buggy
    void extract_local_critical_points(Node_T &n, Box &dom, int level, Mesh &mesh, Spatial_Subdivision &division, flat_areas &fa); ///buggy

    void extract_local_critical_points(Node_V &n, Mesh &mesh, Spatial_Subdivision &division);
    void extract_local_critical_points(Node_T &n, Box &dom, int level, Mesh &mesh, Spatial_Subdivision &division);

    void analyze_leaf(leaf_VV &vvs, leaf_VT &vts, itype v_start, Mesh &mesh, flat_areas &fa); ///buggy
    void analyze_leaf(leaf_VV &vvs, leaf_VT &vts, itype v_start, Mesh &mesh);

    void extract_critical_points_from_flat_areas(flat_areas &fa, Mesh &mesh); ///buggy

    void init_adjacent_vertices(itype v_id, Vertex &v, VT &vt, Mesh &mesh, map<itype, iset> &adj_upper, map<itype, iset> &adj_lower);
    utype get_components_num(map<itype, iset> &adj_map, map<itype, itype> &v_flag);

};

#endif // CRITICAL_POINTS_EXTRACTOR_H
