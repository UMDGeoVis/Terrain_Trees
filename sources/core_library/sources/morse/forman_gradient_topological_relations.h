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

#ifndef FORMAN_GRADIENT_TOPOLOGICAL_RELATIONS_H
#define FORMAN_GRADIENT_TOPOLOGICAL_RELATIONS_H

#include "forman_gradient.h"
#include "dangling_paths_handler.h"

class Forman_Gradient_Topological_Relations
{
public:
    static void set_VTstar(Node_V &n, itype v, itype new_vt_star, leaf_VTstar &all_vtstar, vtstar_cache &cache, Node_V &root, Spatial_Subdivision &division);

    //return a VT relation of one of f vertices, and the position of that vertex in f
    static iset& get_VV(Node_V &n, itype v, leaf_VV &all_vv, asc3_cache &cache,
                        Node_V &root, Spatial_Subdivision &division, Mesh &mesh);
    static itype get_VTstar(Node_V &n, itype v, leaf_VTstar &all_vtstar, vtstar_cache &cache,
                   Node_V &root, Spatial_Subdivision &division, Mesh &mesh, Forman_Gradient &gradient);
    static itype get_VTstar(Node_T &n, itype v_start, itype v_end, itype v, leaf_VTstar &all_vtstar, vtstar_cache &cache,
                   Node_T &root, Spatial_Subdivision &division, Mesh &mesh, Forman_Gradient &gradient);
    static pair<itype,itype> get_ET(Node_V &n, ivect &e, leaf_ET &local_ef, et_cache &cache,
                                    Node_V &root, Spatial_Subdivision &division, Mesh &mesh);
    static pair<itype,itype> get_ET(Node_T &n, itype v_start, itype v_end, ivect &e, leaf_ET &local_ef, et_cache &cache,
                                    Node_T &root, Spatial_Subdivision &division, Mesh &mesh);

    static void get_VTstar_ET(local_VTstar_ET &all_rels, Node_V &n, Mesh& mesh, Forman_Gradient &gradient);
    static void get_VTstar_ET(local_VTstar_ET &all_rels, Node_T &n, itype v_start, itype v_end, Mesh& mesh, Forman_Gradient &gradient);
    static void get_VTstar_ETstar(desc1rels &all_rels, Node_V &n, Mesh& mesh, Forman_Gradient &gradient);

    static void get_VTstar_VV(asc2rels &all_rels, Node_V &n, Mesh& mesh, Forman_Gradient &gradient);

    static itype get_triangle_id(ivect max_tri, VT &vt, Mesh &mesh);
    static void get_VTstar_ET(local_VTstar_ET &all_rels, Node_V &n, Mesh& mesh, Forman_Gradient &gradient, mig_cache& cache);

private:
    Forman_Gradient_Topological_Relations() {}

    static void get_VTstar(leaf_VTstar &vtstars, Node_V &n, Mesh& mesh, Forman_Gradient &gradient);
    static void get_VTstar(leaf_VTstar &vtstars, Node_T &n, itype v_start, itype v_end, Mesh& mesh, Forman_Gradient &gradient);

    static bool check_VTstar(itype v_id, itype t_id, Triangle& t, leaf_VTstar &vtstars, itype v_start, Mesh &mesh, Forman_Gradient &gradient);
    static bool check_ETstar(ivect &e, itype t_id, Triangle &t, leaf_ETstar &etstars, Mesh &mesh, Forman_Gradient &gradient);

};

#endif // FORMAN_GRADIENT_TOPOLOGICAL_RELATIONS_H
