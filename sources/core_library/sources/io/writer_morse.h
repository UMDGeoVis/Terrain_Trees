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

#ifndef _WRITER_MORSE_H
#define	_WRITER_MORSE_H



#include "io/writer.h"
#include "morse/forman_gradient_aux_structure.h"
#include "morse/ig.h"

using namespace std;

///A class that provides an interface for writing to file or standard output some data structures or statistics
class Writer_Morse : public Writer {
public:
    static void write_asc1cells_VTK(string mesh_name, string operation_type, itype vertices_per_leaf, simplices_multimap &triangles, Mesh &mesh,
                                     ivect &original_triangle_indices, ivect &original_vertex_indices, dvect &original_vertex_fields, bool revert_to_original_field);
    static void write_desc2cells_VTK(string mesh_name, string operation_type, itype vertices_per_leaf, ivect &segmentation, Mesh &mesh,
                                     ivect &original_triangle_indices, ivect &original_vertex_indices, dvect &original_vertex_fields, bool revert_to_original_field);
    static void write_desc1cells_VTK(string mesh_name, string operation_type, itype vertices_per_leaf, simplices_multimap &edges, Mesh &mesh,
                                     ivect &original_vertex_indices, dvect &original_vertex_fields, bool revert_to_original_field);
    static void write_asc2cells_VTK(string mesh_name, string operation_type, itype vertices_per_leaf, ivect &segmentation, Mesh &mesh,
                                    ivect &original_vertex_indices, dvect &original_vertex_fields, bool revert_to_original_field);
    static void write_incidence_graph_VTK(string mesh_name, string operation_type, itype vertices_per_leaf, IG &forman_ig, Mesh &mesh,
                                          ivect &original_vertex_indices, dvect &original_vertex_fields, bool revert_to_original_field);
    static void write_critical_clusters(string mesh_name, forman_aux_structures::critical_clusters &cc, Mesh &mesh);

    static void write_asc1cells_OBJ(string mesh_name, string operation_type, itype vertices_per_leaf, simplices_multimap &triangles, Mesh &mesh,
                                     ivect &original_triangle_indices, ivect &original_vertex_indices, dvect &original_vertex_fields, bool revert_to_original_field);
    static void write_asc1cells_CSV(string mesh_name, string operation_type, itype vertices_per_leaf, simplices_multimap &triangles, Mesh &mesh,
                                     ivect &original_triangle_indices, ivect &original_vertex_indices, dvect &original_vertex_fields, bool revert_to_original_field);
    static void write_asc1cells_vertices_CSV(string mesh_name, string operation_type, itype vertices_per_leaf, simplices_multimap &triangles, Mesh &mesh,
                                     ivect &original_triangle_indices, ivect &original_vertex_indices, dvect &original_vertex_fields, bool revert_to_original_field);
private:
    ///A constructor method
    Writer_Morse() {}
    ///A constructor method
    Writer_Morse(const Writer_Morse&) {}
    ///A destructor method
    virtual ~Writer_Morse() {}
};


#endif	/* _WRITER_MORSE_H */

