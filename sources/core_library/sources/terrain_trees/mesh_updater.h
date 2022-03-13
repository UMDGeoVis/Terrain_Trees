#ifndef MESH_UPDATER_H
#define MESH_UPDATER_H

#include "basic_types/mesh.h"

/**
 * @brief The Mesh_Updater class represents an interface for updating a simplicial (or CP) mesh.
 *
 */
class Mesh_Updater
{
public:
    Mesh_Updater() {}

    /**
     * @brief A public method that compresses the vertices array of the mesh
     * The procedure removes the vertices flagged as deleted during a simplification procedure
     * and returns the new position indexes of the vertices not deleted during the simplification.
     *
     * @param mesh a Mesh& representing the indexed mesh
     * @param new_v_positions a vector of integers that at the end of the procedure contains the new position indexes of the vertices
     * @param surviving_vertices a vector of integers containing the non-deleted vertices
     *
     * NOTA: the size of the new_v_positions array is equal to the size of the original array of the mesh, and the
     * entries corresponding to deleted elements are flagged with -1 value.
     */
    void clean_vertices_array(Mesh &mesh, ivect &new_v_positions, ivect &surviving_vertices);
    /**
     * @brief A public method that compresses the triangles arrays of the mesh
     * The procedure removes the triangles flagged as deleted during a simplification procedure
     * and returns, for each top k-cell array, if it has been completely erased.
     *
     * @param mesh a Mesh& representing the indexed mesh
     * @param new_v_positions a vector of integers containing the updated position indexes of the vertices
     * @param new_top_positions a vector of vectors that, at the end of the procedure, contains the new position indexes of the triangles
     * @param simplification_counters a vector referring to the number of triangles removed during a simplification procedure
     * @return a bit-vector sized as the number of triangles encoded, returning the triangles type that have been completely erased by a simplification procedure
     *
     * NOTA: the size of the new_v_positions and new_top_positions arrays is equal to the size of the original arrays of the mesh, and the
     * entries corresponding to deleted elements are flagged with -1 value.
     *
     * NOTA2: during the procedure the boundary relations of the triangles are also updated 
     */
     bool update_and_clean_triangles_arrays(Mesh &mesh, ivect &new_v_positions,
                                                                                         ivect &new_triangle_positions, itype simplification_counters);

    /**
     * @brief A public method that updates the boundary relations of the triangles, by assigning the new position indexes of the corresponding vertices,
     * and resorting the vertices array in the mesh in order to be coherent to the new vertices ordering (calls reorder_vertices_array)
     *
     * @param mesh a Mesh& representing the indexed mesh
     * @param new_v_positions a vector of integers containing the updated position indexes of the vertices
     */
    void update_triangles_boundary(Mesh &mesh, ivect &new_v_positions);

    /**
     * @brief A private method that update the vertices array of the mesh by assigning the new position indexes.
     * The procedure execute the swaps inline an thus it does not require extra storage.
     * The number of swap operation is exactly |V|, where |V| is the number of vertices of the mesh
     *
     * @param mesh a Mesh& representing the indexed mesh
     * @param new_v_positions a vector of integers that at the end of the procedure contains the new position indexes of the vertices
     */
    void reorder_vertices_array(Mesh &mesh, ivect &new_v_positions);
 

private:
    /**
     * @brief A public method that compresses the top d-cells arrays of the mesh.
     * The procedure removes the top d-cells flagged as deleted during a simplification procedure
     * and returns if the corresponding array has been completely erased.
     * @param mesh a Mesh& representing the indexed mesh
     * @param new_v_positions a vector of integers containing the updated position indexes of the vertices
     * @param new_top_positions a vector that, at the end of the procedure, contains the new position indexes of the top d-cells
     * @param counter an integer representing the number of the top d-cells deleted during the simplification procedure
     * @return a boolean, true if all the top d-cells have been deleted, false otherwise
     */
     bool update_and_clean_triangles_array(Mesh &mesh, ivect &new_v_positions, ivect &new_triangle_positions, int counter);
    /**
     * @brief A public method that compresses the top d-cells array of the mesh.
     * The procedure removes the top d-cell deleted during a simplification procedure
     * and returns the new position indexes of the top d-cells not deleted during the simplification.
     *
     * @param mesh a Mesh& representing the indexed mesh

     * @param new_triangle_positions a vector of integers that at the end of the procedure contains the new position indexes of the top d-cells
     *
     * NOTA: the size of the new_triangle_positions array is equal to the size of the original array of the mesh, and the
     * entries corresponding to deleted elements are flagged with -1 value.
     */
     void clean_triangles_array(Mesh &mesh,  ivect &new_triangle_positions);

    /**
     * @brief A public procedure that coherently reorder all the top d-cells of the mesh
     *
     * @param mesh a Mesh& representing the indexed mesh
  
     * @param new_top_positions a vector that, at the end of the procedure, contains the new position indexes of the top d-cells
     */
     void reorder_triangles_array(Mesh &mesh, ivect new_triangles_positions);
};


#endif // MESH_UPDATER_H
