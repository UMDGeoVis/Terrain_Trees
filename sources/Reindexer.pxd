from PRT_Tree cimport PRT_Tree
from libcpp.vector cimport vector
from libcpp cimport bool
cdef extern from "terrain_trees/reindexer.h":
    cdef cppclass Reindexer:
        Reindexer() except +
        void reindex_tree_and_mesh(PRT_Tree& tree, bool save_v_indices, vector[int] &original_vertex_indices, bool save_t_indices, vector[int] &original_triangle_indices)