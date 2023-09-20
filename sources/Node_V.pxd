from libcpp.vector cimport vector
from Mesh cimport Mesh
from libcpp cimport bool
cdef extern from "terrain_trees/node_v.h":
    cdef cppclass Node_V:
        Node_V() except +
        int get_v_start()
        int get_v_end()
        void get_VT(vector[vector[int]] &all_vt, Mesh &mesh)
        bool is_leaf()
        Node_V* get_son(int)
        bool indexes_vertices()

         