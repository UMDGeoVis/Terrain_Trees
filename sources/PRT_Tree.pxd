from libcpp cimport bool
from Mesh cimport Mesh
from Node_V cimport Node_V
cdef extern from "terrain_trees/prt_tree.h":
    cdef cppclass PRT_Tree:
        PRT_Tree(int, int) except +
        void build_tree()
        Mesh& get_mesh()
        Node_V& get_root()
    